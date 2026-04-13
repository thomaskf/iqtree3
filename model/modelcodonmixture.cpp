//
//  modelcodonmixture.cpp
//  iqtree
//
//  Created by Minh Bui on 4/3/2025.
//
//  This file implements the codon mixture models M1a, M2a, M3, M7 and M8
//  using the new parameter syntax described in
//  docs/251022_codon_mix_syntax.pdf .  In short:
//
//      M1a:  GY{kappa}+CMIX1a{p_0, omega_0}                  (2 params)
//      M2a:  GY{kappa}+CMIX2a{p_0, omega_0, p_2, omega_2}    (4 params)
//      M3:   GY{kappa}+CMIX3[.K]{p_0,o_0,p_1,o_1,...}        (2K params)
//      M7:   GY{kappa}+CMIX7[.K]{p, q}                       (2 params, beta-dist shape)
//      M8:   GY{kappa}+CMIX8[.K]{p, q, p_0, omega_0}         (4 params)
//
//  Each value can be a number (fixed), a single '?' (optimise from default)
//  or a number followed by '?' (initial value, then optimise).
//

#include "modelcodonmixture.h"

#include "ratebeta.h"

#include <cfloat>
#include <utility>
#include <vector>

namespace {

// ----------------------------------------------------------------------
// ParamSpec helpers
// ----------------------------------------------------------------------
//
// A ParamSpec captures one user-supplied parameter token in the
// {...} list following +CMIXxxx.  Three forms are supported:
//
//      "?"        ->  has_value = false, optimize = true
//      "0.5"      ->  has_value = true,  optimize = false   (fixed)
//      "0.5?"     ->  has_value = true,  optimize = true    (init then optimise)
//
struct ParamSpec {
    bool   has_value = false;
    double value     = 0.0;
    bool   optimize  = false;
};

static string trim(const string &s) {
    size_t a = 0, b = s.size();
    while (a < b && isspace((unsigned char)s[a])) ++a;
    while (b > a && isspace((unsigned char)s[b-1])) --b;
    return s.substr(a, b - a);
}

static ParamSpec parseParamSpec(const string &raw) {
    ParamSpec ps;
    string t = trim(raw);
    if (t.empty())
        outError("Empty parameter inside CMIX{...}");
    if (t == "?") {
        ps.has_value = false;
        ps.optimize  = true;
        return ps;
    }
    if (t.back() == '?') {
        ps.has_value = true;
        ps.optimize  = true;
        ps.value     = convert_double(t.substr(0, t.size() - 1).c_str());
    } else {
        ps.has_value = true;
        ps.optimize  = false;
        ps.value     = convert_double(t.c_str());
    }
    return ps;
}

// Resolve a ParamSpec into an initial numeric value, falling back to a
// caller-supplied default if the user used a bare "?" without a number.
static inline double resolveValue(const ParamSpec &ps, double default_value) {
    return ps.has_value ? ps.value : default_value;
}

} // anonymous namespace


ModelCodonMixture::ModelCodonMixture(string orig_model_name, string model_name,
                                     ModelsBlock *models_block, StateFreqType freq, string freq_params,
                                     PhyloTree *tree, bool optimize_weights)
: ModelMarkov(tree), ModelMixture(tree) {

    if (tree->aln->seq_type != SEQ_CODON)
        outError("Can't apply codon mixture model as sequence type is not codon");

    // ------------------------------------------------------------------
    // 1. Locate "+CMIX..." in the original model name and decode the
    //    sub-type (1a / 2a / 3 / 7 / 8) and the optional ".K" override
    //    for the number of categories.
    // ------------------------------------------------------------------
    auto cmix_pos = orig_model_name.find("+CMIX");
    ASSERT(cmix_pos != string::npos);
    auto end_pos = orig_model_name.find_first_of("+*{", cmix_pos + 1);

    string cmix_type;
    if (end_pos == string::npos)
        cmix_type = orig_model_name.substr(cmix_pos + 5);
    else
        cmix_type = orig_model_name.substr(cmix_pos + 5, end_pos - cmix_pos - 5);

    int ncat = 0;
    {
        auto ncat_pos = cmix_type.find('.');
        if (ncat_pos != string::npos) {
            string ncat_str = cmix_type.substr(ncat_pos + 1);
            ncat = stoi(ncat_str);
            if (ncat < 1)
                outError("Invalid number of categories in " + orig_model_name.substr(cmix_pos + 1));
            cmix_type = cmix_type.substr(0, ncat_pos);
        }
    }

    cmix_subtype = cmix_type;

    alpha = 1.0;
    beta  = 1.0;
    iteration_num = 0;

    // ------------------------------------------------------------------
    // 2. Parse the parameter list inside +CMIXxxx{...}
    // ------------------------------------------------------------------
    vector<ParamSpec> specs;
    bool user_input_param = false;
    if (end_pos != string::npos && orig_model_name[end_pos] == '{') {
        auto close_br_pos = orig_model_name.find_first_of("}", end_pos + 1);
        if (close_br_pos == string::npos)
            outError("Missing closing '}' in CMIX parameter list of " + orig_model_name);
        string cmix_param_list = orig_model_name.substr(end_pos + 1, close_br_pos - end_pos - 1);
        StrVector vec;
        convert_string_vec(cmix_param_list.c_str(), vec);
        for (auto &s : vec) specs.push_back(parseParamSpec(s));
        user_input_param = true;
    }

    // ------------------------------------------------------------------
    // 3. Parse the kappa specification carried inside model_name
    //    (e.g. "GY{1.0}", "GY{?}" or just "GY").
    // ------------------------------------------------------------------
    ParamSpec kappa_spec; // default: !has_value, !optimize  (= use base-model defaults)
    {
        auto s = model_name.find('{');
        if (s != string::npos) {
            auto t = model_name.find('}', s);
            if (t != string::npos && t > s + 1) {
                kappa_spec = parseParamSpec(model_name.substr(s + 1, t - s - 1));
                model_name = model_name.substr(0, s);
            }
        }
    }
    // String fragment used when building the per-class model_list.  We
    // pass kappa as an explicit value only when the user supplied one;
    // the fix_kappa flag itself is set later (after initMixture).
    string kappa_str;
    if (kappa_spec.has_value)
        kappa_str = "," + convertDoubleToString(kappa_spec.value);

    // ------------------------------------------------------------------
    // 4. Resolve ncat for each model type
    // ------------------------------------------------------------------
    if (cmix_type == "1a")            ncat = 2;
    else if (cmix_type == "2a")       ncat = 3;
    else if (cmix_type == "3" && ncat == 0)  ncat = 3;
    else if (cmix_type == "7" && ncat == 0)  ncat = 10;
    else if (cmix_type == "8" && ncat == 0)  ncat = 11;

    // ------------------------------------------------------------------
    // 5. Validate parameter count when the user supplied a {...} block
    // ------------------------------------------------------------------
    if (user_input_param) {
        size_t expected = 0;
        if (cmix_type == "1a")      expected = 2;
        else if (cmix_type == "2a") expected = 4;
        else if (cmix_type == "3")  expected = 2 * (size_t)ncat;
        else if (cmix_type == "7")  expected = 2;
        else if (cmix_type == "8")  expected = 4;
        else
            outError("Unknown codon mixture " + orig_model_name.substr(cmix_pos));
        if (specs.size() != expected) {
            outError("CMIX" + cmix_type + "{} expects " + convertIntToString((int)expected) +
                     " parameters but got " + convertIntToString((int)specs.size()) +
                     " in " + orig_model_name);
        }
    }

    // ------------------------------------------------------------------
    // 6. Build the model_list (one entry per class)
    // ------------------------------------------------------------------
    string model_list;
    auto appendClass = [&](double init_omega) {
        if (!model_list.empty()) model_list += ",";
        model_list += model_name + "{" + convertDoubleToString(init_omega) + kappa_str + "}";
    };

    if (cmix_type == "1a" || cmix_type == "2a" || cmix_type == "3") {
        vector<double> init_omegas(ncat);

        if (cmix_type == "1a") {
            // class 0: user-controlled omega (default 0.5); class 1: omega = 1
            init_omegas[0] = user_input_param ? resolveValue(specs[1], 0.5) : 0.5;
            init_omegas[1] = 1.0;
        }
        else if (cmix_type == "2a") {
            // class 0 (omega < 1), class 1 (omega = 1), class 2 (omega > 1)
            init_omegas[0] = user_input_param ? resolveValue(specs[1], 0.5) : 0.5;
            init_omegas[1] = 1.0;
            init_omegas[2] = user_input_param ? resolveValue(specs[3], 2.0) : 2.0;
        }
        else { // M3
            if (user_input_param) {
                for (int i = 0; i < ncat; i++)
                    init_omegas[i] = resolveValue(specs[2*i + 1], 0.5);
            } else {
                // Default initial values: spread according to a Beta(1,1)
                // sample (uniform on [0,1]), divided by 0.5 to widen the
                // range -- this matches the previous code's behaviour.
                double *o = RateBeta::SampleOmegas(ncat, 1.0, 1.0);
                for (int i = 0; i < ncat; i++)
                    init_omegas[i] = o[i] / 0.5;
            }
        }
        for (int i = 0; i < ncat; i++) appendClass(init_omegas[i]);
    }
    else if (cmix_type == "7") {
        if (user_input_param) {
            alpha = resolveValue(specs[0], 1.0);
            beta  = resolveValue(specs[1], 1.0);
        }
        double *o = RateBeta::SampleOmegas(ncat, alpha, beta);
        for (int i = 0; i < ncat; i++) appendClass(o[i]);
    }
    else if (cmix_type == "8") {
        double extra_omega_init = 2.0;
        if (user_input_param) {
            alpha            = resolveValue(specs[0], 1.0);
            beta             = resolveValue(specs[1], 1.0);
            extra_omega_init = resolveValue(specs[3], 2.0);
        }
        if (ncat < 2)
            outError("M8 requires at least 2 categories");
        double *o = RateBeta::SampleOmegas(ncat - 1, alpha, beta);
        for (int i = 0; i < ncat - 1; i++) appendClass(o[i]);
        appendClass(extra_omega_init);
    }
    else {
        outError("Unknown codon mixture " + orig_model_name.substr(cmix_pos));
    }

    // ------------------------------------------------------------------
    // 7. Construct the underlying mixture
    // ------------------------------------------------------------------
    initMixture(orig_model_name, model_name, model_list, models_block,
                freq, freq_params, tree, optimize_weights);

    // ------------------------------------------------------------------
    // 8. Override per-class omega bounds and fix flags according to the
    //    parsed user input (or defaults).
    // ------------------------------------------------------------------
    if (cmix_type == "1a") {
        // Class 0: 0 < omega < 1
        ModelCodon *m0 = (ModelCodon*)at(0);
        m0->min_omega = MIN_OMEGA_KAPPA;
        m0->max_omega = 0.999;
        if (m0->omega < m0->min_omega || m0->omega > m0->max_omega)
            m0->omega = 0.5;
        m0->fix_omega = user_input_param ? !specs[1].optimize : false;
        if (m0->fix_omega && m0->omega >= 1.0)
            outError("Omega of class 0 in M1a must be < 1");
        // Class 1: omega is permanently fixed at 1.0
        ModelCodon *m1 = (ModelCodon*)at(1);
        m1->omega = 1.0;
        m1->min_omega = m1->max_omega = 1.0;
        m1->fix_omega = true;
    }
    else if (cmix_type == "2a") {
        ModelCodon *m0 = (ModelCodon*)at(0);
        m0->min_omega = MIN_OMEGA_KAPPA;
        m0->max_omega = 0.999;
        if (m0->omega < m0->min_omega || m0->omega > m0->max_omega)
            m0->omega = 0.5;
        m0->fix_omega = user_input_param ? !specs[1].optimize : false;
        if (m0->fix_omega && m0->omega >= 1.0)
            outError("Omega of class 0 in M2a must be < 1");

        ModelCodon *m1 = (ModelCodon*)at(1);
        m1->omega = 1.0;
        m1->min_omega = m1->max_omega = 1.0;
        m1->fix_omega = true;

        ModelCodon *m2 = (ModelCodon*)at(2);
        m2->min_omega = 1.001;
        m2->max_omega = MAX_OMEGA_KAPPA;
        if (m2->omega < m2->min_omega) m2->omega = 2.0;
        m2->fix_omega = user_input_param ? !specs[3].optimize : false;
        if (m2->fix_omega && m2->omega <= 1.0)
            outError("Omega of class 2 in M2a must be > 1");
    }
    else if (cmix_type == "3") {
        for (int i = 0; i < ncat; i++) {
            ModelCodon *mi = (ModelCodon*)at(i);
            mi->min_omega = MIN_OMEGA_KAPPA;
            mi->max_omega = MAX_OMEGA_KAPPA;
            mi->fix_omega = user_input_param ? !specs[2*i + 1].optimize : false;
        }
    }
    else if (cmix_type == "7") {
        // Per-class omegas are derived from the beta distribution; never
        // optimised independently.  Hide them from BFGS.
        for (int i = 0; i < (int)size(); i++)
            ((ModelCodon*)at(i))->fix_omega = true;
        if (user_input_param) {
            fix_alpha = !specs[0].optimize;
            fix_beta  = !specs[1].optimize;
        }
    }
    else if (cmix_type == "8") {
        for (int i = 0; i < (int)size() - 1; i++)
            ((ModelCodon*)at(i))->fix_omega = true;
        // Extra (positive-selection) class: omega > 1
        ModelCodon *mlast = (ModelCodon*)at(size() - 1);
        mlast->min_omega = 1.001;
        mlast->max_omega = MAX_OMEGA_KAPPA;
        if (mlast->omega < mlast->min_omega) mlast->omega = 2.0;
        // The extra category's omega is managed by getVariables/setVariables
        // through fix_extra_omega -- always mark fix_omega so the per-class
        // BFGS path leaves it alone.
        mlast->fix_omega = true;
        if (user_input_param) {
            fix_alpha        = !specs[0].optimize;
            fix_beta         = !specs[1].optimize;
            fix_extra_weight = !specs[2].optimize;
            fix_extra_omega  = !specs[3].optimize;
        }
    }

    // ------------------------------------------------------------------
    // 9. Apply kappa overrides.  Kappa is shared across all classes
    //    (link_kappa = true) and getVariables() copies it from class 0
    //    onto every class on each iteration, so we set it on every class
    //    here for consistency.
    // ------------------------------------------------------------------
    for (size_t i = 0; i < size(); i++) {
        ModelCodon *mi = (ModelCodon*)at(i);
        if (kappa_spec.has_value) {
            mi->kappa = kappa_spec.value;
            mi->fix_kappa = !kappa_spec.optimize;
        } else if (kappa_spec.optimize) {
            // user wrote GY{?}: optimise from the existing default
            mi->fix_kappa = false;
        }
        // else: leave whatever the underlying ModelCodon defaulted to.
    }

    // ------------------------------------------------------------------
    // 10. Kappa is shared (linked) across all mixture classes —
    //     getVariables() copies it from class 0 to every other class on
    //     each iteration.  If classes 1..K-1 also expose their own kappa
    //     to BFGS, those K-1 slots become degenerate dimensions: BFGS
    //     adjusts them, then getVariables immediately overwrites them,
    //     which can drive BFGS into bogus local optima (the same root
    //     cause that the upstream M7/M8 freeze addresses).
    //
    //     Freeze kappa on every class except class 0 (which holds the
    //     single shared kappa that matters), and refresh num_params so
    //     the BFGS dimension count is correct.  This applies to ALL
    //     codon mixture types, not only M7/M8.
    // ------------------------------------------------------------------
    for (int i = 0; i < (int)size(); i++) {
        ModelCodon *m = (ModelCodon*)at(i);
        if (i > 0)
            m->fix_kappa = true;
        m->updateNumParams();
    }

    // ------------------------------------------------------------------
    // 11. Set initial mixture weights and decide whether to fix them
    // ------------------------------------------------------------------
    if (cmix_type == "1a" || cmix_type == "2a" || cmix_type == "3") {
        vector<double> init_props((size_t)ncat, 1.0 / ncat);
        bool any_weight_optimised = !user_input_param;

        if (user_input_param) {
            if (cmix_type == "1a") {
                init_props[0] = resolveValue(specs[0], 0.5);
                init_props[1] = 1.0 - init_props[0];
                if (specs[0].optimize) any_weight_optimised = true;
            }
            else if (cmix_type == "2a") {
                init_props[0] = resolveValue(specs[0], 0.4);
                init_props[2] = resolveValue(specs[2], 0.2);
                init_props[1] = 1.0 - init_props[0] - init_props[2];
                if (specs[0].optimize || specs[2].optimize)
                    any_weight_optimised = true;
            }
            else { // M3
                for (int i = 0; i < ncat; i++) {
                    init_props[i] = resolveValue(specs[2*i], 1.0 / ncat);
                    if (specs[2*i].optimize) any_weight_optimised = true;
                }
            }
        }

        // Sanity-check and normalise
        double wsum = 0.0;
        for (double p : init_props) {
            if (p < 0.0)
                outError("Mixture weight cannot be negative in " + orig_model_name);
            wsum += p;
        }
        if (wsum <= 0.0)
            outError("Sum of mixture weights must be > 0 in " + orig_model_name);
        for (auto &p : init_props) p /= wsum;
        for (int i = 0; i < ncat; i++) prop[i] = init_props[i];

        fix_prop = !any_weight_optimised;
    }
    else if (cmix_type == "7") {
        // Equal weights, not user-controllable for M7
        for (int i = 0; i < ncat; i++) prop[i] = 1.0 / ncat;
        fix_prop = true;
    }
    else if (cmix_type == "8") {
        // p_0 in the PDF table is the weight of the beta-distributed
        // sub-component (combined over all K-1 classes); the extra
        // positive-selection category therefore takes weight 1 - p_0.
        double p_beta_init = 0.9;  // sensible default: most sites neutral/purifying
        if (user_input_param)
            p_beta_init = resolveValue(specs[2], 0.9);
        if (p_beta_init < 0.0 || p_beta_init > 1.0)
            outError("M8: weight of beta-distributed component must be in [0,1]");
        double p_extra = 1.0 - p_beta_init;
        double share   = p_beta_init / (ncat - 1);
        for (int i = 0; i < ncat - 1; i++) prop[i] = share;
        prop[ncat - 1] = p_extra;
        // Mixture weights for M8 are managed by getVariables/setVariables;
        // the underlying ModelMixture optimisation is bypassed.
        fix_prop = true;
    }

    // Show the initial parameters
    cout << "Initial parameters in the Codon Mixture:" << endl;
    writeInfo(cout);
    cout << endl;
}

ModelCodonMixture::~ModelCodonMixture()
{

}

// ----------------------------------------------------------------------
// M7/M8 extra-parameter table
// ----------------------------------------------------------------------
//
// Returns the list of extra parameters (beta-distribution alpha/beta and
// the M8 extra-category omega/weight) in the order they should be packed
// into the variables array, starting at index ModelMixture::getNDim()+1.
// Fixed parameters are kept in the list (with fixed=true) but excluded
// from the dimension count and from BFGS bounds.
//
// The bounds for alpha/beta were widened from the historical [0.01, 10]
// to [0.01, 50] so the optimizer isn't pinned against 10 on datasets
// whose true beta is close to the old cap (e.g. dataset 170 has
// beta = 7.52 and BFGS used to terminate at beta = 10 at a worse
// log-likelihood).
//
std::vector<ModelCodonMixture::ExtraParam> ModelCodonMixture::getExtraParams() {
    std::vector<ExtraParam> out;
    if (cmix_subtype == "7") {
        out.push_back({&alpha, 0.01, 50.0,  "alpha", fix_alpha});
        out.push_back({&beta,  0.01, 50.0,  "beta",  fix_beta});
    }
    else if (cmix_subtype == "8") {
        // The order here must match the historical layout used by the
        // checkpoint/optimisation code: omega_extra, weight_extra, alpha, beta
        ModelCodon *mlast = (ModelCodon*)at(size() - 1);
        out.push_back({&mlast->omega,   1.01,  10.0,  "omega_extra",  fix_extra_omega});
        out.push_back({&prop[size()-1], 0.001, 0.999, "weight_extra", fix_extra_weight});
        out.push_back({&alpha,          0.01,  50.0,  "alpha",        fix_alpha});
        out.push_back({&beta,           0.01,  50.0,  "beta",         fix_beta});
    }
    return out;
}

bool ModelCodonMixture::getVariables(double *variables) {
    bool changed = ModelMixture::getVariables(variables);
    auto kappa  = ((ModelCodon*)at(0))->kappa;
    auto kappa2 = ((ModelCodon*)at(0))->kappa2;

    if (cmix_subtype == "7" || cmix_subtype == "8") {
        // Walk the extra-parameter list, reading from variables[] for the
        // free parameters and leaving fixed ones at their stored value.
        int idx = ModelMixture::getNDim();   // last index used by the base mixture
        auto extras = getExtraParams();
        for (auto &ep : extras) {
            if (!ep.fixed) {
                idx++;
                *ep.target = variables[idx];
            }
        }
        // Recompute the per-class omega values from the (possibly updated)
        // alpha/beta and (for M8) propagate the extra category's weight.
        if (cmix_subtype == "7") {
            double *omega = RateBeta::SampleOmegas(size(), alpha, beta);
            for (int i = 0; i < (int)size(); i++) {
                ModelCodon *m = (ModelCodon*)at(i);
                m->omega = omega[i];
            }
        } else { // M8
            double *omega = RateBeta::SampleOmegas(size() - 1, alpha, beta);
            double w_extra = prop[size() - 1];
            for (int i = 0; i < (int)size() - 1; i++) {
                ModelCodon *m = (ModelCodon*)at(i);
                m->omega = omega[i];
                prop[i] = (1.0 - w_extra) / (size() - 1);
            }
            ModelCodon *mlast = (ModelCodon*)at(size() - 1);
            cout << "alpha: " << alpha << "\tbeta: " << beta
                 << "\tomega_free: " << mlast->omega
                 << "\tfree_weight: " << w_extra
                 << "\tkappa: " << kappa << endl;
        }
    }

    // Propagate kappa from class 0 onto every other class (link_kappa)
    for (int i = 1; i < (int)size(); i++) {
        ModelCodon *m = (ModelCodon*)at(i);
        m->kappa  = kappa;
        m->kappa2 = kappa2;
    }
    rescale_codon_mix();
    return changed;
}

int ModelCodonMixture::getNDim() {
    int base = ModelMixture::getNDim();
    if (cmix_subtype == "7" || cmix_subtype == "8") {
        int extras = 0;
        // Inline counting (avoids constructing the helper vector here)
        if (cmix_subtype == "7") {
            extras = (fix_alpha ? 0 : 1) + (fix_beta ? 0 : 1);
        } else {
            extras = (fix_extra_omega  ? 0 : 1)
                   + (fix_extra_weight ? 0 : 1)
                   + (fix_alpha        ? 0 : 1)
                   + (fix_beta         ? 0 : 1);
        }
        return base + extras;
    }
    return base;
}

void ModelCodonMixture::setVariables(double *variables) {
    ModelMixture::setVariables(variables);
    if (cmix_subtype == "7" || cmix_subtype == "8") {
        int idx = ModelMixture::getNDim();
        auto extras = getExtraParams();
        for (auto &ep : extras) {
            if (!ep.fixed) {
                idx++;
                variables[idx] = *ep.target;
            }
        }
    }
}

void ModelCodonMixture::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
    ModelMixture::setBounds(lower_bound, upper_bound, bound_check);
    if (cmix_subtype == "7" || cmix_subtype == "8") {
        int idx = ModelMixture::getNDim();
        auto extras = getExtraParams();
        for (auto &ep : extras) {
            if (!ep.fixed) {
                idx++;
                lower_bound[idx] = ep.lo;
                upper_bound[idx] = ep.hi;
                bound_check[idx] = true;
            }
        }
    }
}

// Impose restrictions on the omega values (kept for API compatibility,
// but the constructor now sets bounds directly so this is a no-op).
void ModelCodonMixture::restrict_omega_values(string /*cmix_type*/) {
}

double ModelCodonMixture::optimizeParameters(double gradient_epsilon) {

    int dim = getNDim();
    if (dim == 0)
        return 0.0;
    double score = 0.0;
    IntVector params;
    int i, j, ncategory = size();
    if (dim != 0) {
        score = 1.0;
    }

    if (!phylo_tree->getModelFactory()->unobserved_ptns.empty()) {
        outError("Mixture model +ASC is not supported yet. Contact author if needed.");
    }

    score = phylo_tree->computeLikelihood();
    cout << "before parameter optimization, score = " << score << endl;

    if (Params::getInstance().optimize_alg_qmix == "EM") {
        assert(size() > 0);
        ModelCodon* first_model = (ModelCodon*)at(0);
        if (!first_model->fix_kappa) {
            // fix all the omega values
            bool* orig_fix_omega = new bool[size()];
            for (i = 0; i < size(); i++) {
                orig_fix_omega[i] = ((ModelCodon*)at(i))->fix_omega;
                ((ModelCodon*)at(i))->fix_omega = true;
            }
            // fix the weights
            bool orig_fix_prop = fix_prop;
            fix_prop = true;
            // optimize the shared kappa value using BFGS
            score = ModelMarkov::optimizeParameters(gradient_epsilon);
            cout << "after shared kappa optimization, score = " << score << endl;

            // restore the original values
            for (i = 0; i < size(); i++) {
                ((ModelCodon*)at(i))->fix_omega = orig_fix_omega[i];
            }
            fix_prop = orig_fix_prop;
            delete[] orig_fix_omega;
        }
        // optimize the other parameters using EM algorithm
        // fix the kappa values
        bool orig_fix_kappa = first_model->fix_kappa;
        first_model->fix_kappa = true;
        score = optimizeWithEM(gradient_epsilon);
        // restore the original value
        first_model->fix_kappa = orig_fix_kappa;
        cout << "after EM optimization ";
    } else if (Params::getInstance().optimize_alg_qmix == "BFGS") {
        // BFGS-then-EM-weights coordinate descent can get stuck in local
        // optima where the mixture weights stay near 1/K. The mechanism:
        // when we enter the round with the initial uniform weights (or
        // weights very close to uniform), BFGS — running with the weights
        // fixed — tunes the per-class omegas so the data is best explained
        // *given* uniform weights. After BFGS the omegas have been pushed
        // into a configuration in which every class accounts for ~1/K of
        // the sites by construction, so the subsequent EM weight update
        // has no posterior signal to move the weights anywhere else. Once
        // we are in that basin we stay there for the rest of the run.
        //
        // We break the symmetry by running ONE EM iteration *before* BFGS,
        // on the very first call. With the original (well-spread) initial
        // omegas, the posterior of each pattern under each class differs
        // and the EM update produces non-uniform weights. BFGS then tunes
        // the omegas in that non-uniform-weight regime, which leads to a
        // strictly better optimum on data sets where the true mixture is
        // far from uniform (e.g. dataset3 / GY+M3.3 where the true weights
        // are 0.657, 0.113, 0.230). The change is monotone-safe because
        // it only affects which initial point BFGS sees, not the BFGS
        // itself.
        bool orig_fix_prop = fix_prop;
        if (iteration_num == 0 && !fix_prop) {
            double pre_score = optimizeWeights(1);
            cout << "before parameter optimization, after EM weight init, score = "
                 << pre_score << endl;
        }
        // first optimize the other parameters using BFGS
        fix_prop = true;

        // Multi-start over (alpha, beta) for M7/M8 on the very first call.
        // The (alpha=1, beta=1) default drops some datasets (e.g. 191, with
        // true alpha=0.12, beta=5.75) into a U-shape local optimum
        // hundreds of log-likelihood units worse than the truth. Trying a
        // small grid of starting shape parameters and keeping the best BFGS
        // result reliably escapes that basin while never making good cases
        // worse (we always keep the best score). We also re-optimise branch
        // lengths between starts so that each one is compared on its own
        // best tree, not on a tree fitted to the previous start's model.
        //
        // The multistart is skipped when the user has *fixed* both alpha
        // and beta via the new "+CMIX7/8{...}" syntax, since in that case
        // there is nothing to search over.
        bool can_multistart = (cmix_subtype == "7" || cmix_subtype == "8")
                              && !(fix_alpha && fix_beta);
        if (!multistart_done && can_multistart) {
            multistart_done = true;
            // A compact set of starting points that covers the qualitative
            // shapes of Beta(alpha, beta) on [0,1]: uniform, sharp decay,
            // moderate decay, broad bell, U-shape. Kept small because each
            // start runs a full model+BL alternation and we pay 5x this
            // cost.  The BFGS dimension here is tiny (3 for M7, 5 for M8).
            const std::vector<std::pair<double,double>> start_points = {
                {1.0, 1.0},   // uniform (historical default)
                {0.1, 7.0},   // very sharp decay toward 0
                {2.0, 5.0},   // moderate decay
                {4.0, 2.0},   // broad bell peaked > 0.5
                {0.5, 0.5},   // U-shape
            };
            // Snapshot the initial state so every start runs BFGS from
            // identical kappa / prop / free-omega / branch-lengths and only
            // differs in (alpha, beta). This ensures the starts are
            // independent samples of the landscape rather than a chained
            // sequence whose later starts inherit drift from earlier ones.
            double init_kappa = ((ModelCodon*)at(0))->kappa;
            double init_free_omega = (cmix_subtype == "8")
                ? ((ModelCodon*)at(size()-1))->omega : 0.0;
            DoubleVector init_prop(prop, prop + size());
            DoubleVector init_bl;
            phylo_tree->saveBranchLengths(init_bl);

            double best_score = -DBL_MAX;
            double best_alpha = alpha, best_beta = beta;
            double best_kappa = init_kappa;
            double best_free_omega = init_free_omega;
            double best_free_weight = init_prop[size()-1];
            DoubleVector best_prop = init_prop;
            DoubleVector best_bl = init_bl;
            for (size_t s = 0; s < start_points.size(); s++) {
                // If the user fixed alpha (or beta), keep the user's value
                // for that coordinate and only vary the free one.
                alpha = fix_alpha ? alpha : start_points[s].first;
                beta  = fix_beta  ? beta  : start_points[s].second;
                // restore kappa / prop / free-omega / BL from the snapshot
                // so every start begins from the same baseline
                for (int k = 0; k < size(); k++) {
                    ((ModelCodon*)at(k))->kappa = init_kappa;
                    prop[k] = init_prop[k];
                }
                if (cmix_subtype == "8") {
                    ((ModelCodon*)at(size()-1))->omega = init_free_omega;
                }
                phylo_tree->restoreBranchLengths(init_bl);
                phylo_tree->clearAllPartialLH();
                // Alternate model + BL optimisation a couple of times so
                // each start converges fairly. Without the BL refit, the
                // first start's BL bleeds into all subsequent starts and
                // BFGS sees a stale tree, "moving" the model away from a
                // perfectly good (alpha, beta) basin.
                double s_score = 0.0;
                for (int it = 0; it < 2; it++) {
                    s_score = ModelMarkov::optimizeParameters(gradient_epsilon);
                    s_score = phylo_tree->optimizeAllBranches(1,
                                                              gradient_epsilon);
                }
                cout << "  multistart point (alpha=" << start_points[s].first
                     << ", beta=" << start_points[s].second << ") -> score = "
                     << s_score << " (final alpha=" << alpha
                     << ", beta=" << beta << ")" << endl;
                if (s_score > best_score) {
                    best_score = s_score;
                    best_alpha = alpha;
                    best_beta  = beta;
                    best_kappa = ((ModelCodon*)at(0))->kappa;
                    if (cmix_subtype == "8") {
                        best_free_omega  = ((ModelCodon*)at(size()-1))->omega;
                        best_free_weight = prop[size()-1];
                    }
                    for (int k = 0; k < size(); k++) best_prop[k] = prop[k];
                    phylo_tree->saveBranchLengths(best_bl);
                }
            }
            // Restore the best start's converged state (model + branch
            // lengths) so the outer iteration loop continues from it.
            alpha = best_alpha;
            beta  = best_beta;
            for (int k = 0; k < size(); k++) {
                ((ModelCodon*)at(k))->kappa = best_kappa;
                prop[k] = best_prop[k];
            }
            if (cmix_subtype == "8") {
                ((ModelCodon*)at(size()-1))->omega = best_free_omega;
                prop[size()-1] = best_free_weight;
            }
            // Recompute the per-class omegas from (alpha, beta) and refresh
            // the Q matrices so the model state is consistent with the
            // restored (alpha, beta, kappa). Without this, computeLikelihood
            // would use whichever per-class omegas the LAST start left in
            // the components, not the ones implied by best_alpha/beta.
            {
                double* omega_arr = (cmix_subtype == "7")
                    ? RateBeta::SampleOmegas((int)size(), alpha, beta)
                    : RateBeta::SampleOmegas((int)size()-1, alpha, beta);
                int n_beta = (cmix_subtype == "7") ? (int)size() : (int)size() - 1;
                for (int k = 0; k < n_beta; k++)
                    ((ModelCodon*)at(k))->omega = omega_arr[k];
                if (cmix_subtype == "8")
                    ((ModelCodon*)at(size()-1))->omega = best_free_omega;
            }
            phylo_tree->restoreBranchLengths(best_bl);
            rescale_codon_mix();   // refresh Q matrices, decompose, clear LH
            score = phylo_tree->computeLikelihood();
            cout << "after multistart BFGS, best score = " << score
                 << " (alpha=" << alpha << ", beta=" << beta << ")" << endl;
        } else {
            score = ModelMarkov::optimizeParameters(gradient_epsilon);
            cout << "after parameter optimization, score = " << score << endl;
        }
        // then optimize the weights using EM
        fix_prop = orig_fix_prop;
        if (!fix_prop) {
            // Save total_num_subst before EM. rescale_codon_mix() (called inside
            // optimizeWeights) will change total_num_subst to re-normalise the
            // mixture, which is equivalent to changing the effective branch-length
            // scale by R = total_num_subst_new / total_num_subst_old.
            // Scaling the branch lengths by the same factor R exactly compensates,
            // so the post-EM likelihood is guaranteed >= pre-EM likelihood.
            score = optimizeWeights(++iteration_num);
            cout << "after weight optimization (EM) and rescaling, score = " << score << endl;
            return score;
        }
    } else {
        score = ModelMarkov::optimizeParameters(gradient_epsilon);
        cout << "after all-BFGS parameter optimization ";
    }

    // rescale the Codon Q Matrices (for EM/all-BFGS paths, or BFGS with fixed weights)
    rescale_codon_mix();
    score = phylo_tree->computeLikelihood();
    cout << "and rescaling, score = " << score << endl;


    return score;
}

/**
    write information
    @param out output stream
*/
void ModelCodonMixture::writeInfo(ostream &out) {
    if (cmix_subtype == "7" || cmix_subtype == "8") {
        cout << "Beta distribution -- " << "alpha: " << alpha << ", beta: " << beta << endl;
    }
    ModelMixture::writeInfo(out);
}
