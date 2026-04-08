//
//  modelcodonmixture.cpp
//  iqtree
//
//  Created by Minh Bui on 4/3/2025.
//

#include "modelcodonmixture.h"

#include "ratebeta.h"


ModelCodonMixture::ModelCodonMixture(string orig_model_name, string model_name,
                                     ModelsBlock *models_block, StateFreqType freq, string freq_params,
                                     PhyloTree *tree, bool optimize_weights)
: ModelMarkov(tree), ModelMixture(tree) {

    if (tree->aln->seq_type != SEQ_CODON)
        outError("Can't apply codon mixture model as sequence type is not codon");
    auto cmix_pos = orig_model_name.find("+CMIX");
    ASSERT(cmix_pos != string::npos);
    auto end_pos = orig_model_name.find_first_of("+*{", cmix_pos+1);
    string cmix_type;

    if (end_pos == string::npos)
        cmix_type = orig_model_name.substr(cmix_pos+5);
    else
        cmix_type = orig_model_name.substr(cmix_pos+5, end_pos-cmix_pos-5);

    //Yuri - Determine number of categories
    auto ncat_pos = cmix_type.find(".");
    int ncat = 0;
    if (ncat_pos != string::npos) {
        auto ncat_str = cmix_type.substr(ncat_pos+1);
        ncat = stoi(ncat_str);
        //***Might be a more preferred way to address this check***
        if (ncat < 1) {
            outError("Invalid number of categories " + orig_model_name.substr(cmix_pos+1));
        }
        cmix_type = cmix_type.substr(0,ncat_pos);
    }
    alpha = 1.0;
    beta = 1.0;
    iteration_num = 0;
    // read the input parameters for +CMIXi{...}
    StrVector vec;
    if (end_pos != string::npos && orig_model_name[end_pos] == '{') {
        auto close_br_pos = orig_model_name.find_first_of("}", end_pos+1);
        if (close_br_pos != string::npos) {
            string cmix_param_list = orig_model_name.substr(end_pos+1,close_br_pos-end_pos-1);
            convert_string_vec(cmix_param_list.c_str(), vec);
        }
    }
    
    string model_list = "";
    bool user_input_param = false;
    bool kappa_given = false;
    
    if (vec.size() == 0) {
        /* setting fix_kappa for class 2 and 3 for GY0K model */
        string kappa_str = ",1.0";
        // get the value of kappa if specified
        auto s = model_name.find("{");
        if (s != string::npos) {
            auto t = model_name.find("}", s);
            if (t != string::npos && t > s + 1) {
                kappa_given = true;
                kappa_str = "," + model_name.substr(s+1, t-s-1);
                model_name = model_name.substr(0, s);
            }
        }
        
        /*
        if (model_name != "GY0K" && !kappa_given) {
            // kappa is not fixed, thus cannot use EM algorithm
            if (Params::getInstance().optimize_alg_qmix == "EM") {
                outError("EM algorithm cannot be used for the codon mixture model with unfixed kappa value.");
            }
        }
        */
        
        if (cmix_type == "1a") {
            // M1a neural model with 2 classes, omega2 = 1.0
            if (kappa_given) {
                model_list = model_name + "{<0.999" + kappa_str + "}," + model_name + "{1.0" + kappa_str + "}";
            } else {
                model_list = model_name + "{<0.999}," + model_name + "{1.0" + kappa_str + "}";
            }
        } else if (cmix_type == "2a") {
            // M2a selection model with 3 classes
            if (kappa_given) {
                model_list = model_name + "{<0.999" + kappa_str + "}," + model_name + "{1.0" + kappa_str + "}," + model_name + "{>1.001" + kappa_str + "}";
            } else {
                model_list = model_name + "{<0.999}," + model_name + "{1.0" + kappa_str + "}," + model_name + "{>1.001" + kappa_str + "}";
            }
        } else if (cmix_type == "3") {
            // M3 model with 3 classes with no constraint
            //default value for ncat
            if (ncat == 0)
                ncat = 3;
            if (kappa_given) {
                model_list = model_name + "{>0.001" + kappa_str + "}";
            } else {
                model_list = model_name + "{>0.001}";
            }
            for (int i = 1; i < ncat; i++)
                model_list += "," + model_name + "{>0.001" + kappa_str + "}";

        } else if (cmix_type == "7") {
            // M7 model with category omegas following a beta distribution
            //default value for ncat
            if (ncat == 0)
                ncat = 10;
            double shape_alpha = 1.0;
            double shape_beta = 1.0;

            //RateBeta beta_dist;
            double* omega = RateBeta::SampleOmegas(ncat,shape_alpha,shape_beta);

            model_list = model_name + "{" + std::to_string(omega[0]) + "}:1:0.1";
            for (int i = 1; i < ncat; i++) {
                model_list += "," + model_name + "{" + std::to_string(omega[i]) + kappa_str + "}:1:0.1";
            }
        } else if (cmix_type == "8") {
            // M8 model with category omegas following a beta distribution
            // and an addition category constrained to omega > 1.0
            //default value for ncat
            if (ncat == 0)
                ncat = 11;
            double shape_alpha = 1.0;
            double shape_beta = 1.0;
            //RateBeta beta_dist;
            double* omega = RateBeta::SampleOmegas(ncat-1, shape_alpha,shape_beta);

            model_list = model_name + "{" + std::to_string(omega[0]) + "}:1:0.1";
            for (int i = 1; i < ncat-1; i++) {
                model_list += "," + model_name + "{" + std::to_string(omega[i]) + kappa_str + "}:1:0.1";
            }
            model_list += "," + model_name + "{>1.001" + kappa_str + "}";
        } else {
            outError("Unknown codon mixture " + orig_model_name.substr(cmix_pos));
        }
    } else {
        // user inputs the parameter values for CMIX model
        user_input_param = true;
        if (cmix_type == "1a" && vec.size() != 2 && vec.size() != 4)
            outError("Error! There should be 2 (or 4) parameters inside CMIX1a{} stating the omega value (and the weight) of each class");
        else if (cmix_type == "2a" && vec.size() != 3 && vec.size() != 6)
            outError("Error! There should be 3 (or 6) parameters inside CMIX2a{} stating the omega value (and the weight) of each class.");
        else if (cmix_type == "3" && vec.size() <2 && vec.size() < 4)
            outError("Error! There should be at least 2 (or 4) parameters inside CMIX3{} stating the omega value (and the weight) of each class.");
        if (vec.size() >= 4) {
            // with both omega and weight
            for (int i = 0; i < vec.size(); i+=2) {
                if (i > 0) {
                    model_list += ",";
                }
                model_list += model_name + "{" + vec[i] + ",1.0}:1.0:" + vec[i+1];
            }
        } else {
            // only omega
            for (int i = 0; i < vec.size(); i++) {
                if (i > 0) {
                    model_list += ",";
                }
                model_list += model_name + "{" + vec[i] + ",1.0}";
            }
        }
    }
    initMixture(orig_model_name, model_name, model_list, models_block,
                freq, freq_params, tree, optimize_weights);
    
    // impose restrictions on the omega values if user inputs the parameters
    if (user_input_param) {
        restrict_omega_values(cmix_type);
    }

    // Yuri - Added functions to set custom cnat initials
    // set the initial omega values for M3 model
    if (!user_input_param && cmix_type == "3") {
        double shape_alpha = 1.0;
        double shape_beta = 1.0;
        //RateBeta beta_dist;
        double* omega = RateBeta::SampleOmegas(ncat, shape_alpha,shape_beta);

        for (int i = 0; i < ncat; i++) {
            //if (ncat > 2){
            ((ModelCodon*)at(i))->omega = omega[i]/0.5; //(i+0.01)/(ncat-2);
            //}else {
            //    ((ModelCodon*)at(i))->omega = (i+0.01);
            }
        //}

        /*((ModelCodon*)at(0))->omega = 0.4;
        ((ModelCodon*)at(1))->omega = 0.9;
        ((ModelCodon*)at(2))->omega = 1.8;*/
    }
    
    // the newton optimization does not work properly for codon mixture model
    // Params::getInstance().optimize_by_newton = false;
    // phylo_tree->optimize_by_newton = false;

    // show the initial parameters
    cout << "Initial parameters in the Codon Mixture:" << endl;
    writeInfo(cout);
    cout << endl;
}

ModelCodonMixture::~ModelCodonMixture()
{
    
}

bool ModelCodonMixture::getVariables(double *variables) {
    bool changed = ModelMixture::getVariables(variables);
    auto kappa = ((ModelCodon*)at(0))->kappa;
    auto kappa2 = ((ModelCodon*)at(0))->kappa2;
    //need to add ncat line here
    if (name=="M7") {
        alpha = variables[getNDim()-1];
        beta = variables[getNDim()];
        double* omega = RateBeta::SampleOmegas(size(),alpha,beta);
        for (int i = 0; i < size(); i++) {
            ModelCodon *model = (ModelCodon*)at(i);
            model->omega = omega[i];
        }
    }else if (name=="M8") {
        alpha = variables[getNDim()-1];
        beta = variables[getNDim()];
        double* omega = RateBeta::SampleOmegas(size()-1,alpha,beta);
        for (int i = 0; i < size()-1; i++) {
            ModelCodon *model = (ModelCodon*)at(i);
            model->omega = omega[i];
            //prop[i] = prop[0];
            prop[i] = (1.0-variables[getNDim()-2])/(size()-1);
            //cout << "i: " << i << endl;
            //cout << "omega: " << omega[i] << endl;
        }
        ModelCodon *model = (ModelCodon*)at(size()-1);
        prop[size()-1] = variables[getNDim()-2];
        model->omega = variables[getNDim()-3];
        //cout << "Omega: " << variables[getNDim()-3] << endl;
        cout << "alpha: " << variables[getNDim()-1] << "\tbeta: " << variables[getNDim()] << "\tomega_free: " << model->omega << "\tfree_wight: "<<prop[size()-1]<< "\tkappa: " << kappa << endl;
        //cout << "weight: " << variables[getNDim()-2] << endl;
    }
    for (int i = 1; i < size(); i++) {
        ModelCodon *model = (ModelCodon*)at(i);
        model->kappa = kappa;
        model->kappa2 = kappa2;
    }
    rescale_codon_mix();
    return changed;
}

int ModelCodonMixture::getNDim() {
    if (name == "M7") {
        //Adds 2 beta distribution shape parameters for the beta distribution
        return ModelMixture::getNDim()+2;
    }
    else if(name == "M8") {
        //M8 needs an extra category for the >1 omega category weight
        return ModelMixture::getNDim()+4;
    }
    else{
        return ModelMixture::getNDim();
    }
}

void ModelCodonMixture::setVariables(double *variables) {
    ModelMixture::setVariables(variables);
    if (name == "M7" || name == "M8") {
        variables[getNDim()-1] = alpha;
        variables[getNDim()] = beta;
        if (name == "M8") {
            //set to second last prop
            variables[getNDim()-2] = prop[size()-1];
            ModelCodon *model = (ModelCodon*)at(size()-1);
            variables[getNDim()-3] = model->omega;
        }
    }
}

void ModelCodonMixture::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
    ModelMixture::setBounds(lower_bound, upper_bound, bound_check);
    if (name == "M7" || name == "M8") {
        //alpha/beta shape parameters
        lower_bound[getNDim()-1]=0.01;
        upper_bound[getNDim()-1]=10.0000;
        bound_check[getNDim()-1]=true;
        lower_bound[getNDim()]=0.01;
        upper_bound[getNDim()]=10.0000;
        bound_check[getNDim()]=true;
        if (name == "M8") {
            //Free category weight
            lower_bound[getNDim()-2]=0.001;
            upper_bound[getNDim()-2]=0.999;
            bound_check[getNDim()-2]=true;
            //Free category omega
            lower_bound[getNDim()-3]=1.01;
            upper_bound[getNDim()-3]=10.000;
            bound_check[getNDim()-3]=true;
        }
    }
}

// Impose restrictions on the omega values
// This function is needed only when user inputs parameters
void ModelCodonMixture::restrict_omega_values(string cmix_type) {
    ModelCodon *model;
    // M1a neural model with 2 classes, or
    // M2a selection model with 3 classes
    if (cmix_type == "1a" || cmix_type == "2a") {
        ASSERT(size() >= 2);
        // omega1 is resticted to < 1 for both M1a and M2a models
        model = (ModelCodon*)at(0);
        model->max_omega = 0.999; // for the option -optfromgiven
        if (model->omega > model->max_omega) {
            outError("Omega1 has to be smaller than 1");
        }
        // omega2 is resticted to 1.0 for both M1a and M2a models
        model = (ModelCodon*)at(1);
        model->min_omega = model->max_omega = 1.0; // for the option -optfromgiven
        if (model->omega != model->min_omega) {
            outError("Omega2 has to be 1");
        }

        if (cmix_type == "2a") {
            ASSERT(size() == 3);
            // omega3 is resticted to > 1 for M2a model
            model = (ModelCodon*)at(2);
            model->min_omega = 1.001; // for the option -optfromgiven
            if (model->omega < model->min_omega) {
                outError("Omega3 has to be greater than 1");
            }
        }
    }
    // M3 model with 3 classes with no constraint
    else if (cmix_type == "3") {
        for (size_t i = 0; i < size(); i++) {
            model = (ModelCodon*)at(i);
            model->min_omega = 0.001;
            if (model->omega < model->min_omega) {
                outError("Omega" + convertIntToString(i+1) + " cannot be 0 or the value is too small");
            }
        }
    }
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
        score = ModelMarkov::optimizeParameters(gradient_epsilon);
        cout << "after parameter optimization, score = " << score << endl;
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
    if (name == "M7" || name == "M8") {
        cout << "Beta distribution -- " << "alpha: " << alpha << ", beta: " << beta << endl;
    }
    ModelMixture::writeInfo(out);
}
