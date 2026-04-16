//
//  modelcodonmixture.h
//  iqtree
//
//  Created by Minh Bui on 4/3/2025.
//

#ifndef modelcodonmixture_h
#define modelcodonmixture_h

#include <stdio.h>
#include <vector>
#include "modelcodon.h"
#include "modelmixture.h"

/** Codon Mixture models like M2a, M3, M7, M8, notations following PAML package */
class ModelCodonMixture : public ModelMixture {
public:
    
    /**
        constructor
        @param model_name model name, e.g., CM2a, CM3
        @param freq state frequency type
        @param tree associated phylogenetic tree
    */
    ModelCodonMixture(string orig_model_name, string model_name, ModelsBlock *models_block,
            StateFreqType freq, string freq_params, PhyloTree *tree, bool optimize_weights);

    /**
     * destructor
     */
    virtual ~ModelCodonMixture();

    //shape parameters for the beta distribution (M7 and M8 model)
    double alpha;
    double beta;

    // Fix flags for the M7/M8 beta distribution shape parameters and the
    // M8 extra category. These are honoured by getNDim()/getVariables()/
    // setVariables()/setBounds() so that the user can request fixing some
    // and optimising others via the new "?", "number", "number?" syntax.
    bool fix_alpha = false;
    bool fix_beta = false;
    // M8 extra (positive-selection) category
    bool fix_extra_omega = false;
    bool fix_extra_weight = false;

    // remember the cmix sub-type ("1a", "2a", "3", "7", "8") so other
    // helper functions can branch on it without re-parsing the model name
    string cmix_subtype;

private:

    bool link_kappa = true;

    // iteration #
    int iteration_num;

    // set once the multistart over (alpha, beta) has been performed, so
    // we don't pay its cost on every subsequent call to optimizeParameters
    bool multistart_done = false;

    // Helper used by getVariables/setVariables/setBounds for M7/M8.
    // Returns a list of (target pointer, lower, upper, label, fixed) for
    // the per-mixture-model "extra" parameters in the order they should be
    // packed into the variables array.  Index 0 of the returned vector
    // corresponds to ModelMixture::getNDim()+1 in the variables array.
    struct ExtraParam {
        double *target;     // points to the live value (alpha, beta, omega, weight)
        double lo, hi;
        const char *label;
        bool fixed;
    };
    std::vector<ExtraParam> getExtraParams();

protected:
    /**
        this function is served for the multi-dimension optimization. It should pack the model parameters
        into a vector that is index from 1 (NOTE: not from 0)
        @param variables (OUT) vector of variables, indexed from 1
    */
    virtual void setVariables(double *variables);

    void setBounds(double *lower_bound, double *upper_bound, bool *bound_check);

    /**
        this function is served for the multi-dimension optimization. It should assign the model parameters
        from a vector of variables that is index from 1 (NOTE: not from 0)
        @param variables vector of variables, indexed from 1
        @return TRUE if parameters are changed, FALSE otherwise (2015-10-20)
    */
    virtual bool getVariables(double *variables);

    int getNDim();


    // impose restrictions on the omega values if user inputs the parameters
    // omega1 is resticted to < 1 for both M1a and M2a models
    // omega2 is resticted to 1.0 for both M1a and M2a models
    // omega3 is resticted to > 1 for M2a model
    void restrict_omega_values(string cmix_type);

    double optimizeParameters(double gradient_epsilon);
    
    /**
        write information
        @param out output stream
    */
    virtual void writeInfo(ostream &out);
};

#endif /* modelcodonmixture_h */
