//
//  phylotreebrmodel.h
//  tree
//
//  Created by Thomas Wong on 28/01/25.
//

#ifndef phylotreebrmodel_h
#define phylotreebrmodel_h

#include "iqtree.h"
#include "model/modelbranch.h"

class PhyloTreeBranchModel : public IQTree
{
public:
    
    /**
     * default constructor
     */
    PhyloTreeBranchModel();

    /**
     * Constructor with given alignment
     * @param alignment
     */
    PhyloTreeBranchModel(Alignment *aln);
    
    /**
     * default destructor
     */
    ~PhyloTreeBranchModel();

    virtual void initializeModel(Params &params, string model_name, ModelsBlock *models_block);

    /**
     *      return the associated substitution model
     */

    virtual ModelSubst *getModel() {
        return model;
    }

    virtual ModelSubst* getModel(int hal_id);

    virtual ModelFactory *getModelFactory() {
        return model_factory;
    }

    virtual ModelFactory* getModelFactory(int hal_id);

    /*
     * check how many different branch models from the tree
     */
    int numBranchModels(Node *node = NULL, Node *dad = NULL);
    
    /*
     * obtain the user input model parameters
     */
    void getUserInputModelParams(vector<string> &modelparams, Node *node = NULL, Node *dad = NULL);

    /**
        @return true as this is a branch model (i.e. each branch has different substitution model)
     */
    virtual bool isBranchModel() { return true; }

    /**
        @return number of branch models, default: 1
    */
    virtual int getNumBrModel();

    virtual void computeTipPartialLikelihood();

    virtual double computeLikelihood(double *pattern_lh = NULL, bool save_log_value = true);
    
    virtual void computePtnInvar();

    /**
     model factories
     */
    vector<ModelFactory*> model_facts;
    
    /**
     corresponding branch model
     */
    ModelBranch* br_models;
    
    // has the model been initialized
    bool model_initialized;
    
private:
    // original site rate model
    vector<RateHeterogeneity*> orig_site_rate_models;
    // original sub model
    ModelSubst* orig_sub_model;
};

#endif /* phylotreebrmodel_h */
