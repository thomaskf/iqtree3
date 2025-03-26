//
//  phylotreebrmodel.h
//  tree
//
//  Created by Thomas Wong on 28/01/25.
//

#ifndef phylotreebrmodel_h
#define phylotreebrmodel_h

#include "iqtree.h"

class PhyloTreeBranchModel : public IQTree
{
public:
    
    /**
     * default constructor
     */
    PhyloTreeBranchModel() : IQTree() {}

    /**
     * Constructor with given alignment
     * @param alignment
     */
    PhyloTreeBranchModel(Alignment *aln) : IQTree(aln) {}
    
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
     * check how many different branch models
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

    /**
     models
     */
    vector<ModelSubst*> models;
    vector<ModelFactory*> modelFacts;


};

#endif /* phylotreebrmodel_h */
