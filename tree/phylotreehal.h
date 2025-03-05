//
//  phylotreehal.h
//  tree
//
//  Created by Thomas Wong on 28/01/25.
//

#ifndef phylotreehal_h
#define phylotreehal_h

#include "iqtree.h"

class PhyloTreeHal : public IQTree
{
public:
    
    /**
     * default constructor
     */
    PhyloTreeHal() : IQTree() {}

    /**
     * Constructor with given alignment
     * @param alignment
     */
    PhyloTreeHal(Alignment *aln) : IQTree(aln) {}
    
    /**
     * default destructor
     */
    ~PhyloTreeHal();

    virtual void initializeModel(Params &params, string model_name, ModelsBlock *models_block);

    /**
     *      return the associated substitution model
     */
    virtual ModelSubst* getModel(int hal_id);

    /*
     * check how many different HAL models
     */
    int numHALModels(Node *node = NULL, Node *dad = NULL);
    
    /*
     * obtain the user input model parameters
     */
    void getUserInputModelParams(vector<string> &modelparams, Node *node = NULL, Node *dad = NULL);

    /**
        @return true as this is a HAL model (i.e. heterogeneity across sites model)
     */
    virtual bool isHAL() { return true; }

    /**
     models
     */
    vector<ModelSubst*> models;
    vector<ModelFactory*> modelFacts;


};

#endif /* phylotreehal_h */
