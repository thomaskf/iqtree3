//
//  phylotreehal.cpp
//  tree
//
//  Created by Thomas Wong on 29/1/25.
//

#include "phylotreehal.h"

/**
 * default destructor
 */
PhyloTreeHal::~PhyloTreeHal() {
    models.clear();
    for (int i = 0; i < modelFacts.size(); i++) {
        delete(modelFacts[i]);
    }
    modelFacts.clear();
}

void PhyloTreeHal::initializeModel(Params &params, string model_name, ModelsBlock *models_block) {
    IQTree::initializeModel(params, model_name, models_block);
    
    int numHALs = numHALModels();
    if (numHALs > 0) {
        cout << "Number of HAL models: " << numHALs << endl;
        models.push_back(getModelFactory()->model);
        for (int i = 1; i < numHALs; i++) {
            ModelFactory* mf = new ModelFactory(params, model_name, this, models_block);
            modelFacts.push_back(mf);
            models.push_back(mf->model);
        }
    }
}

/**
 *      return the associated substitution model
 */
ModelSubst* PhyloTreeHal::getModel(int hal_id) {
    if (hal_id == -1 || hal_id >= models.size())
        return model;
    
    return models[hal_id];
}

/*
 * check how many different HAL models
 */
int PhyloTreeHal::numHALModels(Node *node, Node *dad) {
    
    int max = 0;
    if (node == NULL)
        node = root;
    
    FOR_NEIGHBOR_IT(node, dad, it) {
        if ((*it)->hal_id+1 > max)
            max = (*it)->hal_id+1;
        int c = numHALModels((*it)->node, node);
        if (c > max)
            max = c;
    }
    
    return max;
}
