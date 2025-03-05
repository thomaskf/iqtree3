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
    cout << "[PhyloTreeHal::initializeModel] model_name = " << model_name << endl;
    
    int numHALs = numHALModels();
    if (numHALs > 0) {
        cout << "  Number of HAL models: " << numHALs << endl;
        
        // obtain the user input model parameters if exists
        vector<string> modelparams;
        getUserInputModelParams(modelparams);
        
        models.push_back(getModelFactory()->model);
        for (int i = 1; i < numHALs; i++) {
            string curr_model = model_name;
            if (modelparams.size() > i && modelparams[i] != "") {
                // model with user input parameters
                curr_model = modelparams[i];
                cout << "modelparams[" << i << "]=" << curr_model << endl;
            }
            ModelFactory* mf = new ModelFactory(params, curr_model, this, models_block);
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

/*
 * obtain the user input model parameters
 */
void PhyloTreeHal::getUserInputModelParams(vector<string> &modelparams, Node *node, Node *dad) {
    
    if (node == NULL)
        node = root;
    
    FOR_NEIGHBOR_IT(node, dad, it) {
        if (!(*it)->attributes.empty()) {
            // get the input branch_id or clade_id
            int branch_id = (*it)->hal_id;
            int clade_id = -1;
            for (auto attr : (*it)->attributes) {
                if (attr.first == "clade") {
                    clade_id = atoi(attr.second.c_str());
                }
            }
            // get the user input model parameters
            for (auto attr : (*it)->attributes) {
                int curr_id = -1;
                if (attr.first == "branch_model" && branch_id >= 0) {
                    curr_id = branch_id;
                } else if (attr.first == "clade_model" && clade_id >= 0) {
                    curr_id = clade_id;
                }
                if (curr_id >= 0) {
                    if (modelparams.size() < curr_id+1) {
                        modelparams.resize(curr_id+1);
                    }
                    modelparams[curr_id] = attr.second;
                }
            }
        }
        getUserInputModelParams(modelparams, (*it)->node, node);
    }
}
