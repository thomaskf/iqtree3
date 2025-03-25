//
//  phylotreebrmodel.cpp
//  tree
//
//  Created by Thomas Wong on 29/1/25.
//

#include "phylotreebrmodel.h"

/**
 * default destructor
 */
PhyloTreeBranchModel::~PhyloTreeBranchModel() {
    models.clear();
    for (int i = 0; i < modelFacts.size(); i++) {
        delete(modelFacts[i]);
    }
    modelFacts.clear();
}

void PhyloTreeBranchModel::initializeModel(Params &params, string model_name, ModelsBlock *models_block) {
    IQTree::initializeModel(params, model_name, models_block);
    cout << "[PhyloTreeBranchModel::initializeModel] model_name = " << model_name << endl;
    
    int nBranchModels = numBranchModels();
    if (nBranchModels > 0) {
        cout << "  Number of Branch models: " << nBranchModels << endl;
        
        // obtain the user input model parameters if exists
        vector<string> modelparams;
        getUserInputModelParams(modelparams);
        
        models.push_back(IQTree::getModelFactory()->model);
        for (int i = 1; i < nBranchModels; i++) {
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
ModelSubst* PhyloTreeBranchModel::getModel(int branchmodel_id) {
    if (branchmodel_id == -1 || branchmodel_id >= models.size())
        return model;
    
    return models[branchmodel_id];
}

ModelFactory* PhyloTreeBranchModel::getModelFactory(int branchmodel_id) {
    if (branchmodel_id == -1 || branchmodel_id >= modelFacts.size())
        return model_factory;
    
    return modelFacts[branchmodel_id];
}

/*
 * check how many different branch models
 */
int PhyloTreeBranchModel::numBranchModels(Node *node, Node *dad) {
    
    int max = 0;
    if (node == NULL)
        node = root;
    
    FOR_NEIGHBOR_IT(node, dad, it) {
        if ((*it)->branchmodel_id+1 > max)
            max = (*it)->branchmodel_id+1;
        int c = numBranchModels((*it)->node, node);
        if (c > max)
            max = c;
    }
    
    return max;
}

/*
 * obtain the user input model parameters
 */
void PhyloTreeBranchModel::getUserInputModelParams(vector<string> &modelparams, Node *node, Node *dad) {
    
    if (node == NULL)
        node = root;
    
    FOR_NEIGHBOR_IT(node, dad, it) {
        if (!(*it)->attributes.empty()) {
            // get the input branch_id or clade_id
            int branch_id = (*it)->branchmodel_id;
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
