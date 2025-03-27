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
    
    int nBranchModels = numBranchModels();
    if (nBranchModels > 0) {
        cout << endl;
        cout << "Number of branch models: " << nBranchModels << endl;
        
        // obtain the user input model parameters if exists
        vector<string> modelparams;
        getUserInputModelParams(modelparams);
        
        models.push_back(IQTree::getModelFactory()->model);
        cout << "Base model: " << model_name << endl;
        for (int i = 1; i < nBranchModels; i++) {
            string curr_model = model_name;
            if (modelparams.size() > i && modelparams[i] != "") {
                // model with user input parameters
                curr_model = modelparams[i];
                cout << "Model " << i << "   : " << curr_model << endl;
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
    if (branchmodel_id >= models.size())
        return model;
    
    return models[branchmodel_id];
}

ModelFactory* PhyloTreeBranchModel::getModelFactory(int branchmodel_id) {
    if (branchmodel_id >= modelFacts.size())
        return model_factory;
    
    return modelFacts[branchmodel_id];
}

/**
    @return number of branch models, default: 1
*/
int PhyloTreeBranchModel::getNumBrModel() {
    return models.size();
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
            int clade_id = 0;
            for (auto attr : (*it)->attributes) {
                if (attr.first == "clade") {
                    clade_id = atoi(attr.second.c_str());
                }
            }
            // get the user input model parameters
            for (auto attr : (*it)->attributes) {
                int curr_id = 0;
                if (attr.first == "branch_model" && branch_id > 0) {
                    curr_id = branch_id;
                } else if (attr.first == "clade_model" && clade_id > 0) {
                    curr_id = clade_id;
                }
                if (curr_id > 0) {
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

void PhyloTreeBranchModel::computeTipPartialLikelihood() {
    if ((tip_partial_lh_computed & 1) != 0)
        return;
    tip_partial_lh_computed |= 1;
    
    
    //-------------------------------------------------------
    // initialize ptn_freq and ptn_invar
    //-------------------------------------------------------

    computePtnFreq();
    // for +I model
    computePtnInvar();

    if (getModel()->isSiteSpecificModel()) {
        // TODO: THIS NEEDS TO BE CHANGED TO USE ModelSubst::computeTipLikelihood()
//        ModelSet *models = (ModelSet*)model;
        size_t nptn = aln->getNPattern(), max_nptn = ((nptn+vector_size-1)/vector_size)*vector_size, tip_block_size = max_nptn * aln->num_states;
        int nstates = aln->num_states;
        size_t nseq = aln->getNSeq();
        ASSERT(vector_size > 0);
        
        
#ifdef _OPENMP
        #pragma omp parallel for schedule(static)
#endif
        for (int nodeid = 0; nodeid < nseq; nodeid++) {
            auto stateRow = getConvertedSequenceByNumber(nodeid);
            Node* leafnode = findNodeName(aln->getSeqName(nodeid));
            cout << "nodeid = " << nodeid << "; leafnode->id = " << leafnode->id << endl;
            ASSERT(nodeid == leafnode->id);
            // get the first neighbor
            Neighbor* first_nei = leafnode->neighbors[0];
            cout << "first neighbor's branch model: " << first_nei->branchmodel_id << endl;
            ModelSubst* curr_model = getModel(first_nei->branchmodel_id);
            double *partial_lh = tip_partial_lh + tip_block_size*nodeid;
            for (size_t ptn = 0; ptn < nptn; ptn+=vector_size, partial_lh += nstates*vector_size) {
                double *inv_evec = &curr_model->getInverseEigenvectors()[ptn*nstates*nstates];
                for (int v = 0; v < vector_size; v++) {
                    int state = 0;
                    if (ptn+v < nptn) {
                        if (stateRow!=nullptr) {
                            state = stateRow[ptn+v];
                        } else {
                            state = aln->at(ptn+v)[nodeid];
                        }
                    }
                    if (state < nstates) {
                        for (int i = 0; i < nstates; i++)
                            partial_lh[i*vector_size+v] = inv_evec[(i*nstates+state)*vector_size+v];
                    } else if (state == aln->STATE_UNKNOWN) {
                        // special treatment for unknown char
                        for (int i = 0; i < nstates; i++) {
                            double lh_unknown = 0.0;
                            for (int x = 0; x < nstates; x++) {
                                lh_unknown += inv_evec[(i*nstates+x)*vector_size+v];
                            }
                            partial_lh[i*vector_size+v] = lh_unknown;
                        }
                    } else {
                        double lh_ambiguous;
                        // ambiguous characters
                        int ambi_aa[] = {
                            4+8, // B = N or D
                            32+64, // Z = Q or E
                            512+1024 // U = I or L
                            };
                        switch (aln->seq_type) {
                        case SEQ_DNA:
                            {
                                int cstate = state-nstates+1;
                                for (int i = 0; i < nstates; i++) {
                                    lh_ambiguous = 0.0;
                                    for (int x = 0; x < nstates; x++)
                                        if ((cstate) & (1 << x))
                                            lh_ambiguous += inv_evec[(i*nstates+x)*vector_size+v];
                                    partial_lh[i*vector_size+v] = lh_ambiguous;
                                }
                            }
                            break;
                        case SEQ_PROTEIN:
                            //map[(unsigned char)'B'] = 4+8+19; // N or D
                            //map[(unsigned char)'Z'] = 32+64+19; // Q or E
                            {
                                int cstate = state-nstates;
                                for (int i = 0; i < nstates; i++) {
                                    lh_ambiguous = 0.0;
                                    for (int x = 0; x < 11; x++)
                                        if (ambi_aa[cstate] & (1 << x))
                                            lh_ambiguous += inv_evec[(i*nstates+x)*vector_size+v];
                                    partial_lh[i*vector_size+v] = lh_ambiguous;
                                }
                            }
                            break;
                        default:
                            ASSERT(0);
                            break;
                        }
                    }
                    // sanity check
    //                bool all_zero = true;
    //                for (i = 0; i < nstates; i++)
    //                    if (partial_lh[i] != 0) {
    //                        all_zero = false;
    //                        break;
    //                    }
    //                assert(!all_zero && "some tip_partial_lh are all zeros");
                    
                } // FOR v
            } // FOR ptn
            // NO Need to copy dummy anymore
            // dummy values
//            for (ptn = nptn; ptn < max_nptn; ptn++, partial_lh += nstates)
//                memcpy(partial_lh, partial_lh-nstates, nstates*sizeof(double));
        } // FOR nodeid
        return;
    }
    
    // 2020-06-23: refactor to use computeTipLikelihood
    int nmixtures = 1;
    if (getModel()->useRevKernel())
        nmixtures = getModel()->getNMixtures();
    int nstates = getModel()->num_states;
    int state;
    if (aln->seq_type == SEQ_POMO) {
        if (aln->pomo_sampling_method != SAMPLING_WEIGHTED_BINOM &&
            aln->pomo_sampling_method != SAMPLING_WEIGHTED_HYPER)
            outError("Sampling method not supported by PoMo.");
        ASSERT(aln->STATE_UNKNOWN == nstates + aln->pomo_sampled_states.size());
    }

    // assign tip_partial_lh for all admissible states
    int nmodels = models.size();
    for (int modelid = 0; modelid < nmodels; modelid++) {
        int s = modelid * (aln->STATE_UNKNOWN+1) * nstates * nmixtures;
        for (state = 0; state <= aln->STATE_UNKNOWN; state++) {
            double *state_partial_lh = &tip_partial_lh[state*nstates*nmixtures + s];
            getModel(modelid)->computeTipLikelihood(state, state_partial_lh);
            if (getModel(modelid)->useRevKernel()) {
                // transform to inner product of tip likelihood and inverse-eigenvector
                getModel(modelid)->multiplyWithInvEigenvector(state_partial_lh);
            }
        }
    }
}

double PhyloTreeBranchModel::computeLikelihood(double *pattern_lh, bool save_log_value) {
    ASSERT(root->isLeaf());
    if (!current_it) {
        current_it = (PhyloNeighbor*)root->neighbors[0];
        current_it_back = (PhyloNeighbor*)current_it->node->findNeighbor(root);
    }
    return PhyloTree::computeLikelihood(pattern_lh, save_log_value);
}
