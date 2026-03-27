//
//  phylotreebrmodel.cpp
//  tree
//
//  Created by Thomas Wong on 29/1/25.
//

#include "phylotreebrmodel.h"

/**
 * default constructor
 */
PhyloTreeBranchModel::PhyloTreeBranchModel() : IQTree() {
    model_initialized = false;
    br_models = NULL;
    orig_sub_model = NULL;
}

/**
 * Constructor with given alignment
 * @param alignment
 */
PhyloTreeBranchModel::PhyloTreeBranchModel(Alignment *aln) : IQTree(aln) {
    model_initialized = false;
    br_models = NULL;
    orig_sub_model = NULL;
}

/**
 * default destructor
 */
PhyloTreeBranchModel::~PhyloTreeBranchModel() {
    // return the original site rate model to each model factory
    int i;
    for (i = 0; i < model_facts.size(); i++) {
        model_facts[i]->site_rate = orig_site_rate_models[i];
    }
    getModelFactory()->model = orig_sub_model;
    
    // delete all model factories
    for (i = 0; i < model_facts.size(); i++) {
        delete(model_facts[i]);
    }
    model_facts.clear();
}

void PhyloTreeBranchModel::initializeModel(Params &params, string model_name, ModelsBlock *models_block) {
    
    // error if this function has been called more than one time
    ASSERT(model_initialized == false);
    model_initialized = true;
    
    if (aln->ordered_pattern.empty())
        aln->orderPatternByNumChars(PAT_VARIANT);

    ModelBranch* model_branch = new ModelBranch(this);
    br_models = model_branch;

    int nBranchModels = numBranchModels();
    ASSERT(nBranchModels > 0);
//    cout << endl;
//    cout << "Number of branch models: " << nBranchModels << endl;

    // optain the shared rate model
    string shared_rate = "";
    vector<string> modelparams;
    size_t p1 = model_name.find("BR{");
    size_t p2;
    if (p1 != string::npos) {
        p2 = model_name.find_last_of("}");
        ASSERT(p2 != string::npos && p2 > p1 + 3);
        shared_rate = model_name.substr(p2+1);
        model_name = model_name.substr(0, p2+1);
    } else {
        // a branch model with single class (and root freq)
        p2 = model_name.find_last_of("+");
        if (p2 != string::npos && p2 < model_name.length()-1) {
            // check whether it is rate model
            if (model_name[p2+1] == 'R' || model_name[p2+1] == 'G' || model_name[p2+1] == 'I' || model_name[p2+1] == 'H') {
                shared_rate = model_name.substr(p2);
                model_name = model_name.substr(0, p2);
            }
        }
    }
    
    // for model_joint if set
    string orig_model_joint = "";
    if (!params.model_joint.empty()) {
        orig_model_joint = params.model_joint;
        model_name = params.model_joint;
        params.model_joint.clear();
    }
    
    // remove the "BR{" and "}" from model_name
    p1 = model_name.find("BR{");
    if (p1 != string::npos) {
        p2 = model_name.find_last_of("}");
        ASSERT(p2 != string::npos && p2 > p1 + 3);
        model_name = model_name.substr(p1+3, p2-p1-3);
    }
    
    // create a set of models by splitting on commas that are NOT inside {...} blocks
    // (commas inside F{...} frequency specs must not be treated as model separators)
    {
        size_t fr_pos = 0;
        int depth = 0;
        for (size_t i = 0; i <= model_name.size(); i++) {
            char c = (i < model_name.size()) ? model_name[i] : '\0';
            if (c == '{') {
                ++depth;
            } else if (c == '}') {
                if (depth > 0) --depth;
            } else if ((c == ',' && depth == 0) || c == '\0') {
                modelparams.push_back(model_name.substr(fr_pos, i - fr_pos) + shared_rate);
                fr_pos = i + 1;
            }
        }
    }
    ASSERT(modelparams.size() > 0);
    
    if (modelparams.size() != nBranchModels) {
        outError("The number of models on the trees (i.e. " + convertIntToString(nBranchModels) + ") does not match with the number of classes specified in the branch model (i.e. " + convertIntToString(modelparams.size()) + ")");
    }
    
    // a dummy model with the root frequencies and the shared site rate
    ModelFactory *dm = new ModelFactory(params, modelparams[0], this, models_block);
    orig_sub_model = dm->model;
    dm->model = model_branch;
    dm->site_rate->setTree(this);
    setModelFactory(dm);
    setModel(model_branch);
    setRate(dm->site_rate);
    model_branch->init(FREQ_ESTIMATE);

    // load the models
    for (int i = 0; i < nBranchModels; i++) {
        string curr_model = model_name;
        if (i < modelparams.size() && modelparams[i] != "") {
            // model with user input parameters
            curr_model = modelparams[i];
        }
//        if (i == 0) {
//            cout << "Base model: " << curr_model << endl;
//        } else {
//            cout << "Model " << i << "   : " << curr_model << endl;
//        }
        ModelFactory *mf = new ModelFactory(params, curr_model, this, models_block);
        model_facts.push_back(mf);
        model_branch->push_back((ModelMarkov*)mf->model);
        orig_site_rate_models.push_back(mf->site_rate);
        // replace the site rate by the site rate which is shared among all the model factories
        mf->site_rate = getRate();
    }
    
//    if (shared_rate.length() > 0) {
//        cout << "Rate shared among them: " << shared_rate << endl;
//    }
    
    // set the root frequences
    if (params.root_freq_str != "") {
        model_branch->setRootFrequency(params.root_freq_str);
    }
    
    // set checkpoint
    getModelFactory()->setCheckpoint(checkpoint);
    for (int i = 0; i < nBranchModels; i++) {
        model_facts[i]->setCheckpoint(checkpoint);
    }
    
    // set back the params.model_joint to the original value
    if (!orig_model_joint.empty())
        params.model_joint = orig_model_joint;

}

/**
 *      return the associated substitution model
 */
ModelSubst* PhyloTreeBranchModel::getModel(int branchmodel_id) {
    ASSERT (branchmodel_id < br_models->size());
    return br_models->at(branchmodel_id);
}

ModelFactory* PhyloTreeBranchModel::getModelFactory(int branchmodel_id) {
    ASSERT (branchmodel_id < br_models->size());
    return model_facts[branchmodel_id];
}

/**
    obtain the root frequency vector
*/
void PhyloTreeBranchModel::getRootFrequency(double *state_freq) {
    br_models->getRootFrequency(state_freq);
}

/**
    @return number of branch models, default: 1
*/
int PhyloTreeBranchModel::getNumBrModel() {
    return br_models->size();
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
            ASSERT(nodeid == leafnode->id);
            // get the first neighbor
            Neighbor* first_nei = leafnode->neighbors[0];
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
    int nmodels = getNumBrModel();
    for (int modelid = 0; modelid < nmodels; modelid++) {
        int s = modelid * (aln->STATE_UNKNOWN+1) * nstates * nmixtures;
        for (state = 0; state <= aln->STATE_UNKNOWN; state++) {
            double *state_partial_lh = &tip_partial_lh[state*nstates*nmixtures + s];
            getModel(modelid)->computeTipLikelihood(state, state_partial_lh);
            if (getModel()->useRevKernel()) {
                ASSERT(getModel(modelid)->useRevKernel());
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

void PhyloTreeBranchModel::computePtnInvar() {
    // store the current pointers of model and model factory
    ModelSubst* modelorig = model;
    ModelFactory* modelfactorig = model_factory;
    // set all pointers to the base model and model factory
    model = br_models->at(0);
    model_factory = model_facts[0];
    // compute the invar of the patterns according to the base model4
    // TODO: may need to update in the future
    PhyloTree::computePtnInvar();
    // restore the original values for the pointers
    model = modelorig;
    model_factory = modelfactorig;
}

string PhyloTreeBranchModel::optimizeModelParameters(bool printInfo, double logl_epsilon) {
    br_models->logl_epsilon = logl_epsilon;
    if (logl_epsilon == -1)
        br_models->logl_epsilon = params->modelEps;
    return IQTree::optimizeModelParameters(printInfo, logl_epsilon);
}

int PhyloTreeBranchModel::getNParameters() {
    int df = 0;
    
    if (verbose_mode >= VB_MED)
        cout << endl << "Number of parameters:" << endl;
    
    int nmodels = getNumBrModel();
    for (int modelid = 0; modelid < nmodels; modelid++) {
        df += getModel(modelid)->getNDim() + getModel(modelid)->getNDimFreq();
    }
    df += model->getNDimFreq(); // root frequency
    df += site_rate->getNDim(); // site model
    if (params->fixed_branch_length != BRLEN_FIX) {
        df += site_rate->getTree()->getNBranchParameters(BRLEN_OPTIMIZE);
    }
    
    return df;
}
