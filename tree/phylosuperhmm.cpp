//
//  phylosuperhmm.cpp
//
//  Created by Thomas Wong on 5/12/24.
//

#include "phylosuperhmm.h"

// get the trees inside the infile
void getTreesInFile(char* infile, vector<string>& treestrs) {
    ifstream fin;
    string aline;
    string tree_str = "";
    size_t i;
    fin.open(infile);
    while (getline(fin,aline)) {
        for (i=0; i<aline.length(); i++) {
            if (aline[i] >= '!' && aline[i] <= '~') {
                tree_str.append(1, aline[i]);
                if (aline[i] == ';') {
                    treestrs.push_back(tree_str);
                    tree_str = "";
                }
            }
        }
    }
    fin.close();
    if (tree_str.length() > 0)
        treestrs.push_back(tree_str);
}

// get number of trees inside the infile
int getNumTreesInFile(char* infile) {
    ifstream fin;
    string aline;
    size_t i;
    int ntree = 0;
    fin.open(infile);
    while (getline(fin,aline)) {
        for (i=0; i<aline.length(); i++) {
            if (aline[i] == ';') {
                ntree++;
            }
        }
    }
    fin.close();
    return ntree;
}

/**
    constructor
*/
PhyloSuperHmm::PhyloSuperHmm() : PhyloSuperTree() {
}

/**
 constructor
 */
PhyloSuperHmm::PhyloSuperHmm(SuperAlignment *alignment, Params &params) : PhyloSuperTree(alignment,false,false) {
    ASSERT(strlen(params.user_file) > 0);
    int ntree = getNumTreesInFile(params.user_file);
    int npart = alignment->partitions.size();
    size_t i, j;
    
    // Create a set of IQTreeMixHmm
    for (i = 0; i < npart; i++) {
        push_back(new IQTreeMixHmm(params, alignment->partitions[i]));
    }

    for (i = 0; i < ntree; i++) {
        PhyloSuperTree* stree = new PhyloSuperTree(alignment, false, false);
        for (j = 0; j < npart; j++) {
            IQTreeMixHmm* mtreeHmm = (IQTreeMixHmm*) at(j);
            stree->push_back(mtreeHmm->at(i));
        }
        superTreeSet.push_back(stree);
    }
}

/**
 destructor
 */
PhyloSuperHmm::~PhyloSuperHmm() {
    model_factory = nullptr;
    model = nullptr;
    site_rate = nullptr;
    for (reverse_iterator it = rbegin(); it != rend(); it++) {
        delete (*it);
    }
    clear();
}

/**
    set minimum branch length
*/
void PhyloSuperHmm::setMinBranchLen(Params& params) {
    
    if (params.min_branch_length <= 0.0) {
        params.min_branch_length = MAST_MIN_BRANCH_LEN;
    }
    cout << setprecision(7) << "Minimum branch length is set to " << params.min_branch_length << endl;
}

void PhyloSuperHmm::initSettings(Params &params) {
    IQTree::initSettings(params);
    // setLikelihoodKernel(params.SSE);
    // setNumThreads(params.num_threads);
    int ntree = superTreeSet.size();
    int npart = size();
    size_t i, j;
    for (i = 0; i < ntree; i++) {
        superTreeSet[i]->IQTree::initSettings(params);
    }
    for (i = 0; i < npart; i++) {
        IQTreeMixHmm* treeMixHmm = (IQTreeMixHmm*) at(i);
        treeMixHmm->initSettings(params);
    }
}

void PhyloSuperHmm::initializeModel(Params &params, string model_name, ModelsBlock *models_block) {
    if (params.treemix_optimize_methods != "hmm") {
        outError("Tree mixture model on partitioned alignment only support HMM model. Use -hmmster option");
    }
    int npart = size();
    size_t i;
    IQTree::initializeModel(params, model_name, models_block);
    if (model_name == "")
        model_name = params.model_name;
    for (i = 0; i < npart; i++) {
        IQTreeMixHmm* treeMixHmm = (IQTreeMixHmm*) at(i);
        treeMixHmm->initializeModel(params, model_name, models_block);
    }
}

/**
 * Generate the initial tree (usually used for model parameter estimation)
 */
void PhyloSuperHmm::computeInitialTree(LikelihoodKernel kernel, istream* in) {
    if(!params->user_file) {
        outError("HMM model on partitioned alignment must have input trees");
    }
    ifstream fin;
    fin.open(params->user_file);
    for (int i = 0; i < superTreeSet.size(); i++) {
        superTreeSet[i]->computeInitialTree(kernel, &fin);
        superTreeSet[i]->mapTrees();
    }
    fin.close();
    // compute the initial tree individually.
    int npart = size();
    for (int i = 0; i < npart; i++) {
        stringstream ss;
        IQTreeMixHmm* tmixhmm = (IQTreeMixHmm*) at(i);
        for (int j = 0; j < tmixhmm->size(); j++)
            tmixhmm->at(j)->printTree(ss, WT_NEWLINE);
        ss.seekg(ios_base::beg);
        tmixhmm->computeInitialTree(kernel, &ss);
    }
}

void PhyloSuperHmm::setRootNode(const char *my_root, bool multi_taxa) {
    for (int i = 0; i < superTreeSet.size(); i++)
        superTreeSet[i]->setRootNode(my_root, multi_taxa);
}

// show the assignment of the categories along sites with max likelihood
// cat_assign_method:
//  0 - the categories along sites is assigned according to the path with maximum probability (default)
//  1 - the categories along sites is assigned according to the max posterior probability
void PhyloSuperHmm::printResults(string prefix, string ext, int cat_assign_method) {
    if (size() == 0)
        return;
    
    iterator it = begin();
    int ntree = ((IQTreeMixHmm*)(*it))->size();
    int* numSiteCat = new int[ntree];
    int* numSiteCatTot = new int[ntree];
    int numSiteSum = 0;
    int k;
    memset(numSiteCatTot, 0, sizeof(int)*ntree);
    for (iterator it = begin(); it != end(); it++) {
        string filename = prefix + "." + (*it)->aln->name + ext;
        ((IQTreeMixHmm*)(*it))->printResults(filename.c_str(), cat_assign_method, numSiteCat);
        for (k = 0; k < ntree; k++)
            numSiteCatTot[k] += numSiteCat[k];
    }
    for (k = 0; k < ntree; k++)
        numSiteSum += numSiteCatTot[k];
    
    // print out overall percentage of sites over the trees
    cout << "Overall percentage of sites over the trees:";
    for (k = 0; k < ntree; k++)
        cout << " " << (double) numSiteCatTot[k] / numSiteSum;
    cout << endl;
    
    delete[] numSiteCat;
    delete[] numSiteCatTot;
}

// get number of trees
int PhyloSuperHmm::getNumTrees() {
    return superTreeSet.size();
}

/**
    set checkpoint object
    @param checkpoint
*/
void PhyloSuperHmm::setCheckpoint(Checkpoint *checkpoint) {
    size_t i;
    int ntree = superTreeSet.size();
    int npart = size();
    IQTree::setCheckpoint(checkpoint);
    for (i=0; i<ntree; i++) {
        superTreeSet[i]->IQTree::setCheckpoint(checkpoint);
    }
    for (i=0; i<npart; i++) {
        IQTreeMixHmm* treemixhmm = (IQTreeMixHmm*) at(i);
        treemixhmm->setCheckpoint(checkpoint);
    }
}

void PhyloSuperHmm::startCheckpoint() {
    checkpoint->startStruct("PhyloSuperHmm" + convertIntToString(getNumTrees()));
}

void PhyloSuperHmm::saveCheckpoint() {
    size_t i,j;
    int ntree = superTreeSet.size();
    int npart = size();
    startCheckpoint();
    // save the trees
    for (i=0; i<ntree; i++) {
        checkpoint->startStruct("Tree" + convertIntToString(i+1));
        superTreeSet[i]->saveCheckpoint();
        checkpoint->endStruct();
    }
    // save the tree weights
    for (i=0; i<npart; i++) {
        IQTreeMixHmm* treemixhmm = (IQTreeMixHmm*) at(i);
        checkpoint->startStruct("Part" + convertIntToString(i+1));
        CKP_ARRAY_SAVE(ntree, &(treemixhmm->weights[0]));
        checkpoint->endStruct();
    }
    endCheckpoint();
}

void PhyloSuperHmm::restoreCheckpoint() {
    size_t i,j;
    int ntree = superTreeSet.size();
    int npart = size();
    startCheckpoint();
    // load the trees
    for (i=0; i<ntree; i++) {
        checkpoint->startStruct("Tree" + convertIntToString(i+1));
        superTreeSet[i]->restoreCheckpoint();
        checkpoint->endStruct();
    }
    // save the tree weights
    for (i=0; i<npart; i++) {
        IQTreeMixHmm* treemixhmm = (IQTreeMixHmm*) at(i);
        checkpoint->startStruct("Part" + convertIntToString(i+1));
        if (CKP_ARRAY_RESTORE(ntree, &(treemixhmm->weights[0]))) {
            for (j = 0; j < ntree; j++) {
                treemixhmm->weight_logs[j] = log(treemixhmm->weights[j]);
            }
        }
        checkpoint->endStruct();
    }
    endCheckpoint();
    clearAllPartialLH();
}

void PhyloSuperHmm::setParams(Params* params) {
    int ntree = superTreeSet.size();
    int npart = size();
    IQTree::setParams(params);
    size_t i;
    for (i=0; i<ntree; i++)
        superTreeSet[i]->IQTree::setParams(params);
    for (i=0; i<npart; i++)
        at(i)->setParams(params);
}

/**
 * save branch lengths into a vector
 */
void PhyloSuperHmm::saveBranchLengths(DoubleVector &lenvec, int startid, PhyloNode *node, PhyloNode *dad) {
    ASSERT(getMixlen() == 1); // supertree and treemixlen not allowed together
    int totalBranchNum = branchNum * getMixlen();
    int ntree = superTreeSet.size();
    int npart = size();
    size_t i, j;
    for (i = 0; i < ntree; i++) {
        for (j = 0; j < npart; j++) {
            totalBranchNum += superTreeSet[i]->at(j)->branchNum * superTreeSet[i]->at(j)->getMixlen();
        }
    }
    lenvec.resize(startid + totalBranchNum);

    PhyloTree::saveBranchLengths(lenvec, startid);
    startid += branchNum * getMixlen();
    for (i = 0; i < ntree; i++) {
        for (j = 0; j < npart; j++) {
            superTreeSet[i]->at(j)->saveBranchLengths(lenvec, startid);
            startid += superTreeSet[i]->at(j)->branchNum * superTreeSet[i]->at(j)->getMixlen();
        }
    }
}

/**
 * restore branch lengths from a vector previously called with saveBranchLengths
 */
void PhyloSuperHmm::restoreBranchLengths(DoubleVector &lenvec, int startid, PhyloNode *node, PhyloNode *dad) {
    int ntree = superTreeSet.size();
    int npart = size();
    size_t i, j;
    PhyloTree::restoreBranchLengths(lenvec, startid);
    startid += branchNum * getMixlen();
    for (i = 0; i < ntree; i++) {
        for (j = 0; j < npart; j++) {
            superTreeSet[i]->at(j)->restoreBranchLengths(lenvec, startid);
            startid += superTreeSet[i]->at(j)->branchNum * superTreeSet[i]->at(j)->getMixlen();
        }
    }
}

//ModelFactory* PhyloSuperHmm::getModelFactory() {
//    return ((IQTreeMixHmm*)at(0))->getModelFactory();
//}
//
//ModelSubst* PhyloSuperHmm::getModel() {
//    return ((IQTreeMixHmm*)at(0))->getModel();
//}
//
//RateHeterogeneity* PhyloSuperHmm::getRate() {
//    return ((IQTreeMixHmm*)at(0))->getRate();
//}

/**
    test the best number of threads
*/
int PhyloSuperHmm::testNumThreads() {
    int bestNThres = size();
    setNumThreads(bestNThres);
    return bestNThres;
}

void PhyloSuperHmm::saveModelCheckpoint() {
    // save the models
    size_t i;
    int npart = size();
    startCheckpoint();
    for (i=0; i<npart; i++) {
        IQTreeMixHmm* treemixhmm = (IQTreeMixHmm*) at(i);
        checkpoint->startStruct("Part" + convertIntToString(i+1));
        treemixhmm->saveModelCheckpoint();
        checkpoint->endStruct();
    }
    endCheckpoint();
}

void PhyloSuperHmm::restoreModelCheckpoint() {
    // restore the models
    size_t i;
    int npart = size();
    startCheckpoint();
    for (i=0; i<npart; i++) {
        IQTreeMixHmm* treemixhmm = (IQTreeMixHmm*) at(i);
        checkpoint->startStruct("Part" + convertIntToString(i+1));
        treemixhmm->restoreModelCheckpoint();
        checkpoint->endStruct();
    }
    endCheckpoint();
}

/**
    compute the weighted average of branch lengths over the weighted tree of partitions
*/
void PhyloSuperHmm::computeBranchLengths() {
    if (verbose_mode >= VB_DEBUG)
        cout << "Assigning branch lengths for full tree with weighted average..." << endl;
    
    int part = 0, i, k;
    int ntree = superTreeSet.size();
    
    for (k = 0; k < ntree; k++) {
        iterator it;
        NodeVector nodes1, nodes2;
        superTreeSet[k]->getBranches(nodes1, nodes2);
        vector<SuperNeighbor*> neighbors1;
        vector<SuperNeighbor*> neighbors2;
        IntVector occurence;
        occurence.resize(nodes1.size(), 0);
        for (i = 0; i < nodes1.size(); i++) {
            neighbors1.push_back((SuperNeighbor*)nodes1[i]->findNeighbor(nodes2[i]) );
            neighbors2.push_back((SuperNeighbor*)nodes2[i]->findNeighbor(nodes1[i]) );
            neighbors1.back()->length = 0.0;
        }
        for (it = superTreeSet[k]->begin(), part = 0; it != superTreeSet[k]->end(); it++, part++) {
            IntVector brfreq;
            brfreq.resize((*it)->branchNum, 0);
            for (i = 0; i < nodes1.size(); i++) {
                PhyloNeighbor *nei1 = neighbors1[i]->link_neighbors[part];
                if (!nei1) continue;
                brfreq[nei1->id]++;
            }
            double treeweight = ((IQTreeMix*)at(part))->weights[k];
            double inv_weight = 1.0 / treeweight;
            for (i = 0; i < nodes1.size(); i++) {
                PhyloNeighbor *nei1 = neighbors1[i]->link_neighbors[part];
                if (!nei1) continue;
                if ((*it)->aln->seq_type == SEQ_CODON && rescale_codon_brlen) {
                    // rescale branch length by 3
                    neighbors1[i]->length += (nei1->length) * (*it)->aln->getNSite() / brfreq[nei1->id] * treeweight;
                    occurence[i] += (*it)->aln->getNSite()*3*treeweight;
                } else {
                    neighbors1[i]->length += (nei1->length) * (*it)->aln->getNSite() / brfreq[nei1->id] * treeweight;
                    occurence[i] += (*it)->aln->getNSite()*treeweight;
                }
                //cout << neighbors1[i]->id << "  " << nodes1[i]->id << nodes1[i]->name <<"," << nodes2[i]->id << nodes2[i]->name <<": " << (nei1->length) / brfreq[nei1->id] << endl;
            }
            //cout << endl;
        }
        for (i = 0; i < nodes1.size(); i++) {
            if (occurence[i])
                neighbors1[i]->length /= occurence[i];
            neighbors2[i]->length = neighbors1[i]->length;
        }
    }
}



string PhyloSuperHmm::getTreeString() {
    stringstream tree_stream;
    size_t i;
    int ntree = superTreeSet.size();
    for (i=0; i<ntree; i++) {
        tree_stream << superTreeSet[i]->getTreeString() << endl;
    }
    return tree_stream.str();
}

void PhyloSuperHmm::printResultTree(string suffix) {
    if (MPIHelper::getInstance().isWorker()) {
        return;
    }
    if (params->suppress_output_flags & OUT_TREEFILE)
        return;
    
    int ntree = superTreeSet.size();
    int npart = size();
    size_t i, j;

    // print out the main trees
    string main_tree_file_name = params->out_prefix;
    main_tree_file_name += ".treefile";
    if (suffix.compare("") != 0) {
        main_tree_file_name += "." + suffix;
    }
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(main_tree_file_name);
        for (i = 0; i < ntree; i++) {
            superTreeSet[i]->printTree(out, WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
        }
        out.close();
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, main_tree_file_name);
    }

    // print out the part trees
    string part_tree_file_name = params->out_prefix;
    part_tree_file_name += ".parttrees";
    if (suffix.compare("") != 0) {
        part_tree_file_name += "." + suffix;
    }
    try {
        ofstream out;
        out.exceptions(ios::failbit | ios::badbit);
        out.open(part_tree_file_name);
        for (i = 0; i < npart; i++) {
            out << "#part:" << at(i)->aln->name << endl;
            IQTreeMix* treemix = (IQTreeMix*) at(i);
            treemix->setRootNode(params->root, true);
            for (j = 0; j < ntree; j++) {
                treemix->at(j)->printTree(out, WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
            }
        }
        out.close();
    } catch (ios::failure) {
        outError(ERR_WRITE_OUTPUT, part_tree_file_name);
    }
    if (verbose_mode >= VB_MED)
        cout << "Partition trees printed to " << part_tree_file_name << endl;
}
