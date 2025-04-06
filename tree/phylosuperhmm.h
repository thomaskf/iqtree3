//
//  phylosuperhmm.h
//
//  Created by Thomas Wong on 5/12/24.
//

#ifndef phylosuperhmm_h
#define phylosuperhmm_h

#include <stdio.h>
#include "phylosupertree.h"
#include "iqtreemixhmm.h"

class PhyloSuperHmm : public PhyloSuperTree
{
public:
    
    /**
     constructor
     */
    PhyloSuperHmm();
    
    /**
     constructor
     */
    PhyloSuperHmm(SuperAlignment *alignment, Params &params);
    
    /**
     destructor
     */
    ~PhyloSuperHmm();
    
    /**
     @return true if this is a mixture of trees, default: false
     */
    virtual bool isTreeMix() { return true; }
    
    /**
     set minimum branch length
     */
    virtual void setMinBranchLen(Params& params);
    
    virtual void initSettings(Params &params);
    
    virtual void initializeModel(Params &params, string model_name, ModelsBlock *models_block);
    
    /**
     * Generate the initial tree (usually used for model parameter estimation)
     */
    virtual void computeInitialTree(LikelihoodKernel kernel, istream* in);
    
    virtual void setRootNode(const char *my_root, bool multi_taxa);

    virtual void setParams(Params* params);

    // show the assignment of the categories along sites with max likelihood
    // cat_assign_method:
    //  0 - the categories along sites is assigned according to the path with maximum probability (default)
    //  1 - the categories along sites is assigned according to the max posterior probability
    void printResults(string prefix, string ext, int cat_assign_method);
    
    // get number of trees
    int getNumTrees();
    
    /**
        set checkpoint object
        @param checkpoint
    */
    virtual void setCheckpoint(Checkpoint *checkpoint);

    virtual void startCheckpoint();
    
    virtual void saveCheckpoint();
    
    virtual void restoreCheckpoint();

//    virtual ModelFactory* getModelFactory();
//
//    virtual ModelSubst* getModel();

//    virtual RateHeterogeneity* getRate();

    /**
     * save branch lengths into a vector
     */
    void saveBranchLengths(DoubleVector &lenvec, int startid, PhyloNode *node, PhyloNode *dad);

    /**
     * restore branch lengths from a vector previously called with saveBranchLengths
     */
    void restoreBranchLengths(DoubleVector &lenvec, int startid, PhyloNode *node, PhyloNode *dad);
    
    void saveModelCheckpoint();

    void restoreModelCheckpoint();

    /**
        test the best number of threads
    */
    virtual int testNumThreads();

    /**
        compute the weighted average of branch lengths over partitions
    */
    virtual void computeBranchLengths();

    /**
     * Return the tree string contining taxon names and branch lengths
     * @return
     */
    virtual string getTreeString();

    virtual void printResultTree(string suffix);
    
    // a set of phylosupertrees
    vector<PhyloSuperTree*> superTreeSet;
};

#endif
