//
//  iqtreemixhmm.h
//  tree
//
//  Created by Thomas Wong on 19/01/23.
//

#ifndef iqtreemixhmm_h
#define iqtreemixhmm_h

#include <cmath>
#include "iqtreemix.h"
#include "tree/phylohmm.h"
#include "model/modelhmm.h"
#include "model/modelhmmgm.h"

// snapshot of every parameter an optimization step may update
struct HmmParamSnapshot {
    Checkpoint model_ckp;          // substitution models, RHAS models and tree weights
    vector<DoubleVector> brlens;   // branch lengths of every tree
    DoubleVector prob_arr;         // HMM category probabilities
    DoubleVector tran_par;         // HMM transition model parameters
};

class IQTreeMixHmm : public IQTreeMix, public PhyloHmm {
public:
    
    /**
     default constructor
     */
    IQTreeMixHmm();
    
    IQTreeMixHmm(Params &params, Alignment *aln);
    
    /**
     destructor
     */
    ~IQTreeMixHmm() override;
    
    // initialize the model
    void initializeModel(Params &params, string model_name, ModelsBlock *models_block) override;

    // initialize the transition model
    void initializeTransitModel(Params &params);
    
    // initialize the parameters
    void initializeParams();
    
    // set the tree weights according to the marginal probabilities along the sites
    void setWeightToMarginalProb();

    // compute the log of dotproduct of the logorithm arrays
    double logDotProd(double* ln_x, double* ln_y, int n);
    
    // obtain the log-likelihoods for every pattern and every tree
    // output ptn_like_cat[i * ntree + j] : log-likelihood of pattern i and tree j
    void computeLogLikelihoodSiteTree(int updateTree = -1);
    
    // compute backward log-likelihood
    virtual double computeLikelihood(double *pattern_lh = nullptr, bool save_log_value = true) override;
    
    /**
     optimize all branch lengths of one tree
     @param my_iterations number of iterations to loop through all branches
     */
    void optimizeAllBranchesOneTree(int whichtree, int my_iterations = 100, double tolerance = TOL_LIKELIHOOD, int maxNRStep = 100);
    
    /**
     optimize all branch lengths of all trees
     @param my_iterations number of iterations to loop through all branches
     @return the likelihood of the tree
     */
    double optimizeAllBranches(double* pattern_mix_lh = nullptr, int my_iterations = 100, double tolerance = TOL_LIKELIHOOD, int maxNRStep = 100);
    
    double optimizeAllBranchLensByBFGS(double gradient_epsilon, double logl_epsilon, int maxsteps = 3);
    /**
     @return true if this is a HMM model
     */
    virtual bool isHMM() override { return true; }
    
    virtual void startCheckpoint() override;
    
    virtual string optimizeModelParameters(bool printInfo, double logl_epsilon) override;

    // Optimize parameters according to the MAST model
    string optimizeModelParamMAST(bool printInfo, double logl_epsilon);

    // Optimize parameters according to the HMM model
    string optimizeModelParamHMM(bool printInfo, double logl_epsilon);

    virtual void setNumThreads(int num_threads) override;
    
    /**
     test the best number of threads
     */
    virtual int testNumThreads() override;
    
    // number of parameters under a given objective function (0: HMM, 1: MAST)
    int getNParameters(int obj_fun);

    virtual int getNParameters() override { return getNParameters(objFun); }
    
    // print out all the results to a file
    void printResults(const char *filename, int cat_assign_method = 0, int* numSiteCat = nullptr);

    // print out the marginal probabilities to a file
    void printMarginalProb(const char *filename);

    // show the values of the parameters
    void showParameters(ostream& out);
    
    // optimize all substitution models
    double optimizeAllSubstModels(double gradient_epsilon, double* pattern_mix_lh = nullptr);
    
    // optimize all RHAS models
    double optimizeAllRHASModels(double gradient_epsilon, double score = 0.0, double* pattern_mix_lh = nullptr);
    
private:
    
    // indicate which tree's unlinked parameters is under optimization
    // -1   : parameters affecting all trees (default)
    // >= 0 : parameters affecting a specific tree
    int optimTree;
    
    // indicate which branch group is under optimization
    // -1  : no branch group is under optimization
    // >=0 : optimizing a specific branch group
    int optimBranchGrp;
    
    // which objective function (default: MAST)
    // 0: backLikelihood
    // 1: MAST
    int objFun;
    
    // whether using the optimization engine in IQTreeMix
    bool isTMixOptimEngine;
    
    string* objAlgo;
    
    // branch lengths of all the trees
    vector<DoubleVector> allbranchlens;
    
    // compute the log-likelihoods for a single tree t
    void computeLogLikelihoodSingleTree(int t);

    // get the branch lengths of all trees to the variable allbranchlens
    void getAllBranchLengths();

    // set the branch lengths of all trees from the variable allbranchlens
    void setAllBranchLengths();
    
    // show the branch lengths of all trees
    // void showAllBranchLengths();
    
    //--------------------------------------------
    // optimization of branch lengths using BFGS
    //--------------------------------------------

    // the following three functions are for dimension = 1
    double computeFunction(double x) override;
    
    double setSingleVariable();
    
    void getSingleVariable(double x);

    // the following four functions are for dimension > 1
    virtual double targetFunk(double x[]) override;
    
    virtual void setVariables(double *variables);
    
    virtual void getVariables(double *variables);
    
    virtual void setBounds(double *lower_bound, double *upper_bound, bool *bound_check);

    virtual int getNDim() override;
    
    double optimizeBranchGroup(int branchgrp, double gradient_epsilon);
    
    void showBranchGrp();

    /**
             If there are multiple branches belonging to the same group
             set all the branches of the same group to their average
     */
    void setAvgLenEachBranchGrp();
    
    // update the ptn_freq array according to the marginal probabilities along each site for each tree
    void computeFreqArray(double* pattern_mix_lh = nullptr, bool need_computeLike = true, int update_which_tree = -1);
    
    // get marginal probabilities along each site for each tree
    void getMarginalProb(bool need_computeLike = true, int update_which_tree = -1);

    // redirect the checkpoint used by saveModelCheckpoint / restoreModelCheckpoint
    void setModelCheckpoint(Checkpoint* ckp);

    // save / restore every parameter an optimization step may update
    void saveHmmParams(HmmParamSnapshot& snapshot);
    void restoreHmmParams(HmmParamSnapshot& snapshot);
};

#endif /* iqtreemixhmm_h */
