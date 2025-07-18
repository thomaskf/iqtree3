/*
 * modelbranch.h
 *
 *  Created on: April 17, 2025
 *      Author: Thomas Wong
 */

#ifndef MODELBRANCH_H_
#define MODELBRANCH_H_

#include "modelmarkov.h"
#include "tree/phylotree.h"

/**
 * branch model
 */
class ModelBranch: virtual public ModelMarkov, public vector<ModelMarkov*> {
public:
    
    // default constructor
    ModelBranch(PhyloTree *tree);
    
    // destructor
    ~ModelBranch();
    
    virtual void setCheckpoint(Checkpoint *checkpoint);
    
    virtual void startCheckpoint();
    
    virtual void saveCheckpoint();
    
    virtual void restoreCheckpoint();
    
    virtual void decomposeRateMatrix();
    
    virtual double optimizeParameters(double gradient_epsilon);
    
    /**
     @return the number of dimensions
     */
    virtual int getNDim();
    
    /**
     * @return model name
     */
    virtual string getName();
    
    /**
     * @return model name and parameters
     */
    virtual string getNameParams(bool show_fixed_params);

    /*
     * initialization of root frequencies
     */
    void initializeRootFreq();
    
    /**
        @return TRUE if model is time-reversible, FALSE otherwise
    */
    virtual bool isReversible() { return true;  }

    /**
     * optimization of root frequencies
     */
    double optimizeRootFreq(double gradient_epsilon);
    
    // value of logl_epsilson
    double logl_epsilon;
    
    // whether the root frequencies are optimized separately, by default it is true
    bool opt_root_freq;
    
    double* rootfreqs;

private:
    
    // in the progress of root frequencies optimization
    bool is_optimizing_root;
};

#endif
