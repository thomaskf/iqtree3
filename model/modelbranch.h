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

    void getRootFrequency(double* state_freq);
    
    void setRootFrequency(double* state_freq);
    
    void setRootFrequency(string root_freq);
    
    /**
        @return TRUE if model is time-reversible, FALSE otherwise
    */
    virtual bool isReversible() { return false;  }

    /**
        write information
        @param out output stream
    */
    virtual void writeInfo(ostream &out);

    /**
     * optimization of root frequencies
     */
    double optimizeRootFreq(double gradient_epsilon);
    
    // value of logl_epsilson
    double logl_epsilon;
    
    // whether the root frequencies are separated, by default yes
    // if not, then the root frequencies are same as the frequencies of the base class
    bool separate_root_freq;
};

#endif
