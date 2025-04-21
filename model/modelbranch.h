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
     * @return HMM model name
     */
    virtual string getName() { return "BrModel"; }
    
    /**
     * @return HMM model full name
     */
    virtual string getFullName() { return "Branch Model"; }
    
    // value of logl_epsilson
    double logl_epsilon;

};

#endif
