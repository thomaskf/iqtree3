#include "modelbranch.h"

// default constructor
ModelBranch::ModelBranch(PhyloTree *tree) : ModelMarkov(tree) {
    logl_epsilon = 0.01;
    
    /*
    // by default, the root frequency needs to optimize separately
    opt_root_freq = true;
    is_optimizing_root = false;
    rootfreqs = NULL;
    */
}

// destructor
ModelBranch::~ModelBranch() {
    for (int i = 0; i < size(); i++) {
        delete(at(i));
    }
    clear();
    if (rootfreqs != NULL)
        delete(rootfreqs);
}

void ModelBranch::setCheckpoint(Checkpoint *checkpoint) {
    CheckpointFactory::setCheckpoint(checkpoint);
    for (iterator it = begin(); it != end(); it++)
        (*it)->setCheckpoint(checkpoint);
}

void ModelBranch::startCheckpoint() {
    checkpoint->startStruct("BranchModel" + convertIntToString(size()));
}

void ModelBranch::saveCheckpoint() {
    startCheckpoint();
    for (int i = 0; i < size(); i++) {
        checkpoint->startStruct("Component" + convertIntToString(i));
        at(i)->saveCheckpoint();
        checkpoint->endStruct();
    }
    endCheckpoint();
}

void ModelBranch::restoreCheckpoint() {
    startCheckpoint();
    for (int i = 0; i < size(); i++) {
        checkpoint->startStruct("Component" + convertIntToString(i));
        at(i)->restoreCheckpoint();
        checkpoint->endStruct();
    }
    endCheckpoint();
    decomposeRateMatrix();
    if (phylo_tree)
        phylo_tree->clearAllPartialLH();
}

void ModelBranch::decomposeRateMatrix() {
    for (iterator it = begin(); it != end(); it++)
        (*it)->decomposeRateMatrix();
}

double ModelBranch::optimizeParameters(double gradient_epsilon) {
    double prev_score = phylo_tree->computeLikelihood();
    double score = prev_score;
    int optimize_steps = 100;
    int digit_prec = 5;
    int step, k;

    if (verbose_mode >= VB_DEBUG) {
        cout << std::setprecision(digit_prec) << "ModelBranch -- optimizing parameters (gradient_epsilon = " << gradient_epsilon << "; logl_epsilon = " << logl_epsilon << ")" << endl;
    }

    for (step = 0; step < optimize_steps; step++) {
        for (k = 0; k < size(); k++) {
            if (!at(k)->fixed_parameters && at(k)->getNDim() > 0) {
                score = at(k)->optimizeParameters(gradient_epsilon);
                if (verbose_mode >= VB_DEBUG) {
                    cout << std::setprecision(digit_prec) << "step " << step << "; model " << k << ": " << score << endl;
                }
            }
        }
        // optimize the root frequency if necessary
        if (opt_root_freq) {
            optimizeRootFreq(gradient_epsilon);
            // score = optimizeRootFreq(gradient_epsilon);
        }
        if (score < prev_score + logl_epsilon) {
            // converged
            break;
        }
        prev_score = score;
    }
    return score;
}

/**
 @return the number of dimensions
 */
int ModelBranch::getNDim() {
    if (is_optimizing_root) {
        return (rootfreqs != NULL)?(num_states-1):0;
    }
    int totndim = 0;
    for (iterator it = begin(); it != end(); it++)
        totndim += (*it)->getNDim();
    return totndim;
}

/*
 * initialization of root frequencies
 */
void ModelBranch::initializeRootFreq() {
    rootfreqs = new double[num_states];
    if (size() > 0) {
        for (int i = 0; i < num_states; i++) {
            rootfreqs[i] = at(0)->state_freq[i];
        }
    }
}

/**
 * optimization of root frequencies
 */
double ModelBranch::optimizeRootFreq(double gradient_epsilon) {
    is_optimizing_root = true;
    cout << "getNDim() = " << getNDim() << endl;
    is_optimizing_root = false;
    return 0.0;
}

