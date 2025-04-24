#include "modelbranch.h"

// default constructor
ModelBranch::ModelBranch(PhyloTree *tree) : ModelMarkov(tree) {
    logl_epsilon = 0.01;
}

// destructor
ModelBranch::~ModelBranch() {
    for (int i = 0; i < size(); i++) {
        delete(at(i));
    }
    clear();
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
    int totndim = 0;
    for (iterator it = begin(); it != end(); it++)
        totndim += (*it)->getNDim();
    return totndim;
}
