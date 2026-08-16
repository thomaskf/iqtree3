//
//  modelhmmgm.h
//  model
//
//  Created by Thomas Wong on 03/02/23.
//

#include "modelhmmgm.h"

/**
 constructor
 @param ncat : number of categories
 */
ModelHmmGm::ModelHmmGm(int numcat) : ModelHmm(numcat) {
    ndim = (numcat - 1) * numcat;
    // memory allocation of the array
    size_t transit_size = get_safe_upper_limit(numcat * numcat);
    transit = aligned_alloc<double>(transit_size);
    transit_normalize = aligned_alloc<double>(transit_size);
}

/**
 destructor
 */
ModelHmmGm::~ModelHmmGm() {
    aligned_free(transit);
    aligned_free(transit_normalize);
}

// initialize transitLog array
void ModelHmmGm::initialize_transitLog() {
    size_t i, j;

    // memory allocation of the array
    size_t transit_size = get_safe_upper_limit(ncat * ncat);
    if (transitLog != nullptr)
        aligned_free(transitLog);
    transitLog = aligned_alloc<double>(transit_size);

    // the transition probability going to the same state is initialized to INITIAL_PROB_SAME_CAT
    // and the other probabilities are initialized as equally distributed
    double init_same_tran = INITIAL_PROB_SAME_CAT;
    // the start must satisfy -hmm_min_stran, else the first EM step is forced down
    if (init_same_tran < Params::getInstance().HMM_min_stran)
        init_same_tran = Params::getInstance().HMM_min_stran;
    double init_other_tran = (1.0 - init_same_tran) / ((double) ncat - 1.0);
    for (i = 0; i < ncat; i++) {
        for (j = 0; j < i; j++)
            transit[i * ncat + j] = init_other_tran;
        transit[i * ncat + i] = init_same_tran;
        for (j = i+1; j < ncat; j++)
            transit[i * ncat + j] = init_other_tran;
    }
    computeNormalizedTransits();
    computeLogTransits();
}

/**
 Optimize the model parameters
 @param gradient_epsilon: epsilon for optimization
 @return log-likelihood value
 */
double ModelHmmGm::optimizeParameters(double gradient_epsilon) {

    // changed to EM algorithm
    return optimizeParametersByEM();

    // use BFGS algorithm
    int ndim = getNDim();
    double *variables = new double[ndim + 1];
    double *upper_bound = new double[ndim + 1];
    double *lower_bound = new double[ndim + 1];
    bool *bound_check = new bool[ndim + 1];
    double logLike;
    
    // optimization of transition matrix
    setVariables(variables);
    setBounds(lower_bound, upper_bound, bound_check);
    logLike = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound, bound_check, gradient_epsilon);
    getVariables(variables);
    
    delete [] bound_check;
    delete [] lower_bound;
    delete [] upper_bound;
    delete [] variables;

    return logLike;
}

/**
 Optimize the transition matrix by EM algorithm
 @return log-likelihood value
 */
double ModelHmmGm::optimizeParametersByEM() {
    // Optimize the transition matrix by EM algorithm
    double* marginalTran;
    int sq_ncat = ncat * ncat;
    int dim = phylo_hmm->nsite-1;
    int i,j,k;
    
    // compute the expected value of transition probability between the same categories
    phylo_hmm->computeMarginalTransitProb();
    marginalTran = phylo_hmm->marginal_tran;
    // raw counts, row-normalized below: symmetrizing is not the M-step maximizer
    memset(transit, 0, sizeof(double) * sq_ncat);
    for (k=0; k<dim; k++) {
        for (i=0; i<ncat; i++) {
            for (j=0; j<ncat; j++) {
                transit[i * ncat + j] += marginalTran[i * ncat + j];
            }
        }
        marginalTran += sq_ncat;
    }
    for (i=0; i<sq_ncat; i++) {
        transit[i] = transit[i] / (double) dim;
    }
    computeNormalizedTransits();
    // the bound is on the conditional stay probability, hence after normalization
    double lo = MIN_TRAN_PROB;
    if (Params::getInstance().HMM_min_stran > lo)
        lo = Params::getInstance().HMM_min_stran;
    for (i=0; i<ncat; i++) {
        double* row = transit_normalize + i * ncat;
        double p = (row[i] < lo) ? lo : ((row[i] > 1.0 - MIN_TRAN_PROB) ? 1.0 - MIN_TRAN_PROB : row[i]);
        if (p == row[i])
            continue;
        double rest = 1.0 - row[i];
        for (j=0; j<ncat; j++)
            row[j] = (rest > 0.0) ? row[j] * (1.0 - p) / rest : (1.0 - p) / (double) (ncat - 1);
        row[i] = p;
    }
    computeLogTransits();
    return phylo_hmm->computeBackLike();
}

/**
 Show parameters
 */
void ModelHmmGm::showParameters(ostream& out) {
    size_t i, j, k;
    out << "Estimated HMM transition matrix :" << endl;
    k = 0;
    for (i = 0; i < ncat; i++) {
        for (j = 0; j < ncat; j++) {
            if (j > 0)
                out << "\t";
            out << fixed << setprecision(5) << transit_normalize[k];
            k++;
        }
        out << endl;
    }
}

int ModelHmmGm::getNParameters() {
    return ndim;
}

// transit_normalize is clamped independently of transit, so both are saved
void ModelHmmGm::saveParameters(DoubleVector& params) {
    params.assign(transit, transit + sq_ncat);
    params.insert(params.end(), transit_normalize, transit_normalize + sq_ncat);
}

void ModelHmmGm::restoreParameters(DoubleVector& params) {
    memcpy(transit, params.data(), sizeof(double) * sq_ncat);
    memcpy(transit_normalize, params.data() + sq_ncat, sizeof(double) * sq_ncat);
    computeLogTransits();
}

void ModelHmmGm::setVariables(double *variables) {
    double* var;
    double* tr;
    size_t i,k;
    k = 1;
    // copy the values from transition matrix
    for (i = 0; i < ncat; i++) {
        var = variables + k;
        tr = transit + i*ncat;
        memcpy(var, tr, sizeof(double) * (ncat - 1));
        k += (ncat - 1);
    }
}

void ModelHmmGm::getVariables(double *variables) {
    double* var;
    double* tr;
    size_t i, k;
    k = 1;
    // copy the values to transition matrix
    for (i = 0; i < ncat; i++) {
        var = variables + k;
        tr = transit + i*ncat;
        memcpy(tr, var, sizeof(double) * (ncat - 1));
        k += (ncat - 1);
    }
    computeNormalizedTransits();
    computeLogTransits();
}

void ModelHmmGm::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {

    size_t i;

    for (i = 1; i <= ndim; i++) {
        lower_bound[i] = MIN_VALUE;
        upper_bound[i] = MAX_VALUE;
        bound_check[i] = false;
    }
}

// compute the normalized values of transition matrix
void ModelHmmGm::computeNormalizedTransits() {
    double sum;
    size_t i,j,k;
    for (i = 0; i < ncat; i++) {
        sum = 0.0;
        k = i*ncat;
        for (j = 0; j < ncat; j++) {
            sum += transit[k + j];
        }
        for (j = 0; j < ncat; j++) {
            transit_normalize[k + j] = transit[k + j] / sum;
        }
    }
}

// compute the log values of transition matrix
void ModelHmmGm::computeLogTransits() {
    // an unused transition can be driven to zero by the EM; log(0) would poison
    // transitLog, so floor every entry and renormalize the row
    size_t i, j;
    for (i = 0; i < ncat; i++) {
        double* row = transit_normalize + i * ncat;
        double sum = 0.0;
        for (j = 0; j < ncat; j++) {
            if (!(row[j] >= MIN_TRAN_PROB))
                row[j] = MIN_TRAN_PROB;
            sum += row[j];
        }
        for (j = 0; j < ncat; j++) {
            row[j] = row[j] / sum;
            transitLog[i * ncat + j] = log(row[j]);
        }
    }
}
