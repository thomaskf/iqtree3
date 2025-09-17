#include "modelbranch.h"

// default constructor
ModelBranch::ModelBranch(PhyloTree *tree) : ModelMarkov(tree, true, false) {
    logl_epsilon = 0.01;
    
    // whether the root frequencies are separated, by default yes
    // if not, then the root frequencies are same as the frequencies of the base class
    separate_root_freq = Params::getInstance().separate_root_freq;
    
    // only stores the root frequency, the other classes are stored in the array
    num_params = 0;
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
    ModelMarkov::saveCheckpoint();
    endCheckpoint();
}

void ModelBranch::restoreCheckpoint() {
    startCheckpoint();
    for (int i = 0; i < size(); i++) {
        checkpoint->startStruct("Component" + convertIntToString(i));
        at(i)->restoreCheckpoint();
        checkpoint->endStruct();
    }
    ModelMarkov::restoreCheckpoint();
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
    
    int ndim = getNDim();
    if (ndim == 0)
        return 0.0;
    
    // since the model is not reversible, should not call ModelMarkov::optimizeParameters()
    double *variables = new double[ndim+1]; // used for BFGS numerical recipes
    double *variables2 = new double[ndim+1]; // used for L-BFGS-B
    double *upper_bound = new double[ndim+1];
    double *lower_bound = new double[ndim+1];
    bool *bound_check = new bool[ndim+1];
    double score;

    // by BFGS algorithm
    setVariables(variables);
    setBounds(lower_bound, upper_bound, bound_check);
    score = -minimizeMultiDimen(variables, ndim, lower_bound, upper_bound, bound_check, gradient_epsilon);
    
    // scale the frequencies
    scaleStateFreq(true);

    delete [] bound_check;
    delete [] lower_bound;
    delete [] upper_bound;
    delete [] variables2;
    delete [] variables;

    return score;
    
    /*
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
        if (separate_root_freq) {
            score = ModelMarkov::optimizeParameters(gradient_epsilon);
        }
        if (score < prev_score + logl_epsilon) {
            // converged
            break;
        }
        prev_score = score;
    }
    return score;
    */
}

string ModelBranch::getName() {
    if (name != "") return name;
    string retname = "BR";
    retname += OPEN_BRACKET;
    for (iterator it = begin(); it != end(); it++) {
        if (it != begin()) retname += ",";
        retname += (*it)->getName();
    }
    retname += CLOSE_BRACKET;
    return retname;
}

string ModelBranch::getNameParams(bool show_fixed_params) {
    // if (full_name != "")
    //    return full_name;
    string retname = "BR";
    retname += OPEN_BRACKET;
    for (iterator it = begin(); it != end(); it++) {
        if (it != begin()) retname += ",";
        retname += (*it)->getNameParams(show_fixed_params);
    }
    retname += CLOSE_BRACKET;
    return retname;
}

void ModelBranch::getRootFrequency(double* sfreq) {
    if (separate_root_freq) {
        ModelMarkov::getStateFrequency(sfreq);
    } else {
        at(0)->ModelMarkov::getStateFrequency(sfreq);
    }
}

void ModelBranch::setRootFrequency(double* state_freq) {
    if (separate_root_freq) {
        ModelMarkov::setStateFrequency(state_freq);
    }
}

void ModelBranch::setRootFrequency(string root_freq) {
    if (separate_root_freq) {
        double* state_freq = new double[num_states];
        // check the format of rootfreq
        size_t i = 0;
        size_t s = 0;
        size_t p = root_freq.find_first_of(" ,/", s);
        string f_str;
        while (p != string::npos) {
            f_str = root_freq.substr(s, p-s);
            if (i < num_states) {
                state_freq[i] = atof(f_str.c_str());
                i++;
            } else {
                outError("The number of root frequencies inputted is more than " + convertIntToString(num_states));
            }
            s = p+1;
            p = root_freq.find_first_of(" ,/", s);
        }
        f_str = root_freq.substr(s);
        if (i < num_states) {
            state_freq[i] = atof(f_str.c_str());
            i++;
        } else {
            outError("The number of root frequencies inputted is more than " + convertIntToString(num_states));
        }
        if (i < num_states) {
            outError("The number of root frequencies inputted is less than " + convertIntToString(num_states));
        }
        setRootFrequency(state_freq);
        freq_type = FREQ_USER_DEFINED;
        delete[] state_freq;
    }
}

void ModelBranch::showRootFreq() {
    double* f = new double[num_states];
    getRootFrequency(f);
    cout << "Root frequencies:";
    for (size_t i = 0; i < num_states; i++)
        cout << " " << f[i];
    cout << endl;
    delete[] f;
}

void ModelBranch::writeInfo(ostream &out) {
    size_t i;
    for (i=0; i<size(); i++) {
        out << "Branch Model " << i << ":" << endl;
        at(i)->writeInfo(out);
    }
    // report the root frequency
    double* sfreq = new double[num_states];
    getRootFrequency(sfreq);
    out << "Root frequencies:";
    if (num_states == 4) {
        out << "  A: " << sfreq[0];
        out << "  C: " << sfreq[1];
        out << "  G: " << sfreq[2];
        out << "  T: " << sfreq[3];
        out << endl;
    } else if (num_states == 2) {
        out << "  0: " << sfreq[0];
        out << "  1: " << sfreq[1];
        out << endl;
    } else {
        for (i = 0; i < num_states; i++)
            out << " " << sfreq[i];
        out << endl;
    }
    out << endl;
    delete[] sfreq;
}


/**
 @return the number of dimensions
 */
int ModelBranch::getNDim() {
    int totndim = 0;
    for (iterator it = begin(); it != end(); it++)
        totndim += (*it)->getNDim();
    
    if (separate_root_freq) {
        totndim += ModelMarkov::getNDim();
    }
    
    return totndim;
}

int ModelBranch::getNDimFreq() {
    int dim = 0;
    int num_empirical = 0;
    int num_codon_1x4 = 0;
    int num_codon_3x4 = 0;
    for (iterator it = begin(); it != end(); it++) {
        // count empirical freq only once
        switch ((*it)->freq_type) {
        case FREQ_EMPIRICAL:
            num_empirical++;
            if (num_empirical==1)
                dim += (*it)->getNDimFreq();
            break;
        case FREQ_CODON_1x4:
            num_codon_1x4++;
            if (num_codon_1x4==1)
                dim += (*it)->getNDimFreq();
            break;
        case FREQ_CODON_3x4:
        case FREQ_CODON_3x4C:
            num_codon_3x4++;
            if (num_codon_3x4==1)
                dim += (*it)->getNDimFreq();
            break;
        default:
            dim += (*it)->getNDimFreq();
        }
    }
    
    if (separate_root_freq) {
        // consider the root frequency
        switch (freq_type) {
        case FREQ_EMPIRICAL:
            num_empirical++;
            if (num_empirical==1)
                dim += ModelMarkov::getNDimFreq();
            break;
        case FREQ_CODON_1x4:
            num_codon_1x4++;
            if (num_codon_1x4==1)
                dim += ModelMarkov::getNDimFreq();
            break;
        case FREQ_CODON_3x4:
        case FREQ_CODON_3x4C:
            num_codon_3x4++;
            if (num_codon_3x4==1)
                dim += ModelMarkov::getNDimFreq();
            break;
        default:
            dim += ModelMarkov::getNDimFreq();
        }
    }
    return dim;
}

void ModelBranch::setVariables(double *variables) {
    if (getNDim() == 0)
        return;
    int dim = 0;
    for (iterator it = begin(); it != end(); it++) {
        (*it)->ModelMarkov::setVariables(&variables[dim]);
        dim += (*it)->ModelMarkov::getNDim();
    }
    if (separate_root_freq) {
        double* var = &variables[dim+1];
        memcpy(var, ModelMarkov::state_freq, (ModelMarkov::num_states-1)*sizeof(double));
        dim += (ModelMarkov::num_states-1);
    }
}

bool ModelBranch::getVariables(double *variables) {
    if (getNDim() == 0)
        return false;
    int dim = 0;
    bool changed = false;
    for (iterator it = begin(); it != end(); it++) {
        changed |= (*it)->ModelMarkov::getVariables(&variables[dim]);
        dim += (*it)->ModelMarkov::getNDim();
    }
    if (separate_root_freq) {
        double* var = &variables[dim+1];
        for (size_t i = 0; i < ModelMarkov::num_states-1; i++)
            changed |= (ModelMarkov::state_freq[i] != var[i]);
        memcpy(ModelMarkov::state_freq, var, (ModelMarkov::num_states-1)*sizeof(double));
        dim += (ModelMarkov::num_states-1);
    }
    return changed;
}

void ModelBranch::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
    if (getNDim() == 0)
        return;
    int dim = 0;
    for (iterator it = begin(); it != end(); it++) {
        (*it)->ModelMarkov::setBounds(&lower_bound[dim], &upper_bound[dim], &bound_check[dim]);
        dim += (*it)->ModelMarkov::getNDim();
    }
    if (separate_root_freq) {
        for (size_t i = 0; i < ModelMarkov::num_states-1; i++) {
            lower_bound[dim+1+i] = Params::getInstance().min_state_freq;;
            upper_bound[dim+1+i] = 1.0;
            bound_check[dim+1+i] = false;
        }
        dim += (ModelMarkov::num_states-1);
    }
}

double ModelBranch::targetFunk(double x[]) {
    if (verbose_mode >= VB_DEBUG) {
        int ndim = getNDim();
        for(int i=1; i<=ndim; i++){ cout << x[i] << "; "; }
        cout << endl;
    }
    getVariables(x);
    decomposeRateMatrix();
    phylo_tree->clearAllPartialLH();
    
    double score = -phylo_tree->computeLikelihood();
    return score;
}

void ModelBranch::scaleStateFreq(bool sum_one) {
    for (iterator it = begin(); it != end(); it++) {
        (*it)->scaleStateFreq(sum_one);
    }
    if (separate_root_freq) {
        ModelMarkov::scaleStateFreq(sum_one);
    }
}
