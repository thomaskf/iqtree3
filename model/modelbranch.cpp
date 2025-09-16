#include "modelbranch.h"

// default constructor
ModelBranch::ModelBranch(PhyloTree *tree) : ModelMarkov(tree) {
    logl_epsilon = 0.01;
    
    // whether the root frequencies are separated, by default yes
    // if not, then the root frequencies are same as the frequencies of the base class
    separate_root_freq = true;
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
        optimizeRootFreq(gradient_epsilon);
        // score = optimizeRootFreq(gradient_epsilon);

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

void ModelBranch::getRootFrequency(double* state_freq) {
    if (separate_root_freq) {
        ModelMarkov::getStateFrequency(state_freq);
    } else {
        at(0)->getStateFrequency(state_freq);
    }
}

void ModelBranch::setRootFrequency(double* state_freq) {
    if (separate_root_freq) {
        ModelMarkov::setStateFrequency(state_freq);
        freq_type = FREQ_USER_DEFINED;
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
        delete[] state_freq;
    }
}

/**
 * optimization of root frequencies
 */
double ModelBranch::optimizeRootFreq(double gradient_epsilon) {
    cout << "getNDim() = " << getNDim() << endl;
    return 0.0;
}

void ModelBranch::writeInfo(ostream &out) {
    size_t i;
    for (i=0; i<size(); i++) {
        out << "Branch Model " << i << ":" << endl;
        at(i)->writeInfo(out);
    }
    // report the root frequency
    double* state_freq = new double[num_states];
    getRootFrequency(state_freq);
    out << "Root frequencies:";
    if (num_states == 4) {
        out << "  A: " << state_freq[0];
        out << "  C: " << state_freq[1];
        out << "  G: " << state_freq[2];
        out << "  T: " << state_freq[3];
        out << endl;
    } else if (num_states == 2) {
        out << "  0: " << state_freq[0];
        out << "  1: " << state_freq[1];
        out << endl;
    } else {
        for (i = 0; i < num_states; i++)
            out << " " << state_freq[i];
        out << endl;
    }
    out << endl;
    delete[] state_freq;
}
