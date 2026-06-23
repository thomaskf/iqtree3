#ifndef ALPHABET_CLASSIFIER_H
#define ALPHABET_CLASSIFIER_H

#include <vector>
#include <string>
#include "phylo-yaml/statespace.h"

/** Classify amino-acid-letter data as protein, 3Di, or TEA (dipeptide-PCA classifier). */
PML::SeqType classifyProteinLikeAlphabet(std::vector<std::string> &sequences);

#endif
