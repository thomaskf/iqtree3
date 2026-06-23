#ifndef MODELTEA_H
#define MODELTEA_H

#include "modelprotein.h"

extern const char* builtin_tea_models;

/** Empirical models for the TEA structural alphabet, reusing ModelProtein. */
class ModelTea : public ModelProtein
{
public:
	ModelTea(const char *model_name, string model_params, StateFreqType freq, string freq_params, PhyloTree *tree, ModelsBlock *models_block);

	virtual void computeTipLikelihood(PML::StateType state, double *state_lk);
};

#endif
