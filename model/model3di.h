#ifndef MODEL3DI_H
#define MODEL3DI_H

#include "modelprotein.h"

extern const char* builtin_3di_models;

/** Empirical models for the 3Di structural alphabet, reusing ModelProtein. */
class Model3Di : public ModelProtein
{
public:
	Model3Di(const char *model_name, string model_params, StateFreqType freq, string freq_params, PhyloTree *tree, ModelsBlock *models_block);

	virtual void computeTipLikelihood(PML::StateType state, double *state_lk);
};

#endif
