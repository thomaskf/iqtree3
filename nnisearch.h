#ifndef NNISEARCH_H
#define NNISEARCH_H

//#ifdef __cplusplus
//extern "C" {
//#endif

#include <stdio.h>
#include <stdlib.h>
#include "tools.h"
#include "pll/pll.h"
#include "pll/hash.h"
#include <string>
#include <sstream>
#include <set>
#include <algorithm>
#include <map>
#include <vector>
#include <unordered_set>
using namespace std;

const int TOPO_ONLY = 0;
const int NO_BRAN_OPT = 1;
const int ONE_BRAN_OPT = 2;
const int FIVE_BRAN_OPT = 4;

#define IQTREE_NEWZPERCYCLE 10

/* This is the info you need to copy the vector*/
typedef struct
{
  int node_number;
  int num_partitions;
  size_t *partition_sizes;
  double **lh_values;
  int **expVector;
} LH_VECTOR;


typedef struct {
	nodeptr p;
	int nniType;
	char* idString;
	string quartetString;
    double z0[PLL_NUM_BRANCHES]; // p
    double z1[PLL_NUM_BRANCHES]; // p->next
    double z2[PLL_NUM_BRANCHES]; // p->next->next
    double z3[PLL_NUM_BRANCHES]; // q->next
    double z4[PLL_NUM_BRANCHES]; // q->next->next
	double likelihood;
	double loglDelta;
	double negLoglDelta;
} pllNNIMove;

inline bool comparePLLNNIMove(const pllNNIMove &a, const pllNNIMove &b)
{
    return a.likelihood < b.likelihood;
}

static int cmp_nni(const void* nni1, const void* nni2);

int compareDouble(const void * a, const void * b);

pllNNIMove *getNonConflictNNIList(pllInstance* tr);

void _update(pllInstance *tr, partitionList *pr, nodeptr p);

#define MAX_NUM_DELTA 10000

typedef struct {
	double delta[MAX_NUM_DELTA];
	int num_delta;
	double delta_min;
	int doNNICut;
} NNICUT;

typedef struct {
	//unordered_map<string, pllNNIMove> nniList;
	vector<pllNNIMove> nniList;
	bool updateNNIList;
	unordered_set<string> tabuNNIs;
	vector<pllNNIMove> posNNIList; // positive NNI list
	bool tabunni;
	unordered_set<int> affectNodes; // Set of nodes that are affected by the previous NNIs
	unordered_set<string> affectBranches;
	double curLogl;
	int numUnevalQuartet; // number of unevaluated quartet because of tabu constraint
	int evalType;
	int numAppliedNNIs; // total number of applied NNIs sofar
	int curNumAppliedNNIs; // number of applied NNIs at the current step
	int curNumNNISteps;
} SearchInfo;

/**
 * get all the nodes affected by the NNI
 * @param p
 * @return
 */
set<int> getAffectedNodes(pllInstance* tr, nodeptr p);

string getBranString(nodeptr p);

bool containsAffectedNodes(nodeptr p, SearchInfo &searchinfo);

void pllSaveQuartet(nodeptr p, SearchInfo &searchinfo);

void updateBranchLengthForNNI(pllInstance* tr, partitionList *pr, pllNNIMove &nni);

void pllEvalAllNNIs(pllInstance *tr, partitionList *pr, SearchInfo &searchinfo);

double pllDoRandNNIs(pllInstance *tr, partitionList *pr, int numNNI);

/**
 *  Evaluate NNI moves for the current internal branch
 *  @param tr the current tree data structure
 *  @param pr partition data structure
 *  @param p the node representing the current branch
 *  @return number of positive NNIs found
 */
int evalNNIForBran(pllInstance* tr, partitionList *pr, nodeptr p, SearchInfo &searchinfo);

/**
 * Perturb the best tree
 *
 * Given the best tree, apply some NNIs to escape local optimum
 * @param tr
 * @param pr
 * @param nnis list of all NNI to apply
 * @param numNNI size of the array nnis
 * @return
 */
double pllPerturbTree(pllInstance *tr, partitionList *pr, vector<pllNNIMove> &nnis);

/**
 * 	do 1 round of fastNNI
 *  return new tree log-likelihood if found improving NNI otherwise -1.0
 *
 *  @param[in] tr the tree data structure
 *  @param[in] pr partition data structure
 *  @param[out] nniList list containing information about the 2(n-3) evaluated NNIs
 *  @param[in/out] tabuNNIs tabu list
 *  @param[out] nni_count number of NNI that have been applied
 *  @param[out] deltaNNI average improvement made by one NNI
 */
double pllDoNNISearch(pllInstance* tr, partitionList *pr, SearchInfo &searchinfo);

void pllUpdateTabuList(pllInstance *tr, SearchInfo &searchinfo);

void pllSaveQuartetForSubTree(pllInstance* tr, nodeptr p, SearchInfo &searchinfo);



/**
 *  @brief Do 1 NNI move.
 *  @param[in] tr: the tree data structure
 *  @param[in] pr partition data structure
 *  @param[in] swap: represents one of the 2 NNI moves. Could be either 0 or 1
 *  @param[in] evalType: NO_NR, WITH_ONE_NR, WITH_FIVE_NR
 */
double doOneNNI(pllInstance * tr, partitionList *pr, nodeptr p, int swap, int evalType);

void pllGetAllInBran(pllInstance *tr, vector<nodeptr> &branlist);

void pllGetAllInBranForSubtree(pllInstance *tr, nodeptr p, vector<nodeptr> &branlist);


string convertQuartet2String(nodeptr p);
/**
 *  Go through all 2(n-3) internal branches of the tree and
 *  evaluate all possible NNI moves
 */
void evalAllNNI(pllInstance* tr);

/**
 * 	@brief evaluate all NNIs within the subtree specified by node p
 * 	populates the list containing all possible NNI moves
 *
 * 	@param[in] tr: the tree data structure
 * 	@param[in] pr partition data structure
 * 	@param[in] p node pointer that specify the subtree
 */
void evalNNIForSubtree(pllInstance* tr, partitionList *pr, nodeptr p, SearchInfo &searchinfo);

/*
 *  @brief return the array which can be used to store evaluated NNI moves
 *
 *  @param[in] tr: the tree data structure
 */
pllNNIMove *getNNIList(pllInstance* tr);



/*
 * ****************************************************************************
 * pllUFBoot area
 * ****************************************************************************
 */

/**
 * DTH:
 * pllUFBootData struct
 * This one keeps all info necessary to run UFBoot in PLL mode
 */
typedef struct{
	int max_candidate_trees;
	int treels_size;
	int save_all_trees;
	pll_boolean save_all_br_lens;
	double logl_cutoff;
	int duplication_counter;
	int n_patterns;
	struct pllHashTable * treels;
	unsigned int candidate_trees_count; /* counter of trees in pllHashTable */
	double * treels_logl; // maintain size == treels_size
	char ** treels_newick; // maintain size == treels_size
	double ** treels_ptnlh; // maintain size == treels_size
	int ** boot_samples;
	double * boot_logl;
	int * boot_counts;
	int * boot_trees;
//	double * random_doubles;
} pllUFBootData;

/**
 * DTH:
 * The PLL version of saveCurrentTree function
 * @param tr: the tree (a pointer to a pllInstance)
 * @param pr: pointer to a partitionList (this one keeps tons of tree info)
 * @param p: root?
 */
void pllSaveCurrentTree(pllInstance* tr, partitionList *pr, nodeptr p);

/**
 * DTH:
 * Extract the array of site log likelihood to be kept in ptnlh
 * And update *cur_log
 * @param tr: the tree (pointer to an pllInstance)
 * @param ptnlh: to-be-kept array of site log likelihood
 * @param cur_logl: pointer to current tree log likelihood
 */
void pllComputePatternLikelihood(pllInstance* tr, double * ptnlh, double * cur_logl);

/**
 * DTH:
 * Announce the memory allocation error (for debugging)
 */
void pllAlertMemoryError();

/**
 * DTH:
 * Resize some of the arrays in UFBootData if they're full
 * Along with update treels_size (to track the size of these arrays)
 */
void pllResizeUFBootData();

/**
 * DTH:
 * (Based on function Tree2StringREC of PLL)
 * Print out the tree topology with IQTree taxa ID (starts at 0) instead of PLL taxa ID (starts at 1)
 * @param All are the same as in PLL's
 */
static char *pllTree2StringREC(char *treestr, pllInstance *tr, partitionList *pr, nodeptr p, pll_boolean printBranchLengths, pll_boolean printNames,
			    pll_boolean printLikelihood, pll_boolean rellTree, pll_boolean finalPrint, int perGene, pll_boolean branchLabelSupport, pll_boolean printSHSupport);


//#ifdef __cplusplus
//}
//#endif

#endif

