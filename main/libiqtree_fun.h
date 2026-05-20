#ifndef LIBIQTREE2_FUN
#define LIBIQTREE2_FUN

#include "tree/mtree.h"
#include "tree/phylotree.h"
#include "utils/tools.h"
#include "main/alisim.h"
#include "utils/starttree.h"
//#include "suppFunc.h"
#include <vector>
#include <string>
#include <cstring>
#include <fstream>

using namespace std;

#ifdef _MSC_VER
#pragma pack(push, 1)
#else
#pragma pack(1)
#endif

typedef struct {
  const char** strings;
  size_t length;
} StringArray;

typedef struct {
  double* doubles;
  size_t length;
} DoubleArray;

typedef struct {
  int value;
  char* errorStr;
} IntegerResult;

typedef struct {
  char* value;
  char* errorStr;
} StringResult;

typedef struct {
  double* value;
  size_t length;
  char* errorStr;
} DoubleArrayResult;

/*
 * Result of fit_tree_hessian: the fitted rooted tree (Newick), per-branch
 * lengths, log-likelihood gradient (first derivatives w.r.t. branch lengths)
 * and the n x n Hessian matrix (row-major). All pointers are owned by the
 * caller and must be released with iqtree_free.
 */
typedef struct {
  char* tree;              /* Newick tree string (matches branch_lengths order) */
  double* branch_lengths;  /* length n_branches */
  double* gradient;        /* length n_branches */
  double* hessian;         /* length n_branches * n_branches, row-major */
  size_t n_branches;
  double log_likelihood;
  char* errorStr;
} HessianResult;

#ifdef _MSC_VER
#pragma pack(pop)
#else
#pragma pack()
#endif

/*
 * Calculates the robinson fould distance between two trees
 */
extern "C" IntegerResult robinson_fould(const char* ctree1, const char* ctree2);

/*
 * Generates a set of random phylogenetic trees
 * tree_gen_mode allows:"YULE_HARDING", "UNIFORM", "CATERPILLAR", "BALANCED", "BIRTH_DEATH", "STAR_TREE"
 * output: a newick tree (in string format)
 */
extern "C" StringResult random_tree(int num_taxa, const char* tree_gen_mode, int num_trees, int rand_seed = 0);

/*
 * Perform phylogenetic analysis on the input alignment
 * With estimation of the best topology
 * num_thres -- number of cpu threads to be used, default: 1; 0 - auto detection of the optimal number of cpu threads
 * output: results in YAML format with the tree and the details of parameters
 */
extern "C" StringResult build_tree(StringArray& names, StringArray& seqs, const char* model, int rand_seed = 0, int bootstrap_rep = 0, int num_thres = 1, const char* other_options = NULL);

/*
 * Perform phylogenetic analysis on the input alignment
 * With restriction to the input toplogy
 * blfix -- whether to fix the branch length as those on the given tree, default: false
 * num_thres -- number of cpu threads to be used, default: 1; 0 - auto detection of the optimal number of cpu threads
 * output: results in YAML format with the details of parameters
 */
extern "C" StringResult fit_tree(StringArray& names, StringArray& seqs, const char* model, const char* intree, bool blfix = false, int rand_seed = 0, int num_thres = 1, const char* other_options = NULL);

/*
 * Like fit_tree, but in addition to fitting the tree this also computes the
 * log-likelihood gradient (first derivatives w.r.t. branch lengths) and the
 * full n x n Hessian matrix at the fitted point. The branch lengths, gradient
 * and Hessian are returned in the MCMCTree branch ordering and share that
 * order with the returned Newick tree string.
 *
 * Only the non-partitioned (single-tree) case is supported; passing a
 * partitioned model raises an error.
 *
 * blfix      -- whether to fix the branch lengths as those on the given tree
 * num_thres  -- number of threads, default 1; 0 = auto-detect
 * other_options -- additional CLI-style options forwarded to IQ-TREE
 *
 * The caller owns the returned tree/branch_lengths/gradient/hessian pointers
 * and must release each one with iqtree_free.
 */
extern "C" HessianResult fit_tree_hessian(StringArray& names, StringArray& seqs, const char* model, const char* intree, bool blfix = false, int rand_seed = 0, int num_thres = 1, const char* other_options = NULL);

/*
 * Open a persistent Hessian session; returned handle is reused by
 * hessian_session_eval and freed by hessian_session_destroy.
 */
extern "C" void* hessian_session_create(StringArray& names, StringArray& seqs, const char* model, const char* intree, bool blfix, int rand_seed, int num_thres, const char* other_options, char** errorStr);

/* Re-evaluate at new branch lengths (length 0 = use current state). */
extern "C" HessianResult hessian_session_eval(void* handle, DoubleArray& branch_lengths);

/* Log-likelihood only at new branch lengths; skips per-branch derivative work. */
extern "C" double hessian_session_score(void* handle, DoubleArray& branch_lengths, char** errorStr);

extern "C" void hessian_session_destroy(void* handle);

/*
 * Perform phylogenetic analysis with ModelFinder
 * on the input alignment (in string format)
 * model_set -- a set of models to consider
 * freq_set -- a set of frequency types
 * rate_set -- a set of RHAS models
 * rand_seed -- random seed, if 0, then will generate a new random seed
 * num_thres -- number of cpu threads to be used, default: 1; 0 - auto detection of the optimal number of cpu threads
 * output: modelfinder results in YAML format
 */
extern "C" StringResult modelfinder(StringArray& names, StringArray& seqs, int rand_seed = 0,
                   const char* model_set = "", const char* freq_set = "", const char* rate_set = "", int num_thres = 1, const char* other_options = NULL);

/*
 * Build pairwise JC distance matrix
 * output: set of distances
 * (n * i + j)-th element of the list represents the distance between i-th and j-th sequence,
 * where n is the number of sequences
 * num_thres -- number of cpu threads to be used, default: 1; 0 - use all available cpu threads on the machine
 */
extern "C" DoubleArrayResult build_distmatrix(StringArray& names, StringArray& seqs, int num_thres = 1);

/*
 * Using Rapid-NJ to build tree from a distance matrix
 * output: a newick tree (in string format)
 */
extern "C" StringResult build_njtree(StringArray& names, DoubleArray& distances);

/*
 * Compute a consensus tree
 * trees -- a set of input trees
 * minsup -- a threshold to build the majority consensus, default is 0.0
 * output: the consensus tree of the set of input trees
 */
extern "C" StringResult consensus_tree(StringArray& trees, double minsup = 0.0);

/*
 * verion number
 */
extern "C" StringResult version();

/*
 * Execute AliSim Simulation
 * output: results in YAML format that contains the simulated alignment and the content of the log file
 * tree -- the NEWICK tree string
 * subst_model -- the substitution model name
 * seed -- the random seed
 * partition_info -- partition information
 * partition_type -- partition type is either ‘equal’, ‘proportion’, or ‘unlinked’
 * seq_length -- the length of sequences
 * insertion_rate -- the insertion rate
 * deletion_rate -- the deletion rate
 * root_seq -- the root sequence
 * num_threads -- the number of threads
 * insertion_size_distribution -- the insertion size distribution
 * deletion_size_distribution -- the deletion size distribution
 * population_size -- the population size
 */
extern "C" StringResult simulate_alignment(const char* tree, const char* subst_model, int seed, const char* partition_info = "", const char* partition_type = "", int seq_length = 1000, double insertion_rate = 0, double deletion_rate = 0, const char* root_seq = "", int num_threads = 1, const char* insertion_size_distribution = "", const char* deletion_size_distribution = "", int population_size = -1);


/*
 * free the pointer
 */
extern "C" void iqtree_free(void *p);

#endif /* LIBIQTREE2_FUN */
