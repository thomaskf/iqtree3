/***************************************************************************
 *   Copyright (C) 2009 by BUI Quang Minh   *
 *   minh.bui@univie.ac.at   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef PHYLOSUPERTREE_H
#define PHYLOSUPERTREE_H

#include "iqtree.h"
#include "supernode.h"
#include "alignment/superalignment.h"

/**
    Partition information for each partition
 */
struct PartitionInfo {
    double cur_score;    // current log-likelihood
    double part_rate;    // partition heterogeneity rate
    int    evalNNIs;    // number of evaluated NNIs on subtree
    
    //DoubleVector null_score; // log-likelihood of each branch collapsed to zero
    //DoubleVector opt_score;  // optimized log-likelihood for every branch
    //DoubleVector nni1_score; // log-likelihood for 1st NNI for every branch
    //DoubleVector nni2_score; // log-likelihood for 2nd NNI for every branch
    
    vector<DoubleVector> cur_brlen;  // current branch lengths
    //DoubleVector opt_brlen;  // optimized branch lengths for every branch
    vector<DoubleVector> nni1_brlen; // branch length for 1st NNI for every branch
    vector<DoubleVector> nni2_brlen; // branch length for 2nd NNI for every branch
    
    //double *mem_ptnlh; // total memory allocated for all pattern likelihood vectors
    double *cur_ptnlh; // current pattern likelihoods of the tree
    //double *nni1_ptnlh; // pattern likelihoods of 1st NNI tree
    //double *nni2_ptnlh; // pattern likelihoods of 2nd NNI tree
    NNIMove nniMoves[2];
};

/**
Phylogenetic tree for partition model (multi-gene alignment)

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class PhyloSuperTree : public IQTree, public vector<PhyloTree* >
{
public:
	/**
		constructor
	*/
    PhyloSuperTree();

    /**
     constructor
     */
    PhyloSuperTree(SuperAlignment *alignment, bool new_iqtree = false, bool create_tree = true);

    /**
		constructor
	*/
    PhyloSuperTree(SuperAlignment *alignment, PhyloSuperTree *super_tree);

    /**
        destructor
    */
    ~PhyloSuperTree() override;

    /**
        set checkpoint object
        @param checkpoint a checkpoint
    */
    virtual void setCheckpoint(Checkpoint *checkpoint) override;

    /**
        save object into the checkpoint
    */
    virtual void saveCheckpoint() override;

    /**
        restore object from the checkpoint
    */
    virtual void restoreCheckpoint() override;

    /**
            set the model factory
            @param model_fac model factory
     */
    virtual void setModelFactory(ModelFactory *model_fac) override;

    /**
     2019-06-03: copy part_info from tree, taking into account -bsam option
     @param tree input super tree
     */
    void setPartInfo(PhyloSuperTree *tree);
    
    /**
            Set the alignment, important to compute parsimony or likelihood score
            Assing taxa ids according to their position in the alignment
            @param alignment associated alignment
     */
    // virtual void setSuperAlignment(Alignment *alignment);

    /** remove identical sequences from the tree */
    virtual void removeIdenticalSeqs(Params &params) override;

    /** reinsert identical sequences into the tree and reset original alignment */
    virtual void reinsertIdenticalSeqs(Alignment *orig_aln) override;

	virtual void setParams(Params* params) override;

	/**
	 * setup all necessary parameters  (declared as virtual needed for phylosupertree)
	 */
	virtual void initSettings(Params& params) override;
    
    /**
        init sequence instances for each nodes when using AliSim
    */
    virtual void initSequences(Node* node = nullptr, Node* dad = nullptr) override;

    virtual void setLikelihoodKernel(LikelihoodKernel lk) override;

    virtual void changeLikelihoodKernel(LikelihoodKernel lk) override;

    virtual void setParsimonyKernel(LikelihoodKernel lk) override;

    virtual void setNumThreads(int num_threads) override;

	virtual bool isSuperTree() override { return true; }

    /**
     print tree to .treefile
     @param suffix suffix of the output file
     */
    virtual void printResultTree(string suffix = "") override;

    /**
     * Return the tree string contining taxon names and branch lengths
     * @return a tree string
     */
    virtual string getTreeString() override;

    /**
            Read the tree saved with Taxon Names and branch lengths.
            @param tree_string tree string to read from
     */
    virtual void readTreeString(const string &tree_string) override;

    /**
     * save branch lengths into a vector
     */
    virtual void saveBranchLengths(DoubleVector &lenvec, int startid = 0, PhyloNode *node = nullptr, PhyloNode *dad = nullptr) override;
    /**
     * restore branch lengths from a vector previously called with saveBranchLengths
     */
    virtual void restoreBranchLengths(DoubleVector &lenvec, int startid = 0, PhyloNode *node = nullptr, PhyloNode *dad = nullptr) override;

    /**
        Collapse all internal branches with length <= threshold
		@param node the starting node, nullptr to start from the root
		@param dad dad of the node, used to direct the search
        @param threshold branch length threshold
        @return number of branches collapsed
    */
    virtual int collapseInternalBranches(Node *node = nullptr, Node *dad = nullptr, double threshold = 0.0) override;

    /**
            allocate a new node. Override this if you have an inherited Node class.
            @param node_id node ID
            @param node_name node name
            @return a new node
     */
    virtual Node* newNode(int node_id = -1, const char* node_name = nullptr) override;

    /**
            allocate a new node. Override this if you have an inherited Node class.
            @param node_id node ID
            @param node_name node name issued by an interger
            @return a new node
     */
    virtual Node* newNode(int node_id, int node_name) override;

	/**
	 *		@return number of alignment patterns
	*/
	virtual size_t getAlnNPattern() override;

	/**
	 *		@return number of alignment sites
	*/
	virtual size_t getAlnNSite() override;

    /**
            compute the distance between 2 sequences.
            @param seq1 index of sequence 1
            @param seq2 index of sequence 2
            @param initial_dist initial distance
            @return distance between seq1 and seq2
     */
    virtual double computeDist(int seq1, int seq2, double initial_dist, double &var) override;

	/**
		create sub-trees T|Y_1,...,T|Y_k of the current super-tree T
		and map F={f_1,...,f_k} the edges of supertree T to edges of subtrees T|Y_i
	*/
	virtual void mapTrees();

	/*
	 * create one map f_i from supertree T to subtree indexed by part (called by mapTrees)
	 * @param part index of subtree
	 * @param part_taxa vector of taxa of T that are present in subtree
	 * @param node the current node of the post-order tree traversal
	 * @param dad the dad of that node used to direct the traversal
	 */
	void linkTree(int part, NodeVector &part_taxa, SuperNode *node = nullptr, SuperNode *dad = nullptr);

	/**
	 * Given current supertree T and subtrees T|Y_1,...,T|Y_k, build all maps f_1,...,f_k
	 */
	virtual void linkTrees();

	/**
	 * link a branch from supertree to subtree (called by linkTree)
	 * @param part index of subtree
	 * @param nei pointer to branch
	 * @param dad_nei pointer to reverse branch
	 */
	void linkBranch(int part, SuperNeighbor *nei, SuperNeighbor *dad_nei);

    /**
        make the rooting consistent between trees
     */
    void syncRooting();
    
    
    /**
            initialize partial_lh vector of all PhyloNeighbors, allocating central_partial_lh
     */
    virtual void initializeAllPartialLh() override;

    /**
            de-allocate central_partial_lh
     */
    virtual void deleteAllPartialLh() override;

    /**
     NEWLY ADDED (2014-12-04): clear all partial likelihood for a clean computation again
     */
    virtual void clearAllPartialLH(bool make_null = false) override;
    

    /**
            compute the tree likelihood
            @param pattern_lh (OUT) if not nullptr, the function will assign pattern log-likelihoods to this vector
                            assuming pattern_lh has the size of the number of patterns
            @return tree likelihood
     */
    virtual double computeLikelihood(double *pattern_lh = nullptr, bool save_log_value = true) override;

    /**
     * @return number of elements per site lhl entry, used in conjunction with computePatternLhCat
     */
    virtual int getNumLhCat(SiteLoglType wsl) override;

    /**
            compute pattern likelihoods only if the accumulated scaling factor is non-zero.
            Otherwise, copy the pattern_lh attribute
            @param pattern_lh (OUT) pattern log-likelihoods,
                            assuming pattern_lh has the size of the number of patterns
            @param cur_logl current log-likelihood (for sanity check)
            @param pattern_lh_cat (OUT) if not nullptr, store all pattern-likelihood per category
     */
    virtual void computePatternLikelihood(double *pattern_lh, double *cur_logl = nullptr,
    		double *pattern_lh_cat = nullptr, SiteLoglType wsl = WSL_RATECAT) override;

    /**
            compute pattern posterior probabilities per rate/mixture category
            @param pattern_prob_cat (OUT) all pattern-probabilities per category
            @param wsl either WSL_RATECAT, WSL_MIXTURE or WSL_MIXTURE_RATECAT
     */
    virtual void computePatternProbabilityCategory(double *pattern_prob_cat, SiteLoglType wsl) override;

    /**
            optimize all branch lengths of all subtrees, then compute branch lengths
            of supertree as weighted average over all subtrees
            @param my_iterations number of iterations to loop through all branches
            @return the likelihood of the tree
     */
    virtual double optimizeAllBranches(int my_iterations = 100, double tolerance = TOL_LIKELIHOOD, int maxNRStep = 100) override;

    /**
            search the best swap for a branch
            @return NNIMove The best Move/Swap
            @param node1 1 of the 2 nodes on the branch
            @param node2 1 of the 2 nodes on the branch
     */
    virtual NNIMove getBestNNIForBran(PhyloNode *node1, PhyloNode *node2, NNIMove *nniMoves = nullptr) override;

    /**
            Do an NNI on the supertree and synchronize all subtrees respectively
            @param move the single NNI
     */
    virtual void doNNI(NNIMove &move, bool clearLH = true) override;

    /**
     *   Apply 5 new branch lengths stored in the NNI move
     *   @param nnimove the NNI move currently in consideration
     */
    virtual void changeNNIBrans(NNIMove &nnimove) override;

    /**
        OBSOLETE!
     * 	 Restore the branch lengths from the saved values
	 * @param node the current node of the post-order tree traversal
	 * @param dad the dad of that node used to direct the traversal
     */
//    virtual void restoreAllBrans(PhyloNode *node = nullptr, PhyloNode *dad = nullptr);

    /**
            reinsert the whole list of leaves back into the supertree then call mapTrees
            @param del_leaves the list of deleted leaves, returned by deleteLeaves() function
     */
    virtual void reinsertLeaves(PhyloNodeVector &del_leaves) override;

	/**
		compute the weighted average of branch lengths over partitions
	*/
	virtual void computeBranchLengths();

	/**
	 * print debug information about all maps
	 */
	virtual void printMapInfo();

	/**
	 * initialize partition information for super tree
	*/
	virtual void initPartitionInfo() override;

	int getMaxPartNameLength();

	/**
		partition information
	*/
	vector<PartitionInfo> part_info;

    /* partition ID sorted in descending order of computation cost */
    IntVector part_order;
    IntVector part_order_by_nptn;

    /* compute part_order vector */
    void computePartitionOrder();

    /**
            get the name of the model
    */
    virtual string getModelName() override;
	/**
	 * extract subtree containing all taxa from partition IDs
	 * @param ids partitions IDs
	 * @return subtree
	 */
    PhyloTree *extractSubtree(set<int> &ids);

    /**
     * compute the memory size required for storing partial likelihood vectors
     * @return memory size required in bytes
     */
    virtual uint64_t getMemoryRequired(size_t ncategory = 1, bool full_mem = false) override;

    /**
     * compute the memory size for top partitions required for storing partial likelihood vectors
     * @return memory size required in bytes
     */
    virtual uint64_t getMemoryRequiredThreaded(size_t ncategory = 1, bool full_mem = false) override;

    /**
     * count the number of super branches that map to no branches in gene trees
     */
    int countEmptyBranches(PhyloNode *node = nullptr, PhyloNode *dad = nullptr);

    /**
            Neighbor-joining/parsimony tree might contain negative branch length. This
            function will fix this.
            @param node the current node
            @param dad dad of the node, used to direct the search
            @return The number of branches that have no/negative length
     */
    virtual int fixNegativeBranch(bool force = false, Node *node = nullptr, Node *dad = nullptr) override;

    // virtual int computeParsimonyBranchObsolete(PhyloNeighbor *dad_branch, PhyloNode *dad, int *branch_subst = nullptr);

    /****************************************************************************
            ancestral sequence reconstruction
     ****************************************************************************/

    /**
        initialize computing ancestral sequence probability for an internal node by marginal reconstruction
    */
    virtual void initMarginalAncestralState(ostream &out, bool &orig_kernel_nonrev, double* &ptn_ancestral_prob, int* &ptn_ancestral_seq) override;

    /**
        compute ancestral sequence probability for an internal node by marginal reconstruction
        (Yang, Kumar and Nei 1995)
        @param dad_branch branch leading to an internal node where to obtain ancestral sequence
        @param dad dad of the target internal node
        @param[out] ptn_ancestral_prob pattern ancestral probability vector of dad_branch->node
    */
    virtual void computeMarginalAncestralState(PhyloNeighbor *dad_branch, PhyloNode *dad,
        double *ptn_ancestral_prob, int *ptn_ancestral_seq) override;

    /**
        compute ancestral sequence probability for an internal node by partial_lh of that subtree
        @param dad_branch branch leading to an internal node where to obtain ancestral sequence
        @param dad dad of the target internal node
        @param[out] ptn_ancestral_prob pattern ancestral probability vector of dad_branch->node
    */
    virtual void computeSubtreeAncestralState(PhyloNeighbor *dad_branch, PhyloNode *dad,
        double *ptn_ancestral_prob, int *ptn_ancestral_seq) override;

    virtual void writeMarginalAncestralState(ostream &out, PhyloNode *node, double *ptn_ancestral_prob, int *ptn_ancestral_seq) override;

    /**
        end computing ancestral sequence probability for an internal node by marginal reconstruction
    */
    virtual void endMarginalAncestralState(bool orig_kernel_nonrev, double* &ptn_ancestral_prob, int* &ptn_ancestral_seq) override;
    
	/**
		write site-rates to a file in the following format:
		1  rate_1
		2  rate_2
		....
		This function will call computePatternRates()
		@param out output stream to write rates
        @param bayes TRUE to use empirical Bayesian, false for ML method
     */
    virtual void writeSiteRates(ostream &out, bool bayes, int partid = -1) override;

    /**
        write site log likelihood to a output stream
        @param out output stream
        @param wsl write site-loglikelihood type
        @param partid partition ID as first column of the line. -1 to omit it
    */
    virtual void writeSiteLh(ostream &out, SiteLoglType wsl, int partid = -1) override;


    virtual void writeBranch(ostream &out, Node* node1, Node* node2);

    /**
     write branches into a csv file
     Feature requested by Rob Lanfear
     @param out output stream
     */
    virtual void writeBranches(ostream &out) override;

    /**
        print partition file with best model parameters
        @param filename output file name
     */
    void printBestPartitionParams(const char *filename);

    
    /** True when mixed codon with other data type */
    bool rescale_codon_brlen;
    
    int totalNNIs, evalNNIs;

};

#endif
