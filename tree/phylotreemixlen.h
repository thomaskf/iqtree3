//
//  phylotreemixlen.h
//  iqtree
//
//  Created by Minh Bui on 24/08/15.
//
//

#ifndef __iqtree__phylotreemixlen__
#define __iqtree__phylotreemixlen__

#include <stdio.h>
#ifdef USE_CPPOPTLIB
#include "cppoptlib/meta.h"
#include "cppoptlib/boundedproblem.h"
#endif
#include "iqtree.h"



/**
    Phylogenetic tree with mixture of branch lengths
    Started within joint project with Stephen Crotty
*/
#ifdef USE_CPPOPTLIB
class PhyloTreeMixlen : public IQTree, public cppoptlib::BoundedProblem<double>
#else
class PhyloTreeMixlen : public IQTree
#endif
{

    friend class ModelFactoryMixlen;

public:

    /**
            default constructor
     */
    PhyloTreeMixlen();

    PhyloTreeMixlen(Alignment *aln, int mixlen);

    virtual ~PhyloTreeMixlen() override;

    /**
        start structure for checkpointing
    */
    virtual void startCheckpoint() override;

    /** 
        save object into the checkpoint
    */
    virtual void saveCheckpoint() override;

    /** 
        restore object from the checkpoint
    */
    virtual void restoreCheckpoint() override;

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
            refactored 2015-12-22: Taxon IDs instead of Taxon names to save space!
            Read the tree saved with Taxon IDs and branch lengths.
            @param tree_string tree string to read from
     */
    virtual void readTreeString(const string &tree_string) override;

    virtual void initializeModel(Params &params, string model_name, ModelsBlock *models_block) override;

    /**
        @return true if this is a tree with mixture branch lengths, default: false
    */
    virtual bool isMixlen() override { return !initializing_mixlen; }

    /**
        @return number of mixture branch lengths, default: 1
    */
    virtual int getMixlen() override {
        if (initializing_mixlen)
            return 1;
        else
            return mixlen;
    }

    /**
        set number of mixture branch lengths
    */
    void setMixlen(int mixlen);

    /**
            @param[out] lenvec tree lengths for each class in mixlen model
            @param node the starting node, nullptr to start from the root
            @param dad dad of the node, used to direct the search
     */
    virtual void treeLengths(DoubleVector &lenvec, Node *node = nullptr, Node *dad = nullptr) override;

    /**
     * assign branch length as mean over all branch lengths of categories
     */
    void assignMeanMixBranches(Node *node = nullptr, Node *dad = nullptr);


    /**
        parse the string containing branch length(s)
        by default, this will parse just one length
        @param lenstr string containing branch length(s)
        @param[out] branch_len output branch length(s)
    */
//    virtual void parseBranchLength(string &lenstr, DoubleVector &branch_len);

    /**
     *  internal function called by printTree to print branch length
     *  @param out output stream
     *  @param length_nei target Neighbor to print
     */
    virtual void printBranchLength(ostream &out, int brtype, bool print_slash, Neighbor *length_nei) override;

    /**
            print tree to .treefile
            @param suffix suffix of the output file
     */
    virtual void printResultTree(string suffix = "") override;

    /**
        initialize mixture branch lengths
    */
    void initializeMixBranches(PhyloNode *node = nullptr, PhyloNode *dad = nullptr);

    /** initialize parameters if necessary */
    void initializeMixlen(double tolerance, bool write_info);

    /**
        called by fixNegativeBranch to fix one branch
        @param branch_length new branch length
        @param dad_branch dad branch
        @param dad dad node
    */
    virtual void fixOneNegativeBranch(double branch_length, Neighbor *dad_branch, Node *dad) override;

    /**
     * IMPORTANT: semantic change: this function does not return score anymore, for efficiency purpose
            optimize one branch length by ML
            @param node1 1st end node of the branch
            @param node2 2nd end node of the branch
            @param clearLH true to clear the partial likelihood, otherwise false
            @param maxNRStep maximum number of Newton-Raphson steps
     */
    virtual void optimizeOneBranch(PhyloNode *node1, PhyloNode *node2, bool clearLH = true, int maxNRStep = 100) override;

    /**
            optimize all branch lengths of the tree
            @param my_iterations number of iterations to loop through all branches
            @return the likelihood of the tree
     */
    virtual double optimizeAllBranches(int my_iterations = 100, double tolerance = TOL_LIKELIHOOD, int maxNRStep = 100) override;

	/**
		This function calculate f(value), first derivative f'(value) and 2nd derivative f''(value).
		used by Newton raphson method to minimize the function.
		Please always override this function to adapt to likelihood or parsimony score.
		The default is for function f(x) = x^2.
		@param value x-value of the function
		@param df (OUT) first derivative
		@param ddf (OUT) second derivative
	*/
	virtual void computeFuncDervMulti(double *value, double *df, double *ddf) override;

    /**
            Inherited from Optimization class.
            This function calculate f(value), first derivative f'(value) and 2nd derivative f''(value).
            used by Newton raphson method to minimize the function.
            @param value current branch length
            @param df (OUT) first derivative
            @param ddf (OUT) second derivative
     */
    virtual void computeFuncDerv(double value, double &df, double &ddf) override;

	/**
		return the number of dimensions
	*/
	virtual int getNDim() override;


	/**
		the target function which needs to be optimized
		@param x the input vector x
		@return the function value at x
	*/
	virtual double targetFunk(double x[]) override;

	/**
		the approximated derivative function
		@param x the input vector x
		@param dfx the derivative at x
		@return the function value at x
	*/
	virtual double derivativeFunk(double x[], double dfx[]) override;

    /** For Mixlen stuffs */
    virtual int getCurMixture() override { return cur_mixture; }

    /**
     *  Optimize current tree using NNI
     *
     *  @return
     *      <number of NNI steps, number of NNIs> done
     */
    virtual pair<int, int> optimizeNNI(bool speedNNI = true) override;

    /** number of mixture categories */
    int mixlen;

    /** current category, for optimizing branch length */
    int cur_mixture;


/*************** Using cppoptlib for branch length optimization ***********/

#ifdef USE_CPPOPTLIB

//    using typename BoundedProblem<double>::TVector;
//    using typename BoundedProblem<double>::THessian;

    /**
    * @brief returns objective value in x
    * @details [long description]
    *
    * @param x [description]
    * @return [description]
    */
    double value(const TVector &x);

    /**
    * @brief returns gradient in x as reference parameter
    * @details should be overwritten by symbolic gradient
    *
    * @param grad [description]
    */
    void gradient(const TVector &x, TVector &grad);

    /**
    * @brief This computes the hessian
    * @details should be overwritten by symbolic hessian, if solver relies on hessian
    */
    void hessian(const TVector &x, THessian &hessian);

#endif

    /**
     * clear the array "relative_treelen"
     */
    void clear_relative_treelen();

protected:
    
    /** relative rate, used to initialize branch lengths */
    DoubleVector relative_treelen;

    /** true if during initialization phase */
    bool initializing_mixlen;

};

#endif /* defined(__iqtree__phylotreemixlen__) */
