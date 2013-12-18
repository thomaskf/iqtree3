#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "pll/pll.h"
#define GLOBAL_VARIABLES_DEFINITION
#include "nnisearch.h"

#if !defined WIN32 && !defined _WIN32 && !defined __WIN32__
#include <sys/resource.h>
#endif

/* program options */
double TOL_LIKELIHOOD_PHYLOLIB;
int numSmoothTree;
int nni0;
int nni5;
extern Params *globalParam;
/* program options */

/*
 * ****************************************************************************
 * pllUFBoot area
 * ****************************************************************************
 */

pllUFBootData * pllUFBootDataPtr = NULL;


int compareDouble(const void * a, const void * b) {
	if (*(double*) a > *(double*) b)
		return 1;
	else if (*(double*) a < *(double*) b)
		return -1;
	else
		return 0;
}

pllNNIMove *getNNIList(pllInstance* tr) {
	static pllNNIMove* nniList;
	if (nniList == NULL) {
		nniList = (pllNNIMove*) malloc(2 * (tr->mxtips - 3) * sizeof(pllNNIMove));
		assert(nniList != NULL);
	}
	return nniList;
}

pllNNIMove *getNonConflictNNIList(pllInstance* tr) {
	static pllNNIMove* nonConfNNIList;
	if (nonConfNNIList == NULL) {
		nonConfNNIList = (pllNNIMove*) malloc((tr->mxtips - 3) * sizeof(pllNNIMove));
		assert(nonConfNNIList != NULL);
	}
	return nonConfNNIList;
}

double pllDoRandNNIs(pllInstance *tr, partitionList *pr, int numNNI) {
	int numInBrans = tr->mxtips - 3;
	int numNNIinStep = (int) numInBrans / 5;

	// devided in multiple round, each round collect 1/5 of numNNI
	int cnt1 = 0;
	unordered_set<int> selectedNodes;
	vector<nodeptr> selectedBrans;
	vector<nodeptr> brans;
	do {
		int cnt2 = 0;
		selectedNodes.clear();
		selectedBrans.clear();
		brans.clear();
		pllGetAllInBran(tr, brans);
		assert(brans.size() == numInBrans);
		while(cnt2 < numNNIinStep && cnt2 < numNNI) {
			int branIndex = random_int(numInBrans);
			if (selectedNodes.find(brans[branIndex]->number) == selectedNodes.end() &&
					selectedNodes.find(brans[branIndex]->back->number) == selectedNodes.end()) {
				selectedNodes.insert(brans[branIndex]->number);
				selectedNodes.insert(brans[branIndex]->back->number);
				selectedBrans.push_back(brans[branIndex]);
				cnt2++;
			}
		}
		for (vector<nodeptr>::iterator it = selectedBrans.begin(); it != selectedBrans.end(); ++it) {
			int nniType = random_int(2);
			doOneNNI(tr, pr, (*it), nniType, TOPO_ONLY);
		}
		cnt1 += selectedBrans.size();
		if (numNNI - cnt1 < numNNIinStep) {
			numNNIinStep = numNNI - cnt1;
		}
	} while (cnt1 < numNNI);
	pllEvaluateGeneric(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
	pllTreeEvaluate(tr, pr, 1);
	return tr->likelihood;
}

void pllGetAllInBran(pllInstance *tr, vector<nodeptr> &branlist) {
	nodeptr p = tr->start->back;
	nodeptr q = p->next;
	while (q != p) {
		pllGetAllInBranForSubtree(tr, q->back, branlist);
		q = q->next;
	}
}

void pllGetAllInBranForSubtree(pllInstance *tr, nodeptr p, vector<nodeptr> &branlist) {
	if (!isTip(p->number, tr->mxtips) && !isTip(p->back->number, tr->mxtips)) {
		branlist.push_back(p);
		nodeptr q = p->next;
		while (q != p) {
			pllGetAllInBranForSubtree(tr, q->back, branlist);
			q = q->next;
		}
	}
}

double pllPerturbTree(pllInstance *tr, partitionList *pr, vector<pllNNIMove> &nnis) {
	//printf("Perturbing %d NNIs \n", numNNI);
	for (vector<pllNNIMove>::iterator it = nnis.begin(); it != nnis.end(); it++) {
		printf("Do pertubing NNI (%d - %d) with logl = %10.4f \n", (*it).p->number, (*it).p->back->number, (*it).likelihood);
		doOneNNI(tr, pr, (*it).p, (*it).nniType, TOPO_ONLY);
		updateBranchLengthForNNI(tr, pr, (*it));

	}
	pllEvaluateGeneric(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
	pllTreeEvaluate(tr, pr, 1);
	return tr->likelihood;
}

void quicksort_nni(pllNNIMove* arr,int left, int right) {
    int i = left, j = right;
    pllNNIMove tmp, pivot = arr[(left + right) / 2];

    /* partition */
    while (i <= j) {
        while (arr[i].likelihood < pivot.likelihood)
            i++;
        while (pivot.likelihood < arr[j].likelihood)
            j--;
        if (i <= j) {
            tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;

            i++;
            j--;
        }
    };

    /* recursion */
    if (left < j)
        quicksort_nni(arr, left, j);
    if (i < right)
        quicksort_nni(arr, i, right);
}

//TODO: Workaround for memory leak problem when calling setupTopol within doNNISearch
topol *_setupTopol(pllInstance* tr) {
	static topol* tree;
	if (tree == NULL)
		tree = setupTopol(tr->mxtips);
	return tree;
}

vector<string> getAffectedBranches(pllInstance* tr, nodeptr p) {
	vector<string> res;
	res.push_back(getBranString(p));
	nodeptr q = p->back;
	nodeptr p_nei = p->next;
	nodeptr q_nei = q->next;
	while (p_nei != p) {
		res.push_back(getBranString(p_nei));
		if (!isTip(p_nei->back->number, tr->mxtips)) {
			res.push_back(getBranString(p_nei->back->next));
			res.push_back(getBranString(p_nei->back->next->next));
		}
		p_nei = p_nei->next;
	}
	while (q_nei != q) {
		res.push_back(getBranString(q_nei));
		if (!isTip(q_nei->back->number, tr->mxtips)) {
			res.push_back(getBranString(q_nei->back->next));
			res.push_back(getBranString(q_nei->back->next->next));
		}
		q_nei = q_nei->next;
	}
	return res;
}

string getBranString(nodeptr p) {
	stringstream branString;
	nodeptr q = p->back;
	if (p->number < q->number) {
		branString << p->number << "-" << q->number;
	} else {
		branString << q->number << "-" << p->number;
	}
	return branString.str();
}

set<int> getAffectedNodes(pllInstance* tr, nodeptr p) {
	set<int> nodeSet;
	nodeptr q = p->back;
	nodeptr p_nei = p->next;
	nodeptr q_nei = q->next;
	nodeSet.insert(p->number);
	nodeSet.insert(q->number);
	while (p_nei != p) {
		nodeptr nei = p_nei->back;
		if (isTip(nei->number, tr->mxtips)) {
			nodeSet.insert(nei->number);
		} else {
			nodeSet.insert(nei->number);
			nodeSet.insert(nei->next->back->number);
			nodeSet.insert(nei->next->next->back->number);
		}
		p_nei = p_nei->next;
	}
	while (q_nei != q) {
		nodeptr nei = q_nei->back;
		if (isTip(nei->number, tr->mxtips)) {
			nodeSet.insert(nei->number);
		} else {
			nodeSet.insert(nei->number);
			nodeSet.insert(nei->next->back->number);
			nodeSet.insert(nei->next->next->back->number);
		}
		q_nei = q_nei->next;
	}
	return nodeSet;
}

void pllEvalAllNNIs(pllInstance *tr, partitionList *pr, SearchInfo &searchinfo) {
	/* DTH: mimic IQTREE::optimizeNNI 's first call to IQTREE::saveCurrentTree */
	if((globalParam->online_bootstrap == PLL_TRUE) &&
			(globalParam->gbo_replicates > 0)){
		tr->fastScaling = PLL_FALSE;
		pllEvaluateGeneric(tr, pr, tr->start, PLL_FALSE, PLL_TRUE);
		pllSaveCurrentTree(tr, pr, tr->start);
	}

	nodeptr p = tr->start->back;
	nodeptr q = p->next;
	while (q != p) {
		evalNNIForSubtree(tr, pr, q->back, searchinfo);
		q = q->next;
		//		int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;
		//		if (numBranches > 1 && !tr->useRecom) {
		//			pllNewviewGeneric(tr, pr,  q->back, PLL_TRUE);
		//		} else {
		//			pllNewviewGeneric(tr, pr,  q->back, PLL_FALSE);
		//		}
	}
}

/*
void pllSaveQuartetForSubTree(pllInstance *tr, nodeptr p, SearchInfo &searchinfo) {
	if (!isTip(p->number, tr->mxtips) && !isTip(p->back->number, tr->mxtips)) {
			evalNNIForBran(tr, pr, p, searchinfo);
		nodeptr q = p->next;
		while (q != p) {
			pllSaveQuartetForSubTree(tr, q->back, searchinfo);
			q = q->next;
		}
	}
}

void pllSaveAllQuartet(pllInstance *tr, SearchInfo &searchinfo) {
	nodeptr p = tr->start->back;
	nodeptr q = p->next;
	while (q != p) {
		pllSaveQuartetForSubTree(tr, q->back, searchinfo);
	}
}
*/

double pllDoNNISearch(pllInstance* tr, partitionList *pr, SearchInfo &searchinfo) {
	// update tabu list
	if (searchinfo.tabunni) {
		pllUpdateTabuList(tr, searchinfo);
	}
	//cout << "Tabu size: " << searchinfo.tabuNNIs.size() << endl;
	double initLH = tr->likelihood;
	double finalLH = initLH;
	vector<pllNNIMove> selectedNNIs;
	unordered_set<int> selectedNodes;
	int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;
	/* data structure to store the initial tree topology + branch length */
	topol* curTree = _setupTopol(tr);
	saveTree(tr, curTree, numBranches);
	pllEvalAllNNIs(tr, pr, searchinfo);

	if (searchinfo.posNNIList.size() != 0) {
		sort(searchinfo.posNNIList.begin(), searchinfo.posNNIList.end(), comparePLLNNIMove);
		for (vector<pllNNIMove>::reverse_iterator rit = searchinfo.posNNIList.rbegin();
				rit != searchinfo.posNNIList.rend(); ++rit) {
			if (selectedNodes.find((*rit).p->number) == selectedNodes.end()
					&& selectedNodes.find((*rit).p->back->number) == selectedNodes.end()) {
				selectedNNIs.push_back((*rit));
				selectedNodes.insert((*rit).p->number);
				selectedNodes.insert((*rit).p->back->number);
			} else {
				continue;
			}

//			bool conflict = false;
//			for (vector<pllNNIMove>::iterator it = selectedNNIs.begin(); it != selectedNNIs.end(); ++it) {
//				if ((*rit).p->number == (*it).p->number || (*rit).p->number == (*it).p->back->number) {
//					conflict = true;
//					break;
//				}
//				if ((*rit).p->back->number == (*it).p->number || (*rit).p->back->number == (*it).p->back->number) {
//					conflict = true;
//					break;
//				}
//			}
//			if (!conflict) {
//				selectedNNIs.push_back((*rit));
//			}

		}
		//searchinfo.affectNodes.clear();
		//searchinfo.affectBranches.clear();
		/* Applying all independent NNI moves */
		searchinfo.curNumAppliedNNIs = selectedNNIs.size();
		for (vector<pllNNIMove>::iterator it = selectedNNIs.begin(); it != selectedNNIs.end(); it++) {
			/* do the topological change */
			doOneNNI(tr, pr, (*it).p, (*it).nniType, TOPO_ONLY);

			/*
			set<int> affectedNodes = getAffectedNodes(tr, (*it).p);
			searchinfo.affectNodes.insert(affectedNodes.begin(),
					affectedNodes.end());
			vector<string> aBranches = getAffectedBranches(tr, (*it).p);
			searchinfo.affectBranches.insert(aBranches.begin(),
					aBranches.end());
			set<int>::iterator iter2;
			for (iter2 = affectedNodes.begin(); iter2 != affectedNodes.end();
					iter2++) {
				cout << " " << *iter2;
			}
			cout << endl;
			string quartetString = convertQuartet2String((*it).p);
			cout << quartetString << " / old logl: " << initLH
					<< "/ logl delta: " << (*it).loglDelta << endl;
			if (searchinfo.tabuNNIs.find(quartetString)
					!= searchinfo.tabuNNIs.end()) {
				cout << "Quartet " << quartetString
						<< " has already been used / Skipped" << endl;
			} else {
				searchinfo.tabuNNIs.insert(
						unordered_map<string, pllNNIMove>::value_type(
								quartetString, (*it)));
			}
			unordered_map<string, pllNNIMove>::iterator iter1 =
					searchinfo.negativeQuartets.find(quartetString);
			if (iter1 != searchinfo.negativeQuartets.end()) {
				cout << "A negative quartet applied: " << quartetString
						<< " / old delta: " << iter1->second.loglDelta
						<< " / new delta: " << selectedNNIs[i].loglDelta
						<< endl;
				set<int> affectedNodes = getAffectedNodes(tr,
						selectedNNIs[i].p);
				set<int>::iterator iter2;
				for (iter2 = affectedNodes.begin();
						iter2 != affectedNodes.end(); iter2++) {
					cout << " " << *iter2;
				}
				cout << endl;
			}
			*/

			updateBranchLengthForNNI(tr, pr, (*it));
		}
		if (selectedNNIs.size() != 0) {
			pllEvaluateGeneric(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
			pllTreeEvaluate(tr, pr, 1);
			/* new tree likelihood should not be smaller the likelihood of the computed best NNI */
			if (tr->likelihood < selectedNNIs.back().likelihood) {
				if (selectedNNIs.size() == 1) {
					printf("ERROR: new logl=%10.4f after applying only the best NNI < best NNI logl=%10.4f\n",
							tr->likelihood, selectedNNIs[0].likelihood);
					exit(1);
				} else {
					cout << "Roll back tree ... " << endl;
					if (!restoreTree(curTree, tr, pr)) {
						printf("ERROR: failed to roll back tree \n");
						exit(1);
					}
					doOneNNI(tr, pr, selectedNNIs.back().p, selectedNNIs.back().nniType, TOPO_ONLY);
					updateBranchLengthForNNI(tr, pr, selectedNNIs.back());
					pllEvaluateGeneric(tr, pr, tr->start, PLL_TRUE, PLL_FALSE);
					pllTreeEvaluate(tr, pr, 1);
					if (tr->likelihood < selectedNNIs.back().likelihood) {
						printf("ERROR: (After rolling back) new logl=%10.4f after applying only the best NNI < best NNI logl=%10.4f\n",
								tr->likelihood, selectedNNIs.front().likelihood);
						exit(1);
					}

					/* Only apply the best NNI after the tree has been rolled back */
					searchinfo.curNumAppliedNNIs = 1;
				}
			}
			if (tr->likelihood - initLH < 0.1) {
				searchinfo.curNumAppliedNNIs = 0;
				// Here we need to update searchinfo.nniList;
				searchinfo.updateNNIList = true;
			}
			finalLH = tr->likelihood;
		}
	} else {
		searchinfo.curNumAppliedNNIs = 0;
	}
//	if (verbose_mode >= VB_MED) {
//		cout << "Step: " << searchinfo.curNumNNISteps << " / NumNNI: " << searchinfo.curNumAppliedNNIs << " / logl: " << finalLH << endl;
//	}
	return finalLH;
}

void pllUpdateTabuList(pllInstance *tr, SearchInfo &searchinfo) {
	nodeptr p = tr->start->back;
	nodeptr q = p->next;
	while (q != p) {
		pllSaveQuartetForSubTree(tr, q->back, searchinfo);
		q = q->next;
	}
}

void pllSaveQuartetForSubTree(pllInstance* tr, nodeptr p, SearchInfo &searchinfo) {
	if (!isTip(p->number, tr->mxtips) && !isTip(p->back->number, tr->mxtips)) {
		pllSaveQuartet(p, searchinfo);
		nodeptr q = p->next;
		while (q != p) {
			pllSaveQuartetForSubTree(tr, q->back, searchinfo);
			q = q->next;
		}
	}
}

void updateBranchLengthForNNI(pllInstance* tr, partitionList *pr, pllNNIMove &nni) {
	int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;
	/*  apply branch lengths */
	for (int i = 0; i < numBranches; i++) {
		nni.p->z[i] = nni.z0[i];
		nni.p->back->z[i] = nni.z0[i];
		nni.p->next->z[i] = nni.z1[i];
		nni.p->next->back->z[i] = nni.z1[i];
		nni.p->next->next->z[i] = nni.z2[i];
		nni.p->next->next->back->z[i] = nni.z2[i];
		nni.p->back->next->z[i] = nni.z3[i];
		nni.p->back->next->back->z[i] = nni.z3[i];
		nni.p->back->next->next->z[i] = nni.z4[i];
		nni.p->back->next->next->back->z[i] = nni.z4[i];
	}
	/* update partial likelihood */
//	if (numBranches > 1 && !tr->useRecom) {
//		pllNewviewGeneric(tr, pr, nni.p, PLL_TRUE);
//		pllNewviewGeneric(tr, pr, nni.p->back, PLL_TRUE);
//	} else {
//		pllNewviewGeneric(tr, pr, nni.p, PLL_FALSE);
//		pllNewviewGeneric(tr, pr, nni.p->back, PLL_FALSE);
//	}
}

/** @brief Optimize the length of a specific branch with variable number of Newton Raphson iterations

 Optimize the length of the branch connecting \a p and \a p->back
 for each partition (\a tr->numBranches) in library instance \a tr.

 @param tr
 The library instance

 @param pr
 Partition list

 @param p
 Endpoints of branch to be optimized

 @param newzpercycle Maximal number of Newton Raphson iterations

 */
void _update(pllInstance *tr, partitionList *pr, nodeptr p) {
	nodeptr q;
	int i;
	double z[PLL_NUM_BRANCHES], z0[PLL_NUM_BRANCHES];
	int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;

	q = p->back;

	for (i = 0; i < numBranches; i++)
		z0[i] = q->z[i];

	if (numBranches > 1)
		makenewzGeneric(tr, pr, p, q, z0, IQTREE_NEWZPERCYCLE, z, PLL_TRUE);
	else
		makenewzGeneric(tr, pr, p, q, z0, IQTREE_NEWZPERCYCLE, z, PLL_FALSE);

	for (i = 0; i < numBranches; i++) {
		if (!tr->partitionConverged[i]) {
			if (PLL_ABS(z[i] - z0[i]) > PLL_DELTAZ) {
				tr->partitionSmoothed[i] = PLL_FALSE;
			}

			p->z[i] = q->z[i] = z[i];
		}
	}
}

double doOneNNI(pllInstance *tr, partitionList *pr, nodeptr p, int swap, int evalType) {
	assert(swap == 0 || swap == 1);
	nodeptr q;
	nodeptr tmp;
	q = p->back;
	assert(!isTip(q->number, tr->mxtips));
	assert(!isTip(p->number, tr->mxtips));
	int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;

	if (swap == 1) {
		tmp = p->next->back;
		hookup(p->next, q->next->back, q->next->z, numBranches);
		hookup(q->next, tmp, tmp->z, numBranches);
	} else {
		tmp = p->next->next->back;
		hookup(p->next->next, q->next->back, q->next->z, numBranches);
		hookup(q->next, tmp, tmp->z, numBranches);
	}

	if (evalType == TOPO_ONLY) {
		return 0.0;
	} else if (evalType == ONE_BRAN_OPT) {
		if (numBranches > 1 && !tr->useRecom) {
			pllNewviewGeneric(tr, pr, p, PLL_TRUE);
			pllNewviewGeneric(tr, pr, q, PLL_TRUE);

		} else {
			pllNewviewGeneric(tr, pr, p, PLL_FALSE);
			pllNewviewGeneric(tr, pr, q, PLL_FALSE);
		}
		_update(tr, pr, p);
		if((globalParam->online_bootstrap == PLL_TRUE) &&
				(globalParam->gbo_replicates > 0)){
			tr->fastScaling = PLL_FALSE;
			pllEvaluateGeneric(tr, pr, p, PLL_FALSE, PLL_TRUE); // DTH: modified the last arg
			pllSaveCurrentTree(tr, pr, p);
		}else{
			pllEvaluateGeneric(tr, pr, p, PLL_FALSE, PLL_FALSE);
		}
	} else if (evalType == NO_BRAN_OPT) {
		if (numBranches > 1 && !tr->useRecom) {
			pllNewviewGeneric(tr, pr, p, PLL_TRUE);
			pllNewviewGeneric(tr, pr, q, PLL_TRUE);

		} else {
			pllNewviewGeneric(tr, pr, p, PLL_FALSE);
			pllNewviewGeneric(tr, pr, q, PLL_FALSE);
		}
		if((globalParam->online_bootstrap == PLL_TRUE) &&
						(globalParam->gbo_replicates > 0)){
			tr->fastScaling = PLL_FALSE;
			pllEvaluateGeneric(tr, pr, p, PLL_FALSE, PLL_TRUE); // DTH: modified the last arg
			pllSaveCurrentTree(tr, pr, p);
		}else{
			pllEvaluateGeneric(tr, pr, p, PLL_FALSE, PLL_FALSE);
		}
	} else { // 5 branches optimization
		if (numBranches > 1 && !tr->useRecom) {
			pllNewviewGeneric(tr, pr, q, PLL_TRUE);
		} else {
			pllNewviewGeneric(tr, pr, q, PLL_FALSE);
		}
		nodeptr r; // temporary node poiter
		r = p->next;
		if (numBranches > 1 && !tr->useRecom)
			pllNewviewGeneric(tr, pr, r, PLL_TRUE);
		else
			pllNewviewGeneric(tr, pr, r, PLL_FALSE);
		_update(tr, pr, r);
		r = p->next->next;
		if (numBranches > 1 && !tr->useRecom)
			pllNewviewGeneric(tr, pr, r, PLL_TRUE);
		else
			pllNewviewGeneric(tr, pr, r, PLL_FALSE);
		_update(tr, pr, r);
		if (numBranches > 1 && !tr->useRecom)
			pllNewviewGeneric(tr, pr, p, PLL_TRUE);
		else
			pllNewviewGeneric(tr, pr, p, PLL_FALSE);
		_update(tr, pr, p);
		// optimize 2 branches at node q
		r = q->next;
		if (numBranches > 1 && !tr->useRecom)
			pllNewviewGeneric(tr, pr, r, PLL_TRUE);
		else
			pllNewviewGeneric(tr, pr, r, PLL_FALSE);
		_update(tr, pr, r);
		r = q->next->next;
		if (numBranches > 1 && !tr->useRecom)
			pllNewviewGeneric(tr, pr, r, PLL_TRUE);
		else
			pllNewviewGeneric(tr, pr, r, PLL_FALSE);
		_update(tr, pr, r);
		if((globalParam->online_bootstrap == PLL_TRUE) &&
						(globalParam->gbo_replicates > 0)){
			tr->fastScaling = PLL_FALSE;
			pllEvaluateGeneric(tr, pr, r, PLL_FALSE, PLL_TRUE); // DTH: modified the last arg
			pllSaveCurrentTree(tr, pr, r);
		}else{
			pllEvaluateGeneric(tr, pr, r, PLL_FALSE, PLL_FALSE);
		}
	}
	return tr->likelihood;

}

string convertQuartet2String(nodeptr p) {
	nodeptr q = p->back;
	int pNr = p->number;
	int qNr = q->number;
	int pNei1Nr = p->next->back->number;
	int pNei2Nr = p->next->next->back->number;
	int qNei1Nr = q->next->back->number;
	int qNei2Nr = q->next->next->back->number;
	stringstream middle;
	stringstream left;
	stringstream right;
	stringstream res;
	if (pNr < qNr) {
		middle << "-" << pNr << "-" << qNr << "-";
	} else {
		middle << "-" << qNr << "-" << pNr << "-";
	}
	if (pNei1Nr < pNei2Nr) {
		left << pNei1Nr << "-" << pNei2Nr;
	} else {
		left << pNei2Nr << "-" << pNei1Nr;
	}
	if (qNei1Nr < qNei2Nr) {
		right << qNei1Nr << "-" << qNei2Nr;
	} else {
		right << qNei2Nr << "-" << qNei1Nr;
	}
	res << left.str() << middle.str() << right.str();
	return res.str();
}

int evalNNIForBran(pllInstance* tr, partitionList *pr, nodeptr p, SearchInfo &searchinfo) {
	nodeptr q = p->back;
	assert(!isTip(p->number, tr->mxtips));
	assert(!isTip(q->number, tr->mxtips));
	bool recomputePartialLH = false;
	int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;
	int numPosNNI = 0;
	int i;
	pllNNIMove nni0; // dummy NNI to store backup information
	nni0.p = p;
	nni0.nniType = 0;
	nni0.likelihood = searchinfo.curLogl;
	for (i = 0; i < numBranches; i++) {
		nni0.z0[i] = p->z[i];
		nni0.z1[i] = p->next->z[i];
		nni0.z2[i] = p->next->next->z[i];
		nni0.z3[i] = q->next->z[i];
		nni0.z4[i] = q->next->next->z[i];
	}

	doOneNNI(tr, pr, p, 0, TOPO_ONLY);
	string quartetString = convertQuartet2String(p);
	doOneNNI(tr, pr, p, 0, TOPO_ONLY);
	if ( !searchinfo.tabunni || searchinfo.tabuNNIs.find(quartetString) == searchinfo.tabuNNIs.end() ) {
		/* do an NNI move of type 1 */
		double lh1 = doOneNNI(tr, pr, p, 0, searchinfo.evalType);
		pllNNIMove nni1;
		nni1.p = p;
		nni1.nniType = 0;
		// Store the optimized branch lengths
		for (i = 0; i < numBranches; i++) {
			nni1.z0[i] = p->z[i];
			nni1.z1[i] = p->next->z[i];
			nni1.z2[i] = p->next->next->z[i];
			nni1.z3[i] = q->next->z[i];
			nni1.z4[i] = q->next->next->z[i];
		}
		nni1.likelihood = lh1;
		nni1.loglDelta = lh1 - nni0.likelihood;
		nni1.negLoglDelta = -nni1.loglDelta;

		string quartetString = convertQuartet2String(p);
		nni1.quartetString = quartetString;
		searchinfo.nniList.push_back(nni1);

		if (nni1.likelihood > searchinfo.curLogl + 1e-6) {
			numPosNNI++;
			searchinfo.posNNIList.push_back(nni1);
		}

		/* Restore previous NNI move */
		doOneNNI(tr, pr, p, 0, TOPO_ONLY);
		/* Restore the old branch length */
		for (i = 0; i < numBranches; i++) {
			p->z[i] = nni0.z0[i];
			q->z[i] = nni0.z0[i];
			p->next->z[i] = nni0.z1[i];
			p->next->back->z[i] = nni0.z1[i];
			p->next->next->z[i] = nni0.z2[i];
			p->next->next->back->z[i] = nni0.z2[i];
			q->next->z[i] = nni0.z3[i];
			q->next->back->z[i] = nni0.z3[i];
			q->next->next->z[i] = nni0.z4[i];
			q->next->next->back->z[i] = nni0.z4[i];
		}
		recomputePartialLH = true;
	} else {
		searchinfo.numUnevalQuartet++;
	}

	doOneNNI(tr, pr, p, 1, TOPO_ONLY);
	quartetString = convertQuartet2String(p);
	doOneNNI(tr, pr, p, 1, TOPO_ONLY);
	if ( !searchinfo.tabunni || searchinfo.tabuNNIs.find(quartetString) == searchinfo.tabuNNIs.end()) {
		/* do an NNI move of type 2 */
		double lh2 = doOneNNI(tr, pr, p, 1, searchinfo.evalType);
		// Create the nniMove struct to store this move
		pllNNIMove nni2;
		nni2.p = p;
		nni2.nniType = 1;
		// Store the optimized and unoptimized central branch length
		for (i = 0; i < numBranches; i++) {
			nni2.z0[i] = p->z[i];
			nni2.z1[i] = p->next->z[i];
			nni2.z2[i] = p->next->next->z[i];
			nni2.z3[i] = q->next->z[i];
			nni2.z4[i] = q->next->next->z[i];
		}
		nni2.likelihood = lh2;
		nni2.loglDelta = lh2 - nni0.likelihood;
		nni2.negLoglDelta = -nni2.loglDelta;

		quartetString = convertQuartet2String(p);
		nni2.quartetString = quartetString;
		searchinfo.nniList.push_back(nni2);

		if (nni2.likelihood > searchinfo.curLogl + 1e-6) {
			numPosNNI++;
			searchinfo.posNNIList.push_back(nni2);
		}

		/* Restore previous NNI move */
		doOneNNI(tr, pr, p, 1, TOPO_ONLY);
	    /* Restore the old branch length */
		for (i = 0; i < numBranches; i++) {
			p->z[i] = nni0.z0[i];
			q->z[i] = nni0.z0[i];
			p->next->z[i] = nni0.z1[i];
			p->next->back->z[i] = nni0.z1[i];
			p->next->next->z[i] = nni0.z2[i];
			p->next->next->back->z[i] = nni0.z2[i];
			q->next->z[i] = nni0.z3[i];
			q->next->back->z[i] = nni0.z3[i];
			q->next->next->z[i] = nni0.z4[i];
			q->next->next->back->z[i] = nni0.z4[i];
		}
		recomputePartialLH = true;
	} else {
		searchinfo.numUnevalQuartet++;
	}

	if (recomputePartialLH) {
		if (numBranches > 1 && !tr->useRecom) {
			pllNewviewGeneric(tr, pr, p, PLL_TRUE);
			pllNewviewGeneric(tr, pr, p->back, PLL_TRUE);
		} else {
			pllNewviewGeneric(tr, pr, p, PLL_FALSE);
			pllNewviewGeneric(tr, pr, p->back, PLL_FALSE);
		}
	}

	return numPosNNI;
}

void pllSaveQuartet(nodeptr p, SearchInfo &searchinfo) {
	string quartetString = convertQuartet2String(p);
	searchinfo.tabuNNIs.insert(quartetString);
}

bool isAffectedBranch(nodeptr p, SearchInfo &searchinfo) {
	string branString = getBranString(p);
	if (searchinfo.affectBranches.find(branString) != searchinfo.affectBranches.end()) {
		return true;
	} else {
		return false;
	}
}

bool containsAffectedNodes(nodeptr p, SearchInfo &searchinfo) {
	bool res = false;
	nodeptr q = p->back;
	int id_list[6];
	id_list[0] = p->number;
	id_list[1] = q->number;
	id_list[2] = p->next->back->number;
	id_list[3] = p->next->next->back->number;
	id_list[4] = q->next->back->number;
	id_list[5] = q->next->next->back->number;
	for (int i = 0; i < 6; i++) {
		if (searchinfo.affectNodes.find(id_list[i]) != searchinfo.affectNodes.end()) {
			res = true;
			break;
		}
	}
	return res;
}

void evalNNIForSubtree(pllInstance* tr, partitionList *pr, nodeptr p, SearchInfo &searchinfo) {
	if (!isTip(p->number, tr->mxtips) && !isTip(p->back->number, tr->mxtips)) {
		evalNNIForBran(tr, pr, p, searchinfo);
//		int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;
		nodeptr q = p->next;
		while (q != p) {

//			if (numBranches > 1 && !tr->useRecom) {
//				pllNewviewGeneric(tr, pr, p->back, PLL_TRUE);
//			} else {
//				pllNewviewGeneric(tr, pr, p->back, PLL_FALSE);
//			}

			evalNNIForSubtree(tr, pr, q->back, searchinfo);

//			if (numBranches > 1 && !tr->useRecom) {
//				pllNewviewGeneric(tr, pr, p, PLL_TRUE);
//			} else {
//				pllNewviewGeneric(tr, pr, p, PLL_FALSE);
//			}

			q = q->next;
		}
	}
}



/**
 * DTH:
 * Announce the memory allocation error (for debugging)
 */
void pllAlertMemoryError(){
	printf("Memory error!!!!!!\n");;
	exit(1);
}

/**
* DTH:
* The PLL version of saveCurrentTree function
* @param tr: the tree (a pointer to a pllInstance)
* @param pr: pointer to a partitionList (this one keeps tons of tree info)
* @param p: root?
*/
void pllSaveCurrentTree(pllInstance* tr, partitionList *pr, nodeptr p){
	srand(gettime());
	double cur_logl = tr->likelihood;

	struct pllHashItem * item_ptr = (struct pllHashItem *) malloc(sizeof(struct pllHashItem));
	item_ptr->data = (int *) malloc(sizeof(int));
	item_ptr->next = NULL;
	item_ptr->str = NULL;

	unsigned int tree_index = -1;
	char * tree_str = NULL;
	pllTree2StringREC(tr->tree_string, tr, pr, tr->start->back, PLL_FALSE,
			PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_TRUE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
	tree_str = (char *) malloc (strlen(tr->tree_string) + 1);
	strcpy(tree_str, tr->tree_string);

	pll_boolean is_stored = PLL_FALSE;

	if(globalParam->store_candidate_trees){
		is_stored = pllHashSearch(pllUFBootDataPtr->treels, tree_str, &(item_ptr->data));
	}

	if(is_stored){ // if found the tree_str
		pllUFBootDataPtr->duplication_counter++;
		tree_index = *((int *)item_ptr->data);
		if (cur_logl <= pllUFBootDataPtr->treels_logl[tree_index] + 1e-4) {
			if (cur_logl < pllUFBootDataPtr->treels_logl[tree_index] - 5.0)
				if (verbose_mode >= VB_MED)
					printf("Current lh %f is much worse than expected %f\n",
							cur_logl, pllUFBootDataPtr->treels_logl[tree_index]);
/*			free(tree_str);
			free(item_ptr->data);
			free(item_ptr);*/
			return;
		}
		if (verbose_mode >= VB_MAX)
			printf("Updated logl %f to %f\n", pllUFBootDataPtr->treels_logl[tree_index], cur_logl);
		pllUFBootDataPtr->treels_logl[tree_index] = cur_logl;

		if (pllUFBootDataPtr->save_all_br_lens) {
			pllTree2StringREC(tr->tree_string, tr, pr, tr->start->back, PLL_TRUE,
        			PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_TRUE, PLL_SUMMARIZE_LENGTH, PLL_FALSE, PLL_FALSE);
        	char * tree_str_br_lens = (char *) malloc (strlen(tr->tree_string) + 1);
        	strcpy(tree_str_br_lens, tr->tree_string);
        	pllUFBootDataPtr->treels_newick[tree_index] = tree_str_br_lens;
		}
		if (pllUFBootDataPtr->boot_samples == NULL) {
			(pllUFBootDataPtr->treels_ptnlh)[tree_index] =
					(double *) malloc(pllUFBootDataPtr->n_patterns * sizeof(double));
			pllComputePatternLikelihood(tr, (pllUFBootDataPtr->treels_ptnlh)[tree_index], &cur_logl);
/*			free(tree_str);
			free(item_ptr->data);
			free(item_ptr);*/
			return;
		}
		if (verbose_mode >= VB_MAX)
			printf("Update treels_logl[%d] := %f\n", tree_index, cur_logl);

	} else {
       if (pllUFBootDataPtr->logl_cutoff != 0.0 && cur_logl <= pllUFBootDataPtr->logl_cutoff + 1e-4){
/*    		free(tree_str);
   		free(item_ptr->data);
   		free(item_ptr);*/
       	return;
       }

       if(pllUFBootDataPtr->treels_size == pllUFBootDataPtr->candidate_trees_count)
			pllResizeUFBootData();

       tree_index = pllUFBootDataPtr->candidate_trees_count;
       pllUFBootDataPtr->candidate_trees_count++;
       if (globalParam->store_candidate_trees){
           *((int *)item_ptr->data) = tree_index;
           item_ptr->str = tree_str;
       	pllHashAdd(pllUFBootDataPtr->treels, tree_str, item_ptr->data);
       }
       pllUFBootDataPtr->treels_logl[tree_index] = cur_logl;

       if (verbose_mode >= VB_MAX)
       	printf("Add    treels_logl[%d] := %f\n", tree_index, cur_logl);
   }

	//if (write_intermediate_trees)
	//        printTree(out_treels, WT_NEWLINE | WT_BR_LEN);

	double *pattern_lh = (double *) malloc(pllUFBootDataPtr->n_patterns * sizeof(double));
	if(!pattern_lh) pllAlertMemoryError();
	pllComputePatternLikelihood(tr, pattern_lh, &cur_logl);

	if (pllUFBootDataPtr->boot_samples == NULL) {
		// for runGuidedBootstrap
		pllUFBootDataPtr->treels_ptnlh[tree_index] = pattern_lh;
	} else {
		// online bootstrap
		int nptn = pllUFBootDataPtr->n_patterns;
		int updated = 0;
		int nsamples = globalParam->gbo_replicates;
		for (int sample = 0; sample < nsamples; sample++) {
			double rell = 0.0;
			for (int ptn = 0; ptn < nptn; ptn++)
				rell += pattern_lh[ptn] * pllUFBootDataPtr->boot_samples[sample][ptn];

//			int rand_pos = (sample + rand()) % nsamples;

			if (rell > pllUFBootDataPtr->boot_logl[sample] + globalParam->ufboot_epsilon ||
				(rell > pllUFBootDataPtr->boot_logl[sample] - globalParam->ufboot_epsilon &&
					random_double() <= 1.0/(pllUFBootDataPtr->boot_counts[sample]+1))) {
				if (!globalParam->store_candidate_trees){
					is_stored = pllHashSearch(pllUFBootDataPtr->treels, tree_str, &(item_ptr->data));
					if(is_stored)
						tree_index = *((int *)item_ptr->data);
					else{
						*((int *)item_ptr->data) = tree_index = pllUFBootDataPtr->candidate_trees_count - 1;
						item_ptr->str = tree_str;
						pllHashAdd(pllUFBootDataPtr->treels, tree_str, item_ptr->data);
					}
				}
				if (rell <= pllUFBootDataPtr->boot_logl[sample] +
						globalParam->ufboot_epsilon) {
					pllUFBootDataPtr->boot_counts[sample]++;
				} else {
					pllUFBootDataPtr->boot_counts[sample] = 1;
				}
				if(rell > pllUFBootDataPtr->boot_logl[sample])
					pllUFBootDataPtr->boot_logl[sample] = rell;
				pllUFBootDataPtr->boot_trees[sample] = tree_index;
				updated++;
			}
		}
/*		if (updated && verbose_mode >= VB_MAX)
		 printf("%d boot trees updated\n", updated);*/
	}
	if (pllUFBootDataPtr->save_all_br_lens) {
		pllTree2StringREC(tr->tree_string, tr, pr, tr->start->back, PLL_TRUE,
				PLL_FALSE, PLL_FALSE, PLL_FALSE, PLL_TRUE, PLL_SUMMARIZE_LH, PLL_FALSE, PLL_FALSE);
		char * s = (char *) malloc (strlen(tr->tree_string) + 1);
		strcpy(s, tr->tree_string);
		pllUFBootDataPtr->treels_newick[tree_index] = s;
	}

//	if(!globalParam->store_candidate_trees){
//		free(tree_str);
//		free(item_ptr->data);
//		free(item_ptr);
//	}
	if (pllUFBootDataPtr->boot_samples){
		free(pattern_lh);
		pllUFBootDataPtr->treels_ptnlh[tree_index] = NULL;
	}

//	printf("Done freeing: max = %d, count = %d, size = %d\n",
//			pllUFBootDataPtr->max_candidate_trees,
//			pllUFBootDataPtr->candidate_trees_count,
//			pllUFBootDataPtr->treels_size);
}

/**
* DTH:
* Extract the array of site log likelihood to be kept in ptnlh
* And update *cur_log
* @param tr: the tree (pointer to an pllInstance)
* @param ptnlh: to-be-kept array of site log likelihood
* @param cur_logl: pointer to current tree log likelihood
*/
void pllComputePatternLikelihood(pllInstance* tr, double * ptnlh, double * cur_logl){
	int i;
	double tree_logl = 0;
	for(i = 0; i < pllUFBootDataPtr->n_patterns; i++){
		ptnlh[i] = tr->lhs[i];
		tree_logl += tr->lhs[i] * tr->aliaswgt[i];
	}
	*cur_logl = tree_logl;
}

/**
* DTH:
* Resize some of the arrays in UFBootData if they're full
* Along with update treels_size (to track the size of these arrays)
*/
void pllResizeUFBootData(){
	int count = pllUFBootDataPtr->candidate_trees_count;
	pllUFBootDataPtr->treels_size = 2 * count;

	double * rtreels_logl =
			(double *) malloc(2 * count * (sizeof(double)));
	if(!rtreels_logl) pllAlertMemoryError();
//	memset(rtreels_logl, 0, 2 * count * sizeof(double));
	memcpy(rtreels_logl, pllUFBootDataPtr->treels_logl, count * sizeof(double));
	free(pllUFBootDataPtr->treels_logl);
	pllUFBootDataPtr->treels_logl = rtreels_logl;

	char ** rtreels_newick =
			(char **) malloc(2 * count * (sizeof(char *)));
	if(!rtreels_newick) pllAlertMemoryError();
	memset(rtreels_newick, 0, 2 * count * sizeof(char *));
	memcpy(rtreels_newick, pllUFBootDataPtr->treels_newick, count * sizeof(char *));
	free(pllUFBootDataPtr->treels_newick);
	pllUFBootDataPtr->treels_newick = rtreels_newick;

	double ** rtreels_ptnlh =
		(double **) malloc(2 * count * (sizeof(double *)));
	if(!rtreels_ptnlh) pllAlertMemoryError();
	memset(rtreels_ptnlh, 0, 2 * count * sizeof(double *));
	memcpy(rtreels_ptnlh, pllUFBootDataPtr->treels_ptnlh, count * sizeof(double *));
	free(pllUFBootDataPtr->treels_ptnlh);
	pllUFBootDataPtr->treels_ptnlh = rtreels_ptnlh;
}


/**
* DTH:
* (Based on function Tree2StringREC of PLL)
* Print out the tree topology with IQTree taxa ID (starts at 0) instead of PLL taxa ID (starts at 1)
* @param All are the same as in PLL's
*/
static char *pllTree2StringREC(char *treestr, pllInstance *tr, partitionList *pr, nodeptr p, pll_boolean printBranchLengths, pll_boolean printNames,
			    pll_boolean printLikelihood, pll_boolean rellTree, pll_boolean finalPrint, int perGene, pll_boolean branchLabelSupport, pll_boolean printSHSupport)
{
	char * result = treestr; // DTH: added this var to be able to remove the '\n' at the end
 char  *nameptr;

 if(isTip(p->number, tr->mxtips))
   {
     if(printNames)
	{
	  nameptr = tr->nameList[p->number];
	  sprintf(treestr, "%s", nameptr);
	}
     else
	sprintf(treestr, "%d", p->number - 1);

     while (*treestr) treestr++;
   }
 else
   {
     *treestr++ = '(';
     treestr = pllTree2StringREC(treestr, tr, pr, p->next->back, printBranchLengths, printNames, printLikelihood, rellTree,
			       finalPrint, perGene, branchLabelSupport, printSHSupport);
     *treestr++ = ',';
     treestr = pllTree2StringREC(treestr, tr, pr, p->next->next->back, printBranchLengths, printNames, printLikelihood, rellTree,
			       finalPrint, perGene, branchLabelSupport, printSHSupport);
     if(p == tr->start->back)
	{
	  *treestr++ = ',';
	  treestr = pllTree2StringREC(treestr, tr, pr, p->back, printBranchLengths, printNames, printLikelihood, rellTree,
				   finalPrint, perGene, branchLabelSupport, printSHSupport);
	}
     *treestr++ = ')';
   }

 if(p == tr->start->back)
   {
     if(printBranchLengths && !rellTree)
	sprintf(treestr, ":0.0;\n");
     else
	sprintf(treestr, ";\n");
   }
 else
   {
     if(rellTree || branchLabelSupport || printSHSupport)
	{
	  if(( !isTip(p->number, tr->mxtips)) &&
	     ( !isTip(p->back->number, tr->mxtips)))
	    {
	      assert(p->bInf != (branchInfo *)NULL);

	      if(rellTree)
		sprintf(treestr, "%d:%8.20f", p->bInf->support, p->z[0]);
	      if(branchLabelSupport)
		sprintf(treestr, ":%8.20f[%d]", p->z[0], p->bInf->support);
	      if(printSHSupport)
		sprintf(treestr, ":%8.20f[%d]", getBranchLength(tr, pr, perGene, p), p->bInf->support);

	    }
	  else
	    {
	      if(rellTree || branchLabelSupport)
		sprintf(treestr, ":%8.20f", p->z[0]);
	      if(printSHSupport)
		sprintf(treestr, ":%8.20f", getBranchLength(tr, pr, perGene, p));
	    }
	}
     else
	{
	  if(printBranchLengths)
	    sprintf(treestr, ":%8.20f", getBranchLength(tr, pr, perGene, p));
	  else
	    sprintf(treestr, "%s", "\0");
	}
   }

 if(result[strlen(result) - 1] == '\n') result[strlen(result) - 1] = '\0';
 while (*treestr) treestr++;
 return  treestr;
}

