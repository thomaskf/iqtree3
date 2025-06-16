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
#ifndef SUPERNODE_H
#define SUPERNODE_H

#include "phylonode.h"

typedef vector<PhyloNeighbor*> PhyloNeighborVec;

/**
A neighbor in a phylogenetic SUPER tree

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class SuperNeighbor : public PhyloNeighbor {

	friend class SuperNode;
	friend class PhyloSuperTree;

public:
	/**
		construct class with a node and length
		@param anode the other end of the branch
		@param alength length of branch
	*/
	SuperNeighbor(Node *anode, double alength) : PhyloNeighbor(anode, alength) {	
	}

	/**
		construct class with a node and length
		@param anode the other end of the branch
		@param alength length of branch
		@param aid branch ID
	*/
	SuperNeighbor(Node *anode, double alength, int aid) : PhyloNeighbor(anode, alength, aid) {	
	}

    /**
     construct class with another Neighbor
     @param nei another Neighbor
     */
    SuperNeighbor(SuperNeighbor *nei) : PhyloNeighbor(nei) {
    }
    
    
    /**
     allocate a new Neighbor by just copying from this one
     @return pointer to newly created Neighbor
     */
    virtual Neighbor* newNeighbor() override {
        return (new SuperNeighbor(this));
    }

	/**
		vector of size m (m = #partitions)
	*/
	PhyloNeighborVec link_neighbors;

};

/**
Node of a super tree

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/
class SuperNode : public PhyloNode
{
	friend class PhyloSuperTree;

public:
	/**
		constructor 
	*/
    SuperNode();

	/**
		constructor 
		@param aid id of this node
	*/
	SuperNode(int aid);

	/**
		constructor 
		@param aid id of this node
		@param aname name of this node
	*/
	SuperNode(int aid, int aname);

	/**
		constructor 
		@param aid id of this node
		@param aname name of this node
	*/
	SuperNode(int aid, const char *aname);

	/**
		initialization
	*/
	void init();

	/**
		add a neighbor
		@param node the neighbor node
		@param length branch length
		@param id branch ID
	*/
	virtual void addNeighbor(Node *node, double length, int id = -1) override;

    ~SuperNode() override;

};

#endif
