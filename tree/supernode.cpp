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
#include "supernode.h"

SuperNode::SuperNode()
 : PhyloNode()
{
	init();
}


SuperNode::SuperNode(int aid) : PhyloNode(aid)
{
	init();
}

SuperNode::SuperNode(int aid, int aname) : PhyloNode (aid, aname) {
	init();
}


SuperNode::SuperNode(int aid, const char *aname) : PhyloNode(aid, aname) {
	init();
}

void SuperNode::init() {
	//partial_lh = nullptr;
}


void SuperNode::addNeighbor(Node *node, double length, int id) {
	neighbors.push_back(new SuperNeighbor(node, length, id));
}

SuperNode::~SuperNode()
{
}


