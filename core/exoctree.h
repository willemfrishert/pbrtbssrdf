#ifndef PBRT_EXOCTREE_H
#define PBRT_EXOCTREE_H
// octree.h*
#include "pbrt.h"
#include "geometry.h"

// ExOctree Declarations
template <class ExNodeData> struct ExOctNode 
{
	ExOctNode() 
	{
		childLeaves = 0;
		childData = NULL;

		for (int i = 0; i < 8; ++i)
			children[i] = NULL;
	}
	
	~ExOctNode() 
	{
		delete childData;

		for (int i = 0; i < 8; ++i)
			delete children[i];
	}
	
	/**
	 * @description number of child leaves
	 */
	u_int childLeaves;

	/**
	 * @description pointers to the 8 potential children
	 */
	ExOctNode *children[8];
	
	/**
	 * @description stores average data from its children
	 */
	ExNodeData* childData;
	
	/**
	 * @description stores data if it's a leaf node
	 */
	vector<ExNodeData> data;
};
template <class ExNodeData, class LookupProc> class ExOctree {
public:
	// ExOctree Public Methods
	ExOctree(const BBox &b, int md = 16)
		: bound(b) {
			maxDepth = md;
	}
	
	void Add(const ExNodeData &dataItem, const BBox &dataBound, const LookupProc &process) 
	{
		addPrivate(&root, bound, dataItem, dataBound,
			DistanceSquared(dataBound.pMin, dataBound.pMax), process);
	}
	
	void Lookup(const Point &p, const LookupProc &process) 
	{
		if (!bound.Inside(p)) return;
		lookupPrivate(&root, bound, p, process);
	}

	// ExOctree Private Methods
private:
	void addPrivate(ExOctNode<ExNodeData> *node, const BBox &nodeBound,
		const ExNodeData &dataItem, const BBox &dataBound, float diag2,
		const LookupProc &process, int depth = 0);

	void lookupPrivate(ExOctNode<ExNodeData> *node, const BBox &nodeBound, const Point &P,
		const LookupProc &process);

	void addChildInformation(ExOctNode<ExNodeData> *node, const ExNodeData& dataItem, 
		const BBox &nodeBound, const BBox &dataBound, const LookupProc &process);

	// ExOctree Private Data
private:
	int maxDepth;
	BBox bound;
	ExOctNode<ExNodeData> root;
};

/* ******* ExOctree Method Definitions ******* */

/**
 * @param node
 * @param dataItem
 * @param nodeBound the child node's bounding box
 * @param dataBound the data's bounding box
 * @param process
 */
template <class ExNodeData, class LookupProc>
void ExOctree<ExNodeData, LookupProc>::addChildInformation( ExOctNode<ExNodeData> *node, 
														   const ExNodeData& dataItem, const BBox &nodeBound, 
														   const BBox &dataBound, const LookupProc &process )
{
	node->childLeaves++;
	
	// the process will manage the dataItem and add it in the 
	// convenient way to the node
	process.addChildNode(node, dataItem, dataBound, nodeBound);
}

template <class ExNodeData, class LookupProc>
void ExOctree<ExNodeData, LookupProc>::addPrivate(ExOctNode<ExNodeData> *node, const BBox &nodeBound, 
	const ExNodeData &dataItem, const BBox &dataBound, 
	float diag2, const LookupProc &process, int depth) 
{
	// Possibly add data item to current octree node
	if (depth == maxDepth ||
		DistanceSquared(nodeBound.pMin,
		nodeBound.pMax) < diag2) 
	{
		/************************************************************************/
		/* HOW TO PROCEED IF IT IS A LEAF NODE??????? PROCESS SHOULD TAKE CARE OF IT
		/************************************************************************/

		// Adding data to the current child
		node->data.push_back(dataItem);
		return;
	}

	// Otherwise add data item to octree children
	Point pMid = .5 * nodeBound.pMin + .5 * nodeBound.pMax;

	// Determine which children the item overlaps/spans to
	bool over[8];
	over[0] = over[1] =
		over[2] =
		over[3] = (dataBound.pMin.x <= pMid.x);
	over[4] = over[5] =
		over[6] =
		over[7] = (dataBound.pMax.x  > pMid.x);
	over[0] &= (dataBound.pMin.y <= pMid.y);
	over[1] &= (dataBound.pMin.y <= pMid.y);
	over[4] &= (dataBound.pMin.y <= pMid.y);
	over[5] &= (dataBound.pMin.y <= pMid.y);
	over[2] &= (dataBound.pMax.y  > pMid.y);
	over[3] &= (dataBound.pMax.y  > pMid.y);
	over[6] &= (dataBound.pMax.y  > pMid.y);
	over[7] &= (dataBound.pMax.y  > pMid.y);
	over[0] &= (dataBound.pMin.z <= pMid.z);
	over[2] &= (dataBound.pMin.z <= pMid.z);
	over[4] &= (dataBound.pMin.z <= pMid.z);
	over[6] &= (dataBound.pMin.z <= pMid.z);
	over[1] &= (dataBound.pMax.z  > pMid.z);
	over[3] &= (dataBound.pMax.z  > pMid.z);
	over[5] &= (dataBound.pMax.z  > pMid.z);
	over[7] &= (dataBound.pMax.z  > pMid.z);

	for (int child = 0; child < 8; ++child) 
	{
		// if it doesn't overlap the child, skip the node
		if (!over[child]) continue;

		if (!node->children[child])
		{
			node->children[child] = new ExOctNode<ExNodeData>;
		}

		// Compute _childBound_ for octree child _child_
		BBox childBound;
		childBound.pMin.x = (child & 4) ? pMid.x : nodeBound.pMin.x;
		childBound.pMax.x = (child & 4) ? nodeBound.pMax.x : pMid.x;
		childBound.pMin.y = (child & 2) ? pMid.y : nodeBound.pMin.y;
		childBound.pMax.y = (child & 2) ? nodeBound.pMax.y : pMid.y;
		childBound.pMin.z = (child & 1) ? pMid.z : nodeBound.pMin.z;
		childBound.pMax.z = (child & 1) ? nodeBound.pMax.z : pMid.z;

		// otherwise store the information that is going to be passed to this child
		addChildInformation(node, dataItem, childBound, dataBound, process);

		addPrivate(node->children[child], childBound,
			dataItem, dataBound, diag2, process, depth+1);
	}
}
template <class ExNodeData, class LookupProc>
void ExOctree<ExNodeData, LookupProc>::lookupPrivate( 
	ExOctNode<ExNodeData> *node, const BBox &nodeBound,
	const Point &p, const LookupProc &process) 
{

	// if the lookup process decides to stop the search
	// stop it immediatly, otherwise just let the
	// octree continue the recursion
	if ( process(p, node->data, node->childData) )
	{
		return;
	}

	// Determine which octree child node _p_ is inside
	// Midpoint in a cube
	Point pMid = .5f * nodeBound.pMin + .5f * nodeBound.pMax;

	int child = (p.x > pMid.x ? 4 : 0) +
		(p.y > pMid.y ? 2 : 0) + (p.z > pMid.z ? 1 : 0);

	if (node->children[child]) 
	{
		// Compute _childBound_ for octree child _child_
		BBox childBound;
		childBound.pMin.x = (child & 4) ? pMid.x : nodeBound.pMin.x;
		childBound.pMax.x = (child & 4) ? nodeBound.pMax.x : pMid.x;
		childBound.pMin.y = (child & 2) ? pMid.y : nodeBound.pMin.y;
		childBound.pMax.y = (child & 2) ? nodeBound.pMax.y : pMid.y;
		childBound.pMin.z = (child & 1) ? pMid.z : nodeBound.pMin.z;
		childBound.pMax.z = (child & 1) ? nodeBound.pMax.z : pMid.z;
		lookupPrivate(node->children[child], childBound, p,
			process);
	}
}
#endif // PBRT_EXOCTREE_H
