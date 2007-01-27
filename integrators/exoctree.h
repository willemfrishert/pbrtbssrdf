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

		for (int i = 0; i < 8; ++i)
			children[i] = NULL;
	}
	
	~ExOctNode() 
	{
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
	 * @description stores data if it's a leaf node
	 */
	vector<ExNodeData> data;
};
template <class ExNodeData, class LookupProc> class ExOctree {
public:
	// ExOctree Public Methods
	ExOctree(const BBox &b, int md = 5)
		: bound(b) {
			maxDepth = md;
	}
	
	void Add(const ExNodeData &dataItem, const BBox &dataBound, const LookupProc &process) 
	{
		addPrivate(&root, bound, dataItem, dataBound,
			DistanceSquared(dataBound.pMin, dataBound.pMax), process);
	}
	
	void Lookup(const Point &p, LookupProc &process) 
	{
		if (!bound.Inside(p)) return;
		lookupPrivate(&root, bound, p, process);
	}

	void Print(ExOctNode<ExNodeData> *node, FILE *f, int child = 0, int level = 0);

	void Print(FILE *f);

	// ExOctree Private Methods
private:
	void addPrivate(ExOctNode<ExNodeData> *node, const BBox &nodeBound,
		const ExNodeData &dataItem, const BBox &dataBound, float diag2,
		const LookupProc &process, int depth = 0);

	void lookupPrivate(ExOctNode<ExNodeData> *node, const BBox &nodeBound, const Point &P,
		LookupProc &process);

	void lookupPrivate(ExOctNode<ExNodeData> *node, const Point &P,
		LookupProc &process);

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
	float dataBBSquaredDiagonal, const LookupProc &process, int depth) 
{
	// Possibly add data item to current octree node
	if (depth == maxDepth ||
		DistanceSquared(nodeBound.pMin,
		nodeBound.pMax) < dataBBSquaredDiagonal) 
	{
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

	// store the information that is going to be passed to the child nodes
	addChildInformation(node, dataItem, nodeBound, dataBound, process);

	for (int child = 0; child < 8; ++child) 
	{
		// if it doesn't overlap the child, skip the node
		if (!over[child]) continue;

		if ( ! node->children[child] )
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

		addPrivate(node->children[child], childBound,
			dataItem, dataBound, dataBBSquaredDiagonal, process, depth+1);
	}
}

template <class ExNodeData, class LookupProc>
void ExOctree<ExNodeData, LookupProc>::lookupPrivate( 
	ExOctNode<ExNodeData> *node, const Point &p, LookupProc &process) 
{
	//printf("*");
	for (int i = 0; i < 8; i++)
	{
		ExOctNode<ExNodeData>* child = node->children[ i ];
		if( child != NULL )
		{
			// if SO, continue recursion
			if( (child->childLeaves != 0) && 
				process.subdivide(p, child->data, child->childLeaves) )
			{
				lookupPrivate(child, p, process);
			}
			else // just evaluate it
			{
				process.evaluate(p, child->data, child->childLeaves);
			}
			//printf("\n");
		}
	}
}

template <class ExNodeData, class LookupProc>
void ExOctree<ExNodeData, LookupProc>::lookupPrivate( 
	ExOctNode<ExNodeData> *node, const BBox &nodeBound,
	const Point &p, LookupProc &process) 
{
	/*** ******** MAIN BRANCH ********* ***/
	if ( node->childLeaves == 0 )
	{
		process.evaluate(p, node->data, 0);
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

	/*** ******** SECONDARY BRANCHES: for hierarchical evaluation ******** ***/
	for (int i = 0; i < 8; ++i)
	{
		ExOctNode<ExNodeData>* childNode = node->children[ i ];

		// if not the previously followed branch AND it exists
		if( (i != child) && (childNode != NULL) )
		{
			// if leaf: evaluate and step to next branch
			if( childNode->childLeaves == 0 )
			{
				process.evaluate(p, childNode->data, childNode->childLeaves);
				continue;
			} 

			// continue subdivision if process allows that 
			if( process.subdivide(p, childNode->data, childNode->childLeaves) )
			{
				lookupPrivate(childNode, p, process);
				//printf("\n");
			}
			else
			{
				process.evaluate(p, childNode->data, childNode->childLeaves);
			}
		}
	}
}

template <class ExNodeData, class LookupProc>
void ExOctree<ExNodeData, LookupProc>::Print(ExOctNode<ExNodeData> *node, FILE *f, int child, int level)
{
	
	if ( node )
	{
		if (node->childLeaves == 0)
		{
			fprintf(f, "<L>\n", child);
		}
		else
		{
			fprintf(f, "N%d( %d )\n", child, node->childLeaves);
		}
		vector<ExNodeData>::iterator dataIt = node->data.begin();
		for (; dataIt != node->data.end() ;dataIt++)
		{
			// indent data
			for (int i = 0; i < level; i++)
			{
				fprintf(f, "   ");
			}

			dataIt->Print( f );
			fprintf(f, "\n");
		}

		// kill recursion
		if (node->childLeaves == 0)
		{
			return;
		}
		
		level++;

		// start printing the children
		//fprintf(f, "\n");

		// iterate children
		for (int i = 0; i < 8; ++i)
		{
			if ( node->children[ i ] )
			{
				// indent child
				for (int l = 0; l < level; l++)
				{
					fprintf(f, "%d--", level);
				}

				Print(node->children[ i ], f, i, level);
			}
		}
	}
}

template <class ExNodeData, class LookupProc>
void ExOctree<ExNodeData, LookupProc>::Print(FILE *f)
{
	Print(&this->root, f);
}

#endif // PBRT_EXOCTREE_H
