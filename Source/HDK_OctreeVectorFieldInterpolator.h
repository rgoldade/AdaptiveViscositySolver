#ifndef HDK_OCTREEVECTORFIELDINTERPOLATOR_H
#define HDK_OCTREEVECTORFIELDINTERPOLATOR_H

#include <SIM/SIM_RawField.h>
#include <SIM/SIM_RawIndexField.h>
#include <SIM/SIM_ScalarField.h>

#include <UT/UT_Array.h>
#include <UT/UT_ParallelUtil.h>
#include <UT/UT_ThreadedAlgorithm.h>
#include <UT/UT_Vector3.h>
#include <UT/UT_VoxelArray.h>

#include "HDK_Utilities.h"
#include "HDK_OctreeGrid.h"

#include <iomanip>

class HDK_OctreeVectorFieldInterpolator
{
    enum NodeLabel
    {
        INACTIVENODE,
        ACTIVENODE,
        DEPENDENTNODE
    };

public:

    HDK_OctreeVectorFieldInterpolator(const HDK_OctreeGrid &octreeLabels,
					const UT_Array<UT_Array<SIM_RawField>> &velocity,
					const UT_Array<UT_Array<SIM_RawIndexField>> &velocityIndices)
	: myOctreeLabels(octreeLabels),
	myVelocity(velocity),
	myVelocityIndices(velocityIndices)
    {
        const int octreeLevels = myOctreeLabels.getOctreeLevels();

	assert(octreeLevels == velocity.entries() && octreeLevels == velocityIndices.entries());

	for (int level = 0; level < octreeLevels; ++level)
	{
	    assert(velocity[level].entries() == 3 && velocityIndices[level].entries() == 3);
	    for (int axis : {0,1,2})
		velocityIndices[level][axis].isAligned(&velocity[level][axis]);
	}

	myNodeValues.setSize(octreeLevels);
	myNodeLabels.setSize(octreeLevels);

        UT_Vector3 size = myOctreeLabels.getSize(), orig = myOctreeLabels.getOrig();

        UT_Array<bool> isTileOccupiedList;

        UT_Array<UT_Array<SIM_RawField>> nodeWeights;
	nodeWeights.setSize(octreeLevels);

        UT_Array<SIM_RawIndexField> nodeFlags;
        nodeFlags.setSize(octreeLevels);

        for (int level = 0; level < octreeLevels; ++level)
        {
	    UT_Vector3 voxelRes = myOctreeLabels.getVoxelRes(level);

            myNodeLabels[level].init(HDK_NODE, orig, size, voxelRes[0], voxelRes[1], voxelRes[2]);
            myNodeLabels[level].makeConstant(INACTIVENODE);

            myNodeValues[level].setSize(3);
            nodeWeights[level].setSize(3);

	    for (int axis : {0,1,2})
	    {
		myNodeValues[level][axis].init(HDK_NODE, orig, size, voxelRes[0], voxelRes[1], voxelRes[2]);
		myNodeValues[level][axis].makeConstant(0.);

    		nodeWeights[level][axis].init(HDK_NODE, orig, size, voxelRes[0], voxelRes[1], voxelRes[2]);
		nodeWeights[level][axis].makeConstant(0.);
            }

            nodeFlags[level].init(HDK_NODE, orig, size, voxelRes[0], voxelRes[1], voxelRes[2]);
            nodeFlags[level].makeConstant(0);

            isTileOccupiedList.clear();
            isTileOccupiedList.setSize(myNodeLabels[level].field()->numTiles());

#if !defined(NDEBUG)
	    const exint tileCount = myNodeLabels[level].field()->numTiles();

	    assert(nodeFlags[level].field()->numTiles() == tileCount);

	    for (int axis : {0,1,2})
		assert(myNodeValues[level][axis].field()->numTiles() == tileCount &&
			nodeWeights[level][axis].field()->numTiles() == tileCount);
#endif

            isTileOccupiedList.constant(false);

            findOccupiedNodeTiles(isTileOccupiedList, myOctreeLabels.getGridLabels(level), level);

            uncompressTiles(myNodeLabels[level], isTileOccupiedList);

	    for (int axis : {0,1,2})
            {
		uncompressTiles(myNodeValues[level][axis], isTileOccupiedList);
		uncompressTiles(nodeWeights[level][axis], isTileOccupiedList);
            }

            uncompressIndexTiles(nodeFlags[level], isTileOccupiedList);
        }

        // Set node labels. For each node, check if there is an active cell adjacent
	// to it at the same level. If yes, set active. We'll deal with t-junction
	// and split-edge nodes after.
	for (int level = 0; level < octreeLevels; ++level)
            setActiveNodes(level);

        // Sample values from all active adjacent faces at each node at each level.
        for (int level = 0; level < octreeLevels; ++level)
            sampleActiveNodes(nodeWeights[level], nodeFlags[level], level);

        // Pass the node values upwards in level if there is a co-located coarse node.
        for (int level = 0; level < octreeLevels - 1; ++level)
            bubbleActiveNodeValues(nodeWeights, nodeFlags, level);

        // Some nodes will not have a full bit count for their flag. This
        // happens at T-junction nodes and split edges. For each unset bit, we
        // need to sample the region around the node to build values where the
        // missing faces should exist.
        for (int level = 0; level < octreeLevels - 1; ++level)
            finishIncompleteNodes(nodeWeights, nodeFlags, level);

        for (int level = 0; level < octreeLevels; ++level)
            normalizeActiveNodes(nodeWeights[level], nodeFlags, level);

        for (int level = octreeLevels - 2; level >= 0; --level)
            distributeNodeValuesDown(level);

    }

    fpreal interpSPGrid(const UT_Vector3 &pos, const int axis) const;

private:

    THREADED_METHOD3_CONST(HDK_OctreeVectorFieldInterpolator,
			    octreeLabels.shouldMultiThread(),
			    findOccupiedNodeTiles,
			    UT_Array<bool> &, isTileOccupiedList,
			    const SIM_RawField &, octreeLabels,
			    const int, level)

    void findOccupiedNodeTilesPartial(UT_Array<bool> &isTileOccupiedList,
					const SIM_RawField &octreeLabels,
					const int level,
					const UT_JobInfo &info) const;

    THREADED_METHOD2(HDK_OctreeVectorFieldInterpolator,
			isTileOccupiedList.entries() > 20,
			uncompressTiles,
			SIM_RawField &, grid,
			const UT_Array<bool> &, isTileOccupiedList)

    // The macro approach to multithreaded operations over grids makes the
    // following method hard to template. Therefore we have two that do
    // essentially the same thing.

    void uncompressTilesPartial(SIM_RawField &grid,
				const UT_Array<bool> &isTileOccupiedList,
				const UT_JobInfo &info);

    THREADED_METHOD2(HDK_OctreeVectorFieldInterpolator,
			isTileOccupiedList.entries() > 20,
			uncompressIndexTiles,
			SIM_RawIndexField &, grid,
			const UT_Array<bool> &, isTileOccupiedList)

    void uncompressIndexTilesPartial(SIM_RawIndexField &grid,
					const UT_Array<bool> &isTileOccupiedList,
					const UT_JobInfo &info);

    THREADED_METHOD1(HDK_OctreeVectorFieldInterpolator,
			myNodeLabels[level].shouldMultiThread(),
			setActiveNodes,
			const int, level)

    void setActiveNodesPartial(const int level, const UT_JobInfo &info);

    THREADED_METHOD3(HDK_OctreeVectorFieldInterpolator,
			myNodeLabels[level].shouldMultiThread(),
			sampleActiveNodes,
			UT_Array<SIM_RawField> &, nodeWeights,
			SIM_RawIndexField &, nodeFlags,
			const int, level)

    void sampleActiveNodesPartial(UT_Array<SIM_RawField> &nodeWeights,
				    SIM_RawIndexField &nodeFlags,
				    const int level,
				    const UT_JobInfo &info);

    THREADED_METHOD3(HDK_OctreeVectorFieldInterpolator,
			myNodeLabels[level].shouldMultiThread(),
			bubbleActiveNodeValues,
			UT_Array<UT_Array<SIM_RawField>> &, nodeWeights,
			UT_Array<SIM_RawIndexField> &, nodeFlags,
			const int, level)

    void bubbleActiveNodeValuesPartial(UT_Array<UT_Array<SIM_RawField>> &nodeWeights,
					UT_Array<SIM_RawIndexField> &nodeFlags,
					const int level,
					const UT_JobInfo &info);

    THREADED_METHOD3(HDK_OctreeVectorFieldInterpolator,
			myNodeLabels[level].shouldMultiThread(),
			finishIncompleteNodes,
			UT_Array<UT_Array<SIM_RawField>> &, nodeWeights,
			UT_Array<SIM_RawIndexField> &, nodeFlags,
			const int, level)

    void finishIncompleteNodesPartial(UT_Array<UT_Array<SIM_RawField>> &nodeWeights,
					UT_Array<SIM_RawIndexField> &nodeFlags,
					const int level,
					const UT_JobInfo &info);

    THREADED_METHOD3(HDK_OctreeVectorFieldInterpolator,
			myNodeLabels[level].shouldMultiThread(),
			normalizeActiveNodes,
			UT_Array<SIM_RawField> &, nodeWeights,
			UT_Array<SIM_RawIndexField> &, nodeFlags,
			const int, level)

    void normalizeActiveNodesPartial(UT_Array<SIM_RawField> &nodeWeights,
					UT_Array<SIM_RawIndexField> &nodeFlags,
					const int level,
					const UT_JobInfo &info);

    THREADED_METHOD1(HDK_OctreeVectorFieldInterpolator,
			myNodeLabels[level].shouldMultiThread(),
			distributeNodeValuesDown,
			const int, level)

    void distributeNodeValuesDownPartial(const int level, const UT_JobInfo &info);

    UT_Array<SIM_RawField> myNodeLabels;
    UT_Array<UT_Array<SIM_RawField>> myNodeValues;

    const HDK_OctreeGrid &myOctreeLabels;
    const UT_Array<UT_Array<SIM_RawField>> &myVelocity;
    const UT_Array<UT_Array<SIM_RawIndexField>> &myVelocityIndices;
};

#endif