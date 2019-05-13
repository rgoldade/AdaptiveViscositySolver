#ifndef HDK_OCTREEGRID_H
#define HDK_OCTREEGRID_H

#include <SIM/SIM_RawField.h>
#include <SIM/SIM_ScalarField.h>
#include <SIM/SIM_VectorField.h>

#include <UT/UT_ParallelUtil.h>
#include <UT/UT_PerfMonAutoEvent.h>
#include <UT/UT_ThreadedAlgorithm.h>
#include <UT/UT_VoxelArray.h>

#include "HDK_Utilities.h"

// Houdini uses a prefix instead of namespaces, so we'll use HDK to prevent
// and collisions with their internal tools.

class HDK_OctreeGrid
{
public:
    // Every cell in the octree will be either:
    // 1. INACTIVE --> a grid cell at the finest level that has no ACTIVE ancestors.
    //			    Only the finest level should be labelled INACTIVE
    //			    and it corresponds to grid cells that are outside
    //			    out of the region of interest.
    // 2. ACTIVE --> a grid cell that is equivalent to a leaf node in an
    //			actual tree representation of the octree.
    // 3. UP --> a grid cell that is conceptually a descendant of a leaf node
    //		 in an actual tree. UP cells are necessary in this pyramid
    //		 structure to locate ACTIVE cells.
    // 4. DOWN --> a grid cell that is the ancestor of an ACTIVE cell. Conceptually
    //		    equivalent to an internal node in a tree structure.
    enum OctreeCellLabel
    {
        INACTIVE,
        ACTIVE,
        UP,
        DOWN
    };

    HDK_OctreeGrid() {}

    // Refinement mask should be aligned with the reference grid and store > 0
    // values for inactive voxels, 0 values for active voxels and < 0 values for
    // voxels inside the octree (e.g. UP voxels).
    void init(const SIM_RawField &initializerMask, const int desiredLevels)
    {
	UT_Vector3i voxelRes;
        initializerMask.getVoxelRes(voxelRes[0], voxelRes[1], voxelRes[2]);

	// We will stretch the octree grid to be the smallest power of 2
	// that will contain the entire supplied grid size.
	// Because of this, we need to keep the original grid size so
	// we don't set ACTIVE cells outside of the supplied grid region.
        myBaseVoxelRes[0] = voxelRes[0];
        myBaseVoxelRes[1] = voxelRes[1];
        myBaseVoxelRes[2] = voxelRes[2];

        for (int axis : {0,1,2})
        {
            fpreal logSize = std::log2(fpreal(voxelRes[axis]));
	    logSize = std::ceil(logSize);

	    voxelRes[axis] = std::exp2(logSize);
        }

        // Cap levels to be the floor of the log-base-2 of the smallest
        // axis-aligned dimension.
	// We do this because we don't want empty levels. For example, an
        // input grid with that is 64^3 can only represent 6 levels so even
	// if the input is 7 levels, we should only generate 6.

        myOctreeLevels = desiredLevels;
        if (voxelRes[0] > 0 && voxelRes[1] > 0 && voxelRes[2] > 0)
        {
            for (const int axis : {0,1,2})
            {
                if (std::log2(voxelRes[axis]) < myOctreeLevels)
                    myOctreeLevels = std::log2(voxelRes[axis]);
            }
        }

        myOctreeGridLabels.clear();
	myOctreeGridLabels.setSize(myOctreeLevels);

	UT_Vector3 size = initializerMask.getSize();
        UT_Vector3 origin = initializerMask.getOrig();

	UT_Vector3 newSize;

        if (myBaseVoxelRes[0] != voxelRes[0] || myBaseVoxelRes[1] != voxelRes[1] || myBaseVoxelRes[2] != voxelRes[2])
	{
	    assert(myBaseVoxelRes[0] <= voxelRes[0] && myBaseVoxelRes[1] <= voxelRes[1] && myBaseVoxelRes[2] <= voxelRes[2]);
            newSize = initializerMask.getVoxelSize() * UT_Vector3(voxelRes[0], voxelRes[1], voxelRes[2]);
	}
        else
            newSize = size;

	myOctreeGridLabels[0].init(SIM_SAMPLE_CENTER, origin, newSize, voxelRes[0], voxelRes[1], voxelRes[2]);
        myOctreeGridLabels[0].makeConstant(INACTIVE);


	for (int level = 1; level < myOctreeLevels; ++level)
        {
            voxelRes[0] /= 2;
            voxelRes[1] /= 2;
            voxelRes[2] /= 2;

            myOctreeGridLabels[level].init(SIM_SAMPLE_CENTER, origin, newSize, voxelRes[0], voxelRes[1], voxelRes[2]);
            myOctreeGridLabels[level].makeConstant(INACTIVE);
        }

        //
        // Build the octree
        //

	// Pass over the finest resolution and activate cells where the mask grid is 0. Then activate its siblings.
	setBaseGridLabels(initializerMask);

        int threadCount = UT_Thread::getNumProcessors();

        UT_Array<UT_Array<UT_Vector3i>> setParentDownList;
        setParentDownList.setSize(threadCount);

        UT_Array<UT_Array<UT_Vector3i>> setParentActivelist;
        setParentActivelist.setSize(threadCount);

        UT_Array<UT_Array<UT_Vector3i>> setParentUpList;
        setParentUpList.setSize(threadCount);

        UT_Array<UT_Vector3i> newParentList;
        UT_Array<bool> isTileOccupiedList;

        auto vectorCompare = [](const UT_Vector3I &a, const UT_Vector3I &b) -> bool
	{
            if (a[0] < b[0])
                return true;

            else if (a[0] == b[0])
            {
                if (a[1] < b[1])
                    return true;
                else if (a[1] == b[1] && a[2] < b[2])
                    return true;
            }

            return false;
        };

        for (int level = 0; level < myOctreeLevels - 1; ++level)
        {
            // First build pass:
            // If cell is labeled UP and its sibling is labeled ACTIVE, set it to ACTIVE.
            // If a cell is labelled ACTIVE, add parent to list of cells to label DOWN.

            isTileOccupiedList.clear();
	    isTileOccupiedList.setSize(myOctreeGridLabels[level + 1].field()->numTiles());
            isTileOccupiedList.constant(false);

            setActiveCellsAndParentList(setParentDownList, isTileOccupiedList, level);

            uncompressParentTiles(isTileOccupiedList, level);

            newParentList.clear();
            for (int thread = 0; thread < threadCount; ++thread)
            {
                newParentList.concat(setParentDownList[thread]);
                setParentDownList[thread].clear();
            }

            UTparallelSort(newParentList.begin(), newParentList.end(), vectorCompare);
            setParentCellLabel(newParentList, level, DOWN);

            // Second pass:
            // If cell is labeled ACTIVE, any face-adjacent cells that are labeled UP
            // will have their parents listed as ACTIVE. If cell is DOWN, set parent to DOWN.

            isTileOccupiedList.constant(false);
	    setFaceGrading(setParentDownList, setParentActivelist, isTileOccupiedList, level);

            uncompressParentTiles(isTileOccupiedList, level);

            // Set parents to DOWN
            newParentList.clear();
            for (int thread = 0; thread < threadCount; ++thread)
            {
                newParentList.concat(setParentDownList[thread]);
                setParentDownList[thread].clear();
            }

            UTparallelSort(newParentList.begin(), newParentList.end(), vectorCompare);
            setParentCellLabel(newParentList, level, DOWN);

            // Set parents to ACTIVE
            newParentList.clear();
	    for (int thread = 0; thread < threadCount; ++thread)
            {
                newParentList.concat(setParentActivelist[thread]);
                setParentActivelist[thread].clear();
            }

            UTparallelSort(newParentList.begin(), newParentList.end(), vectorCompare);
            setParentCellLabel(newParentList, level, ACTIVE);

            // Third pass:
            // If cell is UP and the parent is INACTIVE, set parent to UP.

            isTileOccupiedList.constant(false);
            setParentsUp(setParentUpList, isTileOccupiedList, level);

            uncompressParentTiles(isTileOccupiedList, level);

            newParentList.clear();
            for (int thread = 0; thread < threadCount; ++thread)
            {
                newParentList.concat(setParentUpList[thread]);
                setParentUpList[thread].clear();
            }

            UTparallelSort(newParentList.begin(), newParentList.end(), vectorCompare);
            setParentCellLabel(newParentList, level, UP);

            myOctreeGridLabels[level + 1].fieldNC()->collapseAllTiles();
        }

        // Clean up pass on the top level
        setTopLevel(myOctreeLevels - 1);

        // Traverse the grid levels and find the lowest level without any ACTIVE
        // voxels. Since there are no active grid cells above this level, we
        // might as well set the internal maximum level to this capped level.

        int cappedLevel = 0;
        for (; cappedLevel < myOctreeLevels; ++cappedLevel)
        {
            bool levelHasActiveCell = false;
            checkActiveCellAtLevel(levelHasActiveCell, cappedLevel);

            if (!levelHasActiveCell)
            {
                assert(cappedLevel != 0);
                break;
            }
        }

        myOctreeLevels = cappedLevel;

#if !defined(NDEBUG)

	// TODO: this could be made parallel
        UT_Interrupt *boss = UTgetInterrupt();
        for (int level = 0; level < myOctreeLevels; ++level)
        {
            UT_VoxelArrayIteratorF vit;
	    vit.setConstArray(myOctreeGridLabels[level].field());
            UT_VoxelTileIteratorF vitt;

            int activeCellCount = 0;

            for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
            {
                if (boss->opInterrupt())
                    break;

                if (!vit.isTileConstant() || vit.getValue() == ACTIVE)
                {
                    vitt.setTile(vit);

                    for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
                    {
                        if (vitt.getValue() == ACTIVE)
                            ++activeCellCount;
                    }
                }
            }
        }
#endif
    }

    // Make sure the tree has one and only one active cell in a vertical column.
    bool unitTest() const;

    void outputOctreeGeometry(GU_Detail &octreeGeometryDetail) const;

    SYS_FORCE_INLINE UT_Vector3I getParentCell(UT_Vector3i cell) const
    {
	// TODO: check if I can use /= with an integer
	cell /= UT_Vector3i(2);
        return cell;
    }

    // Computing parent faces is the same as parent cells.
    SYS_FORCE_INLINE UT_Vector3I getParentFace(UT_Vector3i face) const
    {
	return getParentCell(face);
    }

    SYS_FORCE_INLINE UT_Vector3I getParentNode(UT_Vector3i node) const
    {
	return getParentCell(node);
    }

    SYS_FORCE_INLINE UT_Vector3i getChildCell(UT_Vector3i cell, const int childIndex) const
    {
        assert(childIndex < 8);
	cell *= 2;
	cell = childCellOffset(cell, childIndex);
        return cell;
    }

    SYS_FORCE_INLINE UT_Vector3i childCellOffset(UT_Vector3i cell, const int childIndex) const
    {
	assert(childIndex < 8);

	// Child index is a binary flag indicating a forward offset of 1
	// in each axis direction.
	for (int axis : {0,1,2})
	{
	    if (childIndex & (1 << axis))
		++cell[axis];
	}

	return cell;
    }

    SYS_FORCE_INLINE UT_Vector3i getChildFace(UT_Vector3i face, const int axis, const int childIndex) const
    {
        assert(axis < 3 && childIndex < 4);

	face *= 2;

	if (childIndex & 1)
	    ++face[(axis + 1) % 3];
	if (childIndex & 2)
	    ++face[(axis + 2) % 3];

        return face;
    }

    SYS_FORCE_INLINE UT_Vector3i getChildEdge(UT_Vector3i edge, const int edgeAxis, const int childIndex) const
    {
        assert(edgeAxis < 3 && childIndex < 2);

	edge *= 2;
	if (childIndex > 0)
	    ++edge[edgeAxis];

        return edge;
    }

    SYS_FORCE_INLINE UT_Vector3i getChildNode(UT_Vector3i node) const
    {
	node *= 2;
        return node;
    }

    // Get edges that lay inset in a face from a level lower
    SYS_FORCE_INLINE UT_Vector3i getChildEdgeInFace(UT_Vector3i face, const int faceAxis, const int edgeAxis, const int childIndex) const
    {
        assert(edgeAxis != faceAxis &&
		faceAxis < 3 &&
		edgeAxis < 3 &&
		childIndex < 2);

	face *= 2;

	if (childIndex == 1)
	    ++face[edgeAxis];

	int offsetAxis = 3 - faceAxis - edgeAxis;
	++face[offsetAxis];

        return face;
    }

    SYS_FORCE_INLINE const SIM_RawField& getGridLabels(const int level) const
    {
        assert(level < myOctreeGridLabels.size());
        return myOctreeGridLabels[level];
    }

    SYS_FORCE_INLINE UT_Vector3i getVoxelRes(const int level) const
    {
	UT_Vector3i voxelCount;
        assert(level < myOctreeGridLabels.size());
        myOctreeGridLabels[level].getVoxelRes(voxelCount[0], voxelCount[1], voxelCount[2]);

	return voxelCount;
    }

    SYS_FORCE_INLINE UT_Vector3 getSize() const { return myOctreeGridLabels[0].getSize(); }

    SYS_FORCE_INLINE UT_Vector3 getOrig() const { return myOctreeGridLabels[0].getOrig(); }

    SYS_FORCE_INLINE UT_Vector3 getVoxelSize(int level) const
    {
        assert(level >= 0 && level < myOctreeLevels);
        return myOctreeGridLabels[level].getVoxelSize();
    }

    SYS_FORCE_INLINE UT_Vector3 indexToPos(const UT_Vector3i &indexPosition, const int level) const
    {
	UT_Vector3 worldPosition;
        myOctreeGridLabels[level].indexToPos(indexPosition[0],
						indexPosition[1],
						indexPosition[2],
						worldPosition);
	return worldPosition;
    }

    SYS_FORCE_INLINE UT_Vector3 posToIndex(const UT_Vector3 &worldPosition, const int level) const
    {
	UT_Vector3 indexPosition;
        myOctreeGridLabels[level].posToIndex(worldPosition, indexPosition);
	return indexPosition;
    }

    SYS_FORCE_INLINE int getOctreeLevels() const { return myOctreeLevels; }

    SYS_FORCE_INLINE OctreeCellLabel getCellLabel(const UT_Vector3i &cell, const int level) const
    {
        assert(level >= 0 && level < myOctreeLevels);
	return OctreeCellLabel(HDKgetFieldValue(myOctreeGridLabels[level], cell));
    }

    SYS_FORCE_INLINE bool isCellActive(const UT_Vector3i &cell, const int level) const
    {
        return getCellLabel(cell, level) == ACTIVE;
    }

    void refineGrid();

private:

    //
    // Unit test helpers
    //

    THREADED_METHOD1_CONST(SIM_OctreeGrid, myOctreeGridLabels[0].shouldMultiThread(),
                           activeCountUnitTest, bool &, passed)

    void activeCountUnitTestPartial(bool &passed, const UT_JobInfo &info) const;

    THREADED_METHOD2_CONST(SIM_OctreeGrid,
                           myOctreeGridLabels[level].shouldMultiThread(),
                           upAdjacentUnitTest, bool &, passed, const int, level)

    void upAdjacentUnitTestPartial(bool &passed, const int level,
                                   const UT_JobInfo &info) const;

    THREADED_METHOD2_CONST(SIM_OctreeGrid,
                           myOctreeGridLabels[level].shouldMultiThread(),
                           activeUnitTest, bool &, passed, const int, level)

    void activeUnitTestPartial(bool &passed, const int level,
                               const UT_JobInfo &info) const;

    UT_Array<UT_Vector4i> getFaceAdjacentCells(const UT_Vector3i &cell,
						const int axis,
						const int direction,
						const int level) const;

    //
    // Build octree
    //

    THREADED_METHOD1(SIM_OctreeGrid, myOctreeGridLabels[0].shouldMultiThread(),
                     setBaseGridLabels,
		     const SIM_RawField &, initializerMask)

    void setBaseGridLabelsPartial(const SIM_RawField &initializerMask, const UT_JobInfo &info);

    THREADED_METHOD3(SIM_OctreeGrid, myOctreeGridLabels[level].shouldMultiThread(),
			setActiveCellsAndParentList,
			UT_Array<UT_Array<UT_Vector3i>> &, parallelSetParentDownList,
			UT_Array<bool> &, isTileOccupiedList,
			const int, level)

    void setActiveCellsAndParentListPartial(UT_Array<UT_Array<UT_Vector3i>> &parallelSetParentDownList,
					    UT_Array<bool> &isTileOccupiedList,
					    const int level,
					    const UT_JobInfo &info);

    THREADED_METHOD2(SIM_OctreeGrid, isTileOccupiedList.entries() > 100,
			uncompressParentTiles,
			const UT_Array<bool> &, isTileOccupiedList,
			const int, level)

    void uncompressParentTilesPartial(const UT_Array<bool> &isTileOccupiedList,
					const int level,
					const UT_JobInfo &info);

    THREADED_METHOD3(SIM_OctreeGrid, newParentCellList.entries() > 100,
			setParentCellLabel,
			const UT_Array<UT_Vector3i> &, newParentCellList,
			const int, level,
			const OctreeCellLabel, label)

    void setParentCellLabelPartial(const UT_Array<UT_Vector3i> &newParentCellList,
				    const int level,
				    const OctreeCellLabel label,
				    const UT_JobInfo &info);

    THREADED_METHOD4(SIM_OctreeGrid, myOctreeGridLabels[level].shouldMultiThread(),
			setFaceGrading,
			UT_Array<UT_Array<UT_Vector3i>> &, parallelSetParentDownList,
			UT_Array<UT_Array<UT_Vector3i>> &, parallelSetParentActiveList,
			UT_Array<bool> &, isTileOccupiedList,
			const int, level)

    void setFaceGradingPartial(UT_Array<UT_Array<UT_Vector3i>> &parallelSetParentDownList,
				UT_Array<UT_Array<UT_Vector3i>> &parallelSetParentActiveList,
				UT_Array<bool> &isTileOccupiedList,
				const int level,
				const UT_JobInfo &info);

    THREADED_METHOD3(SIM_OctreeGrid, myOctreeGridLabels[level].shouldMultiThread(),
			setParentsUp,
			UT_Array<UT_Array<UT_Vector3i>> &, parallelNewParentList,
			UT_Array<bool> &, isTileOccupiedList,
			const int, level)

    void setParentsUpPartial(UT_Array<UT_Array<UT_Vector3i>> &parallelNewParentList,
				UT_Array<bool> &isTileOccupiedList,
				const int level,
				const UT_JobInfo &info);

    THREADED_METHOD1(SIM_OctreeGrid, myOctreeGridLabels[level].shouldMultiThread(),
			setTopLevel,
			const int, level);

    void setTopLevelPartial(const int level, const UT_JobInfo &info);

    THREADED_METHOD2_CONST(SIM_OctreeGrid,
			    myOctreeGridLabels[level].shouldMultiThread(),
			    checkActiveCellAtLevel,
			    bool &, hasActiveCell,
			    const int, level);

    void checkActiveCellAtLevelPartial(bool &hasActiveCell,
					const int level,
					const UT_JobInfo &info) const;


    THREADED_METHOD2(SIM_OctreeGrid,
			newGridLabels.shouldMultiThread(),
			setGridFromParent,
			SIM_RawField &, newGridLabels,
			const int, level)

    void setGridFromParentPartial(SIM_RawField &newGridLabels,
				    const int level,
				    const UT_JobInfo &info);

    const UT_Vector3I cell_children[8] = {UT_Vector3I(0, 0, 0), UT_Vector3I(1, 0, 0),
					    UT_Vector3I(1, 0, 1), UT_Vector3I(0, 0, 1),
					    UT_Vector3I(0, 1, 0), UT_Vector3I(1, 1, 0),
					    UT_Vector3I(1, 1, 1), UT_Vector3I(0, 1, 1)};

    // Markers to label active cells, if the cell is "above"
    // or "below" the leaf (aka active cell).
    UT_Array<SIM_RawField> myOctreeGridLabels;
    UT_Vector3i myBaseVoxelRes;
    int myOctreeLevels;
};

#endif