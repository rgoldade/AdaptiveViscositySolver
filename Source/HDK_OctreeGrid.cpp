#include "HDK_OctreeGrid.h"
#include <GU/GU_Detail.h>

void
HDK_OctreeGrid::outputOctreeGeometry(GU_Detail &octreeGeometryDetail) const
{
    // Get pscale
    octreeGeometryDetail.clear();
    GA_RWHandleF gridSizeHandle(&octreeGeometryDetail, GA_ATTRIB_POINT, "dx");
    if (!gridSizeHandle.isValid())
    {
        octreeGeometryDetail.addFloatTuple(GA_ATTRIB_POINT, "dx", 1, GA_Defaults(0));
        gridSizeHandle = GA_RWHandleF(&octreeGeometryDetail, GA_ATTRIB_POINT, "dx");
        gridSizeHandle.bumpDataId();
    }

    GA_RWHandleI octreeLevelHandle(&octreeGeometryDetail, GA_ATTRIB_POINT, "octreeLevel");
    if (!octreeLevelHandle.isValid())
    {
        octreeGeometryDetail.addIntTuple(GA_ATTRIB_POINT, "octreeLevel", 1, GA_Defaults(-1));
        octreeLevelHandle = GA_RWHandleI(&octreeGeometryDetail, GA_ATTRIB_POINT, "octreeLevel");
        octreeLevelHandle.bumpDataId();
    }

    UT_Interrupt *boss = UTgetInterrupt();

    for (int level = 0; level < myOctreeLevels; ++level)
    {
        fpreal dx = myOctreeGridLabels[level].getVoxelSize().maxComponent();

        UT_VoxelArrayIteratorF vit;
        vit.setConstArray(myOctreeGridLabels[level].field());

        UT_VoxelTileIteratorF vitt;
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
                    {
                        UT_Vector3 cellCenter;
                        myOctreeGridLabels[level].indexToPos(vitt.x(),
							vitt.y(),
                                                        vitt.z(),
							cellCenter);

                        GA_Offset pointOffset = octreeGeometryDetail.appendPoint();
                        gridSizeHandle.set(pointOffset, dx);
                        octreeLevelHandle.set(pointOffset, level);
                        octreeGeometryDetail.setPos3(pointOffset, cellCenter);
                    }
                }
            }
        }
    }

    gridSizeHandle.bumpDataId();
    octreeLevelHandle.bumpDataId();
    octreeGeometryDetail.getAttributes().bumpAllDataIds(GA_ATTRIB_POINT);
}

void
HDK_OctreeGrid::setBaseGridLabelsPartial(const SIM_RawField &initializerMask,
						const UT_JobInfo &info)
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorF vit(myOctreeGridLabels[0].fieldNC());
    vit.setCompressOnExit(true);
    vit.splitByTile(info);

    UT_VoxelTileIteratorF vitt;

    UT_Vector3i voxelRes;
    initializerMask.getVoxelRes(voxelRes[0], voxelRes[1], voxelRes[2]);

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

	UT_Vector3i tileIndex;

        int tileNumber = vit.getLinearTileNum();
	myOctreeGridLabels[0].field()->linearTileToXYZ(tileNumber,
						    tileIndex[0],
						    tileIndex[1],
						    tileIndex[2]);

        // Check that tile is within the mask grid's bounds
        if (tileIndex[0] >= initializerMask.field()->getTileRes(0) ||
            tileIndex[1] >= initializerMask.field()->getTileRes(1) ||
            tileIndex[2] >= initializerMask.field()->getTileRes(2))
        {
            assert(vit.getValue() == INACTIVE);
            continue;
        }

        UT_VoxelTile<fpreal32> *maskTile = initializerMask.field()->getTile(tileIndex[0],
									    tileIndex[1],
									    tileIndex[2]);

        // Set constant tiles for marker grid
        if (maskTile->isConstant())
        {
            fpreal32 tileValue = (*maskTile)(0, 0, 0);
            if (tileValue == 0)
                vit.getTile()->makeConstant(ACTIVE);
            else if (tileValue < 0)
                vit.getTile()->makeConstant(UP);
            else
            {
                // Debug check that constant tile is inactive.
		// Everything should be set to inactive by default
		// so let's avoid writing to the tile.
                assert(vit.getValue() == INACTIVE);
            }
        }
        else
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                UT_Vector3i cell(vitt.x(), vitt.y(), vitt.z());

                if (cell[0] >= voxelRes[0] || cell[1] >= voxelRes[1] || cell[2] >= voxelRes[2])
                {
                    assert(vitt.getValue() == INACTIVE);
                    continue;
                }

                fpreal32 maskValue = HDKgetFieldValue(initializerMask, cell);

                if (maskValue == 0)
                    vitt.setValue(ACTIVE);
                else if (maskValue < 0)
                    vitt.setValue(UP);
                else // Value already set to inactive
                    assert(vitt.getValue() == INACTIVE);
            }
        }
    }
}

void
HDK_OctreeGrid::setActiveCellsAndParentListPartial(UT_Array<UT_Array<UT_Vector3i>> &parallelSetParentDownList,
							UT_Array<bool> &isTileOccupiedList,
							const int level,
							const UT_JobInfo &info)
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorF vit(myOctreeGridLabels[level].fieldNC());
    vit.setCompressOnExit(true);
    vit.splitByTile(info);

    UT_Array<UT_Vector3i> &localSetParentDownList = parallelSetParentDownList[info.job()];

    UT_Array<bool> localIsTileOccupiedList;
    localIsTileOccupiedList.setSize(isTileOccupiedList.entries());
    localIsTileOccupiedList.constant(false);

    UT_VoxelTileIteratorF vitt;

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        bool doProcessTile = false;

        if (!vit.isTileConstant() || vit.getValue() == ACTIVE)
            doProcessTile = true;
        else if (vit.getValue() == UP)
        {
            int tileNumber = vit.getLinearTileNum();

	    UT_Vector3i tileIndex;
            myOctreeGridLabels[level].field()->linearTileToXYZ(tileNumber,
								tileIndex[0],
								tileIndex[1],
								tileIndex[2]);

	    UT_Vector3I startTileIndex, endTileIndex;
	    vit.getTileVoxels(startTileIndex, endTileIndex);

	    for (int axis = 0; axis < 3 && !doProcessTile; ++axis)
		for (int direction : {0,1})
		{
		    UT_Vector3i adjacentTileIndex = HDKcellToCell(tileIndex, axis, direction);

		    bool doProcessTileBorder = false;

		    // If the adjacent tile is not constant or out of bounds,
		    // process voxels along the border and check for ACTIVE voxels.
		    // Otherwise, check the adjacent constant tile for ACTIVE labels.

		    if (adjacentTileIndex[axis] >= 0 &&
			adjacentTileIndex[axis] < myOctreeGridLabels[level].field()->getTileRes(axis))
		    {
			UT_VoxelTile<fpreal32> *adjacentTile = myOctreeGridLabels[level].field()->getTile(adjacentTileIndex[0],
													adjacentTileIndex[1],
													adjacentTileIndex[2]);

			if (adjacentTile->isConstant())
			{
			    if ((*adjacentTile)(0, 0, 0) == ACTIVE)
			    {
				doProcessTileBorder = true;
				break;
			    }
			}
			else
			    doProcessTileBorder = true;
		    }
		    else
			doProcessTileBorder = true;

		    // We don't yet know if we need to actually iterate through the
		    // current tile. Our last check is to see if there is an adjacent
		    // ACTIVE cell in the neighbouring tile.
		    if (doProcessTileBorder)
		    {
			auto checkAdjacentCellActivity = [&](const UT_Vector3I &start, const UT_Vector3I &end) -> bool
			{
			    UT_Vector3i cell = start;
			    for (; cell[0] < end[0]; ++cell[0])
				for (; cell[1] < end[1]; ++cell[1])
				    for (; cell[2] < end[2]; ++cell[2])
				    {
					UT_Vector3i adjacentCell = HDKcellToCell(cell, axis, direction);
					if (HDKgetFieldValue(myOctreeGridLabels[level], adjacentCell) == ACTIVE)
					    return true;
				    }

			    return false;
			};
			
			UT_Vector3i localStartIndex(startTileIndex[0], startTileIndex[1], startTileIndex[2]);
			UT_Vector3i localEndIndex(endTileIndex[0], endTileIndex[1], endTileIndex[2]);

			if (direction == 0)
			    endTileIndex[axis] = startTileIndex[axis] + 1;
			else
			    startTileIndex[axis] = endTileIndex[axis] - 1;

			if (checkAdjacentCellActivity(startTileIndex, endTileIndex))
			{
			    doProcessTile = true;
			    break;
			}
		    }
		}
        }

        if (doProcessTile)
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
		// If the current cell is labelled UP but a sibling cell is labelled ACTIVE,
		// then this cell must be re-labelled as ACTIVE.
                if (vitt.getValue() == UP)
                {
                    UT_Vector3i cell(vitt.x(), vitt.y(), vitt.z());
                    UT_Vector3I parentCell = getParentCell(cell);

		    for (int childIndex = 0; childIndex < 8; ++childIndex)
                    {
                        UT_Vector3I siblingCell = getChildCell(parentCell, childIndex);

			if (HDKgetFieldValue(myOctreeGridLabels[level], siblingCell) == ACTIVE)
                        {
                            vitt.setValue(ACTIVE);
                            break;
                        }
                    }
#if !defined(NDEBUG)
		    // There should never be an instance where a group of silbing cells
		    // contain both an UP and INACTIVE label.
                    for (int childIndex = 0; childIndex < 8; ++childIndex)
                    {
                        UT_Vector3I siblingCell = getChildCell(parentCell, childIndex);
                        assert(HDKgetFieldValue(myOctreeGridLabels[level], siblingCell) != INACTIVE);
                    }
#endif
                }
                // If the cell is labelled ACTIVE, the parent must be labelled DOWN.
		// The parent tile may be compressed so we cannot safely write to it
                // (another thread might be trying to do the same thing). Instead, we 
		// defer these operations by setting the parent's entire tile to be
		// uncompressed later and storing the parent's cell to be set DOWN.
                else if (vitt.getValue() == ACTIVE)
                {
                    UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());
                    UT_Vector3I parentCell = getParentCell(cell);

                    int parentTileNumber = myOctreeGridLabels[level + 1].field()->indexToLinearTile(parentCell[0],
												parentCell[1],
												parentCell[2]);

                    localIsTileOccupiedList[parentTileNumber] = true;

                    localSetParentDownList.append(parentCell);
                }
            }
        }
    }

    for (int i = 0; i < isTileOccupiedList.entries(); ++i)
    {
        if (localIsTileOccupiedList[i])
            isTileOccupiedList[i] = true;
    }
}

void
HDK_OctreeGrid::uncompressParentTilesPartial(const UT_Array<bool> &isTileOccupiedList,
                                             const int level,
                                             const UT_JobInfo &info)
{
    UT_Interrupt *boss = UTgetInterrupt();

    int start, end;
    int elements = isTileOccupiedList.entries();

    info.divideWork(elements, start, end);

    if (start == end)
        return;

    for (int i = start, iend = end; i < iend; ++i)
    {
        if (!(i & 127))
        {
            if (boss->opInterrupt())
                break;
        }

        if (isTileOccupiedList[i])
        {
            myOctreeGridLabels[level + 1].field()->getLinearTile(i)->uncompress();
        }
    }
}

void
HDK_OctreeGrid::setParentCellLabelPartial(const UT_Array<UT_Vector3i> &newParentCellList,
						const int level,
						const OctreeCellLabel label,
						const UT_JobInfo &info)
{
    UT_Interrupt *boss = UTgetInterrupt();

    int start, end;
    int elements = newParentCellList.entries();

    info.divideWork(elements, start, end);

    if (start == end)
        return;

    // If the uncompress stage exited early, it's not safe to continue.
    if (boss->opInterrupt())
        return;

    // It's possible to have duplicate entries since multiple children could
    // be trying to set the same label for a parent. Check that the value at
    // the start isn't also included in the previous work group.

    int localStart = start;
    const int localEnd = end;

    if (localStart > 0 && newParentCellList[localStart] == newParentCellList[localStart - 1])
    {
        UT_Vector3i localCell = newParentCellList[localStart];

        while (localStart < localEnd && localCell == newParentCellList[localStart])
            ++localStart;

        if (localStart == localEnd)
            return;
    }

    UT_Vector3i oldCell(-1, -1, -1);

    for (int i = localStart; i < localEnd; ++i)
    {
        if (!(i & 127))
        {
            if (boss->opInterrupt())
                break;
        }

        UT_Vector3i parentCell = newParentCellList[i];

        // We could have duplicates in the list. We only need to set the value
        // once per unique parent index.
        if (parentCell != oldCell)
	    HDKsetFieldValue(myOctreeGridLabels[level + 1], parentCell, label);

        oldCell = parentCell;
    }
}

void
HDK_OctreeGrid::setFaceGradingPartial(UT_Array<UT_Array<UT_Vector3i>> &parallelSetParentDownList,
					    UT_Array<UT_Array<UT_Vector3i>> &parallelSetParentActiveList,
					    UT_Array<bool> &isTileOccupiedList,
					    const int level, const UT_JobInfo &info)
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorF vit;
    vit.setConstArray(myOctreeGridLabels[level].field());
    vit.splitByTile(info);

    UT_VoxelTileIteratorF vitt;

    UT_Array<UT_Vector3i> &localSetParentDownList = parallelSetParentDownList[info.job()];
    UT_Array<UT_Vector3i> &localSetParentActiveList = parallelSetParentActiveList[info.job()];

    UT_Array<bool> localIsTileOccupiedList;
    localIsTileOccupiedList.setSize(isTileOccupiedList.entries());
    localIsTileOccupiedList.constant(false);

    UT_Vector3i voxelRes;
    myOctreeGridLabels[level].getVoxelRes(voxelRes[0], voxelRes[1], voxelRes[2]);

    UT_Vector3i parentGridSize;
    myOctreeGridLabels[level + 1].getVoxelRes(parentGridSize[0],
					    parentGridSize[1],
					    parentGridSize[2]);

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        // We can skip UP and INACTIVE tiles since there's nothing to do forthem.
        if (!vit.isTileConstant() || vit.getValue() != UP || vit.getValue() != INACTIVE)
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

                int activity = vitt.getValue();
                if (activity == ACTIVE)
                {
                    for (int axis : {0,1,2})
			for (int direction : {0,1})
			{
			    UT_Vector3i adjacentCell = HDKcellToCell(cell, axis, direction);

			    if (adjacentCell[axis] < 0 || adjacentCell[axis] >= voxelRes[axis])
				continue;

			    if (HDKgetFieldValue(myOctreeGridLabels[level], adjacentCell) == UP)
			    {
				UT_Vector3i parentCell = getParentCell(adjacentCell);

				// This adjacent cell can't be a sibling because if it was UP,
				// then it should have be set to ACTIVE already.
				assert(parentCell != getParentCell(cell));

				// The parent should be INACTIVE because it should not have been
				// written to yet (meaning the parent has no ACTIVE children).
				assert(HDKgetFieldValue(myOctreeGridLabels[level + 1], parentCell) == INACTIVE);

				int tileNumber = myOctreeGridLabels[level + 1].field()->indexToLinearTile(parentCell[0],
													parentCell[1],
													parentCell[2]);
				localIsTileOccupiedList[tileNumber] = true;

				localSetParentActiveList.append(parentCell);
			    }
                        }
                }
                else if (activity == DOWN)
                {
                    UT_Vector3i parentCell = getParentCell(cell);

#if !defined(NDEBUG)
			for (int childIndex = 0; childIndex < 8; ++childIndex)
			{
			    UT_Vector3i siblingCell = getChildCell(parentCell, childIndex);
			    assert(HDKgetFieldValue(myOctreeGridLabels[level], siblingCell) != UP);
			}
#endif

		    int tileNumber = myOctreeGridLabels[level + 1].field()->indexToLinearTile(parentCell[0],
											    parentCell[1],
											    parentCell[2]);
                    localIsTileOccupiedList[tileNumber] = true;

                    localSetParentDownList.append(parentCell);
                }
            }
        }
    }

    for (int i = 0; i < isTileOccupiedList.entries(); ++i)
    {
        if (localIsTileOccupiedList[i])
	    isTileOccupiedList[i] = true;
    }
}

void
HDK_OctreeGrid::setParentsUpPartial(UT_Array<UT_Array<UT_Vector3i>> &parallelNewParentList,
				    UT_Array<bool> &isTileOccupiedList,
				    const int level,
				    const UT_JobInfo &info)
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorF vit;
    vit.setConstArray(myOctreeGridLabels[level].field());
    vit.splitByTile(info);

    UT_VoxelTileIteratorF vitt;

    UT_Array<UT_Vector3i> &localNewParentList = parallelNewParentList[info.job()];

    UT_Array<bool> localIsTileOccupiedList;
    localIsTileOccupiedList.setSize(isTileOccupiedList.entries());
    localIsTileOccupiedList.constant(false);

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        if (!vit.isTileConstant() || vit.getValue() == UP)
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                if (vitt.getValue() == UP)
                {
                    UT_Vector3i cell(vitt.x(), vitt.y(), vitt.z());
                    UT_Vector3i parentCell = getParentCell(cell);

                    // If the parent is uninitialized (i.e., INACTIVE), then
                    // pass the "UP" label from the child.
		    if (HDKgetFieldValue(myOctreeGridLabels[level + 1], parentCell)  == INACTIVE)
                    {
                        int tileNumber = myOctreeGridLabels[level + 1].field()->indexToLinearTile(parentCell[0],
												parentCell[1],
												parentCell[2]);
                        localIsTileOccupiedList[tileNumber] = true;

                        localNewParentList.append(parentCell);
                    }
                    // Since the child is UP, the parent cannot be DOWN. This
                    // would imply that a sibling is ACTIVE and that would mean
                    // this child should be ACTIVE from the first pass. Also,
                    // the parent cannot be UP since we are only just setting
                    // parents to UP now. The only other option is ACTIVE since
                    // we've capture the INACTIVE event above.
                    else
			assert(HDKgetFieldValue(myOctreeGridLabels[level + 1], parentCell) == ACTIVE);

#if !defined(NDEBUG)
                    // From the first pass, if the sibling is ACTIVE then this
                    // child should already be set to ACTIVE. Because INACTIVE
                    // cells are "OUTSIDE" of the simulation domain, UP cells
                    // are contained inside a boundary of "ACTIVE" cells.
                    // Therefore a UP child should never have an INACTIVE
                    // sibling. Finally, an UP child cannot have a down sibling
                    // because that would violate the face grading rule. If the
                    // child is currently UP then all of its siblings must be UP
                    // as well.

		    for (int childIndex = 0; childIndex < 8; ++childIndex)
                    {
                        UT_Vector3I siblingCell = getChildCell(parentCell, childIndex);

			assert(HDKgetFieldValue(myOctreeGridLabels[level], siblingCell) == UP);
                    }
#endif
                }
            }
        }
    }

    for (int i = 0; i < isTileOccupiedList.entries(); ++i)
    {
        if (localIsTileOccupiedList[i])
            isTileOccupiedList[i] = true;
    }
}

// Any voxels at the top of the pyramid that are still labeled as UP must be set to ACTIVE.
void
HDK_OctreeGrid::setTopLevelPartial(const int level, const UT_JobInfo &info)
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorF vit(myOctreeGridLabels[level].fieldNC());
    vit.setCompressOnExit(true);
    vit.splitByTile(info);

    UT_VoxelTileIteratorF vitt;

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        if (vit.isTileConstant())
        {
            if (vit.getValue() == UP)
                vit.getTile()->makeConstant(ACTIVE);
        }
        else
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                if (vitt.getValue() == UP)
                    vitt.setValue(ACTIVE);
            }
        }
    }
}

void
HDK_OctreeGrid::checkActiveCellAtLevelPartial(bool &hasActiveCell,
						const int level,
						const UT_JobInfo &info) const
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorF vit;
    vit.setConstArray(myOctreeGridLabels[level].field());
    vit.splitByTile(info);

    UT_VoxelTileIteratorF vitt;

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        if (hasActiveCell)
            return;

        if (vit.isTileConstant())
        {
            if (vit.getValue() == ACTIVE)
            {
                hasActiveCell = true;
                return;
            }
        }
        else
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                if (vitt.getValue() == ACTIVE)
                {
                    hasActiveCell = true;
                    return;
                }
            }
        }
    }
}

UT_Array<UT_Vector4i> 
HDK_OctreeGrid::getFaceAdjacentCells(const UT_Vector3i &cell, const int axis,
					const int direction, const int level) const
{
    assert(axis < 3 && direction < 2);
    // TODO: replace with UT_Array
    UT_Array<UT_Vector4i> adjacentCellList;

    UT_Vector3i adjacentCell = HDKcellToCell(cell, axis, direction);

    OctreeCellLabel adjacentLabel = getCellLabel(adjacentCell, level);

    if (adjacentLabel == ACTIVE)
	adjacentCellList.append(UT_Vector4i(adjacentCell[0],
						adjacentCell[1],
						adjacentCell[2],
						level));
    else if (adjacentLabel == UP)
    {
	UT_Vector3i parentCell = getParentCell(adjacentCell);

	assert(isCellActive(parentCell, level + 1));

        adjacentCellList.append(UT_Vector4i(parentCell[0],
						parentCell[1],
						parentCell[2],
						level + 1));
    }
    else if (adjacentLabel == DOWN)
    {
	// Find each child cell that is face-adjacent to the original cell
	for (int secondDirection : {0,1})
	    for (int thirdDirection : {0,1})
	    {
		int childIndex = 0;
		if (direction == 0)
		    childIndex += (1 << axis);
		if (secondDirection == 1)
		    childIndex += (1 << (axis + 1) % 3);
		if (thirdDirection == 1)
		    childIndex += (1 << (axis + 2) % 3);
		
		UT_Vector3i childCell = getChildCell(adjacentCell, childIndex);

		if (isCellActive(childCell, level - 1))
		    adjacentCellList.append(UT_Vector4i(childCell[0],
							    childCell[1],
							    childCell[2],
							    level - 1));
		else
		    assert(HDKgetFieldValue(myOctreeGridLabels[level - 1], childCell) == INACTIVE);

	    }
    }

    return adjacentCellList;
}

//
// Unit test methods
//

// An INACTIVE fine cell can only have INACTIVE or DOWN ancestors.
// An ACTIVE fine cell can only have DOWN ancestors.
// An UP fine cell has one and only one ACTIVE ancestor.
// A DOWN fine cell should never happen.
void
HDK_OctreeGrid::activeCountUnitTestPartial(bool &passed,
					    const UT_JobInfo &info) const
{
    UT_Interrupt *boss = UTgetInterrupt();

    int octreeLevels = getOctreeLevels();

    UT_VoxelArrayIteratorF vit;
    vit.setConstArray(myOctreeGridLabels[0].field());
    vit.splitByTile(info);

    UT_VoxelTileIteratorF vitt;

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        vitt.setTile(vit);

        for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
        {
            if (!passed)
                return;

            UT_Vector3I cell(vitt.x(), vitt.y(), vitt.z());

            if (vitt.getValue() == INACTIVE)
            {
                bool foundDownCell = false;
                for (int level = 1; level < octreeLevels; ++level)
                {
                    cell = getParentCell(cell);
		    OctreeCellLabel cellLabel = getCellLabel(cell, level);

                    if (cellLabel == DOWN)
                        foundDownCell = true;
                    else if (cellLabel == INACTIVE)
                    {
                        if (foundDownCell)
                            passed = false;
                    }
                    else
                        passed = false;
                }
            }
            else if (vitt.getValue() == ACTIVE)
            {
                for (int level = 1; level < octreeLevels; ++level)
                {
                    cell = getParentCell(cell);
                    OctreeCellLabel cellLabel = getCellLabel(cell, level);

                    if (cellLabel != DOWN)
                        passed = false;
                }
            }
            else if (vitt.getValue() == UP)
            {
                bool foundActiveCell = false;

                for (int level = 1; level < octreeLevels; ++level)
                {
                    cell = getParentCell(cell);
                    OctreeCellLabel cellLabel = getCellLabel(cell, level);

                    if (cellLabel == ACTIVE)
                    {
                        if (foundActiveCell)
                            passed = false;

                        foundActiveCell = true;
                    }
                    else if (cellLabel == UP)
                    {
                        if (foundActiveCell)
                            passed = false;
                    }
                    else if (cellLabel == DOWN)
                    {
                        if (!foundActiveCell)
                            passed = false;
                    }
                    else
                        passed = false;
                }
            }
            else // DOWN should not be set at the finest level
                passed = false;
       }
    }
}

// An UP cell can only have an UP sibling and an ACTIVE or UP
// adjacent cell.
void
HDK_OctreeGrid::upAdjacentUnitTestPartial(bool &passed,
					    const int level,
					    const UT_JobInfo &info) const
{
    UT_Interrupt *boss = UTgetInterrupt();

    int octreeLevels = getOctreeLevels();

    UT_Vector3i voxelRes;
    myOctreeGridLabels[level].getVoxelRes(voxelRes[0],
					    voxelRes[1],
					    voxelRes[2]);

    UT_VoxelArrayIteratorF vit;
    vit.setConstArray(myOctreeGridLabels[level].field());
    vit.splitByTile(info);

    UT_VoxelTileIteratorF vitt;

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        if (!vit.isTileConstant() || vit.getValue() == UP)
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                if (vitt.getValue() == UP)
                {
                    if (level == octreeLevels)
                    {
                        passed = false;
                        return;
                    }

                    UT_Vector3i cell(vitt.x(), vitt.y(), vitt.z());

                    UT_Vector3i parentCell = getParentCell(cell);

                    // Check siblings
		    for (int childIndex = 0; childIndex < 8; ++childIndex)
                    {
                        UT_Vector3i siblingCell = getChildCell(parentCell, childIndex);

			if (getCellLabel(siblingCell, level) != UP)
                        {
                            passed = false;
                            return;
                        }
                    }

                    // Check adjacent cells
                    for (int axis : {0,1,2})
			for (int direction : {0,1})
			{
			    UT_Vector3i adjacentCell = HDKcellToCell(cell, axis, direction);
			    
			    if (adjacentCell[axis] < 0 || adjacentCell[axis] >= voxelRes[axis])
				continue;

			    OctreeCellLabel adjacentLabel = getCellLabel(adjacentCell, level);

			    if (adjacentLabel != ACTIVE && adjacentLabel != UP)
			    {
				passed = false;
				return;
			    }
			}
                }
            }
        }
    }
}

// If the adjacent cell is ACTIVE, things are good.
// If the adjacent cell is DOWN, then the adjacent child cells must be ACTIVE.
// It the adjacent cell is UP, then the adjacent parent must be ACTIVE.
// The adjacent cells must also reciprocate.
void
HDK_OctreeGrid::activeUnitTestPartial(bool &passed,
					const int level,
					const UT_JobInfo &info) const
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_Vector3i voxelRes = getVoxelRes(level);

    UT_VoxelArrayIteratorF vit;
    vit.setConstArray(myOctreeGridLabels[level].field());
    vit.splitByTile(info);

    UT_VoxelTileIteratorF vitt;

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        if (!vit.isTileConstant() || vit.getValue() == ACTIVE)
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                if (vitt.getValue() == ACTIVE && !passed)
                {
                    const UT_Vector3i cell(vitt.x(), vitt.y(), vitt.z());

		    for (int axis : {0,1,2})
			for (int direction : {0,1})
			{
			    UT_Vector3i adjacentCell = HDKcellToCell(cell, axis, direction);

			    if (adjacentCell[axis] < 0 || adjacentCell[axis] >= voxelRes[axis])
				continue;

			    auto adjacentFaceCellList = getFaceAdjacentCells(cell, axis, direction, level);
			    OctreeCellLabel adjacentLabel = getCellLabel(adjacentCell, level);

			    if (adjacentLabel == DOWN)
			    {
				if (adjacentFaceCellList.size() != 4)
				{
				    passed = false;
				    return;
				}

				// Check that adjacent cells are ACTIVE. This
				// enforces face grading.
				for (auto cellItem : adjacentFaceCellList)
				{
				    UT_Vector3i localCell(cellItem[0], cellItem[1], cellItem[2]);
				    assert(cellItem[3] == level - 1);

				    OctreeCellLabel nephewLabel = getCellLabel(localCell, level -1);
				    if (nephewLabel != ACTIVE)
				    {
					passed = false;
					return;
				    }
				}
			    }
			    // The parent activity check enforces face grading
			    else if (adjacentLabel == UP)
			    {
				if (level == getOctreeLevels() - 1)
				{
				    passed = false;
				    return;
				}

				UT_Vector3I parentCell = getParentCell(adjacentCell);

				OctreeCellLabel parentLabel = getCellLabel(parentCell, level + 1);

				if (parentLabel != ACTIVE)
				{
				    passed = false;
				    return;
				}
			    }

			    // Check that adjacent active cells also reciprocate
			    for (auto cellItem : adjacentFaceCellList)
			    {
				UT_Vector3i localCell(cellItem[0], cellItem[1], cellItem[2]);

				int otherDirection = (direction + 1) % 2;
				int localLevel = cellItem[3];

				auto reciprocatingCellList = getFaceAdjacentCells(localCell, axis, otherDirection, localLevel);

				auto result = std::find(reciprocatingCellList.begin(),
							reciprocatingCellList.end(),
							UT_Vector4i(cell[0], cell[1], cell[2], level));

				if (result == reciprocatingCellList.end())
				{
				    passed = false;
				    return;
				}
			    }
			}
                }
            }
        }
    }
}

bool
HDK_OctreeGrid::unitTest() const
{
    bool passed = true;
    activeCountUnitTest(passed);

    if (!passed)
        return false;

    int octreeLevels = getOctreeLevels();
    for (int level = 0; level < octreeLevels; ++level)
    {
        upAdjacentUnitTest(passed, level);

        if (!passed)
            return false;
    }

    for (int level = 0; level < octreeLevels; ++level)
    {
        activeUnitTest(passed, level);

        if (!passed)
            return false;
    }

    return true;
}

void
HDK_OctreeGrid::refineGrid()
{
    const UT_Vector3 size = getSize();
    const UT_Vector3 origin = getOrig();

    const int octreeLevels = getOctreeLevels();

    UT_Array<SIM_RawField> localOctreeGridLabels;
    localOctreeGridLabels.setSize(octreeLevels);

    UT_Vector3i voxelRes = getVoxelRes(0);
    voxelRes *= 2;

    myBaseVoxelRes = voxelRes;
    for (int level = 0; level < octreeLevels; ++level)
    {
	localOctreeGridLabels[level].init(SIM_SAMPLE_CENTER, origin, size, voxelRes[0], voxelRes[1], voxelRes[2]);
	localOctreeGridLabels[level].makeConstant(INACTIVE);
	voxelRes /= 2;
    }

    for (int level = 0; level < octreeLevels; ++level)
	setGridFromParent(localOctreeGridLabels[level], level);

    myOctreeGridLabels = std::move(localOctreeGridLabels);
}

void
HDK_OctreeGrid::setGridFromParentPartial(SIM_RawField &newGridLabels,
					    const int level,
					    const UT_JobInfo &info)
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorF vit(newGridLabels.fieldNC());
    vit.setCompressOnExit(true);
    vit.splitByTile(info);

    UT_VoxelTileIteratorF vitt;

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

	vitt.setTile(vit);

	for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
	{
	    UT_Vector3i cell(vitt.x(), vitt.y(), vitt.z());

	    UT_Vector3i parentCell = getParentCell(cell);
	    vitt.setValue(HDKgetFieldValue(myOctreeGridLabels[level], parentCell));
	}
    }
}