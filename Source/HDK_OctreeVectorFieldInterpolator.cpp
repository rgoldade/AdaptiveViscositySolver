#include "HDK_OctreeVectorFieldInterpolator.h"

#include <SYS/SYS_Math.h>

void
HDK_OctreeVectorFieldInterpolator::findOccupiedNodeTilesPartial(UT_Array<bool> &isTileOccupiedList,
								const SIM_RawField &octreeLabels,
								const int level,
								const UT_JobInfo &info) const
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorF vit;
    vit.setConstArray(octreeLabels.field());
    vit.splitByTile(info);

    const exint tileCount = isTileOccupiedList.entries();

    UT_Array<bool> localIsTileOccupiedList;
    localIsTileOccupiedList.setSize(tileCount);
    localIsTileOccupiedList.constant(false);

    UT_VoxelTileIteratorF vitt;

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        if (!vit.isTileConstant() || vit.getValue() == HDK_OctreeGrid::ACTIVE)
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                if (vitt.getValue() == HDK_OctreeGrid::ACTIVE)
                {
                    UT_Vector3i cell(vitt.x(), vitt.y(), vitt.z());

                    for (int nodeIndex = 0; nodeIndex < 8; ++nodeIndex)
                    {
			UT_Vector3i node = HDKcellToNode(cell, nodeIndex);

                        const exint tileNumber = myNodeLabels[level].field()->indexToLinearTile(node[0],
												node[1],
												node[2]);

                        assert(tileNumber == myNodeValues[level][0].field()->indexToLinearTile(node[0], node[1], node[2]) &&
                               tileNumber == myNodeValues[level][1].field()->indexToLinearTile(node[0], node[1], node[2]) &&
                               tileNumber == myNodeValues[level][2].field()->indexToLinearTile(node[0], node[1], node[2]));

                        localIsTileOccupiedList[tileNumber] = true;
                    }
                }
            }
        }
    }

    for (exint tileNumber = 0; tileNumber < tileCount; ++tileNumber)
    {
	if (localIsTileOccupiedList[tileNumber])
	    isTileOccupiedList[tileNumber] = true;
    }
}

void
HDK_OctreeVectorFieldInterpolator::uncompressTilesPartial(SIM_RawField &grid,
							    const UT_Array<bool> &isTileOccupiedList,
							    const UT_JobInfo &info)
{
    UT_Interrupt *boss = UTgetInterrupt();

    exint start, end;
    exint elements = isTileOccupiedList.entries();

    info.divideWork(elements, start, end);

    const exint localEnd = end;
    for (exint i = start; i < localEnd; ++i)
    {
        if (!(i & 127))
        {
            if (boss->opInterrupt())
                break;
        }

        if (isTileOccupiedList[i])
            grid.field()->getLinearTile(i)->uncompress();
    }
}

void
HDK_OctreeVectorFieldInterpolator::uncompressIndexTilesPartial(SIM_RawIndexField &grid,
								const UT_Array<bool> &isTileOccupiedList,
								const UT_JobInfo &info)
{
    UT_Interrupt *boss = UTgetInterrupt();

    exint start, end;
    exint elements = isTileOccupiedList.entries();

    info.divideWork(elements, start, end);

    const exint localEnd = end;
    for (exint i = start; i < localEnd; ++i)
    {
        if (!(i & 127))
        {
            if (boss->opInterrupt())
                break;
        }

        if (isTileOccupiedList[i])
            grid.field()->getLinearTile(i)->uncompress();
    }
}

void
HDK_OctreeVectorFieldInterpolator::setActiveNodesPartial(const int level,
                                                         const UT_JobInfo &info)
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorF vit(myNodeLabels[level].fieldNC());
    vit.setCompressOnExit(true);
    vit.splitByTile(info);

    UT_VoxelTileIteratorF vitt;

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        if (!vit.isTileConstant())
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                UT_Vector3i node(vitt.x(), vitt.y(), vitt.z());

                bool isNodeActive = false;
                bool isNodeInactive = false;

                for (int faceAxis = 0; !isNodeInactive && faceAxis < 3; ++faceAxis)
                {
		    UT_Vector3i faceRes;
                    faceRes[0] = myVelocityIndices[level][faceAxis].getXRes();
                    faceRes[1] = myVelocityIndices[level][faceAxis].getYRes();
                    faceRes[2] = myVelocityIndices[level][faceAxis].getZRes();

		    for (int faceIndex = 0; faceIndex < 4; ++faceIndex)
                    {
			const UT_Vector3i face = HDKnodeToFace(node, faceAxis, faceIndex);

                        if (face[(faceAxis + 1) % 3] < 0 ||
                            face[(faceAxis + 2) % 3] < 0 ||
                            face[(faceAxis + 1) % 3] >= faceRes[(faceAxis + 1) % 3] ||
                            face[(faceAxis + 2) % 3] >= faceRes[(faceAxis + 2) % 3])
                        {
                            isNodeInactive = true;
                            continue;
                        }

                        const exint velocityIndex = HDKgetFieldValue(myVelocityIndices[level][faceAxis], face);

                        if (velocityIndex >= 0)
                            isNodeActive = true;
                        // If the node touches a solid or "outside" face, it
                        // MUST not be near a transition. This means any
                        // adjacent faces that are active must be on the uniform
                        // level and will be referenced directly instead of
                        // through interpolation.
                        else if (velocityIndex == HDK_SOLIDBOUNDARY || velocityIndex == HDK_OUTSIDE)
                        {
                            isNodeInactive = true;
                            break;
                        }
                    }
                }

                if (isNodeActive && !isNodeInactive)
		    vitt.setValue(ACTIVENODE);
            }
        }
    }
}

void
HDK_OctreeVectorFieldInterpolator::sampleActiveNodesPartial(UT_Array<SIM_RawField> &nodeWeights,
							    SIM_RawIndexField &nodeFlags,
							    const int level,
							    const UT_JobInfo &info)
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorF vit;
    vit.setConstArray(myNodeLabels[level].field());
    vit.splitByTile(info);

    UT_VoxelTileIteratorF vitt;

    const int octreeLevels = myOctreeLabels.getOctreeLevels();

    // We need a weight to emphasize contributions from smaller, closer faces.
    const fpreal weight = fpreal(1 << (octreeLevels - level - 1));

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        if (!vit.isTileConstant() || vit.getValue() == ACTIVENODE)
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                if (vitt.getValue() == ACTIVENODE)
                {
                    UT_Vector3i node(vitt.x(), vitt.y(), vitt.z());

                    // In order to keep track of what adjacent faces have been
                    // accounted for at the various levels, we keep a bit-flag
                    // and accumulate all the bits later.
                    exint localFlag = 0;

                    for (int faceAxis : {0,1,2})
                    {
                        UT_Vector3i faceRes;
                        faceRes[0] = myVelocityIndices[level][faceAxis].getXRes();
                        faceRes[1] = myVelocityIndices[level][faceAxis].getYRes();
                        faceRes[2] = myVelocityIndices[level][faceAxis].getZRes();

                        fpreal accumulatedValue = 0.;
                        fpreal accumulatedWeight = 0.;

                        // Grab the adjacent faces of their ancestors to build
                        // the sample at the node
			for (int faceIndex = 0; faceIndex < 4; ++faceIndex)
                        {
			    const UT_Vector3i face = HDKnodeToFace(node, faceAxis, faceIndex);

			    if (face[(faceAxis + 1) % 3] < 0 ||
				face[(faceAxis + 2) % 3] < 0 ||
				face[(faceAxis + 1) % 3] >= faceRes[(faceAxis + 1) % 3] ||
				face[(faceAxis + 2) % 3] >= faceRes[(faceAxis + 2) % 3])
			    {
				localFlag += exint(1 << (faceAxis * 4 + faceIndex));
				accumulatedWeight += weight;

				continue;
			    }

			    assert(face[faceAxis] >= 0 && face[faceAxis] < faceRes[faceAxis]);

                            // Is the directly referenced face active?
			    const exint velocityIndex = HDKgetFieldValue(myVelocityIndices[level][faceAxis], face);
			    
                            if (velocityIndex >= 0)
                            {
				accumulatedValue += weight * HDKgetFieldValue(myVelocity[level][faceAxis], face);
                                accumulatedWeight += weight;
                                localFlag += exint(1 << (faceAxis * 4 + faceIndex));
                            }
                            // This face is neither active or unassigned, it must be specifically
			    // set to be inactive so we don't count its influence but flag it as finished.
                            else if (velocityIndex != HDK_UNASSIGNED)
                            {
                                assert(level == 0);
                                accumulatedWeight += weight;
                                localFlag += exint(1 << (faceAxis * 4 + faceIndex));
                            }
                        }

			HDKsetFieldValue(myNodeValues[level][faceAxis], node, accumulatedValue);
			HDKsetFieldValue(nodeWeights[faceAxis], node, accumulatedWeight);
                    }

		    HDKsetFieldValue(nodeFlags, node, localFlag);
                }
            }
        }
    }
}

void
HDK_OctreeVectorFieldInterpolator::bubbleActiveNodeValuesPartial(UT_Array<UT_Array<SIM_RawField>> &nodeWeights,
								    UT_Array<SIM_RawIndexField> &nodeFlags,
								    const int level,
								    const UT_JobInfo &info)
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorF vit(myNodeLabels[level].fieldNC());
    vit.setCompressOnExit(true);
    vit.splitByTile(info);

    UT_VoxelTileIteratorF vitt;

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        if (!vit.isTileConstant() || vit.getValue() == ACTIVENODE)
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                if (vitt.getValue() == ACTIVENODE)
                {
                    UT_Vector3i node(vitt.x(), vitt.y(), vitt.z());

                    const int oddIndexCount = (node[0] % 2) + (node[1] % 2) + (node[2] % 2);

                    // Co-located parent nodes can only exist at even indices.
                    if (oddIndexCount == 0)
                    {
			UT_Vector3i parentNode = myOctreeLabels.getParentNode(node);

			if (HDKgetFieldValue(myNodeLabels[level + 1], parentNode) == ACTIVENODE)
                        {
			    exint flag = HDKgetFieldValue(nodeFlags[level], node);

			    exint parentFlag = HDKgetFieldValue(nodeFlags[level + 1], parentNode);

			    HDKsetFieldValue(nodeFlags[level + 1], parentNode, flag + parentFlag);

                            assert(!(flag & parentFlag));

                            for (int axis : {0,1,2})
                            {
				fpreal weight = HDKgetFieldValue(nodeWeights[level][axis], node);
				fpreal value = HDKgetFieldValue(myNodeValues[level][axis], node);

				fpreal parentWeight = HDKgetFieldValue(nodeWeights[level + 1][axis], parentNode);
				fpreal parentValue = HDKgetFieldValue(myNodeValues[level + 1][axis], parentNode);

				HDKsetFieldValue(nodeWeights[level + 1][axis], parentNode, weight + parentWeight);
				HDKsetFieldValue(myNodeValues[level + 1][axis], parentNode, value + parentValue);
                            }

                            // Because the parent has accumulated the contribution,
			    // this node needs to be marked as follower.
                            vitt.setValue(DEPENDENTNODE);
                        }
                    }
                }
            }
        }
    }
}

void
HDK_OctreeVectorFieldInterpolator::finishIncompleteNodesPartial(UT_Array<UT_Array<SIM_RawField>> &nodeWeights,
								UT_Array<SIM_RawIndexField> &nodeFlags,
								const int level,
								const UT_JobInfo &info)
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorF vit;
    vit.setConstArray(myNodeLabels[level].field());
    vit.splitByTile(info);

    UT_VoxelTileIteratorF vitt;

    const int octreeLevels = myOctreeLabels.getOctreeLevels();

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        if (!vit.isTileConstant() || vit.getValue() == ACTIVENODE)
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                if (vitt.getValue() == ACTIVENODE)
                {
                    UT_Vector3i node(vitt.x(), vitt.y(), vitt.z());

		    exint localFlag = HDKgetFieldValue(nodeFlags[level], node);

                    // Only nodes in split edges and T-junction nodes should
                    // have incomplete flags after bubbling up the node values.
                    if (localFlag != 0xFFF)
                    {
#if !defined(NDEBUG)
			const int oddIndexCount = (node[0] % 2) + (node[1] % 2) + (node[2] % 2);
                        assert(oddIndexCount != 0 && oddIndexCount != 3);

                        int upCellCount = 0;
			UT_Vector3i voxelRes = myOctreeLabels.getVoxelRes(level);

			for (int cellIndex = 0; cellIndex < 8; ++cellIndex)
                        {
			    UT_Vector3i cell = HDKnodeToCell(node, cellIndex);

                            if (cell[0] < 0 || cell[1] < 0 || cell[2] < 0 ||
				cell[0] >= voxelRes[0] || cell[1] >= voxelRes[1] || cell[2] >= voxelRes[2])
                                continue;

			    if (myOctreeLabels.getCellLabel(cell, level) == HDK_OctreeGrid::UP)
                                ++upCellCount;
                        }

                        if (upCellCount > 0)
                        {
			    // A completely odd node should not have any adjacent up cells
                            assert(oddIndexCount < 3);

			    // It would be weird if there was only one up pointing cell
                            assert(upCellCount >= 2);

                            if (upCellCount != 2 && upCellCount != 4 && upCellCount != 6)
                                assert(false);
                        }
#endif
                        int bitShiftCount = 0;
                        exint tempFlag = localFlag;
			// Only check 12 bits because there are only 12 adjacent faces to a node.
                        while (localFlag != 0xFFF && bitShiftCount < 12)
                        {
                            if (!(tempFlag & 0x1))
                            {
                                int faceAxis = bitShiftCount / 4;
                                exint faceIndex = bitShiftCount % 4;

                                // A corner node might have a valid adjacent
                                // face that is one level higher. We should
                                // check there first before trying to
                                // interpolate from a coarse cell.
                                bool foundAdjacentFace = false;

				if (node[faceAxis] % 2 == 0)
                                {
				    UT_Vector3i face = HDKnodeToFace(node, faceAxis, faceIndex);
				    UT_Vector3i parentFace = myOctreeLabels.getParentFace(face);

				    const exint velocityIndex = HDKgetFieldValue(myVelocityIndices[level + 1][faceAxis], parentFace);

                                    if (velocityIndex >= 0)
                                    {
					fpreal ghostVelocity = HDKgetFieldValue(myVelocity[level + 1][faceAxis], parentFace);

					fpreal localNodeValue = HDKgetFieldValue(myNodeValues[level][faceAxis], node);
					fpreal weight = fpreal(1 << (octreeLevels - level - 1));

					localNodeValue += weight * ghostVelocity;

					HDKsetFieldValue(myNodeValues[level][faceAxis], node, localNodeValue);
                                        
					fpreal localNodeWeight = HDKgetFieldValue(nodeWeights[level][faceAxis], node);
                                        localNodeWeight += weight;

					HDKsetFieldValue(nodeWeights[level][faceAxis], node, localNodeWeight);

                                        localFlag += (1 << bitShiftCount);
                                        foundAdjacentFace = true;
                                    }
                                }

                                if (!foundAdjacentFace)
                                {
                                    // The axis that did not find a valid face MUST not be even.
                                    assert((node[faceAxis] % 2) != 0);

				    UT_Vector3i face = HDKnodeToFace(node, faceAxis, faceIndex);

                                    // Find the active cell that the face falls into
				    UT_Vector3i cell = HDKfaceToCell(face, faceAxis, 1);
				    assert(!myOctreeLabels.isCellActive(cell, level));

#if !defined(NDEBUG)
				    UT_Vector3i tempCell = HDKfaceToCell(face, faceAxis, 0);
                                    assert(!myOctreeLabels.isCellActive(tempCell, level));
#endif

                                    int searchLevel = level;
				    while (!myOctreeLabels.isCellActive(cell, searchLevel))
                                    {
                                        cell = myOctreeLabels.getParentCell(cell);
                                        ++searchLevel;
					assert(searchLevel < octreeLevels);
                                    }

                                    // We will now interpolate to build the ghost value at the unprocessed face.
                                    UT_Vector3 facePoint, indexNodePoint;
				    myVelocity[level][faceAxis].indexToPos(face[0], face[1], face[2], facePoint);
                                    myNodeLabels[searchLevel].posToIndex(facePoint, indexNodePoint);

                                    fpreal interpolationWeight = indexNodePoint[faceAxis] - floor(indexNodePoint[faceAxis]);

				    fpreal ghostVelocity = 0.;

                                    // Average velocities from adjacent big faces
				    for (int direction : {0,1})
                                    {
					UT_Vector3i offsetFace = HDKcellToFace(cell, faceAxis, direction);

					const exint offsetVelocityIndex = HDKgetFieldValue(myVelocityIndices[searchLevel][faceAxis], offsetFace);

					fpreal localInterpolationWeight = 0;
					if (direction == 0)
					    localInterpolationWeight = 1. - interpolationWeight;
					else
					    localInterpolationWeight = interpolationWeight;

                                        // We average the two big faces together for the dangling edge
                                        if (offsetVelocityIndex >= 0)
					    ghostVelocity += localInterpolationWeight * HDKgetFieldValue(myVelocity[searchLevel][faceAxis], offsetFace);
                                        // If the small inset faces are active, we will want to accept all of them
                                        else if (offsetVelocityIndex == HDK_UNASSIGNED)
                                        {
					    for (int childIndex = 0; childIndex < 4; ++childIndex)
                                            {
                                                // The way the face grids are set up is a little counter
                                                // intuitive since the lower level faces must be inserts
                                                // in the big face.
						UT_Vector3i childFace = myOctreeLabels.getChildFace(offsetFace, faceAxis, childIndex);

						const exint childOffsetVelocityIndex = HDKgetFieldValue(myVelocityIndices[searchLevel - 1][faceAxis], childFace);

                                                if (childOffsetVelocityIndex >= 0)
						    ghostVelocity += .25 * localInterpolationWeight * HDKgetFieldValue(myVelocity[searchLevel - 1][faceAxis], childFace);
                                                else
                                                    assert(false);
                                            }
                                        }
                                    }

				    fpreal searchLevelWeight = fpreal(1 << (octreeLevels - level - 1));

				    fpreal localNodeValue = HDKgetFieldValue(myNodeValues[level][faceAxis], node);
                                    //localNodeValue += weight * ghostVelocity;
				    localNodeValue += searchLevelWeight * ghostVelocity;

				    HDKsetFieldValue(myNodeValues[level][faceAxis], node, localNodeValue);

				    fpreal localNodeWeight = HDKgetFieldValue(nodeWeights[level][faceAxis], node);
                                    localNodeWeight += searchLevelWeight;

				    HDKsetFieldValue(nodeWeights[level][faceAxis], node, localNodeWeight);

                                    localFlag += exint(1 << bitShiftCount);
				}
                            }

                            ++bitShiftCount;
                            tempFlag = tempFlag >> 1;
                        }

                        assert(localFlag == 0xFFF);

			HDKsetFieldValue(nodeFlags[level], node, localFlag);
                    }
                }
            }
        }
    }
}

void
HDK_OctreeVectorFieldInterpolator::normalizeActiveNodesPartial(UT_Array<SIM_RawField> &nodeWeights,
								UT_Array<SIM_RawIndexField> &nodeFlags,
								const int level,
								const UT_JobInfo &info)
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorF vit;
    vit.setConstArray(myNodeLabels[level].field());
    vit.splitByTile(info);

    UT_VoxelTileIteratorF vitt;

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        if (!vit.isTileConstant() || vit.getValue() == ACTIVENODE)
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                if (vitt.getValue() == ACTIVENODE)
                {
                    UT_Vector3i node(vitt.x(), vitt.y(), vitt.z());

		    assert(HDKgetFieldValue(nodeFlags[level], node) == 0xFFF);

                    for (int axis : {0,1,2})
                    {
			fpreal localNodeValue = HDKgetFieldValue(myNodeValues[level][axis], node);
			fpreal localNodeWeight = HDKgetFieldValue(nodeWeights[axis], node);

			assert(localNodeWeight > 0);

			HDKsetFieldValue(myNodeValues[level][axis], node, localNodeValue / localNodeWeight);
                    }
                }
            }
        }
    }
}

void
HDK_OctreeVectorFieldInterpolator::distributeNodeValuesDownPartial(const int level, const UT_JobInfo &info)
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorF vit(myNodeLabels[level].fieldNC());
    vit.setCompressOnExit(true);
    vit.splitByTile(info);

    UT_VoxelTileIteratorF vitt;

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        if (!vit.isTileConstant() || vit.getValue() == DEPENDENTNODE)
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                if (vitt.getValue() == DEPENDENTNODE)
                {
                    UT_Vector3i node(vitt.x(), vitt.y(), vitt.z());
		    UT_Vector3i parentNode = myOctreeLabels.getParentNode(node);

		    if (HDKgetFieldValue(myNodeLabels[level + 1], parentNode) == ACTIVENODE)
                    {
                        for (int axis : {0,1,2})
                        {
			    HDKsetFieldValue(myNodeValues[level][axis], node,
						HDKgetFieldValue(myNodeValues[level + 1][axis], parentNode));
                        }
                    }
                    else
                        assert(false);

                    vitt.setValue(ACTIVENODE);
                }
            }
        }
    }
}

fpreal
HDK_OctreeVectorFieldInterpolator::interpSPGrid(const UT_Vector3 &samplePoint,
                                                const int axis) const
{
    const int octreeLevels = myOctreeLabels.getOctreeLevels();

    UT_Vector3 indexPoint;
    myNodeValues[0][axis].posToIndex(samplePoint, indexPoint);

    UT_Vector3i cell(floor(indexPoint[0]), floor(indexPoint[1]), floor(indexPoint[2]));

    for (int level = 0; level < octreeLevels; ++level)
    {
	if (myOctreeLabels.isCellActive(cell, level))
        {
            // Determine if we can linear interpolate
            bool isAtTransition = false;

	    UT_Vector3 indexFacePoint;
	    myVelocity[level][axis].posToIndex(samplePoint, indexFacePoint);

	    UT_Vector3i face(floor(indexFacePoint[0]), floor(indexFacePoint[1]), floor(indexFacePoint[2]));

	    for (int faceIndex = 0; faceIndex < 8; ++faceIndex)
	    {
		// This can be a little confusing. We have applied a floor to the index position relative to
		// a set of faces. Now we want to check the eight faces that would combine to create an interpolation
		// for this sample point. If any of those faces are unassigned, then we have to drop back to node-based
		// interpolation. We can check these eight faces using the cellToNode mapping.
		UT_Vector3i neighbourFace = HDKcellToNode(face, faceIndex);

		const exint velocityIndex = HDKgetFieldValue(myVelocityIndices[level][axis], neighbourFace);

		if (velocityIndex == HDK_UNASSIGNED)
		{
		    isAtTransition = true;
		    break;
		}
	    }

	    if (!isAtTransition)
	    {
		UT_Vector3 interpolationWeight;

		for (int localAxis : {0,1,2})
		{
		    interpolationWeight[localAxis] = indexFacePoint[localAxis] - fpreal(face[localAxis]);
		    interpolationWeight[localAxis] = SYSclamp(interpolationWeight[localAxis], 0., 1.);
		}

		fpreal interpolatedValue = 0;
		for (int faceIndex = 0; faceIndex < 8; ++faceIndex)
		{
		    UT_Vector3i neighbourFace = HDKcellToNode(face, faceIndex);

		    fpreal weight = 1.;
		    for (int localAxis : {0,1,2})
		    {
			if (neighbourFace[localAxis] - face[localAxis] == 0)
			    weight *= (1. - interpolationWeight[localAxis]);
			else
			    weight *= interpolationWeight[localAxis];
		    }

		    interpolatedValue += weight * HDKgetFieldValue(myVelocity[level][axis], neighbourFace);
		}

		return interpolatedValue;
	    }
            else
            {
                myNodeValues[level][axis].posToIndex(samplePoint, indexPoint);

                // The interpolation weight within cell (along the face axis)
                fpreal cellInterpolationWeight = indexPoint[axis] - fpreal(cell[axis]);

                cellInterpolationWeight = SYSclamp(cellInterpolationWeight, 0., 1.);

                fpreal faceInterpolationValue[2] = {0., 0.};

                int faceInterpolationAxis[2];
                faceInterpolationAxis[0] = (axis + 1) % 3;
                faceInterpolationAxis[1] = (axis + 2) % 3;

		for (int direction : {0,1})
                {
                    // Find the face to interpolate within on either end of the cell
		    UT_Vector3i adjacentFace = HDKcellToFace(cell, axis, direction);

		    const exint velocityIndex = HDKgetFieldValue(myVelocityIndices[level][axis], adjacentFace);

                    int faceLevel = level;

                    // If the big face is inactive, find a small face to project onto
                    if (velocityIndex == HDK_UNASSIGNED)
                    {
			assert(level >= 0);
                        UT_Vector3 childIndexPoint;
			myNodeValues[level - 1][axis].posToIndex(samplePoint, childIndexPoint);

                        UT_Vector3i childFace;
#if !defined(NDEBUG)
                        bool foundValidChildFace = false;
#endif
			for (int childIndex = 0; childIndex < 4; ++childIndex)
                        {
                            // The way the face grids are set up is a little
                            // counter intuitive since the lower level faces
                            // must be inserts in the big face.
			    childFace = myOctreeLabels.getChildFace(adjacentFace, axis, childIndex);

			    assert(HDKgetFieldValue(myVelocityIndices[level - 1][axis], childFace) >= 0);

                            if (fpreal(childFace[faceInterpolationAxis[0]]) <= childIndexPoint[faceInterpolationAxis[0]] &&
                                fpreal(childFace[faceInterpolationAxis[1]]) <= childIndexPoint[faceInterpolationAxis[1]] &&
                                fpreal(childFace[faceInterpolationAxis[0]] + 1) >= childIndexPoint[faceInterpolationAxis[0]] &&
                                fpreal(childFace[faceInterpolationAxis[1]] + 1) >= childIndexPoint[faceInterpolationAxis[1]])
                            {
                                faceLevel = level - 1;
                                adjacentFace = childFace;

#if !defined(NDEBUG)
                                foundValidChildFace = true;
#endif

                                break;
                            }
                        }

                        assert(foundValidChildFace);
                    }

		    assert(HDKgetFieldValue(myVelocityIndices[faceLevel][axis], adjacentFace) >= 0);

                    fpreal averageNodeVelocity = 0;
                    UT_Vector3 indexNodePoint;

		    myNodeValues[faceLevel][axis].posToIndex(samplePoint, indexNodePoint);

                    fpreal faceInterpolationWeights[2];
                    faceInterpolationWeights[0] = indexNodePoint[faceInterpolationAxis[0]] - floor(indexNodePoint[faceInterpolationAxis[0]]);
                    faceInterpolationWeights[1] = indexNodePoint[faceInterpolationAxis[1]] - floor(indexNodePoint[faceInterpolationAxis[1]]);

                    assert(faceInterpolationWeights[0] >= 0. && faceInterpolationWeights[0] <= 1. &&
			    faceInterpolationWeights[1] >= 0. && faceInterpolationWeights[1] <= 1.);

		    fpreal faceVelocity = HDKgetFieldValue(myVelocity[faceLevel][axis], adjacentFace);

		    for (int nodeIndex = 0; nodeIndex < 4; ++nodeIndex)
                    {
			UT_Vector3i node = HDKfaceToNode(adjacentFace, axis, nodeIndex);

			fpreal weight = 1.;

                        for (int localAxisIndex : {0,1})
                        {
                            int localAxis = faceInterpolationAxis[localAxisIndex];

                            if (node[localAxis] - adjacentFace[localAxis] == 0)
                                weight *= 1. - faceInterpolationWeights[localAxisIndex];
                            else
                                weight *= faceInterpolationWeights[localAxisIndex];
                        }

			fpreal nodeValue = HDKgetFieldValue(myNodeValues[faceLevel][axis], node);

                        averageNodeVelocity += nodeValue;
			faceInterpolationValue[direction] += nodeValue * weight;
                    }

                    faceInterpolationValue[direction] += 2. * (faceVelocity - .25 * averageNodeVelocity) *
							SYSmin(faceInterpolationWeights[0],
							    SYSmin(faceInterpolationWeights[1],
								SYSmin(1. - faceInterpolationWeights[0], 1 - faceInterpolationWeights[1])));
                }

                return (1. - cellInterpolationWeight) * faceInterpolationValue[0] + cellInterpolationWeight * faceInterpolationValue[1];
            }
        }

	cell = myOctreeLabels.getParentCell(cell);
    }

    assert(false);
    return 0.;
}