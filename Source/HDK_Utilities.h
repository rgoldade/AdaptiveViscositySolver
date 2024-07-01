#ifndef HDK_UTILITIES_H
#define HDK_UTILITIES_H

#ifdef USEEIGEN
    #include "Eigen/Sparse"
#endif

#include <SIM/SIM_RawField.h>
#include <SIM/SIM_RawIndexField.h>
#include <UT/UT_SparseMatrix.h>
#include <UT/UT_Vector3.h>

#define HDK_XEDGE SIM_SAMPLE_EDGEYZ
#define HDK_YEDGE SIM_SAMPLE_EDGEXZ
#define HDK_ZEDGE SIM_SAMPLE_EDGEXY
#define HDK_NODE SIM_SAMPLE_CORNER

constexpr exint HDK_FLUID = 0;
constexpr exint HDK_UNASSIGNED = -1;
constexpr exint HDK_SOLIDBOUNDARY = -2;
constexpr exint HDK_OUTSIDE = -3;

using UT_Vector5I = UT_FixedVector<int, 5>;

#ifdef USESINGLEPRECISION
    using SolveType = fpreal32;
    
    #ifdef USEEIGEN
	using Vector = Eigen::VectorXf;
    #endif
#else
    using SolveType = fpreal64;

    #ifdef USEEIGEN
	using Vector = Eigen::VectorXd;
    #endif
#endif

#ifndef USEEIGEN
    using SparseMatrix = UT_SparseMatrixT<SolveType, false>;
    using SparseMatrixELLT = UT_SparseMatrixELLT<SolveType, false>;
    using Vector = UT_VectorT<SolveType>;
#endif


SYS_FORCE_INLINE UT_Vector3i HDKcellToFace(const UT_Vector3i &cell, const int axis, const int direction)
{
    UT_Vector3i face(cell);
    if (direction == 1)
	++face[axis]; 
    else assert(direction == 0);

    return face;
}

SYS_FORCE_INLINE UT_Vector3i HDKcellToCell(const UT_Vector3i &cell, const int axis, const int direction)
{
    UT_Vector3i adjacentCell(cell);
    if (direction == 0)
	--adjacentCell[axis];
    else
    {
	++adjacentCell[axis];
	assert(direction == 1);
    }

    return adjacentCell;
}

SYS_FORCE_INLINE UT_Vector3i HDKcellToEdge(const UT_Vector3i &cell, const int edgeAxis, const int edgeIndex)
{
    assert(edgeAxis >= 0 && edgeAxis < 3);
    assert(edgeIndex >= 0 && edgeIndex < 4);

    UT_Vector3i edge(cell);
    for (int axisOffset : {0,1})
    {
	if (edgeIndex & (1 << axisOffset))
	{
	    int localAxis = (edgeAxis + 1 + axisOffset) % 3;
	    ++edge[localAxis];
	}
    }

    return edge;
}

SYS_FORCE_INLINE UT_Vector3i HDKcellToNode(const UT_Vector3i &cell, const int nodeIndex)
{
    assert(nodeIndex >= 0 && nodeIndex < 8);

    UT_Vector3i node(cell);
    for (int axis : {0,1,2})
    {
	if (nodeIndex & (1 << axis))
	    ++node[axis];
    }

    return node;
}

SYS_FORCE_INLINE UT_Vector3i HDKfaceToCell(const UT_Vector3i &face, const int axis, const int direction)
{
    assert(axis >= 0 && axis < 3);

    UT_Vector3i cell(face);
    if (direction == 0)
	--cell[axis];
    else
	assert(direction == 1);

    return cell;
}

SYS_FORCE_INLINE UT_Vector3i HDKfaceToEdge(const UT_Vector3i &face, const int faceAxis,
					    const int edgeAxis, const int direction)
{   
    assert(faceAxis >= 0 && faceAxis < 3 && edgeAxis >= 0 && edgeAxis < 3);
    assert(faceAxis != edgeAxis);

    UT_Vector3i edge(face);
    if (direction == 1)
    {
	int offsetAxis = 3 - faceAxis - edgeAxis;
	++edge[offsetAxis];
    }
    else
	assert(direction == 0);

    return edge;
}

SYS_FORCE_INLINE UT_Vector3i HDKfaceToNode(const UT_Vector3i &face, const int faceAxis, const int nodeIndex)
{
    assert(faceAxis >= 0 && faceAxis < 3);
    assert(nodeIndex >= 0 && nodeIndex < 4);

    UT_Vector3i node(face);
    for (int axisOffset : {0,1})
    {
	if (nodeIndex & (1 << axisOffset))
	{
	    int localAxis = (faceAxis + 1 + axisOffset) % 3;
	    ++node[localAxis];
	}
    }

    return node;
}

SYS_FORCE_INLINE UT_Vector3i HDKedgeToFace(const UT_Vector3i &edge, const int edgeAxis,
					    const int faceAxis, const int direction)
{
    assert(faceAxis >= 0 && faceAxis < 3 && edgeAxis >= 0 && edgeAxis < 3);
    assert(faceAxis != edgeAxis);

    UT_Vector3i face(edge);
    if (direction == 0)
    {
	int offsetAxis = 3 - faceAxis - edgeAxis;
	--face[offsetAxis];
    }
    else
	assert(direction == 1);

    return face;
}

SYS_FORCE_INLINE UT_Vector3i HDKedgeToCell(const UT_Vector3i &edge, const int edgeAxis, const int cellIndex)
{
    assert(edgeAxis >= 0 && edgeAxis < 3);
    assert(cellIndex >= 0 && cellIndex < 4);

    UT_Vector3i cell(edge);
    for (int axisOffset : {0,1})
    {
	if (!(cellIndex & (1 << axisOffset)))
	{
	    int localAxis = (edgeAxis + 1 + axisOffset) % 3;
	    --cell[localAxis];
	}
    }

    return cell;
}

SYS_FORCE_INLINE UT_Vector3i HDKnodeToFace(const UT_Vector3i &node, const int faceAxis, const int faceIndex)
{
    assert(faceAxis >= 0 && faceAxis < 3);
    assert(faceIndex >= 0 && faceIndex < 4);

    UT_Vector3i face(node);
    for (int axisOffset : {0,1})
    {
	if (!(faceIndex & (1 << axisOffset)))
	{
	    int localAxis = (faceAxis + 1 + axisOffset) % 3;
	    --face[localAxis];
	}
    }

    return face;
}

SYS_FORCE_INLINE UT_Vector3i HDKnodeToCell(const UT_Vector3i &node, const int cellIndex)
{
    assert(cellIndex >= 0 && cellIndex < 8);

    UT_Vector3i cell(node);
    for (int axis : {0,1,2})
    {
	if (!(cellIndex & (1 << axis)))
	    --cell[axis];
    }

    return cell;
}

SYS_FORCE_INLINE fpreal32 HDKgetFieldValue(const SIM_RawField &field, const UT_Vector3i &cell)
{
    return (*field.field())(cell[0], cell[1], cell[2]);
}

SYS_FORCE_INLINE exint HDKgetFieldValue(const SIM_RawIndexField &field, const UT_Vector3i &cell)
{
    return (*field.field())(cell[0], cell[1], cell[2]);
}

SYS_FORCE_INLINE void HDKsetFieldValue(SIM_RawField &field, const UT_Vector3i &cell, const fpreal32 value)
{
    field.fieldNC()->setValue(cell[0], cell[1], cell[2], value);
}

SYS_FORCE_INLINE void HDKsetFieldValue(SIM_RawIndexField &field, const UT_Vector3i &cell, const exint value)
{
    field.fieldNC()->setValue(cell[0], cell[1], cell[2], value);
}

#endif