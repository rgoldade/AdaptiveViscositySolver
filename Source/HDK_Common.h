#ifndef HDK_COMMON_H
#define HDK_COMMON_H

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

class UT_Vector5I : public UT_FixedVector<int, 5, true>
{
public:
    SYS_FORCE_INLINE UT_Vector5I() : UT_Vector5I(0, 0, 0, 0, 0) {}

    SYS_FORCE_INLINE UT_Vector5I(int v0, int v1, int v2, int v3, int v4)
    {
        vec[0] = v0;
        vec[1] = v1;
        vec[2] = v2;
        vec[3] = v3;
        vec[4] = v4;
    }
};

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


SYS_FORCE_INLINE UT_Vector3i HDKcellToFace(UT_Vector3i cell, const int axis, const int direction)
{
    assert(axis >= 0 && axis < 3);
    if (direction == 1)
	++cell[axis]; 
    else assert(direction == 0);

    return cell;
}

SYS_FORCE_INLINE UT_Vector3i HDKcellToCell(UT_Vector3i cell, const int axis, const int direction)
{
    assert(axis >= 0 && axis < 3);
    if (direction == 0)
	--cell[axis];
    else
    {
	++cell[axis];
	assert(direction == 1);
    }

    return cell;
}

SYS_FORCE_INLINE UT_Vector3i HDKcellToEdge(UT_Vector3i cell, const int edgeAxis, const int edgeIndex)
{
    assert(edgeAxis >= 0 && edgeAxis < 3);
    assert(edgeIndex >= 0 && edgeIndex < 4);

    for (int axisOffset : {0,1})
    {
	if (edgeIndex & (1 << axisOffset))
	{
	    int localAxis = (edgeAxis + 1 + axisOffset) % 3;
	    ++cell[localAxis];
	}
    }

    return cell;
}

SYS_FORCE_INLINE UT_Vector3i HDKcellToNode(UT_Vector3i cell, const int nodeIndex)
{
    assert(nodeIndex >= 0 && nodeIndex < 8);

    for (int axis : {0,1,2})
    {
	if (nodeIndex & (1 << axis))
	    ++cell[axis];
    }

    return cell;
}

SYS_FORCE_INLINE UT_Vector3i HDKfaceToCell(UT_Vector3i face, const int axis, const int direction)
{
    assert(axis >= 0 && axis < 3);
    if (direction == 0)
	--face[axis];
    else
	assert(direction == 1);

    return face;
}

SYS_FORCE_INLINE UT_Vector3i HDKfaceToEdge(UT_Vector3i face, const int faceAxis,
					    const int edgeAxis, const int direction)
{   
    assert(faceAxis >= 0 && faceAxis < 3 && edgeAxis >= 0 && edgeAxis < 3);
    assert(faceAxis != edgeAxis);

    if (direction == 1)
    {
	int offsetAxis = 3 - faceAxis - edgeAxis;
	++face[offsetAxis];
    }
    else
	assert(direction == 0);

    return face;
}

SYS_FORCE_INLINE UT_Vector3i HDKfaceToNode(UT_Vector3i face, const int faceAxis, const int nodeIndex)
{
    assert(faceAxis >= 0 && faceAxis < 3);
    assert(nodeIndex >= 0 && nodeIndex < 4);

    for (int axisOffset : {0,1})
    {
	if (nodeIndex & (1 << axisOffset))
	{
	    int localAxis = (faceAxis + 1 + axisOffset) % 3;
	    ++face[localAxis];
	}
    }

    return face;
}

SYS_FORCE_INLINE UT_Vector3i HDKedgeToFace(UT_Vector3i edge, const int edgeAxis,
					    const int faceAxis, const int direction)
{
    assert(faceAxis >= 0 && faceAxis < 3 && edgeAxis >= 0 && edgeAxis < 3);
    assert(faceAxis != edgeAxis);

    if (direction == 0)
    {
	int offsetAxis = 3 - faceAxis - edgeAxis;
	--edge[offsetAxis];
    }
    else
	assert(direction == 1);

    return edge;
}

SYS_FORCE_INLINE UT_Vector3i HDKedgeToCell(UT_Vector3i edge, const int edgeAxis, const int cellIndex)
{
    assert(edgeAxis >= 0 && edgeAxis < 3);
    assert(cellIndex >= 0 && cellIndex < 4);

    for (int axisOffset : {0,1})
    {
	if (!(cellIndex & (1 << axisOffset)))
	{
	    int localAxis = (edgeAxis + 1 + axisOffset) % 3;
	    --edge[localAxis];
	}
    }

    return edge;
}

SYS_FORCE_INLINE UT_Vector3i HDKnodeToFace(UT_Vector3i node, const int faceAxis, const int faceIndex)
{
    assert(faceAxis >= 0 && faceAxis < 3);
    assert(faceIndex >= 0 && faceIndex < 4);

    for (int axisOffset : {0,1})
    {
	if (!(faceIndex & (1 << axisOffset)))
	{
	    int localAxis = (faceAxis + 1 + axisOffset) % 3;
	    --node[localAxis];
	}
    }

    return node;
}

SYS_FORCE_INLINE UT_Vector3i HDKnodeToCell(UT_Vector3i node, const int cellIndex)
{
    assert(cellIndex >= 0 && cellIndex < 8);

    for (int axis : {0,1,2})
    {
	if (!(cellIndex & (1 << axis)))
	    --node[axis];
    }

    return node;
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