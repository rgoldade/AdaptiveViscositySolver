#ifndef HDK_ADAPTIVEVISCOSITY_H
#define HDK_ADAPTIVEVISCOSITY_H

#include <GAS/GAS_SubSolver.h>
#include <GAS/GAS_Utils.h>

#include "HDK_Utilities.h"
#include "HDK_OctreeGrid.h"
#include "HDK_OctreeVectorFieldInterpolator.h"

class HDK_OctreeVectorFieldInterpolator;

class GAS_API HDK_AdaptiveViscosity : public GAS_SubSolver
{
public:

    struct StressStencilFace
    {
        StressStencilFace(exint faceIndex, fpreal coefficient)
	    : myFaceIndex(faceIndex),
	    myCoefficient(coefficient)
	    {}

        exint myFaceIndex;
        fpreal myCoefficient;
    };

    GET_DATA_FUNC_F(SIM_NAME_TOLERANCE, SolverTolerance);
    GET_DATA_FUNC_I("maxIterations", MaxIterations);
    GET_DATA_FUNC_I("numberSuperSamples", NumberSuperSamples);

    GET_DATA_FUNC_I("octreeLevels", OctreeLevels);
    GET_DATA_FUNC_I("fineBandwidth", FineBandwidth);

    GET_DATA_FUNC_B("useEnhancedGradients", UseEnhancedGradients);

    GET_DATA_FUNC_B("doApplySolidWeights", DoApplySolidWeights);
    GET_DATA_FUNC_B("doPrintOctree", DoPrintOctree);
    GET_DATA_FUNC_B("onlyPrintOctree", OnlyPrintOctree);

    GET_DATA_FUNC_F("extrapolation", Extrapolation);

protected:
    explicit HDK_AdaptiveViscosity(const SIM_DataFactory *factory);
    virtual ~HDK_AdaptiveViscosity();

    // Used to determine if the field is complicated enough to justify
    // the overhead of multithreading.
    bool shouldMultiThread(const SIM_RawField *field) const
    {
        return field->field()->numTiles() > 1;
    }

    // The overloaded callback that GAS_SubSolver will invoke to
    // perform our actual computation.  We are giving a single object
    // at a time to work on.
    virtual bool solveGasSubclass(SIM_Engine &engine, SIM_Object *obj,
                                  SIM_Time time, SIM_Time timestep);

private:
    // We define this to be a DOP_Auto node which means we do not
    // need to implement a DOP_Node derivative for this data.  Instead,
    // this description is used to define the interface.
    static const SIM_DopDescription *getDopDescription();
    /// These macros are necessary to bind our node to the factory and
    /// ensure useful constants like BaseClass are defined.
    DECLARE_STANDARD_GETCASTTOTYPE();
    DECLARE_DATAFACTORY(HDK_AdaptiveViscosity, GAS_SubSolver,
                        "HDK Adaptive Viscosity", getDopDescription());

    void buildIntegrationWeights(SIM_RawField &centerIntegrationWeights,
				    SIM_RawField (&edgeIntegrationWeights)[3],
				    const SIM_RawField &liquidSurface,
				    const SIM_RawField &solidSurface,
				    const fpreal extrapolation,
				    const bool doApplySolidWeights) const;

    void buildOctree(HDK_OctreeGrid &octree,
			const SIM_RawField &liquidSurface,
			const SIM_RawField &solidSurface,
			const int octreeLevels,
			const fpreal extrapolation,
			const fpreal outerFineCellBandwidth,
			const fpreal innerFineCellBandwidth) const;

    //
    // Uncompress velocity and edge grids before labelling
    //

    THREADED_METHOD4_CONST(HDK_AdaptiveViscosity,
			    octreeGridLabels.shouldMultiThread(),
			    findOccupiedOctreeVelocityTiles,
			    UT_Array<bool> &, isisTileOccupiedListList,
			    const SIM_RawField &, octreeGridLabels,
			    const SIM_RawIndexField &, octreeVelocityIndices,
			    const int, axis)

    void findOccupiedOctreeVelocityTilesPartial(UT_Array<bool> &isisTileOccupiedListList,
                                                const SIM_RawField &octreeGridLabels,
                                                const SIM_RawIndexField &octreeVelocityIndices,
                                                const int axis,
                                                const UT_JobInfo &info) const;

    THREADED_METHOD4_CONST(HDK_AdaptiveViscosity,
			    liquidSurface.shouldMultiThread(),
			    findOccupiedRegularVelocityTiles,
			    UT_Array<bool> &, isisTileOccupiedListList,
			    const SIM_RawField &, liquidSurface,
			    const SIM_RawIndexField &, regularVelocityIndices,
			    const int, axis)

    void findOccupiedRegularVelocityTilesPartial(UT_Array<bool> &isTileOccupiedList,
                                                 const SIM_RawField &liquidSurface,
                                                 const SIM_RawIndexField &regularVelocityIndices,
                                                 const int axis,
                                                 const UT_JobInfo &info) const;

    THREADED_METHOD4_CONST(HDK_AdaptiveViscosity,
			    octreeGridLabels.shouldMultiThread(),
			    findOccupiedEdgeStressTiles,
			    UT_Array<bool> &, isTileOccupiedList,
			    const SIM_RawField &, octreeGridLabels,
			    const SIM_RawIndexField &, edgeStressIndices,
			    const int, axis)

    void findOccupiedEdgeStressTilesPartial(UT_Array<bool> &isTileOccupiedList,
					    const SIM_RawField &octreeGridLabels,
					    const SIM_RawIndexField &edgeStressIndices,
					    const int axis,
					    const UT_JobInfo &info) const;

    THREADED_METHOD2_CONST(HDK_AdaptiveViscosity,
			    isTileOccupiedList.entries() > 20,
			    uncompressTiles,
			    SIM_RawIndexField &, grid,
			    const UT_Array<bool> &, isTileOccupiedList)

    void uncompressTilesPartial(SIM_RawIndexField &grid,
				const UT_Array<bool> &isTileOccupiedList,
				const UT_JobInfo &info) const;

    //
    // Build labels for regular and octree velocities and octree edges
    //

    THREADED_METHOD7_CONST(HDK_AdaptiveViscosity,
                           regularVelocityIndices.shouldMultiThread(),
                           classifyRegularVelocityFaces,
			   SIM_RawIndexField &, regularVelocityIndices,
			   const SIM_RawField &, liquidSurface,
                           const SIM_RawField &, solidSurface,
                           const SIM_RawField &, centerIntegrationWeights,
                           const SIM_RawField *, edgeIntegrationWeights,
			   const int, axis,
                           const fpreal, extrapolation);

    void classifyRegularVelocityFacesPartial(SIM_RawIndexField &regularVelocityIndices,
						const SIM_RawField &liquidSurface,
						const SIM_RawField &solidSurface,
						const SIM_RawField &centerIntegrationWeights,
						const SIM_RawField *edgeIntegrationWeights,
						const int axis,
						const fpreal extrapolation, const UT_JobInfo &info) const;

    exint buildRegularVelocityIndices(SIM_RawIndexField (&regularVelocityIndices)[3],
					const SIM_RawField &liquidSurface,
					const SIM_RawField &solidSurface,
					const SIM_RawField &centerIntegrationWeights,
					const SIM_RawField (&edgeIntegrationWeights)[3],
					const fpreal extrapolation) const;

    THREADED_METHOD9_CONST(HDK_AdaptiveViscosity,
			    octreeVelocityIndices.shouldMultiThread(),
			    classifyOctreeVelocityFaces,
			    SIM_RawIndexField &, octreeVelocityIndices,
			    const HDK_OctreeGrid &, octree,
			    const SIM_RawField &, liquidSurface,
			    const SIM_RawField &, solidSurface,
			    const SIM_RawField &, centerIntegrationWeights,
			    const SIM_RawField *, edgeIntegrationWeights,
			    const int, axis, const int, level,
			    const fpreal, extrapolation);

    void classifyOctreeVelocityFacesPartial(SIM_RawIndexField &octreeVelocityIndices,
					    const HDK_OctreeGrid &octree,
					    const SIM_RawField &liquidSurface,
					    const SIM_RawField &solidSurface,
					    const SIM_RawField &centerIntegrationWeights,
					    const SIM_RawField *edgeIntegrationWeights,
					    const int axis, const int level,
					    const fpreal extrapolation,
					    const UT_JobInfo &info) const;

    exint buildOctreeVelocityIndices(UT_Array<UT_Array<SIM_RawIndexField>> &octreeVelocityIndices,
					const SIM_RawField &liquidSurface,
					const SIM_RawField &solidSurface,
					const HDK_OctreeGrid &octreeLabels,
					const SIM_RawField &centerIntegrationWeights,
					const SIM_RawField (&edgeIntegrationWeights)[3],
					const fpreal extrapolation) const;

    THREADED_METHOD6_CONST(HDK_AdaptiveViscosity,
			    edgeStressIndices.shouldMultiThread(),
			    classifyEdgeStresses,
			    SIM_RawIndexField &, edgeStressIndices,
			    const HDK_OctreeGrid &, octreeLabels,
			    const SIM_RawField &, liquidSurface,
			    const SIM_RawField &, edgeIntegrationWeights,
			    const int, axis, const int, level);

    void classifyEdgeStressesPartial(SIM_RawIndexField &edgeStressIndices,
                                     const HDK_OctreeGrid &octree,
                                     const SIM_RawField &liquidSurface,
                                     const SIM_RawField &edgeIntegrationWeights,
                                     const int axis, const int level,
                                     const UT_JobInfo &info) const;

    exint buildEdgeStressIndices(UT_Array<UT_Array<SIM_RawIndexField>> &edgeStressIndices,
				    const SIM_RawField &liquidSurface,
				    const HDK_OctreeGrid &octree,
				    const SIM_RawField (&edgeIntegrationWeights)[3]) const;

    THREADED_METHOD4_CONST(HDK_AdaptiveViscosity,
			    octreeGridLabels.shouldMultiThread(),
			    classifyCenterStresses,
			    SIM_RawIndexField &, centerStressIndices,
			    const SIM_RawField &, octreeGridLabels,
			    const SIM_RawField &, centerIntegrationWeights,
			    const int, level)

    void classifyCenterStressesPartial(SIM_RawIndexField &centerStressIndices,
					const SIM_RawField &octreeGridLabels,
					const SIM_RawField &centerIntegrationWeights,
					const int level,
					const UT_JobInfo &info) const;

    exint buildCenterStressIndices(UT_Array<SIM_RawIndexField> &centerStressIndices,
				    const HDK_OctreeGrid &octree,
				    const SIM_RawField &centerIntegrationWeights) const;

    struct edgeStressStencilParms
    {
        edgeStressStencilParms(const UT_Array<UT_Array<SIM_RawIndexField>> &octreeVelocityIndices,
				const UT_Array<UT_Array<SIM_RawIndexField>> &edgeStressIndices,
				const HDK_OctreeGrid &octreeLabels,
				const SIM_VectorField &solidVelocity,
				const SIM_RawField (&edgeIntegrationWeights)[3],
				const SIM_RawField &viscosity,
				const fpreal dt,
				const bool useEnhancedGradients = false)
            : myOctreeVelocityIndices(octreeVelocityIndices),
              myEdgeStressIndices(edgeStressIndices),
              myOctreeLabels(octreeLabels),
              mySolidVelocity(solidVelocity),
              myEdgeIntegrationWeights(edgeIntegrationWeights),
              myViscosity(viscosity),
              myDt(dt),
              myUseEnhancedGradients(useEnhancedGradients)
        {
        }

        const UT_Array<UT_Array<SIM_RawIndexField>> &myOctreeVelocityIndices;
        const UT_Array<UT_Array<SIM_RawIndexField>> &myEdgeStressIndices;
        const HDK_OctreeGrid &myOctreeLabels;
        const SIM_VectorField &mySolidVelocity;
        const SIM_RawField (&myEdgeIntegrationWeights)[3];
        const SIM_RawField &myViscosity;
        const bool myUseEnhancedGradients;
        const fpreal myDt;
    };

    THREADED_METHOD6_CONST(HDK_AdaptiveViscosity,
                           parms.myEdgeStressIndices[level][axis].shouldMultiThread(),
                           buildEdgeStressStencils,
                           UT_Array<UT_Array<StressStencilFace>> &, edgeStressStencils,
			   UT_Array<UT_Array<fpreal>> &, edgeStressBoundaryStencils,
			   UT_Array<fpreal> &, edgeStressStencilWeights,
			   const edgeStressStencilParms &, parms,
                           const int, axis, const int, level)

    void buildEdgeStressStencilsPartial(UT_Array<UT_Array<StressStencilFace>> &edgeStressStencils,
					UT_Array<UT_Array<fpreal>> &edgeStressBoundaryStencils,
					UT_Array<fpreal> &edgeStressStencilWeights,
					const edgeStressStencilParms &parms,
					const int axis, const int level,
					const UT_JobInfo &info) const;

    struct centerStressStencilParms
    {
        centerStressStencilParms(const UT_Array<UT_Array<SIM_RawIndexField>> &octreeVelocityIndices,
				    const UT_Array<SIM_RawIndexField> &centerStressIndices,
				    const HDK_OctreeGrid &octreeLabels,
				    const SIM_VectorField &solidVelocity)
            : myOctreeVelocityIndices(octreeVelocityIndices),
              myCenterStressIndices(centerStressIndices),
              myOctreeLabels(octreeLabels),
              mySolidVelocity(solidVelocity)
        {
        }

        const UT_Array<UT_Array<SIM_RawIndexField>> &myOctreeVelocityIndices;
        const UT_Array<SIM_RawIndexField> &myCenterStressIndices;
        const HDK_OctreeGrid &myOctreeLabels;
	const SIM_VectorField &mySolidVelocity;
    };

    THREADED_METHOD6_CONST(HDK_AdaptiveViscosity,
			    parms.myCenterStressIndices[level].shouldMultiThread(),
			    buildCenterStressStencils,
			    UT_Array<UT_Array<StressStencilFace>> &, centerStressStencils,
			    UT_Array<UT_Array<fpreal>> &, centerStressBoundaryStencils,
			    const centerStressStencilParms &, parms,
			    const int, axis, const int, level,
			    const int, cellcount)

    void buildCenterStressStencilsPartial(UT_Array<UT_Array<StressStencilFace>> &centerStressStencils,
					    UT_Array<UT_Array<fpreal>> &centerStressBoundaryStencils,
					    const centerStressStencilParms &parms,
					    const int axis, const int level,
					    const exint cellcount,
					    const UT_JobInfo &info) const;

    THREADED_METHOD6_CONST(HDK_AdaptiveViscosity,
                           centerStressIndices.shouldMultiThread(),
			   buildCenterStressWeights,
                           UT_Array<fpreal> &, centerStressStencilWeights,
                           const SIM_RawIndexField &, centerStressIndices,
                           const SIM_RawField &, centerIntegrationWeights,
                           const SIM_RawField &, viscosity,
			   const fpreal, dt,
                           const int, level)

    void buildCenterStressWeightsPartial(UT_Array<fpreal> &centerStressStencilWeights,
					    const SIM_RawIndexField &centerStressIndices,
					    const SIM_RawField &centerIntegrationWeights,
					    const SIM_RawField &viscosity,
					    const fpreal dt,
					    const int level,
					    const UT_JobInfo &info) const;

    THREADED_METHOD8_CONST(HDK_AdaptiveViscosity,
			    octreeVelocityIndices.shouldMultiThread(),
			    buildVelocityMapping,
			    Vector &, initialGuess,
			    const SIM_RawField &, regularVelocity,
			    const SIM_RawIndexField &, regularVelocityIndices,
			    const SIM_RawIndexField &, octreeVelocityIndices,
			    const SIM_VectorField &, solidVelocity,
			    const HDK_OctreeGrid &, octreeLabels,
			    const int, axis, const int, level);

    void buildVelocityMappingPartial(Vector &initialGuess,
					const SIM_RawField &regularVelocity,
					const SIM_RawIndexField &regularVelocityIndices,
					const SIM_RawIndexField &octreeVelocityIndices,
					const SIM_VectorField &solidVelocity,
					const HDK_OctreeGrid &octreeLabels,
					const int axis, const int level,
					const UT_JobInfo &info) const;

    struct octreeSystemStencilParms
    {
        octreeSystemStencilParms(const UT_Array<UT_Array<StressStencilFace>> &edgeStressStencils,
				    const UT_Array<UT_Array<fpreal>> &edgeStressBoundaryStencils,
				    const UT_Array<fpreal> &edgeStressStencilWeights,

				    const UT_Array<UT_Array<StressStencilFace>> &centerStressStencils,
				    const UT_Array<UT_Array<fpreal>> &centerStressBoundaryStencils,
				    const UT_Array<fpreal> &centerStressStencilWeights,

				    const UT_Array<UT_Array<SIM_RawIndexField>> &octreeVelocityIndices,
				    const UT_Array<UT_Array<SIM_RawIndexField>> &edgeStressIndices,
				    const UT_Array<SIM_RawIndexField> &centerStressIndices,

				    const HDK_OctreeGrid &octreeLabels,

				    const SIM_VectorField &faceIntegrationWeights,
				    const SIM_RawField &centerIntegrationWeights,
				    const SIM_RawField &density,
				    
				    const Vector &initialGuess,
				    const bool useEnhancedGradients = false)

		: myEdgeStressStencils(edgeStressStencils),
		myEdgeBoundaryStressStencils(edgeStressBoundaryStencils),
		myEdgeStressStencilWeights(edgeStressStencilWeights),

		myCenterStressStencils(centerStressStencils),
		myCenterBoundaryStressstencils(centerStressBoundaryStencils),
		myCenterStressStencilWeights(centerStressStencilWeights),

		myOctreeVelocityIndices(octreeVelocityIndices),
		myEdgeStressIndices(edgeStressIndices),
		myCenterStressIndices(centerStressIndices),

		myOctreeLabels(octreeLabels),

		myFaceIntegrationWeights(faceIntegrationWeights),
		myCenterIntegrationWeights(centerIntegrationWeights),
		myDensity(density),
		myInitialGuess(initialGuess),

		myUseEnhancedGradients(useEnhancedGradients)
        {
        }

        const UT_Array<UT_Array<StressStencilFace>> &myEdgeStressStencils;
        const UT_Array<UT_Array<fpreal>> &myEdgeBoundaryStressStencils;
        const UT_Array<fpreal> &myEdgeStressStencilWeights;

        const UT_Array<UT_Array<StressStencilFace>> &myCenterStressStencils;
        const UT_Array<UT_Array<fpreal>> &myCenterBoundaryStressstencils;
        const UT_Array<fpreal> &myCenterStressStencilWeights;

        const UT_Array<UT_Array<SIM_RawIndexField>> &myOctreeVelocityIndices;
        const UT_Array<UT_Array<SIM_RawIndexField>> &myEdgeStressIndices;
        const UT_Array<SIM_RawIndexField> &myCenterStressIndices;

        const HDK_OctreeGrid &myOctreeLabels;

        const SIM_VectorField &myFaceIntegrationWeights;
        const SIM_RawField &myCenterIntegrationWeights;
        const SIM_RawField &myDensity;
        const Vector &myInitialGuess;

        const bool myUseEnhancedGradients;
    };

    THREADED_METHOD6_CONST(HDK_AdaptiveViscosity,
			    parms.myOctreeVelocityIndices[level][axis].shouldMultiThread(),
			    buildOctreeSystemFromStencils,
#ifdef USEEIGEN
			    std::vector<std::vector<Eigen::Triplet<SolveType>>> &, parallelSparseMatrixElements,
#else
			    UT_Array<SparseMatrix> &, parallelSparseMatrixElements,
#endif
			    Vector &, rhs,
			    const octreeSystemStencilParms &, parms,
			    const int, axis, const int, level,
			    const exint, cellcount);

    void buildOctreeSystemFromStencilsPartial(
#ifdef USEEIGEN
						std::vector<std::vector<Eigen::Triplet<SolveType>>> &parallelSparseMatrixElements,
#else
						UT_Array<SparseMatrix> &parallelSparseMatrixElements,
#endif
						Vector &rhs,
						const octreeSystemStencilParms &parms,
						const int axis, const int level,
						const exint cellcount,
						const UT_JobInfo &info) const;

    THREADED_METHOD3_CONST(HDK_AdaptiveViscosity,
			    octreeVelocityIndices.shouldMultiThread(),
			    setOctreeVelocity,
			    SIM_RawField &, octreeVelocity,
			    const SIM_RawIndexField &, octreeVelocityIndices,
			    const Vector &, solution)

    void setOctreeVelocityPartial(SIM_RawField &octreeVelocity,
				    const SIM_RawIndexField &octreeVelocityIndices,
				    const Vector &solution,
				    const UT_JobInfo &info) const;

    THREADED_METHOD8_CONST(HDK_AdaptiveViscosity,
			    regularVelocityIndices.shouldMultiThread(),
			    applyVelocitiesToRegularGrid,
			    SIM_RawField &, regularVelocity,
			    const HDK_OctreeVectorFieldInterpolator &, interpolator,
			    const SIM_RawIndexField &, octreeVelocityIndices,
			    const SIM_RawIndexField &, regularVelocityIndices,
			    const SIM_RawField &, solidSurface,
			    const SIM_VectorField &, solidVelocity,
			    const Vector &, velocitySolution,
			    const int, axis)

    void applyVelocitiesToRegularGridPartial(SIM_RawField &regularVelocity,
						const HDK_OctreeVectorFieldInterpolator &interpolator,
						const SIM_RawIndexField &octreeVelocityIndices,
						const SIM_RawIndexField &regularVelocityIndices,
						const SIM_RawField &solidSurface,
						const SIM_VectorField &solidVelocity,
						const Vector &velocitySolution,
						const int axis,
						const UT_JobInfo &info) const;

    //
    // Tests to verify assumptions in the adaptive viscosity system.
    //

    THREADED_METHOD5_CONST(HDK_AdaptiveViscosity,
                           octreeVelocityIndices[level][axis].shouldMultiThread(),
                           octreeVelocityGradingUnitTest,
			   bool &, passed,
                           const UT_Array<UT_Array<SIM_RawIndexField>> &, octreeVelocityIndices,
                           const HDK_OctreeGrid &, octreeLabels,
                           const int, axis, const int, level)

    void octreeVelocityGradingUnitTestPartial(bool &passed,
						const UT_Array<UT_Array<SIM_RawIndexField>> &octreeVelocityIndices,
						const HDK_OctreeGrid &octreeLabels,
						const int axis, const int level, 
						const UT_JobInfo &info) const;

    bool octreeVelocityUnitTest(const UT_Array<UT_Array<SIM_RawIndexField>> &octreeVelocityIndices,
				const HDK_OctreeGrid &octreeLabels) const;

    THREADED_METHOD6_CONST(HDK_AdaptiveViscosity,
			    edgeStressIndices[level][axis].shouldMultiThread(),
			    edgeStressUnitTest,
			    bool &, passed,
			    const UT_Array<UT_Array<SIM_RawIndexField>> &, edgeStressIndices,
			    const UT_Array<UT_Array<SIM_RawIndexField>> &, octreeVelocityIndices,
			    const HDK_OctreeGrid &, octreeLabels,
			    const int, axis, const int, level)

    void edgeStressUnitTestPartial(bool &passed,
				    const UT_Array<UT_Array<SIM_RawIndexField>> &edgeStressIndices,
				    const UT_Array<UT_Array<SIM_RawIndexField>> &octreeVelocityIndices,
				    const HDK_OctreeGrid &octreeLabels,
				    const int axis, const int level,
				    const UT_JobInfo &info) const;

    bool edgeStressUnitTest(const UT_Array<UT_Array<SIM_RawIndexField>> &edgeStressIndices,
			    const UT_Array<UT_Array<SIM_RawIndexField>> &velocityIndices,
			    const HDK_OctreeGrid &octree) const;

    THREADED_METHOD6_CONST(HDK_AdaptiveViscosity,
                           centerStressIndices[level].shouldMultiThread(),
			   centerStressUnitTest,
                           bool &, passed,
			   const UT_Array<SIM_RawIndexField> &, centerStressIndices,
                           const UT_Array<UT_Array<SIM_RawIndexField>> &, edgeStressIndices,
                           const UT_Array<UT_Array<SIM_RawIndexField>> &, velocityIndices,
                           const HDK_OctreeGrid &, octreeLabels,
			   const int, level)

    void centerStressUnitTestPartial(bool &passed,
				    const UT_Array<SIM_RawIndexField> &centerStressIndices,
				    const UT_Array<UT_Array<SIM_RawIndexField>> &edgeStressIndices,
				    const UT_Array<UT_Array<SIM_RawIndexField>> &octreeVelocityIndices,
				    const HDK_OctreeGrid &octreeLabels,
				    const int level,
				    const UT_JobInfo &info) const;

    bool centerStresUnitTest(const UT_Array<SIM_RawIndexField> &centerStressIndices,
				const UT_Array<UT_Array<SIM_RawIndexField>> &edgeStressIndices,
				const UT_Array<UT_Array<SIM_RawIndexField>> &velocityIndices,
				const HDK_OctreeGrid &octreeLabels) const;
};

#endif