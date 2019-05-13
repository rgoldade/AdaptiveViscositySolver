#include <queue>

#include "HDK_AdaptiveViscosity.h"

#include <GU/GU_Detail.h>

#include <PRM/PRM_Include.h>

#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_GeometryCopy.h>
#include <SIM/SIM_Object.h>
#include <SIM/SIM_PRMShared.h>
#include <SIM/SIM_ScalarField.h>
#include <SIM/SIM_VectorField.h>

#include <UT/UT_DSOVersion.h>
#include <UT/UT_FixedVector.h>
#include <UT/UT_PerfMonAutoEvent.h>

void
initializeSIM(void *)
{
    IMPLEMENT_DATAFACTORY(HDK_AdaptiveViscosity);
}

// Standard constructor, note that BaseClass was crated by the
// DECLARE_DATAFACTORY and provides an easy way to chain through
// the class hierarchy.
HDK_AdaptiveViscosity::HDK_AdaptiveViscosity(const SIM_DataFactory *factory)
    : BaseClass(factory)
{
}

HDK_AdaptiveViscosity::~HDK_AdaptiveViscosity() {}

const SIM_DopDescription *
HDK_AdaptiveViscosity::getDopDescription()
{
    static PRM_Name theLiquidSurfaceName(GAS_NAME_SURFACE, "Liquid Surface Field");
    static PRM_Default theLiquidSurfaceDefault(0, "surface");

    static PRM_Name theFaceWeightsName("faceWeights", "Face Weights Field");
    static PRM_Default theFaceWeightsDefault(0, "surfaceweights");

    static PRM_Name theLiquidVelocityName(GAS_NAME_VELOCITY, "Liquid Velocity Field");
    static PRM_Default theLiquidVelocityDefault(0, "vel");

    static PRM_Name theViscosityName("viscosity", "Viscosity Field");
    static PRM_Default theViscosityDefault(0, "viscosity");

    static PRM_Name theDensityName(GAS_NAME_DENSITY, "Density Field");
    static PRM_Default theDensityDefault(0, "massdensity");

    static PRM_Name theSolidSurfaceName(GAS_NAME_COLLISION, "Solid Surface Field");
    static PRM_Default theSolidSurfaceDefault(0, "collision");

    static PRM_Name theSolidVelocityName(GAS_NAME_COLLISIONVELOCITY, "Solid Velocity Field");
    static PRM_Default theSolidVelocityDefault(0, "collisionvel");

    static PRM_Name theApplyCollisionWeights("applySolidWeights", "Apply Solid Weights");

    static PRM_Name theToleranceName(SIM_NAME_TOLERANCE, "Relative Solver Tolerance");
    static PRM_Default theToleranceDefault(1e-3);

    static PRM_Name theMaxIterations("maxIterations", "Max Solver Iterations");
    static PRM_Default theMaxIterationsDefault(2500);

    static PRM_Name theExtrapolationName("extrapolation", "Extrapolation");
    static PRM_Default theExtrapolationDefault(0.5);

    static PRM_Name theSupersampling("numberSuperSamples", "Samples Per Axis");

    static PRM_Name theOctreeLevelName("octreeLevels", "Octree Levels");
    static PRM_Name theFineBandwidthName("fineLayerBandwidth", "Fine Layer Bandwidth");

    static PRM_Name theUseEnhancedGradientsName("useEnhancedGradients", "Use Enhanced Gradients");

    static PRM_Name thePrintOctreeName("doOutputOctree", "Output Octree Geometry");
    static PRM_Name theOnlyPrintOctreeName("onlyOutputOctree", "Only Output Octree");

    static PRM_Name theOctreeGeometryName("octreeGeometry", "Octree Geometry");
    static PRM_Default theOctreeGeometryNameDefault(0, "OctreeGeometry");

    static PRM_Template theTemplates[] = {
        PRM_Template(PRM_STRING, 1, &theLiquidSurfaceName, &theLiquidSurfaceDefault),
        PRM_Template(PRM_STRING, 1, &theFaceWeightsName, &theFaceWeightsDefault),

        PRM_Template(PRM_STRING, 1, &theLiquidVelocityName, &theLiquidVelocityDefault),

        PRM_Template(PRM_STRING, 1, &theSolidSurfaceName, &theSolidSurfaceDefault),
        PRM_Template(PRM_STRING, 1, &theSolidVelocityName, &theSolidVelocityDefault),

        PRM_Template(PRM_TOGGLE, 1, &theApplyCollisionWeights, PRMzeroDefaults),

        PRM_Template(PRM_STRING, 1, &theViscosityName, &theViscosityDefault),
        PRM_Template(PRM_STRING, 1, &theDensityName, &theDensityDefault),

        PRM_Template(PRM_FLT, 1, &theToleranceName, &theToleranceDefault),

        PRM_Template(PRM_FLT, 1, &theMaxIterations, &theMaxIterationsDefault),

        PRM_Template(PRM_FLT, 1, &theExtrapolationName, &theExtrapolationDefault),

        PRM_Template(PRM_INT, 1, &theSupersampling, PRMthreeDefaults),

        PRM_Template(PRM_INT, 1, &theOctreeLevelName, PRMfourDefaults),
        PRM_Template(PRM_INT, 1, &theFineBandwidthName, PRMtwoDefaults),

        PRM_Template(PRM_TOGGLE, 1, &theUseEnhancedGradientsName, PRMoneDefaults),

        PRM_Template(PRM_TOGGLE, 1, &thePrintOctreeName, PRMzeroDefaults),
        PRM_Template(PRM_TOGGLE, 1, &theOnlyPrintOctreeName, PRMzeroDefaults),

        PRM_Template(PRM_STRING, 1, &theOctreeGeometryName, &theOctreeGeometryNameDefault),

        PRM_Template()};

    static SIM_DopDescription theDopDescription(true, "HDK_AdaptiveViscosity", "HDK Adaptive Viscosity",
						"$OS", classname(), theTemplates);

    setGasDescription(theDopDescription);

    return &theDopDescription;
}

bool
HDK_AdaptiveViscosity::solveGasSubclass(SIM_Engine &engine, SIM_Object *obj,
                                        SIM_Time time, SIM_Time timestep)
{
    const fpreal dt = timestep;

    ////////////////////////////////////////////
    //
    // Load in all the fields.
    //
    ////////////////////////////////////////////

    const SIM_ScalarField *liquidSurfaceField = getConstScalarField(obj, GAS_NAME_SURFACE);
    SIM_VectorField *liquidVelocity = getVectorField(obj, GAS_NAME_VELOCITY);

    const SIM_ScalarField *solidSurfaceField = getConstScalarField(obj, GAS_NAME_COLLISION);
    const SIM_VectorField *solidVelocity = getConstVectorField(obj, GAS_NAME_COLLISIONVELOCITY);

    const SIM_VectorField *faceIntegrationWeights = getConstVectorField(obj, "faceWeights");

    ////////////////////////////////////////////
    //
    // Verify fields are structured as expected.
    //
    ////////////////////////////////////////////

    if (liquidVelocity == nullptr)
    {
        addError(obj, SIM_MESSAGE, "Liquid velocity field missing", UT_ERROR_WARNING);
        return false;
    }
    else if (!liquidVelocity->isFaceSampled())
    {
        addError(obj, SIM_MESSAGE, "Liquid velocity field must be a staggered grid", UT_ERROR_WARNING);
        return false;
    }

    if (faceIntegrationWeights == nullptr)
    {
        addError(obj, SIM_MESSAGE, "Face weights field missing",
                 UT_ERROR_WARNING);
        return false;
    }
    else if (!faceIntegrationWeights->isAligned(liquidVelocity))
    {
        addError(obj, SIM_MESSAGE, "Face weights must align with velocity samples", UT_ERROR_WARNING);
        return false;
    }

    if (solidSurfaceField == nullptr)
    {
        addError(obj, SIM_MESSAGE, "Solid surface field missing", UT_ERROR_WARNING);
        return false;
    }

    const SIM_RawField &solidSurface = *solidSurfaceField->getField();

    if (solidVelocity == nullptr)
    {
        addError(obj, SIM_MESSAGE, "Solid velocity field missing", UT_ERROR_WARNING);
        return false;
    }

    if (liquidSurfaceField == nullptr)
    {
        addError(obj, SIM_MESSAGE, "Liquid surface field is missing", UT_ERROR_WARNING);
        return false;
    }

    const SIM_RawField &liquidSurface = *liquidSurfaceField->getField();

    ////////////////////////////////////////////
    //
    // Get fluid parameters
    //
    ////////////////////////////////////////////

    const SIM_ScalarField *viscosityField = getConstScalarField(obj, "viscosity");

    if (viscosityField == nullptr)
    {
        addError(obj, SIM_MESSAGE, "Viscosity field is missing", UT_ERROR_WARNING);
        return false;
    }
    else if (!viscosityField->getField()->isAligned(&liquidSurface))
    {
        addError(obj, SIM_MESSAGE, "Viscosity field must align with the surface volume", UT_ERROR_WARNING);
        return false;
    }

    const SIM_RawField &viscosity = *viscosityField->getField();

    const SIM_ScalarField *densityField = getConstScalarField(obj, GAS_NAME_DENSITY);

    if (densityField == nullptr)
    {
	addError(obj, SIM_MESSAGE, "Density field is missing", UT_ERROR_WARNING);
        return false;
    }
    else if (!densityField->getField()->isAligned(&liquidSurface))
    {
        addError(obj, SIM_MESSAGE, "Density field must align with the surface volume", UT_ERROR_WARNING);
        return false;
    }

    const SIM_RawField &density = *densityField->getField();    

    ////////////////////////////////////////////
    //
    // Build control volume weights for finest grid level
    //
    ////////////////////////////////////////////

    SIM_RawField centerIntegrationWeights;
    SIM_RawField edgeIntegrationWeights[3];

    const fpreal dx = liquidVelocity->getVoxelSize().maxComponent();
    const fpreal extrapolation = dx * getExtrapolation();
    const bool doApplySolidWeights = getDoApplySolidWeights();

    buildIntegrationWeights(centerIntegrationWeights,
			    edgeIntegrationWeights,
			    liquidSurface,
			    solidSurface,
			    extrapolation,
			    doApplySolidWeights);

    ////////////////////////////////////////////
    //
    // Build octree from liquid and solid surfaces
    //
    ////////////////////////////////////////////

    fpreal fineVoxelWidth = SYSmax(2., fpreal(getFineBandwidth()));

    const fpreal innerFineBandwidth = dx * fineVoxelWidth;
    const fpreal outerFineBandwidth = 3. * dx;
    const int desiredOctreeLevels = getOctreeLevels();

    HDK_OctreeGrid octreeLabels;

    buildOctree(octreeLabels,
		liquidSurface,
		solidSurface,
		desiredOctreeLevels,
                extrapolation,
		outerFineBandwidth,
		innerFineBandwidth);

    const int octreeLevels = octreeLabels.getOctreeLevels();

    ////////////////////////////////////////////
    //
    // Dump out geometry of the octree
    //
    ////////////////////////////////////////////

    if (getOnlyPrintOctree())
    {
        SIM_GeometryCopy *octreeGeometry = getOrCreateGeometry(obj, "octreeGeometry");

        SIM_GeometryAutoWriteLock autoLockOctree(octreeGeometry, SIM_DATA_ID_PRESERVE);
        GU_Detail *octreeDetail = &autoLockOctree.getGdp();

        octreeLabels.outputOctreeGeometry(*octreeDetail);

        if (getOnlyPrintOctree())
            return true;
    }

    ////////////////////////////////////////////
    //
    // Build mapping of regular grid velocities
    // faces to octree faces.
    //
    ////////////////////////////////////////////

    SIM_RawIndexField regularVelocityIndices[3];
    exint regularVelocityDOFcount = 0;
    {
        UT_PerfMonAutoSolveEvent event(this, "Build Regular Grid Velocity Labels");

        UT_Vector3 size = liquidSurface.getSize();
        UT_Vector3 orig = liquidSurface.getOrig();

	UT_Vector3i voxelRes;
        liquidSurface.getVoxelRes(voxelRes[0], voxelRes[1], voxelRes[2]);

        regularVelocityIndices[0].init(SIM_SAMPLE_FACEX, orig, size,
					voxelRes[0], voxelRes[1], voxelRes[2]);

        regularVelocityIndices[1].init(SIM_SAMPLE_FACEY, orig, size,
					voxelRes[0], voxelRes[1], voxelRes[2]);

        regularVelocityIndices[2].init(SIM_SAMPLE_FACEZ, orig, size,
					voxelRes[0], voxelRes[1], voxelRes[2]);

        regularVelocityDOFcount = buildRegularVelocityIndices(regularVelocityIndices,
								liquidSurface,
								solidSurface,
								centerIntegrationWeights,
								edgeIntegrationWeights,
								extrapolation);
    }

    ////////////////////////////////////////////
    //
    // Build velocity and edge indices
    //
    ////////////////////////////////////////////

    UT_Array<UT_Array<SIM_RawIndexField>> octreeVelocityIndices;
    octreeVelocityIndices.setSize(octreeLevels);

    UT_Array<UT_Array<SIM_RawIndexField>> edgeStressIndices;
    edgeStressIndices.setSize(octreeLevels);

    UT_Array<SIM_RawIndexField> centerStressIndices;
    centerStressIndices.setSize(octreeLevels);

    for (int level = 0; level < octreeLevels; ++level)
    {
        octreeVelocityIndices[level].setSize(3);
        edgeStressIndices[level].setSize(3);
    }

    // Initialize velocity and edge labels to construct the DOFs and deformation
    // rate labels

    exint octreeVelocityDOFCount = 0;
    exint edgeStressDOFCount = 0;
    exint centerStressDOFCount = 0;

    {
        UT_PerfMonAutoSolveEvent event(this, "Build Octree Velocity and Stress Labels");

        UT_Vector3 size = octreeLabels.getSize();
        UT_Vector3 orig = octreeLabels.getOrig();

        for (int level = 0; level < octreeLevels; ++level)
        {
            UT_Vector3i voxelRes = octreeLabels.getVoxelRes(level);

            // Deal with the proper resolution offsets for sample differences?
            octreeVelocityIndices[level][0].init(SIM_SAMPLE_FACEX, orig, size,
						    voxelRes[0], voxelRes[1], voxelRes[2]);

            octreeVelocityIndices[level][1].init(SIM_SAMPLE_FACEY, orig, size,
						    voxelRes[0], voxelRes[1], voxelRes[2]);

            octreeVelocityIndices[level][2].init(SIM_SAMPLE_FACEZ, orig, size,
						    voxelRes[0], voxelRes[1], voxelRes[2]);

            // This is the most intuitive way based on how I've implemented
            // edge-sampled grids in the past. Houdini indicates edges by the
            // plane they form if they're the normal vector.
            edgeStressIndices[level][0].init(HDK_XEDGE, orig, size,
						voxelRes[0], voxelRes[1], voxelRes[2]); // X-directed edge

            edgeStressIndices[level][1].init(HDK_YEDGE, orig, size,
						voxelRes[0], voxelRes[1], voxelRes[2]); // Y-directed edge

            edgeStressIndices[level][2].init(HDK_ZEDGE, orig, size,
						voxelRes[0], voxelRes[1], voxelRes[2]); // Z-directed edge

            centerStressIndices[level].init(SIM_SAMPLE_CENTER, orig, size,
						voxelRes[0], voxelRes[1], voxelRes[2]);
        }

        octreeVelocityDOFCount = buildOctreeVelocityIndices(octreeVelocityIndices,
							    liquidSurface,
							    solidSurface,
							    octreeLabels,
							    centerIntegrationWeights,
							    edgeIntegrationWeights,
							    extrapolation);

	edgeStressDOFCount = buildEdgeStressIndices(edgeStressIndices,
						    liquidSurface,
						    octreeLabels,
						    edgeIntegrationWeights);

        centerStressDOFCount = buildCenterStressIndices(centerStressIndices, octreeLabels, centerIntegrationWeights);

#if !defined(NDEBUG)
        assert(octreeVelocityUnitTest(octreeVelocityIndices, octreeLabels));
	assert(edgeStressUnitTest(edgeStressIndices, octreeVelocityIndices, octreeLabels));
	assert(centerStresUnitTest(centerStressIndices, edgeStressIndices, octreeVelocityIndices, octreeLabels));
#endif

    }

    ////////////////////////////////////////////
    //
    // Precompute stress gradients.
    //
    ////////////////////////////////////////////

    // The final system is (Mu + 2. * dt * D^T * K * Mtau * U * D) * u^{n+1} = Mu * u^n.
    // However we don't want to perform sparse matrix multiplication to
    // build our linear system. Instead, we will precompute the rows of D now
    // and build the whole system in a single pass after.

    UT_Array<UT_Array<StressStencilFace>> edgeStressStencils;
    edgeStressStencils.setSize(edgeStressDOFCount);

    UT_Array<UT_Array<fpreal>> edgeStressBoundaryStencils;
    edgeStressBoundaryStencils.setSize(edgeStressDOFCount);

    UT_Array<fpreal> edgeStressWeights;
    edgeStressWeights.setSize(edgeStressDOFCount);

    const bool useEnhancedGradients = getUseEnhancedGradients();

    {
        UT_PerfMonAutoSolveEvent event(this, "Build Edge Stress Stencils");

        for (int level = 0; level < octreeLevels; ++level)
            for (int axis : {0,1,2})
            {
                edgeStressStencilParms parms(octreeVelocityIndices,
						edgeStressIndices,
						octreeLabels,
						*solidVelocity,
						edgeIntegrationWeights,
						viscosity,
						dt,
						useEnhancedGradients);

                buildEdgeStressStencils(edgeStressStencils,
					edgeStressBoundaryStencils,
					edgeStressWeights,
					parms,
					axis, level);
            }
    }

    UT_Array<UT_Array<StressStencilFace>> centerStressStencils;
    centerStressStencils.setSize(centerStressDOFCount * 3);

    UT_Array<UT_Array<fpreal>> centerStressBoundaryStencils;
    centerStressBoundaryStencils.setSize(centerStressDOFCount * 3);

    UT_Array<fpreal> centerStressWeights;
    centerStressWeights.setSize(centerStressDOFCount);

    {
        UT_PerfMonAutoSolveEvent event(this, "Build Cell Stress Stencils");

        for (int level = 0; level < octreeLevels; ++level)
        {
            for (int axis : {0,1,2})
            {
                centerStressStencilParms parms(octreeVelocityIndices,
						centerStressIndices,
						octreeLabels,
						*solidVelocity);

                buildCenterStressStencils(centerStressStencils,
					    centerStressBoundaryStencils,
					    parms,
					    axis, level,
					    centerStressDOFCount);
            }

            // Because the cell centered stresses are all co-located, we only
            // need on weight value for all three axes.
            buildCenterStressWeights(centerStressWeights,
					centerStressIndices[level],
					centerIntegrationWeights,
					viscosity,
					dt, level);
        }
    }

    ////////////////////////////////////////////
    //
    // Map regular grid stress to octree system
    //
    ////////////////////////////////////////////

#ifdef USEEIGEN
    Vector viscositySolution = Vector::Zero(octreeVelocityDOFCount);
#else
    Vector viscositySolution;
    viscositySolution.init(0, octreeVelocityDOFCount - 1);
    viscositySolution.zero();
#endif

    {
        UT_PerfMonAutoSolveEvent event(this, "Interpolate Regular Grid Velocities at Octree Velocity Faces");

        for (int level = 0; level < octreeLevels; ++level)
            for (int axis : {0,1,2})
            {
                buildVelocityMapping(viscositySolution,
					*liquidVelocity->getField(axis),
					regularVelocityIndices[axis],
					octreeVelocityIndices[level][axis],
					*solidVelocity,
					octreeLabels,
					axis, level);
            }
    }

    ////////////////////////////////////////////
    //
    // Build linear system for implicit solve
    //
    ////////////////////////////////////////////

#ifdef USEEIGEN
    std::vector<Eigen::Triplet<SolveType>> sparseMatrixElements;
    sparseMatrixElements.reserve(octreeVelocityDOFCount * 15);

    Vector rhs = Vector::Zero(octreeVelocityDOFCount);

#else
    SparseMatrix sparseMatrixElements;
    sparseMatrixElements.init(octreeVelocityDOFCount, octreeVelocityDOFCount);
    sparseMatrixElements.reserve(octreeVelocityDOFCount * 15);

    Vector rhs;
    rhs.init(0, octreeVelocityDOFCount - 1);
    rhs.zero();
#endif

    {
        UT_PerfMonAutoSolveEvent event(this, "Build Octree Linear System");

        const int threadCount = UT_Thread::getNumProcessors();

#ifdef USEEIGEN
        std::vector<std::vector<Eigen::Triplet<SolveType>>> parallelSparseMatrixElements(threadCount);
#else
        UT_Array<SparseMatrix> parallelSparseMatrixElements;
        parallelSparseMatrixElements.setSize(threadCount);

        for (int threads = 0; threads < threadCount; ++threads)
        {
            parallelSparseMatrixElements[threads].init(octreeVelocityDOFCount, octreeVelocityDOFCount);
            parallelSparseMatrixElements[threads].reserve(octreeVelocityDOFCount * 14. / fpreal(threadCount));
        }
#endif

        octreeSystemStencilParms parms(edgeStressStencils, edgeStressBoundaryStencils, edgeStressWeights,
					centerStressStencils, centerStressBoundaryStencils, centerStressWeights,
					octreeVelocityIndices, edgeStressIndices, centerStressIndices,
					octreeLabels, *faceIntegrationWeights, centerIntegrationWeights,
					density, viscositySolution, useEnhancedGradients);

        for (int level = 0; level < octreeLevels; ++level)
            for (int axis : {0,1,2})
            {
                buildOctreeSystemFromStencils(parallelSparseMatrixElements,
						rhs,
						parms,
						axis, level,
						centerStressDOFCount);
            }

#ifdef USEEIGEN
        for (int threads = 0; threads < threadCount; ++threads)
            sparseMatrixElements.insert(sparseMatrixElements.end(), parallelSparseMatrixElements[threads].begin(), parallelSparseMatrixElements[threads].end());
#else
        for (int threads = 0; threads < threadCount; ++threads)
            sparseMatrixElements += parallelSparseMatrixElements[threads];
#endif
    }

    ////////////////////////////////////////////
    //
    // Iterative solve
    //
    ////////////////////////////////////////////

    {
        UT_PerfMonAutoSolveEvent event(this, "Solve linear system");

        const fpreal solverTolerance = getSolverTolerance();
        const int maxSolverIterations = getMaxIterations();
	
	fpreal solverError = 0;
	int numberOfIterations = 0;

#ifdef USEEIGEN

        Eigen::SparseMatrix<SolveType> sparseMatrix(octreeVelocityDOFCount, octreeVelocityDOFCount);
        sparseMatrix.setFromTriplets(sparseMatrixElements.begin(), sparseMatrixElements.end());

        const int eigenThreads = Eigen::nbThreads();

        Eigen::ConjugateGradient<Eigen::SparseMatrix<SolveType>, Eigen::Lower | Eigen::Upper> solver;
        solver.compute(sparseMatrix);

        if (solver.info() != Eigen::Success)
	    return false;

        solver.setTolerance(solverTolerance);
        solver.setMaxIterations(maxSolverIterations);
        
	viscositySolution = solver.solveWithGuess(rhs, viscositySolution);

	numberOfIterations = solver.iterations();
	solverError = solver.error();

#else
        sparseMatrixElements.compile();

        UT_SparseMatrixRowT<SolveType> sparseMatrix;
        sparseMatrix.buildFrom(sparseMatrixElements);

	solverError = sparseMatrix.solveConjugateGradient(viscositySolution,
							    rhs, nullptr,
							    solverTolerance,
							    maxSolverIterations,
							    &numberOfIterations);
#endif

	UT_WorkBuffer extrainfo;
        extrainfo.sprintf("iterations=%d, error=%.6f, octree DOFS=%d, regular DOFs=%d",
			    numberOfIterations,
			    solverError,
			    octreeVelocityDOFCount,
			    regularVelocityDOFcount);
        
	event.setExtraInfo(extrainfo.buffer());
    }

    ////////////////////////////////////////////
    //
    // Apply octree velocites back to reguilar grid
    //
    ////////////////////////////////////////////

    {
        UT_PerfMonAutoSolveEvent event(this, "Apply octree solution to regular grid");

    	UT_Array<UT_Array<SIM_RawField>> octreeVelocity;
	octreeVelocity.setSize(octreeLevels);

	{
	    UT_Vector3 size = octreeLabels.getSize(), orig = octreeLabels.getOrig();

	    for (int level = 0; level < octreeLevels; ++level)
	    {
		octreeVelocity[level].setSize(3);

		UT_Vector3i voxelRes = octreeLabels.getVoxelRes(level);

		// Deal with the proper resolution offsets for sample differences?
		octreeVelocity[level][0].init(SIM_SAMPLE_FACEX, orig, size,
						voxelRes[0], voxelRes[1], voxelRes[2]);

		octreeVelocity[level][1].init(SIM_SAMPLE_FACEY, orig, size,
						voxelRes[0], voxelRes[1], voxelRes[2]);
		
		octreeVelocity[level][2].init(SIM_SAMPLE_FACEZ, orig, size,
						voxelRes[0], voxelRes[1], voxelRes[2]);

		for (int axis : {0,1,2})
		{
		    octreeVelocity[level][axis].makeConstant(0);
		    setOctreeVelocity(octreeVelocity[level][axis], octreeVelocityIndices[level][axis], viscositySolution);
		}
	    }
	}

	HDK_OctreeVectorFieldInterpolator interpolator(octreeLabels, octreeVelocity, octreeVelocityIndices);

	for (int axis : {0,1,2})
        {
            applyVelocitiesToRegularGrid(*liquidVelocity->getField(axis),
					    interpolator,
					    octreeVelocityIndices[0][axis],
					    regularVelocityIndices[axis],
					    solidSurface,
					    *solidVelocity,
					    viscositySolution,
					    axis);
        }
    }

    return true;
}

void
computeIntegrationWeights(SIM_RawField &integrationWeights,
			    const SIM_RawField &liquidSurface,
			    const SIM_FieldSample sample,
			    const int numberOfSamples)
{
    UT_Vector3 size = liquidSurface.getSize(), orig = liquidSurface.getOrig();

    UT_Vector3i voxelRes;
    liquidSurface.getVoxelRes(voxelRes[0], voxelRes[1], voxelRes[2]);

    integrationWeights.init(sample, orig, size, voxelRes[0], voxelRes[1], voxelRes[2]);
    integrationWeights.makeConstant(0);
    integrationWeights.computeSDFWeightsSampled(&liquidSurface, numberOfSamples, false, 0);
}

void
computeSolidIntegrationWeights(SIM_RawField &solidIntegrationWeights,
				const SIM_RawField &liquidSurface,
				const SIM_RawField &solidSurface,
				const SIM_FieldSample sample,
				const int numberOfSamples,
				const fpreal extrapolation)
{
    UT_Vector3 size = liquidSurface.getSize(), orig = liquidSurface.getOrig();
    
    // We want the solid weights to match the surface field
    UT_Vector3i voxelRes;
    liquidSurface.getVoxelRes(voxelRes[0], voxelRes[1], voxelRes[2]);

    solidIntegrationWeights.init(sample, orig, size, voxelRes[0], voxelRes[1], voxelRes[2]);
    solidIntegrationWeights.makeConstant(0);
    solidIntegrationWeights.computeSDFWeightsSampled(&solidSurface, numberOfSamples,
							false, 0, -extrapolation);
}

void
HDK_AdaptiveViscosity::buildIntegrationWeights(SIM_RawField &centerIntegrationWeights,
						SIM_RawField (&edgeIntegrationWeights)[3],
						const SIM_RawField &liquidSurface,
						const SIM_RawField &solidSurface,
						const fpreal extrapolation,
						const bool doApplySolidWeights) const
{
    const int numberOfSamples = getNumberSuperSamples();

    {
        UT_PerfMonAutoSolveEvent event(this, "Compute Surface Weights");

        computeIntegrationWeights(centerIntegrationWeights, liquidSurface, SIM_SAMPLE_CENTER, numberOfSamples);

        computeIntegrationWeights(edgeIntegrationWeights[0], liquidSurface, HDK_XEDGE, numberOfSamples);
        computeIntegrationWeights(edgeIntegrationWeights[1], liquidSurface, HDK_YEDGE, numberOfSamples);
        computeIntegrationWeights(edgeIntegrationWeights[2], liquidSurface, HDK_ZEDGE, numberOfSamples);
    }

    // The collision weights provide the "theta" value for applying ghost fluid
    // collision velocities The weights can then be directly applied to the
    // normal variational weights as they are "1" for any purely air/liquid
    // voxel
    if (doApplySolidWeights)
    {
        SIM_RawField solidCenterWeights, solidEdgeWeights[3];

        UT_PerfMonAutoSolveEvent event(this, "Compute Collision Weights");

        computeSolidIntegrationWeights(solidCenterWeights, liquidSurface, solidSurface,
					SIM_SAMPLE_CENTER, numberOfSamples, extrapolation);

	centerIntegrationWeights.setScaleDivideThreshold(1, nullptr, &solidCenterWeights, 0);

        computeSolidIntegrationWeights(solidEdgeWeights[0], liquidSurface, solidSurface, HDK_XEDGE, numberOfSamples, extrapolation);
        computeSolidIntegrationWeights(solidEdgeWeights[1], liquidSurface, solidSurface, HDK_YEDGE, numberOfSamples, extrapolation);
        computeSolidIntegrationWeights(solidEdgeWeights[2], liquidSurface, solidSurface, HDK_ZEDGE, numberOfSamples, extrapolation);

        edgeIntegrationWeights[0].setScaleDivideThreshold(1, nullptr, &solidEdgeWeights[0], 0);
        edgeIntegrationWeights[1].setScaleDivideThreshold(1, nullptr, &solidEdgeWeights[1], 0);
        edgeIntegrationWeights[2].setScaleDivideThreshold(1, nullptr, &solidEdgeWeights[2], 0);
    }
}

void
HDK_AdaptiveViscosity::buildOctree(HDK_OctreeGrid &octreeLabels,
				   const SIM_RawField &liquidSurface,
				   const SIM_RawField &solidSurface,
                                   const int octreeLevels,
				   const fpreal extrapolation,
                                   const fpreal outerFineBandwidth,
                                   const fpreal innerFineBandwidth) const
{
    // Create a mask of the surface volume to refine the octree
    UT_Vector3i voxelRes;
    liquidSurface.getVoxelRes(voxelRes[0], voxelRes[1], voxelRes[2]);

    SIM_RawField fineCellMask;
    fineCellMask.init(SIM_SAMPLE_CENTER,
			liquidSurface.getOrig(), liquidSurface.getSize(),
			voxelRes[0], voxelRes[1], voxelRes[2]);
    fineCellMask.makeConstant(1.);

    {
        UT_PerfMonAutoSolveEvent event(this, "Build Mask for Octree");

        const auto buildMaskFunctor = [&](const UT_BlockedRange<int64> &range)
	{
            UT_VoxelArrayIteratorF vit;
            vit.setConstArray(liquidSurface.field());
            UT_VoxelTileIteratorF vitt;

            for (int64 i = range.begin(); i != range.end(); ++i)
            {
                vit.myTileStart = i;
                vit.myTileEnd = i + 1;
                vit.rewind();

                if (!vit.atEnd())
                {
                    if (!vit.isTileConstant() || vit.getValue() < outerFineBandwidth)
                    {
                        vitt.setTile(vit);

                        for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
                        {
                            UT_Vector3i cell(vitt.x(), vitt.y(), vitt.z());

                            // Voxels within the band are set to zero, indicating that
			    // these voxels should be active in the octree
                            fpreal sdf = vitt.getValue();

                            fpreal maskLabel = 0;
                            if (sdf > 0 && sdf < outerFineBandwidth)
                                maskLabel = 0;
                            else if (sdf <= 0.)
                            {
                                if (sdf > -innerFineBandwidth)
                                    maskLabel = 0;
                                else
                                {
                                    UT_Vector3 point;
                                    fineCellMask.indexToPos(vitt.x(), vitt.y(), vitt.z(), point);

				    if (solidSurface.getValue(point) > (-innerFineBandwidth - extrapolation))
                                        maskLabel = 0;
                                    else
                                        maskLabel = -1;
                                }
                            }
                            else
                                maskLabel = 1;

			    HDKsetFieldValue(fineCellMask, cell, maskLabel);
                        }
                    }
                }
            }
        };

        int64 tiles = fineCellMask.field()->numTiles();
        UTparallelForEachNumber(tiles, buildMaskFunctor);
    }

    {
        UT_PerfMonAutoSolveEvent event(this, "Build Octree");
        octreeLabels.init(fineCellMask, octreeLevels);
    }

#if !defined(NDEBUG)
    {
        UT_PerfMonAutoSolveEvent event(this, "Unit Test Octree");
	assert(octreeLabels.unitTest());
    }
#endif
}

void
HDK_AdaptiveViscosity::findOccupiedRegularVelocityTilesPartial(UT_Array<bool> &isTileOccupiedList,
								const SIM_RawField &liquidSurface,
								const SIM_RawIndexField &octreeVelocityIndices,
								const int axis,
								const UT_JobInfo &info) const
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorF vit;
    vit.setConstArray(liquidSurface.field());
    vit.splitByTile(info);

    const exint tileCount = isTileOccupiedList.entries();

    UT_Array<bool> localIsTileOccupiedList;
    localIsTileOccupiedList.setSize(tileCount);
    localIsTileOccupiedList.constant(false);

    UT_VoxelTileIteratorF vitt;

    fpreal sdf = 2. * liquidSurface.getVoxelSize().maxComponent();

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        if (!vit.isTileConstant() || vit.getValue() < sdf)
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                if (vitt.getValue() < sdf)
                {
                    UT_Vector3i cell(vitt.x(), vitt.y(), vitt.z());

		    for (int direction : {0,1})
                    {
			UT_Vector3i face = HDKcellToFace(cell, axis, direction);

                        exint tileNumber = octreeVelocityIndices.field()->indexToLinearTile(face[0],
											    face[1],
											    face[2]);
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
HDK_AdaptiveViscosity::findOccupiedOctreeVelocityTilesPartial(UT_Array<bool> &isTileOccupiedList,
								const SIM_RawField &octreeLabels,
								const SIM_RawIndexField &octreeVelocityIndices,
								const int axis,
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

                    for (int direction : {0,1})
                    {
			UT_Vector3i face = HDKcellToFace(cell, axis, direction);

                        exint tileNumber = octreeVelocityIndices.field()->indexToLinearTile(face[0],
											    face[1],
											    face[2]);
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
HDK_AdaptiveViscosity::findOccupiedEdgeStressTilesPartial(UT_Array<bool> &isTileOccupiedList,
							    const SIM_RawField &octreeLabels,
							    const SIM_RawIndexField &edgeStressIndices,
							    const int axis,
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

		    for (int edgeIndex = 0; edgeIndex < 4; ++edgeIndex)
                    {
			UT_Vector3i edge = HDKcellToEdge(cell, axis, edgeIndex);

                        exint tileNumber = edgeStressIndices.field()->indexToLinearTile(edge[0],
											    edge[1],
											    edge[2]);
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
HDK_AdaptiveViscosity::uncompressTilesPartial(SIM_RawIndexField &grid,
					      const UT_Array<bool> &isTileOccupiedList,
                                              const UT_JobInfo &info) const
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
        {
            grid.field()->getLinearTile(i)->uncompress();
        }
    }
}

void
HDK_AdaptiveViscosity::classifyRegularVelocityFacesPartial(SIM_RawIndexField &regularVelocityIndices,
							    const SIM_RawField &liquidSurface,
							    const SIM_RawField &solidSurface,
							    const SIM_RawField &centerIntegrationWeights,
							    const SIM_RawField *edgeIntegrationWeights,
							    const int axis,
							    const fpreal extrapolation,
							    const UT_JobInfo &info) const
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorI vit(regularVelocityIndices.fieldNC());
    vit.splitByTile(info);

    UT_VoxelTileIteratorI vitt;

    UT_Vector3i voxelRes;
    liquidSurface.getVoxelRes(voxelRes[0], voxelRes[1], voxelRes[2]);

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        // Tiles for the index grid must have already been uncompressed
        if (!vit.isTileConstant())
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                UT_Vector3i face(vitt.x(), vitt.y(), vitt.z());

                UT_Vector3i backwardCell = HDKfaceToCell(face, axis, 0);
		UT_Vector3i forwardCell = HDKfaceToCell(face, axis, 1);

		if (backwardCell[axis] < 0 || forwardCell[axis] >= voxelRes[axis])
                    continue;

                bool isActiveVelocity = false;

		if (HDKgetFieldValue(centerIntegrationWeights, backwardCell) > 0. ||
		    HDKgetFieldValue(centerIntegrationWeights, forwardCell) > 0.)
                    isActiveVelocity = true;

                if (!isActiveVelocity)
                {
		    for (int edgeAxis = 0; edgeAxis < 3 && !isActiveVelocity; ++edgeAxis)
		    {
			if (edgeAxis == axis) continue;

			for (int direction : {0,1})
			{
			    UT_Vector3i edge = HDKfaceToEdge(face, axis, edgeAxis, direction);

			    if (HDKgetFieldValue(edgeIntegrationWeights[edgeAxis], edge) > 0)
			    {
				isActiveVelocity = true;
				break;
			    }
			}
                    }
                }

                if (isActiveVelocity)
                {
		    UT_Vector3 point;
                    regularVelocityIndices.indexToPos(face[0], face[1], face[2], point);

                    if (solidSurface.getValue(point) > -extrapolation)
                        vitt.setValue(HDK_SOLIDBOUNDARY);
                    else
                        vitt.setValue(HDK_FLUID);
                }
            }
        }
    }
}

void
HDK_AdaptiveViscosity::classifyOctreeVelocityFacesPartial(SIM_RawIndexField &octreeVelocityIndices,
							    const HDK_OctreeGrid &octreeLabels,
							    const SIM_RawField &liquidSurface,
							    const SIM_RawField &solidSurface,
							    const SIM_RawField(&centerIntegrationWeights),
							    const SIM_RawField *edgeIntegrationWeights,
							    const int axis, const int level,
							    const fpreal extrapolation,
							    const UT_JobInfo &info) const
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorI vit(octreeVelocityIndices.fieldNC());
    vit.splitByTile(info);

    UT_VoxelTileIteratorI vitt;

    UT_Vector3i voxelRes;
    if (level == 0)
        liquidSurface.getVoxelRes(voxelRes[0], voxelRes[1], voxelRes[2]);
    else
        voxelRes = octreeLabels.getVoxelRes(level);

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        // Tiles for the index grid must have already been uncompressed
        if (!vit.isTileConstant())
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                UT_Vector3i face(vitt.x(), vitt.y(), vitt.z());

                UT_Vector3i backwardCell = HDKfaceToCell(face, axis, 0);
		UT_Vector3i forwardCell = HDKfaceToCell(face, axis, 1);

                // Boundary faces and faces outside of the surface can get
                // preemptively set to outside
                if (backwardCell[axis] < 0 || forwardCell[axis] >= voxelRes[axis])
                {
                    if (level == 0)
                        vitt.setValue(HDK_OUTSIDE);
                    continue;
                }

		const int backwardLabel = octreeLabels.getCellLabel(backwardCell, level);
                const int forwardLabel = octreeLabels.getCellLabel(forwardCell, level);

                // The finest level needs to be processed in two pieces.
                // The faces without any grading use the free surface, regular
                // grid, approach for testing labels.
                // The faces with grading are far enough away from the free
                // surface to only only need the octree method.

                if (level == 0)
                {
                    // A face with two ACTIVE adjacent cells at the same level
                    // could be close to the free surface. We need to check if
                    // there are non-zero volumes at stress positions in order
                    // to activate the velocity sample.
                    if (backwardLabel == HDK_OctreeGrid::ACTIVE &&
                        forwardLabel == HDK_OctreeGrid::ACTIVE)
                    {
			bool isActiveVelocity = false;

                        if (HDKgetFieldValue(centerIntegrationWeights, backwardCell) > 0. ||
			    HDKgetFieldValue(centerIntegrationWeights, forwardCell) > 0.)
			    isActiveVelocity = true;

                        if (!isActiveVelocity)
                        {
 			    for (int edgeAxis = 0; edgeAxis < 3 && !isActiveVelocity; ++edgeAxis)
			    {
				if (edgeAxis == axis) continue;

				for (int direction : {0,1})
				{
				    UT_Vector3i edge = HDKfaceToEdge(face, axis, edgeAxis, direction);

				    if (HDKgetFieldValue(edgeIntegrationWeights[edgeAxis], edge) > 0)
				    {
					isActiveVelocity = true;
					break;
				    }
				}
			    }
                        }

                        if (isActiveVelocity)
                        {
                            UT_Vector3 point;
                            octreeVelocityIndices.indexToPos(face[0], face[1], face[2], point);

			    if (solidSurface.getValue(point) > -extrapolation)
				vitt.setValue(HDK_SOLIDBOUNDARY);
                            else
                                vitt.setValue(HDK_FLUID);
                        }
                        else
                            vitt.setValue(HDK_OUTSIDE);
                    }
                    else if (backwardLabel == HDK_OctreeGrid::INACTIVE || forwardLabel == HDK_OctreeGrid::INACTIVE)
                        vitt.setValue(HDK_OUTSIDE);
                    // A face between an ACTIVE and an UP cell must be inside
                    // the fluid domain since it is the transition from a fine
                    // to a coarse grid cell in the simulation.
                    else if ((backwardLabel == HDK_OctreeGrid::UP && forwardLabel == HDK_OctreeGrid::ACTIVE) ||
                             (backwardLabel == HDK_OctreeGrid::ACTIVE && forwardLabel == HDK_OctreeGrid::UP))
                    {
                        vitt.setValue(HDK_FLUID);

#if !defined(NDEBUG)
			UT_Vector3 point;
                        octreeVelocityIndices.indexToPos(face[0], face[1], face[2], point);

			assert(solidSurface.getValue(point) <= -extrapolation);
			assert(liquidSurface.getValue(point) < 0.);
#endif
                    }
                    else
                    {
                        assert(!((backwardLabel == HDK_OctreeGrid::ACTIVE && forwardLabel == HDK_OctreeGrid::DOWN) ||
				(backwardLabel == HDK_OctreeGrid::DOWN && forwardLabel == HDK_OctreeGrid::ACTIVE)) &&
				!((backwardLabel == HDK_OctreeGrid::DOWN && forwardLabel == HDK_OctreeGrid::UP) ||
				(backwardLabel == HDK_OctreeGrid::UP && forwardLabel == HDK_OctreeGrid::DOWN)) &&
				!((backwardLabel == HDK_OctreeGrid::DOWN && forwardLabel == HDK_OctreeGrid::INACTIVE) ||
				(backwardLabel == HDK_OctreeGrid::INACTIVE && forwardLabel == HDK_OctreeGrid::DOWN)));
                    }
                }
                else
                {
                    // The three possible instances of an active face.
                    if ((backwardLabel == HDK_OctreeGrid::ACTIVE && forwardLabel == HDK_OctreeGrid::ACTIVE) ||
			(backwardLabel == HDK_OctreeGrid::UP && forwardLabel == HDK_OctreeGrid::ACTIVE) ||
			(backwardLabel == HDK_OctreeGrid::ACTIVE && forwardLabel == HDK_OctreeGrid::UP))
                    {
#if !defined(NDEBUG)
                        // No octree-level fluid should ever be inside a
                        // collision so we don't need to check for it
                        UT_Vector3 point;
                        octreeVelocityIndices.indexToPos(face[0], face[1], face[2], point);

			assert(solidSurface.getValue(point) <= -extrapolation);
			assert(liquidSurface.getValue(point) < 0);
#endif
                        vitt.setValue(HDK_FLUID);
                    }
                }
            }
        }
    }
}

void
HDK_AdaptiveViscosity::classifyEdgeStressesPartial(SIM_RawIndexField &edgeStressIndices,
						    const HDK_OctreeGrid &octreeLabels,
						    const SIM_RawField &liquidSurface,
						    const SIM_RawField &edgeIntegrationWeights,
						    const int axis, const int level, 
						    const UT_JobInfo &info) const
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorI vit(edgeStressIndices.fieldNC());
    vit.splitByTile(info);

    UT_VoxelTileIteratorI vitt;

    UT_Vector3i voxelRes;
    if (level == 0)
        liquidSurface.getVoxelRes(voxelRes[0], voxelRes[1], voxelRes[2]);
    else
	voxelRes = octreeLabels.getVoxelRes(level);

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        if (!vit.isTileConstant())
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                UT_Vector3i edge(vitt.x(), vitt.y(), vitt.z());

                bool isStressActive = false;

		for (int cellIndex = 0; cellIndex < 4; ++cellIndex)
                {
                    UT_Vector3i cell = HDKedgeToCell(edge, axis, cellIndex);

                    if (cell[0] < 0 || cell[1] < 0 || cell[2] < 0 ||
                        cell[0] >= voxelRes[0] || cell[1] >= voxelRes[1] || cell[2] >= voxelRes[2])
                    {
                        vitt.setValue(HDK_OUTSIDE);
                        break;
                    }

                    // If the grid points down then a descendant edge must be active instead.
                    if (octreeLabels.getCellLabel(cell, level) == HDK_OctreeGrid::DOWN)
                    {
                        isStressActive = false;
                        break;
                    }
                    else if (octreeLabels.getCellLabel(cell, level) == HDK_OctreeGrid::ACTIVE)
                        isStressActive = true;
                }
                if (isStressActive)
                {
                    if (level == 0)
                    {
			if (HDKgetFieldValue(edgeIntegrationWeights, edge) > 0)
                            vitt.setValue(HDK_FLUID);
                        else
                            vitt.setValue(HDK_OUTSIDE);
                    }
                    else
                    {
                        vitt.setValue(HDK_FLUID);

#if !defined(NDEBUG)
                        // The edge must be inside the surface if it is active at a coarse level.
                        UT_Vector3 point;
                        edgeStressIndices.indexToPos(edge[0], edge[1], edge[2], point);
			assert(liquidSurface.getValue(point) <= 0);
#endif
                    }
                }
            }
        }
    }
}

void
HDK_AdaptiveViscosity::classifyCenterStressesPartial(SIM_RawIndexField &centerStressIndices,
							const SIM_RawField &octreeLabels,
							const SIM_RawField &centerIntegrationWeights,
							const int level,
							const UT_JobInfo &info) const
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorF vit;
    vit.setConstArray(octreeLabels.field());
    vit.splitByTile(info);

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

		    if (level != 0 || HDKgetFieldValue(centerIntegrationWeights, cell) > 0.)
			HDKsetFieldValue(centerStressIndices, cell, HDK_FLUID);
                }
            }
        }
    }
}

exint
HDK_AdaptiveViscosity::buildRegularVelocityIndices(SIM_RawIndexField (&regularVelocityIndices)[3],
						    const SIM_RawField &liquidSurface,
						    const SIM_RawField &solidSurface,
						    const SIM_RawField &centerIntegrationWeights,
						    const SIM_RawField (&edgeIntegrationWeights)[3],
						    const fpreal extrapolation) const
{
    // Parallel loop to build regular grid velocity indices"

    UT_Array<bool> isTileOccupiedList;

    for (int axis : {0,1,2})
    {
        regularVelocityIndices[axis].makeConstant(HDK_UNASSIGNED);

        isTileOccupiedList.clear();
        isTileOccupiedList.setSize(regularVelocityIndices[axis].field()->numTiles());
        isTileOccupiedList.constant(false);

        // Uncompress velocity tiles either inside or near the liquid surface
        findOccupiedRegularVelocityTiles(isTileOccupiedList,
					    liquidSurface,
					    regularVelocityIndices[axis],
					    axis);

        uncompressTiles(regularVelocityIndices[axis], isTileOccupiedList);

        // Launch parallel loop to classify velocity faces
	classifyRegularVelocityFaces(regularVelocityIndices[axis],
					liquidSurface,
					solidSurface,
					centerIntegrationWeights,
					edgeIntegrationWeights,
					axis,
					extrapolation);
    }

    // Serial loop to build the actual index labels
    UT_Interrupt *boss = UTgetInterrupt();

    exint velocityIndex = 0;
    for (int axis : {0,1,2})
    {
        UT_VoxelArrayIteratorI vit(regularVelocityIndices[axis].fieldNC());
        vit.setCompressOnExit(true);
        UT_VoxelTileIteratorI vitt;

        for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
        {
            if (boss->opInterrupt())
                break;

            if (!vit.isTileConstant() || vit.getValue() == HDK_FLUID)
            {
                vitt.setTile(vit);

                for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
                {
                    if (vitt.getValue() == HDK_FLUID)
                        vitt.setValue(velocityIndex++);
                }
            }
        }
    }

    return velocityIndex;
}

exint
HDK_AdaptiveViscosity::buildOctreeVelocityIndices(UT_Array<UT_Array<SIM_RawIndexField>> &octreeVelocityIndices,
						    const SIM_RawField &liquidSurface,
						    const SIM_RawField &solidSurface,
						    const HDK_OctreeGrid &octreeLabels,
						    const SIM_RawField &centerIntegrationWeights,
						    const SIM_RawField (&edgeIntegrationWeights)[3],
						    const fpreal extrapolation) const
{
    const int octreeLevels = octreeLabels.getOctreeLevels();

    // Run parallel loops to build the flags for active velocities
    // then clean up with a sweep to set the index.

    UT_Array<bool> isTileOccupiedList;
    for (int level = 0; level < octreeLevels; ++level)
	for (int axis : {0,1,2})
        {
            octreeVelocityIndices[level][axis].makeConstant(HDK_UNASSIGNED);

            isTileOccupiedList.clear();
            isTileOccupiedList.setSize(octreeVelocityIndices[level][axis].field()->numTiles());
            isTileOccupiedList.constant(false);

            // Uncompress velocity tiles either inside or near the liquid surface
            if (level == 0)
                findOccupiedRegularVelocityTiles(isTileOccupiedList,
						    liquidSurface,
						    octreeVelocityIndices[level][axis],
						    axis);
            else
                findOccupiedOctreeVelocityTiles(isTileOccupiedList,
						octreeLabels.getGridLabels(level),
                                                octreeVelocityIndices[level][axis],
						axis);

            uncompressTiles(octreeVelocityIndices[level][axis], isTileOccupiedList);

            // Launch parallel loop to classify velocity faces
            classifyOctreeVelocityFaces(octreeVelocityIndices[level][axis],
					octreeLabels,
					liquidSurface,
					solidSurface,
					centerIntegrationWeights,
					edgeIntegrationWeights,
					axis, level,
					extrapolation);
        }
    
    // Serial loop to build the actual index labels
    UT_Interrupt *boss = UTgetInterrupt();

    exint velocityIndex = 0;

    for (int level = 0; level < octreeLevels; ++level)
        for (int axis : {0,1,2})
        {
            UT_VoxelArrayIteratorI vit(octreeVelocityIndices[level][axis].fieldNC());
	    vit.setCompressOnExit(true);
            UT_VoxelTileIteratorI vitt;

            for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
            {
                if (boss->opInterrupt())
                    break;

                if (!vit.isTileConstant() || vit.getValue() == HDK_FLUID)
                {
                    vitt.setTile(vit);

                    for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
                    {
                        if (vitt.getValue() == HDK_FLUID)
                            vitt.setValue(velocityIndex++);
                    }
                }
            }
        }

    return velocityIndex;
}

exint
HDK_AdaptiveViscosity::buildEdgeStressIndices(UT_Array<UT_Array<SIM_RawIndexField>> &edgeStressIndices,
						const SIM_RawField &liquidSurface,
						const HDK_OctreeGrid &octreeLabels,
						const SIM_RawField (&edgeIntegrationWeights)[3]) const
{
    const int octreeLevels = octreeLabels.getOctreeLevels();

    // Parallel loop to build edge indices

    UT_Array<bool> isTileOccupiedList;
    for (int level = 0; level < octreeLevels; ++level)
        for (int axis : {0,1,2})
        {
            edgeStressIndices[level][axis].makeConstant(HDK_UNASSIGNED);

            isTileOccupiedList.clear();
            isTileOccupiedList.setSize(edgeStressIndices[level][axis].field()->numTiles());
            isTileOccupiedList.constant(false);

            // Uncompress tiles for edges on around each active cell.
            findOccupiedEdgeStressTiles(isTileOccupiedList,
					octreeLabels.getGridLabels(level),
					edgeStressIndices[level][axis],
					axis);

            uncompressTiles(edgeStressIndices[level][axis], isTileOccupiedList);

            // Launch parallel loop to classify edge stress
            classifyEdgeStresses(edgeStressIndices[level][axis],
				    octreeLabels,
				    liquidSurface,
				    edgeIntegrationWeights[axis],
				    axis, level);
        }

    // Serial loop to build the actual index labels
    UT_Interrupt *boss = UTgetInterrupt();

    exint stressIndex = 0;

    for (int level = 0; level < octreeLevels; ++level)
        for (int axis : {0,1,2})
        {
            UT_VoxelArrayIteratorI vit(edgeStressIndices[level][axis].fieldNC());
	    vit.setCompressOnExit(true);
            UT_VoxelTileIteratorI vitt;

            for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
            {
                if (boss->opInterrupt())
                    break;

                if (!vit.isTileConstant() || vit.getValue() == HDK_FLUID)
                {
                    vitt.setTile(vit);

                    for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
                    {
                        if (vitt.getValue() == HDK_FLUID)
                            vitt.setValue(stressIndex++);
                    }
                }
            }
        }

    return stressIndex;
}

exint
HDK_AdaptiveViscosity::buildCenterStressIndices(UT_Array<SIM_RawIndexField> &centerStressIndices,
						const HDK_OctreeGrid &octreeLabels,
						const SIM_RawField &centerIntegrationWeights) const
{
    const int octreeLevels = octreeLabels.getOctreeLevels();

    // Parallel loop to build cell indices

    for (int level = 0; level < octreeLevels; ++level)
    {
        centerStressIndices[level].makeConstant(HDK_UNASSIGNED);

        // Launch parallel loop to classify center stress
        classifyCenterStresses(centerStressIndices[level],
				octreeLabels.getGridLabels(level),
				centerIntegrationWeights,
				level);
    }

    // Serial loop to build the actual index labels
    UT_Interrupt *boss = UTgetInterrupt();

    exint stressIndex = 0;

    for (int level = 0; level < octreeLevels; ++level)
    {
        UT_VoxelArrayIteratorI vit(centerStressIndices[level].fieldNC());
	vit.setCompressOnExit(true);
        UT_VoxelTileIteratorI vitt;

        for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
        {
            if (boss->opInterrupt())
                break;

            if (!vit.isTileConstant() || vit.getValue() == HDK_FLUID)
            {
                vitt.setTile(vit);

                for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
                {
                    if (vitt.getValue() == HDK_FLUID)
                        vitt.setValue(stressIndex++);
                }
            }
        }
    }

    return stressIndex;
}

void
getEdgeStressFaces(UT_Array<HDK_AdaptiveViscosity::StressStencilFace> &stencilFaces,
		    UT_Array<fpreal> &boundaryFaces,
		    const UT_Vector3i &edge,
		    const UT_Array<UT_Array<SIM_RawIndexField>> &octreeVelocityIndices,
		    const UT_Array<UT_Array<SIM_RawIndexField>> &edgeStressIndices,
		    const HDK_OctreeGrid &octreeLabels,
		    const SIM_VectorField &solidVelocity,
		    const int axis, const int level,
		    const bool useEnhancedGradients)
{
    assert(HDKgetFieldValue(edgeStressIndices[level][axis], edge) >= 0);

    stencilFaces.clear();
    boundaryFaces.clear();

    fpreal dx = edgeStressIndices[level][axis].getVoxelSize().maxComponent();

    bool isAtTransition[3] = {false};
    bool isFaceOutside[3] = {false};
    
    UT_Vector3 gradientDx(0.);

    for (int faceAxis : {0,1,2})
    {
	if (faceAxis == axis)
	    continue;

	for (int direction : {0,1})
	{
	    UT_Vector3i face = HDKedgeToFace(edge, axis, faceAxis, direction);

	    UT_Vector3i faceRes;
	    faceRes[0] = octreeVelocityIndices[level][faceAxis].getXRes();
	    faceRes[1] = octreeVelocityIndices[level][faceAxis].getYRes();
	    faceRes[2] = octreeVelocityIndices[level][faceAxis].getZRes();

	    const int gradientAxis = 3 - faceAxis - axis;

	    // Since edges on the boundary are included, we need to account for
	    // ghost faces out of bounds.
	    if (face[gradientAxis] < 0 || face[gradientAxis] >= faceRes[gradientAxis])
	    {
		gradientDx[gradientAxis] += .5 * dx;
		isFaceOutside[gradientAxis] = true;
	    }
	    else
	    {
		const exint velocityIndex = HDKgetFieldValue(octreeVelocityIndices[level][faceAxis], face);

		if (velocityIndex >= 0)
		    gradientDx[gradientAxis] += .5 * dx;
		else if (velocityIndex == HDK_OUTSIDE || velocityIndex == HDK_SOLIDBOUNDARY)
		{
		    gradientDx[gradientAxis] += .5 * dx;
		    isFaceOutside[gradientAxis] = true;
		}
		// The edge must be at the finest level of the adjacent grids.
		// Face grading makes it so that we can assume only
		// two levels are adjacent to the edge.
		else if (velocityIndex == HDK_UNASSIGNED)
		{
		    gradientDx[gradientAxis] += dx;
		    if (useEnhancedGradients)
			isAtTransition[gradientAxis] = true;
		}
		else
		    assert(false);
	    }
        }
    }

    for (int faceAxis : {0,1,2})
    {
	if (faceAxis == axis)
	    continue;

	for (int direction : {0,1})
	{
	    UT_Vector3i face = HDKedgeToFace(edge, axis, faceAxis, direction);

	     UT_Vector3i faceRes;
	    faceRes[0] = octreeVelocityIndices[level][faceAxis].getXRes();
	    faceRes[1] = octreeVelocityIndices[level][faceAxis].getYRes();
	    faceRes[2] = octreeVelocityIndices[level][faceAxis].getZRes();

	    const int gradientAxis = 3 - faceAxis - axis;
	    fpreal sign = (direction == 0) ? -1 : 1;

	    if (face[gradientAxis] < 0 || face[gradientAxis] >= faceRes[gradientAxis])
		continue;
	    const exint velocityIndex = HDKgetFieldValue(octreeVelocityIndices[level][faceAxis], face);

	    // If the face is set, life is easy. Add to gradient.
	    if (velocityIndex >= 0)
	    {
		// If this edge is at a transition AND this is the direction of the small faces
		if (isAtTransition[gradientAxis] && !isFaceOutside[gradientAxis])
		{
		    UT_Vector3i siblingFace = face;
		    siblingFace[axis] += (edge[axis] % 2 == 0) ? 1 : -1;

		    const exint siblingVelocityIndex = HDKgetFieldValue(octreeVelocityIndices[level][faceAxis], siblingFace);
		    assert(siblingVelocityIndex >= 0);

		    stencilFaces.append(HDK_AdaptiveViscosity::StressStencilFace(siblingVelocityIndex, .25 * sign / gradientDx[gradientAxis]));
		    stencilFaces.append(HDK_AdaptiveViscosity::StressStencilFace(velocityIndex, .25 * sign / gradientDx[gradientAxis]));
		}
		else
		    // We use .5 because the D operator has .5 coeffs for the off-diagonal terms
		    stencilFaces.append(HDK_AdaptiveViscosity::StressStencilFace(velocityIndex, .5 * sign / gradientDx[gradientAxis]));
	    }
	    // Face is inactive. We're possibly at a dangling node or an adjacent big face.
	    else if (velocityIndex == HDK_UNASSIGNED)
	    {
		// If the adjacent face is inactive, we need to handle the
		// possibility of a dangling edge. Dangling edges have "odd" index
		// positions in the non-edge aligned axis.
		if (edge[faceAxis] % 2 != 0)
		{
		    // The dangling edge case has a few cases that need to be
		    // handled. If the two axis-aligned faces in the neighbouring
		    // cell are one level higher, we just average them together in
		    // the sparse matrix. However if one is small, then the small
		    // faces must be averaged together to match the big face's
		    // position and then added to the matrix.

		    for (int offset : {-1,1})
		    {
			UT_Vector3i offsetFace = face;
			offsetFace[faceAxis] += offset;

			UT_Vector3i parentFace = octreeLabels.getParentFace(offsetFace);
			exint parentVelocityIndex = HDKgetFieldValue(octreeVelocityIndices[level + 1][faceAxis], parentFace);

			// We average the two big faces together for the dangling edge
			if (parentVelocityIndex >= 0)
			    // We use .25 because it is an average of two faces with a .5 coefficient
			    stencilFaces.append(HDK_AdaptiveViscosity::StressStencilFace(parentVelocityIndex, .25 * sign / gradientDx[gradientAxis]));
			// If the small inset faces are active, we will want to accept all of them
			else if (parentVelocityIndex == HDK_UNASSIGNED)
			{
			    for (int childIndex = 0; childIndex < 4; ++childIndex)
			    {
				// The way the face grids are set up is a little
				// counter intuitive since the lower level faces
				// must be inserts in the big face.
				UT_Vector3i childFace = octreeLabels.getChildFace(parentFace, faceAxis, childIndex);

				// Since we went up a level to get the parent and
				// then got that child, we're back at our starting level.
				const exint childVelocityIndex = HDKgetFieldValue(octreeVelocityIndices[level][faceAxis], childFace);

				if (childVelocityIndex >= 0)
				{
				    // We're averaging the four small faces together
				    // at the big face center, which is then
				    // averaged again at the cell center.
				    stencilFaces.append(HDK_AdaptiveViscosity::StressStencilFace(childVelocityIndex, .0625 * sign / gradientDx[gradientAxis]));
				}
				else
				    assert(false);
			    }
			}
			else
			    assert(parentVelocityIndex != HDK_SOLIDBOUNDARY);
		    }
		}
		// The adjacent face is at the parent level. It can added directly.
		else
		{
		    UT_Vector3i parentFace = octreeLabels.getParentFace(face);
		    
		    exint parentVelocityIndex = HDKgetFieldValue(octreeVelocityIndices[level + 1][faceAxis], parentFace);
		    assert(parentVelocityIndex >= 0);

		    stencilFaces.append(HDK_AdaptiveViscosity::StressStencilFace(parentVelocityIndex, .5 * sign / gradientDx[gradientAxis]));
		}
	    }
	    else if (velocityIndex == HDK_SOLIDBOUNDARY)
	    {
		UT_Vector3 point;
		octreeVelocityIndices[level][faceAxis].indexToPos(face[0], face[1], face[2], point);

		fpreal localVelocity = solidVelocity.getField(axis)->getValue(point);
		boundaryFaces.append(.5 * sign * localVelocity / gradientDx[gradientAxis]);

		assert(level == 0);
	    }
	}
    }
}

void
getCenterStressFaces(UT_Array<HDK_AdaptiveViscosity::StressStencilFace> &stencilFaces,
			UT_Array<fpreal> &boundaryFaces,
			const UT_Vector3i &cell,
			const UT_Array<UT_Array<SIM_RawIndexField>> &octreeVelocityIndices,
			const HDK_OctreeGrid &octreeLabels,
			const SIM_VectorField &solidVelocity,
			const int axis, const int level)
{
    assert(octreeLabels.isCellActive(cell, level));
    stencilFaces.clear();
    boundaryFaces.clear();

    fpreal dx = octreeLabels.getVoxelSize(level).maxComponent();

    for (int direction : {0,1})
    {
        UT_Vector3i face = HDKcellToFace(cell, axis, direction);

        fpreal sign = (direction == 0) ? -1 : 1;

	const exint velocityIndex = HDKgetFieldValue(octreeVelocityIndices[level][axis], face);

        if (velocityIndex >= 0)
            stencilFaces.append(HDK_AdaptiveViscosity::StressStencilFace(velocityIndex, sign / dx));
        // If the adjacent face is not active, the child faces must be active.
        // we average them together to use in our velocity gradient.
        else if (velocityIndex == HDK_UNASSIGNED)
        {
            for (int childIndex = 0; childIndex < 4; ++childIndex)
            {
		assert(level > 0);

		UT_Vector3i childFace = octreeLabels.getChildFace(face, axis, childIndex);

		exint childVelocityIndex = HDKgetFieldValue(octreeVelocityIndices[level - 1][axis], childFace);
                assert(childVelocityIndex >= 0);

                // Note the coefficient of .25 because we're averaging the small faces
                stencilFaces.append(HDK_AdaptiveViscosity::StressStencilFace(childVelocityIndex, .25 * sign / dx));
            }
        }
	else if (velocityIndex == HDK_SOLIDBOUNDARY)
        {
            UT_Vector3 point;
            octreeVelocityIndices[level][axis].indexToPos(face[0], face[1], face[2], point);

	    fpreal localVelocity = solidVelocity.getField(axis)->getValue(point);
            boundaryFaces.append(sign * localVelocity / dx);

            assert(level == 0);
        }
    }
}

fpreal
faceOctreeVolumes(const UT_Vector3i &face,
		    const UT_Array<UT_Array<SIM_RawIndexField>> &octreeVelocityIndices,
		    const HDK_OctreeGrid &octreeLabels,
		    const int axis, const int level)
{
    const UT_Vector3i voxelRes = octreeLabels.getVoxelRes(level);

    // Edge-aligned spacing
    fpreal dx = fpreal(1 << level);
    fpreal gradientDx = 0;

    // Find gradient dx. There are no safety checks here. It should be verified
    // that these gradient samples are active
    for (int direction : {0,1})
    {
        const UT_Vector3i cell = HDKfaceToCell(face, axis, direction);

        if (cell[axis] < 0 || cell[axis] >= voxelRes[axis])
            gradientDx += .5 * dx;
        else
        {
            if (octreeLabels.getCellLabel(cell, level) == HDK_OctreeGrid::ACTIVE ||
                octreeLabels.getCellLabel(cell, level) == HDK_OctreeGrid::INACTIVE)
                gradientDx += .5 * dx;
            else
            {
                const UT_Vector3i parentCell = octreeLabels.getParentCell(cell);
                if (octreeLabels.isCellActive(parentCell, level + 1))
                    gradientDx += dx;
                else
                    assert(false);
            }
        }
    }

    return dx * dx * gradientDx;
}

fpreal
edgeOctreeVolumes(const UT_Vector3i &edge,
		    const UT_Array<UT_Array<SIM_RawIndexField>> &octreeVelocityIndices,
		    const UT_Array<UT_Array<SIM_RawIndexField>> &edgeStressIndices,
		    const HDK_OctreeGrid &octreeLabels,
		    const int axis, const int level)
{
    assert(HDKgetFieldValue(edgeStressIndices[level][axis], edge) >= 0);

    // Build the dx per gradient axis
    fpreal dx = fpreal(1 << level);

    UT_Vector3 volumeDx(0.);
    volumeDx[axis] = dx;

     for (int faceAxis : {0,1,2})
    {
	if (faceAxis == axis)
	    continue;
	for (int direction : {0,1})
	{
	    UT_Vector3i face = HDKedgeToFace(edge, axis, faceAxis, direction);

	    int gradientAxis = 3 - faceAxis - axis;

	    UT_Vector3i faceRes;
	    faceRes[0] = octreeVelocityIndices[level][faceAxis].getXRes();
	    faceRes[1] = octreeVelocityIndices[level][faceAxis].getYRes();
	    faceRes[2] = octreeVelocityIndices[level][faceAxis].getZRes();

	    if (face[gradientAxis] < 0 || face[gradientAxis] >= faceRes[gradientAxis])
		volumeDx[gradientAxis] += .5 * dx;
	    else
	    {
		exint velocityIndex = HDKgetFieldValue(octreeVelocityIndices[level][faceAxis], face);

		if (velocityIndex >= 0 ||
		    velocityIndex == HDK_OUTSIDE ||
		    velocityIndex == HDK_SOLIDBOUNDARY)
			volumeDx[gradientAxis] += .5 * dx;
		// This assumes that the parent face is indeed active
		// because if it isn't, then we have a huge problem.
		// If the face is inactive, it's possible we're at a dangling edge.
		// If so, we could have a overlapping volumes. We need to check.
		else if (velocityIndex == HDK_UNASSIGNED)
		    volumeDx[gradientAxis] += dx;
		else
		    assert(false);
	    }
	}
    }

    return volumeDx[0] * volumeDx[1] * volumeDx[2];
}

void
HDK_AdaptiveViscosity::buildEdgeStressStencilsPartial(UT_Array<UT_Array<StressStencilFace>> &edgeStressStencils,
							UT_Array<UT_Array<fpreal>> &edgeStressBoundaryStencils,
							UT_Array<fpreal> &edgeStressWeights,
							const edgeStressStencilParms &parms,
							const int axis, const int level,
							const UT_JobInfo &info) const
{
    UT_Interrupt *boss = UTgetInterrupt();

    const UT_Array<UT_Array<SIM_RawIndexField>> &octreeVelocityIndices = parms.myOctreeVelocityIndices;
    const UT_Array<UT_Array<SIM_RawIndexField>> &edgeStressIndices = parms.myEdgeStressIndices;

    const HDK_OctreeGrid &octreeLabels = parms.myOctreeLabels;
    const SIM_VectorField &solidVelocity = parms.mySolidVelocity;

    const SIM_RawField(&edgeIntegrationWeights)[3] = parms.myEdgeIntegrationWeights;

    const SIM_RawField &viscosity = parms.myViscosity;

    const bool useEnhancedGradients = parms.myUseEnhancedGradients;

    const fpreal dt = parms.myDt;

    UT_VoxelArrayIteratorI vit;
    vit.setConstArray(edgeStressIndices[level][axis].field());
    vit.splitByTile(info);

    UT_VoxelTileIteratorI vitt;

    fpreal32 constantViscosity = 0.;
    bool isViscosityConsant = viscosity.field()->isConstant(&constantViscosity);

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        // If the tile is constant, it can't possibly be active
        // since the indices should be unique per active edge.
        if (!vit.isTileConstant())
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                exint edgeStressIndex = vitt.getValue();

                if (edgeStressIndex >= 0)
                {
                    auto &stencilFaces = edgeStressStencils[edgeStressIndex];
                    auto &boundaryStencilFaces = edgeStressBoundaryStencils[edgeStressIndex];

                    UT_Vector3i edge(vitt.x(), vitt.y(), vitt.z());

                    getEdgeStressFaces(stencilFaces,
					boundaryStencilFaces,
					edge,
					octreeVelocityIndices,
					edgeStressIndices,
					octreeLabels,
					solidVelocity,
					axis, level,
					useEnhancedGradients);

		    fpreal localEdgeStressWeight;
                    if (level == 0)
                    {
			localEdgeStressWeight = HDKgetFieldValue(edgeIntegrationWeights[axis], edge);

                        // If the control volume is full, it is possible that
                        // the edge is at a transition. In this case, 1. is not
                        // accurate as the volume will be stretched into the
                        // coarse cells.
                        if (localEdgeStressWeight == 1.)
                            localEdgeStressWeight = edgeOctreeVolumes(edge, octreeVelocityIndices,
									edgeStressIndices,
									octreeLabels,
									axis, level);
                    }
                    else
                        localEdgeStressWeight = edgeOctreeVolumes(edge, octreeVelocityIndices,
								    edgeStressIndices,
								    octreeLabels,
								    axis, level);
		    if (isViscosityConsant)
			localEdgeStressWeight *= constantViscosity;
                    else
                    {
                        UT_Vector3 point;
                        edgeStressIndices[level][axis].indexToPos(edge[0], edge[1], edge[2], point);
                        localEdgeStressWeight *= viscosity.getValue(point);
                    }
                    
                    // Scale by 4 here to account for the "2" entry from K and
                    // the "2" leading coefficient in the system.
                    edgeStressWeights[edgeStressIndex] = 4. * dt * localEdgeStressWeight;
                }
            }
        }
    }
}

void
HDK_AdaptiveViscosity::buildCenterStressStencilsPartial(UT_Array<UT_Array<StressStencilFace>> &centerStressStencils,
							UT_Array<UT_Array<fpreal>> &centerStressBoundaryStencils,
							const centerStressStencilParms &parms,
							const int axis, const int level,
							const exint centerStressCount, const UT_JobInfo &info) const
{
    UT_Interrupt *boss = UTgetInterrupt();

    const UT_Array<UT_Array<SIM_RawIndexField>> &octreeVelocityIndices = parms.myOctreeVelocityIndices;
    const UT_Array<SIM_RawIndexField> &centerStressIndices = parms.myCenterStressIndices;

    const HDK_OctreeGrid &octreeLabels = parms.myOctreeLabels;

    const SIM_VectorField &solidVelocity = parms.mySolidVelocity;

    UT_VoxelArrayIteratorI vit;
    vit.setConstArray(centerStressIndices[level].field());
    vit.splitByTile(info);

    UT_VoxelTileIteratorI vitt;

    // We only assigned one index per cell. Therefore each
    // stress stencil axis must be indexed according to an offset.
    const exint offset = centerStressCount * axis;

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        // If the tile is constant, it can't possibly be active
        // since the indices should be unique per active cell.
        if (!vit.isTileConstant())
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                exint cellIndex = vitt.getValue();

                if (cellIndex >= 0)
                {
                    UT_Vector3i cell(vitt.x(), vitt.y(), vitt.z());

                    auto &stencilFaces = centerStressStencils[cellIndex + offset];
                    auto &boundaryFaces = centerStressBoundaryStencils[cellIndex + offset];

                    getCenterStressFaces(stencilFaces,
					    boundaryFaces,
					    cell,
					    octreeVelocityIndices,
					    octreeLabels,
					    solidVelocity,
					    axis, level);
                }
            }
        }
    }
}

void
HDK_AdaptiveViscosity::buildCenterStressWeightsPartial(UT_Array<fpreal> &centerStressWeights,
							const SIM_RawIndexField &centerStressIndices,
							const SIM_RawField &centerIntegrationWeights,
							const SIM_RawField &viscosity,
							const fpreal dt,
							const int level,
							const UT_JobInfo &info) const
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorI vit;
    vit.setConstArray(centerStressIndices.field());
    vit.splitByTile(info);

    UT_VoxelTileIteratorI vitt;

    fpreal cellVolume;
    if (level > 0)
    {
        fpreal dx = fpreal(1 << level);
        cellVolume = dx * dx * dx;
    }

    fpreal32 constantViscosity = 0.;
    bool isViscosityConsant = viscosity.field()->isConstant(&constantViscosity);

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        // If the tile is constant, it can't possibly be active
        // since the indices should be unique per active cell.
        if (!vit.isTileConstant())
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                exint cellIndex = vitt.getValue();

                if (cellIndex >= 0)
                {
                    UT_Vector3i cell(vitt.x(), vitt.y(), vitt.z());

                    fpreal localCenterStressWeight;

                    if (level == 0)
			localCenterStressWeight = HDKgetFieldValue(centerIntegrationWeights, cell);
                    else localCenterStressWeight = cellVolume;

		    if (isViscosityConsant)
			localCenterStressWeight *= constantViscosity;
		    else
		    {
			UT_Vector3 point;
                        centerStressIndices.indexToPos(cell[0], cell[1], cell[2], point);
                        localCenterStressWeight *= viscosity.getValue(point);
		    }

                    centerStressWeights[cellIndex] = 2. * dt * localCenterStressWeight;
                }
            }
        }
    }
}

void
HDK_AdaptiveViscosity::buildVelocityMappingPartial(Vector &initialGuess,
						    const SIM_RawField &regularVelocity,
						    const SIM_RawIndexField &regularVelocityIndices,
						    const SIM_RawIndexField &octreeVelocityIndices,
						    const SIM_VectorField &solidVelocity,
						    const HDK_OctreeGrid &octreeLabels,
						    const int axis, const int level,
						    const UT_JobInfo &info) const
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorI vit;
    vit.setConstArray(octreeVelocityIndices.field());
    vit.splitByTile(info);

    UT_VoxelTileIteratorI vitt;

    struct RestrictionFaceAndWeights
    {
	RestrictionFaceAndWeights() {}
	RestrictionFaceAndWeights(const UT_Vector3i &face, fpreal32 weight, int level = 0)
	    : myFace(face), myWeight(weight), myLevel(level)
	{
	}

	UT_Vector3i myFace;
	fpreal32 myWeight;
	int myLevel;
    };

    // Weights are used to build the upwards interpolation scheme to sample from the lower level cells.
    const fpreal inAxisWeights[3] = {1. / 16., 1. / 8., 1. / 16.};

    std::queue<RestrictionFaceAndWeights> restrictionFaceQueue;

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        // If the tile is constant, it can't possibly be active
        // since the indices should be unique per active face.
        if (!vit.isTileConstant())
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                const exint octreeFaceIndex = vitt.getValue();
                if (octreeFaceIndex >= 0)
                {
                    fpreal restrictedVelocity = 0;
                    UT_Vector3i face(vitt.x(), vitt.y(), vitt.z());

                    assert(restrictionFaceQueue.empty());
                    restrictionFaceQueue.push(RestrictionFaceAndWeights(face, 1., level));

                    while (!restrictionFaceQueue.empty())
                    {
                        // If the current face's level is zero, add it to the
                        // vector. If not, add the children to the queue.

                        const UT_Vector3i localFace = restrictionFaceQueue.front().myFace;
                        const fpreal localWeight = restrictionFaceQueue.front().myWeight;
                        const int localLevel = restrictionFaceQueue.front().myLevel;

                        restrictionFaceQueue.pop();

                        if (localLevel == 0)
                        {
			    exint regularFaceIndex = HDKgetFieldValue(regularVelocityIndices, localFace);
                            if (regularFaceIndex < 0)
                            {
				assert(regularFaceIndex == HDK_SOLIDBOUNDARY);
                                assert(level == 0);
                            }

                            // Even if the regular grid index is a collision, we
                            // should still get values from the velocity field
                            // since that is what we would expect to see in the
                            // uniform grid.
			    restrictedVelocity += localWeight * HDKgetFieldValue(regularVelocity, localFace);
                        }
                        else
                        {
                            // Loop over the 3x4x4 (for x-axis) stencil of adjacent lower-level faces.
			    for (int childIndex = 0; childIndex < 4; ++childIndex)
                            {
				UT_Vector3i childFace = octreeLabels.getChildFace(localFace, axis, childIndex);

                                for (int inAxisOffset = -1; inAxisOffset < 2; ++inAxisOffset)
                                {
                                    UT_Vector3i adjacentFace = childFace;
                                    adjacentFace[axis] += inAxisOffset;

                                    fpreal restrictionWeight = inAxisWeights[inAxisOffset + 1];

				    restrictionFaceQueue.push(RestrictionFaceAndWeights(adjacentFace,
										    restrictionWeight * localWeight,
										    localLevel - 1));
                                }
                            }
                        }
                    }

                    initialGuess(octreeFaceIndex) = restrictedVelocity;
                }
            }
        }
    }
}

SYS_FORCE_INLINE void
applyToMatrix(
#ifdef USEEIGEN
		std::vector<Eigen::Triplet<SolveType>> &sparseMatrixElements,
#else
		SparseMatrix &sparseMatrixElements,
#endif
		Vector &rhs,
		fpreal &diagonalElement,
		fpreal coefficient,
		const exint velocityIndex,
		const UT_Array<HDK_AdaptiveViscosity::StressStencilFace> &stencilFaces,
		const UT_Array<fpreal> &boundaryStencilFaces)
{
#if !defined(NDEBUG)
    bool foundSelf = false;
#endif
    int faceCount = stencilFaces.entries();

    for (int i = 0; i < faceCount; ++i)
    {
        if (stencilFaces[i].myFaceIndex == velocityIndex)
        {
            coefficient *= stencilFaces[i].myCoefficient;

#if !defined(NDEBUG)
            foundSelf = true;
#endif

            break;
        }
    }
    assert(foundSelf);

    // Add edge stress to viscosity stencil
    for (int i = 0; i < faceCount; ++i)
    {
        fpreal element = coefficient * stencilFaces[i].myCoefficient;

        if (stencilFaces[i].myFaceIndex == velocityIndex)
            diagonalElement += element;
        else
#ifdef USEEIGEN
	    sparseMatrixElements.push_back(Eigen::Triplet<SolveType>(velocityIndex, stencilFaces[i].myFaceIndex, element));
#else
            sparseMatrixElements.addToElement(velocityIndex, stencilFaces[i].myFaceIndex, element);
#endif
    }

    int boundaryFaceCount = boundaryStencilFaces.entries();

    for (int i = 0; i < boundaryFaceCount; ++i)
	rhs(velocityIndex) -= coefficient * boundaryStencilFaces[i];
}

void
HDK_AdaptiveViscosity::buildOctreeSystemFromStencilsPartial(
#ifdef USEEIGEN
							    std::vector<std::vector<Eigen::Triplet<SolveType>>> &parallelSparseMatrixElements,
#else
							    UT_Array<SparseMatrix> &parallelSparseMatrixElements,
#endif
							    Vector &rhs,
							    const octreeSystemStencilParms &parms,
							    const int axis, const int level,
							    const exint centerStressCount,
							    const UT_JobInfo &info) const
{
    UT_Interrupt *boss = UTgetInterrupt();

    const UT_Array<UT_Array<StressStencilFace>> &edgeStressStencils = parms.myEdgeStressStencils;
    const UT_Array<UT_Array<fpreal>> &edgeBoundaryStressStencils = parms.myEdgeBoundaryStressStencils;
    const UT_Array<fpreal> &edgeStressStencilWeights = parms.myEdgeStressStencilWeights;

    const UT_Array<UT_Array<StressStencilFace>> &centerStressStencils = parms.myCenterStressStencils;
    const UT_Array<UT_Array<fpreal>> &centerBoundaryStressStencils = parms.myCenterBoundaryStressstencils;
    const UT_Array<fpreal> &centerStressStencilWeights = parms.myCenterStressStencilWeights;

    const UT_Array<UT_Array<SIM_RawIndexField>> &octreeVelocityIndices = parms.myOctreeVelocityIndices;
    const UT_Array<UT_Array<SIM_RawIndexField>> &edgeStressIndices = parms.myEdgeStressIndices;
    const UT_Array<SIM_RawIndexField> &centerStressIndices = parms.myCenterStressIndices;

#if !defined(NDEBUG)
    const SIM_RawField &centerIntegrationWeights = parms.myCenterIntegrationWeights;
#endif

    const SIM_VectorField &faceIntegrationWeights = parms.myFaceIntegrationWeights;

    const HDK_OctreeGrid &octreeLabels = parms.myOctreeLabels;

    const Vector &initialGuess = parms.myInitialGuess;

    const bool useEnhancedGradients = parms.myUseEnhancedGradients;

    const SIM_RawField &density = parms.myDensity;

    fpreal32 constantDensity;
    bool isDensityConstant = density.field()->isConstant(&constantDensity);

    UT_VoxelArrayIteratorI vit;
    vit.setConstArray(octreeVelocityIndices[level][axis].field());
    vit.splitByTile(info);

    UT_VoxelTileIteratorI vitt;

#ifdef USEEIGEN
    std::vector<Eigen::Triplet<SolveType>> &localSparseMatrixElements = parallelSparseMatrixElements[info.job()];
#else
    SparseMatrix &localSparseMatrixElements = parallelSparseMatrixElements[info.job()];
#endif

    UT_Vector3i voxelRes = octreeLabels.getVoxelRes(level);

    UT_Vector3i faceRes;
    faceRes[0] = octreeVelocityIndices[level][axis].getXRes();
    faceRes[1] = octreeVelocityIndices[level][axis].getYRes();
    faceRes[2] = octreeVelocityIndices[level][axis].getZRes();

    // We only assigned one index per cell. Therefore each
    // stress stencil axis must be indexed according to an offset.
    const exint centerStressOffset = centerStressCount * axis;

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        // If the tile is constant, it can't possibly be active
        // since the indices should be unique per active face.
        if (!vit.isTileConstant())
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
		const exint velocityIndex = vitt.getValue();

                if (velocityIndex >= 0)
                {
                    fpreal diagonalElement = 0;
		    UT_Vector3i face(vitt.x(), vitt.y(), vitt.z());

                    // Build cell centered-stress and dangling edge stresses
		    for (int direction : {0,1})
		    {
			const UT_Vector3i cell = HDKfaceToCell(face, axis, direction);

                        if (cell[axis] < 0 || cell[axis] >= voxelRes[axis])
			    continue;
			
			UT_Vector3i stressCell;
			exint stressLevel;

                        if (octreeLabels.isCellActive(cell, level))
                        {
                            stressCell = cell;
                            stressLevel = level;
                        }
                        // Since we're face graded, we can safely assume that if
                        // the cell isn't active then its parent must be.
                        else
                        {
                            stressCell = octreeLabels.getParentCell(cell);
                            stressLevel = level + 1;

                            assert(octreeLabels.getCellLabel(cell, level) == HDK_OctreeGrid::UP);
                            assert(octreeLabels.isCellActive(stressCell, stressLevel));
                            assert(stressLevel < octreeLabels.getOctreeLevels());
                        }

                        // Apply cell stress contribution to matrx A and rhs
                        {
			    const exint cellIndex = HDKgetFieldValue(centerStressIndices[stressLevel], stressCell);

                            if (cellIndex >= 0)
                            {
                                const auto &centerStressFaces = centerStressStencils[cellIndex + centerStressOffset];
                                assert(centerStressFaces.entries() > 0);

				const auto &centerBoundaryFaces = centerBoundaryStressStencils[cellIndex + centerStressOffset];
				
				if (centerBoundaryFaces.entries() > 0)
				    assert(level == 0);
			    
				// The control weight is obviously independent
                                // of axis and therefore doesn't need an offset.
                                fpreal coefficient = centerStressStencilWeights[cellIndex];

                                applyToMatrix(localSparseMatrixElements,
						rhs,
						diagonalElement,
						coefficient,
						velocityIndex,
						centerStressFaces,
						centerBoundaryFaces);
			    }
#if !defined(NDEBUG)
                            else
                            {
                                assert(stressLevel == 0);
				assert(HDKgetFieldValue(centerIntegrationWeights, stressCell) == 0.);
                            }
#endif
                        }

                        // If we're at transition, edge stresses at T-junctions
                        // will create a ghost cell based on adjacent faces in
                        // the coarse cell. This means velocities at transitions
                        // must account for T-junction edges.

			for (int faceAxis : {0,1,2})
			{
			    if (faceAxis == axis)
				continue;

			    for (int faceDirection : {0,1})
			    {
				UT_Vector3i adjacentFace = HDKcellToFace(stressCell, faceAxis, faceDirection);

				if (HDKgetFieldValue(octreeVelocityIndices[stressLevel][faceAxis], adjacentFace) == HDK_UNASSIGNED)
				{
				    int edgeAxis = 3 - faceAxis - axis;

				    for (int insetEdgeIndex : {0,1})
				    {
					UT_Vector3i edge = octreeLabels.getChildEdgeInFace(adjacentFace, faceAxis, edgeAxis, insetEdgeIndex);

					exint edgeStressIndex = HDKgetFieldValue(edgeStressIndices[stressLevel - 1][edgeAxis], edge);
					if (edgeStressIndex >= 0)
					{
					    const auto &edgeStressFaces = edgeStressStencils[edgeStressIndex];
					    const auto &edgeBoundaryFaces = edgeBoundaryStressStencils[edgeStressIndex];
					    fpreal coefficient = edgeStressStencilWeights[edgeStressIndex];

					    applyToMatrix(localSparseMatrixElements,
							    rhs,
							    diagonalElement,
							    coefficient,
							    velocityIndex,
							    edgeStressFaces,
							    edgeBoundaryFaces);
					}
				    }
				}
			    }
			}			    
                    }

		    for (int edgeAxis : {0,1,2})
		    {
			if (edgeAxis == axis)
			    continue;

			for (int direction : {0,1})
			{
			    UT_Vector3i edge = HDKfaceToEdge(face, axis, edgeAxis, direction);
			    const exint edgeStressIndex = HDKgetFieldValue(edgeStressIndices[level][edgeAxis], edge);

			    if (edgeStressIndex >= 0)
			    {
				if (useEnhancedGradients)
				{
				    const int transitionAxis = 3 - edgeAxis - axis;

				    UT_Vector3i adjacentFace = face;
				    adjacentFace[transitionAxis] += (direction == 0) ? -1 : 1;

				    if (adjacentFace[transitionAxis] >= 0 &&
					adjacentFace[transitionAxis] < faceRes[transitionAxis])
				    {
					if (HDKgetFieldValue(octreeVelocityIndices[level][axis], adjacentFace) == HDK_UNASSIGNED)
					{
					    UT_Vector3i siblingEdge = edge;
					    siblingEdge[edgeAxis] += (edge[edgeAxis] % 2 == 0) ? 1 : -1;

					    const exint transitionEdgeStressIndex = HDKgetFieldValue(edgeStressIndices[level][edgeAxis], siblingEdge);
					    assert(transitionEdgeStressIndex >= 0);

					    const auto &edgeStressFaces = edgeStressStencils[transitionEdgeStressIndex];
					    assert(edgeStressFaces.entries() > 0);

					    const auto &edgeBoundaryFaces = edgeBoundaryStressStencils[transitionEdgeStressIndex];
					    fpreal coefficient = edgeStressStencilWeights[transitionEdgeStressIndex];

					    applyToMatrix(localSparseMatrixElements,
							    rhs,
							    diagonalElement,
							    coefficient,
							    velocityIndex,
							    edgeStressFaces,
							    edgeBoundaryFaces);
					}
				    }
				}
				{
				    const auto &edgeStressFaces = edgeStressStencils[edgeStressIndex];
				    assert(edgeStressFaces.entries() > 0);

				    const auto &edgeBoundaryFaces = edgeBoundaryStressStencils[edgeStressIndex];
				    fpreal coefficient = edgeStressStencilWeights[edgeStressIndex];

				    applyToMatrix(localSparseMatrixElements,
						    rhs,
						    diagonalElement,
						    coefficient,
						    velocityIndex,
						    edgeStressFaces,
						    edgeBoundaryFaces);
				}
			    }
			    else if (edgeStressIndex == HDK_UNASSIGNED)
			    {
				assert(level > 0);

				for (int childIndex = 0; childIndex < 2; ++childIndex)
				{
				    UT_Vector3i childEdge = octreeLabels.getChildEdge(edge, edgeAxis, childIndex);

				    const exint childEdgeStressIndex = HDKgetFieldValue(edgeStressIndices[level - 1][edgeAxis], childEdge);

				    if (childEdgeStressIndex >= 0)
				    {
					const auto &edgeStressFaces = edgeStressStencils[childEdgeStressIndex];
					assert(edgeStressFaces.entries() > 0);

					const auto &edgeBoundaryFaces = edgeBoundaryStressStencils[childEdgeStressIndex];
					fpreal coefficient = edgeStressStencilWeights[childEdgeStressIndex];

					applyToMatrix(localSparseMatrixElements,
							rhs,
							diagonalElement,
							coefficient,
							velocityIndex,
							edgeStressFaces,
							edgeBoundaryFaces);
				    }
				    else assert(childEdgeStressIndex == HDK_OUTSIDE);
				}
			    }
			    else assert(edgeStressIndex == HDK_OUTSIDE);
			}
                    }
                    // Build velocity control weights

                    fpreal faceWeight;
                    if (level == 0)
                    {
			faceWeight = HDKgetFieldValue(*faceIntegrationWeights.getField(axis), face);

                        if (faceWeight == 1.)
			    faceWeight = faceOctreeVolumes(face, octreeVelocityIndices, octreeLabels, axis, level);
                    }
                    else
                        faceWeight = faceOctreeVolumes(face, octreeVelocityIndices, octreeLabels, axis, level);

		    if (isDensityConstant)
			faceWeight *= constantDensity;
		    else
		    {
			UT_Vector3 point;
                        octreeVelocityIndices[level][axis].indexToPos(face[0], face[1], face[2], point);
                        faceWeight *= density.getValue(point);
		    }
#ifdef USEEIGEN
		    localSparseMatrixElements.push_back(Eigen::Triplet<SolveType>(velocityIndex, velocityIndex, faceWeight + diagonalElement));
#else
		    localSparseMatrixElements.addToElement(velocityIndex, velocityIndex, faceWeight + diagonalElement);
#endif
		    rhs(velocityIndex) += faceWeight * initialGuess(velocityIndex);
                }
            }
        }
    }
}

void
HDK_AdaptiveViscosity::setOctreeVelocityPartial(SIM_RawField &octreeVelocity,
						const SIM_RawIndexField &octreeVelocityIndices,
						const Vector &solution,
						const UT_JobInfo &info) const
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorI vit;
    vit.setConstArray(octreeVelocityIndices.field());
    vit.splitByTile(info);

    UT_VoxelTileIteratorI vitt;

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

	if (!vit.isTileConstant())
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
		const exint velocityIndex = vitt.getValue();
		if (velocityIndex >= 0)
		{
		    UT_Vector3i face(vitt.x(), vitt.y(), vitt.z());
		    HDKsetFieldValue(octreeVelocity, face, solution(velocityIndex));
		}
	    }
	}
    }   
}

void
HDK_AdaptiveViscosity::applyVelocitiesToRegularGridPartial(SIM_RawField &regularVelocity,
							    const HDK_OctreeVectorFieldInterpolator &interpolator,
							    const SIM_RawIndexField &octreeVelocityIndices,
							    const SIM_RawIndexField &regularVelocityIndices,
							    const SIM_RawField &solidSurface,
							    const SIM_VectorField &solidVelocity,
							    const Vector &solution,
							    const int axis, const UT_JobInfo &info) const
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorI vit;
    vit.setConstArray(regularVelocityIndices.field());
    vit.splitByTile(info);

    UT_VoxelTileIteratorI vitt;

    UT_Vector3i faceRes;
    faceRes[0] = regularVelocityIndices.getXRes();
    faceRes[1] = regularVelocityIndices.getYRes();
    faceRes[2] = regularVelocityIndices.getZRes();

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

	if (!vit.isTileConstant() || vit.getValue() == HDK_SOLIDBOUNDARY)
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                if (vitt.getValue() >= 0)
                {
                    UT_Vector3i face(vitt.x(), vitt.y(), vitt.z());

                    // Check for matching octree faces at the regular grid level
		    const exint octreeFaceIndex = HDKgetFieldValue(octreeVelocityIndices, face);

                    if (octreeFaceIndex >= 0)
			HDKsetFieldValue(regularVelocity, face, solution(octreeFaceIndex));
                    else
                    {
			if (octreeFaceIndex == HDK_SOLIDBOUNDARY)
                        {
			    UT_Vector3 point;
			    octreeVelocityIndices.indexToPos(face[0], face[1], face[2], point);
			    fpreal localVelocity = solidVelocity.getField(axis)->getValue(point);

                            HDKsetFieldValue(regularVelocity, face, localVelocity);
                        }
                        else if (octreeFaceIndex == HDK_UNASSIGNED)
                        {
			    UT_Vector3 point;
			    regularVelocityIndices.indexToPos(face[0], face[1], face[2], point);
			    fpreal localVelocity = interpolator.interpSPGrid(point, axis);
			    //fpreal localVelocity = interpolator.interp(point, axis);

			    HDKsetFieldValue(regularVelocity, face, localVelocity);
                        }
                        else
                            assert(false);
                    }
                }
		else if (vitt.getValue() == HDK_SOLIDBOUNDARY)
                {
                    UT_Vector3i face(vitt.x(), vitt.y(), vitt.z());

		    UT_Vector3 point;
		    regularVelocityIndices.indexToPos(face[0], face[1], face[2], point);
		    fpreal localVelocity = solidVelocity.getField(axis)->getValue(point);

		    HDKsetFieldValue(regularVelocity, face, localVelocity);
                }
            }
        }
    }
}

void
HDK_AdaptiveViscosity::octreeVelocityGradingUnitTestPartial(bool &passed,
							    const UT_Array<UT_Array<SIM_RawIndexField>> &octreeVelocityIndices,
							    const HDK_OctreeGrid &octreeLabels,
							    const int axis, const int level,
							    const UT_JobInfo &info) const
{
    UT_Interrupt *boss = UTgetInterrupt();

    const int octreeLevels = octreeLabels.getOctreeLevels();

    UT_VoxelArrayIteratorI vit;
    vit.setConstArray(octreeVelocityIndices[level][axis].field());
    vit.splitByTile(info);

    UT_VoxelTileIteratorI vitt;

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        // Tiles for the index grid must have already been uncompressed
        if (!vit.isTileConstant())
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                if (vitt.getValue() >= 0)
                {
                    UT_Vector3i face(vitt.x(), vitt.y(), vitt.z());

		    UT_Vector3i backwardCell = HDKfaceToCell(face, axis, 0);
		    UT_Vector3i forwardCell = HDKfaceToCell(face, axis, 1);

		    auto backwardLabel = octreeLabels.getCellLabel(backwardCell, level);
		    auto forwardLabel = octreeLabels.getCellLabel(forwardCell, level);

                    if (backwardLabel == HDK_OctreeGrid::ACTIVE &&
                        forwardLabel == HDK_OctreeGrid::UP)
                    {
                        UT_Vector3i parentCell = octreeLabels.getParentCell(forwardCell);

                        if (level == octreeLevels - 1 ||
                            octreeLabels.getCellLabel(parentCell, level + 1) != HDK_OctreeGrid::ACTIVE)
                        {
                            passed = false;
                            return;
                        }
                    }
                    else if (backwardLabel == HDK_OctreeGrid::UP &&
                             forwardLabel == HDK_OctreeGrid::ACTIVE)
                    {
			UT_Vector3i parentCell = octreeLabels.getParentCell(backwardCell);

                        if (level == octreeLevels - 1 ||
                            octreeLabels.getCellLabel(parentCell, level + 1) != HDK_OctreeGrid::ACTIVE)
                        {
                            passed = false;
                            return;
                        }
                    }
                    else if (!(backwardLabel == HDK_OctreeGrid::ACTIVE &&
                               forwardLabel == HDK_OctreeGrid::ACTIVE))
                    {
                        passed = false;
                        return;
                    }
                }
                else if (vitt.getValue() == HDK_OUTSIDE || vitt.getValue() == HDK_SOLIDBOUNDARY)
                {
                    if (level != 0)
                    {
                        passed = false;
                        return;
                    }
                }
            }
        }
    }
}

// Verify that velocities are face graded
bool
HDK_AdaptiveViscosity::octreeVelocityUnitTest(const UT_Array<UT_Array<SIM_RawIndexField>> &octreeVelocityIndices,
						const HDK_OctreeGrid &octreeLabels) const
{
    bool passed = true;
    const int octreeLevels = octreeLabels.getOctreeLevels();

    for (int level = 0; level < octreeLevels; ++level)
    {
        for (int axis : {0,1,2})
        {
            octreeVelocityGradingUnitTest(passed, octreeVelocityIndices, octreeLabels, axis, level);

            if (!passed)
                return false;
        }
    }

    return true;
}

void
HDK_AdaptiveViscosity::edgeStressUnitTestPartial(bool &passed,
						    const UT_Array<UT_Array<SIM_RawIndexField>> &edgeStressIndices,
						    const UT_Array<UT_Array<SIM_RawIndexField>> &octreeVelocityIndices,
						    const HDK_OctreeGrid &octreeLabels,
						    const int axis, const int level,
						    const UT_JobInfo &info) const
{
    UT_Interrupt *boss = UTgetInterrupt();

    const int octreeLevels = octreeLabels.getOctreeLevels();

    UT_VoxelArrayIteratorI vit;
    vit.setConstArray(edgeStressIndices[level][axis].field());
    vit.splitByTile(info);

    UT_VoxelTileIteratorI vitt;

    UT_Vector3i voxelRes = octreeLabels.getVoxelRes(level);

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        // Tiles for the index grid must have already been uncompressed
        if (!vit.isTileConstant())
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                if (vitt.getValue() >= 0)
                {
                    UT_Vector3i edge(vitt.x(), vitt.y(), vitt.z());

                    // Check that the adjacent cells are either UP or ACTIVE
                    for (int cellIndex = 0; cellIndex < 4; ++cellIndex)
                    {
			UT_Vector3i cell = HDKedgeToCell(edge, axis, cellIndex);

                        if (cell[0] < 0 || cell[1] < 0 || cell[2] < 0 ||
                            cell[0] >= voxelRes[0] || cell[1] >= voxelRes[1] || cell[2] >= voxelRes[2])
                            break;

                        if (octreeLabels.getCellLabel(cell, level) == HDK_OctreeGrid::DOWN ||
                            octreeLabels.getCellLabel(cell, level) == HDK_OctreeGrid::INACTIVE)
                        {
                            passed = false;
                            return;
                        }
		    }

		    // Check that the adjacent faces are ACTIVE or the
		    // parent is ACTIVE. Check if the parent is a dangling
		    // face and if so, check that the coarse cell it falls
		    // into is ACTIVE.

		    for (int faceAxis : {0,1,2})
		    {
			if (faceAxis == axis)
			    continue;

			for (int direction : {0,1})
			{
			    UT_Vector3i face = HDKedgeToFace(edge, axis, faceAxis, direction);

			    UT_Vector3i faceRes;
			    faceRes[0] = octreeVelocityIndices[level][faceAxis].getXRes();
			    faceRes[1] = octreeVelocityIndices[level][faceAxis].getYRes();
			    faceRes[2] = octreeVelocityIndices[level][faceAxis].getZRes();

			    int offsetAxis = 3 - faceAxis - axis;
			    if (face[offsetAxis] < 0 || face[offsetAxis] >= faceRes[offsetAxis])
				continue;

			    exint velocityIndex = HDKgetFieldValue(octreeVelocityIndices[level][faceAxis], face);

			    if (velocityIndex < 0)
			    {
				if (velocityIndex == HDK_SOLIDBOUNDARY || velocityIndex == HDK_OUTSIDE)
                                {
                                    if (level != 0)
                                    {
                                        passed = false;
                                        return;
                                    }
                                }
				else if (velocityIndex == HDK_UNASSIGNED)
				{
				     // Check for t-junction
				    if (edge[faceAxis] % 2 != 0)
				    {
					UT_Vector3i parentCell = octreeLabels.getParentCell(face);
					if (level == octreeLevels - 1 ||
					    octreeLabels.getCellLabel(parentCell, level + 1) != HDK_OctreeGrid::ACTIVE)
					{
					    passed = false;
					    return;
					}
				    }
				    else
				    {
					UT_Vector3i parentFace = octreeLabels.getParentFace(face);

					if (level == octreeLevels - 1 ||
					    (HDKgetFieldValue(octreeVelocityIndices[level + 1][faceAxis], parentFace) == HDK_UNASSIGNED))
					{
					    passed = false;
					    return;
					}
				    }
				}
				else
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
}

// Verify that edges can only have 2 level faces. The levels can only be
// at the active level or a level higher. Verify that edge is at the lowest
// active level of adjacent cells.
bool
HDK_AdaptiveViscosity::edgeStressUnitTest(const UT_Array<UT_Array<SIM_RawIndexField>> &edgeStressIndices,
					    const UT_Array<UT_Array<SIM_RawIndexField>> &octreeVelocityIndices,
					    const HDK_OctreeGrid &octreeLabels) const
{
    bool passed = true;
    const int octreeLevels = octreeLabels.getOctreeLevels();

    for (int level = 0; level < octreeLevels; ++level)
        for (int axis : {0,1,2})
        {
	    edgeStressUnitTest(passed, edgeStressIndices, octreeVelocityIndices, octreeLabels, axis, level);

            if (!passed)
                return false;
        }

    return true;
}

void
HDK_AdaptiveViscosity::centerStressUnitTestPartial(bool &passed,
						    const UT_Array<SIM_RawIndexField> &centerStressIndices,
						    const UT_Array<UT_Array<SIM_RawIndexField>> &edgeStressIndices,
						    const UT_Array<UT_Array<SIM_RawIndexField>> &octreeVelocityIndices,
						    const HDK_OctreeGrid &octreeLabels,
						    const int level,
						    const UT_JobInfo &info) const
{
    UT_Interrupt *boss = UTgetInterrupt();

    UT_VoxelArrayIteratorI vit;
    vit.setConstArray(centerStressIndices[level].field());
    vit.splitByTile(info);

    UT_VoxelTileIteratorI vitt;

    for (vit.rewind(); !vit.atEnd(); vit.advanceTile())
    {
        if (boss->opInterrupt())
            break;

        // Tiles for the index grid must have already been uncompressed
        if (!vit.isTileConstant())
        {
            vitt.setTile(vit);

            for (vitt.rewind(); !vitt.atEnd(); vitt.advance())
            {
                if (vitt.getValue() >= 0)
                {
                    UT_Vector3i cell(vitt.x(), vitt.y(), vitt.z());

                    if (octreeLabels.getCellLabel(cell, level) != HDK_OctreeGrid::ACTIVE)
                    {
                        passed = false;
                        return;
                    }

		    // Check adjacent faces
		    for (int axis : {0,1,2})
			for (int direction : {0,1})
			{
			    UT_Vector3i face = HDKcellToFace(cell, axis, direction);

			    exint velocityIndex = HDKgetFieldValue(octreeVelocityIndices[level][axis], face);

			    if (velocityIndex == HDK_UNASSIGNED)
			    {
				if (level == 0)
				{
				    passed = false;
				    return;
				}

				// Check that the inset faces are active
				for (int childIndex = 0; childIndex < 4; ++childIndex)
				{
				    UT_Vector3i childFace = octreeLabels.getChildFace(face, axis, childIndex);

				    if (HDKgetFieldValue(octreeVelocityIndices[level - 1][axis], childFace) < 0)
				    {
					passed = false;
					return;
				    }
				}
			    }
			    else if (velocityIndex == HDK_OUTSIDE || velocityIndex == HDK_SOLIDBOUNDARY)
			    {
				if (level != 0)
				{
				    passed = false;
				    return;
				}
			    }
			    else if (velocityIndex < 0)
			    {
				passed = false;
				return;
			    }
			}

		    for (int axis : {0,1,2})
			for (int edgeIndex = 0; edgeIndex < 4; ++edgeIndex)
                        {
			    UT_Vector3i edge = HDKcellToEdge(cell, axis, edgeIndex);

			    exint edgeStressIndex = HDKgetFieldValue(edgeStressIndices[level][axis], edge);

			    if (edgeStressIndex == HDK_UNASSIGNED)
                            {
				for (int childIndex = 0; childIndex < 2; ++childIndex)
                                {
				    UT_Vector3i childEdge = octreeLabels.getChildEdge(edge, axis, childIndex);

				    exint childStressIndex = HDKgetFieldValue(edgeStressIndices[level - 1][axis], childEdge);

				    if (childStressIndex < 0)
                                    {
                                        for (int grandchildIndex = 0; grandchildIndex < 2; ++grandchildIndex)
                                        {
					    UT_Vector3i grandchildEdge = octreeLabels.getChildEdge(childEdge, axis, grandchildIndex);

					    exint grandchildStressIndex = HDKgetFieldValue(edgeStressIndices[level - 2][axis], grandchildEdge);

                                            if (grandchildStressIndex < 0)
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
        }
    }
}

// Verify that an active cell here is active in the octree. Verify that all
// faces are active or their children are active. Verify that all edges are
// active or that their children are.
bool
HDK_AdaptiveViscosity::centerStresUnitTest(const UT_Array<SIM_RawIndexField> &centerStressIndices,
					    const UT_Array<UT_Array<SIM_RawIndexField>> &edgeStressIndices,
					    const UT_Array<UT_Array<SIM_RawIndexField>> &octreeVelocityIndices,
					    const HDK_OctreeGrid &octreeLabels) const
{
    bool passed = true;
    const int octreeLevels = octreeLabels.getOctreeLevels();

    for (int level = 0; level < octreeLevels; ++level)
    {
        centerStressUnitTest(passed,
				centerStressIndices,
				edgeStressIndices,
				octreeVelocityIndices,
				octreeLabels,
				level);

        if (!passed)
            return false;
    }
    return true;
}
