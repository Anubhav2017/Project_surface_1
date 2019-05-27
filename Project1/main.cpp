#include <QCoreApplication>
#include "SurfaceTessellation/surfacevoronoigenerator.h"
#include "SurfaceTessellation/quadmeshgenerator.h"
#include "BSpline/bsplinetriangulation.h"
#include "MBSI/viewpointcandidategenerator.h"
#include "MBSI/subdivisionthinplateenergymetric.h"
#include "MBSI/subdivisionexactcurvaturemetric.h"
#include "MBSI/subdivisionsurfaceareametric.h"
#include "MBSI/subdivisionmanualhighlightmetric.h"
#include "MBSI/subdivisionnormaldeviationmetric.h"
#include "MBSI/subdivisioninspectionanglemetric.h"
#include "MBSI/colormapgenerator.h"
#include "MBSI/surfacesubdivider.h"
#include "parameteroptimizer.h"
#include "energymatrixgenerator.h"
#include "splinefitter.h"
#include "splineevaluator.h"
#include "voxelizer.h"
#include "plywriter.h"
#include"stlwriter.h"
#include <iostream>
#include <limits>

#include <QQueue>

#include <QDir>
#include <QFileDialog>
#include <QDebug>

#include "IterativeLinearSolvers"   //includes least squares solver which is only dev. version right now
#include "SparseCholesky"
#include "SparseLU"
#include <Sparse>
#include <Dense>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_integration.h>

#include <omp.h>

#include <nlopt.h>

#include "external/Miniball.hpp"
#include "external/DataIO.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    return a.exec();


}
int main2(int argc, char *argv[])
{
    //int m_startPoint = 0;
    //int m_endPoint = 6;
    //Start timer
    QTime timeTotal, timeStep;
    timeTotal.start();
    timeStep.start();
    //QString meshOpenFilename = QFileDialog::getOpenFileName(0, QString("Load triangulation"), QString("../.."), "*.stl");
    QString meshOpenFilename="/u/a/agarwala/Desktop/AdaptiveSurfaceReconstruction/project_files/stl_files/sphere-d100-sf1-ml.stl";
    qDebug()<<"Filename is"<<meshOpenFilename;
    if (meshOpenFilename.isEmpty())
        return -1;




    m_mesh = Mesh::loadFromFile(meshOpenFilename);
    if (!m_mesh) {
        writeOutput(QString("ERROR: Could not open mesh file: %1").arg(meshOpenFilename));
        exit(-1);
    }
    writeOutput(QString("mesh imported. Time: %1 ms").arg(timeStep.restart()));
    writeOutput(QString("Vertices: %1, Edges: %2, Faces: %3").arg(m_mesh->getNumberOfVertices()).arg(m_mesh->getNumberOfEdges()).arg(m_mesh->getNumberOfFaces()));

//    if (storeImportedMeshAsStl) {
//        m_stlWriter->writeMeshFacesToStl(m_mesh, "importedMesh.stl");
//        writeOutput(QString("imported Mesh stored. Time: %1 ms").arg(timeStep.restart()));
//    }

    //0 - Load or compute the isosurface
    //if (m_startPoint <= 0 && m_endPoint >= 0) {
//        QString filename = QFileDialog::getOpenFileName(0, QString("Load triangulation"), QString("../.."), "*.stl");
//        if (filename.isEmpty()) {
//            writeOutput(QString("No .stl file selected"));
//            exit(-1);
//        }

//        m_mesh = Mesh::loadFromBinaryStl(filename);

//        writeOutput("Loaded triangulation: " + filename);
//        writeOutput(QString("Vertices: %1, Edges: %2, Faces: %3").arg(m_mesh->getNumberOfVertices()).arg(m_mesh->getNumberOfEdges()).arg(m_mesh->getNumberOfFaces()));
//        if (m_mesh->getNumberOfFaces() == 0) {
//            writeOutput("Isosurface Empty!");
//            return;
//        }

//        if (storeImportedMeshAsStl) {
//            m_stlWriter->writeMeshFacesToStl(m_mesh, "exact_isosurface_noSmoothing.stl");
//            writeOutput(QString("Isosurface stored without smoothing. Time: %1 ms").arg(timeStep.restart()));
//        }
    //}
    //1 - Do mesh smoothing

 //   if (m_startPoint <= 1 && m_endPoint >= 1 && doMeshSmoothing) {
        m_mesh->smoothSurface(m_maxSmoothIterations, maxSmoothDistFromOrigin, abortSmoothAverageDist);

        writeOutput(QString("Smoothing done. Time: %1 ms").arg(timeStep.restart()));
        if (storeSmoothedMeshAsStl) {
            m_stlWriter->writeMeshFacesToStl(m_mesh, "exact_isosurface.stl");
            writeOutput(QString("Isosurface stored with smoothing. Time: %1 ms").arg(timeStep.restart()));
        }
  //  }

    if (m_scaleMesh == 0)    //fit mesh into unit cube
        m_mesh->scaleMeshToFitIntoCube(m_keepMeshProportions, m_scalingFactor);
    else if (m_scaleMesh == 1)  //fit mesh into sphere of diameter 1
        m_mesh->scaleMeshToFitIntoSphere(0.5, QVector3D(0.5, 0.5, 0.5));


    //Export Mesh
    if (exportMesh && ((m_startPoint <= 1 && m_endPoint >= 0) || m_scaleMesh > -1)) {
        m_mesh->saveToFile(meshSaveFilename);
        writeOutput(QString("Mesh exported. Time: %1 ms").arg(timeStep.restart()));
    }

    //2 - Compute the quad mesh

        QuadmeshGenerator qmGenerator(m_mesh);
        qmGenerator.setParameterizationParameters(param_bt_forQuadmeshGeneration, param_wt, m_parameterInterval, parameterizationStorageReservation_ExpectedEdgesPerVertex);
        if (m_quadmeshAlgorithm == QMA_DelaunayMatching) {
            m_cellMesh = qmGenerator.computeQuadmeshByDelaunayMatching(m_instantTopologyCheck, m_minNumberOfQuadCells, m_minNumberOfInitialCells);
            writeOutput(QString("Quadmesh computed by Delaunay matching! Number of cells: %1. Time: %2 ms").arg(m_cellMesh->getNumberOfFaces()).arg(timeStep.restart()));
        } else if (m_quadmeshAlgorithm == QMA_VoronoiSubdivision){
            m_cellMesh = qmGenerator.computeQuadmeshBySubdividingVoronoi(m_instantTopologyCheck, m_skipRefinement, m_minNumberOfQuadCells, m_minNumberOfInitialCells);
            writeOutput(QString("Quadmesh computed by Voronoi subdivision! Number of cells: %1. Time: %2 ms").arg(m_cellMesh->getNumberOfFaces()).arg(timeStep.restart()));
        } else if (m_quadmeshAlgorithm == QMA_CuttingPlanes) {
            QVector<QVector3D> basePoints, normals;

            QString filename = QFileDialog::getOpenFileName(0, QString("Load cutting planes"), m_dmFolder, "*.cp.dm");
            if (filename.isEmpty()) {
                writeOutput(QString("No cutting plane file selected"));
                exit(-1);
            }
            Controller::loadBasepointsAndNormalsFromFile(filename, &basePoints, &normals);

            m_cellMesh = qmGenerator.subdivideMeshByCuttingPlanes(normals, basePoints);
            writeOutput(QString("Quadmesh computed by use of cutting planes, loaded from: %1\n Number of cells: %2. Time: %3 ms").arg(filename).arg(m_cellMesh->getNumberOfFaces()).arg(timeStep.restart()));
        } else {
            writeOutput("ERROR: No valid quadmesh algorithm chosen!");
            exit(-1);
        }

        //TODO maybe check debugging mode and do regular output if debuging > 1
        qDebug() << "Vertex Valencies" << m_cellMesh->getValencyHistogramCellVertices() << "total:" << m_cellMesh->getNumberOfVertices();
        qDebug() << "Face Valencies" << m_cellMesh->getValencyHistogramCellFaces() << "total:" << m_cellMesh->getNumberOfFaces();

        m_cellMesh->checkAndEnforceCellOrientation();
        m_cellMesh->computeFullParameterization(param_bt_forQuadmeshGeneration);

        if (storeInitialCellMeshAsStl)
            this->cellMeshToStl(QString("initialQuadmesh"), &timeStep, cellMeshIndividualVerticesAndEdgesToStl, cellMeshContentsToStl);


    //3 - Optimize geometry of surface cells

        for (int i = 0; i < m_numberOfCentralSubdivisions; i++) {
            m_cellMesh->doSubdivisionForFullMesh();
            writeOutput(QString("Cell subdivision %1 / %2 done...").arg(i+1).arg(m_numberOfCentralSubdivisions));
        }

        for (int i = 0; i < m_numberOfVertexOptimizations; i++) {
            m_cellMesh->optimizeAllVertices();
            //m_cellMesh->optimizeAllVerticesSafe();
            writeOutput(QString("Vertex optimization %1 / %2 done...").arg(i+1).arg(m_numberOfVertexOptimizations));
        }

        m_cellMesh->checkAndEnforceCellOrientation();
        writeOutput(QString("Quadmesh optimized. Time: %1 ms").arg(timeStep.restart()));
        if (storeOptimizedCellMeshAsStl)
            this->cellMeshToStl(QString("optimizedSurfaceCellMesh"), &timeStep, cellMeshIndividualVerticesAndEdgesToStl, cellMeshContentsToStl);


    //Invert orientation
    if (m_endPoint >= 2 && m_invertCellmeshOrientation) {
        m_cellMesh->invertTotalOrientation();
        writeOutput(QString("Cellmesh orientation inverted"));
    }

    //Export cellmesh
    if (exportCellMesh && (m_startPoint <= 3 || m_invertCellmeshOrientation) && m_endPoint >= 2) {
        m_cellMesh->checkAndEnforceCellOrientation();
        if (QFile::exists(cellMeshSaveFilename)) {
            if (QFile::exists(cellMeshBackupSaveFilename))
                QFile::remove(cellMeshBackupSaveFilename);
            QFile::copy(cellMeshSaveFilename, cellMeshBackupSaveFilename);
        }

        m_cellMesh->saveToFile(cellMeshSaveFilename);
        writeOutput(QString("Cellmesh exported. Time: %1 ms").arg(timeStep.restart()));
    }

    //Export parameterization
    if (exportParameterization && m_startPoint <= 3 && m_endPoint >= 2) {
        if (saveParameterizationsAsImage)
            m_cellMesh->saveParameterizations(parameterizationSaveFilename, parameterizationImageFilename);
        else
            m_cellMesh->saveParameterizations(parameterizationSaveFilename, QString());
        writeOutput(QString("Parameterizations exported. Time: %1 ms").arg(timeStep.restart()));
    }


    //Construct energy matrix
    Eigen::SparseMatrix<double> energyMatrix(m_numberOfControlPoints * m_numberOfControlPoints, m_numberOfControlPoints * m_numberOfControlPoints);
    EnergyMatrixGenerator energyMatGenerator(m_energyFolder, m_loadEnergyMatricesIfPossible, m_storeEnergyMatricesOnHardDrive);
    if (m_startPoint <= 6 && m_endPoint >= 4 && m_weightEnergyTerm > 0.000001) {
        const int numDerivatives = energyDerivatives.size();
        for (int i = 0; i < numDerivatives; i++) {
            const double intervalStart = m_bSplineProperties->getInnerKnotsStart();
            const double intervalEnd = m_bSplineProperties->getInnerKnotsEnd();
            energyMatrix += energyDerivativeFactors[i] * energyMatGenerator.getEnergyMatrixSurface(m_bSplineProperties,
                                                                                                   energyDerivatives[i].first,
                                                                                                   energyDerivatives[i].second,
                                                                                                   intervalStart, intervalEnd, intervalStart, intervalEnd);
        }
        writeOutput(QString("Energy matrix constructed. Time %1 ms").arg(timeStep.restart()));
    }


    //4 - Do B-Spline Approximation
    if (m_startPoint <= 4 && m_endPoint >= 4) {
        //check for consistent orientation
        if (!m_cellMesh->checkForConsistentOrientation())
            writeOutput(QString("WARNING: Individual cell faces are not orientet consistently"));

        //correct parameterizations
        m_cellMesh->computeFullParameterization(param_bt_forSplineFitting);
        const double paramCorrection = m_cellMesh->limitQuadCellParameterizationsToParameterIntervall();
        if (paramCorrection > 0)
            writeOutput(QString("Parameterizations limited to [%1,%2]. Total error: %3").arg(m_parameterInterval.first).arg(m_parameterInterval.second).arg(paramCorrection));

        Splinefitter bSplineFitter(m_cellMesh, m_bSplineProperties);
        bSplineFitter.setEnergyMatrix(energyMatrix);
        bSplineFitter.setWeightEnergyTerm(m_weightEnergyTerm, true);
        bSplineFitter.setNumberOfG1EdgePoints(m_numberOfG1EdgePoints);
        bSplineFitter.setNumberOfG1VertexPoints(m_numberOfG1VertexPoints);
        bSplineFitter.setSizeOfG1VertexDomain(m_g1VertexParameterRange);
        bSplineFitter.setNloMaxIterations(m_nloptMaxIterations);
        bSplineFitter.setNloXPrec(m_nloptXPrec);
        bSplineFitter.setNloFPrec(m_nloptFPrec);
        bSplineFitter.setNloEqualConstPrec(m_nloptConstrPrec);
        bSplineFitter.setNloG0constPrec(m_nloptConstrPrec);
        bSplineFitter.setNloG1constPrec(m_nloptConstrPrec);
        bSplineFitter.setG1VertexConstraintsAsOpt(m_doVertexFittingWithG1Opt);
        bSplineFitter.setG1VertexConstraintsOptFactor(m_vertexFitG1Weight);
        bSplineFitter.setG1EdgeConstraintsAsOpt(m_doEdgeFittingWithG1Opt);
        bSplineFitter.setG1EdgeConstraintsOptFactor(m_edgeFitG1Weight);
        bSplineFitter.setDetailedOutput(m_detailedFittingOutput);

        if (m_useInitialGuess)
            bSplineFitter.setInitialGuess(m_bSplinePatchNetwork);

        //Do B-spline fitting
        if (m_bSplineAlgorithm == BSA_CurveFittingG0) {
            double wEnergy = m_weightEnergyTerm * m_weightEnergyTerm;
            double wLS = 1 - wEnergy;
            m_bSplinePatchNetwork = bSplineFitter.doApproximationLocalC0WithBoundaryInterpolation(true, wLS, wEnergy);
        } else if (m_bSplineAlgorithm == BSA_CurveFittingG0Iterative) {
            QVector<double> wEnergyCurve, wEnergySurf;
            //wEnergyCurve << 0.5 << 0.95 << 0.75 << 0.5 << 0.25 << 0.1 << 0.01 << 0.001 << 0.0001 << 0;
            wEnergyCurve << 0.5 << 0.95 << 0.66 << 0.33 << 0.1 << 0.001 << m_weightEnergyTerm * m_weightEnergyTerm;// << 0.0001 << 0;
            wEnergySurf << 0.5 << 0.95 << 0.75 << 0.5 << 0.25 << 0.1 << 0.01 << m_weightEnergyTerm;// << 0.001 << 0.0001 << 0;
            m_bSplinePatchNetwork = bSplineFitter.doIterativeCurveAndSurfaceFittingC0BoundaryInterpolation(wEnergyCurve, wEnergySurf, true);
        } else if (m_bSplineAlgorithm == BSA_GlobalG0) {
            m_bSplinePatchNetwork = bSplineFitter.doApproximationGlobalC0G1Constraints(true);
        } else if (m_bSplineAlgorithm == BSA_SimplifiedGlobalG1) {
            m_bSplinePatchNetwork = bSplineFitter.doApproximationGlobalC0G1Constraints(false);
        } else if (m_bSplineAlgorithm == BSA_LocalG1Multiphase) {
            m_bSplinePatchNetwork = bSplineFitter.doApproximationLocalG0ConstraintG1Sample();
        } else {
            writeOutput("ERROR: No valid b-spline algorithm chosen!");
            exit(-1);
        }

        writeOutput(QString("B-Spline Approximation done. Algorithm %1. Time: %2 ms").arg(m_bSplineAlgorithm).arg(timeStep.restart()));

        if (storeComputedBSplinesAsStl)
            this->bSplinesToStl(m_stlBSplineSamplePoints, QString("finalBSpline"), &timeStep, bSplineControlPointsToStl, bSplineKnotLinesToStl, bSplineDataErrorToStl, bSplineDataErrorOfIndividualSurfacesToStl);
    }

//    for (int i = 0; i < m_bSplineSurfaces.size(); i++)
//        qDebug() << calculatePatchEnergy(m_bSplineSurfaces[i], energyMatrix, m_numberOfControlPoints);

    //Export B-Splines
    if (exportBSplines && m_startPoint <= 4 && m_endPoint >= 4) {
        if (QFile::exists(bSplineSaveFilename)) {
            if (QFile::exists(bSplineBackupSaveFilename))
                QFile::remove(bSplineBackupSaveFilename);
            QFile::copy(bSplineSaveFilename, bSplineBackupSaveFilename);
        }
        m_bSplinePatchNetwork->saveToFile(bSplineSaveFilename);
        writeOutput(QString("B-Splines exported. Time: %1 ms").arg(timeStep.restart()));
    }

    //5 - Optimize parameters (projection onto tangential space)
    if (m_startPoint <= 5 && m_endPoint >= 5 && doParamOpt) {
        //optimize parameterizations
        for (int iOpt = 0; iOpt < m_numberOfParameterOptimizations; iOpt++) {
            const double improvement = ParameterOptimizer::optimizeParameterizationFullSurfaceNetwork(m_bSplinePatchNetwork, true, m_useNLOforParameterOpt);

            writeOutput(QString("Parameterizations optimized, error new/old: %1 [%2 / %3]. Time: %4 ms").arg(improvement).arg(iOpt + 1).arg(m_numberOfParameterOptimizations).arg(timeStep.restart()));
        }

        //save parameterizations
        if (exportParameterization) {
            if (saveParameterizationsAsImage)
                m_cellMesh->saveParameterizations(parameterizationSaveFilename, parameterizationImageFilename);
            else
                m_cellMesh->saveParameterizations(parameterizationSaveFilename, QString());
            writeOutput(QString("Parameterizations exported. Time: %1 ms").arg(timeStep.restart()));
        }
    }

    //TODO measure:
    //-energy
    //-lenghts of border curves (individual and total)

    //6 - Error measurement
    if (m_startPoint <= 6 && m_endPoint >= 6) {
        //correct parameterizations (in case of new parameterization or step 5 messed something up)
        const double paramCorrection = m_cellMesh->limitQuadCellParameterizationsToParameterIntervall();
        if (paramCorrection > 0)
            writeOutput(QString("Parameterizations limited to [%1,%2]. Total error: %3").arg(m_parameterInterval.first).arg(m_parameterInterval.second).arg(paramCorrection));

        if (m_scaleForEvaluation == 0) {
            double minX, maxX, minY, maxY, minZ, maxZ;
            m_mesh->calculateMinMaxVertexPositions(minX, maxX, minY, maxY, minZ, maxZ);
            m_mesh->scaleMeshToFitIntoCube(m_keepMeshProportions, m_scalingFactor);

            double difX = maxX - minX;
            double difY = maxY - minY;
            double difZ = maxZ - minZ;
            if (m_keepMeshProportions) {
                double max = difX;
                if (difY > max)
                    max = difY;
                if (difZ > max)
                    max = difZ;
                difX = max;
                difY = max;
                difZ = max;
            }

            const QVector3D translationVector(-minX, -minY, -minZ);

            m_cellMesh->translate(translationVector);
            m_cellMesh->scale(m_scalingFactor / difX, m_scalingFactor / difY, m_scalingFactor / difZ);

            m_bSplinePatchNetwork->translate(translationVector);
            m_bSplinePatchNetwork->scale(m_scalingFactor / difX, m_scalingFactor / difY, m_scalingFactor / difZ);
        } else if (m_scaleForEvaluation == 1) {

            double radius;
            QVector3D center;
            m_mesh->calculateMinimumBoundingSphere(radius, center);

            const double scalingFactor = m_scalingFactor/(2 * radius);
            const QVector3D targetCenter(m_scalingFactor/2, m_scalingFactor/2, m_scalingFactor/2);
            m_mesh->scaleMeshToFitIntoSphere(m_scalingFactor/2, targetCenter);

            m_cellMesh->translate(-center);
            m_cellMesh->scale(scalingFactor, scalingFactor, scalingFactor);
            m_cellMesh->translate(targetCenter);

            m_bSplinePatchNetwork->translate(-center);
            m_bSplinePatchNetwork->scale(scalingFactor, scalingFactor, scalingFactor);
            m_bSplinePatchNetwork->translate(targetCenter);

        }

        //evaluate
        SplineEvaluator evaluator(m_bSplinePatchNetwork, &energyMatrix, numberOfEvalPointsForG1, doOptimizationForDataError);

        evaluator.writeQualityFile(qualityFilename);

        writeOutput(QString("Spline area: %1 Energy: %2 Energy/Area: %3").arg(evaluator.getTotalSurfaceAreaBSplines()).arg(evaluator.getTotalEnergy()).arg(evaluator.getTotalEnergy()/evaluator.getTotalSurfaceAreaBSplines()));
        writeOutput(QString("RMS Error: %1 Max Error: %2").arg(evaluator.getRmsError()).arg(evaluator.getMaxDataError()));
        writeOutput(QString("BBNormalized: %1 BBNMax: %2").arg(evaluator.getBoundingBoxScaledRmsError()).arg(evaluator.getBoundingBoxScaledMaxDataError()));
        writeOutput(QString("BSNormalized [%]: %1 BSNMax: %2").arg(evaluator.getBoundingSphereScaledRmsError()).arg(evaluator.getBoundingSphereScaledMaxDataError()));
        writeOutput(QString("Average G1 Angle Error: %1 Maximum Error: %2").arg(evaluator.getTotalAverageAngleError()).arg(evaluator.getMaxAngleError()));
        writeOutput(QString("B-Spline quality measured. Time: %1 ms").arg(timeStep.restart()));

//        //matlab export
//        if (m_doMatlabOutput) {
//            Controller::removeRecursively(m_dmFolder + matlabSubfolder);
//            if (!QDir(m_dmFolder + matlabSubfolder).exists())
//                QDir().mkpath(m_dmFolder + matlabSubfolder);

//            evaluator.writeMatlabTriangulationStructureFile(matlabTriStructureFilename);
//            evaluator.writeMatlabDataErrorFiles(matlabTriErrorFilename, matlabBSErrorFilename);
//            evaluator.writeMatlabG1ErrorAtPLVFiles(matlabTriBordersFilename, matlabBSBorderFilename);
//            evaluator.writeMatlabRegularSamplingStructureFile(numberOfEvalPointsForCurvature, matlabRegularSamplingStructureFilename);
//            evaluator.writeMatlabRegularSampledCurvature(numberOfEvalPointsForCurvature, matlabRegularSamplingCurvatureFilename);
//            evaluator.writeMatlabRegularSamplingG1Error(matlabRegularSamplingBorderFilename);

//            writeOutput(QString("Matlab export done: %1 ms").arg(timeStep.restart()));
//        }

        if (m_do3DImageOutput) {
            int imageSizeX = m_3DimageSizeX;
            int imageSizeY = m_3DimageSizeY;
            int imageSizeZ = m_3DimageSizeZ;
            QVector3D translationVector(0, 0, 0);
            if (m_3DimageSizeX == -1) {
                int minX, maxX, minY, maxY, minZ, maxZ;
                m_bSplinePatchNetwork->getBoundingBox(minX, maxX, minY, maxY, minZ, maxZ);
                //translate by -(minX, minY, minZ), so that lower corner of the bounding box is (0, 0, 0)
                //add padding in each direction
                translationVector = QVector3D(iassPadding - minX, iassPadding - minY, iassPadding - minZ);
                imageSizeX = 2 * iassPadding + maxX - minX;
                imageSizeY = 2 * iassPadding + maxY - minY;
                imageSizeZ = 2 * iassPadding + maxZ - minZ;
            }

            m_bSplinePatchNetwork->translate(translationVector);
            Voxelizer::splinesToVoxelImage(m_bSplinePatchNetwork, imageSizeX, imageSizeY, imageSizeZ, iassOutputFilename, iassSurfaceSamplePoints);
            m_bSplinePatchNetwork->translate(-translationVector);

            writeOutput(QString("Volume image done. Time: %1 ms").arg(timeStep.restart()));
        }
    }

    if (m_startPoint <= 7 && m_endPoint >= 7) {
        //Generation of view point candidates (VPC) for surface inspection

        //const int m_mbsiMetricType = 1; //TODO replace with enum later and make it a parameter
        const int maxSubdivisionDepth = 50; //TODO make that a parameter

        const int precomputeDepth = 1;  //Only used by thin plate metric
        const int surfaceMetricNumericalSampleRatePerInterval = 2 * m_splineOrder;
        const bool integrateCurvatureIn3D = true;

        QString factorToken("a");
        if (m_mbsiThresholdIsAVGFactor)
            factorToken = "f";
        if (m_mbsiTargetNumberOfVPC > -1) {
            factorToken = "n";
            m_mbsiThreshold = m_mbsiTargetNumberOfVPC;
        }

        //Set up metric and compute average for entire surface network
        SubdivisionAbstractMetric *metric = 0;
        if (m_mbsiMetricType == 0) {
            metric = new SubdivisionThinPlateEnergyMetric(&energyMatGenerator, m_bSplineProperties);
            ((SubdivisionThinPlateEnergyMetric *) metric)->precomputeThinPlateMatrices(precomputeDepth);
        } else if (m_mbsiMetricType == 1) {
            metric = new SubdivisionExactCurvatureMetric(surfaceMetricNumericalSampleRatePerInterval, integrateCurvatureIn3D);
        } else if (m_mbsiMetricType == 2) {
            metric = new SubdivisionSurfaceAreaMetric(surfaceMetricNumericalSampleRatePerInterval);
        } else if (m_mbsiMetricType == 3) {
            metric = new SubdivisionNormalDeviationMetric(surfaceMetricNumericalSampleRatePerInterval);
        } else if (m_mbsiMetricType == 4) {
            metric = new SubdivisionInspectionAngleMetric(m_mbsiCamDistance, surfaceMetricNumericalSampleRatePerInterval);
        } else if (m_mbsiMetricType == 5) {
            //Stretched cylinder
            metric = new SubdivisionManualHighlightMetric(surfaceMetricNumericalSampleRatePerInterval);
            ((SubdivisionManualHighlightMetric *) metric)->addHighlightBoundingSphere(QVector3D(21, -9, 18), 5, 1);
        } else if (m_mbsiMetricType == 6) {
            //Spring
            metric = new SubdivisionManualHighlightMetric(surfaceMetricNumericalSampleRatePerInterval);
            ((SubdivisionManualHighlightMetric *) metric)->addHighlightBoundingSphere(QVector3D(4, 4, 20), 3, 1);
        } else if (m_mbsiMetricType == 7) {
            //Tangle cube
            metric = new SubdivisionManualHighlightMetric(surfaceMetricNumericalSampleRatePerInterval);
            ((SubdivisionManualHighlightMetric *) metric)->addHighlightBoundingSphere(QVector3D(70, 70, 70), 20, 1);
        } else {
            qDebug() << "MBSI: No valid metric specified!";
            exit(-1);
        }
//            srand(time(NULL));
//            const int numberOfRandomHighlights = 5;
//            const double randomHighlighRadius = 0.15;
//            for (int i = 0; i < numberOfRandomHighlights; i++) {
//                const int patchId = rand() % m_bSplinePatchNetwork->getNumberOfBSplineSurfaces();
//                const double randomU = m_parameterInterval.first + (((double) rand() / (RAND_MAX)) * (m_parameterInterval.second - m_parameterInterval.first));
//                const double randomV = m_parameterInterval.first + (((double) rand() / (RAND_MAX)) * (m_parameterInterval.second - m_parameterInterval.first));
//                const QVector3D point = m_bSplinePatchNetwork->getBSplineSurface(patchId)->evaluate(randomU, randomV);

//                qDebug() << patchId << randomU << randomV;
//                qDebug() << point << randomHighlighRadius;

//                ((SubdivisionManualHighlightMetric *) metric)->addHighlightBoundingSphere(point, randomHighlighRadius, 1);
//            }
        const QString metricNameToken = metric->getMetricNameToken();

        //Define threshold
        const int numberOfSurfaces = m_bSplinePatchNetwork->getNumberOfBSplineSurfaces();
        double averageMetric = 1;
        if (m_mbsiThresholdIsAVGFactor && m_mbsiTargetNumberOfVPC < 0) {
            averageMetric = 0;
            //TODO these values get recomputed at depth = 0 for the subdivision
            for (int i = 0; i < numberOfSurfaces; i++) {
                const double value = metric->evaluate(m_bSplinePatchNetwork->getBSplineSurface(i),
                                                      m_parameterInterval.first, m_parameterInterval.second,
                                                      m_parameterInterval.first, m_parameterInterval.second);
                averageMetric += value;
                qDebug() << i << value;
            }
            averageMetric /= (double) numberOfSurfaces;
            qDebug() << "Total:" << averageMetric * numberOfSurfaces << "Average:" << averageMetric;
        }
        const double subdivisionThreshold = m_mbsiThreshold * averageMetric;

        writeOutput(QString("Metric and threshold setup done. Absolute threshold value: %1, time: %2").arg(subdivisionThreshold).arg(timeStep.restart()));


        //Do the subdivision
        QTime timeSubdivision;
        timeSubdivision.start();

        QVector<SubSurfaceTree *> subdivisions;
        if (m_mbsiTargetNumberOfVPC > -1)
            subdivisions = SurfaceSubdivider::subdivideSurfacesUntilFixedNumber(m_bSplinePatchNetwork->getBSplineSurfaces(), metric, 100, maxSubdivisionDepth);
        else
            subdivisions = SurfaceSubdivider::subdivideSurfacesThresholdBased(m_bSplinePatchNetwork->getBSplineSurfaces(), metric, subdivisionThreshold, maxSubdivisionDepth);

        writeOutput(QString("Subdivision done! Time %1 ms").arg(timeSubdivision.restart()));

        //Write subdivisions to stl
        if (writeMBSISubdivisionsToSTL)
            m_stlWriter->writeBSplineSubdivisionsToStl(m_bSplinePatchNetwork->getBSplineSurfaces(), subdivisions, QString("subdivisions_%1_th%2%3_LW%4.stl").arg(metricNameToken).arg(m_mbsiThreshold).arg(factorToken).arg(m_stlLineWidth), m_stlLineWidth, m_stlBSplineSamplePoints);

        //Set up vectors for output
        QVector<QVector3D> pivotPointsCurv;
        QVector<QVector3D> camPositionsCurv;
        QVector<int> depths;

        //Calculate VPC
        ViewpointCandidateGenerator viewpointGenerator(m_bSplinePatchNetwork);
        viewpointGenerator.generateViewPoints_subdivision(camPositionsCurv, pivotPointsCurv, depths, m_mbsiCamDistance, subdivisions);

        //Process result
        const int numSamplePoints = camPositionsCurv.size();
        QVector<QVector3D> normalsCurv(numSamplePoints);
        QVector<Edge *> camDirectionsCurv(numSamplePoints);
        for (int i = 0; i < numSamplePoints; i++) {
            normalsCurv[i] = (camPositionsCurv[i] - pivotPointsCurv[i]).normalized();   //This is not computed by the VPC generator anymore
            Vertex *v0 = new Vertex(pivotPointsCurv[i], -1);
            Vertex *v1 = new Vertex(pivotPointsCurv[i] + m_mbsiCamDirVectorLength * normalsCurv[i], -1);
            camDirectionsCurv[i] = new Edge(v0, v1, -1);
        }

        if (writeMBSIViewpointCandidatesToSTL) {
            m_stlWriter->writePoints(&camPositionsCurv, m_mbsiCamPointSize, QString("vpc_cam_%1_th%2%3_size%4.stl").arg(metricNameToken).arg(m_mbsiThreshold).arg(factorToken).arg(m_mbsiCamPointSize));
            m_stlWriter->writePoints(&pivotPointsCurv, m_mbsiPivotPointSize, QString("vpc_piv_%1_th%2%3_size%4.stl").arg(metricNameToken).arg(m_mbsiThreshold).arg(factorToken).arg(m_mbsiPivotPointSize));
            m_stlWriter->writeEdgesToStl(&camDirectionsCurv, QString("vpc_normals_%1_th%2%3_size%4_dist%5.stl").arg(metricNameToken).arg(m_mbsiThreshold).arg(factorToken).arg(m_mbsiCamDirVectorSize).arg(m_mbsiCamDirVectorLength), m_mbsiCamDirVectorSize);
        }

        std::vector<std::vector<double> > positionList(numSamplePoints);
        for (int i = 0; i < numSamplePoints; i++) {
            std::vector<double> point(10);
            //Cam position
            point[0] = camPositionsCurv[i].x();
            point[1] = camPositionsCurv[i].y();
            point[2] = camPositionsCurv[i].z();
            //Cam direction is the negative of the normal vector
            const QVector3D normal = normalsCurv[i];
            point[3] = -normal.x();
            point[4] = -normal.y();
            point[5] = -normal.z();
            //Orientation (up vector)
            QVector3D up = QVector3D(0, 0, 1) - QVector3D::dotProduct(QVector3D(0, 0, 1), normal) * normal;
            if (up == QVector3D(0, 0, 0))
                up = QVector3D(0, 1, 0);
            up.normalize();
            point[6] = up.x();
            point[7] = up.y();
            point[8] = up.z();

            //depth of the subdivision
            point[9] = depths[i];

            positionList[i] = point;
        }

        QString positionFilename(m_dmFolder + QString("VPC_%1_dist%2_th%3%4.txt").arg(metricNameToken).arg(m_mbsiCamDistance).arg(m_mbsiThreshold).arg(factorToken));
        std::string stdFilename = positionFilename.toStdString();
        writePositions(positionList, stdFilename);

        writeOutput(QString("Subdivision based sample points done. Metric: %1, number of points %2. Time: %3").arg(metricNameToken).arg(numSamplePoints).arg(timeStep.restart()));

        if ((m_mbsiMetricType == 0 || m_mbsiMetricType == 1 || m_mbsiMetricType == 4) && doExtendedMBSIOutput) {
            QTime t;
            t.start();
            qDebug() << "starting timing for new plv output";
            BSplineTriangulation *triangulation = BSplineTriangulation::createRegularSampledTriangulation(m_bSplinePatchNetwork->getNumberOfBSplineSurfaces(), m_stlBSplineSamplePoints,
                                                                                                          m_parameterInterval.first, m_parameterInterval.second, m_parameterInterval.first, m_parameterInterval.second,
                                                                                                          true);    //TODO test with false as well
            qDebug() << "triangulation done" << t.restart();
            ColorMapGenerator colMapGenerator(true);    //TODO test with different fixed color ranges
            colMapGenerator.setColorMapToHeat(false);
            QVector<QColor> colorMap = colMapGenerator.createColorMap(m_bSplinePatchNetwork->getBSplineSurfaces(), triangulation, metric);
            qDebug() << "color map done" << t.restart();
            PlyWriter::writeBSplineTriangulationWithVertexColors(m_bSplinePatchNetwork, triangulation, colorMap, m_dmFolder + QString("metricValues_%1.ply").arg(metricNameToken));
            qDebug() << "ply output done" << t.restart();
            QString matlabPlotCommand = colMapGenerator.generateMatlabPlotCommand(101, QString("./cbar_%1.png").arg(metricNameToken));
            std::cout << matlabPlotCommand.toStdString() << std::endl;
            QFile file(m_dmFolder + QString("colorbar_%1.m").arg(metricNameToken));
            file.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text);
            QTextStream out(&file);
            out << matlabPlotCommand;
            file.close();
            qDebug() << "matlab colorbar done" << t.restart();
    //        for (int i = 0; i < numberOfSurfaces; i++) {
    //            qDebug() << "surface" << i;
    //            this->createMatlab2DOutputForMetric(m_bSplinePatchNetwork->getBSplineSurface(i), QString("../out/metric_%1_%2.csv").arg(metricNameToken).arg(i), metric, subdivisions[i], surfaceMetricNumericalSampleRate, m_parameterInterval);
    //        }
    //        qDebug() << "2D matlab output done" << t.restart();
        }

        delete metric;
        foreach (SubSurfaceTree *subsurface, subdivisions)
            delete subsurface;

        //Sample point generation with reference algorithms
        if (m_mbsiDoReferenceAlgorithms) {
            QVector<QVector3D> camPositionsVertexBased, pivotPointsVertexBased, camPositionsRandomCenter, pivotPointsRandomCenter, camPositionsRandomVertex, pivotPointsRandomVertex;
            viewpointGenerator.generateViewPoints_useVertices(camPositionsVertexBased, pivotPointsVertexBased, m_mbsiCamDistance);
            viewpointGenerator.generateViewPoints_random(camPositionsRandomCenter, pivotPointsRandomCenter, m_mbsiCamDistance, numSamplePoints, false);
            viewpointGenerator.generateViewPoints_random(camPositionsRandomVertex, pivotPointsRandomVertex, m_mbsiCamDistance, numSamplePoints, true);

            QVector<Edge *> camDirectionsVertexBased(camPositionsVertexBased.size());
            for (int i = 0; i < camPositionsVertexBased.size(); i++) {
                Vertex *v0 = new Vertex(camPositionsVertexBased[i], -1);
                Vertex *v1 = new Vertex(pivotPointsVertexBased[i], -1);
                camDirectionsVertexBased[i] = new Edge(v0, v1, -1);
            }
            QVector<Edge *> camDirectionsRandomCenter(camPositionsRandomCenter.size());
            for (int i = 0; i < camPositionsRandomCenter.size(); i++) {
                Vertex *v0 = new Vertex(camPositionsRandomCenter[i], -1);
                Vertex *v1 = new Vertex(pivotPointsRandomCenter[i], -1);
                camDirectionsRandomCenter[i] = new Edge(v0, v1, -1);
            }
            QVector<Edge *> camDirectionsRandomVertex(camPositionsRandomVertex.size());
            for (int i = 0; i < camPositionsRandomVertex.size(); i++) {
                Vertex *v0 = new Vertex(camPositionsRandomVertex[i], -1);
                Vertex *v1 = new Vertex(pivotPointsRandomVertex[i], -1);
                camDirectionsRandomVertex[i] = new Edge(v0, v1, -1);
            }
            m_stlWriter->writePoints(&camPositionsVertexBased, m_mbsiCamPointSize, QString("mbsi_vertex_camPos_%1.stl").arg(m_mbsiCamPointSize));
            m_stlWriter->writePoints(&pivotPointsVertexBased, m_mbsiPivotPointSize, QString("mbsi_vertex_pivPos_%1.stl").arg(m_mbsiPivotPointSize));
            m_stlWriter->writeEdgesToStl(&camDirectionsVertexBased, QString("mbsi_vertex_normals_%1.stl").arg(m_mbsiCamDirVectorSize), m_mbsiCamDirVectorSize);

            m_stlWriter->writePoints(&camPositionsRandomCenter, m_mbsiCamPointSize, QString("mbsi_rndCent_camPos_%1.stl").arg(m_mbsiCamPointSize));
            m_stlWriter->writePoints(&pivotPointsRandomCenter, m_mbsiPivotPointSize, QString("mbsi_rndCent_pivPos_%1.stl").arg(m_mbsiPivotPointSize));
            m_stlWriter->writeEdgesToStl(&camDirectionsRandomCenter, QString("mbsi_rndCent_normals_%1.stl").arg(m_mbsiCamDirVectorSize), m_mbsiCamDirVectorSize);

            m_stlWriter->writePoints(&camPositionsRandomVertex, m_mbsiCamPointSize, QString("mbsi_rndVert_camPos_%1.stl").arg(m_mbsiCamPointSize));
            m_stlWriter->writePoints(&pivotPointsRandomVertex, m_mbsiPivotPointSize, QString("mbsi_rndVert_pivPos_%1.stl").arg(m_mbsiPivotPointSize));
            m_stlWriter->writeEdgesToStl(&camDirectionsRandomVertex, QString("mbsi_rndVert_normals_%1.stl").arg(m_mbsiCamDirVectorSize), m_mbsiCamDirVectorSize);

            writeOutput(QString("Sample points based on reference algorithms done. Time: %2").arg(timeStep.restart()));
        }

    }
    writeOutput(QString("All done! Total Time: %1 ms").arg(timeTotal.elapsed()));

}
