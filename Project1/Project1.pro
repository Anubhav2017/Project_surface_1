QT += core
QT -= gui
QT += widgets
CONFIG += c++11

TARGET = Project1
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
    BSpline/bsplinecurve.cpp \
    BSpline/bsplinepatchnetwork.cpp \
    BSpline/bsplineproperties.cpp \
    BSpline/bsplinesurface.cpp \
    BSpline/bsplinetriangulation.cpp \
    CellMesh/celledge.cpp \
    CellMesh/cellface.cpp \
    CellMesh/cellmesh.cpp \
    CellMesh/cellvertex.cpp \
    CellMesh/parameterization.cpp \
    CellMesh/polylinevertex.cpp \
    external/DataIO.cpp \
    external/IASS.cpp \
    MBSI/colormapgenerator.cpp \
    MBSI/subdivisionexactcurvaturemetric.cpp \
    MBSI/subdivisioninspectionanglemetric.cpp \
    MBSI/subdivisionmanualhighlightmetric.cpp \
    MBSI/subdivisionnormaldeviationmetric.cpp \
    MBSI/subdivisionsurfaceareametric.cpp \
    MBSI/subdivisionthinplateenergymetric.cpp \
    MBSI/subsurfacetree.cpp \
    MBSI/surfacesubdivider.cpp \
    MBSI/viewpointcandidategenerator.cpp \
    Mesh/edge.cpp \
    Mesh/face.cpp \
    Mesh/mesh.cpp \
    Mesh/vertex.cpp \
    SurfaceTessellation/blossomgraphmatching.cpp \
    SurfaceTessellation/quadmeshgenerator.cpp \
    SurfaceTessellation/surfacevoronoigenerator.cpp \
    energymatrixgenerator.cpp \
    parameteroptimizer.cpp \
    plywriter.cpp \
    splineevaluator.cpp \
    splinefitter.cpp \
    stlwriter.cpp \
    voxelizer.cpp


HEADERS += \
    BSpline/bsplinecurve.h \
    BSpline/bsplinepatchnetwork.h \
    BSpline/bsplineproperties.h \
    BSpline/bsplinesurface.h \
    BSpline/bsplinetriangulation.h \
    CellMesh/celledge.h \
    CellMesh/cellface.h \
    CellMesh/cellmesh.h \
    CellMesh/cellvertex.h \
    CellMesh/parameterization.h \
    CellMesh/polylinevertex.h \
    external/DataIO.h \
    external/IASS.h \
    external/Miniball.hpp \
    MBSI/colormapgenerator.h \
    MBSI/subdivisionabstractmetric.h \
    MBSI/subdivisionexactcurvaturemetric.h \
    MBSI/subdivisioninspectionanglemetric.h \
    MBSI/subdivisionmanualhighlightmetric.h \
    MBSI/subdivisionnormaldeviationmetric.h \
    MBSI/subdivisionsurfaceareametric.h \
    MBSI/subdivisionthinplateenergymetric.h \
    MBSI/subsurfacetree.h \
    MBSI/surfacesubdivider.h \
    MBSI/viewpointcandidategenerator.h \
    Mesh/edge.h \
    Mesh/face.h \
    Mesh/mesh.h \
    Mesh/vertex.h \
    SurfaceTessellation/blossomgraphmatching.h \
    SurfaceTessellation/priorityqueue.h \
    SurfaceTessellation/quadmeshgenerator.h \
    SurfaceTessellation/surfacevoronoigenerator.h \
    energymatrixgenerator.h \
    parameteroptimizer.h \
    plywriter.h \
    splineevaluator.h \
    splinefitter.h \
    stlwriter.h \
    voxelizer.h \


# Eigen v.3.3.2 [2017_02_05]
# http://eigen.tuxfamily.org/index.php?title=Main_Page
INCLUDEPATH += "../external/Eigen"

# GNU Scientific Library (GSL) [2017_02_05]
# -> Linux (v.2.3): https://www.gnu.org/software/gsl/
# -> Windows (v.1.8): http://gnuwin32.sourceforge.net/packages/gsl.htm

INCLUDEPATH += "../external/gsl/include"
LIBS += "../external/gsl/lib/libgsl.a"
LIBS += "../external/gsl/lib/libgslcblas.a"

# NLOpt 2.4.2 [2017_05_28]
# http://ab-initio.mit.edu/wiki/index.php/NLopt
INCLUDEPATH += "../external/nlopt-2.4.2/include"
win32 {
   LIBS += -L"../external/nlopt-2.4.2/lib/" -llibnlopt-0
}
unix {
   LIBS += "../external/nlopt-2.4.2/lib/libnlopt.a"
}

