QMAKE_CXXFLAGS += -fpermissive

INCLUDEPATH += $$PWD
INCLUDEPATH += $$PWD/patterns

LPSOLVER_PATH = $$PWD/patterns/lp_solve_5.5
OPENMESH_PATH = $$PWD/patterns/OpenMesh
PATCHGEN_PATH = $$PWD/patterns/patchgen
KT_PATH = $$PWD/patterns/kt84


#Libs
LIBS += -ldl -lm -lmpfr


#Main files
SOURCES += \
    $$PWD/patterns/Patch.cpp \
    $$PWD/patterns/determine_geometry.cpp \
    $$PWD/patterns/generate_patch.cpp

HEADERS += \
    $$PWD/patterns/Patch.h \
    $$PWD/patterns/determine_geometry.h \
    $$PWD/patterns/generate_patch.h


#LPSolver
INCLUDEPATH += $$LPSOLVER_PATH
INCLUDEPATH += $$LPSOLVER_PATH/bfp
INCLUDEPATH += $$LPSOLVER_PATH/bfp/bfp_LUSOL
INCLUDEPATH += $$LPSOLVER_PATH/bfp/bfp_LUSOL/LUSOL
INCLUDEPATH += $$LPSOLVER_PATH/colamd
INCLUDEPATH += $$LPSOLVER_PATH/shared

DEFINES += LoadInverseLib=FALSE
DEFINES += LoadLanguageLib=FALSE
DEFINES += LoadableBlasLib=0
DEFINES += INVERSE_ACTIVE=INVERSE_LUSOL
DEFINES += RoleIsExternalInvEngine
DEFINES += NOFORTIFY
DEFINES += CHECK_SOLUTION
DEFINES += YY_NEVER_INTERACTIVE
DEFINES += PARSER_LP
DEFINES += LP_MAXLINELEN=0


SOURCES += \
    $$LPSOLVER_PATH/bfp/bfp_LUSOL/lp_LUSOL.c \
    $$LPSOLVER_PATH/bfp/bfp_LUSOL/LUSOL/lusol.c \
    $$LPSOLVER_PATH/colamd/colamd.c \
    $$LPSOLVER_PATH/fortify.c \
    $$LPSOLVER_PATH/ini.c \
    $$LPSOLVER_PATH/lp_Hash.c \
    $$LPSOLVER_PATH/lp_MDO.c \
    $$LPSOLVER_PATH/lp_MPS.c \
    $$LPSOLVER_PATH/lp_SOS.c \
    $$LPSOLVER_PATH/lp_crash.c \
    $$LPSOLVER_PATH/lp_lib.c \
    $$LPSOLVER_PATH/lp_matrix.c \
    $$LPSOLVER_PATH/lp_mipbb.c \
    $$LPSOLVER_PATH/lp_params.c \
    $$LPSOLVER_PATH/lp_presolve.c \
    $$LPSOLVER_PATH/lp_price.c \
    $$LPSOLVER_PATH/lp_pricePSE.c \
    $$LPSOLVER_PATH/lp_report.c \
    $$LPSOLVER_PATH/lp_rlp.c \
    $$LPSOLVER_PATH/lp_scale.c \
    $$LPSOLVER_PATH/lp_simplex.c \
    $$LPSOLVER_PATH/lp_utils.c \
    $$LPSOLVER_PATH/lp_wlp.c \
    $$LPSOLVER_PATH/shared/commonlib.c \
    $$LPSOLVER_PATH/shared/mmio.c \
    $$LPSOLVER_PATH/shared/myblas.c \
    $$LPSOLVER_PATH/yacc_read.c
HEADERS += \
    $$LPSOLVER_PATH/colamd/colamd.h \
    $$LPSOLVER_PATH/declare.h \
    $$LPSOLVER_PATH/fortify.h \
    $$LPSOLVER_PATH/ini.h \
    $$LPSOLVER_PATH/lp_Hash.h \
    $$LPSOLVER_PATH/lp_MDO.h \
    $$LPSOLVER_PATH/lp_MPS.h \
    $$LPSOLVER_PATH/lp_SOS.h \
    $$LPSOLVER_PATH/lp_bit.h \
    $$LPSOLVER_PATH/lp_crash.h \
    $$LPSOLVER_PATH/lp_explicit.h \
    $$LPSOLVER_PATH/lp_fortify.h \
    $$LPSOLVER_PATH/lp_lib.h \
    $$LPSOLVER_PATH/lp_matrix.h \
    $$LPSOLVER_PATH/lp_mipbb.h \
    $$LPSOLVER_PATH/lp_presolve.h \
    $$LPSOLVER_PATH/lp_price.h \
    $$LPSOLVER_PATH/lp_pricePSE.h \
    $$LPSOLVER_PATH/lp_report.h \
    $$LPSOLVER_PATH/lp_rlp.h \
    $$LPSOLVER_PATH/lp_rlp.y \
    $$LPSOLVER_PATH/lp_scale.h \
    $$LPSOLVER_PATH/lp_simplex.h \
    $$LPSOLVER_PATH/lp_types.h \
    $$LPSOLVER_PATH/lp_utils.h \
    $$LPSOLVER_PATH/lp_wlp.h \
    $$LPSOLVER_PATH/lpkit.h \
    $$LPSOLVER_PATH/lpsolve.h \
    $$LPSOLVER_PATH/shared/commonlib.h \
    $$LPSOLVER_PATH/shared/mmio.h \
    $$LPSOLVER_PATH/shared/myblas.h \
    $$LPSOLVER_PATH/ufortify.h \
    $$LPSOLVER_PATH/yacc_read.h


#OpenMesh
#DEFINES += _USE_MATH_DEFINES
#DEFINES += NOMINMAX

SOURCES += \
    $$OPENMESH_PATH/Core/Geometry/NormalConeT.cc \
    $$OPENMESH_PATH/Core/IO/BinaryHelper.cc \
    $$OPENMESH_PATH/Core/IO/IOManager.cc \
    $$OPENMESH_PATH/Core/IO/OMFormat.cc \
    $$OPENMESH_PATH/Core/IO/OMFormatT.cc \
    $$OPENMESH_PATH/Core/IO/reader/BaseReader.cc \
    $$OPENMESH_PATH/Core/IO/reader/OBJReader.cc \
    $$OPENMESH_PATH/Core/IO/reader/OFFReader.cc \
    $$OPENMESH_PATH/Core/IO/reader/OMReader.cc \
    $$OPENMESH_PATH/Core/IO/reader/PLYReader.cc \
    $$OPENMESH_PATH/Core/IO/reader/STLReader.cc \
    $$OPENMESH_PATH/Core/IO/writer/BaseWriter.cc \
    $$OPENMESH_PATH/Core/IO/writer/OBJWriter.cc \
    $$OPENMESH_PATH/Core/IO/writer/OFFWriter.cc \
    $$OPENMESH_PATH/Core/IO/writer/OMWriter.cc \
    $$OPENMESH_PATH/Core/IO/writer/PLYWriter.cc \
    $$OPENMESH_PATH/Core/IO/writer/STLWriter.cc \
    $$OPENMESH_PATH/Core/IO/writer/VTKWriter.cc \
    $$OPENMESH_PATH/Core/Mesh/ArrayKernel.cc \
    $$OPENMESH_PATH/Core/Mesh/ArrayKernelT.cc \
    $$OPENMESH_PATH/Core/Mesh/BaseKernel.cc \
    $$OPENMESH_PATH/Core/Mesh/PolyConnectivity.cc \
    $$OPENMESH_PATH/Core/Mesh/PolyMeshT.cc \
    $$OPENMESH_PATH/Core/Mesh/TriConnectivity.cc \
    $$OPENMESH_PATH/Core/Mesh/TriMeshT.cc \
    $$OPENMESH_PATH/Core/System/omstream.cc \
    $$OPENMESH_PATH/Core/Utils/BaseProperty.cc \
    $$OPENMESH_PATH/Core/Utils/Endian.cc \
    $$OPENMESH_PATH/Core/Utils/RandomNumberGenerator.cc \
    $$OPENMESH_PATH/Core/Utils/SingletonT.cc

HEADERS += \
    $$OPENMESH_PATH/Core/Geometry/Config.hh \
    $$OPENMESH_PATH/Core/Geometry/LoopSchemeMaskT.hh \
    $$OPENMESH_PATH/Core/Geometry/MathDefs.hh \
    $$OPENMESH_PATH/Core/Geometry/NormalConeT.hh \
    $$OPENMESH_PATH/Core/Geometry/Plane3d.hh \
    $$OPENMESH_PATH/Core/Geometry/QuadricT.hh \
    $$OPENMESH_PATH/Core/Geometry/VectorT.hh \
    $$OPENMESH_PATH/Core/Geometry/VectorT_inc.hh \
    $$OPENMESH_PATH/Core/IO/BinaryHelper.hh \
    $$OPENMESH_PATH/Core/IO/IOInstances.hh \
    $$OPENMESH_PATH/Core/IO/IOManager.hh \
    $$OPENMESH_PATH/Core/IO/MeshIO.hh \
    $$OPENMESH_PATH/Core/IO/OFFFormat.hh \
    $$OPENMESH_PATH/Core/IO/OMFormat.hh \
    $$OPENMESH_PATH/Core/IO/Options.hh \
    $$OPENMESH_PATH/Core/IO/SR_binary.hh \
    $$OPENMESH_PATH/Core/IO/SR_binary_spec.hh \
    $$OPENMESH_PATH/Core/IO/SR_binary_vector_of_bool.inl \
    $$OPENMESH_PATH/Core/IO/SR_binary_vector_of_fundamentals.inl \
    $$OPENMESH_PATH/Core/IO/SR_binary_vector_of_string.inl \
    $$OPENMESH_PATH/Core/IO/SR_rbo.hh \
    $$OPENMESH_PATH/Core/IO/SR_store.hh \
    $$OPENMESH_PATH/Core/IO/SR_types.hh \
    $$OPENMESH_PATH/Core/IO/StoreRestore.hh \
    $$OPENMESH_PATH/Core/IO/exporter/BaseExporter.hh \
    $$OPENMESH_PATH/Core/IO/exporter/ExporterT.hh \
    $$OPENMESH_PATH/Core/IO/importer/BaseImporter.hh \
    $$OPENMESH_PATH/Core/IO/importer/ImporterT.hh \
    $$OPENMESH_PATH/Core/IO/reader/BaseReader.hh \
    $$OPENMESH_PATH/Core/IO/reader/OBJReader.hh \
    $$OPENMESH_PATH/Core/IO/reader/OFFReader.hh \
    $$OPENMESH_PATH/Core/IO/reader/OMReader.hh \
    $$OPENMESH_PATH/Core/IO/reader/PLYReader.hh \
    $$OPENMESH_PATH/Core/IO/reader/STLReader.hh \
    $$OPENMESH_PATH/Core/IO/writer/BaseWriter.hh \
    $$OPENMESH_PATH/Core/IO/writer/OBJWriter.hh \
    $$OPENMESH_PATH/Core/IO/writer/OFFWriter.hh \
    $$OPENMESH_PATH/Core/IO/writer/OMWriter.hh \
    $$OPENMESH_PATH/Core/IO/writer/PLYWriter.hh \
    $$OPENMESH_PATH/Core/IO/writer/STLWriter.hh \
    $$OPENMESH_PATH/Core/IO/writer/VTKWriter.hh \
    $$OPENMESH_PATH/Core/Mesh/ArrayItems.hh \
    $$OPENMESH_PATH/Core/Mesh/ArrayKernel.hh \
    $$OPENMESH_PATH/Core/Mesh/AttribKernelT.hh \
    $$OPENMESH_PATH/Core/Mesh/Attributes.hh \
    $$OPENMESH_PATH/Core/Mesh/BaseKernel.hh \
    $$OPENMESH_PATH/Core/Mesh/BaseMesh.hh \
    $$OPENMESH_PATH/Core/Mesh/Casts.hh \
    $$OPENMESH_PATH/Core/Mesh/CirculatorsT.hh \
    $$OPENMESH_PATH/Core/Mesh/FinalMeshItemsT.hh \
    $$OPENMESH_PATH/Core/Mesh/Handles.hh \
    $$OPENMESH_PATH/Core/Mesh/IteratorsT.hh \
    $$OPENMESH_PATH/Core/Mesh/PolyConnectivity.hh \
    $$OPENMESH_PATH/Core/Mesh/PolyMeshT.hh \
    $$OPENMESH_PATH/Core/Mesh/PolyMesh_ArrayKernelT.hh \
    $$OPENMESH_PATH/Core/Mesh/Status.hh \
    $$OPENMESH_PATH/Core/Mesh/Traits.hh \
    $$OPENMESH_PATH/Core/Mesh/TriConnectivity.hh \
    $$OPENMESH_PATH/Core/Mesh/TriMeshT.hh \
    $$OPENMESH_PATH/Core/Mesh/TriMesh_ArrayKernelT.hh \
    $$OPENMESH_PATH/Core/System/OpenMeshDLLMacros.hh \
    $$OPENMESH_PATH/Core/System/compiler.hh \
    $$OPENMESH_PATH/Core/System/config.h \
    $$OPENMESH_PATH/Core/System/config.hh \
    $$OPENMESH_PATH/Core/System/mostream.hh \
    $$OPENMESH_PATH/Core/System/omstream.hh \
    $$OPENMESH_PATH/Core/Utils/AutoPropertyHandleT.hh \
    $$OPENMESH_PATH/Core/Utils/BaseProperty.hh \
    $$OPENMESH_PATH/Core/Utils/Endian.hh \
    $$OPENMESH_PATH/Core/Utils/GenProg.hh \
    $$OPENMESH_PATH/Core/Utils/Noncopyable.hh \
    $$OPENMESH_PATH/Core/Utils/Property.hh \
    $$OPENMESH_PATH/Core/Utils/PropertyContainer.hh \
    $$OPENMESH_PATH/Core/Utils/PropertyManager.hh \
    $$OPENMESH_PATH/Core/Utils/RandomNumberGenerator.hh \
    $$OPENMESH_PATH/Core/Utils/SingletonT.hh \
    $$OPENMESH_PATH/Core/Utils/color_cast.hh \
    $$OPENMESH_PATH/Core/Utils/vector_cast.hh \
    $$OPENMESH_PATH/Core/Utils/vector_traits.hh \

#Patchgen
SOURCES += \
    $$PATCHGEN_PATH/adjust_parameter.cpp \
    $$PATCHGEN_PATH/get_boundary_geometry.cpp \
    $$PATCHGEN_PATH/get_constraint_matrix.cpp \
    $$PATCHGEN_PATH/get_constraint_rhs.cpp \
    $$PATCHGEN_PATH/get_default_parameter.cpp \
    $$PATCHGEN_PATH/get_param_str.cpp \
    $$PATCHGEN_PATH/get_variable.cpp \
    $$PATCHGEN_PATH/get_variable_indicators.cpp \
    $$PATCHGEN_PATH/switch_pattern.cpp

HEADERS += \
    $$PATCHGEN_PATH/ILP.h \
    $$PATCHGEN_PATH/PatchParam.h \
    $$PATCHGEN_PATH/PatchVertexTag.h \
    $$PATCHGEN_PATH/PatchVertexTraits.h \
    $$PATCHGEN_PATH/Pattern.h \
    $$PATCHGEN_PATH/Pattern_2_0.h \
    $$PATCHGEN_PATH/Pattern_2_1.h \
    $$PATCHGEN_PATH/Pattern_3_0.h \
    $$PATCHGEN_PATH/Pattern_3_1.h \
    $$PATCHGEN_PATH/Pattern_3_2.h \
    $$PATCHGEN_PATH/Pattern_3_3.h \
    $$PATCHGEN_PATH/Pattern_4_0.h \
    $$PATCHGEN_PATH/Pattern_4_1.h \
    $$PATCHGEN_PATH/Pattern_4_2.h \
    $$PATCHGEN_PATH/Pattern_4_3.h \
    $$PATCHGEN_PATH/Pattern_4_4.h \
    $$PATCHGEN_PATH/Pattern_5_0.h \
    $$PATCHGEN_PATH/Pattern_5_1.h \
    $$PATCHGEN_PATH/Pattern_5_2.h \
    $$PATCHGEN_PATH/Pattern_5_3.h \
    $$PATCHGEN_PATH/Pattern_5_4.h \
    $$PATCHGEN_PATH/Pattern_6_0.h \
    $$PATCHGEN_PATH/Pattern_6_1.h \
    $$PATCHGEN_PATH/Pattern_6_2.h \
    $$PATCHGEN_PATH/Pattern_6_3.h \
    $$PATCHGEN_PATH/Pattern_all.h \
    $$PATCHGEN_PATH/Permutation.h \
    $$PATCHGEN_PATH/add_padding.h \
    $$PATCHGEN_PATH/decl.h \
    $$PATCHGEN_PATH/generate_subtopology.h \
    $$PATCHGEN_PATH/generate_topology.h

