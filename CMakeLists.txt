#--------------------------------------------------
cmake_minimum_required(VERSION 2.8)
IF (ParaView_SOURCE_DIR)
  INCLUDE_DIRECTORIES(
    ${VTK_INCLUDE_DIRS}
    BaseCluster.h
    HierarchicalCluster.h
    SpectralCluster.h
    KmeansCluster.h
    SelfTuningCluster.h
  )
ELSE (ParaView_SOURCE_DIR)
  FIND_PACKAGE(ParaView REQUIRED)
  INCLUDE(${PARAVIEW_USE_FILE})
ENDIF (ParaView_SOURCE_DIR)

ADD_PARAVIEW_PLUGIN(SelfTuningCluster "1.0"
  SERVER_MANAGER_XML SelfTuningCluster.xml
SERVER_MANAGER_SOURCES EntryCluster.cxx
	SOURCES BaseCluster.cpp HierarchicalCluster.cpp SpectralCluster.cpp KmeansCluster.cpp SelfTuningCluster.cpp
	
)

