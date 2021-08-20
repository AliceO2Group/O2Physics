o2physics_add_dpl_workflow(gammaconversions
                    SOURCES gammaConversions.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisDataModel O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)
o2physics_add_dpl_workflow(gammaconversionsmc
                    SOURCES gammaConversionsmc.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2::DetectorsBase O2Physics::AnalysisDataModel O2Physics::AnalysisCore O2::DetectorsVertexing
                    COMPONENT_NAME Analysis)
