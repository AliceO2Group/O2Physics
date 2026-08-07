o2physics_add_dpl_workflow(pid-feature-extractor
                           SOURCES pidFeatureExtractor.cxx
                           PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

# pid-onnx-inference temporarily disabled - it consumed aod::PidFeaturesData/
# aod::PidFeaturesMc, which no longer exist now that pidFeatureExtractor.cxx
# writes TFile/TTree instead of O2 AOD tables. Needs rework before
# re-enabling - see project notes.
# o2physics_add_dpl_workflow(pid-onnx-inference
#                            SOURCES pidOnnxInference.cxx
#                            PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
#                            COMPONENT_NAME Analysis)
