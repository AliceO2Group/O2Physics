o2physics_add_dpl_workflow(pid-feature-extractor
                           SOURCES pidFeatureExtractor.cxx
                           PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pid-onnx-inference
                           SOURCES pidOnnxInference.cxx
                           PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore O2Physics::MLCore
                           COMPONENT_NAME Analysis)
