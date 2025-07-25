# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

o2physics_add_dpl_workflow(pid-ml-producer
                           SOURCES pidMlProducer.cxx
                           JOB_POOL analysis
                           PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pid-ml-batch-eff-and-pur-producer
                           SOURCES pidMlBatchEffAndPurProducer.cxx
                           JOB_POOL analysis
                           PUBLIC_LINK_LIBRARIES O2::Framework ONNXRuntime::ONNXRuntime O2::CCDB O2Physics::DataModel
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(pid-ml-eff-and-pur-producer
                           SOURCES pidMlEffAndPurProducer.cxx
                           JOB_POOL analysis
                           PUBLIC_LINK_LIBRARIES O2::Framework ONNXRuntime::ONNXRuntime O2::CCDB O2Physics::DataModel
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(simple-apply-pid-onnx-model
                           SOURCES simpleApplyPidOnnxModel.cxx
                           JOB_POOL analysis
                           PUBLIC_LINK_LIBRARIES O2::Framework ONNXRuntime::ONNXRuntime O2::CCDB O2Physics::DataModel
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(simple-apply-pid-onnx-interface
                           SOURCES simpleApplyPidOnnxInterface.cxx
                           JOB_POOL analysis
                           PUBLIC_LINK_LIBRARIES O2::Framework ONNXRuntime::ONNXRuntime O2::CCDB O2Physics::DataModel
                           COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-pid
                  SOURCES qaPid.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2::CCDB
                  COMPONENT_NAME Analysis)

o2physics_add_dpl_workflow(qa-pid-ml
                  SOURCES qaPidMl.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore ONNXRuntime::ONNXRuntime
                  COMPONENT_NAME Analysis)
