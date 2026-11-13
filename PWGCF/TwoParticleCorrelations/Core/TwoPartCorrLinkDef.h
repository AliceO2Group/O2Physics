// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file TwoPartCorrLinkDef.h
/// \brief Tow particle correlations core library external interface definition
/// \author victor.gonzalez.sebastian@gmail.com

// NOLINT(build/header_guard)

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class o2::analysis::pwgcf::CutBrick < float> + ;
#pragma link C++ class o2::analysis::pwgcf::CutBrick < int> + ;
#pragma link C++ class o2::analysis::pwgcf::CutBrickLimit < float> + ;
#pragma link C++ class o2::analysis::pwgcf::CutBrickLimit < int> + ;
#pragma link C++ class o2::analysis::pwgcf::CutBrickFnLimit < float> + ;
#pragma link C++ class o2::analysis::pwgcf::CutBrickFnLimit < int> + ;
#pragma link C++ class o2::analysis::pwgcf::CutBrickThreshold < float> + ;
#pragma link C++ class o2::analysis::pwgcf::CutBrickThreshold < int> + ;
#pragma link C++ class o2::analysis::pwgcf::CutBrickFnThreshold < float> + ;
#pragma link C++ class o2::analysis::pwgcf::CutBrickFnThreshold < int> + ;
#pragma link C++ class o2::analysis::pwgcf::CutBrickExtToRange < float> + ;
#pragma link C++ class o2::analysis::pwgcf::CutBrickExtToRange < int> + ;
#pragma link C++ class o2::analysis::pwgcf::CutBrickFnExtToRange < float> + ;
#pragma link C++ class o2::analysis::pwgcf::CutBrickFnExtToRange < int> + ;
#pragma link C++ class o2::analysis::pwgcf::CutBrickRange < float> + ;
#pragma link C++ class o2::analysis::pwgcf::CutBrickRange < int> + ;
#pragma link C++ class o2::analysis::pwgcf::CutBrickFnRange < float> + ;
#pragma link C++ class o2::analysis::pwgcf::CutBrickFnRange < int> + ;
#pragma link C++ class o2::analysis::pwgcf::CutBrickSelectorMultipleRanges < float> + ;
#pragma link C++ class o2::analysis::pwgcf::CutBrickSelectorMultipleRanges < int> + ;
#pragma link C++ class o2::analysis::pwgcf::CutWithVariations < float> + ;
#pragma link C++ class o2::analysis::pwgcf::CutWithVariations < int> + ;
#pragma link C++ class o2::analysis::pwgcf::SpecialCutBrick + ;
#pragma link C++ class o2::analysis::pwgcf::TrackSelectionBrick + ;
#pragma link C++ class o2::analysis::pwgcf::SelectionFilterAndAnalysis + ;
#pragma link C++ class o2::analysis::pwgcf::TrackSelectionConfigurable + ;
#pragma link C++ class o2::analysis::pwgcf::TrackSelectionFilterAndAnalysis + ;
#pragma link C++ class o2::analysis::pwgcf::PIDSelectionConfigurable + ;
#pragma link C++ class o2::analysis::pwgcf::PIDSelectionFilterAndAnalysis + ;
#pragma link C++ class o2::analysis::pwgcf::EventSelectionConfigurable + ;
#pragma link C++ class o2::analysis::pwgcf::EventSelectionFilterAndAnalysis::ComplexBrickHelper + ;
#pragma link C++ class o2::analysis::pwgcf::EventSelectionFilterAndAnalysis::MultiplicityBrick + ;
#pragma link C++ class o2::analysis::pwgcf::EventSelectionFilterAndAnalysis::PileUpRejBrick + ;
#pragma link C++ class o2::analysis::pwgcf::EventSelectionFilterAndAnalysis + ;
#pragma link C++ class o2::analysis::pwgcf::FilterAndAnalysisFramework + ;
#pragma link C++ class DihadronContainer + ;
