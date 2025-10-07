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
//
///
/// \file pseudoEffMFT.cxx
/// \brief  Task for calculating pseudo-efficiency of MFT disks; based on PWGDQ/Tasks/taskMFTTrkEfficiency.cxx
/// \author Gyula Bencedi, gyula.bencedi@cern.ch
/// \since  OCT 2025

#include "Functions.h"
#include "Index.h"
#include "bestCollisionTable.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/MathConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RuntimeError.h"
#include "Framework/runDataProcessing.h"
#include "MathUtils/Utils.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include <DataFormatsITSMFT/TimeDeadMap.h>
#include <Framework/HistogramRegistry.h>
#include <ITSMFTReconstruction/ChipMappingMFT.h>

#include "TPDGCode.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using namespace o2::aod::fwdtrack;
using namespace o2::constants::physics;
using namespace o2::constants::math;
using namespace pwgmm::mult;
using namespace o2::aod::rctsel;

namespace mft_xy
{
constexpr float kMx[936] = {
  -8.8, -8.8, 8.2, 8.2, 9.9, 9.9, -9.9, -9.9, -8.2, -8.2,
  8.8, 8.8, -7.1, -7.1, -7.1, -5.4, -5.4, -5.4, -3.7, -3.7,
  -3.7, -2, -2, -2, -0.3, -0.3, -0.3, 1.4, 1.4, 1.4,
  3.1, 3.1, 3.1, 4.8, 4.8, 4.8, 6.5, 6.5, 6.5, -6.5,
  -6.5, -6.5, -4.8, -4.8, -4.8, -3.1, -3.1, -3.1, -1.4, -1.4,
  -1.4, 0.3, 0.3, 0.3, 2, 2, 2, 3.7, 3.7, 3.7,
  5.4, 5.4, 5.4, 7.1, 7.1, 7.1, -8.8, -8.8, 8.2, 8.2,
  9.9, 9.9, -9.9, -9.9, -8.2, -8.2, 8.8, 8.8, -7.1, -7.1,
  -7.1, -5.4, -5.4, -5.4, -3.7, -3.7, -3.7, -2, -2, -2,
  -0.3, -0.3, -0.3, 1.4, 1.4, 1.4, 3.1, 3.1, 3.1, 4.8,
  4.8, 4.8, 6.5, 6.5, 6.5, -6.5, -6.5, -6.5, -4.8, -4.8,
  -4.8, -3.1, -3.1, -3.1, -1.4, -1.4, -1.4, 0.3, 0.3, 0.3,
  2, 2, 2, 3.7, 3.7, 3.7, 5.4, 5.4, 5.4, 7.1,
  7.1, 7.1, -10.5, -10.5, 9.9, 9.9, -9.9, -9.9, 10.5, 10.5,
  -8.8, -8.8, -8.8, -7.1, -7.1, -7.1, -2, -2, -2, -0.3,
  -0.3, -0.3, 1.4, 1.4, 1.4, 6.5, 6.5, 6.5, 8.2, 8.2,
  8.2, -8.2, -8.2, -8.2, -6.5, -6.5, -6.5, -1.4, -1.4, -1.4,
  0.3, 0.3, 0.3, 2, 2, 2, 7.1, 7.1, 7.1, 8.8,
  8.8, 8.8, -5.4, -5.4, -5.4, -5.4, -3.7, -3.7, -3.7, -3.7,
  3.1, 3.1, 3.1, 3.1, 4.8, 4.8, 4.8, 4.8, -4.8, -4.8,
  -4.8, -4.8, -3.1, -3.1, -3.1, -3.1, 3.7, 3.7, 3.7, 3.7,
  5.4, 5.4, 5.4, 5.4, -12.2, -12.2, -12.2, -10.5, -10.5, -10.5,
  9.9, 9.9, 9.9, 11.6, 11.6, 11.6, 13.3, 13.3, 13.3, -13.3,
  -13.3, -13.3, -11.6, -11.6, -11.6, -9.9, -9.9, -9.9, 10.5, 10.5,
  10.5, 12.2, 12.2, 12.2, -8.8, -8.8, -8.8, -8.8, -7.1, -7.1,
  -7.1, -7.1, -5.4, -5.4, -5.4, -5.4, -3.7, -3.7, -3.7, -3.7,
  -2, -2, -2, -2, -0.3, -0.3, -0.3, -0.3, 1.4, 1.4,
  1.4, 1.4, 3.1, 3.1, 3.1, 3.1, 4.8, 4.8, 4.8, 4.8,
  6.5, 6.5, 6.5, 6.5, 8.2, 8.2, 8.2, 8.2, -8.2, -8.2,
  -8.2, -8.2, -6.5, -6.5, -6.5, -6.5, -4.8, -4.8, -4.8, -4.8,
  -3.1, -3.1, -3.1, -3.1, -1.4, -1.4, -1.4, -1.4, 0.3, 0.3,
  0.3, 0.3, 2, 2, 2, 2, 3.7, 3.7, 3.7, 3.7,
  5.4, 5.4, 5.4, 5.4, 7.1, 7.1, 7.1, 7.1, 8.8, 8.8,
  8.8, 8.8, -13.9, -13.9, -13.9, -12.2, -12.2, -12.2, 11.6, 11.6,
  11.6, 13.3, 13.3, 13.3, -13.3, -13.3, -13.3, -11.6, -11.6, -11.6,
  12.2, 12.2, 12.2, 13.9, 13.9, 13.9, -10.5, -10.5, -10.5, -10.5,
  -8.8, -8.8, -8.8, -8.8, -3.7, -3.7, -3.7, -3.7, -2, -2,
  -2, -2, -0.3, -0.3, -0.3, -0.3, 1.4, 1.4, 1.4, 1.4,
  3.1, 3.1, 3.1, 3.1, 8.2, 8.2, 8.2, 8.2, 9.9, 9.9,
  9.9, 9.9, -9.9, -9.9, -9.9, -9.9, -8.2, -8.2, -8.2, -8.2,
  -3.1, -3.1, -3.1, -3.1, -1.4, -1.4, -1.4, -1.4, 0.3, 0.3,
  0.3, 0.3, 2, 2, 2, 2, 3.7, 3.7, 3.7, 3.7,
  8.8, 8.8, 8.8, 8.8, 10.5, 10.5, 10.5, 10.5, -7.1, -7.1,
  -7.1, -7.1, -7.1, -5.4, -5.4, -5.4, -5.4, -5.4, 4.8, 4.8,
  4.8, 4.8, 4.8, 6.5, 6.5, 6.5, 6.5, 6.5, -6.5, -6.5,
  -6.5, -6.5, -6.5, -4.8, -4.8, -4.8, -4.8, -4.8, 5.4, 5.4,
  5.4, 5.4, 5.4, 7.1, 7.1, 7.1, 7.1, 7.1, 8.8, 8.8,
  -8.2, -8.2, -9.9, -9.9, 9.9, 9.9, 8.2, 8.2, -8.8, -8.8,
  7.1, 7.1, 7.1, 5.4, 5.4, 5.4, 3.7, 3.7, 3.7, 2,
  2, 2, 0.3, 0.3, 0.3, -1.4, -1.4, -1.4, -3.1, -3.1,
  -3.1, -4.8, -4.8, -4.8, -6.5, -6.5, -6.5, 6.5, 6.5, 6.5,
  4.8, 4.8, 4.8, 3.1, 3.1, 3.1, 1.4, 1.4, 1.4, -0.3,
  -0.3, -0.3, -2, -2, -2, -3.7, -3.7, -3.7, -5.4, -5.4,
  -5.4, -7.1, -7.1, -7.1, 8.8, 8.8, -8.2, -8.2, -9.9, -9.9,
  9.9, 9.9, 8.2, 8.2, -8.8, -8.8, 7.1, 7.1, 7.1, 5.4,
  5.4, 5.4, 3.7, 3.7, 3.7, 2, 2, 2, 0.3, 0.3,
  0.3, -1.4, -1.4, -1.4, -3.1, -3.1, -3.1, -4.8, -4.8, -4.8,
  -6.5, -6.5, -6.5, 6.5, 6.5, 6.5, 4.8, 4.8, 4.8, 3.1,
  3.1, 3.1, 1.4, 1.4, 1.4, -0.3, -0.3, -0.3, -2, -2,
  -2, -3.7, -3.7, -3.7, -5.4, -5.4, -5.4, -7.1, -7.1, -7.1,
  10.5, 10.5, -9.9, -9.9, 9.9, 9.9, -10.5, -10.5, 8.8, 8.8,
  8.8, 7.1, 7.1, 7.1, 2, 2, 2, 0.3, 0.3, 0.3,
  -1.4, -1.4, -1.4, -6.5, -6.5, -6.5, -8.2, -8.2, -8.2, 8.2,
  8.2, 8.2, 6.5, 6.5, 6.5, 1.4, 1.4, 1.4, -0.3, -0.3,
  -0.3, -2, -2, -2, -7.1, -7.1, -7.1, -8.8, -8.8, -8.8,
  5.4, 5.4, 5.4, 5.4, 3.7, 3.7, 3.7, 3.7, -3.1, -3.1,
  -3.1, -3.1, -4.8, -4.8, -4.8, -4.8, 4.8, 4.8, 4.8, 4.8,
  3.1, 3.1, 3.1, 3.1, -3.7, -3.7, -3.7, -3.7, -5.4, -5.4,
  -5.4, -5.4, 12.2, 12.2, 12.2, 10.5, 10.5, 10.5, -9.9, -9.9,
  -9.9, -11.6, -11.6, -11.6, -13.3, -13.3, -13.3, 13.3, 13.3, 13.3,
  11.6, 11.6, 11.6, 9.9, 9.9, 9.9, -10.5, -10.5, -10.5, -12.2,
  -12.2, -12.2, 8.8, 8.8, 8.8, 8.8, 7.1, 7.1, 7.1, 7.1,
  5.4, 5.4, 5.4, 5.4, 3.7, 3.7, 3.7, 3.7, 2, 2,
  2, 2, 0.3, 0.3, 0.3, 0.3, -1.4, -1.4, -1.4, -1.4,
  -3.1, -3.1, -3.1, -3.1, -4.8, -4.8, -4.8, -4.8, -6.5, -6.5,
  -6.5, -6.5, -8.2, -8.2, -8.2, -8.2, 8.2, 8.2, 8.2, 8.2,
  6.5, 6.5, 6.5, 6.5, 4.8, 4.8, 4.8, 4.8, 3.1, 3.1,
  3.1, 3.1, 1.4, 1.4, 1.4, 1.4, -0.3, -0.3, -0.3, -0.3,
  -2, -2, -2, -2, -3.7, -3.7, -3.7, -3.7, -5.4, -5.4,
  -5.4, -5.4, -7.1, -7.1, -7.1, -7.1, -8.8, -8.8, -8.8, -8.8,
  13.9, 13.9, 13.9, 12.2, 12.2, 12.2, -11.6, -11.6, -11.6, -13.3,
  -13.3, -13.3, 13.3, 13.3, 13.3, 11.6, 11.6, 11.6, -12.2, -12.2,
  -12.2, -13.9, -13.9, -13.9, 10.5, 10.5, 10.5, 10.5, 8.8, 8.8,
  8.8, 8.8, 3.7, 3.7, 3.7, 3.7, 2, 2, 2, 2,
  0.3, 0.3, 0.3, 0.3, -1.4, -1.4, -1.4, -1.4, -3.1, -3.1,
  -3.1, -3.1, -8.2, -8.2, -8.2, -8.2, -9.9, -9.9, -9.9, -9.9,
  9.9, 9.9, 9.9, 9.9, 8.2, 8.2, 8.2, 8.2, 3.1, 3.1,
  3.1, 3.1, 1.4, 1.4, 1.4, 1.4, -0.3, -0.3, -0.3, -0.3,
  -2, -2, -2, -2, -3.7, -3.7, -3.7, -3.7, -8.8, -8.8,
  -8.8, -8.8, -10.5, -10.5, -10.5, -10.5, 7.1, 7.1, 7.1, 7.1,
  7.1, 5.4, 5.4, 5.4, 5.4, 5.4, -4.8, -4.8, -4.8, -4.8,
  -4.8, -6.5, -6.5, -6.5, -6.5, -6.5, 6.5, 6.5, 6.5, 6.5,
  6.5, 4.8, 4.8, 4.8, 4.8, 4.8, -5.4, -5.4, -5.4, -5.4,
  -5.4, -7.1, -7.1, -7.1, -7.1, -7.1};

constexpr float kMy[936] = {
  -1.7, -4.715, -1.7, -4.715, -1.7, -4.715, -1.7, -4.715, -1.7, -4.715,
  -1.7, -4.715, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -1.7, -4.715,
  -7.73, -3.546, -6.561, -9.576, -3.8, -6.815, -9.83, -3.706, -6.721, -9.736,
  -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -1.7,
  -4.715, -7.73, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -3.706, -6.721,
  -9.736, -3.8, -6.815, -9.83, -3.546, -6.561, -9.576, -1.7, -4.715, -7.73,
  -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -1.7, -4.715, -1.7, -4.715,
  -1.7, -4.715, -1.7, -4.715, -1.7, -4.715, -1.7, -4.715, -1.7, -4.715,
  -7.73, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -3.546, -6.561, -9.576,
  -3.8, -6.815, -9.83, -3.706, -6.721, -9.736, -1.7, -4.715, -7.73, -1.7,
  -4.715, -7.73, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -1.7, -4.715,
  -7.73, -1.7, -4.715, -7.73, -3.706, -6.721, -9.736, -3.8, -6.815, -9.83,
  -3.546, -6.561, -9.576, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -1.7,
  -4.715, -7.73, -1.7, -4.715, -1.7, -4.715, -1.7, -4.715, -1.7, -4.715,
  -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -3.55, -6.565, -9.58, -3.8,
  -6.815, -9.83, -3.71, -6.725, -9.74, -1.7, -4.715, -7.73, -1.7, -4.715,
  -7.73, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -3.71, -6.725, -9.74,
  -3.8, -6.815, -9.83, -3.55, -6.565, -9.58, -1.7, -4.715, -7.73, -1.7,
  -4.715, -7.73, -1.7, -4.715, -7.73, -10.745, -1.7, -4.715, -7.73, -10.745,
  -1.7, -4.715, -7.73, -10.745, -1.7, -4.715, -7.73, -10.745, -1.7, -4.715,
  -7.73, -10.745, -1.7, -4.715, -7.73, -10.745, -1.7, -4.715, -7.73, -10.745,
  -1.7, -4.715, -7.73, -10.745, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73,
  -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -1.7,
  -4.715, -7.73, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -1.7, -4.715,
  -7.73, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -10.745, -1.7, -4.715,
  -7.73, -10.745, -1.7, -4.715, -7.73, -10.745, -3.5, -6.515, -9.53, -12.545,
  -4.735, -7.75, -10.765, -13.78, -4.9, -7.915, -10.93, -13.945, -4.84, -7.855,
  -10.87, -13.885, -3.97, -6.985, -10, -13.015, -1.7, -4.715, -7.73, -10.745,
  -1.7, -4.715, -7.73, -10.745, -1.7, -4.715, -7.73, -10.745, -1.7, -4.715,
  -7.73, -10.745, -1.7, -4.715, -7.73, -10.745, -1.7, -4.715, -7.73, -10.745,
  -3.97, -6.985, -10, -13.015, -4.84, -7.855, -10.87, -13.885, -4.9, -7.915,
  -10.93, -13.945, -4.735, -7.75, -10.765, -13.78, -3.5, -6.515, -9.53, -12.545,
  -1.7, -4.715, -7.73, -10.745, -1.7, -4.715, -7.73, -10.745, -1.7, -4.715,
  -7.73, -10.745, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -1.7, -4.715,
  -7.73, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73,
  -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -1.7, -4.715, -7.73, -10.745,
  -1.7, -4.715, -7.73, -10.745, -4.125, -7.14, -10.155, -13.17, -5.155, -8.17,
  -11.185, -14.2, -5.3, -8.315, -11.33, -14.345, -5.245, -8.26, -11.275, -14.29,
  -4.5, -7.515, -10.53, -13.545, -1.7, -4.715, -7.73, -10.745, -1.7, -4.715,
  -7.73, -10.745, -1.7, -4.715, -7.73, -10.745, -1.7, -4.715, -7.73, -10.745,
  -4.5, -7.515, -10.53, -13.545, -5.245, -8.26, -11.275, -14.29, -5.3, -8.315,
  -11.33, -14.345, -5.155, -8.17, -11.185, -14.2, -4.125, -7.14, -10.155, -13.17,
  -1.7, -4.715, -7.73, -10.745, -1.7, -4.715, -7.73, -10.745, -1.7, -4.715,
  -7.73, -10.745, -13.76, -1.7, -4.715, -7.73, -10.745, -13.76, -1.7, -4.715,
  -7.73, -10.745, -13.76, -1.7, -4.715, -7.73, -10.745, -13.76, -1.7, -4.715,
  -7.73, -10.745, -13.76, -1.7, -4.715, -7.73, -10.745, -13.76, -1.7, -4.715,
  -7.73, -10.745, -13.76, -1.7, -4.715, -7.73, -10.745, -13.76, 1.7, 4.715,
  1.7, 4.715, 1.7, 4.715, 1.7, 4.715, 1.7, 4.715, 1.7, 4.715,
  1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 3.546,
  6.561, 9.576, 3.8, 6.815, 9.83, 3.706, 6.721, 9.736, 1.7, 4.715,
  7.73, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73,
  1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 3.706, 6.721, 9.736, 3.8,
  6.815, 9.83, 3.546, 6.561, 9.576, 1.7, 4.715, 7.73, 1.7, 4.715,
  7.73, 1.7, 4.715, 7.73, 1.7, 4.715, 1.7, 4.715, 1.7, 4.715,
  1.7, 4.715, 1.7, 4.715, 1.7, 4.715, 1.7, 4.715, 7.73, 1.7,
  4.715, 7.73, 1.7, 4.715, 7.73, 3.546, 6.561, 9.576, 3.8, 6.815,
  9.83, 3.706, 6.721, 9.736, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73,
  1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7,
  4.715, 7.73, 3.706, 6.721, 9.736, 3.8, 6.815, 9.83, 3.546, 6.561,
  9.576, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73,
  1.7, 4.715, 1.7, 4.715, 1.7, 4.715, 1.7, 4.715, 1.7, 4.715,
  7.73, 1.7, 4.715, 7.73, 3.55, 6.565, 9.58, 3.8, 6.815, 9.83,
  3.71, 6.725, 9.74, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7,
  4.715, 7.73, 1.7, 4.715, 7.73, 3.71, 6.725, 9.74, 3.8, 6.815,
  9.83, 3.55, 6.565, 9.58, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73,
  1.7, 4.715, 7.73, 10.745, 1.7, 4.715, 7.73, 10.745, 1.7, 4.715,
  7.73, 10.745, 1.7, 4.715, 7.73, 10.745, 1.7, 4.715, 7.73, 10.745,
  1.7, 4.715, 7.73, 10.745, 1.7, 4.715, 7.73, 10.745, 1.7, 4.715,
  7.73, 10.745, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7, 4.715,
  7.73, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73,
  1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7,
  4.715, 7.73, 1.7, 4.715, 7.73, 10.745, 1.7, 4.715, 7.73, 10.745,
  1.7, 4.715, 7.73, 10.745, 3.5, 6.515, 9.53, 12.545, 4.735, 7.75,
  10.765, 13.78, 4.9, 7.915, 10.93, 13.945, 4.84, 7.855, 10.87, 13.885,
  3.97, 6.985, 10, 13.015, 1.7, 4.715, 7.73, 10.745, 1.7, 4.715,
  7.73, 10.745, 1.7, 4.715, 7.73, 10.745, 1.7, 4.715, 7.73, 10.745,
  1.7, 4.715, 7.73, 10.745, 1.7, 4.715, 7.73, 10.745, 3.97, 6.985,
  10, 13.015, 4.84, 7.855, 10.87, 13.885, 4.9, 7.915, 10.93, 13.945,
  4.735, 7.75, 10.765, 13.78, 3.5, 6.515, 9.53, 12.545, 1.7, 4.715,
  7.73, 10.745, 1.7, 4.715, 7.73, 10.745, 1.7, 4.715, 7.73, 10.745,
  1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7,
  4.715, 7.73, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 1.7, 4.715,
  7.73, 1.7, 4.715, 7.73, 1.7, 4.715, 7.73, 10.745, 1.7, 4.715,
  7.73, 10.745, 4.125, 7.14, 10.155, 13.17, 5.155, 8.17, 11.185, 14.2,
  5.3, 8.315, 11.33, 14.345, 5.245, 8.26, 11.275, 14.29, 4.5, 7.515,
  10.53, 13.545, 1.7, 4.715, 7.73, 10.745, 1.7, 4.715, 7.73, 10.745,
  1.7, 4.715, 7.73, 10.745, 1.7, 4.715, 7.73, 10.745, 4.5, 7.515,
  10.53, 13.545, 5.245, 8.26, 11.275, 14.29, 5.3, 8.315, 11.33, 14.345,
  5.155, 8.17, 11.185, 14.2, 4.125, 7.14, 10.155, 13.17, 1.7, 4.715,
  7.73, 10.745, 1.7, 4.715, 7.73, 10.745, 1.7, 4.715, 7.73, 10.745,
  13.76, 1.7, 4.715, 7.73, 10.745, 13.76, 1.7, 4.715, 7.73, 10.745,
  13.76, 1.7, 4.715, 7.73, 10.745, 13.76, 1.7, 4.715, 7.73, 10.745,
  13.76, 1.7, 4.715, 7.73, 10.745, 13.76, 1.7, 4.715, 7.73, 10.745,
  13.76, 1.7, 4.715, 7.73, 10.745, 13.76};
} // namespace mft_xy

auto static constexpr kMinCharge = 3.f;
auto static constexpr kIntZero = 0;
auto static constexpr kZero = 0.f;
auto static constexpr kIntOne = 1;
auto static constexpr kNlayers = 10;
auto static constexpr kNhits = 15;
auto static constexpr kDisks = 5;

enum TrkSel {
  trkSelAll,
  trkSelNCls,
  trkSelChi2Ncl,
  trkSelEta,
  trkSelPhiCut,
  trkSelPt,
  trkSelCA,
  nTrkSel
};

enum TrkBestSel {
  trkBestSelAll,
  trkBestSelCollID,
  trkBestSelDCAxyCut,
  trkBestSelDCAzCut,
  trkBestSelWoAmbiguous,
  nTrkBestSel
};

enum OccupancyEst { TrkITS = 1,
                    Ft0C
};

struct PseudoEffMFT {

  HistogramRegistry registry{
    "registry",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  Configurable<bool> cfgDoIR{"cfgDoIR", false, "Flag to retrieve Interaction rate from CCDB"};
  Configurable<bool> cfgUseIRCut{"cfgUseIRCut", false, "Flag to cut on IR rate"};
  Configurable<bool> cfgIRCrashOnNull{"cfgIRCrashOnNull", false, "Flag to avoid CTP RateFetcher crash"};
  Configurable<std::string> cfgIRSource{"cfgIRSource", "ZNC hadronic", "Estimator of the interaction rate (Pb-Pb: ZNC hadronic)"};
  Configurable<bool> cfgUseEventkSel{"cfgUseEventkSel", false, "Flag to apply event selection"};
  Configurable<bool> cfgUseTrackSel{"cfgUseTrackSel", false, "Flag to apply track selection"};
  Configurable<bool> cfgUseParticleSel{"cfgUseParticleSel", false, "Flag to apply particle selection"};
  Configurable<bool> cfgUseCcdbForRun{"cfgUseCcdbForRun", false, "Get ccdb object based on run number instead of timestamp"};
  Configurable<bool> cfgRejectDeadChips{"cfgRejectDeadChips", true, "Reject tracks passing by dead chips per MFT layer "};
  Configurable<std::string> cfgDeadMapCcdbPath{"cfgDeadMapCcdbPath", "MFT/Calib/TimeDeadMap", "CCDB path for MFT dead map"};

  struct : ConfigurableGroup {
    ConfigurableAxis interactionRateBins{"interactionRateBins", {500, 0, 50}, "Binning for the interaction rate (kHz)"};
    ConfigurableAxis occupancyBins{"occupancyBins", {VARIABLE_WIDTH, 0.0f, 250.0f, 500.0f, 750.0f, 1000.0f, 1500.0f, 2000.0f, 3000.0f, 4500.0f, 6000.0f, 8000.0f, 10000.0f, 50000.0f}, "Occupancy"};
    ConfigurableAxis centralityBins{"centralityBins", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, "Centrality"};
    ConfigurableAxis irBins{"irBins", {500, 0, 50}, "Interaction rate (kHz)"};
    ConfigurableAxis ptBins{"ptBins", {11, -0.5, 10.5}, "pT binning (GeV/c)"};
    ConfigurableAxis phiBins{"phiBins", {63, 0., TwoPI}, "#varphi binning (rad)"};
    ConfigurableAxis etaBins{"etaBins", {20, -4., -2.}, "#eta binning"};
    Configurable<int> nEtaBins{"nEtaBins", 400, "Number of Eta bins"};
    Configurable<int> nPhiBins{"nPhiBins", 400, "Number of Phi bins"};
  } binOpt;

  struct : ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", false, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_fw", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } rctCuts;

  struct : ConfigurableGroup {
    Configurable<bool> usephiCut{"usephiCut", false, "use azimuthal angle cut"};
    Configurable<float> phiCut{"phiCut", 0.1f, "Cut on azimuthal angle of MFT tracks"};
    Configurable<float> minPhi{"minPhi", 0.f, ""};
    Configurable<float> maxPhi{"maxPhi", 6.2832, ""};
    Configurable<float> minEta{"minEta", -3.6f, ""};
    Configurable<float> maxEta{"maxEta", -2.5f, ""};
    Configurable<int> minNclusterMft{"minNclusterMft", 5, "minimum number of MFT clusters"};
    Configurable<bool> useChi2Cut{"useChi2Cut", false, "use track chi2 cut"};
    Configurable<float> maxChi2NCl{"maxChi2NCl", 1000.f, "maximum chi2 per MFT clusters"};
    Configurable<bool> usePtCut{"usePtCut", false, "use track pT cut"};
    Configurable<float> minPt{"minPt", 0., "minimum pT of the MFT tracks"};
    Configurable<bool> requireCA{"requireCA", false, "Use Cellular Automaton track-finding algorithm"};
    Configurable<bool> excludeAmbiguous{"excludeAmbiguous", false, "Exclude Ambiguous tracks"};
    Configurable<float> maxDCAxy{"maxDCAxy", 0.01f, "Cut on dca XY"};
    Configurable<float> maxDCAz{"maxDCAz", 0.01f, "Cut on dca Z"};
  } trackCuts;

  struct : ConfigurableGroup {
    Configurable<float> maxZvtx{"maxZvtx", 20.0f, "maximum cut on z-vtx (cm)"};
    Configurable<float> minZvtx{"minZvtx", -20.0f, "minimum cut on z-vtx (cm)"};
    Configurable<bool> useZDiffCut{"useZDiffCut", false, "use Zvtx reco-mc diff. cut"};
    Configurable<float> maxZvtxDiff{"maxZvtxDiff", 1.0f, "max allowed Z vtx difference for reconstruced collisions (cm)"};
    Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "require events with PV position along z consistent (within 1 cm) between PV reconstructed using tracks and PV using FT0 A-C time difference"};
    Configurable<bool> requireRejectSameBunchPileup{"requireRejectSameBunchPileup", true, "reject collisions in case of pileup with another collision in the same foundBC"};
    Configurable<bool> requireNoCollInTimeRangeStrict{"requireNoCollInTimeRangeStrict", true, " requireNoCollInTimeRangeStrict"};
    Configurable<bool> requireNoCollInRofStrict{"requireNoCollInRofStrict", false, "requireNoCollInRofStrict"};
    Configurable<bool> requireNoCollInRofStandard{"requireNoCollInRofStandard", false, "requireNoCollInRofStandard"};
    Configurable<bool> requireNoHighMultCollInPrevRof{"requireNoHighMultCollInPrevRof", false, "requireNoHighMultCollInPrevRof"};
    Configurable<bool> requireNoCollInTimeRangeStd{"requireNoCollInTimeRangeStd", true, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds"};
    Configurable<bool> requireNoCollInTimeRangeNarrow{"requireNoCollInTimeRangeNarrow", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds"};
    Configurable<uint> occupancyEstimator{"occupancyEstimator", 1, "Occupancy estimator: 1 = trackOccupancyInTimeRange, 2 = ft0cOccupancyInTimeRange"};
    Configurable<float> minOccupancy{"minOccupancy", -1, "minimum occupancy from neighbouring collisions"};
    Configurable<float> maxOccupancy{"maxOccupancy", -1, "maximum occupancy from neighbouring collisions"};
    Configurable<float> minIR{"minIR", -1, "minimum IR (kHz) collisions"};
    Configurable<float> maxIR{"maxIR", -1, "maximum IR (kHz) collisions"};
  } eventCuts;

  Service<o2::framework::O2DatabasePDG> pdg;
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPath{"ccdbPath", "MFT/Calib/TimeDeadMap", "ccdb path to MFT deadmap"};

  int nBCsPerOrbit = 3564;
  uint64_t mOrbit;
  uint64_t mPrevOrbit;
  o2::itsmft::TimeDeadMap* deadmap = nullptr;
  std::array<std::vector<int>, 10> chipsPerLayer{};
  std::array<TH2I*, 10> layerMasks;
  const o2::itsmft::ChipMappingMFT mapping;
  const std::array<o2::itsmft::MFTChipMappingData, 936> chipMap = mapping.getChipMappingData();

  float dX = 1.7;
  float dY = 3.015;
  std::array<float, 10> layersZ = {-45.3, -46.7, -48.6, -50.0, -52.4, -53.8, -67.7, -69.1, -76.1, -77.5};

  int mRunNumber{-1};
  uint64_t mSOR{0};
  float mMinSeconds{-1.};
  std::unordered_map<int, TH2*> gHadronicRate;
  ctpRateFetcher rateFetcher;
  TH2* gCurrentHadronicRate;
  RCTFlagsChecker rctChecker;

  /// @brief init function, definition of histograms
  void init(InitContext&)
  {
    const AxisSpec centralityAxis = {binOpt.centralityBins, "Centrality", "centrality axis"};
    const AxisSpec occupancyAxis = {binOpt.occupancyBins, "Occupancy", "occupancy axis"};
    const AxisSpec irAxis = {binOpt.interactionRateBins, "Interaction Rate", "IR axis"};
    const AxisSpec ptAxis = {binOpt.ptBins, "Pt axis (GeV/c)"};
    const AxisSpec phiAxis = {binOpt.phiBins, "#phi axis"};
    const AxisSpec etaAxis = {binOpt.etaBins, "#eta axis"};

    rctChecker.init(rctCuts.cfgEvtRCTFlagCheckerLabel, rctCuts.cfgEvtRCTFlagCheckerZDCCheck, rctCuts.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(ccdbNoLaterThan.value);
    ccdb->setFatalWhenNull(false);

    auto hev = registry.add<TH1>("Events/hEvtSel", "hEvtSel", HistType::kTH1F,
                                 {{15, -0.5f, +14.5f}});
    hev->GetXaxis()->SetBinLabel(1, "All collisions");
    hev->GetXaxis()->SetBinLabel(2, "Ev. sel.");
    hev->GetXaxis()->SetBinLabel(3, "kIsGoodZvtxFT0vsPV");
    hev->GetXaxis()->SetBinLabel(4, "NoSameBunchPileup");
    hev->GetXaxis()->SetBinLabel(5, "Z-vtx cut");
    hev->GetXaxis()->SetBinLabel(6, "kNoCollInTimeRangeStd");
    hev->GetXaxis()->SetBinLabel(7, "kNoCollInTimeRangeNarrow");
    hev->GetXaxis()->SetBinLabel(8, "kNoCollInTimeRangeStrict");
    hev->GetXaxis()->SetBinLabel(9, "kNoCollInRofStrict");
    hev->GetXaxis()->SetBinLabel(10, "kNoCollInRofStandard");
    hev->GetXaxis()->SetBinLabel(11, "kNoHighMultCollInPrevRof");
    hev->GetXaxis()->SetBinLabel(12, "Below min occup.");
    hev->GetXaxis()->SetBinLabel(13, "Above max occup.");
    hev->GetXaxis()->SetBinLabel(14, "RCT Flag Checker");

    registry.add("Tracks/hBestTrkSel", "Number of best tracks; Cut; #Tracks Passed Cut", {HistType::kTH1F, {{nTrkBestSel, -0.5, +nTrkBestSel - 0.5}}});
    registry.get<TH1>(HIST("Tracks/hBestTrkSel"))->GetXaxis()->SetBinLabel(trkBestSelAll + 1, "All");
    registry.get<TH1>(HIST("Tracks/hBestTrkSel"))->GetXaxis()->SetBinLabel(trkBestSelCollID + 1, "Assigned (ID>=0)");
    registry.get<TH1>(HIST("Tracks/hBestTrkSel"))->GetXaxis()->SetBinLabel(trkBestSelDCAxyCut + 1, "DCA xy cut");
    registry.get<TH1>(HIST("Tracks/hBestTrkSel"))->GetXaxis()->SetBinLabel(trkBestSelDCAzCut + 1, "DCA z cut");
    registry.get<TH1>(HIST("Tracks/hBestTrkSel"))->GetXaxis()->SetBinLabel(trkBestSelWoAmbiguous + 1, "No Ambiguous");

    registry.add("Tracks/hTrkSel", "Number of tracks; Cut; #Tracks Passed Cut", {HistType::kTH1F, {{nTrkSel, -0.5, +nTrkSel - 0.5}}});
    registry.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelAll + 1, "All");
    registry.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelNCls + 1, "Ncl cut");
    registry.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelChi2Ncl + 1, "#chi^{2}/Ncl cut");
    registry.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelEta + 1, "#eta cut");
    registry.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelPhiCut + 1, "#varphi cut");
    registry.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelPt + 1, "#it{p}_{T} cut");
    registry.get<TH1>(HIST("Tracks/hTrkSel"))->GetXaxis()->SetBinLabel(trkSelCA + 1, "Tracking algorithm (CA)");

    auto hBcSel = registry.add<TH1>("hBcSel", "hBcSel", HistType::kTH1F, {{3, -0.5f, +2.5f}});
    hBcSel->GetXaxis()->SetBinLabel(1, "Good BCs");
    hBcSel->GetXaxis()->SetBinLabel(2, "BCs with collisions");
    hBcSel->GetXaxis()->SetBinLabel(3, "BCs with pile-up/splitting");

    if (doprocessPseudoEffInclusive || doprocessPseudoEffMCInlusive ||
        doprocessPseudoEffCentFT0C || doprocessPseudoEffMCCentFT0C) {
      registry.add("Events/hInteractionRate", "; IR (kHz)", kTH1F, {irAxis});
      for (int i = 0; i < kNlayers; i++) {
        layerMasks[i] = new TH2I("", "", binOpt.nEtaBins, trackCuts.minEta, trackCuts.maxEta, binOpt.nPhiBins, -PI, +PI);
      }

      const AxisSpec axisMFtBitMap{1031, -0.5, 1030.5, "axisMFtBitMap"};
      registry.add("hMftBitMap", "hMftBitMap", {HistType::kTH1F, {axisMFtBitMap}});

      const AxisSpec axisNhits{15, -0.5, 15.5, ""};
      const char* elabels[15] = {"N12", "N10", "N02", "N34", "N30", "N04", "N56", "N50", "N06", "N78", "N70", "N08", "N910", "N90", "N010"};
      HistogramConfigSpec defaultNhitsEtaPtPhi({HistType::kTHnF, {{axisNhits}, {etaAxis}, {ptAxis}, {phiAxis}}});
      HistogramConfigSpec defaultNhitsCentEtaPtPhi({HistType::kTHnF, {{axisNhits}, {etaAxis}, {ptAxis}, {phiAxis}, {centralityAxis}}});

      if (doprocessPseudoEffInclusive) {
        registry.add("hEtaPtPhi", "hEtaPtPhi", defaultNhitsEtaPtPhi);
        auto hEtaPtPhi = registry.get<THn>(HIST("hEtaPtPhi"));
        for (int i = 0; i < kNhits; i++) {
          hEtaPtPhi->GetAxis(0)->SetBinLabel(i + 1, elabels[i]);
        }
      }

      if (doprocessPseudoEffCentFT0C) {
        registry.add("Cent/hEtaPtPhi", "hEtaPtPhi", defaultNhitsCentEtaPtPhi);
        auto hEtaPtPhi = registry.get<THn>(HIST("Cent/hEtaPtPhi"));
        for (int i = 0; i < kNhits; i++) {
          hEtaPtPhi->GetAxis(0)->SetBinLabel(i + 1, elabels[i]);
        }
      }

      if (doprocessPseudoEffMCInlusive) {
        registry.add("hPtRecVsGen", "; #it{p}_{T} (GeV/#it{c}); #it{p}_{T} (GeV/#it{c}) Gen", {HistType::kTH2F, {{ptAxis}, {ptAxis}}});
        registry.add("hEtaRecVsGen", "; #eta; #eta Gen", {HistType::kTH2F, {{etaAxis}, {etaAxis}}});
        registry.add("hPhiRecVsGen", "; #varphi; #varphi Gen", {HistType::kTH2F, {{phiAxis}, {phiAxis}}});
        registry.add("hEtaPtPhiGen", "hEtaPtPhiGen", defaultNhitsEtaPtPhi);
        auto hEtaPtPhiGen = registry.get<THn>(HIST("hEtaPtPhiGen"));
        for (int i = 0; i < kNhits; i++) {
          hEtaPtPhiGen->GetAxis(0)->SetBinLabel(i + 1, elabels[i]);
        }
        registry.add("hMftBitMapGen", "hMftBitMapGen", {HistType::kTH1F, {axisMFtBitMap}});
      }

      if (doprocessPseudoEffMCCentFT0C) {
        registry.add("Cent/hPtRecVsGen", "; #it{p}_{T} (GeV/#it{c}); #it{p}_{T} (GeV/#it{c}) Gen", {HistType::kTHnF, {{ptAxis}, {ptAxis}, {centralityAxis}}});
        registry.add("Cent/hEtaRecVsGen", "; #eta; #eta Gen", {HistType::kTHnF, {{etaAxis}, {etaAxis}, {centralityAxis}}});
        registry.add("Cent/hPhiRecVsGen", "; #varphi; #varphi Gen", {HistType::kTHnF, {{phiAxis}, {phiAxis}, {centralityAxis}}});
        registry.add("Cent/hEtaPtPhiGen", "hEtaPtPhiGen", defaultNhitsCentEtaPtPhi);
        auto hEtaPtPhiGen = registry.get<THn>(HIST("Cent/hEtaPtPhiGen"));
        for (int i = 0; i < kNhits; i++) {
          hEtaPtPhiGen->GetAxis(0)->SetBinLabel(i + 1, elabels[i]);
        }
      }
    }
  }

  /// Filters - tracks
  Filter filtTrkEta = (aod::fwdtrack::eta < trackCuts.maxEta) && (aod::fwdtrack::eta > trackCuts.minEta);
  Filter filtATrackID = (aod::fwdtrack::bestCollisionId >= 0);
  Filter filtATrackDCAxy = (nabs(aod::fwdtrack::bestDCAXY) < trackCuts.maxDCAxy);
  Filter filtATrackDCAz = (nabs(aod::fwdtrack::bestDCAZ) < trackCuts.maxDCAz);
  /// Filters - mc particles
  Filter primaries = (aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary && (aod::mcparticle::eta < trackCuts.maxEta) && (aod::mcparticle::eta > trackCuts.minEta);

  /// Joined tables
  using FullBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
  using CollBCs = soa::Join<aod::BCsWithTimestamps, aod::Run3MatchedToBCSparse>;
  using Colls = soa::Join<aod::Collisions, aod::EvSels>;
  using CollsCentFT0C = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::EvSels>;

  using MFTTracksLabeled = soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>;

  /// Filtered tables
  using FiltMftTracks = soa::Filtered<aod::MFTTracks>;
  using FiltMcMftTracks = soa::Filtered<MFTTracksLabeled>;
  using FiltParticles = soa::Filtered<aod::McParticles>;

  template <typename ObjType>
  ObjType* getForTsOrRun(std::string const& fullPath, int64_t timestamp, int runNumber)
  {
    if (cfgUseCcdbForRun) {
      return ccdb->getForRun<ObjType>(fullPath, runNumber);
    } else {
      return ccdb->getForTimeStamp<ObjType>(fullPath, timestamp);
    }
  }

  void initCCDB(CollBCs::iterator const& bc)
  {
    auto timestamp = bc.timestamp();
    auto runnumber = bc.runNumber();
    if (cfgRejectDeadChips) {
      deadmap = getForTsOrRun<o2::itsmft::TimeDeadMap>(cfgDeadMapCcdbPath, timestamp, runnumber);
      if (deadmap != nullptr) {
        LOGF(info, "Using deadmap for run %d", bc.runNumber());
      } else {
        LOGF(fatal, "DeadMap is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
    }
  }

  void decodeChipVector(std::vector<uint16_t>& inChips, std::array<std::vector<int>, 10>& outChipsPerLayer)
  {
    bool prevIsDead = false;
    uint16_t lastDead = 0;
    for (auto const& chip : inChips) {
      if (chip & static_cast<uint16_t>(0x8000)) {
        auto firstChip = chip - 0x8000;
        outChipsPerLayer[chipMap[firstChip].layer].push_back(firstChip);
        prevIsDead = true;
        lastDead = static_cast<uint16_t>(firstChip);
      } else {
        if (prevIsDead) {
          for (int i = 1; i < chip - lastDead; i++) {
            outChipsPerLayer[chipMap[i + lastDead].layer].push_back(i + lastDead);
          }
        }
        outChipsPerLayer[chipMap[chip].layer].push_back(chip);
        prevIsDead = false;
      }
    }
  }

  std::tuple<float, float, float, float> computeEtaPhiCoverage(float posX, float posY, float dX, float dY, float z)
  {
    float x1 = posX - dX / 2.0;
    float x2 = posX + dX / 2.0;
    float y1 = posY - dY / 2.0;
    float y2 = posY + dY / 2.0;

    std::vector<float> etas;
    std::vector<float> phis;
    etas.reserve(4);
    phis.reserve(4);

    for (auto const& x : {x1, x2}) {
      for (auto const& y : {y1, y2}) {
        float r = std::sqrt(x * x + y * y);
        float theta = std::atan2(r, z);
        float eta = -std::log(std::tan(theta / 2.0));
        float phi = std::atan2(y, x);
        etas.push_back(eta);
        phis.push_back(phi);
      }
    }

    float etaMin = *std::min_element(etas.begin(), etas.end());
    float etaMax = *std::max_element(etas.begin(), etas.end());
    float phiMin = *std::min_element(phis.begin(), phis.end());
    float phiMax = *std::max_element(phis.begin(), phis.end());

    RecoDecay::constrainAngle(phiMax - phiMin, 0.f);

    phiMin = *std::min_element(phis.begin(), phis.end());
    phiMax = *std::max_element(phis.begin(), phis.end());

    return std::make_tuple(etaMin, etaMax, phiMin, phiMax);
  }

  void applyChipToMask(TH2* mask, float etaMin, float etaMax, float phiMin, float phiMax)
  {
    if (!mask)
      return;

    int etaBinMin = mask->GetXaxis()->FindBin(etaMin);
    int etaBinMax = mask->GetXaxis()->FindBin(etaMax);
    int phiBinMin = mask->GetYaxis()->FindBin(phiMin);
    int phiBinMax = mask->GetYaxis()->FindBin(phiMax);

    for (int iEta = etaBinMin; iEta <= etaBinMax; ++iEta) {
      for (int iPhi = phiBinMin; iPhi <= phiBinMax; ++iPhi) {
        mask->SetBinContent(iEta, iPhi, 1.0);
      }
    }
  }

  void computeExclusionMap(int layer)
  {
    auto chips = chipsPerLayer[layer];
    float z = layersZ[layer];
    for (const auto& chip : chips) {
      float posX = mft_xy::kMx[chip];
      float posY = mft_xy::kMy[chip];
      auto [etaMin, etaMax, phi_min, phi_max] = computeEtaPhiCoverage(posX, posY, dX, dY, z);
      applyChipToMask(layerMasks[layer], etaMin, etaMax, phi_min, phi_max);
    }
  }

  bool isExcluded(float eta, float phi, int layer)
  {
    int binEta = layerMasks[layer]->GetXaxis()->FindBin(eta);
    int binPhi = layerMasks[layer]->GetYaxis()->FindBin(phi);
    bool isBad = (layerMasks[layer]->GetBinContent(binEta, binPhi) > 0);
    return isBad;
  }

  bool isHitAtDisk(uint16_t map, int ilayer)
  {
    LOGP(debug, " map %i --> %i", map, (map >> (ilayer * 6)) & 0x3F);
    return (map >> (ilayer * 6)) & 0x3F;
  }

  template <bool isCent, bool isMC = false>
  void fillHistos(float cent, float eta, float pt, float phi, uint16_t map)
  {
    if constexpr (!isMC) {
      registry.fill(HIST("hMftBitMap"), map);
    } else {
      registry.fill(HIST("hMftBitMapGen"), map);
    }

    bool iN[10];
    for (int ilayer = 0; ilayer < kNlayers; ilayer++) {
      iN[ilayer] = false;
      bool ishit = isHitAtDisk(map, ilayer);
      if (ishit) {
        iN[ilayer] = true;
      }
    }

    for (int i = 0; i < kDisks; i++) {
      if (iN[2 * i] && iN[2 * i + 1]) {
        if constexpr (!isMC) {
          if constexpr (isCent) {
            registry.get<THn>(HIST("Cent/hEtaPtPhi"))->Fill(3 * i, eta, pt, phi, cent);
          } else {
            registry.get<THn>(HIST("hEtaPtPhi"))->Fill(3 * i, eta, pt, phi);
          }
        } else {
          if constexpr (isCent) {
            registry.get<THn>(HIST("Cent/hEtaPtPhiGen"))->Fill(3 * i, eta, pt, phi, cent);
          } else {
            registry.get<THn>(HIST("hEtaPtPhiGen"))->Fill(3 * i, eta, pt, phi);
          }
        }
      }
      if (iN[2 * i] && (!iN[2 * i + 1]) && !isExcluded(eta, phi, 2 * i + 1)) {
        if constexpr (!isMC) {
          if constexpr (isCent) {
            registry.get<THn>(HIST("Cent/hEtaPtPhi"))->Fill(3 * i + 1, eta, pt, phi, cent);
          } else {
            registry.get<THn>(HIST("hEtaPtPhi"))->Fill(3 * i + 1, eta, pt, phi);
          }
        } else {
          if constexpr (isCent) {
            registry.get<THn>(HIST("Cent/hEtaPtPhiGen"))->Fill(3 * i + 1, eta, pt, phi, cent);
          } else {
            registry.get<THn>(HIST("hEtaPtPhiGen"))->Fill(3 * i + 1, eta, pt, phi);
          }
        }
      }
      if ((!iN[2 * i]) && iN[2 * i + 1] && !isExcluded(eta, phi, 2 * i)) {
        if constexpr (!isMC) {
          if constexpr (isCent) {
            registry.get<THn>(HIST("Cent/hEtaPtPhi"))->Fill(3 * i + 2, eta, pt, phi, cent);
          } else {
            registry.get<THn>(HIST("hEtaPtPhi"))->Fill(3 * i + 2, eta, pt, phi);
          }
        } else {
          if constexpr (isCent) {
            registry.get<THn>(HIST("Cent/hEtaPtPhiGen"))->Fill(3 * i + 2, eta, pt, phi, cent);
          } else {
            registry.get<THn>(HIST("hEtaPtPhiGen"))->Fill(3 * i + 2, eta, pt, phi);
          }
        }
      }
    }
  }

  template <bool fillHis = true, typename B>
  bool isBestTrackSelected(const B& besttrack)
  {
    if (fillHis) {
      registry.fill(HIST("Tracks/hBestTrkSel"), trkBestSelAll);
    }
    if (besttrack.bestCollisionId() < kIntZero) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hBestTrkSel"), trkBestSelCollID);
    }
    if (std::abs(besttrack.bestDCAXY()) >= trackCuts.maxDCAxy) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hBestTrkSel"), trkBestSelDCAxyCut);
    }
    if (std::abs(besttrack.bestDCAZ()) >= trackCuts.maxDCAxy) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hBestTrkSel"), trkBestSelDCAzCut);
    }
    if (trackCuts.excludeAmbiguous && besttrack.ambDegree() > kIntOne) {
      return false;
    }
    if (fillHis) {
      registry.fill(HIST("Tracks/hBestTrkSel"), trkBestSelWoAmbiguous);
    }
    return true;
  }

  template <bool fillHis = true, typename T>
  bool isTrackSelected(const T& track)
  {
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hTrkSel"), trkSelAll);
    }
    if (track.nClusters() < trackCuts.minNclusterMft) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hTrkSel"), trkSelNCls);
    }
    if (trackCuts.useChi2Cut) {
      float nclMft = std::max(2.0f * track.nClusters() - 5.0f, 1.0f);
      float mftChi2NCl = track.chi2() / nclMft;
      if (mftChi2NCl > trackCuts.maxChi2NCl)
        return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hTrkSel"), trkSelChi2Ncl);
    }
    if (track.eta() < trackCuts.minEta || track.eta() > trackCuts.maxEta) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hTrkSel"), trkSelEta);
    }
    if (trackCuts.usephiCut) {
      float phi = track.phi();
      o2::math_utils::bringTo02Pi(phi);
      if (phi < trackCuts.minPhi || trackCuts.maxPhi < phi) {
        return false;
      }
      if ((phi < trackCuts.phiCut) ||
          ((phi > PI - trackCuts.phiCut) && (phi < PI + trackCuts.phiCut)) ||
          (phi > TwoPI - trackCuts.phiCut) ||
          ((phi > ((PIHalf - 0.1) * PI) - trackCuts.phiCut) &&
           (phi < ((PIHalf - 0.1) * PI) + trackCuts.phiCut)))
        return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hTrkSel"), trkSelPhiCut);
    }
    if (trackCuts.usePtCut && track.pt() < trackCuts.minPt) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hTrkSel"), trkSelPt);
    }
    if (trackCuts.requireCA && !track.isCA()) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Tracks/hTrkSel"), trkSelCA);
    }
    return true;
  }

  template <typename P>
  bool isParticleSelected(P const& particle)
  {
    if (particle.eta() < trackCuts.minEta || particle.eta() > trackCuts.maxEta) {
      return false;
    }
    if (trackCuts.usephiCut) {
      float phi = particle.phi();
      o2::math_utils::bringTo02Pi(phi);
      if (phi < trackCuts.minPhi || trackCuts.maxPhi < phi) {
        return false;
      }
      if ((phi < trackCuts.phiCut) ||
          ((phi > PI - trackCuts.phiCut) && (phi < PI + trackCuts.phiCut)) ||
          (phi > TwoPI - trackCuts.phiCut) ||
          ((phi > ((PIHalf - 0.1) * PI) - trackCuts.phiCut) &&
           (phi < ((PIHalf - 0.1) * PI) + trackCuts.phiCut)))
        return false;
    }
    return true;
  }

  template <typename C>
  float getOccupancy(C const& collision, uint occEstimator)
  {
    switch (occEstimator) {
      case OccupancyEst::TrkITS:
        return collision.trackOccupancyInTimeRange();
      case OccupancyEst::Ft0C:
        return collision.ft0cOccupancyInTimeRange();
      default:
        LOG(fatal) << "No valid occupancy estimator ";
        break;
    }
    return -1.f;
  }

  void initHadronicRate(CollBCs::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
    if (gHadronicRate.find(mRunNumber) == gHadronicRate.end()) {
      auto runDuration = ccdb->getRunDuration(mRunNumber);
      mSOR = runDuration.first;
      mMinSeconds = std::floor(mSOR * 1.e-3);               /// round tsSOR to the highest integer lower than tsSOR
      float maxSec = std::ceil(runDuration.second * 1.e-3); /// round tsEOR to the lowest integer higher than tsEOR
      const AxisSpec axisSeconds{static_cast<int>((maxSec - mMinSeconds) / 20.f), 0, maxSec - mMinSeconds, "Seconds since SOR"};
      int hadronicRateBins = static_cast<int>(eventCuts.maxIR - eventCuts.minIR);
      gHadronicRate[mRunNumber] = registry.add<TH2>(Form("HadronicRate/%i", mRunNumber), ";Time since SOR (s);Hadronic rate (kHz)", kTH2D, {axisSeconds, {hadronicRateBins, eventCuts.minIR, eventCuts.maxIR}}).get();
    }
    gCurrentHadronicRate = gHadronicRate[mRunNumber];
  }

  template <bool fillHis = false, typename C>
  bool isGoodEvent(C const& collision)
  {
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 0);
    }
    if (!collision.sel8()) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 1);
    }
    if (eventCuts.requireIsGoodZvtxFT0VsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 2);
    }
    if (eventCuts.requireRejectSameBunchPileup && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 3);
    }
    if (collision.posZ() <= eventCuts.minZvtx || collision.posZ() >= eventCuts.maxZvtx) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 4);
    }
    if (eventCuts.requireNoCollInTimeRangeStd &&
        !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 5);
    }
    if (eventCuts.requireNoCollInTimeRangeNarrow &&
        !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 6);
    }
    if (eventCuts.requireNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 7);
    }
    if (eventCuts.requireNoCollInRofStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 8);
    }
    if (eventCuts.requireNoCollInRofStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 9);
    }
    if (eventCuts.requireNoHighMultCollInPrevRof && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 10);
    }
    if (eventCuts.minOccupancy >= 0 &&
        getOccupancy(collision, eventCuts.occupancyEstimator) <
          eventCuts.minOccupancy) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 11);
    }
    if (eventCuts.maxOccupancy >= 0 &&
        getOccupancy(collision, eventCuts.occupancyEstimator) >
          eventCuts.maxOccupancy) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 12);
    }

    if (rctCuts.requireRCTFlagChecker && !rctChecker(collision)) {
      return false;
    }
    if constexpr (fillHis) {
      registry.fill(HIST("Events/hEvtSel"), 13);
    }
    return true;
  }

  /// @brief Selection of charged particles
  /// @return true: charged; false: not charged
  bool isChrgParticle(int code)
  {
    auto p = pdg->GetParticle(code);
    auto charge = 0.;
    if (p != nullptr) {
      charge = p->Charge();
    }
    return std::abs(charge) >= kMinCharge;
  }

  /// @brief process function for general event statistics
  void processTagging(FullBCs const& bcs, CollsCentFT0C const& collisions)
  {
    std::vector<typename std::decay_t<decltype(collisions)>::iterator> cols;
    for (auto const& bc : bcs) {
      if ((bc.selection_bit(aod::evsel::kIsBBT0A) &&
           bc.selection_bit(aod::evsel::kIsBBT0C)) != 0) {
        registry.fill(HIST("hBcSel"), 0);
        cols.clear();
        for (auto const& collision : collisions) {
          if (collision.has_foundBC()) {
            if (collision.foundBCId() == bc.globalIndex()) {
              cols.emplace_back(collision);
            }
          } else if (collision.bcId() == bc.globalIndex()) {
            cols.emplace_back(collision);
          }
        }
        LOGP(debug, "BC {} has {} collisions", bc.globalBC(), cols.size());
        if (!cols.empty()) {
          registry.fill(HIST("hBcSel"), 1);
          if (cols.size() > 1) {
            registry.fill(HIST("hBcSel"), 2);
          }
        }
      }
    }
  }

  PROCESS_SWITCH(PseudoEffMFT, processTagging, "Collect event sample stats", true);

  template <typename C>
  void processPseudoEff(typename C::iterator const& collision,
                        FiltMftTracks const& /*tracks*/,
                        soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks,
                        CollBCs const& /*bcs*/)
  {
    auto bcFound = collision.template foundBC_as<CollBCs>();

    if (cfgDoIR) {
      initHadronicRate(bcFound);
      float ir = !cfgIRSource.value.empty() ? rateFetcher.fetch(ccdb.service, bcFound.timestamp(), bcFound.runNumber(), cfgIRSource, cfgIRCrashOnNull) * 1.e-3 : -1;
      registry.fill(HIST("Events/hInteractionRate"), ir);
      float seconds = bcFound.timestamp() * 1.e-3 - mMinSeconds;
      if (cfgUseIRCut && (ir < eventCuts.minIR || ir > eventCuts.maxIR)) { // cut on hadronic rate
        return;
      }
      gCurrentHadronicRate->Fill(seconds, ir);
    }

    for (auto const& atrack : besttracks) {
      if (cfgUseTrackSel && !isBestTrackSelected<true>(atrack)) {
        continue;
      }
      auto track = atrack.mfttrack_as<FiltMftTracks>();
      if (!track.has_collision()) {
        continue;
      }
      const auto& coll = track.collision_as<C>();
      if (cfgRejectDeadChips) {
        auto bc = coll.template bc_as<CollBCs>();
        int currentRun = bc.runNumber();
        if (mRunNumber != currentRun) {
          initCCDB(bc);
          mRunNumber = currentRun;
          auto orbits = deadmap->getEvolvingMapKeys();
        }
        if (mOrbit != (bc.globalBC() / nBCsPerOrbit)) {
          mOrbit = (bc.globalBC() / nBCsPerOrbit);
          std::vector<uint16_t> encodeChips;
          auto lowerOrbit = deadmap->getMapAtOrbit(mOrbit, encodeChips);
          if ((mOrbit - lowerOrbit) > mPrevOrbit) {
            for (int i = 0; i < kNlayers; i++) {
              chipsPerLayer[i].clear();
            }
            // for (auto& v : chipsPerLayer) {
            //   v.clear();
            // }
            // for (auto& h : layerMasks) {
            //   if (h)
            //     h->Reset("ICES");
            // }
            decodeChipVector(encodeChips, chipsPerLayer);
            for (int i = 0; i < kNlayers; i++) {
              computeExclusionMap(i);
            }
            mPrevOrbit = mOrbit - lowerOrbit;
          }
        }
      }

      if (cfgUseEventkSel && !isGoodEvent<true>(coll)) {
        continue;
      }
      if (cfgUseTrackSel && !isTrackSelected<true>(track)) {
        continue;
      }

      auto eta = track.eta();
      auto pt = track.pt();
      auto phi = track.phi();
      o2::math_utils::bringTo02Pi(phi);
      if (phi < kZero || TwoPI < phi) {
        continue;
      }
      auto map = track.mftClusterSizesAndTrackFlags();
      fillHistos<has_reco_cent<C>>(getRecoCent(collision), eta, pt, phi, map);
    }
  }

  void processPseudoEffInclusive(Colls::iterator const& collision,
                                 FiltMftTracks const& tracks,
                                 soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks,
                                 CollBCs const& bcs)
  {
    processPseudoEff<Colls>(collision, tracks, besttracks, bcs);
  }

  PROCESS_SWITCH(PseudoEffMFT, processPseudoEffInclusive, "Process PseudoEffMFT data/reco (inclusive)", true);

  void processPseudoEffCentFT0C(CollsCentFT0C::iterator const& collision,
                                FiltMftTracks const& tracks,
                                soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks,
                                CollBCs const& bcs)
  {
    processPseudoEff<CollsCentFT0C>(collision, tracks, besttracks, bcs);
  }

  PROCESS_SWITCH(PseudoEffMFT, processPseudoEffCentFT0C, "Process PseudoEffMFT data/reco (in FT0C centrality bins)", false);

  template <typename C, typename MC>
  void processPseudoEffMC(typename soa::Join<C, aod::McCollisionLabels>::iterator const& collision,
                          MC const& /*mccollisions*/,
                          FiltParticles const& /*particles*/,
                          FiltMcMftTracks const& /*tracks*/,
                          soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks)
  {
    if (cfgUseEventkSel && !isGoodEvent<false>(collision)) {
      return;
    }
    if (!collision.has_mcCollision()) {
      return;
    }
    auto mcCollision = collision.mcCollision();
    if (eventCuts.useZDiffCut) {
      if (std::abs(collision.posZ() - mcCollision.posZ()) > eventCuts.maxZvtxDiff) {
        return;
      }
    }

    for (auto const& atrack : besttracks) {
      if (cfgUseTrackSel && !isBestTrackSelected<false>(atrack)) {
        continue;
      }
      auto track = atrack.mfttrack_as<FiltMcMftTracks>();
      auto eta = track.eta();
      auto pt = track.pt();

      auto phi = track.phi();
      o2::math_utils::bringTo02Pi(phi);
      if (phi < kZero || TwoPI < phi) {
        continue;
      }

      auto map = track.mftClusterSizesAndTrackFlags();
      if (cfgUseTrackSel && !isTrackSelected<false>(track)) {
        continue;
      }
      if (!track.has_mcParticle()) {
        LOGF(warning, "No MC particle for track, skipping...");
        continue;
      }
      auto particle = track.template mcParticle_as<FiltParticles>();

      if (!isChrgParticle(particle.pdgCode())) {
        continue;
      }
      if (cfgUseParticleSel && !isParticleSelected(particle)) {
        continue;
      }

      auto partphi = particle.phi();
      o2::math_utils::bringTo02Pi(partphi);
      if (partphi < 0.f || TwoPI < partphi) {
        continue;
      }

      float cGen = getRecoCent(collision);
      if constexpr (has_reco_cent<C>) {
        registry.fill(HIST("Cent/hPtRecVsGen"), pt, particle.pt(), cGen);
        registry.fill(HIST("Cent/hEtaRecVsGen"), eta, particle.eta(), cGen);
        registry.fill(HIST("Cent/hPhiRecVsGen"), phi, partphi, cGen);
      } else {
        registry.fill(HIST("hPtRecVsGen"), pt, particle.pt());
        registry.fill(HIST("hEtaRecVsGen"), eta, particle.eta());
        registry.fill(HIST("hPhiRecVsGen"), phi, partphi);
      }
      fillHistos<has_reco_cent<C>, true>(cGen, particle.eta(), particle.pt(), partphi, map);
    }
  }

  void processPseudoEffMCInlusive(soa::Join<Colls, aod::McCollisionLabels>::iterator const& collision,
                                  aod::McCollisions const& mccollisions,
                                  FiltParticles const& particles,
                                  FiltMcMftTracks const& tracks,
                                  soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks)
  {
    processPseudoEffMC<Colls, aod::McCollisions>(collision, mccollisions, particles, tracks, besttracks);
  }

  PROCESS_SWITCH(PseudoEffMFT, processPseudoEffMCInlusive, "Process PseudoEffMFT mc generated (inclusive)", false);

  void processPseudoEffMCCentFT0C(soa::Join<CollsCentFT0C, aod::McCollisionLabels>::iterator const& collision,
                                  aod::McCollisions const& mccollisions,
                                  FiltParticles const& particles,
                                  FiltMcMftTracks const& tracks,
                                  soa::SmallGroups<aod::BestCollisionsFwd3d> const& besttracks)
  {
    processPseudoEffMC<CollsCentFT0C, aod::McCollisions>(collision, mccollisions, particles, tracks, besttracks);
  }

  PROCESS_SWITCH(PseudoEffMFT, processPseudoEffMCCentFT0C, "Process PseudoEffMFT mc generated (in FT0C centrality bins)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PseudoEffMFT>(cfgc)};
}
