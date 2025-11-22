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

/// \file skimEmcClusterConverter.cxx
/// \brief Converter task to convert SkimEMCClusters
/// \author M. Hemmer, marvin.hemmer@cern.ch

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

#include <algorithm>
#include <gsl/span>
#include <vector>

using namespace o2;
using namespace o2::framework;

void convertSpanToVector(std::vector<float>& destData, const gsl::span<const float>& srcData)
{
  destData.reserve(srcData.size());
  std::copy(srcData.begin(), srcData.end(), destData.begin());
}

// Converts SkimEMCClusters_000 into SkimEMCClusters_001
struct SkimEmcClusterConverter {
  Produces<aod::SkimEMCClusters_001> tableGammaEMCReco001;

  void process(aod::SkimEMCClusters_000 const& emcClusters)
  {
    std::vector<float> vDummy = {};
    std::vector<float> vPhi, vEta, vPt, vP;
    for (const auto& emcCluster : emcClusters) {
      // using convertSpanToVector is just a temporal solution, since right now tables return gsl::span
      // while filling a table needs std::span which can not be transformed. So going over std::vector
      // as a middle point is the current solution
      convertSpanToVector(vPhi, emcCluster.deltaPhi());
      convertSpanToVector(vEta, emcCluster.deltaEta());
      convertSpanToVector(vP, emcCluster.trackp());
      convertSpanToVector(vPt, emcCluster.trackpt());
      tableGammaEMCReco001(emcCluster.collisionId(), emcCluster.definition(), emcCluster.e(), emcCluster.eta(), emcCluster.phi(), emcCluster.m02(), emcCluster.nCells(), emcCluster.time(), emcCluster.isExotic(), vPhi, vEta, vP, vPt, vDummy, vDummy, vDummy, vDummy);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SkimEmcClusterConverter>(cfgc),
  };
}
