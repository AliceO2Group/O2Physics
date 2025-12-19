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
/// \file selectionStudyTable.cxx
/// \brief Produces tables for centrality selection bias studies
///
/// \author ALICE
///

#include "Common/DataModel/SelectionStudyTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <cmath>
#include <cstdlib>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct SelectionStudyTable {
  Produces<aod::PIDPts> pidpts;

  // could be done in a vector of vectors
  // left for future iteration
  std::vector<float> ptpi;
  std::vector<float> ptka;
  std::vector<float> ptpr;
  std::vector<float> ptk0;
  std::vector<float> ptla;
  std::vector<float> ptxi;
  std::vector<float> ptom;
  std::vector<float> ptph;
  std::vector<float> ptks;
  std::vector<float> ptd;
  std::vector<float> ptlc;
  std::vector<float> ptjp;

  void init(InitContext&)
  {
  }

  void process(aod::McCollision const&, aod::McParticles const& mcParticles)
  {
    ptpi.clear();
    ptka.clear();
    ptpr.clear();
    ptk0.clear();
    ptla.clear();
    ptxi.clear();
    ptom.clear();
    ptph.clear();
    ptks.clear();
    ptd.clear();
    ptlc.clear();
    ptjp.clear();
    for (auto const& mcPart : mcParticles) {
      if (std::fabs(mcPart.y()) > 0.5) {
        continue; // only do midrapidity particles
      }

      // handle resonances first to make sure phys prim crit does not reject them
      if (mcPart.pdgCode() == 333) {
        ptph.push_back(mcPart.pt());
      }
      if (std::abs(mcPart.pdgCode()) == 313) {
        ptks.push_back(mcPart.pt());
      }

      // resonances handled, move to primaries
      if (!mcPart.isPhysicalPrimary()) {
        continue;
      }
      if (std::abs(mcPart.pdgCode()) == 211) {
        ptpi.push_back(mcPart.pt());
      }
      if (std::abs(mcPart.pdgCode()) == 321) {
        ptka.push_back(mcPart.pt());
      }
      if (std::abs(mcPart.pdgCode()) == 2212) {
        ptpr.push_back(mcPart.pt());
      }
      if (std::abs(mcPart.pdgCode()) == 310) {
        ptk0.push_back(mcPart.pt());
      }
      if (std::abs(mcPart.pdgCode()) == 3122) {
        ptla.push_back(mcPart.pt());
      }
      if (std::abs(mcPart.pdgCode()) == 3312) {
        ptxi.push_back(mcPart.pt());
      }
      if (std::abs(mcPart.pdgCode()) == 3334) {
        ptom.push_back(mcPart.pt());
      }
      if (std::abs(mcPart.pdgCode()) == 3334) {
        ptom.push_back(mcPart.pt());
      }
      // inclusive HF for now
      if (std::abs(mcPart.pdgCode()) == 421) {
        ptd.push_back(mcPart.pt());
      }
      if (std::abs(mcPart.pdgCode()) == 4122) {
        ptd.push_back(mcPart.pt());
      }
      if (std::abs(mcPart.pdgCode()) == 443) {
        ptjp.push_back(mcPart.pt());
      }
    }

    pidpts(ptpi, ptka, ptpr, ptk0, ptla, ptxi, ptom, ptph, ptks, ptd, ptlc, ptjp);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<SelectionStudyTable>(cfgc)};
}
