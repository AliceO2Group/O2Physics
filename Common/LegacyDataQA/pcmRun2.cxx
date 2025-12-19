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

/// \file pcmRun2.cxx
/// \brief Analysis using PCM photons from Run2
/// \author M. Hemmer, marvin.hemmer@cern.ch

#include "Common/DataModel/PIDResponseTPC.h"

#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>

#include <array>
#include <cmath>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct PcmRun2 {
  struct : ConfigurableGroup {
    std::string prefix = "trackcuts";
    Configurable<float> tpcFindOverFoundable{"tpcFindOverFoundable", 0.6, "Ratio of number of found TPC clusters to the number of foundable TPC clusters a track must have at least."};
    Configurable<float> minPt{"minPt", 0.05, "Minimum pT cut for the tracks."};
  } trackcuts;

  struct : ConfigurableGroup {
    std::string prefix = "pidcuts";
    Configurable<float> minNSigmaEl{"minNSigmaEl", -3., "Minimum NSgimal electron allowed for V0 leg."};
    Configurable<float> maxNSigmaEl{"maxNSigmaEl", +4., "Maximum NSgimal electron allowed for V0 leg."};
    Configurable<bool> doPionRejection{"doPionRejection", true, "Flag to enable pion rejection based on TPC PID for V0 legs."};
    Configurable<float> minNSigmaPiLowP{"minNSigmaPiLowP", 1., "Minimum NSgimal pion to reject V0 legs for low 0.4 < pT/(GeV/c) < 3.5."};
    Configurable<float> minNSigmaPiHighP{"minNSigmaPiHighP", +0.5, "Minimum NSgimal pion to reject V0 legs for pT > 3.5 GeV/c."};
  } pidcuts;

  struct : ConfigurableGroup {
    std::string prefix = "v0cuts";
    Configurable<float> minPt{"minPt", 0.02, "Minimum pT cut for the V0."};
    Configurable<float> maxEta{"maxEta", 0.8, "Maximum absolut pseudorapidity cut for the V0."};
    Configurable<float> maxZconv{"maxZconv", 0.8, "Maximum absolut z conversion position cut for the V0."};
    Configurable<float> qTFactor{"qTFactor", 0.125, "qT < this * pT"};
    Configurable<float> cosP{"cosP", 0.85, "cos of pointing angle > this value."};
    Configurable<bool> rejectTooCloseV0{"rejectTooCloseV0", true, "."};
  } v0cuts;

  using TracksWithPID = soa::Join<o2::aod::FullTracks, o2::aod::pidTPCPi, o2::aod::pidTPCEl>;

  SliceCache cache;
  Preslice<o2::aod::Run2OTFV0s> perCollisionV0 = aod::run2::oftv0::collisionId;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  void init(InitContext&)
  {

    const AxisSpec thnAxisInvMass{400, 0.0, 0.8, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{100, 0., 20., "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisCent{20, 0., 100., "Centrality (%)"};
    const AxisSpec thnAxisCuts{5, -0.5, 4.5, "cuts"};
    const AxisSpec thnAxisNSgimaE{254, -6.35f, 6.35f, "n#sigma_{e}"};

    registry.add("hSameEvent2D", "hSameEvent2D", HistType::kTH2D, {thnAxisInvMass, thnAxisPt});
    registry.add("hNSigmaEPt", "hNSigmaEPt", HistType::kTH2D, {thnAxisNSgimaE, thnAxisPt});

    registry.add("hCuts", "hCuts", HistType::kTH1D, {thnAxisCuts});
    registry.add("hCutsPair", "hCutsPair", HistType::kTH1D, {thnAxisCuts});

  }; // end init

  template <typename CollisionIter, typename V0Iter>
  float getCosP(CollisionIter const& coll, V0Iter const& v0)
  {
    // V0 momentum
    float px = v0.px();
    float py = v0.py();
    float pz = v0.pz();

    // V0 decay position
    float xV0 = v0.x();
    float yV0 = v0.y();
    float zV0 = v0.z();

    // Primary vertex
    float xPV = coll.posX();
    float yPV = coll.posY();
    float zPV = coll.posZ();

    // Displacement vector r = V0 - PV
    float rx = xV0 - xPV;
    float ry = yV0 - yPV;
    float rz = zV0 - zPV;

    // Dot product p·r
    float dot = px * rx + py * ry + pz * rz;

    // Magnitudes |p| and |r|
    float magP = std::sqrt(px * px + py * py + pz * pz);
    float magR = std::sqrt(rx * rx + ry * ry + rz * rz);

    // Cosine of pointing angle
    return dot / (magP * magR);
  }

  template <typename V0Iter>
  bool checkLegs(V0Iter const& v0)
  {

    auto posLeg = v0.template posTrack_as<TracksWithPID>();
    auto negLeg = v0.template negTrack_as<TracksWithPID>();

    registry.fill(HIST("hNSigmaEPt"), posLeg.tpcNSigmaEl(), posLeg.pt());
    registry.fill(HIST("hNSigmaEPt"), negLeg.tpcNSigmaEl(), negLeg.pt());

    // first let's check the positive leg
    if (posLeg.pt() <= trackcuts.minPt.value) {
      registry.fill(HIST("hCuts"), 1);
      return false;
    }
    if (posLeg.tpcFoundOverFindableCls() <= trackcuts.tpcFindOverFoundable.value) {
      registry.fill(HIST("hCuts"), 1);
      return false;
    }
    if (posLeg.tpcNSigmaEl() <= pidcuts.minNSigmaEl.value || posLeg.tpcNSigmaEl() >= pidcuts.maxNSigmaEl.value) {
      registry.fill(HIST("hCuts"), 2);
      return false;
    }

    float minP = 0.4f; // minimum momentum for legs
    float midP = 3.5f; // momentum where we change NSigma window
    if (pidcuts.doPionRejection.value && posLeg.p() > minP) {
      if (posLeg.p() < midP && posLeg.tpcNSigmaPi() <= pidcuts.minNSigmaPiLowP.value) {
        registry.fill(HIST("hCuts"), 2);
        return false;
      } else if (posLeg.p() >= midP && posLeg.tpcNSigmaPi() <= pidcuts.minNSigmaPiHighP.value) {
        registry.fill(HIST("hCuts"), 2);
        return false;
      }
    }

    // second let's check the negative leg
    if (negLeg.pt() <= trackcuts.minPt.value) {
      registry.fill(HIST("hCuts"), 1);
      return false;
    }
    if (negLeg.tpcFoundOverFindableCls() <= trackcuts.tpcFindOverFoundable.value) {
      registry.fill(HIST("hCuts"), 1);
      return false;
    }
    if (negLeg.tpcNSigmaEl() <= pidcuts.minNSigmaEl.value || negLeg.tpcNSigmaEl() >= pidcuts.maxNSigmaEl.value) {
      registry.fill(HIST("hCuts"), 2);
      return false;
    }

    if (pidcuts.doPionRejection.value && negLeg.p() > minP) {
      if (negLeg.p() < midP && negLeg.tpcNSigmaPi() <= pidcuts.minNSigmaPiLowP.value) {
        registry.fill(HIST("hCuts"), 2);
        return false;
      } else if (negLeg.p() >= midP && negLeg.tpcNSigmaPi() <= pidcuts.minNSigmaPiHighP.value) {
        registry.fill(HIST("hCuts"), 2);
        return false;
      }
    }
    return true;
  }

  // Pi0 from EMCal
  void process(o2::aod::Collisions const& collisions, o2::aod::Run2OTFV0s const& v0s, TracksWithPID const& /*tracks*/)
  {
    // TODO:
    // - Add TooCloseToV0 cut: if two V0s are within dAngle < 0.02 and dR < 6. then remove the V0 with higher Chi2
    // - Everything here!
    // nothing yet

    const float zRslope = std::tan(2.f * std::atan(std::exp(-1.f * v0cuts.maxEta.value)));
    const float z0 = 7.f; // 7 cm
    const float maxZ = 10.f;
    const float minR = 5.f;
    const float maxR = 280.f;
    const float minRExclude = 55.f;
    const float maxRExclude = 72.f;

    for (const auto& collision : collisions) {
      if (std::fabs(collision.posZ()) > maxZ) {
        continue;
      }
      auto photonsPerCollision = v0s.sliceBy(perCollisionV0, collision.globalIndex());

      for (const auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photonsPerCollision, photonsPerCollision))) {
        registry.fill(HIST("hCutsPair"), 0);
        // next V0 cuts
        float cX1 = g1.x();
        float cY1 = g1.y();
        float cZ1 = g1.z();
        float cR1 = std::sqrt(std::pow(cX1, 2.) + std::pow(cY1, 2.));

        float cX2 = g2.x();
        float cY2 = g2.y();
        float cZ2 = g2.z();
        float cR2 = std::sqrt(std::pow(cX2, 2.) + std::pow(cY2, 2.));

        ROOT::Math::PxPyPzEVector v4Photon2(g2.px(), g2.py(), g2.pz(), g2.e());
        ROOT::Math::PxPyPzEVector v4Photon1(g1.px(), g1.py(), g1.pz(), g1.e());

        if (!checkLegs(g1)) {
          registry.fill(HIST("hCutsPair"), 1);
          continue;
        }

        if (!checkLegs(g2)) {
          registry.fill(HIST("hCutsPair"), 1);
          continue;
        }

        if (cR1 < minR || (cR1 > minRExclude && cR1 < maxRExclude) || cR1 > maxR) {
          registry.fill(HIST("hCutsPair"), 2);
          continue;
        }
        if (v4Photon1.Pt() < v0cuts.minPt.value) {
          registry.fill(HIST("hCutsPair"), 2);
          continue;
        }
        if (std::fabs(v4Photon1.Eta()) >= v0cuts.maxEta.value) {
          registry.fill(HIST("hCutsPair"), 2);
          continue;
        }
        if (std::fabs(cZ1) > v0cuts.maxZconv.value) {
          registry.fill(HIST("hCutsPair"), 2);
          continue;
        }
        if (cR1 <= ((std::fabs(cZ1) * zRslope) - z0)) {
          registry.fill(HIST("hCutsPair"), 2);
          continue;
        }
        if (g1.psiPair() >= (0.18f * std::exp(-0.055f * g1.chi2NDF()))) {
          registry.fill(HIST("hCutsPair"), 2);
          continue;
        }
        if (g1.qt() >= (v0cuts.qTFactor.value * v4Photon1.Pt())) {
          registry.fill(HIST("hCutsPair"), 2);
          continue;
        }
        if (getCosP(collision, g1) <= v0cuts.cosP.value) {
          registry.fill(HIST("hCutsPair"), 2);
          continue;
        }

        if (cR2 < minR || (cR2 > minRExclude && cR2 < maxRExclude) || cR2 > maxR) {
          registry.fill(HIST("hCutsPair"), 2);
          continue;
        }
        if (v4Photon2.Pt() < v0cuts.minPt.value) {
          registry.fill(HIST("hCutsPair"), 2);
          continue;
        }
        if (std::fabs(v4Photon2.Eta()) >= v0cuts.maxEta.value) {
          registry.fill(HIST("hCutsPair"), 2);
          continue;
        }
        if (std::fabs(cZ2) > v0cuts.maxZconv.value) {
          registry.fill(HIST("hCutsPair"), 2);
          continue;
        }
        if (cR2 <= ((std::fabs(cZ2) * zRslope) - z0)) {
          registry.fill(HIST("hCutsPair"), 2);
          continue;
        }
        if (g2.psiPair() >= (0.18f * std::exp(-0.055f * g2.chi2NDF()))) {
          registry.fill(HIST("hCutsPair"), 2);
          continue;
        }
        if (g2.qt() >= (v0cuts.qTFactor.value * v4Photon2.Pt())) {
          registry.fill(HIST("hCutsPair"), 2);
          continue;
        }
        if (getCosP(collision, g2) <= v0cuts.cosP.value) {
          registry.fill(HIST("hCutsPair"), 2);
          continue;
        }

        ROOT::Math::PxPyPzEVector vMeson = v4Photon1 + v4Photon2;
        registry.fill(HIST("hCutsPair"), 3);
        registry.fill(HIST("hSameEvent2D"), vMeson.M(), vMeson.Pt());
      } // end of loop over photon pairs
    } // end of loop over collision
  }
  PROCESS_SWITCH(PcmRun2, process, "Default process function", true);
}; // End struct PcmRun2

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PcmRun2>(cfgc)};
}
