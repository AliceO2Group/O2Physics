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
///
/// \brief this is a starting point for the Resonances tutorial
/// \author sourav kundu
/// \since 02/11/2023

#include "PWGLF/DataModel/ReducedDoublePhiTables.h"

#include "Common/Core/trackUtilities.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include <Framework/Configurable.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TVector2.h>

#include <fairlogger/Logger.h>

#include <iostream>
#include <iterator>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct doublephimeson {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> strategyPID1{"strategyPID1", 0, "PID strategy 1"};
  Configurable<int> strategyPID2{"strategyPID2", 0, "PID strategy 2"};
  Configurable<float> daughterDeltaR{"daughterDeltaR", 0.0, "delta R of daughter"};
  Configurable<float> minPhiMass1{"minPhiMass1", 1.01, "Minimum phi mass1"};
  Configurable<float> maxPhiMass1{"maxPhiMass1", 1.03, "Maximum phi mass1"};
  Configurable<float> minPhiPt{"minPhiPt", 0, "Minimum phi Pt"};
  Configurable<float> maxPhiPt{"maxPhiPt", 100, "Maximum phi Pt"};
  Configurable<float> minPhiMass2{"minPhiMass2", 1.01, "Minimum phi mass2"};
  Configurable<float> maxPhiMass2{"maxPhiMass2", 1.03, "Maximum phi mass2"};
  Configurable<float> minExoticMass{"minExoticMass", 2.0, "Minimum Exotic mass"};
  Configurable<float> maxExoticMass{"maxExoticMass", 3.6, "Maximum Exotic mass"};
  Configurable<bool> additionalEvsel{"additionalEvsel", false, "Additional event selection"};
  Configurable<bool> isDeep{"isDeep", true, "Store deep angle"};
  Configurable<float> cutMinNsigmaTPC{"cutMinNsigmaTPC", -2.5, "nsigma cut TPC"};
  Configurable<float> cutNsigmaTPC{"cutNsigmaTPC", 3.0, "nsigma cut TPC"};
  Configurable<float> cutNsigmaTOF{"cutNsigmaTOF", 3.0, "nsigma cut TOF"};
  Configurable<float> momTOFCut{"momTOFCut", 1.8, "minimum pT cut for madnatory TOF"};
  Configurable<float> maxKaonPt{"maxKaonPt", 100.0, "maximum kaon pt cut"};
  // Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 1, "Number of events to mix"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {10, -10, 10}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0, 20.0, 40.0, 60.0, 80.0, 500.0}, "Mixing bins - number of contributor"};

  // THnsparse bining
  ConfigurableAxis configThnAxisPtCorr{"configThnAxisPtCorr", {1000, 0.0, 100}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {1500, 2.0, 3.5}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisInvMassPhi{"configThnAxisInvMassPhi", {20, 1.01, 1.03}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisInvMassDeltaPhi{"configThnAxisInvMassDeltaPhi", {80, 0.0, 0.08}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisDaugherPt{"configThnAxisDaugherPt", {25, 0.0, 50.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {40, 0.0, 20.}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisKstar{"configThnAxisKstar", {200, 0.0, 2.0}, "#it{k}^{*} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisDeltaR{"configThnAxisDeltaR", {200, 0.0, 2.0}, "#it{k}^{*} (GeV/#it{c})"};
  ConfigurableAxis configThnAxisCosTheta{"configThnAxisCosTheta", {160, 0.0, 3.2}, "cos #theta{*}"};
  ConfigurableAxis configThnAxisNumPhi{"configThnAxisNumPhi", {101, -0.5, 100.5}, "cos #theta{*}"};
  ConfigurableAxis configThnAxisDeltaPt{"configThnAxisDeltaPt", {100, 0.0, 1.0}, "delta pt"};
  Configurable<float> maxDeltaMPhi{"maxDeltaMPhi", 0.01f, "Delta-m cut on the two phi masses: sqrt((m1-mPDG)^2 + (m2-mPDG)^2) < maxDeltaMPhi (GeV/c^2)"};

  // Initialize the ananlysis task
  void init(o2::framework::InitContext&)
  {
    // register histograms
    histos.add("hnsigmaTPCKaonPlusBefore", "hnsigmaTPCKaonPlusBefore", kTH2F, {{1000, -3.0, 3.0f}, {100, 0.0f, 10.0f}});
    histos.add("hnsigmaTPCKaonMinusBefore", "hnsigmaTPCKaonMinusBefore", kTH2F, {{1000, -3.0, 3.0f}, {100, 0.0f, 10.0f}});
    histos.add("hnsigmaTPCTOFKaonBefore", "hnsigmaTPCTOFKaonBefore", kTH3F, {{500, -3.0, 3.0f}, {500, -3.0, 3.0f}, {100, 0.0f, 10.0f}});
    histos.add("hnsigmaTPCKaonPlus", "hnsigmaTPCKaonPlus", kTH2F, {{1000, -3.0, 3.0f}, {100, 0.0f, 10.0f}});
    histos.add("hnsigmaTPCKaonMinus", "hnsigmaTPCKaonMinus", kTH2F, {{1000, -3.0, 3.0f}, {100, 0.0f, 10.0f}});
    histos.add("hnsigmaTPCTOFKaon", "hnsigmaTPCTOFKaon", kTH3F, {{500, -3.0, 3.0f}, {500, -3.0, 3.0f}, {100, 0.0f, 10.0f}});
    histos.add("hPhiMass", "hPhiMass", kTH2F, {{40, 1.0, 1.04f}, {100, 0.0f, 10.0f}});
    histos.add("hPhiMass2", "hPhiMass2", kTH2F, {{40, 1.0, 1.04f}, {40, 1.0f, 1.04f}});
    histos.add("hkPlusDeltaetaDeltaPhi", "hkPlusDeltaetaDeltaPhi", kTH2F, {{400, -2.0, 2.0}, {640, -2.0 * TMath::Pi(), 2.0 * TMath::Pi()}});
    histos.add("hkMinusDeltaetaDeltaPhi", "hkMinusDeltaetaDeltaPhi", kTH2F, {{400, -2.0, 2.0}, {640, -2.0 * TMath::Pi(), 2.0 * TMath::Pi()}});
    histos.add("hDeltaRkaonplus", "hDeltaRkaonplus", kTH1F, {{800, 0.0, 8.0}});
    histos.add("hDeltaRkaonminus", "hDeltaRkaonminus", kTH1F, {{800, 0.0, 8.0}});
    histos.add("hPtCorrelation", "hPtCorrelation", kTH2F, {{400, 0.0, 40.0}, {5000, 0.0, 100.0}});
    const AxisSpec thnAxisdeltapt{configThnAxisDeltaPt, "Delta pt"};
    const AxisSpec thnAxisInvMass{configThnAxisInvMass, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{configThnAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisInvMassPhi{configThnAxisInvMassPhi, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisInvMassDeltaPhi{configThnAxisInvMassDeltaPhi, "#it{M} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisDeltaR{configThnAxisDeltaR, "#Delta R)"};
    const AxisSpec thnAxisCosTheta{configThnAxisCosTheta, "cos #theta"};
    const AxisSpec thnAxisNumPhi{configThnAxisNumPhi, "Number of phi meson"};
    const AxisSpec thnAxisPtCorr{configThnAxisPtCorr, "Pt Corr var"};

    histos.add("SEMassUnlike", "SEMassUnlike", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisDeltaR, thnAxisPt, thnAxisDeltaR, thnAxisInvMassDeltaPhi, thnAxisPtCorr});
    // histos.add("SEMassLike", "SEMassLike", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisDeltaR, thnAxisInvMassPhi, thnAxisInvMassPhi, thnAxisNumPhi});
    histos.add("MEMassUnlike", "MEMassUnlike", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisDeltaR, thnAxisInvMassDeltaPhi});
  }

  // get kstar
  TLorentzVector trackSum, PartOneCMS, PartTwoCMS, trackRelK;
  float getkstar(const TLorentzVector part1,
                 const TLorentzVector part2)
  {
    // const TLorentzVector trackSum = part1 + part2;
    trackSum = part1 + part2;
    const float beta = trackSum.Beta();
    const float betax = beta * std::cos(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betay = beta * std::sin(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betaz = beta * std::cos(trackSum.Theta());
    // TLorentzVector PartOneCMS(part1);
    // TLorentzVector PartTwoCMS(part2);
    PartOneCMS.SetXYZM(part1.Px(), part1.Py(), part1.Pz(), part1.M());
    PartTwoCMS.SetXYZM(part2.Px(), part2.Py(), part2.Pz(), part2.M());
    const ROOT::Math::Boost boostPRF = ROOT::Math::Boost(-betax, -betay, -betaz);
    PartOneCMS = boostPRF(PartOneCMS);
    PartTwoCMS = boostPRF(PartTwoCMS);
    // const TLorentzVector trackRelK = PartOneCMS - PartTwoCMS;
    trackRelK = PartOneCMS - PartTwoCMS;
    return 0.5 * trackRelK.P();
  }

  float deepangle2(const ROOT::Math::PtEtaPhiMVector candidate1,
                   const ROOT::Math::PtEtaPhiMVector candidate2)
  {
    double pt1, pt2, pz1, pz2, p1, p2, angle;
    pt1 = candidate1.Pt();
    pt2 = candidate2.Pt();
    pz1 = candidate1.Pz();
    pz2 = candidate2.Pz();
    p1 = candidate1.P();
    p2 = candidate2.P();
    angle = TMath::ACos((pt1 * pt2 + pz1 * pz2) / (p1 * p2));
    return angle;
  }

  float deepangle(const TLorentzVector candidate1,
                  const TLorentzVector candidate2)
  {
    double pt1, pt2, pz1, pz2, p1, p2, angle;
    pt1 = candidate1.Pt();
    pt2 = candidate2.Pt();
    pz1 = candidate1.Pz();
    pz2 = candidate2.Pz();
    p1 = candidate1.P();
    p2 = candidate2.P();
    angle = TMath::ACos((pt1 * pt2 + pz1 * pz2) / (p1 * p2));
    return angle;
  }

  // get cosTheta
  TLorentzVector daughterCMS;
  ROOT::Math::XYZVector threeVecDauCM, threeVecMother;
  float getCosTheta(const TLorentzVector mother,
                    const TLorentzVector daughter)
  {
    threeVecMother = mother.Vect();
    const float beta = mother.Beta();
    const float betax = beta * std::cos(mother.Phi()) * std::sin(mother.Theta());
    const float betay = beta * std::sin(mother.Phi()) * std::sin(mother.Theta());
    const float betaz = beta * std::cos(mother.Theta());
    const ROOT::Math::Boost boostPRF = ROOT::Math::Boost(-betax, -betay, -betaz);
    daughterCMS = boostPRF(daughter);
    threeVecDauCM = daughterCMS.Vect();
    float cosThetaStar = TMath::Abs(threeVecDauCM.Dot(threeVecMother) / std::sqrt(threeVecMother.Mag2()) / std::sqrt(threeVecDauCM.Mag2()));
    return cosThetaStar;
  }

  bool selectionPID(float nsigmaTPC, float nsigmaTOF, int TOFHit, int PIDStrategy, float ptcand)
  {

    if (PIDStrategy == 100) {
      if (ptcand < 1.2) {
        if (TOFHit == 1 && std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.5) {
          return true;
        } else if (TOFHit != 1) {
          if (ptcand < 0.5 && nsigmaTPC > -2.0 && nsigmaTPC < 2.5) {
            return true;
          }
          if (ptcand >= 0.5 && ptcand < 0.6 && nsigmaTPC > -1.5 && nsigmaTPC < 2.5) {
            return true;
          }
          if (ptcand >= 0.6 && ptcand < 0.7 && nsigmaTPC > -1.0 && nsigmaTPC < 2.0) {
            return true;
          }
          if (ptcand >= 0.7 && ptcand < 0.8 && nsigmaTPC > -0.4 && nsigmaTPC < 2.0) {
            return true;
          }
          if (ptcand >= 0.8 && ptcand < 1.0 && nsigmaTPC > 0.0 && nsigmaTPC < 2.0) {
            return true;
          }
          if (ptcand >= 1.0 && ptcand < 1.2 && nsigmaTPC > -2.5 && nsigmaTPC < 0.5) {
            return true;
          }
        }
      } else {
        if (TOFHit == 1 && std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.0) {
          return true;
        } else if (TOFHit != 1 && nsigmaTPC > -2.5 && nsigmaTPC < 1.0) {
          return true;
        }
      }
    }

    if (PIDStrategy == 101) {
      if (ptcand < 1.0) {
        if (TOFHit == 1 && std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.5) {
          return true;
        } else if (TOFHit != 1) {
          if (ptcand < 0.5 && nsigmaTPC > -2.0 && nsigmaTPC < 2.5) {
            return true;
          }
          if (ptcand >= 0.5 && ptcand < 0.6 && nsigmaTPC > -1.5 && nsigmaTPC < 2.5) {
            return true;
          }
          if (ptcand >= 0.6 && ptcand < 0.7 && nsigmaTPC > -1.0 && nsigmaTPC < 2.0) {
            return true;
          }
          if (ptcand >= 0.7 && ptcand < 0.8 && nsigmaTPC > -0.4 && nsigmaTPC < 2.0) {
            return true;
          }
          if (ptcand >= 0.8 && ptcand < 1.0 && nsigmaTPC > 0.0 && nsigmaTPC < 2.0) {
            return true;
          }
        }
      } else if (ptcand >= 1.0 && ptcand < 2.0 && TOFHit == 1 && std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.5) {
        return true;
      } else if (ptcand > 2.0) {
        if (TOFHit == 1 && std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.0) {
          return true;
        } else if (TOFHit != 1 && nsigmaTPC > -2.5 && nsigmaTPC < 1.0) {
          return true;
        }
      }
    }

    if (PIDStrategy == 102) {
      if (TOFHit != 1) {
        if (ptcand < 0.5 && nsigmaTPC > -2.0 && nsigmaTPC < 2.5) {
          return true;
        }
        if (ptcand >= 0.5 && ptcand < 0.6 && nsigmaTPC > -1.5 && nsigmaTPC < 2.5) {
          return true;
        }
        if (ptcand >= 0.6 && ptcand < 0.7 && nsigmaTPC > -1.0 && nsigmaTPC < 2.0) {
          return true;
        }
        if (ptcand >= 2.2 && nsigmaTPC > -2.5 && nsigmaTPC < 1.0) {
          return true;
        }
      }
      if (TOFHit == 1 && ptcand > 0.4 && std::sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.0) {
        return true;
      }
    }

    // optimized TPC TOF
    if (PIDStrategy == 0) {
      if (ptcand < 0.4) {
        if (nsigmaTPC > -3.0 && nsigmaTPC < 3.0) {
          return true;
        }
      } else if (ptcand >= 0.4 && ptcand < 0.5) {
        if (nsigmaTPC > -2.0 && nsigmaTPC < 3.0) {
          return true;
        }
      } else if (ptcand >= 0.5 && ptcand < 5.0 && TOFHit == 1) {
        if (ptcand < 2.0 && TMath::Sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.5) {
          return true;
        }
        if (ptcand >= 2.0 && TMath::Sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.0) {
          return true;
        }
      } else if (ptcand >= 0.5 && ptcand < 5.0 && TOFHit != 1) {
        if (ptcand >= 0.5 && ptcand < 0.6 && nsigmaTPC > -1.5 && nsigmaTPC < 2.0) {
          return true;
        }
        if (ptcand >= 0.6 && ptcand < 0.7 && nsigmaTPC > -1.0 && nsigmaTPC < 2.0) {
          return true;
        }
        if (ptcand >= 0.7 && ptcand < 0.8 && nsigmaTPC > -0.4 && nsigmaTPC < 2.0) {
          return true;
        }
        if (ptcand >= 0.8 && ptcand < 1.0 && nsigmaTPC > -0.0 && nsigmaTPC < 2.0) {
          return true;
        }
        if (ptcand >= 1.0 && ptcand < 1.8 && nsigmaTPC > -2.0 && nsigmaTPC < 2.0) {
          return true;
        }
        if (ptcand >= 1.8 && ptcand < 2.0 && nsigmaTPC > -2.0 && nsigmaTPC < 1.5) {
          return true;
        }
        if (ptcand >= 2.0 && nsigmaTPC > -2.0 && nsigmaTPC < 1.0) {
          return true;
        }
      } else if (ptcand >= 5.0 && nsigmaTPC > -2.0 && nsigmaTPC < 2.0) {
        return true;
      }
    }
    // optimized TPC TOF combined
    if (PIDStrategy == 1) {
      if (ptcand < 0.4) {
        if (nsigmaTPC > cutMinNsigmaTPC && nsigmaTPC < cutNsigmaTPC) {
          return true;
        }
      } else if (ptcand >= 0.4 && ptcand < 0.5) {
        if (nsigmaTPC > -2.0 && nsigmaTPC < cutNsigmaTPC) {
          return true;
        }
      } else if (ptcand >= 0.5 && ptcand < 5.0 && TOFHit == 1) {
        if (ptcand < 2.0 && TMath::Sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.5) {
          return true;
        }
        if (ptcand >= 2.0 && TMath::Sqrt(nsigmaTOF * nsigmaTOF + nsigmaTPC * nsigmaTPC) < 2.0) {
          return true;
        }
      } else if (ptcand >= 5.0 && nsigmaTPC > -2.0 && nsigmaTPC < 2.0) {
        return true;
      }
    }

    if (PIDStrategy == 2) {
      if (ptcand < 0.5) {
        if (nsigmaTPC > cutMinNsigmaTPC && nsigmaTPC < cutNsigmaTPC) {
          return true;
        }
      }
      if (ptcand >= 0.5) {
        if (TOFHit != 1 && ptcand < momTOFCut) {
          if (ptcand >= 0.5 && ptcand < 0.6 && nsigmaTPC > -1.5 && nsigmaTPC < cutNsigmaTPC) {
            return true;
          }
          if (ptcand >= 0.6 && ptcand < 0.7 && nsigmaTPC > -1.0 && nsigmaTPC < cutNsigmaTPC) {
            return true;
          }
          if (ptcand >= 0.7 && ptcand < 0.8 && nsigmaTPC > -0.4 && nsigmaTPC < cutNsigmaTPC) {
            return true;
          }
          if (ptcand >= 0.8 && ptcand < 1.0 && nsigmaTPC > -0.0 && nsigmaTPC < cutNsigmaTPC) {
            return true;
          }
          if (ptcand >= 1.0 && ptcand < 1.8 && nsigmaTPC > -2.0 && nsigmaTPC < 2.0) {
            return true;
          }
          if (ptcand >= 1.8 && ptcand < 2.0 && nsigmaTPC > -2.0 && nsigmaTPC < 1.5) {
            return true;
          }
          if (ptcand >= 2.0 && nsigmaTPC > -2.0 && nsigmaTPC < 1.0) {
            return true;
          }
        }
        if (TOFHit == 1) {
          if (TMath::Sqrt((nsigmaTPC * nsigmaTPC + nsigmaTOF * nsigmaTOF) / 2.0) < cutNsigmaTOF) {
            return true;
          }
        }
      }
    }
    if (PIDStrategy == 3) {
      if (ptcand < 0.5) {
        if (nsigmaTPC > cutMinNsigmaTPC && nsigmaTPC < cutNsigmaTPC) {
          return true;
        }
      }
      if (ptcand >= 0.5) {
        if (TOFHit != 1) {
          if (nsigmaTPC > cutMinNsigmaTPC && nsigmaTPC < cutNsigmaTPC) {
            return true;
          }
        }
        if (TOFHit == 1) {
          if (TMath::Sqrt((nsigmaTPC * nsigmaTPC + nsigmaTOF * nsigmaTOF) / 2.0) < cutNsigmaTOF) {
            return true;
          }
        }
      }
    }
    return false;
  }

  TLorentzVector exotic, Phid1, Phid2;
  TLorentzVector Phi1kaonplus, Phi1kaonminus, Phi2kaonplus, Phi2kaonminus;
  // TLorentzVector exoticRot, Phid1Rot;

  void processSE(aod::RedPhiEvents::iterator const& collision, aod::PhiTracks const& phitracks)
  {
    if (additionalEvsel && (collision.numPos() < 2 || collision.numNeg() < 2)) {
      return;
    }
    int phimult = 0;
    for (auto phitrackd1 : phitracks) {
      if (phitrackd1.phiMass() < minPhiMass1 || phitrackd1.phiMass() > maxPhiMass1) {
        continue;
      }
      auto kaonplusd1pt = TMath::Sqrt(phitrackd1.phid1Px() * phitrackd1.phid1Px() + phitrackd1.phid1Py() * phitrackd1.phid1Py());
      auto kaonminusd1pt = TMath::Sqrt(phitrackd1.phid2Px() * phitrackd1.phid2Px() + phitrackd1.phid2Py() * phitrackd1.phid2Py());
      if (kaonplusd1pt > maxKaonPt) {
        continue;
      }
      if (kaonminusd1pt > maxKaonPt) {
        continue;
      }
      if (!selectionPID(phitrackd1.phid1TPC(), phitrackd1.phid1TOF(), phitrackd1.phid1TOFHit(), strategyPID1, kaonplusd1pt)) {
        continue;
      }
      if (!selectionPID(phitrackd1.phid2TPC(), phitrackd1.phid2TOF(), phitrackd1.phid2TOFHit(), strategyPID2, kaonminusd1pt)) {
        continue;
      }
      phimult = phimult + 1;
    }
    for (auto phitrackd1 : phitracks) {
      auto kaonplusd1pt = TMath::Sqrt(phitrackd1.phid1Px() * phitrackd1.phid1Px() + phitrackd1.phid1Py() * phitrackd1.phid1Py());
      auto kaonminusd1pt = TMath::Sqrt(phitrackd1.phid2Px() * phitrackd1.phid2Px() + phitrackd1.phid2Py() * phitrackd1.phid2Py());
      if (kaonplusd1pt > maxKaonPt) {
        continue;
      }
      if (kaonminusd1pt > maxKaonPt) {
        continue;
      }
      if (!selectionPID(phitrackd1.phid1TPC(), phitrackd1.phid1TOF(), phitrackd1.phid1TOFHit(), strategyPID1, kaonplusd1pt)) {
        continue;
      }
      if (!selectionPID(phitrackd1.phid2TPC(), phitrackd1.phid2TOF(), phitrackd1.phid2TOFHit(), strategyPID1, kaonminusd1pt)) {
        continue;
      }
      histos.fill(HIST("hnsigmaTPCTOFKaon"), phitrackd1.phid1TPC(), phitrackd1.phid1TOF(), kaonplusd1pt);
      histos.fill(HIST("hnsigmaTPCKaonPlus"), phitrackd1.phid1TPC(), kaonplusd1pt);
      histos.fill(HIST("hnsigmaTPCKaonMinus"), phitrackd1.phid2TPC(), kaonminusd1pt);
      histos.fill(HIST("hPhiMass"), Phid1.M(), Phid1.Pt());
      auto phid1id = phitrackd1.index();
      Phid1.SetXYZM(phitrackd1.phiPx(), phitrackd1.phiPy(), phitrackd1.phiPz(), phitrackd1.phiMass());
      Phi1kaonplus.SetXYZM(phitrackd1.phid1Px(), phitrackd1.phid1Py(), phitrackd1.phid1Pz(), 0.493);
      Phi1kaonminus.SetXYZM(phitrackd1.phid2Px(), phitrackd1.phid2Py(), phitrackd1.phid2Pz(), 0.493);
      for (auto phitrackd2 : phitracks) {
        auto phid2id = phitrackd2.index();
        if (phid2id <= phid1id) {
          continue;
        }
        auto kaonplusd2pt = TMath::Sqrt(phitrackd2.phid1Px() * phitrackd2.phid1Px() + phitrackd2.phid1Py() * phitrackd2.phid1Py());
        auto kaonminusd2pt = TMath::Sqrt(phitrackd2.phid2Px() * phitrackd2.phid2Px() + phitrackd2.phid2Py() * phitrackd2.phid2Py());
        if (kaonplusd2pt > maxKaonPt) {
          continue;
        }
        if (kaonminusd2pt > maxKaonPt) {
          continue;
        }
        if (!selectionPID(phitrackd2.phid1TPC(), phitrackd2.phid1TOF(), phitrackd2.phid1TOFHit(), strategyPID2, kaonplusd2pt)) {
          continue;
        }
        if (!selectionPID(phitrackd2.phid2TPC(), phitrackd2.phid2TOF(), phitrackd2.phid2TOFHit(), strategyPID2, kaonminusd2pt)) {
          continue;
        }
        if (phitrackd1.phid1Index() == phitrackd2.phid1Index()) {
          continue;
        }
        if (phitrackd1.phid2Index() == phitrackd2.phid2Index()) {
          continue;
        }
        Phid2.SetXYZM(phitrackd2.phiPx(), phitrackd2.phiPy(), phitrackd2.phiPz(), phitrackd2.phiMass());
        Phi2kaonplus.SetXYZM(phitrackd2.phid1Px(), phitrackd2.phid1Py(), phitrackd2.phid1Pz(), 0.493);
        Phi2kaonminus.SetXYZM(phitrackd2.phid2Px(), phitrackd2.phid2Py(), phitrackd2.phid2Pz(), 0.493);

        /*
              // Like
              Phid1like = Phi1kaonplus + Phi2kaonplus;
              Phid2like = Phi1kaonminus + Phi2kaonminus;
              exoticlike = Phid1like + Phid2like;
              auto deltaRlike = TMath::Sqrt(TMath::Power(Phid1like.Phi() - Phid2like.Phi(), 2.0) + TMath::Power(Phid1like.Eta() - Phid2like.Eta(), 2.0));
              auto costhetalike = (Phid1like.Px() * Phid2like.Px() + Phid1like.Py() * Phid2like.Py() + Phid1like.Pz() * Phid2like.Pz()) / (Phid1like.P() * Phid2like.P());
              auto deltamlike = TMath::Sqrt(TMath::Power(Phid1like.M() - 1.0192, 2.0) + TMath::Power(Phid2like.M() - 1.0192, 2.0));
              if (!isDeep) {
                histos.fill(HIST("SEMassLike"), exoticlike.M(), exoticlike.Pt(), deltaRlike, costhetalike, deltamlike, phimult);
              }
              if (isDeep) {
                histos.fill(HIST("SEMassLike"), exoticlike.M(), exoticlike.Pt(), deltaRlike, deepangle(Phid1like, Phid2like), deltamlike, phimult);
              }
        */

        // Unlike
        histos.fill(HIST("hPhiMass2"), Phid1.M(), Phid2.M());
        if (phitrackd2.phiMass() < minPhiMass2 || phitrackd2.phiMass() > maxPhiMass2) {
          continue;
        }
        if (phitrackd1.phiMass() < minPhiMass1 || phitrackd1.phiMass() > maxPhiMass1) {
          continue;
        }
        exotic = Phid1 + Phid2;
        if (exotic.M() < minExoticMass || exotic.M() > maxExoticMass) {
          continue;
        }
        histos.fill(HIST("hkPlusDeltaetaDeltaPhi"), Phi1kaonplus.Eta() - Phi2kaonplus.Eta(), Phi1kaonplus.Phi() - Phi2kaonplus.Phi());
        histos.fill(HIST("hkMinusDeltaetaDeltaPhi"), Phi1kaonminus.Eta() - Phi2kaonminus.Eta(), Phi1kaonminus.Phi() - Phi2kaonminus.Phi());
        // auto cosThetaStar = getCosTheta(exotic, Phid1);
        // auto kstar = getkstar(Phid1, Phid2);
        auto deltaR = TMath::Sqrt(TMath::Power(Phid1.Phi() - Phid2.Phi(), 2.0) + TMath::Power(Phid1.Eta() - Phid2.Eta(), 2.0));
        auto deltaRd1 = TMath::Sqrt(TMath::Power(Phi1kaonplus.Phi() - Phi2kaonplus.Phi(), 2.0) + TMath::Power(Phi1kaonplus.Eta() - Phi2kaonplus.Eta(), 2.0));
        auto deltaRd2 = TMath::Sqrt(TMath::Power(Phi1kaonminus.Phi() - Phi2kaonminus.Phi(), 2.0) + TMath::Power(Phi1kaonminus.Eta() - Phi2kaonminus.Eta(), 2.0));
        auto deltam = TMath::Sqrt(TMath::Power(Phid1.M() - 1.0192, 2.0) + TMath::Power(Phid2.M() - 1.0192, 2.0));
        if (deltaRd1 < daughterDeltaR) {
          continue;
        }
        if (deltaRd2 < daughterDeltaR) {
          continue;
        }
        if (!isDeep) {
          histos.fill(HIST("SEMassUnlike"), exotic.M(), std::abs(Phid1.Pt() - Phid2.Pt()) / exotic.Pt(), exotic.Pt(), deltaR, deltam, phimult);
        }
        if (isDeep) {
          histos.fill(HIST("SEMassUnlike"), exotic.M(), std::abs(Phid1.Pt() - Phid2.Pt()) / exotic.Pt(), exotic.Pt(), deltaR, deltam, phimult);
        }
      }
    }
  }
  PROCESS_SWITCH(doublephimeson, processSE, "Process Same Event", false);
  void processopti(aod::RedPhiEvents::iterator const& collision, aod::PhiTracks const& phitracks)
  {
    std::vector<ROOT::Math::PtEtaPhiMVector> exoticresonance, phiresonanced1, phiresonanced2, kaonplus1, kaonplus2, kaonminus1, kaonminus2;
    std::vector<int> d1trackid = {};
    std::vector<int> d2trackid = {};
    std::vector<int> d3trackid = {};
    std::vector<int> d4trackid = {};
    if (additionalEvsel && (collision.numPos() < 2 || collision.numNeg() < 2)) {
      return;
    }
    int phimult = 0;

    for (auto phitrackd1 : phitracks) {
      if (phitrackd1.phiMass() < minPhiMass1 || phitrackd1.phiMass() > maxPhiMass1) {
        continue;
      }
      auto kaonplusd1pt = TMath::Sqrt(phitrackd1.phid1Px() * phitrackd1.phid1Px() + phitrackd1.phid1Py() * phitrackd1.phid1Py());
      auto kaonminusd1pt = TMath::Sqrt(phitrackd1.phid2Px() * phitrackd1.phid2Px() + phitrackd1.phid2Py() * phitrackd1.phid2Py());
      if (kaonplusd1pt > maxKaonPt) {
        continue;
      }
      if (kaonminusd1pt > maxKaonPt) {
        continue;
      }
      if (!selectionPID(phitrackd1.phid1TPC(), phitrackd1.phid1TOF(), phitrackd1.phid1TOFHit(), strategyPID1, kaonplusd1pt)) {
        continue;
      }
      if (!selectionPID(phitrackd1.phid2TPC(), phitrackd1.phid2TOF(), phitrackd1.phid2TOFHit(), strategyPID1, kaonminusd1pt)) {
        continue;
      }
      phimult = phimult + 1;
    }
    for (auto phitrackd1 : phitracks) {
      auto kaonplusd1pt = TMath::Sqrt(phitrackd1.phid1Px() * phitrackd1.phid1Px() + phitrackd1.phid1Py() * phitrackd1.phid1Py());
      auto kaonminusd1pt = TMath::Sqrt(phitrackd1.phid2Px() * phitrackd1.phid2Px() + phitrackd1.phid2Py() * phitrackd1.phid2Py());

      histos.fill(HIST("hnsigmaTPCTOFKaonBefore"), phitrackd1.phid1TPC(), phitrackd1.phid1TOF(), kaonplusd1pt);
      histos.fill(HIST("hnsigmaTPCKaonPlusBefore"), phitrackd1.phid1TPC(), kaonplusd1pt);
      histos.fill(HIST("hnsigmaTPCKaonMinusBefore"), phitrackd1.phid2TPC(), kaonminusd1pt);

      if (kaonplusd1pt > maxKaonPt) {
        continue;
      }
      if (kaonminusd1pt > maxKaonPt) {
        continue;
      }
      if (!selectionPID(phitrackd1.phid1TPC(), phitrackd1.phid1TOF(), phitrackd1.phid1TOFHit(), strategyPID1, kaonplusd1pt)) {
        continue;
      }
      if (!selectionPID(phitrackd1.phid2TPC(), phitrackd1.phid2TOF(), phitrackd1.phid2TOFHit(), strategyPID1, kaonminusd1pt)) {
        continue;
      }
      histos.fill(HIST("hnsigmaTPCTOFKaon"), phitrackd1.phid1TPC(), phitrackd1.phid1TOF(), kaonplusd1pt);
      histos.fill(HIST("hnsigmaTPCKaonPlus"), phitrackd1.phid1TPC(), kaonplusd1pt);
      histos.fill(HIST("hnsigmaTPCKaonMinus"), phitrackd1.phid2TPC(), kaonminusd1pt);
      histos.fill(HIST("hPhiMass"), Phid1.M(), Phid1.Pt());
      auto phid1id = phitrackd1.index();
      Phid1.SetXYZM(phitrackd1.phiPx(), phitrackd1.phiPy(), phitrackd1.phiPz(), phitrackd1.phiMass());
      Phi1kaonplus.SetXYZM(phitrackd1.phid1Px(), phitrackd1.phid1Py(), phitrackd1.phid1Pz(), 0.493);
      Phi1kaonminus.SetXYZM(phitrackd1.phid2Px(), phitrackd1.phid2Py(), phitrackd1.phid2Pz(), 0.493);
      for (auto phitrackd2 : phitracks) {
        auto phid2id = phitrackd2.index();
        if (phid2id <= phid1id) {
          continue;
        }
        auto kaonplusd2pt = TMath::Sqrt(phitrackd2.phid1Px() * phitrackd2.phid1Px() + phitrackd2.phid1Py() * phitrackd2.phid1Py());
        auto kaonminusd2pt = TMath::Sqrt(phitrackd2.phid2Px() * phitrackd2.phid2Px() + phitrackd2.phid2Py() * phitrackd2.phid2Py());
        if (kaonplusd2pt > maxKaonPt) {
          continue;
        }
        if (kaonminusd2pt > maxKaonPt) {
          continue;
        }
        if (!selectionPID(phitrackd2.phid1TPC(), phitrackd2.phid1TOF(), phitrackd2.phid1TOFHit(), strategyPID2, kaonplusd2pt)) {
          continue;
        }
        if (!selectionPID(phitrackd2.phid2TPC(), phitrackd2.phid2TOF(), phitrackd2.phid2TOFHit(), strategyPID2, kaonminusd2pt)) {
          continue;
        }
        if ((phitrackd1.phid1Index() == phitrackd2.phid1Index()) || (phitrackd1.phid2Index() == phitrackd2.phid2Index())) {
          continue;
        }
        Phid2.SetXYZM(phitrackd2.phiPx(), phitrackd2.phiPy(), phitrackd2.phiPz(), phitrackd2.phiMass());
        Phi2kaonplus.SetXYZM(phitrackd2.phid1Px(), phitrackd2.phid1Py(), phitrackd2.phid1Pz(), 0.493);
        Phi2kaonminus.SetXYZM(phitrackd2.phid2Px(), phitrackd2.phid2Py(), phitrackd2.phid2Pz(), 0.493);

        // unlike
        if (phitrackd1.phiMass() < minPhiMass1 || phitrackd1.phiMass() > maxPhiMass1) {
          continue;
        }
        if (phitrackd2.phiMass() < minPhiMass2 || phitrackd2.phiMass() > maxPhiMass2) {
          continue;
        }
        exotic = Phid1 + Phid2;
        if (exotic.M() < minExoticMass || exotic.M() > maxExoticMass) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector temp1(exotic.Pt(), exotic.Eta(), exotic.Phi(), exotic.M());
        ROOT::Math::PtEtaPhiMVector temp2(Phid1.Pt(), Phid1.Eta(), Phid1.Phi(), Phid1.M());
        ROOT::Math::PtEtaPhiMVector temp3(Phid2.Pt(), Phid2.Eta(), Phid2.Phi(), Phid2.M());
        exoticresonance.push_back(temp1);
        phiresonanced1.push_back(temp2);
        phiresonanced2.push_back(temp3);
        d1trackid.push_back(phitrackd1.phid1Index());
        d2trackid.push_back(phitrackd2.phid1Index());
        d3trackid.push_back(phitrackd1.phid2Index());
        d4trackid.push_back(phitrackd2.phid2Index());

        ROOT::Math::PtEtaPhiMVector temp4(Phi1kaonplus.Pt(), Phi1kaonplus.Eta(), Phi1kaonplus.Phi(), 0.493);
        ROOT::Math::PtEtaPhiMVector temp5(Phi1kaonminus.Pt(), Phi1kaonminus.Eta(), Phi1kaonminus.Phi(), 0.493);
        ROOT::Math::PtEtaPhiMVector temp6(Phi2kaonplus.Pt(), Phi2kaonplus.Eta(), Phi2kaonplus.Phi(), 0.493);
        ROOT::Math::PtEtaPhiMVector temp7(Phi2kaonminus.Pt(), Phi2kaonminus.Eta(), Phi2kaonminus.Phi(), 0.493);
        kaonplus1.push_back(temp4);
        kaonplus2.push_back(temp6);
        kaonminus1.push_back(temp5);
        kaonminus2.push_back(temp7);
      }
    }
    if (exoticresonance.size() == 0) {
      return;
    }
    // LOGF(info, "Total number of exotic: %d", exoticresonance.size());
    if (exoticresonance.size() == 2) {
      for (auto if1 = exoticresonance.begin(); if1 != exoticresonance.end(); ++if1) {
        auto i5 = std::distance(exoticresonance.begin(), if1);

        auto exotic1phi1 = phiresonanced1.at(i5);
        auto exotic1phi2 = phiresonanced2.at(i5);
        auto exotic1 = exoticresonance.at(i5);

        auto exotic1kaonplus1 = kaonplus1.at(i5);
        auto exotic1kaonminus1 = kaonminus1.at(i5);
        auto exotic1kaonplus2 = kaonplus2.at(i5);
        auto exotic1kaonminus2 = kaonminus2.at(i5);
        auto deltaRkaonplus1 = TMath::Sqrt(TMath::Power(exotic1kaonplus1.Phi() - exotic1kaonplus2.Phi(), 2.0) + TMath::Power(exotic1kaonplus1.Eta() - exotic1kaonplus2.Eta(), 2.0));
        auto deltaRkaonminus1 = TMath::Sqrt(TMath::Power(exotic1kaonminus1.Phi() - exotic1kaonminus2.Phi(), 2.0) + TMath::Power(exotic1kaonminus1.Eta() - exotic1kaonminus2.Eta(), 2.0));
        histos.fill(HIST("hDeltaRkaonplus"), deltaRkaonplus1);
        histos.fill(HIST("hDeltaRkaonminus"), deltaRkaonminus1);

        auto deltam1 = TMath::Sqrt(TMath::Power(exotic1phi1.M() - 1.0192, 2.0) + TMath::Power(exotic1phi2.M() - 1.0192, 2.0));
        auto deltaR1 = TMath::Sqrt(TMath::Power(exotic1phi1.Phi() - exotic1phi2.Phi(), 2.0) + TMath::Power(exotic1phi1.Eta() - exotic1phi2.Eta(), 2.0));

        if (deltaRkaonplus1 < daughterDeltaR) {
          continue;
        }
        if (deltaRkaonminus1 < daughterDeltaR) {
          continue;
        }

        for (auto if2 = if1 + 1; if2 != exoticresonance.end(); ++if2) {
          auto i6 = std::distance(exoticresonance.begin(), if2);
          auto exotic2phi1 = phiresonanced1.at(i6);
          auto exotic2phi2 = phiresonanced2.at(i6);
          auto exotic2 = exoticresonance.at(i6);

          auto exotic2kaonplus1 = kaonplus1.at(i6);
          auto exotic2kaonminus1 = kaonminus1.at(i6);
          auto exotic2kaonplus2 = kaonplus2.at(i6);
          auto exotic2kaonminus2 = kaonminus2.at(i6);
          auto deltaRkaonplus2 = TMath::Sqrt(TMath::Power(exotic2kaonplus1.Phi() - exotic2kaonplus2.Phi(), 2.0) + TMath::Power(exotic2kaonplus1.Eta() - exotic2kaonplus2.Eta(), 2.0));
          auto deltaRkaonminus2 = TMath::Sqrt(TMath::Power(exotic2kaonminus1.Phi() - exotic2kaonminus2.Phi(), 2.0) + TMath::Power(exotic2kaonminus1.Eta() - exotic2kaonminus2.Eta(), 2.0));

          auto deltam2 = TMath::Sqrt(TMath::Power(exotic2phi1.M() - 1.0192, 2.0) + TMath::Power(exotic2phi2.M() - 1.0192, 2.0));
          auto deltaR2 = TMath::Sqrt(TMath::Power(exotic2phi1.Phi() - exotic2phi2.Phi(), 2.0) + TMath::Power(exotic2phi1.Eta() - exotic2phi2.Eta(), 2.0));

          if ((d1trackid.at(i5) == d1trackid.at(i6) || d1trackid.at(i5) == d2trackid.at(i6)) &&
              (d2trackid.at(i5) == d1trackid.at(i6) || d2trackid.at(i5) == d2trackid.at(i6)) &&
              (d3trackid.at(i5) == d3trackid.at(i6) || d3trackid.at(i5) == d4trackid.at(i6)) &&
              (d4trackid.at(i5) == d3trackid.at(i6) || d4trackid.at(i5) == d4trackid.at(i6))) {

            if (deltam2 < deltam1 && deltaRkaonplus2 > daughterDeltaR && deltaRkaonminus2 > daughterDeltaR) {
              histos.fill(HIST("SEMassUnlike"), exotic2.M(), std::abs(exotic2phi1.Pt() - exotic2phi2.Pt()) / exotic2.Pt(), exotic2.Pt(), deltaR2, deltam2, phimult);
              // LOGF(info, "Fill exotic Id %d which is pair of Id %d", i6, i5);
            } else {
              histos.fill(HIST("SEMassUnlike"), exotic1.M(), std::abs(exotic2phi1.Pt() - exotic2phi2.Pt()) / exotic1.Pt(), exotic1.Pt(), deltaR1, deltam1, phimult);
            }
          } else {
            histos.fill(HIST("SEMassUnlike"), exotic1.M(), std::abs(exotic2phi1.Pt() - exotic2phi2.Pt()) / exotic1.Pt(), exotic1.Pt(), deltaR1, deltam1, phimult);
          }
        }
      }
    } else {
      for (auto if1 = exoticresonance.begin(); if1 != exoticresonance.end(); ++if1) {
        auto i5 = std::distance(exoticresonance.begin(), if1);
        auto exotic1phi1 = phiresonanced1.at(i5);
        auto exotic1phi2 = phiresonanced2.at(i5);
        auto exotic1 = exoticresonance.at(i5);

        auto exotic1kaonplus1 = kaonplus1.at(i5);
        auto exotic1kaonminus1 = kaonminus1.at(i5);
        auto exotic1kaonplus2 = kaonplus2.at(i5);
        auto exotic1kaonminus2 = kaonminus2.at(i5);
        auto deltaRkaonplus1 = TMath::Sqrt(TMath::Power(exotic1kaonplus1.Phi() - exotic1kaonplus2.Phi(), 2.0) + TMath::Power(exotic1kaonplus1.Eta() - exotic1kaonplus2.Eta(), 2.0));
        auto deltaRkaonminus1 = TMath::Sqrt(TMath::Power(exotic1kaonminus1.Phi() - exotic1kaonminus2.Phi(), 2.0) + TMath::Power(exotic1kaonminus1.Eta() - exotic1kaonminus2.Eta(), 2.0));
        auto deltam1 = TMath::Sqrt(TMath::Power(exotic1phi1.M() - 1.0192, 2.0) + TMath::Power(exotic1phi2.M() - 1.0192, 2.0));
        auto deltaR1 = TMath::Sqrt(TMath::Power(exotic1phi1.Phi() - exotic1phi2.Phi(), 2.0) + TMath::Power(exotic1phi1.Eta() - exotic1phi2.Eta(), 2.0));

        histos.fill(HIST("hDeltaRkaonplus"), deltaRkaonplus1);
        histos.fill(HIST("hDeltaRkaonminus"), deltaRkaonminus1);

        if (deltaRkaonplus1 < daughterDeltaR) {
          continue;
        }
        if (deltaRkaonminus1 < daughterDeltaR) {
          continue;
        }

        histos.fill(HIST("SEMassUnlike"), exotic1.M(), std::abs(exotic1phi1.Pt() - exotic1phi2.Pt()) / exotic1.Pt(), exotic1.Pt(), deltaR1, deltam1, phimult);
      }
    }
  }
  PROCESS_SWITCH(doublephimeson, processopti, "Process Optimized same event", false);

  void processopti3(aod::RedPhiEvents::iterator const& collision, aod::PhiTracks const& phitracks)
  {
    if (additionalEvsel && (collision.numPos() < 2 || collision.numNeg() < 2))
      return;

    // --- φ multiplicity with PID ---
    int phimult = 0;
    for (auto const& t : phitracks) {
      if (t.phiMass() < minPhiMass1 || t.phiMass() > maxPhiMass1)
        continue;
      const double kpluspt = std::hypot(t.phid1Px(), t.phid1Py());
      const double kminuspt = std::hypot(t.phid2Px(), t.phid2Py());
      if (kpluspt > maxKaonPt || kminuspt > maxKaonPt)
        continue;
      if (!selectionPID(t.phid1TPC(), t.phid1TOF(), t.phid1TOFHit(), strategyPID1, kpluspt))
        continue;
      if (!selectionPID(t.phid2TPC(), t.phid2TOF(), t.phid2TOFHit(), strategyPID2, kminuspt))
        continue;
      ++phimult;
    }
    if (phimult < 2)
      return;

    // --- helpers ---
    constexpr double mPhiPDG = 1.019461; // GeV/c^2

    const auto deltaMPhi = [=](double m1, double m2) {
      const double d1 = m1 - mPhiPDG, d2 = m2 - mPhiPDG;
      return std::sqrt(d1 * d1 + d2 * d2);
    };

    const auto deltaR = [](double phi1, double eta1, double phi2, double eta2) {
      const double dphi = TVector2::Phi_mpi_pi(phi1 - phi2);
      const double deta = eta1 - eta2;
      return std::sqrt(dphi * dphi + deta * deta);
    };

    // minimum ΔR among all kaons in the candidate (4 kaons → 6 combinations)
    const auto minKaonDeltaR = [&](const ROOT::Math::PtEtaPhiMVector& kplusA,
                                   const ROOT::Math::PtEtaPhiMVector& kplusB,
                                   const ROOT::Math::PtEtaPhiMVector& kminusA,
                                   const ROOT::Math::PtEtaPhiMVector& kminusB) {
      // same-sign first (keep your QA histos)
      const double dRkplus = deltaR(kplusA.Phi(), kplusA.Eta(), kplusB.Phi(), kplusB.Eta());
      const double dRkminus = deltaR(kminusA.Phi(), kminusA.Eta(), kminusB.Phi(), kminusB.Eta());
      histos.fill(HIST("hDeltaRkaonplus"), dRkplus);
      histos.fill(HIST("hDeltaRkaonminus"), dRkminus);

      // all other combinations
      const double dR_k1p_k1m = deltaR(kplusA.Phi(), kplusA.Eta(), kminusA.Phi(), kminusA.Eta());
      const double dR_k1p_k2m = deltaR(kplusA.Phi(), kplusA.Eta(), kminusB.Phi(), kminusB.Eta());
      const double dR_k2p_k1m = deltaR(kplusB.Phi(), kplusB.Eta(), kminusA.Phi(), kminusA.Eta());
      const double dR_k2p_k2m = deltaR(kplusB.Phi(), kplusB.Eta(), kminusB.Phi(), kminusB.Eta());

      double minDR = dRkplus;
      minDR = std::min(minDR, dRkminus);
      minDR = std::min(minDR, dR_k1p_k1m);
      minDR = std::min(minDR, dR_k1p_k2m);
      minDR = std::min(minDR, dR_k2p_k1m);
      minDR = std::min(minDR, dR_k2p_k2m);
      return minDR;
    };

    // --- collect candidates once ---
    std::vector<ROOT::Math::PtEtaPhiMVector> pairV, phi1V, phi2V, kplus1V, kplus2V, kminus1V, kminus2V;
    std::vector<double> minDRV; // store minimum ΔR for each pair

    for (auto const& t1 : phitracks) {
      const double kplus1pt = std::hypot(t1.phid1Px(), t1.phid1Py());
      const double kminus1pt = std::hypot(t1.phid2Px(), t1.phid2Py());

      // PID QA before
      histos.fill(HIST("hnsigmaTPCTOFKaonBefore"), t1.phid1TPC(), t1.phid1TOF(), kplus1pt);
      histos.fill(HIST("hnsigmaTPCKaonPlusBefore"), t1.phid1TPC(), kplus1pt);
      histos.fill(HIST("hnsigmaTPCKaonMinusBefore"), t1.phid2TPC(), kminus1pt);

      if (kplus1pt > maxKaonPt || kminus1pt > maxKaonPt)
        continue;
      if (!selectionPID(t1.phid1TPC(), t1.phid1TOF(), t1.phid1TOFHit(), strategyPID1, kplus1pt))
        continue;
      if (!selectionPID(t1.phid2TPC(), t1.phid2TOF(), t1.phid2TOFHit(), strategyPID2, kminus1pt))
        continue;

      TLorentzVector phi1, k1p, k1m;
      phi1.SetXYZM(t1.phiPx(), t1.phiPy(), t1.phiPz(), t1.phiMass());
      k1p.SetXYZM(t1.phid1Px(), t1.phid1Py(), t1.phid1Pz(), 0.493);
      k1m.SetXYZM(t1.phid2Px(), t1.phid2Py(), t1.phid2Pz(), 0.493);

      // φ mass windows
      if (t1.phiMass() < minPhiMass1 || t1.phiMass() > maxPhiMass1)
        continue;
      if (phi1.Pt() < minPhiPt || phi1.Pt() > maxPhiPt)
        continue;
      // PID QA after
      histos.fill(HIST("hnsigmaTPCTOFKaon"), t1.phid1TPC(), t1.phid1TOF(), kplus1pt);
      histos.fill(HIST("hnsigmaTPCKaonPlus"), t1.phid1TPC(), kplus1pt);
      histos.fill(HIST("hnsigmaTPCKaonMinus"), t1.phid2TPC(), kminus1pt);
      histos.fill(HIST("hPhiMass"), phi1.M(), phi1.Pt());

      const auto id1 = t1.index();

      for (auto const& t2 : phitracks) {
        const auto id2 = t2.index();
        if (id2 <= id1)
          continue;

        const double kplus2pt = std::hypot(t2.phid1Px(), t2.phid1Py());
        const double kminus2pt = std::hypot(t2.phid2Px(), t2.phid2Py());
        if (kplus2pt > maxKaonPt || kminus2pt > maxKaonPt)
          continue;
        if (!selectionPID(t2.phid1TPC(), t2.phid1TOF(), t2.phid1TOFHit(), strategyPID1, kplus2pt))
          continue;
        if (!selectionPID(t2.phid2TPC(), t2.phid2TOF(), t2.phid2TOFHit(), strategyPID2, kminus2pt))
          continue;

        // block shared same-sign daughters
        if ((t1.phid1Index() == t2.phid1Index()) || (t1.phid2Index() == t2.phid2Index()))
          continue;

        TLorentzVector phi2, k2p, k2m;
        phi2.SetXYZM(t2.phiPx(), t2.phiPy(), t2.phiPz(), t2.phiMass());
        k2p.SetXYZM(t2.phid1Px(), t2.phid1Py(), t2.phid1Pz(), 0.493);
        k2m.SetXYZM(t2.phid2Px(), t2.phid2Py(), t2.phid2Pz(), 0.493);
        if (t2.phiMass() < minPhiMass2 || t2.phiMass() > maxPhiMass2)
          continue;
        if (phi1.Pt() < minPhiPt || phi1.Pt() > maxPhiPt)
          continue;
        // Δm cut (configurable)
        const double dM = deltaMPhi(phi1.M(), phi2.M());
        if (dM > maxDeltaMPhi)
          continue;

        TLorentzVector pair = phi1 + phi2;
        if (pair.M() < minExoticMass || pair.M() > maxExoticMass)
          continue;

        // daughter ΔR QA and minΔR (NO CUT anymore)
        ROOT::Math::PtEtaPhiMVector k1pV(k1p.Pt(), k1p.Eta(), k1p.Phi(), 0.493);
        ROOT::Math::PtEtaPhiMVector k1mV(k1m.Pt(), k1m.Eta(), k1m.Phi(), 0.493);
        ROOT::Math::PtEtaPhiMVector k2pV(k2p.Pt(), k2p.Eta(), k2p.Phi(), 0.493);
        ROOT::Math::PtEtaPhiMVector k2mV(k2m.Pt(), k2m.Eta(), k2m.Phi(), 0.493);
        const double minDR = minKaonDeltaR(k1pV, k2pV, k1mV, k2mV);

        // store for one-pass fill
        pairV.emplace_back(pair.Pt(), pair.Eta(), pair.Phi(), pair.M());
        phi1V.emplace_back(phi1.Pt(), phi1.Eta(), phi1.Phi(), phi1.M());
        phi2V.emplace_back(phi2.Pt(), phi2.Eta(), phi2.Phi(), phi2.M());
        kplus1V.emplace_back(k1p.Pt(), k1p.Eta(), k1p.Phi(), 0.493);
        kminus1V.emplace_back(k1m.Pt(), k1m.Eta(), k1m.Phi(), 0.493);
        kplus2V.emplace_back(k2p.Pt(), k2p.Eta(), k2p.Phi(), 0.493);
        kminus2V.emplace_back(k2m.Pt(), k2m.Eta(), k2m.Phi(), 0.493);
        minDRV.emplace_back(minDR); // per-candidate minimum ΔR of kaons
      }
    }

    if (pairV.empty())
      return;

    // --- fill the single THnSparse ---
    for (size_t i = 0; i < pairV.size(); ++i) {
      TLorentzVector p1, p2, pair;
      p1.SetPtEtaPhiM(phi1V[i].Pt(), phi1V[i].Eta(), phi1V[i].Phi(), phi1V[i].M());
      p2.SetPtEtaPhiM(phi2V[i].Pt(), phi2V[i].Eta(), phi2V[i].Phi(), phi2V[i].M());
      pair.SetPtEtaPhiM(pairV[i].Pt(), pairV[i].Eta(), pairV[i].Phi(), pairV[i].M());

      const double dM = deltaMPhi(p1.M(), p2.M());
      const double M = pair.M();
      const double dR = deltaR(p1.Phi(), p1.Eta(), p2.Phi(), p2.Eta());
      const double minDR = minDRV[i];
      double ptcorr = p1.Pt() / (pair.Pt() - p1.Pt());
      histos.fill(HIST("hPtCorrelation"), pair.Pt(), ptcorr);
      // NOTE: second axis is now minΔR(all kaons), ΔpT/pT has been removed
      histos.fill(HIST("SEMassUnlike"),
                  M,
                  minDR,
                  pair.Pt(),
                  dR,
                  dM,
                  ptcorr);
    }
  }
  PROCESS_SWITCH(doublephimeson, processopti3, "Process Optimized same event", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<doublephimeson>(cfgc)}; }
