// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file FemtoUniverseDetaDphiStar.h
/// \brief FemtoUniverseDetaDphiStar - Checks particles for the close pair rejection.
/// \author Laura Serksnyte, TU München, laura.serksnyte@tum.de
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch
/// \author Pritam Chakraborty, WUT Warsaw, pritam.chakraborty@cern.ch
/// \author Shirajum Monira, WUT Warsaw, shirajum.monira@cern.ch

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEDETADPHISTAR_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEDETADPHISTAR_H_

#include "PWGCF/FemtoUniverse/Core/FemtoUniverseAngularContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseFemtoContainer.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseTrackSelection.h"
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"

#include "Framework/HistogramRegistry.h"

#include "TMath.h"

#include <memory>
#include <string>
#include <vector>

namespace o2::analysis
{
namespace femto_universe
{

/// \class FemtoUniverseDetaDphiStar
/// \brief Class to check particles for the close pair rejection.
/// \tparam partOne Type of particle 1 (Track/V0/Cascade/...)
/// \tparam partTwo Type of particle 2 (Track/V0/Cascade/...)
template <o2::aod::femtouniverseparticle::ParticleType partOne, o2::aod::femtouniverseparticle::ParticleType partTwo>
class FemtoUniverseDetaDphiStar
{
 public:
  FemtoUniverseTrackSelection trackCuts;
  /// Destructor
  virtual ~FemtoUniverseDetaDphiStar() = default;
  /// Initialization of the histograms and setting required values
  void init(HistogramRegistry* registry, HistogramRegistry* registryQA, float ldeltaphistarcutmin, float ldeltaphistarcutmax, float ldeltaetacutmin, float ldeltaetacutmax, float lchosenradii, bool lplotForEveryRadii, float lPhiMassMin = 1.014, float lPhiMassMax = 1.026, bool lisSameSignCPR = false)
  {
    chosenRadii = lchosenradii;
    cutDeltaPhiStarMax = ldeltaphistarcutmax;
    cutDeltaPhiStarMin = ldeltaphistarcutmin;
    cutDeltaEtaMax = ldeltaetacutmax;
    cutDeltaEtaMin = ldeltaetacutmin;
    plotForEveryRadii = lplotForEveryRadii;
    mHistogramRegistry = registry;
    mHistogramRegistryQA = registryQA;
    cutPhiInvMassLow = lPhiMassMin;
    cutPhiInvMassHigh = lPhiMassMax;
    isSameSignCPR = lisSameSignCPR;

    if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kTrack) {
      std::string dirName = static_cast<std::string>(DirNames[0]);
      histdetadpisame[0][0] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesSame[0][0])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
      histdetadpisame[0][1] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesSame[1][0])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
      histdetadpimixed[0][0] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesMixed[0][0])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
      histdetadpimixed[0][1] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesMixed[1][0])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});

      histdetadpiqlcmssame = mHistogramRegistry->add<TH3>((dirName + static_cast<std::string>(HistNamesSame[1][7])).c_str(), "; #it{q}_{LCMS}; #Delta #eta; #Delta #phi", kTH3F, {{100, 0.0, 0.5}, {100, -0.15, 0.15}, {100, -0.15, 0.15}});
      histdetadpiqlcmsmixed = mHistogramRegistry->add<TH3>((dirName + static_cast<std::string>(HistNamesMixed[1][7])).c_str(), "; #it{q}_{LCMS}; #Delta #eta; #Delta #phi", kTH3F, {{100, 0.0, 0.5}, {100, -0.15, 0.15}, {100, -0.15, 0.15}});

      if (plotForEveryRadii) {
        for (int i = 0; i < 9; i++) {
          histdetadpiRadii[0][i] = mHistogramRegistryQA->add<TH2>((dirName + static_cast<std::string>(HistNamesRadii[0][i])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        }
      }
    }
    if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kV0) {
      for (int i = 0; i < 2; i++) {
        std::string dirName = static_cast<std::string>(DirNames[1]);
        histdetadpisame[i][0] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesSame[0][i])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpisame[i][1] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesSame[1][i])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpimixed[i][0] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesMixed[0][i])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpimixed[i][1] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesMixed[1][i])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});

        if (plotForEveryRadii) {
          for (int j = 0; j < 9; j++) {
            histdetadpiRadii[i][j] = mHistogramRegistryQA->add<TH2>((dirName + static_cast<std::string>(HistNamesRadii[i][j])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
          }
        }
      }
    }
    if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kV0 && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kV0) {
      /// V0-V0 combination
      for (int k = 0; k < 2; k++) {
        std::string dirName = static_cast<std::string>(DirNames[2]);
        histdetadpisame[k][0] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesSame[0][k])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpisame[k][1] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesSame[1][k])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpimixed[k][0] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesMixed[0][k])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpimixed[k][1] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesMixed[1][k])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        if (plotForEveryRadii) {
          for (int l = 0; l < 9; l++) {
            histdetadpiRadii[k][l] = mHistogramRegistryQA->add<TH2>((dirName + static_cast<std::string>(HistNamesRadii[k][l])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
          }
        }
      }
    }
    if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kCascade && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kCascade) {
      /// Cascade-Cascade combination
      for (int k = 0; k < 7; k++) {
        std::string dirName = static_cast<std::string>(DirNames[5]);
        histdetadpisame[k][0] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesSame[0][k])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpisame[k][1] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesSame[1][k])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpimixed[k][0] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesMixed[0][k])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpimixed[k][1] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesMixed[1][k])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        if (plotForEveryRadii) {
          for (int l = 0; l < 9; l++) {
            histdetadpiRadii[k][l] = mHistogramRegistryQA->add<TH2>((dirName + static_cast<std::string>(HistNamesRadii[k][l])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
          }
        }
      }
    }
    if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kCascade) {
      /// Track-Cascade combination
      for (int k = 0; k < 3; k++) {
        std::string dirName = static_cast<std::string>(DirNames[6]);
        histdetadpisame[k][0] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesSame[0][k])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpisame[k][1] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesSame[1][k])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpimixed[k][0] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesMixed[0][k])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpimixed[k][1] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesMixed[1][k])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        if (plotForEveryRadii) {
          for (int l = 0; l < 9; l++) {
            histdetadpiRadii[k][l] = mHistogramRegistryQA->add<TH2>((dirName + static_cast<std::string>(HistNamesRadii[k][l])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
          }
        }
      }
    }
    if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kV0 && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kCascade) {
      /// V0-Cascade combination
      for (int k = 0; k < 3; k++) {
        std::string dirName = static_cast<std::string>(DirNames[7]);
        histdetadpisame[k][0] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesSame[0][k])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpisame[k][1] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesSame[1][k])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpimixed[k][0] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesMixed[0][k])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpimixed[k][1] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesMixed[1][k])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        if (plotForEveryRadii) {
          for (int l = 0; l < 9; l++) {
            histdetadpiRadii[k][l] = mHistogramRegistryQA->add<TH2>((dirName + static_cast<std::string>(HistNamesRadii[k][l])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
          }
        }
      }
    }
    if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kPhi) {
      for (int i = 0; i < 2; i++) {
        std::string dirName = static_cast<std::string>(DirNames[3]);
        histdetadpisame[i][0] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesSame[0][i])).c_str(), "; #Delta #eta; #Delta #varphi*", kTH2F, {{400, -0.30, 0.30}, {400, -0.30, 0.30}});
        histdetadpisame[i][1] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesSame[1][i])).c_str(), "; #Delta #eta; #Delta #varphi*", kTH2F, {{400, -0.30, 0.30}, {400, -0.30, 0.30}});
        histdetadpimixed[i][0] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesMixed[0][i])).c_str(), "; #Delta #eta; #Delta #varphi*", kTH2F, {{400, -0.30, 0.30}, {400, -0.30, 0.30}});
        histdetadpimixed[i][1] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesMixed[1][i])).c_str(), "; #Delta #eta; #Delta #varphi*", kTH2F, {{400, -0.30, 0.30}, {400, -0.30, 0.30}});

        if (plotForEveryRadii) {
          for (int j = 0; j < 9; j++) {
            histdetadpiRadii[i][j] = mHistogramRegistryQA->add<TH2>((dirName + static_cast<std::string>(HistNamesRadii[i][j])).c_str(), "; #Delta #eta; #Delta #varphi*", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
          }
        }
      }
    }
    if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kD0) {
      for (int i = 0; i < 2; i++) {
        std::string dirName = static_cast<std::string>(DirNames[4]);
        histdetadpisame[i][0] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesSame[0][i])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpisame[i][1] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesSame[1][i])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpimixed[i][0] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesMixed[0][i])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpimixed[i][1] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(HistNamesMixed[1][i])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});

        if (plotForEveryRadii) {
          for (int j = 0; j < 9; j++) {
            histdetadpiRadii[i][j] = mHistogramRegistryQA->add<TH2>((dirName + static_cast<std::string>(HistNamesRadii[i][j])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
          }
        }
      }
    }
  }

  template <typename t1>
  void init_kT(HistogramRegistry* registry, t1& ktbins, std::vector<float> ldeltaphistarcutmin, std::vector<float> ldeltaphistarcutmax, std::vector<float> ldeltaetacutmin, std::vector<float> ldeltaetacutmax)
  {
    mHistogramRegistry = registry;
    ktBins = ktbins;

    cutDeltaPhiStarMaxVector = ldeltaphistarcutmax;
    cutDeltaPhiStarMinVector = ldeltaphistarcutmin;
    cutDeltaEtaMaxVector = ldeltaetacutmax;
    cutDeltaEtaMinVector = ldeltaetacutmin;

    if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kTrack) {
      std::string dirName = static_cast<std::string>(DirNames[0]);
      for (int j = 1; j < static_cast<int>(ktBins.size() - 1); j++) {
        std::string histSuffixkT1 = std::to_string(static_cast<int>(ktBins[j] * 100.0));
        std::string histSuffixkT2 = std::to_string(static_cast<int>(ktBins[j + 1] * 100.0));
        std::string histFolderkT = "kT_" + histSuffixkT1 + "_" + histSuffixkT2 + "/";
        histdetadphisamebeforekT[j] = mHistogramRegistry->add<TH2>((dirName + histFolderkT + "detadphidetadphiBeforeSame").c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadphimixedbeforekT[j] = mHistogramRegistry->add<TH2>((dirName + histFolderkT + "detadphidetadphiBeforeMixed").c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadphisameafterkT[j] = mHistogramRegistry->add<TH2>((dirName + histFolderkT + "detadphidetadphiAfterSame").c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadphimixedafterkT[j] = mHistogramRegistry->add<TH2>((dirName + histFolderkT + "detadphidetadphiAfterMixed").c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
      }
    }
  }

  ///  Check if pair is close or not
  template <typename Part, typename Parts>
  bool isClosePair(Part const& part1, Part const& part2, Parts const& particles, float lmagfield, uint8_t ChosenEventType)
  {
    magfield = lmagfield;

    if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kTrack) {
      /// Track-Track combination
      // check if provided particles are in agreement with the class instantiation
      if (part1.partType() != o2::aod::femtouniverseparticle::ParticleType::kTrack || part2.partType() != o2::aod::femtouniverseparticle::ParticleType::kTrack) {
        LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar instantiation! Please provide kTrack,kTrack candidates.";
        return false;
      }
      auto deta = part1.eta() - part2.eta();
      auto dphiAvg = averagePhiStar(part1, part2, 0);
      if (ChosenEventType == femto_universe_container::EventType::same) {
        histdetadpisame[0][0]->Fill(deta, dphiAvg);
      } else if (ChosenEventType == femto_universe_container::EventType::mixed) {
        histdetadpimixed[0][0]->Fill(deta, dphiAvg);
      } else {
        LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar's type of events! Please provide same or mixed.";
      }

      if (std::pow(dphiAvg, 2) / std::pow(cutDeltaPhiStarMax, 2) + std::pow(deta, 2) / std::pow(cutDeltaEtaMax, 2) < 1.) {
        return true;
      } else {
        if (ChosenEventType == femto_universe_container::EventType::same) {
          histdetadpisame[0][1]->Fill(deta, dphiAvg);
        } else if (ChosenEventType == femto_universe_container::EventType::mixed) {
          histdetadpimixed[0][1]->Fill(deta, dphiAvg);
        } else {
          LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar's type of events! Please provide same or mixed.";
        }
        return false;
      }

    } else if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kV0) {
      /// Track-V0 combination
      // check if provided particles are in agreement with the class instantiation
      if (part1.partType() != o2::aod::femtouniverseparticle::ParticleType::kTrack || part2.partType() != o2::aod::femtouniverseparticle::ParticleType::kV0) {
        LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar instantiation! Please provide kTrack,kV0 candidates.";
        return false;
      }

      bool pass = false;
      for (int i = 0; i < 2; i++) {
        auto indexOfDaughter = (ChosenEventType == femto_universe_container::EventType::mixed ? part2.globalIndex() : part2.index()) - 2 + i;
        auto daughter = particles.begin() + indexOfDaughter;
        auto deta = part1.eta() - daughter.eta();
        auto dphiAvg = averagePhiStar(part1, *daughter, i);
        if (ChosenEventType == femto_universe_container::EventType::same) {
          histdetadpisame[i][0]->Fill(deta, dphiAvg);
        } else if (ChosenEventType == femto_universe_container::EventType::mixed) {
          histdetadpimixed[i][0]->Fill(deta, dphiAvg);
        } else {
          LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar's type of events! Please provide same or mixed.";
        }

        if (std::pow(dphiAvg, 2) / std::pow(cutDeltaPhiStarMax, 2) + std::pow(deta, 2) / std::pow(cutDeltaEtaMax, 2) < 1.) {
          pass = true;
        } else {
          if (ChosenEventType == femto_universe_container::EventType::same) {
            histdetadpisame[i][1]->Fill(deta, dphiAvg);
          } else if (ChosenEventType == femto_universe_container::EventType::mixed) {
            histdetadpimixed[i][1]->Fill(deta, dphiAvg);
          } else {
            LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar's type of events! Please provide same or mixed.";
          }
        }
      }
      return pass;

    } else if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kV0 && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kV0) {
      /// V0-V0 combination
      if (part1.partType() != o2::aod::femtouniverseparticle::ParticleType::kV0 || part2.partType() != o2::aod::femtouniverseparticle::ParticleType::kV0) {
        LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar instantiation! Please provide kV0,kV0 candidates.";
        return false;
      }

      bool pass = false;
      for (int i = 0; i < 2; i++) {
        auto indexOfDaughterpart1 = (ChosenEventType == femto_universe_container::EventType::mixed ? part1.globalIndex() : part1.index()) - 2 + i;
        auto indexOfDaughterpart2 = (ChosenEventType == femto_universe_container::EventType::mixed ? part2.globalIndex() : part2.index()) - 2 + i;
        auto daughterpart1 = particles.begin() + indexOfDaughterpart1;
        auto daughterpart2 = particles.begin() + indexOfDaughterpart2;
        auto deta = daughterpart1.eta() - daughterpart2.eta();
        auto dphiAvg = averagePhiStar(*daughterpart1, *daughterpart2, i);
        if (ChosenEventType == femto_universe_container::EventType::same) {
          histdetadpisame[i][0]->Fill(deta, dphiAvg);
        } else if (ChosenEventType == femto_universe_container::EventType::mixed) {
          histdetadpimixed[i][0]->Fill(deta, dphiAvg);
        } else {
          LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar's type of events! Please provide same or mixed.";
        }

        // if (std::pow(dphiAvg, 2) / std::pow(cutDeltaPhiStarMax, 2) + std::pow(deta, 2) / std::pow(cutDeltaEtaMax, 2) < 1.) {
        if ((dphiAvg > cutDeltaPhiStarMin) && (dphiAvg < cutDeltaPhiStarMax) && (deta > cutDeltaEtaMin) && (deta < cutDeltaEtaMax)) {
          pass = true;
        } else {
          if (ChosenEventType == femto_universe_container::EventType::same) {
            histdetadpisame[i][1]->Fill(deta, dphiAvg);
          } else if (ChosenEventType == femto_universe_container::EventType::mixed) {
            histdetadpimixed[i][1]->Fill(deta, dphiAvg);
          } else {
            LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar's type of events! Please provide same or mixed.";
          }
        }
      }
      return pass;

    } else if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kCascade && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kCascade) {
      /// Cascade-Cascade combination
      if (part1.partType() != o2::aod::femtouniverseparticle::ParticleType::kCascade || part2.partType() != o2::aod::femtouniverseparticle::ParticleType::kCascade) {
        LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar instantiation! Please provide kCascade,kCascade candidates.";
        return false;
      }

      bool pass = false;
      static constexpr int CascChildTable[][2] = {{-1, -1}, {-1, -2}, {-1, -3}, {-2, -2}, {-3, -3}, {-2, -1}, {-3, -1}};
      for (int i = 0; i < 7; i++) {
        auto indexOfDaughterpart1 = (ChosenEventType == femto_universe_container::EventType::mixed ? part1.globalIndex() : part1.index()) + CascChildTable[i][0];
        auto indexOfDaughterpart2 = (ChosenEventType == femto_universe_container::EventType::mixed ? part2.globalIndex() : part2.index()) + CascChildTable[i][1];
        auto daughterpart1 = particles.begin() + indexOfDaughterpart1;
        auto daughterpart2 = particles.begin() + indexOfDaughterpart2;
        if (isSameSignCPR && (daughterpart1.mAntiLambda() != daughterpart2.mAntiLambda())) // mAntiLambda() is used here as sign getter
          continue;
        auto deta = daughterpart1.eta() - daughterpart2.eta();
        auto dphiAvg = averagePhiStar(*daughterpart1, *daughterpart2, i);
        if (ChosenEventType == femto_universe_container::EventType::same) {
          histdetadpisame[i][0]->Fill(deta, dphiAvg);
        } else if (ChosenEventType == femto_universe_container::EventType::mixed) {
          histdetadpimixed[i][0]->Fill(deta, dphiAvg);
        } else {
          LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar's type of events! Please provide same or mixed.";
        }

        if ((dphiAvg > cutDeltaPhiStarMin) && (dphiAvg < cutDeltaPhiStarMax) && (deta > cutDeltaEtaMin) && (deta < cutDeltaEtaMax)) {
          pass = true;
        } else {
          if (ChosenEventType == femto_universe_container::EventType::same) {
            histdetadpisame[i][1]->Fill(deta, dphiAvg);
          } else if (ChosenEventType == femto_universe_container::EventType::mixed) {
            histdetadpimixed[i][1]->Fill(deta, dphiAvg);
          } else {
            LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar's type of events! Please provide same or mixed.";
          }
        }
      }
      return pass;

    } else if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kCascade) {
      /// Track-Cascade combination
      if (part1.partType() != o2::aod::femtouniverseparticle::ParticleType::kTrack || part2.partType() != o2::aod::femtouniverseparticle::ParticleType::kCascade) {
        LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar instantiation! Please provide kTrack,kCascade candidates.";
        return false;
      }

      bool pass = false;
      for (int i = 0; i < 3; i++) {
        auto indexOfDaughter = (ChosenEventType == femto_universe_container::EventType::mixed ? part2.globalIndex() : part2.index()) - 3 + i;
        auto daughter = particles.begin() + indexOfDaughter;
        if (isSameSignCPR && (part1.mAntiLambda() != daughter.mAntiLambda())) // mAntiLambda() is used here as sign getter
          continue;
        auto deta = part1.eta() - daughter.eta();
        auto dphiAvg = averagePhiStar(*part1, *daughter, i);
        if (ChosenEventType == femto_universe_container::EventType::same) {
          histdetadpisame[i][0]->Fill(deta, dphiAvg);
        } else if (ChosenEventType == femto_universe_container::EventType::mixed) {
          histdetadpimixed[i][0]->Fill(deta, dphiAvg);
        } else {
          LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar's type of events! Please provide same or mixed.";
        }

        if ((dphiAvg > cutDeltaPhiStarMin) && (dphiAvg < cutDeltaPhiStarMax) && (deta > cutDeltaEtaMin) && (deta < cutDeltaEtaMax)) {
          pass = true;
        } else {
          if (ChosenEventType == femto_universe_container::EventType::same) {
            histdetadpisame[i][1]->Fill(deta, dphiAvg);
          } else if (ChosenEventType == femto_universe_container::EventType::mixed) {
            histdetadpimixed[i][1]->Fill(deta, dphiAvg);
          } else {
            LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar's type of events! Please provide same or mixed.";
          }
        }
      }
      return pass;

    } else if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kV0 && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kCascade) {
      /// V0-Cascade combination
      if (part1.partType() != o2::aod::femtouniverseparticle::ParticleType::kV0 || part2.partType() != o2::aod::femtouniverseparticle::ParticleType::kCascade) {
        LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar instantiation! Please provide kV0,kCascade candidates.";
        return false;
      }

      bool pass = false;
      static constexpr int V0CascChildTable[][2] = {{-1, -1}, {-1, -2}, {-1, -3}, {-2, -1}, {-2, -2}, {-2, -3}};
      for (int i = 0; i < 3; i++) {
        auto indexOfDaughterV0 = (ChosenEventType == femto_universe_container::EventType::mixed ? part1.globalIndex() : part1.index()) + V0CascChildTable[i][0];
        auto indexOfDaughterCasc = (ChosenEventType == femto_universe_container::EventType::mixed ? part2.globalIndex() : part2.index()) + V0CascChildTable[i][1];
        auto daughterV0 = particles.begin() + indexOfDaughterV0;
        auto daughterCasc = particles.begin() + indexOfDaughterCasc;
        if (isSameSignCPR && (daughterV0.mAntiLambda() != daughterCasc.mAntiLambda())) // mAntiLambda() is used here as sign getter
          continue;
        auto deta = daughterV0.eta() - daughterCasc.eta();
        auto dphiAvg = averagePhiStar(*daughterV0, *daughterCasc, i);
        if (ChosenEventType == femto_universe_container::EventType::same) {
          histdetadpisame[i][0]->Fill(deta, dphiAvg);
        } else if (ChosenEventType == femto_universe_container::EventType::mixed) {
          histdetadpimixed[i][0]->Fill(deta, dphiAvg);
        } else {
          LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar's type of events! Please provide same or mixed.";
        }

        if ((dphiAvg > cutDeltaPhiStarMin) && (dphiAvg < cutDeltaPhiStarMax) && (deta > cutDeltaEtaMin) && (deta < cutDeltaEtaMax)) {
          pass = true;
        } else {
          if (ChosenEventType == femto_universe_container::EventType::same) {
            histdetadpisame[i][1]->Fill(deta, dphiAvg);
          } else if (ChosenEventType == femto_universe_container::EventType::mixed) {
            histdetadpimixed[i][1]->Fill(deta, dphiAvg);
          } else {
            LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar's type of events! Please provide same or mixed.";
          }
        }
      }
      return pass;

    } else if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kD0) {
      /// Track-D0 combination
      // check if provided particles are in agreement with the class instantiation
      if (part1.partType() != o2::aod::femtouniverseparticle::ParticleType::kTrack || part2.partType() != o2::aod::femtouniverseparticle::ParticleType::kD0) {
        LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar instantiation! Please provide kTrack, kD0 candidates.";
        return false;
      }

      bool pass = false;
      for (int i = 0; i < 2; i++) {
        auto indexOfDaughter = 0;
        if (ChosenEventType == femto_universe_container::EventType::mixed) {
          indexOfDaughter = part2.globalIndex() - 2 + i;
        } else if (ChosenEventType == femto_universe_container::EventType::same) {
          indexOfDaughter = part2.index() - 2 + i;
        }

        auto daughter = particles.begin() + indexOfDaughter;
        auto deta = part1.eta() - daughter.eta();
        auto dphiAvg = averagePhiStar(part1, *daughter, i); // auto dphiAvg = calculateDphiStar(part1, *daughter);
        dphiAvg = TVector2::Phi_mpi_pi(dphiAvg);
        if (ChosenEventType == femto_universe_container::EventType::same) {
          histdetadpisame[i][0]->Fill(deta, dphiAvg);
        } else if (ChosenEventType == femto_universe_container::EventType::mixed) {
          histdetadpimixed[i][0]->Fill(deta, dphiAvg);
        } else {
          LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar's type of events! Please provide same or mixed.";
        }

        if ((dphiAvg > cutDeltaPhiStarMin) && (dphiAvg < cutDeltaPhiStarMax) && (deta > cutDeltaEtaMin) && (deta < cutDeltaEtaMax)) {
          pass = true; // pair is close
        } else {
          if (ChosenEventType == femto_universe_container::EventType::same) {
            histdetadpisame[i][1]->Fill(deta, dphiAvg);
          } else if (ChosenEventType == femto_universe_container::EventType::mixed) {
            histdetadpimixed[i][1]->Fill(deta, dphiAvg);
          } else {
            LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar's type of events! Please provide same or mixed.";
          }
        }
      }
      return pass;
    } else if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kPhi) {
      /// Track-Phi combination
      // check if provided particles are in agreement with the class instantiation
      if (part1.partType() != o2::aod::femtouniverseparticle::ParticleType::kTrack || part2.partType() != o2::aod::femtouniverseparticle::ParticleType::kPhi) {
        LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar instantiation! Please provide kTrack,kPhi candidates.";
        return false;
      }

      bool pass = false;
      for (int i = 0; i < 2; i++) {
        auto indexOfDaughter = 0;
        if (ChosenEventType == femto_universe_container::EventType::mixed) {
          indexOfDaughter = part2.globalIndex() - 2 + i;
        } else if (ChosenEventType == femto_universe_container::EventType::same) {
          indexOfDaughter = part2.index() - 2 + i;
        }

        auto daughter = particles.begin() + indexOfDaughter;
        auto deta = part1.eta() - daughter.eta();
        auto dphiAvg = averagePhiStar(part1, *daughter, i); // calculateDphiStar(part1, *daughter);
        dphiAvg = TVector2::Phi_mpi_pi(dphiAvg);
        if (ChosenEventType == femto_universe_container::EventType::same) {
          histdetadpisame[i][0]->Fill(deta, dphiAvg);
        } else if (ChosenEventType == femto_universe_container::EventType::mixed) {
          histdetadpimixed[i][0]->Fill(deta, dphiAvg);
        } else {
          LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar's type of events! Please provide same or mixed.";
        }

        // REMOVING THE "RING" -- CALCULATING THE INVARIANT MASS
        TLorentzVector part1Vec;
        TLorentzVector part2Vec;
        float mMassOne = o2::constants::physics::MassKPlus;
        float mMassTwo = o2::constants::physics::MassKMinus;
        part1Vec.SetPtEtaPhiM(part1.pt(), part1.eta(), part1.phi(), mMassOne);
        part2Vec.SetPtEtaPhiM(daughter.pt(), daughter.eta(), daughter.phi(), mMassTwo);
        TLorentzVector sumVec(part1Vec);
        sumVec += part2Vec;
        float phiM = sumVec.M();
        if ((phiM > cutPhiInvMassLow) && (phiM < cutPhiInvMassHigh)) {
          pass = true; // pair comes from Phi meson decay
        }

        // APPLYING THE CUTS
        if ((dphiAvg > cutDeltaPhiStarMin) && (dphiAvg < cutDeltaPhiStarMax) && (deta > cutDeltaEtaMin) && (deta < cutDeltaEtaMax)) {
          pass = true; // pair is close
        } else {
          if (ChosenEventType == femto_universe_container::EventType::same) {
            histdetadpisame[i][1]->Fill(deta, dphiAvg);
          } else if (ChosenEventType == femto_universe_container::EventType::mixed) {
            histdetadpimixed[i][1]->Fill(deta, dphiAvg);
          } else {
            LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar's type of events! Please provide same or mixed.";
          }
        }
      }
      return pass;
    } else {
      LOG(fatal) << "FemtoUniversePairCleaner: Combination of objects not defined - quitting!";
      return false;
    }
  }

  ///  Check if pair is close or not
  template <typename Part>
  bool isClosePairAtITS(Part const& part1, Part const& part2, float lmagfield, uint8_t ChosenEventType)
  {
    magfield = lmagfield;

    if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kTrack) {
      /// Track-Track combination
      // check if provided particles are in agreement with the class instantiation
      if (part1.partType() != o2::aod::femtouniverseparticle::ParticleType::kTrack || part2.partType() != o2::aod::femtouniverseparticle::ParticleType::kTrack) {
        LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar instantiation! Please provide kTrack,kTrack candidates.";
        return false;
      }
      auto deta = part1.eta() - part2.eta();
      auto dphiAvg = part1.phi() - part2.phi();
      if (ChosenEventType == femto_universe_container::EventType::same) {
        histdetadpisame[0][0]->Fill(deta, dphiAvg);
      } else if (ChosenEventType == femto_universe_container::EventType::mixed) {
        histdetadpimixed[0][0]->Fill(deta, dphiAvg);
      } else {
        LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar's type of events! Please provide same or mixed.";
      }

      if (std::pow(dphiAvg, 2) / std::pow(cutDeltaPhiStarMax, 2) + std::pow(deta, 2) / std::pow(cutDeltaEtaMax, 2) < 1.) {
        return true;
      } else {
        if (ChosenEventType == femto_universe_container::EventType::same) {
          histdetadpisame[0][1]->Fill(deta, dphiAvg);
        } else if (ChosenEventType == femto_universe_container::EventType::mixed) {
          histdetadpimixed[0][1]->Fill(deta, dphiAvg);
        } else {
          LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar's type of events! Please provide same or mixed.";
        }
        return false;
      }
    } else {
      LOG(fatal) << "FemtoUniversePairCleaner: Combination of objects not defined - quitting!";
      return false;
    }
  }

  template <typename Part>
  bool isClosePairFrac(Part const& part1, Part const& part2, float lmagfield, uint8_t ChosenEventType, bool IsDphiAvgOrDist, float DistMax, float FracMax, bool CircCut)
  {
    magfield = lmagfield;

    if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kTrack) {
      /// Track-Track combination
      // check if provided particles are in agreement with the class instantiation
      if (part1.partType() != o2::aod::femtouniverseparticle::ParticleType::kTrack || part2.partType() != o2::aod::femtouniverseparticle::ParticleType::kTrack) {
        LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar instantiation! Please provide kTrack,kTrack candidates.";
        return false;
      }
      auto deta = part1.eta() - part2.eta();
      auto dphiAvg = averagePhiStar(part1, part2, 0);
      auto distfrac = averagePhiStarFrac(part1, part2, DistMax);
      if (ChosenEventType == femto_universe_container::EventType::same) {
        histdetadpisame[0][0]->Fill(deta, dphiAvg);
      } else if (ChosenEventType == femto_universe_container::EventType::mixed) {
        histdetadpimixed[0][0]->Fill(deta, dphiAvg);
      } else {
        LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar's type of events! Please provide same or mixed.";
      }

      if (IsDphiAvgOrDist) {
        if (CircCut && (std::pow(dphiAvg, 2) / std::pow(cutDeltaPhiStarMax, 2) + std::pow(deta, 2) / std::pow(cutDeltaEtaMax, 2) < 1.)) {
          return true;
        } else if (!CircCut && (dphiAvg > cutDeltaPhiStarMin) && (dphiAvg < cutDeltaPhiStarMax) && (deta > cutDeltaEtaMin) && (deta < cutDeltaEtaMax)) {
          return true;
        } else {
          if (ChosenEventType == femto_universe_container::EventType::same) {
            histdetadpisame[0][1]->Fill(deta, dphiAvg);
          } else if (ChosenEventType == femto_universe_container::EventType::mixed) {
            histdetadpimixed[0][1]->Fill(deta, dphiAvg);
          } else {
            LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar's type of events! Please provide same or mixed.";
          }
          return false;
        }
      } else {
        if (((deta > cutDeltaEtaMin) && (deta < cutDeltaEtaMax)) && (distfrac > FracMax)) {
          return true;
        } else {
          if (ChosenEventType == femto_universe_container::EventType::same) {
            histdetadpisame[0][1]->Fill(deta, dphiAvg);
          } else if (ChosenEventType == femto_universe_container::EventType::mixed) {
            histdetadpimixed[0][1]->Fill(deta, dphiAvg);
          } else {
            LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar's type of events! Please provide same or mixed.";
          }
          return false;
        }
      }

    } else {
      LOG(fatal) << "FemtoUniversePairCleaner: Combination of objects not defined - quitting!";
      return false;
    }
  }

  ///  Check if pair is close or not
  template <typename Part>
  bool isClosePairkT(Part const& part1, Part const& part2, uint8_t ChosenEventType, float ktval, bool CircCut)
  {
    /// Track-Track combination
    // check if provided particles are in agreement with the class instantiation
    if (part1.partType() != o2::aod::femtouniverseparticle::ParticleType::kTrack || part2.partType() != o2::aod::femtouniverseparticle::ParticleType::kTrack) {
      LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar instantiation! Please provide kTrack,kTrack candidates.";
      return false;
    }

    int ktbinval = 1;
    if (ktval >= ktBins[1] && ktval < ktBins[2]) {
      ktbinval = 1;
    } else if (ktval >= ktBins[2] && ktval < ktBins[3]) {
      ktbinval = 2;
    } else if (ktval >= ktBins[3] && ktval < ktBins[4]) {
      ktbinval = 3;
    } else if (ktval >= ktBins[4] && ktval < ktBins[5]) {
      ktbinval = 4;
    }

    auto deta = part1.eta() - part2.eta();
    auto dphiAvg = averagePhiStar(part1, part2, 0);
    if (ChosenEventType == femto_universe_container::EventType::same) {
      histdetadphisamebeforekT[ktbinval]->Fill(deta, dphiAvg);
    } else if (ChosenEventType == femto_universe_container::EventType::mixed) {
      histdetadphimixedbeforekT[ktbinval]->Fill(deta, dphiAvg);
    } else {
      LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar's type of events! Please provide same or mixed.";
    }

    if (CircCut && (std::pow(dphiAvg, 2) / std::pow(cutDeltaPhiStarMaxVector[ktbinval], 2) + std::pow(deta, 2) / std::pow(cutDeltaEtaMaxVector[ktbinval], 2) < 1.)) {
      return true;
    } else if (!CircCut && (dphiAvg > cutDeltaPhiStarMinVector[ktbinval]) && (dphiAvg < cutDeltaPhiStarMaxVector[ktbinval]) && (deta > cutDeltaEtaMinVector[ktbinval]) && (deta < cutDeltaEtaMaxVector[ktbinval])) {
      return true;
    } else {
      if (ChosenEventType == femto_universe_container::EventType::same) {
        histdetadphisameafterkT[ktbinval]->Fill(deta, dphiAvg);
      } else if (ChosenEventType == femto_universe_container::EventType::mixed) {
        histdetadphimixedafterkT[ktbinval]->Fill(deta, dphiAvg);
      } else {
        LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar's type of events! Please provide same or mixed.";
      }
      return false;
    }
  }

  ///  Check if pair is close or not
  template <typename Part>
  void ClosePairqLCMS(Part const& part1, Part const& part2, float lmagfield, uint8_t ChosenEventType, double qlcms) // add typename Parts and variable parts for adding MClabels
  {
    magfield = lmagfield;
    if constexpr (kPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && kPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kTrack) {
      auto deta = part1.eta() - part2.eta();
      auto dphiAvg = averagePhiStar(part1, part2, 0);

      if (ChosenEventType == femto_universe_container::EventType::same) {
        histdetadpiqlcmssame->Fill(qlcms, deta, dphiAvg);
      } else if (ChosenEventType == femto_universe_container::EventType::mixed) {
        histdetadpiqlcmsmixed->Fill(qlcms, deta, dphiAvg);
      } else {
        LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar's type of events! Please provide same or mixed.";
      }
    }
  }

 private:
  HistogramRegistry* mHistogramRegistry = nullptr;   ///< For main output
  HistogramRegistry* mHistogramRegistryQA = nullptr; ///< For QA output
  static constexpr std::string_view DirNames[8] = {"kTrack_kTrack/", "kTrack_kV0/", "kV0_kV0/", "kTrack_kPhi/", "kTrack_kD0/", "kCascade_kCascade/", "kTrack_kCascade/", "kV0_kCascade/"};

  static constexpr std::string_view HistNamesSame[2][8] = {{"detadphidetadphi0BeforeSame_0", "detadphidetadphi0BeforeSame_1", "detadphidetadphi0BeforeSame_2",
                                                            "detadphidetadphi0BeforeSame_3", "detadphidetadphi0BeforeSame_4", "detadphidetadphi0BeforeSame_5",
                                                            "detadphidetadphi0BeforeSame_6", "detadphidetadphi0BeforeSameqLCMS"},
                                                           {"detadphidetadphi0AfterSame_0", "detadphidetadphi0AfterSame_1", "detadphidetadphi0AfterSame_2",
                                                            "detadphidetadphi0AfterSame_3", "detadphidetadphi0AfterSame_4", "detadphidetadphi0AfterSame_5",
                                                            "detadphidetadphi0AfterSame_6", "detadphidetadphi0AfterSameqLCMS"}};
  static constexpr std::string_view HistNamesMixed[2][8] = {{"detadphidetadphi0BeforeMixed_0", "detadphidetadphi0BeforeMixed_1", "detadphidetadphi0BeforeMixed_2",
                                                             "detadphidetadphi0BeforeMixed_3", "detadphidetadphi0BeforeMixed_4", "detadphidetadphi0BeforeMixed_5",
                                                             "detadphidetadphi0BeforeMixed_6", "detadphidetadphi0BeforeMixedqLCMS"},
                                                            {"detadphidetadphi0AfterMixed_0", "detadphidetadphi0AfterMixed_1", "detadphidetadphi0AfterMixed_2",
                                                             "detadphidetadphi0AfterMixed_3", "detadphidetadphi0AfterMixed_4", "detadphidetadphi0AfterMixed_5",
                                                             "detadphidetadphi0AfterMixed_6", "detadphidetadphi0AfterMixedqLCMS"}};

  static constexpr std::string_view HistNamesRadii[7][9] = {{"detadphidetadphi0Before_0_0", "detadphidetadphi0Before_0_1", "detadphidetadphi0Before_0_2",
                                                             "detadphidetadphi0Before_0_3", "detadphidetadphi0Before_0_4", "detadphidetadphi0Before_0_5",
                                                             "detadphidetadphi0Before_0_6", "detadphidetadphi0Before_0_7", "detadphidetadphi0Before_0_8"},
                                                            {"detadphidetadphi0Before_1_0", "detadphidetadphi0Before_1_1", "detadphidetadphi0Before_1_2",
                                                             "detadphidetadphi0Before_1_3", "detadphidetadphi0Before_1_4", "detadphidetadphi0Before_1_5",
                                                             "detadphidetadphi0Before_1_6", "detadphidetadphi0Before_1_7", "detadphidetadphi0Before_1_8"},
                                                            {"detadphidetadphi0Before_2_0", "detadphidetadphi0Before_2_1", "detadphidetadphi0Before_2_2",
                                                             "detadphidetadphi0Before_2_3", "detadphidetadphi0Before_2_4", "detadphidetadphi0Before_2_5",
                                                             "detadphidetadphi0Before_2_6", "detadphidetadphi0Before_2_7", "detadphidetadphi0Before_2_8"},
                                                            {"detadphidetadphi0Before_3_0", "detadphidetadphi0Before_3_1", "detadphidetadphi0Before_3_2",
                                                             "detadphidetadphi0Before_3_3", "detadphidetadphi0Before_3_4", "detadphidetadphi0Before_3_5",
                                                             "detadphidetadphi0Before_3_6", "detadphidetadphi0Before_3_7", "detadphidetadphi0Before_3_8"},
                                                            {"detadphidetadphi0Before_4_0", "detadphidetadphi0Before_4_1", "detadphidetadphi0Before_4_2",
                                                             "detadphidetadphi0Before_4_3", "detadphidetadphi0Before_4_4", "detadphidetadphi0Before_4_5",
                                                             "detadphidetadphi0Before_4_6", "detadphidetadphi0Before_4_7", "detadphidetadphi0Before_4_8"},
                                                            {"detadphidetadphi0Before_5_0", "detadphidetadphi0Before_5_1", "detadphidetadphi0Before_5_2",
                                                             "detadphidetadphi0Before_5_3", "detadphidetadphi0Before_5_4", "detadphidetadphi0Before_5_5",
                                                             "detadphidetadphi0Before_5_6", "detadphidetadphi0Before_5_7", "detadphidetadphi0Before_5_8"},
                                                            {"detadphidetadphi0Before_6_0", "detadphidetadphi0Before_6_1", "detadphidetadphi0Before_6_2",
                                                             "detadphidetadphi0Before_6_3", "detadphidetadphi0Before_6_4", "detadphidetadphi0Before_6_5",
                                                             "detadphidetadphi0Before_6_6", "detadphidetadphi0Before_6_7", "detadphidetadphi0Before_6_8"}};

  static constexpr o2::aod::femtouniverseparticle::ParticleType kPartOneType = partOne; ///< Type of particle 1
  static constexpr o2::aod::femtouniverseparticle::ParticleType kPartTwoType = partTwo; ///< Type of particle 2

  static constexpr float TmpRadiiTPC[9] = {85., 105., 125., 145., 165., 185., 205., 225., 245.};

  static constexpr uint32_t kSignMinusMask = 1;
  static constexpr uint32_t kSignPlusMask = 1 << 1;
  static constexpr uint32_t kValue0 = 0;

  float chosenRadii;
  float cutDeltaPhiStarMax;
  float cutDeltaPhiStarMin;
  float cutDeltaEtaMax;
  float cutDeltaEtaMin;

  std::vector<float> cutDeltaPhiStarMaxVector;
  std::vector<float> cutDeltaPhiStarMinVector;
  std::vector<float> cutDeltaEtaMaxVector;
  std::vector<float> cutDeltaEtaMinVector;

  float magfield;
  bool plotForEveryRadii = false;
  float cutPhiInvMassLow;
  float cutPhiInvMassHigh;
  bool isSameSignCPR = false;
  std::vector<double> ktBins;

  std::array<std::array<std::shared_ptr<TH2>, 2>, 7> histdetadpisame{};
  std::array<std::array<std::shared_ptr<TH2>, 2>, 7> histdetadpimixed{};
  std::array<std::shared_ptr<TH2>, 4> histdetadphisamebeforekT{};
  std::array<std::shared_ptr<TH2>, 4> histdetadphimixedbeforekT{};
  std::array<std::shared_ptr<TH2>, 4> histdetadphisameafterkT{};
  std::array<std::shared_ptr<TH2>, 4> histdetadphimixedafterkT{};

  std::array<std::array<std::shared_ptr<TH2>, 9>, 7> histdetadpiRadii{};

  std::shared_ptr<TH3> histdetadpiqlcmssame{};
  std::shared_ptr<TH3> histdetadpiqlcmsmixed{};

  ///  Calculate phi at all required radii stored in TmpRadiiTPC
  /// Magnetic field to be provided in Tesla
  template <typename T>
  void phiAtRadiiTPC(const T& part, std::vector<float>& tmpVec)
  {

    float phi0 = part.phi();
    // Start: Get the charge from cutcontainer using masks
    float charge = 0.;
    if ((part.cut() & kSignMinusMask) == kValue0 && (part.cut() & kSignPlusMask) == kValue0) {
      charge = 0;
    } else if ((part.cut() & kSignPlusMask) == kSignPlusMask) {
      charge = 1;
    } else if ((part.cut() & kSignMinusMask) == kSignMinusMask) {
      charge = -1;
    } else {
      LOG(fatal) << "FemtoUniverseDetaDphiStar: Charge bits are set wrong!";
    }
    // End: Get the charge from cutcontainer using masks
    float pt = part.pt();
    for (size_t i = 0; i < 9; i++) {
      double arg = 0.3 * charge * magfield * TmpRadiiTPC[i] * 0.01 / (2. * pt);
      if (std::abs(arg) < 1.0) {
        tmpVec.push_back(phi0 - std::asin(arg));
      } else {
        tmpVec.push_back(999.0);
      }
    }
  }

  ///  Calculate average phi
  template <typename T1, typename T2>
  float averagePhiStar(const T1& part1, const T2& part2, int iHist)
  {
    std::vector<float> tmpVec1;
    std::vector<float> tmpVec2;
    phiAtRadiiTPC(part1, tmpVec1);
    phiAtRadiiTPC(part2, tmpVec2);
    int num = tmpVec1.size();
    float dPhiAvg = 0;
    float dphi = 0;
    int entries = 0;
    for (int i = 0; i < num; i++) {
      if (tmpVec1.at(i) != 999 && tmpVec2.at(i) != 999) {
        dphi = tmpVec1.at(i) - tmpVec2.at(i);
        entries++;
      } else {
        dphi = 0;
      }
      dphi = TVector2::Phi_mpi_pi(dphi);
      dPhiAvg += dphi;
      if (plotForEveryRadii) {
        histdetadpiRadii[iHist][i]->Fill(part1.eta() - part2.eta(), dphi);
      }
    }
    return dPhiAvg / static_cast<float>(entries);
  }

  ///  Calculate average phi
  template <typename T1, typename T2>
  float averagePhiStarFrac(const T1& part1, const T2& part2, float maxdist)
  {
    std::vector<float> tmpVec1;
    std::vector<float> tmpVec2;
    phiAtRadiiTPC(part1, tmpVec1);
    phiAtRadiiTPC(part2, tmpVec2);
    int num = tmpVec1.size();
    float dphi = 0;
    int entries = 0;
    double distance = 0;
    int badpoints = 0;

    for (int i = 0; i < num; i++) {
      if (tmpVec1.at(i) != 999 && tmpVec2.at(i) != 999) {
        dphi = tmpVec1.at(i) - tmpVec2.at(i);
        entries++;
      } else {
        dphi = 0;
      }
      dphi = TVector2::Phi_mpi_pi(dphi);
      distance = 2 * TMath::Sin(TMath::Abs(dphi) * 0.5) * TmpRadiiTPC[i];
      if (distance < maxdist) {
        badpoints++;
      }
    }
    return badpoints / entries;
  }

  // Get particle charge from mask
  template <typename T1>
  float getCharge(const T1& part)
  {
    float charge = 0;
    if ((part.cut() & kSignMinusMask) == kValue0 && (part.cut() & kSignPlusMask) == kValue0) {
      charge = 0;
    } else if ((part.cut() & kSignPlusMask) == kSignPlusMask) {
      charge = 1;
    } else if ((part.cut() & kSignMinusMask) == kSignMinusMask) {
      charge = -1;
    } else {
      LOG(fatal) << "FemtoUniverseDetaDphiStar: Charge bits are set wrong!";
    }
    return charge;
  }

  // Calculate phi* as in https://github.com/alisw/AliPhysics/blob/master/PWGCF/FEMTOSCOPY/AliFemtoUser/AliFemtoPairCutRadialDistance.cxx
  template <typename T1, typename T2>
  double calculateDphiStar(const T1& part1, const T2& part2)
  {
    float charge1 = getCharge(part1);
    float charge2 = getCharge(part2);

    double deltaphiconstFD = 0.3 / 2;
    // double deltaphiconstAF = 0.15;
    double afsi0b = deltaphiconstFD * magfield * charge1 * chosenRadii / part1.pt();
    double afsi1b = deltaphiconstFD * magfield * charge2 * chosenRadii / part2.pt();
    double dphis = 0.0;

    if (std::abs(afsi0b) < 1.0 && std::abs(afsi0b) < 1.0) {
      dphis = part2.phi() - part1.phi() + std::asin(afsi1b) - std::asin(afsi0b);
    }
    return dphis;
  }
};

} /* namespace femto_universe */
} /* namespace o2::analysis */

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEDETADPHISTAR_H_
