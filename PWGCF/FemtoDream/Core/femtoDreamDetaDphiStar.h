// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file FemtoDreamDetaDphiStar.h
/// \brief FemtoDreamDetaDphiStar - Checks particles for the close pair rejection.
/// \author Laura Serksnyte, TU München, laura.serksnyte@tum.de

#ifndef PWGCF_FEMTODREAM_CORE_FEMTODREAMDETADPHISTAR_H_
#define PWGCF_FEMTODREAM_CORE_FEMTODREAMDETADPHISTAR_H_

#include <memory>
#include <string>
#include <vector>
#include "PWGCF/DataModel/FemtoDerived.h"
#include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::framework;

namespace o2::analysis
{
namespace femtoDream
{

/// \class FemtoDreamDetaDphiStar
/// \brief Class to check particles for the close pair rejection.
/// \tparam partOne Type of particle 1 (Track/V0/Cascade/...)
/// \tparam partTwo Type of particle 2 (Track/V0/Cascade/...)

enum ProngCharmHadron {
  Prong0 = 0,
  Prong1 = 1,
  Prong2 = 2,
  Nprongs = 3
};

template <o2::aod::femtodreamparticle::ParticleType partOne, o2::aod::femtodreamparticle::ParticleType partTwo>
class FemtoDreamDetaDphiStar
{
 public:
  /// Destructor
  virtual ~FemtoDreamDetaDphiStar() = default;
  /// Initialization of the histograms and setting required values
  // atWhichRadiiToCut - at which radii apply deta dphi cut; 0 - PV; 1 - average phi at given tpc radii, 2 - at 80 cm
  void init(HistogramRegistry* registry, HistogramRegistry* registryQA, float ldeltaPhiMax, float ldeltaEtaMax, bool lplotForEveryRadii, int meORse = 0, bool oldversion = true, float Q3Limit = 8., bool isMELambda = false, int atWhichRadiiToCut = 1, float radiiTPCtoCut = 85., bool fillTHSparse = false)
  {
    deltaPhiMax = ldeltaPhiMax;
    deltaEtaMax = ldeltaEtaMax;
    plotForEveryRadii = lplotForEveryRadii;
    upperQ3LimitForPlotting = Q3Limit;
    isMixedEventLambda = isMELambda;
    runOldVersion = oldversion;
    mHistogramRegistry = registry;
    mHistogramRegistryQA = registryQA;
    atWhichRadiiToSelect = atWhichRadiiToCut;
    radiiTPC = radiiTPCtoCut;
    fillQA = fillTHSparse;

    if constexpr (mPartOneType == o2::aod::femtodreamparticle::ParticleType::kTrack && (mPartTwoType == o2::aod::femtodreamparticle::ParticleType::kTrack || mPartTwoType == o2::aod::femtodreamparticle::ParticleType::kCascadeV0Child || mPartTwoType == o2::aod::femtodreamparticle::ParticleType::kCascadeBachelor)) {
      std::string dirName = static_cast<std::string>(dirNames[0]);
      histdetadpi[0][0] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(histNames[0][0]) + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
      histdetadpi[0][1] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(histNames[1][0]) + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
      histdetadpi[0][2] = mHistogramRegistry->add<TH2>((dirName + "at_PV_before" + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
      histdetadpi[0][3] = mHistogramRegistry->add<TH2>((dirName + "at_PV_after" + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
      if (plotForEveryRadii) {
        for (int i = 0; i < 9; i++) {
          histdetadpiRadii[0][i] = mHistogramRegistryQA->add<TH2>((dirName + static_cast<std::string>(histNamesRadii[0][i]) + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        }
      }
      if (fillQA) {
        histdetadpi_eta[0] = mHistogramRegistry->add<THnSparse>((dirName + "dEtadPhi_Eta" + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}; #eta_{1}; #eta_{2}", kTHnSparseF, {{100, -0.15, 0.15}, {100, -0.15, 0.15}, {100, -0.8, 0.8}, {100, -0.8, 0.8}});
        histdetadpi_phi[0] = mHistogramRegistry->add<THnSparse>((dirName + "dEtadPhi_Phi" + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}; #phi_{1}; #phi_{2}", kTHnSparseF, {{100, -0.15, 0.15}, {100, -0.15, 0.15}, {100, 0, 6.28}, {100, 0, 6.28}});
      }
    }
    if constexpr (mPartOneType == o2::aod::femtodreamparticle::ParticleType::kTrack && mPartTwoType == o2::aod::femtodreamparticle::ParticleType::kV0) {
      for (int i = 0; i < 2; i++) {
        std::string dirName = static_cast<std::string>(dirNames[1]);
        histdetadpi[i][0] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(histNames[0][i]) + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpi[i][1] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(histNames[1][i]) + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpi[i][2] = mHistogramRegistry->add<TH2>((dirName + "at_PV_" + std::to_string(i) + "_before" + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpi[i][3] = mHistogramRegistry->add<TH2>((dirName + "at_PV_" + std::to_string(i) + "_after" + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        if (plotForEveryRadii) {
          for (int j = 0; j < 9; j++) {
            histdetadpiRadii[i][j] = mHistogramRegistryQA->add<TH2>((dirName + static_cast<std::string>(histNamesRadii[i][j]) + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
          }
        }
        if (fillQA) {
          histdetadpi_eta[i] = mHistogramRegistry->add<THnSparse>((dirName + "dEtadPhi_Eta_" + std::to_string(i) + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}; #eta_{1}; #eta_{2}", kTHnSparseF, {{100, -0.15, 0.15}, {100, -0.15, 0.15}, {100, -0.8, 0.8}, {100, -0.8, 0.8}});
          histdetadpi_phi[i] = mHistogramRegistry->add<THnSparse>((dirName + "dEtadPhi_Phi_" + std::to_string(i) + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}; #phi_{1}; #phi_{2}", kTHnSparseF, {{100, -0.15, 0.15}, {100, -0.15, 0.15}, {100, 0, 6.28}, {100, 0, 6.28}});
        }
      }
    }
    if constexpr (mPartOneType == o2::aod::femtodreamparticle::ParticleType::kTrack && mPartTwoType == o2::aod::femtodreamparticle::ParticleType::kCharmHadron) {
      for (int i = 0; i < Nprongs; i++) {
        std::string dirName = static_cast<std::string>(dirNames[2]);
        histdetadpi[i][0] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(histNames[0][i]) + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpi[i][1] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(histNames[1][i]) + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpi[i][2] = mHistogramRegistry->add<TH2>((dirName + "at_PV_" + std::to_string(i) + "_before" + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpi[i][3] = mHistogramRegistry->add<TH2>((dirName + "at_PV_" + std::to_string(i) + "_after" + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        if (plotForEveryRadii) {
          for (int j = 0; j < 9; j++) {
            histdetadpiRadii[i][j] = mHistogramRegistryQA->add<TH2>((dirName + static_cast<std::string>(histNamesRadii[i][j]) + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
          }
        }
        if (fillQA) {
          histdetadpi_eta[i] = mHistogramRegistry->add<THnSparse>((dirName + "dEtadPhi_Eta_" + std::to_string(i) + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}; #eta_{1}; #eta_{2}", kTHnSparseF, {{100, -0.15, 0.15}, {100, -0.15, 0.15}, {100, -0.8, 0.8}, {100, -0.8, 0.8}});
          histdetadpi_phi[i] = mHistogramRegistry->add<THnSparse>((dirName + "dEtadPhi_Phi_" + std::to_string(i) + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}; #phi_{1}; #phi_{2}", kTHnSparseF, {{100, -0.15, 0.15}, {100, -0.15, 0.15}, {100, 0, 6.28}, {100, 0, 6.28}});
        }
      }
    }
    if constexpr (mPartOneType == o2::aod::femtodreamparticle::ParticleType::kTrack && mPartTwoType == o2::aod::femtodreamparticle::ParticleType::kCascade) {
      for (int i = 0; i < 3; i++) {
        std::string dirName = static_cast<std::string>(dirNames[3]);
        histdetadpi[i][0] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(histNames[0][i]) + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpi[i][1] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(histNames[1][i]) + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpi[i][2] = mHistogramRegistry->add<TH2>((dirName + "at_PV_" + std::to_string(i) + "_before" + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpi[i][3] = mHistogramRegistry->add<TH2>((dirName + "at_PV_" + std::to_string(i) + "_after" + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        if (plotForEveryRadii) {
          for (int j = 0; j < 9; j++) {
            histdetadpiRadii[i][j] = mHistogramRegistryQA->add<TH2>((dirName + static_cast<std::string>(histNamesRadii[i][j]) + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
          }
        }
        if (fillQA) {
          histdetadpi_eta[i] = mHistogramRegistry->add<THnSparse>((dirName + "dEtadPhi_Eta_" + std::to_string(i) + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}; #eta_{1}; #eta_{2}", kTHnSparseF, {{100, -0.15, 0.15}, {100, -0.15, 0.15}, {100, -0.8, 0.8}, {100, -0.8, 0.8}});
          histdetadpi_phi[i] = mHistogramRegistry->add<THnSparse>((dirName + "dEtadPhi_Phi_" + std::to_string(i) + static_cast<std::string>(histNameSEorME[meORse])).c_str(), "; #Delta #eta; #Delta #phi^{*}; #phi_{1}; #phi_{2}", kTHnSparseF, {{100, -0.15, 0.15}, {100, -0.15, 0.15}, {100, 0, 6.28}, {100, 0, 6.28}});
        }
      }
    }
  }
  ///  Check if pair is close or not
  template <typename Part1, typename Part2, typename Parts>
  bool isClosePair(Part1 const& part1, Part2 const& part2, Parts const& particles, float lmagfield, float Q3 = 999.)
  {
    magfield = lmagfield;

    if constexpr (mPartOneType == o2::aod::femtodreamparticle::ParticleType::kTrack && (mPartTwoType == o2::aod::femtodreamparticle::ParticleType::kTrack || mPartTwoType == o2::aod::femtodreamparticle::ParticleType::kCascadeV0Child)) {
      /// Track-Track combination
      // check if provided particles are in agreement with the class instantiation
      if (part1.partType() != o2::aod::femtodreamparticle::ParticleType::kTrack || !(part2.partType() == o2::aod::femtodreamparticle::ParticleType::kTrack || part2.partType() == o2::aod::femtodreamparticle::ParticleType::kCascadeV0Child)) { // hotfix to use the CPR
        // LOG(fatal) << "FemtoDreamDetaDphiStar: passed arguments don't agree with FemtoDreamDetaDphiStar instantiation! Please provide kTrack,kTrack candidates.";
        LOGF(fatal, "FemtoDreamDetaDphiStar: passed arguments don't agree with FemtoDreamDetaDphiStar instantiation! Please provide kTrack,kTrack candidates. Currently: %i", part2.partType());
        return false;
      }
      auto deta = part1.eta() - part2.eta();
      auto dphi_AT_PV = part1.phi() - part2.phi();
      auto dphi_AT_SpecificRadii = PhiAtSpecificRadiiTPC(part1, radiiTPC) - PhiAtSpecificRadiiTPC(part2, radiiTPC);
      bool sameCharge = false;
      auto dphiAvg = AveragePhiStar(part1, part2, 0, &sameCharge);
      if (Q3 == 999) {
        histdetadpi[0][0]->Fill(deta, dphiAvg);
        histdetadpi[0][2]->Fill(deta, dphi_AT_PV);
        if (fillQA) {
          histdetadpi_eta[0]->Fill(deta, dphiAvg, part1.eta(), part2.eta());
          histdetadpi_phi[0]->Fill(deta, dphiAvg, part1.phi(), part2.phi());
        }
      } else if (Q3 < upperQ3LimitForPlotting) {
        histdetadpi[0][0]->Fill(deta, dphiAvg);
        histdetadpi[0][2]->Fill(deta, dphi_AT_PV);
        if (fillQA) {
          histdetadpi_eta[0]->Fill(deta, dphiAvg, part1.eta(), part2.eta());
          histdetadpi_phi[0]->Fill(deta, dphiAvg, part1.phi(), part2.phi());
        }
      }
      if (sameCharge) {
        if (atWhichRadiiToSelect == 1) {
          if (pow(dphiAvg, 2) / pow(deltaPhiMax, 2) + pow(deta, 2) / pow(deltaEtaMax, 2) < 1.) {
            return true;
          } else {
            if (Q3 == 999) {
              histdetadpi[0][1]->Fill(deta, dphiAvg);
              histdetadpi[0][3]->Fill(deta, dphi_AT_PV);
            } else if (Q3 < upperQ3LimitForPlotting) {
              histdetadpi[0][1]->Fill(deta, dphiAvg);
              histdetadpi[0][3]->Fill(deta, dphi_AT_PV);
            }
            return false;
          }
        } else if (atWhichRadiiToSelect == 0) {
          if (pow(dphi_AT_PV, 2) / pow(deltaPhiMax, 2) + pow(deta, 2) / pow(deltaEtaMax, 2) < 1.) {
            return true;
          } else {
            if (Q3 == 999) {
              histdetadpi[0][1]->Fill(deta, dphiAvg);
              histdetadpi[0][3]->Fill(deta, dphi_AT_PV);
            } else if (Q3 < upperQ3LimitForPlotting) {
              histdetadpi[0][1]->Fill(deta, dphiAvg);
              histdetadpi[0][3]->Fill(deta, dphi_AT_PV);
            }
            return false;
          }
        } else if (atWhichRadiiToSelect == 2) {
          if (pow(dphi_AT_SpecificRadii, 2) / pow(deltaPhiMax, 2) + pow(deta, 2) / pow(deltaEtaMax, 2) < 1.) {
            return true;
          } else {
            if (Q3 == 999) {
              histdetadpi[0][1]->Fill(deta, dphiAvg);
              histdetadpi[0][3]->Fill(deta, dphi_AT_PV);
            } else if (Q3 < upperQ3LimitForPlotting) {
              histdetadpi[0][1]->Fill(deta, dphiAvg);
              histdetadpi[0][3]->Fill(deta, dphi_AT_PV);
            }
            return false;
          }
        } else {
          return true;
        }
      } else {
        return false;
      }

    } else if constexpr (mPartOneType == o2::aod::femtodreamparticle::ParticleType::kTrack && mPartTwoType == o2::aod::femtodreamparticle::ParticleType::kV0) {
      /// Track-V0 combination
      // check if provided particles are in agreement with the class instantiation
      if (part1.partType() != o2::aod::femtodreamparticle::ParticleType::kTrack || part2.partType() != o2::aod::femtodreamparticle::ParticleType::kV0) {
        LOG(fatal) << "FemtoDreamDetaDphiStar: passed arguments don't agree with FemtoDreamDetaDphiStar instantiation! Please provide kTrack,kV0 candidates.";
        return false;
      }

      bool pass = false;
      for (int i = 0; i < 2; i++) {
        int indexOfDaughter;
        if (isMixedEventLambda) {
          indexOfDaughter = part2.globalIndex() - 2 + i;
        } else {
          indexOfDaughter = part2.index() - 2 + i;
        }
        auto daughter = particles.begin() + indexOfDaughter;
        auto deta = part1.eta() - daughter.eta();
        auto dphi_AT_PV = part1.phi() - daughter.phi();
        auto dphi_AT_SpecificRadii = PhiAtSpecificRadiiTPC(part1, radiiTPC) - PhiAtSpecificRadiiTPC(daughter, radiiTPC);
        bool sameCharge = false;
        auto dphiAvg = AveragePhiStar(part1, *daughter, i, &sameCharge);
        if (Q3 == 999) {
          histdetadpi[i][0]->Fill(deta, dphiAvg);
          histdetadpi[i][2]->Fill(deta, dphi_AT_PV);
          if (fillQA) {
            histdetadpi_eta[i]->Fill(deta, dphiAvg, part1.eta(), daughter.eta());
            histdetadpi_phi[i]->Fill(deta, dphiAvg, part1.phi(), daughter.phi());
          }
        } else if (Q3 < upperQ3LimitForPlotting) {
          histdetadpi[i][0]->Fill(deta, dphiAvg);
          histdetadpi[i][2]->Fill(deta, dphi_AT_PV);
          if (fillQA) {
            histdetadpi_eta[i]->Fill(deta, dphiAvg, part1.eta(), daughter.eta());
            histdetadpi_phi[i]->Fill(deta, dphiAvg, part1.phi(), daughter.phi());
          }
        }
        if (sameCharge) {
          if (atWhichRadiiToSelect == 1) {
            if (pow(dphiAvg, 2) / pow(deltaPhiMax, 2) + pow(deta, 2) / pow(deltaEtaMax, 2) < 1.) {
              pass = true;
            } else {
              if (Q3 == 999) {
                histdetadpi[i][1]->Fill(deta, dphiAvg);
                histdetadpi[i][3]->Fill(deta, dphi_AT_PV);
              } else if (Q3 < upperQ3LimitForPlotting) {
                histdetadpi[i][1]->Fill(deta, dphiAvg);
                histdetadpi[i][3]->Fill(deta, dphi_AT_PV);
              }
            }
          } else if (atWhichRadiiToSelect == 0) {
            if (pow(dphi_AT_PV, 2) / pow(deltaPhiMax, 2) + pow(deta, 2) / pow(deltaEtaMax, 2) < 1.) {
              pass = true;
            } else {
              if (Q3 == 999) {
                histdetadpi[i][1]->Fill(deta, dphiAvg);
                histdetadpi[i][3]->Fill(deta, dphi_AT_PV);
              } else if (Q3 < upperQ3LimitForPlotting) {
                histdetadpi[i][1]->Fill(deta, dphiAvg);
                histdetadpi[i][3]->Fill(deta, dphi_AT_PV);
              }
            }
          } else if (atWhichRadiiToSelect == 2) {
            if (pow(dphi_AT_SpecificRadii, 2) / pow(deltaPhiMax, 2) + pow(deta, 2) / pow(deltaEtaMax, 2) < 1.) {
              pass = true;
            } else {
              if (Q3 == 999) {
                histdetadpi[i][1]->Fill(deta, dphiAvg);
                histdetadpi[i][3]->Fill(deta, dphi_AT_PV);
              } else if (Q3 < upperQ3LimitForPlotting) {
                histdetadpi[i][1]->Fill(deta, dphiAvg);
                histdetadpi[i][3]->Fill(deta, dphi_AT_PV);
              }
            }
          }
        }
      }
      return pass;
    } else if constexpr (mPartOneType == o2::aod::femtodreamparticle::ParticleType::kTrack && mPartTwoType == o2::aod::femtodreamparticle::ParticleType::kCharmHadron) {
      // check if provided particles are in agreement with the class instantiation
      if (part2.candidateSelFlag() < o2::aod::fdhf::lcToPKPi) {
        LOG(fatal) << "FemtoDreamDetaDphiStar: passed arguments don't agree with FemtoDreamDetaDphiStar instantiation! Please provide Charm Hadron candidates.";
        return false;
      }

      bool pass = false;

      for (int i = 0; i < Nprongs; ++i) {
        double deta, dphiAvg, dphi_AT_PV, dphi_AT_SpecificRadii, daughterEta, daughterPhi;
        bool sameCharge = false;
        daughterEta = -999.;
        daughterPhi = -999.;

        switch (i) {
          case Prong0:
            daughterEta = part2.prong0Eta();
            daughterPhi = part2.prong0Phi();
            deta = part1.eta() - daughterEta;
            dphi_AT_PV = part1.phi() - daughterPhi;
            dphi_AT_SpecificRadii = PhiAtSpecificRadiiTPC(part1, radiiTPC) - PhiAtSpecificRadiiTPC<true, 0>(part2, radiiTPC);
            dphiAvg = AveragePhiStar<true>(part1, part2, 0, &sameCharge);
            // histdetadpi[0][0]->Fill(deta, dphiAvg);
            break;
          case Prong1:
            daughterEta = part2.prong1Eta();
            daughterPhi = part2.prong1Phi();
            deta = part1.eta() - daughterEta;
            dphi_AT_PV = part1.phi() - daughterPhi;
            dphi_AT_SpecificRadii = PhiAtSpecificRadiiTPC(part1, radiiTPC) - PhiAtSpecificRadiiTPC<true, 1>(part2, radiiTPC);
            dphiAvg = AveragePhiStar<true>(part1, part2, 1, &sameCharge);
            // histdetadpi[1][0]->Fill(deta, dphiAvg);
            break;
          case Prong2:
            daughterEta = part2.prong2Eta();
            daughterPhi = part2.prong2Phi();
            deta = part1.eta() - daughterEta;
            dphi_AT_PV = part1.phi() - daughterPhi;
            dphi_AT_SpecificRadii = PhiAtSpecificRadiiTPC(part1, radiiTPC) - PhiAtSpecificRadiiTPC<true, 2>(part2, radiiTPC);
            dphiAvg = AveragePhiStar<true>(part1, part2, 2, &sameCharge);
            // histdetadpi[2][0]->Fill(deta, dphiAvg);
            break;
        }
        if (Q3 == 999) {
          histdetadpi[i][0]->Fill(deta, dphiAvg);
          histdetadpi[i][2]->Fill(deta, dphi_AT_PV);
          if (fillQA) {
            histdetadpi_eta[i]->Fill(deta, dphiAvg, part1.eta(), daughterEta);
            histdetadpi_phi[i]->Fill(deta, dphiAvg, part1.phi(), daughterPhi);
          }
        } else if (Q3 < upperQ3LimitForPlotting) {
          histdetadpi[i][0]->Fill(deta, dphiAvg);
          histdetadpi[i][2]->Fill(deta, dphi_AT_PV);
          if (fillQA) {
            histdetadpi_eta[i]->Fill(deta, dphiAvg, part1.eta(), daughterEta);
            histdetadpi_phi[i]->Fill(deta, dphiAvg, part1.phi(), daughterPhi);
          }
        }

        if (atWhichRadiiToSelect == 1) {
          if (pow(dphiAvg, 2) / pow(deltaPhiMax, 2) + pow(deta, 2) / pow(deltaEtaMax, 2) < 1.) {
            pass = true;
          } else {
            if (Q3 == 999) {
              histdetadpi[i][1]->Fill(deta, dphiAvg);
              histdetadpi[i][3]->Fill(deta, dphi_AT_PV);
            } else if (Q3 < upperQ3LimitForPlotting) {
              histdetadpi[i][1]->Fill(deta, dphiAvg);
              histdetadpi[i][3]->Fill(deta, dphi_AT_PV);
            }
          }
        } else if (atWhichRadiiToSelect == 0) {
          if (pow(dphi_AT_PV, 2) / pow(deltaPhiMax, 2) + pow(deta, 2) / pow(deltaEtaMax, 2) < 1.) {
            pass = true;
          } else {
            if (Q3 == 999) {
              histdetadpi[i][1]->Fill(deta, dphiAvg);
              histdetadpi[i][3]->Fill(deta, dphi_AT_PV);
            } else if (Q3 < upperQ3LimitForPlotting) {
              histdetadpi[i][1]->Fill(deta, dphiAvg);
              histdetadpi[i][3]->Fill(deta, dphi_AT_PV);
            }
          }
        } else if (atWhichRadiiToSelect == 2) {
          if (pow(dphi_AT_SpecificRadii, 2) / pow(deltaPhiMax, 2) + pow(deta, 2) / pow(deltaEtaMax, 2) < 1.) {
            pass = true;
          } else {
            if (Q3 == 999) {
              histdetadpi[i][1]->Fill(deta, dphiAvg);
              histdetadpi[i][3]->Fill(deta, dphi_AT_PV);
            } else if (Q3 < upperQ3LimitForPlotting) {
              histdetadpi[i][1]->Fill(deta, dphiAvg);
              histdetadpi[i][3]->Fill(deta, dphi_AT_PV);
            }
          }
        }
      }

      return pass;

    } else if constexpr (mPartOneType == o2::aod::femtodreamparticle::ParticleType::kTrack && mPartTwoType == o2::aod::femtodreamparticle::ParticleType::kCascade) {
      /// Track-V0 combination
      // check if provided particles are in agreement with the class instantiation
      if (part1.partType() != o2::aod::femtodreamparticle::ParticleType::kTrack || part2.partType() != o2::aod::femtodreamparticle::ParticleType::kCascade) {
        LOG(fatal) << "FemtoDreamDetaDphiStar: passed arguments don't agree with FemtoDreamDetaDphiStar instantiation! Please provide kTrack,kV0 candidates.";
        return false;
      }

      bool pass = false;
      for (int i = 0; i < 3; i++) {
        int indexOfDaughter;
        if (isMixedEventLambda) {
          indexOfDaughter = part2.globalIndex() - 3 + i;
        } else {
          indexOfDaughter = part2.index() - 3 + i;
        }
        auto daughter = particles.begin() + indexOfDaughter;
        auto deta = part1.eta() - daughter.eta();
        auto dphi_AT_PV = part1.phi() - daughter.phi();
        auto dphi_AT_SpecificRadii = PhiAtSpecificRadiiTPC(part1, radiiTPC) - PhiAtSpecificRadiiTPC(daughter, radiiTPC);
        bool sameCharge = false;
        auto dphiAvg = AveragePhiStar(part1, *daughter, i, &sameCharge);
        if (Q3 == 999) {
          histdetadpi[i][0]->Fill(deta, dphiAvg);
          histdetadpi[i][2]->Fill(deta, dphi_AT_PV);
          if (fillQA) {
            histdetadpi_eta[i]->Fill(deta, dphiAvg, part1.eta(), daughter.eta());
            histdetadpi_phi[i]->Fill(deta, dphiAvg, part1.phi(), daughter.phi());
          }
        } else if (Q3 < upperQ3LimitForPlotting) {
          histdetadpi[i][0]->Fill(deta, dphiAvg);
          histdetadpi[i][2]->Fill(deta, dphi_AT_PV);
          if (fillQA) {
            histdetadpi_eta[i]->Fill(deta, dphiAvg, part1.eta(), daughter.eta());
            histdetadpi_phi[i]->Fill(deta, dphiAvg, part1.phi(), daughter.phi());
          }
        }
        if (sameCharge) {
          if (atWhichRadiiToSelect == 1) {
            if (pow(dphiAvg, 2) / pow(deltaPhiMax, 2) + pow(deta, 2) / pow(deltaEtaMax, 2) < 1.) {
              pass = true;
            } else {
              if (Q3 == 999) {
                histdetadpi[i][1]->Fill(deta, dphiAvg);
                histdetadpi[i][3]->Fill(deta, dphi_AT_PV);
              } else if (Q3 < upperQ3LimitForPlotting) {
                histdetadpi[i][1]->Fill(deta, dphiAvg);
                histdetadpi[i][3]->Fill(deta, dphi_AT_PV);
              }
            }

          } else if (atWhichRadiiToSelect == 0) {
            if (pow(dphi_AT_PV, 2) / pow(deltaPhiMax, 2) + pow(deta, 2) / pow(deltaEtaMax, 2) < 1.) {
              pass = true;
            } else {
              if (Q3 == 999) {
                histdetadpi[i][1]->Fill(deta, dphiAvg);
                histdetadpi[i][3]->Fill(deta, dphi_AT_PV);
              } else if (Q3 < upperQ3LimitForPlotting) {
                histdetadpi[i][1]->Fill(deta, dphiAvg);
                histdetadpi[i][3]->Fill(deta, dphi_AT_PV);
              }
            }
          } else if (atWhichRadiiToSelect == 2) {
            if (pow(dphi_AT_SpecificRadii, 2) / pow(deltaPhiMax, 2) + pow(deta, 2) / pow(deltaEtaMax, 2) < 1.) {
              pass = true;
            } else {
              if (Q3 == 999) {
                histdetadpi[i][1]->Fill(deta, dphiAvg);
                histdetadpi[i][3]->Fill(deta, dphi_AT_PV);
              } else if (Q3 < upperQ3LimitForPlotting) {
                histdetadpi[i][1]->Fill(deta, dphiAvg);
                histdetadpi[i][3]->Fill(deta, dphi_AT_PV);
              }
            }
          }
        }
      }
      return pass;
    } else {
      LOG(fatal) << "FemtoDreamPairCleaner: Combination of objects not defined - quitting!";
      return false;
    }
  }

 private:
  HistogramRegistry* mHistogramRegistry = nullptr;   ///< For main output
  HistogramRegistry* mHistogramRegistryQA = nullptr; ///< For QA output
  static constexpr std::string_view dirNames[4] = {"kTrack_kTrack/", "kTrack_kV0/", "kTrack_kCharmHadron/", "kTrack_kCascade/"};

  static constexpr std::string_view histNameSEorME[3] = {"_SEandME", "_SE", "_ME"};

  static constexpr std::string_view histNames[2][4] = {{"detadphidetadphi0Before_0", "detadphidetadphi0Before_1", "detadphidetadphi0Before_2", "detadphidetadphi0Before_3"},
                                                       {"detadphidetadphi0After_0", "detadphidetadphi0After_1", "detadphidetadphi0After_2", "detadphidetadphi0After_3"}};

  static constexpr std::string_view histNamesRadii[4][9] = {
    {"detadphidetadphi0Before_0_0", "detadphidetadphi0Before_0_1", "detadphidetadphi0Before_0_2",
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
     "detadphidetadphi0Before_3_6", "detadphidetadphi0Before_3_7", "detadphidetadphi0Before_3_8"}};

  static constexpr o2::aod::femtodreamparticle::ParticleType mPartOneType = partOne; ///< Type of particle 1
  static constexpr o2::aod::femtodreamparticle::ParticleType mPartTwoType = partTwo; ///< Type of particle 2

  static constexpr float tmpRadiiTPC[9] = {85., 105., 125., 145., 165., 185., 205., 225., 245.};

  static constexpr uint32_t kSignMinusMask = 1;
  static constexpr uint32_t kSignPlusMask = 1 << 1;
  static constexpr uint32_t kValue0 = 0;

  float deltaPhiMax;
  float deltaEtaMax;
  float magfield;
  bool plotForEveryRadii = false;
  bool isMixedEventLambda = false;
  float upperQ3LimitForPlotting = 8.;
  int atWhichRadiiToSelect = 1;
  float radiiTPC = 85.;
  bool fillQA = false;
  // a possible bug was found, but this must be tested on hyperloop with larger statistics
  // possiboility to run old code is turned on so a proper comparison of both code versions can be done
  bool runOldVersion = true;

  std::array<std::array<std::shared_ptr<TH2>, 4>, 3> histdetadpi{};
  std::array<std::array<std::shared_ptr<TH2>, 9>, 3> histdetadpiRadii{};
  std::array<std::shared_ptr<THnSparse>, 3> histdetadpi_eta{};
  std::array<std::shared_ptr<THnSparse>, 3> histdetadpi_phi{};

  ///  Calculate phi at all required radii stored in tmpRadiiTPC
  /// Magnetic field to be provided in Tesla
  template <typename T>
  int PhiAtRadiiTPC(const T& part, std::vector<float>& tmpVec)
  {

    float phi0 = part.phi();
    // Start: Get the charge from cutcontainer using masks
    int charge = 0.;
    if ((part.cut() & kSignMinusMask) == kValue0 && (part.cut() & kSignPlusMask) == kValue0) {
      charge = 0;
    } else if ((part.cut() & kSignPlusMask) == kSignPlusMask) {
      charge = 1;
    } else if ((part.cut() & kSignMinusMask) == kSignMinusMask) {
      charge = -1;
    } else {
      LOG(fatal) << "FemtoDreamDetaDphiStar: Charge bits are set wrong!";
    }
    // End: Get the charge from cutcontainer using masks
    float pt = part.pt();
    for (size_t i = 0; i < 9; i++) {
      if (runOldVersion) {
        tmpVec.push_back(phi0 - std::asin(0.3 * charge * 0.1 * magfield * tmpRadiiTPC[i] * 0.01 / (2. * pt)));
      }
      if (!runOldVersion) {
        auto arg = 0.3 * charge * magfield * tmpRadiiTPC[i] * 0.01 / (2. * pt);
        // for very low pT particles, this value goes outside of range -1 to 1 at at large tpc radius; asin fails
        if (std::fabs(arg) < 1) {
          tmpVec.push_back(phi0 - std::asin(0.3 * charge * magfield * tmpRadiiTPC[i] * 0.01 / (2. * pt)));
        } else {
          tmpVec.push_back(999);
        }
      }
    }
    return charge;
  }

  ///  Calculate phi at specific radii
  /// Magnetic field to be provided in Tesla
  template <bool isHF = false, int prong = 0, typename T>
  float PhiAtSpecificRadiiTPC(const T& part, float radii)
  {
    int charge = 0;
    float phi0, pt;
    if constexpr (isHF) {
      switch (prong) {
        case Prong0:
          charge = part.charge(); // charge calculation according to 3-prong decay, Lc^+ --> P^+ + K^- + pi^+
          phi0 = part.prong0Phi();
          pt = part.prong0Pt();
          break;
        case Prong1:
          charge = -part.charge();
          phi0 = part.prong1Phi();
          pt = part.prong1Pt();
          break;
        case Prong2:
          charge = part.charge();
          phi0 = part.prong2Phi();
          pt = part.prong2Pt();
          break;
        default:
          // Handle invalid prong value if necessary
          break;
      }
    } else {
      phi0 = part.phi();
      // Start: Get the charge from cutcontainer using masks
      if ((part.cut() & kSignMinusMask) == kValue0 && (part.cut() & kSignPlusMask) == kValue0) {
        charge = 0;
      } else if ((part.cut() & kSignPlusMask) == kSignPlusMask) {
        charge = 1;
      } else if ((part.cut() & kSignMinusMask) == kSignMinusMask) {
        charge = -1;
      } else {
        LOG(fatal) << "FemtoDreamDetaDphiStar: Charge bits are set wrong!";
      }
      pt = part.pt();
    }
    // End: Get the charge from cutcontainer using masks
    float phiAtRadii = 0;
    if (runOldVersion) {
      phiAtRadii = phi0 - std::asin(0.3 * charge * 0.1 * magfield * radii * 0.01 / (2. * pt));
    }
    if (!runOldVersion) {
      auto arg = 0.3 * charge * magfield * radii * 0.01 / (2. * pt);
      // for very low pT particles, this value goes outside of range -1 to 1 at at large tpc radius; asin fails
      if (std::fabs(arg) < 1) {
        phiAtRadii = phi0 - std::asin(0.3 * charge * magfield * radii * 0.01 / (2. * pt));
      } else {
        phiAtRadii = 999.;
      }
    }

    return phiAtRadii;
  }

  template <typename T>
  int PhiAtRadiiTPCForHF(const T& part, std::vector<float>& tmpVec, int prong)
  {
    int charge = 0;
    float pt = -999.;
    float phi0 = -999.;
    switch (prong) {
      case Prong0:
        pt = part.prong0Pt();
        phi0 = part.prong0Phi();
        charge = part.charge();
        break;
      case Prong1:
        pt = part.prong1Pt();
        phi0 = part.prong1Phi();
        charge = -part.charge();
        break;
      case Prong2:
        pt = part.prong2Pt();
        phi0 = part.prong2Phi();
        charge = part.charge();
        break;
      default:
        // Handle invalid prong value
        break;
    }
    for (size_t i = 0; i < 9; i++) {
      if (runOldVersion) {
        tmpVec.push_back(phi0 - std::asin(0.3 * charge * 0.1 * magfield * tmpRadiiTPC[i] * 0.01 / (2. * pt)));
      }
      if (!runOldVersion) {
        auto arg = 0.3 * charge * magfield * tmpRadiiTPC[i] * 0.01 / (2. * pt);
        // for very low pT particles, this value goes outside of range -1 to 1 at at large tpc radius; asin fails
        if (std::fabs(arg) < 1) {
          tmpVec.push_back(phi0 - std::asin(0.3 * charge * magfield * tmpRadiiTPC[i] * 0.01 / (2. * pt)));
        } else {
          tmpVec.push_back(999);
        }
      }
    }
    return charge;
  }

  ///  Calculate average phi
  template <bool isHF = false, typename T1, typename T2>
  float AveragePhiStar(const T1& part1, const T2& part2, int iHist, bool* sameCharge)
  {
    std::vector<float> tmpVec1;
    std::vector<float> tmpVec2;
    auto charge1 = PhiAtRadiiTPC(part1, tmpVec1);
    if constexpr (!isHF) {
      auto charge2 = PhiAtRadiiTPC(part2, tmpVec2);
      if (charge1 == charge2) {
        *sameCharge = true;
      }
    } else {
      PhiAtRadiiTPCForHF(part2, tmpVec2, iHist);
      *sameCharge = true; // always true as we checked the condition in the HF task
    }
    int num = tmpVec1.size();
    int meaningfulEntries = num;
    float dPhiAvg = 0;
    float dphi;
    for (int i = 0; i < num; i++) {
      if (tmpVec1.at(i) != 999 && tmpVec2.at(i) != 999) {
        dphi = tmpVec1.at(i) - tmpVec2.at(i);
      } else {
        dphi = 0;
        meaningfulEntries = meaningfulEntries - 1;
      }
      dphi = TVector2::Phi_mpi_pi(dphi);
      dPhiAvg += dphi;
      if (plotForEveryRadii) {
        histdetadpiRadii[iHist][i]->Fill(part1.eta() - part2.eta(), dphi);
      }
    }
    return dPhiAvg / static_cast<float>(meaningfulEntries);
  }
};

} /* namespace femtoDream */
} /* namespace o2::analysis */

#endif // PWGCF_FEMTODREAM_CORE_FEMTODREAMDETADPHISTAR_H_
