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

/// \file FemtoUniverseDetaDphiStar.h
/// \brief FemtoUniverseDetaDphiStar - Checks particles for the close pair rejection.
/// \author Laura Serksnyte, TU München, laura.serksnyte@tum.de
/// \author Zuzanna Chochulska, WUT Warsaw, zuzanna.chochulska.stud@pw.edu.pl

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEDETADPHISTAR_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEDETADPHISTAR_H_

#include <memory>
#include <string>
#include <vector>
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::framework;

namespace o2::analysis
{
namespace femtoUniverse
{

/// \class FemtoUniverseDetaDphiStar
/// \brief Class to check particles for the close pair rejection.
/// \tparam partOne Type of particle 1 (Track/V0/Cascade/...)
/// \tparam partTwo Type of particle 2 (Track/V0/Cascade/...)
template <o2::aod::femtouniverseparticle::ParticleType partOne, o2::aod::femtouniverseparticle::ParticleType partTwo>
class FemtoUniverseDetaDphiStar
{
 public:
  /// Destructor
  virtual ~FemtoUniverseDetaDphiStar() = default;
  /// Initialization of the histograms and setting required values
  void init(HistogramRegistry* registry, HistogramRegistry* registryQA, float ldeltaPhiMax, float ldeltaEtaMax, bool lplotForEveryRadii)
  {
    deltaPhiMax = ldeltaPhiMax;
    deltaEtaMax = ldeltaEtaMax;
    plotForEveryRadii = lplotForEveryRadii;
    mHistogramRegistry = registry;
    mHistogramRegistryQA = registryQA;

    if constexpr (mPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && mPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kTrack) {
      std::string dirName = static_cast<std::string>(dirNames[0]);
      histdetadpi[0][0] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(histNames[0][0])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
      histdetadpi[0][1] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(histNames[1][0])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
      if (plotForEveryRadii) {
        for (int i = 0; i < 9; i++) {
          histdetadpiRadii[0][i] = mHistogramRegistryQA->add<TH2>((dirName + static_cast<std::string>(histNamesRadii[0][i])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        }
      }
    }
    if constexpr (mPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && mPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kV0) {
      for (int i = 0; i < 2; i++) {
        std::string dirName = static_cast<std::string>(dirNames[1]);
        histdetadpi[i][0] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(histNames[0][i])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpi[i][1] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(histNames[1][i])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        if (plotForEveryRadii) {
          for (int j = 0; j < 9; j++) {
            histdetadpiRadii[i][j] = mHistogramRegistryQA->add<TH2>((dirName + static_cast<std::string>(histNamesRadii[i][j])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
          }
        }
      }
    }
    if constexpr (mPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && mPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kPhi) {
      for (int i = 0; i < 2; i++) {
        std::string dirName = static_cast<std::string>(dirNames[2]);
        histdetadpi[i][0] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(histNames[0][i])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        histdetadpi[i][1] = mHistogramRegistry->add<TH2>((dirName + static_cast<std::string>(histNames[1][i])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
        if (plotForEveryRadii) {
          for (int j = 0; j < 9; j++) {
            histdetadpiRadii[i][j] = mHistogramRegistryQA->add<TH2>((dirName + static_cast<std::string>(histNamesRadii[i][j])).c_str(), "; #Delta #eta; #Delta #phi", kTH2F, {{100, -0.15, 0.15}, {100, -0.15, 0.15}});
          }
        }
      }
    }
  }
  ///  Check if pair is close or not
  template <typename Part, typename Parts>
  bool isClosePair(Part const& part1, Part const& part2, Parts const& particles, float lmagfield)
  {
    magfield = lmagfield;

    if constexpr (mPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && mPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kTrack) {
      /// Track-Track combination
      // check if provided particles are in agreement with the class instantiation
      if (part1.partType() != o2::aod::femtouniverseparticle::ParticleType::kTrack || part2.partType() != o2::aod::femtouniverseparticle::ParticleType::kTrack) {
        LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar instantiation! Please provide kTrack,kTrack candidates.";
        return false;
      }
      auto deta = part1.eta() - part2.eta();
      auto dphiAvg = AveragePhiStar(part1, part2, 0);
      histdetadpi[0][0]->Fill(deta, dphiAvg);
      if (pow(dphiAvg, 2) / pow(deltaPhiMax, 2) + pow(deta, 2) / pow(deltaEtaMax, 2) < 1.) {
        return true;
      } else {
        histdetadpi[0][1]->Fill(deta, dphiAvg);
        return false;
      }

    } else if constexpr (mPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && mPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kV0) {
      /// Track-V0 combination
      // check if provided particles are in agreement with the class instantiation
      if (part1.partType() != o2::aod::femtouniverseparticle::ParticleType::kTrack || part2.partType() != o2::aod::femtouniverseparticle::ParticleType::kV0) {
        LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar instantiation! Please provide kTrack,kV0 candidates.";
        return false;
      }

      bool pass = false;
      for (int i = 0; i < 2; i++) {
        auto indexOfDaughter = part2.index() - 2 + i;
        auto daughter = particles.begin() + indexOfDaughter;
        auto deta = part1.eta() - daughter.eta();
        auto dphiAvg = AveragePhiStar(part1, *daughter, i);
        histdetadpi[i][0]->Fill(deta, dphiAvg);
        if (pow(dphiAvg, 2) / pow(deltaPhiMax, 2) + pow(deta, 2) / pow(deltaEtaMax, 2) < 1.) {
          pass = true;
        } else {
          histdetadpi[i][1]->Fill(deta, dphiAvg);
        }
      }
      return pass;
    } else if constexpr (mPartOneType == o2::aod::femtouniverseparticle::ParticleType::kTrack && mPartTwoType == o2::aod::femtouniverseparticle::ParticleType::kPhi) {
      /// Track-Phi combination
      // check if provided particles are in agreement with the class instantiation
      if (part1.partType() != o2::aod::femtouniverseparticle::ParticleType::kTrack || part2.partType() != o2::aod::femtouniverseparticle::ParticleType::kPhi) {
        LOG(fatal) << "FemtoUniverseDetaDphiStar: passed arguments don't agree with FemtoUniverseDetaDphiStar instantiation! Please provide kTrack,kPhi candidates.";
        return false;
      }

      bool pass = false;
      for (int i = 0; i < 2; i++) {
        auto indexOfDaughter = part2.index() - 2 + i;
        auto daughter = particles.begin() + indexOfDaughter;
        auto deta = part1.eta() - daughter.eta();
        auto dphiAvg = AveragePhiStar(part1, *daughter, i);
        histdetadpi[i][0]->Fill(deta, dphiAvg);
        if (pow(dphiAvg, 2) / pow(deltaPhiMax, 2) + pow(deta, 2) / pow(deltaEtaMax, 2) < 1.) {
          pass = true;
        } else {
          histdetadpi[i][1]->Fill(deta, dphiAvg);
        }
      }
      return pass;
    } else {
      LOG(fatal) << "FemtoUniversePairCleaner: Combination of objects not defined - quitting!";
      return false;
    }
  }

 private:
  HistogramRegistry* mHistogramRegistry = nullptr;   ///< For main output
  HistogramRegistry* mHistogramRegistryQA = nullptr; ///< For QA output
  static constexpr std::string_view dirNames[3] = {"kTrack_kTrack/", "kTrack_kV0/", "kTrack_kPhi/"};

  static constexpr std::string_view histNames[2][2] = {{"detadphidetadphi0Before_0", "detadphidetadphi0Before_1"},
                                                       {"detadphidetadphi0After_0", "detadphidetadphi0After_1"}};
  static constexpr std::string_view histNamesRadii[2][9] = {{"detadphidetadphi0Before_0_0", "detadphidetadphi0Before_0_1", "detadphidetadphi0Before_0_2",
                                                             "detadphidetadphi0Before_0_3", "detadphidetadphi0Before_0_4", "detadphidetadphi0Before_0_5",
                                                             "detadphidetadphi0Before_0_6", "detadphidetadphi0Before_0_7", "detadphidetadphi0Before_0_8"},
                                                            {"detadphidetadphi0Before_1_0", "detadphidetadphi0Before_1_1", "detadphidetadphi0Before_1_2",
                                                             "detadphidetadphi0Before_1_3", "detadphidetadphi0Before_1_4", "detadphidetadphi0Before_1_5",
                                                             "detadphidetadphi0Before_1_6", "detadphidetadphi0Before_1_7", "detadphidetadphi0Before_1_8"}};

  static constexpr o2::aod::femtouniverseparticle::ParticleType mPartOneType = partOne; ///< Type of particle 1
  static constexpr o2::aod::femtouniverseparticle::ParticleType mPartTwoType = partTwo; ///< Type of particle 2

  static constexpr float tmpRadiiTPC[9] = {85., 105., 125., 145., 165., 185., 205., 225., 245.};

  static constexpr uint32_t kSignMinusMask = 1;
  static constexpr uint32_t kSignPlusMask = 1 << 1;
  static constexpr uint32_t kValue0 = 0;

  float deltaPhiMax;
  float deltaEtaMax;
  float magfield;
  bool plotForEveryRadii = false;

  std::array<std::array<std::shared_ptr<TH2>, 2>, 2> histdetadpi{};
  std::array<std::array<std::shared_ptr<TH2>, 9>, 2> histdetadpiRadii{};

  ///  Calculate phi at all required radii stored in tmpRadiiTPC
  /// Magnetic field to be provided in Tesla
  template <typename T>
  void PhiAtRadiiTPC(const T& part, std::vector<float>& tmpVec)
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
      tmpVec.push_back(phi0 - std::asin(0.3 * charge * 0.1 * magfield * tmpRadiiTPC[i] * 0.01 / (2. * pt)));
    }
  }

  ///  Calculate average phi
  template <typename T1, typename T2>
  float AveragePhiStar(const T1& part1, const T2& part2, int iHist)
  {
    std::vector<float> tmpVec1;
    std::vector<float> tmpVec2;
    PhiAtRadiiTPC(part1, tmpVec1);
    PhiAtRadiiTPC(part2, tmpVec2);
    int num = tmpVec1.size();
    float dPhiAvg = 0;
    for (int i = 0; i < num; i++) {
      float dphi = tmpVec1.at(i) - tmpVec2.at(i);
      dphi = TVector2::Phi_mpi_pi(dphi);
      dPhiAvg += dphi;
      if (plotForEveryRadii) {
        histdetadpiRadii[iHist][i]->Fill(part1.eta() - part2.eta(), dphi);
      }
    }
    return dPhiAvg / num;
  }
};

} /* namespace femtoUniverse */
} /* namespace o2::analysis */

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEDETADPHISTAR_H_
