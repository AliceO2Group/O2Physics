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

/// \file FemtoWorldParticleHisto.h
/// \brief FemtoWorldParticleHisto - Histogram class for tracks, V0s and cascades
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Zuzanna Chochulska, WUT Warsaw, zchochul@cern.ch

#ifndef FEMTOWORLDPARTICLEHISTO_H_
#define FEMTOWORLDPARTICLEHISTO_H_

#include "PWGCF/FemtoWorld/DataModel/FemtoWorldDerived.h"
#include "Framework/HistogramRegistry.h"

using namespace o2::framework;
// using namespace o2::aod::o2::aod;

namespace o2::analysis::femtoWorld
{

/// \class FemtoWorldParticleHisto
/// \brief Class for histogramming particle properties
/// \tparam particleType Type of the particle (Track/V0/Cascade/...)
/// \tparam suffixType (optional) Takes care of the suffix for the folder name in case of analyses of pairs of the same kind (T-T, V-V, C-C)
template <o2::aod::femtoworldparticle::ParticleType particleType, int suffixType = 0>
class FemtoWorldParticleHisto
{
 public:
  /// Destructor
  virtual ~FemtoWorldParticleHisto() = default;

  /// Initalization of the QA histograms
  /// \param registry HistogramRegistry
  void init(HistogramRegistry* registry)
  {
    if (registry) {
      mHistogramRegistry = registry;
      /// The folder names are defined by the type of the object and the suffix (if applicable)
      std::string folderName = static_cast<std::string>(o2::aod::femtoworldparticle::ParticleTypeName[mParticleType]);
      folderName += static_cast<std::string>(mFolderSuffix[mFolderSuffixType]);

      /// Histograms of the kinematic properties
      mHistogramRegistry->add((folderName + "/hPt").c_str(), "; #it{p}_{T} (GeV/#it{c}); Entries", kTH1F, {{240, 0, 6}});
      mHistogramRegistry->add((folderName + "/hEta").c_str(), "; #eta; Entries", kTH1F, {{200, -1.5, 1.5}});
      mHistogramRegistry->add((folderName + "/hPhi").c_str(), "; #phi; Entries", kTH1F, {{200, 0, 2. * M_PI}});
      mHistogramRegistry->add((folderName + "/dEdxTPCVsMomentum").c_str(), "; #it{p} (GeV/#it{c}); dE/dx (keV/cm)", kTH2F, {{100, 0., 5.}, {250, 0., 500.}});
      mHistogramRegistry->add((folderName + "/TOFBetaVsMomentum").c_str(), "; #it{p} (GeV/#it{c}); TOF #beta", kTH2F, {{100, 0., 5.}, {250, 0.4, 1.1}});

      /// Particle-type specific histograms
      if constexpr (mParticleType == o2::aod::femtoworldparticle::ParticleType::kTrack) {
        /// Track histograms
        LOGF(info, "Tracks============================particle-histo");
        mHistogramRegistry->add((folderName + "/hDCAxy").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {{20, 0.5, 4.05}, {500, -5, 5}});
        mHistogramRegistry->add((folderName + "/hDCAz").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{z} (cm)", kTH2F, {{20, 0.5, 4.05}, {500, -5, 5}});
      } else if constexpr (mParticleType == o2::aod::femtoworldparticle::ParticleType::kV0) {
        /// V0 histograms
        mHistogramRegistry->add((folderName + "/hCPA").c_str(), "; #it{p}_{T} (GeV/#it{c}); cos#alpha", kTH2F, {{8, 0.3, 4.3}, {1000, 0.9, 1}});
      } else if constexpr (mParticleType == o2::aod::femtoworldparticle::ParticleType::kCascade) {
        /// Cascade histograms
      } else if constexpr (mParticleType == o2::aod::femtoworldparticle::ParticleType::kPhi) {
        /// Phi histograms
        LOGF(info, "Phi===================================particle-histo");
      } else {
        LOG(fatal) << "FemtoWorldParticleHisto: Histogramming for requested object not defined - quitting!";
      }
    }
  }

  /// Filling of the histograms
  /// \tparam T Data type of the particle
  /// \param part Particle
  template <typename T>
  void fillQA(T const& part)
  {
    if (mHistogramRegistry) {
      /// Histograms of the kinematic properties
      mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hPt"), part.pt());
      mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hEta"), part.eta());
      mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hPhi"), part.phi());
      mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/dEdxTPCVsMomentum"), part.p(), part.tpcSignal());
      mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/TOFBetaVsMomentum"), part.p(), part.beta());

      /// Particle-type specific histograms
      if constexpr (mParticleType == o2::aod::femtoworldparticle::ParticleType::kTrack) {
        /// Track histograms
        mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hDCAxy"), part.pt(), part.tempFitVar());
        mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hDCAz"), part.pt(), part.dcaZ());
      } else if constexpr (mParticleType == o2::aod::femtoworldparticle::ParticleType::kV0) {
        /// V0 histograms
        mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hCPA"),
                                 part.pt(), part.tempFitVar());
      } else if constexpr (mParticleType == o2::aod::femtoworldparticle::ParticleType::kCascade) {
        /// Cascade histograms
      } else if constexpr (mParticleType == o2::aod::femtoworldparticle::ParticleType::kPhi) {
        /// Phi histograms
      } else {
        LOG(fatal) << "FemtoWorldParticleHisto: Histogramming for requested object not defined - quitting!";
      }
    }
  }

 private:
  HistogramRegistry* mHistogramRegistry;                                                                      ///< For QA output
  static constexpr o2::aod::femtoworldparticle::ParticleType mParticleType = particleType;                    ///< Type of the particle under analysis
  static constexpr int mFolderSuffixType = suffixType;                                                        ///< Counter for the folder suffix specified below
  static constexpr std::string_view mFolderSuffix[5] = {"", "_one", "_one_rejected", "_two", "two_rejected"}; ///< Suffix for the folder name in case of analyses of pairs of the same kind (T-T, V-V, C-C)
};
} // namespace o2::analysis::femtoWorld

#endif /* FEMTOWORLDPARTICLEHISTO_H_ */
