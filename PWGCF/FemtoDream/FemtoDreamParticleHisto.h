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

/// \file FemtoDreamParticleHisto.h
/// \brief FemtoDreamParticleHisto - Histogram class for tracks, V0s and cascades
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#ifndef PWGCF_FEMTODREAM_FEMTODREAMPARTICLEHISTO_H_
#define PWGCF_FEMTODREAM_FEMTODREAMPARTICLEHISTO_H_

#include <string>
#include "PWGCF/DataModel/FemtoDerived.h"
#include "Framework/HistogramRegistry.h"

using namespace o2::framework;

namespace o2::analysis::femtoDream
{

/// \class FemtoDreamParticleHisto
/// \brief Class for histogramming particle properties
/// \tparam particleType Type of the particle (Track/V0/Cascade/...)
/// \tparam suffixType (optional) Takes care of the suffix for the folder name in case of analyses of pairs of the same kind (T-T, V-V, C-C)
template <o2::aod::femtodreamparticle::ParticleType particleType, int suffixType = 0>
class FemtoDreamParticleHisto
{
 public:
  /// Destructor
  virtual ~FemtoDreamParticleHisto() = default;

  /// Initialization of the QA histograms
  /// \param registry HistogramRegistry
  template <typename T>
  void init(HistogramRegistry* registry, T& tempFitVarpTBins, T& tempFitVarBins)
  {
    if (registry) {
      mHistogramRegistry = registry;
      /// The folder names are defined by the type of the object and the suffix (if applicable)
      std::string folderName = static_cast<std::string>(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]);
      folderName += static_cast<std::string>(mFolderSuffix[mFolderSuffixType]);

      /// Histograms of the kinematic properties
      mHistogramRegistry->add((folderName + "/hPt").c_str(), "; #it{p}_{T} (GeV/#it{c}); Entries", kTH1F, {{240, 0, 6}});
      mHistogramRegistry->add((folderName + "/hEta").c_str(), "; #eta; Entries", kTH1F, {{200, -1.5, 1.5}});
      mHistogramRegistry->add((folderName + "/hPhi").c_str(), "; #phi; Entries", kTH1F, {{200, 0, 2. * M_PI}});

      /// particle specific histogramms for the TempFitVar column in FemtoDreamParticles

      std::string tempFitVarName;
      std::string tempFitVarAxisTitle;
      if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kTrack) {
        /// Track histograms
        tempFitVarName = "/hDCAxy";
        tempFitVarAxisTitle = "DCA_{xy} (cm)";
      } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kV0) {
        /// V0 histograms
        tempFitVarName = "/hCPA";
        tempFitVarAxisTitle = "cos#alpha";
      } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascade) {
        /// Cascade histograms
        tempFitVarName = "/hCPA";
        tempFitVarAxisTitle = "cos#alpha";
      } else {
        LOG(fatal) << "FemtoDreamParticleHisto: Histogramming for requested object not defined - quitting!";
      }

      framework::AxisSpec tempFitVarpTAxis = {tempFitVarpTBins, "#it{p}_{T} (GeV/#it{c})"}; // the pT binning may vary
      framework::AxisSpec tempFitVarAxis = {tempFitVarBins, tempFitVarAxisTitle};
      mHistogramRegistry->add((folderName + tempFitVarName).c_str(), ("; #it{p}_{T} (GeV/#it{c}); " + tempFitVarAxisTitle).c_str(), kTH2F, {{tempFitVarpTAxis}, {tempFitVarAxis}});
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
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hPt"), part.pt());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hEta"), part.eta());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hPhi"), part.phi());

      /// Particle-type specific histograms
      if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kTrack) {
        /// Track histograms
        mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hDCAxy"),
                                 part.pt(), part.tempFitVar());
      } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kV0) {
        /// V0 histograms
        mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("/hCPA"),
                                 part.pt(), part.tempFitVar());
      } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascade) {
        /// Cascade histograms
      } else {
        LOG(fatal) << "FemtoDreamParticleHisto: Histogramming for requested object not defined - quitting!";
      }
    }
  }

 private:
  HistogramRegistry* mHistogramRegistry;                                                   ///< For QA output
  static constexpr o2::aod::femtodreamparticle::ParticleType mParticleType = particleType; ///< Type of the particle under analysis
  static constexpr int mFolderSuffixType = suffixType;                                     ///< Counter for the folder suffix specified below
  static constexpr std::string_view mFolderSuffix[3] = {"", "_one", "_two"};               ///< Suffix for the folder name in case of analyses of pairs of the same kind (T-T, V-V, C-C)
};
} // namespace o2::analysis::femtoDream

#endif // PWGCF_FEMTODREAM_FEMTODREAMPARTICLEHISTO_H_
