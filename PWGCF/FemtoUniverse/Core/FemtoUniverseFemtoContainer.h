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

/// \file FemtoUniverseFemtoContainer.h
/// \brief Definition of the FemtoUniverseFemtoContainer
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Valentina Mantovani Sarti, valentina.mantovani-sarti@tum.de
/// \author Georgios Mantzaridis, TU München, georgios.mantzaridis@tum.de
/// \author Anton Riedel, TU München, anton.riedel@tum.de
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEFEMTOCONTAINER_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEFEMTOCONTAINER_H_

#include "PWGCF/FemtoUniverse/Core/FemtoUniverseMath.h"

#include "Framework/HistogramRegistry.h"
#include <Framework/Logger.h>

#include "Math/Vector4D.h"
#include "TDatabasePDG.h"
#include "TMath.h"

#include <string>
#include <vector>

using namespace o2::framework;

namespace o2::analysis::femto_universe
{

namespace femto_universe_femto_container
{
/// Femtoscopic observable to be computed
enum Observable { kstar ///< kstar
};

/// Type of the event processind
enum EventType { same, ///< Pair from same event
                 mixed ///< Pair from mixed event
};
}; // namespace femto_universe_femto_container

/// \class FemtoUniverseFemtoContainer
/// \brief Container for all histogramming related to the correlation function. The two
/// particles of the pair are passed here, and the correlation function and QA histograms
/// are filled according to the specified observable
/// \tparam eventType Type of the event (same/mixed)
/// \tparam obs Observable to be computed (k*/Q_inv/...)
template <femto_universe_femto_container::EventType eventType, femto_universe_femto_container::Observable obs>
class FemtoUniverseFemtoContainer
{
 public:
  /// Destructor
  virtual ~FemtoUniverseFemtoContainer() = default;

  /// Initializes histograms for the task
  /// Called by init both in case of reconstructed data/ Monte Carlo, and for Monte Carlo Truth
  /// \tparam T type of the axis Object
  /// \param folderName Name of the directory in the output file (no suffix for reconstructed data/ Monte Carlo; "_MC" for Monte Carlo Truth)
  /// \param femtoObs Title of the femto observable axis
  /// \param femtoObsAxis axis object for the femto observable axis
  /// \param multAxis axis object for the multiplicity axis
  /// \param kTAxis axis object for the kT axis
  /// \param mTAxis axis object for the mT axis
  template <typename T>
  void initBase(std::string folderName, std::string femtoObs, T femtoObsAxis, T multAxis, T kTAxis, T mTAxis, T multAxis3D, T mTAxis3D, bool use3dplots)
  {
    kHistogramRegistry->add((folderName + "/relPairDist").c_str(), ("; " + femtoObs + "; Entries").c_str(), kTH1F, {femtoObsAxis});
    kHistogramRegistry->add((folderName + "/relPairkT").c_str(), "; #it{k}_{T} (GeV/#it{c}); Entries", kTH1F, {kTAxis});
    kHistogramRegistry->add((folderName + "/relPairkstarkT").c_str(), ("; " + femtoObs + "; #it{k}_{T} (GeV/#it{c})").c_str(), kTH2F, {femtoObsAxis, kTAxis});
    kHistogramRegistry->add((folderName + "/relPairkstarmT").c_str(), ("; " + femtoObs + "; #it{m}_{T} (GeV/#it{c}^{2})").c_str(), kTH2F, {femtoObsAxis, mTAxis});
    kHistogramRegistry->add((folderName + "/relPairkstarMult").c_str(), ("; " + femtoObs + "; Multiplicity").c_str(), kTH2F, {femtoObsAxis, multAxis});
    kHistogramRegistry->add((folderName + "/kstarPtPart1").c_str(), ("; " + femtoObs + "; #it{p} _{T} Particle 1 (GeV/#it{c})").c_str(), kTH2F, {femtoObsAxis, {375, 0., 7.5}});
    kHistogramRegistry->add((folderName + "/kstarPtPart2").c_str(), ("; " + femtoObs + "; #it{p} _{T} Particle 2 (GeV/#it{c})").c_str(), kTH2F, {femtoObsAxis, {375, 0., 7.5}});
    kHistogramRegistry->add((folderName + "/MultPtPart1").c_str(), "; #it{p} _{T} Particle 1 (GeV/#it{c}); Multiplicity", kTH2F, {{375, 0., 7.5}, multAxis});
    kHistogramRegistry->add((folderName + "/MultPtPart2").c_str(), "; #it{p} _{T} Particle 2 (GeV/#it{c}); Multiplicity", kTH2F, {{375, 0., 7.5}, multAxis});
    kHistogramRegistry->add((folderName + "/PtPart1PtPart2").c_str(), "; #it{p} _{T} Particle 1 (GeV/#it{c}); #it{p} _{T} Particle 2 (GeV/#it{c})", kTH2F, {{375, 0., 7.5}, {375, 0., 7.5}});
    if (use3dplots) {
      kHistogramRegistry->add((folderName + "/relPairkstarmTMult").c_str(), ("; " + femtoObs + "; #it{m}_{T} (GeV/#it{c}^{2}); Multiplicity").c_str(), kTH3F, {femtoObsAxis, mTAxis3D, multAxis3D});
    }
  }

  /// Initializes specialized Monte Carlo truth histograms for the task
  /// internal function called by init only in case of Monte Carlo truth
  /// \tparam T type of the xxis Object
  /// \param folderName Name of the directory in the output file (no suffix for reconstructed data/ Monte Carlo; "_MC" for Monte Carlo Truth)
  /// \param femtoObsAxis axis object for the femto observable axis
  template <typename T>
  void initMC(std::string folderName, std::string femtoObs, T femtoObsAxis, T multAxis, T mTAxis)
  {
    kHistogramRegistry->add((folderName + "/relPairDist_ReconNoFake").c_str(), ("; " + femtoObs + "; Entries").c_str(), kTH1F, {femtoObsAxis});
    kHistogramRegistry->add((folderName + "/relPairkstarmT_ReconNoFake").c_str(), ("; " + femtoObs + "; #it{m}_{T} (GeV/#it{c}^{2})").c_str(), kTH2F, {femtoObsAxis, mTAxis});
    kHistogramRegistry->add((folderName + "/relPairkstarMult_ReconNoFake").c_str(), ("; " + femtoObs + "; Multiplicity").c_str(), kTH2F, {femtoObsAxis, multAxis});
    kHistogramRegistry->add((folderName + "/hNoMCtruthPairsCounter").c_str(), "; Counter; Entries", kTH1I, {{1, 0, 1}});
    kHistogramRegistry->add((folderName + "/hFakePairsCounter").c_str(), "; Counter; Entries", kTH1I, {{1, 0, 1}});
    kHistogramRegistry->add((folderName + "/kstar_resolution").c_str(), "; #it{k} _{T} reconstructed (GeV/#it{c}); #it{k} _{T} truth (GeV/#it{c})", kTH2F, {femtoObsAxis, femtoObsAxis});
  }

  /// Templated function to initialize the histograms for the task
  /// Always calls initBase to initialize the histograms for data/ Monte Carlo reconstructed
  /// In case of Monte Carlo, calls initBase again for Monte Carlo truth and the specialized function initMC for additional histogramms
  /// \tparam T type of the configurable for the axis configuration
  /// \param registry Histogram registry to be passed
  /// \param kstarBins k* binning for the histograms
  /// \param multBins multiplicity binning for the histograms
  /// \param kTBins kT binning for the histograms
  /// \param mTBins mT binning for the histograms
  /// \param isMC add Monte Carlo truth histograms to the output file
  template <typename T>
  void init(HistogramRegistry* registry, T& kstarBins, T& multBins, T& kTBins, T& mTBins, T& multBins3D, T& mTBins3D, bool isMC, bool use3dplots)
  {
    kHistogramRegistry = registry;
    std::string femtoObs;
    if constexpr (kFemtoObs == femto_universe_femto_container::Observable::kstar) {
      femtoObs = "#it{k*} (GeV/#it{c})";
    }
    std::vector<double> tmpVecMult = multBins;
    framework::AxisSpec multAxis = {tmpVecMult, "Multiplicity"};
    framework::AxisSpec femtoObsAxis = {kstarBins, femtoObs.c_str()};
    framework::AxisSpec kTAxis = {kTBins, "#it{k}_{T} (GeV/#it{c})"};
    framework::AxisSpec mTAxis = {mTBins, "#it{m}_{T} (GeV/#it{c}^{2})"};

    framework::AxisSpec multAxis3D = {multBins3D, "Multiplicity"};
    framework::AxisSpec mTAxis3D = {mTBins3D, "#it{m}_{T} (GeV/#it{c})"};

    std::string folderName = static_cast<std::string>(kFolderSuffix[kEventType]) + static_cast<std::string>(o2::aod::femtouniverse_mc_particle::MCTypeName[o2::aod::femtouniverse_mc_particle::MCType::kRecon]);

    initBase(folderName, femtoObs, femtoObsAxis, multAxis, kTAxis, mTAxis, multAxis3D, mTAxis3D, use3dplots);
    if (isMC) {
      folderName = static_cast<std::string>(kFolderSuffix[kEventType]) + static_cast<std::string>(o2::aod::femtouniverse_mc_particle::MCTypeName[o2::aod::femtouniverse_mc_particle::MCType::kTruth]);
      initBase(folderName, femtoObs, femtoObsAxis, multAxis, kTAxis, mTAxis, multAxis3D, mTAxis3D, use3dplots);
      initMC(folderName, femtoObs, femtoObsAxis, multAxis, mTAxis);
    }
  }

  /// Set the PDG codes of the two particles involved
  /// \param pdg1 PDG code of particle one
  /// \param pdg2 PDG code of particle two
  void setPDGCodes(const int pdg1, const int pdg2)
  {
    kMassOne = TDatabasePDG::Instance()->GetParticle(pdg1)->Mass();
    kMassTwo = TDatabasePDG::Instance()->GetParticle(pdg2)->Mass();
    kPDGOne = pdg1;
    kPDGTwo = pdg2;
  }

  /// Pass a pair to the container and compute all the relevant observables
  /// Called by setPair both in case of data/ and Monte Carlo reconstructed and for Monte Carlo truth
  /// \tparam T type of the femtouniverseparticle
  /// \param part1 Particle one
  /// \param part2 Particle two
  /// \param mult Multiplicity of the event
  template <o2::aod::femtouniverse_mc_particle::MCType mc, typename T>
  void setPairBase(const float femtoObs, const float mT, T const& part1, T const& part2, const int mult, bool use3dplots)
  {
    const float kT = FemtoUniverseMath::getkT(part1, kMassOne, part2, kMassTwo);

    kHistogramRegistry->fill(HIST(kFolderSuffix[kEventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/relPairDist"), femtoObs);
    kHistogramRegistry->fill(HIST(kFolderSuffix[kEventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/relPairkT"), kT);
    kHistogramRegistry->fill(HIST(kFolderSuffix[kEventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/relPairkstarkT"), femtoObs, kT);
    kHistogramRegistry->fill(HIST(kFolderSuffix[kEventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/relPairkstarmT"), femtoObs, mT);
    kHistogramRegistry->fill(HIST(kFolderSuffix[kEventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/relPairkstarMult"), femtoObs, mult);
    kHistogramRegistry->fill(HIST(kFolderSuffix[kEventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/kstarPtPart1"), femtoObs, part1.pt());
    kHistogramRegistry->fill(HIST(kFolderSuffix[kEventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/kstarPtPart2"), femtoObs, part2.pt());
    kHistogramRegistry->fill(HIST(kFolderSuffix[kEventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/MultPtPart1"), part1.pt(), mult);
    kHistogramRegistry->fill(HIST(kFolderSuffix[kEventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/MultPtPart2"), part2.pt(), mult);
    kHistogramRegistry->fill(HIST(kFolderSuffix[kEventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/PtPart1PtPart2"), part1.pt(), part2.pt());
    if (use3dplots) {
      kHistogramRegistry->fill(HIST(kFolderSuffix[kEventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/relPairkstarmTMult"), femtoObs, mT, mult);
    }
  }

  /// Called by setPair only in case of Monte Carlo truth
  /// Fills MC truth specific histogramms:
  /// - kstar distribution plots with RECONSTRUCTED information but ONLY for non-fake candidates; needed for purity calculations of tracks
  /// - kstar resolution matrix
  /// Note: Standard histogramms with MC truth information are filled with the setPairBase function
  /// \param part1 Particle one
  /// \param part2 Particle two
  /// \param mult Multiplicity of the event
  void setPairMC(const float femtoObsMC, const float femtoObs, const float mT, const int mult)
  {
    if (kHistogramRegistry) {
      // Fill the kstar distributions with the reconstructed information but only for particles with the right PDG code
      kHistogramRegistry->fill(HIST(kFolderSuffix[kEventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[o2::aod::femtouniverse_mc_particle::MCType::kTruth]) + HIST("/relPairDist_ReconNoFake"), femtoObs);
      kHistogramRegistry->fill(HIST(kFolderSuffix[kEventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[o2::aod::femtouniverse_mc_particle::MCType::kTruth]) + HIST("/relPairkstarmT_ReconNoFake"), femtoObs, mT);
      kHistogramRegistry->fill(HIST(kFolderSuffix[kEventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[o2::aod::femtouniverse_mc_particle::MCType::kTruth]) + HIST("/relPairkstarMult_ReconNoFake"), femtoObs, mult);

      kHistogramRegistry->fill(HIST(kFolderSuffix[kEventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[o2::aod::femtouniverse_mc_particle::MCType::kTruth]) + HIST("/kstar_resolution"), femtoObsMC, femtoObs);
    }
  }

  /// Templated function to handle data/ Monte Carlo reconstructed and Monte Carlo truth
  /// Always calls setPairBase to compute the observables with reconstructed data
  /// In case of Monte Carlo, calls setPairBase with MC info and specialized function setPairMC for additional histogramms
  /// \tparam T type of the femtouniverseparticle
  /// \param part1 Particle one
  /// \param part2 Particle two
  /// \param mult Multiplicity of the event
  template <bool isMC, typename T>
  void setPair(T const& part1, T const& part2, const int mult, bool use3dplots)
  {
    float femtoObs, femtoObsMC;
    // Calculate femto observable and the mT with reconstructed information
    if constexpr (kFemtoObs == femto_universe_femto_container::Observable::kstar) {
      femtoObs = FemtoUniverseMath::getkstar(part1, kMassOne, part2, kMassTwo);
    }
    const float mT = FemtoUniverseMath::getmT(part1, kMassOne, part2, kMassTwo);

    if (kHistogramRegistry) {
      setPairBase<o2::aod::femtouniverse_mc_particle::MCType::kRecon>(femtoObs, mT, part1, part2, mult, use3dplots);

      if constexpr (isMC) {
        if (part1.has_fdMCParticle() && part2.has_fdMCParticle()) {
          // calculate the femto observable and the mT with MC truth information
          if constexpr (kFemtoObs == femto_universe_femto_container::Observable::kstar) {
            femtoObsMC = FemtoUniverseMath::getkstar(part1.fdMCParticle(), kMassOne, part2.fdMCParticle(), kMassTwo);
          }
          const float mTMC = FemtoUniverseMath::getmT(part1.fdMCParticle(), kMassOne, part2.fdMCParticle(), kMassTwo);

          if (std::abs(part1.fdMCParticle().pdgMCTruth()) == std::abs(kPDGOne) && std::abs(part2.fdMCParticle().pdgMCTruth()) == std::abs(kPDGTwo)) { // Note: all pair-histogramms are filled with MC truth information ONLY in case of non-fake candidates
            setPairBase<o2::aod::femtouniverse_mc_particle::MCType::kTruth>(femtoObsMC, mTMC, part1.fdMCParticle(), part2.fdMCParticle(), mult, use3dplots);
            setPairMC(femtoObsMC, femtoObs, mT, mult);
          } else {
            kHistogramRegistry->fill(HIST(kFolderSuffix[kEventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[o2::aod::femtouniverse_mc_particle::MCType::kTruth]) + HIST("/hFakePairsCounter"), 0);
          }

        } else {
          kHistogramRegistry->fill(HIST(kFolderSuffix[kEventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[o2::aod::femtouniverse_mc_particle::MCType::kTruth]) + HIST("/hNoMCtruthPairsCounter"), 0);
        }
      }
    }
  }

 protected:
  HistogramRegistry* kHistogramRegistry = nullptr;                                  ///< For QA output
  static constexpr std::string_view kFolderSuffix[2] = {"SameEvent", "MixedEvent"}; ///< Folder naming for the output according to kEventType
  static constexpr femto_universe_femto_container::Observable kFemtoObs = obs;      ///< Femtoscopic observable to be computed (according to femto_universe_femto_container::Observable)
  static constexpr int kEventType = eventType;                                      ///< Type of the event (same/mixed, according to femto_universe_femto_container::EventType)
  float kMassOne = 0.f;                                                             ///< PDG mass of particle 1
  float kMassTwo = 0.f;                                                             ///< PDG mass of particle 2
  int kPDGOne = 0;                                                                  ///< PDG code of particle 1
  int kPDGTwo = 0;                                                                  ///< PDG code of particle 2
};

} // namespace o2::analysis::femto_universe

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEFEMTOCONTAINER_H_
