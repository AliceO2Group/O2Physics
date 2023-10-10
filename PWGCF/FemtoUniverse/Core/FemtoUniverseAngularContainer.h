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

/// \file FemtoUniverseAngularContainer.h
/// \brief Definition of the FemtoUniverseAngularContainer
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Valentina Mantovani Sarti, valentina.mantovani-sarti@tum.de
/// \author Georgios Mantzaridis, TU München, georgios.mantzaridis@tum.de
/// \author Anton Riedel, TU München, anton.riedel@tum.de
/// \author Zuzanna Chochulska, WUT Warsaw, zuzanna.chochulska.stud@pw.edu.pl

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEANGULARCONTAINER_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEANGULARCONTAINER_H_

#include <fairlogger/Logger.h>
#include <vector>
#include <string>

#include "Framework/HistogramRegistry.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseMath.h"

#include "Math/Vector4D.h"
#include "TMath.h"
#include "TDatabasePDG.h"

using namespace o2::framework;

namespace o2::analysis::femtoUniverse
{

namespace femtoUniverseAngularContainer
{
/// Femtoscopic observable to be computed
enum Observable { kstar ///< kstar
};

/// Type of the event processind
enum EventType { same, ///< Pair from same event
                 mixed ///< Pair from mixed event
};
}; // namespace femtoUniverseAngularContainer

/// \class FemtoUniverseAngularContainer
/// \brief Container for all histogramming related to the correlation function. The two
/// particles of the pair are passed here, and the correlation function and QA histograms
/// are filled according to the specified observable
/// \tparam eventType Type of the event (same/mixed)
/// \tparam obs Observable to be computed (k*/Q_inv/...)
template <femtoUniverseAngularContainer::EventType eventType, femtoUniverseAngularContainer::Observable obs>
class FemtoUniverseAngularContainer
{
 public:
  /// Destructor
  virtual ~FemtoUniverseAngularContainer() = default;

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
  void init_base(std::string folderName, std::string femtoObs, T femtoObsAxis, T multAxis, T kTAxis, T mTAxis, T multAxis3D, T mTAxis3D, T etaAxis, T phiAxis, bool use3dplots)
  {
    mHistogramRegistry->add((folderName + "/DeltaEtaDeltaPhi").c_str(), ";  #Delta#varphi (rad); #Delta#eta", kTH2F, {phiAxis, etaAxis});
    if (use3dplots) {
      // use 3d plots
    }
  }

  /// Initializes specialized Monte Carlo truth histograms for the task
  /// internal function called by init only in case of Monte Carlo truth
  /// \tparam T type of the xxis Object
  /// \param folderName Name of the directory in the output file (no suffix for reconstructed data/ Monte Carlo; "_MC" for Monte Carlo Truth)
  /// \param femtoObsAxis axis object for the femto observable axis
  template <typename T>
  void init_MC(std::string folderName, std::string femtoObs, T femtoObsAxis, T multAxis, T mTAxis)
  {
  }

  /// Templated function to initialize the histograms for the task
  /// Always calls init_base to initialize the histograms for data/ Monte Carlo reconstructed
  /// In case of Monte Carlo, calls init_base again for Monte Carlo truth and the specialized function init_MC for additional histogramms
  /// \tparam T type of the configurable for the axis configuration
  /// \param registry Histogram registry to be passed
  /// \param kstarBins k* binning for the histograms
  /// \param multBins multiplicity binning for the histograms
  /// \param kTBins kT binning for the histograms
  /// \param mTBins mT binning for the histograms
  /// \param etaBins eta binning for the histograms
  /// \param phiBins phi binning for the histograms
  /// \param isMC add Monte Carlo truth histograms to the output file
  template <typename T, typename P>
  void init(HistogramRegistry* registry, T& kstarBins, T& multBins, T& kTBins, T& mTBins, T& multBins3D, T& mTBins3D, P& etaBins, P& phiBins, bool isMC, bool use3dplots)
  {
    mHistogramRegistry = registry;
    std::string femtoObs;
    if constexpr (mFemtoObs == femtoUniverseAngularContainer::Observable::kstar) {
      femtoObs = "#it{k*} (GeV/#it{c})";
    }

    std::vector<double> tmpVecMult = multBins;
    framework::AxisSpec multAxis = {tmpVecMult, "Multiplicity"};
    framework::AxisSpec femtoObsAxis = {kstarBins, femtoObs.c_str()};
    framework::AxisSpec kTAxis = {kTBins, "#it{k}_{T} (GeV/#it{c})"};
    framework::AxisSpec mTAxis = {mTBins, "#it{m}_{T} (GeV/#it{c}^{2})"};

    framework::AxisSpec multAxis3D = {multBins3D, "Multiplicity"};
    framework::AxisSpec mTAxis3D = {mTBins3D, "#it{m}_{T} (GeV/#it{c})"};

    // angular correlations
    mPhiLow = (-static_cast<int>(phiBins / 4) + 0.5) * 2. * o2::constants::math::PI / phiBins;
    mPhiHigh = 2 * o2::constants::math::PI + (-static_cast<int>(phiBins / 4) + 0.5) * 2. * o2::constants::math::PI / phiBins;
    framework::AxisSpec phiAxis = {phiBins, mPhiLow, mPhiHigh};
    framework::AxisSpec etaAxis = {etaBins, -2.0, 2.0};

    std::string folderName = static_cast<std::string>(mFolderSuffix[mEventType]) + static_cast<std::string>(o2::aod::femtouniverseMCparticle::MCTypeName[o2::aod::femtouniverseMCparticle::MCType::kRecon]);

    init_base(folderName, femtoObs, femtoObsAxis, multAxis, kTAxis, mTAxis, multAxis3D, mTAxis3D, etaAxis, phiAxis, use3dplots);
    if (isMC) {
      folderName = static_cast<std::string>(mFolderSuffix[mEventType]) + static_cast<std::string>(o2::aod::femtouniverseMCparticle::MCTypeName[o2::aod::femtouniverseMCparticle::MCType::kTruth]);
      init_base(folderName, femtoObs, femtoObsAxis, multAxis, kTAxis, mTAxis, multAxis3D, mTAxis3D, etaAxis, phiAxis, use3dplots);
      init_MC(folderName, femtoObs, femtoObsAxis, multAxis, mTAxis);
    }
  }

  /// Set the PDG codes of the two particles involved
  /// \param pdg1 PDG code of particle one
  /// \param pdg2 PDG code of particle two
  void setPDGCodes(const int pdg1, const int pdg2)
  {
    mMassOne = TDatabasePDG::Instance()->GetParticle(pdg1)->Mass();
    mMassTwo = TDatabasePDG::Instance()->GetParticle(pdg2)->Mass();
    mPDGOne = pdg1;
    mPDGTwo = pdg2;
  }

  /// Pass a pair to the container and compute all the relevant observables
  /// Called by setPair both in case of data/ and Monte Carlo reconstructed and for Monte Carlo truth
  /// \tparam T type of the femtouniverseparticle
  /// \param part1 Particle one
  /// \param part2 Particle two
  /// \param mult Multiplicity of the event
  template <o2::aod::femtouniverseMCparticle::MCType mc, typename T>
  void setPair_base(const float femtoObs, const float mT, T const& part1, T const& part2, const int mult, bool use3dplots)
  {
    delta_eta = part1.eta() - part2.eta();
    delta_phi = part1.phi() - part2.phi();

    while (delta_phi < mPhiLow) {
      delta_phi += o2::constants::math::TwoPI;
    }
    while (delta_phi > mPhiHigh) {
      delta_phi -= o2::constants::math::TwoPI;
    }

    mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/DeltaEtaDeltaPhi"), delta_phi, delta_eta);
    if (use3dplots) {
      // use 3d plots
    }
  }

  /// Called by setPair only in case of Monte Carlo truth
  /// Fills MC truth specific histogramms:
  /// - kstar distribution plots with RECONSTRUCTED information but ONLY for non-fake candidates; needed for purity calculations of tracks
  /// - kstar resolution matrix
  /// Note: Standard histogramms with MC truth information are filled with the setPair_base function
  /// \param part1 Particle one
  /// \param part2 Particle two
  /// \param mult Multiplicity of the event
  void setPair_MC(const float femtoObsMC, const float femtoObs, const float mT, const int mult)
  {
    if (mHistogramRegistry) {
      // Fill the kstar distributions with the reconstructed information but only for particles with the right PDG code
    }
  }

  /// Templated function to handle data/ Monte Carlo reconstructed and Monte Carlo truth
  /// Always calls setPair_base to compute the observables with reconstructed data
  /// In case of Monte Carlo, calls setPair_base with MC info and specialized function setPair_MC for additional histogramms
  /// \tparam T type of the femtouniverseparticle
  /// \param part1 Particle one
  /// \param part2 Particle two
  /// \param mult Multiplicity of the event
  template <bool isMC, typename T>
  void setPair(T const& part1, T const& part2, const int mult, bool use3dplots)
  {
    float femtoObs, femtoObsMC;
    // Calculate femto observable and the mT with reconstructed information
    if constexpr (mFemtoObs == femtoUniverseAngularContainer::Observable::kstar) {
      femtoObs = FemtoUniverseMath::getkstar(part1, mMassOne, part2, mMassTwo);
    }
    const float mT = FemtoUniverseMath::getmT(part1, mMassOne, part2, mMassTwo);

    if (mHistogramRegistry) {
      setPair_base<o2::aod::femtouniverseMCparticle::MCType::kRecon>(femtoObs, mT, part1, part2, mult, use3dplots);

      if constexpr (isMC) {
        if (part1.has_fdMCParticle() && part2.has_fdMCParticle()) {
          // calculate the femto observable and the mT with MC truth information
          if constexpr (mFemtoObs == femtoUniverseAngularContainer::Observable::kstar) {
            femtoObsMC = FemtoUniverseMath::getkstar(part1.fdMCParticle(), mMassOne, part2.fdMCParticle(), mMassTwo);
          }
          const float mTMC = FemtoUniverseMath::getmT(part1.fdMCParticle(), mMassOne, part2.fdMCParticle(), mMassTwo);

          if (abs(part1.fdMCParticle().pdgMCTruth()) == abs(mPDGOne) && abs(part2.fdMCParticle().pdgMCTruth()) == abs(mPDGTwo)) { // Note: all pair-histogramms are filled with MC truth information ONLY in case of non-fake candidates
            setPair_base<o2::aod::femtouniverseMCparticle::MCType::kTruth>(femtoObsMC, mTMC, part1.fdMCParticle(), part2.fdMCParticle(), mult, use3dplots);
            setPair_MC(femtoObsMC, femtoObs, mT, mult);
          } else {
          }

        } else {
        }
      }
    }
  }

 protected:
  HistogramRegistry* mHistogramRegistry = nullptr;                                  ///< For QA output
  static constexpr std::string_view mFolderSuffix[2] = {"SameEvent", "MixedEvent"}; ///< Folder naming for the output according to mEventType
  static constexpr femtoUniverseAngularContainer::Observable mFemtoObs = obs;       ///< Femtoscopic observable to be computed (according to femtoUniverseAngularContainer::Observable)
  static constexpr int mEventType = eventType;                                      ///< Type of the event (same/mixed, according to femtoUniverseAngularContainer::EventType)
  float mMassOne = 0.f;                                                             ///< PDG mass of particle 1
  float mMassTwo = 0.f;                                                             ///< PDG mass of particle 2
  int mPDGOne = 0;                                                                  ///< PDG code of particle 1
  int mPDGTwo = 0;                                                                  ///< PDG code of particle 2
  double mPhiLow;
  double mPhiHigh;
  double delta_eta;
  double delta_phi;
};

} // namespace o2::analysis::femtoUniverse

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEANGULARCONTAINER_H_
