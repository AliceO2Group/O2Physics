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

/// \file FemtoDreamContainer.h
/// \brief Definition of the FemtoDreamContainer
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Valentina Mantovani Sarti, valentina.mantovani-sarti@tum.de
/// \author Georgios Mantzaridis, TU München, georgios.mantzaridis@tum.de
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTODREAM_CORE_FEMTODREAMCONTAINER_H_
#define PWGCF_FEMTODREAM_CORE_FEMTODREAMCONTAINER_H_

#include <fairlogger/Logger.h>
#include <vector>
#include <string>

#include "Framework/HistogramRegistry.h"
#include "PWGCF/FemtoDream/Core/femtoDreamMath.h"
#include "PWGCF/FemtoDream/Core/femtoDreamUtils.h"
#include "PWGCF/DataModel/FemtoDerived.h"

#include "Math/Vector4D.h"
#include "TMath.h"

using namespace o2::framework;

namespace o2::analysis::femtoDream
{

namespace femtoDreamContainer
{
/// Femtoscopic observable to be computed
enum Observable { kstar ///< kstar
};

/// Type of the event processind
enum EventType { same, ///< Pair from same event
                 mixed ///< Pair from mixed event
};
}; // namespace femtoDreamContainer

/// \class FemtoDreamContainer
/// \brief Container for all histogramming related to the correlation function. The two
/// particles of the pair are passed here, and the correlation function and QA histograms
/// are filled according to the specified observable
/// \tparam eventType Type of the event (same/mixed)
/// \tparam obs Observable to be computed (k*/Q_inv/...)
template <femtoDreamContainer::EventType eventType, femtoDreamContainer::Observable obs>
class FemtoDreamContainer
{
 public:
  /// Destructor
  virtual ~FemtoDreamContainer() = default;

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
  void init_base(std::string folderName, std::string femtoObs,
                 T& femtoObsAxis, T& pTAxis, T& kTAxis, T& mTAxis, T& multAxis, T& multPercentileAxis,
                 T& /*kstarAxis4D*/, T& mTAxis4D, T& multAxis4D, T& multPercentileAxis4D,
                 bool use4dplots, bool extendedplots)
  {

    mHistogramRegistry->add((folderName + "/relPairDist").c_str(), ("; " + femtoObs + "; Entries").c_str(), kTH1F, {femtoObsAxis});
    mHistogramRegistry->add((folderName + "/relPairkT").c_str(), "; #it{k}_{T} (GeV/#it{c}); Entries", kTH1F, {kTAxis});
    mHistogramRegistry->add((folderName + "/relPairkstarkT").c_str(), ("; " + femtoObs + "; #it{k}_{T} (GeV/#it{c})").c_str(), kTH2F, {femtoObsAxis, kTAxis});
    mHistogramRegistry->add((folderName + "/relPairkstarmT").c_str(), ("; " + femtoObs + "; #it{m}_{T} (GeV/#it{c}^{2})").c_str(), kTH2F, {femtoObsAxis, mTAxis});
    mHistogramRegistry->add((folderName + "/relPairkstarMult").c_str(), ("; " + femtoObs + "; Multiplicity").c_str(), kTH2F, {femtoObsAxis, multAxis});
    mHistogramRegistry->add((folderName + "/relPairkstarMultPercentile").c_str(), ("; " + femtoObs + "; Multiplicity Percentile").c_str(), kTH2F, {femtoObsAxis, multPercentileAxis4D});
    mHistogramRegistry->add((folderName + "/kstarPtPart1").c_str(), ("; " + femtoObs + "; #it{p} _{T} Particle 1 (GeV/#it{c})").c_str(), kTH2F, {femtoObsAxis, pTAxis});
    mHistogramRegistry->add((folderName + "/kstarPtPart2").c_str(), ("; " + femtoObs + "; #it{p} _{T} Particle 2 (GeV/#it{c})").c_str(), kTH2F, {femtoObsAxis, pTAxis});
    mHistogramRegistry->add((folderName + "/MultPtPart1").c_str(), "; #it{p} _{T} Particle 1 (GeV/#it{c}); Multiplicity", kTH2F, {pTAxis, multAxis});
    mHistogramRegistry->add((folderName + "/MultPtPart2").c_str(), "; #it{p} _{T} Particle 2 (GeV/#it{c}); Multiplicity", kTH2F, {pTAxis, multAxis});
    mHistogramRegistry->add((folderName + "/MultPercentilePtPart1").c_str(), "; #it{p} _{T} Particle 1 (GeV/#it{c}); Multiplicity Percentile", kTH2F, {pTAxis, multPercentileAxis});
    mHistogramRegistry->add((folderName + "/MultPercentilePtPart2").c_str(), "; #it{p} _{T} Particle 2 (GeV/#it{c}); Multiplicity Percentile", kTH2F, {pTAxis, multPercentileAxis});
    mHistogramRegistry->add((folderName + "/PtPart1PtPart2").c_str(), "; #it{p} _{T} Particle 1 (GeV/#it{c}); #it{p} _{T} Particle 2 (GeV/#it{c})", kTH2F, {pTAxis, pTAxis});
    if (use4dplots) {
      mHistogramRegistry->add((folderName + "/relPairkstarmTMultMultPercentile").c_str(), ("; " + femtoObs + "; #it{m}_{T} (GeV/#it{c}^{2}); Multiplicity").c_str(), kTHnSparseF, {femtoObsAxis, mTAxis4D, multAxis4D, multPercentileAxis4D});
    }
    if (extendedplots) {
      mHistogramRegistry->add((folderName + "/mTPtPart1PtPart2MultPercentile").c_str(), "; #it{m}_{T} (GeV/#it{c}^{2}); #it{p} _{T} Particle 1 (GeV/#it{c}); #it{p} _{T} Particle 2 (GeV/#it{c}); Multiplicity Percentile (%)", kTHnSparseF, {mTAxis4D, pTAxis, pTAxis, multPercentileAxis4D});
    }
  }

  /// Initializes specialized Monte Carlo truth histograms for the task
  /// internal function called by init only in case of Monte Carlo truth
  /// \tparam T type of the xxis Object
  /// \param folderName Name of the directory in the output file (no suffix for reconstructed data/ Monte Carlo; "_MC" for Monte Carlo Truth)
  /// \param femtoObsAxis axis object for the femto observable axis
  template <typename T>
  void init_MC(std::string folderName, std::string femtoObs, T femtoObsAxis, T multAxis, T mTAxis, bool smearingByOrigin)
  {
    mHistogramRegistry->add((folderName + "/relPairDist_ReconNoFake").c_str(), ("; " + femtoObs + "; Entries").c_str(), kTH1F, {femtoObsAxis});
    mHistogramRegistry->add((folderName + "/relPairkstarmT_ReconNoFake").c_str(), ("; " + femtoObs + "; #it{m}_{T} (GeV/#it{c}^{2})").c_str(), kTH2F, {femtoObsAxis, mTAxis});
    mHistogramRegistry->add((folderName + "/relPairkstarMult_ReconNoFake").c_str(), ("; " + femtoObs + "; Multiplicity").c_str(), kTH2F, {femtoObsAxis, multAxis});
    mHistogramRegistry->add((folderName + "/hNoMCtruthPairsCounter").c_str(), "; Counter; Entries", kTH1I, {{1, 0, 1}});
    mHistogramRegistry->add((folderName + "/hFakePairsCounter").c_str(), "; Counter; Entries", kTH1I, {{1, 0, 1}});
    if (smearingByOrigin) {
      const int nOriginBins = o2::aod::femtodreamMCparticle::ParticleOriginMCTruth::kNOriginMCTruthTypes;
      // framework::AxisSpec mcOriginAxisPart1 = {{nOriginBins, 0, nOriginBins}, "MC origin particle 1"};
      // framework::AxisSpec mcOriginAxisPart2 = {{nOriginBins, 0, nOriginBins}, "MC origin particle 2"};
      mHistogramRegistry->add((folderName + "/kstar_resolution").c_str(), "; #it{k}_{*} reconstructed (GeV/#it{c}); #it{k}_{*} truth (GeV/#it{c}); MC origin particle 1; MC origin particle 2; ", kTHnSparseF, {femtoObsAxis, femtoObsAxis, {nOriginBins, 0, nOriginBins}, {nOriginBins, 0, nOriginBins}});
    } else {
      mHistogramRegistry->add((folderName + "/kstar_resolution").c_str(), "; #it{k}_{*} reconstructed (GeV/#it{c}); #it{k}_{*} truth (GeV/#it{c})", kTH2F, {femtoObsAxis, femtoObsAxis});
    }
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
  /// \param isMC add Monte Carlo truth histograms to the output file
  template <typename T>
  void init(HistogramRegistry* registry,
            T& kstarBins, T& pTBins, T& kTBins, T& mTBins, T& multBins, T& /*multPercentileBins*/,
            T& kstarBins4D, T& mTBins4D, T& multBins4D, T& multPercentileBins4D,
            bool isMC, bool use4dplots, bool extendedplots,
            float highkstarCut,
            bool smearingByOrigin = false)
  {
    mHistogramRegistry = registry;
    std::string femtoObs;
    if constexpr (mFemtoObs == femtoDreamContainer::Observable::kstar) {
      femtoObs = "#it{k*} (GeV/#it{c})";
    }
    mHighkstarCut = highkstarCut;
    framework::AxisSpec femtoObsAxis = {kstarBins, femtoObs.c_str()};
    framework::AxisSpec pTAxis = {pTBins, "#it{p}_{T} (GeV/#it{c})"};
    framework::AxisSpec kTAxis = {kTBins, "#it{k}_{T} (GeV/#it{c})"};
    framework::AxisSpec mTAxis = {mTBins, "#it{m}_{T} (GeV/#it{c}^{2})"};
    framework::AxisSpec multAxis = {multBins, "Multiplicity"};
    framework::AxisSpec multPercentileAxis = {multBins, "Multiplicity percentile (%)"};

    framework::AxisSpec kstarAxis4D = {kstarBins4D, "Multiplicity"};
    framework::AxisSpec mTAxis4D = {mTBins4D, "#it{m}_{T} (GeV/#it{c})"};
    framework::AxisSpec multAxis4D = {multBins4D, "Multiplicity"};
    framework::AxisSpec multPercentileAxis4D = {multPercentileBins4D, "Multiplicity Percentile (%)"};

    std::string folderName = static_cast<std::string>(mFolderSuffix[mEventType]) + static_cast<std::string>(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kRecon]);

    init_base(folderName, femtoObs,
              femtoObsAxis, pTAxis, kTAxis, mTAxis, multAxis, multPercentileAxis,
              kstarAxis4D, mTAxis4D, multAxis4D, multPercentileAxis4D,
              use4dplots, extendedplots);
    if (isMC) {
      folderName = static_cast<std::string>(mFolderSuffix[mEventType]) + static_cast<std::string>(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kTruth]);
      init_base(folderName, femtoObs,
                femtoObsAxis, pTAxis, kTAxis, mTAxis, multAxis, multPercentileAxis,
                kstarAxis4D, mTAxis4D, multAxis4D, multPercentileAxis4D,
                use4dplots, extendedplots);
      init_MC(folderName, femtoObs, femtoObsAxis, multAxis, mTAxis, smearingByOrigin);
    }
  }

  /// Set the PDG codes of the two particles involved
  /// \param pdg1 PDG code of particle one
  /// \param pdg2 PDG code of particle two
  void setPDGCodes(const int pdg1, const int pdg2)
  {
    mMassOne = o2::analysis::femtoDream::getMass(pdg1);
    mMassTwo = o2::analysis::femtoDream::getMass(pdg2);
    mPDGOne = pdg1;
    mPDGTwo = pdg2;
  }

  /// Pass a pair to the container and compute all the relevant observables
  /// Called by setPair both in case of data/ and Monte Carlo reconstructed and for Monte Carlo truth
  /// \tparam T type of the femtodreamparticle
  /// \param part1 Particle one
  /// \param part2 Particle two
  /// \param mult Multiplicity of the event
  template <o2::aod::femtodreamMCparticle::MCType mc, typename T>
  void setPair_base(const float femtoObs, const float mT, T const& part1, T const& part2, const int mult, const float multPercentile, bool use4dplots, bool extendedplots)
  {
    const float kT = FemtoDreamMath::getkT(part1, mMassOne, part2, mMassTwo);
    mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/relPairDist"), femtoObs);
    mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/relPairkT"), kT);
    mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/relPairkstarkT"), femtoObs, kT);
    mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/relPairkstarmT"), femtoObs, mT);
    mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/relPairkstarMult"), femtoObs, mult);
    mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/relPairkstarMultPercentile"), femtoObs, multPercentile);
    mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/kstarPtPart1"), femtoObs, part1.pt());
    mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/kstarPtPart2"), femtoObs, part2.pt());
    mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/MultPtPart1"), part1.pt(), mult);
    mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/MultPtPart2"), part2.pt(), mult);
    mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/MultPercentilePtPart1"), part1.pt(), multPercentile);
    mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/MultPercentilePtPart2"), part2.pt(), multPercentile);
    mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/PtPart1PtPart2"), part1.pt(), part2.pt());
    if (use4dplots) {
      mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/relPairkstarmTMultMultPercentile"), femtoObs, mT, mult, multPercentile);
    }
    if (extendedplots) {
      mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/mTPtPart1PtPart2MultPercentile"), mT, part1.pt(), part2.pt(), multPercentile);
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
  void setPair_MC(const float femtoObsMC, const float femtoObs, const float mT, const int mult, const int originPart1, const int originPart2, bool smearingByOrigin)
  {
    if (mHistogramRegistry) {
      // Fill the kstar distributions with the reconstructed information but only for particles with the right PDG code
      mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kTruth]) + HIST("/relPairDist_ReconNoFake"), femtoObs);
      mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kTruth]) + HIST("/relPairkstarmT_ReconNoFake"), femtoObs, mT);
      mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kTruth]) + HIST("/relPairkstarMult_ReconNoFake"), femtoObs, mult);
      if (smearingByOrigin) {
        mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kTruth]) + HIST("/kstar_resolution"), femtoObsMC, femtoObs, originPart1, originPart2);
      } else {
        mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kTruth]) + HIST("/kstar_resolution"), femtoObsMC, femtoObs);
      }
    }
  }

  /// Templated function to handle data/ Monte Carlo reconstructed and Monte Carlo truth
  /// Always calls setPair_base to compute the observables with reconstructed data
  /// In case of Monte Carlo, calls setPair_base with MC info and specialized function setPair_MC for additional histogramms
  /// \tparam T type of the femtodreamparticle
  /// \param part1 Particle one
  /// \param part2 Particle two
  /// \param mult Multiplicity of the event
  template <bool isMC, typename T>
  void setPair(T const& part1, T const& part2, const int mult, const float multPercentile, bool use4dplots, bool extendedplots, bool smearingByOrigin = false)
  {
    float femtoObs, femtoObsMC;
    // Calculate femto observable and the mT with reconstructed information
    if constexpr (mFemtoObs == femtoDreamContainer::Observable::kstar) {
      femtoObs = FemtoDreamMath::getkstar(part1, mMassOne, part2, mMassTwo);
    }
    if (mHighkstarCut > 0) {
      if (femtoObs > mHighkstarCut) {
        return;
      }
    }
    const float mT = FemtoDreamMath::getmT(part1, mMassOne, part2, mMassTwo);

    if (mHistogramRegistry) {
      setPair_base<o2::aod::femtodreamMCparticle::MCType::kRecon>(femtoObs, mT, part1, part2, mult, multPercentile, use4dplots, extendedplots);

      if constexpr (isMC) {
        if (part1.has_fdMCParticle() && part2.has_fdMCParticle()) {
          // calculate the femto observable and the mT with MC truth information
          if constexpr (mFemtoObs == femtoDreamContainer::Observable::kstar) {
            femtoObsMC = FemtoDreamMath::getkstar(part1.fdMCParticle(), mMassOne, part2.fdMCParticle(), mMassTwo);
          }
          const float mTMC = FemtoDreamMath::getmT(part1.fdMCParticle(), mMassOne, part2.fdMCParticle(), mMassTwo);

          if (abs(part1.fdMCParticle().pdgMCTruth()) == mPDGOne && abs(part2.fdMCParticle().pdgMCTruth()) == mPDGTwo) { // Note: all pair-histogramms are filled with MC truth information ONLY in case of non-fake candidates
            setPair_base<o2::aod::femtodreamMCparticle::MCType::kTruth>(femtoObsMC, mTMC, part1.fdMCParticle(), part2.fdMCParticle(), mult, multPercentile, use4dplots, extendedplots);
            setPair_MC(femtoObsMC, femtoObs, mT, mult, part1.fdMCParticle().partOriginMCTruth(), part2.fdMCParticle().partOriginMCTruth(), smearingByOrigin);
          } else {
            mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kTruth]) + HIST("/hFakePairsCounter"), 0);
          }

        } else {
          mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kTruth]) + HIST("/hNoMCtruthPairsCounter"), 0);
        }
      }
    }
  }

 protected:
  HistogramRegistry* mHistogramRegistry = nullptr;                                  ///< For QA output
  static constexpr std::string_view mFolderSuffix[2] = {"SameEvent", "MixedEvent"}; ///< Folder naming for the output according to mEventType
  static constexpr femtoDreamContainer::Observable mFemtoObs = obs;                 ///< Femtoscopic observable to be computed (according to femtoDreamContainer::Observable)
  static constexpr int mEventType = eventType;                                      ///< Type of the event (same/mixed, according to femtoDreamContainer::EventType)
  float mMassOne = 0.f;                                                             ///< PDG mass of particle 1
  float mMassTwo = 0.f;                                                             ///< PDG mass of particle 2
  int mPDGOne = 0;                                                                  ///< PDG code of particle 1
  int mPDGTwo = 0;                                                                  ///< PDG code of particle 2
  float mHighkstarCut = 6.;
};

} // namespace o2::analysis::femtoDream

#endif // PWGCF_FEMTODREAM_CORE_FEMTODREAMCONTAINER_H_
