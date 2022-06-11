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

/// \file AliFemtoContainer.h
/// \brief Definition of the AliFemtoContainer
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Valentina Mantovani Sarti, valentina.mantovani-sarti@tum.de

#ifndef ANALYSIS_TASKS_PWGCF_ALIFEMTO_ALIFEMTOCONTAINER_H_
#define ANALYSIS_TASKS_PWGCF_ALIFEMTO_ALIFEMTOCONTAINER_H_

#include "Framework/HistogramRegistry.h"
#include "AliFemtoMath.h"

#include "Math/Vector4D.h"
#include "TMath.h"
#include "TDatabasePDG.h"

using namespace o2::framework;

namespace o2::analysis::aliFemto
{

namespace aliFemtoContainer
{
/// Femtoscopic observable to be computed
enum Observable { kstar ///< kstar
};

/// Type of the event processind
enum EventType { same, ///< Pair from same event
                 mixed ///< Pair from mixed event
};
}; // namespace aliFemtoContainer

/// \class AliFemtoContainer
/// \brief Container for all histogramming related to the correlation function. The two
/// particles of the pair are passed here, and the correlation function and QA histograms
/// are filled according to the specified observable
/// \tparam eventType Type of the event (same/mixed)
/// \tparam obs Observable to be computed (k*/Q_inv/...)
template <aliFemtoContainer::EventType eventType, aliFemtoContainer::Observable obs>
class AliFemtoContainer
{
 public:
  /// Destructor
  virtual ~AliFemtoContainer() = default;

  /// Initializes histograms for the task
  /// \tparam T Type of the configurable for the axis configuration
  /// \param registry Histogram registry to be passed
  /// \param kstarBins k* binning for the histograms
  /// \param multBins multiplicity binning for the histograms
  /// \param kTBins kT binning for the histograms
  /// \param mTBins mT binning for the histograms
  template <typename T>
  void init(HistogramRegistry* registry, T& kstarBins, T& multBins, T& kTBins, T& mTBins)
  {
    mHistogramRegistry = registry;
    std::string femtoObs;
    if constexpr (mFemtoObs == aliFemtoContainer::Observable::kstar) {
      femtoObs = "#it{k*} (GeV/#it{c})";
    }
    std::vector<double> tmpVecMult = multBins;
    framework::AxisSpec multAxis = {tmpVecMult, "Multiplicity"};
    framework::AxisSpec femtoObsAxis = {kstarBins, femtoObs.c_str()};
    framework::AxisSpec kTAxis = {kTBins, "#it{k}_{T} (GeV/#it{c})"};
    framework::AxisSpec mTAxis = {mTBins, "#it{m}_{T} (GeV/#it{c}^{2})"};

    std::string folderName = static_cast<std::string>(mFolderSuffix[mEventType]);
    mHistogramRegistry->add((folderName + "relPairDist").c_str(), ("; " + femtoObs + "; Entries").c_str(), kTH1F, {femtoObsAxis});
    mHistogramRegistry->add((folderName + "relPairkT").c_str(), "; #it{k}_{T} (GeV/#it{c}); Entries", kTH1F, {kTAxis});
    mHistogramRegistry->add((folderName + "relPairkstarkT").c_str(), ("; " + femtoObs + "; #it{k}_{T} (GeV/#it{c})").c_str(), kTH2F, {femtoObsAxis, kTAxis});
    mHistogramRegistry->add((folderName + "relPairkstarmT").c_str(), ("; " + femtoObs + "; #it{m}_{T} (GeV/#it{c}^{2})").c_str(), kTH2F, {femtoObsAxis, mTAxis});
    mHistogramRegistry->add((folderName + "relPairkstarMult").c_str(), ("; " + femtoObs + "; Multiplicity").c_str(), kTH2F, {femtoObsAxis, multAxis});
    mHistogramRegistry->add((folderName + "kstarPtPart1").c_str(), ("; " + femtoObs + "; #it{p} _{T} Particle 1 (GeV/#it{c})").c_str(), kTH2F, {femtoObsAxis, {375, 0., 7.5}});
    mHistogramRegistry->add((folderName + "kstarPtPart2").c_str(), ("; " + femtoObs + "; #it{p} _{T} Particle 2 (GeV/#it{c})").c_str(), kTH2F, {femtoObsAxis, {375, 0., 7.5}});
    mHistogramRegistry->add((folderName + "MultPtPart1").c_str(), "; #it{p} _{T} Particle 1 (GeV/#it{c}); Multiplicity", kTH2F, {{375, 0., 7.5}, multAxis});
    mHistogramRegistry->add((folderName + "MultPtPart2").c_str(), "; #it{p} _{T} Particle 2 (GeV/#it{c}); Multiplicity", kTH2F, {{375, 0., 7.5}, multAxis});
    mHistogramRegistry->add((folderName + "PtPart1PtPart2").c_str(), "; #it{p} _{T} Particle 1 (GeV/#it{c}); #it{p} _{T} Particle 2 (GeV/#it{c})", kTH2F, {{375, 0., 7.5}, {375, 0., 7.5}});
  }

  /// Set the PDG codes of the two particles involved
  /// \param pdg1 PDG code of particle one
  /// \param pdg2 PDG code of particle two
  void setPDGCodes(const int pdg1, const int pdg2)
  {
    mMassOne = TDatabasePDG::Instance()->GetParticle(pdg1)->Mass();
    mMassTwo = TDatabasePDG::Instance()->GetParticle(pdg2)->Mass();
  }

  /// Pass a pair to the container and compute all the relevant observables
  /// \tparam T type of the alifemtoparticle
  /// \param part1 Particle one
  /// \param part2 Particle two
  /// \param mult Multiplicity of the event
  template <typename T>
  void setPair(T const& part1, T const& part2, const int mult)
  {
    float femtoObs;
    if constexpr (mFemtoObs == aliFemtoContainer::Observable::kstar) {
      femtoObs = AliFemtoMath::getkstar(part1, mMassOne, part2, mMassTwo);
    }
    const float kT = AliFemtoMath::getkT(part1, mMassOne, part2, mMassTwo);
    const float mT = AliFemtoMath::getmT(part1, mMassOne, part2, mMassTwo);

    if (mHistogramRegistry) {
      mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST("relPairDist"), femtoObs);
      mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST("relPairkT"), kT);
      mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST("relPairkstarkT"), femtoObs, kT);
      mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST("relPairkstarmT"), femtoObs, mT);
      mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST("relPairkstarMult"), femtoObs, mult);
      mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST("kstarPtPart1"), femtoObs, part1.pt());
      mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST("kstarPtPart2"), femtoObs, part2.pt());
      mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST("MultPtPart1"), part1.pt(), mult);
      mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST("MultPtPart2"), part2.pt(), mult);
      mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST("PtPart1PtPart2"), part1.pt(), part2.pt());
    }
  }

 protected:
  HistogramRegistry* mHistogramRegistry = nullptr;                                    ///< For QA output
  static constexpr std::string_view mFolderSuffix[2] = {"SameEvent/", "MixedEvent/"}; ///< Folder naming for the output according to mEventType
  static constexpr aliFemtoContainer::Observable mFemtoObs = obs;                   ///< Femtoscopic observable to be computed (according to aliFemtoContainer::Observable)
  static constexpr int mEventType = eventType;                                        ///< Type of the event (same/mixed, according to aliFemtoContainer::EventType)
  float mMassOne = 0.f;                                                               ///< PDG mass of particle 1
  float mMassTwo = 0.f;                                                               ///< PDG mass of particle 2
};

} // namespace o2::analysis::aliFemto

#endif /* ANALYSIS_TASKS_PWGCF_ALIFEMTO_ALIFEMTOCONTAINER_H_ */
