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

/// \file FemtoWorldContainer.h
/// \brief Definition of the FemtoWorldContainer
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Valentina Mantovani Sarti, valentina.mantovani-sarti@tum.de
/// \author Zuzanna Chochulska, WUT Warsaw, zchochul@cern.ch NOWE

#ifndef PWGCF_FEMTOWORLD_CORE_FEMTOWORLDCONTAINER_H_
#define PWGCF_FEMTOWORLD_CORE_FEMTOWORLDCONTAINER_H_

#include <vector>
#include <string>
#include "Framework/HistogramRegistry.h"
#include "PWGCF/FemtoWorld/Core/FemtoWorldMath.h"

#include "Math/Vector4D.h"
#include "TMath.h"
#include "TDatabasePDG.h"

#include "TLorentzVector.h"
#include "CommonConstants/MathConstants.h"
#include "TRandom.h"

using namespace o2::framework;

namespace o2::analysis::femtoWorld
{

namespace femtoWorldContainer
{
/// Femtoscopic observable to be computed
enum Observable { kstar ///< kstar
};

/// Type of the event processind
enum EventType { same, ///< Pair from same event
                 mixed ///< Pair from mixed event
};
}; // namespace femtoWorldContainer

/// \class FemtoWorldContainer
/// \brief Container for all histogramming related to the correlation function. The two
/// particles of the pair are passed here, and the correlation function and QA histograms
/// are filled according to the specified observable
/// \tparam eventType Type of the event (same/mixed)
/// \tparam obs Observable to be computed (k*/Q_inv/...)
template <femtoWorldContainer::EventType eventType, femtoWorldContainer::Observable obs>
class FemtoWorldContainer
{
 public:
  /// Destructor
  virtual ~FemtoWorldContainer() = default;

  /// Initializes histograms for the task
  /// \tparam T Type of the configurable for the axis configuration
  /// \param registry Histogram registry to be passed
  /// \param kstarBins k* binning for the histograms
  /// \param multBins multiplicity binning for the histograms
  /// \param kTBins kT binning for the histograms
  /// \param mTBins mT binning for the histograms
  /// \param etaBins eta binning for the histograms
  /// \param phiBins phi binning for the histograms
  /// \param mInvBins invariant mass binning for the histograms

  template <typename T1, typename T2>
  void init(HistogramRegistry* registry, T1& kstarBins, T1& multBins, T1& kTBins, T1& mTBins, T2& phiBins, T2& etaBins, T2& mInvBins)
  {
    mHistogramRegistry = registry;
    std::string femtoObs;
    if constexpr (mFemtoObs == femtoWorldContainer::Observable::kstar) {
      femtoObs = "#it{k*} (GeV/#it{c})";
    }
    std::vector<double> tmpVecMult = multBins;
    framework::AxisSpec multAxis = {tmpVecMult, "Multiplicity"};
    framework::AxisSpec femtoObsAxis = {kstarBins, femtoObs.c_str()};
    framework::AxisSpec kTAxis = {kTBins, "#it{k}_{T} (GeV/#it{c})"};
    framework::AxisSpec mTAxis = {mTBins, "#it{m}_{T} (GeV/#it{c}^{2})"};

    mPhiLow = (-static_cast<int>(phiBins / 4) + 0.5) * 2. * o2::constants::math::PI / phiBins;
    mPhiHigh = 2 * o2::constants::math::PI + (-static_cast<int>(phiBins / 4) + 0.5) * 2. * o2::constants::math::PI / phiBins;

    framework::AxisSpec phiAxis = {phiBins, mPhiLow, mPhiHigh};
    framework::AxisSpec etaAxis = {etaBins, -2.0, 2.0};
    framework::AxisSpec mInvAxis = {mInvBins, 0.0, 10.0};

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
    mHistogramRegistry->add((folderName + "relPairDetaDphi").c_str(), ";  #Delta#varphi (rad); #Delta#eta", kTH2D, {phiAxis, etaAxis});
    mHistogramRegistry->add((folderName + "relPairInvariantMass").c_str(), ";M_{K^{+}K^{-}} (GeV/#it{c}^{2});", kTH1D, {mInvAxis});
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
  /// \tparam T type of the femtoworldparticle
  /// \param part1 Particle one
  /// \param part2 Particle two
  /// \param mult Multiplicity of the event
  template <typename T>
  void setPair(T const& part1, T const& part2, const int mult)
  {
    float femtoObs;
    if constexpr (mFemtoObs == femtoWorldContainer::Observable::kstar) {
      femtoObs = FemtoWorldMath::getkstar(part1, mMassOne, part2, mMassTwo);
    }
    const float kT = FemtoWorldMath::getkT(part1, mMassOne, part2, mMassTwo);
    const float mT = FemtoWorldMath::getmT(part1, mMassOne, part2, mMassTwo);

    // // HERE RANDOMLY CHOOSE PARTICLES
    // TRandom* rndm = new TRandom(0);
    // double ran = rndm->Rndm();

    // double delta_eta;
    // double delta_phi;

    // if (ran < 0.5) {
    //   delta_eta = part1.eta() - part2.eta();
    //   delta_phi = part1.phi() - part2.phi();
    // } else if (ran >= 0.5) {
    //   delta_eta = part2.eta() - part1.eta();
    //   delta_phi = part2.phi() - part1.phi();
    // } else {
    //   LOGF(error, "-------Error in FemtoWorldContainer (randomly picking particles in a pair) - value %i, not 0 or 1", ran);
    // }

    double delta_eta = part1.eta() - part2.eta();
    double delta_phi = part1.phi() - part2.phi();

    while (delta_phi < mPhiLow) {
      delta_phi += o2::constants::math::TwoPI;
    }
    while (delta_phi > mPhiHigh) {
      delta_phi -= o2::constants::math::TwoPI;
    }
    TLorentzVector part1Vec;
    part1Vec.SetPtEtaPhiM(part1.pt(), part1.eta(), part1.phi(), mMassOne);
    TLorentzVector part2Vec;
    part2Vec.SetPtEtaPhiM(part2.pt(), part2.eta(), part2.phi(), mMassTwo);

    TLorentzVector sumVec(part1Vec);
    sumVec += part2Vec;

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
      mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST("relPairDetaDphi"), delta_phi, delta_eta);
      mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST("relPairInvariantMass"), sumVec.M());
    }
  }

 protected:
  HistogramRegistry* mHistogramRegistry = nullptr;                                    ///< For QA output
  static constexpr std::string_view mFolderSuffix[2] = {"SameEvent/", "MixedEvent/"}; ///< Folder naming for the output according to mEventType
  static constexpr femtoWorldContainer::Observable mFemtoObs = obs;                   ///< Femtoscopic observable to be computed (according to femtoWorldContainer::Observable)
  static constexpr int mEventType = eventType;                                        ///< Type of the event (same/mixed, according to femtoWorldContainer::EventType)
  float mMassOne = 0.f;                                                               ///< PDG mass of particle 1
  float mMassTwo = 0.f;                                                               ///< PDG mass of particle 2
  double mPhiLow;
  double mPhiHigh;
};

} // namespace o2::analysis::femtoWorld

#endif // PWGCF_FEMTOWORLD_CORE_FEMTOWORLDCONTAINER_H_
