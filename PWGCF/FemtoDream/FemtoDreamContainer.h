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

/// \file FemtoDreamContainer.h
/// \brief Definition of the FemtoDreamContainer
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Valentina Mantovani Sarti, valentina.mantovani-sarti@tum.de

#ifndef PWGCF_FEMTODREAM_FEMTODREAMCONTAINER_H_
#define PWGCF_FEMTODREAM_FEMTODREAMCONTAINER_H_

#include <vector>
#include <string>

#include "Framework/HistogramRegistry.h"
#include "FemtoDreamMath.h"

#include "Math/Vector4D.h"
#include "TMath.h"
#include "TDatabasePDG.h"

using namespace o2::framework;

namespace o2::analysis::femtoDream
{

namespace femtoDreamContainer
{
/// Femtoscopic observable to be computed
enum Observable { kstar ///< kstar
};

/// Type of the event processind
enum EventType { same,  ///< Pair from same event
                 mixed, ///< Pair from mixed event
                 sameMC,
                 mixedMC
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
    if constexpr (mFemtoObs == femtoDreamContainer::Observable::kstar) {
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
  /// \tparam T type of the femtodreamparticle
  /// \param part1 Particle one
  /// \param part2 Particle two
  /// \param mult Multiplicity of the event
  template <typename T>
  void setPair(T const& part1, T const& part2, const int mult)
  {
    float femtoObs;
    if constexpr (mFemtoObs == femtoDreamContainer::Observable::kstar) {
      femtoObs = FemtoDreamMath::getkstar(part1, mMassOne, part2, mMassTwo);
    }
    const float kT = FemtoDreamMath::getkT(part1, mMassOne, part2, mMassTwo);
    const float mT = FemtoDreamMath::getmT(part1, mMassOne, part2, mMassTwo);

    /*
    if(part1.pt()==part2.pt()){
    std::cout<<"New particle";
    std::cout<<"particle 1";
    std::cout<<"pt = "<<part1.pt();
    std::cout<<"eta = "<<part1.eta();
    std::cout<<"phi = "<<part1.phi();
    std::cout<<"particle 2";
    std::cout<<"pt = "<<part2.pt();
    std::cout<<"eta = "<<part2.eta();
    std::cout<<"phi = "<<part2.phi();
    std::cout<<" ";
    std::cout<<"kstar = "<<femtoObs;
    }
    */

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

  template <typename T>
  void initMC(HistogramRegistry* registry, T& kstarBins)
  {
    mHistogramRegistry = registry;
    std::string femtoObs;
    if constexpr (mFemtoObs == femtoDreamContainer::Observable::kstar) {
      femtoObs = "#it{k*} (GeV/#it{c})";
    }
    framework::AxisSpec femtoObsAxis = {kstarBins, femtoObs.c_str()};
    std::string folderName = static_cast<std::string>(mFolderSuffix[mEventType]);
    mHistogramRegistry->add((folderName + "kstar_resolution").c_str(), "; #it{k} _{T} reconstructed (GeV/#it{c}); #it{k} _{T} truth (GeV/#it{c})", kTH2F, {femtoObsAxis, femtoObsAxis});
  }

  /// Pass a pair to the container and compute all the relevant observables
  /// \tparam T type of the femtodreamparticle
  /// \param part1 Particle one
  /// \param part2 Particle two
  /// \param mult Multiplicity of the event
  template <typename T, typename TMC>
  void setPairMC(T const& part1, T const& part2, TMC const& partMC1, TMC const& partMC2, const int mult)
  {
    float femtoObs, femtoObsMC;
    if constexpr (mFemtoObs == femtoDreamContainer::Observable::kstar) {
      femtoObs = FemtoDreamMath::getkstar(part1, mMassOne, part2, mMassTwo);
      femtoObsMC = FemtoDreamMath::getkstar(partMC1, mMassOne, partMC2, mMassTwo);
    }

    std::cout << " " << std::endl;
    std::cout << "New particle" << std::endl;
    std::cout << "particle 1" << std::endl;
    std::cout << "reco     -----     sim" << std::endl;
    std::cout << "pt = " << part1.pt() << " -- pt = " << partMC1.pt() << std::endl;
    std::cout << "eta = " << part1.eta() << " -- eta = " << partMC1.eta() << std::endl;
    std::cout << "phi = " << part1.phi() << " -- phi = " << partMC1.phi() << std::endl;
    std::cout << "CollisionID = " << part1.femtoDreamCollisionId() << " -- CollisionID = " << partMC1.femtoDreamCollisionId() << std::endl;
    std::cout << "Particle ID= " << part1.index() << " -- Particle ID= " << partMC1.index() << std::endl;
    std::cout << "Particle ID global= " << part1.globalIndex() << " -- Particle ID global= " << partMC1.globalIndex() << std::endl;
    if (part1.index() == partMC1.index()) {
      std::cout << " equal index 1" << std::endl;
    }
    std::cout << "  mc PDG: " << partMC1.pdgMCTruth() << std::endl;
    std::cout << "  mc origin: " << static_cast<int>(partMC1.partOriginMCTruth()) << std::endl;

    std::cout << "### particle 2 ###" << std::endl;
    std::cout << "reco     -----     sim" << std::endl;
    std::cout << "pt = " << part2.pt() << " -- pt = " << partMC2.pt() << std::endl;
    std::cout << "eta = " << part2.eta() << " -- eta = " << partMC2.eta() << std::endl;
    std::cout << "phi = " << part2.phi() << " -- phi = " << partMC2.phi() << std::endl;
    std::cout << "CollisionID = " << part2.femtoDreamCollisionId() << " -- CollisionID = " << partMC2.femtoDreamCollisionId() << std::endl;
    std::cout << "Particle ID= " << part2.index() << " -- Particle ID= " << partMC2.index() << std::endl;
    std::cout << "Particle ID global= " << part2.globalIndex() << " -- Particle ID global= " << partMC2.globalIndex() << std::endl;
    if (part2.index() == partMC2.index()) {
      std::cout << " equal index 2" << std::endl;
    }
    std::cout << "  mc PDG: " << partMC2.pdgMCTruth() << std::endl;
    std::cout << "  mc origin: " << static_cast<int>(partMC2.partOriginMCTruth()) << std::endl;
    std::cout << "--> kstar = " << femtoObs << " -- kstar truth = " << femtoObsMC << std::endl;
    if (femtoObs < 0.05 || femtoObsMC < 0.05) {
      std::cout << "==> zks";
    }
    if (part1.femtoDreamCollisionId() != partMC1.femtoDreamCollisionId() || part2.femtoDreamCollisionId() != partMC2.femtoDreamCollisionId()) {
      std::cout << " necID";
    }
    if (part1.index() != partMC1.index() || part2.index() != partMC2.index()) {
      std::cout << " nei";
    }
    if (part1.globalIndex() != partMC1.globalIndex() || part2.globalIndex() != partMC2.globalIndex()) {
      std::cout << " negi";
    }
    std::cout << " " << std::endl;
    // if( !( (femtoObs<0.01 || femtoObsMC<0.01) && (partMC1.femtoDreamCollisionId()==0 || partMC1.femtoDreamCollisionId()==0 || partMC2.femtoDreamCollisionId()==1 || partMC2.femtoDreamCollisionId()==1) )){
    //    std::cout<<"zerokstar"<<std::endl;
    // }

    if (mHistogramRegistry) {
      mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST("kstar_resolution"), femtoObsMC, femtoObs);
    }
  }

 protected:
  HistogramRegistry* mHistogramRegistry = nullptr;                                                                     ///< For QA output
  static constexpr std::string_view mFolderSuffix[4] = {"SameEvent/", "MixedEvent/", "SameEventMC/", "MixedEventMC/"}; ///< Folder naming for the output according to mEventType
  static constexpr femtoDreamContainer::Observable mFemtoObs = obs;                                                    ///< Femtoscopic observable to be computed (according to femtoDreamContainer::Observable)
  static constexpr int mEventType = eventType;                                                                         ///< Type of the event (same/mixed, according to femtoDreamContainer::EventType)
  float mMassOne = 0.f;                                                                                                ///< PDG mass of particle 1
  float mMassTwo = 0.f;                                                                                                ///< PDG mass of particle 2
};

} // namespace o2::analysis::femtoDream

#endif // PWGCF_FEMTODREAM_FEMTODREAMCONTAINER_H_
