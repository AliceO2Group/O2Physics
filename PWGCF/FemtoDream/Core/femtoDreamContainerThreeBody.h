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
/// \author Laura Serksnyte, TU München, laura.serksnyte@tum.de

#ifndef PWGCF_FEMTODREAM_CORE_FEMTODREAMCONTAINERTHREEBODY_H_
#define PWGCF_FEMTODREAM_CORE_FEMTODREAMCONTAINERTHREEBODY_H_

#include <fairlogger/Logger.h>
#include <vector>
#include <string>

#include "Framework/HistogramRegistry.h"
#include "PWGCF/FemtoDream/Core/femtoDreamMath.h"

#include "Math/Vector4D.h"
#include "TMath.h"
#include "TDatabasePDG.h"

using namespace o2::framework;

namespace o2::analysis::femtoDream
{

namespace femtoDreamContainerThreeBody
{
/// Femtoscopic observable to be computed
enum Observable { Q3 ///< Q3
};

/// Type of the event processind
enum EventType { same, ///< Pair from same event
                 mixed ///< Pair from mixed event
};
}; // namespace femtoDreamContainerThreeBody

/// \class FemtoDreamContainerThreeBody
/// \brief Container for all histogramming related to the correlation function. The two
/// particles of the pair are passed here, and the correlation function and QA histograms
/// are filled according to the specified observable
/// \tparam eventType Type of the event (same/mixed)
/// \tparam obs Observable to be computed (Q_3/...)
template <femtoDreamContainerThreeBody::EventType eventType, femtoDreamContainerThreeBody::Observable obs>
class FemtoDreamContainerThreeBody
{
 public:
  /// Destructor
  virtual ~FemtoDreamContainerThreeBody() = default;

  /// Initializes histograms for the task in case of three-body femtoscopy
  /// Called by init both in case of reconstructed data/ Monte Carlo, and for Monte Carlo Truth
  /// \tparam T type of the axis Object
  /// \param folderName Name of the directory in the output file (no suffix for reconstructed data/ Monte Carlo; "_MC" for Monte Carlo Truth)
  /// \param femtoObs Title of the femto observable axis
  /// \param femtoObsAxis axis object for the femto observable axis
  /// \param multAxis axis object for the multiplicity axis
  template <typename T>
  void init_base(std::string folderName, std::string femtoObs, T femtoObsAxis, T multAxis)
  {

    mHistogramRegistry->add((folderName + "/relTripletDist").c_str(), ("; " + femtoObs + "; Entries").c_str(), kTH1F, {femtoObsAxis});
    mHistogramRegistry->add((folderName + "/relTripletQ3Mult").c_str(), ("; " + femtoObs + "; Multiplicity").c_str(), kTH2F, {femtoObsAxis, multAxis});
    mHistogramRegistry->add((folderName + "/mT1").c_str(), ";mT; Q3", kTH2F, {{1000, 0, 25}, femtoObsAxis});
    mHistogramRegistry->add((folderName + "/mT2").c_str(), ";mT; Q3", kTH2F, {{1000, 0, 25}, femtoObsAxis});
    mHistogramRegistry->add((folderName + "/mT3").c_str(), ";mT; Q3", kTH2F, {{1000, 0, 25}, femtoObsAxis});
    mHistogramRegistry->add((folderName + "/mTAverage").c_str(), ";mT; Q3", kTH2F, {{1000, 0, 25}, femtoObsAxis});
  }

  /// Initializes specialized Monte Carlo truth histograms for the task in case of three-body femtoscopy
  /// internal function called by init only in case of Monte Carlo truth
  /// \tparam T type of the xxis Object
  /// \param folderName Name of the directory in the output file (no suffix for reconstructed data/ Monte Carlo; "_MC" for Monte Carlo Truth)
  /// \param femtoObsAxis axis object for the femto observable axis
  template <typename T>
  void init_MC(std::string folderName, std::string femtoObs, T femtoObsAxis, T multAxis)
  {
    mHistogramRegistry->add((folderName + "/relTripletDist_ReconNoFake").c_str(), ("; " + femtoObs + "; Entries").c_str(), kTH1F, {femtoObsAxis});
    mHistogramRegistry->add((folderName + "/relTripletQ3Mult_ReconNoFake").c_str(), ("; " + femtoObs + "; Multiplicity").c_str(), kTH2F, {femtoObsAxis, multAxis});
    mHistogramRegistry->add((folderName + "/hNoMCtruthTripletCounter").c_str(), "; Counter; Entries", kTH1I, {{1, 0, 1}});
    mHistogramRegistry->add((folderName + "/hFakeTripletCounter").c_str(), "; Counter; Entries", kTH1I, {{1, 0, 1}});
    mHistogramRegistry->add((folderName + "/Q3_resolution").c_str(), "; #it{Q}_{3} reconstructed (GeV/#it{c}); #it{Q}_{3} truth (GeV/#it{c})", kTH2F, {femtoObsAxis, femtoObsAxis});
  }

  /// Templated function to initialize the histograms for the task in case of three-body femtoscopy
  /// Always calls init_base to initialize the histograms for data/ Monte Carlo reconstructed
  /// In case of Monte Carlo, calls init_base again for Monte Carlo truth and the specialized function init_MC for additional histogramms
  /// \tparam T type of the configurable for the axis configuration
  /// \param registry Histogram registry to be passed
  /// \param Q3Bins Q3 binning for the histograms
  /// \param multBins multiplicity binning for the histograms
  /// \param isMC add Monte Carlo truth histograms to the output file
  template <typename T>
  void init(HistogramRegistry* registry, T& Q3Bins, T& multBins, bool isMC)
  {
    mHistogramRegistry = registry;
    std::string femtoObs;
    if constexpr (mFemtoObs == femtoDreamContainerThreeBody::Observable::Q3) {
      femtoObs = "#it{Q}_{3} (GeV/#it{c})";
    }
    std::vector<double> tmpVecMult = multBins;
    framework::AxisSpec multAxis = {tmpVecMult, "Multiplicity"};
    framework::AxisSpec femtoObsAxis = {Q3Bins, femtoObs.c_str()};

    std::string folderName = static_cast<std::string>(mFolderSuffix[mEventType]) + static_cast<std::string>(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kRecon]);

    init_base(folderName, femtoObs, femtoObsAxis, multAxis);
    if (isMC) {
      folderName = static_cast<std::string>(mFolderSuffix[mEventType]) + static_cast<std::string>(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kTruth]);
      init_base(folderName, femtoObs, femtoObsAxis, multAxis);
      init_MC(folderName, femtoObs, femtoObsAxis, multAxis);
    }
  }

  /// Set the PDG codes of the three particles involved
  /// \param pdg1 PDG code of particle one
  /// \param pdg2 PDG code of particle two
  /// \param pdg3 PDG code of particle three
  void setPDGCodes(const int pdg1, const int pdg2, const int pdg3)
  {
    mMassOne = TDatabasePDG::Instance()->GetParticle(pdg1)->Mass();
    mMassTwo = TDatabasePDG::Instance()->GetParticle(pdg2)->Mass();
    mMassThree = TDatabasePDG::Instance()->GetParticle(pdg3)->Mass();
    mPDGOne = pdg1;
    mPDGTwo = pdg2;
    mPDGThree = pdg3;
  }

  /// Pass a triplet to the container and compute all the relevant observables
  /// Called by setTriplet both in case of data/ and Monte Carlo reconstructed and for Monte Carlo truth
  /// \tparam T type of the femtodreamparticle
  /// \param part1 Particle one
  /// \param part2 Particle two
  /// \param part3 Particle three
  /// \param mult Multiplicity of the event
  template <o2::aod::femtodreamMCparticle::MCType mc, typename T>
  void setTriplet_base(const float femtoObs, T const& /*part1*/, T const& /*part2*/, T const& /*part3*/, const int mult)
  {
    // FILL MAIN Q3 info
    mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/relTripletDist"), femtoObs);
    mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/relTripletQ3Mult"), femtoObs, mult);
  }

  /// Called by setTriplet only in case of Monte Carlo truth
  /// Fills MC truth specific histogramms:
  /// - Q3 distribution plots with RECONSTRUCTED information but ONLY for non-fake candidates; needed for purity calculations of tracks
  /// - Q3 resolution matrix
  /// Note: Standard histograms with MC truth information are filled with the setPair_base function
  /// \param part1 Particle one
  /// \param part2 Particle two
  /// \param part3 Particle three
  /// \param mult Multiplicity of the event
  void setTriplet_MC(const float femtoObsMC, const float femtoObs, const int mult)
  {
    if (mHistogramRegistry) {
      // Fill the Q3 distributions with the reconstructed information but only for particles with the right PDG code
      mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kTruth]) + HIST("/relTripletDist_ReconNoFake"), femtoObs);
      mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kTruth]) + HIST("/relTripletQ3Mult_ReconNoFake"), femtoObs, mult);
      mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kTruth]) + HIST("/Q3_resolution"), femtoObsMC, femtoObs);
    }
  }

  /// Templated function to handle data/ Monte Carlo reconstructed and Monte Carlo truth
  /// Always calls setTriplet_base to compute the observables with reconstructed data
  /// In case of Monte Carlo, calls setTriplet_base with MC info and specialized function setTriplet_MC for additional histogramms
  /// \tparam T type of the femtodreamparticle
  /// \param part1 Particle one
  /// \param part2 Particle two
  /// \param part3 Particle three
  /// \param mult Multiplicity of the event
  template <bool isMC, typename T>
  void setTriplet(T const& part1, T const& part2, T const& part3, const int mult, const float femtoObs)
  {
    float femtoObsMC;
    float mT1 = FemtoDreamMath::getmT(part1, mMassOne, part2, mMassTwo);
    float mT2 = FemtoDreamMath::getmT(part2, mMassTwo, part3, mMassThree);
    float mT3 = FemtoDreamMath::getmT(part1, mMassOne, part3, mMassThree);
    float mTAverage = (mT1 + mT2 + mT3) / 3.;
    // FILL mT info
    mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kRecon]) + HIST("/mT1"), mT1, femtoObs);
    mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kRecon]) + HIST("/mT2"), mT2, femtoObs);
    mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kRecon]) + HIST("/mT3"), mT3, femtoObs);
    mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kRecon]) + HIST("/mTAverage"), mTAverage, femtoObs);

    if (mHistogramRegistry) {
      setTriplet_base<o2::aod::femtodreamMCparticle::MCType::kRecon>(femtoObs, part1, part2, part3, mult);

      if constexpr (isMC) {
        if (part1.has_fdMCParticle() && part2.has_fdMCParticle() && part3.has_fdMCParticle()) {
          // calculate the femto observable with MC truth information
          if constexpr (mFemtoObs == femtoDreamContainerThreeBody::Observable::Q3) {
            femtoObsMC = FemtoDreamMath::getQ3(part1.fdMCParticle(), mMassOne, part2.fdMCParticle(), mMassTwo, part3.fdMCParticle(), mMassThree);
          }

          if (abs(part1.fdMCParticle().pdgMCTruth()) == mPDGOne && abs(part2.fdMCParticle().pdgMCTruth()) == mPDGTwo && abs(part3.fdMCParticle().pdgMCTruth()) == mPDGThree) { // Note: all triplet-histogramms are filled with MC truth information ONLY in case of non-fake candidates
            setTriplet_base<o2::aod::femtodreamMCparticle::MCType::kTruth>(femtoObsMC, part1.fdMCParticle(), part2.fdMCParticle(), part3.fdMCParticle(), mult);
            setTriplet_MC(femtoObsMC, femtoObs, mult);
          } else {
            mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kTruth]) + HIST("/hFakeTripletCounter"), 0);
          }

        } else {
          mHistogramRegistry->fill(HIST(mFolderSuffix[mEventType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kTruth]) + HIST("/hNoMCtruthTripletCounter"), 0);
        }
      }
    }
  }

 protected:
  HistogramRegistry* mHistogramRegistry = nullptr;                                  ///< For QA output
  static constexpr std::string_view mFolderSuffix[2] = {"SameEvent", "MixedEvent"}; ///< Folder naming for the output according to mEventType
  static constexpr femtoDreamContainerThreeBody::Observable mFemtoObs = obs;        ///< Femtoscopic observable to be computed (according to femtoDreamContainerThreeBody::Observable)
  static constexpr int mEventType = eventType;                                      ///< Type of the event (same/mixed, according to femtoDreamContainerThreeBody::EventType)
  float mMassOne = 0.f;                                                             ///< PDG mass of particle 1
  float mMassTwo = 0.f;                                                             ///< PDG mass of particle 1
  float mMassThree = 0.f;                                                           ///< PDG mass of particle 3; if relevant
  int mPDGOne = 0;                                                                  ///< PDG code of particle 1
  int mPDGTwo = 0;                                                                  ///< PDG code of particle 2
  int mPDGThree = 0;                                                                ///< PDG code of particle 3; if relevant
};

} // namespace o2::analysis::femtoDream

#endif // PWGCF_FEMTODREAM_CORE_FEMTODREAMCONTAINERTHREEBODY_H_
