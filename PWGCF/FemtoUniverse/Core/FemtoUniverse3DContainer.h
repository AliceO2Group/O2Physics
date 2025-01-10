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

/// \file FemtoUniverse3DContainer.h
/// \brief Definition of the FemtoUniverse3DContainer
/// \remark This file is inherited from ~/FemtoUniverse/Core/FemtoUniverseFemtoContainer.h on 10/01/2024
/// \author Pritam Chakraborty, WUT Warsaw, pritam.chakraborty@pw.edu.pl8

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSE3DCONTAINER_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSE3DCONTAINER_H_

#include <fairlogger/Logger.h>
#include <vector>
#include <string>

#include "Framework/HistogramRegistry.h"
#include "PWGCF/FemtoUniverse/Core/FemtoUniverseMath.h"
#include "Math/Vector4D.h"
#include "TMath.h"
#include "TDatabasePDG.h"

using namespace o2::framework;

namespace o2::analysis::femto_universe
{

namespace femto_universe3d_container
{
/// Femtoscopic observable to be computed
enum Observable { kstar ///< kstar
};

/// Type of the event processind
enum EventType { same, ///< Pair from same event
                 mixed ///< Pair from mixed event
};
}; // namespace femto_universe3d_container

/// \class FemtoUniverse3DContainer
/// \brief Container for all histogramming related to the correlation function. The two
/// particles of the pair are passed here, and the correlation function and QA histograms
/// are filled according to the specified observable
/// \tparam eventType Type of the event (same/mixed)
/// \tparam obs Observable to be computed (k*/Q_inv/...)
template <femto_universe3d_container::EventType eventType, femto_universe3d_container::Observable obs>
class FemtoUniverse3DContainer
{
 public:
  /// Destructor
  virtual ~FemtoUniverse3DContainer() = default;

  /// Initializes histograms for the task
  /// Called by init both in case of reconstructed data/ Monte Carlo, and for Monte Carlo Truth
  /// \tparam T type of the axis Object
  /// \param folderName Name of the directory in the output file (no suffix for reconstructed data/ Monte Carlo; "_MC" for Monte Carlo Truth)
  /// \param femtoObs1D Title of the femto observable kstar
  /// \param femtoObsKout Title of the femto observable kout
  /// \param femtoObsKside Title of the femto observable kside
  /// \param femtoObsKlong Title of the femto observable klong
  /// \param femtoObsAxis1D axis object for the femto observable kstar
  /// \param femtoObsAxisOut axis object for the femto observable kout
  /// \param femtoObsAxisSide axis object for the femto observable kside
  /// \param femtoObsAxisLong axis object for the femto observable klong
  /// \param multAxis axis object for the multiplicity axis
  /// \param kTAxis axis object for the kT axis
  /// \param mTAxis axis object for the mT axis
  /// \param use3dplots Flag to fill 3D plots
  /// \param isiden Identical or non-identical particle pair
  template <typename T>
  void initBase(std::string folderName, std::string femtoObs1D, std::string femtoObsKout, std::string femtoObsKside, std::string femtoObsKlong, T femtoObsAxis1D, T femtoObsAxisOut, T femtoObsAxisSide, T femtoObsAxisLong, T multAxis, T kTAxis, T mTAxis, T multAxis3D, T mTAxis3D, bool use3dplots, bool isiden)
  {
    mHistogramRegistry->add((folderName + "/relPairMom3D").c_str(), ("; " + femtoObsKout + "; " + femtoObsKside + "; " + femtoObsKlong).c_str(), kTH3F, {femtoObsAxisOut, femtoObsAxisSide, femtoObsAxisLong});
    mHistogramRegistry->add((folderName + "/relPairMomOut").c_str(), ("; " + femtoObsKout + "; Entries").c_str(), kTH1F, {femtoObsAxisOut});
    mHistogramRegistry->add((folderName + "/relPairMomSide").c_str(), ("; " + femtoObsKside + "; Entries").c_str(), kTH1F, {femtoObsAxisSide});
    mHistogramRegistry->add((folderName + "/relPairMomLong").c_str(), ("; " + femtoObsKlong + "; Entries").c_str(), kTH1F, {femtoObsAxisLong});
    mHistogramRegistry->add((folderName + "/relPairMom1D").c_str(), ("; " + femtoObs1D + "; Entries").c_str(), kTH1F, {femtoObsAxis1D});
    if (!isiden) {
      mHistogramRegistry->add((folderName + "/KStarOutP").c_str(), ("; " + femtoObs1D + "; Entries").c_str(), kTH1F, {femtoObsAxis1D});
      mHistogramRegistry->add((folderName + "/KStarSideP").c_str(), ("; " + femtoObs1D + "; Entries").c_str(), kTH1F, {femtoObsAxis1D});
      mHistogramRegistry->add((folderName + "/KStarLongP").c_str(), ("; " + femtoObs1D + "; Entries").c_str(), kTH1F, {femtoObsAxis1D});
      mHistogramRegistry->add((folderName + "/KStarOutN").c_str(), ("; " + femtoObs1D + "; Entries").c_str(), kTH1F, {femtoObsAxis1D});
      mHistogramRegistry->add((folderName + "/KStarSideN").c_str(), ("; " + femtoObs1D + "; Entries").c_str(), kTH1F, {femtoObsAxis1D});
      mHistogramRegistry->add((folderName + "/KStarLongN").c_str(), ("; " + femtoObs1D + "; Entries").c_str(), kTH1F, {femtoObsAxis1D});
    }
    mHistogramRegistry->add((folderName + "/relPairkT").c_str(), "; #it{k}_{T} (GeV/#it{c}); Entries", kTH1F, {kTAxis});
    mHistogramRegistry->add((folderName + "/relPairkstarkT").c_str(), ("; " + femtoObs1D + "; #it{k}_{T} (GeV/#it{c})").c_str(), kTH2F, {femtoObsAxis1D, kTAxis});
    mHistogramRegistry->add((folderName + "/relPairkstarmT").c_str(), ("; " + femtoObs1D + "; #it{m}_{T} (GeV/#it{c}^{2})").c_str(), kTH2F, {femtoObsAxis1D, mTAxis});
    mHistogramRegistry->add((folderName + "/relPairkstarMult").c_str(), ("; " + femtoObs1D + "; Multiplicity").c_str(), kTH2F, {femtoObsAxis1D, multAxis});
    mHistogramRegistry->add((folderName + "/kstarPtPart1").c_str(), ("; " + femtoObs1D + "; #it{p} _{T} Particle 1 (GeV/#it{c})").c_str(), kTH2F, {femtoObsAxis1D, {375, 0., 7.5}});
    mHistogramRegistry->add((folderName + "/kstarPtPart2").c_str(), ("; " + femtoObs1D + "; #it{p} _{T} Particle 2 (GeV/#it{c})").c_str(), kTH2F, {femtoObsAxis1D, {375, 0., 7.5}});
    mHistogramRegistry->add((folderName + "/MultPtPart1").c_str(), "; #it{p} _{T} Particle 1 (GeV/#it{c}); Multiplicity", kTH2F, {{375, 0., 7.5}, multAxis});
    mHistogramRegistry->add((folderName + "/MultPtPart2").c_str(), "; #it{p} _{T} Particle 2 (GeV/#it{c}); Multiplicity", kTH2F, {{375, 0., 7.5}, multAxis});
    mHistogramRegistry->add((folderName + "/PtPart1PtPart2").c_str(), "; #it{p} _{T} Particle 1 (GeV/#it{c}); #it{p} _{T} Particle 2 (GeV/#it{c})", kTH2F, {{375, 0., 7.5}, {375, 0., 7.5}});
    if (use3dplots) {
      mHistogramRegistry->add((folderName + "/relPairkstarmTMult").c_str(), ("; " + femtoObs1D + "; #it{m}_{T} (GeV/#it{c}^{2}); Multiplicity").c_str(), kTH3F, {femtoObsAxis1D, mTAxis3D, multAxis3D});
    }
  }

  /// Templated function to initialize the histograms for the task
  /// Always calls initBase to initialize the histograms for data/ Monte Carlo reconstructed
  /// \tparam T type of the configurable for the axis configuration
  /// \param registry Histogram registry to be passed
  /// \param kstarBins k* binning for the histograms
  /// \param multBins multiplicity binning for the histograms
  /// \param kTBins kT binning for the histograms
  /// \param mTBins mT binning for the histograms
  /// \param isMC add Monte Carlo truth histograms to the output file
  /// \param use3dplots Flag to fill 3D plots
  /// \param isiden Identical or non-identical particle pair
  template <typename T>
  void init(HistogramRegistry* registry, T& kstarBins, T& multBins, T& kTBins, T& mTBins, T& multBins3D, T& mTBins3D, bool /*isMC*/, bool use3dplots, bool isiden)
  {
    mHistogramRegistry = registry;
    std::string femtoObs1D, femtoObsKout, femtoObsKside, femtoObsKlong;

    if (isiden) {
      femtoObs1D = "#it{q} (GeV/#it{c})";
      femtoObsKout = "#it{q}_{out} (GeV/#it{c})";
      femtoObsKside = "#it{q}_{side} (GeV/#it{c})";
      femtoObsKlong = "#it{q}_{long} (GeV/#it{c})";
    } else {
      femtoObs1D = "#it{k*} (GeV/#it{c})";
      femtoObsKout = "#it{k*}_{out} (GeV/#it{c})";
      femtoObsKside = "#it{k*}_{side} (GeV/#it{c})";
      femtoObsKlong = "#it{k*}_{long} (GeV/#it{c})";
    }
    framework::AxisSpec femtoObsAxis1D = {kstarBins, femtoObs1D.c_str()};
    framework::AxisSpec femtoObsAxisOut = {kstarBins, femtoObsKout.c_str()};
    framework::AxisSpec femtoObsAxisSide = {kstarBins, femtoObsKside.c_str()};
    framework::AxisSpec femtoObsAxisLong = {kstarBins, femtoObsKlong.c_str()};

    std::vector<double> tmpVecMult = multBins;
    framework::AxisSpec multAxis = {tmpVecMult, "Multiplicity"};
    framework::AxisSpec kTAxis = {kTBins, "#it{k}_{T} (GeV/#it{c})"};
    framework::AxisSpec mTAxis = {mTBins, "#it{m}_{T} (GeV/#it{c}^{2})"};

    framework::AxisSpec multAxis3D = {multBins3D, "Multiplicity"};
    framework::AxisSpec mTAxis3D = {mTBins3D, "#it{m}_{T} (GeV/#it{c})"};

    std::string folderName = static_cast<std::string>(FolderSuffix[EventType]) + static_cast<std::string>(o2::aod::femtouniverse_mc_particle::MCTypeName[o2::aod::femtouniverse_mc_particle::MCType::kRecon]);

    initBase(folderName, femtoObs1D, femtoObsKout, femtoObsKside, femtoObsKlong, femtoObsAxis1D, femtoObsAxisOut, femtoObsAxisSide, femtoObsAxisLong, multAxis, kTAxis, mTAxis, multAxis3D, mTAxis3D, use3dplots, isiden);
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
  /// Called by setPair
  /// \tparam T type of the femtouniverseparticle
  /// \param part1 Particle one
  /// \param part2 Particle two
  /// \param mult Multiplicity of the event
  template <o2::aod::femtouniverse_mc_particle::MCType mc, typename T>
  void setPairBase(const float femtoObsKout, const float femtoObsKside, const float femtoObsKlong, const float femtoObs1D, const float kT, const float mT, T const& part1, T const& part2, const int mult, bool use3dplots, const float isiden)
  {
    mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/relPairMomOut"), femtoObsKout);
    mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/relPairMomSide"), femtoObsKside);
    mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/relPairMomLong"), femtoObsKlong);
    mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/relPairMom3D"), femtoObsKout, femtoObsKside, femtoObsKlong);
    mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/MultPtPart1"), part1.pt(), mult);
    mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/MultPtPart2"), part2.pt(), mult);
    mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/PtPart1PtPart2"), part1.pt(), part2.pt());

    if (isiden) {
      mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/relPairMom1D"), (2.0 * femtoObs1D));
      mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/relPairkT"), kT);
      mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/relPairkstarkT"), (2.0 * femtoObs1D), kT);
      mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/relPairkstarmT"), (2.0 * femtoObs1D), mT);
      mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/relPairkstarMult"), (2.0 * femtoObs1D), mult);
      mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/kstarPtPart1"), (2.0 * femtoObs1D), part1.pt());
      mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/kstarPtPart2"), (2.0 * femtoObs1D), part2.pt());
      if (use3dplots) {
        mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/relPairkstarmTMult"), (2.0 * femtoObs1D), mT, mult);
      }
    } else {
      mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/relPairMom1D"), femtoObs1D);
      mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/relPairkT"), kT);
      mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/relPairkstarkT"), femtoObs1D, kT);
      mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/relPairkstarmT"), femtoObs1D, mT);
      mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/relPairkstarMult"), femtoObs1D, mult);
      mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/kstarPtPart1"), femtoObs1D, part1.pt());
      mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/kstarPtPart2"), femtoObs1D, part2.pt());
      if (use3dplots) {
        mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[mc]) + HIST("/relPairkstarmTMult"), femtoObs1D, mT, mult);
      }
    }
  }

  /// Templated function to compute the necessary observables and fill the respective histograms
  /// Always calls setPairBase to compute the observables with reconstructed data
  /// \tparam T type of the femtouniverseparticle
  /// \param part1 Particle one
  /// \param part2 Particle two
  /// \param mult Multiplicity of the event
  /// \param use3dplots Flag to fill 3D plots
  /// \param isiden Choosing identical or non-identical pairs
  /// \param islcm Choosing LCMS or PRF
  template <bool isMC, typename T>
  void setPair(T const& part1, T const& part2, const int mult, bool use3dplots, bool isiden)
  {
    std::vector<double> f3d;
    const float kT = FemtoUniverseMath::getkT(part1, mMassOne, part2, mMassTwo);
    const float mT = FemtoUniverseMath::getmT(part1, mMassOne, part2, mMassTwo);

    f3d = FemtoUniverseMath::newpairfunc(part1, mMassOne, part2, mMassTwo, isiden);

    const float femtoObs1D = f3d[0];
    const float femtoObsKout = f3d[1];
    const float femtoObsKside = f3d[2];
    const float femtoObsKlong = f3d[3];

    if (mHistogramRegistry) {
      setPairBase<o2::aod::femtouniverse_mc_particle::MCType::kRecon>(femtoObsKout, femtoObsKside, femtoObsKlong, femtoObs1D, kT, mT, part1, part2, mult, use3dplots, isiden);
      if (!isiden) {
        if (femtoObsKout > 0.0) {
          mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[o2::aod::femtouniverse_mc_particle::MCType::kRecon]) + HIST("/KStarOutP"), femtoObs1D);
        } else {
          mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[o2::aod::femtouniverse_mc_particle::MCType::kRecon]) + HIST("/KStarOutN"), femtoObs1D);
        }
        if (femtoObsKside > 0.0) {
          mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[o2::aod::femtouniverse_mc_particle::MCType::kRecon]) + HIST("/KStarSideP"), femtoObs1D);
        } else {
          mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[o2::aod::femtouniverse_mc_particle::MCType::kRecon]) + HIST("/KStarSideN"), femtoObs1D);
        }
        if (femtoObsKlong > 0.0) {
          mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[o2::aod::femtouniverse_mc_particle::MCType::kRecon]) + HIST("/KStarLongP"), femtoObs1D);
        } else {
          mHistogramRegistry->fill(HIST(FolderSuffix[EventType]) + HIST(o2::aod::femtouniverse_mc_particle::MCTypeName[o2::aod::femtouniverse_mc_particle::MCType::kRecon]) + HIST("/KStarLongN"), femtoObs1D);
        }
      }
    }
  }

 protected:
  HistogramRegistry* mHistogramRegistry = nullptr;                                 ///< For QA output
  static constexpr std::string_view FolderSuffix[2] = {"SameEvent", "MixedEvent"}; ///< Folder naming for the output according to EventType
  static constexpr int EventType = eventType;                                      ///< Type of the event (same/mixed, according to femto_universe3d_container::EventType)
  float mMassOne = 0.f;                                                            ///< PDG mass of particle 1
  float mMassTwo = 0.f;                                                            ///< PDG mass of particle 2
  int mPDGOne = 0;                                                                 ///< PDG code of particle 1
  int mPDGTwo = 0;                                                                 ///< PDG code of particle 2
};

} // namespace o2::analysis::femto_universe

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSE3DCONTAINER_H_
