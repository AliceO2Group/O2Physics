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

/// \file FemtoFlowFemtoContainer.h
/// \brief Definition of the FemtoFlowFemtoContainer
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Valentina Mantovani Sarti, valentina.mantovani-sarti@tum.de
/// \author Georgios Mantzaridis, TU München, georgios.mantzaridis@tum.de
/// \author Anton Riedel, TU München, anton.riedel@tum.de
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch

#ifndef PWGCF_FEMTO_CORE_FEMTOFLOWFEMTOCONTAINER_H_
#define PWGCF_FEMTO_CORE_FEMTOFLOWFEMTOCONTAINER_H_

#include "PWGCF/Femto/Core/FemtoFlowMath.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "Math/Vector4D.h"
#include "TMath.h"

#include <fairlogger/Logger.h>

#include <string>
#include <vector>

namespace o2::analysis::femto_flow
{

namespace femto_flow_femto_container
{
/// Femtoscopic observable to be computed
enum Observable { kstar ///< kstar
};

/// Type of the event processind
enum EventType { same, ///< Pair from same event
                 mixed ///< Pair from mixed event
};
}; // namespace femto_flow_femto_container

/// \class FemtoFlowFemtoContainer
/// \brief Container for all histogramming related to the correlation function. The two
/// particles of the pair are passed here, and the correlation function and QA histograms
/// are filled according to the specified observable
/// \tparam eventType Type of the event (same/mixed)
/// \tparam obs Observable to be computed (k*/Q_inv/...)
template <femto_flow_femto_container::EventType eventType, femto_flow_femto_container::Observable obs>
class FemtoFlowFemtoContainer
{
 public:
  /// Destructor
  virtual ~FemtoFlowFemtoContainer() = default;

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
  void initBase(std::string folderName, std::string femtoObs, T femtoObsAxis, T multAxis, T kTAxis, T mTAxis, T multAxis3D, T mTAxis3D, bool use3dplots, bool useqnDivide)
  {
    using namespace o2::framework;

    mHistogramRegistry->add((folderName + "/relPairDist").c_str(), ("; " + femtoObs + "; Entries").c_str(), kTH1F, {femtoObsAxis});
    mHistogramRegistry->add((folderName + "/relPairkT").c_str(), "; #it{k}_{T} (GeV/#it{c}); Entries", kTH1F, {kTAxis});
    mHistogramRegistry->add((folderName + "/relPairkstarkT").c_str(), ("; " + femtoObs + "; #it{k}_{T} (GeV/#it{c})").c_str(), kTH2F, {femtoObsAxis, kTAxis});
    mHistogramRegistry->add((folderName + "/relPairkstarmT").c_str(), ("; " + femtoObs + "; #it{m}_{T} (GeV/#it{c}^{2})").c_str(), kTH2F, {femtoObsAxis, mTAxis});
    mHistogramRegistry->add((folderName + "/relPairkstarMult").c_str(), ("; " + femtoObs + "; Multiplicity").c_str(), kTH2F, {femtoObsAxis, multAxis});
    mHistogramRegistry->add((folderName + "/kstarPtPart1").c_str(), ("; " + femtoObs + "; #it{p} _{T} Particle 1 (GeV/#it{c})").c_str(), kTH2F, {femtoObsAxis, {375, 0., 7.5}});
    mHistogramRegistry->add((folderName + "/kstarPtPart2").c_str(), ("; " + femtoObs + "; #it{p} _{T} Particle 2 (GeV/#it{c})").c_str(), kTH2F, {femtoObsAxis, {375, 0., 7.5}});
    mHistogramRegistry->add((folderName + "/MultPtPart1").c_str(), "; #it{p} _{T} Particle 1 (GeV/#it{c}); Multiplicity", kTH2F, {{375, 0., 7.5}, multAxis});
    mHistogramRegistry->add((folderName + "/MultPtPart2").c_str(), "; #it{p} _{T} Particle 2 (GeV/#it{c}); Multiplicity", kTH2F, {{375, 0., 7.5}, multAxis});
    mHistogramRegistry->add((folderName + "/PtPart1PtPart2").c_str(), "; #it{p} _{T} Particle 1 (GeV/#it{c}); #it{p} _{T} Particle 2 (GeV/#it{c})", kTH2F, {{375, 0., 7.5}, {375, 0., 7.5}});
    if (use3dplots) {
      mHistogramRegistry->add((folderName + "/relPairkstarmTMult").c_str(), ("; " + femtoObs + "; #it{m}_{T} (GeV/#it{c}^{2}); Multiplicity").c_str(), kTH3F, {femtoObsAxis, mTAxis3D, multAxis3D});
    }
    if (useqnDivide) {
      for (int iqn(0); iqn < numqnBins; ++iqn) {
        mHistogramRegistry->add((folderName + std::to_string(iqn) + "/relPairDist").c_str(), ("; " + femtoObs + "; Entries").c_str(), kTH1F, {femtoObsAxis});
        mHistogramRegistry->add((folderName + std::to_string(iqn) + "/relPairkstarmT").c_str(), ("; " + femtoObs + "; #it{m}_{T} (GeV/#it{c}^{2})").c_str(), kTH2F, {femtoObsAxis, mTAxis});
        mHistogramRegistry->add((folderName + std::to_string(iqn) + "/relPairkstarMult").c_str(), ("; " + femtoObs + "; Multiplicity").c_str(), kTH2F, {femtoObsAxis, multAxis});
      }
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
    using namespace o2::framework;

    mHistogramRegistry->add((folderName + "/relPairDist_ReconNoFake").c_str(), ("; " + femtoObs + "; Entries").c_str(), kTH1F, {femtoObsAxis});
    mHistogramRegistry->add((folderName + "/relPairkstarmT_ReconNoFake").c_str(), ("; " + femtoObs + "; #it{m}_{T} (GeV/#it{c}^{2})").c_str(), kTH2F, {femtoObsAxis, mTAxis});
    mHistogramRegistry->add((folderName + "/relPairkstarMult_ReconNoFake").c_str(), ("; " + femtoObs + "; Multiplicity").c_str(), kTH2F, {femtoObsAxis, multAxis});
    mHistogramRegistry->add((folderName + "/hNoMCtruthPairsCounter").c_str(), "; Counter; Entries", kTH1I, {{1, 0, 1}});
    mHistogramRegistry->add((folderName + "/hFakePairsCounter").c_str(), "; Counter; Entries", kTH1I, {{1, 0, 1}});
    mHistogramRegistry->add((folderName + "/kstar_resolution").c_str(), "; #it{k} _{T} reconstructed (GeV/#it{c}); #it{k} _{T} truth (GeV/#it{c})", kTH2F, {femtoObsAxis, femtoObsAxis});
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
  void init(o2::framework::HistogramRegistry* registry, T& kstarBins, T& multBins, T& kTBins, T& mTBins, T& multBins3D, T& mTBins3D, bool isMC, bool use3dplots, bool useqnDivide)
  {
    using namespace o2::framework;

    mHistogramRegistry = registry;
    std::string femtoObs;
    if constexpr (FemtoObs == femto_flow_femto_container::Observable::kstar) {
      femtoObs = "#it{k*} (GeV/#it{c})";
    }
    std::vector<double> tmpVecMult = multBins;
    o2::framework::AxisSpec multAxis = {tmpVecMult, "Multiplicity"};
    o2::framework::AxisSpec femtoObsAxis = {kstarBins, femtoObs.c_str()};
    o2::framework::AxisSpec kTAxis = {kTBins, "#it{k}_{T} (GeV/#it{c})"};
    o2::framework::AxisSpec mTAxis = {mTBins, "#it{m}_{T} (GeV/#it{c}^{2})"};

    o2::framework::AxisSpec multAxis3D = {multBins3D, "Multiplicity"};
    o2::framework::AxisSpec mTAxis3D = {mTBins3D, "#it{m}_{T} (GeV/#it{c})"};

    std::string folderName = static_cast<std::string>(FolderSuffix[kEventType]) + static_cast<std::string>(o2::aod::femtoflow_mc_particle::MCTypeName[o2::aod::femtoflow_mc_particle::MCType::kRecon]);

    initBase(folderName, femtoObs, femtoObsAxis, multAxis, kTAxis, mTAxis, multAxis3D, mTAxis3D, use3dplots, useqnDivide);
    if (isMC) {
      folderName = static_cast<std::string>(FolderSuffix[kEventType]) + static_cast<std::string>(o2::aod::femtoflow_mc_particle::MCTypeName[o2::aod::femtoflow_mc_particle::MCType::kTruth]);
      initBase(folderName, femtoObs, femtoObsAxis, multAxis, kTAxis, mTAxis, multAxis3D, mTAxis3D, use3dplots, useqnDivide);
      initMC(folderName, femtoObs, femtoObsAxis, multAxis, mTAxis);
    }
  }

  /// Set the PDG codes of the two particles involved
  /// \param pdg1 PDG code of particle one
  /// \param pdg2 PDG code of particle two
  void setPDGCodesMass(const int pdg1, const int pdg2, const double mass1, const int mass2)
  {
    kMassOne = mass1;
    kMassTwo = mass2;
    kPDGOne = pdg1;
    kPDGTwo = pdg2;
  }

  /// Pass a pair to the container and compute all the relevant observables
  /// Called by setPair both in case of data/ and Monte Carlo reconstructed and for Monte Carlo truth
  /// \tparam T type of the femtoflowparticle
  /// \param part1 Particle one
  /// \param part2 Particle two
  /// \param mult Multiplicity of the event
  template <o2::aod::femtoflow_mc_particle::MCType mc, typename T>
  void setPairBase(const float femtoObs, const float mT, T const& part1, T const& part2, const int mult, bool use3dplots, bool useqnDivide, int mybinNum)
  {
    using namespace o2::framework;

    const float kT = FemtoFlowMath::getkT(part1, kMassOne, part2, kMassTwo);

    mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("/relPairDist"), femtoObs);
    mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("/relPairkT"), kT);
    mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("/relPairkstarkT"), femtoObs, kT);
    mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("/relPairkstarmT"), femtoObs, mT);
    mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("/relPairkstarMult"), femtoObs, mult);
    mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("/kstarPtPart1"), femtoObs, part1.pt());
    mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("/kstarPtPart2"), femtoObs, part2.pt());
    mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("/MultPtPart1"), part1.pt(), mult);
    mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("/MultPtPart2"), part2.pt(), mult);
    mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("/PtPart1PtPart2"), part1.pt(), part2.pt());
    if (use3dplots) {
      mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("/relPairkstarmTMult"), femtoObs, mT, mult);
    }
    if (useqnDivide && mybinNum >= 0 && mybinNum < numqnBins) {
      switch (mybinNum) {
        case 0:
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("0") + HIST("/relPairDist"), femtoObs);
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("0") + HIST("/relPairkstarmT"), femtoObs, mT);
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("0") + HIST("/relPairkstarMult"), femtoObs, mult);
          break;
        case 1:
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("1") + HIST("/relPairDist"), femtoObs);
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("1") + HIST("/relPairkstarmT"), femtoObs, mT);
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("1") + HIST("/relPairkstarMult"), femtoObs, mult);
          break;
        case 2:
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("2") + HIST("/relPairDist"), femtoObs);
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("2") + HIST("/relPairkstarmT"), femtoObs, mT);
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("2") + HIST("/relPairkstarMult"), femtoObs, mult);
          break;
        case 3:
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("3") + HIST("/relPairDist"), femtoObs);
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("3") + HIST("/relPairkstarmT"), femtoObs, mT);
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("3") + HIST("/relPairkstarMult"), femtoObs, mult);
          break;
        case 4:
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("4") + HIST("/relPairDist"), femtoObs);
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("4") + HIST("/relPairkstarmT"), femtoObs, mT);
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("4") + HIST("/relPairkstarMult"), femtoObs, mult);
          break;
        case 5:
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("5") + HIST("/relPairDist"), femtoObs);
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("5") + HIST("/relPairkstarmT"), femtoObs, mT);
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("5") + HIST("/relPairkstarMult"), femtoObs, mult);
          break;
        case 6:
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("6") + HIST("/relPairDist"), femtoObs);
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("6") + HIST("/relPairkstarmT"), femtoObs, mT);
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("6") + HIST("/relPairkstarMult"), femtoObs, mult);
          break;
        case 7:
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("7") + HIST("/relPairDist"), femtoObs);
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("7") + HIST("/relPairkstarmT"), femtoObs, mT);
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("7") + HIST("/relPairkstarMult"), femtoObs, mult);
          break;
        case 8:
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("8") + HIST("/relPairDist"), femtoObs);
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("8") + HIST("/relPairkstarmT"), femtoObs, mT);
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("8") + HIST("/relPairkstarMult"), femtoObs, mult);
          break;
        case 9:
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("9") + HIST("/relPairDist"), femtoObs);
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("9") + HIST("/relPairkstarmT"), femtoObs, mT);
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[mc]) + HIST("9") + HIST("/relPairkstarMult"), femtoObs, mult);
          break;
        default:
          return; // invalid qn bin
      }
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
    using namespace o2::framework;

    if (mHistogramRegistry) {
      // Fill the kstar distributions with the reconstructed information but only for particles with the right PDG code
      mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[o2::aod::femtoflow_mc_particle::MCType::kTruth]) + HIST("/relPairDist_ReconNoFake"), femtoObs);
      mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[o2::aod::femtoflow_mc_particle::MCType::kTruth]) + HIST("/relPairkstarmT_ReconNoFake"), femtoObs, mT);
      mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[o2::aod::femtoflow_mc_particle::MCType::kTruth]) + HIST("/relPairkstarMult_ReconNoFake"), femtoObs, mult);

      mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[o2::aod::femtoflow_mc_particle::MCType::kTruth]) + HIST("/kstar_resolution"), femtoObsMC, femtoObs);
    }
  }

  /// Templated function to handle data/ Monte Carlo reconstructed and Monte Carlo truth
  /// Always calls setPairBase to compute the observables with reconstructed data
  /// In case of Monte Carlo, calls setPairBase with MC info and specialized function setPairMC for additional histogramms
  /// \tparam T type of the femtoflowparticle
  /// \param part1 Particle one
  /// \param part2 Particle two
  /// \param mult Multiplicity of the event
  template <bool isMC, typename T>
  void setPair(T const& part1, T const& part2, const int mult, bool use3dplots, bool useqnDivide, int mybinNum)
  {
    using namespace o2::framework;

    float femtoObs, femtoObsMC;
    // Calculate femto observable and the mT with reconstructed information
    if constexpr (FemtoObs == femto_flow_femto_container::Observable::kstar) {
      femtoObs = FemtoFlowMath::getkstar(part1, kMassOne, part2, kMassTwo);
    }
    const float mT = FemtoFlowMath::getmT(part1, kMassOne, part2, kMassTwo);

    if (mHistogramRegistry) {
      setPairBase<o2::aod::femtoflow_mc_particle::MCType::kRecon>(femtoObs, mT, part1, part2, mult, use3dplots, useqnDivide, mybinNum);

      if constexpr (isMC) {
        if (part1.has_fDMCParticle() && part2.has_fDMCParticle()) {
          // calculate the femto observable and the mT with MC truth information
          if constexpr (FemtoObs == femto_flow_femto_container::Observable::kstar) {
            femtoObsMC = FemtoFlowMath::getkstar(part1.fDMCParticle(), kMassOne, part2.fDMCParticle(), kMassTwo);
          }
          const float mTMC = FemtoFlowMath::getmT(part1.fDMCParticle(), kMassOne, part2.fDMCParticle(), kMassTwo);

          if (std::abs(part1.fDMCParticle().pdgMCTruth()) == std::abs(kPDGOne) && std::abs(part2.fDMCParticle().pdgMCTruth()) == std::abs(kPDGTwo)) { // Note: all pair-histogramms are filled with MC truth information ONLY in case of non-fake candidates
            setPairBase<o2::aod::femtoflow_mc_particle::MCType::kTruth>(femtoObsMC, mTMC, part1.fDMCParticle(), part2.fDMCParticle(), mult, use3dplots, useqnDivide, mybinNum);
            setPairMC(femtoObsMC, femtoObs, mT, mult);
          } else {
            mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[o2::aod::femtoflow_mc_particle::MCType::kTruth]) + HIST("/hFakePairsCounter"), 0);
          }

        } else {
          mHistogramRegistry->fill(HIST(FolderSuffix[kEventType]) + HIST(o2::aod::femtoflow_mc_particle::MCTypeName[o2::aod::femtoflow_mc_particle::MCType::kTruth]) + HIST("/hNoMCtruthPairsCounter"), 0);
        }
      }
    }
  }

 protected:
  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;                   ///< For QA output
  static constexpr std::string_view FolderSuffix[2] = {"SameEvent", "MixedEvent"}; ///< Folder naming for the output according to kEventType
  static constexpr femto_flow_femto_container::Observable FemtoObs = obs;          ///< Femtoscopic observable to be computed (according to femto_flow_femto_container::Observable)
  static constexpr int kEventType = eventType;                                     ///< Type of the event (same/mixed, according to femto_flow_femto_container::EventType)
  float kMassOne = 0.f;                                                            ///< PDG mass of particle 1
  float kMassTwo = 0.f;                                                            ///< PDG mass of particle 2
  int kPDGOne = 0;                                                                 ///< PDG code of particle 1
  int kPDGTwo = 0;
  int numqnBins = 10; ///< Max num of devided qn bins
};

} // namespace o2::analysis::femto_flow

#endif // PWGCF_FEMTO_CORE_FEMTOFLOWFEMTOCONTAINER_H_
