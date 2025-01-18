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

/// \file FemtoDreamParticleHisto.h
/// \brief FemtoDreamParticleHisto - Histogram class for tracks, V0s and cascades
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Georgios Mantzaridis, TU München, georgios.mantzaridis@tum.de
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTODREAM_CORE_FEMTODREAMPARTICLEHISTO_H_
#define PWGCF_FEMTODREAM_CORE_FEMTODREAMPARTICLEHISTO_H_

#include <TMath.h>
#include <string>
#include "PWGCF/DataModel/FemtoDerived.h"
#include "Framework/HistogramRegistry.h"
#include "CommonConstants/PhysicsConstants.h"

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

  /// Initializes particle histograms
  /// Called by init both in case of reconstructed data/ Monte Carlo, and for Monte Carlo Truth
  /// \tparam T type of the axis Object
  /// \tparam mc enum object to get the suffix ("" for data/ Monte Cartlo reconstructed, "_MC" for Monte Carlo truth) for the folder in the output file
  /// \param folderName base path of the directory in the output file, in which to store the histograms
  /// \param tempFitVarAxisTitle  Title of the axis of the tempFitVar (DCA_xy in case of tracks, CPA in case of V0s, etc.)
  /// \param tempFitVarpTAxis axis object for the pT axis in the pT vs. tempFitVar plots
  /// \param tempFitVarAxis axis object for the tempFitVar axis
  template <o2::aod::femtodreamMCparticle::MCType mc, typename T>
  void init_base(std::string folderName, std::string tempFitVarAxisTitle, T& pTAxis, T& tempFitVarAxis, T& InvMassAxis, T& /*multAxis*/)
  {
    std::string folderSuffix = static_cast<std::string>(o2::aod::femtodreamMCparticle::MCTypeName[mc]).c_str();
    /// Histograms of the kinematic properties
    mHistogramRegistry->add((folderName + folderSuffix + "/hPt").c_str(), "; #it{p}_{T} (GeV/#it{c}); Entries", kTH1F, {{240, 0, 6}});
    mHistogramRegistry->add((folderName + folderSuffix + "/hEta").c_str(), "; #eta; Entries", kTH1F, {{200, -1.5, 1.5}});
    mHistogramRegistry->add((folderName + folderSuffix + "/hPhi").c_str(), "; #phi; Entries", kTH1F, {{200, 0, TMath::TwoPi()}});

    /// particle specific histogramms for the TempFitVar column in FemtoDreamParticles
    if constexpr (o2::aod::femtodreamMCparticle::MCType::kRecon == mc) {
      mHistogramRegistry->add((folderName + folderSuffix + static_cast<std::string>(o2::aod::femtodreamparticle::TempFitVarName[mParticleType])).c_str(), ("; #it{p}_{T} (GeV/#it{c}); " + tempFitVarAxisTitle).c_str(), kTH2F, {{pTAxis}, {tempFitVarAxis}});
    }
    if constexpr ((mParticleType == o2::aod::femtodreamparticle::ParticleType::kV0 || mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascadeV0) && mc == o2::aod::femtodreamMCparticle::MCType::kRecon) {
      mHistogramRegistry->add((folderName + folderSuffix + "/hInvMassLambda").c_str(), "; M_{#Lambda}; Entries", kTH1F, {InvMassAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/hpTInvMassLambda").c_str(), "; p_{T} (GeV/#it{c{}); M_{#Lambda}", kTH2F, {pTAxis, InvMassAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/hInvMassAntiLambda").c_str(), "; M_{#bar{#Lambda}}; Entries", kTH1F, {InvMassAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/hpTInvMassAntiLambda").c_str(), "; M_{#bar{#Lambda}}; Entries", kTH2F, {pTAxis, InvMassAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/hInvMassLambdaAntiLambda").c_str(), "; M_{#Lambda}; M_{#bar{#Lambda}}", kTH2F, {InvMassAxis, InvMassAxis});
    }
    if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascade) {
      mHistogramRegistry->add((folderName + folderSuffix + "/hInvMassCascade").c_str(), "; M_{Cascade}; Entries", kTH1F, {InvMassAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/hpTInvMassCascade").c_str(), "; p_{T} (GeV/#it{c{}); M_{Cascade}", kTH2F, {pTAxis, InvMassAxis});
    }
  }

  // comment
  template <o2::aod::femtodreamMCparticle::MCType mc, typename T>
  void init_debug(std::string folderName, T& multAxis, T& multPercentileAxis, T& pTAxis, T& etaAxis, T& phiAxis, T& tempFitVarAxis, T& dcazAxis, T& NsigmaTPCAxis, T& NsigmaTOFAxis, T& NsigmaTPCTOFAxis, T& NsigmaITSAxis, bool correlatedPlots)
  {

    std::string folderSuffix = static_cast<std::string>(o2::aod::femtodreamMCparticle::MCTypeName[mc]).c_str();

    mHistogramRegistry->add((folderName + folderSuffix + "/hMomentumVsEta").c_str(), "; #it{p} (GeV/#it{c}); #eta", kTH2F, {{500, 0, 10}, {300, -1.5, 1.5}});
    mHistogramRegistry->add((folderName + folderSuffix + "/hMomentumVsPhi").c_str(), "; #it{p} (GeV/#it{c}); #phi", kTH2F, {{500, 0, 10}, {360, 0., TMath::TwoPi()}});
    mHistogramRegistry->add((folderName + folderSuffix + "/hEtaVsPhi").c_str(), "; #eta; #phi", kTH2F, {{300, -1.5, 1.5}, {360, 0., TMath::TwoPi()}});

    if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kTrack || mParticleType == o2::aod::femtodreamparticle::ParticleType::kV0Child || mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascadeV0Child || mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascadeBachelor) {
      mHistogramRegistry->add((folderName + folderSuffix + "/hCharge").c_str(), "; Charge; Entries", kTH1F, {{5, -2.5, 2.5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hTPCfindable").c_str(), "; TPC findable clusters; Entries", kTH1F, {{163, -0.5, 162.5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hTPCfound").c_str(), "; TPC found clusters; Entries", kTH1F, {{163, -0.5, 162.5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hTPCcrossedOverFindable").c_str(), "; TPC ratio findable over crossed; Entries", kTH1F, {{100, 0.5, 1.5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hTPCcrossedRows").c_str(), "; TPC crossed rows; Entries", kTH1F, {{163, -0.5, 162.5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hTPCshared").c_str(), "; TPC shared clusters; Entries", kTH1F, {{163, -0.5, 162.5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hTPCsharedOverFound").c_str(), "; TPC ratio shared over found; Entries", kTH1F, {{1000, 0, 1}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hTPCfindableVsCrossed").c_str(), ";TPC findable clusters ; TPC crossed rows;", kTH2F, {{163, -0.5, 162.5}, {163, -0.5, 162.5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hTPCfoundVsShared").c_str(), ";TPC found clusters ; TPC shared clusters;", kTH2F, {{163, -0.5, 162.5}, {163, -0.5, 162.5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hITSclusters").c_str(), "; ITS clusters; Entries", kTH1F, {{10, -0.5, 9.5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hITSclustersIB").c_str(), "; ITS clusters in IB; Entries", kTH1F, {{10, -0.5, 9.5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hDCAz").c_str(), "; #it{p} (GeV/#it{c}); DCA_{z} (cm)", kTH2F, {pTAxis, dcazAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/hDCA").c_str(), "; #it{p} (GeV/#it{c}); DCA (cm)", kTH2F, {pTAxis, {300, 0., 1.5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hTPCdEdX").c_str(), "; #it{p} (GeV/#it{c}); TPC Signal", kTH2F, {{100, 0, 10}, {1000, 0, 1000}});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaTPC_el").c_str(), "n#sigma_{TPC}^{e}", kTH2F, {pTAxis, NsigmaTPCAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaTPC_pi").c_str(), "n#sigma_{TPC}^{#pi}", kTH2F, {pTAxis, NsigmaTPCAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaTPC_K").c_str(), "n#sigma_{TPC}^{K}", kTH2F, {pTAxis, NsigmaTPCAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaTPC_p").c_str(), "n#sigma_{TPC}^{p}", kTH2F, {pTAxis, NsigmaTPCAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaTPC_d").c_str(), "n#sigma_{TPC}^{d}", kTH2F, {pTAxis, NsigmaTPCAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaTPC_tr").c_str(), "n#sigma_{TPC}^{tr}", kTH2F, {pTAxis, NsigmaTPCAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaTPC_he3").c_str(), "n#sigma_{TPC}^{he3}", kTH2F, {pTAxis, NsigmaTPCAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaTOF_el").c_str(), "n#sigma_{TOF}^{e}", kTH2F, {pTAxis, NsigmaTOFAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaTOF_pi").c_str(), "n#sigma_{TOF}^{#pi}", kTH2F, {pTAxis, NsigmaTOFAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaTOF_K").c_str(), "n#sigma_{TOF}^{K}", kTH2F, {pTAxis, NsigmaTOFAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaTOF_p").c_str(), "n#sigma_{TOF}^{p}", kTH2F, {pTAxis, NsigmaTOFAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaTOF_d").c_str(), "n#sigma_{TOF}^{d}", kTH2F, {pTAxis, NsigmaTOFAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaTOF_tr").c_str(), "n#sigma_{TOF}^{tr}", kTH2F, {pTAxis, NsigmaTOFAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaTOF_he3").c_str(), "n#sigma_{TOF}^{he3}", kTH2F, {pTAxis, NsigmaTOFAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaComb_el").c_str(), "n#sigma_{comb}^{e}", kTH2F, {pTAxis, NsigmaTPCTOFAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaComb_pi").c_str(), "n#sigma_{comb}^{#pi}", kTH2F, {pTAxis, NsigmaTPCTOFAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaComb_K").c_str(), "n#sigma_{comb}^{K}", kTH2F, {pTAxis, NsigmaTPCTOFAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaComb_p").c_str(), "n#sigma_{comb}^{p}", kTH2F, {pTAxis, NsigmaTPCTOFAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaComb_d").c_str(), "n#sigma_{comb}^{d}", kTH2F, {pTAxis, NsigmaTPCTOFAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaComb_tr").c_str(), "n#sigma_{comb}^{tr}", kTH2F, {pTAxis, NsigmaTPCTOFAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaComb_he3").c_str(), "n#sigma_{comb}^{he3}", kTH2F, {pTAxis, NsigmaTPCTOFAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/ITSSignal").c_str(), "<cluster size>x<cos#lambda>", kTH2F, {pTAxis, NsigmaITSAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaITS_el").c_str(), "n#sigma_{ITS}^{e}", kTH2F, {pTAxis, NsigmaITSAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaITS_pi").c_str(), "n#sigma_{ITS}^{#pi}", kTH2F, {pTAxis, NsigmaITSAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaITS_K").c_str(), "n#sigma_{ITS}^{K}", kTH2F, {pTAxis, NsigmaITSAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaITS_p").c_str(), "n#sigma_{ITS}^{p}", kTH2F, {pTAxis, NsigmaITSAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaITS_d").c_str(), "n#sigma_{ITS}^{d}", kTH2F, {pTAxis, NsigmaITSAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaITS_tr").c_str(), "n#sigma_{ITS}^{tr}", kTH2F, {pTAxis, NsigmaITSAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaITS_he3").c_str(), "n#sigma_{ITS}^{he3}", kTH2F, {pTAxis, NsigmaITSAxis});
      if (correlatedPlots) {
        mHistogramRegistry->add((folderName + folderSuffix + "/HighDcorrelator").c_str(), "", kTHnSparseF, {multAxis, multPercentileAxis, pTAxis, etaAxis, phiAxis, tempFitVarAxis, dcazAxis, NsigmaTPCAxis, NsigmaTOFAxis});
      }
    } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kV0 || mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascadeV0) {
      mHistogramRegistry->add((folderName + folderSuffix + "/hDaughDCA").c_str(), "; DCA^{daugh} (cm); Entries", kTH1F, {{1000, 0, 10}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hTransRadius").c_str(), "; #it{r}_{xy} (cm); Entries", kTH1F, {{1500, 0, 150}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hDecayVtxX").c_str(), "; #it{Vtx}_{x} (cm); Entries", kTH1F, {{2000, 0, 200}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hDecayVtxY").c_str(), "; #it{Vtx}_{y} (cm)); Entries", kTH1F, {{2000, 0, 200}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hDecayVtxZ").c_str(), "; #it{Vtx}_{z} (cm); Entries", kTH1F, {{2000, 0, 200}});
    } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascade) {
      mHistogramRegistry->add((folderName + folderSuffix + "/hDaughDCA").c_str(), "; DCA^{daugh} (cm); Entries", kTH1F, {{1000, 0, 10}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hTransRadius").c_str(), "; #it{r}_{xy} (cm); Entries", kTH1F, {{1500, 0, 150}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hDecayVtxX").c_str(), "; #it{Vtx}_{x} (cm); Entries", kTH1F, {{2000, 0, 200}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hDecayVtxY").c_str(), "; #it{Vtx}_{y} (cm)); Entries", kTH1F, {{2000, 0, 200}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hDecayVtxZ").c_str(), "; #it{Vtx}_{z} (cm); Entries", kTH1F, {{2000, 0, 200}});
    }
  }

  /// Initializes specialized Monte Carlo particle histograms
  /// internal function called by init only in case of Monte Carlo truth
  /// \tparam T type of the axis Object
  /// \param folderName base path of the directory in the output file, in which to store the histograms
  /// \param tempFitVarAxisTitle  Title of the axis of the tempFitVar (DCA_xy in case of tracks, CPA in case of V0s, etc.)
  /// \param tempFitVarpTAxis axis object for the pT axis in the pT vs. tempFitVar plots
  /// \param tempFitVarAxis axis object for the tempFitVar axis
  template <typename T>
  void init_MC(std::string folderName, std::string /*tempFitVarAxisTitle*/, T& tempFitVarpTAxis, T& tempFitVarAxis, T& dcazAxis, T& multAxis, bool isDebug)
  {
    /// Particle-type specific histograms
    std::string folderSuffix = static_cast<std::string>(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kTruth]).c_str();

    mHistogramRegistry->add((folderName + folderSuffix + "/hPt_ReconNoFake").c_str(), "; #it{p}_{T} (GeV/#it{c}); Entries", kTH1F, {{240, 0, 6}});
    mHistogramRegistry->add((folderName + folderSuffix + "/hPDG").c_str(), "; PDG; Entries", kTH1I, {{6001, -3000.5, 3000.5}});
    mHistogramRegistry->add((folderName + folderSuffix + "/hOrigin_MC").c_str(), "; Origin; Entries", kTH1I, {{7, -0.5, 6.5}});
    mHistogramRegistry->add((folderName + folderSuffix + "/hNoMCtruthCounter").c_str(), "; Counter; Entries", kTH1I, {{1, -0.5, 0.5}});
    mHistogramRegistry->add((folderName + folderSuffix + "/hPt_DiffTruthReco").c_str(), "; p^{truth}_{T}; (p^{reco}_{T} - p^{truth}_{T}) / p^{truth}_{T}", kTH2F, {tempFitVarpTAxis, {200, -1, 1}});
    mHistogramRegistry->add((folderName + folderSuffix + "/hEta_DiffTruthReco").c_str(), "; #eta^{truth}; #eta^{reco} - #eta^{truth}", kTH2F, {{200, -1, 1}, {200, -1, 1}});
    mHistogramRegistry->add((folderName + folderSuffix + "/hPhi_DiffTruthReco").c_str(), "; #varphi^{truth}; #varphi^{reco} - #varphi^{truth}", kTH2F, {{720, 0, TMath::TwoPi()}, {200, -1, 1}});

    if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kTrack || mParticleType == o2::aod::femtodreamparticle::ParticleType::kV0Child || mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascadeV0Child || mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascadeBachelor) {
      /// Track histograms
      if (isDebug) {
        mHistogramRegistry->add((folderName + folderSuffix + "/Debug/hPDGmother_Primary").c_str(), "; PDG mother; Entries", kTH1I, {{6001, -3000.5, 3000.5}});
        mHistogramRegistry->add((folderName + folderSuffix + "/Debug/hPDGmother_Secondary").c_str(), "; PDG mother; Entries", kTH1I, {{6001, -3000.5, 3000.5}});
        mHistogramRegistry->add((folderName + folderSuffix + "/Debug/hPDGmother_Material").c_str(), "; PDG mother; Entries", kTH1I, {{6001, -3000.5, 3000.5}});
        mHistogramRegistry->add((folderName + folderSuffix + "/Debug/hPDGmother_WrongCollision").c_str(), "; PDG mother; Entries", kTH1I, {{6001, -3000.5, 3000.5}});
        mHistogramRegistry->add((folderName + folderSuffix + "/Debug/hPDGmother_Fake").c_str(), "; PDG mother; Entries", kTH1I, {{6001, -3000.5, 3000.5}});
        mHistogramRegistry->add((folderName + folderSuffix + "/Debug/hPDGmother_Else").c_str(), "; PDG mother; Entries", kTH1I, {{6001, -3000.5, 3000.5}});
        mHistogramRegistry->add((folderName + folderSuffix + "/Debug/hPDGmother_SecondaryDaughterLambda").c_str(), "; PDG mother; Entries", kTH1I, {{12001, -6000.5, 6000.5}});
        mHistogramRegistry->add((folderName + folderSuffix + "/Debug/hPDGmother_SecondaryDaughterSigmaplus").c_str(), "; PDG mother; Entries", kTH1I, {{12001, -6000.5, 6000.5}});

        mHistogramRegistry->add((folderName + folderSuffix + "/hDCAxy_Primary").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTHnSparseF, {tempFitVarpTAxis, tempFitVarAxis, dcazAxis, multAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hDCAxy_Secondary").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTHnSparseF, {tempFitVarpTAxis, tempFitVarAxis, dcazAxis, multAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hDCAxy_Material").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTHnSparseF, {tempFitVarpTAxis, tempFitVarAxis, dcazAxis, multAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hDCAxy_WrongCollision").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTHnSparseF, {tempFitVarpTAxis, tempFitVarAxis, dcazAxis, multAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hDCAxy_Fake").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTHnSparseF, {tempFitVarpTAxis, tempFitVarAxis, dcazAxis, multAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hDCAxy_SecondaryDaughterLambda").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTHnSparseF, {tempFitVarpTAxis, tempFitVarAxis, dcazAxis, multAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hDCAxy_SecondaryDaughterSigmaplus").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTHnSparseF, {tempFitVarpTAxis, tempFitVarAxis, dcazAxis, multAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hDCAxy_Else").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTHnSparseF, {tempFitVarpTAxis, tempFitVarAxis, dcazAxis, multAxis});
      } else {
        mHistogramRegistry->add((folderName + folderSuffix + "/hDCAxy_Primary").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hDCAxy_Secondary").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hDCAxy_Material").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hDCAxy_WrongCollision").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hDCAxy_Fake").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hDCAxy_SecondaryDaughterLambda").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hDCAxy_SecondaryDaughterSigmaplus").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hDCAxy_Else").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
      }

      // DCA plots
    } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kV0 || mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascadeV0) {
      if (isDebug) {
        mHistogramRegistry->add((folderName + folderSuffix + "/hCPA_Primary").c_str(), "; #it{p}_{T} (GeV/#it{c}); CPA", kTHnSparseF, {tempFitVarpTAxis, tempFitVarAxis, multAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hCPA_Secondary").c_str(), "; #it{p}_{T} (GeV/#it{c}); CPA", kTHnSparseF, {tempFitVarpTAxis, tempFitVarAxis, multAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hCPA_Material").c_str(), "; #it{p}_{T} (GeV/#it{c}); CPA", kTHnSparseF, {tempFitVarpTAxis, tempFitVarAxis, multAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hCPA_WrongCollision").c_str(), "; #it{p}_{T} (GeV/#it{c}); CPA", kTHnSparseF, {tempFitVarpTAxis, tempFitVarAxis, multAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hCPA_Fake").c_str(), "; #it{p}_{T} (GeV/#it{c}); CPA", kTHnSparseF, {tempFitVarpTAxis, tempFitVarAxis, multAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hCPA_Else").c_str(), "; #it{p}_{T} (GeV/#it{c}); CPA", kTHnSparseF, {tempFitVarpTAxis, tempFitVarAxis, multAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hCPA_Secondary_SIGMA0").c_str(), "; #it{p}_{T} (GeV/#it{c}); CPA", kTHnSparseF, {tempFitVarpTAxis, tempFitVarAxis, multAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hCPA_Secondary_XIMinus").c_str(), "; #it{p}_{T} (GeV/#it{c}); CPA", kTHnSparseF, {tempFitVarpTAxis, tempFitVarAxis, multAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hCPA_Secondary_XI0").c_str(), "; #it{p}_{T} (GeV/#it{c}); CPA", kTHnSparseF, {tempFitVarpTAxis, tempFitVarAxis, multAxis});
      } else {
        mHistogramRegistry->add((folderName + folderSuffix + "/hCPA_Primary").c_str(), "; #it{p}_{T} (GeV/#it{c}); CPA", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hCPA_Secondary").c_str(), "; #it{p}_{T} (GeV/#it{c}); CPA", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hCPA_Material").c_str(), "; #it{p}_{T} (GeV/#it{c}); CPA", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hCPA_WrongCollision").c_str(), "; #it{p}_{T} (GeV/#it{c}); CPA", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hCPA_Fake").c_str(), "; #it{p}_{T} (GeV/#it{c}); CPA", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hCPA_Else").c_str(), "; #it{p}_{T} (GeV/#it{c}); CPA", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hCPA_Secondary_SIGMA0").c_str(), "; #it{p}_{T} (GeV/#it{c}); CPA", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hCPA_Secondary_XIMinus").c_str(), "; #it{p}_{T} (GeV/#it{c}); CPA", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
        mHistogramRegistry->add((folderName + folderSuffix + "/hCPA_Secondary_XI0").c_str(), "; #it{p}_{T} (GeV/#it{c}); CPA", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
      }
      /// V0 histograms
      ///  to be implemented
    } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascade) {
      /// Cascade histograms
      /// to be implemented
    } else {
      LOG(fatal) << "FemtoDreamParticleHisto: Histogramming for requested object not defined - quitting!";
    }
  }

  /// Templated function for the initialization of the QA histograms
  /// Always calls init_base to initialize the histograms with data/ Monte Carlo reconstructed
  /// In case of Monte Carlo, calls init_base again for Monte Carlo truth and the specialized function init_MC for additional Monte Carlo histogramms
  /// \tparam T type of the axis binning
  /// \param registry Histogram registry to be passed
  /// \param tempFitVarpTBins binning of the pT axis in the pT vs. tempFitVar
  /// \param tempFitVarBins binning of the tempFitVar (DCA_xy in case of tracks, CPA in case of V0s, etc.)
  /// \param isMC add Monte Carlo truth histograms to the output file
  template <typename T>
  void init(HistogramRegistry* registry, T& MultBins, T& PercentileBins, T& pTBins, T& etaBins, T& phiBins, T& tempFitVarBins, T& NsigmaTPCBins, T& NsigmaTOFBins, T& NsigmaTPCTOFBins, T& NsigmaITSBins, T& InvMassBins, bool isMC, int pdgCode, bool isDebug = false, bool correlatedPlots = false)
  {
    mPDG = pdgCode;
    if (registry) {
      mHistogramRegistry = registry;
      /// The folder names are defined by the type of the object and the suffix (if applicable)
      std::string tempFitVarAxisTitle;
      if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kTrack || mParticleType == o2::aod::femtodreamparticle::ParticleType::kV0Child || mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascadeV0Child || mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascadeBachelor) {
        /// Track histograms
        tempFitVarAxisTitle = "DCA_{xy} (cm)";
      } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kV0 || mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascadeV0) {
        /// V0 histograms
        tempFitVarAxisTitle = "cos#alpha";
      } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascade) {
        /// Cascade histograms
        tempFitVarAxisTitle = "cos#alpha";
      } else {
        LOG(fatal) << "FemtoDreamParticleHisto: Histogramming for requested object not defined - quitting!";
      }

      framework::AxisSpec multAxis = {MultBins, " Multiplicity"};
      framework::AxisSpec multPercentileAxis = {PercentileBins, "Multiplicity Percentile (%)"};
      framework::AxisSpec pTAxis = {pTBins, "#it{p}_{T} (GeV/#it{c})"};
      framework::AxisSpec etaAxis = {etaBins, "#eta"};
      framework::AxisSpec phiAxis = {phiBins, "#phi"};
      framework::AxisSpec tempFitVarAxis = {tempFitVarBins, tempFitVarAxisTitle};
      framework::AxisSpec dcazAxis = {tempFitVarBins, "DCA_{z} (cm)"};
      framework::AxisSpec NsigmaTPCAxis = {NsigmaTPCBins, "n#sigma_{TPC}"};
      framework::AxisSpec NsigmaTOFAxis = {NsigmaTOFBins, "n#sigma_{TOF}"};
      framework::AxisSpec NsigmaTPCTOFAxis = {NsigmaTPCTOFBins, "n#sigma_{TPC+TOF}"};
      framework::AxisSpec NsigmaITSAxis = {NsigmaITSBins, "n#sigma_{ITS}"};
      framework::AxisSpec InvMassAxis = {InvMassBins, "M_{inv} (GeV/#it{c}^{2})"};

      std::string folderName = (static_cast<std::string>(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]).c_str() + static_cast<std::string>(mFolderSuffix[mFolderSuffixType])).c_str();

      // Fill here the actual histogramms by calling init_base and init_MC
      init_base<o2::aod::femtodreamMCparticle::MCType::kRecon>(folderName, tempFitVarAxisTitle, pTAxis, tempFitVarAxis, InvMassAxis, multAxis);
      if (isDebug) {
        init_debug<o2::aod::femtodreamMCparticle::MCType::kRecon>(folderName, multAxis, multPercentileAxis, pTAxis, etaAxis, phiAxis, tempFitVarAxis, dcazAxis, NsigmaTPCAxis, NsigmaTOFAxis, NsigmaTPCTOFAxis, NsigmaITSAxis, correlatedPlots);
      }
      if (isMC) {
        init_base<o2::aod::femtodreamMCparticle::MCType::kTruth>(folderName, tempFitVarAxisTitle, pTAxis, tempFitVarAxis, InvMassAxis, multAxis);
        init_MC(folderName, tempFitVarAxisTitle, pTAxis, tempFitVarAxis, dcazAxis, multAxis, isDebug);
      }
    }
  }

  /// Filling of the histograms
  /// Called by init both in case of reconstructed data/ Monte Carlo, and for Monte Carlo Truth
  /// \tparam T Data type of the particle
  /// \param part Particle
  template <o2::aod::femtodreamMCparticle::MCType mc, typename T>
  void fillQA_base(T const& part, const int /*mult*/)
  {
    /// Histograms of the kinematic properties
    mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hPt"), part.pt());
    mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hEta"), part.eta());
    mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hPhi"), part.phi());

    /// particle specific histogramms for the TempFitVar column in FemtoDreamParticles
    if constexpr (mc == o2::aod::femtodreamMCparticle::MCType::kRecon) {
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST(o2::aod::femtodreamparticle::TempFitVarName[mParticleType]), part.pt(), part.tempFitVar());
    }
    if constexpr ((mParticleType == o2::aod::femtodreamparticle::ParticleType::kV0 || mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascadeV0) && mc == o2::aod::femtodreamMCparticle::MCType::kRecon) {
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hInvMassLambda"), part.mLambda());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hpTInvMassLambda"), part.pt(), part.mLambda());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hInvMassAntiLambda"), part.mAntiLambda());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hpTInvMassAntiLambda"), part.pt(), part.mAntiLambda());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hInvMassLambdaAntiLambda"), part.mLambda(), part.mAntiLambda());
    }
    if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascade) {
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hInvMassCascade"), part.mLambda());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hpTInvMassCascade"), part.pt(), part.mLambda());
      // mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hInvMassCascade"), part.mLambda());
      // mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hpTInvMassCascade"), part.pt(), part.mLambda());
    }
  }

  template <o2::aod::femtodreamMCparticle::MCType mc, typename Tpart>
  void fillQA_debug(Tpart const& part, aod::femtodreamparticle::MomentumType MomentumType, const int mult, const float cent, bool correlatedPlots)
  {
    float momentum;
    switch (MomentumType) {
      case aod::femtodreamparticle::kPt:
        momentum = part.pt();
        break;
      case aod::femtodreamparticle::kPreco:
        momentum = part.p();
        break;
      case aod::femtodreamparticle::kPtpc:
        LOG(warn) << "Momentum at TPC inner wall is only available if the Producer Task generates debug information.";
        momentum = part.tpcInnerParam();
        break;
      default:
        LOG(warn) << "MomentumType " << MomentumType << " not implemented. Use pT of the Track for debug Histograms.";
        momentum = part.pt();
    }

    mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hMomentumVsEta"), momentum, part.eta());
    mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hMomentumVsPhi"), momentum, part.phi());
    mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hEtaVsPhi"), part.eta(), part.phi());

    // Histograms holding further debug information
    if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kTrack || mParticleType == o2::aod::femtodreamparticle::ParticleType::kV0Child || mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascadeV0Child || mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascadeBachelor) {
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hCharge"), part.sign());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hTPCfindable"), part.tpcNClsFindable());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hTPCfound"), part.tpcNClsFound());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hTPCcrossedOverFindable"), part.tpcCrossedRowsOverFindableCls());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hTPCcrossedRows"), part.tpcNClsCrossedRows());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hTPCfindableVsCrossed"), part.tpcNClsFindable(), part.tpcNClsCrossedRows());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hTPCshared"), part.tpcNClsShared());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hTPCsharedOverFound"), static_cast<float>(part.tpcNClsShared()) / static_cast<float>(part.tpcNClsFound()));
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hTPCfoundVsShared"), part.tpcNClsFound(), part.tpcNClsShared());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hITSclusters"), part.itsNCls());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hITSclustersIB"), part.itsNClsInnerBarrel());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hDCAz"), momentum, part.dcaZ());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hDCA"), momentum, std::sqrt(std::pow(part.dcaXY(), 2.) + std::pow(part.dcaZ(), 2.)));
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hTPCdEdX"), momentum, part.tpcSignal());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaTPC_el"), momentum, part.tpcNSigmaEl());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaTPC_pi"), momentum, part.tpcNSigmaPi());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaTPC_K"), momentum, part.tpcNSigmaKa());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaTPC_p"), momentum, part.tpcNSigmaPr());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaTPC_d"), momentum, part.tpcNSigmaDe());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaTPC_tr"), momentum, part.tpcNSigmaTr());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaTPC_he3"), momentum, part.tpcNSigmaHe());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaTOF_el"), momentum, part.tofNSigmaEl());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaTOF_pi"), momentum, part.tofNSigmaPi());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaTOF_K"), momentum, part.tofNSigmaKa());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaTOF_p"), momentum, part.tofNSigmaPr());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaTOF_d"), momentum, part.tofNSigmaDe());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaTOF_tr"), momentum, part.tofNSigmaTr());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaTOF_he3"), momentum, part.tofNSigmaHe());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaComb_el"), momentum, std::sqrt(part.tpcNSigmaEl() * part.tpcNSigmaEl() + part.tofNSigmaEl() * part.tofNSigmaEl()));
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaComb_pi"), momentum, std::sqrt(part.tpcNSigmaPi() * part.tpcNSigmaPi() + part.tofNSigmaPi() * part.tofNSigmaPi()));
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaComb_K"), momentum, std::sqrt(part.tpcNSigmaKa() * part.tpcNSigmaKa() + part.tofNSigmaKa() * part.tofNSigmaKa()));
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaComb_p"), momentum, std::sqrt(part.tpcNSigmaPr() * part.tpcNSigmaPr() + part.tofNSigmaPr() * part.tofNSigmaPr()));
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaComb_d"), momentum, std::sqrt(part.tpcNSigmaDe() * part.tpcNSigmaDe() + part.tofNSigmaDe() * part.tofNSigmaDe()));
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaComb_tr"), momentum, std::sqrt(part.tpcNSigmaTr() * part.tpcNSigmaTr() + part.tofNSigmaTr() * part.tofNSigmaTr()));
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaComb_he3"), momentum, std::sqrt(part.tpcNSigmaHe() * part.tpcNSigmaHe() + part.tofNSigmaHe() * part.tofNSigmaHe()));
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/ITSSignal"), momentum, part.itsSignal());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaITS_el"), momentum, part.itsNSigmaEl());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaITS_pi"), momentum, part.itsNSigmaPi());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaITS_K"), momentum, part.itsNSigmaKa());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaITS_p"), momentum, part.itsNSigmaPr());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaITS_d"), momentum, part.itsNSigmaDe());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaITS_tr"), momentum, part.itsNSigmaTr());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/nSigmaITS_he3"), momentum, part.itsNSigmaHe());

      if (correlatedPlots) {

        float pidTPC = 0.;
        float pidTOF = 0.;

        switch (abs(mPDG)) {
          case kElectron:
            pidTPC = part.tpcNSigmaEl();
            pidTOF = part.tofNSigmaEl();
            break;
          case kPiPlus:
            pidTPC = part.tpcNSigmaPi();
            pidTOF = part.tofNSigmaPi();
            break;
          case kKPlus:
            pidTPC = part.tpcNSigmaKa();
            pidTOF = part.tofNSigmaKa();
            break;
          case kProton:
            pidTPC = part.tpcNSigmaPr();
            pidTOF = part.tofNSigmaPr();
            break;
          case constants::physics::kDeuteron:
            pidTPC = part.tpcNSigmaDe();
            pidTOF = part.tofNSigmaDe();
            break;
          case constants::physics::kTriton:
            pidTPC = part.tpcNSigmaTr();
            pidTOF = part.tofNSigmaTr();
            break;
          case constants::physics::kHelium3:
            pidTPC = part.tpcNSigmaHe();
            pidTOF = part.tofNSigmaHe();
            break;
          default:
            LOG(warn) << "PDG code " << mPDG << " not supported. No PID information will be used.";
        }

        mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/HighDcorrelator"), mult,
                                 cent,
                                 momentum,
                                 part.eta(),
                                 part.phi(),
                                 part.dcaXY(),
                                 part.dcaZ(),
                                 pidTPC,
                                 pidTOF);
      }
    } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kV0 || mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascadeV0) {
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hDaughDCA"), part.daughDCA());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hTransRadius"), part.transRadius());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hDecayVtxX"), part.decayVtxX());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hDecayVtxY"), part.decayVtxY());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hDecayVtxZ"), part.decayVtxZ());
    } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascade) {
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hDaughDCA"), part.daughDCA());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hTransRadius"), part.transRadius());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hDecayVtxX"), part.decayVtxX());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hDecayVtxY"), part.decayVtxY());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hDecayVtxZ"), part.decayVtxZ());
    }
  }

  /// Filling specialized histograms for Monte Carlo truth
  /// internal function called by init only in case of Monte Carlo truth
  /// \tparam T Data type of the particle
  /// \tparam TMC Data typ of the Monte Carlo Particle
  /// \param part Particle
  /// \param mctruthorigin Origin of the associated mc Truth particle
  /// \param pdgcode PDG of the associated mc Truth particle associated to the reconstructed particle part
  template <bool isDebug, typename T>
  void fillQA_MC(T const& part, int mctruthorigin, int pdgcode, const int mult)
  {
    if (mHistogramRegistry) {
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hPDG"), pdgcode);
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hOrigin_MC"), mctruthorigin);

      auto MCpart = part.fdMCParticle();

      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hPt_DiffTruthReco"), MCpart.pt(), (part.pt() - MCpart.pt()) / MCpart.pt());
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hEta_DiffTruthReco"), MCpart.eta(), (part.eta() - MCpart.eta()));
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hPhi_DiffTruthReco"), MCpart.phi(), (part.phi() - MCpart.phi()));

      if (abs(pdgcode) == mPDG) { // fill this histogramm only for TRUE protons (independently of their origin) for the track purity estimation
        mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hPt_ReconNoFake"), part.pt());
      }
      if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kTrack || mParticleType == o2::aod::femtodreamparticle::ParticleType::kV0Child) {
        /// Track histograms
        if constexpr (isDebug) {
          switch (mctruthorigin) {
            case (o2::aod::femtodreamMCparticle::kPrimary):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/Debug/hPDGmother_Primary"), part.fdExtMCParticle().motherPDG());
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_Primary"),
                                       part.pt(), part.tempFitVar(), part.dcaZ(), mult);
              break;
            case (o2::aod::femtodreamMCparticle::kSecondary):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/Debug/hPDGmother_Secondary"), part.fdExtMCParticle().motherPDG());
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_Secondary"),
                                       part.pt(), part.tempFitVar(), part.dcaZ(), mult);
              break;
            case (o2::aod::femtodreamMCparticle::kMaterial):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/Debug/hPDGmother_Material"), part.fdExtMCParticle().motherPDG());
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_Material"),
                                       part.pt(), part.tempFitVar(), part.dcaZ(), mult);
              break;
            case (o2::aod::femtodreamMCparticle::kWrongCollision):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/Debug/hPDGmother_WrongCollision"), part.fdExtMCParticle().motherPDG());
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_WrongCollision"),
                                       part.pt(), part.tempFitVar(), part.dcaZ(), mult);
              break;
            case (o2::aod::femtodreamMCparticle::kFake):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/Debug/hPDGmother_Fake"), part.fdExtMCParticle().motherPDG());
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_Fake"),
                                       part.pt(), part.tempFitVar(), part.dcaZ(), mult);
              break;
            case (o2::aod::femtodreamMCparticle::kSecondaryDaughterLambda):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/Debug/hPDGmother_SecondaryDaughterLambda"), part.fdExtMCParticle().motherPDG());
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_SecondaryDaughterLambda"),
                                       part.pt(), part.tempFitVar(), part.dcaZ(), mult);
              break;
            case (o2::aod::femtodreamMCparticle::kSecondaryDaughterSigmaplus):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/Debug/hPDGmother_SecondaryDaughterSigmaplus"), part.fdExtMCParticle().motherPDG());
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_SecondaryDaughterSigmaplus"),
                                       part.pt(), part.tempFitVar(), part.dcaZ(), mult);
              break;
            case (o2::aod::femtodreamMCparticle::kElse):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/Debug/hPDGmother_Else"), part.fdExtMCParticle().motherPDG());
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_Else"),
                                       part.pt(), part.tempFitVar(), part.dcaZ(), mult);
              break;
            default:
              LOG(fatal) << "femtodreamparticleMC: not known value for ParticleOriginMCTruth - please check. Quitting!";
          }
        } else {
          switch (mctruthorigin) {
            case (o2::aod::femtodreamMCparticle::kPrimary):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_Primary"),
                                       part.pt(), part.tempFitVar());
              break;
            case (o2::aod::femtodreamMCparticle::kSecondary):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_Secondary"),
                                       part.pt(), part.tempFitVar());
              break;
            case (o2::aod::femtodreamMCparticle::kMaterial):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_Material"),
                                       part.pt(), part.tempFitVar());
              break;
            case (o2::aod::femtodreamMCparticle::kWrongCollision):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_WrongCollision"),
                                       part.pt(), part.tempFitVar());
              break;
            case (o2::aod::femtodreamMCparticle::kFake):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_Fake"),
                                       part.pt(), part.tempFitVar());
              break;
            case (o2::aod::femtodreamMCparticle::kSecondaryDaughterLambda):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_SecondaryDaughterLambda"),
                                       part.pt(), part.tempFitVar());
              break;
            case (o2::aod::femtodreamMCparticle::kSecondaryDaughterSigmaplus):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_SecondaryDaughterSigmaplus"),
                                       part.pt(), part.tempFitVar());
              break;
            case (o2::aod::femtodreamMCparticle::kElse):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_Else"),
                                       part.pt(), part.tempFitVar());
              break;
            default:
              LOG(fatal) << "femtodreamparticleMC: not known value for ParticleOriginMCTruth - please check. Quitting!";
          }
        }
      } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kV0) {
        /// V0 histograms
        if constexpr (isDebug) {
          switch (mctruthorigin) {
            case (o2::aod::femtodreamMCparticle::kPrimary):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hCPA_Primary"),
                                       part.pt(), part.tempFitVar(), mult);
              break;
            case (o2::aod::femtodreamMCparticle::kSecondary):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hCPA_Secondary"),
                                       part.pt(), part.tempFitVar(), mult);
              break;
            case (o2::aod::femtodreamMCparticle::kMaterial):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hCPA_Material"),
                                       part.pt(), part.tempFitVar(), mult);
              break;
            case (o2::aod::femtodreamMCparticle::kWrongCollision):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hCPA_WrongCollision"),
                                       part.pt(), part.tempFitVar(), mult);
              break;
            case (o2::aod::femtodreamMCparticle::kFake):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hCPA_Fake"),
                                       part.pt(), part.tempFitVar(), mult);
              break;
            case (o2::aod::femtodreamMCparticle::kElse):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hCPA_Else"),
                                       part.pt(), part.tempFitVar(), mult);
              break;
            case (o2::aod::femtodreamMCparticle::kSecondaryDaughterSigma0):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hCPA_Secondary_SIGMA0"),
                                       part.pt(), part.tempFitVar(), mult);
              break;
            case (o2::aod::femtodreamMCparticle::kSecondaryDaughterXiMinus):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hCPA_Secondary_XIMinus"),
                                       part.pt(), part.tempFitVar(), mult);
              break;
            case (o2::aod::femtodreamMCparticle::kSecondaryDaughterXi0):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hCPA_Secondary_XI0"),
                                       part.pt(), part.tempFitVar(), mult);
              break;
            default:
              LOG(fatal) << "femtodreamparticleMC: not known value for ParticleOriginMCTruth - please check. Quitting!";
          }
        } else {
          switch (mctruthorigin) {
            case (o2::aod::femtodreamMCparticle::kPrimary):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hCPA_Primary"),
                                       part.pt(), part.tempFitVar());
              break;
            case (o2::aod::femtodreamMCparticle::kSecondary):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hCPA_Secondary"),
                                       part.pt(), part.tempFitVar());
              break;
            case (o2::aod::femtodreamMCparticle::kMaterial):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hCPA_Material"),
                                       part.pt(), part.tempFitVar());
              break;
            case (o2::aod::femtodreamMCparticle::kWrongCollision):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hCPA_WrongCollision"),
                                       part.pt(), part.tempFitVar());
              break;
            case (o2::aod::femtodreamMCparticle::kFake):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hCPA_Fake"),
                                       part.pt(), part.tempFitVar());
              break;
            case (o2::aod::femtodreamMCparticle::kElse):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hCPA_Else"),
                                       part.pt(), part.tempFitVar());
              break;
            case (o2::aod::femtodreamMCparticle::kSecondaryDaughterSigma0):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hCPA_Secondary_SIGMA0"),
                                       part.pt(), part.tempFitVar());
              break;
            case (o2::aod::femtodreamMCparticle::kSecondaryDaughterXiMinus):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hCPA_Secondary_XIMinus"),
                                       part.pt(), part.tempFitVar());
              break;
            case (o2::aod::femtodreamMCparticle::kSecondaryDaughterXi0):
              mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hCPA_Secondary_XI0"),
                                       part.pt(), part.tempFitVar());
              break;
            default:
              LOG(fatal) << "femtodreamparticleMC: not known value for ParticleOriginMCTruth - please check. Quitting!";
          }
        }
      } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascade) {
        /// Cascade histograms
      } else {
        LOG(fatal) << "FemtoDreamParticleHisto: Histogramming for requested object not defined - quitting!";
      }
    }
  }

  /// Templated function to fill particle histograms for data/ Monte Carlo reconstructed and Monte Carlo truth
  /// Always calls fillQA_base fill histogramms with data/ Monte Carlo reconstructed
  /// In case of Monte Carlo, calls fillQA_base with Monte Carlo truth info and specialized function fillQA_MC for additional histogramms
  /// \tparam T particle type
  /// \tparam isMC fills the additional histograms for Monte Carlo truth
  /// \param part particle for which the histograms should be filled
  template <bool isMC, bool isDebug, typename Tpart>
  void fillQA(Tpart const& part, aod::femtodreamparticle::MomentumType MomemtumType, const int mult, const float cent, bool correlatedPlots = false)
  {
    if (mHistogramRegistry) {
      fillQA_base<o2::aod::femtodreamMCparticle::MCType::kRecon>(part, mult);
      if constexpr (isDebug) {
        fillQA_debug<o2::aod::femtodreamMCparticle::MCType::kRecon>(part, MomemtumType, mult, cent, correlatedPlots);
      }
      if constexpr (isMC) {
        if (part.has_fdMCParticle()) {
          fillQA_base<o2::aod::femtodreamMCparticle::MCType::kTruth>(part.fdMCParticle(), mult);
          fillQA_MC<isDebug>(part, (part.fdMCParticle()).partOriginMCTruth(), (part.fdMCParticle()).pdgMCTruth(), mult);
        } else {
          mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hNoMCtruthCounter"), 0);
        }
      }
    }
  }

 private:
  HistogramRegistry* mHistogramRegistry;                                                                                                                      ///< For QA output
  static constexpr o2::aod::femtodreamparticle::ParticleType mParticleType = particleType;                                                                    ///< Type of the particle under analysis
  static constexpr int mFolderSuffixType = suffixType;                                                                                                        ///< Counter for the folder suffix specified below
  static constexpr std::string_view mFolderSuffix[9] = {"", "_one", "_two", "_pos", "_neg", "_allSelected", "_allSelected_pos", "_allSelected_neg", "_bach"}; ///< Suffix for the folder name in case of analyses of pairs of the same kind (T-T, V-V, C-C)
  int mPDG = 0;                                                                                                                                               ///< PDG code of the selected particle
};
} // namespace o2::analysis::femtoDream

#endif // PWGCF_FEMTODREAM_CORE_FEMTODREAMPARTICLEHISTO_H_
