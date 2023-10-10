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

/// \file FemtoUniverseParticleHisto.h
/// \brief FemtoUniverseParticleHisto - Histogram class for tracks, V0s and cascades
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Georgios Mantzaridis, TU München, georgios.mantzaridis@tum.de
/// \author Anton Riedel, TU München, anton.riedel@tum.de
/// \author Zuzanna Chochulska, WUT Warsaw, zuzanna.chochulska.stud@pw.edu.pl

#ifndef PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEPARTICLEHISTO_H_
#define PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEPARTICLEHISTO_H_

#include <string>
#include "PWGCF/FemtoUniverse/DataModel/FemtoDerived.h"
#include "Framework/HistogramRegistry.h"

using namespace o2::framework;

namespace o2::analysis::femtoUniverse
{

/// \class FemtoUniverseParticleHisto
/// \brief Class for histogramming particle properties
/// \tparam particleType Type of the particle (Track/V0/Cascade/Phi/...)
/// \tparam suffixType (optional) Takes care of the suffix for the folder name in case of analyses of pairs of the same kind (T-T, V-V, C-C)
template <o2::aod::femtouniverseparticle::ParticleType particleType, int suffixType = 0>
class FemtoUniverseParticleHisto
{
 public:
  /// Destructor
  virtual ~FemtoUniverseParticleHisto() = default;

  /// Initializes particle histograms
  /// Called by init both in case of reconstructed data/ Monte Carlo, and for Monte Carlo Truth
  /// \tparam T type of the axis Object
  /// \tparam mc enum object to get the suffix ("" for data/ Monte Cartlo reconstructed, "_MC" for Monte Carlo truth) for the folder in the output file
  /// \param folderName base path of the directory in the output file, in which to store the histograms
  /// \param tempFitVarAxisTitle  Title of the axis of the tempFitVar (DCA_xy in case of tracks, CPA in case of V0s, etc.)
  /// \param tempFitVarpTAxis axis object for the pT axis in the pT vs. tempFitVar plots
  /// \param tempFitVarAxis axis object for the tempFitVar axis
  template <o2::aod::femtouniverseMCparticle::MCType mc, typename T>
  void init_base(std::string folderName, std::string tempFitVarAxisTitle, T& tempFitVarpTAxis, T& tempFitVarAxis)
  {
    std::string folderSuffix = static_cast<std::string>(o2::aod::femtouniverseMCparticle::MCTypeName[mc]).c_str();
    /// Histograms of the kinematic properties
    mHistogramRegistry->add((folderName + folderSuffix + "/hPt").c_str(), "; #it{p}_{T} (GeV/#it{c}); Entries", kTH1F, {{240, 0, 6}});
    mHistogramRegistry->add((folderName + folderSuffix + "/hEta").c_str(), "; #eta; Entries", kTH1F, {{200, -1.5, 1.5}});
    mHistogramRegistry->add((folderName + folderSuffix + "/hPhi").c_str(), "; #phi; Entries", kTH1F, {{200, 0, 2. * M_PI}});
    mHistogramRegistry->add((folderName + folderSuffix + "/hPhiEta").c_str(), "; #phi; #eta", kTH2F, {{200, 0, 2. * M_PI}, {200, -1.5, 1.5}});

    /// particle specific histogramms for the TempFitVar column in FemtoUniverseParticles
    if constexpr (o2::aod::femtouniverseMCparticle::MCType::kRecon == mc) {
      mHistogramRegistry->add((folderName + folderSuffix + static_cast<std::string>(o2::aod::femtouniverseparticle::TempFitVarName[mParticleType])).c_str(), ("; #it{p}_{T} (GeV/#it{c}); " + tempFitVarAxisTitle).c_str(), kTH2F, {{tempFitVarpTAxis}, {tempFitVarAxis}});
    }
  }

  // comment
  template <o2::aod::femtouniverseMCparticle::MCType mc>
  void init_debug(std::string folderName)
  {
    std::string folderSuffix = static_cast<std::string>(o2::aod::femtouniverseMCparticle::MCTypeName[mc]).c_str();
    if constexpr (mParticleType == o2::aod::femtouniverseparticle::ParticleType::kTrack || mParticleType == o2::aod::femtouniverseparticle::ParticleType::kV0Child || mParticleType == o2::aod::femtouniverseparticle::ParticleType::kMCTruthTrack) {
      mHistogramRegistry->add((folderName + folderSuffix + "/hCharge").c_str(), "; Charge; Entries", kTH1F, {{5, -2.5, 2.5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hTPCfindable").c_str(), "; TPC findable clusters; Entries", kTH1F, {{163, -0.5, 162.5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hTPCfound").c_str(), "; TPC found clusters; Entries", kTH1F, {{163, -0.5, 162.5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hTPCcrossedOverFindable").c_str(), "; TPC ratio findable; Entries", kTH1F, {{100, 0.5, 1.5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hTPCcrossedRows").c_str(), "; TPC crossed rows; Entries", kTH1F, {{163, -0.5, 162.5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hTPCfindableVsCrossed").c_str(), ";TPC findable clusters ; TPC crossed rows;", kTH2F, {{163, -0.5, 162.5}, {163, -0.5, 162.5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hTPCshared").c_str(), "; TPC shared clusters; Entries", kTH1F, {{163, -0.5, 162.5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hITSclusters").c_str(), "; ITS clusters; Entries", kTH1F, {{10, -0.5, 9.5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hITSclustersIB").c_str(), "; ITS clusters in IB; Entries", kTH1F, {{10, -0.5, 9.5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hDCAz").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{z} (cm)", kTH2F, {{100, 0, 10}, {500, -5, 5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hDCA").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA (cm)", kTH2F, {{100, 0, 10}, {301, 0., 1.5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hTPCdEdX").c_str(), "; #it{p} (GeV/#it{c}); TPC Signal", kTH2F, {{100, 0, 10}, {1000, 0, 1000}});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaTPC_el").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{e}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaTPC_pi").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{#pi}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaTPC_K").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{K}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaTPC_p").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{p}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaTPC_d").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{d}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaTOF_el").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{e}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaTOF_pi").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{#pi}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaTOF_K").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{K}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaTOF_p").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{p}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaTOF_d").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{d}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaComb_el").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{e}", kTH2F, {{100, 0, 10}, {100, 0, 5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaComb_pi").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{#pi}", kTH2F, {{100, 0, 10}, {100, 0, 5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaComb_K").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{K}", kTH2F, {{100, 0, 10}, {100, 0, 5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaComb_p").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{p}", kTH2F, {{100, 0, 10}, {100, 0, 5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/nSigmaComb_d").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{d}", kTH2F, {{100, 0, 10}, {100, 0, 5}});
    } else if constexpr (mParticleType == o2::aod::femtouniverseparticle::ParticleType::kV0) {
      mHistogramRegistry->add((folderName + folderSuffix + "/hDaughDCA").c_str(), "; DCA^{daugh} (cm); Entries", kTH1F, {{1000, 0, 10}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hTransRadius").c_str(), "; #it{r}_{xy} (cm); Entries", kTH1F, {{1500, 0, 150}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hDecayVtxX").c_str(), "; #it{Vtx}_{x} (cm); Entries", kTH1F, {{2000, 0, 200}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hDecayVtxY").c_str(), "; #it{Vtx}_{y} (cm)); Entries", kTH1F, {{2000, 0, 200}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hDecayVtxZ").c_str(), "; #it{Vtx}_{z} (cm); Entries", kTH1F, {{2000, 0, 200}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hInvMassLambda").c_str(), "; M_{#Lambda}; Entries", kTH1F, {{2000, 1.f, 3.f}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hInvMassAntiLambda").c_str(), "; M_{#bar{#Lambda}}; Entries", kTH1F, {{2000, 1.f, 3.f}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hInvMassLambdaAntiLambda").c_str(), "; M_{#Lambda}; M_{#bar{#Lambda}}", kTH2F, {{2000, 1.f, 3.f}, {2000, 1.f, 3.f}});
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
  void init_MC(std::string folderName, std::string tempFitVarAxisTitle, T& tempFitVarpTAxis, T& tempFitVarAxis)
  {
    /// Particle-type specific histograms
    std::string folderSuffix = static_cast<std::string>(o2::aod::femtouniverseMCparticle::MCTypeName[o2::aod::femtouniverseMCparticle::MCType::kTruth]).c_str();

    mHistogramRegistry->add((folderName + folderSuffix + "/hPt_ReconNoFake").c_str(), "; #it{p}_{T} (GeV/#it{c}); Entries", kTH1F, {{240, 0, 6}});

    if constexpr (mParticleType == o2::aod::femtouniverseparticle::ParticleType::kTrack || mParticleType == o2::aod::femtouniverseparticle::ParticleType::kV0Child || mParticleType == o2::aod::femtouniverseparticle::ParticleType::kMCTruthTrack) {
      /// Track histograms
      mHistogramRegistry->add((folderName + folderSuffix + "/hPDG").c_str(), "; PDG; Entries", kTH1I, {{6001, -3000, 3000}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hOrigin_MC").c_str(), "; Origin; Entries", kTH1I, {{7, 0, 7}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hNoMCtruthCounter").c_str(), "; Counter; Entries", kTH1I, {{1, 0, 1}});
      // DCA plots
      mHistogramRegistry->add((folderName + folderSuffix + "/hDCAxy_Material").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/hDCAxy_Fake").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/hDCAxy_DaughterLambda").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/hDCAxy_DaughterSigmaplus").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/hDCAxy_Primary").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
      mHistogramRegistry->add((folderName + folderSuffix + "/hDCAxy_Daughter").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
    } else if constexpr (mParticleType == o2::aod::femtouniverseparticle::ParticleType::kV0) {
      /// V0 histograms
      ///  to be implemented
    } else if constexpr (mParticleType == o2::aod::femtouniverseparticle::ParticleType::kCascade) {
      /// Cascade histograms
      /// to be implemented
    } else if constexpr (mParticleType == o2::aod::femtouniverseparticle::ParticleType::kPhi) {
      // Phi histograms
    } else {
      LOG(fatal) << "FemtoUniverseParticleHisto: Histogramming for requested object not defined - quitting!";
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
  void init(HistogramRegistry* registry, T& tempFitVarpTBins, T& tempFitVarBins, bool isMC, int pdgCode, bool isDebug = false)
  {
    mPDG = pdgCode;
    if (registry) {
      mHistogramRegistry = registry;
      /// The folder names are defined by the type of the object and the suffix (if applicable)
      std::string tempFitVarAxisTitle;
      if constexpr (mParticleType == o2::aod::femtouniverseparticle::ParticleType::kTrack || mParticleType == o2::aod::femtouniverseparticle::ParticleType::kV0Child) {
        /// Track histograms
        tempFitVarAxisTitle = "DCA_{xy} (cm)";
      } else if constexpr (mParticleType == o2::aod::femtouniverseparticle::ParticleType::kMCTruthTrack) {
        /// MC Truth Track histograms
        tempFitVarAxisTitle = "PDG code";
      } else if constexpr (mParticleType == o2::aod::femtouniverseparticle::ParticleType::kV0) {
        /// V0 histograms
        tempFitVarAxisTitle = "cos#alpha";
      } else if constexpr (mParticleType == o2::aod::femtouniverseparticle::ParticleType::kCascade) {
        /// Cascade histograms
        tempFitVarAxisTitle = "cos#alpha";
      } else if constexpr (mParticleType == o2::aod::femtouniverseparticle::ParticleType::kPhi) {
        // Phi histograms
        tempFitVarAxisTitle = "#Phi invariant mass";
      } else {
        LOG(fatal) << "FemtoUniverseParticleHisto: Histogramming for requested object not defined - quitting!";
      }

      framework::AxisSpec tempFitVarpTAxis = {tempFitVarpTBins, "#it{p}_{T} (GeV/#it{c})"}; // the pT binning may vary
      framework::AxisSpec tempFitVarAxis = {tempFitVarBins, tempFitVarAxisTitle};

      std::string folderName = (static_cast<std::string>(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]).c_str() + static_cast<std::string>(mFolderSuffix[mFolderSuffixType])).c_str();

      // Fill here the actual histogramms by calling init_base and init_MC
      init_base<o2::aod::femtouniverseMCparticle::MCType::kRecon>(folderName, tempFitVarAxisTitle, tempFitVarpTAxis, tempFitVarAxis);
      if (isDebug) {
        init_debug<o2::aod::femtouniverseMCparticle::MCType::kRecon>(folderName);
      }
      if (isMC) {
        init_base<o2::aod::femtouniverseMCparticle::MCType::kTruth>(folderName, tempFitVarAxisTitle, tempFitVarpTAxis, tempFitVarAxis);
        init_MC(folderName, tempFitVarAxisTitle, tempFitVarpTAxis, tempFitVarAxis);
      }
    }
  }

  /// Filling of the histograms
  /// Called by init both in case of reconstructed data/ Monte Carlo, and for Monte Carlo Truth
  /// \tparam T Data type of the particle
  /// \param part Particle
  template <o2::aod::femtouniverseMCparticle::MCType mc, typename T>
  void fillQA_base(T const& part)
  {
    /// Histograms of the kinematic properties
    mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/hPt"), part.pt());
    mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/hEta"), part.eta());
    mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/hPhi"), part.phi());
    mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/hPhiEta"), part.phi(), part.eta());

    /// particle specific histogramms for the TempFitVar column in FemtoUniverseParticles
    if constexpr (mc == o2::aod::femtouniverseMCparticle::MCType::kRecon) {
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST(o2::aod::femtouniverseparticle::TempFitVarName[mParticleType]), part.pt(), part.tempFitVar());
    }
  }

  template <o2::aod::femtouniverseMCparticle::MCType mc, typename T>
  void fillQA_debug(T const& part)
  {
    // Histograms holding further debug information
    if constexpr (mParticleType == o2::aod::femtouniverseparticle::ParticleType::kTrack || mParticleType == o2::aod::femtouniverseparticle::ParticleType::kV0Child || mParticleType == o2::aod::femtouniverseparticle::ParticleType::kMCTruthTrack) {
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/hCharge"), part.sign());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/hTPCfindable"), part.tpcNClsFindable());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/hTPCfound"), part.tpcNClsFound());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/hTPCcrossedOverFindable"), part.tpcCrossedRowsOverFindableCls());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/hTPCcrossedRows"), part.tpcNClsCrossedRows());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/hTPCfindableVsCrossed"), part.tpcNClsFindable(), part.tpcNClsCrossedRows());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/hTPCshared"), part.tpcNClsShared());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/hITSclusters"), part.itsNCls());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/hITSclustersIB"), part.itsNClsInnerBarrel());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/hDCAz"), part.pt(), part.dcaZ());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/hDCA"), part.pt(), std::sqrt(std::pow(part.dcaXY(), 2.) + std::pow(part.dcaZ(), 2.)));
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/hTPCdEdX"), part.p(), part.tpcSignal());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/nSigmaTPC_el"), part.p(), part.tpcNSigmaEl());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/nSigmaTPC_pi"), part.p(), part.tpcNSigmaPi());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/nSigmaTPC_K"), part.p(), part.tpcNSigmaKa());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/nSigmaTPC_p"), part.p(), part.tpcNSigmaPr());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/nSigmaTPC_d"), part.p(), part.tpcNSigmaDe());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/nSigmaTOF_el"), part.p(), part.tofNSigmaEl());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/nSigmaTOF_pi"), part.p(), part.tofNSigmaPi());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/nSigmaTOF_K"), part.p(), part.tofNSigmaKa());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/nSigmaTOF_p"), part.p(), part.tofNSigmaPr());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/nSigmaTOF_d"), part.p(), part.tofNSigmaDe());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/nSigmaComb_el"), part.p(), std::sqrt(part.tpcNSigmaEl() * part.tpcNSigmaEl() + part.tofNSigmaEl() * part.tofNSigmaEl()));
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/nSigmaComb_pi"), part.p(), std::sqrt(part.tpcNSigmaPi() * part.tpcNSigmaPi() + part.tofNSigmaPi() * part.tofNSigmaPi()));
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/nSigmaComb_K"), part.p(), std::sqrt(part.tpcNSigmaKa() * part.tpcNSigmaKa() + part.tofNSigmaKa() * part.tofNSigmaKa()));
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/nSigmaComb_p"), part.p(), std::sqrt(part.tpcNSigmaPr() * part.tpcNSigmaPr() + part.tofNSigmaPr() * part.tofNSigmaPr()));
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/nSigmaComb_d"), part.p(), std::sqrt(part.tpcNSigmaDe() * part.tpcNSigmaDe() + part.tofNSigmaDe() * part.tofNSigmaDe()));
    } else if constexpr (mParticleType == o2::aod::femtouniverseparticle::ParticleType::kV0) {
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/hDaughDCA"), part.daughDCA());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/hTransRadius"), part.transRadius());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/hDecayVtxX"), part.decayVtxX());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/hDecayVtxY"), part.decayVtxY());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/hDecayVtxZ"), part.decayVtxZ());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/hInvMassLambda"), part.mLambda());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/hInvMassAntiLambda"), part.mAntiLambda());
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtouniverseMCparticle::MCTypeName[mc]) + HIST("/hInvMassLambdaAntiLambda"), part.mLambda(), part.mAntiLambda());
    }
  }

  /// Filling specialized histograms for Monte Carlo truth
  /// internal function called by init only in case of Monte Carlo truth
  /// \tparam T Data type of the particle
  /// \tparam TMC Data typ of the Monte Carlo Particle
  /// \param part Particle
  /// \param mctruthorigin Origin of the associated mc Truth particle
  /// \param pdgcode PDG of the associated mc Truth particle associated to the reconstructed particle part
  template <typename T>
  void fillQA_MC(T const& part, int mctruthorigin, int pdgcode)
  {
    if (mHistogramRegistry) {
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hPDG"), pdgcode);
      mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hOrigin_MC"), mctruthorigin);

      if (abs(pdgcode) == mPDG) { // fill this histogramm only for TRUE protons (independently of their origin) for the track purity estimation
        mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hPt_ReconNoFake"), part.pt());
      }

      if constexpr (mParticleType == o2::aod::femtouniverseparticle::ParticleType::kTrack || mParticleType == o2::aod::femtouniverseparticle::ParticleType::kV0Child || mParticleType == o2::aod::femtouniverseparticle::ParticleType::kMCTruthTrack) {
        /// Track histograms
        switch (mctruthorigin) {
          case (o2::aod::femtouniverseMCparticle::kPrimary):
            mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_Primary"),
                                     part.pt(), part.tempFitVar());
            break;
          case (o2::aod::femtouniverseMCparticle::kDaughter):
            mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_Daughter"),
                                     part.pt(), part.tempFitVar());
            break;
          case (o2::aod::femtouniverseMCparticle::kMaterial):
            mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_Material"),
                                     part.pt(), part.tempFitVar());
            break;
          case (o2::aod::femtouniverseMCparticle::kFake):
            mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_Fake"),
                                     part.pt(), part.tempFitVar());
            break;
          case (o2::aod::femtouniverseMCparticle::kDaughterLambda):
            mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_DaughterLambda"),
                                     part.pt(), part.tempFitVar());
            break;
          case (o2::aod::femtouniverseMCparticle::kDaughterSigmaplus):
            mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_DaughterSigmaplus"),
                                     part.pt(), part.tempFitVar());
            break;
          default:
            LOG(fatal) << "femtouniverseparticleMC: not known value for ParticleOriginMCTruth - please check. Quitting!";
        }
      } else if constexpr (mParticleType == o2::aod::femtouniverseparticle::ParticleType::kV0) {
        /// V0 histograms
      } else if constexpr (mParticleType == o2::aod::femtouniverseparticle::ParticleType::kCascade) {
        /// Cascade histograms
      } else if constexpr (mParticleType == o2::aod::femtouniverseparticle::ParticleType::kPhi) {
        // Phi histograms
      } else {
        LOG(fatal) << "FemtoUniverseParticleHisto: Histogramming for requested object not defined - quitting!";
      }
    }
  }

  /// Templated function to fill particle histograms for data/ Monte Carlo reconstructed and Monte Carlo truth
  /// Always calls fillQA_base fill histogramms with data/ Monte Carlo reconstructed
  /// In case of Monte Carlo, calls fillQA_base with Monte Carlo truth info and specialized function fillQA_MC for additional histogramms
  /// \tparam T particle type
  /// \tparam isMC fills the additional histograms for Monte Carlo truth
  /// \param part particle for which the histograms should be filled
  template <bool isMC, bool isDebug, typename T>
  void fillQA(T const& part)
  {
    std::string tempFitVarName;
    if (mHistogramRegistry) {
      fillQA_base<o2::aod::femtouniverseMCparticle::MCType::kRecon>(part);
      if constexpr (isDebug) {
        fillQA_debug<o2::aod::femtouniverseMCparticle::MCType::kRecon>(part);
      }
      if constexpr (isMC) {
        if (part.has_fdMCParticle()) {
          fillQA_base<o2::aod::femtouniverseMCparticle::MCType::kTruth>(part.fdMCParticle());
          fillQA_MC(part, (part.fdMCParticle()).partOriginMCTruth(), (part.fdMCParticle()).pdgMCTruth());
        } else {
          mHistogramRegistry->fill(HIST(o2::aod::femtouniverseparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hNoMCtruthCounter"), 0);
        }
      }
    }
  }

 private:
  HistogramRegistry* mHistogramRegistry;                                                      ///< For QA output
  static constexpr o2::aod::femtouniverseparticle::ParticleType mParticleType = particleType; ///< Type of the particle under analysis
  static constexpr int mFolderSuffixType = suffixType;                                        ///< Counter for the folder suffix specified below
  static constexpr std::string_view mFolderSuffix[5] = {"", "_one", "_two", "_pos", "_neg"};  ///< Suffix for the folder name in case of analyses of pairs of the same kind (T-T, V-V, C-C)
  int mPDG = 0;                                                                               ///< PDG code of the selected particle
};
} // namespace o2::analysis::femtoUniverse

#endif // PWGCF_FEMTOUNIVERSE_CORE_FEMTOUNIVERSEPARTICLEHISTO_H_
