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

/// \file dptDptEfficiencyAndQc.cxx
/// \brief Provides efficiency extraction and QC for track cuts and PID
/// \author victor.gonzalez.sebastian@gmail.com

#include <TH2F.h>
#include <TProfile2D.h>
#include <TPDGCode.h>
#include <CCDB/BasicCCDBManager.h>
#include <vector>
#include <cstdio>
#include <memory>
#include <string>
#include "ReconstructionDataFormats/PID.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TableHelper.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/Expressions.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Math/MatrixFunctions.h"
#include "Math/SMatrix.h"

#include "PWGCF/Core/AnalysisConfigurableCuts.h"
#include "PWGCF/DataModel/DptDptFiltered.h"
#include "PWGCF/TableProducer/dptdptfilter.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::framework::expressions;

#define ADDHISTOGRAM(thetype, thedirectory, thename, thetitle, thekind, thebinning...) \
  registry.add<thetype>(TString::Format("%s/%s", thedirectory, thename).Data(), thetitle, thekind, thebinning)
#define FORMATSTRING(theformat, theparams...) TString::Format(theformat, theparams).Data()
#define DIRECTORYSTRING(thedirectoryfmt, thedirectorypars...) FORMATSTRING(thedirectoryfmt, thedirectorypars)
#define HNAMESTRING(thehnamefmt, thehnamepars...) FORMATSTRING(thehnamefmt, thehnamepars)
#define HTITLESTRING(thehtitlefmt, thehtitlepars...) FORMATSTRING(thehtitlefmt, thehtitlepars)

namespace o2::analysis::dptdptfilter
{
TpcExcludeTrack tpcExcluder; ///< the TPC excluder object instance
} // namespace o2::analysis::dptdptfilter

namespace efficiencyandqatask
{
/// \enum KindOfData
/// \brief The kind of data for templating the procedures
enum KindOfData {
  kReco = 0, ///< processing over reconstructed particles/tracks
  kGen       ///< processing over generated particles
};

/// \enum KindOfProcess
/// \brief The kind of processing for templating the procedures and produce histograms
enum KindOfProcess {
  kBASIC,   ///< produce the basic histograms
  kEXTRA,   ///< produce the extra pair based histograms
  kPID,     ///< produce the basic PID histograms
  kPIDEXTRA ///< produce the extra PID histograms
};

/// \enum BeforeAfter
/// \brief The kind of filling, before or after track selection
enum BeforeAfter {
  kBefore = 0, ///< filling before track selection
  kAfter       ///< filling after track selection
};

/* the structures for checking the TPC sector borders impact */
constexpr int kNoOfTpcSectors = 18;
constexpr float kTpcPhiSectorWidth = (constants::math::TwoPI) / kNoOfTpcSectors;

/* the configuration of the nsigma axis */
float minNSigma = -4.05f;
float maxNSigma = 4.05f;
float widthNSigmaBin = 0.1f;
int noOfNSigmaBins = static_cast<int>((maxNSigma - minNSigma) / widthNSigmaBin);

/* the PID selector object to help with the configuration and the id of the selected particles */
o2::analysis::dptdptfilter::PIDSpeciesSelection pidselector;

// initialized during self configuration
std::vector<std::string> poinames; ///< the species of interest names
std::vector<std::string> tnames;   ///< the track names

static const std::vector<o2::track::PID::ID> allmainspecies{o2::track::PID::Electron, o2::track::PID::Muon, o2::track::PID::Pion, o2::track::PID::Kaon, o2::track::PID::Proton};
static const std::vector<std::string> allmainspnames{"ElectronP", "ElectronM", "MuonP", "MuonM", "PionP", "PionM", "KaonP", "KaonM", "ProtonP", "ProtonM"};
static const std::vector<std::string> allmainsptitles{"e^{#plus}", "e^{#minus}", "#mu^{#plus}", "#mu^{#minus}", "#pi^{#plus}", "#pi^{#minus}", "K^{#plus}", "K^{#minus}", "p", "#bar{p}"};
static const std::vector<o2::track::PID::ID> mainspecies{o2::track::PID::Pion, o2::track::PID::Kaon, o2::track::PID::Proton};
static const std::vector<std::string> mainspnames{"PionP", "PionM", "KaonP", "KaonM", "ProtonP", "ProtonM"};
static const std::vector<std::string> mainsptitles{"#pi^{#plus}", "#pi^{#minus}", "K^{#plus}", "K^{#minus}", "p", "#bar{p}"};
static const std::vector<int> pdgcodes = {kElectron, kMuonPlus, kPiPlus, kKPlus, kProton};
} // namespace efficiencyandqatask

/* the QA data collecting engine */
struct QADataCollectingEngine {
  uint nsp = static_cast<uint>(efficiencyandqatask::tnames.size());
  uint nmainsp = static_cast<uint>(efficiencyandqatask::mainspnames.size());
  uint nallmainsp = static_cast<uint>(efficiencyandqatask::allmainspnames.size());

  //===================================================
  // The QA output objects
  //===================================================
  /* momentum histograms */
  std::shared_ptr<TH2> fhPvsInnerP = nullptr;
  std::shared_ptr<TH2> fhTruePvsP = nullptr;
  std::shared_ptr<TH2> fhTruePvsInnerP = nullptr;
  /* efficiency histograms histograms */
  /* when two indexes, first index reco and detector level, second index generator level */
  /* when no indexes, reco and detector level */
  std::vector<std::shared_ptr<TH1>> fhPtB{2, nullptr};
  std::vector<std::shared_ptr<TH2>> fhPtVsEtaB{2, nullptr};
  std::vector<std::shared_ptr<TH2>> fhPtVsZvtxB{2, nullptr};
  std::shared_ptr<TH2> fhPhiVsPtPosB{nullptr};
  std::shared_ptr<TH3> fhNchVsPhiVsPtPosB{nullptr};
  TH2* fhPerColNchVsPhiVsPtPosB{nullptr};
  std::shared_ptr<TH2> fhPhiVsInnerWallMomPosB{nullptr};
  std::shared_ptr<TH3> fhNchVsPhiVsInnerWallMomPosB{nullptr};
  TH2* fhPerColNchVsPhiVsInnerWallMomPosB{nullptr};
  std::shared_ptr<TH2> fhPhiVsPtNegB{nullptr};
  std::shared_ptr<TH3> fhNchVsPhiVsPtNegB{nullptr};
  TH2* fhPerColNchVsPhiVsPtNegB{nullptr};
  std::shared_ptr<TH2> fhPhiVsInnerWallMomNegB{nullptr};
  std::shared_ptr<TH3> fhNchVsPhiVsInnerWallMomNegB{nullptr};
  TH2* fhPerColNchVsPhiVsInnerWallMomNegB{nullptr};
  std::vector<std::vector<std::shared_ptr<TH1>>> fhPtA{2, {nsp, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPtVsEtaA{2, {nsp, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPtVsZvtxA{2, {nsp, nullptr}};
  std::vector<std::shared_ptr<TH2>> fhPhiVsPtA{nsp, nullptr};
  std::vector<std::shared_ptr<TH3>> fhNchVsPhiVsPtA{nsp, nullptr};
  std::vector<TH2*> fhPerColNchVsPhiVsPtA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhPhiVsInnerWallMomA{nsp, nullptr};
  std::vector<std::shared_ptr<TH3>> fhNchVsPhiVsInnerWallMomA{nsp, nullptr};
  std::vector<TH2*> fhPerColNchVsPhiVsInnerWallMomA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhPhiShiftedVsPtA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhPhiShiftedVsInnerWallMomA{nsp, nullptr};
  std::shared_ptr<TH2> fhPtVsEtaItsAcc{nullptr};
  std::shared_ptr<TH2> fhPtVsEtaTpcAcc{nullptr};
  std::shared_ptr<TH2> fhPtVsEtaItsTpcAcc{nullptr};
  std::shared_ptr<TH2> fhPtVsEtaItsTofAcc{nullptr};
  std::shared_ptr<TH2> fhPtVsEtaTpcTofAcc{nullptr};
  std::shared_ptr<TH2> fhPtVsEtaItsTpcTofAcc{nullptr};
  std::vector<std::shared_ptr<TH2>> fhPtVsEtaItsA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhPtVsEtaTpcA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhPtVsEtaItsTpcA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhPtVsEtaItsTofA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhPtVsEtaTpcTofA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhPtVsEtaItsTpcTofA{nsp, nullptr};
  /* primaries and secondaries */
  /* overall, first index detector level second index generator level */
  /* detailed, first index detector level, second index associated particle */
  std::shared_ptr<TH3> fhPtPurityPosPrimA{nullptr};
  std::shared_ptr<TH3> fhPtPurityNegPrimA{nullptr};
  std::vector<std::shared_ptr<TH2>> fhPtVsEtaPrimA{nsp, nullptr};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPtVsEtaPrimItsA{2, {nsp, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPtVsEtaPrimItsTpcA{2, {nsp, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPtVsEtaPrimItsTpcTofA{2, {nsp, nullptr}};
  std::shared_ptr<TH3> fhPtPurityPosSecA{nullptr};
  std::shared_ptr<TH3> fhPtPurityNegSecA{nullptr};
  std::vector<std::shared_ptr<TH2>> fhPtVsEtaSecA{nsp, nullptr};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPtVsEtaSecItsA{2, {nsp, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPtVsEtaSecItsTpcA{2, {nsp, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPtVsEtaSecItsTpcTofA{2, {nsp, nullptr}};
  std::shared_ptr<TH3> fhPtPurityPosMatA{nullptr};
  std::shared_ptr<TH3> fhPtPurityNegMatA{nullptr};
  std::vector<std::shared_ptr<TH2>> fhPtVsEtaMatA{nsp, nullptr};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPtVsEtaMatItsA{2, {nsp, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPtVsEtaMatItsTpcA{2, {nsp, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPtVsEtaMatItsTpcTofA{2, {nsp, nullptr}};
  /* QC histograms */
  std::shared_ptr<TH2> fhItsNClsVsPtB{nullptr};
  std::shared_ptr<TH2> fhItsChi2NClsVsPtB{nullptr};
  std::shared_ptr<TH2> fhTpcFindableNClsVsPtB{nullptr};
  std::shared_ptr<TH2> fhTpcFoundNClsVsPtB{nullptr};
  std::shared_ptr<TH2> fhTpcSharedNClsVsPtB{nullptr};
  std::shared_ptr<TH2> fhTpcFractionSharedClsVsPtB{nullptr};
  std::shared_ptr<TH2> fhTpcCrossedRowsVsPtB{nullptr};
  std::shared_ptr<TH2> fhTpcCrossedRowsOverFindableClsVsPtB{nullptr};
  std::shared_ptr<TH2> fhTpcChi2NClsVsPtB{nullptr};
  std::vector<std::shared_ptr<TH2>> fhItsNClsVsPtA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhItsChi2NClsVsPtA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhTpcFindableNClsVsPtA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhTpcFoundNClsVsPtA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhTpcSharedNClsVsPtA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhTpcFractionSharedClsVsPtA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhTpcCrossedRowsVsPtA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhTpcCrossedRowsOverFindableClsVsPtA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhTpcChi2NClsVsPtA{nsp, nullptr};

  template <efficiencyandqatask::KindOfData kindOfData>
  void init(HistogramRegistry& registry, const char* dirname)
  {
    using namespace efficiencyandqatask;
    using namespace analysis::dptdptfilter;

    AxisSpec pidPtAxis{150, 0.1, 5.0, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pidPtAxisReduced{50, 0.1, 5.0, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pidPAxis{150, 0.1, 5.0, "#it{p} (GeV/#it{c})"};
    AxisSpec pidPAxisReduced{50, 0.1, 5.0, "#it{p} (GeV/#it{c})"};
    pidPtAxis.makeLogarithmic();
    pidPAxis.makeLogarithmic();
    const AxisSpec ptAxis{ptbins, ptlow, ptup, "#it{p}_{T} (GeV/c)"};
    const AxisSpec etaAxis{etabins, etalow, etaup, "#eta"};
    const AxisSpec phiAxis{360, 0.0f, constants::math::TwoPI, "#varphi (rad)"};
    const AxisSpec phiSectorAxis{144, 0.0f, 0.36f, "#varphi (mod(2#pi/18) (rad))"};
    const AxisSpec phiSectorAxisReduced{36, 0.0f, 0.36f, "#varphi (mod(2#pi/18) (rad))"};
    const AxisSpec nChargeAxis{100, 0.0f, 100.0f, "#it{N}_{ch}"};
    const AxisSpec phiShiftedSectorAxis{220, -55.0f, 55.0f, "% of the sector"};
    const AxisSpec zvtxAxis{zvtxbins, zvtxlow, zvtxup, "#it{z}_{vtx}"};
    const AxisSpec itsNClsAxis{8, -0.5, 7.5, "ITS n clusters"};
    const AxisSpec itsCh2Axis{100, 0, 40, "#Chi^{2}/Cls ITS"};
    const AxisSpec tpcNClsAxis{165, -0.5, 164.5, "TPC n clusters"};
    const AxisSpec tpcNRowsAxis{165, -0.5, 164.5, "TPC n rows"};
    const AxisSpec tpcFractionAxis{100, 0, 1, "fraction"};
    const AxisSpec tpcXRowsOverFindClsAxis{60, 0.7, 1.3, "fraction"};
    const AxisSpec tpcCh2Axis{100, 0, 10, "#Chi^{2}/Cls TPC"};

    /* the reconstructed and generated levels histograms */
    std::string recogen = (kindOfData == kReco) ? "Reco" : "Gen";
    fhPtB[kindOfData] = ADDHISTOGRAM(TH1, DIRECTORYSTRING("%s/%s/%s", dirname, recogen.c_str(), "Before"), "Pt", "#it{p}_{T}", kTH1F, {ptAxis});
    fhPtVsEtaB[kindOfData] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, recogen.c_str(), "Before"), "PtVsEta", "#it{p}_T vs #eta", kTH2F, {etaAxis, ptAxis});
    fhPtVsZvtxB[kindOfData] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, recogen.c_str(), "Before"), "PtVsZvtx", "#it{p}_T vs #it{z}_{vtx}", kTH2F, {zvtxAxis, ptAxis});
    for (uint isp = 0; isp < nsp; ++isp) {
      fhPtA[kindOfData][isp] = ADDHISTOGRAM(TH1, DIRECTORYSTRING("%s/%s/%s", dirname, recogen.c_str(), "After"), HNAMESTRING("Pt_%s", tnames[isp].c_str()), HTITLESTRING("#it{p}_{T} %s", tnames[isp].c_str()), kTH1F, {ptAxis});
      fhPtVsEtaA[kindOfData][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, recogen.c_str(), "After"), HNAMESTRING("PtVsEta_%s", tnames[isp].c_str()), HTITLESTRING("#it{p}_{T} vs #eta %s", tnames[isp].c_str()), kTH2F, {etaAxis, ptAxis});
      fhPtVsZvtxA[kindOfData][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, recogen.c_str(), "After"), HNAMESTRING("PtVsZvtx_%s", tnames[isp].c_str()), HTITLESTRING("#it{p}_{T} vs #it{z}_{zvtx} %s", tnames[isp].c_str()), kTH2F, {zvtxAxis, ptAxis});
    }

    if constexpr (kindOfData == kReco) {
      /* only the reconstructed level histograms*/
      fhPhiVsPtPosB = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "PhiVsPtPos", "#varphi (mod(2#pi/18))", kTH2F, {pidPtAxis, phiSectorAxis});
      fhNchVsPhiVsPtPosB = ADDHISTOGRAM(TH3, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "NchVsPhiVsPtPos", "#it{N}_{ch}^{#plus} #varphi (mod(2#pi/18))", kTH3F, {pidPtAxisReduced, phiSectorAxisReduced, nChargeAxis});
      fhPhiVsInnerWallMomPosB = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "PhiVsIwMomPos", "#varphi (mod(2#pi/18)) TPC_{iw} #it{p}", kTH2F, {pidPAxis, phiSectorAxis});
      fhNchVsPhiVsInnerWallMomPosB = ADDHISTOGRAM(TH3, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "NchVsPhiVsIwMomPos", "#it{N}_{ch}^{#plus} #varphi (mod(2#pi/18)) TPC_{iw} #it{p}", kTH3F, {pidPAxisReduced, phiSectorAxisReduced, nChargeAxis});
      fhPhiVsPtNegB = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "PhiVsPtNeg", "#varphi (mod(2#pi/18))", kTH2F, {pidPtAxis, phiSectorAxis});
      fhNchVsPhiVsPtNegB = ADDHISTOGRAM(TH3, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "NchVsPhiVsPtNeg", "#it{N}_{ch}^{#minus} #varphi (mod(2#pi/18))", kTH3F, {pidPtAxisReduced, phiSectorAxisReduced, nChargeAxis});
      fhPhiVsInnerWallMomNegB = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "PhiVsIwMomNeg", "#varphi (mod(2#pi/18)) TPC_{iw} #it{p}", kTH2F, {pidPAxis, phiSectorAxis});
      fhNchVsPhiVsInnerWallMomNegB = ADDHISTOGRAM(TH3, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "NchVsPhiVsIwMomNeg", "#it{N}_{ch}^{#minus} #varphi (mod(2#pi/18)) TPC_{iw} #it{p}", kTH3F, {pidPAxisReduced, phiSectorAxisReduced, nChargeAxis});
      fhItsNClsVsPtB = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "ITSNCls", "ITS clusters", kTH2F, {ptAxis, itsNClsAxis});
      fhItsChi2NClsVsPtB = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "ITSChi2NCls", "ITS #Chi^{2}", kTH2F, {ptAxis, itsCh2Axis});
      fhTpcFindableNClsVsPtB = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "TPCFindableNCls", "TPC findable clusters", kTH2F, {ptAxis, tpcNClsAxis});
      fhTpcFoundNClsVsPtB = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "TPCFoundNCls", "TPC found clusters", kTH2F, {ptAxis, tpcNClsAxis});
      fhTpcSharedNClsVsPtB = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "TPCSharedNCls", "TPC shared clusters", kTH2F, {ptAxis, tpcNClsAxis});
      fhTpcFractionSharedClsVsPtB = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "TPCFractionSharedCls", "TPC fraction shared clusters", kTH2F, {ptAxis, tpcFractionAxis});
      fhTpcCrossedRowsVsPtB = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "TPCXrows", "TPC crossed rows", kTH2F, {ptAxis, tpcNRowsAxis});
      fhTpcCrossedRowsOverFindableClsVsPtB = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "XRowsOverFindableCls", "TPC xrows over findable clusters", kTH2F, {ptAxis, tpcXRowsOverFindClsAxis});
      fhTpcChi2NClsVsPtB = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "TPCChi2NCls", "TPC #Chi^{2}", kTH2F, {ptAxis, tpcCh2Axis});
      /* efficiency histograms */
      fhPvsInnerP = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), "pVsInnerP", "#it{p} versus TPC inner wall #it{p}", kTH2F, {pidPAxis, pidPAxis});
      fhPtVsEtaItsAcc = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), "ptItsAcc", "ITS tracks within the acceptance", kTH2F, {etaAxis, ptAxis});
      fhPtVsEtaTpcAcc = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), "ptTpcAcc", "TPC tracks within the acceptance", kTH2F, {etaAxis, ptAxis});
      fhPtVsEtaItsTpcAcc = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), "ptItsTpcAcc", "ITS&TPC tracks within the acceptance", kTH2F, {etaAxis, ptAxis});
      fhPtVsEtaItsTofAcc = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), "ptItsTofAcc", "ITS&TOF tracks within the acceptance", kTH2F, {etaAxis, ptAxis});
      fhPtVsEtaTpcTofAcc = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), "ptTpcTofAcc", "TPC&TOF tracks within the acceptance", kTH2F, {etaAxis, ptAxis});
      fhPtVsEtaItsTpcTofAcc = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), "ptItsTpcTofAcc", "ITS&TPC&TOF tracks within the acceptance", kTH2F, {etaAxis, ptAxis});
      /* per collision histograms not going to the results file */
      int nPtBins = fhNchVsPhiVsPtPosB->GetNbinsX();
      float ptLow = fhNchVsPhiVsPtNegB->GetXaxis()->GetBinLowEdge(1);
      float ptHigh = fhNchVsPhiVsPtNegB->GetXaxis()->GetBinUpEdge(nPtBins);
      int nTpcIwMomBins = fhNchVsPhiVsInnerWallMomNegB->GetNbinsX();
      float tpcIwMomLow = fhNchVsPhiVsInnerWallMomNegB->GetXaxis()->GetBinLowEdge(1);
      float tpcIwMomHigh = fhNchVsPhiVsInnerWallMomNegB->GetXaxis()->GetBinUpEdge(nTpcIwMomBins);
      int nPhiSectorBins = fhNchVsPhiVsPtPosB->GetNbinsY();
      float phiSectorLow = fhNchVsPhiVsPtNegB->GetYaxis()->GetBinLowEdge(1);
      float phiSectorHigh = fhNchVsPhiVsPtNegB->GetYaxis()->GetBinUpEdge(nPhiSectorBins);
      fhPerColNchVsPhiVsPtPosB = new TH2F(TString::Format("%s_PerColNchVsPhiVsPtPosB", dirname), "", nPtBins, ptLow, ptHigh, nPhiSectorBins, phiSectorLow, phiSectorHigh);
      fhPerColNchVsPhiVsInnerWallMomPosB = new TH2F(TString::Format("%s_PerColNchVsPhiVsInnerWallMomPosB", dirname), "", nTpcIwMomBins, tpcIwMomLow, tpcIwMomHigh, nPhiSectorBins, phiSectorLow, phiSectorHigh);
      fhPerColNchVsPhiVsPtNegB = new TH2F(TString::Format("%s_PerColNchVsPhiVsPtNegB", dirname), "", nPtBins, ptLow, ptHigh, nPhiSectorBins, phiSectorLow, phiSectorHigh);
      fhPerColNchVsPhiVsInnerWallMomNegB = new TH2F(TString::Format("%s_PerColNchVsPhiVsInnerWallMomNegB", dirname), "", nTpcIwMomBins, tpcIwMomLow, tpcIwMomHigh, nPhiSectorBins, phiSectorLow, phiSectorHigh);
      for (uint isp = 0; isp < nsp; ++isp) {
        fhPhiVsPtA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "After"), HNAMESTRING("PhiVsPt_%s", tnames[isp].c_str()), HTITLESTRING("#varphi %s (mod(2#pi/18))", tnames[isp].c_str()), kTH2F, {pidPtAxis, phiSectorAxis});
        fhNchVsPhiVsPtA[isp] = ADDHISTOGRAM(TH3, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "After"), HNAMESTRING("NchVsPhiVsPt_%s", tnames[isp].c_str()), HTITLESTRING("#it{N}_{ch}^{%s} #varphi (mod(2#pi/18))", tnames[isp].c_str()), kTH3F, {pidPtAxisReduced, phiSectorAxisReduced, nChargeAxis});
        fhPhiVsInnerWallMomA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "After"), HNAMESTRING("PhiVsIwMom_%s", tnames[isp].c_str()), HTITLESTRING("#varphi %s (mod(2#pi/18)) TPC_{iw} #it{p}", tnames[isp].c_str()), kTH2F, {pidPAxis, phiSectorAxis});
        fhNchVsPhiVsInnerWallMomA[isp] = ADDHISTOGRAM(TH3, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "After"), HNAMESTRING("NchVsPhiVsIwMom_%s", tnames[isp].c_str()), HTITLESTRING("#it{N}_{ch}^{%s} #varphi (mod(2#pi/18)) TPC_{iw} #it{p}", tnames[isp].c_str()), kTH3F, {pidPAxisReduced, phiSectorAxisReduced, nChargeAxis});
        fhPhiShiftedVsPtA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "After"), HNAMESTRING("PhiShiftedVsPt_%s", tnames[isp].c_str()), HTITLESTRING("%s TPC sector %%", tnames[isp].c_str()), kTH2F, {pidPtAxis, phiShiftedSectorAxis});
        fhPhiShiftedVsInnerWallMomA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "After"), HNAMESTRING("PhiShiftedVsIwMom_%s", tnames[isp].c_str()), HTITLESTRING("%s TPC sector %% TPC_{iw} #it{p}", tnames[isp].c_str()), kTH2F, {pidPAxis, phiShiftedSectorAxis});
        fhItsNClsVsPtA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "After"), HNAMESTRING("ITSNCls_%s", tnames[isp].c_str()), HTITLESTRING("ITS clusters %s", tnames[isp].c_str()), kTH2F, {ptAxis, itsNClsAxis});
        fhItsChi2NClsVsPtA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "After"), HNAMESTRING("ITSChi2NCls_%s", tnames[isp].c_str()), HTITLESTRING("ITS #Chi^{2} %s", tnames[isp].c_str()), kTH2F, {ptAxis, itsCh2Axis});
        fhTpcFindableNClsVsPtA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "After"), HNAMESTRING("TPCFindableNCls_%s", tnames[isp].c_str()), HTITLESTRING("TPC findable clusters %s", tnames[isp].c_str()), kTH2F, {ptAxis, tpcNClsAxis});
        fhTpcFoundNClsVsPtA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "After"), HNAMESTRING("TPCFoundNCls_%s", tnames[isp].c_str()), HTITLESTRING("TPC found clusters %s", tnames[isp].c_str()), kTH2F, {ptAxis, tpcNClsAxis});
        fhTpcSharedNClsVsPtA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "After"), HNAMESTRING("TPCSharedNCls_%s", tnames[isp].c_str()), HTITLESTRING("TPC shared clusters %s", tnames[isp].c_str()), kTH2F, {ptAxis, tpcNClsAxis});
        fhTpcFractionSharedClsVsPtA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "After"), HNAMESTRING("TPCFractionSharedCls_%s", tnames[isp].c_str()), HTITLESTRING("TPC fraction shared clusters %s", tnames[isp].c_str()), kTH2F, {ptAxis, tpcFractionAxis});
        fhTpcCrossedRowsVsPtA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "After"), HNAMESTRING("TPCXrows_%s", tnames[isp].c_str()), HTITLESTRING("TPC crossed rows %s", tnames[isp].c_str()), kTH2F, {ptAxis, tpcNRowsAxis});
        fhTpcCrossedRowsOverFindableClsVsPtA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "After"), HNAMESTRING("XRowsOverFindableCls_%s", tnames[isp].c_str()), HTITLESTRING("TPC xrows over findable clusters %s", tnames[isp].c_str()), kTH2F, {ptAxis, tpcXRowsOverFindClsAxis});
        fhTpcChi2NClsVsPtA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "After"), HNAMESTRING("TPCChi2NCls_%s", tnames[isp].c_str()), HTITLESTRING("TPC #Chi^{2} %s", tnames[isp].c_str()), kTH2F, {ptAxis, tpcCh2Axis});
        /* efficiency histograms */
        fhPtVsEtaItsA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), HNAMESTRING("ptIts%s", tnames[isp].c_str()), HTITLESTRING("ITS %s tracks", tnames[isp].c_str()), kTH2F, {etaAxis, ptAxis});
        fhPtVsEtaTpcA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), HNAMESTRING("ptTpc%s", tnames[isp].c_str()), HTITLESTRING("TPC %s tracks", tnames[isp].c_str()), kTH2F, {etaAxis, ptAxis});
        fhPtVsEtaItsTpcA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), HNAMESTRING("ptItsTpc%s", tnames[isp].c_str()), HTITLESTRING("ITS&TPC %s tracks", tnames[isp].c_str()), kTH2F, {etaAxis, ptAxis});
        fhPtVsEtaItsTofA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), HNAMESTRING("ptItsTof_%s", tnames[isp].c_str()), HTITLESTRING("ITS&TOF %s tracks", tnames[isp].c_str()), kTH2F, {etaAxis, ptAxis});
        fhPtVsEtaTpcTofA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), HNAMESTRING("ptTpcTof_%s", tnames[isp].c_str()), HTITLESTRING("TPC&TOF %s tracks", tnames[isp].c_str()), kTH2F, {etaAxis, ptAxis});
        fhPtVsEtaItsTpcTofA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), HNAMESTRING("ptItsTpcTof_%s", tnames[isp].c_str()), HTITLESTRING("ITS&TPC&TOF %s tracks", tnames[isp].c_str()), kTH2F, {etaAxis, ptAxis});
        /* per collision histograms not going to the results file */
        fhPerColNchVsPhiVsPtA[isp] = new TH2F(HNAMESTRING("%s_PerColNchVsPhiVsPt_%s", dirname, tnames[isp].c_str()), "", nPtBins, ptLow, ptHigh, nPhiSectorBins, phiSectorLow, phiSectorHigh);
        fhPerColNchVsPhiVsInnerWallMomA[isp] = new TH2F(HNAMESTRING("%s_PerColNchVsPhiVsInnerWallMom_%s", dirname, tnames[isp].c_str()), "", nTpcIwMomBins, tpcIwMomLow, tpcIwMomHigh, nPhiSectorBins, phiSectorLow, phiSectorHigh);
      }
    } else {
      AxisSpec recoSpecies{static_cast<int>(nsp) + 1, -0.5, nsp - 0.5, "reco species"};
      AxisSpec trueSpecies{static_cast<int>(nmainsp) + 1, -0.5, nmainsp + 0.5, "true species"};
      fhTruePvsP = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Gen"), "truePVsP", "#it{p} gen versus reco #it{p}", kTH2F, {pidPAxis, pidPAxis});
      fhTruePvsInnerP = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Gen"), "truePVsInnerP", "#it{p} gen versus reco TPC inner wall #it{p}", kTH2F, {pidPAxis, pidPAxis});
      fhPtPurityPosPrimA = ADDHISTOGRAM(TH3, DIRECTORYSTRING("%s/%s", dirname, "Purity"), "ptPurityPosPrim", "Primaries for reconstructed positive", kTH3F, {recoSpecies, trueSpecies, ptAxis});
      fhPtPurityNegPrimA = ADDHISTOGRAM(TH3, DIRECTORYSTRING("%s/%s", dirname, "Purity"), "ptPurityNegPrim", "Primaries for reconstructed negative", kTH3F, {recoSpecies, trueSpecies, ptAxis});
      fhPtPurityPosSecA = ADDHISTOGRAM(TH3, DIRECTORYSTRING("%s/%s", dirname, "Purity"), "ptPurityPosSec", "Secondaries for reconstructed positive", kTH3F, {recoSpecies, trueSpecies, ptAxis});
      fhPtPurityNegSecA = ADDHISTOGRAM(TH3, DIRECTORYSTRING("%s/%s", dirname, "Purity"), "ptPurityNegSec", "Secondaries for reconstructed negative", kTH3F, {recoSpecies, trueSpecies, ptAxis});
      fhPtPurityPosMatA = ADDHISTOGRAM(TH3, DIRECTORYSTRING("%s/%s", dirname, "Purity"), "ptPurityPosMat", "Secondaries from material for reconstructed positive", kTH3F, {recoSpecies, trueSpecies, ptAxis});
      fhPtPurityNegMatA = ADDHISTOGRAM(TH3, DIRECTORYSTRING("%s/%s", dirname, "Purity"), "ptPurityNegMat", "Secondaries from material for reconstructed negative", kTH3F, {recoSpecies, trueSpecies, ptAxis});
      for (uint isp = 0; isp < nsp; ++isp) {
        /* detector level and generator level histograms */
        fhPtVsEtaPrimA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Gen"),
                                           HNAMESTRING("ptPrim%s", tnames[isp].c_str()),
                                           HTITLESTRING("ITS  %s tracks (primaries)", tnames[isp].c_str()),
                                           kTH2F, {etaAxis, ptAxis});
        fhPtVsEtaSecA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Gen"),
                                          HNAMESTRING("ptSec%s", tnames[isp].c_str()),
                                          HTITLESTRING("ITS %s tracks (secondaries)", tnames[isp].c_str()),
                                          kTH2F, {etaAxis, ptAxis});
        fhPtVsEtaMatA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Gen"),
                                          HNAMESTRING("ptMat%s", tnames[isp].c_str()),
                                          HTITLESTRING("ITS %s tracks (from material)", tnames[isp].c_str()),
                                          kTH2F, {etaAxis, ptAxis});

        const std::vector<std::string> detectedorigin = {"DetReco", "DetAssoc"};
        for (uint ix = 0; ix < detectedorigin.size(); ++ix) {
          fhPtVsEtaPrimItsA[ix][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", detectedorigin[ix].c_str()),
                                                    HNAMESTRING("ptItsPrim_%s", tnames[isp].c_str()),
                                                    HTITLESTRING("ITS %s tracks (primaries)", tnames[isp].c_str()),
                                                    kTH2F, {etaAxis, ptAxis});
          fhPtVsEtaPrimItsTpcA[ix][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", detectedorigin[ix].c_str()),
                                                       HNAMESTRING("ptItsTpcPrim_%s", tnames[isp].c_str()),
                                                       HTITLESTRING("ITS&TPC %s tracks (primaries)", tnames[isp].c_str()),
                                                       kTH2F, {etaAxis, ptAxis});
          fhPtVsEtaPrimItsTpcTofA[ix][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", detectedorigin[ix].c_str()),
                                                          HNAMESTRING("ptItsTpcTofPrim_%s", tnames[isp].c_str()),
                                                          HTITLESTRING("ITS&TPC&TOF %s tracks (primaries)", tnames[isp].c_str()),
                                                          kTH2F, {etaAxis, ptAxis});
          fhPtVsEtaSecItsA[ix][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", detectedorigin[ix].c_str()),
                                                   HNAMESTRING("ptItsSec_%s", tnames[isp].c_str()),
                                                   HTITLESTRING("ITS %s tracks (secondaries)", tnames[isp].c_str()),
                                                   kTH2F, {etaAxis, ptAxis});
          fhPtVsEtaSecItsTpcA[ix][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", detectedorigin[ix].c_str()),
                                                      HNAMESTRING("ptItsTpcSec_%s", tnames[isp].c_str()),
                                                      HTITLESTRING("ITS&TPC %s tracks (secondaries)", tnames[isp].c_str()),
                                                      kTH2F, {etaAxis, ptAxis});
          fhPtVsEtaSecItsTpcTofA[ix][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", detectedorigin[ix].c_str()),
                                                         HNAMESTRING("ptItsTpcTofSec_%s", tnames[isp].c_str()),
                                                         HTITLESTRING("ITS&TPC&TOF %s tracks (secondaries)", tnames[isp].c_str()),
                                                         kTH2F, {etaAxis, ptAxis});
          fhPtVsEtaMatItsA[ix][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", detectedorigin[ix].c_str()),
                                                   HNAMESTRING("ptItsMat_%s", tnames[isp].c_str()),
                                                   HTITLESTRING("ITS %s tracks (from material)", tnames[isp].c_str()),
                                                   kTH2F, {etaAxis, ptAxis});
          fhPtVsEtaMatItsTpcA[ix][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", detectedorigin[ix].c_str()),
                                                      HNAMESTRING("ptItsTpcMat_%s", tnames[isp].c_str()),
                                                      HTITLESTRING("ITS&TPC %s tracks (from material)", tnames[isp].c_str()),
                                                      kTH2F, {etaAxis, ptAxis});
          fhPtVsEtaMatItsTpcTofA[ix][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", detectedorigin[ix].c_str()),
                                                         HNAMESTRING("ptItsTpcTofMat_%s", tnames[isp].c_str()),
                                                         HTITLESTRING("ITS&TPC&TOF %s tracks (from material)", tnames[isp].c_str()),
                                                         kTH2F, {etaAxis, ptAxis});
        }
      }
    }
  }

  template <efficiencyandqatask::KindOfData kindOfData, typename CollisionsObject, typename TrackObject>
  void processTrack(float zvtx, TrackObject const& track)
  {
    using namespace efficiencyandqatask;
    using namespace analysis::dptdptfilter;
    using namespace o2::aod::track;

    fhPtB[kindOfData]->Fill(track.pt());
    fhPtVsEtaB[kindOfData]->Fill(track.eta(), track.pt());
    fhPtVsZvtxB[kindOfData]->Fill(zvtx, track.pt());
    if (!(track.trackacceptedid() < 0)) {
      fhPtA[kindOfData][track.trackacceptedid()]->Fill(track.pt());
      fhPtVsEtaA[kindOfData][track.trackacceptedid()]->Fill(track.eta(), track.pt());
      fhPtVsZvtxA[kindOfData][track.trackacceptedid()]->Fill(zvtx, track.pt());
    }
    if constexpr (kindOfData == kReco) {
      auto fillhisto = [&track](auto& h, bool cond) {
        if (cond) {
          h->Fill(track.eta(), track.pt());
        }
      };
      bool hasits = track.hasITS() && TrackSelectionFlags::checkFlag(track.trackCutFlag(), TrackSelectionITS);
      bool hastpc = track.hasTPC() && TrackSelectionFlags::checkFlag(track.trackCutFlag(), TrackSelectionTPC);
      bool hastof = track.hasTOF();

      float phiInTpcSector = std::fmod(track.phi(), kTpcPhiSectorWidth);
      float phiShiftedPercentInTpcSector = phiInTpcSector * 100 / kTpcPhiSectorWidth;
      phiShiftedPercentInTpcSector = (phiShiftedPercentInTpcSector > 50.0f) ? (phiShiftedPercentInTpcSector - 100.0f) : phiShiftedPercentInTpcSector;
      if (track.sign() > 0) {
        fhPhiVsPtPosB->Fill(track.pt(), phiInTpcSector);
        fhPerColNchVsPhiVsPtPosB->Fill(track.pt(), phiInTpcSector);
        fhPhiVsInnerWallMomPosB->Fill(track.tpcInnerParam(), phiInTpcSector);
        fhPerColNchVsPhiVsInnerWallMomPosB->Fill(track.tpcInnerParam(), phiInTpcSector);
      } else {
        fhPhiVsPtNegB->Fill(track.pt(), phiInTpcSector);
        fhPerColNchVsPhiVsPtNegB->Fill(track.pt(), phiInTpcSector);
        fhPhiVsInnerWallMomNegB->Fill(track.tpcInnerParam(), phiInTpcSector);
        fhPerColNchVsPhiVsInnerWallMomNegB->Fill(track.tpcInnerParam(), phiInTpcSector);
      }
      fhItsNClsVsPtB->Fill(track.pt(), track.itsNCls());
      fhItsChi2NClsVsPtB->Fill(track.pt(), track.itsChi2NCl());
      fhTpcFindableNClsVsPtB->Fill(track.pt(), track.tpcNClsFindable());
      fhTpcFoundNClsVsPtB->Fill(track.pt(), track.tpcNClsFound());
      fhTpcSharedNClsVsPtB->Fill(track.pt(), track.tpcNClsShared());
      fhTpcFractionSharedClsVsPtB->Fill(track.pt(), track.tpcFractionSharedCls());
      fhTpcCrossedRowsVsPtB->Fill(track.pt(), track.tpcNClsCrossedRows());
      fhTpcCrossedRowsOverFindableClsVsPtB->Fill(track.pt(), track.tpcCrossedRowsOverFindableCls());
      fhTpcChi2NClsVsPtB->Fill(track.pt(), track.tpcChi2NCl());
      if (inTheAcceptance<CollisionsObject>(track)) {
        /* efficiency histograms */
        fillhisto(fhPtVsEtaItsAcc, hasits);
        fillhisto(fhPtVsEtaTpcAcc, hastpc);
        fillhisto(fhPtVsEtaItsTpcAcc, hasits && hastpc);
        fillhisto(fhPtVsEtaItsTofAcc, hasits && hastof);
        fillhisto(fhPtVsEtaTpcTofAcc, hastpc && hastof);
        fillhisto(fhPtVsEtaItsTpcTofAcc, hasits && hastpc && hastof);
      }
      if (!(track.trackacceptedid() < 0)) {
        fhPhiVsPtA[track.trackacceptedid()]->Fill(track.pt(), phiInTpcSector);
        fhPerColNchVsPhiVsPtA[track.trackacceptedid()]->Fill(track.pt(), phiInTpcSector);
        fhPhiVsInnerWallMomA[track.trackacceptedid()]->Fill(track.tpcInnerParam(), phiInTpcSector);
        fhPerColNchVsPhiVsInnerWallMomA[track.trackacceptedid()]->Fill(track.tpcInnerParam(), phiInTpcSector);
        fhPhiShiftedVsPtA[track.trackacceptedid()]->Fill(track.pt(), phiShiftedPercentInTpcSector);
        fhPhiShiftedVsInnerWallMomA[track.trackacceptedid()]->Fill(track.tpcInnerParam(), phiShiftedPercentInTpcSector);
        fhItsNClsVsPtA[track.trackacceptedid()]->Fill(track.pt(), track.itsNCls());
        fhItsChi2NClsVsPtA[track.trackacceptedid()]->Fill(track.pt(), track.itsChi2NCl());
        fhTpcFindableNClsVsPtA[track.trackacceptedid()]->Fill(track.pt(), track.tpcNClsFindable());
        fhTpcFoundNClsVsPtA[track.trackacceptedid()]->Fill(track.pt(), track.tpcNClsFound());
        fhTpcSharedNClsVsPtA[track.trackacceptedid()]->Fill(track.pt(), track.tpcNClsShared());
        fhTpcFractionSharedClsVsPtA[track.trackacceptedid()]->Fill(track.pt(), track.tpcFractionSharedCls());
        fhTpcCrossedRowsVsPtA[track.trackacceptedid()]->Fill(track.pt(), track.tpcNClsCrossedRows());
        fhTpcCrossedRowsOverFindableClsVsPtA[track.trackacceptedid()]->Fill(track.pt(), track.tpcCrossedRowsOverFindableCls());
        fhTpcChi2NClsVsPtA[track.trackacceptedid()]->Fill(track.pt(), track.tpcChi2NCl());
        /* efficiency histograms */
        fhPvsInnerP->Fill(track.tpcInnerParam(), track.p());
        fillhisto(fhPtVsEtaItsA[track.trackacceptedid()], hasits);
        fillhisto(fhPtVsEtaTpcA[track.trackacceptedid()], hastpc);
        fillhisto(fhPtVsEtaItsTpcA[track.trackacceptedid()], hasits && hastpc);
        fillhisto(fhPtVsEtaItsTofA[track.trackacceptedid()], hasits && hastof);
        fillhisto(fhPtVsEtaTpcTofA[track.trackacceptedid()], hastpc && hastof);
        fillhisto(fhPtVsEtaItsTpcTofA[track.trackacceptedid()], hasits && hastpc && hastof);
        /* the detector / generator combined level */
        if constexpr (framework::has_type_v<aod::mctracklabel::McParticleId, typename TrackObject::all_columns>) {
          auto findgenid = [&](auto& part) {
            int pdgcode = std::abs(part.pdgCode());
            for (uint ix = 0; ix < pdgcodes.size(); ++ix) {
              if (pdgcode == pdgcodes[ix]) {
                return ix;
              }
            }
            return static_cast<uint>(pdgcodes.size());
          };
          auto fillpurityhistos = [](auto& hpos, auto& hneg, auto& genid, auto& track, bool cond) {
            if (cond) {
              if (track.sign() > 0) {
                hpos->Fill(static_cast<int>(track.trackacceptedid() / 2), genid, track.pt());
              } else {
                hneg->Fill(static_cast<int>(track.trackacceptedid() / 2), genid, track.pt());
              }
            }
          };
          /* get the associated MC particle we are sure it does exist because the track was accepted */
          const auto& mcparticle = track.template mcParticle_as<soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo>>();
          float genid = findgenid(mcparticle);

          bool isprimary = mcparticle.isPhysicalPrimary();
          bool issecdecay = !isprimary && (mcparticle.getProcess() == 4);
          bool isfrommaterial = !isprimary && !issecdecay;
          fillpurityhistos(fhPtPurityPosPrimA, fhPtPurityNegPrimA, genid, track, isprimary);
          fillpurityhistos(fhPtPurityPosSecA, fhPtPurityNegSecA, genid, track, issecdecay);
          fillpurityhistos(fhPtPurityPosMatA, fhPtPurityNegMatA, genid, track, isfrommaterial);
          fhTruePvsP->Fill(track.p(), mcparticle.p());
          fhTruePvsInnerP->Fill(track.tpcInnerParam(), mcparticle.p());

          auto fillhisto = [](auto& h, float pt, float eta, bool cond1, bool cond2) {
            if (cond1 && cond2) {
              h->Fill(eta, pt);
            }
          };
          std::vector<float> tPt = {track.pt(), mcparticle.pt()};
          std::vector<float> tEta = {track.eta(), mcparticle.eta()};
          for (uint ix = 0; ix < tPt.size(); ++ix) {
            fillhisto(fhPtVsEtaPrimItsA[ix][track.trackacceptedid()], tPt[ix], tEta[ix], hasits, isprimary);
            fillhisto(fhPtVsEtaPrimItsTpcA[ix][track.trackacceptedid()], tPt[ix], tEta[ix], hasits && hastpc, isprimary);
            fillhisto(fhPtVsEtaPrimItsTpcTofA[ix][track.trackacceptedid()], tPt[ix], tEta[ix], hasits && hastof, isprimary);
            fillhisto(fhPtVsEtaSecItsA[ix][track.trackacceptedid()], tPt[ix], tEta[ix], hasits, issecdecay);
            fillhisto(fhPtVsEtaSecItsTpcA[ix][track.trackacceptedid()], tPt[ix], tEta[ix], hasits && hastpc, issecdecay);
            fillhisto(fhPtVsEtaSecItsTpcTofA[ix][track.trackacceptedid()], tPt[ix], tEta[ix], hasits && hastof, issecdecay);
            fillhisto(fhPtVsEtaMatItsA[ix][track.trackacceptedid()], tPt[ix], tEta[ix], hasits, isfrommaterial);
            fillhisto(fhPtVsEtaMatItsTpcA[ix][track.trackacceptedid()], tPt[ix], tEta[ix], hasits && hastpc, isfrommaterial);
            fillhisto(fhPtVsEtaMatItsTpcTofA[ix][track.trackacceptedid()], tPt[ix], tEta[ix], hasits && hastof, isfrommaterial);
          }
        }
      }
    }
    if constexpr (kindOfData == kGen) {
      if (!(track.trackacceptedid() < 0)) {
        /* pure generator level */
        if (track.isPhysicalPrimary()) {
          fhPtVsEtaPrimA[track.trackacceptedid()]->Fill(track.eta(), track.pt());
        } else if (track.getProcess() == 4) {
          fhPtVsEtaSecA[track.trackacceptedid()]->Fill(track.eta(), track.pt());
        } else {
          fhPtVsEtaMatA[track.trackacceptedid()]->Fill(track.eta(), track.pt());
        }
      }
    }
  }

  template <efficiencyandqatask::KindOfData kindOfData>
  void newCollision()
  {
    using namespace efficiencyandqatask;
    if constexpr (kindOfData == kReco) {
      fhPerColNchVsPhiVsPtPosB->Reset();
      fhPerColNchVsPhiVsPtNegB->Reset();
      fhPerColNchVsPhiVsInnerWallMomPosB->Reset();
      fhPerColNchVsPhiVsInnerWallMomNegB->Reset();
      for (uint isp = 0; isp < nsp; ++isp) {
        fhPerColNchVsPhiVsPtA[isp]->Reset();
        fhPerColNchVsPhiVsInnerWallMomA[isp]->Reset();
      }
    }
  }

  template <efficiencyandqatask::KindOfData kindOfData>
  void finishedCollision()
  {
    using namespace efficiencyandqatask;
    if constexpr (kindOfData == kReco) {
      auto fillHistogram = [](auto& th, const TH2* sh) {
        int nBinsX = sh->GetNbinsX();
        int nBinsY = sh->GetNbinsY();
        for (int ix = 0; ix < nBinsX; ++ix) {
          for (int iy = 0; iy < nBinsY; ++iy) {
            th->Fill(sh->GetXaxis()->GetBinCenter(ix + 1), sh->GetYaxis()->GetBinCenter(iy + 1), sh->GetBinContent(ix + 1, iy + 1));
          }
        }
      };
      fillHistogram(fhNchVsPhiVsPtPosB, fhPerColNchVsPhiVsPtPosB);
      fillHistogram(fhNchVsPhiVsPtNegB, fhPerColNchVsPhiVsPtNegB);
      fillHistogram(fhNchVsPhiVsInnerWallMomPosB, fhPerColNchVsPhiVsInnerWallMomPosB);
      fillHistogram(fhNchVsPhiVsInnerWallMomNegB, fhPerColNchVsPhiVsInnerWallMomNegB);
      for (uint isp = 0; isp < nsp; ++isp) {
        fillHistogram(fhNchVsPhiVsPtA[isp], fhPerColNchVsPhiVsPtA[isp]);
        fillHistogram(fhNchVsPhiVsInnerWallMomA[isp], fhPerColNchVsPhiVsInnerWallMomA[isp]);
      }
    }
  }
};

/* the QA extra data, pairs, collecting engine */
struct QAExtraDataCollectingEngine {
  uint nsp = static_cast<uint>(efficiencyandqatask::tnames.size());
  uint nmainsp = static_cast<uint>(efficiencyandqatask::mainspnames.size());
  uint nallmainsp = static_cast<uint>(efficiencyandqatask::allmainspnames.size());

  //===================================================
  // The QA output objects
  //===================================================
  /* pairs histograms */
  std::vector<std::vector<std::vector<std::shared_ptr<TH2>>>> fhPhiPhiA{2, {nsp, {nsp, nullptr}}};
  std::vector<std::vector<std::vector<std::shared_ptr<TH3>>>> fhDeltaPhiVsPhiPhiA{2, {nsp, {nsp, nullptr}}};
  std::vector<std::vector<std::vector<std::shared_ptr<TH3>>>> fhDeltaPhiVsEtaEtaA{2, {nsp, {nsp, nullptr}}};
  std::vector<std::vector<std::vector<std::shared_ptr<TH2>>>> fhEtaEtaA{2, {nsp, {nsp, nullptr}}};
  std::vector<std::vector<std::vector<std::shared_ptr<TH3>>>> fhDeltaEtaVsEtaEtaA{2, {nsp, {nsp, nullptr}}};
  std::vector<std::vector<std::vector<std::shared_ptr<TH3>>>> fhDeltaEtaVsPhiPhiA{2, {nsp, {nsp, nullptr}}};

  template <efficiencyandqatask::KindOfData kindOfData>
  void init(HistogramRegistry& registry, const char* dirname)
  {
    using namespace efficiencyandqatask;
    using namespace analysis::dptdptfilter;

    AxisSpec phiAxis = {phibins, 0.0f, constants::math::TwoPI, "#varphi"};
    AxisSpec deltaPhiAxis = {phibins, 0.0f, constants::math::TwoPI, "#Delta#varphi"};
    AxisSpec etaAxis = {etabins, etalow, etaup, "#eta"};
    AxisSpec deltaEtaAxis = {2 * etabins - 1, etalow - etaup, etaup - etalow, "#DeltaEta"};

    /* the reconstructed and generated levels histograms */
    std::string recogen = (kindOfData == kReco) ? "Reco" : "Gen";
    for (uint isp = 0; isp < nsp; ++isp) {
      for (uint jsp = 0; jsp < nsp; ++jsp) {
        fhPhiPhiA[kindOfData][isp][jsp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, recogen.c_str(), "After"), HNAMESTRING("PhiPhi_%s%s", tnames[isp].c_str(), tnames[jsp].c_str()),
                                                       HTITLESTRING("%s%s pairs", tnames[isp].c_str(), tnames[jsp].c_str()), kTH2F, {phiAxis, phiAxis});
        fhDeltaPhiVsPhiPhiA[kindOfData][isp][jsp] = ADDHISTOGRAM(TH3, DIRECTORYSTRING("%s/%s/%s", dirname, recogen.c_str(), "After"), HNAMESTRING("DeltaPhiVsPhiPhi_%s%s", tnames[isp].c_str(), tnames[jsp].c_str()),
                                                                 HTITLESTRING("%s%s pairs", tnames[isp].c_str(), tnames[jsp].c_str()), kTH3F, {phiAxis, phiAxis, deltaPhiAxis});
        fhDeltaEtaVsPhiPhiA[kindOfData][isp][jsp] = ADDHISTOGRAM(TH3, DIRECTORYSTRING("%s/%s/%s", dirname, recogen.c_str(), "After"), HNAMESTRING("DeltaEtaVsPhiPhi_%s%s", tnames[isp].c_str(), tnames[jsp].c_str()),
                                                                 HTITLESTRING("%s%s pairs", tnames[isp].c_str(), tnames[jsp].c_str()), kTH3F, {phiAxis, phiAxis, deltaEtaAxis});
        fhEtaEtaA[kindOfData][isp][jsp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, recogen.c_str(), "After"), HNAMESTRING("EtaEta_%s%s", tnames[isp].c_str(), tnames[jsp].c_str()),
                                                       HTITLESTRING("%s%s pairs", tnames[isp].c_str(), tnames[jsp].c_str()), kTH2F, {etaAxis, etaAxis});
        fhDeltaEtaVsEtaEtaA[kindOfData][isp][jsp] = ADDHISTOGRAM(TH3, DIRECTORYSTRING("%s/%s/%s", dirname, recogen.c_str(), "After"), HNAMESTRING("DeltaEtaVsEtaEta_%s%s", tnames[isp].c_str(), tnames[jsp].c_str()),
                                                                 HTITLESTRING("%s%s pairs", tnames[isp].c_str(), tnames[jsp].c_str()), kTH3F, {etaAxis, etaAxis, deltaEtaAxis});
        fhDeltaPhiVsEtaEtaA[kindOfData][isp][jsp] = ADDHISTOGRAM(TH3, DIRECTORYSTRING("%s/%s/%s", dirname, recogen.c_str(), "After"), HNAMESTRING("DeltaPhiVsEtaEta_%s%s", tnames[isp].c_str(), tnames[jsp].c_str()),
                                                                 HTITLESTRING("%s%s pairs", tnames[isp].c_str(), tnames[jsp].c_str()), kTH3F, {etaAxis, etaAxis, deltaPhiAxis});
      }
    }
  }

  template <efficiencyandqatask::KindOfData kindOfData, typename CollisionsObject, typename TracksObject>
  void processTrackPairs(TracksObject const& tracks1, TracksObject const& tracks2)
  {
    using namespace efficiencyandqatask;
    using namespace analysis::dptdptfilter;

    /* we should only receive accepted tracks */
    for (auto const& track1 : tracks1) {
      for (auto const& track2 : tracks2) {
        /* checking the same track id condition */
        if (track1 == track2) {
          /* exclude autocorrelations */
          continue;
        }
        float deltaPhi = RecoDecay::constrainAngle(track1.phi() - track2.phi());
        float deltaEta = track1.eta() - track2.eta();
        fhPhiPhiA[kindOfData][track1.trackacceptedid()][track2.trackacceptedid()]->Fill(track1.phi(), track2.phi());
        fhDeltaPhiVsPhiPhiA[kindOfData][track1.trackacceptedid()][track2.trackacceptedid()]->Fill(track1.phi(), track2.phi(), deltaPhi);
        fhDeltaEtaVsPhiPhiA[kindOfData][track1.trackacceptedid()][track2.trackacceptedid()]->Fill(track1.phi(), track2.phi(), deltaEta);
        fhEtaEtaA[kindOfData][track1.trackacceptedid()][track2.trackacceptedid()]->Fill(track1.eta(), track2.eta());
        fhDeltaEtaVsEtaEtaA[kindOfData][track1.trackacceptedid()][track2.trackacceptedid()]->Fill(track1.eta(), track2.eta(), deltaEta);
        fhDeltaPhiVsEtaEtaA[kindOfData][track1.trackacceptedid()][track2.trackacceptedid()]->Fill(track1.eta(), track2.eta(), deltaPhi);
      }
    }
  }
};

/* the PID data collecting engine */
struct PidDataCollectingEngine {
  uint nsp = static_cast<uint>(efficiencyandqatask::tnames.size());
  uint nmainsp = static_cast<uint>(efficiencyandqatask::mainspnames.size());
  uint nallmainsp = static_cast<uint>(efficiencyandqatask::allmainspnames.size());

  /* PID histograms */
  /* before and after */
  std::vector<std::shared_ptr<TH2>> fhTPCdEdxSignalVsP{2, nullptr};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhTPCdEdxSignalDiffVsP{2, {nmainsp, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhTPCnSigmasVsP{2, {nallmainsp, nullptr}};
  std::vector<std::shared_ptr<TH2>> fhTOFSignalVsP{2, nullptr};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhTOFSignalDiffVsP{2, {nmainsp, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhTOFnSigmasVsP{2, {nallmainsp, nullptr}};
  std::vector<std::shared_ptr<TH2>> fhPvsTOFSqMass{2, nullptr};
  std::vector<std::vector<std::shared_ptr<TH3>>> fhTPCTOFSigmaVsP{2, {nmainsp, nullptr}};

  template <efficiencyandqatask::KindOfData kindOfData>
  void init(HistogramRegistry& registry, const char* dirname)
  {
    using namespace efficiencyandqatask;

    const AxisSpec dEdxAxis{200, 0.0, 200.0, "dE/dx (au)"};
    AxisSpec pidPAxis{150, 0.1, 5.0, "#it{p} (GeV/#it{c})"};
    pidPAxis.makeLogarithmic();

    if constexpr (kindOfData == kReco) {
      /* PID histograms */
      std::vector<std::string> whenname{"Before", "After"};
      char whenprefix[2]{'B', 'A'};
      std::vector<std::string> whentitle{"before", ""};
      for (uint ix = 0; ix < whenname.size(); ++ix) {
        fhTPCdEdxSignalVsP[ix] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "PID", whenname[ix].c_str()),
                                              HNAMESTRING("tpcSignalVsP%c", whenprefix[ix]),
                                              HTITLESTRING("TPC dE/dx signal %s", whentitle[ix].c_str()), kTH2F, {pidPAxis, dEdxAxis});
        fhTOFSignalVsP[ix] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "PID", whenname[ix].c_str()),
                                          HNAMESTRING("tofSignalVsP%c", whenprefix[ix]),
                                          HTITLESTRING("TOF signal %s", whentitle[ix].c_str()),
                                          kTH2F, {pidPAxis, {200, 0.0, 1.1, "#beta"}});
        fhPvsTOFSqMass[ix] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "PID", whenname[ix].c_str()),
                                          HNAMESTRING("tofPvsMassSq%c", whenprefix[ix]),
                                          HTITLESTRING("Momentum versus #it{m}^{2} %s", whentitle[ix].c_str()),
                                          kTH2F, {{140, 0.0, 1.4, "#it{m}^{2} ((GeV/c^{2})^{2})"}, pidPAxis});
        for (uint isp = 0; isp < nmainsp; ++isp) {
          fhTPCdEdxSignalDiffVsP[ix][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "PID", whenname[ix].c_str()),
                                                         HNAMESTRING("tpcSignalDiffVsP%c_%s", whenprefix[ix], mainspnames[isp].c_str()),
                                                         HTITLESTRING("TPC dE/dx to the %s line %s", mainsptitles[isp].c_str(), whentitle[ix].c_str()),
                                                         kTH2F, {pidPAxis, {400, -200.0, 200.0, FORMATSTRING("dE/dx - <dE/dx>_{%s}", mainsptitles[isp].c_str())}});
          fhTOFSignalDiffVsP[ix][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "PID", whenname[ix].c_str()),
                                                     HNAMESTRING("tofSignalDiffVsP%c_%s", whenprefix[ix], mainspnames[isp].c_str()),
                                                     HTITLESTRING("#Delta^{TOF_{%s}} %s", mainsptitles[isp].c_str(), whentitle[ix].c_str()),
                                                     kTH2F, {pidPAxis, {200, -1000.0, 1000.0, FORMATSTRING("t-t_{ev}-t_{exp_{%s}} (ps)", mainsptitles[isp].c_str())}});
          fhTPCTOFSigmaVsP[ix][isp] = ADDHISTOGRAM(TH3, DIRECTORYSTRING("%s/%s/%s", dirname, "PID", whenname[ix].c_str()),
                                                   HNAMESTRING("toftpcNSigmasVsP%c_%s", whenprefix[ix], mainspnames[isp].c_str()),
                                                   HTITLESTRING("n#sigma to the %s line %s", mainsptitles[isp].c_str(), whentitle[ix].c_str()),
                                                   kTH3F, {pidPAxis, {noOfNSigmaBins, minNSigma, maxNSigma, FORMATSTRING("n#sigma_{TPC}^{%s}", mainsptitles[isp].c_str())}, {120, -6.0, 6.0, FORMATSTRING("n#sigma_{TOF}^{%s}", mainsptitles[isp].c_str())}});
        }
        for (uint isp = 0; isp < nallmainsp; ++isp) {
          fhTPCnSigmasVsP[ix][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "PID", whenname[ix].c_str()),
                                                  HNAMESTRING("tpcNSigmasVsP%c_%s", whenprefix[ix], allmainspnames[isp].c_str()),
                                                  HTITLESTRING("TPC n#sigma to the %s line %s", allmainsptitles[isp].c_str(), whentitle[ix].c_str()),
                                                  kTH2F, {pidPAxis, {noOfNSigmaBins, minNSigma, maxNSigma, FORMATSTRING("n#sigma_{TPC}^{%s}", allmainsptitles[isp].c_str())}});
          fhTOFnSigmasVsP[ix][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "PID", whenname[ix].c_str()),
                                                  HNAMESTRING("tofNSigmasVsP%c_%s", whenprefix[ix], allmainspnames[isp].c_str()),
                                                  HTITLESTRING("TOF n#sigma to the %s line %s", allmainsptitles[isp].c_str(), whentitle[ix].c_str()),
                                                  kTH2F, {pidPAxis, {noOfNSigmaBins, minNSigma, maxNSigma, FORMATSTRING("n#sigma_{TOF}^{%s}", allmainsptitles[isp].c_str())}});
        }
      }
    }
  }

  template <o2::track::PID::ID id, typename TrackObject>
  void fillAllSpeciesPID(uint ix, TrackObject const& track, float tpcmom, float tofmom)
  {
    if (track.sign() < 0) {
      ix = 2 * ix + 1;
    } else {
      ix = 2 * ix;
    }
    for (uint when = 0; when < 2; ++when) {
      fhTPCnSigmasVsP[when][ix]->Fill(tpcmom, o2::aod::pidutils::tpcNSigma<id>(track));
      fhTOFnSigmasVsP[when][ix]->Fill(tofmom, o2::aod::pidutils::tofNSigma<id>(track));
      if (track.trackacceptedid() < 0) {
        /* track not accepted */
        return;
      }
    }
  }

  template <o2::track::PID::ID id, typename TrackObject>
  void fillSpeciesPID(uint ix, TrackObject const& track, float tpcmom, float tofmom)
  {
    if (track.sign() < 0) {
      ix = 2 * ix + 1;
    } else {
      ix = 2 * ix;
    }
    for (uint when = 0; when < 2; ++when) {
      fhTPCdEdxSignalDiffVsP[when][ix]->Fill(tpcmom, o2::aod::pidutils::tpcExpSignalDiff<id>(track));
      fhTOFSignalDiffVsP[when][ix]->Fill(tofmom, o2::aod::pidutils::tofExpSignalDiff<id>(track));
      fhTPCTOFSigmaVsP[when][ix]->Fill(tpcmom, o2::aod::pidutils::tpcNSigma<id>(track), o2::aod::pidutils::tofNSigma<id>(track));
      if (track.trackacceptedid() < 0) {
        /* track not accepted */
        return;
      }
    }
  }

  template <typename TrackObject>
  void fillPID(TrackObject const& track, float tpcmom, float tofmom)
  {
    for (uint when = 0; when < 2; ++when) {
      if constexpr (framework::has_type_v<o2::aod::mcpidtpc::DeDxTunedMc, typename TrackObject::all_columns>) {
        fhTPCdEdxSignalVsP[when]->Fill(tpcmom, track.mcTunedTPCSignal());
      } else {
        fhTPCdEdxSignalVsP[when]->Fill(tpcmom, track.tpcSignal());
      }
      fhTOFSignalVsP[when]->Fill(tofmom, track.beta());
      fhPvsTOFSqMass[when]->Fill(track.mass() * track.mass(), tofmom);
      if (track.trackacceptedid() < 0) {
        /* track not accepted */
        return;
      }
    }
  }

  template <efficiencyandqatask::KindOfData kindOfData, typename TrackObject>
  void processTrack(TrackObject const& track, float tpcmom, float tofmom)
  {
    using namespace efficiencyandqatask;

    if constexpr (kindOfData == kReco) {
      fillPID(track, tpcmom, tofmom);
      fillSpeciesPID<o2::track::PID::Pion>(0, track, tpcmom, tofmom);
      fillSpeciesPID<o2::track::PID::Kaon>(1, track, tpcmom, tofmom);
      fillSpeciesPID<o2::track::PID::Proton>(2, track, tpcmom, tofmom);
      fillAllSpeciesPID<o2::track::PID::Electron>(0, track, tpcmom, tofmom);
      fillAllSpeciesPID<o2::track::PID::Muon>(1, track, tpcmom, tofmom);
      fillAllSpeciesPID<o2::track::PID::Pion>(2, track, tpcmom, tofmom);
      fillAllSpeciesPID<o2::track::PID::Kaon>(3, track, tpcmom, tofmom);
      fillAllSpeciesPID<o2::track::PID::Proton>(4, track, tpcmom, tofmom);
    }
  }
};

/* the PID extra data collecting engine */
struct PidExtraDataCollectingEngine {
  uint nsp = static_cast<uint>(efficiencyandqatask::tnames.size());
  uint nmainsp = static_cast<uint>(efficiencyandqatask::mainspnames.size());
  uint nallmainsp = static_cast<uint>(efficiencyandqatask::allmainspnames.size());

  /* PID histograms */
  /* only after track selection */
  std::vector<std::shared_ptr<TH2>> fhIdTPCdEdxSignalVsP{nsp, nullptr};
  std::vector<std::shared_ptr<TProfile2D>> fpIdTPCdEdxSignalVsPSigmas{nsp, nullptr};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhIdTPCdEdxSignalDiffVsP{nsp, {nmainsp, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhIdTPCnSigmasVsP{nsp, {nallmainsp, nullptr}};
  std::vector<std::shared_ptr<TH2>> fhIdTOFSignalVsP{nsp, nullptr};
  std::vector<std::shared_ptr<TProfile2D>> fpIdTOFSignalVsPSigmas{nsp, nullptr};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhIdTOFSignalDiffVsP{nsp, {nmainsp, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhIdTOFnSigmasVsP{nsp, {nallmainsp, nullptr}};
  std::vector<std::shared_ptr<TH2>> fhIdPvsTOFSqMass{nsp, nullptr};

  template <efficiencyandqatask::KindOfData kindOfData>
  void init(HistogramRegistry& registry, const char* dirname)
  {
    using namespace efficiencyandqatask;

    const AxisSpec dEdxAxis{200, 0.0, 200.0, "dE/dx (au)"};
    AxisSpec pidPAxis{150, 0.1, 5.0, "#it{p} (GeV/#it{c})"};
    pidPAxis.makeLogarithmic();

    if constexpr (kindOfData == kReco) {
      /* PID histograms */
      for (uint isp = 0; isp < nsp; ++isp) {
        fhIdTPCdEdxSignalVsP[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "PID", "Selected"),
                                                 HNAMESTRING("tpcSignalVsPSelected_%s", tnames[isp].c_str()),
                                                 HTITLESTRING("TPC dE/dx for selected %s", tnames[isp].c_str()),
                                                 kTH2F, {pidPAxis, dEdxAxis});
        fpIdTPCdEdxSignalVsPSigmas[isp] = ADDHISTOGRAM(TProfile2D, DIRECTORYSTRING("%s/%s/%s", dirname, "PID", "Selected"),
                                                       HNAMESTRING("tpcSignalSigmasVsPSelected_%s", tnames[isp].c_str()),
                                                       HTITLESTRING("TPC dE/dx and n#sigma for selected %s", tnames[isp].c_str()),
                                                       kTProfile2D, {pidPAxis, dEdxAxis});
        fhIdTOFSignalVsP[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "PID", "Selected"),
                                             HNAMESTRING("tofSignalVsPSelected_%s", tnames[isp].c_str()),
                                             HTITLESTRING("TOF signal for selected %s", tnames[isp].c_str()),
                                             kTH2F, {pidPAxis, {200, 0.0, 1.1, "#beta"}});
        fpIdTOFSignalVsPSigmas[isp] = ADDHISTOGRAM(TProfile2D, DIRECTORYSTRING("%s/%s/%s", dirname, "PID", "Selected"),
                                                   HNAMESTRING("tofSignalSigmasVsPSelected_%s", tnames[isp].c_str()),
                                                   HTITLESTRING("TOF signal and n#sigma for selected %s", tnames[isp].c_str()),
                                                   kTProfile2D, {pidPAxis, {200, 0.0, 1.1, "#beta"}});
        for (uint imainsp = 0; imainsp < nallmainsp; ++imainsp) {
          /* only the same charge makes any sense */
          if (isp % 2 == imainsp % 2) {
            fhIdTPCnSigmasVsP[isp][imainsp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "PID", "Selected"),
                                                           HNAMESTRING("tpcNSigmasVsPSelected_%s_to%s", tnames[isp].c_str(), allmainspnames[imainsp].c_str()),
                                                           HTITLESTRING("TPC n#sigma for selected %s to the %s line", tnames[isp].c_str(), allmainsptitles[imainsp].c_str()),
                                                           kTH2F, {pidPAxis, {noOfNSigmaBins, minNSigma, maxNSigma, FORMATSTRING("n#sigma_{TPC}^{%s}", mainsptitles[isp].c_str())}});
            fhIdTOFnSigmasVsP[isp][imainsp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "PID", "Selected"),
                                                           HNAMESTRING("tofNSigmasVsPSelected_%s_to%s", tnames[isp].c_str(), allmainspnames[imainsp].c_str()),
                                                           HTITLESTRING("TOF n#sigma for selected %s to the %s line", tnames[isp].c_str(), allmainsptitles[imainsp].c_str()),
                                                           kTH2F, {pidPAxis, {noOfNSigmaBins, minNSigma, maxNSigma, FORMATSTRING("n#sigma_{TOF}^{%s}", mainsptitles[isp].c_str())}});
          }
        }
      }
    }
  }

  template <o2::track::PID::ID id, typename TrackObject>
  void fillAllSpeciesPID(uint ix, TrackObject const& track, float tpcmom, float tofmom)
  {
    if (track.trackacceptedid() < 0) {
      /* track not accepted */
      return;
    }
    if (track.sign() < 0) {
      ix = 2 * ix + 1;
    } else {
      ix = 2 * ix;
    }
    fhIdTPCnSigmasVsP[track.trackacceptedid()][ix]->Fill(tpcmom, o2::aod::pidutils::tpcNSigma<id>(track));
    fhIdTOFnSigmasVsP[track.trackacceptedid()][ix]->Fill(tofmom, o2::aod::pidutils::tofNSigma<id>(track));
    if (efficiencyandqatask::pidselector.isGlobalSpecies(track.trackacceptedid() / 2, id)) {
      /* only if the species of the selected track matches the target of the number of sigmas */
      if constexpr (framework::has_type_v<o2::aod::mcpidtpc::DeDxTunedMc, typename TrackObject::all_columns>) {
        fpIdTPCdEdxSignalVsPSigmas[track.trackacceptedid()]->Fill(tpcmom, track.mcTunedTPCSignal(), o2::aod::pidutils::tpcNSigma<id>(track));
      } else {
        fpIdTPCdEdxSignalVsPSigmas[track.trackacceptedid()]->Fill(tpcmom, track.tpcSignal(), o2::aod::pidutils::tpcNSigma<id>(track));
      }
      fpIdTOFSignalVsPSigmas[track.trackacceptedid()]->Fill(tofmom, track.beta(), o2::aod::pidutils::tofNSigma<id>(track));
    }
  }

  template <o2::track::PID::ID id, typename TrackObject>
  void fillSpeciesPID(uint, TrackObject const&, float, float)
  {
  }

  template <typename TrackObject>
  void fillPID(TrackObject const& track, float tpcmom, float tofmom)
  {
    if (track.trackacceptedid() < 0) {
      /* track not accepted */
      return;
    }
    if constexpr (framework::has_type_v<o2::aod::mcpidtpc::DeDxTunedMc, typename TrackObject::all_columns>) {
      fhIdTPCdEdxSignalVsP[track.trackacceptedid()]->Fill(tpcmom, track.mcTunedTPCSignal());
    } else {
      fhIdTPCdEdxSignalVsP[track.trackacceptedid()]->Fill(tpcmom, track.tpcSignal());
    }
    fhIdTOFSignalVsP[track.trackacceptedid()]->Fill(tofmom, track.beta());
  }

  template <efficiencyandqatask::KindOfData kindOfData, typename TrackObject>
  void processTrack(TrackObject const& track, float tpcmom, float tofmom)
  {
    using namespace efficiencyandqatask;

    if constexpr (kindOfData == kReco) {
      fillPID(track, tpcmom, tofmom);
      fillSpeciesPID<o2::track::PID::Pion>(0, track, tpcmom, tofmom);
      fillSpeciesPID<o2::track::PID::Kaon>(1, track, tpcmom, tofmom);
      fillSpeciesPID<o2::track::PID::Proton>(2, track, tpcmom, tofmom);
      fillAllSpeciesPID<o2::track::PID::Electron>(0, track, tpcmom, tofmom);
      fillAllSpeciesPID<o2::track::PID::Muon>(1, track, tpcmom, tofmom);
      fillAllSpeciesPID<o2::track::PID::Pion>(2, track, tpcmom, tofmom);
      fillAllSpeciesPID<o2::track::PID::Kaon>(3, track, tpcmom, tofmom);
      fillAllSpeciesPID<o2::track::PID::Proton>(4, track, tpcmom, tofmom);
    }
  }
};

struct DptDptEfficiencyAndQc {
  /* the data memebers for this task */
  /* the centrality / multiplicity limits for collecting data in this task instance */
  uint ncmranges = 0;
  float* fCentMultMin = nullptr;
  float* fCentMultMax = nullptr;

  /* the data collecting engine instances */
  QADataCollectingEngine** qaDataCE;
  QAExtraDataCollectingEngine** qaExtraDataCE;
  PidDataCollectingEngine** pidDataCE;
  PidExtraDataCollectingEngine** pidExtraDataCE;

  /* the histogram registries */
  HistogramRegistry registryOne{"registryOne", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryTwo{"registryTwo", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryThree{"registryThree", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryFour{"registryFour", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryFive{"registryFive", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registrySix{"registrySix", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registrySeven{"registrySeven", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryEight{"registryEight", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryNine{"registryNine", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryTen{"registryTen", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryExtraOne{"extraregistryOne", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryExtraTwo{"extraregistryTwo", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryExtraThree{"extraregistryThree", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryExtraFour{"extraregistryFour", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryExtraFive{"extraregistryFive", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryExtraSix{"extraregistrySix", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryExtraSeven{"extraregistrySeven", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryExtraEight{"extraregistryEight", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryExtraNine{"extraregistryNine", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryExtraTen{"extraregistryTen", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryPidOne{"pidregistryOne", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryPidTwo{"pidregistryTwo", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryPidThree{"pidregistryThree", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryPidFour{"pidregistryFour", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryPidFive{"pidregistryFive", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryPidSix{"pidregistrySix", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryPidSeven{"pidregistrySeven", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryPidEight{"pidregistryEight", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryPidNine{"pidregistryNine", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registryPidTen{"pidregistryTen", {}, OutputObjHandlingPolicy::AnalysisObject};
  std::vector<HistogramRegistry*> registryBank{&registryOne, &registryTwo, &registryThree, &registryFour, &registryFive,
                                               &registrySix, &registrySeven, &registryEight, &registryNine, &registryTen};
  std::vector<HistogramRegistry*> extraRegistryBank{&registryExtraOne, &registryExtraTwo, &registryExtraThree, &registryExtraFour, &registryExtraFive,
                                                    &registryExtraSix, &registryExtraSeven, &registryExtraEight, &registryExtraNine, &registryExtraTen};
  std::vector<HistogramRegistry*> pidRegistryBank{&registryPidOne, &registryPidTwo, &registryPidThree, &registryPidFour, &registryPidFive,
                                                  &registryPidSix, &registryPidSeven, &registryPidEight, &registryPidNine, &registryPidTen};

  Configurable<bool> useCentrality{"useCentrality", false, "Perform the task using centrality/multiplicity classes. Default value: false"};
  Configurable<bool> useTPCInnerWallMomentum{"useTPCInnerWallMomentum", false, "Use the TPC inner wall momentum. Default: false"};
  Configurable<float> cfgMinNSigma{"cfgMinNSigma", -4.05f, "nsigma axes lowest value. Default: -4.05"};
  Configurable<float> cfgMaxNSigma{"cfgMaxNSigma", 4.05f, "nsigma axes highest value. Default: 4.05"};
  Configurable<float> cfgWidthNSigmaBin{"cfgWidthNSigmaBin", 0.1, "nsigma axes bin width. Deafault: 0.1"};

  void init(o2::framework::InitContext& initContext)
  {
    using namespace efficiencyandqatask;
    using namespace analysis::dptdptfilter;

    /* do nothing if not active */
    if (!doprocessDetectorLevelNotStored &&
        !doprocessExtraDetectorLevelNotStored &&
        !doprocessDetectorLevelNotStoredPID &&
        !doprocessDetectorLevelNotStoredTunedOnDataPID &&
        !doprocessDetectorLevelNotStoredPIDExtra &&
        !doprocessDetectorLevelNotStoredTunedOnDataPIDExtra &&
        !doprocessGeneratorLevelNotStored &&
        !doprocessExtraGeneratorLevelNotStored &&
        !doprocessReconstructedNotStored &&
        !doprocessReconstructedNotStoredPID &&
        !doprocessReconstructedNotStoredPIDExtra) {
      return;
    }

    /* Self configuration: requires dptdptfilter task in the workflow */
    {
      /* the binning */
      getTaskOptionValue(initContext, "dpt-dpt-filter", "overallminp", overallminp, false);
      getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mZVtxbins", zvtxbins, false);
      getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mZVtxmin", zvtxlow, false);
      getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mZVtxmax", zvtxup, false);
      getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mPTbins", ptbins, false);
      getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mPTmin", ptlow, false);
      getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mPTmax", ptup, false);
      getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mEtabins", etabins, false);
      getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mEtamin", etalow, false);
      getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mEtamax", etaup, false);
      getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mPhibins", phibins, false);

      /* configuring the involved species */
      std::vector<std::string> cfgnames = {"elpidsel", "mupidsel", "pipidsel", "kapidsel", "prpidsel"};
      std::vector<uint8_t> spids = {0, 1, 2, 3, 4};
      for (uint i = 0; i < cfgnames.size(); ++i) {
        auto includeIt = [&initContext](int spid, auto name) {
          bool mUseIt = false;
          bool mExcludeIt = false;
          if (getTaskOptionValue(initContext, "dpt-dpt-filter-tracks", TString::Format("%s.mUseIt", name.c_str()).Data(), mUseIt, false) &&
              getTaskOptionValue(initContext, "dpt-dpt-filter-tracks", TString::Format("%s.mExclude", name.c_str()).Data(), mExcludeIt, false)) {
            if (mUseIt && !mExcludeIt) {
              auto cfg = new o2::analysis::TrackSelectionPIDCfg();
              cfg->mUseIt = true;
              cfg->mExclude = false;
              pidselector.addSpecies(spid, cfg);
            }
          }
        };
        includeIt(spids[i], cfgnames[i]);
      }
      uint nspecies = pidselector.getNSpecies();
      if (nspecies == 0) {
        /* unidentified analysis */
        poinames.push_back(pidselector.getHadFName());
        tnames.push_back(std::string(TString::Format("%sP", pidselector.getHadFName()).Data()));
        tnames.push_back(std::string(TString::Format("%sM", pidselector.getHadFName()).Data()));
        LOGF(info, "Incorporated species name %s to the analysis", poinames[0].c_str());
      } else {
        for (uint8_t ix = 0; ix < nspecies; ++ix) {
          poinames.push_back(std::string(pidselector.getSpeciesFName(ix)));
          tnames.push_back(std::string(TString::Format("%sP", pidselector.getSpeciesFName(ix)).Data()));
          tnames.push_back(std::string(TString::Format("%sM", pidselector.getSpeciesFName(ix)).Data()));
          LOGF(info, "Incorporated species name %s to the analysis", poinames[ix].c_str());
        }
      }

      /* create the data collecting engine instances according to the configured centrality/multiplicity ranges */
      std::string centspec;
      if (useCentrality.value && getTaskOptionValue(initContext, "dpt-dpt-filter", "centralities", centspec, false)) {
        LOGF(info, "Got the centralities specification: %s", centspec.c_str());
        auto tokens = TString(centspec.c_str()).Tokenize(",");
        ncmranges = tokens->GetEntries();
        fCentMultMin = new float[ncmranges];
        fCentMultMax = new float[ncmranges];
        for (uint i = 0; i < ncmranges; ++i) {
          float cmmin = 0.0f;
          float cmmax = 0.0f;
          sscanf(tokens->At(i)->GetName(), "%f-%f", &cmmin, &cmmax);
          fCentMultMin[i] = cmmin;
          fCentMultMax[i] = cmmax;
        }
        delete tokens;
      } else {
        LOGF(info, "No centralities specification. Setting it to: 0-100");
        ncmranges = 1;
        fCentMultMin = new float[ncmranges];
        fCentMultMax = new float[ncmranges];
        fCentMultMin[0] = 0.0f;
        fCentMultMax[0] = 100.0f;
      }
      /* configure nsigma axes */
      minNSigma = cfgMinNSigma.value;
      maxNSigma = cfgMaxNSigma.value;
      widthNSigmaBin = cfgWidthNSigmaBin.value;
      noOfNSigmaBins = static_cast<int>((maxNSigma - minNSigma) / widthNSigmaBin);

      bool doBasicAnalysis = doprocessDetectorLevelNotStored || doprocessReconstructedNotStored || doprocessGeneratorLevelNotStored;
      bool doExtraAnalysis = doprocessExtraDetectorLevelNotStored || doprocessExtraReconstructedNotStored || doprocessExtraGeneratorLevelNotStored;
      bool doPidAnalysis = doprocessDetectorLevelNotStoredPID || doprocessDetectorLevelNotStoredTunedOnDataPID || doprocessReconstructedNotStoredPID;
      bool doPidExtraAnalysis = doprocessDetectorLevelNotStoredPIDExtra || doprocessDetectorLevelNotStoredTunedOnDataPIDExtra || doprocessReconstructedNotStoredPIDExtra;

      if (doBasicAnalysis) {
        qaDataCE = new QADataCollectingEngine*[ncmranges];
      }
      if (doExtraAnalysis) {
        qaExtraDataCE = new QAExtraDataCollectingEngine*[ncmranges];
      }
      if (doPidAnalysis) {
        pidDataCE = new PidDataCollectingEngine*[ncmranges];
      }
      if (doPidExtraAnalysis) {
        pidExtraDataCE = new PidExtraDataCollectingEngine*[ncmranges];
      }
      std::string recogen;
      if (ncmranges > registryBank.size()) {
        LOGF(fatal, "There are more centrality ranges configured than registries in the bank. Please fix it!");
      }
      /* in reverse order for proper order in results file */
      for (uint i = 0; i < ncmranges; ++i) {
        auto initializeCEInstance = [&](auto dce, auto name, auto& registry, bool genlevel) {
          /* crete the output list for the passed centrality/multiplicity range */
          /* init the data collection instance */
          dce->template init<kReco>(registry, name.Data());
          if (genlevel) {
            dce->template init<kGen>(registry, name.Data());
          }
        };
        auto buildQACEInstance = [&](float min, float max) {
          auto* dce = new QADataCollectingEngine();
          initializeCEInstance(dce, TString::Format("EfficiencyAndQaData-%d-%d", static_cast<int>(min), static_cast<int>(max)), *registryBank[i], doprocessGeneratorLevelNotStored);
          return dce;
        };
        auto buildQACEExtraInstance = [&](float min, float max) {
          auto* dce = new QAExtraDataCollectingEngine();
          initializeCEInstance(dce, TString::Format("EfficiencyAndQaExtraData-%d-%d", static_cast<int>(min), static_cast<int>(max)), *extraRegistryBank[i], doprocessExtraGeneratorLevelNotStored);
          return dce;
        };
        auto buildPidCEInstance = [&](float min, float max) {
          auto* dce = new PidDataCollectingEngine();
          initializeCEInstance(dce, TString::Format("EfficiencyAndPidData-%d-%d", static_cast<int>(min), static_cast<int>(max)), *pidRegistryBank[i], doprocessGeneratorLevelNotStored);
          return dce;
        };
        auto buildPidExtraCEInstance = [&](float min, float max) {
          auto* dce = new PidExtraDataCollectingEngine();
          initializeCEInstance(dce, TString::Format("EfficiencyAndPidData-%d-%d", static_cast<int>(min), static_cast<int>(max)), *pidRegistryBank[i], doprocessGeneratorLevelNotStored);
          return dce;
        };
        /* in reverse order for proper order in results file */
        if (doBasicAnalysis) {
          qaDataCE[ncmranges - i - 1] = buildQACEInstance(fCentMultMin[ncmranges - i - 1], fCentMultMax[ncmranges - i - 1]);
        }
        if (doExtraAnalysis) {
          qaExtraDataCE[ncmranges - i - 1] = buildQACEExtraInstance(fCentMultMin[ncmranges - i - 1], fCentMultMax[ncmranges - i - 1]);
        }
        if (doPidAnalysis) {
          pidDataCE[ncmranges - i - 1] = buildPidCEInstance(fCentMultMin[ncmranges - i - 1], fCentMultMax[ncmranges - i - 1]);
        }
        if (doPidExtraAnalysis) {
          pidExtraDataCE[ncmranges - i - 1] = buildPidExtraCEInstance(fCentMultMin[ncmranges - i - 1], fCentMultMax[ncmranges - i - 1]);
        }
      }
      for (uint i = 0; i < ncmranges; ++i) {
        LOGF(info, " centrality/multipliicty range: %d, low limit: %0.2f, up limit: %0.2f", i, fCentMultMin[i], fCentMultMax[i]);
      }
    }
  }

  /// \brief Get the data collecting engine index corresponding to the passed collision
  template <typename FilteredCollision>
  int getDCEindex(FilteredCollision collision)
  {
    if (!useCentrality.value) {
      return 0;
    } else {
      int ixDCE = -1;
      float cm = collision.centmult();
      for (uint i = 0; i < ncmranges; ++i) {
        if (cm < fCentMultMax[i]) {
          ixDCE = i;
          break;
        }
      }
      if (!(ixDCE < 0)) {
        if (cm < fCentMultMin[ixDCE]) {
          ixDCE = -1;
        }
      }
      return ixDCE;
    }
  }

  template <typename FilteredCollisions, efficiencyandqatask::KindOfProcess kindOfProcess, efficiencyandqatask::KindOfData kindOfData, typename PassedTracks>
  void processTracks(FilteredCollisions::iterator const& collision, PassedTracks const& tracks)
  {
    using namespace efficiencyandqatask;

    int ixDCE = getDCEindex(collision);
    if (!(ixDCE < 0)) {
      if constexpr (kindOfProcess == kBASIC) {
        qaDataCE[ixDCE]->newCollision<kindOfData>();
      }
      if constexpr (kindOfProcess == kEXTRA) {
        qaExtraDataCE[ixDCE]->processTrackPairs<kindOfData, FilteredCollisions>(tracks, tracks);
      }
      for (auto const& track : tracks) {
        float tpcmom = track.p();
        float tofmom = track.p();
        if (useTPCInnerWallMomentum.value) {
          if constexpr (!framework::has_type_v<aod::mcparticle::PdgCode, typename PassedTracks::iterator::all_columns>) {
            tpcmom = track.tpcInnerParam();
          }
        }
        if constexpr (kindOfProcess == kBASIC) {
          qaDataCE[ixDCE]->processTrack<kindOfData, FilteredCollisions>(collision.posZ(), track);
        }
        if constexpr (kindOfProcess == kPID) {
          pidDataCE[ixDCE]->processTrack<kindOfData>(track, tpcmom, tofmom);
        }
        if constexpr (kindOfProcess == kPIDEXTRA) {
          pidExtraDataCE[ixDCE]->processTrack<kindOfData>(track, tpcmom, tofmom);
        }
      }
      if constexpr (kindOfProcess == kBASIC) {
        qaDataCE[ixDCE]->finishedCollision<kindOfData>();
      }
    }
  }

  void process(aod::Collisions const& collisions)
  {
    /* void function for alow processing the timestamp check task on faulty productions */
    LOGF(debug, "Received %d collisions", collisions.size());
  }

  using TpcPID = soa::Join<aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr>;
  using TofPID = soa::Join<aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::TOFSignal, aod::pidTOFbeta, aod::pidTOFmass>;

  Filter onlyacceptedcollisions = (aod::dptdptfilter::collisionaccepted == uint8_t(true));
  Filter onlyacceptedtracks = (aod::dptdptfilter::trackacceptedid >= int8_t(0));

  void processReconstructedNotStored(soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>::iterator const& collision,
                                     soa::Join<aod::FullTracks, aod::DptDptCFTracksInfo, aod::TrackSelection> const& tracks)
  {
    using namespace efficiencyandqatask;

    processTracks<soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>, kBASIC, kReco>(collision, tracks);
  }
  PROCESS_SWITCH(DptDptEfficiencyAndQc, processReconstructedNotStored, "Process reconstructed efficiency and QA for not stored derived data", false);

  void processDetectorLevelNotStored(soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>::iterator const& collision,
                                     soa::Join<aod::FullTracks, aod::DptDptCFTracksInfo, aod::TrackSelection, aod::McTrackLabels> const& tracks,
                                     soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo> const&)
  {
    using namespace efficiencyandqatask;

    processTracks<soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>, kBASIC, kReco>(collision, tracks);
  }
  PROCESS_SWITCH(DptDptEfficiencyAndQc, processDetectorLevelNotStored, "Process MC detector level efficiency and QA for not stored derived data", false);

  void processGeneratorLevelNotStored(soa::Filtered<soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo>>::iterator const& collision,
                                      soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo> const& particles)
  {
    using namespace efficiencyandqatask;

    processTracks<soa::Filtered<soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo>>, kBASIC, kGen>(collision, particles);
  }
  PROCESS_SWITCH(DptDptEfficiencyAndQc, processGeneratorLevelNotStored, "Process MC generator level efficiency and QA for not stored derived data", true);

  void processExtraReconstructedNotStored(soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>::iterator const& collision,
                                          soa::Filtered<soa::Join<aod::FullTracks, aod::DptDptCFTracksInfo, aod::TrackSelection>> const& tracks)
  {
    using namespace efficiencyandqatask;

    processTracks<soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>, kEXTRA, kReco>(collision, tracks);
  }
  PROCESS_SWITCH(DptDptEfficiencyAndQc, processExtraReconstructedNotStored, "Process reconstructed extra efficiency and QA for not stored derived data", false);

  void processExtraDetectorLevelNotStored(soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>::iterator const& collision,
                                          soa::Filtered<soa::Join<aod::FullTracks, aod::DptDptCFTracksInfo, aod::TrackSelection, aod::McTrackLabels>> const& tracks,
                                          soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo> const&)
  {
    using namespace efficiencyandqatask;

    processTracks<soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>, kEXTRA, kReco>(collision, tracks);
  }
  PROCESS_SWITCH(DptDptEfficiencyAndQc, processExtraDetectorLevelNotStored, "Process MC detector level extra efficiency and QA for not stored derived data", false);

  void processExtraGeneratorLevelNotStored(soa::Filtered<soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo>>::iterator const& collision,
                                           soa::Filtered<soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo>> const& particles)
  {
    using namespace efficiencyandqatask;

    processTracks<soa::Filtered<soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo>>, kEXTRA, kGen>(collision, particles);
  }
  PROCESS_SWITCH(DptDptEfficiencyAndQc, processExtraGeneratorLevelNotStored, "Process MC generator level extra efficiency and QA for not stored derived data", true);

  void processReconstructedNotStoredPID(soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>::iterator const& collision,
                                        soa::Join<aod::FullTracks, aod::DptDptCFTracksInfo, aod::TrackSelection, TpcPID, TofPID> const& tracks)
  {
    using namespace efficiencyandqatask;

    processTracks<soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>, kPID, kReco>(collision, tracks);
  }
  PROCESS_SWITCH(DptDptEfficiencyAndQc, processReconstructedNotStoredPID, "Process reconstructed PID QA for not stored derived data", false);

  void processReconstructedNotStoredPIDExtra(soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>::iterator const& collision,
                                             soa::Join<aod::FullTracks, aod::DptDptCFTracksInfo, aod::TrackSelection, TpcPID, TofPID> const& tracks)
  {
    using namespace efficiencyandqatask;

    processTracks<soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>, kPIDEXTRA, kReco>(collision, tracks);
  }
  PROCESS_SWITCH(DptDptEfficiencyAndQc, processReconstructedNotStoredPIDExtra, "Process reconstructed PID extra QA for not stored derived data", false);

  void processDetectorLevelNotStoredPID(soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>::iterator const& collision,
                                        soa::Join<aod::FullTracks, aod::DptDptCFTracksInfo, aod::TrackSelection, aod::McTrackLabels, TpcPID, TofPID> const& tracks,
                                        soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo> const&)
  {
    using namespace efficiencyandqatask;

    processTracks<soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>, kPID, kReco>(collision, tracks);
  }
  PROCESS_SWITCH(DptDptEfficiencyAndQc, processDetectorLevelNotStoredPID, "Process MC detector level PID QA for not stored derived data", false);

  void processDetectorLevelNotStoredTunedOnDataPID(soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>::iterator const& collision,
                                                   soa::Join<aod::FullTracks, aod::DptDptCFTracksInfo, aod::TrackSelection, aod::McTrackLabels, TpcPID, TofPID, aod::mcTPCTuneOnData> const& tracks,
                                                   soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo> const&)
  {
    using namespace efficiencyandqatask;

    processTracks<soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>, kPID, kReco>(collision, tracks);
  }
  PROCESS_SWITCH(DptDptEfficiencyAndQc, processDetectorLevelNotStoredTunedOnDataPID, "Process MC detector level tuned on data PID QA for not stored derived data", true);

  void processDetectorLevelNotStoredPIDExtra(soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>::iterator const& collision,
                                             soa::Join<aod::FullTracks, aod::DptDptCFTracksInfo, aod::TrackSelection, aod::McTrackLabels, TpcPID, TofPID> const& tracks,
                                             soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo> const&)
  {
    using namespace efficiencyandqatask;

    processTracks<soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>, kPIDEXTRA, kReco>(collision, tracks);
  }
  PROCESS_SWITCH(DptDptEfficiencyAndQc, processDetectorLevelNotStoredPIDExtra, "Process MC detector level PID extra QA for not stored derived data", false);

  void processDetectorLevelNotStoredTunedOnDataPIDExtra(soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>::iterator const& collision,
                                                        soa::Join<aod::FullTracks, aod::DptDptCFTracksInfo, aod::TrackSelection, aod::McTrackLabels, TpcPID, TofPID, aod::mcTPCTuneOnData> const& tracks,
                                                        soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo> const&)
  {
    using namespace efficiencyandqatask;

    processTracks<soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>, kPIDEXTRA, kReco>(collision, tracks);
  }
  PROCESS_SWITCH(DptDptEfficiencyAndQc, processDetectorLevelNotStoredTunedOnDataPIDExtra, "Process MC detector level tuned on data PID extra QA for not stored derived data", true);
};

using BCsWithTimestamps = soa::Join<aod::BCs, aod::Timestamps>;

struct CheckTimestamp {

  o2::ccdb::CcdbApi ccdbApi;
  int mRunNumber;
  uint64_t runsor = 0;
  uint64_t runeor = 0;
  std::shared_ptr<TH2> hTimeStampDiffNegative = nullptr;
  std::shared_ptr<TH2> hTimeStampDiffPositive = nullptr;

  Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    LOG(info) << "Initializing check timestamp task";
    AxisSpec diffAxis{10000, 0.0001, 1000000, "time (s)"};
    diffAxis.makeLogarithmic();

    hTimeStampDiffNegative = registry.add<TH2>("DiffNegative", "Time before SOR (s)", kTH2F, {{100, 0.5, 100.5, "Run number"}, diffAxis});
    hTimeStampDiffPositive = registry.add<TH2>("DiffPositive", "Time after EOR (s)", kTH2F, {{100, 0.5, 100.5, "Run number"}, diffAxis});
    ccdbApi.init(ccdburl);
  }

  void process(aod::Collisions const& collisions, BCsWithTimestamps const&)
  {
    for (auto const& collision : collisions) {
      /* check the previous run number */
      auto bc = collision.bc_as<BCsWithTimestamps>();
      if (bc.runNumber() != mRunNumber) {
        LOGF(info, "timestamp=%llu, run number=%d", bc.timestamp(), bc.runNumber());
        mRunNumber = bc.runNumber();

        // read SOR and EOR timestamps from RCT CCDB via utility function
        auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdbApi, mRunNumber, false);
        runeor = soreor.second;
        runsor = soreor.first;
      }
      if (bc.timestamp() < runsor || runeor < bc.timestamp()) {
        /* we got a wrong timestamp let's convert the out of run time to seconds */
        if (bc.timestamp() < runsor) {
          hTimeStampDiffNegative->Fill(Form("%d", mRunNumber), (runsor - bc.timestamp()) / 1000, 1);
        } else {
          hTimeStampDiffPositive->Fill(Form("%d", mRunNumber), (bc.timestamp() - runeor) / 1000, 1);
        }
      }
    }
  }
};

//****************************************************************************************
/**
 * Workflow definition.
 */
//****************************************************************************************
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<DptDptEfficiencyAndQc>(cfgc),
    adaptAnalysisTask<CheckTimestamp>(cfgc)};
  return workflow;
}
