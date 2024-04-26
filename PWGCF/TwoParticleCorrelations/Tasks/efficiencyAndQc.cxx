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

#include <TH2F.h>
#include <TProfile2D.h>
#include <CCDB/BasicCCDBManager.h>
#include "ReconstructionDataFormats/PID.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TableHelper.h"
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

namespace efficiencyandqatask
{
/// \enum KindOfProcessQA
/// \brief The kind of processing for templating the procedures
enum KindOfProcess {
  kReco = 0, ///< processing over reconstructed particles/tracks
  kGen       ///< processing over generated particles
};

/// \enum BeforeAfter
/// \brief The kind of filling, before or after track selection
enum BeforeAfter {
  kBefore = 0, ///< filling before track selection
  kAfter       ///< filling after track selection
};

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
static const std::vector<int> pdgcodes = {11, 13, 211, 321, 2212};
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
  std::vector<std::shared_ptr<TH2>> fhPt_vs_EtaB{2, nullptr};
  std::vector<std::shared_ptr<TH2>> fhPt_vs_ZvtxB{2, nullptr};
  std::vector<std::vector<std::shared_ptr<TH1>>> fhPtA{2, {nsp, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPt_vs_EtaA{2, {nsp, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPt_vs_ZvtxA{2, {nsp, nullptr}};
  std::shared_ptr<TH2> fhPt_vs_EtaItsAcc{nullptr};
  std::shared_ptr<TH2> fhPt_vs_EtaTpcAcc{nullptr};
  std::shared_ptr<TH2> fhPt_vs_EtaItsTpcAcc{nullptr};
  std::shared_ptr<TH2> fhPt_vs_EtaItsTofAcc{nullptr};
  std::shared_ptr<TH2> fhPt_vs_EtaTpcTofAcc{nullptr};
  std::shared_ptr<TH2> fhPt_vs_EtaItsTpcTofAcc{nullptr};
  std::vector<std::shared_ptr<TH2>> fhPt_vs_EtaItsA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhPt_vs_EtaTpcA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhPt_vs_EtaItsTpcA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhPt_vs_EtaItsTofA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhPt_vs_EtaTpcTofA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhPt_vs_EtaItsTpcTofA{nsp, nullptr};
  /* primaries and secondaries */
  /* overall, first index detector level second index generator level */
  /* detailed, first index detector level, second index associated particle */
  std::shared_ptr<TH3> fhPtPurityPosPrimA{nullptr};
  std::shared_ptr<TH3> fhPtPurityNegPrimA{nullptr};
  std::vector<std::shared_ptr<TH2>> fhPt_vs_EtaPrimA{nsp, nullptr};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPt_vs_EtaPrimItsA{2, {nsp, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPt_vs_EtaPrimItsTpcA{2, {nsp, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPt_vs_EtaPrimItsTpcTofA{2, {nsp, nullptr}};
  std::shared_ptr<TH3> fhPtPurityPosSecA{nullptr};
  std::shared_ptr<TH3> fhPtPurityNegSecA{nullptr};
  std::vector<std::shared_ptr<TH2>> fhPt_vs_EtaSecA{nsp, nullptr};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPt_vs_EtaSecItsA{2, {nsp, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPt_vs_EtaSecItsTpcA{2, {nsp, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPt_vs_EtaSecItsTpcTofA{2, {nsp, nullptr}};
  std::shared_ptr<TH3> fhPtPurityPosMatA{nullptr};
  std::shared_ptr<TH3> fhPtPurityNegMatA{nullptr};
  std::vector<std::shared_ptr<TH2>> fhPt_vs_EtaMatA{nsp, nullptr};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPt_vs_EtaMatItsA{2, {nsp, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPt_vs_EtaMatItsTpcA{2, {nsp, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPt_vs_EtaMatItsTpcTofA{2, {nsp, nullptr}};
  /* QC histograms */
  std::shared_ptr<TH2> fhITS_NCls_vs_PtB{nullptr};
  std::shared_ptr<TH2> fhITS_Chi2NCls_vs_PtB{nullptr};
  std::shared_ptr<TH2> fhTPC_FindableNCls_vs_PtB{nullptr};
  std::shared_ptr<TH2> fhTPC_FoundNCls_vs_PtB{nullptr};
  std::shared_ptr<TH2> fhTPC_SharedNCls_vs_PtB{nullptr};
  std::shared_ptr<TH2> fhTPC_FractionSharedCls_vs_PtB{nullptr};
  std::shared_ptr<TH2> fhTPC_CrossedRows_vs_PtB{nullptr};
  std::shared_ptr<TH2> fhTPC_CrossedRowsOverFindableCls_vs_PtB{nullptr};
  std::shared_ptr<TH2> fhTPC_Chi2NCls_vs_PtB{nullptr};
  std::vector<std::shared_ptr<TH2>> fhITS_NCls_vs_PtA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhITS_Chi2NCls_vs_PtA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhTPC_FindableNCls_vs_PtA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhTPC_FoundNCls_vs_PtA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhTPC_SharedNCls_vs_PtA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhTPC_FractionSharedCls_vs_PtA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhTPC_CrossedRows_vs_PtA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhTPC_CrossedRowsOverFindableCls_vs_PtA{nsp, nullptr};
  std::vector<std::shared_ptr<TH2>> fhTPC_Chi2NCls_vs_PtA{nsp, nullptr};

  template <efficiencyandqatask::KindOfProcess kind>
  void init(HistogramRegistry& registry, const char* dirname)
  {
    using namespace efficiencyandqatask;
    using namespace analysis::dptdptfilter;

    AxisSpec pidPAxis{150, 0.1, 5.0, "#it{p} (GeV/#it{c})"};
    pidPAxis.makeLogarithmic();
    const AxisSpec ptAxis{ptbins, ptlow, ptup, "#it{p}_{T} (GeV/c)"};
    const AxisSpec etaAxis{etabins, etalow, etaup, "#eta"};
    const AxisSpec phiAxis{360, 0.0f, constants::math::TwoPI, "#varphi"};
    const AxisSpec zvtxAxis{zvtxbins, zvtxlow, zvtxup, "#it{z}_{vtx}"};
    const AxisSpec itsNClsAxis{8, -0.5, 7.5, "ITS n clusters"};
    const AxisSpec itsCh2Axis{100, 0, 40, "#Chi^{2}/Cls ITS"};
    const AxisSpec tpcNClsAxis{165, -0.5, 164.5, "TPC n clusters"};
    const AxisSpec tpcNRowsAxis{165, -0.5, 164.5, "TPC n rows"};
    const AxisSpec tpcFractionAxis{100, 0, 1, "fraction"};
    const AxisSpec tpcXRowsOverFindClsAxis{60, 0.7, 1.3, "fraction"};
    const AxisSpec tpcCh2Axis{100, 0, 10, "#Chi^{2}/Cls TPC"};

    /* the reconstructed and generated levels histograms */
    std::string recogen = (kind == kReco) ? "Reco" : "Gen";
    fhPtB[kind] = ADDHISTOGRAM(TH1, DIRECTORYSTRING("%s/%s/%s", dirname, recogen.c_str(), "Before"), "Pt", "#it{p}_{T}", kTH1F, {ptAxis});
    fhPt_vs_EtaB[kind] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, recogen.c_str(), "Before"), "PtVsEta", "#it{p}_T vs #eta", kTH2F, {etaAxis, ptAxis});
    fhPt_vs_ZvtxB[kind] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, recogen.c_str(), "Before"), "PtVsZvtx", "#it{p}_T vs #it{z}_{vtx}", kTH2F, {zvtxAxis, ptAxis});
    for (uint isp = 0; isp < nsp; ++isp) {
      fhPtA[kind][isp] = ADDHISTOGRAM(TH1, DIRECTORYSTRING("%s/%s/%s", dirname, recogen.c_str(), "After"), HNAMESTRING("Pt_%s", tnames[isp].c_str()), HTITLESTRING("#it{p}_{T} %s", tnames[isp].c_str()), kTH1F, {ptAxis});
      fhPt_vs_EtaA[kind][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, recogen.c_str(), "After"), HNAMESTRING("PtVsEta_%s", tnames[isp].c_str()), HTITLESTRING("#it{p}_{T} vs #eta %s", tnames[isp].c_str()), kTH2F, {etaAxis, ptAxis});
      fhPt_vs_ZvtxA[kind][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, recogen.c_str(), "After"), HNAMESTRING("PtVsZvtx_%s", tnames[isp].c_str()), HTITLESTRING("#it{p}_{T} vs #it{z}_{zvtx} %s", tnames[isp].c_str()), kTH2F, {zvtxAxis, ptAxis});
    }

    if constexpr (kind == kReco) {
      /* only the reconstructed level histograms*/
      fhITS_NCls_vs_PtB = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "ITSNCls", "ITS clusters", kTH2F, {ptAxis, itsNClsAxis});
      fhITS_Chi2NCls_vs_PtB = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "ITSChi2NCls", "ITS #Chi^{2}", kTH2F, {ptAxis, itsCh2Axis});
      fhTPC_FindableNCls_vs_PtB = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "TPCFindableNCls", "TPC findable clusters", kTH2F, {ptAxis, tpcNClsAxis});
      fhTPC_FoundNCls_vs_PtB = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "TPCFoundNCls", "TPC found clusters", kTH2F, {ptAxis, tpcNClsAxis});
      fhTPC_SharedNCls_vs_PtB = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "TPCSharedNCls", "TPC shared clusters", kTH2F, {ptAxis, tpcNClsAxis});
      fhTPC_FractionSharedCls_vs_PtB = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "TPCFractionSharedCls", "TPC fraction shared clusters", kTH2F, {ptAxis, tpcFractionAxis});
      fhTPC_CrossedRows_vs_PtB = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "TPCXrows", "TPC crossed rows", kTH2F, {ptAxis, tpcNRowsAxis});
      fhTPC_CrossedRowsOverFindableCls_vs_PtB = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "XRowsOverFindableCls", "TPC xrows over findable clusters", kTH2F, {ptAxis, tpcXRowsOverFindClsAxis});
      fhTPC_Chi2NCls_vs_PtB = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "Before"), "TPCChi2NCls", "TPC #Chi^{2}", kTH2F, {ptAxis, tpcCh2Axis});
      /* efficiency histograms */
      fhPvsInnerP = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), "pVsInnerP", "#it{p} versus TPC inner wall #it{p}", kTH2F, {pidPAxis, pidPAxis});
      fhPt_vs_EtaItsAcc = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), "ptItsAcc", "ITS tracks within the acceptance", kTH2F, {etaAxis, ptAxis});
      fhPt_vs_EtaTpcAcc = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), "ptTpcAcc", "TPC tracks within the acceptance", kTH2F, {etaAxis, ptAxis});
      fhPt_vs_EtaItsTpcAcc = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), "ptItsTpcAcc", "ITS&TPC tracks within the acceptance", kTH2F, {etaAxis, ptAxis});
      fhPt_vs_EtaItsTofAcc = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), "ptItsTofAcc", "ITS&TOF tracks within the acceptance", kTH2F, {etaAxis, ptAxis});
      fhPt_vs_EtaTpcTofAcc = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), "ptTpcTofAcc", "TPC&TOF tracks within the acceptance", kTH2F, {etaAxis, ptAxis});
      fhPt_vs_EtaItsTpcTofAcc = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), "ptItsTpcTofAcc", "ITS&TPC&TOF tracks within the acceptance", kTH2F, {etaAxis, ptAxis});
      for (uint isp = 0; isp < nsp; ++isp) {
        fhITS_NCls_vs_PtA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "After"), HNAMESTRING("ITSNCls_%s", tnames[isp].c_str()), HTITLESTRING("ITS clusters %s", tnames[isp].c_str()), kTH2F, {ptAxis, itsNClsAxis});
        fhITS_Chi2NCls_vs_PtA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "After"), HNAMESTRING("ITSChi2NCls_%s", tnames[isp].c_str()), HTITLESTRING("ITS #Chi^{2} %s", tnames[isp].c_str()), kTH2F, {ptAxis, itsCh2Axis});
        fhTPC_FindableNCls_vs_PtA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "After"), HNAMESTRING("TPCFindableNCls_%s", tnames[isp].c_str()), HTITLESTRING("TPC findable clusters %s", tnames[isp].c_str()), kTH2F, {ptAxis, tpcNClsAxis});
        fhTPC_FoundNCls_vs_PtA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "After"), HNAMESTRING("TPCFoundNCls_%s", tnames[isp].c_str()), HTITLESTRING("TPC found clusters %s", tnames[isp].c_str()), kTH2F, {ptAxis, tpcNClsAxis});
        fhTPC_SharedNCls_vs_PtA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "After"), HNAMESTRING("TPCSharedNCls_%s", tnames[isp].c_str()), HTITLESTRING("TPC shared clusters %s", tnames[isp].c_str()), kTH2F, {ptAxis, tpcNClsAxis});
        fhTPC_FractionSharedCls_vs_PtA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "After"), HNAMESTRING("TPCFractionSharedCls_%s", tnames[isp].c_str()), HTITLESTRING("TPC fraction shared clusters %s", tnames[isp].c_str()), kTH2F, {ptAxis, tpcFractionAxis});
        fhTPC_CrossedRows_vs_PtA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "After"), HNAMESTRING("TPCXrows_%s", tnames[isp].c_str()), HTITLESTRING("TPC crossed rows %s", tnames[isp].c_str()), kTH2F, {ptAxis, tpcNRowsAxis});
        fhTPC_CrossedRowsOverFindableCls_vs_PtA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "After"), HNAMESTRING("XRowsOverFindableCls_%s", tnames[isp].c_str()), HTITLESTRING("TPC xrows over findable clusters %s", tnames[isp].c_str()), kTH2F, {ptAxis, tpcXRowsOverFindClsAxis});
        fhTPC_Chi2NCls_vs_PtA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Reco", "After"), HNAMESTRING("TPCChi2NCls_%s", tnames[isp].c_str()), HTITLESTRING("TPC #Chi^{2} %s", tnames[isp].c_str()), kTH2F, {ptAxis, tpcCh2Axis});
        /* efficiency histograms */
        fhPt_vs_EtaItsA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), HNAMESTRING("ptIts_%s", tnames[isp].c_str()), HTITLESTRING("ITS %s tracks", tnames[isp].c_str()), kTH2F, {etaAxis, ptAxis});
        fhPt_vs_EtaTpcA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), HNAMESTRING("ptTpc_%s", tnames[isp].c_str()), HTITLESTRING("TPC %s tracks", tnames[isp].c_str()), kTH2F, {etaAxis, ptAxis});
        fhPt_vs_EtaItsTpcA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), HNAMESTRING("ptItsTpc_%s", tnames[isp].c_str()), HTITLESTRING("ITS&TPC %s tracks", tnames[isp].c_str()), kTH2F, {etaAxis, ptAxis});
        fhPt_vs_EtaItsTofA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), HNAMESTRING("ptItsTof_%s", tnames[isp].c_str()), HTITLESTRING("ITS&TOF %s tracks", tnames[isp].c_str()), kTH2F, {etaAxis, ptAxis});
        fhPt_vs_EtaTpcTofA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), HNAMESTRING("ptTpcTof_%s", tnames[isp].c_str()), HTITLESTRING("TPC&TOF %s tracks", tnames[isp].c_str()), kTH2F, {etaAxis, ptAxis});
        fhPt_vs_EtaItsTpcTofA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Reco"), HNAMESTRING("ptItsTpcTof_%s", tnames[isp].c_str()), HTITLESTRING("ITS&TPC&TOF %s tracks", tnames[isp].c_str()), kTH2F, {etaAxis, ptAxis});
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
        fhPt_vs_EtaPrimA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Gen"),
                                             HNAMESTRING("ptPrim%s", tnames[isp].c_str()),
                                             HTITLESTRING("ITS  %s tracks (primaries)", tnames[isp].c_str()),
                                             kTH2F, {etaAxis, ptAxis});
        fhPt_vs_EtaSecA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Gen"),
                                            HNAMESTRING("ptSec%s", tnames[isp].c_str()),
                                            HTITLESTRING("ITS %s tracks (secondaries)", tnames[isp].c_str()),
                                            kTH2F, {etaAxis, ptAxis});
        fhPt_vs_EtaMatA[isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", "Gen"),
                                            HNAMESTRING("ptMat%s", tnames[isp].c_str()),
                                            HTITLESTRING("ITS %s tracks (from material)", tnames[isp].c_str()),
                                            kTH2F, {etaAxis, ptAxis});

        const std::vector<std::string> detectedorigin = {"DetReco", "DetAssoc"};
        for (uint ix = 0; ix < detectedorigin.size(); ++ix) {
          fhPt_vs_EtaPrimItsA[ix][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", detectedorigin[ix].c_str()),
                                                      HNAMESTRING("ptItsPrim_%s", tnames[isp].c_str()),
                                                      HTITLESTRING("ITS %s tracks (primaries)", tnames[isp].c_str()),
                                                      kTH2F, {etaAxis, ptAxis});
          fhPt_vs_EtaPrimItsTpcA[ix][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", detectedorigin[ix].c_str()),
                                                         HNAMESTRING("ptItsTpcPrim_%s", tnames[isp].c_str()),
                                                         HTITLESTRING("ITS&TPC %s tracks (primaries)", tnames[isp].c_str()),
                                                         kTH2F, {etaAxis, ptAxis});
          fhPt_vs_EtaPrimItsTpcTofA[ix][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", detectedorigin[ix].c_str()),
                                                            HNAMESTRING("ptItsTpcTofPrim_%s", tnames[isp].c_str()),
                                                            HTITLESTRING("ITS&TPC&TOF %s tracks (primaries)", tnames[isp].c_str()),
                                                            kTH2F, {etaAxis, ptAxis});
          fhPt_vs_EtaSecItsA[ix][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", detectedorigin[ix].c_str()),
                                                     HNAMESTRING("ptItsSec_%s", tnames[isp].c_str()),
                                                     HTITLESTRING("ITS %s tracks (secondaries)", tnames[isp].c_str()),
                                                     kTH2F, {etaAxis, ptAxis});
          fhPt_vs_EtaSecItsTpcA[ix][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", detectedorigin[ix].c_str()),
                                                        HNAMESTRING("ptItsTpcSec_%s", tnames[isp].c_str()),
                                                        HTITLESTRING("ITS&TPC %s tracks (secondaries)", tnames[isp].c_str()),
                                                        kTH2F, {etaAxis, ptAxis});
          fhPt_vs_EtaSecItsTpcTofA[ix][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", detectedorigin[ix].c_str()),
                                                           HNAMESTRING("ptItsTpcTofSec_%s", tnames[isp].c_str()),
                                                           HTITLESTRING("ITS&TPC&TOF %s tracks (secondaries)", tnames[isp].c_str()),
                                                           kTH2F, {etaAxis, ptAxis});
          fhPt_vs_EtaMatItsA[ix][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", detectedorigin[ix].c_str()),
                                                     HNAMESTRING("ptItsMat_%s", tnames[isp].c_str()),
                                                     HTITLESTRING("ITS %s tracks (from material)", tnames[isp].c_str()),
                                                     kTH2F, {etaAxis, ptAxis});
          fhPt_vs_EtaMatItsTpcA[ix][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", detectedorigin[ix].c_str()),
                                                        HNAMESTRING("ptItsTpcMat_%s", tnames[isp].c_str()),
                                                        HTITLESTRING("ITS&TPC %s tracks (from material)", tnames[isp].c_str()),
                                                        kTH2F, {etaAxis, ptAxis});
          fhPt_vs_EtaMatItsTpcTofA[ix][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "Efficiency", detectedorigin[ix].c_str()),
                                                           HNAMESTRING("ptItsTpcTofMat_%s", tnames[isp].c_str()),
                                                           HTITLESTRING("ITS&TPC&TOF %s tracks (from material)", tnames[isp].c_str()),
                                                           kTH2F, {etaAxis, ptAxis});
        }
      }
    }
  }

  template <efficiencyandqatask::KindOfProcess kind, typename TrackObject>
  void processTrack(float zvtx, TrackObject const& track)
  {
    using namespace efficiencyandqatask;
    using namespace analysis::dptdptfilter;
    using namespace o2::aod::track;

    fhPtB[kind]->Fill(track.pt());
    fhPt_vs_EtaB[kind]->Fill(track.eta(), track.pt());
    fhPt_vs_ZvtxB[kind]->Fill(zvtx, track.pt());
    if (!(track.trackacceptedid() < 0)) {
      fhPtA[kind][track.trackacceptedid()]->Fill(track.pt());
      fhPt_vs_EtaA[kind][track.trackacceptedid()]->Fill(track.eta(), track.pt());
      fhPt_vs_ZvtxA[kind][track.trackacceptedid()]->Fill(zvtx, track.pt());
    }
    if constexpr (kind == kReco) {
      auto fillhisto = [&track](auto& h, bool cond) {
        if (cond) {
          h->Fill(track.eta(), track.pt());
        }
      };
      bool hasits = track.hasITS() && TrackSelectionFlags::checkFlag(track.trackCutFlag(), trackSelectionITS);
      bool hastpc = track.hasTPC() && TrackSelectionFlags::checkFlag(track.trackCutFlag(), trackSelectionTPC);
      bool hastof = track.hasTOF();

      fhITS_NCls_vs_PtB->Fill(track.pt(), track.itsNCls());
      fhITS_Chi2NCls_vs_PtB->Fill(track.pt(), track.itsChi2NCl());
      fhTPC_FindableNCls_vs_PtB->Fill(track.pt(), track.tpcNClsFindable());
      fhTPC_FoundNCls_vs_PtB->Fill(track.pt(), track.tpcNClsFound());
      fhTPC_SharedNCls_vs_PtB->Fill(track.pt(), track.tpcNClsShared());
      fhTPC_FractionSharedCls_vs_PtB->Fill(track.pt(), track.tpcFractionSharedCls());
      fhTPC_CrossedRows_vs_PtB->Fill(track.pt(), track.tpcNClsCrossedRows());
      fhTPC_CrossedRowsOverFindableCls_vs_PtB->Fill(track.pt(), track.tpcCrossedRowsOverFindableCls());
      fhTPC_Chi2NCls_vs_PtB->Fill(track.pt(), track.tpcChi2NCl());
      if (InTheAcceptance(track)) {
        /* efficiency histograms */
        fillhisto(fhPt_vs_EtaItsAcc, hasits);
        fillhisto(fhPt_vs_EtaTpcAcc, hastpc);
        fillhisto(fhPt_vs_EtaItsTpcAcc, hasits && hastpc);
        fillhisto(fhPt_vs_EtaItsTofAcc, hasits && hastof);
        fillhisto(fhPt_vs_EtaTpcTofAcc, hastpc && hastof);
        fillhisto(fhPt_vs_EtaItsTpcTofAcc, hasits && hastpc && hastof);
      }
      if (!(track.trackacceptedid() < 0)) {
        fhITS_NCls_vs_PtA[track.trackacceptedid()]->Fill(track.pt(), track.itsNCls());
        fhITS_Chi2NCls_vs_PtA[track.trackacceptedid()]->Fill(track.pt(), track.itsChi2NCl());
        fhTPC_FindableNCls_vs_PtA[track.trackacceptedid()]->Fill(track.pt(), track.tpcNClsFindable());
        fhTPC_FoundNCls_vs_PtA[track.trackacceptedid()]->Fill(track.pt(), track.tpcNClsFound());
        fhTPC_SharedNCls_vs_PtA[track.trackacceptedid()]->Fill(track.pt(), track.tpcNClsShared());
        fhTPC_FractionSharedCls_vs_PtA[track.trackacceptedid()]->Fill(track.pt(), track.tpcFractionSharedCls());
        fhTPC_CrossedRows_vs_PtA[track.trackacceptedid()]->Fill(track.pt(), track.tpcNClsCrossedRows());
        fhTPC_CrossedRowsOverFindableCls_vs_PtA[track.trackacceptedid()]->Fill(track.pt(), track.tpcCrossedRowsOverFindableCls());
        fhTPC_Chi2NCls_vs_PtA[track.trackacceptedid()]->Fill(track.pt(), track.tpcChi2NCl());
        /* efficiency histograms */
        fhPvsInnerP->Fill(track.tpcInnerParam(), track.p());
        fillhisto(fhPt_vs_EtaItsA[track.trackacceptedid()], hasits);
        fillhisto(fhPt_vs_EtaTpcA[track.trackacceptedid()], hastpc);
        fillhisto(fhPt_vs_EtaItsTpcA[track.trackacceptedid()], hasits && hastpc);
        fillhisto(fhPt_vs_EtaItsTofA[track.trackacceptedid()], hasits && hastof);
        fillhisto(fhPt_vs_EtaTpcTofA[track.trackacceptedid()], hastpc && hastof);
        fillhisto(fhPt_vs_EtaItsTpcTofA[track.trackacceptedid()], hasits && hastpc && hastof);
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
          std::vector<float> t_pt = {track.pt(), mcparticle.pt()};
          std::vector<float> t_eta = {track.eta(), mcparticle.eta()};
          for (uint ix = 0; ix < t_pt.size(); ++ix) {
            fillhisto(fhPt_vs_EtaPrimItsA[ix][track.trackacceptedid()], t_pt[ix], t_eta[ix], hasits, isprimary);
            fillhisto(fhPt_vs_EtaPrimItsTpcA[ix][track.trackacceptedid()], t_pt[ix], t_eta[ix], hasits && hastpc, isprimary);
            fillhisto(fhPt_vs_EtaPrimItsTpcTofA[ix][track.trackacceptedid()], t_pt[ix], t_eta[ix], hasits && hastof, isprimary);
            fillhisto(fhPt_vs_EtaSecItsA[ix][track.trackacceptedid()], t_pt[ix], t_eta[ix], hasits, issecdecay);
            fillhisto(fhPt_vs_EtaSecItsTpcA[ix][track.trackacceptedid()], t_pt[ix], t_eta[ix], hasits && hastpc, issecdecay);
            fillhisto(fhPt_vs_EtaSecItsTpcTofA[ix][track.trackacceptedid()], t_pt[ix], t_eta[ix], hasits && hastof, issecdecay);
            fillhisto(fhPt_vs_EtaMatItsA[ix][track.trackacceptedid()], t_pt[ix], t_eta[ix], hasits, isfrommaterial);
            fillhisto(fhPt_vs_EtaMatItsTpcA[ix][track.trackacceptedid()], t_pt[ix], t_eta[ix], hasits && hastpc, isfrommaterial);
            fillhisto(fhPt_vs_EtaMatItsTpcTofA[ix][track.trackacceptedid()], t_pt[ix], t_eta[ix], hasits && hastof, isfrommaterial);
          }
        }
      }
    }
    if constexpr (kind == kGen) {
      if (!(track.trackacceptedid() < 0)) {
        /* pure generator level */
        if (track.isPhysicalPrimary()) {
          fhPt_vs_EtaPrimA[track.trackacceptedid()]->Fill(track.eta(), track.pt());
        } else if (track.getProcess() == 4) {
          fhPt_vs_EtaSecA[track.trackacceptedid()]->Fill(track.eta(), track.pt());
        } else {
          fhPt_vs_EtaMatA[track.trackacceptedid()]->Fill(track.eta(), track.pt());
        }
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

  template <efficiencyandqatask::KindOfProcess kind>
  void init(HistogramRegistry& registry, const char* dirname)
  {
    using namespace efficiencyandqatask;

    const AxisSpec dEdxAxis{200, 0.0, 200.0, "dE/dx (au)"};
    AxisSpec pidPAxis{150, 0.1, 5.0, "#it{p} (GeV/#it{c})"};
    pidPAxis.makeLogarithmic();

    if constexpr (kind == kReco) {
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
                                                   kTH3F, {pidPAxis, {120, -6.0, 6.0, FORMATSTRING("n#sigma_{TPC}^{%s}", mainsptitles[isp].c_str())}, {120, -6.0, 6.0, FORMATSTRING("n#sigma_{TOF}^{%s}", mainsptitles[isp].c_str())}});
        }
        for (uint isp = 0; isp < nallmainsp; ++isp) {
          fhTPCnSigmasVsP[ix][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "PID", whenname[ix].c_str()),
                                                  HNAMESTRING("tpcNSigmasVsP%c_%s", whenprefix[ix], allmainspnames[isp].c_str()),
                                                  HTITLESTRING("TPC n#sigma to the %s line %s", allmainsptitles[isp].c_str(), whentitle[ix].c_str()),
                                                  kTH2F, {pidPAxis, {120, -6.0, 6.0, FORMATSTRING("n#sigma_{TPC}^{%s}", allmainsptitles[isp].c_str())}});
          fhTOFnSigmasVsP[ix][isp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "PID", whenname[ix].c_str()),
                                                  HNAMESTRING("tofNSigmasVsP%c_%s", whenprefix[ix], allmainspnames[isp].c_str()),
                                                  HTITLESTRING("TOF n#sigma to the %s line %s", allmainsptitles[isp].c_str(), whentitle[ix].c_str()),
                                                  kTH2F, {pidPAxis, {120, -6.0, 6.0, FORMATSTRING("n#sigma_{TOF}^{%s}", allmainsptitles[isp].c_str())}});
        }
      }
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
                                                           kTH2F, {pidPAxis, {120, -6.0, 6.0, FORMATSTRING("n#sigma_{TPC}^{%s}", allmainsptitles[isp].c_str())}});
            fhIdTOFnSigmasVsP[isp][imainsp] = ADDHISTOGRAM(TH2, DIRECTORYSTRING("%s/%s/%s", dirname, "PID", "Selected"),
                                                           HNAMESTRING("tofNSigmasVsPSelected_%s_to%s", tnames[isp].c_str(), allmainspnames[imainsp].c_str()),
                                                           HTITLESTRING("TOF n#sigma for selected %s to the %s line", tnames[isp].c_str(), allmainsptitles[imainsp].c_str()),
                                                           kTH2F, {pidPAxis, {120, -6.0, 6.0, FORMATSTRING("n#sigma_{TOF}^{%s}", allmainsptitles[isp].c_str())}});
          }
        }
      }
    }
  }

  template <o2::track::PID::ID id, typename TrackObject>
  void fillAllSpeciesPID(uint ix, TrackObject const& track, float mom)
  {
    if (track.sign() < 0) {
      ix = 2 * ix + 1;
    } else {
      ix = 2 * ix;
    }
    for (uint when = 0; when < 2; ++when) {
      fhTPCnSigmasVsP[when][ix]->Fill(mom, o2::aod::pidutils::tpcNSigma<id>(track));
      fhTOFnSigmasVsP[when][ix]->Fill(mom, o2::aod::pidutils::tofNSigma<id>(track));
      if (track.trackacceptedid() < 0) {
        /* track not accepted */
        return;
      }
    }
    fhIdTPCnSigmasVsP[track.trackacceptedid()][ix]->Fill(mom, o2::aod::pidutils::tpcNSigma<id>(track));
    fhIdTOFnSigmasVsP[track.trackacceptedid()][ix]->Fill(mom, o2::aod::pidutils::tofNSigma<id>(track));
    if (efficiencyandqatask::pidselector.isGlobalSpecies(track.trackacceptedid() / 2, id)) {
      /* only if the species of the selected track matches the target of the number of sigmas */
      fpIdTPCdEdxSignalVsPSigmas[track.trackacceptedid()]->Fill(mom, track.tpcSignal(), o2::aod::pidutils::tpcNSigma<id>(track));
      fpIdTOFSignalVsPSigmas[track.trackacceptedid()]->Fill(mom, track.beta(), o2::aod::pidutils::tofNSigma<id>(track));
    }
  }

  template <o2::track::PID::ID id, typename TrackObject>
  void fillSpeciesPID(uint ix, TrackObject const& track, float mom)
  {
    if (track.sign() < 0) {
      ix = 2 * ix + 1;
    } else {
      ix = 2 * ix;
    }
    for (uint when = 0; when < 2; ++when) {
      fhTPCdEdxSignalDiffVsP[when][ix]->Fill(mom, o2::aod::pidutils::tpcExpSignalDiff<id>(track));
      fhTOFSignalDiffVsP[when][ix]->Fill(mom, o2::aod::pidutils::tofExpSignalDiff<id>(track));
      fhTPCTOFSigmaVsP[when][ix]->Fill(mom, o2::aod::pidutils::tpcNSigma<id>(track), o2::aod::pidutils::tofNSigma<id>(track));
      if (track.trackacceptedid() < 0) {
        /* track not accepted */
        return;
      }
    }
  }

  template <typename TrackObject>
  void fillPID(TrackObject const& track, float mom)
  {
    for (uint when = 0; when < 2; ++when) {
      fhTPCdEdxSignalVsP[when]->Fill(mom, track.tpcSignal());
      fhTOFSignalVsP[when]->Fill(mom, track.beta());
      fhPvsTOFSqMass[when]->Fill(track.mass() * track.mass(), mom);
      if (track.trackacceptedid() < 0) {
        /* track not accepted */
        return;
      }
    }
    fhIdTPCdEdxSignalVsP[track.trackacceptedid()]->Fill(mom, track.tpcSignal());
    fhIdTOFSignalVsP[track.trackacceptedid()]->Fill(mom, track.beta());
  }

  template <efficiencyandqatask::KindOfProcess kind, typename TrackObject>
  void processTrack(TrackObject const& track, float mom)
  {
    using namespace efficiencyandqatask;

    if constexpr (kind == kReco) {
      fillPID(track, mom);
      fillSpeciesPID<o2::track::PID::Pion>(0, track, mom);
      fillSpeciesPID<o2::track::PID::Kaon>(1, track, mom);
      fillSpeciesPID<o2::track::PID::Proton>(2, track, mom);
      fillAllSpeciesPID<o2::track::PID::Electron>(0, track, mom);
      fillAllSpeciesPID<o2::track::PID::Muon>(1, track, mom);
      fillAllSpeciesPID<o2::track::PID::Pion>(2, track, mom);
      fillAllSpeciesPID<o2::track::PID::Kaon>(3, track, mom);
      fillAllSpeciesPID<o2::track::PID::Proton>(4, track, mom);
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
  PidDataCollectingEngine** pidDataCE;

  /* the histogram registries */
  HistogramRegistry registry_one{"registry_one", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registry_two{"registry_two", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registry_three{"registry_three", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registry_four{"registry_four", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registry_five{"registry_five", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registry_six{"registry_six", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registry_seven{"registry_seven", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registry_eight{"registry_eight", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registry_nine{"registry_nine", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registry_ten{"registry_ten", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registry_pidone{"pidregistry_one", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registry_pidtwo{"pidregistry_two", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registry_pidthree{"pidregistry_three", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registry_pidfour{"pidregistry_four", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registry_pidfive{"pidregistry_five", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registry_pidsix{"pidregistry_six", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registry_pidseven{"pidregistry_seven", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registry_pideight{"pidregistry_eight", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registry_pidnine{"pidregistry_nine", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry registry_pidten{"pidregistry_ten", {}, OutputObjHandlingPolicy::AnalysisObject};
  std::vector<HistogramRegistry*> registrybank{&registry_one, &registry_two, &registry_three, &registry_four, &registry_five,
                                               &registry_six, &registry_seven, &registry_eight, &registry_nine, &registry_ten};
  std::vector<HistogramRegistry*> pidregistrybank{&registry_pidone, &registry_pidtwo, &registry_pidthree, &registry_pidfour, &registry_pidfive,
                                                  &registry_pidsix, &registry_pidseven, &registry_pideight, &registry_pidnine, &registry_pidten};

  Configurable<bool> inCentralityClasses{"usecentrality", false, "Perform the task using centrality/multiplicity classes. Default value: false"};
  Configurable<bool> useTPCInnerWallMomentum{"useinnerwallmomm", false, "Use the TPC inner wall momentum. Default: false"};

  void init(o2::framework::InitContext& initContext)
  {
    using namespace efficiencyandqatask;
    using namespace analysis::dptdptfilter;

    /* do nothing if not active */
    if (!doprocessDetectorLevelNotStored && !doprocessDetectorLevelNotStoredNoPID && !doprocessGeneratorLevelNotStored && !doprocessReconstructedNotStored && !doprocessReconstructedNotStoredNoPID) {
      return;
    }

    fPDG = TDatabasePDG::Instance();
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
              pidselector.Add(spid, cfg);
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
      if (inCentralityClasses.value && getTaskOptionValue(initContext, "dpt-dpt-filter", "centralities", centspec, false)) {
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
      bool doPidAnalysis = doprocessDetectorLevelNotStored || doprocessReconstructedNotStored;
      qaDataCE = new QADataCollectingEngine*[ncmranges];
      if (doPidAnalysis) {
        pidDataCE = new PidDataCollectingEngine*[ncmranges];
      }
      std::string recogen;
      if (!(doprocessReconstructedNotStored || doprocessReconstructedNotStoredNoPID) && !(doprocessDetectorLevelNotStored || doprocessDetectorLevelNotStoredNoPID)) {
        LOGF(fatal, "Neither reco nor detector level not configured. Please, fix it!");
      }
      if (ncmranges > registrybank.size()) {
        LOGF(fatal, "There are more centrality ranges configured than registries in the bank. Please fix it!");
      }
      /* in reverse order for proper order in results file */
      for (uint i = 0; i < ncmranges; ++i) {
        auto initializeCEInstance = [&](auto dce, auto name, auto& registry) {
          /* crete the output list for the passed centrality/multiplicity range */
          /* init the data collection instance */
          dce->template init<kReco>(registry, name.Data());
          if (doprocessGeneratorLevelNotStored) {
            dce->template init<kGen>(registry, name.Data());
          }
        };
        auto buildQACEInstance = [&](float min, float max) {
          auto* dce = new QADataCollectingEngine();
          initializeCEInstance(dce, TString::Format("EfficiencyAndQaData-%d-%d", static_cast<int>(min), static_cast<int>(max)), *registrybank[i]);
          return dce;
        };
        auto buildPidCEInstance = [&](float min, float max) {
          auto* dce = new PidDataCollectingEngine();
          initializeCEInstance(dce, TString::Format("EfficiencyAndPidData-%d-%d", static_cast<int>(min), static_cast<int>(max)), *pidregistrybank[i]);
          return dce;
        };
        /* in reverse order for proper order in results file */
        qaDataCE[ncmranges - i - 1] = buildQACEInstance(fCentMultMin[ncmranges - i - 1], fCentMultMax[ncmranges - i - 1]);
        if (doPidAnalysis) {
          pidDataCE[ncmranges - i - 1] = buildPidCEInstance(fCentMultMin[ncmranges - i - 1], fCentMultMax[ncmranges - i - 1]);
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
    if (!inCentralityClasses.value) {
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

  template <bool dopid, efficiencyandqatask::KindOfProcess kind, typename FilterdCollision, typename PassedTracks>
  void processTracks(FilterdCollision const& collision, PassedTracks const& tracks)
  {
    using namespace efficiencyandqatask;

    int ixDCE = getDCEindex(collision);
    if (!(ixDCE < 0)) {
      for (auto& track : tracks) {
        float mom = track.p();
        if (useTPCInnerWallMomentum.value) {
          if constexpr (!framework::has_type_v<aod::mcparticle::PdgCode, typename PassedTracks::iterator::all_columns>) {
            mom = track.tpcInnerParam();
          }
        }
        qaDataCE[ixDCE]->processTrack<kind>(collision.posZ(), track);
        if constexpr (dopid) {
          pidDataCE[ixDCE]->processTrack<kind>(track, mom);
        }
      }
    }
  }

  void process(aod::Collisions const& collisions)
  {
    /* void function for alow processing the timestamp check task on faulty productions */
    LOGF(debug, "Received %d collisions", collisions.size());
  }

  using tpcPID = soa::Join<aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr>;
  using tofPID = soa::Join<aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::TOFSignal, aod::pidTOFbeta, aod::pidTOFmass>;

  Filter onlyacceptedcollisions = (aod::dptdptfilter::collisionaccepted == uint8_t(true));

  void processReconstructedNotStored(soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>::iterator const& collision,
                                     soa::Join<aod::FullTracks, aod::DptDptCFTracksInfo, aod::TrackSelection, tpcPID, tofPID>& tracks)
  {
    using namespace efficiencyandqatask;

    processTracks<true, kReco>(collision, tracks);
  }
  PROCESS_SWITCH(DptDptEfficiencyAndQc, processReconstructedNotStored, "Process reconstructed efficiency and QA for not stored derived data", false);

  void processDetectorLevelNotStored(soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>::iterator const& collision,
                                     soa::Join<aod::FullTracks, aod::DptDptCFTracksInfo, aod::TrackSelection, tpcPID, tofPID, aod::McTrackLabels>& tracks,
                                     soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo> const&)
  {
    using namespace efficiencyandqatask;

    processTracks<true, kReco>(collision, tracks);
  }
  PROCESS_SWITCH(DptDptEfficiencyAndQc, processDetectorLevelNotStored, "Process MC detector level efficiency and QA for not stored derived data", false);

  void processGeneratorLevelNotStored(soa::Filtered<soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo>>::iterator const& collision,
                                      soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo>& particles)
  {
    using namespace efficiencyandqatask;

    processTracks<false, kGen>(collision, particles);
  }
  PROCESS_SWITCH(DptDptEfficiencyAndQc, processGeneratorLevelNotStored, "Process MC generator level efficiency and QA for not stored derived data", true);

  void processReconstructedNotStoredNoPID(soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>::iterator const& collision,
                                          soa::Join<aod::FullTracks, aod::DptDptCFTracksInfo, aod::TrackSelection>& tracks)
  {
    using namespace efficiencyandqatask;

    processTracks<false, kReco>(collision, tracks);
  }
  PROCESS_SWITCH(DptDptEfficiencyAndQc, processReconstructedNotStoredNoPID, "Process reconstructed efficiency and QA for not stored derived data", false);

  void processDetectorLevelNotStoredNoPID(soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>::iterator const& collision,
                                          soa::Join<aod::FullTracks, aod::DptDptCFTracksInfo, aod::TrackSelection, aod::McTrackLabels>& tracks,
                                          soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo> const&)
  {
    using namespace efficiencyandqatask;

    processTracks<false, kReco>(collision, tracks);
  }
  PROCESS_SWITCH(DptDptEfficiencyAndQc, processDetectorLevelNotStoredNoPID, "Process MC detector level efficiency and QA for not stored derived data", true);
};

using BCsWithTimestamps = soa::Join<aod::BCs, aod::Timestamps>;

struct checkTimestamp {

  o2::ccdb::CcdbApi ccdb_api;
  int mRunNumber;
  ULong64_t runsor = 0;
  ULong64_t runeor = 0;
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
    ccdb_api.init(ccdburl);
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
        auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdb_api, mRunNumber, false);
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
    adaptAnalysisTask<checkTimestamp>(cfgc)};
  return workflow;
}
