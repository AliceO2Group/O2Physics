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

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TableHelper.h"
#include "Common/DataModel/TrackSelectionTables.h"
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

namespace efficiencyandqatask
{
// initialized during self configuration
int ptbins = 0;
float ptlow = 0.0f;
float ptup = 0.0f;
int etabins = 0;
float etalow = 0.0f;
float etaup = 0.0f;
int zvtxbins = 0;
float zvtxlow = 0.0f;
float zvtxup = 0.0f;

/// \enum KindOfProcess
/// \brief The kind of processing for templating the procedures
enum KindOfProcess {
  kReco = 0, ///< processing over reconstructed particles/tracks
  kGen       ///< processing over generated particles
};

// initialized during self configuration
std::vector<std::string> poinames; ///< the species of interest names
std::vector<std::string> tnames;   ///< the track names
} // namespace efficiencyandqatask

/* the QA data collecting engine */
struct QADataCollectingEngine {
  size_t nch = efficiencyandqatask::tnames.size();

  //===================================================
  // The QA output objects
  //===================================================
  /* efficiency histograms histograms */
  /* when two indexes, first index reco and detector level, second index generator level */
  /* when no indexes, reco and detector level */
  std::vector<std::shared_ptr<TH1>> fhPtB{2, nullptr};
  std::vector<std::shared_ptr<TH2>> fhPt_vs_EtaB{2, nullptr};
  std::vector<std::shared_ptr<TH2>> fhPt_vs_ZvtxB{2, nullptr};
  std::vector<std::vector<std::shared_ptr<TH1>>> fhPtA{2, {nch, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPt_vs_EtaA{2, {nch, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPt_vs_ZvtxA{2, {nch, nullptr}};
  std::vector<std::shared_ptr<TH2>> fhPt_vs_EtaItsA{nch, nullptr};
  std::vector<std::shared_ptr<TH2>> fhPt_vs_EtaTpcA{nch, nullptr};
  std::vector<std::shared_ptr<TH2>> fhPt_vs_EtaItsTpcA{nch, nullptr};
  std::vector<std::shared_ptr<TH2>> fhPt_vs_EtaItsTofA{nch, nullptr};
  std::vector<std::shared_ptr<TH2>> fhPt_vs_EtaTpcTofA{nch, nullptr};
  std::vector<std::shared_ptr<TH2>> fhPt_vs_EtaItsTpcTofA{nch, nullptr};
  /* primaries and secondaries */
  /* overall, first index detector level second index generator level */
  /* detailed, first index detector level, second index associated particle */
  std::vector<std::shared_ptr<TH2>> fhPt_vs_EtaPrimA{nch, nullptr};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPt_vs_EtaPrimItsA{2, {nch, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPt_vs_EtaPrimItsTpcA{2, {nch, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPt_vs_EtaPrimItsTpcTofA{2, {nch, nullptr}};
  std::vector<std::shared_ptr<TH2>> fhPt_vs_EtaSecA{nch, nullptr};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPt_vs_EtaSecItsA{2, {nch, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPt_vs_EtaSecItsTpcA{2, {nch, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPt_vs_EtaSecItsTpcTofA{2, {nch, nullptr}};
  std::vector<std::shared_ptr<TH2>> fhPt_vs_EtaMatA{nch, nullptr};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPt_vs_EtaMatItsA{2, {nch, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPt_vs_EtaMatItsTpcA{2, {nch, nullptr}};
  std::vector<std::vector<std::shared_ptr<TH2>>> fhPt_vs_EtaMatItsTpcTofA{2, {nch, nullptr}};
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
  std::vector<std::shared_ptr<TH2>> fhITS_NCls_vs_PtA{nch, nullptr};
  std::vector<std::shared_ptr<TH2>> fhITS_Chi2NCls_vs_PtA{nch, nullptr};
  std::vector<std::shared_ptr<TH2>> fhTPC_FindableNCls_vs_PtA{nch, nullptr};
  std::vector<std::shared_ptr<TH2>> fhTPC_FoundNCls_vs_PtA{nch, nullptr};
  std::vector<std::shared_ptr<TH2>> fhTPC_SharedNCls_vs_PtA{nch, nullptr};
  std::vector<std::shared_ptr<TH2>> fhTPC_FractionSharedCls_vs_PtA{nch, nullptr};
  std::vector<std::shared_ptr<TH2>> fhTPC_CrossedRows_vs_PtA{nch, nullptr};
  std::vector<std::shared_ptr<TH2>> fhTPC_CrossedRowsOverFindableCls_vs_PtA{nch, nullptr};
  std::vector<std::shared_ptr<TH2>> fhTPC_Chi2NCls_vs_PtA{nch, nullptr};

  template <efficiencyandqatask::KindOfProcess kind>
  void init(HistogramRegistry& registry, const char* dirname)
  {
    using namespace efficiencyandqatask;

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

#define ADDHISTOGRAM(thetype, thedir, therecogen, thewhen, thename, thetitle, thekind, thebinning...) \
  registry.add<thetype>(TString::Format("%s/%s/%s/%s", thedir, therecogen, thewhen, thename).Data(), thetitle, thekind, thebinning)
#define ADDHISTOGRAMID(thetype, thedir, therecogen, thewhen, thename, thetitle, theid, thekind, thebinning...) \
  registry.add<thetype>(TString::Format("%s/%s/%s/%s_%s", thedir, therecogen, thewhen, thename, theid.c_str()).Data(), TString::Format("%s %s", thetitle, theid.c_str()).Data(), thekind, {thebinning})

    /* the reconstructed and generated levels histograms */
    std::string recogen = (kind == kReco) ? "Reco" : "Gen";
    fhPtB[kind] = ADDHISTOGRAM(TH1, dirname, recogen.c_str(), "Before", "Pt", "#it{p}_{T}", kTH1F, {ptAxis});
    fhPt_vs_EtaB[kind] = ADDHISTOGRAM(TH2, dirname, recogen.c_str(), "Before", "PtVsEta", "#it{p}_T vs #eta", kTH2F, {etaAxis, ptAxis});
    fhPt_vs_ZvtxB[kind] = ADDHISTOGRAM(TH2, dirname, recogen.c_str(), "Before", "PtVsZvtx", "#it{p}_T vs #it{z}_{vtx}", kTH2F, {zvtxAxis, ptAxis});
    for (uint isp = 0; isp < nch; ++isp) {
      fhPtA[kind][isp] = ADDHISTOGRAMID(TH1, dirname, recogen.c_str(), "After", "Pt", "#it{p}_{T}", tnames[isp], kTH1F, {ptAxis});
      fhPt_vs_EtaA[kind][isp] = ADDHISTOGRAMID(TH2, dirname, recogen.c_str(), "After", "PtVsEta", "#it{p}_{T} vs #eta", tnames[isp], kTH2F, {etaAxis, ptAxis});
      fhPt_vs_ZvtxA[kind][isp] = ADDHISTOGRAMID(TH2, dirname, recogen.c_str(), "After", "PtVsZvtx", "#it{p}_{T} vs #it{z}_{zvtx}", tnames[isp], kTH2F, {zvtxAxis, ptAxis});
    }

    if constexpr (kind == kReco) {
      /* only the reconstructed level histograms*/
      fhITS_NCls_vs_PtB = ADDHISTOGRAM(TH2, dirname, "Reco", "Before", "ITSNCls", "ITS clusters", kTH2F, {ptAxis, itsNClsAxis});
      fhITS_Chi2NCls_vs_PtB = ADDHISTOGRAM(TH2, dirname, "Reco", "Before", "ITSChi2NCls", "ITS #Chi^{2}", kTH2F, {ptAxis, itsCh2Axis});
      fhTPC_FindableNCls_vs_PtB = ADDHISTOGRAM(TH2, dirname, "Reco", "Before", "TPCFindableNCls", "TPC findable clusters", kTH2F, {ptAxis, tpcNClsAxis});
      fhTPC_FoundNCls_vs_PtB = ADDHISTOGRAM(TH2, dirname, "Reco", "Before", "TPCFoundNCls", "TPC found clusters", kTH2F, {ptAxis, tpcNClsAxis});
      fhTPC_SharedNCls_vs_PtB = ADDHISTOGRAM(TH2, dirname, "Reco", "Before", "TPCSharedNCls", "TPC shared clusters", kTH2F, {ptAxis, tpcNClsAxis});
      fhTPC_FractionSharedCls_vs_PtB = ADDHISTOGRAM(TH2, dirname, "Reco", "Before", "TPCFractionSharedCls", "TPC fraction shared clusters", kTH2F, {ptAxis, tpcFractionAxis});
      fhTPC_CrossedRows_vs_PtB = ADDHISTOGRAM(TH2, dirname, "Reco", "Before", "TPCXrows", "TPC crossed rows", kTH2F, {ptAxis, tpcNRowsAxis});
      fhTPC_CrossedRowsOverFindableCls_vs_PtB = ADDHISTOGRAM(TH2, dirname, "Reco", "Before", "XRowsOverFindableCls", "TPC xrows over findable clusters", kTH2F, {ptAxis, tpcXRowsOverFindClsAxis});
      fhTPC_Chi2NCls_vs_PtB = ADDHISTOGRAM(TH2, dirname, "Reco", "Before", "TPCChi2NCls", "TPC #Chi^{2}", kTH2F, {ptAxis, tpcCh2Axis});
      for (uint isp = 0; isp < nch; ++isp) {
        fhITS_NCls_vs_PtA[isp] = ADDHISTOGRAMID(TH2, dirname, "Reco", "After", "ITSNCls", "ITS clusters", tnames[isp], kTH2F, {ptAxis, itsNClsAxis});
        fhITS_Chi2NCls_vs_PtA[isp] = ADDHISTOGRAMID(TH2, dirname, "Reco", "After", "ITSChi2NCls", "ITS #Chi^{2}", tnames[isp], kTH2F, {ptAxis, itsCh2Axis});
        fhTPC_FindableNCls_vs_PtA[isp] = ADDHISTOGRAMID(TH2, dirname, "Reco", "After", "TPCFindableNCls", "TPC findable clusters", tnames[isp], kTH2F, {ptAxis, tpcNClsAxis});
        fhTPC_FoundNCls_vs_PtA[isp] = ADDHISTOGRAMID(TH2, dirname, "Reco", "After", "TPCFoundNCls", "TPC found clusters", tnames[isp], kTH2F, {ptAxis, tpcNClsAxis});
        fhTPC_SharedNCls_vs_PtA[isp] = ADDHISTOGRAMID(TH2, dirname, "Reco", "After", "TPCSharedNCls", "TPC shared clusters", tnames[isp], kTH2F, {ptAxis, tpcNClsAxis});
        fhTPC_FractionSharedCls_vs_PtA[isp] = ADDHISTOGRAMID(TH2, dirname, "Reco", "After", "TPCFractionSharedCls", "TPC fraction shared clusters", tnames[isp], kTH2F, {ptAxis, tpcFractionAxis});
        fhTPC_CrossedRows_vs_PtA[isp] = ADDHISTOGRAMID(TH2, dirname, "Reco", "After", "TPCXrows", "TPC crossed rows", tnames[isp], kTH2F, {ptAxis, tpcNRowsAxis});
        fhTPC_CrossedRowsOverFindableCls_vs_PtA[isp] = ADDHISTOGRAMID(TH2, dirname, "Reco", "After", "XRowsOverFindableCls", "TPC xrows over findable clusters", tnames[isp], kTH2F, {ptAxis, tpcXRowsOverFindClsAxis});
        fhTPC_Chi2NCls_vs_PtA[isp] = ADDHISTOGRAMID(TH2, dirname, "Reco", "After", "TPCChi2NCls", "TPC #Chi^{2}", tnames[isp], kTH2F, {ptAxis, tpcCh2Axis});
        /* efficiency histograms */
        fhPt_vs_EtaItsA[isp] = ADDHISTOGRAMID(TH2, dirname, "Efficiency", "Reco", "ptIts", "ITS tracks", tnames[isp], kTH2F, {etaAxis, ptAxis});
        fhPt_vs_EtaTpcA[isp] = ADDHISTOGRAMID(TH2, dirname, "Efficiency", "Reco", "ptTpc", "TPC tracks", tnames[isp], kTH2F, {etaAxis, ptAxis});
        fhPt_vs_EtaItsTpcA[isp] = ADDHISTOGRAMID(TH2, dirname, "Efficiency", "Reco", "ptItsTpc", "ITS&TPC tracks", tnames[isp], kTH2F, {etaAxis, ptAxis});
        fhPt_vs_EtaItsTofA[isp] = ADDHISTOGRAMID(TH2, dirname, "Efficiency", "Reco", "ptItsTof", "ITS&TOF tracks", tnames[isp], kTH2F, {etaAxis, ptAxis});
        fhPt_vs_EtaTpcTofA[isp] = ADDHISTOGRAMID(TH2, dirname, "Efficiency", "Reco", "ptTpcTof", "TPC&TOF tracks", tnames[isp], kTH2F, {etaAxis, ptAxis});
        fhPt_vs_EtaItsTpcTofA[isp] = ADDHISTOGRAMID(TH2, dirname, "Efficiency", "Reco", "ptItsTpcTof", "ITS&TPC&TOF tracks", tnames[isp], kTH2F, {etaAxis, ptAxis});
      }
    } else {
      for (uint isp = 0; isp < nch; ++isp) {
        /* detector level and generator level histograms */
        fhPt_vs_EtaPrimA[isp] = ADDHISTOGRAMID(TH2, dirname, "Efficiency", "Gen", "ptPrim", "ITS tracks (primaries)", tnames[isp], kTH2F, {etaAxis, ptAxis});
        fhPt_vs_EtaSecA[isp] = ADDHISTOGRAMID(TH2, dirname, "Efficiency", "Gen", "ptSec", "ITS tracks (secondaries)", tnames[isp], kTH2F, {etaAxis, ptAxis});
        fhPt_vs_EtaMatA[isp] = ADDHISTOGRAMID(TH2, dirname, "Efficiency", "Gen", "ptMat", "ITS tracks (from material)", tnames[isp], kTH2F, {etaAxis, ptAxis});
        const std::vector<std::string> detectedorigin = {"DetReco", "DetAssoc"};
        for (uint ix = 0; ix < detectedorigin.size(); ++ix) {
          fhPt_vs_EtaPrimItsA[ix][isp] = ADDHISTOGRAMID(TH2, dirname, "Efficiency", detectedorigin[ix].c_str(), "ptItsPrim", "ITS tracks (primaries)", tnames[isp], kTH2F, {etaAxis, ptAxis});
          fhPt_vs_EtaPrimItsTpcA[ix][isp] = ADDHISTOGRAMID(TH2, dirname, "Efficiency", detectedorigin[ix].c_str(), "ptItsTpcPrim", "ITS&TPC tracks (primaries)", tnames[isp], kTH2F, {etaAxis, ptAxis});
          fhPt_vs_EtaPrimItsTpcTofA[ix][isp] = ADDHISTOGRAMID(TH2, dirname, "Efficiency", detectedorigin[ix].c_str(), "ptItsTpcTofPrim", "ITS&TPC&TOF tracks (primaries)", tnames[isp], kTH2F, {etaAxis, ptAxis});
          fhPt_vs_EtaSecItsA[ix][isp] = ADDHISTOGRAMID(TH2, dirname, "Efficiency", detectedorigin[ix].c_str(), "ptItsSec", "ITS tracks (secondaries)", tnames[isp], kTH2F, {etaAxis, ptAxis});
          fhPt_vs_EtaSecItsTpcA[ix][isp] = ADDHISTOGRAMID(TH2, dirname, "Efficiency", detectedorigin[ix].c_str(), "ptItsTpcSec", "ITS&TPC tracks (secondaries)", tnames[isp], kTH2F, {etaAxis, ptAxis});
          fhPt_vs_EtaSecItsTpcTofA[ix][isp] = ADDHISTOGRAMID(TH2, dirname, "Efficiency", detectedorigin[ix].c_str(), "ptItsTpcTofSec", "ITS&TPC&TOF tracks (secondaries)", tnames[isp], kTH2F, {etaAxis, ptAxis});
          fhPt_vs_EtaMatItsA[ix][isp] = ADDHISTOGRAMID(TH2, dirname, "Efficiency", detectedorigin[ix].c_str(), "ptItsMat", "ITS tracks (from material)", tnames[isp], kTH2F, {etaAxis, ptAxis});
          fhPt_vs_EtaMatItsTpcA[ix][isp] = ADDHISTOGRAMID(TH2, dirname, "Efficiency", detectedorigin[ix].c_str(), "ptItsTpcMat", "ITS&TPC tracks (from material)", tnames[isp], kTH2F, {etaAxis, ptAxis});
          fhPt_vs_EtaMatItsTpcTofA[ix][isp] = ADDHISTOGRAMID(TH2, dirname, "Efficiency", detectedorigin[ix].c_str(), "ptItsTpcTofMat", "ITS&TPC&TOF tracks (from material)", tnames[isp], kTH2F, {etaAxis, ptAxis});
        }
      }
    }
  }

  template <efficiencyandqatask::KindOfProcess kind, typename TrackListObject>
  void processTracks(float zvtx, TrackListObject const& passedtracks)
  {
    using namespace efficiencyandqatask;

    for (auto& track : passedtracks) {
      fhPtB[kind]->Fill(track.pt());
      fhPt_vs_EtaB[kind]->Fill(track.eta(), track.pt());
      fhPt_vs_ZvtxB[kind]->Fill(zvtx, track.pt());
      if (!(track.trackacceptedid() < 0)) {
        fhPtA[kind][track.trackacceptedid()]->Fill(track.pt());
        fhPt_vs_EtaA[kind][track.trackacceptedid()]->Fill(track.eta(), track.pt());
        fhPt_vs_ZvtxA[kind][track.trackacceptedid()]->Fill(zvtx, track.pt());
      }
      if constexpr (kind == kReco) {
        bool hasits = track.hasITS();
        bool hastpc = track.hasTPC();
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
          auto fillhisto = [&track](auto& h, bool cond) {
            if (cond) {
              h->Fill(track.eta(), track.pt());
            }
          };
          fillhisto(fhPt_vs_EtaItsA[track.trackacceptedid()], hasits);
          fillhisto(fhPt_vs_EtaTpcA[track.trackacceptedid()], hastpc);
          fillhisto(fhPt_vs_EtaItsTpcA[track.trackacceptedid()], hasits && hastpc);
          fillhisto(fhPt_vs_EtaItsTofA[track.trackacceptedid()], hasits && hastof);
          fillhisto(fhPt_vs_EtaTpcTofA[track.trackacceptedid()], hastpc && hastof);
          fillhisto(fhPt_vs_EtaItsTpcTofA[track.trackacceptedid()], hasits && hastpc && hastof);
          /* the detector / generator combined level */
          if constexpr (framework::has_type_v<aod::mctracklabel::McParticleId, typename TrackListObject::iterator::all_columns>) {
            /* get the associated MC particle we are sure it does exist because the track was accepted */
            const auto& mcparticle = track.template mcParticle_as<soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo>>();

            /* TODO: what if the id of the generated is not the same as the id of the reconstructed */
            bool isprimary = mcparticle.isPhysicalPrimary();
            bool issecdecay = !isprimary && (mcparticle.getProcess() == 4);
            bool isfrommaterial = !isprimary && !issecdecay;
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
  }
};

struct EfficiencyCollectingEngine {
  size_t nch = efficiencyandqatask::tnames.size();

  std::vector<TH1F*> fhN1_pt{nch, nullptr};
};

struct DptDptEfficiencyAndQc {
  /* the data memebers for this task */
  /* the centrality / multiplicity limits for collecting data in this task instance */
  uint ncmranges = 0;
  float* fCentMultMin = nullptr;
  float* fCentMultMax = nullptr;

  /* the data collecting engine instances */
  QADataCollectingEngine** dataCE;

  /* this is a nightmare but not other way found to overcome the histogram registry limit */
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
  std::vector<HistogramRegistry*> registrybank{&registry_one, &registry_two, &registry_three, &registry_four, &registry_five,
                                               &registry_six, &registry_seven, &registry_eight, &registry_nine, &registry_ten};

  void init(o2::framework::InitContext& initContext)
  {
    using namespace efficiencyandqatask;
    /* Self configuration: requires dptdptfilter task in the workflow */
    {
      /* the binning */
      getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mZVtxbins", zvtxbins);
      getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mZVtxmin", zvtxlow);
      getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mZVtxmax", zvtxup);
      getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mPTbins", ptbins);
      getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mPTmin", ptlow);
      getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mPTmax", ptup);
      getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mEtabins", etabins);
      getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mEtamin", etalow);
      getTaskOptionValue(initContext, "dpt-dpt-filter", "binning.mEtamax", etaup);

      /* configuring the involved species */
      o2::analysis::dptdptfilter::PIDSpeciesSelection pidselector;
      std::vector<std::string> cfgnames = {"pipidsel", "kapidsel", "prpidsel"};
      std::vector<uint8_t> spids = {2, 3, 4};
      for (uint i = 0; i < cfgnames.size(); ++i) {
        auto includeIt = [&pidselector, &initContext](int spid, auto name) {
          bool mUseIt = false;
          bool mExcludeIt = false;
          if (getTaskOptionValue(initContext, "dpt-dpt-filter-tracks", TString::Format("%s.mUseIt", name.c_str()).Data(), mUseIt) &&
              getTaskOptionValue(initContext, "dpt-dpt-filter-tracks", TString::Format("%s.mExclude", name.c_str()).Data(), mExcludeIt)) {
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
      if (getTaskOptionValue(initContext, "dpt-dpt-filter", "centralities", centspec)) {
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
      dataCE = new QADataCollectingEngine*[ncmranges];
      std::string recogen;
      if (!doprocessReconstructedNotStored && !doprocessDetectorLevelNotStored) {
        LOGF(fatal, "Neither reco nor detector level not configured. Please, fix it!");
      }
      if (ncmranges > registrybank.size()) {
        LOGF(fatal, "There are more centrality ranges configured than registries in the bank. Please fix it!");
      }
      /* in reverse order for proper order in results file */
      for (uint i = 0; i < ncmranges; ++i) {
        auto initializeCEInstance = [&](auto dce, auto name) {
          /* crete the output list for the passed centrality/multiplicity range */
          /* init the data collection instance */
          dce->template init<kReco>(*registrybank[i], name.Data());
          if (doprocessGeneratorLevelNotStored) {
            dce->template init<kGen>(*registrybank[i], name.Data());
          }
        };
        auto buildCEInstance = [&initializeCEInstance](float min, float max) {
          QADataCollectingEngine* dce = new QADataCollectingEngine();
          initializeCEInstance(dce, TString::Format("EfficiencyAndQaData-%d-%d", int(min), int(max)));
          return dce;
        };
        /* in reverse order for proper order in results file */
        dataCE[ncmranges - i - 1] = buildCEInstance(fCentMultMin[ncmranges - i - 1], fCentMultMax[ncmranges - i - 1]);
      }
      for (uint i = 0; i < ncmranges; ++i) {
        LOGF(info, " centrality/multipliicty range: %d, low limit: %d, up limit: %d", i, int(fCentMultMin[i]), int(fCentMultMax[i]));
      }
    }
  }

  /// \brief Get the data collecting engine index corresponding to the passed collision
  template <typename FilteredCollision>
  int getDCEindex(FilteredCollision collision)
  {
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

  template <efficiencyandqatask::KindOfProcess kind, typename FilterdCollision, typename PassedTracks>
  void processTracks(FilterdCollision const& collision, PassedTracks const& tracks)
  {
    using namespace efficiencyandqatask;

    int ixDCE = getDCEindex(collision);
    if (!(ixDCE < 0)) {
      dataCE[ixDCE]->processTracks<kind>(collision.posZ(), tracks);
    }
  }

  Filter onlyacceptedcollisions = (aod::dptdptfilter::collisionaccepted == uint8_t(true));

  void processReconstructedNotStored(soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>::iterator const& collision, soa::Join<aod::FullTracks, aod::DptDptCFTracksInfo>& tracks)
  {
    using namespace efficiencyandqatask;

    processTracks<kReco>(collision, tracks);
  }
  PROCESS_SWITCH(DptDptEfficiencyAndQc, processReconstructedNotStored, "Process reconstructed efficiency and QA for not stored derived data", false);

  void processDetectorLevelNotStored(soa::Filtered<soa::Join<aod::Collisions, aod::DptDptCFCollisionsInfo>>::iterator const& collision,
                                     soa::Join<aod::FullTracks, aod::DptDptCFTracksInfo, aod::McTrackLabels>& tracks,
                                     soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo> const&)
  {
    using namespace efficiencyandqatask;

    processTracks<kReco>(collision, tracks);
  }
  PROCESS_SWITCH(DptDptEfficiencyAndQc, processDetectorLevelNotStored, "Process MC detector level efficiency and QA for not stored derived data", true);

  void processGeneratorLevelNotStored(soa::Filtered<soa::Join<aod::McCollisions, aod::DptDptCFGenCollisionsInfo>>::iterator const& collision,
                                      soa::Join<aod::McParticles, aod::DptDptCFGenTracksInfo>& particles)
  {
    using namespace efficiencyandqatask;

    processTracks<kGen>(collision, particles);
  }
  PROCESS_SWITCH(DptDptEfficiencyAndQc, processGeneratorLevelNotStored, "Process MC generator level efficiency and QA for not stored derived data", true);
};

//****************************************************************************************
/**
 * Workflow definition.
 */
//****************************************************************************************
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<DptDptEfficiencyAndQc>(cfgc)};
  return workflow;
}