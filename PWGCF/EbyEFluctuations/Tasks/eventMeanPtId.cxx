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
//
/// \file eventMeanPtId.cxx
/// \brief Analysis task to study Mean pT Fluctuations using two particle correlator using Cumulant Method
/// \author Sweta Singh (sweta.singh@cern.ch)

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"

#include <TPDGCode.h>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

double massPi = o2::constants::physics::MassPionCharged;
double massKa = o2::constants::physics::MassKaonCharged;
double massPr = o2::constants::physics::MassProton;

using namespace o2::constants::physics;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;
using o2::constants::physics::Pdg;

namespace o2::aod
{
using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFT0Ms>;
using MyTracks = soa::Join<aod::FullTracks,
                           aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                           aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
                           aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::StoredTracks,
                           aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFbeta, aod::TOFSignal, aod::TracksExtra, aod::TracksIU, aod::TracksDCA, aod::TrackSelection>;

using MyMCRecoCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::McCollisionLabels>;
using MyMCRecoCollision = MyMCRecoCollisions::iterator;

using MyMCRecoTracks = soa::Join<aod::FullTracks,
                                 aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                 aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
                                 aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::StoredTracks,
                                 aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFbeta, aod::TOFSignal, aod::TracksExtra, aod::TracksIU, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>;
using MyMCRecoTrack = MyMCRecoTracks::iterator;

using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Cs, aod::CentFT0Ms, aod::Mults>;
using MyCollision = MyCollisions::iterator;
using MyTrack = MyTracks::iterator;
} // namespace o2::aod

struct EventMeanPtId {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;

  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> cfgUrlCCDB{"cfgUrlCCDB", "http://alice-ccdb.cern.ch", "url of ccdb"};
  Configurable<std::string> cfgPathCCDB{"cfgPathCCDB", "Users/s/swsingh/My/Object/eff_Pb", "Path for ccdb-object"};
  Configurable<bool> cfgLoadEff{"cfgLoadEff", true, "Load efficiency"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  TH1D* ptHistogramAllchargeRec = nullptr;
  TH1D* ptHistogramPionrec = nullptr;
  TH1D* ptHistogramKaonrec = nullptr;
  TH1D* ptHistogramProtonrec = nullptr;
  TH1D* hRecoPi = nullptr;
  TH1D* hRecoKa = nullptr;
  TH1D* hRecoPr = nullptr;
  TH2D* hPtyPion = nullptr;
  TH2D* hPtyKaon = nullptr;
  TH2D* hPtyProton = nullptr;

  Configurable<float> ptMax{"ptMax", 2.0, "maximum pT"};
  Configurable<float> ptMin{"ptMin", 0.15, "minimum pT"};
  Configurable<std::vector<double>> ptBins{"ptBins", {0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00}, "p_{T} bins"};
  Configurable<float> piluprejection{"piluprejection", false, "Pileup rejection"};

  void init(o2::framework::InitContext&)
  {
    if (cfgLoadEff) {
      // Set CCDB url
      ccdb->setURL(cfgUrlCCDB.value);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      // ccdb->setCreatedNotAfter(ccdbNoLaterThan.value);
      // LOGF(info, "Getting object %s", ccdbPath.value.data());

      TList* lst = ccdb->getForTimeStamp<TList>(cfgPathCCDB.value, -1);
      ptHistogramAllchargeRec = reinterpret_cast<TH1D*>(lst->FindObject("ptHistogramAllchargeRec"));
      ptHistogramPionrec = reinterpret_cast<TH1D*>(lst->FindObject("ptHistogramPionrec"));
      ptHistogramKaonrec = reinterpret_cast<TH1D*>(lst->FindObject("ptHistogramKaonrec"));
      ptHistogramProtonrec = reinterpret_cast<TH1D*>(lst->FindObject("ptHistogramProtonrec"));
      hRecoPi = reinterpret_cast<TH1D*>(lst->FindObject("hRecoPi"));
      hRecoKa = reinterpret_cast<TH1D*>(lst->FindObject("hRecoKa"));
      hRecoPr = reinterpret_cast<TH1D*>(lst->FindObject("hRecoPr"));
      hPtyPion = reinterpret_cast<TH2D*>(lst->FindObject("hPtyPion"));
      hPtyKaon = reinterpret_cast<TH2D*>(lst->FindObject("hPtyKaon"));
      hPtyProton = reinterpret_cast<TH2D*>(lst->FindObject("hPtyProton"));

      if (!ptHistogramAllchargeRec || !ptHistogramPionrec || !ptHistogramKaonrec || !ptHistogramProtonrec || !hRecoPi || !hRecoKa || !hRecoPr || !hPtyPion || !hPtyKaon || !hPtyProton) {
        LOGF(info, "FATAL!! Could not find required histograms in CCDB");
      }
    }

    std::vector<double> ptBinning = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0};
    //  AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec vtxZAxis = {100, -20.0, 20.0, "Z (cm)"};
    AxisSpec dcaAxis = {1002, -5.01, 5.01, "DCA_{xy} (cm)"};
    AxisSpec dcazAxis = {1002, -5.01, 5.01, "DCA_{z} (cm)"};
    AxisSpec ptAxis = {600, 0.0, 6.0, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pAxis = {400, 0.0, 4.0, "#it{p} (GeV/#it{c})"};
    AxisSpec betaAxis = {200, 0.0, 2.0, "TOF_{#beta} (GeV/#it{c})"};
    AxisSpec dEdxAxis = {2000, 0.0, 200.0, "dE/dx (GeV/#it{c})"};
    AxisSpec etaAxis = {300, -1.5, 1.5, "#eta"}; // 300, -1.5, 1.5
    AxisSpec nSigmaTPCAxis = {170, -8.5, 8.5, "n#sigma_{TPC}^{proton}"};
    AxisSpec nSigmaTPCAxispid = {170, -8.5, 8.5, "n#sigma_{TPC}"};
    AxisSpec nSigmaTOFAxispid = {170, -8.5, 8.5, "n#sigma_{TOF}"};
    AxisSpec centAxis = {100, 0., 100., "centrality"};
    AxisSpec subAxis = {30, 0., 30., "sample"};
    AxisSpec nchAxis = {4000, 0., 4000., "nch"};
    AxisSpec varAxis1 = {400, 0., 4., "var1"};
    AxisSpec varAxis2 = {400, 0., 4., "var2"};
    AxisSpec chi2Axis = {100, 0., 100., "Chi2"};
    AxisSpec crossedRowTpcAxis = {600, 0., 600., "TPC Crossed rows"};
    AxisSpec counter = {10, 0., 10., "events"};

    // QA Plots
    histos.add("hEventcounter", "event counts", kTH1D, {counter});
    auto h = histos.add<TH1>("tracksel", "tracksel", HistType::kTH1D, {{10, 0.5, 10.5}});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, "Global track passed");
    h->GetXaxis()->SetBinLabel(3, "DCAxy passed");
    h->GetXaxis()->SetBinLabel(4, "DCAz passed");
    h->GetXaxis()->SetBinLabel(5, "Eta-cut passed");
    h->GetXaxis()->SetBinLabel(6, "pT-cut passed");
    h->GetXaxis()->SetBinLabel(7, "TPC crossed rows passed");
    h->GetXaxis()->SetBinLabel(8, "TPC Chai2cluster passed");
    h->GetXaxis()->SetBinLabel(9, "ITS Chai2cluster passed");

    histos.add("hEventcounter_recMC", "event counts rec MC", kTH1D, {counter});
    auto hRec = histos.add<TH1>("trackSelRec", "trackSelRec", HistType::kTH1D, {{10, 0.5, 10.5}});
    hRec->GetXaxis()->SetBinLabel(1, "has_mcCollision() read");
    hRec->GetXaxis()->SetBinLabel(2, "Vertex Z > 10cm passed");
    hRec->GetXaxis()->SetBinLabel(3, "sel 8 passed");
    hRec->GetXaxis()->SetBinLabel(4, "kNoSameBunchPileup passed");
    hRec->GetXaxis()->SetBinLabel(5, "kNoITSROFrameBorder passed");
    hRec->GetXaxis()->SetBinLabel(6, "klsGoodZvtxFT0vsPV passed");
    hRec->GetXaxis()->SetBinLabel(7, "klsVertexITSTPC passed");

    histos.add("Data/hZvtx_before_sel", "hZvtx_before_sel", kTH1D, {vtxZAxis});
    histos.add("Data/hZvtx_after_sel", "hZvtx_after_sel", kTH1D, {vtxZAxis});
    histos.add("Data/hZvtx_after_sel8", "hZvtx_after_sel8", kTH1D, {vtxZAxis});
    histos.add("Data/hP", "hP", kTH1D, {pAxis});
    histos.add("Data/hEta", ";hEta", kTH1D, {etaAxis});
    histos.add("Data/hPt", ";#it{p}_{T} (GeV/#it{c})", kTH1D, {ptAxis});
    histos.add("Data/hNsigmaTPC", "hNsigmaTPC", kTH2D, {pAxis, nSigmaTPCAxis});
    histos.add("Data/hDCAxy", "hDCAxy", kTH1D, {dcaAxis});
    histos.add("Data/hDCAz", "hDCAz", kTH1D, {dcazAxis});
    histos.add("Data/hPtDCAxy", "hPtDCAxy", kTH2D, {ptAxis, dcaAxis});
    histos.add("Data/hPtDCAz", "hPtDCAz", kTH2D, {ptAxis, dcazAxis});
    histos.add("Data/NSigamaTPCpion", "NSigamaTPCpion", kTH2D, {ptAxis, nSigmaTPCAxispid});
    histos.add("Data/NSigamaTPCkaon", "NSigamaTPCkaon", kTH2D, {ptAxis, nSigmaTPCAxispid});
    histos.add("Data/NSigamaTPCproton", "NSigamaTPCproton", kTH2D, {ptAxis, nSigmaTPCAxispid});
    histos.add("Data/NSigamaTOFpion", "NSigamaTOFpion", kTH2D, {ptAxis, nSigmaTOFAxispid});
    histos.add("Data/NSigamaTOFkaon", "NSigamaTOFkaon", kTH2D, {ptAxis, nSigmaTOFAxispid});
    histos.add("Data/NSigamaTOFproton", "NSigamaTOFproton", kTH2D, {ptAxis, nSigmaTOFAxispid});
    histos.add("Data/NSigamaTPCTOFpion", "NSigamaTPCTOFpion", kTH2D, {nSigmaTPCAxispid, nSigmaTOFAxispid});
    histos.add("Data/NSigamaTPCTOFkaon", "NSigamaTPCTOFkaon", kTH2D, {nSigmaTPCAxispid, nSigmaTOFAxispid});
    histos.add("Data/NSigamaTPCTOFproton", "NSigamaTPCTOFproton", kTH2D, {nSigmaTPCAxispid, nSigmaTOFAxispid});
    histos.add("Data/hPtPion", ";#it{p}_{T} (GeV/#it{c})", kTH1D, {ptAxis});
    histos.add("Data/hPtKaon", ";#it{p}_{T} (GeV/#it{c})", kTH1D, {ptAxis});
    histos.add("Data/hPtProton", ";#it{p}_{T} (GeV/#it{c})", kTH1D, {ptAxis});
    histos.add("Data/hEtaPion", ";hEta", kTH1D, {etaAxis});
    histos.add("Data/hEtaKaon", ";hEta", kTH1D, {etaAxis});
    histos.add("Data/hEtaProton", ";hEta", kTH1D, {etaAxis});
    histos.add("Data/hyPion", ";hyPion", kTH1D, {etaAxis});
    histos.add("Data/hyKaon", ";hyKaon", kTH1D, {etaAxis});
    histos.add("Data/hyProton", ";hyProton", kTH1D, {etaAxis});
    histos.add("Data/hPtCh", "hPtCh", kTH2D, {nchAxis, ptAxis});
    histos.add("Data/hPtChPion", "hPtChPion", kTH2D, {nchAxis, ptAxis});
    histos.add("Data/hPtChKaon", "hPtChKaon", kTH2D, {nchAxis, ptAxis});
    histos.add("Data/hPtChProton", "hPtChProton", kTH2D, {nchAxis, ptAxis});
    histos.add("Data/hPtCent", "hPtCent", kTH2D, {centAxis, ptAxis});
    histos.add("Data/hPtCentPion", "hPtCentPion", kTH2D, {centAxis, ptAxis});
    histos.add("Data/hPtCentKaon", "hPtCentKaon", kTH2D, {centAxis, ptAxis});
    histos.add("Data/hPtCentProton", "hPtCentProton", kTH2D, {centAxis, ptAxis});
    histos.add("Data/hMeanPtCh", "hMeanPtCh", kTH2D, {nchAxis, ptAxis});
    histos.add("Data/hCent", "hCent", kTH2D, {nchAxis, centAxis});

    histos.add("Data/hVar1", "hVar1", kTH2D, {subAxis, centAxis});
    histos.add("Data/hVar2", "hVar2", kTH2D, {subAxis, centAxis});
    histos.add("Data/hVar2meanpt", "hVar2meanpt", kTH2D, {centAxis, varAxis2});
    histos.add("Data/hVar", "hVar", kTH2D, {subAxis, centAxis});
    histos.add("Data/hVarc", "hVarc", kTH2D, {subAxis, centAxis});
    histos.add("Data/hVar1pi", "hVar1pi", kTH2D, {subAxis, centAxis});
    histos.add("Data/hVar2pi", "hVar2pi", kTH2D, {subAxis, centAxis});
    histos.add("Data/hVarpi", "hVarpi", kTH2D, {subAxis, centAxis});
    histos.add("Data/hVar2meanptpi", "hVar2meanptpi", kTH2D, {centAxis, varAxis2});
    histos.add("Data/hVar1k", "hVar1k", kTH2D, {subAxis, centAxis});
    histos.add("Data/hVar2k", "hVar2k", kTH2D, {subAxis, centAxis});
    histos.add("Data/hVark", "hVark", kTH2D, {subAxis, centAxis});
    histos.add("Data/hVar2meanptk", "hVar2meanptk", kTH2D, {centAxis, varAxis2});
    histos.add("Data/hVar1p", "hVar1p", kTH2D, {subAxis, centAxis});
    histos.add("Data/hVar2p", "hVar2p", kTH2D, {subAxis, centAxis});
    histos.add("Data/hVarp", "hVarp", kTH2D, {subAxis, centAxis});
    histos.add("Data/hVar2meanptp", "hVar2meanptp", kTH2D, {centAxis, varAxis2});

    histos.add("Data/hnchAll", ";hnchAll", kTH1D, {nchAxis});
    histos.add("Data/hnchAll_bf_cut", ";hnchAll_bf_cut", kTH1D, {nchAxis});
    histos.add("Data/hnch", ";hnch", kTH1D, {nchAxis});
    histos.add("Data/hnchTrue", ";hnchTrue", kTH1D, {nchAxis});
    histos.add("Data/hnchTrue_pt", ";hnchTrue_pt", kTH1D, {nchAxis});

    histos.add("Data/hVar1x", "hVar1x", kTH2D, {subAxis, nchAxis});
    histos.add("Data/hVar2x", "hVar2x", kTH2D, {subAxis, nchAxis});
    histos.add("Data/hVarx", "hVarx", kTH2D, {subAxis, nchAxis});
    histos.add("Data/hVar2meanptx", "hVar2meanptx", kTH2D, {nchAxis, varAxis2});
    histos.add("Data/hVar1pix", "hVar1pix", kTH2D, {subAxis, nchAxis});
    histos.add("Data/hVar2pix", "hVar2pix", kTH2D, {subAxis, nchAxis});
    histos.add("Data/hVarpix", "hVarpix", kTH2D, {subAxis, nchAxis});
    histos.add("Data/hVar2meanptpix", "hVar2meanptpix", kTH2D, {nchAxis, varAxis2});
    histos.add("Data/hVar1kx", "hVar1kx", kTH2D, {subAxis, nchAxis});
    histos.add("Data/hVar2kx", "hVar2kx", kTH2D, {subAxis, nchAxis});
    histos.add("Data/hVarkx", "hVarkx", kTH2D, {subAxis, nchAxis});
    histos.add("Data/hVar2meanptkx", "hVar2meanptkx", kTH2D, {nchAxis, varAxis2});
    histos.add("Data/hVar1px", "hVar1px", kTH2D, {subAxis, nchAxis});
    histos.add("Data/hVar2px", "hVar2px", kTH2D, {subAxis, nchAxis});
    histos.add("Data/hVarpx", "hVarpx", kTH2D, {subAxis, nchAxis});
    histos.add("Data/hVar2meanptpx", "hVar2meanptpx", kTH2D, {nchAxis, varAxis2});
    histos.add("Data/ht", "ht", kTH1D, {centAxis});
    histos.add("Data/hCentrality", "hCentrality", kTH1D, {centAxis});
    histos.add("Data/hPEta", "hPEta", kTH2D, {pAxis, etaAxis});
    histos.add("Data/hPtEta", "hPtEta", kTH2D, {ptAxis, etaAxis});
    histos.add("Data/hPy", "hPy", kTH2D, {pAxis, etaAxis});
    histos.add("Data/hPty", "hPty", kTH2D, {ptAxis, etaAxis});
    histos.add("Data/hPtyPion", "hPtyPion", kTH2D, {ptAxis, etaAxis});
    histos.add("Data/hPtyKaon", "hPtyKaon", kTH2D, {ptAxis, etaAxis});
    histos.add("Data/hPtyProton", "hPtyProton", kTH2D, {ptAxis, etaAxis});
    histos.add("Data/hTOFbeta", "hTOFbeta", kTH2D, {pAxis, betaAxis});
    histos.add("Data/hdEdx", "hdEdx", kTH2D, {pAxis, dEdxAxis});
    histos.add("Data/hTOFbeta_afterselection", "hTOFbeta_afterselection", kTH2D, {pAxis, betaAxis});
    histos.add("Data/hdEdx_afterselection", "hdEdx_afterselection", kTH2D, {pAxis, dEdxAxis});
    histos.add("Data/hTOFbeta_afterselection1", "hTOFbeta_afterselection1", kTH2D, {pAxis, betaAxis});
    histos.add("Data/hdEdx_afterselection1", "hdEdx_afterselection1", kTH2D, {pAxis, dEdxAxis});
    histos.add("Data/hTPCchi2perCluster_before", "TPC #Chi^{2}/Cluster", kTH1D, {chi2Axis});
    histos.add("Data/hITSchi2perCluster_before", "ITS #Chi^{2}/Cluster", kTH1D, {chi2Axis});
    histos.add("Data/hTPCCrossedrows_before", "Crossed TPC rows", kTH1D, {crossedRowTpcAxis});
    histos.add("Data/hTPCchi2perCluster_after", "TPC #Chi^{2}/Cluster", kTH1D, {chi2Axis});
    histos.add("Data/hITSchi2perCluster_after", "ITS #Chi^{2}/Cluster", kTH1D, {chi2Axis});
    histos.add("Data/hTPCCrossedrows_after", "Crossed TPC rows", kTH1D, {crossedRowTpcAxis});
    histos.add("Data/hdEdx_rec_bf_anycut", "hdEdx_rec_bf_anycut", kTH2D, {pAxis, dEdxAxis});
    histos.add("Data/hcent_nacc", "hcent_nacc", kTH2D, {centAxis, nchAxis});

    histos.addClone("Data/", "Rec/");
    // rec histograms
    histos.add("NSigamaTPCpion_rec", "NSigamaTPCpion_rec", kTH2D, {pAxis, nSigmaTPCAxispid});
    histos.add("NSigamaTPCkaon_rec", "NSigamaTPCkaon_rec", kTH2D, {pAxis, nSigmaTPCAxispid});
    histos.add("NSigamaTPCproton_rec", "NSigamaTPCproton_rec", kTH2D, {pAxis, nSigmaTPCAxispid});
    histos.add("NSigamaTOFpion_rec", "NSigamaTOFpion_rec", kTH2D, {pAxis, nSigmaTOFAxispid});
    histos.add("NSigamaTOFkaon_rec", "NSigamaTOFkaon_rec", kTH2D, {pAxis, nSigmaTOFAxispid});
    histos.add("NSigamaTOFproton_rec", "NSigamaTOFproton_rec", kTH2D, {pAxis, nSigmaTOFAxispid});
    histos.add("NSigamaTPCTOFpion_rec", "NSigamaTPCTOFpion_rec", kTH2D, {nSigmaTPCAxispid, nSigmaTOFAxispid});
    histos.add("NSigamaTPCTOFkaon_rec", "NSigamaTPCTOFkaon_rec", kTH2D, {nSigmaTPCAxispid, nSigmaTOFAxispid});
    histos.add("NSigamaTPCTOFproton_rec", "NSigamaTPCTOFproton_rec", kTH2D, {nSigmaTPCAxispid, nSigmaTOFAxispid});
    histos.add("NSigamaTPCpion_rec_bf_sel", "NSigamaTPCpion_rec_bf_sel", kTH2D, {pAxis, nSigmaTPCAxispid});
    histos.add("NSigamaTPCkaon_rec_bf_sel", "NSigamaTPCkaon_rec_bf_sel", kTH2D, {pAxis, nSigmaTPCAxispid});
    histos.add("NSigamaTPCproton_rec_bf_sel", "NSigamaTPCproton_rec_bf_sel", kTH2D, {pAxis, nSigmaTPCAxispid});
    histos.add("NSigamaTOFpion_rec_bf_sel", "NSigamaTOFpion_rec_bf_sel", kTH2D, {pAxis, nSigmaTOFAxispid});
    histos.add("NSigamaTOFkaon_rec_bf_sel", "NSigamaTOFkaon_rec_bf_sel", kTH2D, {pAxis, nSigmaTOFAxispid});
    histos.add("NSigamaTOFproton_rec_bf_sel", "NSigamaTOFproton_rec_bf_sel", kTH2D, {pAxis, nSigmaTOFAxispid});
    histos.add("NSigamaTPCTOFpion_rec_bf_sel", "NSigamaTPCTOFpion_rec_bf_sel", kTH2D, {nSigmaTPCAxispid, nSigmaTOFAxispid});
    histos.add("NSigamaTPCTOFkaon_rec_bf_sel", "NSigamaTPCTOFkaon_rec_bf_sel", kTH2D, {nSigmaTPCAxispid, nSigmaTOFAxispid});
    histos.add("NSigamaTPCTOFproton_rec_bf_sel", "NSigamaTPCTOFproton_rec_bf_sel", kTH2D, {nSigmaTPCAxispid, nSigmaTOFAxispid});
    histos.add("hPtyPion_rec", "hPtyPion_rec", kTH2D, {ptAxis, etaAxis});
    histos.add("hPtyKaon_rec", "hPtyKaon_rec", kTH2D, {ptAxis, etaAxis});
    histos.add("hPtyProton_rec", "hPtyProton_rec", kTH2D, {ptAxis, etaAxis});
    histos.add("hPyPion_rec", "hPyPion_rec", kTH2D, {pAxis, etaAxis});
    histos.add("hPyKaon_rec", "hPyKaon_rec", kTH2D, {pAxis, etaAxis});
    histos.add("hPyProton_rec", "hPyProton_rec", kTH2D, {pAxis, etaAxis});
    histos.add("hTOFbeta_afterselection_rec_afterpidcut", "hTOFbeta_afterselection_rec_afterpidcut", kTH2D, {pAxis, betaAxis});
    histos.add("hdEdx_afterselection_rec_afterpidcut", "hdEdx_afterselection_rec_afterpidcut", kTH2D, {pAxis, dEdxAxis});
    histos.add("hTOFbeta_afterselection_rec_beforepidcut", "hTOFbeta_afterselection_rec_beforepidcut", kTH2D, {pAxis, betaAxis});
    histos.add("hdEdx_afterselection_rec_beforepidcut", "hdEdx_afterselection_rec_beforepidcut", kTH2D, {pAxis, dEdxAxis});

    histos.add("heffVar1x", "heffVar1x", kTH2D, {subAxis, nchAxis});
    histos.add("heffVar2x", "heffVar2x", kTH2D, {subAxis, nchAxis});
    histos.add("heffVarx", "heffVarx", kTH2D, {subAxis, nchAxis});
    histos.add("heffVar2meanptx", "heffVar2meanptx", kTH2D, {nchAxis, varAxis2});
    histos.add("hnchRec_all", ";hnchRec_all", kTH1D, {nchAxis});
    histos.add("hnchRec", ";hnchRec", kTH1D, {nchAxis});
    histos.add("hnchRec_true", ";hnchRec_true", kTH1D, {nchAxis});
    histos.add("hVar1x_rec_old", "hVar1x_rec_old", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2x_rec_old", "hVar2x_rec_old", kTH2D, {subAxis, nchAxis});
    histos.add("hVarx_rec_old", "hVarx_rec_old", kTH2D, {subAxis, nchAxis});
    histos.add("hVar1x_rec", "hVar1x_rec", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2x_rec", "hVar2x_rec", kTH2D, {subAxis, nchAxis});
    histos.add("hVarx_rec", "hVarx_rec", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2meanptx_rec", "hVar2meanptx_rec", kTH2D, {nchAxis, varAxis2});
    histos.add("hVar1pix_rec", "hVar1pix_rec", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2pix_rec", "hVar2pix_rec", kTH2D, {subAxis, nchAxis});
    histos.add("hVarpix_rec", "hVarpix_rec", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2meanptpix_rec", "hVar2meanptpix_rec", kTH2D, {nchAxis, varAxis2});
    histos.add("hVar1kx_rec", "hVar1kx_rec", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2kx_rec", "hVar2kx_rec", kTH2D, {subAxis, nchAxis});
    histos.add("hVarkx_rec", "hVarkx_rec", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2meanptkx_rec", "hVar2meanptkx_rec", kTH2D, {nchAxis, varAxis2});
    histos.add("hVar1px_rec", "hVar1px_rec", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2px_rec", "hVar2px_rec", kTH2D, {subAxis, nchAxis});
    histos.add("hVarpx_rec", "hVarpx_rec", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2meanptpx_rec", "hVar2meanptpx_rec", kTH2D, {nchAxis, varAxis2});
    histos.add("hZvtx_after_sel_rec", "hZvtx_after_sel_rec", kTH1D, {vtxZAxis});
    histos.add("hZvtx_after_sel8_rec", "hZvtx_after_sel8_rec", kTH1D, {vtxZAxis});
    histos.add("etaHistogram_allcharge_rec", "etaHistogram_allcharge_rec", kTH1D, {etaAxis});
    histos.add("ptHistogram_allcharge_bfptcut_rec", "ptHistogram_allcharge_bfptcut_rec", kTH1D, {ptAxis});
    histos.add("ptHistogramAllchargeRec", "ptHistogramAllchargeRec", kTH1D, {ptAxis});
    histos.add("ptHistogramPionrec", "ptHistogramPionrec", kTH1D, {ptAxis});
    histos.add("ptHistogramKaonrec", "ptHistogramKaonrec", kTH1D, {ptAxis});
    histos.add("ptHistogramProtonrec", "ptHistogramProtonrec", kTH1D, {ptAxis});
    histos.add("ptHistogramPionrec_purity", "ptHistogramPionrec_purity", kTH1D, {ptAxis});
    histos.add("ptHistogramKaonrec_purity", "ptHistogramKaonrec_purity", kTH1D, {ptAxis});
    histos.add("ptHistogramProtonrec_purity", "ptHistogramProtonrec_purity", kTH1D, {ptAxis});
    histos.add("ptHistogramPionrec_pdg", "ptHistogramPionrec_pdg", kTH1D, {ptAxis});
    histos.add("ptHistogramKaonrec_pdg", "ptHistogramKaonrec_pdg", kTH1D, {ptAxis});
    histos.add("ptHistogramProtonrec_pdg", "ptHistogramProtonrec_pdg", kTH1D, {ptAxis});
    histos.add("Histogram_mass2_p_rec_beforesel", "Histogram_mass2_p_rec_beforesel", kTH1D, {ptAxis});
    histos.add("Histogram_mass2_p_rec_aftersel", "Histogram_mass2_p_rec_aftersel", kTH1D, {ptAxis});
    histos.add("hEffVar1x", "hEffVar1x", kTH2D, {subAxis, nchAxis});
    histos.add("hEffVar2x", "hEffVar2x", kTH2D, {subAxis, nchAxis});
    histos.add("hEffVarx", "hEffVarx", kTH2D, {subAxis, nchAxis});
    histos.add("hEffVar1pix", "hEffVar1pix", kTH2D, {subAxis, nchAxis});
    histos.add("hEffVar2pix", "hEffVar2pix", kTH2D, {subAxis, nchAxis});
    histos.add("hEffVarpix", "hEffVarpix", kTH2D, {subAxis, nchAxis});
    histos.add("hEffVar1kx", "hEffVar1kx", kTH2D, {subAxis, nchAxis});
    histos.add("hEffVar2kx", "hEffVar2kx", kTH2D, {subAxis, nchAxis});
    histos.add("hEffVarkx", "hEffVarkx", kTH2D, {subAxis, nchAxis});
    histos.add("hEffVar1px", "hEffVar1px", kTH2D, {subAxis, nchAxis});
    histos.add("hEffVar2px", "hEffVar2px", kTH2D, {subAxis, nchAxis});
    histos.add("hEffVarpx", "hEffVarpx", kTH2D, {subAxis, nchAxis});
    histos.add("hEffVar2Meanptx", "hEffVar2Meanptx", kTH2D, {nchAxis, varAxis2});
    histos.add("hEffVar2Meanptpix", "hEffVar2Meanptpix", kTH2D, {nchAxis, varAxis2});
    histos.add("hEffVar2Meanptkx", "hEffVar2Meanptkx", kTH2D, {nchAxis, varAxis2});
    histos.add("hEffVar2Meanptpx", "hEffVar2Meanptpx", kTH2D, {nchAxis, varAxis2});
    //=======================MC histograms Generated ================================================
    histos.add("ptHistogram_allcharge_gen", "ptHistogram_allcharge_gen", kTH1D, {ptAxis});
    histos.add("ptHistogramPion", "ptHistogramPion", kTH1D, {ptAxis});
    histos.add("ptHistogramKaon", "ptHistogramKaon", kTH1D, {ptAxis});
    histos.add("ptHistogramProton", "ptHistogramProton", kTH1D, {ptAxis});
    histos.add("hMC_Pt", ";#it{p}_{T} (GeV/#it{c})", kTH1D, {ptAxis});
    histos.add("MC_hZvtx_after_sel", ";#it{p}_{T} (GeV/#it{c})", kTH1D, {vtxZAxis});
    histos.add("hTOFbeta_gen_pion", "hTOFbeta_gen_pion", kTH2D, {pAxis, betaAxis});
    histos.add("hdEdx_gen_pion", "hdEdx_gen_pion", kTH2D, {pAxis, dEdxAxis});
    histos.add("hnch_gen_all", ";hnch_gen_all", kTH1D, {nchAxis});
    histos.add("hnch_gen", ";hnch_gen", kTH1D, {nchAxis});
    histos.add("hnch_gen_true", ";hnch_gen_true", kTH1D, {nchAxis});
    histos.add("hnch_gen_eta", ";hnch_gen_eta", kTH1D, {etaAxis});
    histos.add("hnch1", ";hnch1", kTH1D, {nchAxis});
    histos.add("hnch2", ";hnch2", kTH1D, {nchAxis});
    histos.add("hnch3", ";hnch3", kTH1D, {nchAxis});
    histos.add("hnch_pi", ";hnch_pi", kTH1D, {nchAxis});
    histos.add("hnch_ka", ";hnch_ka", kTH1D, {nchAxis});
    histos.add("hnch_pr", ";hnch_pr", kTH1D, {nchAxis});

    histos.add("hVar1x_gen_old", "hVar1x_gen_old", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2x_gen_old", "hVar2x_gen_old", kTH2D, {subAxis, nchAxis});
    histos.add("hVarx_gen_old", "hVarx_gen_old", kTH2D, {subAxis, nchAxis});
    histos.add("hVar1x_gen", "hVar1x_gen", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2x_gen", "hVar2x_gen", kTH2D, {subAxis, nchAxis});
    histos.add("hVarx_gen", "hVarx_gen", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2meanptx_gen", "hVar2meanptx_gen", kTH2D, {nchAxis, varAxis2});
    histos.add("hVar1pix_gen", "hVar1pix_gen", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2pix_gen", "hVar2pix_gen", kTH2D, {subAxis, nchAxis});
    histos.add("hVarpix_gen", "hVarpix_gen", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2meanptpix_gen", "hVar2meanptpix_gen", kTH2D, {nchAxis, varAxis2});
    histos.add("hVar1kx_gen", "hVar1kx_gen", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2kx_gen", "hVar2kx_gen", kTH2D, {subAxis, nchAxis});
    histos.add("hVarkx_gen", "hVarkx_gen", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2meanptkx_gen", "hVar2meanptkx_gen", kTH2D, {nchAxis, varAxis2});
    histos.add("hVar1px_gen", "hVar1px_gen", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2px_gen", "hVar2px_gen", kTH2D, {subAxis, nchAxis});
    histos.add("hVarpx_gen", "hVarpx_gen", kTH2D, {subAxis, nchAxis});
    histos.add("hVar2meanptpx_gen", "hVar2meanptpx_gen", kTH2D, {nchAxis, varAxis2});
    histos.add("hcent_nacc_rec", "hcent_nacc_rec", kTH2D, {centAxis, nchAxis});
    histos.add("hcent_nacc_gen", "hcent_nacc_gen", kTH2D, {centAxis, nchAxis});
    histos.add("hGenCentrality", "hGenCentrality", kTH1D, {centAxis});
    histos.add("hVtxZ_before_gen", "", kTH1F, {vtxZAxis});
    histos.add("hVtxZ_after_gen", "", kTH1F, {vtxZAxis});
    histos.add("hEta_gen", "", kTH1F, {etaAxis});
    histos.add("hEta_rec", "", kTH1F, {etaAxis});
    histos.add("hPt_gen", "", kTH1F, {ptAxis});
    histos.add("hPt_rec", "", kTH1F, {ptAxis});
  }
  // Configurables
  Configurable<float> cVtxZcut{"cVtxZcut", 10.f, "Vertex Z"};
  Configurable<float> cEtacut{"cEtacut", 0.8, "Eta cut"};
  Configurable<float> cPtmincut{"cPtmincut", 0.2, "Pt min cut"};
  Configurable<float> cPtmaxcut{"cPtmaxcut", 2.0, "Pt max cut"};
  Configurable<float> cDcaXYcut{"cDcaXYcut", 0.12, "DCA XY cut"};
  Configurable<float> cDcaZcut{"cDcaZcut", 0.3, "DCA Z cut"};
  Configurable<float> cCentmincut{"cCentmincut", 0.0, "Min cent cut"};
  Configurable<float> cCentmaxcut{"cCentmaxcut", 90.0, "Max cent cut"};
  Configurable<int> cTPCcrosscut{"cTPCcrosscut", 70, "TPC crossrows cut"};
  Configurable<int> cItsChiCut{"cItsChiCut", 70, "ITS chi2 cluster cut"};
  Configurable<int> cTpcChiCut{"cTpcChiCut", 70, "TPC chi2 cluster cut"};

  // Event selections
  Configurable<bool> cSel8Trig{"cSel8Trig", true, "Sel8 (T0A + T0C) Selection Run3"};
  Configurable<bool> cTFBorder{"cTFBorder", true, "Timeframe Border Selection"};
  Configurable<bool> cNoItsROBorder{"cNoItsROBorder", true, "No ITSRO Border Cut"};
  Configurable<bool> cItsTpcVtx{"cItsTpcVtx", true, "ITS+TPC Vertex Selection"};
  Configurable<bool> cPileupReject{"cPileupReject", true, "Pileup rejection"};
  Configurable<bool> cZVtxTimeDiff{"cZVtxTimeDiff", true, "z-vtx time diff selection"};
  Configurable<bool> cIsGoodITSLayers{"cIsGoodITSLayers", true, "Good ITS Layers All"};
  Configurable<bool> cItslayerall{"cItslayerall", true, "dead staves of ITS removed"};
  Configurable<bool> cvtxtofmatched{"cvtxtofmatched", true, "TOF vertex matched"};
  Configurable<bool> cfgRejEl{"cfgRejEl", true, "Rejected electrons"};

  // PID selection configurables
  Configurable<float> cPionPmincut{"cPionPmincut", 0.2, "pion min cut of pion"};
  Configurable<float> cKaonPmincut{"cKaonPmincut", 0.2, "kaon min cut of kaon"};
  Configurable<float> cProtonPmincut{"cProtonPmincut", 0.2, "proton min cut of proton"};
  Configurable<float> cPionPmaxcut{"cPionPmaxcut", 2.0, "pion min cut of pion"};
  Configurable<float> cKaonPmaxcut{"cKaonPmaxcut", 2.0, "kaon min cut of kaon"};
  Configurable<float> cProtonPmaxcut{"cProtonPmaxcut", 2.0, "proton min cut of proton"};
  Configurable<float> cPionPthcut{"cPionPthcut", 0.65, "pion threshold cut of pion"};
  Configurable<float> cKaonPthcut{"cKaonPthcut", 0.65, "kaon threshold cut of kaon"};
  Configurable<float> cProtonPthcut{"cProtonPthcut", 1.0, "proton threshold cut of proton"};
  Configurable<float> cNSigCut2{"cNSigCut2", 2.0, "nSigma cut (2)"};
  Configurable<float> cNSigCut3{"cNSigCut3", 3.0, "nSigma cut (3)"};
  Configurable<float> cElMinCut{"cElMinCut", -3.0, "electron min cut"};
  Configurable<float> cElMaxCut{"cElMaxCut", 5.0, "electron max cut"};
  Configurable<float> cTwoPtlCut2{"cTwoPtlCut2", 2.0, "n2ptl cut"};
  Configurable<float> cRapidityCut05{"cRapidityCut05", 0.5, "rapidity cut"};

  template <typename C>
  bool selCollision(C const& coll)
  {

    if (std::abs(coll.posZ()) > cVtxZcut) {
      return false;
    } // Reject the collisions with large vertex-z
    histos.fill(HIST("hEventcounter"), 2.);

    //  cent = coll.centFT0M(); //centrality for run3
    if (cSel8Trig && !coll.sel8()) {
      return false;
    } // require min bias trigger
    histos.fill(HIST("hEventcounter"), 3.);

    if (cTFBorder && !coll.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (cNoItsROBorder && !coll.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    histos.fill(HIST("trackSelRec"), 4);

    if (cPileupReject && !coll.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    histos.fill(HIST("trackSelRec"), 5);

    if (cZVtxTimeDiff && !coll.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    histos.fill(HIST("trackSelRec"), 6);

    if (cItsTpcVtx && !coll.selection_bit(aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    histos.fill(HIST("trackSelRec"), 7);

    //  if (cItslayerall && !coll.selection_bit(aod::evsel::kIsGoodITSLayersAll))         {return false;}
    histos.fill(HIST("trackSelRec"), 8);

    if (cvtxtofmatched && !coll.selection_bit(aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    histos.fill(HIST("trackSelRec"), 9);

    return true; // if all checks pass, accept the collision
  }

  template <typename T>
  bool selTrack(T const& track)
  {
    if (!track.isGlobalTrack()) {
      return false;
    } // accept only global tracks
    histos.fill(HIST("tracksel"), 2);

    //        if (std::fabs(track.dcaXY()) > cDcaXYcut)             {return false;}
    histos.fill(HIST("tracksel"), 3);

    //        if (std::fabs(track.dcaZ()) > cDcaZcut)               {return false;}
    histos.fill(HIST("tracksel"), 4);

    if (std::fabs(track.eta()) >= cEtacut) {
      return false;
    }
    histos.fill(HIST("tracksel"), 5);

    if (track.pt() < cPtmincut) {
      return false;
    }
    if (track.pt() > cPtmaxcut) {
      return false;
    }
    histos.fill(HIST("tracksel"), 6);

    //        if (track.tpcNClsCrossedRows() < cTPCcrosscut)        {return false;}
    histos.fill(HIST("tracksel"), 7);

    //        if (track.itsChi2NCl() > cItsChiCut)                  {return false;}
    histos.fill(HIST("tracksel"), 8);

    //        if (track.tpcChi2NCl() > cTpcChiCut)                  {return false;}
    histos.fill(HIST("tracksel"), 9);

    if (track.sign() == 0)
      return false;

    return true; // if all checks pass, accept the collision
  }

  template <typename T>
  bool rejEl(T const& track)
  {
    if (track.tpcNSigmaEl() > cElMinCut && track.tpcNSigmaEl() < cElMaxCut && std::fabs(track.tpcNSigmaPi()) > cNSigCut3 && std::fabs(track.tpcNSigmaKa()) > cNSigCut3 && std::fabs(track.tpcNSigmaPr()) > cNSigCut3) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool selProton(T const& track)
  {
    //! if pt < threshold (For tracks without TOF information)
    if (track.p() > cProtonPmincut && track.p() <= cProtonPthcut) {
      if (track.hasTPC() && std::fabs(track.tpcNSigmaPr()) < cNSigCut2 && std::fabs(track.tpcNSigmaPi()) > cNSigCut2 && std::fabs(track.tpcNSigmaKa()) > cNSigCut2) {
        return true;
      }
    }

    //! if pt < threshold (For tracks with TOF information)
    if (track.p() > cProtonPmincut && track.p() <= cProtonPthcut) {
      if (track.hasTOF() && std::fabs(track.tpcNSigmaPr()) < cNSigCut2 && std::fabs(track.tofNSigmaPr()) < cNSigCut2 && std::fabs(track.tpcNSigmaPi()) > cNSigCut2 && std::fabs(track.tpcNSigmaKa()) > cNSigCut2) {
        return true;
      }
    }

    //! if pt > threshold (For tracks with TOF information)
    if (track.p() > cProtonPthcut && track.p() <= cProtonPmaxcut) {
      if (track.hasTPC() && track.hasTOF() && std::fabs(track.tpcNSigmaPr()) < cNSigCut2 && std::fabs(track.tofNSigmaPr()) < cNSigCut2 && std::hypot(track.tofNSigmaPi(), track.tpcNSigmaPi()) > cNSigCut2 && std::hypot(track.tofNSigmaKa(), track.tpcNSigmaKa()) > cNSigCut2) {
        return true;
      }
    }

    return false;
  }

  template <typename T>
  bool selKaon(T const& track)
  {
    //! if pt < threshold (For tracks without TOF information)
    if (track.p() > cKaonPmincut && track.p() <= cKaonPthcut) {
      if (track.hasTPC() && std::fabs(track.tpcNSigmaKa()) < cNSigCut2 && std::fabs(track.tpcNSigmaPi()) > cNSigCut2 && std::fabs(track.tpcNSigmaPr()) > cNSigCut2) {
        return true;
      }
    }

    //! if pt < threshold (For tracks with TOF information)
    if (track.p() > cKaonPmincut && track.p() <= cKaonPthcut) {
      if (track.hasTOF() && std::fabs(track.tpcNSigmaKa()) < cNSigCut2 && std::fabs(track.tofNSigmaKa()) < cNSigCut2 && std::fabs(track.tpcNSigmaPi()) > cNSigCut2 && std::fabs(track.tpcNSigmaPr()) > cNSigCut2) {
        return true;
      }
    }

    //! if pt > threshold (For tracks with TOF information)
    if (track.p() > cKaonPthcut && track.p() <= cKaonPmaxcut) {
      if (track.hasTPC() && track.hasTOF() && std::fabs(track.tpcNSigmaKa()) < cNSigCut2 && std::fabs(track.tofNSigmaKa()) < cNSigCut2 && std::hypot(track.tofNSigmaPi(), track.tpcNSigmaPi()) > cNSigCut2 && std::hypot(track.tofNSigmaPr(), track.tpcNSigmaPr()) > cNSigCut2) {
        return true;
      }
    }

    return false;
  }

  template <typename T>
  bool selPion(T const& track)
  {
    //! if pt < threshold (For tracks without TOF information)
    if (track.p() > cPionPmincut && track.p() <= cPionPthcut) {
      if (track.hasTPC() && std::fabs(track.tpcNSigmaPi()) < cNSigCut2 && std::fabs(track.tpcNSigmaKa()) > cNSigCut2 && std::fabs(track.tpcNSigmaPr()) > cNSigCut2) {
        return true;
      }
    }

    //! if pt < threshold (For tracks with TOF information)
    if (track.p() > cPionPmincut && track.p() <= cPionPthcut) {
      if (track.hasTOF() && std::fabs(track.tpcNSigmaPi()) < cNSigCut2 && std::fabs(track.tofNSigmaPi()) < cNSigCut2 && std::fabs(track.tpcNSigmaKa()) > cNSigCut2 && std::fabs(track.tpcNSigmaPr()) > cNSigCut2) {
        return true;
      }
    }

    //! if pt > threshold (For tracks with TOF information)
    if (track.p() > cPionPthcut && track.p() <= cPionPmaxcut) {
      if (track.hasTPC() && track.hasTOF() && std::fabs(track.tpcNSigmaPi()) < cNSigCut2 && std::fabs(track.tofNSigmaPi()) < cNSigCut2 && std::hypot(track.tofNSigmaKa(), track.tpcNSigmaKa()) > cNSigCut2 && std::hypot(track.tofNSigmaPr(), track.tpcNSigmaPr()) > cNSigCut2) {
        return true;
      }
    }

    return false;
  }

  double getEfficiency(double pt, TH1D* ptHistogramAllchargeRec)
  {
    int bin = ptHistogramAllchargeRec->FindBin(pt);
    double eff = ptHistogramAllchargeRec->GetBinContent(bin);
    return (eff > 0) ? eff : 1e-6; // Avoid division by zero
  }

  //++++++++++++++++++++++++++++++++++++DATA CALCULATION +++++++++++++++++++++++++++++++++++++++++++++++++++++//

  void process(aod::MyCollision const& coll, aod::MyTracks const& inputTracks)
  {
    histos.fill(HIST("hEventcounter"), 1.);
    histos.fill(HIST("Data/hZvtx_before_sel"), coll.posZ());

    if (!selCollision(coll))
      return;
    {
      histos.fill(HIST("Data/hZvtx_after_sel8"), coll.posZ());
    }

    const auto cent = coll.centFT0C();
    histos.fill(HIST("Data/hCentrality"), cent);

    double nch = 0., nchPi = 0., nchKa = 0., nchPr = 0., nchAll = 0., nchAllBfCut = 0., nchEta = 0., nchPt = 0.;
    double q1 = 0., q2 = 0.;
    double q1Pi = 0., q2Pi = 0., q1Ka = 0., q2Ka = 0., q1Pr = 0., q2Pr = 0.;
    double var1 = 0., var2 = 0., twoParAllCharge = 0.;
    double var1Pi = 0., var2Pi = 0.;
    double var1Ka = 0., var2Ka = 0.;
    double var1Pr = 0., var2Pr = 0.;

    int sample = histos.get<TH1>(HIST("Data/hZvtx_after_sel8"))->GetEntries();
    sample = sample % 30; // subsample error estimation
    for (const auto& track : inputTracks) {
      nchAllBfCut += 1.;
      histos.fill(HIST("Data/hnchAll_bf_cut"), nchAllBfCut);

      histos.fill(HIST("tracksel"), 1);
      histos.fill(HIST("Data/hTPCchi2perCluster_before"), track.tpcChi2NCl());
      histos.fill(HIST("Data/hITSchi2perCluster_before"), track.itsChi2NCl());
      histos.fill(HIST("Data/hTPCCrossedrows_before"), track.tpcNClsCrossedRows());

      if (std::fabs(track.eta()) <= cEtacut) {
        nchEta++;
        histos.fill(HIST("Data/hnchTrue"), nchEta);
      }
      if (track.pt() >= cPtmincut && track.pt() <= cPtmaxcut) {
        nchPt += 1.;
        histos.fill(HIST("Data/hnchTrue_pt"), nchPt);
      }

      if (track.sign() == 0)
        continue;
      if (!selTrack(track))
        continue;

      nchAll += 1.;
      histos.fill(HIST("Data/hnchAll"), nchAll);
      histos.fill(HIST("Data/hDCAxy"), track.dcaXY());
      histos.fill(HIST("Data/hDCAz"), track.dcaZ());
      histos.fill(HIST("Data/hTPCCrossedrows_after"), track.tpcNClsCrossedRows());
      histos.fill(HIST("Data/hTPCchi2perCluster_after"), track.tpcChi2NCl());
      histos.fill(HIST("Data/hITSchi2perCluster_after"), track.itsChi2NCl());
      histos.fill(HIST("Data/hP"), track.p());
      histos.fill(HIST("Data/hPt"), track.pt());
      histos.fill(HIST("Data/hEta"), track.eta());
      histos.fill(HIST("Data/hPtDCAxy"), track.pt(), track.dcaXY());
      histos.fill(HIST("Data/hPtDCAz"), track.pt(), track.dcaZ());
      histos.fill(HIST("Data/hPtEta"), track.pt(), track.eta());
      histos.fill(HIST("Data/hPEta"), track.p(), track.eta());
      histos.fill(HIST("Data/hNsigmaTPC"), track.p(), track.tpcNSigmaPr());

      if (track.pt() >= cPtmincut || track.pt() <= cPtmaxcut) // do not change this (it is for different pt work)
      {
        nch += 1.;
        histos.fill(HIST("Data/hnch"), nch);
      }

      q1 += track.pt();
      q2 += (track.pt() * track.pt());

      // only TPC tracks: Pion, Kaon, Proton
      if (track.hasTPC() && std::abs(track.tpcNSigmaPi()) < cNSigCut3)
        histos.fill(HIST("Data/NSigamaTPCpion"), track.pt(), track.tpcNSigmaPi());
      if (track.hasTPC() && std::abs(track.tpcNSigmaKa()) < cNSigCut3)
        histos.fill(HIST("Data/NSigamaTPCkaon"), track.pt(), track.tpcNSigmaKa());
      if (track.hasTPC() && std::abs(track.tpcNSigmaPr()) < cNSigCut3)
        histos.fill(HIST("Data/NSigamaTPCproton"), track.pt(), track.tpcNSigmaPr());

      // only TOF tracks: Pion, Kaon, Proton
      if (track.hasTOF() && std::abs(track.tofNSigmaPi()) < cNSigCut3)
        histos.fill(HIST("Data/NSigamaTOFpion"), track.pt(), track.tofNSigmaPi());
      if (track.hasTOF() && std::abs(track.tofNSigmaKa()) < cNSigCut3)
        histos.fill(HIST("Data/NSigamaTOFkaon"), track.pt(), track.tofNSigmaKa());
      if (track.hasTOF() && std::abs(track.tofNSigmaPr()) < cNSigCut3)
        histos.fill(HIST("Data/NSigamaTOFproton"), track.pt(), track.tofNSigmaPr());

      if (track.hasTPC())
        histos.fill(HIST("Data/hdEdx"), track.p(), track.tpcSignal());
      if (track.hasTOF())
        histos.fill(HIST("Data/hTOFbeta"), track.p(), track.beta());

      //===================================pion==============================================================
      // only TPC+TOF tracks: Pion, Kaon, Proton
      if ((track.hasTPC() && std::abs(track.tpcNSigmaPi()) < cNSigCut3) && (track.hasTOF() && std::abs(track.tofNSigmaPi()) < cNSigCut3)) {
        histos.fill(HIST("Data/NSigamaTPCTOFpion"), track.tpcNSigmaPi(), track.tofNSigmaPi());

        histos.fill(HIST("Data/hdEdx_afterselection"), track.p(), track.tpcSignal());
        histos.fill(HIST("Data/hTOFbeta_afterselection"), track.p(), track.beta());
      }

      if (selPion(track)) {
        histos.fill(HIST("Data/hPtPion"), track.pt());
        histos.fill(HIST("Data/hEtaPion"), track.eta());
        histos.fill(HIST("Data/hyPion"), track.rapidity(massPi));
        histos.fill(HIST("Data/hPtyPion"), track.pt(), track.rapidity(massPi));
        nchPi += 1.;
        q1Pi += track.pt();
        q2Pi += (track.pt() * track.pt());

        if (track.beta() > 1)
          continue;
        histos.fill(HIST("Data/hdEdx_afterselection1"), track.p(), track.tpcSignal());
        histos.fill(HIST("Data/hTOFbeta_afterselection1"), track.p(), track.beta());
      }

      //===========================kaon===============================================================
      if ((track.hasTPC() && std::abs(track.tpcNSigmaKa()) < cNSigCut3) && (track.hasTOF() && std::abs(track.tofNSigmaKa()) < cNSigCut3)) {
        histos.fill(HIST("Data/NSigamaTPCTOFkaon"), track.tpcNSigmaKa(), track.tofNSigmaKa());
        histos.fill(HIST("Data/hdEdx_afterselection"), track.p(), track.tpcSignal());
        histos.fill(HIST("Data/hTOFbeta_afterselection"), track.p(), track.beta());
      }

      if (selKaon(track)) {
        histos.fill(HIST("Data/hPtKaon"), track.pt());
        histos.fill(HIST("Data/hEtaKaon"), track.eta());
        histos.fill(HIST("Data/hyKaon"), track.rapidity(massKa));
        histos.fill(HIST("Data/hPtyKaon"), track.pt(), track.rapidity(massKa));
        nchKa += 1.;
        q1Ka += track.pt();
        q2Ka += (track.pt() * track.pt());

        if (track.beta() > 1)
          continue;
        histos.fill(HIST("Data/hdEdx_afterselection1"), track.p(), track.tpcSignal());
        histos.fill(HIST("Data/hTOFbeta_afterselection1"), track.p(), track.beta());
      }

      //============================proton===========================================================
      if ((track.hasTPC() && std::abs(track.tpcNSigmaPr()) < cNSigCut3) && (track.hasTOF() && std::abs(track.tofNSigmaPr()) < cNSigCut3)) {
        histos.fill(HIST("Data/NSigamaTPCTOFproton"), track.tpcNSigmaPr(), track.tofNSigmaPr());
        histos.fill(HIST("Data/hdEdx_afterselection"), track.p(), track.tpcSignal());
        histos.fill(HIST("Data/hTOFbeta_afterselection"), track.p(), track.beta());
      }

      if (selProton(track)) {
        histos.fill(HIST("Data/hPtProton"), track.pt());
        histos.fill(HIST("Data/hEtaProton"), track.eta());
        histos.fill(HIST("Data/hyProton"), track.rapidity(massPr));
        histos.fill(HIST("Data/hPtyProton"), track.pt(), track.rapidity(massPr));
        nchPr += 1.;
        q1Pr += track.pt();
        q2Pr += (track.pt() * track.pt());

        if (track.beta() > 1)
          continue;
        histos.fill(HIST("Data/hdEdx_afterselection1"), track.p(), track.tpcSignal());
        histos.fill(HIST("Data/hTOFbeta_afterselection1"), track.p(), track.beta());
      }
    } // Track loop ends!
    histos.fill(HIST("Data/hcent_nacc"), cent, nchAll);

    if (nchAll < cTwoPtlCut2)
      return;
    var1 = (q1 * q1 - q2) / (nchAll * (nchAll - 1));
    var2 = (q1 / nchAll);

    //------------------ all charges-------------------------------------
    histos.fill(HIST("Data/hVar1"), sample, cent, var1);
    histos.fill(HIST("Data/hVar2"), sample, cent, var2);
    histos.fill(HIST("Data/hVarc"), sample, cent);
    histos.fill(HIST("Data/hVar2meanpt"), cent, var2);
    twoParAllCharge = (var1 - var2);
    histos.fill(HIST("Data/hVar"), nchAll, twoParAllCharge);

    //---------------------- pions ----------------------------------------
    if (nchPi >= cTwoPtlCut2) {
      var1Pi = (q1Pi * q1Pi - q2Pi) / (nchPi * (nchPi - 1));
      var2Pi = (q1Pi / nchPi);
    }

    //----------------------- kaons ---------------------------------------
    if (nchKa >= cTwoPtlCut2) {
      var1Ka = (q1Ka * q1Ka - q2Ka) / (nchKa * (nchKa - 1));
      var2Ka = (q1Ka / nchKa);
    }

    //---------------------------- protons ----------------------------------
    if (nchPr >= cTwoPtlCut2) {
      var1Pr = (q1Pr * q1Pr - q2Pr) / (nchPr * (nchPr - 1));
      var2Pr = (q1Pr / nchPr);
    }

    //========================centrality==========================================
    histos.fill(HIST("Data/hVar1pi"), sample, cent, var1Pi);
    histos.fill(HIST("Data/hVar2pi"), sample, cent, var2Pi);
    histos.fill(HIST("Data/hVar2meanptpi"), cent, var2Pi);

    histos.fill(HIST("Data/hVar1k"), sample, cent, var1Ka);
    histos.fill(HIST("Data/hVar2k"), sample, cent, var2Ka);
    histos.fill(HIST("Data/hVar2meanptk"), cent, var2Ka);

    histos.fill(HIST("Data/hVar1p"), sample, cent, var1Pr);
    histos.fill(HIST("Data/hVar2p"), sample, cent, var2Pr);
    histos.fill(HIST("Data/hVar2meanptp"), cent, var2Pr);

    //-----------------------nch-------------------------------------
    histos.fill(HIST("Data/hVar1x"), sample, nchAll, var1);
    histos.fill(HIST("Data/hVar2x"), sample, nchAll, var2);
    histos.fill(HIST("Data/hVarx"), sample, nchAll);
    histos.fill(HIST("Data/hVar2meanptx"), nchAll, var2);

    histos.fill(HIST("Data/hVar1pix"), sample, nchAll, var1Pi);
    histos.fill(HIST("Data/hVar2pix"), sample, nchAll, var2Pi);
    histos.fill(HIST("Data/hVarpix"), sample, nchPi);
    histos.fill(HIST("Data/hVar2meanptpix"), nchAll, var2Pi);

    histos.fill(HIST("Data/hVar1kx"), sample, nchAll, var1Ka);
    histos.fill(HIST("Data/hVar2kx"), sample, nchAll, var2Ka);
    histos.fill(HIST("Data/hVarkx"), sample, nchKa);
    histos.fill(HIST("Data/hVar2meanptkx"), nchAll, var2Ka);

    histos.fill(HIST("Data/hVar1px"), sample, nchAll, var1Pr);
    histos.fill(HIST("Data/hVar2px"), sample, nchAll, var2Pr);
    histos.fill(HIST("Data/hVarpx"), sample, nchPr);
    histos.fill(HIST("Data/hVar2meanptpx"), nchAll, var2Pr);

  } // event loop ends!

  PROCESS_SWITCH(EventMeanPtId, process, "process real data information", false);

  //++++++++++++++++++++++++++++++++++++MC Reconstructed +++++++++++++++++++++++++++++++++++++++++++++++++++++//

  SliceCache cache;
  Preslice<aod::McParticles> mcTrack = o2::aod::mcparticle::mcCollisionId;
  void processMcReco(aod::MyMCRecoCollision const& coll, aod::MyMCRecoTracks const& inputTracks, aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    (void)mcCollisions;
    if (!coll.has_mcCollision()) {
      return;
    }
    histos.fill(HIST("Rec/hZvtx_before_sel"), coll.posZ());
    histos.fill(HIST("hVtxZ_before_gen"), coll.mcCollision().posZ());

    if (cTFBorder && !coll.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return;
    }
    if (cNoItsROBorder && !coll.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return;
    }
    if (cPileupReject && !coll.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    if (cZVtxTimeDiff && !coll.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }
    if (cItsTpcVtx && !coll.selection_bit(aod::evsel::kIsVertexITSTPC)) {
      return;
    }
    if (cvtxtofmatched && !coll.selection_bit(aod::evsel::kIsVertexTOFmatched)) {
      return;
    }
    if (std::abs(coll.posZ()) > cVtxZcut) {
      return;
    }

    if (!coll.sel8()) {
      return;
    }
    float cent = coll.centFT0C();
    histos.fill(HIST("Rec/hZvtx_after_sel8"), coll.posZ());

    double nch = 0., nchPi = 0., nchKa = 0., nchPr = 0., nchAll = 0., nchAllBfCut = 0., nchEta = 0., nchPt = 0.;
    double q1 = 0., q2 = 0.;
    double q1Pi = 0., q2Pi = 0., q1Ka = 0., q2Ka = 0., q1Pr = 0., q2Pr = 0.;
    double var1 = 0., var2 = 0., twoParAllCharge = 0.;
    double var1Pi = 0., var2Pi = 0., var1Ka = 0., var2Ka = 0., var1Pr = 0., var2Pr = 0.;
    double sumPtWeight = 0., sumWeight = 0., sumPtPtWeight = 0., var1Eff = 0., var2Eff = 0.;
    double sumPtWeightPi = 0., sumWeightPi = 0., sumPtPtWeightPi = 0., var1EffPi = 0., var2EffPi = 0.;
    double sumPtWeightKa = 0., sumWeightKa = 0., sumPtPtWeightKa = 0., var1EffKa = 0., var2EffKa = 0.;
    double sumPtWeightPr = 0., sumWeightPr = 0., sumPtPtWeightPr = 0., var1EffPr = 0., var2EffPr = 0.;

    int sample = histos.get<TH1>(HIST("Rec/hZvtx_after_sel8"))->GetEntries();
    sample = sample % 30;

    for (const auto& track : inputTracks) {
      nchAllBfCut += 1.;
      histos.fill(HIST("Rec/hnchAll_bf_cut"), nchAllBfCut);
      histos.fill(HIST("Rec/hTPCchi2perCluster_before"), track.tpcChi2NCl());
      histos.fill(HIST("Rec/hITSchi2perCluster_before"), track.itsChi2NCl());
      histos.fill(HIST("Rec/hTPCCrossedrows_before"), track.tpcNClsCrossedRows());

      if (std::fabs(track.eta()) <= cEtacut) {
        nchEta++;
        histos.fill(HIST("Rec/hnchTrue"), nchEta);
      }
      if (track.pt() >= cPtmincut && track.pt() <= cPtmaxcut) {
        nchPt += 1.;
        histos.fill(HIST("Rec/hnchTrue_pt"), nchPt);
      }
      if (!track.isGlobalTrack())
        continue;
      if (std::fabs(track.eta()) > cEtacut)
        continue;
      if ((track.pt() <= cPtmincut) || (track.pt() >= cPtmaxcut))
        continue;
      if (track.sign() == 0)
        continue;
      //          if (std::fabs(track.y()) > 0.5) continue;
      histos.fill(HIST("hPt_rec"), track.pt());
      histos.fill(HIST("hEta_rec"), track.eta());

      auto mcParticle = track.mcParticle();
      nchAll += 1.;
      histos.fill(HIST("Rec/hnchAll"), nchAll);
      histos.fill(HIST("ptHistogramAllchargeRec"), track.pt());
      histos.fill(HIST("Rec/hDCAxy"), track.dcaXY());
      histos.fill(HIST("Rec/hDCAz"), track.dcaZ());
      histos.fill(HIST("Rec/hTPCCrossedrows_after"), track.tpcNClsCrossedRows());
      histos.fill(HIST("Rec/hTPCchi2perCluster_after"), track.tpcChi2NCl());
      histos.fill(HIST("Rec/hITSchi2perCluster_after"), track.itsChi2NCl());
      histos.fill(HIST("Rec/hP"), track.p());
      histos.fill(HIST("Rec/hPt"), track.pt());
      histos.fill(HIST("Rec/hEta"), track.eta());
      histos.fill(HIST("Rec/hPtDCAxy"), track.pt(), track.dcaXY());
      histos.fill(HIST("Rec/hPtDCAz"), track.pt(), track.dcaZ());
      histos.fill(HIST("Rec/hPtEta"), track.pt(), track.eta());
      histos.fill(HIST("Rec/hPEta"), track.p(), track.eta());
      histos.fill(HIST("Rec/hNsigmaTPC"), track.p(), track.tpcNSigmaPr());

      if (track.pt() >= cPtmincut || track.pt() <= cPtmaxcut) // do not change this (it is for different pt work)
      {
        nch += 1.;
        histos.fill(HIST("Rec/hnch"), nch);
      }
      q1 += track.pt();
      q2 += (track.pt() * track.pt());

      double eff = getEfficiency(track.pt(), ptHistogramAllchargeRec);
      //  LOGF(info, " with value %.2f",  eff);
      sumPtWeight += track.pt() / eff;
      sumPtPtWeight += (track.pt() * track.pt()) / (eff * eff);
      sumWeight += 1. / eff;

      if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus)
        histos.fill(HIST("ptHistogramPionrec_pdg"), track.pt());
      if (std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus)
        histos.fill(HIST("ptHistogramKaonrec_pdg"), track.pt());
      if (std::abs(mcParticle.pdgCode()) == PDG_t::kProton)
        histos.fill(HIST("ptHistogramProtonrec_pdg"), track.pt());

      if (cfgRejEl == false && rejEl(track)) {
        return;
      }

      // only TPC tracks: Pion, Kaon, Proton
      if (track.hasTPC() && std::abs(track.tpcNSigmaPi()) < cNSigCut3)
        histos.fill(HIST("Rec/NSigamaTPCpion"), track.pt(), track.tpcNSigmaPi());
      if (track.hasTPC() && std::abs(track.tpcNSigmaKa()) < cNSigCut3)
        histos.fill(HIST("Rec/NSigamaTPCkaon"), track.pt(), track.tpcNSigmaKa());
      if (track.hasTPC() && std::abs(track.tpcNSigmaPr()) < cNSigCut3)
        histos.fill(HIST("Rec/NSigamaTPCproton"), track.pt(), track.tpcNSigmaPr());

      // only TOF tracks: Pion, Kaon, Proton
      if (track.hasTOF() && std::abs(track.tofNSigmaPi()) < cNSigCut3)
        histos.fill(HIST("Rec/NSigamaTOFpion"), track.pt(), track.tofNSigmaPi());
      if (track.hasTOF() && std::abs(track.tofNSigmaKa()) < cNSigCut3)
        histos.fill(HIST("Rec/NSigamaTOFkaon"), track.pt(), track.tofNSigmaKa());
      if (track.hasTOF() && std::abs(track.tofNSigmaPr()) < cNSigCut3)
        histos.fill(HIST("Rec/NSigamaTOFproton"), track.pt(), track.tofNSigmaPr());

      if (track.hasTPC())
        histos.fill(HIST("Rec/hdEdx"), track.p(), track.tpcSignal());
      if (track.hasTOF())
        histos.fill(HIST("Rec/hTOFbeta"), track.p(), track.beta());
      if (track.hasTPC())
        histos.fill(HIST("hdEdx_afterselection_rec_beforepidcut"), track.p(), track.tpcSignal());
      if (track.hasTOF())
        histos.fill(HIST("hTOFbeta_afterselection_rec_beforepidcut"), track.p(), track.beta());

      //===================================pion==============================================================
      if ((track.hasTPC() && std::abs(track.tpcNSigmaPi()) < cNSigCut3) && (track.hasTOF() && std::abs(track.tofNSigmaPi()) < cNSigCut3)) {
        histos.fill(HIST("Rec/NSigamaTPCTOFpion"), track.tpcNSigmaPi(), track.tofNSigmaPi());

        histos.fill(HIST("Rec/hdEdx_afterselection"), track.p(), track.tpcSignal());
        histos.fill(HIST("Rec/hTOFbeta_afterselection"), track.p(), track.beta());
      }

      if (selPion(track)) {
        if (std::fabs(track.y()) > cRapidityCut05)
          continue;
        if (track.beta() > 1)
          continue;
        histos.fill(HIST("ptHistogramPionrec"), track.pt());
        histos.fill(HIST("Rec/hPtPion"), track.pt());
        histos.fill(HIST("Rec/hEtaPion"), track.eta());
        histos.fill(HIST("Rec/hyPion"), track.rapidity(massPi));
        histos.fill(HIST("Rec/hPtyPion"), track.pt(), track.rapidity(massPi));
        histos.fill(HIST("NSigamaTPCpion_rec"), track.p(), track.tpcNSigmaPi());
        histos.fill(HIST("NSigamaTOFpion_rec"), track.p(), track.tofNSigmaPi());
        histos.fill(HIST("NSigamaTPCTOFpion_rec"), track.tpcNSigmaPi(), track.tofNSigmaPi());
        histos.fill(HIST("Rec/hdEdx_afterselection1"), track.p(), track.tpcSignal());
        histos.fill(HIST("Rec/hTOFbeta_afterselection1"), track.p(), track.beta());
        if (std::abs(track.mcParticle().pdgCode()) == PDG_t::kPiPlus) {
          histos.fill(HIST("ptHistogramPionrec_purity"), track.pt());
        }
        nchPi += 1.;
        q1Pi += track.pt();
        q2Pi += (track.pt() * track.pt());

        double effPi = getEfficiency(track.pt(), ptHistogramPionrec);
        //  LOGF(info, " with value %.2f",  eff);
        sumPtWeightPi += track.pt() / effPi;
        sumPtPtWeightPi += (track.pt() * track.pt()) / (effPi * effPi);
        sumWeightPi += 1. / effPi;

        histos.fill(HIST("hPyPion_rec"), track.p(), track.rapidity(massPi));
        histos.fill(HIST("hPtyPion_rec"), track.pt(), track.rapidity(massPi));
      }

      //===========================kaon===============================================================

      if ((track.hasTPC() && std::abs(track.tpcNSigmaKa()) < cNSigCut3) && (track.hasTOF() && std::abs(track.tofNSigmaKa()) < cNSigCut3)) {
        histos.fill(HIST("Rec/NSigamaTPCTOFkaon"), track.tpcNSigmaKa(), track.tofNSigmaKa());
        histos.fill(HIST("Rec/hdEdx_afterselection"), track.p(), track.tpcSignal());
        histos.fill(HIST("Rec/hTOFbeta_afterselection"), track.p(), track.beta());
      }

      if (selKaon(track)) {
        if (std::fabs(track.y()) > cRapidityCut05)
          continue;
        if (track.beta() > 1)
          continue;
        histos.fill(HIST("ptHistogramKaonrec"), track.pt());
        histos.fill(HIST("Rec/hPtKaon"), track.pt());
        histos.fill(HIST("Rec/hEtaKaon"), track.eta());
        histos.fill(HIST("Rec/hyKaon"), track.rapidity(massKa));
        histos.fill(HIST("Rec/hPtyKaon"), track.pt(), track.rapidity(massKa));
        histos.fill(HIST("NSigamaTPCkaon_rec"), track.p(), track.tpcNSigmaKa());
        histos.fill(HIST("NSigamaTOFkaon_rec"), track.p(), track.tofNSigmaKa());
        histos.fill(HIST("NSigamaTPCTOFkaon_rec"), track.tpcNSigmaKa(), track.tofNSigmaKa());
        histos.fill(HIST("Rec/hdEdx_afterselection1"), track.p(), track.tpcSignal());
        histos.fill(HIST("Rec/hTOFbeta_afterselection1"), track.p(), track.beta());
        if (std::abs(track.mcParticle().pdgCode()) == PDG_t::kKPlus) {
          histos.fill(HIST("ptHistogramKaonrec_purity"), track.pt());
        }
        nchKa += 1.;
        q1Ka += track.pt();
        q2Ka += (track.pt() * track.pt());

        double effKa = getEfficiency(track.pt(), ptHistogramKaonrec);
        //  LOGF(info, " with value %.2f",  eff);
        sumPtWeightKa += track.pt() / effKa;
        sumPtPtWeightKa += (track.pt() * track.pt()) / (effKa * effKa);
        sumWeightKa += 1. / effKa;

        histos.fill(HIST("hPyKaon_rec"), track.p(), track.rapidity(massKa));
        histos.fill(HIST("hPtyKaon_rec"), track.pt(), track.rapidity(massKa));
      }

      //============================proton===========================================================

      if ((track.hasTPC() && std::abs(track.tpcNSigmaPr()) < cNSigCut3) && (track.hasTOF() && std::abs(track.tofNSigmaPr()) < cNSigCut3)) {
        histos.fill(HIST("Rec/NSigamaTPCTOFproton"), track.tpcNSigmaPr(), track.tofNSigmaPr());
        histos.fill(HIST("Rec/hdEdx_afterselection"), track.p(), track.tpcSignal());
        histos.fill(HIST("Rec/hTOFbeta_afterselection"), track.p(), track.beta());
      }

      if (selProton(track)) {
        if (std::fabs(track.y()) > cRapidityCut05)
          continue;
        if (track.beta() > 1)
          continue;
        histos.fill(HIST("ptHistogramProtonrec"), track.pt());
        histos.fill(HIST("Rec/hPtProton"), track.pt());
        histos.fill(HIST("Rec/hEtaProton"), track.eta());
        histos.fill(HIST("Rec/hyProton"), track.rapidity(massPr));
        histos.fill(HIST("Rec/hPtyProton"), track.pt(), track.rapidity(massPr));
        histos.fill(HIST("NSigamaTPCproton_rec"), track.p(), track.tpcNSigmaPr());
        histos.fill(HIST("NSigamaTOFproton_rec"), track.p(), track.tofNSigmaPr());
        histos.fill(HIST("NSigamaTPCTOFproton_rec"), track.tpcNSigmaPr(), track.tofNSigmaPr());
        histos.fill(HIST("Rec/hdEdx_afterselection1"), track.p(), track.tpcSignal());
        histos.fill(HIST("Rec/hTOFbeta_afterselection1"), track.p(), track.beta());
        if (std::abs(track.mcParticle().pdgCode()) == PDG_t::kProton) {
          histos.fill(HIST("ptHistogramProtonrec_purity"), track.pt());
        }
        nchPr += 1.;
        q1Pr += track.pt();
        q2Pr += (track.pt() * track.pt());

        double effPr = getEfficiency(track.pt(), ptHistogramProtonrec);
        //  LOGF(info, " with value %.2f",  eff);
        sumPtWeightPr += track.pt() / effPr;
        sumPtPtWeightPr += (track.pt() * track.pt()) / (effPr * effPr);
        sumWeightPr += 1. / effPr;

        histos.fill(HIST("hPyProton_rec"), track.p(), track.rapidity(massPr));
        histos.fill(HIST("hPtyProton_rec"), track.pt(), track.rapidity(massPr));
      }

    } // loop over tracks
    histos.fill(HIST("Rec/hcent_nacc"), cent, nchAll);

    if (nchAll < cTwoPtlCut2)
      return;
    var1 = (q1 * q1 - q2) / (nchAll * (nchAll - 1));
    var2 = (q1 / nchAll);

    //------------------ Efficiency corrected histograms ---------------

    var1Eff = (sumPtWeight * sumPtWeight - sumPtPtWeight) / (sumWeight * (sumWeight - 1));
    var2Eff = (sumPtWeight / sumWeight);

    histos.fill(HIST("Rec/hVar1"), sample, cent, var1);
    histos.fill(HIST("Rec/hVar2"), sample, cent, var2);
    histos.fill(HIST("Rec/hVarc"), sample, cent);
    histos.fill(HIST("Rec/hVar2meanpt"), cent, var2);
    twoParAllCharge = (var1 - var2);
    histos.fill(HIST("Rec/hVar"), nchAll, twoParAllCharge);

    //---------------------- pions ----------------------------------------
    if (nchPi >= cTwoPtlCut2) {
      var1Pi = (q1Pi * q1Pi - q2Pi) / (nchPi * (nchPi - 1));
      var2Pi = (q1Pi / nchPi);

      var1EffPi = (sumPtWeightPi * sumPtWeightPi - sumPtPtWeightPi) / (sumWeightPi * (sumWeightPi - 1));
      var2EffPi = (sumPtWeightPi / sumWeightPi);
    }
    //----------------------- kaons ---------------------------------------
    if (nchKa >= cTwoPtlCut2) {
      var1Ka = (q1Ka * q1Ka - q2Ka) / (nchKa * (nchKa - 1));
      var2Ka = (q1Ka / nchKa);

      var1EffKa = (sumPtWeightKa * sumPtWeightKa - sumPtPtWeightKa) / (sumWeightKa * (sumWeightKa - 1));
      var2EffKa = (sumPtWeightKa / sumWeightKa);
    }
    //---------------------------- protons ----------------------------------
    if (nchPr >= cTwoPtlCut2) {
      var1Pr = (q1Pr * q1Pr - q2Pr) / (nchPr * (nchPr - 1));
      var2Pr = (q1Pr / nchPr);

      var1EffPr = (sumPtWeightPr * sumPtWeightPr - sumPtPtWeightPr) / (sumWeightPr * (sumWeightPr - 1));
      var2EffPr = (sumPtWeightPr / sumWeightPr);
    }
    //========================centrality==========================================

    histos.fill(HIST("Rec/hVar1pi"), sample, cent, var1Pi);
    histos.fill(HIST("Rec/hVar2pi"), sample, cent, var2Pi);
    histos.fill(HIST("Rec/hVar2meanptpi"), cent, var2Pi);
    histos.fill(HIST("Rec/hVar1k"), sample, cent, var1Ka);
    histos.fill(HIST("Rec/hVar2k"), sample, cent, var2Ka);
    histos.fill(HIST("Rec/hVar2meanptk"), cent, var2Ka);
    histos.fill(HIST("Rec/hVar1p"), sample, cent, var1Pr);
    histos.fill(HIST("Rec/hVar2p"), sample, cent, var2Pr);
    histos.fill(HIST("Rec/hVar2meanptp"), cent, var2Pr);

    //-----------------------nch-------------------------------------
    histos.fill(HIST("Rec/hVar1x"), sample, nchAll, var1);
    histos.fill(HIST("Rec/hVar2x"), sample, nchAll, var2);
    histos.fill(HIST("Rec/hVarx"), sample, nchAll);
    histos.fill(HIST("Rec/hVar2meanptx"), nchAll, var2);
    histos.fill(HIST("Rec/hVar1pix"), sample, nchAll, var1Pi);
    histos.fill(HIST("Rec/hVar2pix"), sample, nchAll, var2Pi);
    histos.fill(HIST("Rec/hVarpix"), sample, nchPi);
    histos.fill(HIST("Rec/hVar2meanptpix"), nchAll, var2Pi);
    histos.fill(HIST("Rec/hVar1kx"), sample, nchAll, var1Ka);
    histos.fill(HIST("Rec/hVar2kx"), sample, nchAll, var2Ka);
    histos.fill(HIST("Rec/hVarkx"), sample, nchKa);
    histos.fill(HIST("Rec/hVar2meanptkx"), nchAll, var2Ka);
    histos.fill(HIST("Rec/hVar1px"), sample, nchAll, var1Pr);
    histos.fill(HIST("Rec/hVar2px"), sample, nchAll, var2Pr);
    histos.fill(HIST("Rec/hVarpx"), sample, nchPr);
    histos.fill(HIST("Rec/hVar2meanptpx"), nchAll, var2Pr);

    histos.fill(HIST("hEffVar1x"), sample, nchAll, var1Eff);
    histos.fill(HIST("hEffVar2x"), sample, nchAll, var2Eff);
    histos.fill(HIST("hEffVarx"), sample, nchAll);
    histos.fill(HIST("hEffVar2Meanptx"), nchAll, var2Eff);

    histos.fill(HIST("hEffVar1pix"), sample, nchAll, var1EffPi);
    histos.fill(HIST("hEffVar2pix"), sample, nchAll, var2EffPi);
    histos.fill(HIST("hEffVarpix"), sample, nchAll);
    histos.fill(HIST("hEffVar2Meanptpix"), nchAll, var2EffPi);

    histos.fill(HIST("hEffVar1kx"), sample, nchAll, var1EffKa);
    histos.fill(HIST("hEffVar2kx"), sample, nchAll, var2EffKa);
    histos.fill(HIST("hEffVarkx"), sample, nchAll);
    histos.fill(HIST("hEffVar2Meanptkx"), nchAll, var2EffKa);

    histos.fill(HIST("hEffVar1px"), sample, nchAll, var1EffPr);
    histos.fill(HIST("hEffVar2px"), sample, nchAll, var2EffPr);
    histos.fill(HIST("hEffVarpx"), sample, nchAll);
    histos.fill(HIST("hEffVar2Meanptpx"), nchAll, var2EffPr);

    //================= generated level==============================

    const auto& mccolgen = coll.mcCollision_as<aod::McCollisions>();
    if (std::abs(mccolgen.posZ()) > cVtxZcut) {
      return;
    }
    const auto& mcpartgen = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mccolgen.globalIndex(), cache);
    histos.fill(HIST("hVtxZ_after_gen"), mccolgen.posZ());

    double nchGen = 0., nchGenAll = 0., nchGenTrue = 0.;
    double nchPiGen = 0., nchKaGen = 0., nchPrGen = 0.;
    double nch1 = 0., nch2 = 0., nch3 = 0.;
    double q1AllGen = 0, q2AllGen = 0.;
    double q1PiGen = 0, q2PiGen = 0, q1KaGen = 0, q2KaGen = 0, q1PrGen = 0, q2PrGen = 0;
    double var1AllGen = 0, var2AllGen = 0.;
    double var1PiGen = 0, var2PiGen = 0, var1KaGen = 0, var2KaGen = 0, var1PrGen = 0, var2PrGen = 0;

    int sampleGen = histos.get<TH1>(HIST("hVtxZ_after_gen"))->GetEntries();
    sampleGen = sampleGen % 30;

    for (const auto& mcpart : mcpartgen) {
      // auto  pdgcode = std::abs(mcpart.pdgCode());
      if (!mcpart.isPhysicalPrimary()) {
        continue;
      }
      nch1++;
      histos.fill(HIST("hnch1"), nch1);
      nch2++;
      histos.fill(HIST("hnch2"), nch2);
      nch3++;
      histos.fill(HIST("hnch3"), nch3);

      int pid = mcpart.pdgCode();
      auto sign = 0;
      auto* pd = pdg->GetParticle(pid);
      if (pd != nullptr) {
        sign = pd->Charge() / 3.;
      }
      if (sign == 0) {
        continue;
      }
      //    histos.fill(HIST("gen_hSign"), sign);
      if (std::fabs(mcpart.eta()) > cEtacut)
        continue;
      nchGenTrue++;
      histos.fill(HIST("hnch_gen_true"), nchGenTrue);
      if ((mcpart.pt() <= cPtmincut) || (mcpart.pt() >= cPtmaxcut))
        continue;
      histos.fill(HIST("hPt_gen"), mcpart.pt());
      histos.fill(HIST("hEta_gen"), mcpart.eta());
      histos.fill(HIST("ptHistogram_allcharge_gen"), mcpart.pt());
      nchGenAll += 1.;
      q1AllGen += mcpart.pt();
      q2AllGen += (mcpart.pt() * mcpart.pt());
      histos.fill(HIST("hnch_gen_all"), nchGenAll);
      if (std::fabs(mcpart.y()) < cRapidityCut05) {

        if (mcpart.pdgCode() == PDG_t::kPiPlus || mcpart.pdgCode() == PDG_t::kPiMinus) {
          histos.fill(HIST("ptHistogramPion"), mcpart.pt());
          nchPiGen += 1.;
          q1PiGen += mcpart.pt();
          q2PiGen += (mcpart.pt() * mcpart.pt());
          histos.fill(HIST("hnch_pi"), nchPiGen);
        }

        if (mcpart.pdgCode() == PDG_t::kKPlus || mcpart.pdgCode() == PDG_t::kKMinus) {
          histos.fill(HIST("ptHistogramKaon"), mcpart.pt());
          nchKaGen += 1.;
          q1KaGen += mcpart.pt();
          q2KaGen += (mcpart.pt() * mcpart.pt());
          histos.fill(HIST("hnch_ka"), nchKaGen);
        }

        if (mcpart.pdgCode() == PDG_t::kProton || mcpart.pdgCode() == PDG_t::kProtonBar) {
          histos.fill(HIST("ptHistogramProton"), mcpart.pt());
          nchPrGen += 1.;
          q1PrGen += mcpart.pt();
          q2PrGen += (mcpart.pt() * mcpart.pt());
          histos.fill(HIST("hnch_pr"), nchPrGen);
        }

      } //|y| < 0.5 cut ends!

    } // particle
    histos.fill(HIST("hcent_nacc_gen"), cent, nchGen);

    if (nchGenAll < cTwoPtlCut2)
      return;
    var1AllGen = (q1AllGen * q1AllGen - q2AllGen) / (nchGenAll * (nchGenAll - 1));
    var2AllGen = (q1AllGen / nchGenAll);

    if (nchPiGen >= cTwoPtlCut2) {
      var1PiGen = (q1PiGen * q1PiGen - q2PiGen) / (nchPiGen * (nchPiGen - 1));
      var2PiGen = (q1PiGen / nchPiGen);
    }

    //----------------------- kaons ---------------------------------------
    if (nchKaGen >= cTwoPtlCut2) {
      var1KaGen = (q1KaGen * q1KaGen - q2KaGen) / (nchKaGen * (nchKaGen - 1));
      var2KaGen = (q1KaGen / nchKaGen);
    }
    //---------------------------- protons ----------------------------------
    if (nchPrGen >= cTwoPtlCut2) {
      var1PrGen = (q1PrGen * q1PrGen - q2PrGen) / (nchPrGen * (nchPrGen - 1));
      var2PrGen = (q1PrGen / nchPrGen);
    }
    //-----------------------nch-------------------------------------
    histos.fill(HIST("hVar1x_gen"), sampleGen, nchGenAll, var1AllGen);
    histos.fill(HIST("hVar2x_gen"), sampleGen, nchGenAll, var2AllGen);
    histos.fill(HIST("hVarx_gen"), sampleGen, nchGenAll);
    histos.fill(HIST("hVar2meanptx_gen"), nchGenAll, var2AllGen);
    histos.fill(HIST("hVar1pix_gen"), sampleGen, nchGenAll, var1PiGen);
    histos.fill(HIST("hVar2pix_gen"), sampleGen, nchGenAll, var2PiGen);
    histos.fill(HIST("hVarpix_gen"), sampleGen, nchPiGen);
    histos.fill(HIST("hVar2meanptpix_gen"), nchGenAll, var2PiGen);
    histos.fill(HIST("hVar1kx_gen"), sampleGen, nchGenAll, var1KaGen);
    histos.fill(HIST("hVar2kx_gen"), sampleGen, nchGenAll, var2KaGen);
    histos.fill(HIST("hVarkx_gen"), sampleGen, nchKaGen);
    histos.fill(HIST("hVar2meanptkx_gen"), nchGenAll, var2KaGen);
    histos.fill(HIST("hVar1px_gen"), sampleGen, nchGenAll, var1PrGen);
    histos.fill(HIST("hVar2px_gen"), sampleGen, nchGenAll, var2PrGen);
    histos.fill(HIST("hVarpx_gen"), sampleGen, nchPrGen);
    histos.fill(HIST("hVar2meanptpx_gen"), nchGenAll, var2PrGen);

  } // void process
  PROCESS_SWITCH(EventMeanPtId, processMcReco, "Process reconstructed", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<EventMeanPtId>(cfgc)};
  return workflow;
}
