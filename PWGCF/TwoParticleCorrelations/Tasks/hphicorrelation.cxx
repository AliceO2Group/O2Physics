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
///
/// \brief
///
/// \author Rahul Verma (rahul.verma@iitb.ac.in) :: Sadhana Dash (sadhana@phy.iitb.ac.in)

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include <TLorentzVector.h>
#include <typeinfo>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Delta phi in [-pi/2,3pi/2] calculator
Double_t ComputeDeltaPhi(Double_t particlePhi, Double_t leadPhi)
{
  // Make inputs lie in range [0,2pi]
  if (leadPhi < 0.0) {
    leadPhi = leadPhi + 2.0 * TMath::Pi();
  }
  if (particlePhi < 0.0) {
    particlePhi = particlePhi + 2.0 * TMath::Pi();
  }

  Double_t dphi = 0; // default value for particlePhi = leadPhi
  if (particlePhi > leadPhi && particlePhi <= (TMath::Pi() + leadPhi)) {
    dphi = particlePhi - leadPhi;
  }
  if (particlePhi > leadPhi && particlePhi > (TMath::Pi() + leadPhi)) {
    dphi = particlePhi - leadPhi - 2 * TMath::Pi();
  }
  if (leadPhi > particlePhi && leadPhi <= (TMath::Pi() + particlePhi)) {
    dphi = particlePhi - leadPhi;
  }
  if (leadPhi > particlePhi && leadPhi > (TMath::Pi() + particlePhi)) {
    dphi = particlePhi - leadPhi + 2 * TMath::Pi();
  }

  // changing range of dphi from [-pi, pi] to [-pi/2,3pi/2]
  if (dphi < -TMath::Pi() / 2) {
    dphi = dphi + 2.0 * TMath::Pi();
  }
  return dphi;
}

// Start--Mixing Structures ***************************************************************
struct StoredColl {
  int Bin = -1;
  double PosZ;
  double CentFT0C;
  double GlobalIndex = -1;
  int RunNumber = -1;
  std::string DfName = "DF_-1";
  double DfNo = -1;

  int bin() const { return Bin; }
  double posZ() const { return PosZ; }
  double centFT0C() const { return CentFT0C; }
  double globalIndex() const { return GlobalIndex; }
  int runNumber() const { return RunNumber; }
  std::string dfName() const { return DfName; }
  double dfNo() const { return DfNo; }
};

struct myStoredTrack {
  int64_t GlobalIndex;
  int64_t globalIndex() const { return GlobalIndex; }
  int64_t CollId;
  int64_t collisionId() const { return CollId; }
  float Px;
  float px() const { return Px; }
  float Py;
  float py() const { return Py; }
  float Pz;
  float pz() const { return Pz; }
  float P;
  float p() const { return P; }
  float Pt;
  float pt() const { return Pt; }
  float Signed1Pt;
  float signed1Pt() const { return Signed1Pt; }
  float Eta;
  float eta() const { return Eta; }
  float Phi;
  float phi() const { return Phi; }
  float DcaXY;
  float dcaXY() const { return DcaXY; }
  float DcaZ;
  float dcaZ() const { return DcaZ; }
  bool HasTOF;
  bool hasTOF() const { return HasTOF; }
  float Beta;
  float beta() const { return Beta; }
  float TpcSignal;
  float tpcSignal() const { return TpcSignal; }
  float TpcInnerParam;
  float tpcInnerParam() const { return TpcInnerParam; }
  float TofExpMom;
  float tofExpMom() const { return TofExpMom; }
  int16_t TpcNClsCrossedRows;
  int16_t tpcNClsCrossedRows() const { return TpcNClsCrossedRows; }
  bool IsGlobalTrack;
  bool isGlobalTrack() const { return IsGlobalTrack; }
  int Sign;
  int sign() const { return Sign; }
  float TpcNSigmaPi;
  float tpcNSigmaPi() const { return TpcNSigmaPi; }
  float TpcNSigmaKa;
  float tpcNSigmaKa() const { return TpcNSigmaKa; }
  float TpcNSigmaPr;
  float tpcNSigmaPr() const { return TpcNSigmaPr; }
  float TpcNSigmaEl;
  float tpcNSigmaEl() const { return TpcNSigmaEl; }
  float TpcNSigmaDe;
  float tpcNSigmaDe() const { return TpcNSigmaDe; }
  float TofNSigmaPi;
  float tofNSigmaPi() const { return TofNSigmaPi; }
  float TofNSigmaKa;
  float tofNSigmaKa() const { return TofNSigmaKa; }
  float TofNSigmaPr;
  float tofNSigmaPr() const { return TofNSigmaPr; }
  float TofNSigmaEl;
  float tofNSigmaEl() const { return TofNSigmaEl; }
  float TofNSigmaDe;
  float tofNSigmaDe() const { return TofNSigmaDe; }
  int DfNo;
  int dfNo() const { return DfNo; }
  std::string DfName;
  std::string dfName() const { return DfName; }
  int RunNumber;
  int runNumber() const { return RunNumber; }
};

const int nMixCases = 15;
const int nMixBin = 40;
const int nMixEvt = 5;

int mixingSlotPos[nMixCases][nMixBin]; // For a particular case -> in a particular mixing bin, It gives slot Position To Be Filled For Current Colllision
int collRoll[nMixCases][nMixBin];      // For a particular case -> in a particular mixing bin, It gives collRoll (rollNo of current collision in that bin)

std::vector<std::vector<int>> mixSlotPosList;

std::vector<StoredColl> myStoredMixingCollisions[nMixCases][nMixBin];

// Store in the processsing format
std::vector<myStoredTrack> PerColl_triggerTracks[nMixCases][nMixBin][nMixEvt + 1];
std::vector<myStoredTrack> PerColl_posTracks[nMixCases][nMixBin][nMixEvt + 1];
std::vector<myStoredTrack> PerColl_negTracks[nMixCases][nMixBin][nMixEvt + 1];
std::vector<myStoredTrack> PerColl_associatedTracks_0To2[nMixCases][nMixBin][nMixEvt + 1];
std::vector<myStoredTrack> PerColl_associatedTracks_2To4[nMixCases][nMixBin][nMixEvt + 1];
// End--Mixing Structures

struct hphicorrelation {
  // Hisogram redistry:
  HistogramRegistry SE_recoEvent{"SE_recoEvent", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry SE_recoTracks{"SE_recoTracks", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry SE_recoKaon{"SE_recoKaon", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry SE_recoPhi{"SE_recoPhi", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry SE_recoTrigger{"SE_recoTrigger", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry SE_recoAnalysis{"SE_recoAnalysis", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  HistogramRegistry ME_recoEvent{"ME_recoEvent", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry ME_recoTracks{"ME_recoTracks", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry ME_recoKaon{"ME_recoKaon", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry ME_recoPhi{"ME_recoPhi", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry ME_recoTrigger{"ME_recoTrigger", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry ME_recoAnalysis{"ME_recoAnalysis", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // PDG data base
  Service<o2::framework::O2DatabasePDG> pdgDB;

  // Configurables
  // Event Selection
  Configurable<float> cutZvertex{"cutZvertex", 10.0f, "Accepted z-vertex range (cm)"};

  // track Selection
  Configurable<float> cfgCutPt{"cfgCutPt", 0.15, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};

  // phi meson cuts
  Configurable<float> cfgPhiMassLow{"cfgPhiMassLow", 1.013, "Min Phi invarient Mass"};
  Configurable<float> cfgPhiMassUp{"cfgPhiMassUp", 1.026, "Max Phi invarient Mass"};
  Configurable<float> cfgLSBMassLow{"cfgLSBMassLow", 0.0995, "LSB Min Phi invarient Mass"};
  Configurable<float> cfgLSBMassUp{"cfgLSBMassUp", 1.005, "LSB Max Phi invarient Mass"};
  Configurable<float> cfgRSBMassLow{"cfgRSBMassLow", 1.040, "RSB Min Phi invarient Mass"};
  Configurable<float> cfgRSBMassUp{"cfgRSBMassUp", 1.060, "RSB Max Phi invarient Mass"};

  // Trigger Track Cuts
  Configurable<float> cfgTriggerPtLow{"cfgTriggerPtLow", 4.0, "Trigger Track Low Pt Cut"};
  Configurable<float> cfgTriggerPtHigh{"cfgTriggerPtHigh", 8.0, "Trigger Track High Pt Cut"};

  /// Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 5, "Number of events to mix"};
  // ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  // ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0., 1., 5., 10., 30., 50., 70., 100., 110.}, "Mixing bins - multiplicity"};
  // ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, -1.0f, 20.0f, 50., 70., 100., 110.}, "Mixing bins - multiplicity"};

  std::vector<double> CfgVtxBins{VARIABLE_WIDTH, -10.0, -8., -6., -4., -2., 0., 2., 4., 6., 8., 10.};
  std::vector<double> CfgMultBins{VARIABLE_WIDTH, -1.0f, 20.0, 50., 80., 101.0};

  void init(InitContext const&)
  {
    LOGF(info, "Starting init");
    // Axes
    AxisSpec Axis_vertexZ = {30, -15., 15., "vrtx_{Z} [cm]"};
    AxisSpec Axis_centFT0C = {1200, -10.0, 110.0, "centFT0C(percentile)"};
    AxisSpec Axis_Mult = {150, -1.0, 149.0};

    AxisSpec Axis_p = {200, 0.0f, 10.0f, "#it{p} (GeV/#it{c})"};
    AxisSpec Axis_pt = {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec Axis_tpcInnerParam = {200, 0.0f, 10.0f, "#it{p}_{tpcInnerParam} (GeV/#it{c})"};
    AxisSpec Axis_tofExpMom = {200, 0.0f, 10.0f, "#it{p}_{tofExpMom} (GeV/#it{c})"};

    AxisSpec Axis_eta = {40, -2, 2, "#eta"};
    AxisSpec Axis_phi = {110, -1, 10, "#phi (radians)"};

    AxisSpec Axis_rapidity = {100, -5, 5, "Rapidity (y)"};
    AxisSpec Axis_dcaXY = {80, -4.0, 4.0, "dcaXY"};
    AxisSpec Axis_dcaZ = {80, -4.0, 4.0, "dcaZ"};
    AxisSpec Axis_Sign = {10, -5, 5, "track.sign"};

    AxisSpec Axis_tpcSignal = {10010, -1, 1000, "tpcSignal"};
    AxisSpec Axis_tofBeta = {400, -2.0, 2.0, "tofBeta"};

    AxisSpec Axis_tpcNSigmaPi = {200, -10.0, 10.0, "n#sigma_{TPC}^{Pi}"};
    AxisSpec Axis_tofNSigmaPi = {200, -10.0, 10.0, "n#sigma_{TOF}^{Pi}"};
    AxisSpec Axis_tpcNSigmaKa = {200, -10.0, 10.0, "n#sigma_{TPC}^{Ka}"};
    AxisSpec Axis_tofNSigmaKa = {200, -10.0, 10.0, "n#sigma_{TOF}^{Ka}"};
    AxisSpec Axis_tpcNSigmaPr = {200, -10.0, 10.0, "n#sigma_{TPC}^{Pr}"};
    AxisSpec Axis_tofNSigmaPr = {200, -10.0, 10.0, "n#sigma_{TOF}^{Pr}"};
    AxisSpec Axis_tpcNSigmaEl = {200, -10.0, 10.0, "n#sigma_{TPC}^{El}"};
    AxisSpec Axis_tofNSigmaEl = {200, -10.0, 10.0, "n#sigma_{TOF}^{El}"};
    AxisSpec Axis_tpcNSigmaDe = {200, -10.0, 10.0, "n#sigma_{TPC}^{De}"};
    AxisSpec Axis_tofNSigmaDe = {200, -10.0, 10.0, "n#sigma_{TOF}^{De}"};

    AxisSpec Axis_00 = Axis_p;             // p              //0
    AxisSpec Axis_01 = Axis_pt;            // pt             //1
    AxisSpec Axis_02 = Axis_tpcInnerParam; // tpcInnerParam  //2
    AxisSpec Axis_03 = Axis_tofExpMom;     // tofExpMom      //3

    AxisSpec Axis_05 = Axis_tpcSignal; // tpcSignal      //5
    AxisSpec Axis_06 = Axis_tofBeta;   // Axis_tofBeta   //

    AxisSpec Axis_20 = Axis_tpcNSigmaPi; // tpcNSigmaPi    //5
    AxisSpec Axis_21 = Axis_tofNSigmaPi; // tofNSigmaPi
    AxisSpec Axis_22 = Axis_tpcNSigmaKa; // tpcNSigmaKa    //5
    AxisSpec Axis_23 = Axis_tofNSigmaKa; // tofNSigmaKa
    AxisSpec Axis_24 = Axis_tpcNSigmaPr; // tpcNSigmaPr    //5
    AxisSpec Axis_25 = Axis_tofNSigmaPr; // tofNSigmaPr
    AxisSpec Axis_26 = Axis_tpcNSigmaEl; // tpcNSigmaEl    //5
    AxisSpec Axis_27 = Axis_tofNSigmaEl; // tofNSigmaEl
    AxisSpec Axis_28 = Axis_tpcNSigmaDe; // tpcNSigmaDe    //5
    AxisSpec Axis_29 = Axis_tofNSigmaDe; // tofNSigmaDe

    AxisSpec Axis_PidTag = {34 * 2, -1.0, 33.0, "PID TAG"};
    AxisSpec Axis_paircharge = {6 * 2, -3, 3};
    AxisSpec Axis_PhiMass = {800, 0.99f, 1.07f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};

    AxisSpec Axis_dPhi = {80, -2.0f, 6.0f, "#Delta#phi"};
    AxisSpec Axis_dEta = {40, -2, 2, "#dEta"};

    std::vector<double> centBinning = {0., 20., 50., 80., 101.0};
    AxisSpec Axis_CentBins = {centBinning, "centBins"};
    AxisSpec Axis_vtxZbins = {14, -14, 14, "vtxZbins (Width = 2cm)"};

    LOG(info) << "DEBUG :: SIZE VtxBins  :: " << CfgVtxBins.size() - 2;
    LOG(info) << "DEBUG :: SIZE MultBins :: " << CfgMultBins.size() - 2;
    int BinCount = (CfgVtxBins.size() - 2) * (CfgMultBins.size() - 2);

    AxisSpec Axis_MixingBin = {(int(BinCount * 1.2) - (-1)) * 4 + 2, -1.25, int(BinCount * 1.2) + 0.25, "Mixing Bin"};
    AxisSpec Axis_EventCount = {101, -1, 100, "EventCount"};

    // Histograms
    // Event Selection
    SE_recoEvent.add("SE_hFullCollisionCount", "SE_hFullCollisionCount", {HistType::kTH1D, {{1, 0, 1}}});

    SE_recoEvent.add("SE_hCollisionCount", "SE_hCollisionCount", {HistType::kTH1D, {{1, 0, 1}}});
    SE_recoEvent.add("SE_hVertexXRec", "SE_hVertexXRec", {HistType::kTH1D, {{10000, -0.2, 0.2}}});
    SE_recoEvent.add("SE_hVertexYRec", "SE_hVertexYRec", {HistType::kTH1D, {{10000, -0.2, 0.2}}});
    SE_recoEvent.add("SE_hVertexZRec", "SE_Vertex distribution in Z;Z (cm);Counts", {HistType::kTH1F, {Axis_vertexZ}});
    SE_recoEvent.add("SE_hCentrality", "SE_hCentrality", {HistType::kTH1F, {Axis_centFT0C}});
    SE_recoEvent.add("SE_hCentrality_vs_vtxZ", "SE_hCentrality_vs_vtxZ", {HistType::kTH2F, {Axis_centFT0C, Axis_vertexZ}});

    SE_recoEvent.add("SE_hEvent_TrackSize", "SE_hEvent_TrackSize;track.size();Counts", kTH1F, {Axis_Mult});
    SE_recoEvent.add("SE_hEvent_nTrack", "SE_hEvent_nTrack;N_{tracks};Counts", kTH1F, {Axis_Mult});
    SE_recoEvent.add("SE_hEvent_nTrigger", "SE_hEvent_nTrigger;N_{Trigger};Counts", kTH1F, {Axis_Mult});
    SE_recoEvent.add("SE_hEvent_nPhi", "SE_hEvent_nPhi;N_{#phi};Counts", kTH1F, {Axis_Mult});
    SE_recoEvent.add("SE_hEvent_nPhiPhi", "SE_hEvent_nPhiPhi;N_{#phi#phi pairs};Counts", kTH1F, {Axis_Mult});

    SE_recoEvent.add("SE_hEvent_nPhi_0_2", "SE_hEvent_nPhi_0_2;N_{nPhi_0_2};Counts", kTH1F, {Axis_Mult});
    SE_recoEvent.add("SE_hEvent_nPhi_2_4", "SE_hEvent_nPhi_2_4;N_{nPhi_2_4};Counts", kTH1F, {Axis_Mult});
    SE_recoEvent.add("SE_hEvent_nPhi_4_8", "SE_hEvent_nPhi_4_8;N_{nPhi_4_8};Counts", kTH1F, {Axis_Mult});
    SE_recoEvent.add("SE_hEvent_nPhi_8_i", "SE_hEvent_nPhi_8_i;N_{nPhi_8_i};Counts", kTH1F, {Axis_Mult});
    SE_recoEvent.add("SE_hEvent_nLSB_0_2", "SE_hEvent_nLSB_0_2;N_{nLSB_0_2};Counts", kTH1F, {Axis_Mult});
    SE_recoEvent.add("SE_hEvent_nLSB_2_4", "SE_hEvent_nLSB_2_4;N_{nLSB_2_4};Counts", kTH1F, {Axis_Mult});
    SE_recoEvent.add("SE_hEvent_nLSB_4_8", "SE_hEvent_nLSB_4_8;N_{nLSB_4_8};Counts", kTH1F, {Axis_Mult});
    SE_recoEvent.add("SE_hEvent_nLSB_8_i", "SE_hEvent_nLSB_8_i;N_{nLSB_8_i};Counts", kTH1F, {Axis_Mult});
    SE_recoEvent.add("SE_hEvent_nRSB_0_2", "SE_hEvent_nRSB_0_2;N_{nRSB_0_2};Counts", kTH1F, {Axis_Mult});
    SE_recoEvent.add("SE_hEvent_nRSB_2_4", "SE_hEvent_nRSB_2_4;N_{nRSB_2_4};Counts", kTH1F, {Axis_Mult});
    SE_recoEvent.add("SE_hEvent_nRSB_4_8", "SE_hEvent_nRSB_4_8;N_{nRSB_4_8};Counts", kTH1F, {Axis_Mult});
    SE_recoEvent.add("SE_hEvent_nRSB_8_i", "SE_hEvent_nRSB_8_i;N_{nRSB_8_i};Counts", kTH1F, {Axis_Mult});

    SE_recoEvent.add("SE_hEvent_nLeadPhi", "SE_hEvent_nLeadPhi;N_{nLeadPhi};Counts", kTH1F, {Axis_Mult});
    SE_recoEvent.add("SE_hEvent_nAssoPhi_0_2", "SE_hEvent_nAssoPhi_0_2;N_{nAssoPhi_0_2};Counts", kTH1F, {Axis_Mult});
    SE_recoEvent.add("SE_hEvent_nAssoPhi_2_4", "SE_hEvent_nAssoPhi_2_4;N_{nAssoPhi_2_4};Counts", kTH1F, {Axis_Mult});
    SE_recoEvent.add("SE_hEvent_nAssoHad_0_2", "SE_hEvent_nAssoHad_0_2;N_{nAssoHad_0_2};Counts", kTH1F, {Axis_Mult});
    SE_recoEvent.add("SE_hEvent_nAssoHad_2_4", "SE_hEvent_nAssoHad_2_4;N_{nAssoHad_2_4};Counts", kTH1F, {Axis_Mult});
    //

    // Tracks reconstruction
    // FullTrack
    SE_recoTracks.add("SE_htracks_00_1_0_FullTrack_P", "SE_hTracks_SelectedTrack_P", {HistType::kTH1F, {Axis_p}});
    SE_recoTracks.add("SE_htracks_00_1_0_FullTrack_tpcInnerParam", "SE_hTracks_SelectedTrack_tpcInnerParam", {HistType::kTH1F, {Axis_tpcInnerParam}});
    SE_recoTracks.add("SE_htracks_00_1_0_FullTrack_tofExpMom", "SE_hTracks_SelectedTrack_tofExpMom", {HistType::kTH1F, {Axis_tofExpMom}});
    SE_recoTracks.add("SE_hTracks_00_1_1_FullTrack_Pt", "SE_hTracks_FullTrack_Pt", {HistType::kTH1F, {Axis_pt}});
    SE_recoTracks.add("SE_hTracks_00_1_2_FullTrack_Eta", "SE_hTracks_FullTrack_Eta", {HistType::kTH1F, {Axis_eta}});
    SE_recoTracks.add("SE_hTracks_00_1_3_FullTrack_Phi", "SE_hTracks_FullTrack_Phi", {HistType::kTH1F, {Axis_phi}});
    SE_recoTracks.add("SE_hTracks_00_1_4_FullTrack_DcaXY", "SE_hTracks_FullTrack_DcaXY", {HistType::kTH1F, {Axis_dcaXY}});
    SE_recoTracks.add("SE_hTracks_00_1_5_FullTrack_DcaZ", "SE_hTracks_FullTrack_DcaZ", {HistType::kTH1F, {Axis_dcaZ}});
    SE_recoTracks.add("SE_hTracks_00_1_6_FullTrack_Sign", "SE_hTracks_FullTrack_Sign", {HistType::kTH1D, {Axis_Sign}});

    // DcaXY
    SE_recoTracks.add("SE_htracks_00_1_7_FullTrack_00_DcaXY", "SE_htracks_FullTrack_00_DcaXY", kTH2F, {Axis_00, Axis_dcaXY});
    SE_recoTracks.add("SE_htracks_00_1_7_FullTrack_01_DcaXY", "SE_htracks_FullTrack_01_DcaXY", kTH2F, {Axis_01, Axis_dcaXY});
    SE_recoTracks.add("SE_htracks_00_1_7_FullTrack_02_DcaXY", "SE_htracks_FullTrack_02_DcaXY", kTH2F, {Axis_02, Axis_dcaXY});
    SE_recoTracks.add("SE_htracks_00_1_7_FullTrack_03_DcaXY", "SE_htracks_FullTrack_03_DcaXY", kTH2F, {Axis_03, Axis_dcaXY});

    // DcaZ
    SE_recoTracks.add("SE_htracks_00_1_7_FullTrack_00_DcaZ", "SE_htracks_FullTrack_00_DcaZ", kTH2F, {Axis_00, Axis_dcaZ});
    SE_recoTracks.add("SE_htracks_00_1_7_FullTrack_01_DcaZ", "SE_htracks_FullTrack_01_DcaZ", kTH2F, {Axis_01, Axis_dcaZ});
    SE_recoTracks.add("SE_htracks_00_1_7_FullTrack_02_DcaZ", "SE_htracks_FullTrack_02_DcaZ", kTH2F, {Axis_02, Axis_dcaZ});
    SE_recoTracks.add("SE_htracks_00_1_7_FullTrack_03_DcaZ", "SE_htracks_FullTrack_03_DcaZ", kTH2F, {Axis_03, Axis_dcaZ});

    // momemtum
    SE_recoTracks.add("SE_hTracks_00_2_1_FullTrack_Axis_00_01", "SE_htrack_Axis_00_01", kTH2F, {Axis_00, Axis_01}); // Axis_pt            ;
    SE_recoTracks.add("SE_hTracks_00_2_2_FullTrack_Axis_00_02", "SE_htrack_Axis_00_02", kTH2F, {Axis_00, Axis_02}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_hTracks_00_2_3_FullTrack_Axis_00_03", "SE_htrack_Axis_00_03", kTH2F, {Axis_00, Axis_03}); // Axis_tofExpMom     ;

    // tpcSignal
    SE_recoTracks.add("SE_hTracks_00_3_1_FullTrack_Axis_00_05", "SE_htrack_Axis_00_05", kTH2F, {Axis_00, Axis_05}); // Axis_p             ;
    SE_recoTracks.add("SE_hTracks_00_3_2_FullTrack_Axis_02_05", "SE_htrack_Axis_02_05", kTH2F, {Axis_02, Axis_05}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_hTracks_00_3_3_FullTrack_Axis_03_05", "SE_htrack_Axis_03_05", kTH2F, {Axis_03, Axis_05}); // Axis_tofExpMom     ;

    // tofBeta
    SE_recoTracks.add("SE_hTracks_00_4_1_FullTrack_Axis_00_06", "SE_htrack_Axis_00_06", kTH2F, {Axis_00, Axis_06}); // Axis_p             ;
    SE_recoTracks.add("SE_hTracks_00_4_2_FullTrack_Axis_02_06", "SE_htrack_Axis_02_06", kTH2F, {Axis_02, Axis_06}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_hTracks_00_4_3_FullTrack_Axis_03_06", "SE_htrack_Axis_03_06", kTH2F, {Axis_03, Axis_06}); // Axis_tofExpMom     ;

    // Look at Pion
    SE_recoTracks.add("SE_hTracks_00_5_1_FullTrack_Axis_00_20", "SE_htrack_Axis_00_20", kTH2F, {Axis_00, Axis_20}); // Axis_p             ;
    SE_recoTracks.add("SE_hTracks_00_5_2_FullTrack_Axis_01_20", "SE_htrack_Axis_01_20", kTH2F, {Axis_01, Axis_20}); // Axis_pt            ;
    SE_recoTracks.add("SE_hTracks_00_5_3_FullTrack_Axis_02_20", "SE_htrack_Axis_02_20", kTH2F, {Axis_02, Axis_20}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_hTracks_00_5_4_FullTrack_Axis_03_20", "SE_htrack_Axis_03_20", kTH2F, {Axis_03, Axis_20}); // Axis_tofExpMom     ;
    SE_recoTracks.add("SE_hTracks_00_5_5_FullTrack_Axis_00_21", "SE_htrack_Axis_00_21", kTH2F, {Axis_00, Axis_21}); // Axis_p             ;
    SE_recoTracks.add("SE_hTracks_00_5_6_FullTrack_Axis_01_21", "SE_htrack_Axis_01_21", kTH2F, {Axis_01, Axis_21}); // Axis_pt            ;
    SE_recoTracks.add("SE_hTracks_00_5_7_FullTrack_Axis_02_21", "SE_htrack_Axis_02_21", kTH2F, {Axis_02, Axis_21}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_hTracks_00_5_8_FullTrack_Axis_03_21", "SE_htrack_Axis_03_21", kTH2F, {Axis_03, Axis_21}); // Axis_tofExpMom     ;
    SE_recoTracks.add("SE_hTracks_00_5_9_FullTrack_Axis_20_21", "SE_htrack_Axis_20_21", kTH2F, {Axis_20, Axis_21}); // Axis_tpcInnerParam ;
    // Pion
    // Look at Kaon
    SE_recoTracks.add("SE_hTracks_00_6_1_FullTrack_Axis_00_22", "SE_htrack_Axis_00_22", kTH2F, {Axis_00, Axis_22}); // Axis_p             ;
    SE_recoTracks.add("SE_hTracks_00_6_2_FullTrack_Axis_01_22", "SE_htrack_Axis_01_22", kTH2F, {Axis_01, Axis_22}); // Axis_pt            ;
    SE_recoTracks.add("SE_hTracks_00_6_3_FullTrack_Axis_02_22", "SE_htrack_Axis_02_22", kTH2F, {Axis_02, Axis_22}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_hTracks_00_6_4_FullTrack_Axis_03_22", "SE_htrack_Axis_03_22", kTH2F, {Axis_03, Axis_22}); // Axis_tofExpMom     ;
    SE_recoTracks.add("SE_hTracks_00_6_5_FullTrack_Axis_00_23", "SE_htrack_Axis_00_23", kTH2F, {Axis_00, Axis_23}); // Axis_p             ;
    SE_recoTracks.add("SE_hTracks_00_6_6_FullTrack_Axis_01_23", "SE_htrack_Axis_01_23", kTH2F, {Axis_01, Axis_23}); // Axis_pt            ;
    SE_recoTracks.add("SE_hTracks_00_6_7_FullTrack_Axis_02_23", "SE_htrack_Axis_02_23", kTH2F, {Axis_02, Axis_23}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_hTracks_00_6_8_FullTrack_Axis_03_23", "SE_htrack_Axis_03_23", kTH2F, {Axis_03, Axis_23}); // Axis_tofExpMom     ;
    SE_recoTracks.add("SE_hTracks_00_6_9_FullTrack_Axis_22_23", "SE_htrack_Axis_22_23", kTH2F, {Axis_22, Axis_23}); // Axis_tpcInnerParam ;
    // Kaon
    // Look at Proton
    SE_recoTracks.add("SE_hTracks_00_7_1_FullTrack_Axis_00_24", "SE_htrack_Axis_00_24", kTH2F, {Axis_00, Axis_24}); // Axis_p             ;
    SE_recoTracks.add("SE_hTracks_00_7_2_FullTrack_Axis_01_24", "SE_htrack_Axis_01_24", kTH2F, {Axis_01, Axis_24}); // Axis_pt            ;
    SE_recoTracks.add("SE_hTracks_00_7_3_FullTrack_Axis_02_24", "SE_htrack_Axis_02_24", kTH2F, {Axis_02, Axis_24}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_hTracks_00_7_4_FullTrack_Axis_03_24", "SE_htrack_Axis_03_24", kTH2F, {Axis_03, Axis_24}); // Axis_tofExpMom     ;
    SE_recoTracks.add("SE_hTracks_00_7_5_FullTrack_Axis_00_25", "SE_htrack_Axis_00_25", kTH2F, {Axis_00, Axis_25}); // Axis_p             ;
    SE_recoTracks.add("SE_hTracks_00_7_6_FullTrack_Axis_01_25", "SE_htrack_Axis_01_25", kTH2F, {Axis_01, Axis_25}); // Axis_pt            ;
    SE_recoTracks.add("SE_hTracks_00_7_7_FullTrack_Axis_02_25", "SE_htrack_Axis_02_25", kTH2F, {Axis_02, Axis_25}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_hTracks_00_7_8_FullTrack_Axis_03_25", "SE_htrack_Axis_03_25", kTH2F, {Axis_03, Axis_25}); // Axis_tofExpMom     ;
    SE_recoTracks.add("SE_hTracks_00_7_9_FullTrack_Axis_24_25", "SE_htrack_Axis_24_25", kTH2F, {Axis_24, Axis_25}); // Axis_tpcInnerParam ;
    // Proton
    // Look at Electron
    SE_recoTracks.add("SE_hTracks_00_8_1_FullTrack_Axis_00_26", "SE_htrack_Axis_00_26", kTH2F, {Axis_00, Axis_26}); // Axis_p             ;
    SE_recoTracks.add("SE_hTracks_00_8_2_FullTrack_Axis_01_26", "SE_htrack_Axis_01_26", kTH2F, {Axis_01, Axis_26}); // Axis_pt            ;
    SE_recoTracks.add("SE_hTracks_00_8_3_FullTrack_Axis_02_26", "SE_htrack_Axis_02_26", kTH2F, {Axis_02, Axis_26}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_hTracks_00_8_4_FullTrack_Axis_03_26", "SE_htrack_Axis_03_26", kTH2F, {Axis_03, Axis_26}); // Axis_tofExpMom     ;
    SE_recoTracks.add("SE_hTracks_00_8_5_FullTrack_Axis_00_27", "SE_htrack_Axis_00_27", kTH2F, {Axis_00, Axis_27}); // Axis_p             ;
    SE_recoTracks.add("SE_hTracks_00_8_6_FullTrack_Axis_01_27", "SE_htrack_Axis_01_27", kTH2F, {Axis_01, Axis_27}); // Axis_pt            ;
    SE_recoTracks.add("SE_hTracks_00_8_7_FullTrack_Axis_02_27", "SE_htrack_Axis_02_27", kTH2F, {Axis_02, Axis_27}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_hTracks_00_8_8_FullTrack_Axis_03_27", "SE_htrack_Axis_03_27", kTH2F, {Axis_03, Axis_27}); // Axis_tofExpMom     ;
    SE_recoTracks.add("SE_hTracks_00_8_9_FullTrack_Axis_26_27", "SE_htrack_Axis_26_27", kTH2F, {Axis_26, Axis_27}); // Axis_tpcInnerParam ;
    // Electron
    // Look at Deuteron
    SE_recoTracks.add("SE_hTracks_00_9_1_FullTrack_Axis_00_28", "SE_htrack_Axis_00_28", kTH2F, {Axis_00, Axis_28}); // Axis_p             ;
    SE_recoTracks.add("SE_hTracks_00_9_2_FullTrack_Axis_01_28", "SE_htrack_Axis_01_28", kTH2F, {Axis_01, Axis_28}); // Axis_pt            ;
    SE_recoTracks.add("SE_hTracks_00_9_3_FullTrack_Axis_02_28", "SE_htrack_Axis_02_28", kTH2F, {Axis_02, Axis_28}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_hTracks_00_9_4_FullTrack_Axis_03_28", "SE_htrack_Axis_03_28", kTH2F, {Axis_03, Axis_28}); // Axis_tofExpMom     ;
    SE_recoTracks.add("SE_hTracks_00_9_5_FullTrack_Axis_00_29", "SE_htrack_Axis_00_29", kTH2F, {Axis_00, Axis_29}); // Axis_p             ;
    SE_recoTracks.add("SE_hTracks_00_9_6_FullTrack_Axis_01_29", "SE_htrack_Axis_01_29", kTH2F, {Axis_01, Axis_29}); // Axis_pt            ;
    SE_recoTracks.add("SE_hTracks_00_9_7_FullTrack_Axis_02_29", "SE_htrack_Axis_02_29", kTH2F, {Axis_02, Axis_29}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_hTracks_00_9_8_FullTrack_Axis_03_29", "SE_htrack_Axis_03_29", kTH2F, {Axis_03, Axis_29}); // Axis_tofExpMom     ;
    SE_recoTracks.add("SE_hTracks_00_9_9_FullTrack_Axis_28_29", "SE_htrack_Axis_28_29", kTH2F, {Axis_28, Axis_29}); // Axis_tpcInnerParam ;
    // Deuteron
    // FullTrack

    // SelectedTrack
    SE_recoTracks.add("SE_htracks_11_1_0_SelectedTrack_P", "SE_hTracks_SelectedTrack_P", {HistType::kTH1F, {Axis_p}});
    SE_recoTracks.add("SE_htracks_11_1_0_SelectedTrack_tpcInnerParam", "SE_hTracks_SelectedTrack_tpcInnerParam", {HistType::kTH1F, {Axis_tpcInnerParam}});
    SE_recoTracks.add("SE_htracks_11_1_0_SelectedTrack_tofExpMom", "SE_hTracks_SelectedTrack_tofExpMom", {HistType::kTH1F, {Axis_tofExpMom}});

    SE_recoTracks.add("SE_htracks_11_1_1_SelectedTrack_Pt", "SE_hTracks_SelectedTrack_Pt", {HistType::kTH1F, {Axis_pt}});
    SE_recoTracks.add("SE_htracks_11_1_2_SelectedTrack_Eta", "SE_hTracks_SelectedTrack_Eta", {HistType::kTH1F, {Axis_eta}});
    SE_recoTracks.add("SE_htracks_11_1_3_SelectedTrack_Phi", "SE_hTracks_SelectedTrack_Phi", {HistType::kTH1F, {Axis_phi}});
    SE_recoTracks.add("SE_htracks_11_1_4_SelectedTrack_DcaXY", "SE_hTracks_SelectedTrack_DcaXY", {HistType::kTH1F, {Axis_dcaXY}});
    SE_recoTracks.add("SE_htracks_11_1_5_SelectedTrack_DcaZ", "SE_hTracks_SelectedTrack_DcaZ", {HistType::kTH1F, {Axis_dcaZ}});
    SE_recoTracks.add("SE_htracks_11_1_6_SelectedTrack_Sign", "SE_hTracks_SelectedTrack_Sign", {HistType::kTH1D, {Axis_Sign}});

    // DcaXY
    SE_recoTracks.add("SE_htracks_11_1_7_SelectedTrack_00_DcaXY", "SE_htracks_SelectedTrack_00_DcaXY", kTH2F, {Axis_00, Axis_dcaXY});
    SE_recoTracks.add("SE_htracks_11_1_7_SelectedTrack_01_DcaXY", "SE_htracks_SelectedTrack_01_DcaXY", kTH2F, {Axis_01, Axis_dcaXY});
    SE_recoTracks.add("SE_htracks_11_1_7_SelectedTrack_02_DcaXY", "SE_htracks_SelectedTrack_02_DcaXY", kTH2F, {Axis_02, Axis_dcaXY});
    SE_recoTracks.add("SE_htracks_11_1_7_SelectedTrack_03_DcaXY", "SE_htracks_SelectedTrack_03_DcaXY", kTH2F, {Axis_03, Axis_dcaXY});

    // DcaZ
    SE_recoTracks.add("SE_htracks_11_1_7_SelectedTrack_00_DcaZ", "SE_htracks_SelectedTrack_00_DcaZ", kTH2F, {Axis_00, Axis_dcaZ});
    SE_recoTracks.add("SE_htracks_11_1_7_SelectedTrack_01_DcaZ", "SE_htracks_SelectedTrack_01_DcaZ", kTH2F, {Axis_01, Axis_dcaZ});
    SE_recoTracks.add("SE_htracks_11_1_7_SelectedTrack_02_DcaZ", "SE_htracks_SelectedTrack_02_DcaZ", kTH2F, {Axis_02, Axis_dcaZ});
    SE_recoTracks.add("SE_htracks_11_1_7_SelectedTrack_03_DcaZ", "SE_htracks_SelectedTrack_03_DcaZ", kTH2F, {Axis_03, Axis_dcaZ});

    // momemtum
    SE_recoTracks.add("SE_htracks_11_2_1_SelectedTrack_Axis_00_01", "SE_htrack_Axis_00_01", kTH2F, {Axis_00, Axis_01}); // Axis_pt            ;
    SE_recoTracks.add("SE_htracks_11_2_2_SelectedTrack_Axis_00_02", "SE_htrack_Axis_00_02", kTH2F, {Axis_00, Axis_02}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_htracks_11_2_3_SelectedTrack_Axis_00_03", "SE_htrack_Axis_00_03", kTH2F, {Axis_00, Axis_03}); // Axis_tofExpMom     ;

    // tpcSignal
    SE_recoTracks.add("SE_htracks_11_3_1_SelectedTrack_Axis_00_05", "SE_htrack_Axis_00_05", kTH2F, {Axis_00, Axis_05}); // Axis_p             ;
    SE_recoTracks.add("SE_htracks_11_3_2_SelectedTrack_Axis_02_05", "SE_htrack_Axis_02_05", kTH2F, {Axis_02, Axis_05}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_htracks_11_3_3_SelectedTrack_Axis_03_05", "SE_htrack_Axis_03_05", kTH2F, {Axis_03, Axis_05}); // Axis_tofExpMom     ;

    // tofBeta
    SE_recoTracks.add("SE_htracks_11_4_1_SelectedTrack_Axis_00_06", "SE_htrack_Axis_00_06", kTH2F, {Axis_00, Axis_06}); // Axis_p             ;
    SE_recoTracks.add("SE_htracks_11_4_2_SelectedTrack_Axis_02_06", "SE_htrack_Axis_02_06", kTH2F, {Axis_02, Axis_06}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_htracks_11_4_3_SelectedTrack_Axis_03_06", "SE_htrack_Axis_03_06", kTH2F, {Axis_03, Axis_06}); // Axis_tofExpMom     ;

    // Look at Pion
    SE_recoTracks.add("SE_htracks_11_5_1_SelectedTrack_Axis_00_20", "SE_htrack_Axis_00_20", kTH2F, {Axis_00, Axis_20}); // Axis_p             ;
    SE_recoTracks.add("SE_htracks_11_5_2_SelectedTrack_Axis_01_20", "SE_htrack_Axis_01_20", kTH2F, {Axis_01, Axis_20}); // Axis_pt            ;
    SE_recoTracks.add("SE_htracks_11_5_3_SelectedTrack_Axis_02_20", "SE_htrack_Axis_02_20", kTH2F, {Axis_02, Axis_20}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_htracks_11_5_4_SelectedTrack_Axis_03_20", "SE_htrack_Axis_03_20", kTH2F, {Axis_03, Axis_20}); // Axis_tofExpMom     ;
    SE_recoTracks.add("SE_htracks_11_5_5_SelectedTrack_Axis_00_21", "SE_htrack_Axis_00_21", kTH2F, {Axis_00, Axis_21}); // Axis_p             ;
    SE_recoTracks.add("SE_htracks_11_5_6_SelectedTrack_Axis_01_21", "SE_htrack_Axis_01_21", kTH2F, {Axis_01, Axis_21}); // Axis_pt            ;
    SE_recoTracks.add("SE_htracks_11_5_7_SelectedTrack_Axis_02_21", "SE_htrack_Axis_02_21", kTH2F, {Axis_02, Axis_21}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_htracks_11_5_8_SelectedTrack_Axis_03_21", "SE_htrack_Axis_03_21", kTH2F, {Axis_03, Axis_21}); // Axis_tofExpMom     ;
    SE_recoTracks.add("SE_htracks_11_5_9_SelectedTrack_Axis_20_21", "SE_htrack_Axis_20_21", kTH2F, {Axis_20, Axis_21}); // Axis_tpcInnerParam ;
    // Pion
    // Look at Kaon
    SE_recoTracks.add("SE_htracks_11_6_1_SelectedTrack_Axis_00_22", "SE_htrack_Axis_00_22", kTH2F, {Axis_00, Axis_22}); // Axis_p             ;
    SE_recoTracks.add("SE_htracks_11_6_2_SelectedTrack_Axis_01_22", "SE_htrack_Axis_01_22", kTH2F, {Axis_01, Axis_22}); // Axis_pt            ;
    SE_recoTracks.add("SE_htracks_11_6_3_SelectedTrack_Axis_02_22", "SE_htrack_Axis_02_22", kTH2F, {Axis_02, Axis_22}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_htracks_11_6_4_SelectedTrack_Axis_03_22", "SE_htrack_Axis_03_22", kTH2F, {Axis_03, Axis_22}); // Axis_tofExpMom     ;
    SE_recoTracks.add("SE_htracks_11_6_5_SelectedTrack_Axis_00_23", "SE_htrack_Axis_00_23", kTH2F, {Axis_00, Axis_23}); // Axis_p             ;
    SE_recoTracks.add("SE_htracks_11_6_6_SelectedTrack_Axis_01_23", "SE_htrack_Axis_01_23", kTH2F, {Axis_01, Axis_23}); // Axis_pt            ;
    SE_recoTracks.add("SE_htracks_11_6_7_SelectedTrack_Axis_02_23", "SE_htrack_Axis_02_23", kTH2F, {Axis_02, Axis_23}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_htracks_11_6_8_SelectedTrack_Axis_03_23", "SE_htrack_Axis_03_23", kTH2F, {Axis_03, Axis_23}); // Axis_tofExpMom     ;
    SE_recoTracks.add("SE_htracks_11_6_9_SelectedTrack_Axis_22_23", "SE_htrack_Axis_22_23", kTH2F, {Axis_22, Axis_23}); // Axis_tpcInnerParam ;
    // Kaon
    // Look at Proton
    SE_recoTracks.add("SE_htracks_11_7_1_SelectedTrack_Axis_00_24", "SE_htrack_Axis_00_24", kTH2F, {Axis_00, Axis_24}); // Axis_p             ;
    SE_recoTracks.add("SE_htracks_11_7_2_SelectedTrack_Axis_01_24", "SE_htrack_Axis_01_24", kTH2F, {Axis_01, Axis_24}); // Axis_pt            ;
    SE_recoTracks.add("SE_htracks_11_7_3_SelectedTrack_Axis_02_24", "SE_htrack_Axis_02_24", kTH2F, {Axis_02, Axis_24}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_htracks_11_7_4_SelectedTrack_Axis_03_24", "SE_htrack_Axis_03_24", kTH2F, {Axis_03, Axis_24}); // Axis_tofExpMom     ;
    SE_recoTracks.add("SE_htracks_11_7_5_SelectedTrack_Axis_00_25", "SE_htrack_Axis_00_25", kTH2F, {Axis_00, Axis_25}); // Axis_p             ;
    SE_recoTracks.add("SE_htracks_11_7_6_SelectedTrack_Axis_01_25", "SE_htrack_Axis_01_25", kTH2F, {Axis_01, Axis_25}); // Axis_pt            ;
    SE_recoTracks.add("SE_htracks_11_7_7_SelectedTrack_Axis_02_25", "SE_htrack_Axis_02_25", kTH2F, {Axis_02, Axis_25}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_htracks_11_7_8_SelectedTrack_Axis_03_25", "SE_htrack_Axis_03_25", kTH2F, {Axis_03, Axis_25}); // Axis_tofExpMom     ;
    SE_recoTracks.add("SE_htracks_11_7_9_SelectedTrack_Axis_24_25", "SE_htrack_Axis_24_25", kTH2F, {Axis_24, Axis_25}); // Axis_tpcInnerParam ;
    // Proton
    // Look at Electron
    SE_recoTracks.add("SE_htracks_11_8_1_SelectedTrack_Axis_00_26", "SE_htrack_Axis_00_26", kTH2F, {Axis_00, Axis_26}); // Axis_p             ;
    SE_recoTracks.add("SE_htracks_11_8_2_SelectedTrack_Axis_01_26", "SE_htrack_Axis_01_26", kTH2F, {Axis_01, Axis_26}); // Axis_pt            ;
    SE_recoTracks.add("SE_htracks_11_8_3_SelectedTrack_Axis_02_26", "SE_htrack_Axis_02_26", kTH2F, {Axis_02, Axis_26}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_htracks_11_8_4_SelectedTrack_Axis_03_26", "SE_htrack_Axis_03_26", kTH2F, {Axis_03, Axis_26}); // Axis_tofExpMom     ;
    SE_recoTracks.add("SE_htracks_11_8_5_SelectedTrack_Axis_00_27", "SE_htrack_Axis_00_27", kTH2F, {Axis_00, Axis_27}); // Axis_p             ;
    SE_recoTracks.add("SE_htracks_11_8_6_SelectedTrack_Axis_01_27", "SE_htrack_Axis_01_27", kTH2F, {Axis_01, Axis_27}); // Axis_pt            ;
    SE_recoTracks.add("SE_htracks_11_8_7_SelectedTrack_Axis_02_27", "SE_htrack_Axis_02_27", kTH2F, {Axis_02, Axis_27}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_htracks_11_8_8_SelectedTrack_Axis_03_27", "SE_htrack_Axis_03_27", kTH2F, {Axis_03, Axis_27}); // Axis_tofExpMom     ;
    SE_recoTracks.add("SE_htracks_11_8_9_SelectedTrack_Axis_26_27", "SE_htrack_Axis_26_27", kTH2F, {Axis_26, Axis_27}); // Axis_tpcInnerParam ;
    // Electron
    // Look at Deuteron
    SE_recoTracks.add("SE_htracks_11_9_1_SelectedTrack_Axis_00_28", "SE_htrack_Axis_00_28", kTH2F, {Axis_00, Axis_28}); // Axis_p             ;
    SE_recoTracks.add("SE_htracks_11_9_2_SelectedTrack_Axis_01_28", "SE_htrack_Axis_01_28", kTH2F, {Axis_01, Axis_28}); // Axis_pt            ;
    SE_recoTracks.add("SE_htracks_11_9_3_SelectedTrack_Axis_02_28", "SE_htrack_Axis_02_28", kTH2F, {Axis_02, Axis_28}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_htracks_11_9_4_SelectedTrack_Axis_03_28", "SE_htrack_Axis_03_28", kTH2F, {Axis_03, Axis_28}); // Axis_tofExpMom     ;
    SE_recoTracks.add("SE_htracks_11_9_5_SelectedTrack_Axis_00_29", "SE_htrack_Axis_00_29", kTH2F, {Axis_00, Axis_29}); // Axis_p             ;
    SE_recoTracks.add("SE_htracks_11_9_6_SelectedTrack_Axis_01_29", "SE_htrack_Axis_01_29", kTH2F, {Axis_01, Axis_29}); // Axis_pt            ;
    SE_recoTracks.add("SE_htracks_11_9_7_SelectedTrack_Axis_02_29", "SE_htrack_Axis_02_29", kTH2F, {Axis_02, Axis_29}); // Axis_tpcInnerParam ;
    SE_recoTracks.add("SE_htracks_11_9_8_SelectedTrack_Axis_03_29", "SE_htrack_Axis_03_29", kTH2F, {Axis_03, Axis_29}); // Axis_tofExpMom     ;
    SE_recoTracks.add("SE_htracks_11_9_9_SelectedTrack_Axis_28_29", "SE_htrack_Axis_28_29", kTH2F, {Axis_28, Axis_29}); // Axis_tpcInnerParam ;
    // Deuteron
    // SelectedTrack
    //

    // Kaon identification
    // tpcSignal
    SE_recoKaon.add("SE_hKaon_12_Ka_Id0_1_Axis_00_05", "SE_hKaon_12_Ka_Id0_1_Axis_00_05", kTH2F, {Axis_00, Axis_05}); // track.p            (),track.tpcSignal()) ;//Axis_p             ;
    SE_recoKaon.add("SE_hKaon_12_Ka_Id0_1_Axis_02_05", "SE_hKaon_12_Ka_Id0_1_Axis_02_05", kTH2F, {Axis_02, Axis_05}); // track.tpcInnerParam(),track.tpcSignal()) ;//Axis_tpcInnerParam ;
    SE_recoKaon.add("SE_hKaon_12_Ka_Id0_1_Axis_03_05", "SE_hKaon_12_Ka_Id0_1_Axis_03_05", kTH2F, {Axis_03, Axis_05}); // track.tofExpMom    (),track.tpcSignal()) ;//Axis_tofExpMom     ;
    // tofBeta
    SE_recoKaon.add("SE_hKaon_12_Ka_Id0_2_Axis_00_06", "SE_hKaon_12_Ka_Id0_2_Axis_00_06", kTH2F, {Axis_00, Axis_06}); // track.p            (),track.beta()) ;//Axis_p             ;
    SE_recoKaon.add("SE_hKaon_12_Ka_Id0_2_Axis_02_06", "SE_hKaon_12_Ka_Id0_2_Axis_02_06", kTH2F, {Axis_02, Axis_06}); // track.tpcInnerParam(),track.beta()) ;//Axis_tpcInnerParam ;
    SE_recoKaon.add("SE_hKaon_12_Ka_Id0_2_Axis_03_06", "SE_hKaon_12_Ka_Id0_2_Axis_03_06", kTH2F, {Axis_03, Axis_06}); // track.tofExpMom    (),track.beta()) ;//Axis_tofExpMom     ;
    // Kaon
    SE_recoKaon.add("SE_hKaon_12_Ka_Id0_3_1_Axis_00_22", "SE_hKaon_12_Ka_Id0_3_1_Axis_00_22", kTH2F, {Axis_00, Axis_22}); // track.p            (),track.tpcNSigmaKa()) ;//Axis_p             ;
    SE_recoKaon.add("SE_hKaon_12_Ka_Id0_3_2_Axis_01_22", "SE_hKaon_12_Ka_Id0_3_2_Axis_01_22", kTH2F, {Axis_01, Axis_22}); // track.pt           (),track.tpcNSigmaKa()) ;//Axis_pt            ;
    SE_recoKaon.add("SE_hKaon_12_Ka_Id0_3_3_Axis_02_22", "SE_hKaon_12_Ka_Id0_3_3_Axis_02_22", kTH2F, {Axis_02, Axis_22}); // track.tpcInnerParam(),track.tpcNSigmaKa()) ;//Axis_tpcInnerParam ;
    SE_recoKaon.add("SE_hKaon_12_Ka_Id0_3_4_Axis_03_22", "SE_hKaon_12_Ka_Id0_3_4_Axis_03_22", kTH2F, {Axis_03, Axis_22}); // track.tofExpMom    (),track.tpcNSigmaKa()) ;//Axis_tofExpMom     ;
    SE_recoKaon.add("SE_hKaon_12_Ka_Id0_3_5_Axis_00_23", "SE_hKaon_12_Ka_Id0_3_5_Axis_00_23", kTH2F, {Axis_00, Axis_23}); // track.p            (),track.tofNSigmaKa()) ;//Axis_p             ;
    SE_recoKaon.add("SE_hKaon_12_Ka_Id0_3_6_Axis_01_23", "SE_hKaon_12_Ka_Id0_3_6_Axis_01_23", kTH2F, {Axis_01, Axis_23}); // track.pt           (),track.tofNSigmaKa()) ;//Axis_pt            ;
    SE_recoKaon.add("SE_hKaon_12_Ka_Id0_3_7_Axis_02_23", "SE_hKaon_12_Ka_Id0_3_7_Axis_02_23", kTH2F, {Axis_02, Axis_23}); // track.tpcInnerParam(),track.tofNSigmaKa()) ;//Axis_tpcInnerParam ;
    SE_recoKaon.add("SE_hKaon_12_Ka_Id0_3_8_Axis_03_23", "SE_hKaon_12_Ka_Id0_3_8_Axis_03_23", kTH2F, {Axis_03, Axis_23}); // track.tofExpMom    (),track.tofNSigmaKa()) ;//Axis_tofExpMom     ;
    SE_recoKaon.add("SE_hKaon_12_Ka_Id0_3_9_Axis_22_23", "SE_hKaon_12_Ka_Id0_3_9_Axis_22_23", kTH2F, {Axis_22, Axis_23}); // track.tpcNSigmaKa  (),track.tofNSigmaKa()) ;//Axis_tpcInnerParam ;

    // tpcSignal
    SE_recoKaon.add("SE_hKaon_12_Ka_Id1_1_Axis_00_05", "SE_hKaon_12_Ka_Id1_1_Axis_00_05", kTH2F, {Axis_00, Axis_05}); // track.p            (),track.tpcSignal()) ;//Axis_p             ;
    SE_recoKaon.add("SE_hKaon_12_Ka_Id1_1_Axis_02_05", "SE_hKaon_12_Ka_Id1_1_Axis_02_05", kTH2F, {Axis_02, Axis_05}); // track.tpcInnerParam(),track.tpcSignal()) ;//Axis_tpcInnerParam ;
    SE_recoKaon.add("SE_hKaon_12_Ka_Id1_1_Axis_03_05", "SE_hKaon_12_Ka_Id1_1_Axis_03_05", kTH2F, {Axis_03, Axis_05}); // track.tofExpMom    (),track.tpcSignal()) ;//Axis_tofExpMom     ;
    // tofBeta
    SE_recoKaon.add("SE_hKaon_12_Ka_Id1_2_Axis_00_06", "SE_hKaon_12_Ka_Id1_2_Axis_00_06", kTH2F, {Axis_00, Axis_06}); // track.p            (),track.beta()) ;//Axis_p             ;
    SE_recoKaon.add("SE_hKaon_12_Ka_Id1_2_Axis_02_06", "SE_hKaon_12_Ka_Id1_2_Axis_02_06", kTH2F, {Axis_02, Axis_06}); // track.tpcInnerParam(),track.beta()) ;//Axis_tpcInnerParam ;
    SE_recoKaon.add("SE_hKaon_12_Ka_Id1_2_Axis_03_06", "SE_hKaon_12_Ka_Id1_2_Axis_03_06", kTH2F, {Axis_03, Axis_06}); // track.tofExpMom    (),track.beta()) ;//Axis_tofExpMom     ;
    // Kaon
    SE_recoKaon.add("SE_hKaon_12_Ka_Id1_3_1_Axis_00_22", "SE_hKaon_12_Ka_Id1_3_1_Axis_00_22", kTH2F, {Axis_00, Axis_22}); // track.p            (),track.tpcNSigmaKa()) ;//Axis_p             ;
    SE_recoKaon.add("SE_hKaon_12_Ka_Id1_3_2_Axis_01_22", "SE_hKaon_12_Ka_Id1_3_2_Axis_01_22", kTH2F, {Axis_01, Axis_22}); // track.pt           (),track.tpcNSigmaKa()) ;//Axis_pt            ;
    SE_recoKaon.add("SE_hKaon_12_Ka_Id1_3_3_Axis_02_22", "SE_hKaon_12_Ka_Id1_3_3_Axis_02_22", kTH2F, {Axis_02, Axis_22}); // track.tpcInnerParam(),track.tpcNSigmaKa()) ;//Axis_tpcInnerParam ;
    SE_recoKaon.add("SE_hKaon_12_Ka_Id1_3_4_Axis_03_22", "SE_hKaon_12_Ka_Id1_3_4_Axis_03_22", kTH2F, {Axis_03, Axis_22}); // track.tofExpMom    (),track.tpcNSigmaKa()) ;//Axis_tofExpMom     ;
    SE_recoKaon.add("SE_hKaon_12_Ka_Id1_3_5_Axis_00_23", "SE_hKaon_12_Ka_Id1_3_5_Axis_00_23", kTH2F, {Axis_00, Axis_23}); // track.p            (),track.tofNSigmaKa()) ;//Axis_p             ;
    SE_recoKaon.add("SE_hKaon_12_Ka_Id1_3_6_Axis_01_23", "SE_hKaon_12_Ka_Id1_3_6_Axis_01_23", kTH2F, {Axis_01, Axis_23}); // track.pt           (),track.tofNSigmaKa()) ;//Axis_pt            ;
    SE_recoKaon.add("SE_hKaon_12_Ka_Id1_3_7_Axis_02_23", "SE_hKaon_12_Ka_Id1_3_7_Axis_02_23", kTH2F, {Axis_02, Axis_23}); // track.tpcInnerParam(),track.tofNSigmaKa()) ;//Axis_tpcInnerParam ;
    SE_recoKaon.add("SE_hKaon_12_Ka_Id1_3_8_Axis_03_23", "SE_hKaon_12_Ka_Id1_3_8_Axis_03_23", kTH2F, {Axis_03, Axis_23}); // track.tofExpMom    (),track.tofNSigmaKa()) ;//Axis_tofExpMom     ;
    SE_recoKaon.add("SE_hKaon_12_Ka_Id1_3_9_Axis_22_23", "SE_hKaon_12_Ka_Id1_3_9_Axis_22_23", kTH2F, {Axis_22, Axis_23}); // track.tpcNSigmaKa  (),track.tofNSigmaKa()) ;//Axis_tpcInnerParam ;
    // Kaon

    // trigger Hadrons
    // SE
    SE_recoTrigger.add("SE_hTrigger_0_TriggerTrack_P", "SE_hTrigger_TriggerTrack_P", {HistType::kTH1F, {Axis_p}});
    SE_recoTrigger.add("SE_hTrigger_0_TriggerTrack_tpcInnerParam", "SE_hTrigger_TriggerTrack_tpcInnerParam", {HistType::kTH1F, {Axis_tpcInnerParam}});
    SE_recoTrigger.add("SE_hTrigger_0_TriggerTrack_tofExpMom", "SE_hTrigger_TriggerTrack_tofExpMom", {HistType::kTH1F, {Axis_tofExpMom}});

    SE_recoTrigger.add("SE_hTrigger_1_TriggerTrack_Pt", "SE_hTrigger_TriggerTrack_Pt", {HistType::kTH1F, {Axis_pt}});
    SE_recoTrigger.add("SE_hTrigger_2_TriggerTrack_Eta", "SE_hTrigger_TriggerTrack_Eta", {HistType::kTH1F, {Axis_eta}});
    SE_recoTrigger.add("SE_hTrigger_3_TriggerTrack_Phi", "SE_hTrigger_TriggerTrack_Phi", {HistType::kTH1F, {Axis_phi}});
    SE_recoTrigger.add("SE_hTrigger_4_TriggerTrack_DcaXY", "SE_hTrigger_TriggerTrack_DcaXY", {HistType::kTH1F, {Axis_dcaXY}});
    SE_recoTrigger.add("SE_hTrigger_5_TriggerTrack_DcaZ", "SE_hTrigger_TriggerTrack_DcaZ", {HistType::kTH1F, {Axis_dcaZ}});
    SE_recoTrigger.add("SE_hTrigger_6_TriggerTrack_Sign", "SE_hTrigger_TriggerTrack_Sign", {HistType::kTH1D, {Axis_Sign}});
    SE_recoTrigger.add("SE_hTrigger_11_TriggerTrack_IdentificationTag", "SE_hTrigger_TriggerTrack_IdentificationTag", {HistType::kTH1D, {Axis_PidTag}});
    //
    // ME
    ME_recoTrigger.add("ME_hTrigger_0_TriggerTrack_P", "ME_hTrigger_TriggerTrack_P", {HistType::kTH1F, {Axis_p}});
    ME_recoTrigger.add("ME_hTrigger_0_TriggerTrack_tpcInnerParam", "ME_hTrigger_TriggerTrack_tpcInnerParam", {HistType::kTH1F, {Axis_tpcInnerParam}});
    ME_recoTrigger.add("ME_hTrigger_0_TriggerTrack_tofExpMom", "ME_hTrigger_TriggerTrack_tofExpMom", {HistType::kTH1F, {Axis_tofExpMom}});

    ME_recoTrigger.add("ME_hTrigger_1_TriggerTrack_Pt", "ME_hTrigger_TriggerTrack_Pt", {HistType::kTH1F, {Axis_pt}});
    ME_recoTrigger.add("ME_hTrigger_2_TriggerTrack_Eta", "ME_hTrigger_TriggerTrack_Eta", {HistType::kTH1F, {Axis_eta}});
    ME_recoTrigger.add("ME_hTrigger_3_TriggerTrack_Phi", "ME_hTrigger_TriggerTrack_Phi", {HistType::kTH1F, {Axis_phi}});
    ME_recoTrigger.add("ME_hTrigger_4_TriggerTrack_DcaXY", "ME_hTrigger_TriggerTrack_DcaXY", {HistType::kTH1F, {Axis_dcaXY}});
    ME_recoTrigger.add("ME_hTrigger_5_TriggerTrack_DcaZ", "ME_hTrigger_TriggerTrack_DcaZ", {HistType::kTH1F, {Axis_dcaZ}});
    ME_recoTrigger.add("ME_hTrigger_6_TriggerTrack_Sign", "ME_hTrigger_TriggerTrack_Sign", {HistType::kTH1D, {Axis_Sign}});
    ME_recoTrigger.add("ME_hTrigger_11_TriggerTrack_IdentificationTag", "ME_hTrigger_TriggerTrack_IdentificationTag", {HistType::kTH1D, {Axis_PidTag}});
    //
    // trigger Hadrons

    // Analysis
    // h-Phi correlations
    // Same Event
    // Unlike Sign Peak Region
    // phi_pT in [0.0-inf]
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_AllPhi_Phi_p", "SE_hAnalysis_hPhi_US_Peak_AllPhi_Phi_p", kTH1F, {Axis_p});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_AllPhi_Phi_pT", "SE_hAnalysis_hPhi_US_Peak_AllPhi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_AllPhi_Phi_eta", "SE_hAnalysis_hPhi_US_Peak_AllPhi_Phi_eta;#eta", kTH1F, {Axis_eta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_AllPhi_Phi_phi", "SE_hAnalysis_hPhi_US_Peak_AllPhi_Phi_phi", kTH1F, {Axis_phi});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_AllPhi_Phi_rapidity", "SE_hAnalysis_hPhi_US_Peak_AllPhi_Phi_rapidity", kTH1F, {Axis_rapidity});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_AllPhi_Phi_paircharge", "SE_hAnalysis_hPhi_US_Peak_AllPhi_Phi_paircharge", kTH1F, {Axis_paircharge});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_AllPhi_Phi_InvMass", "SE_hAnalysis_hPhi_US_Peak_AllPhi_Phi_InvMass;Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_AllPhi_dPhi_dEta", "SE_hAnalysis_hPhi_US_Peak_AllPhi_dPhi_dEta;#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_AllPhi_Analysis", "SE_hAnalysis_hPhi_US_Peak_AllPhi_Analysis", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // Phi_pT in [0.0-2.0](Bulk)
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_p", "SE_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_p", kTH1F, {Axis_p});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_pT", "SE_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_eta", "SE_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_eta;#eta", kTH1F, {Axis_eta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_phi", "SE_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_phi", kTH1F, {Axis_phi});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_rapidity", "SE_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_rapidity", kTH1F, {Axis_rapidity});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_paircharge", "SE_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_paircharge", kTH1F, {Axis_paircharge});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_InvMass", "SE_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_InvMass(Bulk);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_0To2Phi_dPhi_dEta", "SE_hAnalysis_hPhi_US_Peak_0To2Phi_dPhi_dEta(Bulk);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_0To2Phi_Analysis", "SE_hAnalysis_hPhi_US_Peak_0To2Phi_Analysis(Bulk)", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // Phi_pT in [2.0-4.0](Required)
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_p", "SE_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_p", kTH1F, {Axis_p});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_pT", "SE_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_eta", "SE_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_eta;#eta", kTH1F, {Axis_eta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_phi", "SE_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_phi", kTH1F, {Axis_phi});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_rapidity", "SE_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_rapidity", kTH1F, {Axis_rapidity});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_paircharge", "SE_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_paircharge", kTH1F, {Axis_paircharge});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_InvMass", "SE_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_InvMass(Required);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_2To4Phi_dPhi_dEta", "SE_hAnalysis_hPhi_US_Peak_2To4Phi_dPhi_dEta(Required);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_2To4Phi_Analysis", "SE_hAnalysis_hPhi_US_Peak_2To4Phi_Analysis(Required)", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // Phi_pT in [4.0-8.0](Hard)
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_p", "SE_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_p", kTH1F, {Axis_p});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_pT", "SE_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_eta", "SE_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_eta;#eta", kTH1F, {Axis_eta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_phi", "SE_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_phi", kTH1F, {Axis_phi});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_rapidity", "SE_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_rapidity", kTH1F, {Axis_rapidity});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_paircharge", "SE_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_paircharge", kTH1F, {Axis_paircharge});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_InvMass", "SE_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_InvMass(Hard);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_4To8Phi_dPhi_dEta", "SE_hAnalysis_hPhi_US_Peak_4To8Phi_dPhi_dEta(Hard);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_4To8Phi_Analysis", "SE_hAnalysis_hPhi_US_Peak_4To8Phi_Analysis(Hard)", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // Phi_pT in [8.0-inf](Very Hard)
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_p", "SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_p", kTH1F, {Axis_p});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_pT", "SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_eta", "SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_eta;#eta", kTH1F, {Axis_eta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_phi", "SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_phi", kTH1F, {Axis_phi});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_rapidity", "SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_rapidity", kTH1F, {Axis_rapidity});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_paircharge", "SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_paircharge", kTH1F, {Axis_paircharge});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_InvMass", "SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_InvMass(VeryHard);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_dPhi_dEta", "SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_dPhi_dEta(VeryHard);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_Analysis", "SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_Analysis(VeryHard)", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    //
    // Unlike Sign LSB Region
    // phi_pT in [0.0-inf]
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_AllPhi_Phi_p", "SE_hAnalysis_hPhi_US_LSB_AllPhi_Phi_p", kTH1F, {Axis_p});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_AllPhi_Phi_pT", "SE_hAnalysis_hPhi_US_LSB_AllPhi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_AllPhi_Phi_eta", "SE_hAnalysis_hPhi_US_LSB_AllPhi_Phi_eta;#eta", kTH1F, {Axis_eta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_AllPhi_Phi_phi", "SE_hAnalysis_hPhi_US_LSB_AllPhi_Phi_phi", kTH1F, {Axis_phi});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_AllPhi_Phi_rapidity", "SE_hAnalysis_hPhi_US_LSB_AllPhi_Phi_rapidity", kTH1F, {Axis_rapidity});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_AllPhi_Phi_paircharge", "SE_hAnalysis_hPhi_US_LSB_AllPhi_Phi_paircharge", kTH1F, {Axis_paircharge});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_AllPhi_Phi_InvMass", "SE_hAnalysis_hPhi_US_LSB_AllPhi_Phi_InvMass;Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_AllPhi_dPhi_dEta", "SE_hAnalysis_hPhi_US_LSB_AllPhi_dPhi_dEta;#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_AllPhi_Analysis", "SE_hAnalysis_hPhi_US_LSB_AllPhi_Analysis", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // Phi_pT in [0.0-2.0](Bulk)
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_p", "SE_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_p", kTH1F, {Axis_p});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_pT", "SE_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_eta", "SE_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_eta;#eta", kTH1F, {Axis_eta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_phi", "SE_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_phi", kTH1F, {Axis_phi});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_rapidity", "SE_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_rapidity", kTH1F, {Axis_rapidity});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_paircharge", "SE_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_paircharge", kTH1F, {Axis_paircharge});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_InvMass", "SE_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_InvMass(Bulk);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_0To2Phi_dPhi_dEta", "SE_hAnalysis_hPhi_US_LSB_0To2Phi_dPhi_dEta(Bulk);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_0To2Phi_Analysis", "SE_hAnalysis_hPhi_US_LSB_0To2Phi_Analysis(Bulk)", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // Phi_pT in [2.0-4.0](Required)
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_p", "SE_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_p", kTH1F, {Axis_p});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_pT", "SE_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_eta", "SE_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_eta;#eta", kTH1F, {Axis_eta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_phi", "SE_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_phi", kTH1F, {Axis_phi});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_rapidity", "SE_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_rapidity", kTH1F, {Axis_rapidity});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_paircharge", "SE_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_paircharge", kTH1F, {Axis_paircharge});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_InvMass", "SE_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_InvMass(Required);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_2To4Phi_dPhi_dEta", "SE_hAnalysis_hPhi_US_LSB_2To4Phi_dPhi_dEta(Required);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_2To4Phi_Analysis", "SE_hAnalysis_hPhi_US_LSB_2To4Phi_Analysis(Required)", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // Phi_pT in [4.0-8.0](Hard)
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_p", "SE_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_p", kTH1F, {Axis_p});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_pT", "SE_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_eta", "SE_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_eta;#eta", kTH1F, {Axis_eta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_phi", "SE_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_phi", kTH1F, {Axis_phi});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_rapidity", "SE_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_rapidity", kTH1F, {Axis_rapidity});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_paircharge", "SE_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_paircharge", kTH1F, {Axis_paircharge});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_InvMass", "SE_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_InvMass(Hard);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_4To8Phi_dPhi_dEta", "SE_hAnalysis_hPhi_US_LSB_4To8Phi_dPhi_dEta(Hard);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_4To8Phi_Analysis", "SE_hAnalysis_hPhi_US_LSB_4To8Phi_Analysis(Hard)", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // Phi_pT in [8.0-inf](Very Hard)
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_p", "SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_p", kTH1F, {Axis_p});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_pT", "SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_eta", "SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_eta;#eta", kTH1F, {Axis_eta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_phi", "SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_phi", kTH1F, {Axis_phi});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_rapidity", "SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_rapidity", kTH1F, {Axis_rapidity});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_paircharge", "SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_paircharge", kTH1F, {Axis_paircharge});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_InvMass", "SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_InvMass(VeryHard);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_dPhi_dEta", "SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_dPhi_dEta(VeryHard);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_Analysis", "SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_Analysis(VeryHard)", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    //
    // Unlike Sign RSB Region
    // phi_pT in [0.0-inf]
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_AllPhi_Phi_p", "SE_hAnalysis_hPhi_US_RSB_AllPhi_Phi_p", kTH1F, {Axis_p});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_AllPhi_Phi_pT", "SE_hAnalysis_hPhi_US_RSB_AllPhi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_AllPhi_Phi_eta", "SE_hAnalysis_hPhi_US_RSB_AllPhi_Phi_eta;#eta", kTH1F, {Axis_eta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_AllPhi_Phi_phi", "SE_hAnalysis_hPhi_US_RSB_AllPhi_Phi_phi", kTH1F, {Axis_phi});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_AllPhi_Phi_rapidity", "SE_hAnalysis_hPhi_US_RSB_AllPhi_Phi_rapidity", kTH1F, {Axis_rapidity});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_AllPhi_Phi_paircharge", "SE_hAnalysis_hPhi_US_RSB_AllPhi_Phi_paircharge", kTH1F, {Axis_paircharge});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_AllPhi_Phi_InvMass", "SE_hAnalysis_hPhi_US_RSB_AllPhi_Phi_InvMass;Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_AllPhi_dPhi_dEta", "SE_hAnalysis_hPhi_US_RSB_AllPhi_dPhi_dEta;#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_AllPhi_Analysis", "SE_hAnalysis_hPhi_US_RSB_AllPhi_Analysis", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // Phi_pT in [0.0-2.0](Bulk)
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_p", "SE_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_p", kTH1F, {Axis_p});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_pT", "SE_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_eta", "SE_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_eta;#eta", kTH1F, {Axis_eta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_phi", "SE_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_phi", kTH1F, {Axis_phi});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_rapidity", "SE_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_rapidity", kTH1F, {Axis_rapidity});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_paircharge", "SE_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_paircharge", kTH1F, {Axis_paircharge});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_InvMass", "SE_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_InvMass(Bulk);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_0To2Phi_dPhi_dEta", "SE_hAnalysis_hPhi_US_RSB_0To2Phi_dPhi_dEta(Bulk);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_0To2Phi_Analysis", "SE_hAnalysis_hPhi_US_RSB_0To2Phi_Analysis(Bulk)", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // Phi_pT in [2.0-4.0](Required)
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_p", "SE_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_p", kTH1F, {Axis_p});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_pT", "SE_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_eta", "SE_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_eta;#eta", kTH1F, {Axis_eta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_phi", "SE_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_phi", kTH1F, {Axis_phi});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_rapidity", "SE_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_rapidity", kTH1F, {Axis_rapidity});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_paircharge", "SE_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_paircharge", kTH1F, {Axis_paircharge});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_InvMass", "SE_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_InvMass(Required);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_2To4Phi_dPhi_dEta", "SE_hAnalysis_hPhi_US_RSB_2To4Phi_dPhi_dEta(Required);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_2To4Phi_Analysis", "SE_hAnalysis_hPhi_US_RSB_2To4Phi_Analysis(Required)", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // Phi_pT in [4.0-8.0](Hard)
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_p", "SE_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_p", kTH1F, {Axis_p});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_pT", "SE_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_eta", "SE_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_eta;#eta", kTH1F, {Axis_eta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_phi", "SE_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_phi", kTH1F, {Axis_phi});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_rapidity", "SE_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_rapidity", kTH1F, {Axis_rapidity});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_paircharge", "SE_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_paircharge", kTH1F, {Axis_paircharge});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_InvMass", "SE_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_InvMass(Hard);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_4To8Phi_dPhi_dEta", "SE_hAnalysis_hPhi_US_RSB_4To8Phi_dPhi_dEta(Hard);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_4To8Phi_Analysis", "SE_hAnalysis_hPhi_US_RSB_4To8Phi_Analysis(Hard)", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // Phi_pT in [8.0-inf](Very Hard)
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_p", "SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_p", kTH1F, {Axis_p});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_pT", "SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_eta", "SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_eta;#eta", kTH1F, {Axis_eta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_phi", "SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_phi", kTH1F, {Axis_phi});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_rapidity", "SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_rapidity", kTH1F, {Axis_rapidity});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_paircharge", "SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_paircharge", kTH1F, {Axis_paircharge});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_InvMass", "SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_InvMass(VeryHard);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_dPhi_dEta", "SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_dPhi_dEta(VeryHard);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_Analysis", "SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_Analysis(VeryHard)", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    //
    //
    // Mixed Event
    // Unlike Sign Peak Region
    // phi_pT in [0.0-inf]
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_AllPhi_Phi_p", "ME_hAnalysis_hPhi_US_Peak_AllPhi_Phi_p", kTH1F, {Axis_p});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_AllPhi_Phi_pT", "ME_hAnalysis_hPhi_US_Peak_AllPhi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_AllPhi_Phi_eta", "ME_hAnalysis_hPhi_US_Peak_AllPhi_Phi_eta;#eta", kTH1F, {Axis_eta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_AllPhi_Phi_phi", "ME_hAnalysis_hPhi_US_Peak_AllPhi_Phi_phi", kTH1F, {Axis_phi});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_AllPhi_Phi_rapidity", "ME_hAnalysis_hPhi_US_Peak_AllPhi_Phi_rapidity", kTH1F, {Axis_rapidity});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_AllPhi_Phi_paircharge", "ME_hAnalysis_hPhi_US_Peak_AllPhi_Phi_paircharge", kTH1F, {Axis_paircharge});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_AllPhi_Phi_InvMass", "ME_hAnalysis_hPhi_US_Peak_AllPhi_Phi_InvMass;Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_AllPhi_dPhi_dEta", "ME_hAnalysis_hPhi_US_Peak_AllPhi_dPhi_dEta;#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_AllPhi_Analysis", "ME_hAnalysis_hPhi_US_Peak_AllPhi_Analysis", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // Phi_pT in [0.0-2.0](Bulk)
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_p", "ME_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_p", kTH1F, {Axis_p});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_pT", "ME_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_eta", "ME_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_eta;#eta", kTH1F, {Axis_eta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_phi", "ME_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_phi", kTH1F, {Axis_phi});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_rapidity", "ME_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_rapidity", kTH1F, {Axis_rapidity});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_paircharge", "ME_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_paircharge", kTH1F, {Axis_paircharge});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_InvMass", "ME_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_InvMass(Bulk);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_0To2Phi_dPhi_dEta", "ME_hAnalysis_hPhi_US_Peak_0To2Phi_dPhi_dEta(Bulk);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_0To2Phi_Analysis", "ME_hAnalysis_hPhi_US_Peak_0To2Phi_Analysis(Bulk)", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // Phi_pT in [2.0-4.0](Required)
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_p", "ME_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_p", kTH1F, {Axis_p});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_pT", "ME_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_eta", "ME_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_eta;#eta", kTH1F, {Axis_eta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_phi", "ME_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_phi", kTH1F, {Axis_phi});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_rapidity", "ME_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_rapidity", kTH1F, {Axis_rapidity});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_paircharge", "ME_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_paircharge", kTH1F, {Axis_paircharge});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_InvMass", "ME_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_InvMass(Required);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_2To4Phi_dPhi_dEta", "ME_hAnalysis_hPhi_US_Peak_2To4Phi_dPhi_dEta(Required);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_2To4Phi_Analysis", "ME_hAnalysis_hPhi_US_Peak_2To4Phi_Analysis(Required)", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // Phi_pT in [4.0-8.0](Hard)
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_p", "ME_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_p", kTH1F, {Axis_p});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_pT", "ME_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_eta", "ME_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_eta;#eta", kTH1F, {Axis_eta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_phi", "ME_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_phi", kTH1F, {Axis_phi});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_rapidity", "ME_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_rapidity", kTH1F, {Axis_rapidity});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_paircharge", "ME_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_paircharge", kTH1F, {Axis_paircharge});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_InvMass", "ME_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_InvMass(Hard);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_4To8Phi_dPhi_dEta", "ME_hAnalysis_hPhi_US_Peak_4To8Phi_dPhi_dEta(Hard);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_4To8Phi_Analysis", "ME_hAnalysis_hPhi_US_Peak_4To8Phi_Analysis(Hard)", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // Phi_pT in [8.0-inf](Very Hard)
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_p", "ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_p", kTH1F, {Axis_p});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_pT", "ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_eta", "ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_eta;#eta", kTH1F, {Axis_eta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_phi", "ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_phi", kTH1F, {Axis_phi});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_rapidity", "ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_rapidity", kTH1F, {Axis_rapidity});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_paircharge", "ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_paircharge", kTH1F, {Axis_paircharge});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_InvMass", "ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_InvMass(VeryHard);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_dPhi_dEta", "ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_dPhi_dEta(VeryHard);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_Analysis", "ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_Analysis(VeryHard)", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    //
    // Unlike Sign LSB Region
    // phi_pT in [0.0-inf]
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_AllPhi_Phi_p", "ME_hAnalysis_hPhi_US_LSB_AllPhi_Phi_p", kTH1F, {Axis_p});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_AllPhi_Phi_pT", "ME_hAnalysis_hPhi_US_LSB_AllPhi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_AllPhi_Phi_eta", "ME_hAnalysis_hPhi_US_LSB_AllPhi_Phi_eta;#eta", kTH1F, {Axis_eta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_AllPhi_Phi_phi", "ME_hAnalysis_hPhi_US_LSB_AllPhi_Phi_phi", kTH1F, {Axis_phi});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_AllPhi_Phi_rapidity", "ME_hAnalysis_hPhi_US_LSB_AllPhi_Phi_rapidity", kTH1F, {Axis_rapidity});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_AllPhi_Phi_paircharge", "ME_hAnalysis_hPhi_US_LSB_AllPhi_Phi_paircharge", kTH1F, {Axis_paircharge});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_AllPhi_Phi_InvMass", "ME_hAnalysis_hPhi_US_LSB_AllPhi_Phi_InvMass;Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_AllPhi_dPhi_dEta", "ME_hAnalysis_hPhi_US_LSB_AllPhi_dPhi_dEta;#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_AllPhi_Analysis", "ME_hAnalysis_hPhi_US_LSB_AllPhi_Analysis", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // Phi_pT in [0.0-2.0](Bulk)
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_p", "ME_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_p", kTH1F, {Axis_p});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_pT", "ME_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_eta", "ME_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_eta;#eta", kTH1F, {Axis_eta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_phi", "ME_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_phi", kTH1F, {Axis_phi});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_rapidity", "ME_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_rapidity", kTH1F, {Axis_rapidity});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_paircharge", "ME_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_paircharge", kTH1F, {Axis_paircharge});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_InvMass", "ME_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_InvMass(Bulk);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_0To2Phi_dPhi_dEta", "ME_hAnalysis_hPhi_US_LSB_0To2Phi_dPhi_dEta(Bulk);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_0To2Phi_Analysis", "ME_hAnalysis_hPhi_US_LSB_0To2Phi_Analysis(Bulk)", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // Phi_pT in [2.0-4.0](Required)
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_p", "ME_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_p", kTH1F, {Axis_p});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_pT", "ME_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_eta", "ME_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_eta;#eta", kTH1F, {Axis_eta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_phi", "ME_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_phi", kTH1F, {Axis_phi});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_rapidity", "ME_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_rapidity", kTH1F, {Axis_rapidity});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_paircharge", "ME_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_paircharge", kTH1F, {Axis_paircharge});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_InvMass", "ME_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_InvMass(Required);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_2To4Phi_dPhi_dEta", "ME_hAnalysis_hPhi_US_LSB_2To4Phi_dPhi_dEta(Required);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_2To4Phi_Analysis", "ME_hAnalysis_hPhi_US_LSB_2To4Phi_Analysis(Required)", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // Phi_pT in [4.0-8.0](Hard)
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_p", "ME_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_p", kTH1F, {Axis_p});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_pT", "ME_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_eta", "ME_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_eta;#eta", kTH1F, {Axis_eta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_phi", "ME_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_phi", kTH1F, {Axis_phi});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_rapidity", "ME_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_rapidity", kTH1F, {Axis_rapidity});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_paircharge", "ME_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_paircharge", kTH1F, {Axis_paircharge});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_InvMass", "ME_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_InvMass(Hard);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_4To8Phi_dPhi_dEta", "ME_hAnalysis_hPhi_US_LSB_4To8Phi_dPhi_dEta(Hard);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_4To8Phi_Analysis", "ME_hAnalysis_hPhi_US_LSB_4To8Phi_Analysis(Hard)", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // Phi_pT in [8.0-inf](Very Hard)
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_p", "ME_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_p", kTH1F, {Axis_p});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_pT", "ME_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_eta", "ME_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_eta;#eta", kTH1F, {Axis_eta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_phi", "ME_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_phi", kTH1F, {Axis_phi});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_rapidity", "ME_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_rapidity", kTH1F, {Axis_rapidity});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_paircharge", "ME_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_paircharge", kTH1F, {Axis_paircharge});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_InvMass", "ME_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_InvMass(VeryHard);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_8ToInfPhi_dPhi_dEta", "ME_hAnalysis_hPhi_US_LSB_8ToInfPhi_dPhi_dEta(VeryHard);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_LSB_8ToInfPhi_Analysis", "ME_hAnalysis_hPhi_US_LSB_8ToInfPhi_Analysis(VeryHard)", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    //
    // Unlike Sign RSB Region
    // phi_pT in [0.0-inf]
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_AllPhi_Phi_p", "ME_hAnalysis_hPhi_US_RSB_AllPhi_Phi_p", kTH1F, {Axis_p});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_AllPhi_Phi_pT", "ME_hAnalysis_hPhi_US_RSB_AllPhi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_AllPhi_Phi_eta", "ME_hAnalysis_hPhi_US_RSB_AllPhi_Phi_eta;#eta", kTH1F, {Axis_eta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_AllPhi_Phi_phi", "ME_hAnalysis_hPhi_US_RSB_AllPhi_Phi_phi", kTH1F, {Axis_phi});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_AllPhi_Phi_rapidity", "ME_hAnalysis_hPhi_US_RSB_AllPhi_Phi_rapidity", kTH1F, {Axis_rapidity});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_AllPhi_Phi_paircharge", "ME_hAnalysis_hPhi_US_RSB_AllPhi_Phi_paircharge", kTH1F, {Axis_paircharge});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_AllPhi_Phi_InvMass", "ME_hAnalysis_hPhi_US_RSB_AllPhi_Phi_InvMass;Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_AllPhi_dPhi_dEta", "ME_hAnalysis_hPhi_US_RSB_AllPhi_dPhi_dEta;#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_AllPhi_Analysis", "ME_hAnalysis_hPhi_US_RSB_AllPhi_Analysis", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // Phi_pT in [0.0-2.0](Bulk)
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_p", "ME_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_p", kTH1F, {Axis_p});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_pT", "ME_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_eta", "ME_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_eta;#eta", kTH1F, {Axis_eta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_phi", "ME_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_phi", kTH1F, {Axis_phi});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_rapidity", "ME_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_rapidity", kTH1F, {Axis_rapidity});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_paircharge", "ME_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_paircharge", kTH1F, {Axis_paircharge});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_InvMass", "ME_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_InvMass(Bulk);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_0To2Phi_dPhi_dEta", "ME_hAnalysis_hPhi_US_RSB_0To2Phi_dPhi_dEta(Bulk);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_0To2Phi_Analysis", "ME_hAnalysis_hPhi_US_RSB_0To2Phi_Analysis(Bulk)", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // Phi_pT in [2.0-4.0](Required)
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_p", "ME_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_p", kTH1F, {Axis_p});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_pT", "ME_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_eta", "ME_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_eta;#eta", kTH1F, {Axis_eta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_phi", "ME_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_phi", kTH1F, {Axis_phi});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_rapidity", "ME_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_rapidity", kTH1F, {Axis_rapidity});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_paircharge", "ME_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_paircharge", kTH1F, {Axis_paircharge});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_InvMass", "ME_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_InvMass(Required);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_2To4Phi_dPhi_dEta", "ME_hAnalysis_hPhi_US_RSB_2To4Phi_dPhi_dEta(Required);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_2To4Phi_Analysis", "ME_hAnalysis_hPhi_US_RSB_2To4Phi_Analysis(Required)", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // Phi_pT in [4.0-8.0](Hard)
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_p", "ME_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_p", kTH1F, {Axis_p});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_pT", "ME_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_eta", "ME_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_eta;#eta", kTH1F, {Axis_eta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_phi", "ME_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_phi", kTH1F, {Axis_phi});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_rapidity", "ME_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_rapidity", kTH1F, {Axis_rapidity});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_paircharge", "ME_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_paircharge", kTH1F, {Axis_paircharge});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_InvMass", "ME_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_InvMass(Hard);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_4To8Phi_dPhi_dEta", "ME_hAnalysis_hPhi_US_RSB_4To8Phi_dPhi_dEta(Hard);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_4To8Phi_Analysis", "ME_hAnalysis_hPhi_US_RSB_4To8Phi_Analysis(Hard)", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // Phi_pT in [8.0-inf](Very Hard)
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_p", "ME_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_p", kTH1F, {Axis_p});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_pT", "ME_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_eta", "ME_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_eta;#eta", kTH1F, {Axis_eta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_phi", "ME_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_phi", kTH1F, {Axis_phi});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_rapidity", "ME_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_rapidity", kTH1F, {Axis_rapidity});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_paircharge", "ME_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_paircharge", kTH1F, {Axis_paircharge});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_InvMass", "ME_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_InvMass(VeryHard);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_8ToInfPhi_dPhi_dEta", "ME_hAnalysis_hPhi_US_RSB_8ToInfPhi_dPhi_dEta(VeryHard);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    ME_recoAnalysis.add("ME_hAnalysis_hPhi_US_RSB_8ToInfPhi_Analysis", "ME_hAnalysis_hPhi_US_RSB_8ToInfPhi_Analysis(VeryHard)", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    //
    //
    //
    //

    // PhiMesons
    SE_recoPhi.add("SE_hPhi_KK_allPairs_01_Mass", "SE_hPhi_KK_allPairs_01_Mass", kTH1F, {Axis_PhiMass});
    SE_recoPhi.add("SE_hPhi_KK_allPairs_02_P", "SE_hPhi_KK_allPairs_02_P", kTH1F, {Axis_p});
    SE_recoPhi.add("SE_hPhi_KK_allPairs_03_Pt", "SE_hPhi_KK_allPairs_03_Pt", kTH1F, {Axis_pt});
    SE_recoPhi.add("SE_hPhi_KK_allPairs_04_Eta", "SE_hPhi_KK_allPairs_04_Eta", kTH1F, {Axis_eta});
    SE_recoPhi.add("SE_hPhi_KK_allPairs_05_Phi", "SE_hPhi_KK_allPairs_05_Phi", kTH1F, {Axis_phi});
    SE_recoPhi.add("SE_hPhi_KK_allPairs_06_Rapidity", "SE_hPhi_KK_allPairs_06_Rapidity", kTH1F, {Axis_rapidity});

    SE_recoPhi.add("SE_hPhi_KK_USpairs_01_Mass", "SE_hPhi_KK_USpairs_01_Mass", kTH1F, {Axis_PhiMass});
    SE_recoPhi.add("SE_hPhi_KK_USpairs_02_P", "SE_hPhi_KK_USpairs_02_P", kTH1F, {Axis_p});
    SE_recoPhi.add("SE_hPhi_KK_USpairs_03_Pt", "SE_hPhi_KK_USpairs_03_Pt", kTH1F, {Axis_pt});
    SE_recoPhi.add("SE_hPhi_KK_USpairs_04_Eta", "SE_hPhi_KK_USpairs_04_Eta", kTH1F, {Axis_eta});
    SE_recoPhi.add("SE_hPhi_KK_USpairs_05_Phi", "SE_hPhi_KK_USpairs_05_Phi", kTH1F, {Axis_phi});
    SE_recoPhi.add("SE_hPhi_KK_USpairs_06_Rapidity", "SE_hPhi_KK_USpairs_06_Rapidity", kTH1F, {Axis_rapidity});

    SE_recoPhi.add("SE_hPhi_KK_LSpairs_01_Mass", "SE_hPhi_KK_LSpairs_01_Mass", kTH1F, {Axis_PhiMass});
    SE_recoPhi.add("SE_hPhi_KK_LSpairs_02_P", "SE_hPhi_KK_LSpairs_02_P", kTH1F, {Axis_p});
    SE_recoPhi.add("SE_hPhi_KK_LSpairs_03_Pt", "SE_hPhi_KK_LSpairs_03_Pt", kTH1F, {Axis_pt});
    SE_recoPhi.add("SE_hPhi_KK_LSpairs_04_Eta", "SE_hPhi_KK_LSpairs_04_Eta", kTH1F, {Axis_eta});
    SE_recoPhi.add("SE_hPhi_KK_LSpairs_05_Phi", "SE_hPhi_KK_LSpairs_05_Phi", kTH1F, {Axis_phi});
    SE_recoPhi.add("SE_hPhi_KK_LSpairs_06_Rapidity", "SE_hPhi_KK_LSpairs_06_Rapidity", kTH1F, {Axis_rapidity});

    // Phi Count and Common Kaons
    SE_recoPhi.add("SE_hPhi_nPhiCandidates_US_Peak", "SE_hPhi_nPhiCandidates_US_Peak", kTH1F, {{12 * 2, -1.5, 10.5}});
    SE_recoPhi.add("SE_hPhi_nPhiCandidates_PeakLS", "SE_hPhi_nPhiCandidates_PeakLS", kTH1F, {{12 * 2, -1.5, 10.5}});

    SE_recoPhi.add("SE_hPhi_nCommonKaonInDifferentPhi_US_Peak", "SE_hPhi_nCommonKaonInDifferentPhi_US_Peak", kTH1F, {{12 * 2, -1.5, 10.5}});
    SE_recoPhi.add("SE_hPhi_nCommonKaonInDifferentPhi_PeakLS", "SE_hPhi_nCommonKaonInDifferentPhi_PeakLS", kTH1F, {{12 * 2, -1.5, 10.5}});
    // PhiMesons

    // Like Sign Invariant Mass Line
    // phi_pT in [0.0-inf]
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_AllPhi_Phi_p", "SE_hAnalysis_hPhi_LS_KK_AllPhi_Phi_p", kTH1F, {Axis_p});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_AllPhi_Phi_pT", "SE_hAnalysis_hPhi_LS_KK_AllPhi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_AllPhi_Phi_eta", "SE_hAnalysis_hPhi_LS_KK_AllPhi_Phi_eta;#eta", kTH1F, {Axis_eta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_AllPhi_Phi_phi", "SE_hAnalysis_hPhi_LS_KK_AllPhi_Phi_phi", kTH1F, {Axis_phi});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_AllPhi_Phi_rapidity", "SE_hAnalysis_hPhi_LS_KK_AllPhi_Phi_rapidity", kTH1F, {Axis_rapidity});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_AllPhi_Phi_paircharge", "SE_hAnalysis_hPhi_LS_KK_AllPhi_Phi_paircharge", kTH1F, {Axis_paircharge});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_AllPhi_Phi_InvMass", "SE_hAnalysis_hPhi_LS_KK_AllPhi_Phi_InvMass;Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    // Phi_pT in [0.0-2.0](Bulk)
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_0To2Phi_Phi_p", "SE_hAnalysis_hPhi_LS_KK_0To2Phi_Phi_p", kTH1F, {Axis_p});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_0To2Phi_Phi_pT", "SE_hAnalysis_hPhi_LS_KK_0To2Phi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_0To2Phi_Phi_eta", "SE_hAnalysis_hPhi_LS_KK_0To2Phi_Phi_eta;#eta", kTH1F, {Axis_eta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_0To2Phi_Phi_phi", "SE_hAnalysis_hPhi_LS_KK_0To2Phi_Phi_phi", kTH1F, {Axis_phi});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_0To2Phi_Phi_rapidity", "SE_hAnalysis_hPhi_LS_KK_0To2Phi_Phi_rapidity", kTH1F, {Axis_rapidity});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_0To2Phi_Phi_paircharge", "SE_hAnalysis_hPhi_LS_KK_0To2Phi_Phi_paircharge", kTH1F, {Axis_paircharge});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_0To2Phi_Phi_InvMass", "SE_hAnalysis_hPhi_LS_KK_0To2Phi_Phi_InvMass(Bulk);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    // Phi_pT in [2.0-4.0](Required)
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_2To4Phi_Phi_p", "SE_hAnalysis_hPhi_LS_KK_2To4Phi_Phi_p", kTH1F, {Axis_p});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_2To4Phi_Phi_pT", "SE_hAnalysis_hPhi_LS_KK_2To4Phi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_2To4Phi_Phi_eta", "SE_hAnalysis_hPhi_LS_KK_2To4Phi_Phi_eta;#eta", kTH1F, {Axis_eta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_2To4Phi_Phi_phi", "SE_hAnalysis_hPhi_LS_KK_2To4Phi_Phi_phi", kTH1F, {Axis_phi});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_2To4Phi_Phi_rapidity", "SE_hAnalysis_hPhi_LS_KK_2To4Phi_Phi_rapidity", kTH1F, {Axis_rapidity});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_2To4Phi_Phi_paircharge", "SE_hAnalysis_hPhi_LS_KK_2To4Phi_Phi_paircharge", kTH1F, {Axis_paircharge});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_2To4Phi_Phi_InvMass", "SE_hAnalysis_hPhi_LS_KK_2To4Phi_Phi_InvMass(Required);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    // Phi_pT in [4.0-8.0](Hard)
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_4To8Phi_Phi_p", "SE_hAnalysis_hPhi_LS_KK_4To8Phi_Phi_p", kTH1F, {Axis_p});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_4To8Phi_Phi_pT", "SE_hAnalysis_hPhi_LS_KK_4To8Phi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_4To8Phi_Phi_eta", "SE_hAnalysis_hPhi_LS_KK_4To8Phi_Phi_eta;#eta", kTH1F, {Axis_eta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_4To8Phi_Phi_phi", "SE_hAnalysis_hPhi_LS_KK_4To8Phi_Phi_phi", kTH1F, {Axis_phi});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_4To8Phi_Phi_rapidity", "SE_hAnalysis_hPhi_LS_KK_4To8Phi_Phi_rapidity", kTH1F, {Axis_rapidity});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_4To8Phi_Phi_paircharge", "SE_hAnalysis_hPhi_LS_KK_4To8Phi_Phi_paircharge", kTH1F, {Axis_paircharge});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_4To8Phi_Phi_InvMass", "SE_hAnalysis_hPhi_LS_KK_4To8Phi_Phi_InvMass(Hard);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    // Phi_pT in [8.0-inf](Very Hard)
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_8ToInfPhi_Phi_p", "SE_hAnalysis_hPhi_LS_KK_8ToInfPhi_Phi_p", kTH1F, {Axis_p});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_8ToInfPhi_Phi_pT", "SE_hAnalysis_hPhi_LS_KK_8ToInfPhi_Phi_pT;p_{T}", kTH1F, {Axis_pt});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_8ToInfPhi_Phi_eta", "SE_hAnalysis_hPhi_LS_KK_8ToInfPhi_Phi_eta;#eta", kTH1F, {Axis_eta});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_8ToInfPhi_Phi_phi", "SE_hAnalysis_hPhi_LS_KK_8ToInfPhi_Phi_phi", kTH1F, {Axis_phi});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_8ToInfPhi_Phi_rapidity", "SE_hAnalysis_hPhi_LS_KK_8ToInfPhi_Phi_rapidity", kTH1F, {Axis_rapidity});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_8ToInfPhi_Phi_paircharge", "SE_hAnalysis_hPhi_LS_KK_8ToInfPhi_Phi_paircharge", kTH1F, {Axis_paircharge});
    SE_recoAnalysis.add("SE_hAnalysis_hPhi_LS_KK_8ToInfPhi_Phi_InvMass", "SE_hAnalysis_hPhi_LS_KK_8ToInfPhi_Phi_InvMass(VeryHard);Mass(K^{+} K^{-})", kTH1F, {Axis_PhiMass});
    //
    //

    // LeadingPhi-Phi Correlations
    //
    SE_recoPhi.add("SE_hPhiPhi_Full_LeadPhi_01_Mass", "SE_hPhiPhi_Full_LeadPhi_01_Mass", kTH1F, {Axis_PhiMass});
    SE_recoPhi.add("SE_hPhiPhi_Full_LeadPhi_02_P", "SE_hPhiPhi_Full_LeadPhi_02_P", kTH1F, {Axis_p});
    SE_recoPhi.add("SE_hPhiPhi_Full_LeadPhi_03_Pt", "SE_hPhiPhi_Full_LeadPhi_03_Pt", kTH1F, {Axis_pt});
    SE_recoPhi.add("SE_hPhiPhi_Full_LeadPhi_04_Eta", "SE_hPhiPhi_Full_LeadPhi_04_Eta", kTH1F, {Axis_eta});
    SE_recoPhi.add("SE_hPhiPhi_Full_LeadPhi_05_Phi", "SE_hPhiPhi_Full_LeadPhi_05_Phi", kTH1F, {Axis_phi});
    SE_recoPhi.add("SE_hPhiPhi_Full_LeadPhi_06_Rapidity", "SE_hPhiPhi_Full_LeadPhi_06_Rapidity", kTH1F, {Axis_rapidity});
    SE_recoPhi.add("SE_hPhiPhi_Full_AssoPhi_01_Mass", "SE_hPhiPhi_Full_AssoPhi_01_Mass", kTH1F, {Axis_PhiMass});
    SE_recoPhi.add("SE_hPhiPhi_Full_AssoPhi_02_P", "SE_hPhiPhi_Full_AssoPhi_02_P", kTH1F, {Axis_p});
    SE_recoPhi.add("SE_hPhiPhi_Full_AssoPhi_03_Pt", "SE_hPhiPhi_Full_AssoPhi_03_Pt", kTH1F, {Axis_pt});
    SE_recoPhi.add("SE_hPhiPhi_Full_AssoPhi_04_Eta", "SE_hPhiPhi_Full_AssoPhi_04_Eta", kTH1F, {Axis_eta});
    SE_recoPhi.add("SE_hPhiPhi_Full_AssoPhi_05_Phi", "SE_hPhiPhi_Full_AssoPhi_05_Phi", kTH1F, {Axis_phi});
    SE_recoPhi.add("SE_hPhiPhi_Full_AssoPhi_06_Rapidity", "SE_hPhiPhi_Full_AssoPhi_06_Rapidity", kTH1F, {Axis_rapidity});

    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_Full_LeadPhi_01_Mass", "SE_hPhiPhi_4to8LeadPhi_Full_LeadPhi_01_Mass", kTH1F, {Axis_PhiMass});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_Full_LeadPhi_02_P", "SE_hPhiPhi_4to8LeadPhi_Full_LeadPhi_02_P", kTH1F, {Axis_p});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_Full_LeadPhi_03_Pt", "SE_hPhiPhi_4to8LeadPhi_Full_LeadPhi_03_Pt", kTH1F, {Axis_pt});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_Full_LeadPhi_04_Eta", "SE_hPhiPhi_4to8LeadPhi_Full_LeadPhi_04_Eta", kTH1F, {Axis_eta});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_Full_LeadPhi_05_Phi", "SE_hPhiPhi_4to8LeadPhi_Full_LeadPhi_05_Phi", kTH1F, {Axis_phi});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_Full_LeadPhi_06_Rapidity", "SE_hPhiPhi_4to8LeadPhi_Full_LeadPhi_06_Rapidity", kTH1F, {Axis_rapidity});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_Full_AssoPhi_01_Mass", "SE_hPhiPhi_4to8LeadPhi_Full_AssoPhi_01_Mass", kTH1F, {Axis_PhiMass});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_Full_AssoPhi_02_P", "SE_hPhiPhi_4to8LeadPhi_Full_AssoPhi_02_P", kTH1F, {Axis_p});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_Full_AssoPhi_03_Pt", "SE_hPhiPhi_4to8LeadPhi_Full_AssoPhi_03_Pt", kTH1F, {Axis_pt});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_Full_AssoPhi_04_Eta", "SE_hPhiPhi_4to8LeadPhi_Full_AssoPhi_04_Eta", kTH1F, {Axis_eta});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_Full_AssoPhi_05_Phi", "SE_hPhiPhi_4to8LeadPhi_Full_AssoPhi_05_Phi", kTH1F, {Axis_phi});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_Full_AssoPhi_06_Rapidity", "SE_hPhiPhi_4to8LeadPhi_Full_AssoPhi_06_Rapidity", kTH1F, {Axis_rapidity});

    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_01_Mass", "SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_01_Mass", kTH1F, {Axis_PhiMass});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_02_P", "SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_02_P", kTH1F, {Axis_p});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_03_Pt", "SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_03_Pt", kTH1F, {Axis_pt});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_04_Eta", "SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_04_Eta", kTH1F, {Axis_eta});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_05_Phi", "SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_05_Phi", kTH1F, {Axis_phi});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_06_Rapidity", "SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_06_Rapidity", kTH1F, {Axis_rapidity});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_01_Mass", "SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_01_Mass", kTH1F, {Axis_PhiMass});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_02_P", "SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_02_P", kTH1F, {Axis_p});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_03_Pt", "SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_03_Pt", kTH1F, {Axis_pt});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_04_Eta", "SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_04_Eta", kTH1F, {Axis_eta});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_05_Phi", "SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_05_Phi", kTH1F, {Axis_phi});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_06_Rapidity", "SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_06_Rapidity", kTH1F, {Axis_rapidity});

    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_01_Mass", "SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_01_Mass", kTH1F, {Axis_PhiMass});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_02_P", "SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_02_P", kTH1F, {Axis_p});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_03_Pt", "SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_03_Pt", kTH1F, {Axis_pt});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_04_Eta", "SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_04_Eta", kTH1F, {Axis_eta});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_05_Phi", "SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_05_Phi", kTH1F, {Axis_phi});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_06_Rapidity", "SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_06_Rapidity", kTH1F, {Axis_rapidity});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_01_Mass", "SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_01_Mass", kTH1F, {Axis_PhiMass});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_02_P", "SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_02_P", kTH1F, {Axis_p});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_03_Pt", "SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_03_Pt", kTH1F, {Axis_pt});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_04_Eta", "SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_04_Eta", kTH1F, {Axis_eta});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_05_Phi", "SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_05_Phi", kTH1F, {Axis_phi});
    SE_recoPhi.add("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_06_Rapidity", "SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_06_Rapidity", kTH1F, {Axis_rapidity});
    //
    //
    ME_recoPhi.add("ME_hPhiPhi_Full_LeadPhi_01_Mass", "ME_hPhiPhi_Full_LeadPhi_01_Mass", kTH1F, {Axis_PhiMass});
    ME_recoPhi.add("ME_hPhiPhi_Full_LeadPhi_02_P", "ME_hPhiPhi_Full_LeadPhi_02_P", kTH1F, {Axis_p});
    ME_recoPhi.add("ME_hPhiPhi_Full_LeadPhi_03_Pt", "ME_hPhiPhi_Full_LeadPhi_03_Pt", kTH1F, {Axis_pt});
    ME_recoPhi.add("ME_hPhiPhi_Full_LeadPhi_04_Eta", "ME_hPhiPhi_Full_LeadPhi_04_Eta", kTH1F, {Axis_eta});
    ME_recoPhi.add("ME_hPhiPhi_Full_LeadPhi_05_Phi", "ME_hPhiPhi_Full_LeadPhi_05_Phi", kTH1F, {Axis_phi});
    ME_recoPhi.add("ME_hPhiPhi_Full_LeadPhi_06_Rapidity", "ME_hPhiPhi_Full_LeadPhi_06_Rapidity", kTH1F, {Axis_rapidity});
    ME_recoPhi.add("ME_hPhiPhi_Full_AssoPhi_01_Mass", "ME_hPhiPhi_Full_AssoPhi_01_Mass", kTH1F, {Axis_PhiMass});
    ME_recoPhi.add("ME_hPhiPhi_Full_AssoPhi_02_P", "ME_hPhiPhi_Full_AssoPhi_02_P", kTH1F, {Axis_p});
    ME_recoPhi.add("ME_hPhiPhi_Full_AssoPhi_03_Pt", "ME_hPhiPhi_Full_AssoPhi_03_Pt", kTH1F, {Axis_pt});
    ME_recoPhi.add("ME_hPhiPhi_Full_AssoPhi_04_Eta", "ME_hPhiPhi_Full_AssoPhi_04_Eta", kTH1F, {Axis_eta});
    ME_recoPhi.add("ME_hPhiPhi_Full_AssoPhi_05_Phi", "ME_hPhiPhi_Full_AssoPhi_05_Phi", kTH1F, {Axis_phi});
    ME_recoPhi.add("ME_hPhiPhi_Full_AssoPhi_06_Rapidity", "ME_hPhiPhi_Full_AssoPhi_06_Rapidity", kTH1F, {Axis_rapidity});

    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_Full_LeadPhi_01_Mass", "ME_hPhiPhi_4to8LeadPhi_Full_LeadPhi_01_Mass", kTH1F, {Axis_PhiMass});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_Full_LeadPhi_02_P", "ME_hPhiPhi_4to8LeadPhi_Full_LeadPhi_02_P", kTH1F, {Axis_p});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_Full_LeadPhi_03_Pt", "ME_hPhiPhi_4to8LeadPhi_Full_LeadPhi_03_Pt", kTH1F, {Axis_pt});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_Full_LeadPhi_04_Eta", "ME_hPhiPhi_4to8LeadPhi_Full_LeadPhi_04_Eta", kTH1F, {Axis_eta});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_Full_LeadPhi_05_Phi", "ME_hPhiPhi_4to8LeadPhi_Full_LeadPhi_05_Phi", kTH1F, {Axis_phi});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_Full_LeadPhi_06_Rapidity", "ME_hPhiPhi_4to8LeadPhi_Full_LeadPhi_06_Rapidity", kTH1F, {Axis_rapidity});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_Full_AssoPhi_01_Mass", "ME_hPhiPhi_4to8LeadPhi_Full_AssoPhi_01_Mass", kTH1F, {Axis_PhiMass});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_Full_AssoPhi_02_P", "ME_hPhiPhi_4to8LeadPhi_Full_AssoPhi_02_P", kTH1F, {Axis_p});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_Full_AssoPhi_03_Pt", "ME_hPhiPhi_4to8LeadPhi_Full_AssoPhi_03_Pt", kTH1F, {Axis_pt});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_Full_AssoPhi_04_Eta", "ME_hPhiPhi_4to8LeadPhi_Full_AssoPhi_04_Eta", kTH1F, {Axis_eta});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_Full_AssoPhi_05_Phi", "ME_hPhiPhi_4to8LeadPhi_Full_AssoPhi_05_Phi", kTH1F, {Axis_phi});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_Full_AssoPhi_06_Rapidity", "ME_hPhiPhi_4to8LeadPhi_Full_AssoPhi_06_Rapidity", kTH1F, {Axis_rapidity});

    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_01_Mass", "ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_01_Mass", kTH1F, {Axis_PhiMass});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_02_P", "ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_02_P", kTH1F, {Axis_p});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_03_Pt", "ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_03_Pt", kTH1F, {Axis_pt});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_04_Eta", "ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_04_Eta", kTH1F, {Axis_eta});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_05_Phi", "ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_05_Phi", kTH1F, {Axis_phi});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_06_Rapidity", "ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_06_Rapidity", kTH1F, {Axis_rapidity});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_01_Mass", "ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_01_Mass", kTH1F, {Axis_PhiMass});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_02_P", "ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_02_P", kTH1F, {Axis_p});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_03_Pt", "ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_03_Pt", kTH1F, {Axis_pt});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_04_Eta", "ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_04_Eta", kTH1F, {Axis_eta});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_05_Phi", "ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_05_Phi", kTH1F, {Axis_phi});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_06_Rapidity", "ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_06_Rapidity", kTH1F, {Axis_rapidity});

    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_01_Mass", "ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_01_Mass", kTH1F, {Axis_PhiMass});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_02_P", "ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_02_P", kTH1F, {Axis_p});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_03_Pt", "ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_03_Pt", kTH1F, {Axis_pt});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_04_Eta", "ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_04_Eta", kTH1F, {Axis_eta});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_05_Phi", "ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_05_Phi", kTH1F, {Axis_phi});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_06_Rapidity", "ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_06_Rapidity", kTH1F, {Axis_rapidity});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_01_Mass", "ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_01_Mass", kTH1F, {Axis_PhiMass});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_02_P", "ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_02_P", kTH1F, {Axis_p});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_03_Pt", "ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_03_Pt", kTH1F, {Axis_pt});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_04_Eta", "ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_04_Eta", kTH1F, {Axis_eta});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_05_Phi", "ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_05_Phi", kTH1F, {Axis_phi});
    ME_recoPhi.add("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_06_Rapidity", "ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_06_Rapidity", kTH1F, {Axis_rapidity});
    //

    // LeadingPhi-Phi Correlations
    SE_recoAnalysis.add("SE_hPhiPhi_Full_dPhi_dEta", "SE_hPhiPhi_Full_dPhi_dEta", kTH2F, {Axis_dPhi, Axis_dEta});
    SE_recoAnalysis.add("SE_hPhiPhi_Full_Analysis", "SE_hPhiPhi_Full_Analysis", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    SE_recoAnalysis.add("SE_hPhiPhi_4to8LeadPhi_Full_dPhi_dEta", "SE_hPhiPhi_4to8LeadPhi_Full_dPhi_dEta", kTH2F, {Axis_dPhi, Axis_dEta});
    SE_recoAnalysis.add("SE_hPhiPhi_4to8LeadPhi_Full_Analysis", "SE_hPhiPhi_4to8LeadPhi_Full_Analysis", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    SE_recoAnalysis.add("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_dPhi_dEta", "SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_dPhi_dEta", kTH2F, {Axis_dPhi, Axis_dEta});
    SE_recoAnalysis.add("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_Analysis", "SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_Analysis ", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    SE_recoAnalysis.add("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_dPhi_dEta", "SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_dPhi_dEta", kTH2F, {Axis_dPhi, Axis_dEta});
    SE_recoAnalysis.add("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_Analysis", "SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_Analysis", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});

    ME_recoAnalysis.add("ME_hPhiPhi_Full_dPhi_dEta", "ME_hPhiPhi_Full_dPhi_dEta", kTH2F, {Axis_dPhi, Axis_dEta});
    ME_recoAnalysis.add("ME_hPhiPhi_Full_Analysis", "ME_hPhiPhi_Full_Analysis", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    ME_recoAnalysis.add("ME_hPhiPhi_4to8LeadPhi_Full_dPhi_dEta", "ME_hPhiPhi_4to8LeadPhi_Full_dPhi_dEta", kTH2F, {Axis_dPhi, Axis_dEta});
    ME_recoAnalysis.add("ME_hPhiPhi_4to8LeadPhi_Full_Analysis", "ME_hPhiPhi_4to8LeadPhi_Full_Analysis", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    ME_recoAnalysis.add("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_dPhi_dEta", "ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_dPhi_dEta", kTH2F, {Axis_dPhi, Axis_dEta});
    ME_recoAnalysis.add("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_Analysis", "ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_Analysis ", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    ME_recoAnalysis.add("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_dPhi_dEta", "ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_dPhi_dEta", kTH2F, {Axis_dPhi, Axis_dEta});
    ME_recoAnalysis.add("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_Analysis", "ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_Analysis", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    //

    // di-hadron correlations
    // SE
    // hadron_pT in [0.0-2.0](Bulk)
    SE_recoAnalysis.add("SE_hAnalysis_hh_0To2AssoHadron_LeadH_02_P", "SE_hAnalysis_hh_0To2AssoHadron_LeadH_02_P", kTH1F, {Axis_p});
    SE_recoAnalysis.add("SE_hAnalysis_hh_0To2AssoHadron_LeadH_03_Pt", "SE_hAnalysis_hh_0To2AssoHadron_LeadH_03_Pt", kTH1F, {Axis_pt});
    SE_recoAnalysis.add("SE_hAnalysis_hh_0To2AssoHadron_LeadH_04_Eta", "SE_hAnalysis_hh_0To2AssoHadron_LeadH_04_Eta", kTH1F, {Axis_eta});
    SE_recoAnalysis.add("SE_hAnalysis_hh_0To2AssoHadron_LeadH_05_Phi", "SE_hAnalysis_hh_0To2AssoHadron_LeadH_05_Phi", kTH1F, {Axis_phi});
    SE_recoAnalysis.add("SE_hAnalysis_hh_0To2AssoHadron_AssoH_02_P", "SE_hAnalysis_hh_0To2AssoHadron_AssoH_02_P", kTH1F, {Axis_p});
    SE_recoAnalysis.add("SE_hAnalysis_hh_0To2AssoHadron_AssoH_03_Pt", "SE_hAnalysis_hh_0To2AssoHadron_AssoH_03_Pt", kTH1F, {Axis_pt});
    SE_recoAnalysis.add("SE_hAnalysis_hh_0To2AssoHadron_AssoH_04_Eta", "SE_hAnalysis_hh_0To2AssoHadron_AssoH_04_Eta", kTH1F, {Axis_eta});
    SE_recoAnalysis.add("SE_hAnalysis_hh_0To2AssoHadron_AssoH_05_Phi", "SE_hAnalysis_hh_0To2AssoHadron_AssoH_05_Phi", kTH1F, {Axis_phi});
    SE_recoAnalysis.add("SE_hAnalysis_hh_0To2AssoHadron_dPhi_dEta", "SE_hAnalysis_hh_0To2AssoHadron_dPhi_dEta(Bulk);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    SE_recoAnalysis.add("SE_hAnalysis_hh_0To2AssoHadron_Analysis", "SE_hAnalysis_hh_0To2AssoHadron_Analysis", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // hadron_pT in [2.0-4.0](Required)
    SE_recoAnalysis.add("SE_hAnalysis_hh_2To4AssoHadron_LeadH_02_P", "SE_hAnalysis_hh_2To4AssoHadron_LeadH_02_P", kTH1F, {Axis_p});
    SE_recoAnalysis.add("SE_hAnalysis_hh_2To4AssoHadron_LeadH_03_Pt", "SE_hAnalysis_hh_2To4AssoHadron_LeadH_03_Pt", kTH1F, {Axis_pt});
    SE_recoAnalysis.add("SE_hAnalysis_hh_2To4AssoHadron_LeadH_04_Eta", "SE_hAnalysis_hh_2To4AssoHadron_LeadH_04_Eta", kTH1F, {Axis_eta});
    SE_recoAnalysis.add("SE_hAnalysis_hh_2To4AssoHadron_LeadH_05_Phi", "SE_hAnalysis_hh_2To4AssoHadron_LeadH_05_Phi", kTH1F, {Axis_phi});
    SE_recoAnalysis.add("SE_hAnalysis_hh_2To4AssoHadron_AssoH_02_P", "SE_hAnalysis_hh_2To4AssoHadron_AssoH_02_P", kTH1F, {Axis_p});
    SE_recoAnalysis.add("SE_hAnalysis_hh_2To4AssoHadron_AssoH_03_Pt", "SE_hAnalysis_hh_2To4AssoHadron_AssoH_03_Pt", kTH1F, {Axis_pt});
    SE_recoAnalysis.add("SE_hAnalysis_hh_2To4AssoHadron_AssoH_04_Eta", "SE_hAnalysis_hh_2To4AssoHadron_AssoH_04_Eta", kTH1F, {Axis_eta});
    SE_recoAnalysis.add("SE_hAnalysis_hh_2To4AssoHadron_AssoH_05_Phi", "SE_hAnalysis_hh_2To4AssoHadron_AssoH_05_Phi", kTH1F, {Axis_phi});
    SE_recoAnalysis.add("SE_hAnalysis_hh_2To4AssoHadron_dPhi_dEta", "SE_hAnalysis_hh_2To4AssoHadron_dPhi_dEta(Required);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    SE_recoAnalysis.add("SE_hAnalysis_hh_2To4AssoHadron_Analysis", "SE_hAnalysis_hh_2To4AssoHadron_Analysis", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // SE
    // ME
    // hadron_pT in [0.0-2.0](Bulk)
    ME_recoAnalysis.add("ME_hAnalysis_hh_0To2AssoHadron_LeadH_02_P", "ME_hAnalysis_hh_0To2AssoHadron_LeadH_02_P", kTH1F, {Axis_p});
    ME_recoAnalysis.add("ME_hAnalysis_hh_0To2AssoHadron_LeadH_03_Pt", "ME_hAnalysis_hh_0To2AssoHadron_LeadH_03_Pt", kTH1F, {Axis_pt});
    ME_recoAnalysis.add("ME_hAnalysis_hh_0To2AssoHadron_LeadH_04_Eta", "ME_hAnalysis_hh_0To2AssoHadron_LeadH_04_Eta", kTH1F, {Axis_eta});
    ME_recoAnalysis.add("ME_hAnalysis_hh_0To2AssoHadron_LeadH_05_Phi", "ME_hAnalysis_hh_0To2AssoHadron_LeadH_05_Phi", kTH1F, {Axis_phi});
    ME_recoAnalysis.add("ME_hAnalysis_hh_0To2AssoHadron_AssoH_02_P", "ME_hAnalysis_hh_0To2AssoHadron_AssoH_02_P", kTH1F, {Axis_p});
    ME_recoAnalysis.add("ME_hAnalysis_hh_0To2AssoHadron_AssoH_03_Pt", "ME_hAnalysis_hh_0To2AssoHadron_AssoH_03_Pt", kTH1F, {Axis_pt});
    ME_recoAnalysis.add("ME_hAnalysis_hh_0To2AssoHadron_AssoH_04_Eta", "ME_hAnalysis_hh_0To2AssoHadron_AssoH_04_Eta", kTH1F, {Axis_eta});
    ME_recoAnalysis.add("ME_hAnalysis_hh_0To2AssoHadron_AssoH_05_Phi", "ME_hAnalysis_hh_0To2AssoHadron_AssoH_05_Phi", kTH1F, {Axis_phi});
    ME_recoAnalysis.add("ME_hAnalysis_hh_0To2AssoHadron_dPhi_dEta", "ME_hAnalysis_hh_0To2AssoHadron_dPhi_dEta(Bulk);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    ME_recoAnalysis.add("ME_hAnalysis_hh_0To2AssoHadron_Analysis", "ME_hAnalysis_hh_0To2AssoHadron_Analysis", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    // hadron_pT in [2.0-4.0](Required)
    ME_recoAnalysis.add("ME_hAnalysis_hh_2To4AssoHadron_LeadH_02_P", "ME_hAnalysis_hh_2To4AssoHadron_LeadH_02_P", kTH1F, {Axis_p});
    ME_recoAnalysis.add("ME_hAnalysis_hh_2To4AssoHadron_LeadH_03_Pt", "ME_hAnalysis_hh_2To4AssoHadron_LeadH_03_Pt", kTH1F, {Axis_pt});
    ME_recoAnalysis.add("ME_hAnalysis_hh_2To4AssoHadron_LeadH_04_Eta", "ME_hAnalysis_hh_2To4AssoHadron_LeadH_04_Eta", kTH1F, {Axis_eta});
    ME_recoAnalysis.add("ME_hAnalysis_hh_2To4AssoHadron_LeadH_05_Phi", "ME_hAnalysis_hh_2To4AssoHadron_LeadH_05_Phi", kTH1F, {Axis_phi});
    ME_recoAnalysis.add("ME_hAnalysis_hh_2To4AssoHadron_AssoH_02_P", "ME_hAnalysis_hh_2To4AssoHadron_AssoH_02_P", kTH1F, {Axis_p});
    ME_recoAnalysis.add("ME_hAnalysis_hh_2To4AssoHadron_AssoH_03_Pt", "ME_hAnalysis_hh_2To4AssoHadron_AssoH_03_Pt", kTH1F, {Axis_pt});
    ME_recoAnalysis.add("ME_hAnalysis_hh_2To4AssoHadron_AssoH_04_Eta", "ME_hAnalysis_hh_2To4AssoHadron_AssoH_04_Eta", kTH1F, {Axis_eta});
    ME_recoAnalysis.add("ME_hAnalysis_hh_2To4AssoHadron_AssoH_05_Phi", "ME_hAnalysis_hh_2To4AssoHadron_AssoH_05_Phi", kTH1F, {Axis_phi});
    ME_recoAnalysis.add("ME_hAnalysis_hh_2To4AssoHadron_dPhi_dEta", "ME_hAnalysis_hh_2To4AssoHadron_dPhi_dEta(Required);#Delta#phi;#Delta#eta", kTH2F, {Axis_dPhi, Axis_dEta});
    ME_recoAnalysis.add("ME_hAnalysis_hh_2To4AssoHadron_Analysis", "ME_hAnalysis_hh_2To4AssoHadron_Analysis", kTHnSparseF, {Axis_CentBins, Axis_vtxZbins, Axis_dPhi, Axis_dEta});
    //

    // Mixing DataFrame information
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBinPerDF_Trig", "ME_hMixingEventsAvailabePerBinPerDF_Trig", kTH2F, {Axis_MixingBin, Axis_EventCount});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBinPerDF_nPhi", "ME_hMixingEventsAvailabePerBinPerDF_nPhi", kTH2F, {Axis_MixingBin, Axis_EventCount});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBin00", "ME_hMixingEventsAvailabePerBin00", kTH1F, {Axis_MixingBin});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBin01", "ME_hMixingEventsAvailabePerBin01", kTH1F, {Axis_MixingBin});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBin02", "ME_hMixingEventsAvailabePerBin02", kTH1F, {Axis_MixingBin});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBin03", "ME_hMixingEventsAvailabePerBin03", kTH1F, {Axis_MixingBin});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBin04", "ME_hMixingEventsAvailabePerBin04", kTH1F, {Axis_MixingBin});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBin05", "ME_hMixingEventsAvailabePerBin05", kTH1F, {Axis_MixingBin});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBin06", "ME_hMixingEventsAvailabePerBin06", kTH1F, {Axis_MixingBin});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBin07", "ME_hMixingEventsAvailabePerBin07", kTH1F, {Axis_MixingBin});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBin08", "ME_hMixingEventsAvailabePerBin08", kTH1F, {Axis_MixingBin});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBin09", "ME_hMixingEventsAvailabePerBin09", kTH1F, {Axis_MixingBin});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBin10", "ME_hMixingEventsAvailabePerBin10", kTH1F, {Axis_MixingBin});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBin11", "ME_hMixingEventsAvailabePerBin11", kTH1F, {Axis_MixingBin});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBin12", "ME_hMixingEventsAvailabePerBin12", kTH1F, {Axis_MixingBin});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBin13", "ME_hMixingEventsAvailabePerBin13", kTH1F, {Axis_MixingBin});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBin14", "ME_hMixingEventsAvailabePerBin14", kTH1F, {Axis_MixingBin});

    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBinPerDF00", "ME_hMixingEventsAvailabePerBinPerDF00", kTH2F, {Axis_MixingBin, Axis_EventCount});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBinPerDF01", "ME_hMixingEventsAvailabePerBinPerDF01", kTH2F, {Axis_MixingBin, Axis_EventCount});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBinPerDF02", "ME_hMixingEventsAvailabePerBinPerDF02", kTH2F, {Axis_MixingBin, Axis_EventCount});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBinPerDF03", "ME_hMixingEventsAvailabePerBinPerDF03", kTH2F, {Axis_MixingBin, Axis_EventCount});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBinPerDF04", "ME_hMixingEventsAvailabePerBinPerDF04", kTH2F, {Axis_MixingBin, Axis_EventCount});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBinPerDF05", "ME_hMixingEventsAvailabePerBinPerDF05", kTH2F, {Axis_MixingBin, Axis_EventCount});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBinPerDF06", "ME_hMixingEventsAvailabePerBinPerDF06", kTH2F, {Axis_MixingBin, Axis_EventCount});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBinPerDF07", "ME_hMixingEventsAvailabePerBinPerDF07", kTH2F, {Axis_MixingBin, Axis_EventCount});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBinPerDF08", "ME_hMixingEventsAvailabePerBinPerDF08", kTH2F, {Axis_MixingBin, Axis_EventCount});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBinPerDF09", "ME_hMixingEventsAvailabePerBinPerDF09", kTH2F, {Axis_MixingBin, Axis_EventCount});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBinPerDF10", "ME_hMixingEventsAvailabePerBinPerDF10", kTH2F, {Axis_MixingBin, Axis_EventCount});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBinPerDF11", "ME_hMixingEventsAvailabePerBinPerDF11", kTH2F, {Axis_MixingBin, Axis_EventCount});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBinPerDF12", "ME_hMixingEventsAvailabePerBinPerDF12", kTH2F, {Axis_MixingBin, Axis_EventCount});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBinPerDF13", "ME_hMixingEventsAvailabePerBinPerDF13", kTH2F, {Axis_MixingBin, Axis_EventCount});
    ME_recoAnalysis.add("ME_hMixingEventsAvailabePerBinPerDF14", "ME_hMixingEventsAvailabePerBinPerDF14", kTH2F, {Axis_MixingBin, Axis_EventCount});
    //

    // Mixing Structures
    for (int iMixCase = 0; iMixCase < nMixCases; iMixCase++) {
      for (int iBin = 0; iBin < nMixBin; iBin++) {
        mixingSlotPos[iMixCase][iBin] = -1;
        collRoll[iMixCase][iBin] = -1;
      }
    }

    std::vector<int> tempList;
    for (int iMixPosition = 0; iMixPosition < (nMixEvt + 1); iMixPosition++) {
      tempList.clear();
      for (int j = 0; j < (nMixEvt + 1); j++) {
        if (j == iMixPosition) {
          continue;
        }
        tempList.push_back(j);
      }
      mixSlotPosList.push_back(tempList);
    }
    // mixSlotPosList
    //

    // Printing the Stored Registry information
    LOG(info) << "Printing Stored Registry Information";
    LOG(info) << "Printing SE_recoEvent ";
    SE_recoEvent.print();
    LOG(info) << "Printing SE_recoTracks ";
    SE_recoTracks.print();
    LOG(info) << "Printing SE_recoKaon ";
    SE_recoKaon.print();
    LOG(info) << "Printing SE_recoPhi ";
    SE_recoPhi.print();
    LOG(info) << "Printing SE_recoAnalysis ";
    SE_recoAnalysis.print();
  }

  template <typename T>
  bool selectionTrack(const T& track)
  {
    if (track.tpcNClsCrossedRows() < 70) {
      return false;
    }
    if (!track.isGlobalTrack()) {
      return false;
    }
    return true;
  }

  // tpc Selections
  template <typename T>
  bool selPionTPCInnerParam(T track)
  {
    if (TMath::Abs(track.tpcNSigmaEl()) > 3.0 && TMath::Abs(track.tpcNSigmaKa()) > 3.0 && TMath::Abs(track.tpcNSigmaPr()) > 3.0 && TMath::Abs(track.tpcNSigmaDe()) > 3.0) {
      if (0.05 <= track.tpcInnerParam() && track.tpcInnerParam() < 0.70 && TMath::Abs(track.tpcNSigmaPi()) < 3.0) {
        return true;
      }
      if (0.70 <= track.tpcInnerParam() && TMath::Abs(track.tpcNSigmaPi()) < 2.0) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selKaonTPCInnerParam(T track)
  {
    // p dependent cuts
    if (TMath::Abs(track.tpcNSigmaEl()) > 3.0 && TMath::Abs(track.tpcNSigmaPi()) > 3.0 && TMath::Abs(track.tpcNSigmaPr()) > 3.0 && TMath::Abs(track.tpcNSigmaDe()) > 3.0) {
      if (0.05 <= track.tpcInnerParam() && track.tpcInnerParam() < 0.70 && TMath::Abs(track.tpcNSigmaKa()) < 3.0) {
        return true;
      }
      if (0.70 <= track.tpcInnerParam() && TMath::Abs(track.tpcNSigmaKa()) < 2.0) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selProtonTPCInnerParam(T track)
  {
    if (TMath::Abs(track.tpcNSigmaEl()) > 3.0 && TMath::Abs(track.tpcNSigmaPi()) > 3.0 && TMath::Abs(track.tpcNSigmaKa()) > 3.0 && TMath::Abs(track.tpcNSigmaDe()) > 3.0) {
      if (0.05 <= track.tpcInnerParam() && track.tpcInnerParam() < 1.60 && TMath::Abs(track.tpcNSigmaPr()) < 3.0) {
        return true;
      }
      if (1.60 <= track.tpcInnerParam() && TMath::Abs(track.tpcNSigmaPr()) < 2.0) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selDeuteronTPCInnerParam(T track)
  {
    if (TMath::Abs(track.tpcNSigmaEl()) > 3.0 && TMath::Abs(track.tpcNSigmaPi()) > 3.0 && TMath::Abs(track.tpcNSigmaKa()) > 3.0 && TMath::Abs(track.tpcNSigmaPr()) > 3.0) {
      if (0.05 <= track.tpcInnerParam() && track.tpcInnerParam() < 1.80 && TMath::Abs(track.tpcNSigmaDe()) < 3.0) {
        return true;
      }
      if (1.80 <= track.tpcInnerParam() && TMath::Abs(track.tpcNSigmaDe()) < 2.0) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selElectronTPCInnerParam(T track)
  {
    if (track.tpcNSigmaEl() < 3.0 && track.tpcNSigmaPi() > 3.0 && track.tpcNSigmaKa() > 3.0 && track.tpcNSigmaPr() > 3.0 && track.tpcNSigmaDe() > 3.0) {
      return true;
    }
    return false;
  }
  //

  // TOF Selections
  // Pion
  template <typename T>
  bool selPionTOF(T track)
  {
    if (track.p() <= 0.75
        // && (TMath::Power(track.tpcNSigmaPi(),2)+TMath::Power(track.tofNSigmaPi(),2)) < 9.0
        && TMath::Abs(track.tpcNSigmaPi()) < 3.0 && TMath::Abs(track.tofNSigmaPi()) < 3.0 && TMath::Abs(track.tofNSigmaEl()) > 3.0) {
      return true;
    } else if (0.75 < track.p() // after p = 0.75, Pi and Ka lines of nSigma 3.0 will start intersecting
                                //  && (TMath::Power(track.tpcNSigmaPi(),2)+TMath::Power(track.tofNSigmaPi(),2)) < 8.0
               && TMath::Abs(track.tpcNSigmaPi()) < 2.0 && TMath::Abs(track.tofNSigmaPi()) < 2.0 && TMath::Abs(track.tofNSigmaEl()) > 3.0 && TMath::Abs(track.tofNSigmaKa()) > 3.0 && TMath::Abs(track.tofNSigmaPr()) > 3.0 && TMath::Abs(track.tofNSigmaDe()) > 3.0) {
      return true;
    }
    return false;
  }

  // Kaon
  template <typename T>
  bool selKaonTOF(T track)
  {
    if (track.p() <= 0.75
        // && (TMath::Power(track.tpcNSigmaKa(),2)+TMath::Power(track.tofNSigmaKa(),2)) < 9.0
        && TMath::Abs(track.tpcNSigmaKa()) < 3.0 && TMath::Abs(track.tofNSigmaKa()) < 3.0) {
      return true;
    }
    if (0.75 < track.p() && track.p() <= 1.30 // after 0.75 Pi and Ka lines of nSigma 3.0 will start intersecting
        && TMath::Abs(track.tpcNSigmaKa()) < 3.0 && TMath::Abs(track.tofNSigmaKa()) < 3.0 && TMath::Abs(track.tofNSigmaPi()) > 3.0 && TMath::Abs(track.tofNSigmaEl()) > 3.0) {
      return true;
    }
    if (1.30 < track.p() // after 1.30 Pr and Ka lines of nSigma 3.0 will start intersecting
                         //  && (TMath::Power(track.tpcNSigmaKa(),2)+TMath::Power(track.tofNSigmaKa(),2)) < 8.0
        && TMath::Abs(track.tpcNSigmaKa()) < 2.0 && TMath::Abs(track.tofNSigmaKa()) < 2.0 && TMath::Abs(track.tofNSigmaEl()) > 3.0 && TMath::Abs(track.tofNSigmaPi()) > 3.0 && TMath::Abs(track.tofNSigmaPr()) > 3.0 && TMath::Abs(track.tofNSigmaDe()) > 3.0) {
      return true;
    }
    return false;
  }

  // Proton
  template <typename T>
  bool selProtonTOF(T track)
  {
    if (track.p() <= 1.30
        // && (TMath::Power(track.tpcNSigmaPr(),2)+TMath::Power(track.tofNSigmaPr(),2)) < 9.0
        && TMath::Abs(track.tpcNSigmaPr()) < 3.0 && TMath::Abs(track.tofNSigmaPr()) < 3.0) {
      return true;
    }
    if (1.30 < track.p() && track.p() <= 3.10                                                                                                                                                                                                                 // after 1.30 Pr and Ka lines of nSigma 3.0 will start intersecting
        && TMath::Abs(track.tpcNSigmaPr()) < 3.0 && TMath::Abs(track.tofNSigmaPr()) < 3.0 && TMath::Abs(track.tofNSigmaEl()) > 3.0 && TMath::Abs(track.tofNSigmaPi()) > 3.0 && TMath::Abs(track.tofNSigmaKa()) > 3.0 && TMath::Abs(track.tofNSigmaDe()) > 3.0 // Some Deuteron contamination is still coming in p dependent cuts
    ) {
      return true;
    }
    if (3.10 < track.p() // after 3.10 Pr and De lines of nSigma 3.0 will start intersecting
                         //  && (TMath::Power(track.tpcNSigmaPr(),2)+TMath::Power(track.tofNSigmaPr(),2)) < 8.0
        && TMath::Abs(track.tpcNSigmaPr()) < 2.0 && TMath::Abs(track.tofNSigmaPr()) < 2.0 && TMath::Abs(track.tofNSigmaEl()) > 3.0 && TMath::Abs(track.tofNSigmaPi()) > 3.0 && TMath::Abs(track.tofNSigmaKa()) > 3.0 && TMath::Abs(track.tofNSigmaDe()) > 3.0) {
      return true;
    }
    return false;
  }

  // Deuteron
  template <typename T>
  bool selDeuteronTOF(T track)
  {
    if (track.p() <= 3.10
        // && (TMath::Power(track.tpcNSigmaDe(),2)+TMath::Power(track.tofNSigmaDe(),2)) < 9.0
        && TMath::Abs(track.tpcNSigmaDe()) < 3.0 && TMath::Abs(track.tofNSigmaDe()) < 3.0 && TMath::Abs(track.tofNSigmaEl()) > 3.0 && TMath::Abs(track.tofNSigmaPi()) > 3.0 && TMath::Abs(track.tofNSigmaKa()) > 3.0 && TMath::Abs(track.tofNSigmaPr()) > 3.0) {
      return true;
    }
    if (3.10 < track.p()                                                                                                                                                                                                                                      // after 3.10 De and Pr lines of nSigma 3.0 will start intersecting
        && TMath::Abs(track.tpcNSigmaDe()) < 2.0 && TMath::Abs(track.tofNSigmaDe()) < 2.0 && TMath::Abs(track.tofNSigmaEl()) > 3.0 && TMath::Abs(track.tofNSigmaPi()) > 3.0 && TMath::Abs(track.tofNSigmaKa()) > 3.0 && TMath::Abs(track.tofNSigmaPr()) > 3.0 // Some Deuteron contamination is still coming in p dependent cuts
    ) {
      return true;
    }
    return false;
  }

  // Electron
  template <typename T>
  bool selElectronTOF(T track)
  {
    if (
      (TMath::Power(track.tpcNSigmaEl(), 2) + TMath::Power(track.tofNSigmaEl(), 2)) < 9.00 && TMath::Abs(track.tofNSigmaPi()) > 3.0 && TMath::Abs(track.tofNSigmaKa()) > 3.0 && TMath::Abs(track.tofNSigmaPr()) > 3.0 && TMath::Abs(track.tofNSigmaDe()) > 3.0) {
      return true;
    }
    return false;
  }
  //

  // SelectionFunctions
  // Pion
  template <typename T>
  bool selPion(T track, int& IdMethod)
  {
    if (!track.hasTOF() && selPionTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() < 0.0 && selPionTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && selPionTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && !selPionTPCInnerParam(track)) {
      IdMethod = 1;
      return selPionTOF(track);
    }
    return false;
  }

  // Kaon
  template <typename T>
  bool selKaon(T track, int& IdMethod)
  {
    if (!track.hasTOF() && selKaonTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() < 0.0 && selKaonTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && selKaonTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && !selKaonTPCInnerParam(track)) {
      IdMethod = 1;
      return selKaonTOF(track);
    }
    return false;
  }

  // Proton
  template <typename T>
  bool selProton(T track, int& IdMethod)
  {
    if (!track.hasTOF() && selProtonTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() < 0.0 && selProtonTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && selProtonTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && !selProtonTPCInnerParam(track)) {
      IdMethod = 1;
      return selProtonTOF(track);
    }
    return false;
  }

  // Deuteron
  template <typename T>
  bool selDeuteron(T track, int& IdMethod)
  {
    if (!track.hasTOF() && selDeuteronTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() < 0.0 && selDeuteronTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && selDeuteronTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && !selDeuteronTPCInnerParam(track)) {
      IdMethod = 1;
      return selDeuteronTOF(track);
    }
    return false;
  }

  // Electron
  template <typename T>
  bool selElectron(T track, int& IdMethod)
  {
    if (!track.hasTOF() && selElectronTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() < 0.0 && selElectronTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && selElectronTPCInnerParam(track)) {
      IdMethod = 0;
      return true;
    }
    if (track.hasTOF() && track.beta() > 0.0 && !selElectronTPCInnerParam(track)) {
      IdMethod = 1;
      return selElectronTOF(track);
    }
    return false;
  }
  //

  enum pidTagValue {
    tagOther = 0,
    tagPion = 1,
    tagKaon = 2,
    tagProton = 4,
    tagElectron = 8,
    tagDeuteron = 16
  };

  template <typename T>
  int FindTrackTag(T track)
  {
    int tempPID = tagOther;
    // 0- Some other particle
    // 1- Pion
    // 2- Kaon
    // 4- Proton
    // 8- Electron
    int PiIdMethod = -1;
    int KaIdMethod = -1;
    int PrIdMethod = -1;
    int ElIdMethod = -1;
    int DeIdMethod = -1;
    if (selPion(track, PiIdMethod)) {
      tempPID += tagPion;
    }
    if (selKaon(track, KaIdMethod)) {
      tempPID += tagKaon;
    }
    if (selProton(track, PrIdMethod)) {
      tempPID += tagProton;
    }
    if (selElectron(track, ElIdMethod)) {
      tempPID += tagElectron;
    }
    if (selDeuteron(track, DeIdMethod)) {
      tempPID += tagDeuteron;
    }

    return tempPID;
  }

  template <typename T>
  bool matchIndices(T track1, T track2)
  {
    if (track1.globalIndex() == track2.globalIndex()) {
      return true;
    }
    return false;
  }

  void InsertionSortVector(std::vector<int64_t>& UnsortedVector)
  {
    for (long unsigned int i = 1; i < UnsortedVector.size(); i++) {
      int currentElement = UnsortedVector[i]; // Element to be Inserted at correct position
      int j;                                  //(j+1) is the correct position of current element
      for (j = i - 1; j >= 0 && (UnsortedVector[j] > currentElement); j--) {
        UnsortedVector[j + 1] = UnsortedVector[j];
      }
      UnsortedVector[j + 1] = currentElement;
    }
  }

  template <typename int64_t>
  int BinarySearchVector(int64_t Key, std::vector<int64_t> List, int low, int high)
  {
    while (low <= high) {
      int mid = low + (high - low) / 2;
      if (Key == List[mid]) {
        return mid;
      }

      if (Key > List[mid]) {
        low = mid + 1;
      } // If Key is greater, ignore left  half, update the low
      else {
        high = mid - 1;
      } // If Key is smaller, ignore right half, update the high
    }
    return -1; // Element is not present
  }
  // do a fastest search in an sorted array. ==> Binary Search in an array.
  template <typename T>
  bool checkTrackInList(T track, std::vector<int64_t> ParticleList)
  {
    // for( auto indexVal : ParticleList){ if (track.globalIndex() == indexVal) {return true;}}
    // int posInVector = BinarySearchVector(track.globalIndex(), ParticleList, 0, ParticleList.size()-1);
    if (BinarySearchVector(track.globalIndex(), ParticleList, 0, ParticleList.size() - 1) != -1) {
      return true;
    }
    return false;
  }

  template <typename T>
  void FindRepeatEntries(std::vector<int64_t> ParticleList, T hist)
  {
    for (long unsigned int ii = 0; ii < ParticleList.size(); ii++) {
      int nCommonCount = 0; // checking the repeat number of track
      for (long unsigned int jj = 0; jj < ParticleList.size(); jj++) {
        if (ParticleList[jj] == ParticleList[ii]) {
          if (jj < ii) {
            break;
          } // break if it was already counted
          nCommonCount++; // To Calculate no of times the entry was repeated
        }
      }
      hist->Fill(nCommonCount);
    }
  }

  void FillNewListFromOldList(std::vector<int64_t>& NewList, std::vector<int64_t> OldList)
  {
    for (long unsigned int ii = 0; ii < OldList.size(); ii++) {
      bool RepeatEntry = false;
      for (long unsigned int jj = 0; jj < NewList.size(); jj++) {
        if (OldList[ii] == NewList[jj]) {
          RepeatEntry = true;
        }
      }
      if (!RepeatEntry) {
        NewList.push_back(OldList[ii]);
      }
    }
  }

  template <typename T>
  void FillFullTrackQA(T track)
  {
    // FullTrack Information
    SE_recoTracks.fill(HIST("SE_htracks_00_1_0_FullTrack_P"), track.p());
    SE_recoTracks.fill(HIST("SE_htracks_00_1_0_FullTrack_tpcInnerParam"), track.tpcInnerParam());
    SE_recoTracks.fill(HIST("SE_htracks_00_1_0_FullTrack_tofExpMom"), track.tofExpMom());

    SE_recoTracks.fill(HIST("SE_hTracks_00_1_1_FullTrack_Pt"), track.pt());
    SE_recoTracks.fill(HIST("SE_hTracks_00_1_2_FullTrack_Eta"), track.eta());
    SE_recoTracks.fill(HIST("SE_hTracks_00_1_3_FullTrack_Phi"), track.phi());
    SE_recoTracks.fill(HIST("SE_hTracks_00_1_4_FullTrack_DcaXY"), track.dcaXY());
    SE_recoTracks.fill(HIST("SE_hTracks_00_1_5_FullTrack_DcaZ"), track.dcaZ());
    SE_recoTracks.fill(HIST("SE_hTracks_00_1_6_FullTrack_Sign"), track.sign());

    // DcaXY
    SE_recoTracks.fill(HIST("SE_htracks_00_1_7_FullTrack_00_DcaXY"), track.p(), track.dcaXY());
    SE_recoTracks.fill(HIST("SE_htracks_00_1_7_FullTrack_01_DcaXY"), track.pt(), track.dcaXY());
    SE_recoTracks.fill(HIST("SE_htracks_00_1_7_FullTrack_02_DcaXY"), track.tpcInnerParam(), track.dcaXY());
    SE_recoTracks.fill(HIST("SE_htracks_00_1_7_FullTrack_03_DcaXY"), track.tofExpMom(), track.dcaXY());

    // DcaZ
    SE_recoTracks.fill(HIST("SE_htracks_00_1_7_FullTrack_00_DcaZ"), track.p(), track.dcaZ());
    SE_recoTracks.fill(HIST("SE_htracks_00_1_7_FullTrack_01_DcaZ"), track.pt(), track.dcaZ());
    SE_recoTracks.fill(HIST("SE_htracks_00_1_7_FullTrack_02_DcaZ"), track.tpcInnerParam(), track.dcaZ());
    SE_recoTracks.fill(HIST("SE_htracks_00_1_7_FullTrack_03_DcaZ"), track.tofExpMom(), track.dcaZ());

    // momemtum
    SE_recoTracks.fill(HIST("SE_hTracks_00_2_1_FullTrack_Axis_00_01"), track.p(), track.pt());            // Axis_pt            ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_2_2_FullTrack_Axis_00_02"), track.p(), track.tpcInnerParam()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_2_3_FullTrack_Axis_00_03"), track.p(), track.tofExpMom());     // Axis_tofExpMom     ;

    // tpcSignal
    SE_recoTracks.fill(HIST("SE_hTracks_00_3_1_FullTrack_Axis_00_05"), track.p(), track.tpcSignal());             // Axis_p             ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_3_2_FullTrack_Axis_02_05"), track.tpcInnerParam(), track.tpcSignal()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_3_3_FullTrack_Axis_03_05"), track.tofExpMom(), track.tpcSignal());     // Axis_tofExpMom     ;

    // tofBeta
    SE_recoTracks.fill(HIST("SE_hTracks_00_4_1_FullTrack_Axis_00_06"), track.p(), track.beta());             // Axis_p             ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_4_2_FullTrack_Axis_02_06"), track.tpcInnerParam(), track.beta()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_4_3_FullTrack_Axis_03_06"), track.tofExpMom(), track.beta());     // Axis_tofExpMom     ;

    // Look at Pion
    SE_recoTracks.fill(HIST("SE_hTracks_00_5_1_FullTrack_Axis_00_20"), track.p(), track.tpcNSigmaPi());             // Axis_p             ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_5_2_FullTrack_Axis_01_20"), track.pt(), track.tpcNSigmaPi());            // Axis_pt            ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_5_3_FullTrack_Axis_02_20"), track.tpcInnerParam(), track.tpcNSigmaPi()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_5_4_FullTrack_Axis_03_20"), track.tofExpMom(), track.tpcNSigmaPi());     // Axis_tofExpMom     ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_5_5_FullTrack_Axis_00_21"), track.p(), track.tofNSigmaPi());             // Axis_p             ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_5_6_FullTrack_Axis_01_21"), track.pt(), track.tofNSigmaPi());            // Axis_pt            ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_5_7_FullTrack_Axis_02_21"), track.tpcInnerParam(), track.tofNSigmaPi()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_5_8_FullTrack_Axis_03_21"), track.tofExpMom(), track.tofNSigmaPi());     // Axis_tofExpMom     ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_5_9_FullTrack_Axis_20_21"), track.tpcNSigmaPi(), track.tofNSigmaPi());   // Axis_tpcInnerParam ;
    // Pion
    // Look at Kaon
    SE_recoTracks.fill(HIST("SE_hTracks_00_6_1_FullTrack_Axis_00_22"), track.p(), track.tpcNSigmaKa());             // Axis_p             ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_6_2_FullTrack_Axis_01_22"), track.pt(), track.tpcNSigmaKa());            // Axis_pt            ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_6_3_FullTrack_Axis_02_22"), track.tpcInnerParam(), track.tpcNSigmaKa()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_6_4_FullTrack_Axis_03_22"), track.tofExpMom(), track.tpcNSigmaKa());     // Axis_tofExpMom     ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_6_5_FullTrack_Axis_00_23"), track.p(), track.tofNSigmaKa());             // Axis_p             ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_6_6_FullTrack_Axis_01_23"), track.pt(), track.tofNSigmaKa());            // Axis_pt            ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_6_7_FullTrack_Axis_02_23"), track.tpcInnerParam(), track.tofNSigmaKa()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_6_8_FullTrack_Axis_03_23"), track.tofExpMom(), track.tofNSigmaKa());     // Axis_tofExpMom     ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_6_9_FullTrack_Axis_22_23"), track.tpcNSigmaKa(), track.tofNSigmaKa());   // Axis_tpcInnerParam ;
    // Kaon
    // Look at Proton
    SE_recoTracks.fill(HIST("SE_hTracks_00_7_1_FullTrack_Axis_00_24"), track.p(), track.tpcNSigmaPr());             // Axis_p             ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_7_2_FullTrack_Axis_01_24"), track.pt(), track.tpcNSigmaPr());            // Axis_pt            ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_7_3_FullTrack_Axis_02_24"), track.tpcInnerParam(), track.tpcNSigmaPr()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_7_4_FullTrack_Axis_03_24"), track.tofExpMom(), track.tpcNSigmaPr());     // Axis_tofExpMom     ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_7_5_FullTrack_Axis_00_25"), track.p(), track.tofNSigmaPr());             // Axis_p             ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_7_6_FullTrack_Axis_01_25"), track.pt(), track.tofNSigmaPr());            // Axis_pt            ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_7_7_FullTrack_Axis_02_25"), track.tpcInnerParam(), track.tofNSigmaPr()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_7_8_FullTrack_Axis_03_25"), track.tofExpMom(), track.tofNSigmaPr());     // Axis_tofExpMom     ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_7_9_FullTrack_Axis_24_25"), track.tpcNSigmaPr(), track.tofNSigmaPr());   // Axis_tpcInnerParam ;
    // Proton
    // Look at Electron
    SE_recoTracks.fill(HIST("SE_hTracks_00_8_1_FullTrack_Axis_00_26"), track.p(), track.tpcNSigmaEl());             // Axis_p             ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_8_2_FullTrack_Axis_01_26"), track.pt(), track.tpcNSigmaEl());            // Axis_pt            ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_8_3_FullTrack_Axis_02_26"), track.tpcInnerParam(), track.tpcNSigmaEl()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_8_4_FullTrack_Axis_03_26"), track.tofExpMom(), track.tpcNSigmaEl());     // Axis_tofExpMom     ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_8_5_FullTrack_Axis_00_27"), track.p(), track.tofNSigmaEl());             // Axis_p             ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_8_6_FullTrack_Axis_01_27"), track.pt(), track.tofNSigmaEl());            // Axis_pt            ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_8_7_FullTrack_Axis_02_27"), track.tpcInnerParam(), track.tofNSigmaEl()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_8_8_FullTrack_Axis_03_27"), track.tofExpMom(), track.tofNSigmaEl());     // Axis_tofExpMom     ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_8_9_FullTrack_Axis_26_27"), track.tpcNSigmaEl(), track.tofNSigmaEl());   // Axis_tpcInnerParam ;
    // Electron
    // Look at Deuteron
    SE_recoTracks.fill(HIST("SE_hTracks_00_9_1_FullTrack_Axis_00_28"), track.p(), track.tpcNSigmaDe());             // Axis_p             ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_9_2_FullTrack_Axis_01_28"), track.pt(), track.tpcNSigmaDe());            // Axis_pt            ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_9_3_FullTrack_Axis_02_28"), track.tpcInnerParam(), track.tpcNSigmaDe()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_9_4_FullTrack_Axis_03_28"), track.tofExpMom(), track.tpcNSigmaDe());     // Axis_tofExpMom     ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_9_5_FullTrack_Axis_00_29"), track.p(), track.tofNSigmaDe());             // Axis_p             ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_9_6_FullTrack_Axis_01_29"), track.pt(), track.tofNSigmaDe());            // Axis_pt            ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_9_7_FullTrack_Axis_02_29"), track.tpcInnerParam(), track.tofNSigmaDe()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_9_8_FullTrack_Axis_03_29"), track.tofExpMom(), track.tofNSigmaDe());     // Axis_tofExpMom     ;
    SE_recoTracks.fill(HIST("SE_hTracks_00_9_9_FullTrack_Axis_28_29"), track.tpcNSigmaDe(), track.tofNSigmaDe());   // Axis_tpcInnerParam ;
    // Deuteron
    //
  }

  template <typename T>
  void FillSelectedTrackQA(T track)
  {
    // Full Track Information
    SE_recoTracks.fill(HIST("SE_htracks_11_1_0_SelectedTrack_P"), track.p());
    SE_recoTracks.fill(HIST("SE_htracks_11_1_0_SelectedTrack_tpcInnerParam"), track.tpcInnerParam());
    SE_recoTracks.fill(HIST("SE_htracks_11_1_0_SelectedTrack_tofExpMom"), track.tofExpMom());

    SE_recoTracks.fill(HIST("SE_htracks_11_1_1_SelectedTrack_Pt"), track.pt());
    SE_recoTracks.fill(HIST("SE_htracks_11_1_2_SelectedTrack_Eta"), track.eta());
    SE_recoTracks.fill(HIST("SE_htracks_11_1_3_SelectedTrack_Phi"), track.phi());
    SE_recoTracks.fill(HIST("SE_htracks_11_1_4_SelectedTrack_DcaXY"), track.dcaXY());
    SE_recoTracks.fill(HIST("SE_htracks_11_1_5_SelectedTrack_DcaZ"), track.dcaZ());
    SE_recoTracks.fill(HIST("SE_htracks_11_1_6_SelectedTrack_Sign"), track.sign());

    // DcaXY
    SE_recoTracks.fill(HIST("SE_htracks_11_1_7_SelectedTrack_00_DcaXY"), track.p(), track.dcaXY());
    SE_recoTracks.fill(HIST("SE_htracks_11_1_7_SelectedTrack_01_DcaXY"), track.pt(), track.dcaXY());
    SE_recoTracks.fill(HIST("SE_htracks_11_1_7_SelectedTrack_02_DcaXY"), track.tpcInnerParam(), track.dcaXY());
    SE_recoTracks.fill(HIST("SE_htracks_11_1_7_SelectedTrack_03_DcaXY"), track.tofExpMom(), track.dcaXY());

    // DcaZ
    SE_recoTracks.fill(HIST("SE_htracks_11_1_7_SelectedTrack_00_DcaZ"), track.p(), track.dcaZ());
    SE_recoTracks.fill(HIST("SE_htracks_11_1_7_SelectedTrack_01_DcaZ"), track.pt(), track.dcaZ());
    SE_recoTracks.fill(HIST("SE_htracks_11_1_7_SelectedTrack_02_DcaZ"), track.tpcInnerParam(), track.dcaZ());
    SE_recoTracks.fill(HIST("SE_htracks_11_1_7_SelectedTrack_03_DcaZ"), track.tofExpMom(), track.dcaZ());

    // momemtum
    SE_recoTracks.fill(HIST("SE_htracks_11_2_1_SelectedTrack_Axis_00_01"), track.p(), track.pt());            // Axis_pt            ;
    SE_recoTracks.fill(HIST("SE_htracks_11_2_2_SelectedTrack_Axis_00_02"), track.p(), track.tpcInnerParam()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_htracks_11_2_3_SelectedTrack_Axis_00_03"), track.p(), track.tofExpMom());     // Axis_tofExpMom     ;

    // tpcSignal
    SE_recoTracks.fill(HIST("SE_htracks_11_3_1_SelectedTrack_Axis_00_05"), track.p(), track.tpcSignal());             // Axis_p             ;
    SE_recoTracks.fill(HIST("SE_htracks_11_3_2_SelectedTrack_Axis_02_05"), track.tpcInnerParam(), track.tpcSignal()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_htracks_11_3_3_SelectedTrack_Axis_03_05"), track.tofExpMom(), track.tpcSignal());     // Axis_tofExpMom     ;

    // tofBeta
    SE_recoTracks.fill(HIST("SE_htracks_11_4_1_SelectedTrack_Axis_00_06"), track.p(), track.beta());             // Axis_p             ;
    SE_recoTracks.fill(HIST("SE_htracks_11_4_2_SelectedTrack_Axis_02_06"), track.tpcInnerParam(), track.beta()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_htracks_11_4_3_SelectedTrack_Axis_03_06"), track.tofExpMom(), track.beta());     // Axis_tofExpMom     ;

    // Look at Pion
    SE_recoTracks.fill(HIST("SE_htracks_11_5_1_SelectedTrack_Axis_00_20"), track.p(), track.tpcNSigmaPi());             // Axis_p             ;
    SE_recoTracks.fill(HIST("SE_htracks_11_5_2_SelectedTrack_Axis_01_20"), track.pt(), track.tpcNSigmaPi());            // Axis_pt            ;
    SE_recoTracks.fill(HIST("SE_htracks_11_5_3_SelectedTrack_Axis_02_20"), track.tpcInnerParam(), track.tpcNSigmaPi()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_htracks_11_5_4_SelectedTrack_Axis_03_20"), track.tofExpMom(), track.tpcNSigmaPi());     // Axis_tofExpMom     ;
    SE_recoTracks.fill(HIST("SE_htracks_11_5_5_SelectedTrack_Axis_00_21"), track.p(), track.tofNSigmaPi());             // Axis_p             ;
    SE_recoTracks.fill(HIST("SE_htracks_11_5_6_SelectedTrack_Axis_01_21"), track.pt(), track.tofNSigmaPi());            // Axis_pt            ;
    SE_recoTracks.fill(HIST("SE_htracks_11_5_7_SelectedTrack_Axis_02_21"), track.tpcInnerParam(), track.tofNSigmaPi()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_htracks_11_5_8_SelectedTrack_Axis_03_21"), track.tofExpMom(), track.tofNSigmaPi());     // Axis_tofExpMom     ;
    SE_recoTracks.fill(HIST("SE_htracks_11_5_9_SelectedTrack_Axis_20_21"), track.tpcNSigmaPi(), track.tofNSigmaPi());   // Axis_tpcInnerParam ;
    // Pion
    // Look at Kaon
    SE_recoTracks.fill(HIST("SE_htracks_11_6_1_SelectedTrack_Axis_00_22"), track.p(), track.tpcNSigmaKa());             // Axis_p             ;
    SE_recoTracks.fill(HIST("SE_htracks_11_6_2_SelectedTrack_Axis_01_22"), track.pt(), track.tpcNSigmaKa());            // Axis_pt            ;
    SE_recoTracks.fill(HIST("SE_htracks_11_6_3_SelectedTrack_Axis_02_22"), track.tpcInnerParam(), track.tpcNSigmaKa()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_htracks_11_6_4_SelectedTrack_Axis_03_22"), track.tofExpMom(), track.tpcNSigmaKa());     // Axis_tofExpMom     ;
    SE_recoTracks.fill(HIST("SE_htracks_11_6_5_SelectedTrack_Axis_00_23"), track.p(), track.tofNSigmaKa());             // Axis_p             ;
    SE_recoTracks.fill(HIST("SE_htracks_11_6_6_SelectedTrack_Axis_01_23"), track.pt(), track.tofNSigmaKa());            // Axis_pt            ;
    SE_recoTracks.fill(HIST("SE_htracks_11_6_7_SelectedTrack_Axis_02_23"), track.tpcInnerParam(), track.tofNSigmaKa()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_htracks_11_6_8_SelectedTrack_Axis_03_23"), track.tofExpMom(), track.tofNSigmaKa());     // Axis_tofExpMom     ;
    SE_recoTracks.fill(HIST("SE_htracks_11_6_9_SelectedTrack_Axis_22_23"), track.tpcNSigmaKa(), track.tofNSigmaKa());   // Axis_tpcInnerParam ;
    // Kaon
    // Look at Proton
    SE_recoTracks.fill(HIST("SE_htracks_11_7_1_SelectedTrack_Axis_00_24"), track.p(), track.tpcNSigmaPr());             // Axis_p             ;
    SE_recoTracks.fill(HIST("SE_htracks_11_7_2_SelectedTrack_Axis_01_24"), track.pt(), track.tpcNSigmaPr());            // Axis_pt            ;
    SE_recoTracks.fill(HIST("SE_htracks_11_7_3_SelectedTrack_Axis_02_24"), track.tpcInnerParam(), track.tpcNSigmaPr()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_htracks_11_7_4_SelectedTrack_Axis_03_24"), track.tofExpMom(), track.tpcNSigmaPr());     // Axis_tofExpMom     ;
    SE_recoTracks.fill(HIST("SE_htracks_11_7_5_SelectedTrack_Axis_00_25"), track.p(), track.tofNSigmaPr());             // Axis_p             ;
    SE_recoTracks.fill(HIST("SE_htracks_11_7_6_SelectedTrack_Axis_01_25"), track.pt(), track.tofNSigmaPr());            // Axis_pt            ;
    SE_recoTracks.fill(HIST("SE_htracks_11_7_7_SelectedTrack_Axis_02_25"), track.tpcInnerParam(), track.tofNSigmaPr()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_htracks_11_7_8_SelectedTrack_Axis_03_25"), track.tofExpMom(), track.tofNSigmaPr());     // Axis_tofExpMom     ;
    SE_recoTracks.fill(HIST("SE_htracks_11_7_9_SelectedTrack_Axis_24_25"), track.tpcNSigmaPr(), track.tofNSigmaPr());   // Axis_tpcInnerParam ;
    // Proton
    // Look at Electron
    SE_recoTracks.fill(HIST("SE_htracks_11_8_1_SelectedTrack_Axis_00_26"), track.p(), track.tpcNSigmaEl());             // Axis_p             ;
    SE_recoTracks.fill(HIST("SE_htracks_11_8_2_SelectedTrack_Axis_01_26"), track.pt(), track.tpcNSigmaEl());            // Axis_pt            ;
    SE_recoTracks.fill(HIST("SE_htracks_11_8_3_SelectedTrack_Axis_02_26"), track.tpcInnerParam(), track.tpcNSigmaEl()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_htracks_11_8_4_SelectedTrack_Axis_03_26"), track.tofExpMom(), track.tpcNSigmaEl());     // Axis_tofExpMom     ;
    SE_recoTracks.fill(HIST("SE_htracks_11_8_5_SelectedTrack_Axis_00_27"), track.p(), track.tofNSigmaEl());             // Axis_p             ;
    SE_recoTracks.fill(HIST("SE_htracks_11_8_6_SelectedTrack_Axis_01_27"), track.pt(), track.tofNSigmaEl());            // Axis_pt            ;
    SE_recoTracks.fill(HIST("SE_htracks_11_8_7_SelectedTrack_Axis_02_27"), track.tpcInnerParam(), track.tofNSigmaEl()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_htracks_11_8_8_SelectedTrack_Axis_03_27"), track.tofExpMom(), track.tofNSigmaEl());     // Axis_tofExpMom     ;
    SE_recoTracks.fill(HIST("SE_htracks_11_8_9_SelectedTrack_Axis_26_27"), track.tpcNSigmaEl(), track.tofNSigmaEl());   // Axis_tpcInnerParam ;
    // Electron
    // Look at Deuteron
    SE_recoTracks.fill(HIST("SE_htracks_11_9_1_SelectedTrack_Axis_00_28"), track.p(), track.tpcNSigmaDe());             // Axis_p             ;
    SE_recoTracks.fill(HIST("SE_htracks_11_9_2_SelectedTrack_Axis_01_28"), track.pt(), track.tpcNSigmaDe());            // Axis_pt            ;
    SE_recoTracks.fill(HIST("SE_htracks_11_9_3_SelectedTrack_Axis_02_28"), track.tpcInnerParam(), track.tpcNSigmaDe()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_htracks_11_9_4_SelectedTrack_Axis_03_28"), track.tofExpMom(), track.tpcNSigmaDe());     // Axis_tofExpMom     ;
    SE_recoTracks.fill(HIST("SE_htracks_11_9_5_SelectedTrack_Axis_00_29"), track.p(), track.tofNSigmaDe());             // Axis_p             ;
    SE_recoTracks.fill(HIST("SE_htracks_11_9_6_SelectedTrack_Axis_01_29"), track.pt(), track.tofNSigmaDe());            // Axis_pt            ;
    SE_recoTracks.fill(HIST("SE_htracks_11_9_7_SelectedTrack_Axis_02_29"), track.tpcInnerParam(), track.tofNSigmaDe()); // Axis_tpcInnerParam ;
    SE_recoTracks.fill(HIST("SE_htracks_11_9_8_SelectedTrack_Axis_03_29"), track.tofExpMom(), track.tofNSigmaDe());     // Axis_tofExpMom     ;
    SE_recoTracks.fill(HIST("SE_htracks_11_9_9_SelectedTrack_Axis_28_29"), track.tpcNSigmaDe(), track.tofNSigmaDe());   // Axis_tpcInnerParam ;
    // Deuteron
    //
  }

  template <typename T>
  void SE_FillTriggerHadronQA(T track)
  {
    SE_recoTrigger.fill(HIST("SE_hTrigger_0_TriggerTrack_P"), track.p());
    SE_recoTrigger.fill(HIST("SE_hTrigger_0_TriggerTrack_tpcInnerParam"), track.tpcInnerParam());
    SE_recoTrigger.fill(HIST("SE_hTrigger_0_TriggerTrack_tofExpMom"), track.tofExpMom());

    SE_recoTrigger.fill(HIST("SE_hTrigger_1_TriggerTrack_Pt"), track.pt());
    SE_recoTrigger.fill(HIST("SE_hTrigger_2_TriggerTrack_Eta"), track.eta());
    SE_recoTrigger.fill(HIST("SE_hTrigger_3_TriggerTrack_Phi"), track.phi());
    SE_recoTrigger.fill(HIST("SE_hTrigger_4_TriggerTrack_DcaXY"), track.dcaXY());
    SE_recoTrigger.fill(HIST("SE_hTrigger_5_TriggerTrack_DcaZ"), track.dcaZ());
    SE_recoTrigger.fill(HIST("SE_hTrigger_6_TriggerTrack_Sign"), track.sign());
    SE_recoTrigger.fill(HIST("SE_hTrigger_11_TriggerTrack_IdentificationTag"), FindTrackTag(track));
  }

  template <typename T>
  void fillKaonQA(T track, int IdMethod)
  {
    if (IdMethod == 0) {
      // tpcSignal
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id0_1_Axis_00_05"), track.p(), track.tpcSignal());             // Axis_p             ;
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id0_1_Axis_02_05"), track.tpcInnerParam(), track.tpcSignal()); // Axis_tpcInnerParam ;
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id0_1_Axis_03_05"), track.tofExpMom(), track.tpcSignal());     // Axis_tofExpMom     ;
      // tofBeta
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id0_2_Axis_00_06"), track.p(), track.beta());             // Axis_p             ;
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id0_2_Axis_02_06"), track.tpcInnerParam(), track.beta()); // Axis_tpcInnerParam ;
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id0_2_Axis_03_06"), track.tofExpMom(), track.beta());     // Axis_tofExpMom     ;

      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id0_3_1_Axis_00_22"), track.p(), track.tpcNSigmaKa());             // Axis_p             ;
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id0_3_2_Axis_01_22"), track.pt(), track.tpcNSigmaKa());            // Axis_pt            ;
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id0_3_3_Axis_02_22"), track.tpcInnerParam(), track.tpcNSigmaKa()); // Axis_tpcInnerParam ;
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id0_3_4_Axis_03_22"), track.tofExpMom(), track.tpcNSigmaKa());     // Axis_tofExpMom     ;
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id0_3_5_Axis_00_23"), track.p(), track.tofNSigmaKa());             // Axis_p             ;
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id0_3_6_Axis_01_23"), track.pt(), track.tofNSigmaKa());            // Axis_pt            ;
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id0_3_7_Axis_02_23"), track.tpcInnerParam(), track.tofNSigmaKa()); // Axis_tpcInnerParam ;
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id0_3_8_Axis_03_23"), track.tofExpMom(), track.tofNSigmaKa());     // Axis_tofExpMom     ;
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id0_3_9_Axis_22_23"), track.tpcNSigmaKa(), track.tofNSigmaKa());   // Axis_tpcInnerParam ;
    } else if (IdMethod == 1) {
      // tpcSignal
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id1_1_Axis_00_05"), track.p(), track.tpcSignal());             // Axis_p             ;
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id1_1_Axis_02_05"), track.tpcInnerParam(), track.tpcSignal()); // Axis_tpcInnerParam ;
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id1_1_Axis_03_05"), track.tofExpMom(), track.tpcSignal());     // Axis_tofExpMom     ;
      // tofBeta
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id1_2_Axis_00_06"), track.p(), track.beta());             // Axis_p             ;
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id1_2_Axis_02_06"), track.tpcInnerParam(), track.beta()); // Axis_tpcInnerParam ;
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id1_2_Axis_03_06"), track.tofExpMom(), track.beta());     // Axis_tofExpMom     ;

      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id1_3_1_Axis_00_22"), track.p(), track.tpcNSigmaKa());             // Axis_p             ;
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id1_3_2_Axis_01_22"), track.pt(), track.tpcNSigmaKa());            // Axis_pt            ;
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id1_3_3_Axis_02_22"), track.tpcInnerParam(), track.tpcNSigmaKa()); // Axis_tpcInnerParam ;
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id1_3_4_Axis_03_22"), track.tofExpMom(), track.tpcNSigmaKa());     // Axis_tofExpMom     ;
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id1_3_5_Axis_00_23"), track.p(), track.tofNSigmaKa());             // Axis_p             ;
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id1_3_6_Axis_01_23"), track.pt(), track.tofNSigmaKa());            // Axis_pt            ;
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id1_3_7_Axis_02_23"), track.tpcInnerParam(), track.tofNSigmaKa()); // Axis_tpcInnerParam ;
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id1_3_8_Axis_03_23"), track.tofExpMom(), track.tofNSigmaKa());     // Axis_tofExpMom     ;
      SE_recoKaon.fill(HIST("SE_hKaon_12_Ka_Id1_3_9_Axis_22_23"), track.tpcNSigmaKa(), track.tofNSigmaKa());   // Axis_tpcInnerParam ;
    }
  }

  template <typename T>
  bool selPhiMeson(T track1, T track2, double& PhiPt)
  {
    daughter1.SetXYZM(track1.px(), track1.py(), track1.pz(), pdgDB->Mass(321)); // set the daughter1 4-momentum
    daughter2.SetXYZM(track2.px(), track2.py(), track2.pz(), pdgDB->Mass(321)); // set the daughter2
    mother = daughter1 + daughter2;                                             // calculate the mother 4-momentum;
    mass = mother.M();
    PhiPt = mother.Pt();
    if (cfgPhiMassLow <= mass && mass <= cfgPhiMassUp) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool selLSBMeson(T track1, T track2, double& LSBPt)
  {
    daughter1.SetXYZM(track1.px(), track1.py(), track1.pz(), pdgDB->Mass(321)); // set the daughter1 4-momentum
    daughter2.SetXYZM(track2.px(), track2.py(), track2.pz(), pdgDB->Mass(321)); // set the daughter2
    mother = daughter1 + daughter2;                                             // calculate the mother 4-momentum;
    mass = mother.M();
    LSBPt = mother.Pt();
    if (cfgLSBMassLow <= mass && mass <= cfgLSBMassUp) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool selRSBMeson(T track1, T track2, double& RSBPt)
  {
    daughter1.SetXYZM(track1.px(), track1.py(), track1.pz(), pdgDB->Mass(321)); // set the daughter1 4-momentum
    daughter2.SetXYZM(track2.px(), track2.py(), track2.pz(), pdgDB->Mass(321)); // set the daughter2
    mother = daughter1 + daughter2;                                             // calculate the mother 4-momentum;
    mass = mother.M();
    RSBPt = mother.Pt();
    if (cfgRSBMassLow <= mass && mass <= cfgRSBMassUp) {
      return true;
    }
    return false;
  }

  template <typename T>
  void FillPhiQA(T track1, T track2)
  {

    daughter1.SetXYZM(track1.px(), track1.py(), track1.pz(), pdgDB->Mass(321)); // set the daughter1 4-momentum
    daughter2.SetXYZM(track2.px(), track2.py(), track2.pz(), pdgDB->Mass(321)); // set the daughter2
    mother = daughter1 + daughter2;                                             // calculate the mother 4-momentum;
    paircharge = track1.sign() * track2.sign();

    SE_recoPhi.fill(HIST("SE_hPhi_KK_allPairs_01_Mass"), mother.M());
    SE_recoPhi.fill(HIST("SE_hPhi_KK_allPairs_02_P"), mother.P());
    SE_recoPhi.fill(HIST("SE_hPhi_KK_allPairs_03_Pt"), mother.Pt());
    SE_recoPhi.fill(HIST("SE_hPhi_KK_allPairs_04_Eta"), mother.Eta());
    SE_recoPhi.fill(HIST("SE_hPhi_KK_allPairs_05_Phi"), mother.Phi());
    SE_recoPhi.fill(HIST("SE_hPhi_KK_allPairs_06_Rapidity"), mother.Rapidity());

    if (paircharge < 0) {
      SE_recoPhi.fill(HIST("SE_hPhi_KK_USpairs_01_Mass"), mother.M());
      SE_recoPhi.fill(HIST("SE_hPhi_KK_USpairs_02_P"), mother.P());
      SE_recoPhi.fill(HIST("SE_hPhi_KK_USpairs_03_Pt"), mother.Pt());
      SE_recoPhi.fill(HIST("SE_hPhi_KK_USpairs_04_Eta"), mother.Eta());
      SE_recoPhi.fill(HIST("SE_hPhi_KK_USpairs_05_Phi"), mother.Phi());
      SE_recoPhi.fill(HIST("SE_hPhi_KK_USpairs_06_Rapidity"), mother.Rapidity());
    }
    if (paircharge > 0) {
      SE_recoPhi.fill(HIST("SE_hPhi_KK_LSpairs_01_Mass"), mother.M());
      SE_recoPhi.fill(HIST("SE_hPhi_KK_LSpairs_02_P"), mother.P());
      SE_recoPhi.fill(HIST("SE_hPhi_KK_LSpairs_03_Pt"), mother.Pt());
      SE_recoPhi.fill(HIST("SE_hPhi_KK_LSpairs_04_Eta"), mother.Eta());
      SE_recoPhi.fill(HIST("SE_hPhi_KK_LSpairs_05_Phi"), mother.Phi());
      SE_recoPhi.fill(HIST("SE_hPhi_KK_LSpairs_06_Rapidity"), mother.Rapidity());
    }
  }

  template <typename T> // Check and Remove it// for debugging purpose
  void FillLS_Analysis(const T& track1, const T& track2)
  {
    daughter1.SetXYZM(track1.px(), track1.py(), track1.pz(), pdgDB->Mass(321)); // set the daughter1 4-momentum
    daughter2.SetXYZM(track2.px(), track2.py(), track2.pz(), pdgDB->Mass(321)); // set the daughter2
    mother = daughter1 + daughter2;                                             // calculate the mother 4-momentum;
    mass = mother.M();

    // //Unlike Sign Peak Region
    // if ( cfgPhiMassLow <= mass &&  mass <= cfgPhiMassUp ){
    //   //Unlike Sign Peak Region
    // phi_pT in [0.0-inf]
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_AllPhi_Phi_p"), mother.P());
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_AllPhi_Phi_pT"), mother.Pt());
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_AllPhi_Phi_eta"), mother.Eta());
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_AllPhi_Phi_phi"), mother.Phi());
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_AllPhi_Phi_rapidity"), mother.Rapidity());
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_AllPhi_Phi_paircharge"), paircharge);
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_AllPhi_Phi_InvMass"), mother.M());

    if (pT < 2.0) { // Phi_pT in [0.0-2.0](Bulk)
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_0To2Phi_Phi_p"), mother.P());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_0To2Phi_Phi_pT"), mother.Pt());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_0To2Phi_Phi_eta"), mother.Eta());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_0To2Phi_Phi_phi"), mother.Phi());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_0To2Phi_Phi_rapidity"), mother.Rapidity());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_0To2Phi_Phi_paircharge"), paircharge);
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_0To2Phi_Phi_InvMass"), mother.M());
    }
    if (2.0 <= pT && pT < 4.0) { // Phi_pT in [2.0-4.0](Required)
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_2To4Phi_Phi_p"), mother.P());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_2To4Phi_Phi_pT"), mother.Pt());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_2To4Phi_Phi_eta"), mother.Eta());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_2To4Phi_Phi_phi"), mother.Phi());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_2To4Phi_Phi_rapidity"), mother.Rapidity());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_2To4Phi_Phi_paircharge"), paircharge);
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_2To4Phi_Phi_InvMass"), mother.M());
    }
    if (4.0 <= pT && pT < 8.0) { // Phi_pT in [4.0-8.0](Hard)
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_4To8Phi_Phi_p"), mother.P());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_4To8Phi_Phi_pT"), mother.Pt());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_4To8Phi_Phi_eta"), mother.Eta());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_4To8Phi_Phi_phi"), mother.Phi());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_4To8Phi_Phi_rapidity"), mother.Rapidity());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_4To8Phi_Phi_paircharge"), paircharge);
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_4To8Phi_Phi_InvMass"), mother.M());
    }
    if (8.0 <= pT) { // Phi_pT in [8.0-inf](Very Hard)
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_8ToInfPhi_Phi_p"), mother.P());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_8ToInfPhi_Phi_pT"), mother.Pt());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_8ToInfPhi_Phi_eta"), mother.Eta());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_8ToInfPhi_Phi_phi"), mother.Phi());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_8ToInfPhi_Phi_rapidity"), mother.Rapidity());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_8ToInfPhi_Phi_paircharge"), paircharge);
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_LS_KK_8ToInfPhi_Phi_InvMass"), mother.M());
    }
    //
    // }//Unlike Sign Peak Region
  }

  double rapidity, mass, pT, p, phi, eta, paircharge;
  TLorentzVector daughter1, daughter2, mother;
  template <typename U, typename T>
  void SE_FillRHCorrelationUS(const U& collision, const T& triggerTrack, const T& track1, const T& track2, int& nCR_Phi, int& nCR_Phi_0_2, int& nCR_Phi_2_4, int& nCR_Phi_4_8, int& nCR_Phi_8_i)
  {
    daughter1.SetXYZM(track1.px(), track1.py(), track1.pz(), pdgDB->Mass(321)); // set the daughter1 4-momentum
    daughter2.SetXYZM(track2.px(), track2.py(), track2.pz(), pdgDB->Mass(321)); // set the daughter2
    mother = daughter1 + daughter2;                                             // calculate the mother 4-momentum;
    mass = mother.M();
    pT = mother.Pt();
    // rapidity = mother.Rapidity();
    // p = mother.P();
    // eta = mother.Eta();
    // phi = mother.Phi();
    paircharge = track1.sign() * track2.sign();

    // Unlike Sign Peak Region
    if (cfgPhiMassLow <= mass && mass <= cfgPhiMassUp) {
      // Unlike Sign Peak Region
      // phi_pT in [0.0-inf]
      nCR_Phi++;
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_AllPhi_Phi_p"), mother.P());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_AllPhi_Phi_pT"), mother.Pt());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_AllPhi_Phi_eta"), mother.Eta());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_AllPhi_Phi_phi"), mother.Phi());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_AllPhi_Phi_rapidity"), mother.Rapidity());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_AllPhi_Phi_paircharge"), paircharge);
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_AllPhi_Phi_InvMass"), mother.M());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_AllPhi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_AllPhi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
      if (pT < 2.0) { // Phi_pT in [0.0-2.0](Bulk)
        nCR_Phi_0_2++;
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_p"), mother.P());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_pT"), mother.Pt());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_eta"), mother.Eta());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_phi"), mother.Phi());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_rapidity"), mother.Rapidity());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_paircharge"), paircharge);
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_InvMass"), mother.M());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_0To2Phi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_0To2Phi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
      }
      if (2.0 <= pT && pT < 4.0) { // Phi_pT in [2.0-4.0](Required)
        nCR_Phi_2_4++;
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_p"), mother.P());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_pT"), mother.Pt());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_eta"), mother.Eta());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_phi"), mother.Phi());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_rapidity"), mother.Rapidity());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_paircharge"), paircharge);
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_InvMass"), mother.M());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_2To4Phi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_2To4Phi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
      }
      if (4.0 <= pT && pT < 8.0) { // Phi_pT in [4.0-8.0](Hard)
        nCR_Phi_4_8++;
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_p"), mother.P());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_pT"), mother.Pt());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_eta"), mother.Eta());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_phi"), mother.Phi());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_rapidity"), mother.Rapidity());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_paircharge"), paircharge);
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_InvMass"), mother.M());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_4To8Phi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_4To8Phi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
      }
      if (8.0 <= pT) { // Phi_pT in [8.0-inf](Very Hard)
        nCR_Phi_8_i++;
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_p"), mother.P());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_pT"), mother.Pt());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_eta"), mother.Eta());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_phi"), mother.Phi());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_rapidity"), mother.Rapidity());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_paircharge"), paircharge);
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_InvMass"), mother.M());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_Peak_8ToInfPhi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
      }
      //
    } // Unlike Sign Peak Region

    // Unlike Sign LSB Region
    if (cfgLSBMassLow <= mass && mass <= cfgLSBMassUp) {
      // Unlike Sign LSB Region
      // phi_pT in [0.0-inf]
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_AllPhi_Phi_p"), mother.P());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_AllPhi_Phi_pT"), mother.Pt());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_AllPhi_Phi_eta"), mother.Eta());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_AllPhi_Phi_phi"), mother.Phi());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_AllPhi_Phi_rapidity"), mother.Rapidity());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_AllPhi_Phi_paircharge"), paircharge);
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_AllPhi_Phi_InvMass"), mother.M());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_AllPhi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_AllPhi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
      if (pT < 2.0) { // Phi_pT in [0.0-2.0](Bulk)
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_p"), mother.P());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_pT"), mother.Pt());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_eta"), mother.Eta());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_phi"), mother.Phi());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_rapidity"), mother.Rapidity());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_paircharge"), paircharge);
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_InvMass"), mother.M());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_0To2Phi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_0To2Phi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
      }
      if (2.0 <= pT && pT < 4.0) { // Phi_pT in [2.0-4.0](Required)
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_p"), mother.P());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_pT"), mother.Pt());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_eta"), mother.Eta());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_phi"), mother.Phi());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_rapidity"), mother.Rapidity());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_paircharge"), paircharge);
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_InvMass"), mother.M());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_2To4Phi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_2To4Phi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
      }
      if (4.0 <= pT && pT < 8.0) { // Phi_pT in [4.0-8.0](Hard)
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_p"), mother.P());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_pT"), mother.Pt());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_eta"), mother.Eta());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_phi"), mother.Phi());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_rapidity"), mother.Rapidity());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_paircharge"), paircharge);
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_4To8Phi_Phi_InvMass"), mother.M());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_4To8Phi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_4To8Phi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
      }
      if (8.0 <= pT) { // Phi_pT in [8.0-inf](Very Hard)
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_p"), mother.P());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_pT"), mother.Pt());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_eta"), mother.Eta());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_phi"), mother.Phi());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_rapidity"), mother.Rapidity());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_paircharge"), paircharge);
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_Phi_InvMass"), mother.M());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_LSB_8ToInfPhi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
      }
      //
    } // Unlike Sign LSB Region

    // Unlike Sign RSB Region
    if (cfgRSBMassLow <= mass && mass <= cfgRSBMassUp) {
      // Unlike Sign RSB Region
      // phi_pT in [0.0-inf]
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_AllPhi_Phi_p"), mother.P());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_AllPhi_Phi_pT"), mother.Pt());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_AllPhi_Phi_eta"), mother.Eta());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_AllPhi_Phi_phi"), mother.Phi());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_AllPhi_Phi_rapidity"), mother.Rapidity());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_AllPhi_Phi_paircharge"), paircharge);
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_AllPhi_Phi_InvMass"), mother.M());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_AllPhi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
      SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_AllPhi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
      if (pT < 2.0) { // Phi_pT in [0.0-2.0](Bulk)
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_p"), mother.P());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_pT"), mother.Pt());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_eta"), mother.Eta());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_phi"), mother.Phi());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_rapidity"), mother.Rapidity());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_paircharge"), paircharge);
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_InvMass"), mother.M());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_0To2Phi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_0To2Phi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
      }
      if (2.0 <= pT && pT < 4.0) { // Phi_pT in [2.0-4.0](Required)
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_p"), mother.P());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_pT"), mother.Pt());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_eta"), mother.Eta());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_phi"), mother.Phi());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_rapidity"), mother.Rapidity());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_paircharge"), paircharge);
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_InvMass"), mother.M());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_2To4Phi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_2To4Phi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
      }
      if (4.0 <= pT && pT < 8.0) { // Phi_pT in [4.0-8.0](Hard)
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_p"), mother.P());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_pT"), mother.Pt());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_eta"), mother.Eta());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_phi"), mother.Phi());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_rapidity"), mother.Rapidity());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_paircharge"), paircharge);
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_4To8Phi_Phi_InvMass"), mother.M());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_4To8Phi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_4To8Phi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
      }
      if (8.0 <= pT) { // Phi_pT in [8.0-inf](Very Hard)
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_p"), mother.P());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_pT"), mother.Pt());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_eta"), mother.Eta());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_phi"), mother.Phi());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_rapidity"), mother.Rapidity());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_paircharge"), paircharge);
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_Phi_InvMass"), mother.M());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
        SE_recoAnalysis.fill(HIST("SE_hAnalysis_hPhi_US_RSB_8ToInfPhi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
      }
      //
    } // Unlike Sign RSB Region
  }

  TLorentzVector phi1mother, daughter3, daughter4, phi2mother;
  TLorentzVector leadMother, assocMother;
  template <typename U, typename T>
  void SE_FillPhiPhiCorrelation(U collision, T track1, T track2, T track3, T track4, int& nLeadPhi, int& nAssoPhi_0_2, int& nAssoPhi_2_4)
  {

    daughter1.SetXYZM(track1.px(), track1.py(), track1.pz(), pdgDB->Mass(321)); // set the daughter1 4-momentum
    daughter2.SetXYZM(track2.px(), track2.py(), track2.pz(), pdgDB->Mass(321)); // set the daughter2
    phi1mother = daughter1 + daughter2;                                         // calculate the mother 4-momentum;

    daughter1.SetXYZM(track3.px(), track3.py(), track3.pz(), pdgDB->Mass(321)); // set the daughter1 4-momentum
    daughter2.SetXYZM(track4.px(), track4.py(), track4.pz(), pdgDB->Mass(321)); // set the daughter2
    phi2mother = daughter1 + daughter2;                                         // calculate the mother 4-momentum;

    if (phi1mother.Pt() > phi2mother.Pt()) {
      leadMother = phi1mother;
      assocMother = phi2mother;
    } else {
      leadMother = phi2mother;
      assocMother = phi1mother;
    }

    double dphi = ComputeDeltaPhi(assocMother.Phi(), leadMother.Phi());
    double deta = assocMother.Eta() - leadMother.Eta();

    SE_recoPhi.fill(HIST("SE_hPhiPhi_Full_LeadPhi_01_Mass"), leadMother.M());
    SE_recoPhi.fill(HIST("SE_hPhiPhi_Full_LeadPhi_02_P"), leadMother.P());
    SE_recoPhi.fill(HIST("SE_hPhiPhi_Full_LeadPhi_03_Pt"), leadMother.Pt());
    SE_recoPhi.fill(HIST("SE_hPhiPhi_Full_LeadPhi_04_Eta"), leadMother.Eta());
    SE_recoPhi.fill(HIST("SE_hPhiPhi_Full_LeadPhi_05_Phi"), leadMother.Phi());
    SE_recoPhi.fill(HIST("SE_hPhiPhi_Full_LeadPhi_06_Rapidity"), leadMother.Rapidity());
    SE_recoPhi.fill(HIST("SE_hPhiPhi_Full_AssoPhi_01_Mass"), assocMother.M());
    SE_recoPhi.fill(HIST("SE_hPhiPhi_Full_AssoPhi_02_P"), assocMother.P());
    SE_recoPhi.fill(HIST("SE_hPhiPhi_Full_AssoPhi_03_Pt"), assocMother.Pt());
    SE_recoPhi.fill(HIST("SE_hPhiPhi_Full_AssoPhi_04_Eta"), assocMother.Eta());
    SE_recoPhi.fill(HIST("SE_hPhiPhi_Full_AssoPhi_05_Phi"), assocMother.Phi());
    SE_recoPhi.fill(HIST("SE_hPhiPhi_Full_AssoPhi_06_Rapidity"), assocMother.Rapidity());
    SE_recoAnalysis.fill(HIST("SE_hPhiPhi_Full_dPhi_dEta"), dphi, deta);
    SE_recoAnalysis.fill(HIST("SE_hPhiPhi_Full_Analysis"), collision.centFT0C(), collision.posZ(), dphi, deta);
    if (4.0 < leadMother.Pt() && leadMother.Pt() < 8.0) {
      nLeadPhi++;
      SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_Full_LeadPhi_01_Mass"), leadMother.M());
      SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_Full_LeadPhi_02_P"), leadMother.P());
      SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_Full_LeadPhi_03_Pt"), leadMother.Pt());
      SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_Full_LeadPhi_04_Eta"), leadMother.Eta());
      SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_Full_LeadPhi_05_Phi"), leadMother.Phi());
      SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_Full_LeadPhi_06_Rapidity"), leadMother.Rapidity());
      SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_Full_AssoPhi_01_Mass"), assocMother.M());
      SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_Full_AssoPhi_02_P"), assocMother.P());
      SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_Full_AssoPhi_03_Pt"), assocMother.Pt());
      SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_Full_AssoPhi_04_Eta"), assocMother.Eta());
      SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_Full_AssoPhi_05_Phi"), assocMother.Phi());
      SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_Full_AssoPhi_06_Rapidity"), assocMother.Rapidity());
      SE_recoAnalysis.fill(HIST("SE_hPhiPhi_4to8LeadPhi_Full_dPhi_dEta"), dphi, deta);
      SE_recoAnalysis.fill(HIST("SE_hPhiPhi_4to8LeadPhi_Full_Analysis"), collision.centFT0C(), collision.posZ(), dphi, deta);
      if (assocMother.Pt() < 2.0) {
        nAssoPhi_0_2++;
        SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_01_Mass"), leadMother.M());
        SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_02_P"), leadMother.P());
        SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_03_Pt"), leadMother.Pt());
        SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_04_Eta"), leadMother.Eta());
        SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_05_Phi"), leadMother.Phi());
        SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_06_Rapidity"), leadMother.Rapidity());
        SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_01_Mass"), assocMother.M());
        SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_02_P"), assocMother.P());
        SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_03_Pt"), assocMother.Pt());
        SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_04_Eta"), assocMother.Eta());
        SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_05_Phi"), assocMother.Phi());
        SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_06_Rapidity"), assocMother.Rapidity());
        SE_recoAnalysis.fill(HIST("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_dPhi_dEta"), dphi, deta);
        SE_recoAnalysis.fill(HIST("SE_hPhiPhi_4to8LeadPhi_0to2AssoPhi_Analysis"), collision.centFT0C(), collision.posZ(), dphi, deta);
      } else if (2.0 < assocMother.Pt() && assocMother.Pt() < 4.0) {
        nAssoPhi_2_4++;
        SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_01_Mass"), leadMother.M());
        SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_02_P"), leadMother.P());
        SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_03_Pt"), leadMother.Pt());
        SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_04_Eta"), leadMother.Eta());
        SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_05_Phi"), leadMother.Phi());
        SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_06_Rapidity"), leadMother.Rapidity());
        SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_01_Mass"), assocMother.M());
        SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_02_P"), assocMother.P());
        SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_03_Pt"), assocMother.Pt());
        SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_04_Eta"), assocMother.Eta());
        SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_05_Phi"), assocMother.Phi());
        SE_recoPhi.fill(HIST("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_06_Rapidity"), assocMother.Rapidity());
        SE_recoAnalysis.fill(HIST("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_dPhi_dEta"), dphi, deta);
        SE_recoAnalysis.fill(HIST("SE_hPhiPhi_4to8LeadPhi_2to4AssoPhi_Analysis"), collision.centFT0C(), collision.posZ(), dphi, deta);
      }
    }
  }

  template <typename U, typename T>
  void SE_Fill_hh_0To2AssoHadron(U collision, T triggerTrack, T assocTrack)
  {
    // hadron_pT in [0.0-2.0](Bulk)
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hh_0To2AssoHadron_LeadH_02_P"), triggerTrack.p());
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hh_0To2AssoHadron_LeadH_03_Pt"), triggerTrack.pt());
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hh_0To2AssoHadron_LeadH_04_Eta"), triggerTrack.eta());
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hh_0To2AssoHadron_LeadH_05_Phi"), triggerTrack.phi());
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hh_0To2AssoHadron_AssoH_02_P"), assocTrack.p());
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hh_0To2AssoHadron_AssoH_03_Pt"), assocTrack.pt());
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hh_0To2AssoHadron_AssoH_04_Eta"), assocTrack.eta());
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hh_0To2AssoHadron_AssoH_05_Phi"), assocTrack.phi());
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hh_0To2AssoHadron_dPhi_dEta"), ComputeDeltaPhi(assocTrack.phi(), triggerTrack.phi()), assocTrack.eta() - triggerTrack.eta());
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hh_0To2AssoHadron_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(assocTrack.phi(), triggerTrack.phi()), assocTrack.eta() - triggerTrack.eta());
  }

  template <typename U, typename T>
  void SE_Fill_hh_2To4AssoHadron(U collision, T triggerTrack, T assocTrack)
  {
    // hadron_pT in [0.0-2.0](Bulk)
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hh_2To4AssoHadron_LeadH_02_P"), triggerTrack.p());
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hh_2To4AssoHadron_LeadH_03_Pt"), triggerTrack.pt());
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hh_2To4AssoHadron_LeadH_04_Eta"), triggerTrack.eta());
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hh_2To4AssoHadron_LeadH_05_Phi"), triggerTrack.phi());
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hh_2To4AssoHadron_AssoH_02_P"), assocTrack.p());
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hh_2To4AssoHadron_AssoH_03_Pt"), assocTrack.pt());
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hh_2To4AssoHadron_AssoH_04_Eta"), assocTrack.eta());
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hh_2To4AssoHadron_AssoH_05_Phi"), assocTrack.phi());
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hh_2To4AssoHadron_dPhi_dEta"), ComputeDeltaPhi(assocTrack.phi(), triggerTrack.phi()), assocTrack.eta() - triggerTrack.eta());
    SE_recoAnalysis.fill(HIST("SE_hAnalysis_hh_2To4AssoHadron_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(assocTrack.phi(), triggerTrack.phi()), assocTrack.eta() - triggerTrack.eta());
  }

  template <typename T>
  void ME_FillTriggerHadronQA(T track)
  {
    ME_recoTrigger.fill(HIST("ME_hTrigger_0_TriggerTrack_P"), track.p());
    ME_recoTrigger.fill(HIST("ME_hTrigger_0_TriggerTrack_tpcInnerParam"), track.tpcInnerParam());
    ME_recoTrigger.fill(HIST("ME_hTrigger_0_TriggerTrack_tofExpMom"), track.tofExpMom());

    ME_recoTrigger.fill(HIST("ME_hTrigger_1_TriggerTrack_Pt"), track.pt());
    ME_recoTrigger.fill(HIST("ME_hTrigger_2_TriggerTrack_Eta"), track.eta());
    ME_recoTrigger.fill(HIST("ME_hTrigger_3_TriggerTrack_Phi"), track.phi());
    ME_recoTrigger.fill(HIST("ME_hTrigger_4_TriggerTrack_DcaXY"), track.dcaXY());
    ME_recoTrigger.fill(HIST("ME_hTrigger_5_TriggerTrack_DcaZ"), track.dcaZ());
    ME_recoTrigger.fill(HIST("ME_hTrigger_6_TriggerTrack_Sign"), track.sign());
    ME_recoTrigger.fill(HIST("ME_hTrigger_11_TriggerTrack_IdentificationTag"), FindTrackTag(track));
  }

  // double rapidity, mass, pT, p, phi, eta, paircharge;
  // TLorentzVector daughter1, daughter2, mother;
  template <typename U, typename T>
  void ME_FillRHCorrelationUS(const int& caseNo, const U& collision, const T& triggerTrack, const T& track1, const T& track2)
  {
    daughter1.SetXYZM(track1.px(), track1.py(), track1.pz(), pdgDB->Mass(321)); // set the daughter1 4-momentum
    daughter2.SetXYZM(track2.px(), track2.py(), track2.pz(), pdgDB->Mass(321)); // set the daughter2
    mother = daughter1 + daughter2;                                             // calculate the mother 4-momentum;
    mass = mother.M();
    pT = mother.Pt();
    paircharge = track1.sign() * track2.sign();

    // Unlike Sign Peak Region
    if (caseNo == 0 || caseNo == 1 || caseNo == 2 || caseNo == 3 || caseNo == 4) {
      if (cfgPhiMassLow <= mass && mass <= cfgPhiMassUp) {
        // Unlike Sign Peak Region
        // phi_pT in [0.0-inf]
        if (caseNo == 0) { // Fill 0-i
          ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_AllPhi_Phi_p"), mother.P());
          ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_AllPhi_Phi_pT"), mother.Pt());
          ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_AllPhi_Phi_eta"), mother.Eta());
          ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_AllPhi_Phi_phi"), mother.Phi());
          ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_AllPhi_Phi_rapidity"), mother.Rapidity());
          ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_AllPhi_Phi_paircharge"), paircharge);
          ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_AllPhi_Phi_InvMass"), mother.M());
          ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_AllPhi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
          ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_AllPhi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
        } else if (caseNo == 1) { //
          if (pT < 2.0) {         // Phi_pT in [0.0-2.0](Bulk)
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_p"), mother.P());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_pT"), mother.Pt());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_eta"), mother.Eta());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_phi"), mother.Phi());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_rapidity"), mother.Rapidity());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_paircharge"), paircharge);
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_0To2Phi_Phi_InvMass"), mother.M());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_0To2Phi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_0To2Phi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
          }
        } else if (caseNo == 2) {      //
          if (2.0 <= pT && pT < 4.0) { // Phi_pT in [2.0-4.0](Required)
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_p"), mother.P());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_pT"), mother.Pt());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_eta"), mother.Eta());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_phi"), mother.Phi());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_rapidity"), mother.Rapidity());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_paircharge"), paircharge);
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_2To4Phi_Phi_InvMass"), mother.M());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_2To4Phi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_2To4Phi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
          }
        } else if (caseNo == 3) {      //
          if (4.0 <= pT && pT < 8.0) { // Phi_pT in [4.0-8.0](Hard)
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_p"), mother.P());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_pT"), mother.Pt());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_eta"), mother.Eta());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_phi"), mother.Phi());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_rapidity"), mother.Rapidity());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_paircharge"), paircharge);
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_4To8Phi_Phi_InvMass"), mother.M());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_4To8Phi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_4To8Phi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
          }
        } else if (caseNo == 4) { //
          if (8.0 <= pT) {        // Phi_pT in [8.0-inf](Very Hard)
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_p"), mother.P());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_pT"), mother.Pt());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_eta"), mother.Eta());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_phi"), mother.Phi());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_rapidity"), mother.Rapidity());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_paircharge"), paircharge);
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_Phi_InvMass"), mother.M());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_Peak_8ToInfPhi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
          }
        }
        //
      } // Unlike Sign Peak Region
    }
    // Unlike Sign LSB Region
    if (caseNo == 11 || caseNo == 12) {
      if (cfgLSBMassLow <= mass && mass <= cfgLSBMassUp) {
        // Unlike Sign LSB Region
        // phi_pT in [0.0-inf]
        ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_AllPhi_Phi_p"), mother.P());
        ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_AllPhi_Phi_pT"), mother.Pt());
        ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_AllPhi_Phi_eta"), mother.Eta());
        ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_AllPhi_Phi_phi"), mother.Phi());
        ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_AllPhi_Phi_rapidity"), mother.Rapidity());
        ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_AllPhi_Phi_paircharge"), paircharge);
        ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_AllPhi_Phi_InvMass"), mother.M());
        ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_AllPhi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
        ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_AllPhi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
        if (caseNo == 11) {
          if (pT < 2.0) { // Phi_pT in [0.0-2.0](Bulk)
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_p"), mother.P());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_pT"), mother.Pt());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_eta"), mother.Eta());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_phi"), mother.Phi());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_rapidity"), mother.Rapidity());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_paircharge"), paircharge);
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_0To2Phi_Phi_InvMass"), mother.M());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_0To2Phi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_0To2Phi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
          }
        }
        if (caseNo == 12) {
          if (2.0 <= pT && pT < 4.0) { // Phi_pT in [2.0-4.0](Required)
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_p"), mother.P());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_pT"), mother.Pt());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_eta"), mother.Eta());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_phi"), mother.Phi());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_rapidity"), mother.Rapidity());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_paircharge"), paircharge);
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_2To4Phi_Phi_InvMass"), mother.M());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_2To4Phi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_LSB_2To4Phi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
          }
        }
        //
      } // Unlike Sign LSB Region
    }
    // Unlike Sign RSB Region
    if (caseNo == 13 || caseNo == 14) {
      if (cfgRSBMassLow <= mass && mass <= cfgRSBMassUp) {
        // Unlike Sign RSB Region
        // phi_pT in [0.0-inf]
        ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_AllPhi_Phi_p"), mother.P());
        ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_AllPhi_Phi_pT"), mother.Pt());
        ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_AllPhi_Phi_eta"), mother.Eta());
        ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_AllPhi_Phi_phi"), mother.Phi());
        ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_AllPhi_Phi_rapidity"), mother.Rapidity());
        ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_AllPhi_Phi_paircharge"), paircharge);
        ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_AllPhi_Phi_InvMass"), mother.M());
        ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_AllPhi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
        ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_AllPhi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
        if (caseNo == 13) {
          if (pT < 2.0) { // Phi_pT in [0.0-2.0](Bulk)
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_p"), mother.P());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_pT"), mother.Pt());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_eta"), mother.Eta());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_phi"), mother.Phi());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_rapidity"), mother.Rapidity());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_paircharge"), paircharge);
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_0To2Phi_Phi_InvMass"), mother.M());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_0To2Phi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_0To2Phi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
          }
        }
        if (caseNo == 14) {
          if (2.0 <= pT && pT < 4.0) { // Phi_pT in [2.0-4.0](Required)
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_p"), mother.P());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_pT"), mother.Pt());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_eta"), mother.Eta());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_phi"), mother.Phi());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_rapidity"), mother.Rapidity());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_paircharge"), paircharge);
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_2To4Phi_Phi_InvMass"), mother.M());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_2To4Phi_dPhi_dEta"), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
            ME_recoAnalysis.fill(HIST("ME_hAnalysis_hPhi_US_RSB_2To4Phi_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(mother.Phi(), triggerTrack.phi()), mother.Eta() - triggerTrack.eta());
          }
        }
        //
      } // Unlike Sign RSB Region
    }
  }

  void MIX_Collisions_h_Phi(int caseNo, int iMixBin, int iMixSlot1, int iMixSlot2)
  {
    LOG(info) << "DEBUG :: Mixing hPhi   :: caseNo = " << caseNo << " :: iMixBin = " << iMixBin << " :: iMixSlot_ij = (" << iMixSlot1 << "," << iMixSlot2 << ")";

    auto const& c1 = myStoredMixingCollisions[caseNo][iMixBin][iMixSlot1];

    auto const& PerColl_triggerTrack_c1 = PerColl_triggerTracks[caseNo][iMixBin][iMixSlot1];

    auto const& PerColl_triggerTrack_c2 = PerColl_triggerTracks[caseNo][iMixBin][iMixSlot2];
    auto const& PerColl_posTracks_c2 = PerColl_posTracks[caseNo][iMixBin][iMixSlot2];
    auto const& PerColl_negTracks_c2 = PerColl_negTracks[caseNo][iMixBin][iMixSlot2];

    std::vector<int64_t> triggerTrackIndexList_c2;
    for (auto& triggerTrack : PerColl_triggerTrack_c2) {
      if (!selectionTrack(triggerTrack)) {
        continue;
      }
      triggerTrackIndexList_c2.push_back(triggerTrack.globalIndex());
    }

    // 01-h-phi All Three Region Unlike Sign correlation
    for (auto& triggerTrack : PerColl_triggerTrack_c1) {
      if (!selectionTrack(triggerTrack)) {
        continue;
      }
      // nTrigger++;
      ME_FillTriggerHadronQA(triggerTrack);
      // 01-Start-obtaining h-Phi correlation
      for (auto posTrack : PerColl_posTracks_c2) {
        if (!selectionTrack(posTrack)) {
          continue;
        }
        int KaIdMethod = -1;
        if (!selKaon(posTrack, KaIdMethod)) {
          continue;
        } // Kaon PID Check
        if (checkTrackInList(posTrack, triggerTrackIndexList_c2)) {
          continue;
        } // make sure that posTrack is not a Trigger
        // nPosTrack++;
        for (auto negTrack : PerColl_negTracks_c2) {
          if (!selectionTrack(negTrack)) {
            continue;
          }
          KaIdMethod = -1;
          if (!selKaon(negTrack, KaIdMethod)) {
            continue;
          }
          if (checkTrackInList(negTrack, triggerTrackIndexList_c2)) {
            continue;
          } // make sure that negTrack is not a Trigger
          // both tracks are Kaon now;
          ME_FillRHCorrelationUS(caseNo, c1, triggerTrack, posTrack, negTrack);
        } // Second Track
      } // first Track
      // //01-End  -obtaining h-Phi correlation
    } // tigger Loop
    // 01-h-phi All Three Region Unlike Sign correlation
  }

  template <typename T>
  bool selPhiMeson(T track1, T track2)
  {
    daughter1.SetXYZM(track1.px(), track1.py(), track1.pz(), pdgDB->Mass(321)); // set the daughter1 4-momentum
    daughter2.SetXYZM(track2.px(), track2.py(), track2.pz(), pdgDB->Mass(321)); // set the daughter2
    mother = daughter1 + daughter2;                                             // calculate the mother 4-momentum;
    mass = mother.M();
    if (cfgPhiMassLow <= mass && mass <= cfgPhiMassUp) {
      return true;
    }
    return false;
  }

  template <typename U, typename T>
  void ME_FillPhiPhiCorrelation(const int& caseNo, const U& collision, T track1, T track2, T track3, T track4)
  {
    daughter1.SetXYZM(track1.px(), track1.py(), track1.pz(), pdgDB->Mass(321)); // set the daughter1 4-momentum
    daughter2.SetXYZM(track2.px(), track2.py(), track2.pz(), pdgDB->Mass(321)); // set the daughter2
    phi1mother = daughter1 + daughter2;                                         // calculate the mother 4-momentum;

    daughter1.SetXYZM(track3.px(), track3.py(), track3.pz(), pdgDB->Mass(321)); // set the daughter1 4-momentum
    daughter2.SetXYZM(track4.px(), track4.py(), track4.pz(), pdgDB->Mass(321)); // set the daughter2
    phi2mother = daughter1 + daughter2;                                         // calculate the mother 4-momentum;

    if (phi1mother.Pt() > phi2mother.Pt()) {
      leadMother = phi1mother;
      assocMother = phi2mother;
    } else {
      leadMother = phi2mother;
      assocMother = phi1mother;
    }

    double dphi = ComputeDeltaPhi(assocMother.Phi(), leadMother.Phi());
    double deta = assocMother.Eta() - leadMother.Eta();

    if (caseNo == 5) {
      ME_recoPhi.fill(HIST("ME_hPhiPhi_Full_LeadPhi_01_Mass"), leadMother.M());
      ME_recoPhi.fill(HIST("ME_hPhiPhi_Full_LeadPhi_02_P"), leadMother.P());
      ME_recoPhi.fill(HIST("ME_hPhiPhi_Full_LeadPhi_03_Pt"), leadMother.Pt());
      ME_recoPhi.fill(HIST("ME_hPhiPhi_Full_LeadPhi_04_Eta"), leadMother.Eta());
      ME_recoPhi.fill(HIST("ME_hPhiPhi_Full_LeadPhi_05_Phi"), leadMother.Phi());
      ME_recoPhi.fill(HIST("ME_hPhiPhi_Full_LeadPhi_06_Rapidity"), leadMother.Rapidity());
      ME_recoPhi.fill(HIST("ME_hPhiPhi_Full_AssoPhi_01_Mass"), assocMother.M());
      ME_recoPhi.fill(HIST("ME_hPhiPhi_Full_AssoPhi_02_P"), assocMother.P());
      ME_recoPhi.fill(HIST("ME_hPhiPhi_Full_AssoPhi_03_Pt"), assocMother.Pt());
      ME_recoPhi.fill(HIST("ME_hPhiPhi_Full_AssoPhi_04_Eta"), assocMother.Eta());
      ME_recoPhi.fill(HIST("ME_hPhiPhi_Full_AssoPhi_05_Phi"), assocMother.Phi());
      ME_recoPhi.fill(HIST("ME_hPhiPhi_Full_AssoPhi_06_Rapidity"), assocMother.Rapidity());
      ME_recoAnalysis.fill(HIST("ME_hPhiPhi_Full_dPhi_dEta"), dphi, deta);
      ME_recoAnalysis.fill(HIST("ME_hPhiPhi_Full_Analysis"), collision.centFT0C(), collision.posZ(), dphi, deta);
    }
    if (caseNo == 6 || caseNo == 7 || caseNo == 8) {
      if (4.0 < leadMother.Pt() && leadMother.Pt() < 8.0) {
        ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_Full_LeadPhi_01_Mass"), leadMother.M());
        ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_Full_LeadPhi_02_P"), leadMother.P());
        ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_Full_LeadPhi_03_Pt"), leadMother.Pt());
        ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_Full_LeadPhi_04_Eta"), leadMother.Eta());
        ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_Full_LeadPhi_05_Phi"), leadMother.Phi());
        ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_Full_LeadPhi_06_Rapidity"), leadMother.Rapidity());
        ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_Full_AssoPhi_01_Mass"), assocMother.M());
        ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_Full_AssoPhi_02_P"), assocMother.P());
        ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_Full_AssoPhi_03_Pt"), assocMother.Pt());
        ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_Full_AssoPhi_04_Eta"), assocMother.Eta());
        ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_Full_AssoPhi_05_Phi"), assocMother.Phi());
        ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_Full_AssoPhi_06_Rapidity"), assocMother.Rapidity());
        ME_recoAnalysis.fill(HIST("ME_hPhiPhi_4to8LeadPhi_Full_dPhi_dEta"), dphi, deta);
        ME_recoAnalysis.fill(HIST("ME_hPhiPhi_4to8LeadPhi_Full_Analysis"), collision.centFT0C(), collision.posZ(), dphi, deta);
        if (caseNo == 7) {
          if (assocMother.Pt() < 2.0) {
            ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_01_Mass"), leadMother.M());
            ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_02_P"), leadMother.P());
            ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_03_Pt"), leadMother.Pt());
            ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_04_Eta"), leadMother.Eta());
            ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_05_Phi"), leadMother.Phi());
            ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_LeadPhi_06_Rapidity"), leadMother.Rapidity());
            ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_01_Mass"), assocMother.M());
            ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_02_P"), assocMother.P());
            ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_03_Pt"), assocMother.Pt());
            ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_04_Eta"), assocMother.Eta());
            ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_05_Phi"), assocMother.Phi());
            ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_AssoPhi_06_Rapidity"), assocMother.Rapidity());
            ME_recoAnalysis.fill(HIST("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_dPhi_dEta"), dphi, deta);
            ME_recoAnalysis.fill(HIST("ME_hPhiPhi_4to8LeadPhi_0to2AssoPhi_Analysis"), collision.centFT0C(), collision.posZ(), dphi, deta);
          }
        } else if (caseNo == 8) {
          if (2.0 < assocMother.Pt() && assocMother.Pt() < 4.0) {
            ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_01_Mass"), leadMother.M());
            ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_02_P"), leadMother.P());
            ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_03_Pt"), leadMother.Pt());
            ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_04_Eta"), leadMother.Eta());
            ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_05_Phi"), leadMother.Phi());
            ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_LeadPhi_06_Rapidity"), leadMother.Rapidity());
            ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_01_Mass"), assocMother.M());
            ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_02_P"), assocMother.P());
            ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_03_Pt"), assocMother.Pt());
            ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_04_Eta"), assocMother.Eta());
            ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_05_Phi"), assocMother.Phi());
            ME_recoPhi.fill(HIST("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_AssoPhi_06_Rapidity"), assocMother.Rapidity());
            ME_recoAnalysis.fill(HIST("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_dPhi_dEta"), dphi, deta);
            ME_recoAnalysis.fill(HIST("ME_hPhiPhi_4to8LeadPhi_2to4AssoPhi_Analysis"), collision.centFT0C(), collision.posZ(), dphi, deta);
          }
        }
      }
    } // case6, case 7 or case 8
  }

  void MIX_Collisions_Phi_Phi(int caseNo, int iMixBin, int iMixSlot1, int iMixSlot2)
  {
    LOG(info) << "DEBUG :: Mixing PhiPhi :: caseNo = " << caseNo << " :: iMixBin = " << iMixBin << " :: iMixSlot_ij = (" << iMixSlot1 << "," << iMixSlot2 << ")";

    auto const& c1 = myStoredMixingCollisions[caseNo][iMixBin][iMixSlot1];

    auto const& PerColl_triggerTrack_c1 = PerColl_triggerTracks[caseNo][iMixBin][iMixSlot1];
    auto const& PerColl_posTracks_c1 = PerColl_posTracks[caseNo][iMixBin][iMixSlot1];
    auto const& PerColl_negTracks_c1 = PerColl_negTracks[caseNo][iMixBin][iMixSlot1];

    auto const& PerColl_triggerTrack_c2 = PerColl_triggerTracks[caseNo][iMixBin][iMixSlot2];
    auto const& PerColl_posTracks_c2 = PerColl_posTracks[caseNo][iMixBin][iMixSlot2];
    auto const& PerColl_negTracks_c2 = PerColl_negTracks[caseNo][iMixBin][iMixSlot2];

    std::vector<int64_t> triggerTrackIndexList_c1;
    std::vector<int64_t> triggerTrackIndexList_c2;

    for (auto& triggerTrack : PerColl_triggerTrack_c1) {
      if (!selectionTrack(triggerTrack)) {
        continue;
      }
      triggerTrackIndexList_c1.push_back(triggerTrack.globalIndex());
    }

    for (auto& triggerTrack : PerColl_triggerTrack_c2) {
      if (!selectionTrack(triggerTrack)) {
        continue;
      }
      triggerTrackIndexList_c2.push_back(triggerTrack.globalIndex());
    }

    // 03-phi-phi correlation
    for (auto posTrack : PerColl_posTracks_c1) {
      int KaIdMethod = -1;
      if (!selectionTrack(posTrack)) {
        continue;
      } // FillSelectedTrackQAMIX(posTrack);
      if (!selKaon(posTrack, KaIdMethod)) {
        continue;
      } // fillKaonQAMIX(posTrack, KaIdMethod);
      if (checkTrackInList(posTrack, triggerTrackIndexList_c1)) {
        continue;
      } // make sure that posTrack is not a Trigger

      for (auto negTrack : PerColl_negTracks_c1) {
        KaIdMethod = -1;
        if (!selectionTrack(negTrack)) {
          continue;
        } // if( FillDebugger1 == 0 ){ FillDebugger1++;}//FillSelectedTrackQAMIX(negTrack)   ;}
        if (!selKaon(negTrack, KaIdMethod)) {
          continue;
        } // if( FillDebugger2 == 0 ){ FillDebugger2++;}//fillKaonQAMIX(negTrack, KaIdMethod);}
        if (checkTrackInList(posTrack, triggerTrackIndexList_c1)) {
          continue;
        } // make sure that posTrack is not a Trigger
        // Two Daughters are Kaon now;

        // FillPhiQA(posTrack, negTrack);//To Get US Results

        if (!selPhiMeson(posTrack, negTrack)) {
          continue;
        }

        // 03-Start-PhiPhi Correlation
        // Two phi mesons already identified
        for (auto posTrack2 : PerColl_posTracks_c2) {
          int KaIdMethod = -1;
          // if( posTrack2.globalIndex() <= posTrack.globalIndex()) {continue;}
          if (!selectionTrack(posTrack2)) {
            continue;
          }
          if (!selKaon(posTrack2, KaIdMethod)) {
            continue;
          }
          if (checkTrackInList(posTrack2, triggerTrackIndexList_c2)) {
            continue;
          } // make sure that posTrack2 is not a Trigger

          for (auto negTrack2 : PerColl_negTracks_c2) {
            KaIdMethod = -1;
            // if( negTrack2.globalIndex() <= negTrack.globalIndex()) {continue;}
            if (!selectionTrack(negTrack2)) {
              continue;
            }
            if (!selKaon(negTrack2, KaIdMethod)) {
              continue;
            }
            if (checkTrackInList(negTrack2, triggerTrackIndexList_c2)) {
              continue;
            } // make sure that negTrack2 is not a Trigger

            if (!selPhiMeson(posTrack2, negTrack2)) {
              continue;
            }
            // nPhiPhi++;
            ME_FillPhiPhiCorrelation(caseNo, c1, posTrack, negTrack, posTrack2, negTrack2);
          } // negTrack2
        } // posTrack2
        // 03-End  -PhiPhi Correlation
      } // negTrack
    } // posTrack
    // phi-phi correlation is over
  }

  template <typename U, typename T>
  void ME_Fill_hh_0To2AssoHadron(U collision, T triggerTrack, T assocTrack)
  {
    // hadron_pT in [0.0-2.0](Bulk)
    ME_recoAnalysis.fill(HIST("ME_hAnalysis_hh_0To2AssoHadron_LeadH_02_P"), triggerTrack.p());
    ME_recoAnalysis.fill(HIST("ME_hAnalysis_hh_0To2AssoHadron_LeadH_03_Pt"), triggerTrack.pt());
    ME_recoAnalysis.fill(HIST("ME_hAnalysis_hh_0To2AssoHadron_LeadH_04_Eta"), triggerTrack.eta());
    ME_recoAnalysis.fill(HIST("ME_hAnalysis_hh_0To2AssoHadron_LeadH_05_Phi"), triggerTrack.phi());
    ME_recoAnalysis.fill(HIST("ME_hAnalysis_hh_0To2AssoHadron_AssoH_02_P"), assocTrack.p());
    ME_recoAnalysis.fill(HIST("ME_hAnalysis_hh_0To2AssoHadron_AssoH_03_Pt"), assocTrack.pt());
    ME_recoAnalysis.fill(HIST("ME_hAnalysis_hh_0To2AssoHadron_AssoH_04_Eta"), assocTrack.eta());
    ME_recoAnalysis.fill(HIST("ME_hAnalysis_hh_0To2AssoHadron_AssoH_05_Phi"), assocTrack.phi());
    ME_recoAnalysis.fill(HIST("ME_hAnalysis_hh_0To2AssoHadron_dPhi_dEta"), ComputeDeltaPhi(assocTrack.phi(), triggerTrack.phi()), assocTrack.eta() - triggerTrack.eta());
    ME_recoAnalysis.fill(HIST("ME_hAnalysis_hh_0To2AssoHadron_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(assocTrack.phi(), triggerTrack.phi()), assocTrack.eta() - triggerTrack.eta());
  }

  template <typename U, typename T>
  void ME_Fill_hh_2To4AssoHadron(U collision, T triggerTrack, T assocTrack)
  {
    // hadron_pT in [2.0-4.0](Bulk)
    ME_recoAnalysis.fill(HIST("ME_hAnalysis_hh_2To4AssoHadron_LeadH_02_P"), triggerTrack.p());
    ME_recoAnalysis.fill(HIST("ME_hAnalysis_hh_2To4AssoHadron_LeadH_03_Pt"), triggerTrack.pt());
    ME_recoAnalysis.fill(HIST("ME_hAnalysis_hh_2To4AssoHadron_LeadH_04_Eta"), triggerTrack.eta());
    ME_recoAnalysis.fill(HIST("ME_hAnalysis_hh_2To4AssoHadron_LeadH_05_Phi"), triggerTrack.phi());
    ME_recoAnalysis.fill(HIST("ME_hAnalysis_hh_2To4AssoHadron_AssoH_02_P"), assocTrack.p());
    ME_recoAnalysis.fill(HIST("ME_hAnalysis_hh_2To4AssoHadron_AssoH_03_Pt"), assocTrack.pt());
    ME_recoAnalysis.fill(HIST("ME_hAnalysis_hh_2To4AssoHadron_AssoH_04_Eta"), assocTrack.eta());
    ME_recoAnalysis.fill(HIST("ME_hAnalysis_hh_2To4AssoHadron_AssoH_05_Phi"), assocTrack.phi());
    ME_recoAnalysis.fill(HIST("ME_hAnalysis_hh_2To4AssoHadron_dPhi_dEta"), ComputeDeltaPhi(assocTrack.phi(), triggerTrack.phi()), assocTrack.eta() - triggerTrack.eta());
    ME_recoAnalysis.fill(HIST("ME_hAnalysis_hh_2To4AssoHadron_Analysis"), collision.centFT0C(), collision.posZ(), ComputeDeltaPhi(assocTrack.phi(), triggerTrack.phi()), assocTrack.eta() - triggerTrack.eta());
  }

  void Mix_Collisions_h_h(int caseNo, int iMixBin, int iMixSlot1, int iMixSlot2)
  {
    // LOG(info)<<"DEBUG :: Mixing hh     :: caseNo = "<<caseNo<<" :: iMixBin = "<<iMixBin<<" :: iMixSlot_ij = ("<<iMixSlot1<<","<<iMixSlot2<<")";
    auto const& c1 = myStoredMixingCollisions[caseNo][iMixBin][iMixSlot1];

    auto const& PerColl_triggerTrack_c1 = PerColl_triggerTracks[caseNo][iMixBin][iMixSlot1];

    auto const& PerColl_triggerTrack_c2 = PerColl_triggerTracks[caseNo][iMixBin][iMixSlot2];
    auto const& PerColl_posTracks_c2 = PerColl_posTracks[caseNo][iMixBin][iMixSlot2];
    auto const& PerColl_negTracks_c2 = PerColl_negTracks[caseNo][iMixBin][iMixSlot2];

    auto const& PerColl_associatedTracks_0To2_c2 = PerColl_associatedTracks_0To2[caseNo][iMixBin][iMixSlot2];
    auto const& PerColl_associatedTracks_2To4_c2 = PerColl_associatedTracks_2To4[caseNo][iMixBin][iMixSlot2];

    std::vector<int64_t> triggerTrackIndexList_c2;
    for (auto& triggerTrack : PerColl_triggerTrack_c2) {
      if (!selectionTrack(triggerTrack)) {
        continue;
      }
      triggerTrackIndexList_c2.push_back(triggerTrack.globalIndex());
    }

    std::vector<int64_t> PhiPosDauKaonList_c2;
    std::vector<int64_t> PhiNegDauKaonList_c2;

    // GetPhiDaughters of second Event // to reject them in h-h correlation
    for (auto posTrack : PerColl_posTracks_c2) {
      int KaIdMethod = -1;
      if (!selectionTrack(posTrack)) {
        continue;
      }
      if (!selKaon(posTrack, KaIdMethod)) {
        continue;
      }
      if (checkTrackInList(posTrack, triggerTrackIndexList_c2)) {
        continue;
      } // make sure that posTrack is not a Trigger
      for (auto negTrack : PerColl_negTracks_c2) {
        KaIdMethod = -1;
        if (!selectionTrack(negTrack)) {
          continue;
        }
        if (!selKaon(negTrack, KaIdMethod)) {
          continue;
        }
        if (checkTrackInList(negTrack, triggerTrackIndexList_c2)) {
          continue;
        } // make sure that posTrack is not a Trigger

        if (!selPhiMeson(posTrack, negTrack)) {
          continue;
        }
        PhiPosDauKaonList_c2.push_back(posTrack.globalIndex());
        PhiNegDauKaonList_c2.push_back(negTrack.globalIndex());
      } // negTrack
    } // posTrack

    std::vector<int64_t> PhiDauKaonList_c2;
    FillNewListFromOldList(PhiDauKaonList_c2, PhiPosDauKaonList_c2);
    FillNewListFromOldList(PhiDauKaonList_c2, PhiNegDauKaonList_c2);
    InsertionSortVector(PhiDauKaonList_c2);

    // 05-h-h correlation
    for (auto& triggerTrack : PerColl_triggerTrack_c1) {
      if (!selectionTrack(triggerTrack)) {
        continue;
      }

      for (auto assocTrack : PerColl_associatedTracks_0To2_c2) {
        if (!selectionTrack(assocTrack)) {
          continue;
        }
        if (checkTrackInList(assocTrack, PhiDauKaonList_c2)) {
          continue;
        } // reject phi daughters
        ME_Fill_hh_0To2AssoHadron(c1, triggerTrack, assocTrack);
      } // Low  pT associated track

      for (auto assocTrack : PerColl_associatedTracks_2To4_c2) {
        if (!selectionTrack(assocTrack)) {
          continue;
        }
        if (checkTrackInList(assocTrack, PhiDauKaonList_c2)) {
          continue;
        } // reject phi daughters
        ME_Fill_hh_2To4AssoHadron(c1, triggerTrack, assocTrack);
      } // High pT associated track
    } // tigger Loop
    // 05-h-h correlation
  }

  void MIX_Collisions_Structure(int caseNo, int iMixBin, int iMixSlot1, int iMixSlot2)
  {
    // LOG(info)<<"DEBUG :: In Mixing Structure :: MixingOf bin = "<<iMixBin<<" :: ("<<iMixSlot1<<","<<iMixSlot2<<")";
    if (caseNo == 0 || caseNo == 1 || caseNo == 2 || caseNo == 3 || caseNo == 4 || caseNo == 11 || caseNo == 12 || caseNo == 13 || caseNo == 14) {
      MIX_Collisions_h_Phi(caseNo, iMixBin, iMixSlot1, iMixSlot2);
    } else if (caseNo == 5 || caseNo == 6 || caseNo == 7 || caseNo == 8) {
      MIX_Collisions_Phi_Phi(caseNo, iMixBin, iMixSlot1, iMixSlot2);
    } else if (caseNo == 9 || caseNo == 10) {
      Mix_Collisions_h_h(caseNo, iMixBin, iMixSlot1, iMixSlot2);
    }
  }

  template <typename T, typename U>
  void SetCollision(int dfNumber, std::string df_Name, int run_Number, int collBin, T& myColl, const U& myOriginalColl)
  {
    myColl.Bin = collBin;
    myColl.PosZ = myOriginalColl.posZ();
    myColl.CentFT0C = myOriginalColl.centFT0C();
    myColl.GlobalIndex = myOriginalColl.globalIndex();
    myColl.RunNumber = run_Number;
    myColl.DfName = df_Name;
    myColl.DfNo = dfNumber;
  }

  template <typename T, typename U>
  void SetTrack(int dfNumber, std::string df_Name, int run_Number, T& myTrack, const U& TableTrack)
  {
    myTrack.GlobalIndex = TableTrack.globalIndex();
    myTrack.CollId = TableTrack.collisionId();
    myTrack.Px = TableTrack.px();
    myTrack.Py = TableTrack.py();
    myTrack.Pz = TableTrack.pz();
    myTrack.P = TableTrack.p();
    myTrack.Pt = TableTrack.pt();
    myTrack.Signed1Pt = TableTrack.signed1Pt();
    myTrack.Eta = TableTrack.eta();
    myTrack.Phi = TableTrack.phi();
    myTrack.DcaXY = TableTrack.dcaXY();
    myTrack.DcaZ = TableTrack.dcaZ();
    myTrack.HasTOF = TableTrack.hasTOF();
    myTrack.Beta = TableTrack.beta();
    myTrack.TpcSignal = TableTrack.tpcSignal();
    myTrack.TpcInnerParam = TableTrack.tpcInnerParam();
    myTrack.TofExpMom = TableTrack.tofExpMom();
    myTrack.TpcNClsCrossedRows = TableTrack.tpcNClsCrossedRows();
    myTrack.IsGlobalTrack = TableTrack.isGlobalTrack();
    myTrack.Sign = TableTrack.sign();
    myTrack.TpcNSigmaPi = TableTrack.tpcNSigmaPi();
    myTrack.TpcNSigmaKa = TableTrack.tpcNSigmaKa();
    myTrack.TpcNSigmaPr = TableTrack.tpcNSigmaPr();
    myTrack.TpcNSigmaEl = TableTrack.tpcNSigmaEl();
    myTrack.TpcNSigmaDe = TableTrack.tpcNSigmaDe();
    myTrack.TofNSigmaPi = TableTrack.tofNSigmaPi();
    myTrack.TofNSigmaKa = TableTrack.tofNSigmaKa();
    myTrack.TofNSigmaPr = TableTrack.tofNSigmaPr();
    myTrack.TofNSigmaEl = TableTrack.tofNSigmaEl();
    myTrack.TofNSigmaDe = TableTrack.tofNSigmaDe();
    myTrack.DfNo = dfNumber;
    myTrack.DfName = df_Name;
    myTrack.RunNumber = run_Number;
  }

  template <typename U, typename V>
  void MixColl_Case(const int& caseNo, int dfNumber, std::string df_Name, int run_Number, int collBin, int mixingSlotPosition, U myOriginalColl, V PerColl_triggerTrack_c1, V PerColl_posTracks_c1, V PerColl_negTracks_c1, V PerColl_associatedTracks_0To2_c1, V PerColl_associatedTracks_2To4_c1)
  {

    StoredColl myColl;
    SetCollision(dfNumber, df_Name, run_Number, collBin, myColl, myOriginalColl);

    // Storing Data for Mixing - Begin
    // Clear Tracks to update vector
    PerColl_triggerTracks[caseNo][collBin][mixingSlotPosition].clear();
    PerColl_posTracks[caseNo][collBin][mixingSlotPosition].clear();
    PerColl_negTracks[caseNo][collBin][mixingSlotPosition].clear();
    PerColl_negTracks[caseNo][collBin][mixingSlotPosition].clear();
    PerColl_associatedTracks_0To2[caseNo][collBin][mixingSlotPosition].clear();
    PerColl_associatedTracks_2To4[caseNo][collBin][mixingSlotPosition].clear();

    // Store the collision at proper position slot
    if (myStoredMixingCollisions[caseNo][collBin].size() <= nMixEvt) {
      // Fill untill 6 events are filled for mixing
      LOG(info) << "DEBUG :: Mixing :: caseNumber = " << caseNo << " :: Only Fill Slot"
                << " :: bin = " << collBin;

      myStoredMixingCollisions[caseNo][collBin].push_back(myColl);

      myStoredTrack trigTrack;
      for (auto const& triggerTrack : PerColl_triggerTrack_c1) {
        SetTrack(dfNumber, df_Name, run_Number, trigTrack, triggerTrack);
        PerColl_triggerTracks[caseNo][collBin][mixingSlotPosition].push_back(trigTrack);
      }
      myStoredTrack pTrack;
      for (auto const& posTrack : PerColl_posTracks_c1) {
        SetTrack(dfNumber, df_Name, run_Number, pTrack, posTrack);
        PerColl_posTracks[caseNo][collBin][mixingSlotPosition].push_back(pTrack);
      }
      myStoredTrack nTrack;
      for (auto const& negTrack : PerColl_negTracks_c1) {
        SetTrack(dfNumber, df_Name, run_Number, nTrack, negTrack);
        PerColl_negTracks[caseNo][collBin][mixingSlotPosition].push_back(nTrack);
      }
      myStoredTrack assoTracks_0To2;
      for (auto const& associatedTracks_0To2 : PerColl_associatedTracks_0To2_c1) {
        SetTrack(dfNumber, df_Name, run_Number, assoTracks_0To2, associatedTracks_0To2);
        PerColl_associatedTracks_0To2[caseNo][collBin][mixingSlotPosition].push_back(assoTracks_0To2);
      }
      myStoredTrack assoTracks_2To4;
      for (auto const& associatedTracks_2To4 : PerColl_associatedTracks_2To4_c1) {
        SetTrack(dfNumber, df_Name, run_Number, assoTracks_2To4, associatedTracks_2To4);
        PerColl_associatedTracks_2To4[caseNo][collBin][mixingSlotPosition].push_back(assoTracks_2To4);
      }
    } else {
      // Replace the Events
      //  LOG(info)<<"DEBUG :: Mixing :: caseNumber = "<<caseNo<<" :: Replace Fill Slot"<<" :: bin = "<<collBin;
      myStoredMixingCollisions[caseNo][collBin][mixingSlotPosition] = myColl;

      myStoredTrack trigTrack;
      for (auto const& triggerTrack : PerColl_triggerTrack_c1) {
        SetTrack(dfNumber, df_Name, run_Number, trigTrack, triggerTrack);
        PerColl_triggerTracks[caseNo][collBin][mixingSlotPosition].push_back(trigTrack);
      }
      myStoredTrack pTrack;
      for (auto const& posTrack : PerColl_posTracks_c1) {
        SetTrack(dfNumber, df_Name, run_Number, pTrack, posTrack);
        PerColl_posTracks[caseNo][collBin][mixingSlotPosition].push_back(pTrack);
      }
      myStoredTrack nTrack;
      for (auto const& negTrack : PerColl_negTracks_c1) {
        SetTrack(dfNumber, df_Name, run_Number, nTrack, negTrack);
        PerColl_negTracks[caseNo][collBin][mixingSlotPosition].push_back(nTrack);
      }
      myStoredTrack assoTracks_0To2;
      for (auto const& associatedTracks_0To2 : PerColl_associatedTracks_0To2_c1) {
        SetTrack(dfNumber, df_Name, run_Number, assoTracks_0To2, associatedTracks_0To2);
        PerColl_associatedTracks_0To2[caseNo][collBin][mixingSlotPosition].push_back(assoTracks_0To2);
      }
      myStoredTrack assoTracks_2To4;
      for (auto const& associatedTracks_2To4 : PerColl_associatedTracks_2To4_c1) {
        SetTrack(dfNumber, df_Name, run_Number, assoTracks_2To4, associatedTracks_2To4);
        PerColl_associatedTracks_2To4[caseNo][collBin][mixingSlotPosition].push_back(assoTracks_2To4);
      }
    }

    if (int(PerColl_triggerTracks[caseNo][collBin][mixingSlotPosition].size()) != int(PerColl_triggerTrack_c1.size())) {
      LOG(info) << "DEBUG :: ERROR :: ERROR :: ERROR :: ERROR in triggerStoring";
    }
    if (int(PerColl_posTracks[caseNo][collBin][mixingSlotPosition].size()) != int(PerColl_posTracks_c1.size())) {
      LOG(info) << "DEBUG :: ERROR :: ERROR :: ERROR :: ERROR in posTracksStoring";
    }
    if (int(PerColl_negTracks[caseNo][collBin][mixingSlotPosition].size()) != int(PerColl_negTracks_c1.size())) {
      LOG(info) << "DEBUG :: ERROR :: ERROR :: ERROR :: ERROR in negTracksStoring";
    }
    if (int(PerColl_associatedTracks_0To2[caseNo][collBin][mixingSlotPosition].size()) != int(PerColl_associatedTracks_0To2_c1.size())) {
      LOG(info) << "DEBUG :: ERROR :: ERROR :: ERROR :: ERROR in associatedTracks_0To2Storing";
    }
    if (int(PerColl_associatedTracks_2To4[caseNo][collBin][mixingSlotPosition].size()) != int(PerColl_associatedTracks_2To4_c1.size())) {
      LOG(info) << "DEBUG :: ERROR :: ERROR :: ERROR :: ERROR in associatedTracks_2To4Storing";
    }
    // Storing Data for Mixing - End

    // Implement Collision Mixing - Begin
    if (myStoredMixingCollisions[caseNo][collBin].size() <= nMixEvt) {
      for (int i = 0; i < int(myStoredMixingCollisions[caseNo][collBin].size()) - 1; i++) {
        MIX_Collisions_Structure(caseNo, collBin, i, mixingSlotPosition);
      }
    } else {
      // for(int i = 0 ; i <= nMixEvt; i++){//   if(i == mixingSlotPosition){ continue;}
      for (int i : mixSlotPosList[mixingSlotPosition]) {
        MIX_Collisions_Structure(caseNo, collBin, i, mixingSlotPosition);
      }
    }
    // Implement Collision Mixing - End
  } // InsideMixingStructure

  // Event Filter
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutZvertex);

  // Track Filter
  Filter PtFilter = (o2::aod::track::pt) > cfgCutPt;
  Filter etaFilter = (nabs(aod::track::eta)) < cfgCutEta;
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using myCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels,
                                               aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>>; // aod::CentFV0As,

  using myTracks = soa::Filtered<soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TOFSignal, aod::pidTOFbeta, aod::pidTOFmass, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullDe, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullDe>>;

  // For manual sliceBy
  Preslice<myTracks> TracksPerCollisionPreslice = o2::aod::track::collisionId;

  // // definition of partitions
  SliceCache cache;
  Partition<myTracks> triggerTracks = 4.0f < aod::track::pt && aod::track::pt < 8.0f;
  Partition<myTracks> associatedTracks_0To2 = 0.0f < aod::track::pt && aod::track::pt < 2.0f;
  Partition<myTracks> associatedTracks_2To4 = 2.0f < aod::track::pt && aod::track::pt < 4.0f;
  Partition<myTracks> posTracks = aod::track::signed1Pt > 0.0f;
  Partition<myTracks> negTracks = aod::track::signed1Pt < 0.0f;

  using BinningTypeVtxZFT0C = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  BinningTypeVtxZFT0C colBinning{{CfgVtxBins, CfgMultBins}, true};

  int dfNumber = -1;
  // For DataFrameQA dont use iterator
  void processSameEvent(myCollisions const& collisions, myTracks const& fullTracks, o2::aod::Origins const& Origins, aod::BCsWithTimestamps const&)
  {
    dfNumber++;
    LOG(info) << "DEBUG :: df_" << dfNumber;

    int nTrack = 0;
    int nTrigger = 0;
    int nPhi = 0;
    int nPhiPhi = 0;
    double Phi1Pt = -1.0;
    double LSB1Pt = -1.0;
    double RSB1Pt = -1.0;
    double Phi2Pt = -1.0;

    int nPhi_0_2 = 0, nPhi_2_4 = 0, nPhi_4_8 = 0, nPhi_8_i = 0;
    int nLSB_0_2 = 0, nLSB_2_4 = 0, nLSB_4_8 = 0, nLSB_8_i = 0;
    int nRSB_0_2 = 0, nRSB_2_4 = 0, nRSB_4_8 = 0, nRSB_8_i = 0;
    int nLeadPhi = 0, nAssoPhi_0_2 = 0, nAssoPhi_2_4 = 0;
    int nAssoHad_0_2 = 0, nAssoHad_2_4 = 0;

    int nCR_Phi = 0;
    int nCR_Phi_0_2 = 0, nCR_Phi_2_4 = 0, nCR_Phi_4_8 = 0, nCR_Phi_8_i = 0;

    std::string df_Name = "DF_" + std::to_string(Origins.iteratorAt(0).dataframeID());
    auto bc = collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>();
    int currentRunNumber = bc.runNumber();
    LOG(info) << "DEBUG :: df_" << dfNumber << " :: ME :: currentRunNumber = " << currentRunNumber << " :: " << collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>().runNumber();

    TH1D* hBin00 = new TH1D("hBin00", "hBin00", 40, 0, 40);
    TH1D* hBin01 = new TH1D("hBin01", "hBin01", 40, 0, 40);
    TH1D* hBin02 = new TH1D("hBin02", "hBin02", 40, 0, 40);
    TH1D* hBin03 = new TH1D("hBin03", "hBin03", 40, 0, 40);
    TH1D* hBin04 = new TH1D("hBin04", "hBin04", 40, 0, 40);
    TH1D* hBin05 = new TH1D("hBin05", "hBin05", 40, 0, 40);
    TH1D* hBin06 = new TH1D("hBin06", "hBin06", 40, 0, 40);
    TH1D* hBin07 = new TH1D("hBin07", "hBin07", 40, 0, 40);
    TH1D* hBin08 = new TH1D("hBin08", "hBin08", 40, 0, 40);
    TH1D* hBin09 = new TH1D("hBin09", "hBin09", 40, 0, 40);
    TH1D* hBin10 = new TH1D("hBin10", "hBin10", 40, 0, 40);

    TH1D* hBin11 = new TH1D("hBin11", "hBin11", 40, 0, 40);
    TH1D* hBin12 = new TH1D("hBin12", "hBin12", 40, 0, 40);
    TH1D* hBin13 = new TH1D("hBin13", "hBin13", 40, 0, 40);
    TH1D* hBin14 = new TH1D("hBin14", "hBin14", 40, 0, 40);

    int bin = -1;
    std::vector<int> nthCaseVector;
    int mixingCounts = 0;

    for (auto collision : collisions) { // CollisionLoop-Start
      const auto tracks = fullTracks.sliceBy(TracksPerCollisionPreslice, collision.globalIndex());
      SE_recoEvent.fill(HIST("SE_hFullCollisionCount"), 0.5);

      nTrack = 0;
      for (auto track : tracks) {
        FillFullTrackQA(track);
        if (!selectionTrack(track)) {
          continue;
        }
        nTrack++;
        // FillSelectedTrackQA(track); //cant fill here, event might not have 3 tracks
      }
      if (nTrack < 3) {
        continue;
      } // {return;}

      nTrigger = 0;
      nPhi = 0;
      nPhiPhi = 0;

      Phi1Pt = -1.0;
      LSB1Pt = -1.0;
      RSB1Pt = -1.0;
      Phi2Pt = -1.0;

      nPhi_0_2 = 0;
      nPhi_2_4 = 0;
      nPhi_4_8 = 0;
      nPhi_8_i = 0;
      nLSB_0_2 = 0;
      nLSB_2_4 = 0;
      nLSB_4_8 = 0;
      nLSB_8_i = 0;
      nRSB_0_2 = 0;
      nRSB_2_4 = 0;
      nRSB_4_8 = 0;
      nRSB_8_i = 0;

      nLeadPhi = 0;
      nAssoPhi_0_2 = 0;
      nAssoPhi_2_4 = 0;
      nAssoHad_0_2 = 0;
      nAssoHad_2_4 = 0;

      nCR_Phi = 0;
      nCR_Phi_0_2 = 0;
      nCR_Phi_2_4 = 0;
      nCR_Phi_4_8 = 0;
      nCR_Phi_8_i = 0;

      // Get the Partitions
      auto triggerTracks_perColl = triggerTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      auto associatedTracks_0To2_perColl = associatedTracks_0To2->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      auto associatedTracks_2To4_perColl = associatedTracks_2To4->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      auto posTracks_perColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      auto negTracks_perColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

      // Option 1 :: We wont take phi meson if their daughter is one of the triggers
      // Option 2 :: Cheking if trigger is daughter of required phi and then skipping that trigger from correlation => if (PhiDaughterCheck(track)) { continue;}
      // We will go with option 1 for now and will check option2 later

      // 01-h-phi All Three Region Unlike Sign correlation
      std::vector<int64_t> triggerTrackIndexList; // long int
      nTrigger = 0;
      for (auto& triggerTrack : triggerTracks_perColl) {
        if (!selectionTrack(triggerTrack)) {
          continue;
        }
        SE_FillTriggerHadronQA(triggerTrack);
        nTrigger++;
        triggerTrackIndexList.push_back(triggerTrack.globalIndex());
        // 01-Start-obtaining h-Phi correlation
        for (auto posTrack : posTracks_perColl) {
          if (!selectionTrack(posTrack)) {
            continue;
          }
          int KaIdMethod = -1;
          if (!selKaon(posTrack, KaIdMethod)) {
            continue;
          } // Kaon PID Check
          if (matchIndices(posTrack, triggerTrack)) {
            continue;
          } // check posTrack is not a Trigger

          for (auto negTrack : negTracks_perColl) {
            if (!selectionTrack(negTrack)) {
              continue;
            }
            KaIdMethod = -1;
            if (!selKaon(negTrack, KaIdMethod)) {
              continue;
            }
            if (matchIndices(negTrack, triggerTrack)) {
              continue;
            } // check negTrack is not a Trigger

            // both tracks are Kaon now;
            SE_FillRHCorrelationUS(collision, triggerTrack, posTrack, negTrack, nCR_Phi, nCR_Phi_0_2, nCR_Phi_2_4, nCR_Phi_4_8, nCR_Phi_8_i);
          } // Second Track
        } // first Track
        // 01-End  -obtaining h-Phi correlation
      } // trigger Loop-End
      // 01-h-phi All Three Region Unlike Sign correlation
      InsertionSortVector(triggerTrackIndexList);

      // 02-Phi-meson reconstruction , US and LS PhiQA,  and  finding kaon daughter list
      // 03-phi-phi correlation are also inserted in this same loop
      std::vector<int64_t> PhiPosDauKaonList; // long int
      std::vector<int64_t> PhiNegDauKaonList; // long int
      nPhi = 0;
      nPhiPhi = 0;
      int fillSelectedTrack = 0;
      int fillKaonTrack = 0;
      for (auto posTrack : posTracks_perColl) {
        int KaIdMethod = -1;
        if (!selectionTrack(posTrack)) {
          continue;
        }
        FillSelectedTrackQA(posTrack);
        if (!selKaon(posTrack, KaIdMethod)) {
          continue;
        }
        fillKaonQA(posTrack, KaIdMethod);
        if (checkTrackInList(posTrack, triggerTrackIndexList)) {
          continue;
        } // make sure that posTrack is not a Trigger

        for (auto negTrack : negTracks_perColl) {
          KaIdMethod = -1;
          if (!selectionTrack(negTrack)) {
            continue;
          }
          if (fillSelectedTrack == 0) {
            FillSelectedTrackQA(negTrack);
          }
          if (!selKaon(negTrack, KaIdMethod)) {
            continue;
          }
          if (fillKaonTrack == 0) {
            fillKaonQA(negTrack, KaIdMethod);
          }
          if (checkTrackInList(negTrack, triggerTrackIndexList)) {
            continue;
          } // make sure that negTrack is not a Trigger

          // Two Daughters are Kaon now;
          FillPhiQA(posTrack, negTrack); // To Get US Results

          LSB1Pt = -1.0;
          RSB1Pt = -1.0;
          if (selLSBMeson(posTrack, negTrack, LSB1Pt)) {
            // nLSB++;
            if (LSB1Pt < 2.0) {
              nLSB_0_2++;
            }
            if (2.0 <= LSB1Pt && LSB1Pt < 4.0) {
              nLSB_2_4++;
            }
            if (4.0 <= LSB1Pt && LSB1Pt < 8.0) {
              nLSB_4_8++;
            }
            if (8.0 <= LSB1Pt) {
              nLSB_8_i++;
            }
          }

          if (selRSBMeson(posTrack, negTrack, RSB1Pt)) {
            // nRSB++;
            if (RSB1Pt < 2.0) {
              nRSB_0_2++;
            }
            if (2.0 <= RSB1Pt && RSB1Pt < 4.0) {
              nRSB_2_4++;
            }
            if (4.0 <= RSB1Pt && RSB1Pt < 8.0) {
              nRSB_4_8++;
            }
            if (8.0 <= RSB1Pt) {
              nRSB_8_i++;
            }
          }

          Phi1Pt = -1.0;
          if (!selPhiMeson(posTrack, negTrack, Phi1Pt)) {
            continue;
          }
          nPhi++;
          if (Phi1Pt < 2.0) {
            nPhi_0_2++;
          }
          if (2.0 <= Phi1Pt && Phi1Pt < 4.0) {
            nPhi_2_4++;
          }
          if (4.0 <= Phi1Pt && Phi1Pt < 8.0) {
            nPhi_4_8++;
          }
          if (8.0 <= Phi1Pt) {
            nPhi_8_i++;
          }

          PhiPosDauKaonList.push_back(posTrack.globalIndex()); // For repeat Kaon Counts
          PhiNegDauKaonList.push_back(negTrack.globalIndex()); // For repeat Kaon Counts

          // 03-Start-PhiPhi Correlation
          // Find second Phi meson - One phi meson already identified
          for (auto posTrack2 : posTracks_perColl) {
            int KaIdMethod = -1;
            if (posTrack2.globalIndex() <= posTrack.globalIndex()) {
              continue;
            }
            if (!selectionTrack(posTrack2)) {
              continue;
            }
            if (!selKaon(posTrack2, KaIdMethod)) {
              continue;
            }
            if (checkTrackInList(posTrack2, triggerTrackIndexList)) {
              continue;
            } // make sure that posTrack2 is not a Trigger

            for (auto negTrack2 : negTracks_perColl) {
              KaIdMethod = -1;
              if (negTrack2.globalIndex() <= negTrack.globalIndex()) {
                continue;
              }
              if (!selectionTrack(negTrack2)) {
                continue;
              }
              if (!selKaon(negTrack2, KaIdMethod)) {
                continue;
              }
              if (checkTrackInList(negTrack2, triggerTrackIndexList)) {
                continue;
              } // make sure that negTrack2 is not a Trigger

              Phi2Pt = -1.0;
              if (!selPhiMeson(posTrack2, negTrack2, Phi2Pt)) {
                continue;
              }
              nPhiPhi++;
              SE_FillPhiPhiCorrelation(collision, posTrack, negTrack, posTrack2, negTrack2, nLeadPhi, nAssoPhi_0_2, nAssoPhi_2_4);
            } // negTrack2
          } // posTrack2
          // 03-End  -PhiPhi Correlation
        } // negTrack
        fillSelectedTrack = 1;
        fillKaonTrack = 1;
      } // posTrack
      //
      FindRepeatEntries(PhiPosDauKaonList, SE_recoPhi.get<TH1>(HIST("SE_hPhi_nCommonKaonInDifferentPhi_US_Peak")));
      FindRepeatEntries(PhiNegDauKaonList, SE_recoPhi.get<TH1>(HIST("SE_hPhi_nCommonKaonInDifferentPhi_US_Peak")));

      // 04-LS KK Pair Plots
      for (auto negTrack1 : negTracks_perColl) {
        int KaIdMethod = -1;
        if (!selectionTrack(negTrack1)) {
          continue;
        }
        if (!selKaon(negTrack1, KaIdMethod)) {
          continue;
        }
        if (checkTrackInList(negTrack1, triggerTrackIndexList)) {
          continue;
        } // make sure that negTrack1 is not a Trigger
        auto negTrack1Id = negTrack1.globalIndex();

        for (auto negTrack2 : negTracks_perColl) {
          auto negTrack2Id = negTrack2.globalIndex();
          if (negTrack2Id <= negTrack1Id) {
            continue;
          }
          KaIdMethod = -1;
          if (!selectionTrack(negTrack2)) {
            continue;
          }
          if (!selKaon(negTrack2, KaIdMethod)) {
            continue;
          }
          if (checkTrackInList(negTrack2, triggerTrackIndexList)) {
            continue;
          } // make sure that negTrack2 is not a Trigger

          // Two Daughters are Kaon now;
          FillPhiQA(negTrack1, negTrack2);
          FillLS_Analysis(negTrack1, negTrack2);
        }
      }
      for (auto posTrack1 : posTracks_perColl) {
        int KaIdMethod = -1;
        if (!selectionTrack(posTrack1)) {
          continue;
        }
        if (!selKaon(posTrack1, KaIdMethod)) {
          continue;
        }
        if (checkTrackInList(posTrack1, triggerTrackIndexList)) {
          continue;
        } // make sure that posTrack1 is not a Trigger
        auto posTrack1Id = posTrack1.globalIndex();

        for (auto posTrack2 : posTracks_perColl) {
          auto posTrack2Id = posTrack2.globalIndex();
          if (posTrack2Id <= posTrack1Id) {
            continue;
          }
          KaIdMethod = -1;
          if (!selectionTrack(posTrack2)) {
            continue;
          }
          if (!selKaon(posTrack2, KaIdMethod)) {
            continue;
          }
          if (checkTrackInList(posTrack2, triggerTrackIndexList)) {
            continue;
          } // make sure that posTrack2 is not a Trigger

          // Two Daughters are Kaon now;
          FillPhiQA(posTrack1, posTrack2);
          FillLS_Analysis(posTrack1, posTrack2);
        }
      }
      // 04-LS KK Pair Plots

      std::vector<int64_t> PhiDauKaonList;
      FillNewListFromOldList(PhiDauKaonList, PhiPosDauKaonList);
      FillNewListFromOldList(PhiDauKaonList, PhiNegDauKaonList);
      InsertionSortVector(PhiDauKaonList);

      // 05-h-h correlation
      for (auto& triggerTrack : triggerTracks_perColl) {
        if (!selectionTrack(triggerTrack)) {
          continue;
        }

        for (auto assocTrack : associatedTracks_0To2_perColl) {
          if (!selectionTrack(assocTrack)) {
            continue;
          }
          if (checkTrackInList(assocTrack, PhiDauKaonList)) {
            continue;
          } // reject phi daughters
          SE_Fill_hh_0To2AssoHadron(collision, triggerTrack, assocTrack);
          nAssoHad_0_2++;
        } // Low  pT associated track

        for (auto assocTrack : associatedTracks_2To4_perColl) {
          if (!selectionTrack(assocTrack)) {
            continue;
          }
          if (checkTrackInList(assocTrack, PhiDauKaonList)) {
            continue;
          } // reject phi daughters
          SE_Fill_hh_2To4AssoHadron(collision, triggerTrack, assocTrack);
          nAssoHad_2_4++;
        } // High pT associated track

      } // tigger Loop
      // 05-h-h correlation-End

      SE_recoEvent.fill(HIST("SE_hCollisionCount"), 0.5);
      SE_recoEvent.fill(HIST("SE_hEvent_TrackSize"), tracks.size());
      SE_recoEvent.fill(HIST("SE_hVertexXRec"), collision.posX());
      SE_recoEvent.fill(HIST("SE_hVertexYRec"), collision.posY());
      SE_recoEvent.fill(HIST("SE_hVertexZRec"), collision.posZ());
      SE_recoEvent.fill(HIST("SE_hCentrality"), collision.centFT0C());
      SE_recoEvent.fill(HIST("SE_hCentrality_vs_vtxZ"), collision.centFT0C(), collision.posZ());

      SE_recoEvent.fill(HIST("SE_hEvent_nTrack"), nTrack);
      SE_recoEvent.fill(HIST("SE_hEvent_nTrigger"), nTrigger);
      SE_recoEvent.fill(HIST("SE_hEvent_nPhi"), nPhi);
      SE_recoEvent.fill(HIST("SE_hEvent_nPhiPhi"), nPhiPhi);

      SE_recoEvent.fill(HIST("SE_hEvent_nPhi_0_2"), nPhi_0_2);
      SE_recoEvent.fill(HIST("SE_hEvent_nPhi_2_4"), nPhi_2_4);
      SE_recoEvent.fill(HIST("SE_hEvent_nPhi_4_8"), nPhi_4_8);
      SE_recoEvent.fill(HIST("SE_hEvent_nPhi_8_i"), nPhi_8_i);
      SE_recoEvent.fill(HIST("SE_hEvent_nLSB_0_2"), nLSB_0_2);
      SE_recoEvent.fill(HIST("SE_hEvent_nLSB_2_4"), nLSB_2_4);
      SE_recoEvent.fill(HIST("SE_hEvent_nLSB_4_8"), nLSB_4_8);
      SE_recoEvent.fill(HIST("SE_hEvent_nLSB_8_i"), nLSB_8_i);
      SE_recoEvent.fill(HIST("SE_hEvent_nRSB_0_2"), nRSB_0_2);
      SE_recoEvent.fill(HIST("SE_hEvent_nRSB_2_4"), nRSB_2_4);
      SE_recoEvent.fill(HIST("SE_hEvent_nRSB_4_8"), nRSB_4_8);
      SE_recoEvent.fill(HIST("SE_hEvent_nRSB_8_i"), nRSB_8_i);

      SE_recoEvent.fill(HIST("SE_hEvent_nLeadPhi"), nLeadPhi);
      SE_recoEvent.fill(HIST("SE_hEvent_nAssoPhi_0_2"), nAssoPhi_0_2);
      SE_recoEvent.fill(HIST("SE_hEvent_nAssoPhi_2_4"), nAssoPhi_2_4);
      SE_recoEvent.fill(HIST("SE_hEvent_nAssoHad_0_2"), nAssoHad_0_2);
      SE_recoEvent.fill(HIST("SE_hEvent_nAssoHad_2_4"), nAssoHad_2_4);

      // ERROR Checks
      if (nTrigger * nPhi != nCR_Phi) {
        LOG(info) << "DEBUG :: SE :: ERROR :: nTrigger*nPhi     != nCR_Phi    ";
      }
      if (nTrigger * nPhi_0_2 != nCR_Phi_0_2) {
        LOG(info) << "DEBUG :: SE :: ERROR :: nTrigger*nPhi_0_2 != nCR_Phi_0_2";
      }
      if (nTrigger * nPhi_2_4 != nCR_Phi_2_4) {
        LOG(info) << "DEBUG :: SE :: ERROR :: nTrigger*nPhi_2_4 != nCR_Phi_2_4";
      }
      if (nTrigger * nPhi_4_8 != nCR_Phi_4_8) {
        LOG(info) << "DEBUG :: SE :: ERROR :: nTrigger*nPhi_4_8 != nCR_Phi_4_8";
      }
      if (nTrigger * nPhi_8_i != nCR_Phi_8_i) {
        LOG(info) << "DEBUG :: SE :: ERROR :: nTrigger*nPhi_8_i != nCR_Phi_8_i";
      }

      // Start Mixing Part

      // Find caseNo of the collision
      nthCaseVector.clear();

      if (nTrigger > 0) {
        // if      (nPhi > 0)                          { nthCaseVector.push_back( 0);} //caseNo =  0;
        if (nPhi_0_2 > 0) {
          nthCaseVector.push_back(1);
        } // caseNo =  1;
        if (nPhi_2_4 > 0) {
          nthCaseVector.push_back(2);
        } // caseNo =  2;
        // if  (nPhi_4_8 > 0)                          { nthCaseVector.push_back( 3);} //caseNo =  3;
        // if  (nPhi_8_i > 0)                          { nthCaseVector.push_back( 4);} //caseNo =  4;

        if (nLSB_0_2 > 0) {
          nthCaseVector.push_back(11);
        } // caseNo =  1;
        if (nLSB_2_4 > 0) {
          nthCaseVector.push_back(12);
        } // caseNo =  2;

        if (nRSB_0_2 > 0) {
          nthCaseVector.push_back(13);
        } // caseNo =  1;
        if (nRSB_2_4 > 0) {
          nthCaseVector.push_back(14);
        } // caseNo =  2;
      }
      if (nPhiPhi > 0) {
        nthCaseVector.push_back(5);
      } // caseNo =  5;
      // if(nLeadPhi > 0 && nPhiPhi > 0)             { nthCaseVector.push_back( 6);} //caseNo =  6;
      if (nLeadPhi > 0 && nAssoPhi_0_2 > 0) {
        nthCaseVector.push_back(7);
      } // caseNo =  7;
      if (nLeadPhi > 0 && nAssoPhi_2_4 > 0) {
        nthCaseVector.push_back(8);
      } // caseNo =  8;

      if (nTrigger > 0 && nAssoHad_0_2 > 0) {
        nthCaseVector.push_back(9);
      } // caseNo =  9;
      if (nTrigger > 0 && nAssoHad_2_4 > 0) {
        nthCaseVector.push_back(10);
      } // caseNo = 10;

      if (nthCaseVector.size() > 0) {

        // Find its mixing bin number;
        bin = colBinning.getBin({collision.posZ(), collision.centFT0C()});

        for (auto caseNo : nthCaseVector) { // Mixing case loop
          // if(caseNo ==  0) { hBin00->Fill(bin); ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBin00"),bin);}
          if (caseNo == 1) {
            hBin01->Fill(bin);
            ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBin01"), bin);
          }
          if (caseNo == 2) {
            hBin02->Fill(bin);
            ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBin02"), bin);
          }
          if (caseNo == 3) {
            hBin03->Fill(bin);
            ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBin03"), bin);
          }
          if (caseNo == 4) {
            hBin04->Fill(bin);
            ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBin04"), bin);
          }
          if (caseNo == 5) {
            hBin05->Fill(bin);
            ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBin05"), bin);
          }
          if (caseNo == 6) {
            hBin06->Fill(bin);
            ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBin06"), bin);
          }
          if (caseNo == 7) {
            hBin07->Fill(bin);
            ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBin07"), bin);
          }
          if (caseNo == 8) {
            hBin08->Fill(bin);
            ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBin08"), bin);
          }
          if (caseNo == 9) {
            hBin09->Fill(bin);
            ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBin09"), bin);
          }
          if (caseNo == 10) {
            hBin10->Fill(bin);
            ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBin10"), bin);
          }
          if (caseNo == 11) {
            hBin11->Fill(bin);
            ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBin11"), bin);
          }
          if (caseNo == 12) {
            hBin12->Fill(bin);
            ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBin12"), bin);
          }
          if (caseNo == 13) {
            hBin13->Fill(bin);
            ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBin13"), bin);
          }
          if (caseNo == 14) {
            hBin14->Fill(bin);
            ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBin14"), bin);
          }

          collRoll[caseNo][bin]++;
          mixingSlotPos[caseNo][bin]++;
          if (mixingSlotPos[caseNo][bin] > nMixEvt) {
            mixingSlotPos[caseNo][bin] = 0;
          }

          // You got the collisions, Start mixing now
          if (caseNo != 9 && caseNo != 10) {
            LOG(info) << "DEBUG :: df_" << dfNumber << " :: caseNumber = " << caseNo << " :: bin = " << bin << " :: collRoll = " << collRoll[caseNo][bin] << " :: mixingSlotPos = " << mixingSlotPos[caseNo][bin];
          }
          MixColl_Case(caseNo, dfNumber, df_Name, currentRunNumber, bin, mixingSlotPos[caseNo][bin], collision, triggerTracks_perColl, posTracks_perColl, negTracks_perColl, associatedTracks_0To2_perColl, associatedTracks_2To4_perColl);
          mixingCounts++;
        }
      } // Mixing Case Loop
    } // CollisionLoop-End

    // Store Dataframe information
    for (int iBin = 1; iBin <= 40; iBin++) {

      ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBinPerDF00"), iBin - 1, hBin00->GetBinContent(iBin));
      ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBinPerDF01"), iBin - 1, hBin01->GetBinContent(iBin));
      ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBinPerDF02"), iBin - 1, hBin02->GetBinContent(iBin));
      ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBinPerDF03"), iBin - 1, hBin03->GetBinContent(iBin));
      ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBinPerDF04"), iBin - 1, hBin04->GetBinContent(iBin));
      ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBinPerDF05"), iBin - 1, hBin05->GetBinContent(iBin));
      ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBinPerDF06"), iBin - 1, hBin06->GetBinContent(iBin));
      ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBinPerDF07"), iBin - 1, hBin07->GetBinContent(iBin));
      ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBinPerDF08"), iBin - 1, hBin08->GetBinContent(iBin));
      ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBinPerDF09"), iBin - 1, hBin09->GetBinContent(iBin));
      ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBinPerDF10"), iBin - 1, hBin10->GetBinContent(iBin));
      ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBinPerDF11"), iBin - 1, hBin11->GetBinContent(iBin));
      ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBinPerDF12"), iBin - 1, hBin12->GetBinContent(iBin));
      ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBinPerDF13"), iBin - 1, hBin13->GetBinContent(iBin));
      ME_recoAnalysis.fill(HIST("ME_hMixingEventsAvailabePerBinPerDF14"), iBin - 1, hBin14->GetBinContent(iBin));
    }

    delete hBin00;
    delete hBin01;
    delete hBin02;
    delete hBin03;
    delete hBin04;
    delete hBin05;
    delete hBin06;
    delete hBin07;
    delete hBin08;
    delete hBin09;
    delete hBin10;

    delete hBin11;
    delete hBin12;
    delete hBin13;
    delete hBin14;
  }
  PROCESS_SWITCH(hphicorrelation, processSameEvent, "Process Same event", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<hphicorrelation>(cfgc)};
}
