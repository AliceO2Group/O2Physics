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

/// \file electronSelectionWithTpcEmcal.cxx
/// \brief Task used to electron selection with tpc and emcal.
/// \author Rashi Gupta <rashi.gupta@cern.ch>, IIT Indore
/// \author Ravindra Singh <ravindra.singh@cern.ch>, IIT Indore

#include "PWGHF/HFL/DataModel/ElectronSelectionTable.h"
#include "PWGJE/DataModel/EMCALClusters.h"

#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Tools/KFparticle/KFUtilities.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TPDGCode.h>

#include <KFPTrack.h>
#include <KFParticle.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

using namespace o2;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

const int kEtaLocal = 221;

struct HfElectronSelectionWithTpcEmcal {

  Produces<aod::HfSelEl> electronSel;
  Produces<aod::HfCorrSelEl> hfElectronSelection;
  Produces<aod::HfMcGenSelEl> hfGenElectronSel;

  // select the emcal or dcal acceptance
  enum EMCalRegion {
    NoAcceptance = 0,
    EMCalAcceptance = 1,
    DCalAcceptance = 2
  };
  // Configurables
  // EMCal Cluster information
  KFParticle kfNonHfe;
  Configurable<bool> fillEmcClusterInfo{"fillEmcClusterInfo", true, "Fill histograms with EMCal cluster info before and after track match"};
  Configurable<bool> fillTrackInfo{"fillTrackInfo", true, "Fill histograms with Track Information info before track match"};
  Configurable<bool> skipNoEmcClusters{"skipNoEmcClusters", false, "Skip events with no EMCal clusters"};

  // select the emcal or dcal acceptance
  Configurable<int> emcalRegion{"emcalRegion", 0, "Select EMCal region for filling histograms (see EMCalRegion enum)"};

  // Event Selection
  Configurable<float> zPvPosMax{"zPvPosMax", 10., "Maximum z of the primary vertex (cm)"};
  Configurable<bool> isRun3{"isRun3", true, "Data is from Run3 or Run2"};

  // Track selection
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 0.5f, "DCA XY cut"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 1.0f, "DCA Z cut"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.6f, "Eta range for electron tracks"};
  Configurable<float> etaTrackMin{"etaTrackMin", -0.6f, "Eta range for electron tracks"};
  Configurable<float> ptTrackMin{"ptTrackMin", 3.0f, "Transverse MOmentum range for electron tracks"};

  // Associated electron selection cut
  Configurable<float> etaAssoTrackMax{"etaAssoTrackMax", 0.9f, "Eta range for Associatred electron tracks"};
  Configurable<float> etaAssoTrackMin{"etaAssoTrackMin", -0.9f, "Eta range for  Associatred  electron tracks"};
  Configurable<float> ptAssoTrackMin{"ptAssoTrackMin", 0.2f, "Transverse MOmentum range for  Associatred electron tracks"};
  Configurable<float> tpcNsigmaAssoElectronMin{"tpcNsigmaAssoElectronMin", -3.0f, "min Associated Electron TPCnsigma"};
  Configurable<float> tpcNsigmaAssoElectronMax{"tpcNsigmaAssoElectronMax", 3.0f, "max Associated Electron TPCnsigma"};
  Configurable<float> invariantMass{"invariantMass", 0.14f, "max Invariant Mass for Photonic electron"};
  Configurable<float> chiSquareMax{"chiSquareMax", 3.0f, "chiSquare on the reconstructed parent particle"};

  // EMcal and Dcal selection cut
  Configurable<float> etaTrackDCalNegativeMax{"etaTrackDCalNegativeMax", -0.22f, "Eta range for electron Dcal tracks"};
  Configurable<float> etaTrackDCalNegativeMin{"etaTrackDCalNegativeMin", -0.6f, "Eta range for electron tracks"};
  Configurable<float> etaTrackDCalPositiveMax{"etaTrackDCalPositiveMax", 0.6f, "Eta range for electron Dcal tracks"};
  Configurable<float> etaTrackDCalPositiveMin{"etaTrackDCalPositiveMin", 0.22f, "Eta range for electron tracks"};
  Configurable<float> phiTrackDCalMax{"phiTrackDCalMax", 5.708f, "phi range for electron tracks associated Dcal"};
  Configurable<float> phiTrackDCalMin{"phiTrackDCalMin", 4.5355f, "phi range for electron tracks associated Dcal"};
  Configurable<float> phiTrackEMCalMax{"phiTrackEMCalMax", 3.3621f, "phi range for electron tracks associated Emcal"};
  Configurable<float> phiTrackEMCalMin{"phiTrackEMCalMin", 1.3955f, "phi range for electron tracks associated Emcal"};

  // Track and  EMCal Cluster matching cut
  Configurable<float> deltaEtaMatchMin{"deltaEtaMatchMin", -0.013f, "Min Eta distance of EMCAL cluster to its closest track"};
  Configurable<float> deltaEtaMatchMax{"deltaEtaMatchMax", 0.0171f, "Max Eta distance of EMCAL cluster to its closest track"};
  Configurable<float> deltaPhiMatchMin{"deltaPhiMatchMin", -0.022f, "Min Phi distance of EMCAL cluster to its closest track"};
  Configurable<float> deltaPhiMatchMax{"deltaPhiMatchMax", 0.028f, "Max Phi distance of EMCAL cluster to its closest track"};
  Configurable<float> timeEmcClusterMax{"timeEmcClusterMax", 50.f, "EMCal Cluster time"};

  // Inclusive electron selection cut
  Configurable<float> eopElectronMin{"eopElectronMin", 0.8f, "Minimum E/p for electron tracks"};
  Configurable<float> eopElectronMax{"eopElectronMax", 1.2f, "Maximum E/p for electron tracks"};
  Configurable<float> m02EmcClusterElectronMax{"m02EmcClusterElectronMax", 0.9f, "max Electron  EMCal Cluster M02"};
  Configurable<float> m02EmcClusterElectronMin{"m02EmcClusterElectronMin", 0.02f, "min Electron  EMCal Cluster M02"};
  Configurable<float> m20EmcClusterElectronMax{"m20EmcClusterElectronMax", 1000.f, "max Electron  EMCal Cluster M20"};
  Configurable<float> m20EmcClusterElectronMin{"m20EmcClusterElectronMin", 0.0f, "min Electron  EMCal Cluster M20"};
  Configurable<float> tpcNsigmaElectronMin{"tpcNsigmaElectronMin", -0.5f, "min Electron TPCnsigma"};
  Configurable<float> tpcNsigmaElectronMax{"tpcNsigmaElectronMax", 3.0f, "max Electron TPCnsigma"};
  Configurable<float> tofNSigmaEl{"tofNSigmaEl", 3.0, "Sigma cut for electrons not in EMCal"};
  Configurable<int> pdgCodeCharmMin{"pdgCodeCharmMin", 400, "Min Charm Hadron PdgCode"};
  Configurable<int> pdgCodeCharmMax{"pdgCodeCharmMax", 600, "Max Charm Hadron PdgCode"};
  Configurable<int> pdgCodeBeautyMin{"pdgCodeBeautyMin", 4000, "Min beauty Hadron PdgCode"};
  Configurable<int> pdgCodeBeautyMax{"pdgCodeBeautyMax", 6000, "Max beauty Hadron PdgCode"};

  using TableCollisions = o2::soa::Filtered<o2::soa::Join<aod::Collisions, aod::Mults, aod::EvSels>>;
  using TableCollision = TableCollisions::iterator;
  using TableTracks = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksCov, o2::aod::TracksExtra, o2::aod::pidTPCFullEl, o2::aod::pidTOFFullEl, o2::aod::TracksDCA, o2::aod::TrackSelection, o2::aod::TrackSelectionExtension>;

  using McTableCollisions = o2::soa::Filtered<o2::soa::Join<TableCollisions, aod::McCollisionLabels>>;
  using McTableCollision = McTableCollisions::iterator;
  using McGenTableCollisions = soa::Join<aod::McCollisions, aod::MultsExtraMC>;
  using McGenTableCollision = McGenTableCollisions::iterator;
  using McTableTracks = soa::Join<TableTracks, aod::McTrackLabels>;
  using McTableEmcals = soa::Join<o2::aod::EMCALClusters, aod::EMCALMCClusters>;

  Filter collisionFilter = nabs(aod::collision::posZ) < zPvPosMax && aod::collision::numContrib > static_cast<uint16_t>(1);
  PresliceUnsorted<o2::aod::EMCALMatchedTracks> perClusterMatchedTracks = o2::aod::emcalmatchedtrack::trackId;

  // configurable axis
  ConfigurableAxis binsPosZ{"binsPosZ", {100, -10., 10.}, "primary vertex z coordinate"};
  ConfigurableAxis binsEta{"binsEta", {100, -2.0, 2.}, "#it{#eta}"};
  ConfigurableAxis binsPhi{"binsPhi", {32, 0.0, o2::constants::math::TwoPI}, "#it{#varphi}"};
  ConfigurableAxis binsPt{"binsPt", {50, 0.0, 50}, "#it{p_{T}}(GeV/#it{c})"};
  ConfigurableAxis binsdEdx{"binsdEdx", {160, 0., 160.}, "dE/dX"};
  ConfigurableAxis binsnSigma{"binsnSigma", {30, -15., 15.}, "#it{#sigma_{TPC}}"};
  ConfigurableAxis binsM02{"binsM02", {50, 0., 2.0}, "M02; entries"};
  ConfigurableAxis binsM20{"binsM20", {50, 0., 2.0}, "M20; entries"};
  ConfigurableAxis binsEoP{"binsEoP", {30, 0., 3.}, "e/p"};
  ConfigurableAxis binsEmcEnergy{"binsEmcEnergy", {50, 0., 50.}, "Cluster Energy (GeV/#it{c}^{2})"};
  ConfigurableAxis binsEmcClsNCells{"binsEmcClsNCells", {50, 0., 50.}, "nCells"};
  ConfigurableAxis binsEmcClsTime{"binsEmcClsTime", {1800, -900.0, 900.}, "Cluster Time"};
  ConfigurableAxis binsPassEMcal{"binsPassEMcal", {3, 0.0, 3.}, "Pass EMcal"};

  ConfigurableAxis binsDeltaEta{"binsDeltaEta", {20, -0.2, 0.2}, "Track Cluser Match #Delta #eta"};
  ConfigurableAxis binsDeltaPhi{"binsDeltaPhi", {20, -0.2, 0.2}, "Track Cluser Match #Delta #varphi"};
  ConfigurableAxis binsMass{"binsMass", {100, 0.0, 2.0}, "Mass (GeV/#it{c}^{2}); entries"};

  HistogramRegistry registry{
    "registry",
    {}};

  void init(o2::framework::InitContext&)
  {
    AxisSpec const axisPosZ = {binsPosZ, "Pos Z"};
    AxisSpec axisMass = {binsMass, "Mass (GeV/#it{c}^{2}); entries"};
    AxisSpec axisPt = {binsPt, "#it{p_{T}}(GeV/#it{c})"};
    AxisSpec axisEta = {binsEta, "#it{#eta}"};
    AxisSpec axisPhi = {binsPhi, "#it{#varphi}"};
    AxisSpec axisdEdx = {binsdEdx, "dE/dX"};
    AxisSpec axisnSigma = {binsnSigma, "it{#sigma_{TPC}}"};
    AxisSpec axisM02 = {binsM02, "M02; entries"};
    AxisSpec axisM20 = {binsM20, "M20; entries"};
    AxisSpec axisEoP = {binsEoP, "E/p"};
    AxisSpec axisEmcEnergy = {binsEmcEnergy, "Cluster Energy (GeV/#it{c}^{2})"};
    AxisSpec axisEmcClsNCells = {binsEmcClsNCells, "nCell"};
    AxisSpec axisEmcClsTime = {binsEmcClsTime, "Cluster Time"};
    AxisSpec axisPassEMcal = {binsPassEMcal, "Pass EMcal"};
    AxisSpec axisDeltaEta = {binsDeltaEta, "#Delta #eta = #eta_{trk}- #eta_{cluster}"};
    AxisSpec axisDeltaPhi = {binsDeltaPhi, "#Delta #varphi = #varphi_{trk}- #varphi_{cluster}"};

    registry.add("hZvertex", "z vertex", {HistType::kTH1D, {axisPosZ}});
    registry.add("hNeventsAfterPassEmcal", "No of events pass the Emcal", {HistType::kTH1D, {{3, 1, 4}}});
    registry.add("hNevents", "No of events", {HistType::kTH1D, {{3, 1, 4}}});
    registry.add("hLikeMass", "Like mass", {HistType::kTH1D, {{axisMass}}});
    registry.add("hUnLikeMass", "unLike mass", {HistType::kTH1D, {{axisMass}}});
    registry.add("hLikeSignPt", "Like sign Momentum ", {HistType::kTH1D, {{axisPt}}});
    registry.add("hUnLikeSignPt", "UnLike sign Momentum", {HistType::kTH1D, {{axisPt}}});
    registry.add("hMcgenInElectron", "Mc Gen Inclusive Electron", {HistType::kTH1D, {{axisPt}}});
    registry.add("hMcRecInElectron", "Mc Rec Inclusive Electron", {HistType::kTH1D, {{axisPt}}});
    registry.add("hMcRecwithoutEMCalInElectron", "Mc Rec Inclusive Electron without Emcal", {HistType::kTH1D, {{axisPt}}});

    registry.add("hMcgenAllNonHfeElectron", "Mc Gen All NonHf Electron", {HistType::kTH1D, {{axisPt}}});
    registry.add("hMcgenNonHfeElectron", "Mc Gen NonHf  Electron with mother", {HistType::kTH1D, {{axisPt}}});
    registry.add("hPi0eEmbTrkPt", "Mc Gen  Pi0 mother NonHf Electron", {HistType::kTH1D, {{axisPt}}});

    registry.add("hEtaeEmbTrkPt", "Mc Gen  Eta mother  NonHf Electron", {HistType::kTH1D, {{axisPt}}});
    registry.add("hEmcClusterM02", "m02", {HistType::kTH1D, {{axisM02}}});
    registry.add("hEmcClusterM20", "m20", {HistType::kTH1D, {{axisM20}}});
    registry.add("hTrackEtaPhi", "TPC EtaPhi Info; #eta;#varphi;passEMcal;", {HistType::kTH3F, {{axisEta}, {axisPhi}, {axisPassEMcal}}});
    registry.add("hTrackEnergyLossVsP", " TPC Energy loss info vs P; dE/dx;#it{p} (GeV#it{/c});passEMcal;", {HistType::kTH3F, {{axisdEdx}, {axisPt}, {axisPassEMcal}}});
    registry.add("hTrackEnergyLossVsPt", "TPC Energy loss info vs Pt; dE/dx;#it{p}_{T} (GeV#it{/c});passEMcal;", {HistType::kTH3F, {{axisdEdx}, {axisPt}, {axisPassEMcal}}});
    registry.add("hTracknSigmaVsP", " TPC nSigma info vs P; n#sigma;#it{p} (GeV#it{/c});passEMcal;", {HistType::kTH3F, {{axisnSigma}, {axisPt}, {axisPassEMcal}}});
    registry.add("hTracknSigmaVsPt", "  TPC nSigma info vs Pt; n#sigma;#it{p}_{T} (GeV#it{/c});passEMcal;", {HistType::kTH3F, {{axisnSigma}, {axisPt}, {axisPassEMcal}}});
    registry.add("hEmcClusterEnergy", "EMCal Cluster Info before match Energy; Energy (GeV); entries;", {HistType::kTH1D, {{axisEmcEnergy}}});
    registry.add("hEmcClusterEtaPhi", "EMCal Cluster Info before match Eta  and Phi; #eta;#varphi;", {HistType::kTH2F, {{axisEta}, {axisPhi}}});
    registry.add("hEmcClusterEnergyCell", "EMCal Cluster Info before match Energy vs nCells; Energy (GeV);ncell;", {HistType::kTH2F, {{axisEmcEnergy}, {axisEmcClsNCells}}});
    registry.add("hEmcClusterEnergyTime", "EMCal Cluster Info before match Energy vs time; Energy (GeV); sec;", {HistType::kTH2F, {{axisEmcEnergy}, {axisEmcClsTime}}});
    registry.add("hEmcClusterAfterMatchEnergy", "EMCal Cluster Info After match Energy; Energy (GeV); entries;", {HistType::kTH1D, {{axisEmcEnergy}}});
    registry.add("hEmcClusterAfterMatchEtaPhi", "EMCal Cluster Info After match Eta  and Phi; #eta;#varphi;", {HistType::kTH2F, {{axisEta}, {axisPhi}}});
    registry.add("hEmcClusterAfterMatchEnergyCells", "EMCal Cluster Info After match Energy vs nCells; Energy (GeV);ncell;", {HistType::kTH2F, {{axisEmcEnergy}, {axisEmcClsNCells}}});
    registry.add("hEmcClusterAfterMatchEnergyTime", "EMCal Cluster Info After match Energy vs time; Energy (GeV); sec;", {HistType::kTH2F, {{axisEmcEnergy}, {axisEmcClsTime}}});
    registry.add("hAfterMatchSigmaVsEoP", "PID Info after  match EoP vs Sigma ; E/P;#it{p}_{T} (GeV#it{/c});n#sigma; m02; m20;", {HistType::kTHnSparseF, {{axisEoP}, {axisPt}, {axisnSigma}, {axisM02}, {axisM20}}});
    registry.add("hAfterMatchEoPVsP", "PID Info after match  EoP vs P; E/P;#it{p} (GeV#it{/c});", {HistType::kTH2F, {{axisEoP}, {axisPt}}});
    registry.add("hAfterMatchSigmaVsP", "PID Info after match Sigma vs Momentum ; n#sigma; #it{p} (GeV#it{/c}; ", {HistType::kTH2F, {{axisnSigma}, {axisPt}}});
    registry.add("hAfterMatchEtaPhi", "PID Info after match Eta vs Phi ; #eta; #varphi; ", {HistType::kTH2F, {{axisEta}, {axisPhi}}});
    registry.add("hAfterMatchEnergyLossVsP", "PID Info after match Energy loss info vs P ; dE/dx;#it{p} (GeV#it{/c});; ", {HistType::kTH2F, {{axisdEdx}, {axisPt}}});
    registry.add("hAfterMatchEnergyLossVsPt", "PID Info after match Energy loss info vs Pt ;dE/dx;#it{p}_{T} (GeV#it{/c}); ", {HistType::kTH2F, {{axisdEdx}, {axisPt}}});

    registry.add("hAfterPIDEtaPhi", "PID Info after PID Cuts Eta vs Phi ; #eta; #varphi; ", {HistType::kTH2F, {{axisEta}, {axisPhi}}});
    registry.add("hEPRatioAfterPID", "E/P Ratio after PID Cuts apply only trackwodca filter", {HistType::kTH2F, {{axisPt}, {axisEmcEnergy}}});

    registry.add("hPIDAfterPIDCuts", "PID Info after PID cuts; E/P;#it{p}_{T} (GeV#it{/c});n#sigma;m02; m20;", {HistType::kTHnSparseF, {{axisEoP}, {axisPt}, {axisnSigma}, {axisM02}, {axisM20}}});
    registry.add("hEmcClsTrkEtaPhiDiffTime", "EmcClsTrkEtaPhiDiffTime;#Delta#eta;#Delta#varphi;Sec;", {HistType::kTH3F, {{axisDeltaEta}, {axisDeltaPhi}, {axisEmcClsTime}}});
    registry.add("hTofNSigmaVsPt", " TOF nSigma vs pt; n#sigma;#it{pt} (GeV/#it{c});", {HistType::kTH2F, {{axisnSigma}, {axisPt}}});
    registry.add("hTpcNSigmaVsPt", " TPC nSigma vs pt; n#sigma;#it{pt} (GeV/#it{c});", {HistType::kTH2F, {{axisnSigma}, {axisPt}}});
  }
  // Track Selection Cut
  template <typename T>
  bool selTracks(T const& track)
  {
    if (!track.isGlobalTrackWoDCA()) {
      return false;
    }
    if (std::abs(track.dcaXY()) > dcaXYTrackMax || std::abs(track.dcaZ()) > dcaZTrackMax) {
      return false;
    }
    if (track.eta() < etaTrackMin || track.eta() > etaTrackMax) {
      return false;
    }
    if ((track.phi() < phiTrackEMCalMin || track.phi() > phiTrackEMCalMax) && (track.phi() < phiTrackDCalMin || track.phi() > phiTrackDCalMax)) {
      return false;
    }
    if (track.pt() < ptTrackMin) {
      return false;
    }
    return true;
  }
  // Associated electron Selection Cut
  template <typename T>
  bool selAssoTracks(T const& track)
  {
    if (!track.isGlobalTrackWoDCA()) {
      return false;
    }
    if (std::abs(track.dcaXY()) > dcaXYTrackMax || std::abs(track.dcaZ()) > dcaZTrackMax) {
      return false;
    }
    if (track.eta() < etaAssoTrackMin || track.eta() > etaAssoTrackMax) {
      return false;
    }

    if (track.pt() < ptAssoTrackMin) {
      return false;
    }
    if (track.tpcNSigmaEl() < tpcNsigmaAssoElectronMin || track.tpcNSigmaEl() > tpcNsigmaAssoElectronMax) {
      return false;
    }

    return true;
  }

  // mc gen particle  selection cut
  template <typename T>
  bool mcGensel(T const& track)
  {
    if (track.eta() < etaTrackMin || track.eta() > etaTrackMax) {
      return false;
    }
    if ((track.phi() < phiTrackEMCalMin || track.phi() > phiTrackEMCalMax) && (track.phi() < phiTrackDCalMin || track.phi() > phiTrackDCalMax)) {
      return false;
    }
    if (track.pt() < ptTrackMin) {
      return false;
    }
    return true;
  }
  // nonHfe Identification

  template <typename ElectronType, typename TracksType>
  void nonHfe(ElectronType const& electron, TracksType const& tracks, bool isEMcal)
  {
    int nElPairsLS = 0;
    int nElPairsUS = 0;
    float invMassElectron = 0.;
    float massLike = 0;
    float massUnLike = 0;
    std::vector<float> vecLSMass;
    std::vector<float> vecULSMass;
    for (const auto& pTrack : tracks) {
      if (pTrack.globalIndex() == electron.globalIndex()) {
        continue;
      }
      // Apply partner electron selection

      if (!selAssoTracks(pTrack)) {
        continue;
      }
      if (electron.pt() <= pTrack.pt()) {
        continue;
      }
      int pdgE1 = kElectron;
      int pdgE2 = kElectron;
      if (electron.sign() > 0) {
        pdgE1 = kPositron;
      }

      if (pTrack.sign() > 0) {
        pdgE2 = kPositron;
      }

      KFPTrack const kfpTrack = createKFPTrackFromTrack(electron);
      KFPTrack const kfpAssociatedTrack = createKFPTrackFromTrack(pTrack);
      KFParticle const kfTrack(kfpTrack, pdgE1);
      KFParticle const kfAssociatedTrack(kfpAssociatedTrack, pdgE2);
      const KFParticle* electronPairs[2] = {&kfTrack, &kfAssociatedTrack};
      kfNonHfe.SetConstructMethod(2);
      kfNonHfe.Construct(electronPairs, 2);

      int const ndf = kfNonHfe.GetNDF();
      double const chi2recg = kfNonHfe.GetChi2() / ndf;
      if (ndf < 1.0) {
        continue;
      }

      if (std::sqrt(std::abs(chi2recg)) > chiSquareMax) {
        continue;
      }

      invMassElectron = RecoDecay::m(std::array{pTrack.pVector(), electron.pVector()}, std::array{MassElectron, MassElectron});
      bool isLSElectron = false;
      bool isULSElectron = false;
      // for like charge
      if (pTrack.sign() == electron.sign()) {
        massLike = invMassElectron;
        vecLSMass.push_back(massLike);
        isLSElectron = true;
        if (isEMcal) {
          registry.fill(HIST("hLikeMass"), massLike);
        }
      }
      // for unlike charge
      if (pTrack.sign() != electron.sign()) {
        massUnLike = invMassElectron;
        vecULSMass.push_back(massUnLike);
        isULSElectron = true;
        if (isEMcal) {
          registry.fill(HIST("hUnLikeMass"), massUnLike);
        }
      }

      // for like charge
      if (isLSElectron && (invMassElectron <= invariantMass)) {
        massLike = invMassElectron;
        ++nElPairsLS;
        if (isEMcal) {
          registry.fill(HIST("hLikeSignPt"), electron.pt());
        }
      }
      // for unlike charge
      if (isULSElectron && (invMassElectron <= invariantMass)) {
        massUnLike = invMassElectron;
        ++nElPairsUS;
        if (isEMcal) {
          registry.fill(HIST("hUnLikeSignPt"), electron.pt());
        }
      }
    }
    // Pass multiplicities and other required parameters for this electron
    // Pass multiplicities and other required parameters for this electron
    hfElectronSelection(electron.collisionId(), electron.globalIndex(), electron.eta(), electron.phi(), electron.pt(), electron.tpcNSigmaEl(), electron.tofNSigmaEl(), vecLSMass, vecULSMass, nElPairsLS, nElPairsUS, isEMcal);
  }
  // Electron Identification
  template <bool IsMc, typename TracksType, typename EmcClusterType, typename MatchType, typename CollisionType, typename ParticleType>
  void fillElectronTrack(CollisionType const& collision, TracksType const& tracks, EmcClusterType const& emcClusters, MatchType const& matchedTracks, ParticleType const&)
  {
    if (!(isRun3 ? collision.sel8() : (collision.sel7() && collision.alias_bit(kINT7)))) {
      return;
    }

    registry.fill(HIST("hNevents"), emcalRegion.value);

    // skip events with no clusters
    if (emcClusters.size() == 0 && skipNoEmcClusters) {
      return;
    }
    registry.fill(HIST("hZvertex"), collision.posZ());
    registry.fill(HIST("hNeventsAfterPassEmcal"), static_cast<int>(emcalRegion));
    /////////////////////////////////
    // EMCal cluster info before match ///
    ///////////////////////////////
    if (fillEmcClusterInfo) {
      for (const auto& emcClusterBefore : emcClusters) {
        registry.fill(HIST("hEmcClusterEnergy"), emcClusterBefore.energy());                                // track etaphi infor after filter bit
        registry.fill(HIST("hEmcClusterEtaPhi"), emcClusterBefore.eta(), emcClusterBefore.phi());           // track etaphi infor after filter bit
        registry.fill(HIST("hEmcClusterEnergyCell"), emcClusterBefore.energy(), emcClusterBefore.nCells()); // track etaphi infor after filter bit
        registry.fill(HIST("hEmcClusterEnergyTime"), emcClusterBefore.energy(), emcClusterBefore.time());   // track etaphi infor after filter bit
        registry.fill(HIST("hEmcClusterM02"), emcClusterBefore.m02());
        registry.fill(HIST("hEmcClusterM20"), emcClusterBefore.m20());
      }
    }
    EMCalRegion passEMCal = NoAcceptance;
    float phiTrack = -999;
    float etaTrack = -999;
    float pTrack = -999;
    float ptTrack = -999;
    float dcaxyTrack = -999;
    float dcazTrack = -999;
    float tpcNsigmaTrack = -999;

    for (const auto& track : tracks) {
      phiTrack = track.phi();
      etaTrack = track.eta();
      pTrack = track.p();
      ptTrack = track.pt();
      dcaxyTrack = track.dcaXY();
      dcazTrack = track.dcaZ();
      tpcNsigmaTrack = track.tpcNSigmaEl();
      // Apply Track Selection
      if (!selTracks(track)) {
        continue;
      }
      if ((phiTrack > phiTrackEMCalMin && phiTrack < phiTrackEMCalMax) && (etaTrack > etaTrackMin && etaTrack < etaTrackMax)) {
        passEMCal = EMCalAcceptance; // EMcal acceptance passed
      }
      if ((phiTrack > phiTrackDCalMin && phiTrack < phiTrackDCalMax) && ((etaTrack > etaTrackDCalPositiveMin && etaTrack < etaTrackDCalPositiveMax) || (etaTrack > etaTrackDCalNegativeMin && etaTrack < etaTrackDCalNegativeMax))) {
        passEMCal = DCalAcceptance; // Dcal acceptance passed
      }

      if (fillTrackInfo) {
        registry.fill(HIST("hTrackEtaPhi"), etaTrack, phiTrack, passEMCal);                 // track etaphi infor after filter bit
        registry.fill(HIST("hTrackEnergyLossVsP"), track.tpcSignal(), pTrack, passEMCal);   // track etaphi infor after filter bit
        registry.fill(HIST("hTrackEnergyLossVsPt"), track.tpcSignal(), ptTrack, passEMCal); // track etaphi infor after filter bit
        registry.fill(HIST("hTracknSigmaVsP"), tpcNsigmaTrack, pTrack, passEMCal);          // track etaphi infor after filter bit
        registry.fill(HIST("hTracknSigmaVsPt"), tpcNsigmaTrack, ptTrack, passEMCal);        // track etaphi infor after filter bit
      }
      auto tracksofcluster = matchedTracks.sliceBy(perClusterMatchedTracks, track.globalIndex());
      float phiMatchTrack = -999;
      float etaMatchTrack = -999;
      float pMatchTrack = -999;
      float ptMatchTrack = -999;
      float tpcNsigmaMatchTrack = -999;
      float phiMatchEmcCluster = -999;
      float etaMatchEmcCluster = -999;
      float eMatchEmcCluster = -999;
      float m02MatchEmcCluster = -999;
      float m20MatchEmcCluster = -999;
      float timeEmcCluster = -999;
      float cellEmcCluster = -999;
      float deltaPhiMatch = -999.;
      float deltaEtaMatch = -999.;
      float eop = -999;
      bool const isEMcal = false;

      float trackRapidity = track.rapidity(MassElectron);

      for (const auto& ematchTrack : tracksofcluster) {

        auto matchTrack = ematchTrack.template track_as<TracksType>();

        auto emcCluster = ematchTrack.template emcalcluster_as<EmcClusterType>();

        phiMatchTrack = matchTrack.phi();
        etaMatchTrack = matchTrack.eta();
        pMatchTrack = matchTrack.p();
        ptMatchTrack = matchTrack.pt();
        tpcNsigmaMatchTrack = matchTrack.tpcNSigmaEl();
        phiMatchEmcCluster = emcCluster.phi();
        etaMatchEmcCluster = emcCluster.eta();
        eMatchEmcCluster = emcCluster.energy();
        m02MatchEmcCluster = emcCluster.m02();
        m20MatchEmcCluster = emcCluster.m20();
        timeEmcCluster = emcCluster.time();
        cellEmcCluster = emcCluster.nCells();

        deltaPhiMatch = ematchTrack.deltaPhi();
        deltaEtaMatch = ematchTrack.deltaEta();

        // Track and EMCal cluster Matching
        if (std::abs(timeEmcCluster) > timeEmcClusterMax) {
          continue;
        }
        if (deltaPhiMatch < deltaPhiMatchMin || deltaPhiMatch > deltaPhiMatchMax || deltaEtaMatch < deltaEtaMatchMin || deltaEtaMatch > deltaEtaMatchMax) {
          continue;
        }

        registry.fill(HIST("hEmcClsTrkEtaPhiDiffTime"), deltaEtaMatch, deltaPhiMatch, timeEmcCluster);

        if (fillEmcClusterInfo) {
          registry.fill(HIST("hEmcClusterAfterMatchEnergy"), emcCluster.energy());                           // track etaphi infor after filter bit
          registry.fill(HIST("hEmcClusterAfterMatchEtaPhi"), emcCluster.eta(), emcCluster.phi());            // track etaphi infor after filter bit
          registry.fill(HIST("hEmcClusterAfterMatchEnergyCells"), emcCluster.energy(), emcCluster.nCells()); // track etaphi infor after filter bit
          registry.fill(HIST("hEmcClusterAfterMatchEnergyTime"), emcCluster.energy(), emcCluster.time());    // track etaphi infor after filter bit
        }

        eop = eMatchEmcCluster / pMatchTrack;

        registry.fill(HIST("hAfterMatchSigmaVsEoP"), eop, ptMatchTrack, tpcNsigmaMatchTrack, m02MatchEmcCluster, m20MatchEmcCluster);
        registry.fill(HIST("hAfterMatchEoPVsP"), eop, pMatchTrack);
        registry.fill(HIST("hAfterMatchSigmaVsP"), tpcNsigmaMatchTrack, pMatchTrack);
        registry.fill(HIST("hAfterMatchEtaPhi"), etaMatchTrack, phiMatchTrack);
        registry.fill(HIST("hAfterMatchEnergyLossVsP"), matchTrack.tpcSignal(), pMatchTrack);
        registry.fill(HIST("hAfterMatchEnergyLossVsPt"), matchTrack.tpcSignal(), ptMatchTrack);
        // Apply Electron Identification cuts

        if ((tpcNsigmaMatchTrack < tpcNsigmaElectronMin || tpcNsigmaMatchTrack > tpcNsigmaElectronMax) || (m02MatchEmcCluster < m02EmcClusterElectronMin || m02MatchEmcCluster > m02EmcClusterElectronMax) || (m20MatchEmcCluster < m20EmcClusterElectronMin || m20MatchEmcCluster > m20EmcClusterElectronMax)) {
          continue;
        }

        registry.fill(HIST("hPIDAfterPIDCuts"), eop, ptMatchTrack, tpcNsigmaMatchTrack, m02MatchEmcCluster, m20MatchEmcCluster);
        registry.fill(HIST("hEPRatioAfterPID"), pMatchTrack, eMatchEmcCluster);
        registry.fill(HIST("hAfterPIDEtaPhi"), etaMatchTrack, phiMatchTrack);
        if (eop < eopElectronMin || eop > eopElectronMax) {
          continue;
        }

        /////////////////          NonHf electron Selection with Emcal       ////////////////////////
        if constexpr (IsMc) {
          if (matchTrack.has_mcParticle()) {
            auto mcParticle = matchTrack.template mcParticle_as<aod::McParticles>();
            if (std::abs(mcParticle.pdgCode()) == kElectron) {

              registry.fill(HIST("hMcRecInElectron"), mcParticle.pt());
              bool isEmbEta = false;
              bool isEmbPi0 = false;

              // Check first mother
              if (mcParticle.has_mothers()) {
                auto const& mother = mcParticle.template mothers_first_as<aod::McParticles>();

                if (std::abs(mother.pdgCode()) == kEtaLocal || std::abs(mother.pdgCode()) == kPi0 || std::abs(mother.pdgCode()) == kGamma) {

                  auto const& gmother = mother.template mothers_first_as<aod::McParticles>();
                  // cases to consider: eta->e, eta->pi0->e, eta->gamma->e, eta->pi0->gamma->e, pi0->e, pi0->gamma->e

                  //=================  eta->e ======================================
                  if (std::abs(mother.pdgCode()) == kEtaLocal) {

                    if (mother.isPhysicalPrimary()) {
                      if ((std::abs(gmother.pdgCode()) >= pdgCodeCharmMin && std::abs(gmother.pdgCode()) < pdgCodeCharmMax) ||
                          (std::abs(gmother.pdgCode()) >= pdgCodeBeautyMin && std::abs(gmother.pdgCode()) < pdgCodeBeautyMax)) {
                        continue;
                      }
                      isEmbEta = true;
                    }
                  }

                  //=================  eta->pi0->e ======================================

                  if (std::abs(mother.pdgCode()) == kPi0) {
                    if (mother.isPhysicalPrimary()) {
                      if ((std::abs(gmother.pdgCode()) >= pdgCodeCharmMin && std::abs(gmother.pdgCode()) < pdgCodeCharmMax) ||
                          (std::abs(gmother.pdgCode()) >= pdgCodeBeautyMin && std::abs(gmother.pdgCode()) < pdgCodeBeautyMax)) {
                        continue;
                      }
                      isEmbPi0 = true; // pi0 -> e
                    }
                    if (std::abs(gmother.pdgCode()) == kEtaLocal) {
                      if (gmother.isPhysicalPrimary() || gmother.has_mothers()) {
                        auto const& ggmother = gmother.template mothers_first_as<aod::McParticles>();
                        if ((std::abs(ggmother.pdgCode()) >= pdgCodeCharmMin && std::abs(ggmother.pdgCode()) < pdgCodeCharmMax) ||
                            (std::abs(ggmother.pdgCode()) >= pdgCodeBeautyMin && std::abs(ggmother.pdgCode()) < pdgCodeBeautyMax)) {
                          continue;
                        }
                        isEmbEta = true; // eta->pi0-> e
                      }
                    }
                  }

                  /// ====================================  eta->gamma->e  and eta->pi0->gamma->e============
                  if (std::abs(mother.pdgCode()) == kGamma) {

                    if (std::abs(gmother.pdgCode()) == kEtaLocal) {
                      if (gmother.isPhysicalPrimary() || gmother.has_mothers()) {
                        auto const& ggmother = gmother.template mothers_first_as<aod::McParticles>();
                        if ((std::abs(ggmother.pdgCode()) >= pdgCodeCharmMin && std::abs(ggmother.pdgCode()) < pdgCodeCharmMax) ||
                            (std::abs(ggmother.pdgCode()) >= pdgCodeBeautyMin && std::abs(ggmother.pdgCode()) < pdgCodeBeautyMax)) {
                          continue;
                        }
                        isEmbEta = true; // eta->gamma-> e
                      }
                    }
                    if (std::abs(gmother.pdgCode()) == kPi0) {
                      if (gmother.isPhysicalPrimary() || gmother.has_mothers()) {
                        auto const& ggmother = gmother.template mothers_first_as<aod::McParticles>();
                        if ((std::abs(ggmother.pdgCode()) >= pdgCodeCharmMin && std::abs(ggmother.pdgCode()) < pdgCodeCharmMax) ||
                            (std::abs(ggmother.pdgCode()) >= pdgCodeBeautyMin && std::abs(ggmother.pdgCode()) < pdgCodeBeautyMax)) {
                          continue;
                        }
                        isEmbPi0 = true; // pi0-> gamma-> e
                      }
                      if (gmother.has_mothers()) {
                        auto const& ggmother = gmother.template mothers_first_as<aod::McParticles>();
                        if (std::abs(ggmother.pdgCode()) == kEtaLocal) {
                          if (ggmother.isPhysicalPrimary() || ggmother.has_mothers()) {
                            auto const& gggmother = ggmother.template mothers_first_as<aod::McParticles>();
                            if ((std::abs(gggmother.pdgCode()) >= pdgCodeCharmMin && std::abs(gggmother.pdgCode()) < pdgCodeCharmMax) ||
                                (std::abs(gggmother.pdgCode()) >= pdgCodeBeautyMin && std::abs(gggmother.pdgCode()) < pdgCodeBeautyMax)) {
                              continue;
                            }
                            isEmbEta = true; // eta->pi0->gamma-> e
                          }
                        }
                      }
                    }
                  }
                  if (!(isEmbPi0 && isEmbEta)) {
                    continue;
                  }
                }
              }
            }
          }
        }
        nonHfe(matchTrack, tracks, true);

        /////////////////          NonHf electron Selection without Emcal       ////////////////////////
        electronSel(track.collisionId(), track.globalIndex(), etaTrack, phiTrack, ptTrack, pTrack, trackRapidity, dcaxyTrack, dcazTrack, track.tpcNSigmaEl(), track.tofNSigmaEl(),
                    eMatchEmcCluster, etaMatchEmcCluster, phiMatchEmcCluster, m02MatchEmcCluster, m20MatchEmcCluster, cellEmcCluster, timeEmcCluster, deltaEtaMatch, deltaPhiMatch, isEMcal);
      }
      /// Electron information without Emcal and use TPC and TOF
      if (isEMcal) {
        continue;
      }
      if (std::abs(track.tofNSigmaEl()) > tofNSigmaEl) {
        continue;
      }
      registry.fill(HIST("hTofNSigmaVsPt"), track.tofNSigmaEl(), track.pt());
      registry.fill(HIST("hTpcNSigmaVsPt"), track.tpcNSigmaEl(), track.pt());

      if ((track.tpcNSigmaEl() < tpcNsigmaElectronMin || track.tpcNSigmaEl() > tpcNsigmaElectronMax)) {
        continue;
      }
      if constexpr (IsMc) {
        if (track.has_mcParticle()) {
          auto mcParticle = track.template mcParticle_as<aod::McParticles>();
          if (std::abs(mcParticle.pdgCode()) == kElectron) {

            registry.fill(HIST("hMcRecwithoutEMCalInElectron"), mcParticle.pt());
            bool isEmbEta = false;
            bool isEmbPi0 = false;

            // Check first mother
            if (mcParticle.has_mothers()) {
              auto const& mother = mcParticle.template mothers_first_as<aod::McParticles>();

              if (std::abs(mother.pdgCode()) == kEtaLocal || std::abs(mother.pdgCode()) == kPi0 || std::abs(mother.pdgCode()) == kGamma) {

                auto const& gmother = mother.template mothers_first_as<aod::McParticles>();
                // cases to consider: eta->e, eta->pi0->e, eta->gamma->e, eta->pi0->gamma->e, pi0->e, pi0->gamma->e

                //=================  eta->e ======================================
                if (std::abs(mother.pdgCode()) == kEtaLocal) {

                  if (mother.isPhysicalPrimary()) {
                    if ((std::abs(gmother.pdgCode()) >= pdgCodeCharmMin && std::abs(gmother.pdgCode()) < pdgCodeCharmMax) ||
                        (std::abs(gmother.pdgCode()) >= pdgCodeBeautyMin && std::abs(gmother.pdgCode()) < pdgCodeBeautyMax)) {
                      continue;
                    }
                    isEmbEta = true;
                  }
                }

                //=================  eta->pi0->e ======================================

                if (std::abs(mother.pdgCode()) == kPi0) {
                  if (mother.isPhysicalPrimary()) {
                    if ((std::abs(gmother.pdgCode()) >= pdgCodeCharmMin && std::abs(gmother.pdgCode()) < pdgCodeCharmMax) ||
                        (std::abs(gmother.pdgCode()) >= pdgCodeBeautyMin && std::abs(gmother.pdgCode()) < pdgCodeBeautyMax)) {
                      continue;
                    }
                    isEmbPi0 = true; // pi0 -> e
                  }
                  if (std::abs(gmother.pdgCode()) == kEtaLocal) {
                    if (gmother.isPhysicalPrimary() || gmother.has_mothers()) {
                      auto const& ggmother = gmother.template mothers_first_as<aod::McParticles>();
                      if ((std::abs(ggmother.pdgCode()) >= pdgCodeCharmMin && std::abs(ggmother.pdgCode()) < pdgCodeCharmMax) ||
                          (std::abs(ggmother.pdgCode()) >= pdgCodeBeautyMin && std::abs(ggmother.pdgCode()) < pdgCodeBeautyMax)) {
                        continue;
                      }
                      isEmbEta = true; // eta->pi0-> e
                    }
                  }
                }

                /// ====================================  eta->gamma->e  and eta->pi0->gamma->e============
                if (std::abs(mother.pdgCode()) == kGamma) {

                  if (std::abs(gmother.pdgCode()) == kEtaLocal) {
                    if (gmother.isPhysicalPrimary() || gmother.has_mothers()) {
                      auto const& ggmother = gmother.template mothers_first_as<aod::McParticles>();
                      if ((std::abs(ggmother.pdgCode()) >= pdgCodeCharmMin && std::abs(ggmother.pdgCode()) < pdgCodeCharmMax) ||
                          (std::abs(ggmother.pdgCode()) >= pdgCodeBeautyMin && std::abs(ggmother.pdgCode()) < pdgCodeBeautyMax)) {
                        continue;
                      }
                      isEmbEta = true; // eta->gamma-> e
                    }
                  }
                  if (std::abs(gmother.pdgCode()) == kPi0) {
                    if (gmother.isPhysicalPrimary() || gmother.has_mothers()) {
                      auto const& ggmother = gmother.template mothers_first_as<aod::McParticles>();
                      if ((std::abs(ggmother.pdgCode()) >= pdgCodeCharmMin && std::abs(ggmother.pdgCode()) < pdgCodeCharmMax) ||
                          (std::abs(ggmother.pdgCode()) >= pdgCodeBeautyMin && std::abs(ggmother.pdgCode()) < pdgCodeBeautyMax)) {
                        continue;
                      }
                      isEmbPi0 = true; // pi0-> gamma-> e
                    }
                    if (gmother.has_mothers()) {
                      auto const& ggmother = gmother.template mothers_first_as<aod::McParticles>();
                      if (std::abs(ggmother.pdgCode()) == kEtaLocal) {
                        if (ggmother.isPhysicalPrimary() || ggmother.has_mothers()) {
                          auto const& gggmother = ggmother.template mothers_first_as<aod::McParticles>();
                          if ((std::abs(gggmother.pdgCode()) >= pdgCodeCharmMin && std::abs(gggmother.pdgCode()) < pdgCodeCharmMax) ||
                              (std::abs(gggmother.pdgCode()) >= pdgCodeBeautyMin && std::abs(gggmother.pdgCode()) < pdgCodeBeautyMax)) {
                            continue;
                          }
                          isEmbEta = true; // eta->pi0->gamma-> e
                        }
                      }
                    }
                  }
                }
                if (!(isEmbPi0 && isEmbEta)) {
                  continue;
                }
              }
            }
          }
        }
      }
      nonHfe(track, tracks, false);
      /////////////////          NonHf electron Selection without Emcal       ////////////////////////
      electronSel(track.collisionId(), track.globalIndex(), etaTrack, phiTrack, ptTrack, pTrack, trackRapidity, dcaxyTrack, dcazTrack, track.tpcNSigmaEl(), track.tofNSigmaEl(),
                  eMatchEmcCluster, etaMatchEmcCluster, phiMatchEmcCluster, m02MatchEmcCluster, m20MatchEmcCluster, cellEmcCluster, timeEmcCluster, deltaEtaMatch, deltaPhiMatch, isEMcal);
    }
  }

  ///  Electron selection - for real data and data-like analysis
  void processData(TableCollision const& collision,
                   TableTracks const& tracks,
                   aod::EMCALClusters const& emcClusters,
                   o2::aod::EMCALMatchedTracks const& matchedTracks)
  {
    fillElectronTrack<false>(collision, tracks, emcClusters, matchedTracks, 0);
  }
  PROCESS_SWITCH(HfElectronSelectionWithTpcEmcal, processData, "process Data info only", true);
  ///  Electron selection - for MC reco-level analysis
  void processMcRec(McTableCollision const& mcCollision,
                    McTableTracks const& mcTracks,
                    McTableEmcals const& mcEmcClusters,
                    o2::aod::EMCALMatchedTracks const& matchedTracks,
                    aod::McParticles const& mcParticles)
  {
    fillElectronTrack<true>(mcCollision, mcTracks, mcEmcClusters, matchedTracks, mcParticles);
  }
  PROCESS_SWITCH(HfElectronSelectionWithTpcEmcal, processMcRec, "Process MC Reco mode", false);

  void processMcGen(McGenTableCollision const& mcCollision, aod::McParticles const& mcParticles)
  {

    bool isNonHfe = false;
    for (const auto& particleMc : mcParticles) {

      if (!mcGensel(particleMc)) {
        continue;
      }
      if (std::abs(particleMc.pdgCode()) == kElectron) {

        registry.fill(HIST("hMcgenInElectron"), particleMc.pt());
        bool isEmbEta = false;
        bool isEmbPi0 = false;

        // Check first mother
        if (particleMc.has_mothers()) {
          auto const& mother = particleMc.mothers_first_as<aod::McParticles>();

          if (std::abs(mother.pdgCode()) == kEtaLocal || std::abs(mother.pdgCode()) == kPi0 || std::abs(mother.pdgCode()) == kGamma) {
            registry.fill(HIST("hMcgenAllNonHfeElectron"), particleMc.pt());

            auto const& gmother = mother.mothers_first_as<aod::McParticles>();
            // cases to consider: eta->e, eta->pi0->e, eta->gamma->e, eta->pi0->gamma->e, pi0->e, pi0->gamma->e

            //=================  eta->e ======================================
            if (std::abs(mother.pdgCode()) == kEtaLocal) {

              if (mother.isPhysicalPrimary()) {
                if ((std::abs(gmother.pdgCode()) >= pdgCodeCharmMin && std::abs(gmother.pdgCode()) < pdgCodeCharmMax) ||
                    (std::abs(gmother.pdgCode()) >= pdgCodeBeautyMin && std::abs(gmother.pdgCode()) < pdgCodeBeautyMax)) {
                  continue;
                }
                isEmbEta = true;
              }
            }

            //=================  eta->pi0->e ======================================

            if (std::abs(mother.pdgCode()) == kPi0) {
              if (mother.isPhysicalPrimary()) {
                if ((std::abs(gmother.pdgCode()) >= pdgCodeCharmMin && std::abs(gmother.pdgCode()) < pdgCodeCharmMax) ||
                    (std::abs(gmother.pdgCode()) >= pdgCodeBeautyMin && std::abs(gmother.pdgCode()) < pdgCodeBeautyMax)) {
                  continue;
                }
                isEmbPi0 = true; // pi0 -> e
              }
              if (std::abs(gmother.pdgCode()) == kEtaLocal) {
                if (gmother.isPhysicalPrimary() || gmother.has_mothers()) {
                  auto const& ggmother = gmother.mothers_first_as<aod::McParticles>();
                  if ((std::abs(ggmother.pdgCode()) >= pdgCodeCharmMin && std::abs(ggmother.pdgCode()) < pdgCodeCharmMax) ||
                      (std::abs(ggmother.pdgCode()) >= pdgCodeBeautyMin && std::abs(ggmother.pdgCode()) < pdgCodeBeautyMax)) {
                    continue;
                  }
                  isEmbEta = true; // eta->pi0-> e
                }
              }
            }

            /// ====================================  eta->gamma->e  and eta->pi0->gamma->e============
            if (std::abs(mother.pdgCode()) == kGamma) {

              if (std::abs(gmother.pdgCode()) == kEtaLocal) {
                if (gmother.isPhysicalPrimary() || gmother.has_mothers()) {
                  auto const& ggmother = gmother.mothers_first_as<aod::McParticles>();
                  if ((std::abs(ggmother.pdgCode()) >= pdgCodeCharmMin && std::abs(ggmother.pdgCode()) < pdgCodeCharmMax) ||
                      (std::abs(ggmother.pdgCode()) >= pdgCodeBeautyMin && std::abs(ggmother.pdgCode()) < pdgCodeBeautyMax)) {
                    continue;
                  }
                  isEmbEta = true; // eta->gamma-> e
                }
              }
              if (std::abs(gmother.pdgCode()) == kPi0) {
                if (gmother.isPhysicalPrimary() || gmother.has_mothers()) {
                  auto const& ggmother = gmother.mothers_first_as<aod::McParticles>();
                  if ((std::abs(ggmother.pdgCode()) >= pdgCodeCharmMin && std::abs(ggmother.pdgCode()) < pdgCodeCharmMax) ||
                      (std::abs(ggmother.pdgCode()) >= pdgCodeBeautyMin && std::abs(ggmother.pdgCode()) < pdgCodeBeautyMax)) {
                    continue;
                  }
                  isEmbPi0 = true; // pi0-> gamma-> e
                }
                if (gmother.has_mothers()) {
                  auto const& ggmother = gmother.mothers_first_as<aod::McParticles>();
                  if (std::abs(ggmother.pdgCode()) == kEtaLocal) {
                    if (ggmother.isPhysicalPrimary() || ggmother.has_mothers()) {
                      auto const& gggmother = ggmother.mothers_first_as<aod::McParticles>();
                      if ((std::abs(gggmother.pdgCode()) >= pdgCodeCharmMin && std::abs(gggmother.pdgCode()) < pdgCodeCharmMax) ||
                          (std::abs(gggmother.pdgCode()) >= pdgCodeBeautyMin && std::abs(gggmother.pdgCode()) < pdgCodeBeautyMax)) {
                        continue;
                      }
                      isEmbEta = true; // eta->pi0->gamma-> e
                    }
                  }
                }
              }
            }
            if (isEmbPi0 || isEmbEta) {
              registry.fill(HIST("hMcgenNonHfeElectron"), particleMc.pt());
              isNonHfe = true;
              if (isEmbPi0) {

                registry.fill(HIST("hPi0eEmbTrkPt"), particleMc.pt());
              }
              if (isEmbEta) {
                registry.fill(HIST("hEtaeEmbTrkPt"), particleMc.pt());
              }
            }
          }
        }
        hfGenElectronSel(mcCollision.globalIndex(), particleMc.globalIndex(), particleMc.eta(), particleMc.phi(), particleMc.pt(), isNonHfe);
      }
    }
  }
  PROCESS_SWITCH(HfElectronSelectionWithTpcEmcal, processMcGen, "Process MC Gen mode", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfElectronSelectionWithTpcEmcal>(cfgc)};
}
