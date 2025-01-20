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

/// \file jetSpectraEseTask.cxx
/// \brief jet spectra analysis framework with ESE (19/08/2024)
///
/// \author Joachim C. K. B. Hansen, Lund University

#include <string>
#include <vector>
#include <map>

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "Common/DataModel/EseTable.h"
#include "Common/DataModel/Qvectors.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"
#include "Framework/StaticFor.h"

struct JetSpectraEseTask {
  ConfigurableAxis binJetPt{"binJetPt", {200, 0., 200.}, ""};
  ConfigurableAxis bindPhi{"bindPhi", {100, -o2::constants::math::PI - 1, o2::constants::math::PI + 1}, ""};
  ConfigurableAxis binESE{"binESE", {100, 0, 100}, ""};
  ConfigurableAxis binCos{"binCos", {100, -1.05, 1.05}, ""};
  ConfigurableAxis binOccupancy{"binOccupancy", {5000, 0, 25000}, ""};
  ConfigurableAxis binQVec{"binQVec", {100, -3, 3}, ""};

  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> jetR{"jetR", 0.2, "jet resolution parameter"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0, "vertex z cut"};
  Configurable<std::vector<float>> centRange{"centRange", {30, 50}, "centrality region of interest"};
  Configurable<double> leadingJetPtCut{"leadingJetPtCut", 5.0, "leading jet pT cut"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  Configurable<int> fColSwitch{"fColSwitch", 0, "collision switch"};

  Configurable<bool> cfgEvSelOccupancy{"cfgEvSelOccupancy", false, "Flag for occupancy cut"};
  Configurable<std::vector<int>> cfgCutOccupancy{"cfgCutOccupancy", {0, 500}, "Occupancy cut"};
  Configurable<std::vector<float>> cfgOccupancyPtCut{"cfgOccupancyPtCut", {0, 100}, "pT cut"};

  Configurable<int> cfgnTotalSystem{"cfgnTotalSystem", 7, "total qvector number // look in Qvector table for this number"};
  Configurable<int> cfgnCorrLevel{"cfgnCorrLevel", 3, "QVector step: 0 = no corr, 1 = rect, 2 = twist, 3 = full"};

  Configurable<std::string> cfgEPRefA{"cfgEPRefA", "FT0A", "EP reference A"};
  Configurable<std::string> cfgEPRefB{"cfgEPRefB", "TPCpos", "EP reference B"};
  Configurable<std::string> cfgEPRefC{"cfgEPRefC", "TPCneg", "EP reference C"};

  AxisSpec jetPtAxis = {binJetPt, "#it{p}_{T,jet}"};
  AxisSpec dPhiAxis = {bindPhi, "#Delta#phi"};
  AxisSpec eseAxis = {binESE, "#it{q}_{2}"};
  AxisSpec cosAxis = {binCos, ""};
  AxisSpec occAxis = {binOccupancy, "Occupancy"};
  AxisSpec qvecAxis = {binQVec, "Q-vector"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  int eventSelection = -1;
  int trackSelection = -1;

  enum class DetID { FT0C,
                     FT0A,
                     FT0M,
                     FV0A,
                     TPCpos,
                     TPCneg,
                     TPCall };

  void init(o2::framework::InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    LOGF(info, "jetSpectraEse::init()");

    switch (fColSwitch) {
      case 0:
        LOGF(info, "JetSpectraEseTask::init() - using data");
        registry.add("hEventCounter", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});
        registry.add("hJetPt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{jetPtAxis}}});
        registry.add("hJetPt_bkgsub", "jet pT background sub;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{jetPtAxis}}});
        registry.add("hJetEta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
        registry.add("hJetPhi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}});
        registry.add("hRho", ";#rho;entries", {HistType::kTH1F, {{100, 0, 200.}}});
        registry.add("hJetArea", ";area_{jet};entries", {HistType::kTH1F, {{100, 0, 10.}}});
        registry.add("hdPhi", "#Delta#phi;entries;", {HistType::kTH1F, {{dPhiAxis}}});
        registry.add("hCentJetPtdPhiq2", "", {HistType::kTHnSparseF, {{100, 0, 100}, {jetPtAxis}, {dPhiAxis}, {eseAxis}}});
        registry.add("hPsi2FT0C", ";Centrality; #Psi_{2}", {HistType::kTH2F, {{100, 0, 100}, {150, -2.5, 2.5}}});
        registry.addClone("hPsi2FT0C", "hPsi2FT0A");
        registry.addClone("hPsi2FT0C", "hPsi2FV0A");
        registry.addClone("hPsi2FT0C", "hPsi2TPCpos");
        registry.addClone("hPsi2FT0C", "hPsi2TPCneg");

        registry.add("hCosPsi2AmC", ";Centrality;cos(2(#Psi_{2}^{A}-#Psi_{2}^{B}));#it{q}_{2}", {HistType::kTH3F, {{100, 0, 100}, {cosAxis}, {eseAxis}}});
        registry.addClone("hCosPsi2AmC", "hCosPsi2AmB");
        registry.addClone("hCosPsi2AmC", "hCosPsi2BmC");

        registry.add("hQvecUncorV2", ";Centrality;Q_x;Q_y", {HistType::kTH3F, {{100, 0, 100}, {qvecAxis}, {qvecAxis}}});
        registry.addClone("hQvecUncorV2", "hQvecRectrV2");
        registry.addClone("hQvecUncorV2", "hQvecTwistV2");
        registry.addClone("hQvecUncorV2", "hQvecFinalV2");

        registry.addClone("hPsi2FT0C", "hEPUncorV2");
        registry.addClone("hPsi2FT0C", "hEPRectrV2");
        registry.addClone("hPsi2FT0C", "hEPTwistV2");
        break;
      case 1:
        LOGF(info, "JetSpectraEseTask::init() - using MC");
        registry.add("hMCEventCounter", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});
        registry.add("hPartJetPt", "particle level jet pT;#it{p}_{T,jet part} (GeV/#it{c});entries", {HistType::kTH1F, {{jetPtAxis}}});
        registry.add("hPartJetEta", "particle level jet #eta;#eta_{jet part};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
        registry.add("hPartJetPhi", "particle level jet #phi;#phi_{jet part};entries", {HistType::kTH1F, {{80, -1.0, 7.}}});
        registry.add("hPartJetPtMatch", "particle level jet pT;#it{p}_{T,jet part} (GeV/#it{c});entries", {HistType::kTH1F, {{jetPtAxis}}});
        registry.add("hPartJetEtaMatch", "particle level jet #eta;#eta_{jet part};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
        registry.add("hPartJetPhiMatch", "particle level jet #phi;#phi_{jet part};entries", {HistType::kTH1F, {{80, -1.0, 7.}}});
        registry.add("hDetectorJetPt", "detector level jet pT;#it{p}_{T,jet det} (GeV/#it{c});entries", {HistType::kTH1F, {{jetPtAxis}}});
        registry.add("hDetectorJetEta", "detector level jet #eta;#eta_{jet det};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
        registry.add("hDetectorJetPhi", "detector level jet #phi;#phi_{jet det};entries", {HistType::kTH1F, {{80, -1.0, 7.}}});
        registry.add("hMatchedJetsPtDelta", "#it{p}_{T,jet part}; det - part", {HistType::kTH2F, {{jetPtAxis}, {200, -20., 20.0}}});
        registry.add("hMatchedJetsEtaDelta", "#eta_{jet part}; det - part", {HistType::kTH2F, {{100, -1.0, 1.0}, {200, -20.0, 20.0}}});
        registry.add("hMatchedJetsPhiDelta", "#phi_{jet part}; det - part", {HistType::kTH2F, {{80, -1.0, 7.}, {200, -20.0, 20.}}});
        registry.add("hResponseMatrixMatch", "#it{p}_{T, jet det}; #it{p}_{T, jet part}", HistType::kTH2F, {jetPtAxis, jetPtAxis});
        break;
      case 2:
        LOGF(info, "JetSpectraEseTask::init() - using Occupancy processing");
        registry.add("hEventCounterOcc", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});
        registry.add("hTrackPt", "track pT;#it{p}_{T,track} (GeV/#it{c});entries", {HistType::kTHnSparseF, {{100, 0, 100}, {100, 0, 100}, {eseAxis}, {occAxis}}});
        registry.add("hTrackEta", "track #eta;#eta_{track};entries", {HistType::kTH3F, {{100, 0, 100}, {100, -1.0, 1.0}, {occAxis}}});
        registry.add("hTrackPhi", "track #phi;#phi_{track};entries", {HistType::kTH3F, {{100, 0, 100}, {80, -1.0, 7.}, {occAxis}}});
        registry.add("hOccupancy", "Occupancy;Occupancy;entries", {HistType::kTH1F, {{occAxis}}});
        registry.add("hPsiOccupancy", "Occupancy;#Psi_{2};entries", {HistType::kTH3F, {{100, 0, 100}, {150, -2.5, 2.5}, {occAxis}}});
        break;
    }
  }

  Filter jetCuts = aod::jet::pt > jetPtMin&& aod::jet::r == nround(jetR.node() * 100.0f) && nabs(aod::jet::eta) < 0.9f - jetR;
  Filter colFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter mcCollisionFilter = nabs(aod::jmccollision::posZ) < vertexZCut;

  void processESEDataCharged(soa::Join<aod::JetCollisions, aod::BkgChargedRhos, aod::Qvectors, aod::QPercentileFT0Cs>::iterator const& collision,
                             soa::Filtered<aod::ChargedJets> const& jets,
                             aod::JetTracks const& tracks)
  {
    float counter{0.5f};
    registry.fill(HIST("hEventCounter"), counter++);
    registry.fill(HIST("hEventCounter"), counter++);

    const auto vPsi2 = procEP<true>(collision);
    const auto qPerc = collision.qPERCFT0C();
    if (qPerc[0] < 0)
      return;
    registry.fill(HIST("hEventCounter"), counter++);

    if (!jetderiveddatautilities::selectCollision(collision, eventSelection))
      return;

    registry.fill(HIST("hEventCounter"), counter++);

    if (!isAcceptedLeadTrack(tracks))
      return;

    if (cfgEvSelOccupancy) {
      auto occupancy = collision.trackOccupancyInTimeRange();
      if (occupancy < cfgCutOccupancy->at(0) || occupancy > cfgCutOccupancy->at(1))
        registry.fill(HIST("hEventCounter"), counter++);
      return;
    }

    registry.fill(HIST("hEventCounter"), counter++);
    registry.fill(HIST("hRho"), collision.rho());
    for (auto const& jet : jets) {
      float jetpTbkgsub = jet.pt() - (collision.rho() * jet.area());
      registry.fill(HIST("hJetPt"), jet.pt());
      registry.fill(HIST("hJetPt_bkgsub"), jetpTbkgsub);
      registry.fill(HIST("hJetEta"), jet.eta());
      registry.fill(HIST("hJetPhi"), jet.phi());
      registry.fill(HIST("hJetArea"), jet.area());

      float dPhi = RecoDecay::constrainAngle(jet.phi() - vPsi2, -o2::constants::math::PI);
      registry.fill(HIST("hdPhi"), dPhi);
      registry.fill(HIST("hCentJetPtdPhiq2"), collision.centrality(), jetpTbkgsub, dPhi, qPerc[0]);
    }
    registry.fill(HIST("hEventCounter"), counter++);

    if (collision.centrality() < centRange->at(0) || collision.centrality() > centRange->at(1)) /* for counting */
      return;
    registry.fill(HIST("hEventCounter"), counter++);
  }
  PROCESS_SWITCH(JetSpectraEseTask, processESEDataCharged, "process ese collisions", true);

  void processESEOccupancy(soa::Join<aod::JetCollisions, aod::Qvectors, aod::QPercentileFT0Cs>::iterator const& collision,
                           soa::Join<aod::JetTracks, aod::JTrackPIs> const& tracks)
  {
    float count{0.5};
    registry.fill(HIST("hEventCounterOcc"), count++);
    const auto vPsi2 = procEP<false>(collision);
    const auto qPerc = collision.qPERCFT0C();

    auto occupancy = collision.trackOccupancyInTimeRange();
    registry.fill(HIST("hPsiOccupancy"), collision.centrality(), vPsi2, occupancy);
    registry.fill(HIST("hOccupancy"), occupancy);

    if (!jetderiveddatautilities::selectCollision(collision, eventSelection))
      return;
    registry.fill(HIST("hEventCounterOcc"), count++);

    for (auto const& track : tracks) {
      if (!jetderiveddatautilities::selectTrack(track, trackSelection))
        continue;

      registry.fill(HIST("hTrackPt"), collision.centrality(), track.pt(), qPerc[0], occupancy);
      if (track.pt() < cfgOccupancyPtCut->at(0) || track.pt() > cfgOccupancyPtCut->at(1))
        continue;
      registry.fill(HIST("hTrackEta"), collision.centrality(), track.eta(), occupancy);
      registry.fill(HIST("hTrackPhi"), collision.centrality(), track.phi(), occupancy);
    }
  }
  PROCESS_SWITCH(JetSpectraEseTask, processESEOccupancy, "process occupancy", false);

  void processMCParticleLevel(soa::Filtered<aod::ChargedMCParticleLevelJets>::iterator const& jet)
  {
    registry.fill(HIST("hPartJetPt"), jet.pt());
    registry.fill(HIST("hPartJetEta"), jet.eta());
    registry.fill(HIST("hPartJetPhi"), jet.phi());
  }
  PROCESS_SWITCH(JetSpectraEseTask, processMCParticleLevel, "jets on particle level MC", false);

  using JetMCPTable = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>;
  void processMCChargedMatched(soa::Filtered<aod::JetCollisionsMCD>::iterator const& collision,
                               soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>> const& mcdjets,
                               JetMCPTable const&,
                               aod::JetTracks const&,
                               aod::JetParticles const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection))
      return;

    float counter{0.5f};
    registry.fill(HIST("hMCEventCounter"), counter++);
    for (const auto& mcdjet : mcdjets) {

      registry.fill(HIST("hDetectorJetPt"), mcdjet.pt());
      registry.fill(HIST("hDetectorJetEta"), mcdjet.eta());
      registry.fill(HIST("hDetectorJetPhi"), mcdjet.phi());
      for (const auto& mcpjet : mcdjet.template matchedJetGeo_as<JetMCPTable>()) {

        registry.fill(HIST("hPartJetPtMatch"), mcpjet.pt());
        registry.fill(HIST("hPartJetEtaMatch"), mcpjet.eta());
        registry.fill(HIST("hPartJetPhiMatch"), mcpjet.phi());

        registry.fill(HIST("hMatchedJetsPtDelta"), mcpjet.pt(), mcdjet.pt() - mcpjet.pt());
        registry.fill(HIST("hMatchedJetsPhiDelta"), mcpjet.phi(), mcdjet.phi() - mcpjet.phi());
        registry.fill(HIST("hMatchedJetsEtaDelta"), mcpjet.eta(), mcdjet.eta() - mcpjet.eta());

        registry.fill(HIST("hResponseMatrixMatch"), mcdjet.pt(), mcpjet.pt());
      }
    }
    registry.fill(HIST("hMCEventCounter"), counter++);
  }
  PROCESS_SWITCH(JetSpectraEseTask, processMCChargedMatched, "jet matched mcp and mcd", false);

  template <typename T>
  bool isAcceptedLeadTrack(T const& tracks)
  {
    double leadingTrackpT{0.0};
    for (const auto& track : tracks) {
      if (track.pt() > leadingJetPtCut) {
        if (track.pt() > leadingTrackpT) {
          leadingTrackpT = track.pt();
        }
      }
    }
    if (leadingTrackpT == 0.0)
      return false;
    else
      return true;
  }
  template <bool Fill, typename EPCol>
  float procEP(EPCol const& vec)
  {
    constexpr std::array<float, 2> AmpCut{1e-8, 0.0};
    auto computeEP = [&AmpCut](std::vector<float> vec, auto det) { return vec[2] > AmpCut[det] ? 0.5 * std::atan2(vec[1], vec[0]) : 999.; };
    std::map<std::string, float> epMap;
    epMap["FT0A"] = computeEP(qVecNoESE<DetID::FT0A, Fill>(vec), 0);
    if constexpr (Fill) {
      epMap["FT0C"] = computeEP(qVecNoESE<DetID::FT0C, false>(vec), 0);
      epMap["FV0A"] = computeEP(qVecNoESE<DetID::FV0A, false>(vec), 0);
      epMap["TPCpos"] = computeEP(qVecNoESE<DetID::TPCpos, false>(vec), 1);
      epMap["TPCneg"] = computeEP(qVecNoESE<DetID::TPCneg, false>(vec), 1);
      fillEP(/*std::make_index_sequence<5>{},*/ vec, epMap);
      auto cosPsi = [](float psiX, float psiY) { return (static_cast<double>(psiX) == 999. || static_cast<double>(psiY) == 999.) ? 999. : std::cos(2.0 * (psiX - psiY)); };
      std::array<float, 3> epCorrContainer{};
      epCorrContainer[0] = cosPsi(epMap[cfgEPRefA], epMap[cfgEPRefC]);
      epCorrContainer[1] = cosPsi(epMap[cfgEPRefA], epMap[cfgEPRefB]);
      epCorrContainer[2] = cosPsi(epMap[cfgEPRefB], epMap[cfgEPRefC]);
      fillEPCos(/*std::make_index_sequence<3>{},*/ vec, epCorrContainer);
    }
    return epMap["FT0A"];
  }
  template </*std::size_t... Idx,*/ typename collision>
  void fillEPCos(/*const std::index_sequence<Idx...>&,*/ const collision& col, const std::array<float, 3>& Corr)
  {
    // static constexpr std::string CosList[] = {"hCosPsi2AmC", "hCosPsi2AmB", "hCosPsi2BmC"};
    // (registry.fill(HIST(CosList[Idx]), col.centrality(), Corr[Idx], col.qPERCFT0C()[0]), ...);
    registry.fill(HIST("hCosPsi2AmC"), col.centrality(), Corr[0], col.qPERCFT0C()[0]);
    registry.fill(HIST("hCosPsi2AmB"), col.centrality(), Corr[1], col.qPERCFT0C()[0]);
    registry.fill(HIST("hCosPsi2BmC"), col.centrality(), Corr[2], col.qPERCFT0C()[0]);
  }

  template </*std::size_t... Idx,*/ typename collision>
  void fillEP(/*const std::index_sequence<Idx...>&,*/ const collision& col, const std::map<std::string, float>& epMap)
  {
    // static constexpr std::string_view EpList[] = {"hPsi2FT0A", "hPsi2FV0A", "hPsi2FT0C", "hPsi2TPCpos", "hPsi2TPCneg"};
    // (registry.fill(HIST(EpList[Idx]), col.centrality(), epMap.at(std::string(RemovePrefix(EpList[Idx])))), ...);
    registry.fill(HIST("hPsi2FT0A"), col.centrality(), epMap.at("FT0A"));
    registry.fill(HIST("hPsi2FV0A"), col.centrality(), epMap.at("FV0A"));
    registry.fill(HIST("hPsi2FT0C"), col.centrality(), epMap.at("FT0C"));
    registry.fill(HIST("hPsi2TPCpos"), col.centrality(), epMap.at("TPCpos"));
    registry.fill(HIST("hPsi2TPCneg"), col.centrality(), epMap.at("TPCneg"));
  }
  constexpr std::string_view RemovePrefix(std::string_view str)
  {
    constexpr std::string_view Prefix = "hPsi2";
    return str.substr(Prefix.size());
  }

  constexpr int detIDN(const DetID id)
  {
    switch (id) {
      case DetID::FT0C:
        return 0;
      case DetID::FT0A:
        return 1;
      case DetID::FT0M:
        return 2;
      case DetID::FV0A:
        return 3;
      case DetID::TPCpos:
        return 4;
      case DetID::TPCneg:
        return 5;
      case DetID::TPCall:
        return 6;
    }
    return -1;
  }

  template <DetID id, bool fill, typename Col>
  std::vector<float> qVecNoESE(Col collision)
  {
    const int nmode = 2;
    int detId = detIDN(id);
    int detInd = detId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
    if constexpr (fill) {
      if (collision.qvecAmp()[detInd] > 1e-8) {
        registry.fill(HIST("hQvecUncorV2"), collision.centrality(), collision.qvecRe()[detInd], collision.qvecIm()[detInd]);
        registry.fill(HIST("hQvecRectrV2"), collision.centrality(), collision.qvecRe()[detInd + 1], collision.qvecIm()[detInd + 1]);
        registry.fill(HIST("hQvecTwistV2"), collision.centrality(), collision.qvecRe()[detInd + 2], collision.qvecIm()[detInd + 2]);
        registry.fill(HIST("hQvecFinalV2"), collision.centrality(), collision.qvecRe()[detInd + 3], collision.qvecIm()[detInd + 3]);
        registry.fill(HIST("hEPUncorV2"), collision.centrality(), 0.5 * std::atan2(collision.qvecIm()[detInd], collision.qvecRe()[detInd]));
        registry.fill(HIST("hEPRectrV2"), collision.centrality(), 0.5 * std::atan2(collision.qvecIm()[detInd + 1], collision.qvecRe()[detInd + 1]));
        registry.fill(HIST("hEPTwistV2"), collision.centrality(), 0.5 * std::atan2(collision.qvecIm()[detInd + 2], collision.qvecRe()[detInd + 2]));
      }
    }
    std::vector<float> qVec{};
    qVec.push_back(collision.qvecRe()[detInd + cfgnCorrLevel]);
    qVec.push_back(collision.qvecIm()[detInd + cfgnCorrLevel]);
    qVec.push_back(collision.qvecAmp()[detId]);
    return qVec;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetSpectraEseTask>(cfgc)}; }
