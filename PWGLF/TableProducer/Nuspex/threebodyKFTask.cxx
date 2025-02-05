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
/// \brief Analysis task for KFVtx3BodyDatas (3body candidates reconstructed with KF)
/// \author Carolina Reetz <c.reetz@cern.ch> --> partly copied from threebodyRecoTask.cxx
// ========================

#include <cmath>
#include <array>
#include <cstdlib>
#include <vector>
#include <TLorentzVector.h>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/Vtx3BodyTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "CommonConstants/PhysicsConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using MCLabeledTracksIU = soa::Join<aod::TracksIU, aod::McTrackLabels>;

struct threebodyKFTask {

  Produces<aod::McKFVtx3BodyDatas> outputMCTable;
  std::vector<unsigned int> filledMothers;
  std::vector<bool> isGoodCollision;

  // Configurables
  Configurable<int> bachelorPdgCode{"bachelorPdgCode", 1000010020, "pdgCode of bachelor daughter"};
  Configurable<int> motherPdgCode{"motherPdgCode", 1010010030, "pdgCode of mother"};

  ConfigurableAxis m2PrPiBins{"m2PrPiBins", {60, 1., 3.3}, "Binning for m2(p,pi) axis"};
  ConfigurableAxis m2PiDeBins{"m2PiDeBins", {120, 4.0, 8.0}, "Binning for m2(pi,d) axis"};

  // collision filter and preslice
  Filter collisionFilter = (aod::evsel::sel8 == true && nabs(aod::collision::posZ) < 10.f);
  Preslice<soa::Join<aod::KFVtx3BodyDatas, aod::McKFVtx3BodyLabels>> perCollisionVtx3BodyDatas = o2::aod::vtx3body::collisionId;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    const AxisSpec axisM2PrPi{m2PrPiBins, "#it{m}^{2}(p,#pi^{-}) ((GeV/#it{c}^{2})^{2})"};
    const AxisSpec axisM2PiDe{m2PiDeBins, "#it{m}^{2}(#pi^{+},#bar{d}) ((GeV/#it{c}^{2})^{2})"};

    registry.add<TH1>("hCentFT0C", "hCentFT0C", HistType::kTH1F, {{100, 0.0f, 100.0f, "FT0C Centrality"}});

    // mass spectrum of reco candidates
    registry.add<TH1>("hMassHypertriton", "Mass hypertriton", HistType::kTH1F, {{80, 2.96f, 3.04f, "#it{m}(p,#pi^{-},d) (GeV/#it{c}^{2})"}});
    registry.add<TH1>("hMassAntiHypertriton", "Mass anti-hypertriton", HistType::kTH1F, {{80, 2.96f, 3.04f, "#it{m}(#bar{p},#pi^{+},#bar{d}) (GeV/#it{c}^{2})"}});

    // Dalitz diagrams of reco candidates
    registry.add<TH2>("hDalitzHypertriton", "Dalitz diagram", HistType::kTH2F, {axisM2PrPi, axisM2PiDe})->GetYaxis()->SetTitle("#it{m}^{2}(#pi^{-},d) ((GeV/#it{c}^{2})^{2})");
    registry.add<TH2>("hDalitzAntiHypertriton", "Dalitz diagram", HistType::kTH2F, {axisM2PrPi, axisM2PiDe})->GetXaxis()->SetTitle("#it{m}^{2}(#bar{p},#pi^{+}) ((GeV/#it{c}^{2})^{2})");

    // bachelor histgrams
    registry.add<TH1>("hAverageITSClusterSizeBachelor", "Average ITS cluster size bachelor track", HistType::kTH1F, {{15, 0.5, 15.5, "#langle ITS cluster size #rangle"}});
    registry.add<TH1>("hdEdxBachelor", "TPC dE/dx bachelor track", HistType::kTH1F, {{200, 0.0f, 200.0f, "Average ITS cluster size"}});
    registry.add<TH1>("hPIDTrackingBachelor", "Tracking PID bachelor track", HistType::kTH1F, {{20, 0.5, 20.5, "Tracking PID identifier"}});

    // for gen information of reco candidates
    auto LabelHist = registry.add<TH1>("hLabelCounter", "Reco MC candidate counter", HistType::kTH1F, {{3, 0.0f, 3.0f}});
    LabelHist->GetXaxis()->SetBinLabel(1, "Total");
    LabelHist->GetXaxis()->SetBinLabel(2, "Have Same MotherTrack");
    LabelHist->GetXaxis()->SetBinLabel(3, "True H3L/Anti-H3L");
    registry.add<TH1>("hTrueHypertritonMCPt", "pT gen. of reco. H3L", HistType::kTH1F, {{100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"}});
    registry.add<TH1>("hTrueAntiHypertritonMCPt", "pT gen. of reco. Anti-H3L", HistType::kTH1F, {{100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"}});
    registry.add<TH1>("hTrueHypertritonMCMass", "mass gen. of reco. H3L", HistType::kTH1F, {{40, 2.96f, 3.04f, "#it{m}(p,#pi^{-},d) (GeV/#it{c}^{2})"}});
    registry.add<TH1>("hTrueAntiHypertritonMCMass", "mass gen. of reco. Anti-H3L", HistType::kTH1F, {{40, 2.96f, 3.04f, "#it{m}(#bar{p},#pi^{+},#bar{d}) (GeV/#it{c}^{2})"}});
    registry.add<TH1>("hTrueHypertritonMCCTau", "#it{c}#tau gen. of reco. H3L", HistType::kTH1F, {{50, 0.0f, 50.0f, "#it{c}#tau(cm)"}});
    registry.add<TH1>("hTrueAntiHypertritonMCCTau", "#it{c}#tau gen. of reco. Anti-H3L", HistType::kTH1F, {{50, 0.0f, 50.0f, "#it{c}#tau(cm)"}});
    registry.add<TH1>("hTrueHypertritonMCMassPrPi", "inv. mass gen. of reco. V0 pair (H3L)", HistType::kTH1F, {{100, 0.0f, 6.0f, "#it{m}(p,#pi^{-}) (GeV/#it{c}^{2})"}});
    registry.add<TH1>("hTrueAntiHypertritonMCMassPrPi", "inv. mass gen. of reco. V0 pair (Anti-H3L)", HistType::kTH1F, {{100, 0.0f, 6.0f, "#it{m}(#bar{p},#pi^{+}) (GeV/#it{c}^{2})"}});

    // for gen information of non reco candidates
    registry.add<TH1>("hTrueHypertritonMCMassPrPi_nonReco", "inv. mass gen. of non-reco. V0 pair (H3L)", HistType::kTH1F, {{100, 0.0f, 6.0f, "#it{m}(p,#pi^{-}) (GeV/#it{c}^{2})"}});
    registry.add<TH1>("hTrueAntiHypertritonMCMassPrPi_nonReco", "inv. mass gen. of non-reco. V0 pair (Anti-H3L)", HistType::kTH1F, {{100, 0.0f, 6.0f, "#it{m}(#bar{p},#pi^{+}) (GeV/#it{c}^{2})"}});
    registry.add<TH1>("hTrueHypertritonMCPtProton_nonReco", "Proton #it{p}_{T} gen. of non-reco. H3L", HistType::kTH1F, {{100, 0.0f, 6.0f, "#it{p}_{T}(p) (GeV/#it{c})"}});
    registry.add<TH1>("hTrueHypertritonMCPtPion_nonReco", "Pion #it{p}_{T} gen. of non-reco. H3L", HistType::kTH1F, {{100, 0.0f, 6.0f, "#it{p}_{T}(#pi) (GeV/#it{c})"}});
    registry.add<TH1>("hTrueAntiHypertritonMCPtProton_nonReco", "Proton #it{p}_{T} gen. of non-reco. Anti-H3L", HistType::kTH1F, {{100, 0.0f, 6.0f, "#it{p}_{T}(p) (GeV/#it{c})"}});
    registry.add<TH1>("hTrueAntiHypertritonMCPtPion_nonReco", "Pion #it{p}_{T} gen. of non-reco. Anti- H3L", HistType::kTH1F, {{100, 0.0f, 6.0f, "#it{p}_{T}(#pi) (GeV/#it{c})"}});
  }

  template <typename TCand>
  void fillQAPlots(TCand const& vtx3body)
  {
    // Mass plot
    if (vtx3body.track2sign() > 0) { // hypertriton
      registry.fill(HIST("hMassHypertriton"), vtx3body.mass());
    } else if (vtx3body.track2sign() < 0) { // anti-hypertriton
      registry.fill(HIST("hMassAntiHypertriton"), vtx3body.mass());
    }

    // Dalitz plot
    auto m2prpi = RecoDecay::m2(array{array{vtx3body.pxtrack0(), vtx3body.pytrack0(), vtx3body.pztrack0()}, array{vtx3body.pxtrack1(), vtx3body.pytrack1(), vtx3body.pztrack1()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged});
    auto m2pide = RecoDecay::m2(array{array{vtx3body.pxtrack1(), vtx3body.pytrack1(), vtx3body.pztrack1()}, array{vtx3body.pxtrack2(), vtx3body.pytrack2(), vtx3body.pztrack2()}}, array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassDeuteron});
    if (std::abs(vtx3body.mass() - o2::constants::physics::MassHyperTriton) <= 0.005) {
      if (vtx3body.track2sign() > 0) { // hypertriton
        registry.fill(HIST("hDalitzHypertriton"), m2prpi, m2pide);
      } else if (vtx3body.track2sign() < 0) { // anti-hypertriton
        registry.fill(HIST("hDalitzAntiHypertriton"), m2prpi, m2pide);
      }
    }

    // ITS cluster sizes
    registry.fill(HIST("hAverageITSClusterSizeBachelor"), vtx3body.itsclussizedeuteron());
    registry.fill(HIST("hdEdxBachelor"), vtx3body.tpcdedxdeuteron());
    registry.fill(HIST("hPIDTrackingBachelor"), vtx3body.pidtrackingdeuteron());
  }

  //------------------------------------------------------------------
  // process real data analysis
  void processData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>::iterator const& collision,
                   aod::KFVtx3BodyDatas const& vtx3bodydatas)
  {
    registry.fill(HIST("hCentFT0C"), collision.centFT0C());

    for (auto& vtx3bodydata : vtx3bodydatas) {
      // QA histograms
      fillQAPlots(vtx3bodydata);
    }
  }
  PROCESS_SWITCH(threebodyKFTask, processData, "Data analysis", true);

  //------------------------------------------------------------------
  // process mc analysis
  void processMC(soa::Join<aod::Collisions, o2::aod::McCollisionLabels, aod::EvSels> const& collisions,
                 soa::Join<aod::KFVtx3BodyDatas, aod::McKFVtx3BodyLabels> const& vtx3bodydatas,
                 aod::McParticles const& particlesMC,
                 MCLabeledTracksIU const&,
                 aod::McCollisions const& mcCollisions)
  {
    filledMothers.clear();
    isGoodCollision.resize(mcCollisions.size(), false);

    // loop over collisions
    for (const auto& collision : collisions) {
      // event selection
      if (!collision.sel8() || abs(collision.posZ()) > 10.f) {
        continue;
      }
      // reco collision survived event selection filter --> fill value for MC collision if collision is "true" MC collision
      if (collision.mcCollisionId() >= 0) {
        isGoodCollision[collision.mcCollisionId()] = true;
      }

      // fill MC table with reco MC candidate information and gen information if matched to MC particle
      auto Decay3BodyTable_thisCollision = vtx3bodydatas.sliceBy(perCollisionVtx3BodyDatas, collision.globalIndex());
      for (auto& vtx3bodydata : Decay3BodyTable_thisCollision) {
        registry.fill(HIST("hLabelCounter"), 0.5);

        // fill QA histograms for all reco candidates
        fillQAPlots(vtx3bodydata);

        double MClifetime = -1.;
        bool isTrueH3L = false;
        bool isTrueAntiH3L = false;
        float genPhi = -1.;
        float genEta = -1.;
        float genRap = -1.;
        float genP = -1.;
        float genPt = -1.;
        std::array<float, 3> genDecVtx{-1.f};
        int vtx3bodyPDGcode = -1;
        double MCmassPrPi = -1.;
        float genPosP = -1.;
        float genPosPt = -1.;
        float genNegP = -1.;
        float genNegPt = -1.;
        float genBachP = -1.;
        float genBachPt = -1.;

        auto track0 = vtx3bodydata.track0_as<MCLabeledTracksIU>();
        auto track1 = vtx3bodydata.track1_as<MCLabeledTracksIU>();
        auto track2 = vtx3bodydata.track2_as<MCLabeledTracksIU>();

        if (vtx3bodydata.has_mcParticle() && vtx3bodydata.mcParticleId() > -1 && vtx3bodydata.mcParticleId() <= particlesMC.size()) { // mother to daughter association already checked in decay3bodybuilder
          auto MCvtx3body = vtx3bodydata.mcParticle();
          registry.fill(HIST("hLabelCounter"), 1.5);
          if (MCvtx3body.has_daughters()) {
            auto lMCTrack0 = track0.mcParticle_as<aod::McParticles>();
            auto lMCTrack1 = track1.mcParticle_as<aod::McParticles>();
            auto lMCTrack2 = track2.mcParticle_as<aod::McParticles>();
            // check PDG codes
            if ((MCvtx3body.pdgCode() == motherPdgCode && lMCTrack0.pdgCode() == 2212 && lMCTrack1.pdgCode() == -211 && lMCTrack2.pdgCode() == bachelorPdgCode) ||
                (MCvtx3body.pdgCode() == -motherPdgCode && lMCTrack0.pdgCode() == 211 && lMCTrack1.pdgCode() == -2212 && lMCTrack2.pdgCode() == -bachelorPdgCode)) {
              vtx3bodyPDGcode = MCvtx3body.pdgCode();
              genDecVtx = {lMCTrack0.vx(), lMCTrack0.vy(), lMCTrack0.vz()};
              MClifetime = RecoDecay::sqrtSumOfSquares(lMCTrack2.vx() - MCvtx3body.vx(), lMCTrack2.vy() - MCvtx3body.vy(), lMCTrack2.vz() - MCvtx3body.vz()) * o2::constants::physics::MassHyperTriton / MCvtx3body.p();
              genPhi = MCvtx3body.phi();
              genEta = MCvtx3body.eta();
              genRap = MCvtx3body.y();
              genP = MCvtx3body.p();
              genPt = MCvtx3body.pt();
              genPosP = lMCTrack0.p();
              genPosPt = lMCTrack0.pt();
              genNegP = lMCTrack1.p();
              genNegPt = lMCTrack1.pt();
              genBachP = lMCTrack2.p();
              genBachPt = lMCTrack2.pt();
              filledMothers.push_back(MCvtx3body.globalIndex());
            } // end is H3L or Anti-H3L
            if (MCvtx3body.pdgCode() == motherPdgCode && lMCTrack0.pdgCode() == 2212 && lMCTrack1.pdgCode() == -211 && lMCTrack2.pdgCode() == bachelorPdgCode) {
              isTrueH3L = true;
              double hypertritonMCMass = RecoDecay::m(array{array{lMCTrack0.px(), lMCTrack0.py(), lMCTrack0.pz()}, array{lMCTrack1.px(), lMCTrack1.py(), lMCTrack1.pz()}, array{lMCTrack2.px(), lMCTrack2.py(), lMCTrack2.pz()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged, o2::constants::physics::MassDeuteron});
              MCmassPrPi = RecoDecay::m(array{array{lMCTrack0.px(), lMCTrack0.py(), lMCTrack0.pz()}, array{lMCTrack1.px(), lMCTrack1.py(), lMCTrack1.pz()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged});
              registry.fill(HIST("hLabelCounter"), 2.5);
              registry.fill(HIST("hTrueHypertritonMCPt"), MCvtx3body.pt());
              registry.fill(HIST("hTrueHypertritonMCCTau"), MClifetime);
              registry.fill(HIST("hTrueHypertritonMCMass"), hypertritonMCMass);
              registry.fill(HIST("hTrueHypertritonMCMassPrPi"), MCmassPrPi);
            } // end is H3L
            if (MCvtx3body.pdgCode() == -motherPdgCode && lMCTrack0.pdgCode() == 211 && lMCTrack1.pdgCode() == -2212 && lMCTrack2.pdgCode() == -bachelorPdgCode) {
              isTrueAntiH3L = true;
              double antiHypertritonMCMass = RecoDecay::m(array{array{lMCTrack0.px(), lMCTrack0.py(), lMCTrack0.pz()}, array{lMCTrack1.px(), lMCTrack1.py(), lMCTrack1.pz()}, array{lMCTrack2.px(), lMCTrack2.py(), lMCTrack2.pz()}}, array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron});
              MCmassPrPi = RecoDecay::m(array{array{lMCTrack0.px(), lMCTrack0.py(), lMCTrack0.pz()}, array{lMCTrack1.px(), lMCTrack1.py(), lMCTrack1.pz()}}, array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton});
              registry.fill(HIST("hLabelCounter"), 2.5);
              registry.fill(HIST("hTrueAntiHypertritonMCPt"), MCvtx3body.pt());
              registry.fill(HIST("hTrueAntiHypertritonMCCTau"), MClifetime);
              registry.fill(HIST("hTrueAntiHypertritonMCMass"), antiHypertritonMCMass);
              registry.fill(HIST("hTrueAntiHypertritonMCMassPrPi"), MCmassPrPi);
            } // end is Anti-H3L
          } // end has daughters
        } // end has matched MC particle

        outputMCTable( // filled for each reconstructed candidate (in KFVtx3BodyDatas)
          vtx3bodydata.mass(),
          vtx3bodydata.x(), vtx3bodydata.y(), vtx3bodydata.z(),
          vtx3bodydata.xerr(), vtx3bodydata.yerr(), vtx3bodydata.zerr(),
          vtx3bodydata.px(), vtx3bodydata.py(), vtx3bodydata.pz(), vtx3bodydata.pt(),
          vtx3bodydata.pxerr(), vtx3bodydata.pyerr(), vtx3bodydata.pzerr(), vtx3bodydata.pterr(),
          vtx3bodydata.sign(),
          vtx3bodydata.dcavtxtopvkf(), vtx3bodydata.dcaxyvtxtopvkf(),
          vtx3bodydata.vtxcospakf(), vtx3bodydata.vtxcosxypakf(),
          vtx3bodydata.vtxcospakftopo(), vtx3bodydata.vtxcosxypakftopo(),
          vtx3bodydata.decaylkf(), vtx3bodydata.decaylxykf(), vtx3bodydata.decayldeltal(),
          vtx3bodydata.chi2geondf(), vtx3bodydata.chi2topondf(),
          vtx3bodydata.ctaukftopo(),
          vtx3bodydata.trackedclsize(),
          vtx3bodydata.massv0(), vtx3bodydata.chi2massv0(),
          vtx3bodydata.cospav0(),
          vtx3bodydata.pxtrack0(), vtx3bodydata.pytrack0(), vtx3bodydata.pztrack0(),                                  // proton
          vtx3bodydata.pxtrack1(), vtx3bodydata.pytrack1(), vtx3bodydata.pztrack1(),                                  // pion
          vtx3bodydata.pxtrack2(), vtx3bodydata.pytrack2(), vtx3bodydata.pztrack2(),                                  // deuteron
          vtx3bodydata.tpcinnerparamtrack0(), vtx3bodydata.tpcinnerparamtrack1(), vtx3bodydata.tpcinnerparamtrack2(), // proton, pion, deuteron
          vtx3bodydata.dcatrack0topvkf(), vtx3bodydata.dcatrack1topvkf(), vtx3bodydata.dcatrack2topvkf(),             // proton, pion, deuteron
          vtx3bodydata.dcaxytrack0topvkf(), vtx3bodydata.dcaxytrack1topvkf(), vtx3bodydata.dcaxytrack2topvkf(),       // proton, pion, deuteron
          vtx3bodydata.dcaxytrack0tosvkf(), vtx3bodydata.dcaxytrack1tosvkf(), vtx3bodydata.dcaxytrack2tosvkf(),       // proton, pion, deuteron
          vtx3bodydata.dcatrack0totrack1kf(), vtx3bodydata.dcatrack0totrack2kf(), vtx3bodydata.dcatrack1totrack2kf(),
          vtx3bodydata.dcavtxdaughterskf(),
          vtx3bodydata.dcaxytrackpostopv(), vtx3bodydata.dcaxytracknegtopv(), vtx3bodydata.dcaxytrackbachtopv(),
          vtx3bodydata.dcatrackpostopv(), vtx3bodydata.dcatracknegtopv(), vtx3bodydata.dcatrackbachtopv(),
          vtx3bodydata.track0sign(), vtx3bodydata.track1sign(), vtx3bodydata.track2sign(), // proton, pion, deuteron
          vtx3bodydata.tpcnsigmaproton(), vtx3bodydata.tpcnsigmapion(), vtx3bodydata.tpcnsigmadeuteron(), vtx3bodydata.tpcnsigmapionbach(),
          vtx3bodydata.tpcdedxproton(), vtx3bodydata.tpcdedxpion(), vtx3bodydata.tpcdedxdeuteron(),
          vtx3bodydata.tofnsigmadeuteron(),
          vtx3bodydata.itsclussizedeuteron(),
          vtx3bodydata.pidtrackingdeuteron(),
          // MC info (-1 if not matched to MC particle)
          genP,
          genPt,
          genDecVtx[0], genDecVtx[1], genDecVtx[2],
          MClifetime,
          genPhi,
          genEta,
          genRap,
          genPosP, genPosPt, genNegP, genNegPt, genBachP, genBachPt,
          isTrueH3L, isTrueAntiH3L,
          vtx3bodyPDGcode,
          true,  // is reconstructed
          true); // reco event passed event selection
      } // end vtx3bodydatas loop
    } // end collision loop

    // generated MC particle analysis
    // fill MC table with gen information for all generated but not reconstructed particles
    for (auto& mcparticle : particlesMC) {

      double genMCmassPrPi = -1.;
      bool isTrueGenH3L = false;
      bool isTrueGenAntiH3L = false;
      float genPBach = -1.;
      float genPtBach = -1.;
      float genPPos = -1.;
      float genPtPos = -1.;
      float genPNeg = -1.;
      float genPtNeg = -1.;

      // check if mcparticle was reconstructed and already filled in the table
      if (std::find(filledMothers.begin(), filledMothers.end(), mcparticle.globalIndex()) != std::end(filledMothers)) {
        continue;
      }

      // set flag if corresponding reco collision survived event selection
      bool survEvSel = isGoodCollision[mcparticle.mcCollisionId()];

      // check if MC particle is hypertriton with 3-body decay
      if (std::abs(mcparticle.pdgCode()) != motherPdgCode) {
        continue;
      }
      bool haveProton = false, havePion = false, haveBachelor = false;
      bool haveAntiProton = false, haveAntiPion = false, haveAntiBachelor = false;
      for (auto& mcparticleDaughter : mcparticle.template daughters_as<aod::McParticles>()) {
        if (mcparticleDaughter.pdgCode() == 2212)
          haveProton = true;
        if (mcparticleDaughter.pdgCode() == -2212)
          haveAntiProton = true;
        if (mcparticleDaughter.pdgCode() == 211)
          havePion = true;
        if (mcparticleDaughter.pdgCode() == -211)
          haveAntiPion = true;
        if (mcparticleDaughter.pdgCode() == bachelorPdgCode)
          haveBachelor = true;
        if (mcparticleDaughter.pdgCode() == -bachelorPdgCode)
          haveAntiBachelor = true;
      }

      // check if particle or anti-particle
      if (haveProton && haveAntiPion && haveBachelor && mcparticle.pdgCode() > 0) {
        isTrueGenH3L = true;
        // get proton and pion daughter
        std::array<float, 3> protonMom{0.f};
        std::array<float, 3> piMinusMom{0.f};
        for (auto& mcparticleDaughter : mcparticle.template daughters_as<aod::McParticles>()) {
          if (mcparticleDaughter.pdgCode() == 2212) {
            protonMom = {mcparticleDaughter.px(), mcparticleDaughter.py(), mcparticleDaughter.pz()};
            genPPos = mcparticleDaughter.p();
            genPtPos = mcparticleDaughter.pt();
          } else if (mcparticleDaughter.pdgCode() == -211) {
            piMinusMom = {mcparticleDaughter.px(), mcparticleDaughter.py(), mcparticleDaughter.pz()};
            genPNeg = mcparticleDaughter.p();
            genPtNeg = mcparticleDaughter.pt();
          }
        }
        genMCmassPrPi = RecoDecay::m(array{protonMom, piMinusMom}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged});
        registry.fill(HIST("hTrueHypertritonMCMassPrPi_nonReco"), genMCmassPrPi);
        registry.fill(HIST("hTrueHypertritonMCPtProton_nonReco"), RecoDecay::sqrtSumOfSquares(protonMom[0], protonMom[1]));
        registry.fill(HIST("hTrueHypertritonMCPtPion_nonReco"), RecoDecay::sqrtSumOfSquares(piMinusMom[0], piMinusMom[1]));
      } else if (haveAntiProton && havePion && haveAntiBachelor && mcparticle.pdgCode() < 0) {
        isTrueGenAntiH3L = true;
        // get anti-proton and pion daughter
        std::array<float, 3> antiProtonMom{0.f};
        std::array<float, 3> piPlusMom{0.f};
        for (auto& mcparticleDaughter : mcparticle.template daughters_as<aod::McParticles>()) {
          if (mcparticleDaughter.pdgCode() == -2212) {
            antiProtonMom = {mcparticleDaughter.px(), mcparticleDaughter.py(), mcparticleDaughter.pz()};
            genPNeg = mcparticleDaughter.p();
            genPtNeg = mcparticleDaughter.pt();
          } else if (mcparticleDaughter.pdgCode() == 211) {
            piPlusMom = {mcparticleDaughter.px(), mcparticleDaughter.py(), mcparticleDaughter.pz()};
            genPPos = mcparticleDaughter.p();
            genPtPos = mcparticleDaughter.pt();
          }
        }
        genMCmassPrPi = RecoDecay::m(array{antiProtonMom, piPlusMom}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged});
        registry.fill(HIST("hTrueAntiHypertritonMCMassPrPi_nonReco"), genMCmassPrPi);
        registry.fill(HIST("hTrueAntiHypertritonMCPtProton_nonReco"), RecoDecay::sqrtSumOfSquares(antiProtonMom[0], antiProtonMom[1]));
        registry.fill(HIST("hTrueAntiHypertritonMCPtPion_nonReco"), RecoDecay::sqrtSumOfSquares(piPlusMom[0], piPlusMom[1]));
      } else {
        continue; // stop if particle is no true H3L or Anti-H3L
      }

      // get gen decay vertex and calculate ctau
      std::array<float, 3> genDecayVtx{0.f};
      for (auto& mcDaughter : mcparticle.daughters_as<aod::McParticles>()) {
        if (std::abs(mcDaughter.pdgCode()) == bachelorPdgCode) {
          genDecayVtx = {mcDaughter.vx(), mcDaughter.vy(), mcDaughter.vz()};
          genPBach = mcDaughter.p();
          genPtBach = mcDaughter.pt();
        }
      }
      double genMClifetime = RecoDecay::sqrtSumOfSquares(genDecayVtx[0] - mcparticle.vx(), genDecayVtx[1] - mcparticle.vy(), genDecayVtx[2] - mcparticle.vz()) * o2::constants::physics::MassHyperTriton / mcparticle.p();
      outputMCTable( // reco information (-1)
        -1,
        -1, -1, -1,
        -1, -1, -1,
        -1, -1, -1, -1,
        -1, -1, -1, -1,
        -1,
        -1, -1,
        -1, -1,
        -1, -1,
        -1, -1, -1,
        -1, -1,
        -1,
        -1,
        -1, -1,
        -1,
        -1, -1, -1,
        -1, -1, -1,
        -1, -1, -1,
        -1, -1, -1,
        -1, -1, -1,
        -1, -1, -1,
        -1, -1, -1,
        -1, -1, -1,
        -1,
        -1, -1, -1,
        -1, -1, -1,
        -1, -1, -1,
        -1, -1, -1, -1,
        -1, -1, -1,
        -1,
        -1,
        -1,
        // gen information
        mcparticle.p(),
        mcparticle.pt(),
        genDecayVtx[0], genDecayVtx[1], genDecayVtx[2],
        genMClifetime,
        mcparticle.phi(),
        mcparticle.eta(),
        mcparticle.y(),
        genPPos, genPtPos, genPNeg, genPtNeg, genPBach, genPtBach,
        isTrueGenH3L, isTrueGenAntiH3L,
        mcparticle.pdgCode(),
        false, // is reconstructed
        survEvSel);
    } // end mcparticles loop
  }
  PROCESS_SWITCH(threebodyKFTask, processMC, "MC analysis", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<threebodyKFTask>(cfgc),
  };
}
