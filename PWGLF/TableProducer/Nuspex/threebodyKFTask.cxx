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
/// \brief Analysis task for KFVtx3BodyDatas (3body candidates reconstructed with KF and (optionally) strangeness tracking)
/// \author Carolina Reetz <c.reetz@cern.ch> --> in large parts copied from threebodyRecoTask.cxx
// ========================

#include <cmath>
#include <array>
#include <cstdlib>
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
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "CommonConstants/PhysicsConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullPi, aod::pidTPCFullDe>;
using MCLabeledTracksIU = soa::Join<FullTracksExtIU, aod::McTrackLabels>;

struct Candidate3body {
  int mcmotherId;
  int track0Id;
  int track1Id;
  int track2Id;
  TLorentzVector lcand;
  TLorentzVector lproton;
  TLorentzVector lpion;
  TLorentzVector lbachelor;
  // 0 - proton, 1 - pion, 2 - bachelor
  uint8_t dautpcNclusters[3];
  uint8_t dauitsclussize[3];
  uint8_t daudcaxytopv[3];
  uint8_t daudcatopv[3];
  float dautpcNsigma[3];
  bool isMatter;
  float invmass;
  float ct;
  float cosPA;
  float dcadaughters;
  float dcacandtopv;
  float bachelortofNsigma;
  TLorentzVector lgencand = {0, 0, 0, 0};
  float genct = -1;
  float genrapidity = -999;
  bool isSignal = false;
  bool isReco = false;
  int pdgCode = -1;
  bool SurvivedEventSelection = false;
};

struct threebodyKFTask {

  Produces<aod::Hyp3BodyCands> outputDataTable;
  Produces<aod::MCHyp3BodyCands> outputMCTable;
  std::vector<Candidate3body> Candidates3body;
  std::vector<unsigned int> filledMothers;
  std::vector<bool> isGoodCollision;

  // Configurables
  Configurable<int> bachelorPdgCode{"bachelorPdgCode", 1000010020, "pdgCode of bachelor daughter"};
  Configurable<int> motherPdgCode{"motherPdgCode", 1010010030, "pdgCode of mother"};

  // collision filter
  Filter collisionFilter = (aod::evsel::sel8 == true && nabs(aod::collision::posZ) < 10.f);

  HistogramRegistry registry{
    "registry",
    {
      {"hCentFT0C", "hCentFT0C", {HistType::kTH1F, {{100, 0.0f, 100.0f, "FT0C Centrality"}}}},
      {"hMassHypertriton", "hMassHypertriton", {HistType::kTH1F, {{80, 2.96f, 3.04f}}}},
      {"hMassAntiHypertriton", "hMassAntiHypertriton", {HistType::kTH1F, {{80, 2.96f, 3.04f}}}},

      {"hDalitz", "hDalitz", {HistType::kTH2F, {{60, 1.1, 1.4, "#it{M}^{2}(p#pi) ((GeV/c^{2})^{2})"}, {120, 7.85, 8.45, "#it{M}^{2}(#pi d) ((GeV/c^{2})^{2})"}}}},

      // for mcparticles information
      {"hLabelCounter", "hLabelCounter", {HistType::kTH1F, {{3, 0.0f, 3.0f}}}},
      {"hTrueHypertritonMCPt", "hTrueHypertritonMCPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
      {"hTrueAntiHypertritonMCPt", "hTrueAntiHypertritonMCPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
      {"hTrueHypertritonMCMass", "hTrueHypertritonMCMass", {HistType::kTH1F, {{40, 2.95f, 3.05f, "Inv. Mass (GeV/c^{2})"}}}},
      {"hTrueAntiHypertritonMCMass", "hTrueAntiHypertritonMCMass", {HistType::kTH1F, {{40, 2.95f, 3.05f, "Inv. Mass (GeV/c^{2})"}}}},
      {"hTrueHypertritonMCLifetime", "hTrueHypertritonMCLifetime", {HistType::kTH1F, {{50, 0.0f, 50.0f, "ct(cm)"}}}},
      {"hTrueAntiHypertritonMCLifetime", "hTrueAntiHypertritonMCLifetime", {HistType::kTH1F, {{50, 0.0f, 50.0f, "ct(cm)"}}}},

      {"hTrueHypertritonCounter", "hTrueHypertritonCounter", {HistType::kTH1F, {{12, 0.0f, 12.0f}}}},
      {"hGeneratedHypertritonCounter", "hGeneratedHypertritonCounter", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
      {"hctGeneratedHypertriton", "hctGeneratedHypertriton", {HistType::kTH1F, {{50, 0, 50, "ct(cm)"}}}},
      {"hEtaGeneratedHypertriton", "hEtaGeneratedHypertriton", {HistType::kTH1F, {{40, -2.0f, 2.0f}}}},
      {"hRapidityGeneratedHypertriton", "hRapidityGeneratedHypertriton", {HistType::kTH1F, {{40, -2.0f, 2.0f}}}},
      {"hPtGeneratedAntiHypertriton", "hPtGeneratedAntiHypertriton", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hctGeneratedAntiHypertriton", "hctGeneratedAntiHypertriton", {HistType::kTH1F, {{50, 0, 50, "ct(cm)"}}}},
      {"hEtaGeneratedAntiHypertriton", "hEtaGeneratedAntiHypertriton", {HistType::kTH1F, {{40, -2.0f, 2.0f}}}},
      {"hRapidityGeneratedAntiHypertriton", "hRapidityGeneratedAntiHypertriton", {HistType::kTH1F, {{40, -2.0f, 2.0f}}}},
    },
  };

  void init(InitContext const&)
  {
    registry.get<TH1>(HIST("hLabelCounter"))->GetXaxis()->SetBinLabel(1, "Total");
    registry.get<TH1>(HIST("hLabelCounter"))->GetXaxis()->SetBinLabel(2, "Have Same MotherTrack");
    registry.get<TH1>(HIST("hLabelCounter"))->GetXaxis()->SetBinLabel(3, "True H3L");
  }

  //------------------------------------------------------------------
  // Fill stats histograms

  struct {
    std::array<int32_t, kNCandSteps> candstats;
    std::array<int32_t, kNCandSteps> truecandstats;
  } statisticsRegistry;

  void resetHistos()
  {
    for (Int_t ii = 0; ii < kNCandSteps; ii++) {
      statisticsRegistry.candstats[ii] = 0;
      statisticsRegistry.truecandstats[ii] = 0;
    }
  }
  void FillCandCounter(int kn, bool istrue = false)
  {
    statisticsRegistry.candstats[kn]++;
    if (istrue) {
      statisticsRegistry.truecandstats[kn]++;
    }
  }
  void fillHistos()
  {
    for (Int_t ii = 0; ii < kNCandSteps; ii++) {
      registry.fill(HIST("hCandidatesCounter"), ii, statisticsRegistry.candstats[ii]);
      registry.fill(HIST("hTrueHypertritonCounter"), ii, statisticsRegistry.truecandstats[ii]);
    }
  }

  ConfigurableAxis dcaBinning{"dca-binning", {200, 0.0f, 1.0f}, ""};
  ConfigurableAxis ptBinning{"pt-binning", {200, 0.0f, 10.0f}, ""};

  void init(InitContext const&)
  {
    AxisSpec dcaAxis = {dcaBinning, "DCA (cm)"};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/c)"};
    AxisSpec massAxisHypertriton = {80, 2.96f, 3.04f, "Inv. Mass (GeV/c^{2})"};

    TString CandCounterbinLabel[kNCandSteps] = {"Total", "VtxCosPA", "TrackEta", "MomRapidity", "Lifetime", "VtxDcaDau", "d TOFPID", "TPCPID", "TPCNcls", "DauPt", "PionDcatoPV", "InvMass"};
    for (int i{0}; i < kNCandSteps; i++) {
      registry.get<TH1>(HIST("hCandidatesCounter"))->GetXaxis()->SetBinLabel(i + 1, CandCounterbinLabel[i]);
      registry.get<TH1>(HIST("hTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(i + 1, CandCounterbinLabel[i]);
    }

    registry.get<TH1>(HIST("hGeneratedHypertritonCounter"))->GetXaxis()->SetBinLabel(1, "Total");
    registry.get<TH1>(HIST("hGeneratedHypertritonCounter"))->GetXaxis()->SetBinLabel(2, "3-body decay");
  }

  //------------------------------------------------------------------
  Preslice<aod::Vtx3BodyDatas> perCollisionVtx3BodyDatas = o2::aod::vtx3body::collisionId;
  //------------------------------------------------------------------
  template <class TMCTrackTo, typename TMCParticle>
  bool is3bodyDecayed(TMCParticle const& particle)
  {
    if (std::abs(particle.pdgCode()) != motherPdgCode) {
      return false;
    }
    bool haveProton = false, havePion = false, haveBachelor = false;
    bool haveAntiProton = false, haveAntiPion = false, haveAntiBachelor = false;
    for (auto& mcparticleDaughter : particle.template daughters_as<TMCTrackTo>()) {
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
    if (haveProton && haveAntiPion && haveBachelor && particle.pdgCode() > 0) {
      return true;
    } else if (haveAntiProton && havePion && haveAntiBachelor && particle.pdgCode() < 0) {
      return true;
    }
    return false;
  }

  //------------------------------------------------------------------
  // Analysis process for a single candidate
  template <class TTrackClass, typename TCollisionTable, typename TCandTable>
  void CandidateAnalysis(TCollisionTable const& dCollision, TCandTable const& candData, bool& if_hasvtx, bool isTrueCand = false, int lLabel = -1, TLorentzVector lmother = {0, 0, 0, 0}, double MClifetime = -1)
  {

    FillCandCounter(kCandAll, isTrueCand);

    // auto track0 = vtx3bodyData.template track0_as<TTrackClass>();
    // auto track1 = vtx3bodyData.template track1_as<TTrackClass>();
    // auto track2 = vtx3bodyData.template track2_as<TTrackClass>();

    Candidate3body cand3body;
    // Hypertriton
    if ((track2.sign() > 0 && candData.mHypertriton() > h3LMassLowerlimit && candData.mHypertriton() < h3LMassUpperlimit)) {
      FillCandCounter(kCandInvMass, isTrueCand);

      registry.fill(HIST("hMassHypertriton"), candData.mHypertriton());
      registry.fill(HIST("hMassHypertritonTotal"), candData.mHypertriton());

      cand3body.isMatter = true;
      cand3body.lproton.SetXYZM(candData.pxtrack0(), candData.pytrack0(), candData.pztrack0(), o2::constants::physics::MassProton);
      cand3body.lpion.SetXYZM(candData.pxtrack1(), candData.pytrack1(), candData.pztrack1(), o2::constants::physics::MassPionCharged);

      if (candData.mHypertriton() > lowersignallimit && candData.mHypertriton() < uppersignallimit) {
        registry.fill(HIST("hDalitz"), RecoDecay::m2(array{array{candData.pxtrack0(), candData.pytrack0(), candData.pztrack0()}, array{candData.pxtrack2(), candData.pytrack2(), candData.pztrack2()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron}), RecoDecay::m2(array{array{candData.pxtrack0(), candData.pytrack0(), candData.pztrack0()}, array{candData.pxtrack1(), candData.pytrack1(), candData.pztrack1()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged}));
      }
    } else if ((track2.sign() < 0 && candData.mAntiHypertriton() > h3LMassLowerlimit && candData.mAntiHypertriton() < h3LMassUpperlimit)) {
      // AntiHypertriton
      FillCandCounter(kCandInvMass, isTrueCand);
      cand3body.isMatter = false;
      cand3body.lproton.SetXYZM(candData.pxtrack1(), candData.pytrack1(), candData.pztrack1(), o2::constants::physics::MassPionCharged);
      cand3body.lpion.SetXYZM(candData.pxtrack0(), candData.pytrack0(), candData.pztrack0(), o2::constants::physics::MassProton);

      registry.fill(HIST("hMassAntiHypertriton"), candData.mAntiHypertriton());
      registry.fill(HIST("hMassHypertritonTotal"), candData.mAntiHypertriton());
      if (candData.mAntiHypertriton() > lowersignallimit && candData.mAntiHypertriton() < uppersignallimit) {
        registry.fill(HIST("hDalitz"), RecoDecay::m2(array{array{candData.pxtrack1(), candData.pytrack1(), candData.pztrack1()}, array{candData.pxtrack2(), candData.pytrack2(), candData.pztrack2()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron}), RecoDecay::m2(array{array{candData.pxtrack1(), candData.pytrack1(), candData.pztrack1()}, array{candData.pxtrack0(), candData.pytrack0(), candData.pztrack0()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged}));
      }
    } else {
      return;
    }

    if_hasvtx = true;
    cand3body.mcmotherId = lLabel;
    cand3body.track0Id = candData.track0Id();
    cand3body.track1Id = candData.track1Id();
    cand3body.track2Id = candData.track2Id();
    cand3body.invmass = cand3body.isMatter ? candData.mHypertriton() : candData.mAntiHypertriton();
    cand3body.lbachelor.SetXYZM(candData.pxtrack2(), candData.pytrack2(), candData.pztrack2(), o2::constants::physics::MassDeuteron);
    cand3body.dautpcNclusters[0] = trackProton.tpcNClsFound();
    cand3body.dautpcNclusters[1] = trackPion.tpcNClsFound();
    cand3body.dautpcNclusters[2] = trackDeuteron.tpcNClsFound();
    cand3body.dauitsclussize[0] = trackProton.itsClusterSizes();
    cand3body.dauitsclussize[1] = trackPion.itsClusterSizes();
    cand3body.dauitsclussize[2] = trackDeuteron.itsClusterSizes();
    cand3body.dautpcNsigma[0] = trackProton.tpcNSigmaPr();
    cand3body.dautpcNsigma[1] = trackPion.tpcNSigmaPi();
    cand3body.dautpcNsigma[2] = trackDeuteron.tpcNSigmaDe();
    cand3body.daudcaxytopv[0] = cand3body.isMatter ? candData.dcaXYtrack0topv() : candData.dcaXYtrack1topv();
    cand3body.daudcaxytopv[1] = cand3body.isMatter ? candData.dcaXYtrack1topv() : candData.dcaXYtrack0topv();
    cand3body.daudcaxytopv[2] = candData.dcaXYtrack2topv();
    cand3body.daudcatopv[0] = cand3body.isMatter ? candData.dcatrack0topv() : candData.dcatrack1topv();
    cand3body.daudcatopv[1] = cand3body.isMatter ? candData.dcatrack1topv() : candData.dcatrack0topv();
    cand3body.daudcatopv[2] = candData.dcatrack2topv();
    cand3body.lcand.SetXYZM(candData.px(), candData.py(), candData.pz(), o2::constants::physics::MassHyperTriton);
    cand3body.ct = ct;
    cand3body.cosPA = cospa;
    cand3body.dcadaughters = candData.dcaVtxdaughters();
    cand3body.dcacandtopv = candData.dcavtxtopv(dCollision.posX(), dCollision.posY(), dCollision.posZ());
    cand3body.bachelortofNsigma = candData.tofNSigmaBachDe();
    if (isTrueCand) {
      cand3body.mcmotherId = lLabel;
      cand3body.lgencand = lmother;
      cand3body.genct = MClifetime;
      cand3body.genrapidity = lmother.Rapidity();
      cand3body.isSignal = true;
      cand3body.isReco = true;
      cand3body.pdgCode = cand3body.isMatter ? motherPdgCode : -motherPdgCode;
      cand3body.SurvivedEventSelection = true;
      filledMothers.push_back(lLabel);
    }

    Candidates3body.push_back(cand3body);

    registry.fill(HIST("hTPCPIDProton"), trackProton.tpcNSigmaPr());
    registry.fill(HIST("hTPCPIDPion"), trackPion.tpcNSigmaPi());
    registry.fill(HIST("hTPCPIDDeuteron"), trackDeuteron.tpcNSigmaDe());
    registry.fill(HIST("hProtonTPCBB"), trackProton.sign() * trackProton.p(), trackProton.tpcSignal());
    registry.fill(HIST("hPionTPCBB"), trackPion.sign() * trackPion.p(), trackPion.tpcSignal());
    registry.fill(HIST("hDeuteronTPCBB"), trackDeuteron.sign() * trackDeuteron.p(), trackDeuteron.tpcSignal());
    registry.fill(HIST("hProtonTPCVsPt"), trackProton.pt(), trackProton.tpcNSigmaPr());
    registry.fill(HIST("hPionTPCVsPt"), trackProton.pt(), trackPion.tpcNSigmaPi());
    registry.fill(HIST("hDeuteronTPCVsPt"), trackDeuteron.pt(), trackDeuteron.tpcNSigmaDe());
    registry.fill(HIST("hTOFPIDDeuteron"), candData.tofNSigmaBachDe());
  }

  //------------------------------------------------------------------
  // collect information for generated hypertriton (should be called after event selection)
  void GetGeneratedH3LInfo(aod::McParticles const& particlesMC)
  {
    for (auto& mcparticle : particlesMC) {
      if (mcparticle.pdgCode() != motherPdgCode && mcparticle.pdgCode() != -motherPdgCode) {
        continue;
      }
      registry.fill(HIST("hGeneratedHypertritonCounter"), 0.5);

      bool haveProton = false, havePionPlus = false, haveDeuteron = false;
      bool haveAntiProton = false, havePionMinus = false, haveAntiDeuteron = false;
      double MClifetime = -1;
      for (auto& mcparticleDaughter : mcparticle.template daughters_as<aod::McParticles>()) {
        if (mcparticleDaughter.pdgCode() == 2212)
          haveProton = true;
        if (mcparticleDaughter.pdgCode() == -2212)
          haveAntiProton = true;
        if (mcparticleDaughter.pdgCode() == 211)
          havePionPlus = true;
        if (mcparticleDaughter.pdgCode() == -211)
          havePionMinus = true;
        if (mcparticleDaughter.pdgCode() == bachelorPdgCode) {
          haveDeuteron = true;
          MClifetime = RecoDecay::sqrtSumOfSquares(mcparticleDaughter.vx() - mcparticle.vx(), mcparticleDaughter.vy() - mcparticle.vy(), mcparticleDaughter.vz() - mcparticle.vz()) * o2::constants::physics::MassHyperTriton / mcparticle.p();
        }
        if (mcparticleDaughter.pdgCode() == -bachelorPdgCode) {
          haveAntiDeuteron = true;
          MClifetime = RecoDecay::sqrtSumOfSquares(mcparticleDaughter.vx() - mcparticle.vx(), mcparticleDaughter.vy() - mcparticle.vy(), mcparticleDaughter.vz() - mcparticle.vz()) * o2::constants::physics::MassHyperTriton / mcparticle.p();
        }
      }
      if (haveProton && havePionMinus && haveDeuteron && mcparticle.pdgCode() == motherPdgCode) {
        registry.fill(HIST("hGeneratedHypertritonCounter"), 1.5);
        registry.fill(HIST("hctGeneratedHypertriton"), MClifetime);
        registry.fill(HIST("hEtaGeneratedHypertriton"), mcparticle.eta());
        registry.fill(HIST("hRapidityGeneratedHypertriton"), mcparticle.y());
      } else if (haveAntiProton && havePionPlus && haveAntiDeuteron && mcparticle.pdgCode() == -motherPdgCode) {
        registry.fill(HIST("hGeneratedHypertritonCounter"), 1.5);
        registry.fill(HIST("hPtGeneratedAntiHypertriton"), mcparticle.pt());
        registry.fill(HIST("hctGeneratedAntiHypertriton"), MClifetime);
        registry.fill(HIST("hEtaGeneratedAntiHypertriton"), mcparticle.eta());
        registry.fill(HIST("hRapidityGeneratedAntiHypertriton"), mcparticle.y());
      }
    }
  }

  //------------------------------------------------------------------
  // process real data analysis
  void processData(soa::Filtered<soa::Join<aod::Collisions, aod::CentFT0Cs>>::iterator const& collision, aod::KFVtx3BodyDatas const& vtx3bodydatas)
  {
    registry.fill(HIST("hCentFT0C"), collision.centFT0C());

    for (auto& vtx3bodydata : vtx3bodydatas) {
      auto isMatter = (vtx3bodydata.track2sign() > 0) ? true : false;

      // mass
      if (isMatter) {
        registry.fill(HIST("hMassHypertriton"), vtx3bodydata.mass());
      } else {
        registry.fill(HIST("hMassAntiHypertriton"), vtx3bodydata.mass());
      }

      // Dalitz plot
      auto m2prpi = RecoDecay::m2(array{array{vtx3bodydata.pxtrack0(), vtx3bodydata.pytrack0(), vtx3bodydata.pztrack0()}, array{vtx3bodydata.pxtrack1(), vtx3bodydata.pytrack1(), vtx3bodydata.pztrack1()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged});
      auto m2pide = RecoDecay::m2(array{array{vtx3bodydata.pxtrack1(), vtx3bodydata.pytrack1(), vtx3bodydata.pztrack1()}, array{vtx3bodydata.pxtrack2(), vtx3bodydata.pytrack2(), vtx3bodydata.pztrack2()}}, array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassDeuteron});
      if (std::abs(vtx3bodydata.mass() - o2::constants::physics::MassHyperTriton) <= 0.005) {
        registry.fill(HIST("hDalitz"), m2prpi, m2pide);
      }
    }
  }
  PROCESS_SWITCH(threebodyKFTask, processData, "Data analysis", true);

  //------------------------------------------------------------------
  // process mc analysis
  void processMC(soa::Filtered<soa::Join<aod::Collisions, o2::aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs>>::iterator const& collision, 
                 soa::Join<aod::KFVtx3BodyDatas, aod::McKFVtx3BodyLabels> const& vtx3bodydatas, 
                 aod::McParticles const& particlesMC, 
                 MCLabeledTracksIU const&, 
                 aod::McCollisions const& mcCollisions)
  {
    // Candidates3body.clear();
    // filledMothers.clear();
    // GetGeneratedH3LInfo(particlesMC);
    // isGoodCollision.resize(mcCollisions.size(), false);

    for (auto& vtx3bodydata : vtx3bodydatas) {
      registry.fill(HIST("hLabelCounter"), 0.5);

      int motherLabel = -1;
      double MClifetime = -1;
      TLorentzVector lmother;
      bool isTrueH3L = false;
      bool isTrueAntiH3L = false;

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
          if (MCvtx3body.pdgCode() == motherPdgCode && lMCTrack0.pdgCode() == 2212 && lMCTrack1.pdgCode() == -211 && lMCTrack2.pdgCode() == bachelorPdgCode) {
            isTrueH3L = true;
            motherLabel = MCvtx3body.mcParticleId();
            double hypertritonMCMass = RecoDecay::m(array{array{lMCTrack0.px(), lMCTrack0.py(), lMCTrack0.pz()}, array{lMCTrack1.px(), lMCTrack1.py(), lMCTrack1.pz()}, array{lMCTrack2.px(), lMCTrack2.py(), lMCTrack2.pz()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged, o2::constants::physics::MassDeuteron});
            MClifetime = RecoDecay::sqrtSumOfSquares(lMCTrack2.vx() - lMother2.vx(), lMCTrack2.vy() - lMother2.vy(), lMCTrack2.vz() - lMother2.vz()) * o2::constants::physics::MassHyperTriton / lMother2.p();
            registry.fill(HIST("hLabelCounter"), 2.5);
            registry.fill(HIST("hTrueHypertritonMCPt"), MCvtx3body.pt());
            registry.fill(HIST("hTrueHypertritonMCLifetime"), MClifetime);
            registry.fill(HIST("hTrueHypertritonMCMass"), hypertritonMCMass);
            lmother.SetXYZM(lMother0.px(), lMother0.py(), lMother0.pz(), o2::constants::physics::MassHyperTriton);
          } // end is H3L
          if (MCvtx3body.pdgCode() == -motherPdgCode && lMCTrack0.pdgCode() == 211 && lMCTrack1.pdgCode() == -2212 && lMCTrack2.pdgCode() == -bachelorPdgCode) {
            isTrueAntiH3L = true;
            motherLabel = MCvtx3body.mcParticleId();
            double antiHypertritonMCMass = RecoDecay::m(array{array{lMCTrack0.px(), lMCTrack0.py(), lMCTrack0.pz()}, array{lMCTrack1.px(), lMCTrack1.py(), lMCTrack1.pz()}, array{lMCTrack2.px(), lMCTrack2.py(), lMCTrack2.pz()}}, array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron});
            MClifetime = RecoDecay::sqrtSumOfSquares(lMCTrack2.vx() - lMother2.vx(), lMCTrack2.vy() - lMother2.vy(), lMCTrack2.vz() - lMother2.vz()) * o2::constants::physics::MassHyperTriton / lMother2.p();
            registry.fill(HIST("hLabelCounter"), 2.5);
            registry.fill(HIST("hTrueAntiHypertritonMCPt"), MCvtx3body.pt());
            registry.fill(HIST("hTrueAntiHypertritonMCLifetime"), MClifetime);
            registry.fill(HIST("hTrueAntiHypertritonMCMass"), antiHypertritonMCMass);
            lmother.SetXYZM(lMother0.px(), lMother0.py(), lMother0.pz(), o2::constants::physics::MassHyperTriton);
          } // end is Anti-H3L
        } // end has daughters
      } // end has matched MC particle

      

      /// TODO: CONTINUE HERE!! 
      /// BRIEF: fill MC table and move rest of CandidateAnalysis here
      /// BRIEF: find out for what GetGeneratedH3LInfois needed
      /// BRIEF: make function to creat mass and Dalitz histograms to also include it here in processMC
      /// BRIEF: understand and potentially add las part about signal candidates that were not reconstructed here





      CandidateAnalysis<MCLabeledTracksIU>(collision, vtx, if_hasvtx, isTrueCand, lLabel, lmother, MClifetime);
    }

    if (if_hasvtx)
      registry.fill(HIST("hEventCounter"), 3.5);
    fillHistos();
    resetHistos();

    for (auto& cand3body : Candidates3body) {
      outputMCTable(collision.centFT0C(),
                    cand3body.isMatter, cand3body.invmass, cand3body.lcand.P(), cand3body.lcand.Pt(), cand3body.ct,
                    cand3body.cosPA, cand3body.dcadaughters, cand3body.dcacandtopv,
                    cand3body.lproton.Pt(), cand3body.lproton.Eta(), cand3body.lproton.Phi(),
                    cand3body.lpion.Pt(), cand3body.lpion.Eta(), cand3body.lpion.Phi(),
                    cand3body.lbachelor.Pt(), cand3body.lbachelor.Eta(), cand3body.lbachelor.Phi(),
                    cand3body.dautpcNclusters[0], cand3body.dautpcNclusters[1], cand3body.dautpcNclusters[2],
                    cand3body.dauitsclussize[0], cand3body.dauitsclussize[1], cand3body.dauitsclussize[2],
                    cand3body.dautpcNsigma[0], cand3body.dautpcNsigma[1], cand3body.dautpcNsigma[2], cand3body.bachelortofNsigma,
                    cand3body.daudcaxytopv[0], cand3body.daudcaxytopv[1], cand3body.daudcaxytopv[2],
                    cand3body.daudcatopv[0], cand3body.daudcatopv[1], cand3body.daudcatopv[2],
                    cand3body.lgencand.P(), cand3body.lgencand.Pt(), cand3body.genct, cand3body.lgencand.Phi(), cand3body.lgencand.Eta(), cand3body.lgencand.Rapidity(),
                    cand3body.isSignal, cand3body.isReco, cand3body.pdgCode, cand3body.SurvivedEventSelection);
    }

    // now we fill only the signal candidates that were not reconstructed
    for (auto& mcparticle : particlesMC) {
      if (!is3bodyDecayed<aod::McParticles>(mcparticle)) {
        continue;
      }
      if (std::find(filledMothers.begin(), filledMothers.end(), mcparticle.globalIndex()) != std::end(filledMothers)) {
        continue;
      }
      bool isSurEvSelection = isGoodCollision[mcparticle.mcCollisionId()];
      std::array<float, 3> posSV{0.f};
      for (auto& mcDaughter : mcparticle.daughters_as<aod::McParticles>()) {
        if (std::abs(mcDaughter.pdgCode()) == bachelorPdgCode) {
          posSV = {mcDaughter.vx(), mcDaughter.vy(), mcDaughter.vz()};
        }
      }
      double MClifetime = RecoDecay::sqrtSumOfSquares(posSV[0] - mcparticle.vx(), posSV[1] - mcparticle.vy(), posSV[2] - mcparticle.vz()) * o2::constants::physics::MassHyperTriton / mcparticle.p();
      outputMCTable(-1,
                    -1, -1, -1, -1, -1,
                    -1, -1, -1,
                    -1, -1, -1,
                    -1, -1, -1,
                    -1, -1, -1,
                    -1, -1, -1,
                    -1, -1, -1,
                    -1, -1, -1, -1,
                    -1, -1, -1,
                    -1, -1, -1,
                    mcparticle.p(), mcparticle.pt(), MClifetime, mcparticle.phi(), mcparticle.eta(), mcparticle.y(),
                    true, false, mcparticle.pdgCode(), isSurEvSelection);
    }
  }
  PROCESS_SWITCH(threebodyKFTask, processMC, "MC analysis", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<threebodyKFTask>(cfgc),
  };
}
