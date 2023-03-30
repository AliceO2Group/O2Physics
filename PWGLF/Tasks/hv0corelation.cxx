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
/// \brief This task is to do the reconstruction of V0s (cascade may be involved later).
 //  The yield will be got by corelation method.
 //  Trigger particle : Hadrons   Associated Particle : V0s
 //  
/// \author Kai Cui
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/Multiplicity.h"

namespace o2::aod
{
namespace assocK0shorts
{
// Needed to have shorter table that does not rely on existing one (filtering!)
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, Tracks, "_Pos"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, Tracks, "_Neg"); //!
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                         //!
DECLARE_SOA_INDEX_COLUMN(V0, v0);                                       //!

DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(MK0short, m, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
}
DECLARE_SOA_TABLE(AssocK0shorts,"AOD","K0SHORTS",o2::soa::Index<>, assocK0shorts::PosTrackId, assocK0shorts::NegTrackId, assocK0shorts::CollisionId, assocK0shorts::V0Id, assocK0shorts::Phi, assocK0shorts::Pt, assocK0shorts::MK0short, assocK0shorts::Eta);

namespace assocLambdas
{
// Needed to have shorter table that does not rely on existing one (filtering!)
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, Tracks, "_Pos"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, Tracks, "_Neg"); //!
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                         //!
DECLARE_SOA_INDEX_COLUMN(V0, v0);                                       //!

DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(MLambda, m, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
}
DECLARE_SOA_TABLE(AssocLambdas,"AOD","MYLAMBDAS",o2::soa::Index<>, assocLambdas::PosTrackId, assocLambdas::NegTrackId, assocLambdas::CollisionId, assocLambdas::V0Id, assocLambdas::Phi, assocLambdas::Pt, assocLambdas::MLambda, assocLambdas::Eta);

namespace assocV0s
{
// Needed to have shorter table that does not rely on existing one (filtering!)
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, Tracks, "_Pos"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, Tracks, "_Neg"); //!
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                         //!
DECLARE_SOA_INDEX_COLUMN(V0, v0);                                       //!

DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(MV0, m, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
}
DECLARE_SOA_TABLE(AssocV0s,"AOD","MYV0S",o2::soa::Index<>, assocV0s::PosTrackId, assocV0s::NegTrackId, assocV0s::CollisionId, assocV0s::V0Id, assocV0s::Phi, assocV0s::Pt, assocV0s::MV0, assocV0s::Eta);

namespace triggerTracks
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                         //!
DECLARE_SOA_COLUMN(TrackType, trackType, uint8_t);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Sign, sign, float);   //! (sign of charge)/Pt in c/GeV. Use pt() and sign() instead
}
DECLARE_SOA_TABLE(TriggerTracks,"AOD","MYTRACKS",o2::soa::Index<>, triggerTracks::CollisionId, triggerTracks::TrackType, triggerTracks::Phi, triggerTracks::Pt, triggerTracks::Eta, triggerTracks::Sign);
}

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;


Double_t ComputeDeltaPhi( Double_t phi1, Double_t phi2) {
//To be completely sure, use inner products
Double_t x1, y1, x2, y2;
   x1 = TMath::Cos( phi1 );
   y1 = TMath::Sin( phi1 );
   x2 = TMath::Cos( phi2 );
   y2 = TMath::Sin( phi2 );
Double_t lInnerProd  = x1*x2 + y1*y2;
Double_t lVectorProd = x1*y2 - x2*y1;

Double_t lReturnVal = 0;
   if( lVectorProd > 1e-8 ) {
       lReturnVal = TMath::ACos(lInnerProd);
   }
   if( lVectorProd < -1e-8 ) {
       lReturnVal = -TMath::ACos(lInnerProd);
   }

   if( lReturnVal < -TMath::Pi()/2. ) {
       lReturnVal += 2.*TMath::Pi();
   }

      return lReturnVal;
  }
 

struct vzerofilter {

  // Selection criteria: 5 basic V0 selection criteria!
  Configurable<double> v0Cospa{"v0cospa", 0.97, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcaV0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> dcaNegtopv{"dcanegtopv", 0.06, "DCA Neg To PV"};
  Configurable<float> dcaPostopv{"dcapostopv", 0.06, "DCA Pos To PV"};
  Configurable<float> v0RadiusMin{"v0radiusmin", 0.5, "v0radius"};
  Configurable<float> v0RadiusMax{"v0radiusmax", 200, "v0radius"};
  Configurable<float> triggerEtaMin{"triggerEtaCutMin", -0.8, "triggeretamin"};
  Configurable<float> triggerEtaMax{"triggerEtaCutMax", 0.8, "triggeretamax"};
  Configurable<float> assocEtaMin{"assocEtaCutMin", -0.8, "triggeretamin"};
  Configurable<float> assocEtaMax{"assocEtaCutMax", 0.8, "triggeretamax"};
  Configurable<float> triggerPtCutMin{"triggerPtCutMin", 1, "triggerptmin"};
  Configurable<float> triggerPtCutMax{"triggerPtCutMax", 3, "triggerptmax"};
  Configurable<float> assocPtCutMin{"assocPtCutMin", 1, "assocptmin"};
  Configurable<float> assocPtCutMax{"assocPtCutMax", 3, "assocptmax"};

  //Cannot filter on dynamic columns, so we cut on DCA to PV and DCA between daus only!
 Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > dcaPostopv&&
                       nabs(aod::v0data::dcanegtopv) > dcaNegtopv&&
                       aod::v0data::dcaV0daughters < dcaV0dau;
  
  // histogram defined with HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {      
      {"correlationHadronHadron", "correlationHadronHadron", {HistType::kTH1F, {{40,-0.5*M_PI, 1.5*M_PI,"#Phi"}}}},
      {"correlationHadronV0", "correlationHadronV0", {HistType::kTH1F, {{40,-0.5*M_PI, 1.5*M_PI,"#Phi"}}}},
      {"hVertexZ", "hVertexZ", {HistType::kTH1F, {{100, -15., 15.}}}},
      {"hV0Radius", "hV0Radius", {HistType::kTH1F, {{250, 0, 250}}}},
      {"hV0Eta", "hV0Eta", {HistType::kTH1F, {{200, -1, 1,"#Eta"}}}},
      {"hTrackEta", "hTrackEta", {HistType::kTH1F, {{200, -1, 1,"#Eta"}}}},
      {"hTrackSign", "hTrackSign", {HistType::kTH1F, {{5, -2, 2}}}},
      {"hV0dauDCA", "hV0dauDCA", {HistType::kTH1F, {{200, -1, 1}}}},
      {"hID", "hID", {HistType::kTH1F, {{20000, 0, 20000}}}},
      {"hV0CPA", "hV0CPA", {HistType::kTH1F, {{100, 0, 1}}}},
      {"hPosDCAtoPV", "hPosDCAtoPV", {HistType::kTH1F, {{400, 0.05, 0.45}}}},
      {"hNegDCAtoPV", "hNegDCAtoPV", {HistType::kTH1F, {{400, 0.05, 0.45}}}},
      {"hMassK0Short", "hMassK0Short", {HistType::kTH1F, {{200, 0.450f, 0.550f}}}},
      {"hMassLambda", "hMassLambda", {HistType::kTH1F, {{200, 1.0f, 1.550f}}}}
    }
  };
  using DauTracks = soa::Join<aod::Tracks,aod::TracksExtra,aod::pidTPCPi,aod::pidTPCPr>;

  Produces<aod::AssocLambdas> assocLambda;
  Produces<aod::AssocK0shorts> assocK0short;
  Produces<aod::AssocV0s> assocV0;
  Produces<aod::TriggerTracks> triggerTrack;

  void process(soa::Join<aod::Collisions, aod::EvSels,aod::Mults>::iterator const& collision,DauTracks const& tracks ,soa::Filtered<aod::V0Datas> const& V0s)
  {

    if (!collision.sel8()) {
      return;
    }
    for(auto& track : tracks)
      {
	if (track.eta()> triggerEtaMax || track.eta()< triggerEtaMin )  {continue;}
//	if (track.sign()= 1 ) {continue;}
	if (track.pt()> triggerPtCutMax || track.pt()< triggerPtCutMin ) {continue;}

        triggerTrack(
        track.collisionId(),
        track.trackType(),
        track.phi(),track.pt(),track.eta(),track.sign());

        registry.fill(HIST("hTrackEta"),track.eta());
        registry.fill(HIST("hTrackSign"),track.sign());
        registry.fill(HIST("hID"),track.collisionId());
      }
    //Basic event selection (all helper tasks are now included!)
    //check getter here: https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html
    registry.get<TH1>(HIST("hVertexZ"))->Fill(collision.posZ());
    for (auto& v0 : V0s) {
	 
	 auto posdau = v0.posTrack_as<DauTracks>();
	 auto negdau = v0.negTrack_as<DauTracks>();

	 if (v0.v0radius() < v0RadiusMin ||v0.v0radius() > v0RadiusMax||v0.eta() > assocEtaMax ||v0.eta() < assocEtaMin|| v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0Cospa){
      continue;
	 }
         registry.fill(HIST("hV0Radius"), v0.v0radius());
         registry.fill(HIST("hV0Eta"), v0.eta());
         registry.fill(HIST("hV0dauDCA"), v0.dcaV0daughters());
         registry.fill(HIST("hPosDCAtoPV"), v0.dcapostopv());
         registry.fill(HIST("hNegDCAtoPV"), v0.dcanegtopv());
         registry.fill(HIST("hV0CPA"), v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));

	 if (TMath::Abs(posdau.tpcNSigmaPi())<5 && TMath::Abs(negdau.tpcNSigmaPi())<5){
	    registry.fill(HIST("hMassK0Short"), v0.mK0Short());

        assocK0short(    
		v0.posTrackId(),
        v0.negTrackId(),
        v0.collisionId(),v0.globalIndex(),
		v0.phi(),v0.pt(),v0.mK0Short(),v0.eta());

	 }
	 if (TMath::Abs(posdau.tpcNSigmaPr())<5 && TMath::Abs(negdau.tpcNSigmaPi())<5){ 
	 registry.fill(HIST("hMassLambda"), v0.mLambda());
		assocLambda(        
	    v0.posTrackId(),
        v0.negTrackId(),
        v0.collisionId(),
		v0.globalIndex(),v0.phi(),v0.pt(),v0.mLambda(),v0.eta());

	 }
		assocV0(
        v0.posTrackId(),
        v0.negTrackId(),
        v0.collisionId(),
        v0.globalIndex(),
        v0.phi(),v0.pt(),v0.mLambda(),v0.eta());

    }
      for (auto& [trackTrigger, trackAssoc] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, tracks))) {
      registry.get<TH1>(HIST("correlationHadronHadron"))->Fill( ComputeDeltaPhi(trackTrigger.phi(), trackAssoc.phi() ));
    }
      for (auto& [trackTrigger, trackAssoc] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, V0s))) {
      registry.get<TH1>(HIST("correlationHadronV0"))->Fill( ComputeDeltaPhi(trackTrigger.phi(), trackAssoc.phi() ));
    }
  }
};

struct MixedEvents {


   HistogramRegistry registry{
      "registry",
      {
        {"correlationHadronV0", "correlationHadronV0", {HistType::kTH2F, {{40,-0.5*M_PI, 1.5*M_PI,"#deltaPhi"},{400,-4,4,"#deltaEta"}}}},
        {"deltaPhiDistribution", "deltaPhiDistribution", {HistType::kTH1F, {{40,-0.5*M_PI, 1.5*M_PI,"#deltaPhi"}}}},
        {"deltaEtaDistribution", "deltaEtaDistribution", {HistType::kTH1F, {{400,-4, 4,"#deltaEta"}}}},
        {"hEta1", "hV0Eta", {HistType::kTH1F, {{200, -1, 1,"#Eta"}}}},
        {"hEta2", "hV0Eta", {HistType::kTH1F, {{200, -1, 1,"#Eta"}}}},
        {"hPhi1", "hPhi", {HistType::kTH1F, {{200, 0, 8,"#Phi"}}}},
        {"hPhi2", "hPhi", {HistType::kTH1F, {{200, 0, 8,"#Phi"}}}},
        {"hId1", "hId1", {HistType::kTH1F, {{2000, 0, 200}}}},
        {"hId2", "hId2", {HistType::kTH1F, {{2000, 0, 200}}}},
      }
    };
 // ConfigurableAxis ConfMultBins{"ConfMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis ConfMultBins{"ConfMultBins", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1}, "Mixing bins - multiplicity"};
  ConfigurableAxis ConfVtxBins{"ConfVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

 // using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFT0A>;
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultFV0M<aod::mult::MultFV0A, aod::mult::MultFV0C>>;
  BinningType colBinning{{ConfVtxBins, ConfMultBins}, true};
  Pair< soa::Join<aod::Collisions, aod::EvSels, aod::Mults> , aod::TriggerTracks, aod::AssocV0s, BinningType> pair{colBinning, 5, -1}; // indicates that 5 events should be mixed and under/overflow (-1) to be ignored

  void process(soa::Join<aod::Collisions, aod::EvSels, aod::Mults> const& collisions,aod::TriggerTracks const& tracks,aod::AssocV0s const& assocv0s)
  {
    int count = 0;
    LOGF(info, "Input data Collisions %d, Tracks %d V0s %d", collisions.size(), tracks.size(), assocv0s.size());
    for (auto& [c1, tracks1, c2, tracks2] : pair) {
      LOGF(info, "Mixed event collisions: (%d, %d)", c1.globalIndex(), c2.globalIndex());
       count++;
     // if (count == 1000) { break;}
	int trackCount = 0;
	for(auto& track1 : tracks1)
	{
	        registry.get<TH1>(HIST("hId1"))->Fill( track1.collisionId());
	}
        for(auto& track2 : tracks2)
        {
                registry.get<TH1>(HIST("hId2"))->Fill( track2.collisionId());
        }

      // Example of using tracks from mixed events -- iterate over all track pairs from the two collisions
      for (auto& [t1, t2] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d), track event: (%d, %d)", t1.collisionId(), t2.collisionId(), c1.index(), c2.index(), t1.collision().index(), t2.collision().index());
	registry.get<TH1>(HIST("deltaPhiDistribution"))->Fill( ComputeDeltaPhi(t1.phi(), t2.phi() ));
	registry.get<TH1>(HIST("deltaEtaDistribution"))->Fill(t1.eta() - t2.eta());
	registry.fill(HIST("correlationHadronV0"),ComputeDeltaPhi(t1.phi(), t2.phi() ),t1.eta() - t2.eta());
        trackCount++;
       // if (trackCount == 100)
        //  break;
      }
    }
  }
};

struct SameEvents {
	

   HistogramRegistry registry{
      "registry",
      {
      	{"correlationHadronK0short", "correlationFunction", {HistType::kTH1F, {{40,-0.5*M_PI, 1.5*M_PI,"#Phi"}}}},
      	{"correlationHadronLambda", "correlationFunction2", {HistType::kTH1F, {{40,-0.5*M_PI, 1.5*M_PI,"#Phi"}}}},
      	{"correlationHadronV0", "correlationFunction3", {HistType::kTH1F, {{40,-0.5*M_PI, 1.5*M_PI,"#Phi"}}}},
        {"hMassK0Short", "hMassK0Short", {HistType::kTH1F, {{200, 0.450f, 0.550f}}}},
	    {"hV0Eta", "hV0Eta", {HistType::kTH1F, {{200, -1, 1,"#Eta"}}}},
	    {"hTrackEta", "hTrackEta", {HistType::kTH1F, {{200, -1, 1,"#Eta"}}}},
	    {"hID", "hID", {HistType::kTH1F, {{2000, 0, 2000}}}},
        {"hPtK0Short", "hPtK0Short", {HistType::kTH1F, {{200, 0, 15,"#it{p}_{T} (GeV/c)"}}}},
        {"hPtLambda", "hPtLambda", {HistType::kTH1F, {{200, 0, 15,"#it{p}_{T} (GeV/c)"}}}},
        {"hMassLambda", "hMassLambda", {HistType::kTH1F, {{200, 0.85f, 1.5}}}}
      }
    };

    void process(soa::Join<aod::Collisions, aod::EvSels,aod::Mults>::iterator const& collision,aod::TriggerTracks const& tracks,aod::AssocLambdas const& lambdas, aod::AssocK0shorts const& k0shorts, aod::AssocV0s const& v0s)
    {

	if (!collision.sel8()) {
      return;
    }
        for (auto& lambda : lambdas)
        {

            registry.fill(HIST("hMassLambda"),lambda.m());
            registry.fill(HIST("hPtLambda"),lambda.pt());
	    }
        for (auto& k0short : k0shorts)
        {

            registry.fill(HIST("hMassK0Short"),k0short.m());
            registry.fill(HIST("hPtK0Short"),k0short.pt());
        }
	    for (auto& v0 : v0s)
	    {
		    registry.fill(HIST("hV0Eta"),v0.eta());
	    }
        for (auto& track : tracks)
        {
                registry.fill(HIST("hTrackEta"),track.eta());
        }

        for (auto& [trackTrigger, trackAssoc] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, lambdas))) 
        {
            registry.get<TH1>(HIST("correlationHadronLambda"))->Fill( ComputeDeltaPhi(trackTrigger.phi(), trackAssoc.phi() ));
            registry.get<TH1>(HIST("hID"))->Fill(trackTrigger.collisionId());
//        LOGF(info, "corelation event tracks pair: (%d, %d)", trackTrigger.index(), trackAssoc.index());
        }

        for (auto& [trackTrigger, trackAssoc] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, k0shorts))) 
        {
            registry.get<TH1>(HIST("correlationHadronK0short"))->Fill( ComputeDeltaPhi(trackTrigger.phi(), trackAssoc.phi() ));
            registry.get<TH1>(HIST("hID"))->Fill(trackTrigger.collisionId());
        }

        for (auto& [trackTrigger, trackAssoc] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, v0s)))
        {
            registry.get<TH1>(HIST("correlationHadronV0"))->Fill( ComputeDeltaPhi(trackTrigger.phi(), trackAssoc.phi() ));
        }

   }

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<vzerofilter>(cfgc),
    adaptAnalysisTask<MixedEvents>(cfgc),
    adaptAnalysisTask<SameEvents>(cfgc)
  };
}
