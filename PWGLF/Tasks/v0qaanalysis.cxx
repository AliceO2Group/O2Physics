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
/// \brief QA task for V0 analysis using derived data
///
/// \author Francesca Ercolessi (francesca.ercolessi@cern.ch)

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using DauTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFPr>;

namespace o2::aod
{

namespace myv0candidates
{

DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(V0Pt, v0pt, float);
DECLARE_SOA_COLUMN(RapLambda, raplambda, float);
DECLARE_SOA_COLUMN(RapK0Short, rapk0short, float);
DECLARE_SOA_COLUMN(MassLambda, masslambda, float);
DECLARE_SOA_COLUMN(MassAntiLambda, massantilambda, float);
DECLARE_SOA_COLUMN(MassK0Short, massk0short, float);
DECLARE_SOA_COLUMN(V0Radius, v0radius, float);
DECLARE_SOA_COLUMN(V0CosPA, v0cospa, float);
DECLARE_SOA_COLUMN(V0DCAPosToPV, v0dcapostopv, float);
DECLARE_SOA_COLUMN(V0DCANegToPV, v0dcanegtopv, float);
DECLARE_SOA_COLUMN(V0DCAV0Daughters, v0dcav0daughters, float);
DECLARE_SOA_COLUMN(V0PosEta, v0poseta, float);
DECLARE_SOA_COLUMN(V0NegEta, v0negeta, float);
DECLARE_SOA_COLUMN(V0PosPhi, v0posphi, float);
DECLARE_SOA_COLUMN(V0NegPhi, v0negphi, float);
DECLARE_SOA_COLUMN(V0PosITSHits, v0positshits, float);
DECLARE_SOA_COLUMN(V0NegITSHits, v0negitshits, float);
DECLARE_SOA_COLUMN(CtauLambda, ctaulambda, float);
DECLARE_SOA_COLUMN(CtauAntiLambda, ctauantilambda, float);
DECLARE_SOA_COLUMN(CtauK0Short, ctauk0short, float);
DECLARE_SOA_COLUMN(NTPCSigmaNegPr, ntpcsigmanegpr, float);
DECLARE_SOA_COLUMN(NTPCSigmaPosPr, ntpcsigmapospr, float);
DECLARE_SOA_COLUMN(NTPCSigmaNegPi, ntpcsigmanegpi, float);
DECLARE_SOA_COLUMN(NTPCSigmaPosPi, ntpcsigmapospi, float);
DECLARE_SOA_COLUMN(NTOFSigmaNegPr, ntofsigmanegpr, float);
DECLARE_SOA_COLUMN(NTOFSigmaPosPr, ntofsigmapospr, float);
DECLARE_SOA_COLUMN(NTOFSigmaNegPi, ntofsigmanegpi, float);
DECLARE_SOA_COLUMN(NTOFSigmaPosPi, ntofsigmapospi, float);

} // namespace myv0candidates

DECLARE_SOA_TABLE(MyV0Candidates, "AOD", "MYV0CANDIDATES", o2::soa::Index<>,
                  myv0candidates::CollisionId, myv0candidates::V0Pt, myv0candidates::RapLambda, myv0candidates::RapK0Short,
                  myv0candidates::MassLambda, myv0candidates::MassAntiLambda, myv0candidates::MassK0Short,
                  myv0candidates::V0Radius, myv0candidates::V0CosPA, myv0candidates::V0DCAPosToPV,
                  myv0candidates::V0DCANegToPV, myv0candidates::V0DCAV0Daughters,
                  myv0candidates::V0PosEta, myv0candidates::V0NegEta, myv0candidates::V0PosPhi, myv0candidates::V0NegPhi,
                  myv0candidates::V0PosITSHits, myv0candidates::V0NegITSHits, myv0candidates::CtauLambda, myv0candidates::CtauAntiLambda, myv0candidates::CtauK0Short,
                  myv0candidates::NTPCSigmaNegPr, myv0candidates::NTPCSigmaPosPr, myv0candidates::NTPCSigmaNegPi, myv0candidates::NTPCSigmaPosPi,
                  myv0candidates::NTOFSigmaNegPr, myv0candidates::NTOFSigmaPosPr, myv0candidates::NTOFSigmaNegPi, myv0candidates::NTOFSigmaPosPi);

} // namespace o2::aod

struct v0qaanalysis {

  // Produces
  Produces<aod::MyV0Candidates> myv0s;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    registry.add("hNEvents", "hNEvents", {HistType::kTH1I, {{1, 0.f, 1.f}}});
  }

  // Event selection criteria
  Configurable<float> cutzvertex{"cutzvertex", 15.0f, "Accepted z-vertex range (cm)"};
  Configurable<bool> sel8{"sel8", 1, "Apply sel8 event selection"};

  // V0 selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.97, "V0 CosPA"};
  Configurable<float> dcav0dau{"dcav0dau", 1.5, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", 0.06, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.06, "DCA Pos To PV"};
  Configurable<float> v0radius{"v0radius", 0.5, "Radius"};
  Configurable<float> rapidity{"rapidity", 0.5, "Rapidity"};
  Configurable<float> etadau{"etadau", 0.8, "Eta Daughters"};

  Filter preFilterV0 = nabs(aod::v0data::dcapostopv) > dcapostopv&&
                                                         nabs(aod::v0data::dcanegtopv) > dcanegtopv&& aod::v0data::dcaV0daughters < dcav0dau;

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Filtered<aod::V0Datas> const& V0s, DauTracks const& tracks)
  {
    // Event selection
    if (sel8 && !collision.sel8()) {
      return;
    }
    if (TMath::Abs(collision.posZ()) > cutzvertex) {
      return;
    }

    registry.fill(HIST("hNEvents"), 0.5);

    for (auto& v0 : V0s) { // loop over V0s

      float ctauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::getMassPDG(3122);
      float ctauAntiLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::getMassPDG(-3122);
      float ctauK0s = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::getMassPDG(310);

      int posITSNhits = 0, negITSNhits = 0;
      for (unsigned int i = 0; i < 7; i++) {
        if (v0.posTrack_as<DauTracks>().itsClusterMap() & (1 << i)) {
          posITSNhits++;
        }
        if (v0.negTrack_as<DauTracks>().itsClusterMap() & (1 << i)) {
          negITSNhits++;
        }
      }

      if (v0.v0radius() > v0radius &&
          v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa &&
          TMath::Abs(v0.posTrack_as<DauTracks>().eta()) < etadau &&
          TMath::Abs(v0.negTrack_as<DauTracks>().eta()) < etadau) {

        // Fill table
        myv0s(v0.globalIndex(), v0.pt(), v0.yLambda(), v0.yK0Short(),
              v0.mLambda(), v0.mAntiLambda(), v0.mK0Short(),
              v0.v0radius(), v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()),
              v0.dcapostopv(), v0.dcanegtopv(), v0.dcaV0daughters(),
              v0.posTrack_as<DauTracks>().eta(), v0.negTrack_as<DauTracks>().eta(),
              v0.posTrack_as<DauTracks>().phi(), v0.negTrack_as<DauTracks>().phi(),
              posITSNhits, negITSNhits, ctauLambda, ctauAntiLambda, ctauK0s,
              v0.negTrack_as<DauTracks>().tpcNSigmaPr(), v0.posTrack_as<DauTracks>().tpcNSigmaPr(),
              v0.negTrack_as<DauTracks>().tpcNSigmaPi(), v0.posTrack_as<DauTracks>().tpcNSigmaPi(),
              v0.negTrack_as<DauTracks>().tofNSigmaPr(), v0.posTrack_as<DauTracks>().tofNSigmaPr(),
              v0.negTrack_as<DauTracks>().tofNSigmaPi(), v0.posTrack_as<DauTracks>().tofNSigmaPi());
      }
    }
  }
};

struct myV0s {

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    registry.add("hMassLambda", "hMassLambda", {HistType::kTH1F, {{200, 1.016f, 1.216f}}});
    registry.add("hPt", "hPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}});
    registry.add("hMassVsPtLambda", "hMassVsPtLambda", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {200, 1.016f, 1.216f}}});
    registry.add("hMassAntiLambda", "hMassAntiLambda", {HistType::kTH1F, {{200, 1.016f, 1.216f}}});
    registry.add("hMassVsPtAntiLambda", "hMassVsPtAntiLambda", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {200, 1.016f, 1.216f}}});
    registry.add("hMassK0Short", "hMassK0Short", {HistType::kTH1F, {{200, 0.4f, 0.6f}}});
    registry.add("hMassVsPtK0Short", "hMassVsPtK0Short", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {200, 0.4f, 0.6f}}});
    registry.add("V0Radius", "V0Radius", {HistType::kTH1D, {{100, 0.0f, 20.0f}}});
    registry.add("CosPA", "CosPA", {HistType::kTH1F, {{100, 0.9f, 1.0f}}});
    registry.add("V0DCANegToPV", "V0DCANegToPV", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
    registry.add("V0DCAPosToPV", "V0DCAPosToPV", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
    registry.add("V0DCAV0Daughters", "V0DCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.20f}}});
    registry.add("CtauK0s", "CtauK0s", {HistType::kTH1F, {{150, 0.0f, 30.0f}}});
    registry.add("CtauLambda", "CtauLambda", {HistType::kTH1F, {{200, 0.0f, 40.0f}}});
    registry.add("CtauAntiLambda", "CtauAntiLambda", {HistType::kTH1F, {{200, 0.0f, 40.0f}}});
    registry.add("TPCNSigmaPosPi", "TPCNSigmaPosPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("TPCNSigmaNegPi", "TPCNSigmaNegPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("TPCNSigmaPosPr", "TPCNSigmaPosPr", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("TPCNSigmaNegPr", "TPCNSigmaNegPr", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("PosITSHits", "PosITSHits", {HistType::kTH1F, {{8, -0.5f, 7.5f}}});
    registry.add("NegITSHits", "NegITSHits", {HistType::kTH1F, {{8, -0.5f, 7.5f}}});
  }

  void process(aod::MyV0Candidates const& myv0s)
  {
    for (auto& candidate : myv0s) {

      registry.fill(HIST("hMassLambda"), candidate.masslambda());
      registry.fill(HIST("hPt"), candidate.v0pt());
      registry.fill(HIST("hMassVsPtLambda"), candidate.v0pt(), candidate.masslambda());
      registry.fill(HIST("hMassAntiLambda"), candidate.massantilambda());
      registry.fill(HIST("hMassVsPtAntiLambda"), candidate.v0pt(), candidate.massantilambda());
      registry.fill(HIST("hMassK0Short"), candidate.massk0short());
      registry.fill(HIST("hMassVsPtK0Short"), candidate.v0pt(), candidate.massk0short());
      registry.fill(HIST("V0Radius"), candidate.v0radius());
      registry.fill(HIST("CosPA"), candidate.v0cospa());
      registry.fill(HIST("V0DCANegToPV"), candidate.v0dcanegtopv());
      registry.fill(HIST("V0DCAPosToPV"), candidate.v0dcapostopv());
      registry.fill(HIST("V0DCAV0Daughters"), candidate.v0dcav0daughters());
      registry.fill(HIST("CtauK0s"), candidate.ctauk0short());
      registry.fill(HIST("CtauLambda"), candidate.ctaulambda());
      registry.fill(HIST("CtauAntiLambda"), candidate.ctauantilambda());
      registry.fill(HIST("TPCNSigmaPosPi"), candidate.ntpcsigmapospi());
      registry.fill(HIST("TPCNSigmaNegPi"), candidate.ntpcsigmanegpi());
      registry.fill(HIST("TPCNSigmaPosPr"), candidate.ntpcsigmapospr());
      registry.fill(HIST("PosITSHits"), candidate.v0positshits());
      registry.fill(HIST("NegITSHits"), candidate.v0negitshits());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<v0qaanalysis>(cfgc, TaskName{"lf-v0qaanalysis"}),
    adaptAnalysisTask<myV0s>(cfgc, TaskName{"lf-myv0s"})};
}
