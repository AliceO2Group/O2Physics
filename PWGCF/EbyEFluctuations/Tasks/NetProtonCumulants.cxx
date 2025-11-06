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
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include <CCDB/BasicCCDBManager.h>

#include "TList.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TRandom3.h"

#include <cmath>
#include <vector>

namespace o2::aod
{
namespace netprotonNum
{
DECLARE_SOA_COLUMN(NetProtNo, net_prot_no, float); //! net proton no. in an event
DECLARE_SOA_COLUMN(N_ch, n_ch, float);             //! no of charged particles/multiplicity in an event
DECLARE_SOA_COLUMN(Centrality, centrality, float); //! Centrality of event
} // namespace netprotonNum
DECLARE_SOA_TABLE(NetProton, "AOD", "NETPROTONNUM", netprotonNum::NetProtNo, netprotonNum::N_ch, netprotonNum::Centrality); //! table to store e-by-e net-proton numbers, multiplicity and centrality
} // namespace o2::aod

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct NetProtonCumulants_Table_QA {

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutPtLower{"cfgCutPtLower", 0.2f, "Lower pT cut"};
  Configurable<float> cfgCutPtUpper{"cfgCutPtUpper", 3.0f, "Higher pT cut"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "absolute Eta cut"};
  Configurable<float> cfgCutPtUpperTPC{"cfgCutPtUpperTPC", 0.6f, "Upper pT cut for PID using TPC only"};
  Configurable<float> cfgCutTpcChi2NCl{"cfgCutTpcChi2NCl", 2.5f, "Maximum TPCchi2NCl"};
  Configurable<float> cfgnSigmaCut{"cfgnSigmaCut", 2.0f, "PID nSigma cut"};

  // Filter command***********
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < 0.8f) && (aod::track::pt > cfgCutPtLower) && (aod::track::pt < 5.0f) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutTpcChi2NCl);

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // filtering collisions and tracks***********
  using aodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFDDMs>>;
  using aodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::pidTPCPr, aod::pidTOFPr, aod::pidTPCKa, aod::pidTOFKa, aod::pidTPCPi, aod::pidTOFPi>>;

  // Equivalent of the AliRoot task UserCreateOutputObjects
  void init(o2::framework::InitContext&)
  {
    // Define your axes
    // Constant bin width axis
    AxisSpec vtxZAxis = {100, -20, 20};
    // Variable bin width axis
    std::vector<double> ptBinning = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4.};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};
    std::vector<double> centBining = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90};
    AxisSpec centAxis = {centBining, "centrality (%)"};
    AxisSpec netProtonAxis = {2001, -1000.5, 1000.5, "net-proton number"};

    // Add histograms to histogram manager (as in the output object of in AliPhysics)
    histos.add("hZvtx_after_sel", ";Z (cm)", kTH1F, {vtxZAxis});
    histos.add("hPtAll", ";#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
    histos.add("hPtProton", ";#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
    histos.add("hPtAntiproton", ";#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
    histos.add("hPhiAll", ";#phi", kTH1F, {{100, 0., 2. * M_PI}});
    histos.add("hPhiProton", ";#phi", kTH1F, {{100, 0., 2. * M_PI}});
    histos.add("hPhiAntiproton", ";#phi", kTH1F, {{100, 0., 2. * M_PI}});
    histos.add("hEtaAll", ";#eta", kTH1F, {{100, -2.01, 2.01}});
    histos.add("hEtaProton", ";#eta", kTH1F, {{100, -2.01, 2.01}});
    histos.add("hEtaAntiproton", ";#eta", kTH1F, {{100, -2.01, 2.01}});
    histos.add("hCentrality", ";centrality (%)", kTH1F, {{90, 0, 90}});
    histos.add("hNetProtonVsCentrality", "", kTH2F, {netProtonAxis, centAxis});
    histos.add("hProfileTotalProton", "", kTProfile, {centAxis});
  }
  Produces<aod::NetProton> net_proton_num; //! table creation

  // void process(aod::Collision const& coll, aod::Tracks const& inputTracks)
  void process(aodCollisions::iterator const& coll, aod::BCsWithTimestamps const&, aodTracks const& inputTracks)
  {
    histos.fill(HIST("hZvtx_after_sel"), coll.posZ());
    histos.fill(HIST("hCentrality"), coll.centFT0C());
    // variables
    float cent = coll.centFT0C();
    float n_ch = 0.0;
    float n_prot = 0.0;
    float n_antiprot = 0.0;

    //! centrality cut
    if (cent > 0.0f && cent < 90.0f) {

      for (auto track : inputTracks) { //! Loop over tracks
        histos.fill(HIST("hPtAll"), track.pt());
        histos.fill(HIST("hEtaAll"), track.eta());
        histos.fill(HIST("hPhiAll"), track.phi());

        //! PID checking
        int flag = 0;
        if (track.pt() > 0.2f && track.pt() <= cfgCutPtUpperTPC) {
          if (track.tpcNSigmaPr() < cfgnSigmaCut) {
            flag = 1;
          }
        }
        if (track.pt() > cfgCutPtUpperTPC && track.pt() < 5.0f) {
          const float combNSigmaPr = std::sqrt(pow(track.tpcNSigmaPr(), 2.0) + pow(track.tofNSigmaPr(), 2.0));
          const float combNSigmaPi = std::sqrt(pow(track.tpcNSigmaPi(), 2.0) + pow(track.tofNSigmaPi(), 2.0));
          const float combNSigmaKa = std::sqrt(pow(track.tpcNSigmaKa(), 2.0) + pow(track.tofNSigmaKa(), 2.0));

          int flag2 = 0;
          if (combNSigmaPr < 3.0)
            flag2 += 1;
          if (combNSigmaPi < 3.0)
            flag2 += 1;
          if (combNSigmaKa < 3.0)
            flag2 += 1;
          if (!(flag2 > 1) && !(combNSigmaPr > combNSigmaPi) && !(combNSigmaPr > combNSigmaKa)) {
            if (combNSigmaPr < cfgnSigmaCut) {
              flag = 1;
            }
          }
        }

        if (track.pt() > cfgCutPtLower && track.pt() < cfgCutPtUpper && track.sign() != 0 && TMath::Abs(track.eta()) < cfgCutEta) {
          if (flag == 1) {
            if (track.sign() > 0) {
              histos.fill(HIST("hPtProton"), track.pt());
              histos.fill(HIST("hEtaProton"), track.eta());
              histos.fill(HIST("hPhiProton"), track.phi());
              n_prot = n_prot + 1.0; //! calculating no. of proton
            }
            if (track.sign() < 0) {
              histos.fill(HIST("hPtAntiproton"), track.pt());
              histos.fill(HIST("hEtaAntiproton"), track.eta());
              histos.fill(HIST("hPhiAntiproton"), track.phi());
              n_antiprot = n_antiprot + 1.0; //! calculating no. of anti-proton
            }
          }
          n_ch = n_ch + 1; //! calculating no. of charged particles
        }
      } //! end loop on tracks

      float net_prot = n_prot - n_antiprot;
      net_proton_num(net_prot, n_ch, cent);
      histos.fill(HIST("hNetProtonVsCentrality"), net_prot, cent);
      histos.fill(HIST("hProfileTotalProton"), cent, (n_prot + n_antiprot));
    }

  } //! end process loop
};

struct NetProtonCumulants_analysis {

  Configurable<int> cfgNSubsample{"cfgNSubsample", 10, "Number of subsamples"};
  ConfigurableAxis centAxis{"centAxis", {90, 0, 90}, ""};
  ConfigurableAxis multAxis{"multAxis", {5000, 0.5, 5000.5}, ""};

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};

  // Define output
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  std::vector<std::vector<std::shared_ptr<TProfile2D>>> Subsample2D;
  std::vector<std::vector<std::shared_ptr<TProfile>>> Subsample;
  TRandom3* fRndm = new TRandom3(0);

  void init(o2::framework::InitContext&)
  {
    // AxisSpec centAxis = {90, 0, 90, "centrality (%)"};
    // AxisSpec multAxis = {5000, 0.5, 5000.5, "#it{N}_{ch,acc}"};

    registry.add("Prof_mu1_netproton", "", {HistType::kTProfile, {centAxis}});
    registry.add("Prof_mu2_netproton", "", {HistType::kTProfile, {centAxis}});
    registry.add("Prof_mu3_netproton", "", {HistType::kTProfile, {centAxis}});
    registry.add("Prof_mu4_netproton", "", {HistType::kTProfile, {centAxis}});
    registry.add("Prof_mu5_netproton", "", {HistType::kTProfile, {centAxis}});
    registry.add("Prof_mu6_netproton", "", {HistType::kTProfile, {centAxis}});
    registry.add("Prof_mu7_netproton", "", {HistType::kTProfile, {centAxis}});
    registry.add("Prof_mu8_netproton", "", {HistType::kTProfile, {centAxis}});
    registry.add("Prof2D_mu1_netproton", "", {HistType::kTProfile2D, {centAxis, multAxis}});
    registry.add("Prof2D_mu2_netproton", "", {HistType::kTProfile2D, {centAxis, multAxis}});
    registry.add("Prof2D_mu3_netproton", "", {HistType::kTProfile2D, {centAxis, multAxis}});
    registry.add("Prof2D_mu4_netproton", "", {HistType::kTProfile2D, {centAxis, multAxis}});
    registry.add("Prof2D_mu5_netproton", "", {HistType::kTProfile2D, {centAxis, multAxis}});
    registry.add("Prof2D_mu6_netproton", "", {HistType::kTProfile2D, {centAxis, multAxis}});
    registry.add("Prof2D_mu7_netproton", "", {HistType::kTProfile2D, {centAxis, multAxis}});
    registry.add("Prof2D_mu8_netproton", "", {HistType::kTProfile2D, {centAxis, multAxis}});

    // initial array
    Subsample2D.resize(cfgNSubsample);
    Subsample.resize(cfgNSubsample);
    for (int i = 0; i < cfgNSubsample; i++) {
      Subsample2D[i].resize(8);
      Subsample[i].resize(8);
    }
    for (int i = 0; i < cfgNSubsample; i++) {
      //! 2D profiles of moments
      Subsample2D[i][0] = std::get<std::shared_ptr<TProfile2D>>(registry.add(Form("Subsample_%d/Prof2D_mu1_netproton", i), "", {HistType::kTProfile2D, {centAxis, multAxis}}));
      Subsample2D[i][1] = std::get<std::shared_ptr<TProfile2D>>(registry.add(Form("Subsample_%d/Prof2D_mu2_netproton", i), "", {HistType::kTProfile2D, {centAxis, multAxis}}));
      Subsample2D[i][2] = std::get<std::shared_ptr<TProfile2D>>(registry.add(Form("Subsample_%d/Prof2D_mu3_netproton", i), "", {HistType::kTProfile2D, {centAxis, multAxis}}));
      Subsample2D[i][3] = std::get<std::shared_ptr<TProfile2D>>(registry.add(Form("Subsample_%d/Prof2D_mu4_netproton", i), "", {HistType::kTProfile2D, {centAxis, multAxis}}));
      Subsample2D[i][4] = std::get<std::shared_ptr<TProfile2D>>(registry.add(Form("Subsample_%d/Prof2D_mu5_netproton", i), "", {HistType::kTProfile2D, {centAxis, multAxis}}));
      Subsample2D[i][5] = std::get<std::shared_ptr<TProfile2D>>(registry.add(Form("Subsample_%d/Prof2D_mu6_netproton", i), "", {HistType::kTProfile2D, {centAxis, multAxis}}));
      Subsample2D[i][6] = std::get<std::shared_ptr<TProfile2D>>(registry.add(Form("Subsample_%d/Prof2D_mu7_netproton", i), "", {HistType::kTProfile2D, {centAxis, multAxis}}));
      Subsample2D[i][7] = std::get<std::shared_ptr<TProfile2D>>(registry.add(Form("Subsample_%d/Prof2D_mu8_netproton", i), "", {HistType::kTProfile2D, {centAxis, multAxis}}));
      //! 1D profiles of moments
      Subsample[i][0] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("Subsample_%d/Prof_mu1_netproton", i), "", {HistType::kTProfile, {centAxis}}));
      Subsample[i][1] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("Subsample_%d/Prof_mu2_netproton", i), "", {HistType::kTProfile, {centAxis}}));
      Subsample[i][2] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("Subsample_%d/Prof_mu3_netproton", i), "", {HistType::kTProfile, {centAxis}}));
      Subsample[i][3] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("Subsample_%d/Prof_mu4_netproton", i), "", {HistType::kTProfile, {centAxis}}));
      Subsample[i][4] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("Subsample_%d/Prof_mu5_netproton", i), "", {HistType::kTProfile, {centAxis}}));
      Subsample[i][5] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("Subsample_%d/Prof_mu6_netproton", i), "", {HistType::kTProfile, {centAxis}}));
      Subsample[i][6] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("Subsample_%d/Prof_mu7_netproton", i), "", {HistType::kTProfile, {centAxis}}));
      Subsample[i][7] = std::get<std::shared_ptr<TProfile>>(registry.add(Form("Subsample_%d/Prof_mu8_netproton", i), "", {HistType::kTProfile, {centAxis}}));
    }
  }

  void process(aod::NetProton::iterator const& event_netproton)
  {
    // LOGF(info, "Centrality= %f Nch= %f net-proton no. = %f", event_netproton.centrality(), event_netproton.n_ch(), event_netproton.net_prot_no());

    // filling profiles for central values
    registry.get<TProfile2D>(HIST("Prof2D_mu1_netproton"))->Fill(event_netproton.centrality(), event_netproton.n_ch(), pow(event_netproton.net_prot_no(), 1.0));
    registry.get<TProfile2D>(HIST("Prof2D_mu2_netproton"))->Fill(event_netproton.centrality(), event_netproton.n_ch(), pow(event_netproton.net_prot_no(), 2.0));
    registry.get<TProfile2D>(HIST("Prof2D_mu3_netproton"))->Fill(event_netproton.centrality(), event_netproton.n_ch(), pow(event_netproton.net_prot_no(), 3.0));
    registry.get<TProfile2D>(HIST("Prof2D_mu4_netproton"))->Fill(event_netproton.centrality(), event_netproton.n_ch(), pow(event_netproton.net_prot_no(), 4.0));
    registry.get<TProfile2D>(HIST("Prof2D_mu5_netproton"))->Fill(event_netproton.centrality(), event_netproton.n_ch(), pow(event_netproton.net_prot_no(), 5.0));
    registry.get<TProfile2D>(HIST("Prof2D_mu6_netproton"))->Fill(event_netproton.centrality(), event_netproton.n_ch(), pow(event_netproton.net_prot_no(), 6.0));
    registry.get<TProfile2D>(HIST("Prof2D_mu7_netproton"))->Fill(event_netproton.centrality(), event_netproton.n_ch(), pow(event_netproton.net_prot_no(), 7.0));
    registry.get<TProfile2D>(HIST("Prof2D_mu8_netproton"))->Fill(event_netproton.centrality(), event_netproton.n_ch(), pow(event_netproton.net_prot_no(), 8.0));

    registry.get<TProfile>(HIST("Prof_mu1_netproton"))->Fill(event_netproton.centrality(), pow(event_netproton.net_prot_no(), 1.0));
    registry.get<TProfile>(HIST("Prof_mu2_netproton"))->Fill(event_netproton.centrality(), pow(event_netproton.net_prot_no(), 2.0));
    registry.get<TProfile>(HIST("Prof_mu3_netproton"))->Fill(event_netproton.centrality(), pow(event_netproton.net_prot_no(), 3.0));
    registry.get<TProfile>(HIST("Prof_mu4_netproton"))->Fill(event_netproton.centrality(), pow(event_netproton.net_prot_no(), 4.0));
    registry.get<TProfile>(HIST("Prof_mu5_netproton"))->Fill(event_netproton.centrality(), pow(event_netproton.net_prot_no(), 5.0));
    registry.get<TProfile>(HIST("Prof_mu6_netproton"))->Fill(event_netproton.centrality(), pow(event_netproton.net_prot_no(), 6.0));
    registry.get<TProfile>(HIST("Prof_mu7_netproton"))->Fill(event_netproton.centrality(), pow(event_netproton.net_prot_no(), 7.0));
    registry.get<TProfile>(HIST("Prof_mu8_netproton"))->Fill(event_netproton.centrality(), pow(event_netproton.net_prot_no(), 8.0));

    // selecting subsample and filling profiles
    float l_Random = fRndm->Rndm();
    int SampleIndex = static_cast<int>(cfgNSubsample * l_Random);
    Subsample2D[SampleIndex][0]->Fill(event_netproton.centrality(), event_netproton.n_ch(), pow(event_netproton.net_prot_no(), 1.0));
    Subsample2D[SampleIndex][1]->Fill(event_netproton.centrality(), event_netproton.n_ch(), pow(event_netproton.net_prot_no(), 2.0));
    Subsample2D[SampleIndex][2]->Fill(event_netproton.centrality(), event_netproton.n_ch(), pow(event_netproton.net_prot_no(), 3.0));
    Subsample2D[SampleIndex][3]->Fill(event_netproton.centrality(), event_netproton.n_ch(), pow(event_netproton.net_prot_no(), 4.0));
    Subsample2D[SampleIndex][4]->Fill(event_netproton.centrality(), event_netproton.n_ch(), pow(event_netproton.net_prot_no(), 5.0));
    Subsample2D[SampleIndex][5]->Fill(event_netproton.centrality(), event_netproton.n_ch(), pow(event_netproton.net_prot_no(), 6.0));
    Subsample2D[SampleIndex][6]->Fill(event_netproton.centrality(), event_netproton.n_ch(), pow(event_netproton.net_prot_no(), 7.0));
    Subsample2D[SampleIndex][7]->Fill(event_netproton.centrality(), event_netproton.n_ch(), pow(event_netproton.net_prot_no(), 8.0));

    Subsample[SampleIndex][0]->Fill(event_netproton.centrality(), pow(event_netproton.net_prot_no(), 1.0));
    Subsample[SampleIndex][1]->Fill(event_netproton.centrality(), pow(event_netproton.net_prot_no(), 2.0));
    Subsample[SampleIndex][2]->Fill(event_netproton.centrality(), pow(event_netproton.net_prot_no(), 3.0));
    Subsample[SampleIndex][3]->Fill(event_netproton.centrality(), pow(event_netproton.net_prot_no(), 4.0));
    Subsample[SampleIndex][4]->Fill(event_netproton.centrality(), pow(event_netproton.net_prot_no(), 5.0));
    Subsample[SampleIndex][5]->Fill(event_netproton.centrality(), pow(event_netproton.net_prot_no(), 6.0));
    Subsample[SampleIndex][6]->Fill(event_netproton.centrality(), pow(event_netproton.net_prot_no(), 7.0));
    Subsample[SampleIndex][7]->Fill(event_netproton.centrality(), pow(event_netproton.net_prot_no(), 8.0));
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Equivalent to the AddTask in AliPhysics
  return WorkflowSpec{
    adaptAnalysisTask<NetProtonCumulants_Table_QA>(cfgc),
    adaptAnalysisTask<NetProtonCumulants_analysis>(cfgc),
  };
}
