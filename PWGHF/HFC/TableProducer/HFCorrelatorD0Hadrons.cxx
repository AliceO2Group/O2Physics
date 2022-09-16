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

/// \file hfCorrelatorD0Hadrons.cxx
/// \brief D0-Hadron correlator task - data-like, MC-reco and MC-kine analyses.
///
/// \author Samrangy Sadhu <samrangy.sadhu@cern.ch>, INFN Bari

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/runDataProcessing.h"
#include <iostream>
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_prong2;
using namespace o2::aod::hf_correlation_d0hadron;
using namespace o2::analysis::hf_cuts_d0_topik;
using namespace o2::constants::math;



///
/// Returns deltaPhi value in range [-pi/2., 3.*pi/2], typically used for correlation studies
///
double getDeltaPhi(double phiD, double phiDbar)
{
  return RecoDecay::constrainAngle(phiDbar - phiD, -o2::constants::math::PI / 2.);
}


const int npTBinsMassAndEfficiency = o2::analysis::hf_cuts_d0_topik::npTBins;
const double efficiencyDmesonDefault[npTBinsMassAndEfficiency] = {};
auto efficiencyDmeson_v = std::vector<double>{efficiencyDmesonDefault, efficiencyDmesonDefault + npTBinsMassAndEfficiency};

// histogram binning definition
const int massAxisBins = 120;
const double massAxisMin = 1.5848;
const double massAxisMax = 2.1848;
const int phiAxisBins = 32;
const double phiAxisMin = 0.;
const double phiAxisMax = 2. * o2::constants::math::PI;
const int yAxisBins = 100;
const double yAxisMin = -5.;
const double yAxisMax = 5.;
const int ptDAxisBins = 180;
const double ptDAxisMin = 0.;
const double ptDAxisMax = 36.;

//using MCParticlesPlus = soa::Join<aod::McParticles, aod::HfCandProng2MCGen>;

struct hfCorrelatorD0Hadrons {
  Produces<aod::DHadronPair> entryD0HadronPair;
  Produces<aod::DHadronRecoInfo> entryD0HadronRecoInfo;

  HistogramRegistry registry{
    "registry",
    // NOTE: use hMassD0 for trigger normalisation (S*0.955), and hMass2DCorrelationPairs (in final task) for 2D-sideband-subtraction purposes
    {{"hPtCand", "D0,D0bar candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng0", "D0,D0bar candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng1", "D0,D0bar candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hSelectionStatus", "D0,D0bar candidates;selection status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}}},
     {"hEta", "D0,D0bar candidates;candidate #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhi", "D0,D0bar candidates;candidate #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hY", "D0,D0bar candidates;candidate #it{y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hMultiplicityPreSelection", "multiplicity prior to selection;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hMultiplicity", "multiplicity;multiplicity;entries", {HistType::kTH1F, {{10000, 0., 10000.}}}},
     {"hPtCandMCRec", "D0,D0bar candidates - MC reco;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng0MCRec", "D0,D0bar candidates - MC reco;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hPtProng1MCRec", "D0,D0bar candidates - MC reco;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hSelectionStatusMCRec", "D0,D0bar candidates - MC reco;selection status;entries", {HistType::kTH1F, {{4, -0.5, 3.5}}}},
     {"hEtaMCRec", "D0,D0bar candidates - MC reco;candidate #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhiMCRec", "D0,D0bar candidates - MC reco;candidate #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hYMCRec", "D0,D0bar candidates - MC reco;candidate #it{y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hMCEvtCount", "Event counter - MC gen;;entries", {HistType::kTH1F, {{1, -0.5, 0.5}}}},
     {"hPtCandMCGen", "D0,D0bar particles - MC gen;particle #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{ptDAxisBins, ptDAxisMin, ptDAxisMax}}}},
     {"hEtaMCGen", "D0,D0bar particles - MC gen;particle #it{#eta};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hPhiMCGen", "D0,D0bar particles - MC gen;particle #it{#varphi};entries", {HistType::kTH1F, {{phiAxisBins, phiAxisMin, phiAxisMax}}}},
     {"hYMCGen", "D0,D0bar candidates - MC gen;candidate #it{y};entries", {HistType::kTH1F, {{yAxisBins, yAxisMin, yAxisMax}}}},
     {"hCountSoftPiData", "soft pion counter -  Data", {HistType::kTH1F, {{20, 0., 20.}}}},
     {"hCountSoftPiMCRec", "soft pion counter - MC rec", {HistType::kTH1F, {{20, 0., 20.}}}},
     {"hCountSoftPiMCGen", "soft pion counter - MC gen", {HistType::kTH1F, {{20, 0., 20.}}}}}};
   

  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<double> cutYCandMax{"cutYCandMax", -1., "max. cand. rapidity"};
  Configurable<double> cutTrackEtaMax{"cutTrackEtaMax", 4., "max. eta of tracks"};
  Configurable<double> cutDCAxyMax{"cutDCAxyMax", 0.0025, "max. DCAxy of tracks"};
  Configurable<double> cutDCAzMax{"cutDCAzMax", 0.0025, "max. DCAz of tracks"};
  Configurable<double> cutPtCandMin{"cutPtCandMin", -1., "min. cand. pT"};
  Configurable<double> cutPtTrackMin{"cutPtTrackMin", -1., "min. track pT"};
  Configurable<double> cutPtCandMax{"cutPtCandMax", -1., "max. cand. pT"};
  Configurable<std::vector<double>> bins{"ptBinsForMassAndEfficiency", std::vector<double>{o2::analysis::hf_cuts_d0_topik::pTBins_v}, "pT bin limits for candidate mass plots and efficiency"};
  Configurable<std::vector<double>> efficiencyDmeson{"efficiencyDmeson", std::vector<double>{efficiencyDmeson_v}, "Efficiency values for D0 meson"};
  Configurable<int> flagApplyEfficiency{"efficiencyFlagD", 1, "Flag for applying D-meson efficiency weights"};
  Configurable<double> multMin{"multMin", 0., "minimum multiplicity accepted"};
  Configurable<double> multMax{"multMax", 10000., "maximum multiplicity accepted"};

  void init(o2::framework::InitContext&)
  {
    auto vbins = (std::vector<double>)bins;
    registry.add("hMass", "D0,D0bar candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMass1D", "D0,D0bar candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{massAxisBins, massAxisMin, massAxisMax}}});
    registry.add("hMassD01D", "D0,D0bar candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{massAxisBins, massAxisMin, massAxisMax}}});
    registry.add("hMassD0bar1D", "D0,D0bar candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{massAxisBins, massAxisMin, massAxisMax}}});    
   
    registry.add("hMassD0MCRecSig", "D0 signal candidates - MC reco;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

registry.add("hMassD0MCRecRefl", "D0 reflection candidates - MC reco;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0MCRecBkg", "D0 background candidates - MC reco;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0barMCRecSig", "D0bar signal candidates - MC reco;inv. mass D0bar only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0barMCRecRefl", "D0bar reflection candidates - MC reco;inv. mass D0bar only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassD0barMCRecBkg", "D0bar background candidates - MC reco;inv. mass D0bar only (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{massAxisBins, massAxisMin, massAxisMax}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCountD0triggersMCGen", "D0 trigger particles - MC gen;;N of trigger D0", {HistType::kTH2F, {{1, -0.5, 0.5}, {vbins, "#it{p}_{T} (GeV/#it{c})"}}});

  }

  Filter filterSelectCandidates = (aod::hf_selcandidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_selcandidate_d0::isSelD0bar >= selectionFlagD0bar);
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // =============================================================================  Process starts for Data ==================================================================================
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /// D0-h correlation pair builder - for real data and data-like analysis (i.e. reco-level w/o matching request via MC truth)
  void processData (aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksDCA>& tracks, soa::Filtered<soa::Join<aod::HfCandProng2, aod::HFSelD0Candidate>> const& candidates)
  {
    int nTracks = 0;
    if (collision.numContrib() > 1) {
      for (const auto& track : tracks) {
        if (std::abs(track.eta()) > cutTrackEtaMax) {
          continue;
        }
        if (abs(track.dcaXY()) > cutDCAxyMax || abs(track.dcaZ()) > cutDCAzMax) {
          continue;
        }
        nTracks++;
      }
    }
    registry.fill(HIST("hMultiplicityPreSelection"), nTracks);
    if (nTracks < multMin || nTracks > multMax) {
      return;
    }
    registry.fill(HIST("hMultiplicity"), nTracks);

    for (auto& candidate1 : candidates) {
      if (cutYCandMax >= 0. && std::abs(YD0(candidate1)) > cutYCandMax) {
        continue;
      }
      if (cutPtCandMin >= 0. && candidate1.pt() < cutPtCandMin) {
        continue;
      }
      // check decay channel flag for candidate1
      if (!(candidate1.hfflag() & 1 << DecayType::D0ToPiK)) {
        continue;
      }

	// ========================== Define parameters for soft pion removal ================================
	double p2prong0 = (candidate1.pxProng0()*candidate1.pxProng0())+(candidate1.pyProng0()*candidate1.pyProng0())+(candidate1.pzProng0()*candidate1.pzProng0());
	double p2prong1 = (candidate1.pxProng1()*candidate1.pxProng1())+(candidate1.pyProng1()*candidate1.pyProng1())+(candidate1.pzProng1()*candidate1.pzProng1());
	double mpi = RecoDecay::getMassPDG(kPiPlus);
	double mk = RecoDecay::getMassPDG(kKPlus);
	double e0pi = TMath::Sqrt(p2prong0+(mpi*mpi));
	double e1pi = TMath::Sqrt(p2prong1+(mpi*mpi));
	double e0k = TMath::Sqrt(p2prong0+(mk*mk));
	double e1k = TMath::Sqrt(p2prong1+(mk*mk));
	
	// ========================== trigger efficiency ================================
      double efficiencyWeight = 1.;
      if (flagApplyEfficiency) {
        efficiencyWeight = 1. / efficiencyDmeson->at(o2::analysis::findBin(bins, candidate1.pt()));
      }
	// ========================== Fill mass histo  ================================
      if (candidate1.isSelD0() >= selectionFlagD0) {
      registry.fill(HIST("hMass"), InvMassD0(candidate1), candidate1.pt(), efficiencyWeight);
      registry.fill(HIST("hMass1D"), InvMassD0(candidate1), efficiencyWeight);
      registry.fill(HIST("hMassD01D"), InvMassD0(candidate1), efficiencyWeight);
      } 
      if (candidate1.isSelD0bar() >= selectionFlagD0bar) {
      registry.fill(HIST("hMass"), InvMassD0bar(candidate1), candidate1.pt(), efficiencyWeight);
      registry.fill(HIST("hMass1D"), InvMassD0bar(candidate1), efficiencyWeight);
      registry.fill(HIST("hMassD0bar1D"), InvMassD0bar(candidate1), efficiencyWeight);
      } 
 	// ========================== Fill general histos ================================     
      registry.fill(HIST("hPtCand"), candidate1.pt());
      registry.fill(HIST("hPtProng0"), candidate1.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate1.ptProng1());
      registry.fill(HIST("hEta"), candidate1.eta());
      registry.fill(HIST("hPhi"), candidate1.phi());
      registry.fill(HIST("hY"), YD0(candidate1));
      registry.fill(HIST("hSelectionStatus"), candidate1.isSelD0bar() + (candidate1.isSelD0() * 2));
    
      // ================================================================================= D-h correlation dedicated section =====================================================

      
      	// ========================== track loop starts here ================================
      for (const auto& track : tracks) {
      //Remove D0 daughters by checking track indices
        if ((candidate1.index0Id() == track.mRowIndex) || (candidate1.index1Id() == track.mRowIndex)) {
          continue;
        }
        if (abs(track.dcaXY()) >= 1.0 || abs(track.dcaZ()) >= 1.0) continue; // Remove secondary tracks
        
                              registry.fill(HIST("hCountSoftPiData"),2);    //fill no. of tracks before soft pion removal    
        
        // ===== soft pion removal ===================================================
        double invmassDstar1 = 0, invmassDstar2 = 0;
        bool isSoftpiD0 = false, isSoftpiD0bar = false;
        double psum2 = ((candidate1.px()+track.px())*(candidate1.px()+track.px()))+((candidate1.py()+track.py())*(candidate1.py()+track.py()))+((candidate1.pz()+track.pz())*(candidate1.pz()+track.pz()));
        double Epion = TMath::Sqrt((track.px()*track.px())+(track.py()*track.py())+(track.pz()*track.pz())+(mpi*mpi));

        invmassDstar1 = TMath::Sqrt(pow(e0pi+e1k+Epion,2.)-psum2);
        invmassDstar2 = TMath::Sqrt(pow(e0k+e1pi+Epion,2.)-psum2);
        if (candidate1.isSelD0() >= selectionFlagD0){
        if ((TMath::Abs(invmassDstar1-InvMassD0(candidate1))-0.14543) < 3*800.*pow(10.,-6.)) {isSoftpiD0 = true; continue;}}
        
        if (candidate1.isSelD0bar() >= selectionFlagD0bar){
        if ((TMath::Abs(invmassDstar2-InvMassD0bar(candidate1))-0.14543) < 3*800.*pow(10.,-6.)) {isSoftpiD0bar = true; continue;}}
                              registry.fill(HIST("hCountSoftPiData"),4);    //fill no. of tracks after soft pion removal    
        
        int SignalStatus = 0;
        if ((candidate1.isSelD0() >= selectionFlagD0) && (isSoftpiD0 = false)) {SignalStatus += 1;}
        if ((candidate1.isSelD0bar() >= selectionFlagD0bar) && (isSoftpiD0bar = false)) {SignalStatus += 2;}           

          
        
        entryD0HadronPair(getDeltaPhi(track.phi(), candidate1.phi()),
                         track.eta() - candidate1.eta(),
                         candidate1.pt(),
                         track.pt());
        entryD0HadronRecoInfo(InvMassD0(candidate1),InvMassD0bar(candidate1),SignalStatus); 
          
      
      } // end inner loop (tracks)

    } // end outer loop
    
  }
  PROCESS_SWITCH(hfCorrelatorD0Hadrons, processData, "Process data", false);
  
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // =============================================================================  Process starts for MCRec ==================================================================================
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
    //Partition<soa::Join<aod::HfCandProng2, aod::HFSelD0Candidate, aod::HfCandProng2MCRec>> recoFlag2Prong = aod::hf_selcandidate_d0::isRecoHFFlag >= 0;
  void processMcRec(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksDCA>& tracks, soa::Filtered<soa::Join<aod::HfCandProng2, aod::HFSelD0Candidate, aod::HfCandProng2MCRec>> const& candidates)
  {
  
  int nTracks = 0;
    if (collision.numContrib() > 1) {
      for (const auto& track : tracks) {
        if (std::abs(track.eta()) > cutTrackEtaMax) {
          continue;
        }
        if (abs(track.dcaXY()) > cutDCAxyMax || abs(track.dcaZ()) > cutDCAzMax) {
          continue;
        }
        nTracks++;
      }
    }
    registry.fill(HIST("hMultiplicityPreSelection"), nTracks);
    if (nTracks < multMin || nTracks > multMax) {
      return;
    }
    registry.fill(HIST("hMultiplicity"), nTracks);
  // MC reco level
    bool flagD0 = false;
    bool flagD0bar = false;
    
    for (auto& candidate1 : candidates) {
      // check decay channel flag for candidate1
      if (!(candidate1.hfflag() & 1 << DecayType::D0ToPiK)) {
        continue;
      }
      if (cutYCandMax >= 0. && std::abs(YD0(candidate1)) > cutYCandMax) {
        continue;
      }
      if (cutPtCandMin >= 0. && candidate1.pt() < cutPtCandMin) {
        continue;
      }

      double efficiencyWeight = 1.;
      if (flagApplyEfficiency) {
        efficiencyWeight = 1. / efficiencyDmeson->at(o2::analysis::findBin(bins, candidate1.pt()));
      }
      
      if (std::abs(candidate1.flagMCMatchRec()) == 1 << DecayType::D0ToPiK) {
        // fill per-candidate distributions from D0/D0bar true candidates
        registry.fill(HIST("hPtCandMCRec"), candidate1.pt());
        registry.fill(HIST("hPtProng0MCRec"), candidate1.ptProng0());
        registry.fill(HIST("hPtProng1MCRec"), candidate1.ptProng1());
        registry.fill(HIST("hEtaMCRec"), candidate1.eta());
        registry.fill(HIST("hPhiMCRec"), candidate1.phi());
        registry.fill(HIST("hYMCRec"), YD0(candidate1));
        registry.fill(HIST("hSelectionStatusMCRec"), candidate1.isSelD0bar() + (candidate1.isSelD0() * 2));
      }
     // fill invariant mass plots from D0/D0bar signal and background candidates
      if (candidate1.isSelD0() >= selectionFlagD0) {                  // only reco as D0
        if (candidate1.flagMCMatchRec() == 1 << DecayType::D0ToPiK) { // also matched as D0
          registry.fill(HIST("hMassD0MCRecSig"), InvMassD0(candidate1), candidate1.pt(), efficiencyWeight);
        } else if (candidate1.flagMCMatchRec() == -(1 << DecayType::D0ToPiK)) {
          registry.fill(HIST("hMassD0MCRecRefl"), InvMassD0(candidate1), candidate1.pt(), efficiencyWeight);
        } else {
          registry.fill(HIST("hMassD0MCRecBkg"), InvMassD0(candidate1), candidate1.pt(), efficiencyWeight);
        }
      }
      if (candidate1.isSelD0bar() >= selectionFlagD0bar) {               // only reco as D0bar
        if (candidate1.flagMCMatchRec() == -(1 << DecayType::D0ToPiK)) { // also matched as D0bar
          registry.fill(HIST("hMassD0barMCRecSig"), InvMassD0bar(candidate1), candidate1.pt(), efficiencyWeight);
        } else if (candidate1.flagMCMatchRec() == 1 << DecayType::D0ToPiK) {
          registry.fill(HIST("hMassD0barMCRecRefl"), InvMassD0bar(candidate1), candidate1.pt(), efficiencyWeight);
        } else {
          registry.fill(HIST("hMassD0barMCRecBkg"), InvMassD0bar(candidate1), candidate1.pt(), efficiencyWeight);
        }
      } 
                      registry.fill(HIST("hCountSoftPiMCRec"),2);    //fill no. of tracks before soft pion removal    
      // ========================== Define parameters for soft pion removal ================================
	double p2prong0 = (candidate1.pxProng0()*candidate1.pxProng0())+(candidate1.pyProng0()*candidate1.pyProng0())+(candidate1.pzProng0()*candidate1.pzProng0());
	double p2prong1 = (candidate1.pxProng1()*candidate1.pxProng1())+(candidate1.pyProng1()*candidate1.pyProng1())+(candidate1.pzProng1()*candidate1.pzProng1());
	double mpi = RecoDecay::getMassPDG(kPiPlus);
	double mk = RecoDecay::getMassPDG(kKPlus);
	double e0pi = TMath::Sqrt(p2prong0+(mpi*mpi));
	double e1pi = TMath::Sqrt(p2prong1+(mpi*mpi));
	double e0k = TMath::Sqrt(p2prong0+(mk*mk));
	double e1k = TMath::Sqrt(p2prong1+(mk*mk));
      
      // D0-h correlation dedicated section
     
      flagD0 = candidate1.flagMCMatchRec() == (1 << DecayType::D0ToPiK);        // flagD0Signal 'true' if candidate1 matched to D0 (particle)
      flagD0bar = candidate1.flagMCMatchRec() == -(1 << DecayType::D0ToPiK); // flagD0Reflection 'true' if candidate1, selected as D0 (particle), is matched to D0bar (antiparticle)


      
      for (const auto& track : tracks) {
        if (std::abs(track.eta()) > cutTrackEtaMax) {
          continue;
        }
        if (track.pt() < cutPtTrackMin) {
          continue;
        }
        // Removing D0 daughters by checking track indices
        if ((candidate1.index0Id() == track.mRowIndex) || (candidate1.index1Id() == track.mRowIndex)) {
          continue;
        }
        if (std::abs(track.dcaXY()) >= 1.0 || std::abs(track.dcaZ()) >= 1.0) {
          continue; // Remove secondary tracks
        }
                registry.fill(HIST("hCountSoftPiMCRec"),2);    //fill no. of tracks before soft pion removal    
        // ===== soft pion removal ===================================================
        double invmassDstar1 = 0, invmassDstar2 = 0;
        bool isSoftpiD0 = false, isSoftpiD0bar = false;
        double psum2 = ((candidate1.px()+track.px())*(candidate1.px()+track.px()))+((candidate1.py()+track.py())*(candidate1.py()+track.py()))+((candidate1.pz()+track.pz())*(candidate1.pz()+track.pz()));
        double Epion = TMath::Sqrt((track.px()*track.px())+(track.py()*track.py())+(track.pz()*track.pz())+(mpi*mpi));

        invmassDstar1 = TMath::Sqrt(pow(e0pi+e1k+Epion,2.)-psum2);
        invmassDstar2 = TMath::Sqrt(pow(e0k+e1pi+Epion,2.)-psum2);

        if (candidate1.isSelD0() >= selectionFlagD0){
        if ((TMath::Abs(invmassDstar1-InvMassD0(candidate1))-0.14543) < 3*800.*pow(10.,-6.)) {
        isSoftpiD0 = true; 
        continue; }}
        
        if (candidate1.isSelD0bar() >= selectionFlagD0bar){
        if ((TMath::Abs(invmassDstar2-InvMassD0bar(candidate1))-0.14543) < 3*800.*pow(10.,-6.)) { 
        isSoftpiD0bar = true; 
        continue;} }
        

        registry.fill(HIST("hCountSoftPiMCRec"),4); // fill no. of tracks after soft pion removal
        
        int SignalStatus = 0;
        if ((flagD0 = true) && (candidate1.isSelD0() >= selectionFlagD0) && (isSoftpiD0 = false))      
                {SignalStatus += 1;} // signal case D0
        if ((flagD0bar = true) && (candidate1.isSelD0() >= selectionFlagD0) && (isSoftpiD0 = false))      
                {SignalStatus += 2;} // reflection case D0
        if ((flagD0 = false) && (flagD0bar = false) && (candidate1.isSelD0() >= selectionFlagD0) && (isSoftpiD0 = false))      
                {SignalStatus += 4;} // background case D0        
                
        if ((flagD0bar = true) && (candidate1.isSelD0bar() >= selectionFlagD0bar) && (isSoftpiD0bar = false))      
                {SignalStatus += 8;} //signal case D0bar     
        if ((flagD0 = true) && (candidate1.isSelD0bar() >= selectionFlagD0bar) && (isSoftpiD0bar = false))      
                {SignalStatus += 16;} //reflection case D0bar       
        if ((flagD0 = false) && (flagD0bar = false ) && (candidate1.isSelD0bar() >= selectionFlagD0bar) && (isSoftpiD0bar = false))      
                {SignalStatus += 32;} //background case D0bar  
                

        
        entryD0HadronPair(getDeltaPhi(track.phi(), candidate1.phi()),
                             track.eta() - candidate1.eta(),
                             candidate1.pt(),
                             track.pt());
        entryD0HadronRecoInfo(InvMassD0(candidate1),InvMassD0bar(candidate1),SignalStatus); 
      } // end inner loop (Tracks)
      
  } // end of outer loop (D0)
 }

  PROCESS_SWITCH(hfCorrelatorD0Hadrons, processMcRec, "Process MC Reco mode", true);
  
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // =============================================================================  Process starts for MCGen ==================================================================================
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  void processMcGen(aod::McCollision const& mccollision, soa::Join<aod::McParticles, aod::HfCandProng2MCGen> const& particlesMC)
  {

    registry.fill(HIST("hMCEvtCount"), 0);
    // MC gen level
    for (auto& particle1 : particlesMC) {
      // check if the particle is D0 or D0bar (for general plot filling and selection, so both cases are fine) - NOTE: decay channel is not probed!
      if (std::abs(particle1.pdgCode()) != pdg::Code::kD0) {
        continue;
      }
      double yD = RecoDecay::y(array{particle1.px(), particle1.py(), particle1.pz()}, RecoDecay::getMassPDG(particle1.pdgCode()));
      if (cutYCandMax >= 0. && std::abs(yD) > cutYCandMax) {
        continue;
      }
      if (cutPtCandMin >= 0. && particle1.pt() < cutPtCandMin) {
        continue;
      }
      registry.fill(HIST("hPtCandMCGen"), particle1.pt());
      registry.fill(HIST("hEtaMCGen"), particle1.eta());
      registry.fill(HIST("hPhiMCGen"), particle1.phi());
      registry.fill(HIST("hYMCGen"), yD);

// D-h correlation dedicated section

      if (std::abs(particle1.pdgCode()) != pdg::Code::kD0) // just checking the particle PDG, not the decay channel (differently from Reco: you have a BR factor btw such levels!)
	continue;
      registry.fill(HIST("hCountD0triggersMCGen"), 0, particle1.pt()); // to count trigger D0 (for normalisation)

      for (auto& particle2 : particlesMC) {

        if (std::abs(particle2.eta()) > cutTrackEtaMax) continue;

        if (particle2.pt() < cutPtTrackMin) continue;

        if ((std::abs(particle2.pdgCode()) != 11) && (std::abs(particle2.pdgCode()) != 13) && (std::abs(particle2.pdgCode()) != 211) && (std::abs(particle2.pdgCode()) != 321) && (std::abs(particle2.pdgCode()) != 2212))  continue;

        // ==============================soft pion removal================================
        registry.fill(HIST("hCountSoftPiMCGen"),2); //fill before soft pi removal
        // method used: indexMother = -1 by default if the mother doesn't match with given PID of the mother. We find mother of pion if it is D* and mother of D0 if it is D*. If they are both positive and they both match each other, then it is detected as a soft pion
        
        auto indexMotherPi = RecoDecay::getMother(particlesMC, particle2, 413, true, nullptr, 1); // last arguement 1 is written to consider immediate decay mother only
        auto indexMotherD0 = RecoDecay::getMother(particlesMC, particle1, 413, true, nullptr, 1);
        if(std::abs(particle2.pdgCode()) == 211 && indexMotherPi >= 0 && indexMotherD0 >= 0 && indexMotherPi==indexMotherD0)  continue;
        
        registry.fill(HIST("hCountSoftPiMCGen"),4); //fill after soft pion removal
        entryD0HadronPair(getDeltaPhi(particle2.phi(), particle1.phi()),
                             particle2.eta() - particle1.eta(),
                             particle1.pt(),
                             particle2.pt());
        entryD0HadronRecoInfo(1.864,
                             1.864,
                             0); //dummy info
      } // end inner loop

    }   // end outer loop
  }

  PROCESS_SWITCH(hfCorrelatorD0Hadrons, processMcGen, "Process MC Gen mode", false);


};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<hfCorrelatorD0Hadrons>(cfgc)};
}
