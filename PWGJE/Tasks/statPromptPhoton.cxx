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

/// \file statPromptPhoton.cxx
/// \brief Reconstruction of Phi yield through track-track Minv correlations for resonance hadrochemistry analysis.
///
///
/// \author Adrian Fereydon Nassirpour <adrian.fereydon.nassirpour@cern.ch>

#include <TLorentzVector.h>
#include <TVector2.h>

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "CommonConstants/PhysicsConstants.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "EMCALBase/Geometry.h"
#include "EMCALCalib/BadChannelMap.h"

#include "DataFormatsEMCAL/Cell.h"
#include "DataFormatsEMCAL/Constants.h"
#include "DataFormatsEMCAL/AnalysisCluster.h"

#include "CommonDataFormat/InteractionRecord.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct statPromptPhoton {
  SliceCache cache;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<double> cfgMaxDCArToPVcut{"cfgMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  Configurable<double> cfgMaxDCAzToPVcut{"cfgMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};
  Configurable<bool> cfgConnectedToPV{"cfgConnectedToPV", true, "PV contributor track selection"};
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"};
  Configurable<double> cfgnFindableTPCClusters{"cfgnFindableTPCClusters", 50, "nFindable TPC Clusters"};
  Configurable<double> cfgnTPCCrossedRows{"cfgnTPCCrossedRows", 70, "nCrossed TPC Rows"};
  Configurable<double> cfgnRowsOverFindable{"cfgnRowsOverFindable", 1.2, "nRowsOverFindable TPC CLusters"};
  Configurable<double> cfgnTPCChi2{"cfgnTPChi2", 4.0, "nTPC Chi2 per Cluster"};
  Configurable<double> cfgnITSChi2{"cfgnITShi2", 36.0, "nITS Chi2 per Cluster"};
  Configurable<int> cfgClusterDefinition{"clusterDefinition", 10, "cluster definition to be selected, e.g. 10=kV3Default, 0 = kV1Default"};
  Configurable<float> cfgMinTime{"MinTime", -30., "Minimum cluster time for time cut"};
  Configurable<float> cfgMaxTime{"MaxTime", +35., "Maximum cluster time for time cut"};
  Configurable<float> cfgMinClusterEnergy{"MinClusterEnergy", 0.7f, "Minimal cluster energy"};
  Configurable<int> cfgMinNCells{"MinNCelss", 2, "Minimal amount of cells per cluster"};
  Configurable<int> cfgMaxNLM{"MaxNLM", 2, "Maximal amount of local Maxima per cluster"};
  Configurable<bool> cfgExoticContribution{"ExoticContribution", false, "Exotic cluster in the data"};
  Configurable<double> cfgtrkMinPt{"cfgtrkMinPt", 0.15, "set track min pT"};
  Configurable<double> cfgtrkMaxEta{"cfgtrkMaxEta", 0.6, "set track max Eta"};
  Configurable<float> cfgMinR{"MinR", 0.1, "Min. Radii of Delta R cone around photon trigger"};
  Configurable<float> cfgMaxR{"MaxR", 0.4, "Max. Radii of Delta R cone around photon trigger"};
  Configurable<float> cfgMinTrig{"MinTrig", 1, "Min. Trigger energy/momentum"};
  Configurable<float> cfgMaxTrig{"MaxTrig", 5, "Max. Trigger energy/momentum"};
  Configurable<float> cfgVtxCut{"cfgVtxCut", 10.0, "V_z cut selection"};

  // INIT
  void init(InitContext const&)
  {
    std::vector<double> ptBinning = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 12.0, 16.0, 20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0, 300.0, 500.0};
    AxisSpec pthadAxis = {ptBinning, "#it{p}_{T}^{had sum} [GeV/c]"};

    histos.add("REC_nEvents", "REC_nEvents", kTH1F, {{4, 0.0, 4.0}});
    histos.add("REC_Cluster_QA", "REC_Cluster_QA", kTH1F, {{10, -0.5, 9.5}});
    histos.add("REC_PtHadSum_Photon", "REC_PtHadSum_Photon", kTH1F, {pthadAxis});
    histos.add("REC_TrackPhi_photontrigger", "REC_TrackPhi_photontrigger", kTH1F, {{64, 0, 2 * TMath::Pi()}});
    histos.add("REC_TrackEta_photontrigger", "REC_TrackEta_photontrigger", kTH1F, {{100, -1, 1}});
    histos.add("REC_True_v_Cluster_Phi", "REC_True_v_Cluster_Phi", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
    histos.add("REC_True_v_Cluster_Eta", "REC_True_v_Cluster_Eta", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
    histos.add("REC_Track_v_Cluster_Phi", "REC_Track_v_Cluster_Phi", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
    histos.add("REC_Track_v_Cluster_Eta", "REC_Track_v_Cluster_Eta", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
    histos.add("REC_Track_v_Cluster_Phi_AC", "REC_Track_v_Cluster_Phi_AC", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
    histos.add("REC_Track_v_Cluster_Eta_AC", "REC_Track_v_Cluster_Eta_AC", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
    histos.add("REC_SumPt_BC", "REC_SumPt_BC", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
    histos.add("REC_SumPt_AC", "REC_SumPt_AC", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
    histos.add("REC_M02_BC", "REC_M02_BC", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});
    histos.add("REC_M02_AC", "REC_M02_AC", kTH1F, {{628, -2 * TMath::Pi(), 2 * TMath::Pi()}});

    histos.add("REC_Trigger_Purity", "REC_Trigger_Purity", kTH1F, {{4, 0.0, 4.0}});
    histos.add("REC_Trigger_Energy", "REC_Trigger_Energy", kTH1F, {{82, -1.0, 40.0}});

    histos.add("REC_Trigger_Energy_GOOD", "REC_Trigger_Energy_GOOD", kTH1F, {{82, -1.0, 40.0}});
    histos.add("REC_Trigger_Energy_MISS", "REC_Trigger_Energy_MISS", kTH1F, {{82, -1.0, 40.0}});
    histos.add("REC_Trigger_Energy_FAKE", "REC_Trigger_Energy_FAKE", kTH1F, {{82, -1.0, 40.0}});

    histos.add("REC_All_Energy", "REC_All_Energy", kTH1F, {{82, -1.0, 40.0}});
    histos.add("REC_True_Trigger_Energy", "REC_True_Trigger_Energy", kTH1F, {{82, -1.0, 40.0}});
    histos.add("REC_True_Prompt_Trigger_Energy", "REC_True_Prompt_Trigger_Energy", kTH1F, {{82, -1.0, 40.0}});

    histos.add("REC_Trigger_V_PtHadSum_Stern", "REC_Trigger_V_PtHadSum_Stern", kTH2F, {{100, 0, 100}, pthadAxis});
    histos.add("REC_Trigger_V_PtHadSum_Photon", "REC_Trigger_V_PtHadSum_Photon", kTH2F, {{100, 0, 100}, pthadAxis});
    histos.add("REC_TrueTrigger_V_PtHadSum_Photon", "REC_Trigger_V_PtHadSum_Photon", kTH2F, {{100, 0, 100}, pthadAxis});
    histos.add("REC_dR_Photon", "REC_dR_Photon", kTH1F, {{628, 0.0, 2 * TMath::Pi()}});
    histos.add("REC_dR_Stern", "REC_dR_Stern", kTH1F, {{628, 0.0, 2 * TMath::Pi()}});

    histos.add("GEN_nEvents", "GEN_nEvents", kTH1F, {{4, 0.0, 4.0}});
    histos.add("GEN_True_Photon_Energy", "GEN_True_Photon_Energy", kTH1F, {{8200, -1.0, 40.0}});
    histos.add("GEN_True_Prompt_Photon_Energy", "GEN_True_Prompt_Photon_Energy", kTH1F, {{8200, -1.0, 40.0}});
    histos.add("GEN_Trigger_V_PtHadSum_Stern", "GEN_Trigger_V_PtHadSum_Stern", kTH2F, {{100, 0, 100}, pthadAxis});
    histos.add("GEN_Trigger_V_PtHadSum_Photon", "GEN_Trigger_V_PtHadSum_Photon", kTH2F, {{100, 0, 100}, pthadAxis});
    histos.add("GEN_TrueTrigger_V_PtHadSum_Photon", "GEN_Trigger_V_PtHadSum_Photon", kTH2F, {{100, 0, 100}, pthadAxis});
    histos.add("GEN_dR_Photon", "GEN_dR_Photon", kTH1F, {{628, 0.0, 2 * TMath::Pi()}});
    histos.add("GEN_dR_Stern", "GEN_dR_Stern", kTH1F, {{628, 0.0, 2 * TMath::Pi()}});

  } // end of init

  Filter clusterDefinitionSelection = (o2::aod::emcalcluster::definition == cfgClusterDefinition) && (o2::aod::emcalcluster::time >= cfgMinTime) && (o2::aod::emcalcluster::time <= cfgMaxTime) && (o2::aod::emcalcluster::energy > cfgMinClusterEnergy) && (o2::aod::emcalcluster::nCells >= cfgMinNCells) && (o2::aod::emcalcluster::nlm <= cfgMaxNLM) && (o2::aod::emcalcluster::isExotic == cfgExoticContribution);
  Filter emccellfilter = aod::calo::caloType == 1; // mc emcal cell
  Filter PosZFilter = nabs(aod::collision::posZ) < cfgVtxCut;
  Filter mcPosZFilter = nabs(aod::mccollision::posZ) < cfgVtxCut;

  using MCCells = o2::soa::Join<aod::Calos, aod::McCaloLabels_001>;

  using MCClusters = o2::soa::Join<o2::aod::EMCALMCClusters, o2::aod::EMCALClusters>;
  using selectedCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using selectedMCCollisions = aod::McCollisions;

  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;

  using filteredMCCells = o2::soa::Filtered<MCCells>;
  using filteredMCClusters = soa::Filtered<MCClusters>;
  using filteredCollisions = soa::Filtered<selectedCollisions>;
  using filteredMCCollisions = soa::Filtered<selectedMCCollisions>;

  Preslice<o2::aod::EMCALMatchedTracks> perClusterMatchedTracks = o2::aod::emcalclustercell::emcalclusterId;

  // Helper functions
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <typename Tracks, typename Trigger>

  double GetPtHadSum(const Tracks& tracks, const Trigger& trigger, double MinR, double MaxR, bool IsStern, bool IsParticle, bool DodR)
  {
    double eta_trigger, phi_trigger;
    if constexpr (requires { trigger.eta(); }) {
      eta_trigger = trigger.eta();
      phi_trigger = trigger.phi();
    } else if constexpr (requires { trigger.Eta(); }) {
      eta_trigger = trigger.Eta();
      phi_trigger = trigger.Phi();
    }
    double pthadsum = 0;

    for (auto& track : tracks) {
      double phi_track = track.phi();
      double eta_track = track.eta();
      double pt_track = track.pt();

      if (!IsParticle) {
        if constexpr (requires { track.trackId(); }) {
          auto originaltrack = track.template track_as<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>>();
          if (!trackSelection(originaltrack)) {
            continue;
          } // reject track
        } else if constexpr (requires { track.sign(); }) { // checking for JTrack
          // done checking for JTrack, now default to normal tracks
          if (!trackSelection(track)) {
            continue;
          } // reject track
        } // done checking for JTrack
      } else {
        if constexpr (requires { track.isPhysicalPrimary(); }) {
          if (track.pt() < 0.15) {
            continue;
          }
          if (std::abs(track.eta()) > cfgtrkMaxEta) {
            continue;
          }
          if (track.getGenStatusCode() < 20) {
            continue;
          }
          if (!track.isPhysicalPrimary()) {
            continue;
          }
          int pdg = std::abs(track.pdgCode());
          if (pdg != 211 && pdg != 321 && pdg != 2212) {
            continue;
          }
        }
      }
      if (IsStern || IsParticle) {
        if constexpr (requires { trigger.globalIndex(); }) {
          if (trigger.globalIndex() == track.globalIndex())
            continue;
        }
      }
      double phidiff = TVector2::Phi_mpi_pi(phi_track - phi_trigger);
      double etadiff = std::abs(eta_track - eta_trigger);
      double dR = TMath::Sqrt((etadiff * etadiff) + (phidiff * phidiff));

      if (DodR) {
        if (dR > MinR && dR < MaxR) {
          if (!IsParticle) {
            if (IsStern) {
              histos.fill(HIST("REC_dR_Stern"), dR);
            }
            if (!IsStern) {
              histos.fill(HIST("REC_dR_Photon"), dR);
            }
          } else {
            if (IsStern) {
              histos.fill(HIST("GEN_dR_Stern"), dR);
            }
            if (!IsStern) {
              histos.fill(HIST("GEN_dR_Photon"), dR);
            }
          }
        }
      }
      if (dR > MinR && dR < MaxR) {
        pthadsum += pt_track;
      }

    } // end of track loop
    return pthadsum;
  } // end of GetPtHadSum
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <typename TrackType>
  bool trackSelection(const TrackType track)
  {
    // basic track cuts
    if (track.pt() < cfgtrkMinPt)
      return false;

    if (std::abs(track.eta()) > cfgtrkMaxEta)
      return false;

    if (std::abs(track.dcaXY()) > cfgMaxDCArToPVcut)
      return false;

    if (std::abs(track.dcaZ()) > cfgMaxDCAzToPVcut)
      return false;

    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;

    if (track.tpcNClsFindable() < cfgnFindableTPCClusters)
      return false;

    if (track.tpcNClsCrossedRows() < cfgnTPCCrossedRows)
      return false;

    if (track.tpcCrossedRowsOverFindableCls() > cfgnRowsOverFindable)
      return false;

    if (track.tpcChi2NCl() > cfgnTPCChi2)
      return false;

    if (track.itsChi2NCl() > cfgnITSChi2)
      return false;

    if (cfgConnectedToPV && !track.isPVContributor())
      return false;

    return true;
  }; // end of track selection
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////
  // PROCESS
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  int nEventsGenMC = 0;

  void processMCGen(filteredMCCollisions::iterator const& collision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, filteredCollisions>> const& recocolls, aod::McParticles const& mcParticles, filteredMCClusters const&)
  {
    nEventsGenMC++;
    if ((nEventsGenMC + 1) % 10000 == 0) {
      std::cout << "Processed Gen MC Events: " << nEventsGenMC << std::endl;
    }
    histos.fill(HIST("GEN_nEvents"), 0.5);
    if (fabs(collision.posZ()) > cfgVtxCut)
      return;

    if (recocolls.size() <= 0) // not reconstructed
      return;
    for (auto& recocoll : recocolls) { // poorly reconstructed
      if (!recocoll.sel8())
        return;
      if (fabs(recocoll.posZ()) > cfgVtxCut)

        return;
      histos.fill(HIST("GEN_nEvents"), 1.5);
      if (!recocoll.alias_bit(kTVXinEMC))
        return;
      histos.fill(HIST("GEN_nEvents"), 2.5);
    }

    for (auto& mcPhoton : mcParticles) {
      bool photontrigger = false;
      if (mcPhoton.pt() < 0.15)
        continue;
      if (std::abs(mcPhoton.eta()) > cfgtrkMaxEta)
        continue;
      if (mcPhoton.getGenStatusCode() < 20)
        continue;

      // first we check for pthadsums for all charged particles a la sternheimer
      if (mcPhoton.isPhysicalPrimary()) {
        int pdg = std::abs(mcPhoton.pdgCode());
        if (pdg == 211 || pdg == 321 || pdg == 2212) {
          bool sterntrigger = false;
          double sternPt = 0.0;
          if (mcPhoton.pt() > cfgMinTrig && mcPhoton.pt() < cfgMaxTrig) {
            if (fabs(mcPhoton.eta()) <= cfgtrkMaxEta) {
              sterntrigger = true;
              sternPt = mcPhoton.pt();
            }
          }
          // stern trigger
          if (sterntrigger) {
            bool doStern = true;
            double sterncount = 1.0;
            while (doStern) {
              TLorentzVector lParticleTrigger;
              lParticleTrigger.SetPxPyPzE(mcPhoton.px(), mcPhoton.py(), mcPhoton.pz(), mcPhoton.e());
              double pthadsum = GetPtHadSum(mcParticles, lParticleTrigger, cfgMinR, cfgMaxR, true, true, true);
              histos.fill(HIST("GEN_Trigger_V_PtHadSum_Stern"), sterncount, pthadsum, 2.0 / sternPt);
              if (sterncount < sternPt) {
                sterncount++;
              } else {
                doStern = false;
              }
            } // While sternin'
          } // stern trigger loop
        } // check if charged pikp
      } // check for primary particles

      // now we do all photons
      if (mcPhoton.pdgCode() == 22) {
        histos.fill(HIST("GEN_True_Photon_Energy"), mcPhoton.e());
        if (mcPhoton.pt() > cfgMinTrig && mcPhoton.pt() < cfgMaxTrig) {
          if (fabs(mcPhoton.eta()) <= cfgtrkMaxEta) {
            photontrigger = true;
          }
        } // check for photon trigger
        if (photontrigger) {
          TLorentzVector lRealPhoton;
          lRealPhoton.SetPxPyPzE(mcPhoton.px(), mcPhoton.py(), mcPhoton.pz(), mcPhoton.e());
          double truepthadsum = GetPtHadSum(mcParticles, lRealPhoton, cfgMinR, cfgMaxR, false, true, false);
          histos.fill(HIST("GEN_Trigger_V_PtHadSum_Photon"), mcPhoton.e(), truepthadsum);
        }
        // now we do all PROMPT photons
        if (std::abs(mcPhoton.getGenStatusCode()) > 19 && std::abs(mcPhoton.getGenStatusCode()) < 70) {
          if (mcPhoton.isPhysicalPrimary()) {
            histos.fill(HIST("GEN_True_Prompt_Photon_Energy"), mcPhoton.e());
            if (photontrigger) {
              TLorentzVector lRealPromptPhoton;
              lRealPromptPhoton.SetPxPyPzE(mcPhoton.px(), mcPhoton.py(), mcPhoton.pz(), mcPhoton.e());
              double truepthadsum = GetPtHadSum(mcParticles, lRealPromptPhoton, cfgMinR, cfgMaxR, false, true, true);
              histos.fill(HIST("GEN_TrueTrigger_V_PtHadSum_Photon"), mcPhoton.e(), truepthadsum);
            } // photontrigger
          } // check for primary photons
        } // prompt photon check

      } // photon check

    } // loop over mc particles

  } // end of process

  PROCESS_SWITCH(statPromptPhoton, processMCGen, "process MC Gen", true);

  Filter PosZFilter_JE = nabs(aod::jcollision::posZ) < cfgVtxCut;
  Filter clusterDefinitionSelection_JE = (o2::aod::jcluster::definition == cfgClusterDefinition) && (o2::aod::jcluster::time >= cfgMinTime) && (o2::aod::jcluster::time <= cfgMaxTime) && (o2::aod::jcluster::energy > cfgMinClusterEnergy) && (o2::aod::jcluster::nCells >= cfgMinNCells) && (o2::aod::jcluster::nlm <= cfgMaxNLM) && (o2::aod::jcluster::isExotic == cfgExoticContribution);

  using jTrackCandidates = soa::Join<aod::JTracks, aod::JTrackPIs, aod::McTrackLabels>;
  using jMCClusters = o2::soa::Join<o2::aod::JMcClusterLbs, o2::aod::JClusters, o2::aod::JClusterTracks>;
  using jselectedCollisions = soa::Join<aod::JCollisions, aod::EvSels, aod::JEMCCollisionLbs>;
  using jfilteredCollisions = soa::Filtered<jselectedCollisions>;
  using jfilteredMCClusters = soa::Filtered<jMCClusters>;

  int nEventsRecMC_JE = 0;
  void processMCRec_JE(jfilteredCollisions::iterator const& collision, jfilteredMCClusters const& mcclusters, jTrackCandidates const& tracks, soa::Join<aod::JTracks, aod::JTrackExtras, aod::JTrackPIs> const& /*caltracks*/, aod::JMcParticles const&, TrackCandidates const&)
  {

    nEventsRecMC_JE++;
    if ((nEventsRecMC_JE + 1) % 10000 == 0) {
      std::cout << "Processed JE Rec MC Events: " << nEventsRecMC_JE << std::endl;
    }

    histos.fill(HIST("REC_nEvents"), 0.5);

    if (fabs(collision.posZ()) > cfgVtxCut)
      return;
    if (!collision.sel8())
      return;

    histos.fill(HIST("REC_nEvents"), 1.5);

    if (!collision.isEmcalReadout())
      return;

    histos.fill(HIST("REC_nEvents"), 2.5);

    // now we do clusters
    bool clustertrigger = false;
    for (auto& mccluster : mcclusters) {
      histos.fill(HIST("REC_M02_BC"), mccluster.energy());
      if (mccluster.m02() < 0.1)
        continue;
      histos.fill(HIST("REC_Cluster_QA"), 0.5);
      histos.fill(HIST("REC_M02_AC"), mccluster.energy());
      histos.fill(HIST("REC_All_Energy"), mccluster.energy());
      bool photontrigger = false;  // is a neutral cluster
      bool vetotrigger = false;    // might be a neutral cluster
      bool chargetrigger = false;  // is definitely not a neutral cluster
      double photonPt = 0.0;
      double truephotonPt = 0.0;
      auto tracksofcluster = mccluster.matchedTracks_as<soa::Join<aod::JTracks, aod::JTrackExtras, aod::JTrackPIs>>();

      // first, we check if veto is required
      double sumptT = 0;
      double sumptTAC = 0;
      for (auto& ctrack : tracksofcluster) {
        auto ogtrack = ctrack.track_as<TrackCandidates>();
        if (!trackSelection(ogtrack)) {
          continue;
        }
        if (!ogtrack.isGlobalTrack()) {
          continue;
        }

        double ptT = ctrack.pt();
        sumptT += ptT;
        double phidiff = TVector2::Phi_mpi_pi(mccluster.phi() - ctrack.phi());
        double etadiff = mccluster.eta() - ctrack.eta();
        histos.fill(HIST("REC_Track_v_Cluster_Phi"), phidiff);
        histos.fill(HIST("REC_Track_v_Cluster_Eta"), etadiff);
      }

      if (sumptT > 0) {
        double mccluster_over_sumptT = mccluster.energy() / sumptT;
        histos.fill(HIST("REC_SumPt_BC"), mccluster_over_sumptT);
        if (mccluster_over_sumptT < 1.7) {
          histos.fill(HIST("REC_Cluster_QA"), 2.5);
          vetotrigger = true;
        } else {
          photontrigger = true;
          histos.fill(HIST("REC_Cluster_QA"), 1.5);
        }
      } else {
        photontrigger = true;
      }

      // veto is required
      if (vetotrigger) {
        for (auto& ctrack : tracksofcluster) {
          auto ogtrack = ctrack.track_as<TrackCandidates>();
          if (!trackSelection(ogtrack)) {
            continue;
          }
          if (!ogtrack.isGlobalTrack()) {
            continue;
          }

          double etaT = ctrack.eta();
          double etaC = mccluster.eta();
          double phiT = ctrack.phi();
          double phiC = mccluster.phi();
          double ptT = ctrack.pt();

          if (fabs(etaT - etaC) < (0.010 + pow(ptT + 4.07, -2.5))) {
            chargetrigger = true;
          }
          if (fabs(TVector2::Phi_mpi_pi(phiT - phiC)) < (0.015 + pow(ptT + 3.65, -2.0))) {
            chargetrigger = true;
          }
          if (chargetrigger) {
            histos.fill(HIST("REC_Cluster_QA"), 3.5);
            break;
          }
        } // tracks
      } // veto

      // check if cluster is good
      if (!chargetrigger) {
        for (auto& ctrack : tracksofcluster) {
          auto ogtrack = ctrack.track_as<TrackCandidates>();
          if (!trackSelection(ogtrack)) {
            continue;
          }
          if (!ogtrack.isGlobalTrack()) {
            continue;
          }

          double ptT = ctrack.pt();
          sumptTAC += ptT;
          double phidiff = TVector2::Phi_mpi_pi(mccluster.phi() - ctrack.phi());
          double etadiff = mccluster.eta() - ctrack.eta();
          histos.fill(HIST("REC_Track_v_Cluster_Phi_AC"), phidiff);
          histos.fill(HIST("REC_Track_v_Cluster_Eta_AC"), etadiff);
        } // tracks

        if (sumptTAC > 0) {
          double mccluster_over_sumptTAC = mccluster.energy() / sumptTAC;
          histos.fill(HIST("REC_SumPt_AC"), mccluster_over_sumptTAC);
        }

        photontrigger = true;
      } // check if there is no charge trigger

      if (photontrigger) { // if no charge trigger, cluster is good!
        histos.fill(HIST("REC_Cluster_QA"), 4.5);
        clustertrigger = true;
        double pthadsum = GetPtHadSum(tracks, mccluster, cfgMinR, cfgMaxR, false, false, true);
        histos.fill(HIST("REC_Trigger_V_PtHadSum_Photon"), photonPt, pthadsum);
        histos.fill(HIST("REC_PtHadSum_Photon"), pthadsum);
        histos.fill(HIST("REC_Trigger_Energy"), mccluster.energy());
      }

      // now we check the realness of our prompt photons
      auto ClusterParticles = mccluster.mcParticles_as<aod::JMcParticles>();
      bool goodgentrigger = true;
      for (auto& clusterparticle : ClusterParticles) {
        if (clusterparticle.pdgCode() == 211 || clusterparticle.pdgCode() == 321 || clusterparticle.pdgCode() == 2212) {
          goodgentrigger = false;
        }

        double phidiff = TVector2::Phi_mpi_pi(mccluster.phi() - clusterparticle.phi());
        double etadiff = mccluster.eta() - clusterparticle.eta();
        histos.fill(HIST("REC_True_v_Cluster_Phi"), phidiff);
        histos.fill(HIST("REC_True_v_Cluster_Eta"), etadiff);

        if (!photontrigger) {
          continue;
        }
        if (clusterparticle.pdgCode() == 22) {
          histos.fill(HIST("REC_True_Trigger_Energy"), clusterparticle.e());
          if (std::abs(clusterparticle.getGenStatusCode()) > 19 && std::abs(clusterparticle.getGenStatusCode()) < 70) {
            histos.fill(HIST("REC_True_Prompt_Trigger_Energy"), clusterparticle.e());
            TLorentzVector lRealPhoton;
            lRealPhoton.SetPxPyPzE(clusterparticle.px(), clusterparticle.py(), clusterparticle.pz(), clusterparticle.e());
            double truepthadsum = GetPtHadSum(tracks, lRealPhoton, cfgMinR, cfgMaxR, false, false, false);
            truephotonPt = clusterparticle.e();
            histos.fill(HIST("REC_TrueTrigger_V_PtHadSum_Photon"), truephotonPt, truepthadsum);
          }
        } // photon check
      } // clusterparticle loop
      if (goodgentrigger && photontrigger) {
        histos.fill(HIST("REC_Trigger_Purity"), 0.5);
        histos.fill(HIST("REC_Trigger_Energy_GOOD"), mccluster.energy());
      }
      if (goodgentrigger && !photontrigger) {
        histos.fill(HIST("REC_Trigger_Purity"), 1.5);
        histos.fill(HIST("REC_Trigger_Energy_MISS"), mccluster.energy());
      }
      if (!goodgentrigger && photontrigger) {
        histos.fill(HIST("REC_Trigger_Purity"), 2.5);
        histos.fill(HIST("REC_Trigger_Energy_FAKE"), mccluster.energy());
      }
    } // cluster loop

    // clusters done, now we do the sternheimer tracks
    for (auto& track : tracks) {
      bool sterntrigger = false;
      double sternPt = 0.0;
      auto ogtrack = track.track_as<TrackCandidates>();
      if (!trackSelection(ogtrack)) {
        continue;
      }
      if (clustertrigger) {
        histos.fill(HIST("REC_TrackPhi_photontrigger"), track.phi());
        histos.fill(HIST("REC_TrackEta_photontrigger"), track.eta());
      }
      if (track.pt() > cfgMinTrig && track.pt() < cfgMaxTrig) {
        if (fabs(track.eta()) <= cfgtrkMaxEta) {
          sterntrigger = true;
          sternPt = track.pt();
        }
      }

      if (sterntrigger) {
        bool doStern = true;
        double sterncount = 1.0;
        while (doStern) {
          double pthadsum = GetPtHadSum(tracks, track, cfgMinR, cfgMaxR, true, false, true);
          histos.fill(HIST("REC_Trigger_V_PtHadSum_Stern"), sterncount, pthadsum, 2.0 / sternPt);
          if (sterncount < sternPt) {
            sterncount++;
          } else {
            doStern = false;
          }
        } // While sternin'
      } // stern trigger loop
    } // track loop

  } // end of process

  PROCESS_SWITCH(statPromptPhoton, processMCRec_JE, "processJE  MC data", false);

}; // end of main struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<statPromptPhoton>(cfgc)};
};
