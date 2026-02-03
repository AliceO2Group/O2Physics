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
// Task for analysing (anti)deuteron production in jets using pT-triggered data - update: 03-02-2026
//
// Executable : o2-analysis-lf-deuteron-in-jet-trg-pt

#include "PWGJE/Core/JetBkgSubUtils.h"
#include "PWGJE/Core/JetUtilities.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"

#include "Framework/AnalysisHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/PseudoJet.hh>

#include "CCDB/BasicCCDBManager.h"

#include <TRandom3.h>
#include <TVector2.h>
#include <TVector3.h>

#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::constants::math;
using namespace o2::constants::physics;

// Define conventional aliases for commoly used table joins
using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels>;
using SelectedTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTOFFullPr, aod::pidTOFFullDe>;


struct DeuteronInJetsTrgPt{
    // Histogram registry for data
    HistogramRegistry registryData{"registryData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

    // Random generator for subsample assignment
    TRandom3 mRand;

    // Setting default selection criteria to select tracks. May be changed when configuring the analysis.
    struct: o2::framework::ConfigurableGroup{
        std::string prefix{"cfgTrackCut"};
        // General specific
        Configurable<bool> requirePvContributor{"requirePvContributor", false, "Require that the track is a PV contributor"};
        Configurable<float> EtaMax{"EtaMax", 0.9f, "Max Eta for track acceptance"};
        Configurable<double> minPt{"minPt", 0.3, "Minimum pt of the tracks"};
        Configurable<double> maxDcaxy{"maxDcaxy", 0.05, "Maximum DCAxy"};
        Configurable<double> maxDcaz{"maxDcaz", 0.05, "Maximum DCAz"};
        // Part relative to ITS
        Configurable<int> ITSnClusMin{"ITSnClsMin", 6, "Minimum number of ITS clusters"};
        Configurable<float> ITSchi2ClusMax{"ITSchi2ClusMax", 36.f, "Max ITS Chi2 per cluster"};
        Configurable<bool> applyItsPid{"applyItsPid", false, "apply ITS PID"};
        Configurable<bool> setMCDefaultItsParams{"setMCDefaultItsParams", true, "Set MC default parameters for ITS PID"};
        Configurable<double> nSigmaItsMin{"nSigmaItsMin", -3.0, "nSigmaITS min"};
        Configurable<double> nSigmaItsMax{"nSigmaItsMax", +3.0, "nSigmaITS max"};
        Configurable<double> ptMaxItsPidProt{"ptMaxItsPidProt", 1.0, "maximum pt for ITS PID for protons"};
        Configurable<double> ptMaxItsPidDeut{"ptMaxItsPidDeut", 1.0, "maximum pt for ITS PID for deuterons"};
        // Part relative to TPC
        Configurable<int> TPCnClsMin{"TPCnClsMin", 110, "Minimum number of TPC clusters"};
        Configurable<float> TPCchi2ClusMin{"TPCchi2ClusMin", 0.f, "Min TPC Chi2 per cluster"};
        Configurable<float> TPCchi2ClusMax{"TPCchi2ClusMax", 4.f, "Max TPC Chi2 per cluster"};
        Configurable<int> TPCnCrossedRowsMin{"TPCnCrossedRowsMin", 100, "Minimum number of TPC crossed rows"};
        Configurable<double> Rtpc{"minRtpc", 0.8, "Minimum value of TPC crossed rows/TPC n cluster findable"};
        Configurable<float> TPCrigidityMin{"TPCrigidityMin", 0.3f, "Minimum TPC rigidity (p/Z) for track"};
        Configurable<double> minNsigmaTpc{"minNsigmaTpc", -3.0, "Minimum nsigma TPC"};
        Configurable<double> maxNsigmaTpc{"maxNsigmaTpc", +3.0, "Maximum nsigma TPC"};
        // Part relatuive to TOF
        Configurable<double> minNsigmaTof{"minNsigmaTof", -3.0, "Minimum nsigma TOF"};
        Configurable<double> maxNsigmaTof{"maxNsigmaTof", +3.5, "Maximum nsigma TOF"};
    } cfgTrackCut;

    // Setting default selection criteria for events. May be changes when configuring the analysis.
    struct: o2::framework::ConfigurableGroup{
        std::string prefix{"cgfEventCut"};
        Configurable<double> zVtx{"zVtx", 10.0, "Maximum z vertex"};
    } cfgEvCut;

    // Skimmed data flag and list of active triggers for processing
    Configurable<bool> cfgSkimmedProcessing{"cfgSkimmedProcessing", false, "Skimmed dataset processing"};
    Configurable<std::string> triggerList{"triggerList", "fJetFullLowPt", "Trigger list"};

    // Setting default selection criteria fr jet identification. May be changes when configuring the analysis.
    struct: o2::framework::ConfigurableGroup{
        std::string prefix{"cgfJetCut"};
        Configurable<double> minJetPt{"minJetPt", 10.0, "Minimum pt of the jet after bkg subtraction"};
        Configurable<double> rJet{"rJet", 0.3, "Jet parameter R"};
        Configurable<double> deltaEtaEdge{"deltaEtaEdge", 0.05, "eta gap from the edge"};
    } cfgJetCut;

    // Setting the number of bins and min and max value for the nsigma distribution
    Configurable<int> cfgNbins{"Nbins", 120, "Number of pT-bins"};
    Configurable<double> cfgpt_min{"pt_min", 0.0, "Min pT value of pT-axis"};
    Configurable<double> cfgpt_max{"pt_max", 6.0, "Max pT value of pT-axis"};

    // CCDB manager service for accessing condition data
    Service<o2::ccdb::BasicCCDBManager> ccdb;

    // Instantiate the main Zorro processing object and define an output to store summary information
    Zorro zorro;
    OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

    // Utility object for jet background subtraction methods
    JetBkgSubUtils backgroundSub;

    // Initiliaze ITS PID Rensponse object
    o2::aod::ITSResponse itsResponse;

    void initCCDB(aod::BCsWithTimestamps::iterator const& bc){
        if (cfgSkimmedProcessing){
            zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), triggerList);
            zorro.populateHistRegistry(registryData, bc.runNumber());
        }
    }

    // Set all configuration to allow task working
    void init(InitContext const&){
        // Set summary object if processing skimmed data
        if (cfgSkimmedProcessing) zorroSummary.setObject(zorro.getZorroSummary());

        // Set default MC parametrization for ITS response
        if (cfgTrackCut.setMCDefaultItsParams) itsResponse.setMCDefaultParameters();

        // Initialize random seed using high-resolution clock to ensure unique sequences across parallel Grid jobs
        auto time_seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        mRand.SetSeed(time_seed);

        // Histrograms for real data
        if (doprocessData){
            registryData.add("number_of_events_data", "number of events in data", HistType::kTH1F, {{4, 0, 4, "counter"}});   // Event counters
            registryData.add("settingData", "settingData", HistType::kTH2F, {{100, 0.0, 50.0, "min #it{p}^{jet}_{T} [GeV/#it{c}]"}, {20, 0.0, 1.0, "#it{R}_{jet}"}});   // Configuration
            registryData.add("jetEffectiveAreaOverPiR2", "jet effective area / piR^2", HistType::kTH1F, {{2000, 0, 2, "Area/#piR^{2}"}});   // Jet effective area over piR^2

            // Antiprotons
            registryData.add("antiproton_jet_tpc", "antiproton_jet_tpc", HistType::kTH2F, {{cfgNbins, cfgpt_min, cfgpt_max, "#it{p}_{T} [GeV/#it{c}]"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
            registryData.add("antiproton_jet_tof", "antiproton_jet_tof", HistType::kTH2F, {{cfgNbins, cfgpt_min, cfgpt_max, "#it{p}_{T} [GeV/#it{c}]"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});
            registryData.add("antiproton_ue_tpc", "antiproton_ue_tpc", HistType::kTH2F, {{cfgNbins, cfgpt_min, cfgpt_max, "#it{p}_{T} [GeV/#it{c}]"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
            registryData.add("antiproton_ue_tof", "antiproton_ue_tof", HistType::kTH2F, {{cfgNbins, cfgpt_min, cfgpt_max, "#it{p}_{T} [GeV/#it{c}]"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});
            registryData.add("antiproton_dca_jet", "antiproton_dca_jet", HistType::kTH2F, {{cfgNbins, cfgpt_min, cfgpt_max, "#it{p}_{T} [GeV/#it{c}]"}, {200, -1.0, 1.0, "DCA_{xy} [cm]"}});
            registryData.add("antiproton_dca_ue", "antiproton_dca_ue", HistType::kTH2F, {{cfgNbins, cfgpt_min, cfgpt_max, "#it{p}_{T} [GeV/#it{c}]"}, {200, -1.0, 1.0, "DCA_{xy} [cm]"}});

            // protons
            registryData.add("proton_jet_tpc", "proton_jet_tpc", HistType::kTH2F, {{cfgNbins, cfgpt_min, cfgpt_max, "#it{p}_{T} [GeV/#it{c}]"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
            registryData.add("proton_jet_tof", "proton_jet_tof", HistType::kTH2F, {{cfgNbins, cfgpt_min, cfgpt_max, "#it{p}_{T} [GeV/#it{c}]"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});
            registryData.add("proton_ue_tpc", "proton_ue_tpc", HistType::kTH2F, {{cfgNbins, cfgpt_min, cfgpt_max, "#it{p}_{T} [GeV/#it{c}]"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
            registryData.add("proton_ue_tof", "proton_ue_tof", HistType::kTH2F, {{cfgNbins, cfgpt_min, cfgpt_max, "#it{p}_{T} [GeV/#it{c}]"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});
            registryData.add("proton_dca_jet", "proton_dca_jet", HistType::kTH2F, {{cfgNbins, cfgpt_min, cfgpt_max, "#it{p}_{T} [GeV/#it{c}]"}, {200, -1.0, 1.0, "DCA_{xy} [cm]"}});
            registryData.add("proton_dca_ue", "proton_dca_ue", HistType::kTH2F, {{cfgNbins, cfgpt_min, cfgpt_max, "#it{p}_{T} [GeV/#it{c}]"}, {200, -1.0, 1.0, "DCA_{xy} [cm]"}});

            // Antideuterons
            registryData.add("antideuteron_jet_tpc", "antideuteron_jet_tpc", HistType::kTH2F, {{cfgNbins, 2 * cfgpt_min, 2 * cfgpt_max, "#it{p}_{T} [GeV/#it{c}]"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
            registryData.add("antideuteron_jet_tof", "antideuteron_jet_tof", HistType::kTH2F, {{cfgNbins, 2 * cfgpt_min, 2 * cfgpt_max, "#it{p}_{T} [GeV/#it{c}]"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});
            registryData.add("antideuteron_ue_tpc", "antideuteron_ue_tpc", HistType::kTH2F, {{cfgNbins, 2 * cfgpt_min, 2 * cfgpt_max, "#it{p}_{T} [GeV/#it{c}]"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
            registryData.add("antideuteron_ue_tof", "antideuteron_ue_tof", HistType::kTH2F, {{cfgNbins, 2 * cfgpt_min, 2 * cfgpt_max, "#it{p}_{T} [GeV/#it{c}]"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});

            // Deuterons
            registryData.add("deuteron_jet_tpc", "deuteron_jet_tpc", HistType::kTH2F, {{cfgNbins, 2 * cfgpt_min, 2 * cfgpt_max, "#it{p}_{T} [GeV/#it{c}]"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
            registryData.add("deuteron_jet_tof", "deuteron_jet_tof", HistType::kTH2F, {{cfgNbins, 2 * cfgpt_min, 2 * cfgpt_max, "#it{p}_{T} [GeV/#it{c}]"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});
            registryData.add("deuteron_ue_tpc", "deuteron_ue_tpc", HistType::kTH2F, {{cfgNbins, 2 * cfgpt_min, 2 * cfgpt_max, "#it{p}_{T} [GeV/#it{c}]"}, {400, -20.0, 20.0, "n#sigma_{TPC}"}});
            registryData.add("deuteron_ue_tof", "deuteron_ue_tof", HistType::kTH2F, {{cfgNbins, 2 * cfgpt_min, 2 * cfgpt_max, "#it{p}_{T} [GeV/#it{c}]"}, {400, -20.0, 20.0, "n#sigma_{TOF}"}});

            // nsigmaITS for antiproton candidates
            registryData.add("antiproton_nsigma_its_data", "antiproton_nsigma_its_data", HistType::kTH2F, {{cfgNbins, cfgpt_min, cfgpt_max, "#it{p}_{T} [GeV/#it{c}]"}, {400, -20.0, 20.0, "n#sigma_{ITS}"}});
        }
    }

    // Compute two transverse directions orthogonal to vector p (useful for UE evaluation)
    void getPerpendicularDirections(const TVector3& p, TVector3& u1, TVector3& u2){
        // Get momentum components
        double px = p.X();
        double py = p.Y();
        double pz = p.Z();

        // Precompute squared terms
        double px2 = px * px;
        double py2 = py * py;
        double pz2 = pz * pz;
        double pz4 = pz2 * pz2;

        // Case 1: vector along z-axis -> undefined perpendiculars
        if (px == 0 && py == 0) {
        u1.SetXYZ(0, 0, 0);
        u2.SetXYZ(0, 0, 0);
        return;
        }

        // Case 2: px = 0 -> avoid division by zero
        if (px == 0 && py != 0) {
        double ux = std::sqrt(py2 - pz4 / py2);
        double uy = -pz2 / py;
        u1.SetXYZ(ux, uy, pz);
        u2.SetXYZ(-ux, uy, pz);
        return;
        }

        // Case 3: py = 0 -> avoid division by zero
        if (py == 0 && px != 0) {
        double ux = -pz2 / px;
        double uy = std::sqrt(px2 - pz4 / px2);
        u1.SetXYZ(ux, uy, pz);
        u2.SetXYZ(ux, -uy, pz);
        return;
        }

        // General case: solve quadratic for perpendicular vectors
        double a = px2 + py2;
        double b = 2.0 * px * pz2;
        double c = pz4 - py2 * py2 - px2 * py2;
        double delta = b * b - 4.0 * a * c;

        // Invalid or degenerate solutions
        if (delta < 0 || a == 0) {
        u1.SetXYZ(0, 0, 0);
        u2.SetXYZ(0, 0, 0);
        return;
        }

        // Solution 1
        double u1x = (-b + std::sqrt(delta)) / (2.0 * a);
        double u1y = (-pz2 - px * u1x) / py;
        u1.SetXYZ(u1x, u1y, pz);

        // Solution 2
        double u2x = (-b - std::sqrt(delta)) / (2.0 * a);
        double u2y = (-pz2 - px * u2x) / py;
        u2.SetXYZ(u2x, u2y, pz);
    }

    // Compute delta phi
    double getDeltaPhi(double a1, double a2){
        double deltaPhi(0);
        double phi1 = TVector2::Phi_0_2pi(a1);
        double phi2 = TVector2::Phi_0_2pi(a2);
        double diff = std::fabs(phi1 - phi2);

        if (diff <= PI)
        deltaPhi = diff;
        if (diff > PI)
        deltaPhi = TwoPI - diff;

        return deltaPhi;
    }

    // Find ITS hit
    template <typename TrackIts>
    bool hasITShit(const TrackIts& track, int layer){
        int ibit = layer - 1;
        return (track.itsClusterMap() & (1 << ibit));
    }

    // Track selection criteria for creating jets
    template <typename JetTrack>
    bool passedTrackSelectionForJetReconstruction(const JetTrack& track){
        static constexpr int MinTpcCr = 70;
        static constexpr double MaxChi2Tpc = 4.0;
        static constexpr double MaxChi2Its = 36.0;
        static constexpr double MinPtTrack = 0.1;
        static constexpr double DcaxyMaxTrackPar0 = 0.0105;
        static constexpr double DcaxyMaxTrackPar1 = 0.035;
        static constexpr double DcaxyMaxTrackPar2 = 1.1;
        static constexpr double DcazMaxTrack = 2.0;

        // General part
        if (std::fabs(track.eta()) > cfgTrackCut.EtaMax) return false;
        if (track.pt() < MinPtTrack) return false;
        if (std::fabs(track.dcaXY()) > (DcaxyMaxTrackPar0 + DcaxyMaxTrackPar1 / std::pow(track.pt(), DcaxyMaxTrackPar2))) return false;     // DCAxy cut
        if (std::fabs(track.dcaZ()) > DcazMaxTrack) return false;       // DCAz cut
        // Part relative to ITS
        if (!track.hasITS()) return false;
        if ((!hasITShit(track, 1)) && (!hasITShit(track, 2)) && (!hasITShit(track, 3))) return false;       // Has Inner Barrel hit
        if (track.itsChi2NCl() >= MaxChi2Its) return false;
        // Part relative to TPC
        if (!track.hasTPC()) return false;
        if (track.tpcNClsCrossedRows() < MinTpcCr) return false;
        if (track.tpcChi2NCl() >= MaxChi2Tpc) return false;

        return true;
    }


    //Track selection for antinuclei
    template <typename AntinucleusTrack>
    bool passedTrackSelection(const AntinucleusTrack& track){

        // General part
        if (cfgTrackCut.requirePvContributor && !(track.isPVContributor())) return false;   // Flag to check if the track contributed to the collision vertex fit
        if (std::fabs(track.eta()) > cfgTrackCut.EtaMax) return false;     // Eta
        if (track.pt() < cfgTrackCut.minPt) return false;
        // Part relative to ITS
        if (!track.hasITS()) return false;                              // Flag to check if track has a ITS match
        if ((!hasITShit(track, 1)) && (!hasITShit(track, 2)) && (!hasITShit(track, 3))) return false;       // Require IB hit
        if (track.itsNCls() < cfgTrackCut.ITSnClusMin) return false;                        // Minimum number of ITS cluster
        if (track.itsChi2NCl() > cfgTrackCut.ITSchi2ClusMax) return false;                  // Minimum chi2 per cluster in ITS
        // Part relative to TPC
        if (!track.hasTPC()) return false;                                                   // Flag to check if track has a ITS match
        if (track.tpcNClsFound() < cfgTrackCut.TPCnClsMin) return false;                     // Minimum number of TPC cluster
        if (track.tpcNClsCrossedRows() < cfgTrackCut.TPCnCrossedRowsMin) return false;       // Minimum number of crossed rows in TPC
        if (track.tpcChi2NCl() < cfgTrackCut.TPCchi2ClusMin) return false;                   // Minimum chi2 per cluster in TPC
        if (track.tpcChi2NCl() > cfgTrackCut.TPCchi2ClusMax) return false;                   // Maximum chi2 per cluster in TPC
        if (track.tpcCrossedRowsOverFindableCls() < cfgTrackCut.Rtpc) return false;          // R_{TPC} > 0.8

        return true;
    }


    void processData(SelectedCollisions::iterator const& collision, SelectedTracks const& tracks, aod::BCsWithTimestamps const&){
        // Event counter before event selection
        registryData.fill(HIST("number_of_events_data"),0.5);
        registryData.fill(HIST("settingData"), cfgJetCut.minJetPt.value, cfgJetCut.rJet.value);

        // Retrieve the bunch crossing information with timestamps from the collision
        auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
        initCCDB(bc);

        // If skimmed processing is enabled, aplly Zorro trigger selection
        if (cfgSkimmedProcessing && !zorro.isSelected(collision.template bc_as<aod::BCsWithTimestamps>().globalBC())) return;
        registryData.fill(HIST("number_of_events_data"), 1.5);

        // Apply standard event selection
        if (!collision.sel8() || std::fabs(collision.posZ()) >= cfgEvCut.zVtx) return;
        registryData.fill(HIST("number_of_events_data"), 2.5);  // Save number of collisions that passed standard selections

        // loop over reco tracks
        int id(-1);
        std::vector<fastjet::PseudoJet> fjParticles;
        for (auto const& track : tracks){
            id++;
            if (!passedTrackSelectionForJetReconstruction(track)) continue;     // Skip tracks that not satisfy tracking selection criteria

            // 4-momentum representation of a particle
            fastjet::PseudoJet fourMomentum(track.px(),track.py(), track.pz(), track.energy(MassPionCharged));
            fourMomentum.set_user_index(id);
            fjParticles.emplace_back(fourMomentum);
        }

        if (fjParticles.empty()) return;    // Reject empty events

        // Cluster particles using the the anti-kT algorithm
        fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, cfgJetCut.rJet);     // Defining the algorithm to cluster, and the jet radius
        fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0));     // Activate area evaluation, and set area evaluation method
        fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);      // Declare the will of applying clustering algorithm with area evaluation to the selected candidates
        std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
        auto [rhoPerp, rhoMPerp] = jetutilities::estimateRhoPerpCone(fjParticles, jets[0], cfgJetCut.rJet);

        // Loop over reco jets
        bool isAtLeastOneJetSelected = false;
        for (auto const& jet : jets){

            if ((std::fabs(jet.eta()) + cfgJetCut.rJet) > (cfgTrackCut.EtaMax - cfgJetCut.deltaEtaEdge)) continue;  // Jet must be fully contained in the acceptance

            // Jet pt must be larger than threshold
            auto jetForSub = jet;
            fastjet::PseudoJet jetMinusBkg = backgroundSub.doRhoAreaSub(jetForSub, rhoPerp, rhoMPerp);
            if (jetMinusBkg.pt() < cfgJetCut.minJetPt) continue;      // Skip jets with pT < pT threshold

            double normalizedJetArea = jet.area() / (PI * cfgJetCut.rJet * cfgJetCut.rJet);
            isAtLeastOneJetSelected = true;

            // Perpendicular cones
            double coneRadius = std::sqrt(jet.area() / PI);
            TVector3 jetAxis(jet.px(), jet.py(), jet.pz());
            TVector3 ueAxis1(0, 0, 0), ueAxis2(0, 0, 0);
            getPerpendicularDirections(jetAxis, ueAxis1, ueAxis2);
            if (ueAxis1.Mag() == 0 || ueAxis2.Mag() == 0) continue;   // Skip not valid orthogonal cones

            registryData.fill(HIST("jetEffectiveAreaOverPiR2"), normalizedJetArea);   // Fill histogram with jet effective area / piR^2
            std::vector<fastjet::PseudoJet> jetConstituents = jet.constituents();       // Get jet constituents

            // Loop over constituents
            for (const auto& particle : jetConstituents){
                auto const& track = tracks.iteratorAt(particle.user_index());       // Get the corresponding track
                if (!passedTrackSelection(track)) continue;     // Apply track selection criteria

                // Define variables
                double nsigmaTPCPr = track.tpcNSigmaPr();
                double nsigmaTOFPr = track.tofNSigmaPr();
                double nsigmaTPCDe = track.tpcNSigmaDe();
                double nsigmaTOFDe = track.tofNSigmaDe();
                double pt = track.pt();
                double dcaxy = track.dcaXY();
                double dcaz = track.dcaZ();

                // Fill DCA distribution for (anti)protons
                if (track.sign() < 0 && std::fabs(dcaz) < cfgTrackCut.maxDcaz) registryData.fill(HIST("antiproton_dca_jet"), pt, dcaxy);
                if (track.sign() > 0 && std::fabs(dcaz) < cfgTrackCut.maxDcaz) registryData.fill(HIST("proton_dca_jet"), pt, dcaxy);

                // Apply DCA selections
                if (std::fabs(dcaxy) > cfgTrackCut.maxDcaxy || std::fabs(dcaz) > cfgTrackCut.maxDcaz) continue;

                // Particle identification using the ITS cluster size
                bool passedItsPidProt(true), passedItsPidDeut(true);
                double nSigmaITSprot = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Proton>(track));
                double nSigmaITSdeut = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Deuteron>(track));

                if (cfgTrackCut.applyItsPid && pt < cfgTrackCut.ptMaxItsPidProt && (nSigmaITSprot < cfgTrackCut.nSigmaItsMin || nSigmaITSprot > cfgTrackCut.nSigmaItsMax)) passedItsPidProt = false;
                if (cfgTrackCut.applyItsPid && pt < cfgTrackCut.ptMaxItsPidDeut && (nSigmaITSdeut < cfgTrackCut.nSigmaItsMin || nSigmaITSdeut > cfgTrackCut.nSigmaItsMax)) passedItsPidDeut = false;

                // Fill histograms for antimatter
                if (track.sign()<0){
                    if (passedItsPidProt){
                        registryData.fill(HIST("antiproton_jet_tpc"), pt, nsigmaTPCPr);
                        if (nsigmaTPCPr > cfgTrackCut.minNsigmaTpc && nsigmaTPCPr < cfgTrackCut.maxNsigmaTpc && track.hasTOF()) registryData.fill(HIST("antiproton_jet_tof"), pt, nsigmaTOFPr);        // requiring that track candidate in TOF have nisgma in TPC < threshold
                    }
                    if (passedItsPidDeut){
                        registryData.fill(HIST("antideuteron_jet_tpc"), pt, nsigmaTPCDe);
                        if (nsigmaTPCDe > cfgTrackCut.minNsigmaTpc && nsigmaTPCDe < cfgTrackCut.maxNsigmaTpc && track.hasTOF()) registryData.fill(HIST("antideuteron_jet_tof"), pt, nsigmaTOFDe);        // requiring that track candidate in TOF have nisgma in TPC < threshold
                    }
                }
                // Fill histograms for matter
                if (track.sign()>0){
                    if (passedItsPidProt){
                        registryData.fill(HIST("proton_jet_tpc"), pt, nsigmaTPCPr);
                        if (nsigmaTPCPr > cfgTrackCut.minNsigmaTpc && nsigmaTPCPr < cfgTrackCut.maxNsigmaTpc && track.hasTOF()) registryData.fill(HIST("proton_jet_tof"), pt, nsigmaTOFPr);        // requiring that track candidate in TOF have nisgma in TPC < threshold
                    }
                    if (passedItsPidDeut){
                        registryData.fill(HIST("deuteron_jet_tpc"), pt, nsigmaTPCDe);
                        if (nsigmaTPCDe > cfgTrackCut.minNsigmaTpc && nsigmaTPCDe < cfgTrackCut.maxNsigmaTpc && track.hasTOF()) registryData.fill(HIST("deuteron_jet_tof"), pt, nsigmaTOFDe);        // requiring that track candidate in TOF have nisgma in TPC < threshold
                    }
                }
            }   // End of loop over jet constituents

            // Loop over tracks in the UE
            for (auto const& track : tracks){
                if (!passedTrackSelection(track)) continue; // Apply track selection criteria

                // Calculate the angular distance between the track and the UE axes in eta-phi space
                double deltaEtaUe1 = track.eta() - ueAxis1.Eta();
                double deltaPhiUe1 = getDeltaPhi(track.phi(), ueAxis1.Phi());
                double deltaRUe1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
                double deltaEtaUe2 = track.eta() - ueAxis2.Eta();
                double deltaPhiUe2 = getDeltaPhi(track.phi(), ueAxis2.Phi());
                double deltaRUe2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

                double maxConeRadius = coneRadius;
                if (deltaRUe1 > maxConeRadius && deltaRUe2 > maxConeRadius) continue;       // Reject tracks that lie outside the maxConeRadius from both UE axes

                // Define variables
                double nsigmaTPCPr = track.tpcNSigmaPr();
                double nsigmaTOFPr = track.tofNSigmaPr();
                double nsigmaTPCDe = track.tpcNSigmaDe();
                double nsigmaTOFDe = track.tofNSigmaDe();
                double pt = track.pt();
                double dcaxy = track.dcaXY();
                double dcaz = track.dcaZ();

                // Fill DCA distribution for (anti)protons
                if (track.sign() < 0 && std::fabs(dcaz) < cfgTrackCut.maxDcaz) registryData.fill(HIST("antiproton_dca_ue"), pt, dcaxy);
                if (track.sign() > 0 && std::fabs(dcaz) < cfgTrackCut.maxDcaz) registryData.fill(HIST("proton_dca_ue"), pt, dcaxy);

                // Apply DCA selections
                if (std::fabs(dcaxy) > cfgTrackCut.maxDcaxy || std::fabs(dcaz) > cfgTrackCut.maxDcaz) continue;

                // Particle identification using the ITS cluster size
                bool passedItsPidProt(true), passedItsPidDeut(true);
                double nSigmaITSprot = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Proton>(track));
                double nSigmaITSdeut = static_cast<double>(itsResponse.nSigmaITS<o2::track::PID::Deuteron>(track));

                if (cfgTrackCut.applyItsPid && pt < cfgTrackCut.ptMaxItsPidProt && (nSigmaITSprot < cfgTrackCut.nSigmaItsMin || nSigmaITSprot > cfgTrackCut.nSigmaItsMax)) passedItsPidProt = false;
                if (cfgTrackCut.applyItsPid && pt < cfgTrackCut.ptMaxItsPidDeut && (nSigmaITSdeut < cfgTrackCut.nSigmaItsMin || nSigmaITSdeut > cfgTrackCut.nSigmaItsMax)) passedItsPidDeut = false;

                // Fill histograms for antimatter
                if (track.sign()<0){
                    if (passedItsPidProt){
                        registryData.fill(HIST("antiproton_ue_tpc"), pt, nsigmaTPCPr);
                        if (nsigmaTPCPr > cfgTrackCut.minNsigmaTpc && nsigmaTPCPr < cfgTrackCut.maxNsigmaTpc && track.hasTOF()) registryData.fill(HIST("antiproton_ue_tof"), pt, nsigmaTOFPr);        // requiring that track candidate in TOF have nisgma in TPC < threshold
                    }
                    if (passedItsPidDeut){
                        registryData.fill(HIST("antideuteron_ue_tpc"), pt, nsigmaTPCDe);
                        if (nsigmaTPCDe > cfgTrackCut.minNsigmaTpc && nsigmaTPCDe < cfgTrackCut.maxNsigmaTpc && track.hasTOF()) registryData.fill(HIST("antideuteron_ue_tof"), pt, nsigmaTOFDe);        // requiring that track candidate in TOF have nisgma in TPC < threshold
                    }
                }
                // Fill histograms for matter
                if (track.sign()>0){
                    if (passedItsPidProt){
                        registryData.fill(HIST("proton_ue_tpc"), pt, nsigmaTPCPr);
                        if (nsigmaTPCPr > cfgTrackCut.minNsigmaTpc && nsigmaTPCPr < cfgTrackCut.maxNsigmaTpc && track.hasTOF()) registryData.fill(HIST("proton_ue_tof"), pt, nsigmaTOFPr);        // requiring that track candidate in TOF have nisgma in TPC < threshold
                    }
                    if (passedItsPidDeut){
                        registryData.fill(HIST("deuteron_ue_tpc"), pt, nsigmaTPCDe);
                        if (nsigmaTPCDe > cfgTrackCut.minNsigmaTpc && nsigmaTPCDe < cfgTrackCut.maxNsigmaTpc && track.hasTOF()) registryData.fill(HIST("deuteron_ue_tof"), pt, nsigmaTOFDe);        // requiring that track candidate in TOF have nisgma in TPC < threshold
                    }
                }
            }
        }

        // Event counter: events with at least one jet selected
        if (isAtLeastOneJetSelected) registryData.fill(HIST("number_of_events_data"), 3.5);
    }
    PROCESS_SWITCH(DeuteronInJetsTrgPt, processData, "Process Data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc){
    return WorkflowSpec{adaptAnalysisTask<DeuteronInJetsTrgPt>(cfgc)};
}
