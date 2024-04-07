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
/// \brief  task for the analysis of rho photoproduction in UPCs
/// \author Jakub Juracka, jakub.juracka@cern.ch

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"

#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h" // has some useful funtions for stuff not available from the tables
#include "PWGUD/Core/DGPIDSelector.h" // possibly useful

// ROOT headers
//#include "TLorentzVector.h" // legacy class
#include <Math/Vector4D.h> // this should apparently be used instead of TLorentzVector, e.g. "ROOT::Math::PxPyPzMVector vector;"
//#include "TEfficiency.h" // for eventual MC studies

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using FullUDCollision = soa::Join<aod::UDCollisions, aod::UDCollisionsSels>::iterator;
using FullUDTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksPID, aod::UDTracksFlags>;

struct upcRhoAnalysis {

    // configurables
    double PcEtaCut = 0.90; // cut on track eta as per physics coordination recommendation
    Configurable<bool> tracksRequireTOF{"tracksRequireTOF", false, "requireTOF"};
    Configurable<double> tracksTpcNSigmaPiCut{"treacksTpcNSigmaPiCut", 3.0, "tpcNSigmaPiCut"};
    Configurable<double> tracksTofNSigmaPiCut{"treacksTofNSigmaPiCut", 3.0, "tofNSigmaPiCut"};
    Configurable<double> tracksPtMaxCut{"tracksPtMaxCut", 2.0, "ptMaxCut"};
    Configurable<double> tracksDcaMaxCut{"tracksDcaMaxCut", 1.0, "dcaMaxCut"};

    Configurable<double> systemYMaxCut{"systemYMaxCut", 0.9, "yMaxCut"};

    ConfigurableAxis mAxis{"mAxis", {250, 0.0, 2.5}, "m (GeV/#it{c}^{2})"};
    ConfigurableAxis ptAxis{"ptAxis", {300, 0.0, 3.0}, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis etaAxis{"etaAxis", {180, -0.9, 0.9}, "#eta"};
    ConfigurableAxis yAxis{"yAxis", {180, -0.9, 0.9}, "y"};
    ConfigurableAxis phiAxis{"phiAxis", {180, 0.0, 2.0*o2::constants::math::PI}, "#phi"};
    ConfigurableAxis nTracksAxis{"nTracksAxis", {101, -0.5, 100.5}, "N_{tracks}"};
    ConfigurableAxis tpcNSigmaPiAxis{"tpcNSigmaPiAxis", {400, -10.0, 30.0}, "TPC n#sigma_{#pi}"};
    ConfigurableAxis tofNSigmaPiAxis{"tofNSigmaPiAxis", {400, -20.0, 20.0}, "TOF n#sigma_{#pi}"};
    ConfigurableAxis dcaAxis{"dcaXYAxis", {1000, -5.0, 5.0}, "DCA (cm)"};

    HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

    void init(o2::framework::InitContext&){
        // selection counter
        std::vector<std::string> selectionNames = {"All", "PV", "Tracks", "TOF", "TPC", "DCA", "ITS"};
        int nSelections = selectionNames.size();
        registry.add("hSelectionCounter", ";;counts", kTH1D, {{nSelections, 0.5, nSelections + 0.5}});
        for (int i = 0; i < nSelections; i++) registry.get<TH1>(HIST("hSelectionCounter"))->GetXaxis()->SetBinLabel(i + 1, selectionNames[i].c_str());
        // collisions
        registry.add("QC/collisions/hNetCharge", ";net charge;counts", kTH1D, {{11, -5.5, 5.5}});
        registry.add("QC/collisions/hNumContributors", ";N_{contributors};counts", kTH1D, {{11, -0.5, 10.5}});
        registry.add("QC/collisions/hRgtrwTOF", ";fraction of tracks with TOF;counts", kTH1D, {{101, -0.005, 1.005}});
        registry.add("QC/collisions/hFt0Amplitude", ";A;C;counts", kTH2D, {{201, -1.5, 200.5}, {201, -1.5, 200.5}});
        registry.add("QC/collisions/hFt0Time", ";t_{A};t_{C};counts", kTH2D, {{1050, -999.5, 50.5}, {1050, -999.5, 50.5}});
        registry.add("QC/collisions/hFddAmplitude", ";A;C;counts", kTH2D, {{300, 0.0, 3000.0}, {300, 0.0, 3000.0}});
        registry.add("QC/collisions/hFddTime", ";t_{A};t_{C};counts", kTH2D, {{1050, -999.5, 50.5}, {1050, -999.5, 50.5}});
        registry.add("QC/collisions/hFv0Amplitude", ";A;counts", kTH1D, {{72, -1.5, 70.5}});
        registry.add("QC/collisions/hFv0Time", ";t;counts", kTH1D, {{1050, -999.5, 50.5}});
        // all tracks
        registry.add("QC/allTracks/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
        registry.add("QC/allTracks/hEta", ";#eta;counts", kTH1D, {etaAxis});
        registry.add("QC/allTracks/hPhi", ";#phi;counts", kTH1D, {phiAxis});
        registry.add("QC/allTracks/hCharge", ";charge;counts", kTH1D, {{5, -2.5, 2.5}});
        registry.add("QC/allTracks/hNTracks", ";N_{tracks};counts", kTH1D, {nTracksAxis});
        registry.add("QC/allTracks/hTpcNSigmaPi", ";TPC n#sigma_{#pi};counts", kTH1D, {tpcNSigmaPiAxis});
        registry.add("QC/allTracks/hTofNSigmaPi", ";TOF n#sigma_{#pi};counts", kTH1D, {tofNSigmaPiAxis});
        registry.add("QC/allTracks/hDcaXYZ", ";DCA_{z} (cm);DCA_{xy} (cm);counts", kTH2D, {dcaAxis, dcaAxis});
        registry.add("QC/allTracks/hTpcSignalVsPt", ";p_{T} (GeV/#it{c});TPC signal;counts", kTH2D, {ptAxis, {100, 0.0, 250.0}});
        registry.add("QC/allTracks/hTpcNsigmaPi2D", ";TPC n#sigma_{#pi}^{+};TPC n#sigma_{#pi}^{-};counts", kTH2D, {tpcNSigmaPiAxis, tpcNSigmaPiAxis});
        registry.add("QC/allTracks/hTpcNClsFindable", ";N_{findable};counts", kTH1D, {{171, -0.5, 170.5}});
        registry.add("QC/allTracks/hTpcNClsFindableMinusFound", ";N_{findable} - N_{found};counts", kTH1D, {{31, -10.5, 20.5}});
        // tracks passing selections
        registry.add("QC/cutTracks/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
        registry.add("QC/cutTracks/hEta", ";#eta;counts", kTH1D, {etaAxis});
        registry.add("QC/cutTracks/hPhi", ";#phi;counts", kTH1D, {phiAxis});
        registry.add("QC/cutTracks/hCharge", ";charge;counts", kTH1D, {{5, -2.5, 2.5}});
        registry.add("QC/cutTracks/hNTracks", ";N_{tracks};counts", kTH1D, {nTracksAxis});
        registry.add("QC/cutTracks/hTpcNSigmaPi", ";TPC n#sigma_{#pi};counts", kTH1D, {tpcNSigmaPiAxis});
        registry.add("QC/cutTracks/hTofNSigmaPi", ";TOF n#sigma_{#pi};counts", kTH1D, {tofNSigmaPiAxis});
        registry.add("QC/cutTracks/hDcaXYZ", ";DCA_{z} (cm);DCA_{xy} (cm);counts", kTH2D, {dcaAxis, dcaAxis});
        registry.add("QC/cutTracks/hTpcSignalVsPt", ";p_{T} (GeV/#it{c});TPC signal;counts", kTH2D, {ptAxis, {100, 0.0, 250.0}});
        registry.add("QC/cutTracks/hTpcNsigmaPi2D", ";TPC n#sigma_{#pi}^{+};TPC n#sigma_{#pi}^{-};counts", kTH2D, {tpcNSigmaPiAxis, tpcNSigmaPiAxis});
        registry.add("QC/cutTracks/hTpcNClsFindable", ";N_{findable};counts", kTH1D, {{171, -0.5, 170.5}});
        registry.add("QC/cutTracks/hTpcNClsFindableMinusFound", ";N_{findable} - N_{found};counts", kTH1D, {{31, -10.5, 20.5}});
        // reco pions
        registry.add("reco/pions/hPt", ";p_{T}(#pi^{+}) (GeV/#it{c});p_{T}(#pi^{-}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
        registry.add("reco/pions/hEta", ";#eta(#pi^{+});#eta(#pi^{-});counts", kTH2D, {etaAxis, etaAxis});
        registry.add("reco/pions/hPhi", ";#phi(#pi^{+});#phi(#pi^{-});counts", kTH2D, {phiAxis, phiAxis});
        // reco rhos
        registry.add("reco/rho/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
        registry.add("reco/rho/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
        registry.add("reco/rho/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
        registry.add("reco/rho/hY", ";y;counts", kTH1D, {yAxis});
        registry.add("reco/rho/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    }

    template <typename T>
    bool trackPassesCuts(T const& track) { // request ITS hits and check DCA
        if (!track.isPVContributor()) return false;
        if (track.dcaZ() > tracksDcaMaxCut) return false;
        if (track.dcaXY() > tracksDcaMaxCut) return false;
        if (!track.hasITS()) return false;
        if (tracksRequireTOF && !track.hasTOF()) return false;
        if (std::abs(eta(track.px(), track.py(), track.pz())) > PcEtaCut) return false;
        if (std::abs(track.pt()) > tracksPtMaxCut) return false;
        // if (std::abs(track.tpcNSigmaPi()) > tracksTpcNSigmaPiCut) return false; // will be done directly in the analysis
        // if (track.hasTOF() && std::abs(track.tofNSigmaPi()) > tracksTofNSigmaPiCut) return false; // will not use TOF (bias)
        else return true;
    }

    template <typename T>
    bool systemPassesCuts(T const& vec) {
        if (std::abs(vec.Rapidity()) > systemYMaxCut) return false;
        else return true;
    }

    void processReco(FullUDCollision const& collision, FullUDTracks const& tracks) {
        //if (collision.netCharge() != 0) return; // only consider events with net charge 0
        registry.fill(HIST("QC/collisions/hNetCharge"), collision.netCharge());
        registry.fill(HIST("QC/collisions/hNumContributors"), collision.numContrib());
        registry.fill(HIST("QC/collisions/hRgtrwTOF"), collision.rgtrwTOF());
        registry.fill(HIST("QC/collisions/hFt0Amplitude"), collision.totalFT0AmplitudeA(), collision.totalFT0AmplitudeC());
        registry.fill(HIST("QC/collisions/hFt0Time"), collision.timeFT0A(), collision.timeFT0C());
        registry.fill(HIST("QC/collisions/hFddAmplitude"), collision.totalFDDAmplitudeA(), collision.totalFDDAmplitudeC());
        registry.fill(HIST("QC/collisions/hFddTime"), collision.timeFDDA(), collision.timeFDDC());
        registry.fill(HIST("QC/collisions/hFv0Amplitude"), collision.totalFV0AmplitudeA());
        registry.fill(HIST("QC/collisions/hFv0Time"), collision.timeFV0A());

        // create some vectors for storing track info
        std::vector<ROOT::Math::PxPyPzMVector> track4Vecs;
        ROOT::Math::PxPyPzMVector rho;
        std::vector<decltype(tracks.begin())> allTracks, cutTracks;
        decltype(tracks.begin()) posTrack, negTrack;
        ROOT::Math::PxPyPzMVector piPos, piNeg;

        for (const auto& track : tracks) {
            allTracks.push_back(track);
            // fill histograms for all tracks
            registry.fill(HIST("QC/allTracks/hPt"), track.pt());
            registry.fill(HIST("QC/allTracks/hEta"), eta(track.px(), track.py(), track.pz()));
            registry.fill(HIST("QC/allTracks/hPhi"), phi(track.px(), track.py()));
            registry.fill(HIST("QC/allTracks/hCharge"), track.sign());
            registry.fill(HIST("QC/allTracks/hTpcNSigmaPi"), track.tpcNSigmaPi());
            registry.fill(HIST("QC/allTracks/hTofNSigmaPi"), track.tofNSigmaPi());
            registry.fill(HIST("QC/allTracks/hDcaXYZ"), track.dcaZ(), track.dcaXY());
            registry.fill(HIST("QC/allTracks/hTpcSignalVsPt"), track.pt(), track.tpcSignal());
            registry.fill(HIST("QC/allTracks/hTpcNClsFindable"), track.tpcNClsFindable());
            registry.fill(HIST("QC/allTracks/hTpcNClsFindableMinusFound"), track.tpcNClsFindableMinusFound());

            if (!trackPassesCuts(track)) continue;
            cutTracks.push_back(track);
            // fill histograms for cut tracks
            registry.fill(HIST("QC/cutTracks/hPt"), track.pt());
            registry.fill(HIST("QC/cutTracks/hEta"), eta(track.px(), track.py(), track.pz()));
            registry.fill(HIST("QC/cutTracks/hPhi"), phi(track.px(), track.py()));
            registry.fill(HIST("QC/cutTracks/hCharge"), track.sign());
            registry.fill(HIST("QC/cutTracks/hTpcNSigmaPi"), track.tpcNSigmaPi());
            registry.fill(HIST("QC/cutTracks/hTofNSigmaPi"), track.tofNSigmaPi());
            registry.fill(HIST("QC/cutTracks/hDcaXYZ"), track.dcaZ(), track.dcaXY());
            registry.fill(HIST("QC/cutTracks/hTpcSignalVsPt"), track.pt(), track.tpcSignal());
            registry.fill(HIST("QC/cutTracks/hTpcNClsFindable"), track.tpcNClsFindable());
            registry.fill(HIST("QC/cutTracks/hTpcNClsFindableMinusFound"), track.tpcNClsFindableMinusFound());

            track4Vecs.push_back(ROOT::Math::PxPyPzMVector(track.px(), track.py(), track.pz(), o2::constants::physics::MassPionCharged)); // assume pion mass
        }
        registry.fill(HIST("QC/allTracks/hNTracks"), allTracks.size());
        registry.fill(HIST("QC/cutTracks/hNTracks"), cutTracks.size());

        if (cutTracks.size() != 2 || cutTracks[0].sign() + cutTracks[1].sign() != 0) return; // only consider events with exactly 2 tracks of opposite charge// unlike-signed pion pair

        if (cutTracks[0].sign() > 0) {
            piPos = track4Vecs[0];
            posTrack = cutTracks[0];
            piNeg = track4Vecs[1];
            negTrack = cutTracks[1];
        } else {
            piNeg = track4Vecs[0];
            negTrack = cutTracks[0];
            piPos = track4Vecs[1];
            posTrack = cutTracks[1];
        }

        registry.fill(HIST("QC/allTracks/hTpcNsigmaPi2D"), posTrack.tpcNSigmaPi(), negTrack.tpcNSigmaPi());

        // 2D PID check
        if (std::pow(posTrack.tpcNSigmaPi(), 2) + std::pow(negTrack.tpcNSigmaPi(), 2) > std::pow(tracksTpcNSigmaPiCut, 2)) return;
        registry.fill(HIST("QC/cutTracks/hTpcNsigmaPi2D"), posTrack.tpcNSigmaPi(), negTrack.tpcNSigmaPi());
        registry.fill(HIST("reco/pions/hPt"), piPos.Pt(), piNeg.Pt());
        registry.fill(HIST("reco/pions/hEta"), piPos.Eta(), piNeg.Eta());
        registry.fill(HIST("reco/pions/hPhi"), piPos.Phi() + o2::constants::math::PI, piNeg.Phi() + o2::constants::math::PI); // shift by pi to get to the range [0, 2pi]
        // reconstruct rho
        rho = piPos + piNeg;
        registry.fill(HIST("reco/rho/hM"), rho.M());
        registry.fill(HIST("reco/rho/hPt"), rho.Pt());
        registry.fill(HIST("reco/rho/hPtVsM"), rho.M(), rho.Pt());
        registry.fill(HIST("reco/rho/hY"), rho.Rapidity());
        registry.fill(HIST("reco/rho/hPhi"), rho.Phi() + o2::constants::math::PI); // shift by pi to get to the range [0, 2pi]
    } PROCESS_SWITCH(upcRhoAnalysis, processReco, "Analyse reco tracks", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
    return WorkflowSpec{
        o2::framework::adaptAnalysisTask<upcRhoAnalysis>(cfgc)
    };
}