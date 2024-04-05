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

using FullUDCollision = aod::UDCollisions::iterator;
using FullUDTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksPID, aod::UDTracksFlags>;

struct upcRhoAnalysis {

    double PcEtaCut = 0.90; // cut on track eta as per physics coordination recommendation
    // configurables
    //Configurable<bool> verbosity{"verbosity", true, "verbosity"};
    Configurable<bool> tracksRequireTOF{"tracksRequireTOF", false, "requireTOF"};
    Configurable<double> tracksTpcNSigmaPiCut{"treacksTpcNSigmaPiCut", 3.0, "tpcNSigmaPiCut"};
    Configurable<double> tracksTofNSigmaPiCut{"treacksTofNSigmaPiCut", 3.0, "tofNSigmaPiCut"};
    Configurable<double> tracksPtMaxCut{"tracksPtMaxCut", 2.0, "ptMaxCut"};
    Configurable<double> tracksDcaMaxCut{"tracksDcaMaxCut", 1.0, "dcaMaxCut"};

    Configurable<double> systemYMaxCut{"systemYMaxCut", 0.9, "yMaxCut"};

    ConfigurableAxis mAxis{"mAxis", {250, 0.0, 2.5}, "m (GeV/#it{c}^{2})"};
    ConfigurableAxis ptAxis{"ptAxis", {250, 0.0, 2.5}, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis etaAxis{"etaAxis", {180, -0.9, 0.9}, "#eta"};
    ConfigurableAxis phiAxis{"phiAxis", {360, 0.0, 2.0*o2::constants::math::PI}, "#phi"};
    ConfigurableAxis nTracksAxis{"nTracksAxis", {101, -0.5, 100.5}, "N_{tracks}"};
    ConfigurableAxis tpcNSigmaPiAxis{"tpcNSigmaPiAxis", {500, -10.0, 10.0}, "TPC n#sigma_{#pi}"};
    ConfigurableAxis tofNSigmaPiAxis{"tofNSigmaPiAxis", {500, -10.0, 10.0}, "TOF n#sigma_{#pi}"};
    ConfigurableAxis dcaAxis{"dcaXYAxis", {100, -5.0, 5.0}, "DCA (cm)"};

    // the following is something stolen from UD_2b tutorial, not entirely sure if/how it works
    DGAnaparHolder anaPars = DGAnaparHolder(); // analysis cuts
    Configurable<DGAnaparHolder> DGPars{"anaPars", {}, "Analysis parameters"};
    DGPIDSelector pidsel = DGPIDSelector(); // PID selector

    HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

    void init(o2::framework::InitContext&){
        // all tracks
        registry.add("QC/allTracks/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
        registry.add("QC/allTracks/hEta", ";#eta;counts", kTH1D, {etaAxis});
        registry.add("QC/allTracks/hPhi", ";#phi;counts", kTH1D, {phiAxis});
        registry.add("QC/allTracks/hCharge", ";charge;counts", kTH1D, {{5, -2.5, 2.5}});
        // registry.add("QC/allTracks/hNTracks", ";N_{tracks};counts", kTH1D, {nTracksAxis});
        registry.add("QC/allTracks/hTpcNSigmaPi", ";TPC n#sigma_{#pi};counts", kTH1D, {tpcNSigmaPiAxis});
        registry.add("QC/allTracks/hTofNSigmaPi", ";TOF n#sigma_{#pi};counts", kTH1D, {tofNSigmaPiAxis});
        registry.add("QC/allTracks/hDcaXYZ", ";DCA_{z} (cm);DCA_{xy} (cm);counts", kTH2D, {dcaAxis, dcaAxis});
        registry.add("QC/allTracks/hTpcSignalVsPt", ";p_{T} (GeV/#it{c});TPC signal;counts", kTH2D, {ptAxis, {100, 0.0, 250.0}});
        registry.add("QC/allTracks/hTpcNsigmaPi2D", ";TPC n#sigma_{#pi}^{+};TPC n#sigma_{#pi}^{-};counts", kTH2D, {tpcNSigmaPiAxis, tpcNSigmaPiAxis});
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
        // reco pions
        registry.add("reco/pions/hNRecoPions", ";N_{#pi};counts", kTH1D, {{11, -0.5, 10.5}});
        registry.add("reco/pions/hPt", ";p_{T}(#pi^{+}) (GeV/#it{c});p_{T}(#pi^{-}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
        registry.add("reco/pions/hEta", ";#eta(#pi^{+});#eta(#pi^{-});counts", kTH2D, {etaAxis, etaAxis});
        registry.add("reco/pions/hPhi", ";#phi(#pi^{+});#phi(#pi^{-});counts", kTH2D, {phiAxis, phiAxis});
        // reco rhos
        registry.add("reco/rho/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
        registry.add("reco/rho/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
        registry.add("reco/rho/hEta", ";#eta;counts", kTH1D, {etaAxis});
        registry.add("reco/rho/hPhi", ";#phi;counts", kTH1D, {phiAxis});

        anaPars = (DGAnaparHolder)DGPars;
        pidsel.init(anaPars);
        LOGF(info, "<UDTutorial02>");
        pidsel.Print();
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
        // if (std::abs(track.tpcNSigmaPi()) > tracksTpcNSigmaPiCut) return false;
        // if (track.hasTOF() && std::abs(track.tofNSigmaPi()) > tracksTofNSigmaPiCut) return false;
        else return true;
    }

    template <typename T>
    bool systemPassesCuts(T const& vec) {
        if (std::abs(vec.Rapidity()) > systemYMaxCut) return false;
        else return true;
    }

    void processMine(FullUDCollision const& collision, FullUDTracks const& tracks) {
        //if (collision.netCharge() != 0) return; // only consider events with net charge 0

        // create some vectors for storing track info
        std::vector<ROOT::Math::PxPyPzMVector> piPos;
        std::vector<ROOT::Math::PxPyPzMVector> piNeg;
        std::vector<decltype(tracks.begin())> posTracks;
        std::vector<decltype(tracks.begin())> negTracks;
        ROOT::Math::PxPyPzMVector rho;
        double posTrackTpcNSigmaPi = 0.0, negTrackTpcNSigmaPi = 0.0, posCutTrackNSigmaPi = 0.0, negCutTrackNSigmaPi = 0.0;

        int nTracks = 0;
        for(const auto& track : tracks) {
            // fill histograms for all tracks
            registry.get<TH1>(HIST("QC/allTracks/hPt"))->Fill(track.pt());
            registry.get<TH1>(HIST("QC/allTracks/hEta"))->Fill(eta(track.px(), track.py(), track.pz()));
            registry.get<TH1>(HIST("QC/allTracks/hPhi"))->Fill(phi(track.px(), track.py()));
            registry.get<TH1>(HIST("QC/allTracks/hCharge"))->Fill(track.sign());
            registry.get<TH1>(HIST("QC/allTracks/hTpcNSigmaPi"))->Fill(track.tpcNSigmaPi());
            if (track.sign() > 0) posTrackTpcNSigmaPi = track.tpcNSigmaPi();
            if (track.sign() < 0) negTrackTpcNSigmaPi = track.tpcNSigmaPi();
            registry.get<TH1>(HIST("QC/allTracks/hTofNSigmaPi"))->Fill(track.tofNSigmaPi());
            registry.get<TH2>(HIST("QC/allTracks/hDcaXYZ"))->Fill(track.dcaZ(), track.dcaXY());
            registry.get<TH2>(HIST("QC/allTracks/hTpcSignalVsPt"))->Fill(track.pt(), track.tpcSignal());

            if (!trackPassesCuts(track)) continue;
            nTracks++;
            if (track.sign() > 0) posCutTrackNSigmaPi = track.tpcNSigmaPi();
            if (track.sign() < 0) negCutTrackNSigmaPi = track.tpcNSigmaPi();
            if (std::abs(track.tpcNSigmaPi()) > tracksTpcNSigmaPiCut) continue;
            // fill histograms for cut tracks
            registry.get<TH1>(HIST("QC/cutTracks/hPt"))->Fill(track.pt());
            registry.get<TH1>(HIST("QC/cutTracks/hEta"))->Fill(eta(track.px(), track.py(), track.pz()));
            registry.get<TH1>(HIST("QC/cutTracks/hPhi"))->Fill(phi(track.px(), track.py()));
            registry.get<TH1>(HIST("QC/cutTracks/hCharge"))->Fill(track.sign());
            registry.get<TH1>(HIST("QC/cutTracks/hTpcNSigmaPi"))->Fill(track.tpcNSigmaPi());
            registry.get<TH1>(HIST("QC/cutTracks/hTofNSigmaPi"))->Fill(track.tofNSigmaPi());
            registry.get<TH2>(HIST("QC/cutTracks/hDcaXYZ"))->Fill(track.dcaZ(), track.dcaXY());
            registry.get<TH2>(HIST("QC/cutTracks/hTpcSignalVsPt"))->Fill(track.pt(), track.tpcSignal());

            if (track.sign() == 1) piPos.push_back(ROOT::Math::PxPyPzMVector(track.px(), track.py(), track.pz(), o2::constants::physics::MassPionCharged)); // assume pion mass
            if (track.sign() == -1) piNeg.push_back(ROOT::Math::PxPyPzMVector(track.px(), track.py(), track.pz(), o2::constants::physics::MassPionCharged)); // assume pion mass
        }
        // registry.get<TH1>(HIST("QC/allTracks/hNTracks"))->Fill(tracks.size());
        registry.get<TH1>(HIST("QC/cutTracks/hNTracks"))->Fill(nTracks);


        registry.get<TH1>(HIST("reco/pions/hNRecoPions"))->Fill(piPos.size()+piNeg.size());

        if (piNeg.size() == 1 && piPos.size() == 1) {
            registry.get<TH2>(HIST("QC/allTracks/hTpcNsigmaPi2D"))->Fill(posTrackTpcNSigmaPi, negTrackTpcNSigmaPi);
            registry.get<TH2>(HIST("QC/cutTracks/hTpcNsigmaPi2D"))->Fill(posCutTrackNSigmaPi, negCutTrackNSigmaPi);

            rho = piPos[0] + piNeg[0];
            registry.get<TH2>(HIST("reco/pions/hPt"))->Fill(piPos[0].Pt(), piNeg[0].Pt());
            registry.get<TH2>(HIST("reco/pions/hEta"))->Fill(piPos[0].Rapidity(), piNeg[0].Rapidity());
            registry.get<TH2>(HIST("reco/pions/hPhi"))->Fill(piPos[0].Phi()+o2::constants::math::PI, piNeg[0].Phi()+o2::constants::math::PI); // shift by pi to get to the range [0, 2pi]

            registry.get<TH1>(HIST("reco/rho/hM"))->Fill(rho.M());
            registry.get<TH1>(HIST("reco/rho/hPt"))->Fill(rho.Pt());
            registry.get<TH1>(HIST("reco/rho/hEta"))->Fill(rho.Rapidity());
            registry.get<TH1>(HIST("reco/rho/hPhi"))->Fill(rho.Phi()+o2::constants::math::PI); // shift by pi to get to the range [0, 2pi]
        }
    } PROCESS_SWITCH(upcRhoAnalysis, processMine, "Analyse using sticks and stones technique", false);

    void processDG(FullUDCollision const& collision, FullUDTracks const& tracks) {
        Partition<FullUDTracks> PVContributors = aod::udtrack::isPVContributor == true;
        PVContributors.bindTable(tracks);
        if (collision.numContrib() < anaPars.minNTracks() || collision.numContrib() > anaPars.maxNTracks()) return;
        auto nIVMs = pidsel.computeIVMs(PVContributors);

        for (auto ivm : pidsel.unlikeIVMs()) {
            registry.get<TH1>(HIST("reco/rho/hM"))->Fill(ivm.M());
            registry.get<TH1>(HIST("reco/rho/hPt"))->Fill(ivm.Perp());
            for (auto ind : ivm.trkinds()) {
                auto track = PVContributors.begin() + ind;
                registry.get<TH1>(HIST("QC/cutTracks/hPt"))->Fill(track.pt());
                registry.get<TH1>(HIST("QC/cutTracks/hEta"))->Fill(eta(track.px(), track.py(), track.pz()));
                registry.get<TH1>(HIST("QC/cutTracks/hPhi"))->Fill(phi(track.px(), track.py()));
                registry.get<TH1>(HIST("QC/cutTracks/hCharge"))->Fill(track.sign());
                registry.get<TH1>(HIST("QC/cutTracks/hTpcNSigmaPi"))->Fill(track.tpcNSigmaPi());
                registry.get<TH1>(HIST("QC/cutTracks/hTofNSigmaPi"))->Fill(track.tofNSigmaPi());
                registry.get<TH2>(HIST("QC/cutTracks/hDcaXYZ"))->Fill(track.dcaZ(), track.dcaXY());
                registry.get<TH2>(HIST("QC/cutTracks/hTpcSignalVsPt"))->Fill(track.pt(), track.tpcSignal());
            }
        }
    } PROCESS_SWITCH(upcRhoAnalysis, processDG, "Analyse using Paul's technique", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
    return WorkflowSpec{
        o2::framework::adaptAnalysisTask<upcRhoAnalysis>(cfgc)
    };
}