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
/// \date   25.3.2024

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"

#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h" // has some useful funtions for stuff not available from the tables

// ROOT headers
#include "TLorentzVector.h"
#include "TEfficiency.h"
//#include <Math/Vector4D.h> // this should apparently be used instead of TLorentzVector, e.g. "ROOT::Math::PxPyPzMVector vector;"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using FullUDCollision = aod::UDCollisions::iterator;
using FullUDTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksPID>;

struct upcRhoAnalysis {

    double PcEtaCut = 0.9; // cut on track eta as per physics coordination recommendation
    // configurables
    //Configurable<bool> verbosity{"verbosity", true, "verbosity"};
    Configurable<bool> requireTOF{"requireTOF", false, "requireTOF"};
    Configurable<double> tpcNSigmaPiCut{"tpcNSigmaPiCut", 3.0, "tpcNSigmaPiCut"};
    Configurable<double> tofNSigmaPiCut{"tofNSigmaPiCut", 3.0, "tofNSigmaPiCut"};
    Configurable<double> ptMaxCut{"ptMaxCut", 10.0, "ptMaxCut"};

    ConfigurableAxis mAxis{"mAxis", {500, 0.0, 5.0}, "m (GeV/#it{c}^{2})"};
    ConfigurableAxis ptAxis{"ptAxis", {100, 0.0, 10.0}, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis etaAxis{"etaAxis", {90, -0.9, 0.9}, "#eta"};
    ConfigurableAxis phiAxis{"phiAxis", {360, 0.0, 2.0*o2::constants::math::PI}, "#phi"};

    HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

    void init(o2::framework::InitContext&){
        registry.add("hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
        registry.add("hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
        registry.add("hEta", ";#eta;counts", kTH1D, {etaAxis});
        registry.add("hPhi", ";#phi;counts", kTH1D, {phiAxis});
        registry.add("hCharge", ";charge;counts", kTH1D, {{5, -2.0, 2.0}});
        registry.add("hNTracks", ";N_{tracks};counts", kTH1D, {{11, 0.0, 10.0}});
        registry.add("hTpcNSigmaPi", ";TPC n#sigma_{#pi};counts", kTH1D, {{100, -5.0, 5.0}});
        registry.add("hTofNSigmaPi", ";TOF n#sigma_{#pi};counts", kTH1D, {{100, -5.0, 5.0}});
    }

    template <typename T>
    bool passesCuts(T const& track) {
        if (track.pt() > ptMaxCut) return false;
        if (eta(track.px(), track.py(), track.pz()) > PcEtaCut) return false;
        if (track.tpcNSigmaPi() > tpcNSigmaPiCut) return false;
        if (requireTOF && !track.hasTOF()) return false;
        if (track.hasTOF() && track.tofNSigmaPi() > 3.0) return false;
        return true;
    }

    void process(FullUDCollision const& collision, FullUDTracks const& tracks) {
        for(const auto& track : tracks) {
            if (!passesCuts(track)) continue;
            // fill histograms
            registry.get<TH1>(HIST("hPt"))->Fill(track.pt());
            registry.get<TH1>(HIST("hEta"))->Fill(eta(track.px(), track.py(), track.pz()));
            registry.get<TH1>(HIST("hPhi"))->Fill(phi(track.px(), track.py()));
            registry.get<TH1>(HIST("hCharge"))->Fill(track.sign());
            registry.get<TH1>(HIST("hNTracks"))->Fill(tracks.size());
            registry.get<TH1>(HIST("hTpcNSigmaPi"))->Fill(track.tpcNSigmaPi());
            registry.get<TH1>(HIST("hTofNSigmaPi"))->Fill(track.tofNSigmaPi());
        }
    }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
    return WorkflowSpec{
        o2::framework::adaptAnalysisTask<upcRhoAnalysis>(cfgc)
    };
}