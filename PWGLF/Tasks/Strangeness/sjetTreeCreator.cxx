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

/// \file sjetTreeCreator.cxx
/// \brief Task for building a TTree with information about the jet, its consistuents,
/// to be used as input for machine learning s-jet identification.
/// \author Lorenzo Bernardinis (lorenzo.bernardinis@cern.ch)
///
/// Inspired by PWGJE/Tasks/bjetTreeCreator.cxx

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetTaggingUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetTagging.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <string>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace jet_info
{
// DECLARE_SOA_INDEX_COLUMN(JetIndex, jetindex); //! The jet index
DECLARE_SOA_COLUMN(JetpT, jetpT, float);           //! jet pT
DECLARE_SOA_COLUMN(JetEta, jetEta, float);         //! jet eta
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);         //! jet phi
DECLARE_SOA_COLUMN(NTracks, nTracks, int16_t);     //! number of charged tracks inside the jet
DECLARE_SOA_COLUMN(NSV, nSV, int16_t);             //! Number of secondary vertices in the jet
DECLARE_SOA_COLUMN(JetMass, jetMass, float);       //! The jet mass
DECLARE_SOA_COLUMN(JetFlavor, jetFlavor, int16_t); //! The jet flavor (b, c, s or udg)
DECLARE_SOA_COLUMN(JetR, jetR, int16_t);           //! The jet radius
} // namespace jet_info

DECLARE_SOA_TABLE(sjetParams, "AOD", "SJETPARAM",
                  o2::soa::Index<>,
                  jet_info::JetpT,
                  jet_info::JetEta,
                  jet_info::JetPhi,
                  jet_info::NTracks,
                  jet_info::JetMass,
                  jet_info::JetFlavor,
                  jet_info::JetR);

using sjetParam = sjetParams::iterator;

namespace track_info
{
DECLARE_SOA_INDEX_COLUMN(sjetParam, jetindex);                             //! The jet index
DECLARE_SOA_COLUMN(TrackpT, trackpT, float);                               //! The track pT
DECLARE_SOA_COLUMN(TrackEta, trackEta, float);                             //! The track eta
DECLARE_SOA_COLUMN(DotProdTrackJet, dotProdTrackJet, float);               //! The dot product between the track and the jet
DECLARE_SOA_COLUMN(DotProdTrackJetOverJet, dotProdTrackJetOverJet, float); //! The dot product between the track and the jet over the jet momentum
DECLARE_SOA_COLUMN(DeltaRJetTrack, deltaRJetTrack, float);                 //! The DR jet-track
DECLARE_SOA_COLUMN(SignedIP2D, signedIP2D, float);                         //! The track signed 2D IP
DECLARE_SOA_COLUMN(SignedIP2DSign, signedIP2DSign, float);                 //! The track signed 2D IP significance
DECLARE_SOA_COLUMN(SignedIPz, signedIPz, float);                           //! The track signed z IP
DECLARE_SOA_COLUMN(SignedIPzSign, signedIPzSign, float);                   //! The track signed z IP significance
DECLARE_SOA_COLUMN(SignedIP3DSign, signedIP3DSign, float);                 //! The track signed 3D IP significance
DECLARE_SOA_COLUMN(MomFraction, momFraction, float);                       //! The track momentum fraction of the jets
DECLARE_SOA_COLUMN(DeltaRTrackVertex, deltaRTrackVertex, float);           //! DR between the track and the closest SV, to be decided whether to add to or not
DECLARE_SOA_COLUMN(TrackPhi, trackPhi, float);                             //! The track phi
DECLARE_SOA_COLUMN(TrackCharge, trackCharge, float);                       //! The track sign (charge)
DECLARE_SOA_COLUMN(TrackITSChi2NCl, trackITSChi2NCl, float);               //! The track ITS Chi2NCl
DECLARE_SOA_COLUMN(TrackTPCChi2NCl, trackTPCChi2NCl, float);               //! The track TPC Chi2NCl
DECLARE_SOA_COLUMN(TrackITSNCls, trackITSNCls, float);                     //! The track ITS NCls
DECLARE_SOA_COLUMN(TrackTPCNCls, trackTPCNCls, float);                     //! The track TPC NCls (Found)
DECLARE_SOA_COLUMN(TrackTPCNCrossedRows, trackTPCNCrossedRows, float);     //! The track TPC NCrossedRows
DECLARE_SOA_COLUMN(TrackOrigin, trackOrigin, int);                         //! The track origin label for GNN track origin predictions
DECLARE_SOA_COLUMN(TrackVtxIndex, trackVtxIndex, int);                     //! The track vertex index for GNN vertex predictions
// DECLARE_SOA_COLUMN(DCATrackJet, dcaTrackJet, float);                   //! The distance between track and jet, unfortunately it cannot be calculated in O2
} // namespace track_info

DECLARE_SOA_TABLE(sjetTracksParams, "AOD", "SJETTRACKSPARAM",
                  o2::soa::Index<>,
                  track_info::sjetParamId,
                  track_info::TrackpT,
                  track_info::TrackEta,
                  track_info::DotProdTrackJet,
                  track_info::DotProdTrackJetOverJet,
                  track_info::DeltaRJetTrack,
                  track_info::SignedIP2D,
                  track_info::SignedIP2DSign,
                  track_info::SignedIPz,
                  track_info::SignedIPzSign,
                  track_info::SignedIP3DSign,
                  track_info::MomFraction,
                  track_info::DeltaRTrackVertex);

using sjetTracksParam = sjetTracksParams::iterator;

DECLARE_SOA_TABLE(sjetTracksParamsExtra, "AOD", "SJETTRACKSEXTRA",
                  // o2::soa::Index<>,
                  track_info::TrackPhi,
                  track_info::TrackCharge,
                  track_info::TrackITSChi2NCl,
                  track_info::TrackTPCChi2NCl,
                  track_info::TrackITSNCls,
                  track_info::TrackTPCNCls,
                  track_info::TrackTPCNCrossedRows);

using sjetTracksParamExtra = sjetTracksParamsExtra::iterator;

namespace constituents
{
DECLARE_SOA_INDEX_COLUMN(sjetParam, jetindex);
DECLARE_SOA_ARRAY_INDEX_COLUMN(sjetTracksParam, tracks);
} // namespace constituents

DECLARE_SOA_TABLE(sjetConstituents, "AOD", "SJETCONSTIT",
                  constituents::sjetParamId,
                  constituents::sjetTracksParamIds);

} // namespace o2::aod

struct SjetTreeCreator {

  Produces<aod::sjetParams> sjetParamsTable;
  Produces<aod::sjetTracksParams> sjetTracksParamsTable;
  Produces<aod::sjetTracksParamsExtra> sjetTracksExtraTable;
  Produces<aod::sjetConstituents> sjetConstituentsTable;

  HistogramRegistry registry;

  // event level configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};

  Configurable<std::vector<double>> jetPtBins{"jetPtBins", std::vector<double>{5, 1000}, "jet pT bins for reduction"};
  Configurable<std::vector<double>> jetReductionFactors{"jetReductionFactors", std::vector<double>{0.0}, "jet reduction factors"};

  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};

  // track level configurables
  Configurable<float> trackPtMin{"trackPtMin", 0.5, "minimum track pT"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track eta"};

  Configurable<float> maxIPxy{"maxIPxy", 10, "maximum track DCA in xy plane"};
  Configurable<float> maxIPz{"maxIPz", 10, "maximum track DCA in z direction"};

  // jet level configurables
  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT"};
  Configurable<float> jetPtMax{"jetPtMax", 1000.0, "maximum jet pT"};
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};

  Configurable<float> eventReductionFactor{"eventReductionFactor", 0.0, "Percentage of events to be removed"};

  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};

  Configurable<bool> produceTree{"produceTree", true, "produce the jet TTree"};

  Configurable<float> vtxRes{"vtxRes", 0.01, "Vertex position resolution (cluster size) for GNN vertex predictions (cm)"};

  std::vector<int> eventSelectionBits;

  std::vector<double> jetRadiiValues;
  std::vector<double> jetPtBinsReduction;
  std::vector<double> jetReductionFactorsPt;

  void init(InitContext const&)
  {
    // Seed the random number generator using current time
    std::srand(static_cast<unsigned int>(std::time(nullptr)));

    jetRadiiValues = (std::vector<double>)jetRadii;
    jetPtBinsReduction = (std::vector<double>)jetPtBins;
    jetReductionFactorsPt = (std::vector<double>)jetReductionFactors;

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));

    registry.add("h_vertexZ", "Vertex Z;#it{Z} (cm)", {HistType::kTH1F, {{40, -20.0, 20.0}}});

    registry.add("h2_nTracks_jetpT", "Number of tracks;#it{p}_{T,jet} (GeV/#it{c});nTracks", {HistType::kTH2F, {{200, 0., 200.}, {100, 0, 100.0}}});

    registry.add("h2_SIPs2D_jetpT", "2D IP significance;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
    registry.add("h2_SIPs3D_jetpT", "3D IP significance;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
    registry.add("h2_jetMass_jetpT", "Jet mass;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{jet} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 50.0}}});

    if (doprocessMCJets) {
      registry.add("h2_SIPs2D_jetpT_hfjet", "2D IP significance hf-jets;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_SIPs3D_jetpT_hfjet", "3D IP significance hf-jets;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_jetMass_jetpT_hfjet", "Jet mass hf-jets;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{jet} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 50.0}}});
      registry.add("h_jetpT_detector_hfjet", "Jet transverse momentum hf-jets;#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});

      registry.add("h2_SIPs2D_jetpT_sjet", "2D IP significance s-jets;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_SIPs3D_jetpT_sjet", "3D IP significance s-jets;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_jetMass_jetpT_sjet", "Jet mass s-jets;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{jet} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 50.0}}});
      registry.add("h_jetpT_detector_sjet", "Jet transverse momentum s-jets;#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});

      registry.add("h2_SIPs2D_jetpT_udgjet", "2D IP significance udg-jets;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_SIPs3D_jetpT_udgjet", "3D IP significance udg-jets;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_jetMass_jetpT_udgjet", "Jet mass udg-jets;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{jet} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 50.0}}});
      registry.add("h_jetpT_detector_udgjet", "Jet transverse momentum udg-jets;#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});

      // Jet information
      registry.add("h_jet_pt", "jet_pt;#it{p}_{T}^{ch jet} (GeV/#it{c});Entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_jet_eta", "jet_eta;#it{#eta}_{ch jet};Entries", {HistType::kTH1F, {{200, -2., 2.}}});
      registry.add("h_jet_phi", "jet_phi;#it{#phi}_{ch jet};Entries", {HistType::kTH1F, {{200, 0., o2::constants::math::TwoPI}}});
      registry.add("h_n_trks", "n_trks;#it{n}_{tracks};Entries", {HistType::kTH1F, {{50, 0., 50.}}});
      registry.add("h_jet_mass", "jet_mass;#it{m}_{jet} (GeV/#it{c}^2);Entries", {HistType::kTH1F, {{200, 0., 50.}}});

      registry.add("h_jet_flav", "jet_flav;jet flavor;Entries", {HistType::kTH1F, {{4, -0.5, 3.5}}});
      auto hJetFlavor = registry.get<TH1>(HIST("h_jet_flav"));
      hJetFlavor->GetXaxis()->SetBinLabel(1, "no mcparticle"); // bin 1
      hJetFlavor->GetXaxis()->SetBinLabel(2, "hf-jet");        // bin 2 --> flavour number 1+2
      hJetFlavor->GetXaxis()->SetBinLabel(3, "s-jet");         // bin 3 --> flavour number 7
      hJetFlavor->GetXaxis()->SetBinLabel(4, "udg-jet");       // bin 4 --> flavour number 6

      // Track information
      registry.add("h_trk_pt", "trk_pt;#it{p}_{T} (GeV/#it{c});Entries", {HistType::kTH1F, {{200, 0., 100.}}});
      registry.add("h_trk_eta", "trk_eta;#it{#eta};Entries", {HistType::kTH1F, {{200, -2., 2.}}});
      registry.add("h_trk_phi", "trk_phi;#it{#phi};Entries", {HistType::kTH1F, {{200, 0., o2::constants::math::TwoPI}}});
      registry.add("h_trk_charge", "trk_charge;#it{q};Entries", {HistType::kTH1F, {{3, -1.5, 1.5}}});
      registry.add("h_trk_dcaxy", "trk_dcaxy;#it{DCA}_{xy} (cm);Entries", {HistType::kTH1F, {{200, -0.1, 0.1}}});
      registry.add("h_trk_dcaxyz", "trk_dcaxyz;#it{DCA}_{xyz} (cm);Entries", {HistType::kTH1F, {{200, -0.1, 0.1}}});
      registry.add("h_trk_sigmadcaxy", "trk_sigmadcaxy;#it{#sigma}_{#it{DCA}_{xy}} (cm);Entries", {HistType::kTH1F, {{200, 0., 0.1}}});
      registry.add("h_trk_sigmadcaxyz", "trk_sigmadcaxyz;#it{#sigma}_{#it{DCA}_{xyz}} (cm);Entries", {HistType::kTH1F, {{200, 0., 0.1}}});
      registry.add("h_trk_itsncls", "trk_itsncls;ITS NCls;Entries", {HistType::kTH1F, {{10, 0., 10.}}});
      registry.add("h_trk_tpcncls", "trk_tpcncls;TPC NCls (Found);Entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_trk_tpcncrs", "trk_tpcncrs;TPC NCrossedRows;Entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_trk_itschi2ncl", "trk_itschi2ncl;ITS #it{#chi}^{2}/ndf;Entries", {HistType::kTH1F, {{200, 0., 20.}}});
      registry.add("h_trk_tpcchi2ncl", "trk_tpcchi2ncl;TPC #it{#chi}^{2}/ndf;Entries", {HistType::kTH1F, {{200, 0., 10.}}});
      registry.add("h2_trk_jtrackpt_vs_origtrackpt", "JTracks::pt vs Tracks::pt", {HistType::kTH2F, {{200, 0., 100.}, {200, 0., 100.}}});
    }
  }

  // FIXME filtering only works when you loop directly over the list, but if you loop over it as a constituent they will not be filtered
  Filter collisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter trackCuts = (aod::jtrack::pt > trackPtMin && aod::jtrack::pt < trackPtMax &&
                      aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter partCuts = (aod::jmcparticle::pt >= trackPtMin && aod::jmcparticle::pt < trackPtMax);
  Filter jetFilter = (aod::jet::pt >= jetPtMin && aod::jet::pt <= jetPtMax &&
                      aod::jet::eta < jetEtaMax - aod::jet::r / 100.f && aod::jet::eta > jetEtaMin + aod::jet::r / 100.f);

  using FilteredCollision = soa::Filtered<soa::Join<aod::JetCollisions, aod::JCollisionPIs>>;
  using JetTrackswID = soa::Filtered<soa::Join<aod::JetTracks, aod::JTrackExtras, aod::JTrackPIs>>;
  using JetTracksMCDwID = soa::Filtered<soa::Join<aod::JetTracksMCD, aod::JTrackExtras, aod::JTrackPIs>>;
  using DataJets = soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents, aod::DataSecondaryVertex3ProngIndices>>;

  using OriginalTracks = soa::Join<aod::Tracks, aod::TracksCov, aod::TrackSelection, aod::TracksDCA, aod::TracksDCACov, aod::TracksExtra>;

  // Function to get the reduction factor based on jet pT
  double getReductionFactor(double jetPT)
  {
    // Loop through the jetPtBins vector
    for (size_t ibin = 0; ibin < jetPtBinsReduction.size() - 1; ++ibin) {
      if (jetPT >= jetPtBinsReduction[ibin] && jetPT < jetPtBinsReduction[ibin + 1]) {
        return jetReductionFactorsPt[ibin];
      }
    }

    // If jetPT is above the last bin, use the last reduction factor
    if (jetPT >= jetPtBinsReduction.back()) {
      return jetReductionFactorsPt.back();
    }

    // If jetPT is below the first bin, return the first reduction factor
    return jetReductionFactorsPt.front();
  }

  template <typename AnyCollision, typename AnalysisJet, typename AnyTracks, typename AnyOriginalTracks>
  void analyzeJetTrackInfo(AnyCollision const& /*collision*/,
                           AnalysisJet const& analysisJet,
                           AnyTracks const& /*allTracks*/,
                           AnyOriginalTracks const&,
                           std::vector<int>& trackIndices,
                           int jetFlavor = 0,
                           double eventweight = 1.0)
  {
    for (const auto& constituent : analysisJet.template tracks_as<AnyTracks>()) {

      if (constituent.pt() < trackPtMin) {
        continue;
      }

      double deltaRJetTrack = jetutilities::deltaR(analysisJet, constituent);
      double dotProduct = RecoDecay::dotProd(std::array<float, 3>{analysisJet.px(), analysisJet.py(), analysisJet.pz()},
                                             std::array<float, 3>{constituent.px(), constituent.py(), constituent.pz()});
      int sign = jettaggingutilities::getGeoSign(analysisJet, constituent);

      registry.fill(HIST("h2_SIPs2D_jetpT"), analysisJet.pt(), sign * std::abs(constituent.dcaXY()) / constituent.sigmadcaXY(), eventweight);
      registry.fill(HIST("h2_SIPs3D_jetpT"), analysisJet.pt(), sign * std::abs(constituent.dcaXYZ()) / constituent.sigmadcaXYZ(), eventweight);

      if (doprocessMCJets) {
        if (jetFlavor == JetTaggingSpecies::beauty || jetFlavor == JetTaggingSpecies::charm) {
          registry.fill(HIST("h2_SIPs2D_jetpT_hfjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXY()) / constituent.sigmadcaXY(), eventweight);
          registry.fill(HIST("h2_SIPs3D_jetpT_hfjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXYZ()) / constituent.sigmadcaXYZ(), eventweight);
        } else if (jetFlavor == JetTaggingSpecies::strange) {
          registry.fill(HIST("h2_SIPs2D_jetpT_sjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXY()) / constituent.sigmadcaXY(), eventweight);
          registry.fill(HIST("h2_SIPs3D_jetpT_sjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXYZ()) / constituent.sigmadcaXYZ(), eventweight);
        } else {
          registry.fill(HIST("h2_SIPs2D_jetpT_udgjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXY()) / constituent.sigmadcaXY(), eventweight);
          registry.fill(HIST("h2_SIPs3D_jetpT_udgjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXYZ()) / constituent.sigmadcaXYZ(), eventweight);
        }

        auto origConstit = constituent.template track_as<AnyOriginalTracks>();

        // Track information
        registry.fill(HIST("h_trk_pt"), constituent.pt(), eventweight);
        registry.fill(HIST("h_trk_eta"), constituent.eta(), eventweight);
        registry.fill(HIST("h_trk_phi"), origConstit.phi(), eventweight);
        registry.fill(HIST("h_trk_charge"), constituent.sign(), eventweight);
        registry.fill(HIST("h_trk_dcaxy"), std::abs(constituent.dcaXY()) * sign, eventweight);
        registry.fill(HIST("h_trk_dcaxyz"), std::abs(constituent.dcaXYZ()) * sign, eventweight);
        registry.fill(HIST("h_trk_sigmadcaxy"), constituent.sigmadcaXY(), eventweight);
        registry.fill(HIST("h_trk_sigmadcaxyz"), constituent.sigmadcaXYZ(), eventweight);
        registry.fill(HIST("h_trk_itsncls"), origConstit.itsNCls(), eventweight);
        registry.fill(HIST("h_trk_tpcncls"), origConstit.tpcNClsFound(), eventweight);
        registry.fill(HIST("h_trk_tpcncrs"), origConstit.tpcNClsCrossedRows(), eventweight);
        registry.fill(HIST("h_trk_itschi2ncl"), origConstit.itsChi2NCl(), eventweight);
        registry.fill(HIST("h_trk_tpcchi2ncl"), origConstit.tpcChi2NCl(), eventweight);
        registry.fill(HIST("h2_trk_jtrackpt_vs_origtrackpt"), constituent.pt(), origConstit.pt(), eventweight);

        if (produceTree) {
          sjetTracksExtraTable(/*sjetParamsTable.lastIndex() + 1, */
                               origConstit.phi(),
                               constituent.sign(),
                               origConstit.itsChi2NCl(),
                               origConstit.tpcChi2NCl(),
                               origConstit.itsNCls(),
                               origConstit.tpcNClsFound(),
                               origConstit.tpcNClsCrossedRows());

          sjetTracksParamsTable(sjetParamsTable.lastIndex() + 1,
                                constituent.pt(),
                                constituent.eta(),
                                dotProduct,
                                dotProduct / analysisJet.p(),
                                deltaRJetTrack,
                                std::abs(constituent.dcaXY()) * sign,
                                constituent.sigmadcaXY(),
                                std::abs(constituent.dcaZ()) * sign,
                                constituent.sigmadcaZ(),
                                constituent.sigmadcaXYZ(),
                                constituent.p() / analysisJet.p(),
                                0.);
        }
        trackIndices.push_back(sjetTracksParamsTable.lastIndex());
      }
    }
  }

  void processDummy(FilteredCollision::iterator const& /*collision*/)
  {
  }
  PROCESS_SWITCH(SjetTreeCreator, processDummy, "Dummy process function turned on by default", true);

  using MCDJetTable = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets,
                                              aod::ChargedMCDetectorLevelJetConstituents,
                                              aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets,
                                              aod::MCDSecondaryVertex3ProngIndices,
                                              aod::ChargedMCDetectorLevelJetEventWeights>>;

  using MCPJetTable = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets,
                                              aod::ChargedMCParticleLevelJetConstituents,
                                              aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets,
                                              aod::ChargedMCParticleLevelJetEventWeights>>;

  using FilteredCollisionMCD = soa::Filtered<soa::Join<aod::JetCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>>;

  Preslice<aod::JMcParticles> mcParticlesPerCollision = aod::jmcparticle::mcCollisionId;
  Preslice<MCPJetTable> mcpJetsPerCollision = aod::jet::mcCollisionId;

  using MCDJetTableNoSV = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets,
                                                  aod::ChargedMCDetectorLevelJetConstituents,
                                                  aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets,
                                                  aod::ChargedMCDetectorLevelJetEventWeights>>;

  using JetParticleswID = soa::Join<aod::JetParticles, aod::JMcParticlePIs>;

  void processMCJets(FilteredCollisionMCD::iterator const& collision,
                     aod::JMcCollisions const&,
                     MCDJetTableNoSV const& MCDjets,
                     MCPJetTable const& MCPjets,
                     JetTracksMCDwID const& allTracks,
                     JetParticleswID const& MCParticles,
                     OriginalTracks const& origTracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) ||
        (static_cast<double>(std::rand()) / RAND_MAX < eventReductionFactor)) {
      return;
    }

    registry.fill(HIST("h_vertexZ"), collision.posZ());

    auto const mcParticlesPerColl = MCParticles.sliceBy(mcParticlesPerCollision, collision.mcCollisionId());
    auto const mcPJetsPerColl = MCPjets.sliceBy(mcpJetsPerCollision, collision.mcCollisionId());

    for (const auto& analysisJet : MCDjets) {

      bool jetIncluded = false;
      for (const auto& jetR : jetRadiiValues) {
        if (analysisJet.r() == static_cast<int>(jetR * 100)) {
          jetIncluded = true;
          break;
        }
      }

      if (!jetIncluded) {
        continue;
      }

      std::vector<int> indicesTracks;

      int16_t jetFlavorIdx = 0;
      int16_t jetFlavor = 0;

      for (const auto& mcpjet : analysisJet.template matchedJetGeo_as<MCPJetTable>()) {
        jetFlavor = jettaggingutilities::getSJetFlavor(mcpjet, mcParticlesPerColl);
      }

      if (jetFlavor == JetTaggingSpecies::strange) {
        jetFlavorIdx = 2;
      } else if (jetFlavor == JetTaggingSpecies::udg) {
        jetFlavorIdx = 3;
      } else if (jetFlavor == JetTaggingSpecies::charm || jetFlavor == JetTaggingSpecies::beauty) {
        jetFlavorIdx = 1;
      } else {
        jetFlavorIdx = jetFlavor;
      }

      if ((jetFlavor != JetTaggingSpecies::strange) &&
          (static_cast<double>(std::rand()) / RAND_MAX < getReductionFactor(analysisJet.pt()))) {
        continue;
      }

      float eventWeight = analysisJet.eventWeight();

      analyzeJetTrackInfo(collision, analysisJet, allTracks, origTracks, indicesTracks, jetFlavor, eventWeight);

      registry.fill(HIST("h2_jetMass_jetpT"), analysisJet.pt(), analysisJet.mass(), eventWeight);
      registry.fill(HIST("h2_nTracks_jetpT"), analysisJet.pt(), indicesTracks.size());

      // Jet info
      registry.fill(HIST("h_jet_pt"), analysisJet.pt());
      registry.fill(HIST("h_jet_eta"), analysisJet.eta());
      registry.fill(HIST("h_jet_phi"), analysisJet.phi());

      registry.fill(HIST("h_jet_flav"), jetFlavorIdx);
      registry.fill(HIST("h_n_trks"), indicesTracks.size());
      registry.fill(HIST("h_jet_mass"), analysisJet.mass());

      if (jetFlavor == JetTaggingSpecies::beauty || jetFlavor == JetTaggingSpecies::charm) {
        registry.fill(HIST("h2_jetMass_jetpT_hfjet"), analysisJet.pt(), analysisJet.mass(), eventWeight);
        registry.fill(HIST("h_jetpT_detector_hfjet"), analysisJet.pt(), eventWeight);
      } else if (jetFlavor == JetTaggingSpecies::strange) {
        registry.fill(HIST("h2_jetMass_jetpT_sjet"), analysisJet.pt(), analysisJet.mass(), eventWeight);
        registry.fill(HIST("h_jetpT_detector_sjet"), analysisJet.pt(), eventWeight);
      } else {
        registry.fill(HIST("h2_jetMass_jetpT_udgjet"), analysisJet.pt(), analysisJet.mass(), eventWeight);
        registry.fill(HIST("h_jetpT_detector_udgjet"), analysisJet.pt(), eventWeight);
      }

      if (produceTree) {
        sjetConstituentsTable(sjetParamsTable.lastIndex() + 1, indicesTracks);
        sjetParamsTable(analysisJet.pt(),
                        analysisJet.eta(),
                        analysisJet.phi(),
                        indicesTracks.size(),
                        analysisJet.mass(),
                        jetFlavor,
                        analysisJet.r());
      }
    }
  }
  PROCESS_SWITCH(SjetTreeCreator, processMCJets, "jet information in MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<SjetTreeCreator>(cfgc)};
}
