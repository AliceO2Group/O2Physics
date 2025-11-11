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

/// \file bjetTreeCreator.cxx
/// \brief Task for building a TTree with information about the jet, its consistuents, and the secondary vertices inside
/// to be used for as input for machine learning b-jet identification.
///
/// \author Hadi Hassan <hadi.hassan@cern.ch>, University of Jyväskylä

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
namespace jetInfo
{
// DECLARE_SOA_INDEX_COLUMN(JetIndex, jetindex); //! The jet index
DECLARE_SOA_COLUMN(JetpT, jetpt, float);                   //! jet pT
DECLARE_SOA_COLUMN(JetEta, jeteta, float);                 //! jet eta
DECLARE_SOA_COLUMN(JetPhi, jetphi, float);                 //! jet phi
DECLARE_SOA_COLUMN(NTracks, nTracks, int16_t);             //! number of charged tracks inside the jet
DECLARE_SOA_COLUMN(NSV, nSV, int16_t);                     //! Number of secondary vertices in the jet
DECLARE_SOA_COLUMN(JetMass, mass, float);                  //! The jet mass
DECLARE_SOA_COLUMN(JetFlavor, jetFl, int16_t);             //! The jet flavor (b, c, or lf)
DECLARE_SOA_COLUMN(JetR, jetR, int16_t);                   //! The jet radius
DECLARE_SOA_COLUMN(JetEventWeight, jetEventWeight, float); //! The jet event weight for pTHat weighting
} // namespace jetInfo

DECLARE_SOA_TABLE(bjetParams, "AOD", "BJETPARAM",
                  o2::soa::Index<>,
                  jetInfo::JetpT,
                  jetInfo::JetEta,
                  jetInfo::JetPhi,
                  jetInfo::NTracks,
                  jetInfo::NSV,
                  jetInfo::JetMass,
                  jetInfo::JetFlavor,
                  jetInfo::JetR);

using bjetParam = bjetParams::iterator;

DECLARE_SOA_TABLE(bjetParamsExtra, "AOD", "BJETPARAMSEXTRA",
                  // o2::soa::Index<>,
                  jetInfo::JetEventWeight);

using bjetParamExtra = bjetParamsExtra::iterator;

namespace trackInfo
{
DECLARE_SOA_INDEX_COLUMN(bjetParam, jetindex);                         //! The jet index
DECLARE_SOA_COLUMN(TrackpT, trackpt, float);                           //! The track pT
DECLARE_SOA_COLUMN(TrackEta, tracketa, float);                         //! The track eta
DECLARE_SOA_COLUMN(DotProdTrackJet, trackdotjet, float);               //! The dot product between the track and the jet
DECLARE_SOA_COLUMN(DotProdTrackJetOverJet, trackdotjetoverjet, float); //! The dot product between the track and the jet over the jet momentum
DECLARE_SOA_COLUMN(DeltaRJetTrack, rjettrack, float);                  //! The DR jet-track
DECLARE_SOA_COLUMN(SignedIP2D, ip2d, float);                           //! The track signed 2D IP
DECLARE_SOA_COLUMN(SignedIP2DSign, ip2dsigma, float);                  //! The track signed 2D IP significance
DECLARE_SOA_COLUMN(SignedIPz, ipz, float);                             //! The track signed z IP
DECLARE_SOA_COLUMN(SignedIPzSign, ipzsigma, float);                    //! The track signed z IP significance
DECLARE_SOA_COLUMN(SignedIP3DSign, ip3dsigma, float);                  //! The track signed 3D IP significance
DECLARE_SOA_COLUMN(MomFraction, momfraction, float);                   //! The track momentum fraction of the jets
DECLARE_SOA_COLUMN(DeltaRTrackVertex, rtrackvertex, float);            //! DR between the track and the closest SV, to be decided whether to add to or not
DECLARE_SOA_COLUMN(TrackPhi, trackphi, float);                         //! The track phi
DECLARE_SOA_COLUMN(TrackCharge, trackcharge, float);                   //! The track sign (charge)
DECLARE_SOA_COLUMN(TrackITSChi2NCl, trackitschi2ncl, float);           //! The track ITS Chi2NCl
DECLARE_SOA_COLUMN(TrackTPCChi2NCl, tracktpcchi2ncl, float);           //! The track TPC Chi2NCl
DECLARE_SOA_COLUMN(TrackITSNCls, trackitsncls, float);                 //! The track ITS NCls
DECLARE_SOA_COLUMN(TrackTPCNCls, tracktpcncls, float);                 //! The track TPC NCls (Found)
DECLARE_SOA_COLUMN(TrackTPCNCrossedRows, tracktpcncrossedrows, float); //! The track TPC NCrossedRows
DECLARE_SOA_COLUMN(TrackOrigin, trk_origin, int);                      //! The track origin label for GNN track origin predictions
DECLARE_SOA_COLUMN(TrackVtxIndex, trk_vtx_index, int);                 //! The track vertex index for GNN vertex predictions
// DECLARE_SOA_COLUMN(DCATrackJet, dcatrackjet, float);                              //! The distance between track and jet, unfortunately it cannot be calculated in O2
} // namespace trackInfo

DECLARE_SOA_TABLE(bjetTracksParams, "AOD", "BJETTRACKSPARAM",
                  o2::soa::Index<>,
                  trackInfo::bjetParamId,
                  trackInfo::TrackpT,
                  trackInfo::TrackEta,
                  trackInfo::DotProdTrackJet,
                  trackInfo::DotProdTrackJetOverJet,
                  trackInfo::DeltaRJetTrack,
                  trackInfo::SignedIP2D,
                  trackInfo::SignedIP2DSign,
                  trackInfo::SignedIPz,
                  trackInfo::SignedIPzSign,
                  trackInfo::SignedIP3DSign,
                  trackInfo::MomFraction,
                  trackInfo::DeltaRTrackVertex);

using bjetTracksParam = bjetTracksParams::iterator;

DECLARE_SOA_TABLE(bjetTracksParamsExtra, "AOD", "BJETTRACKSEXTRA",
                  // o2::soa::Index<>,
                  trackInfo::TrackPhi,
                  trackInfo::TrackCharge,
                  trackInfo::TrackITSChi2NCl,
                  trackInfo::TrackTPCChi2NCl,
                  trackInfo::TrackITSNCls,
                  trackInfo::TrackTPCNCls,
                  trackInfo::TrackTPCNCrossedRows,
                  trackInfo::TrackOrigin,
                  trackInfo::TrackVtxIndex);

using bjetTracksParamExtra = bjetTracksParamsExtra::iterator;

namespace SVInfo
{
DECLARE_SOA_INDEX_COLUMN(bjetParam, jetindex);            //! The jet index
DECLARE_SOA_COLUMN(SVpT, svpt, float);                    //! The SV pT
DECLARE_SOA_COLUMN(DeltaRSVJet, rsvjet, float);           //! The DR jet-SV
DECLARE_SOA_COLUMN(SVMass, mass, float);                  //! The SV mass
DECLARE_SOA_COLUMN(SVfE, svfe, float);                    //! The SV energy fraction
DECLARE_SOA_COLUMN(IPXY, ipxy, float);                    //! The SV 2D IP
DECLARE_SOA_COLUMN(CPA, cpa, float);                      //! Cosine pointing angle between the SV direction and momentum
DECLARE_SOA_COLUMN(Chi2PCA, chi2pca, float);              //! Sum of (non-weighted) distances of the secondary vertex to its prongsm
DECLARE_SOA_COLUMN(Dispersion, dispersion, float);        //! The SV dispersion
DECLARE_SOA_COLUMN(DecayLength2D, lxy, float);            //! The decay length of the SV in XY
DECLARE_SOA_COLUMN(DecayLength2DError, lxysigma, float);  //! The decay length of the SV in XY significance
DECLARE_SOA_COLUMN(DecayLength3D, lxyz, float);           //! The decay length of the SV in 3D
DECLARE_SOA_COLUMN(DecayLength3DError, lxyzsigma, float); //! The decay length of the SV in 3d significance
} // namespace SVInfo

DECLARE_SOA_TABLE(bjetSVParams, "AOD", "BJETSVPARAM",
                  o2::soa::Index<>,
                  SVInfo::bjetParamId,
                  SVInfo::SVpT,
                  SVInfo::DeltaRSVJet,
                  SVInfo::SVMass,
                  SVInfo::SVfE,
                  SVInfo::IPXY,
                  SVInfo::CPA,
                  SVInfo::Chi2PCA,
                  SVInfo::Dispersion,
                  SVInfo::DecayLength2D,
                  SVInfo::DecayLength2DError,
                  SVInfo::DecayLength3D,
                  SVInfo::DecayLength3DError);

using bjetSVParam = bjetSVParams::iterator;

namespace constituents
{
DECLARE_SOA_INDEX_COLUMN(bjetParam, jetindex);
DECLARE_SOA_ARRAY_INDEX_COLUMN(bjetTracksParam, tracks);
DECLARE_SOA_ARRAY_INDEX_COLUMN(bjetSVParam, svs);
} // namespace constituents

DECLARE_SOA_TABLE(bjetConstituents, "AOD", "BJETCONSTIT",
                  constituents::bjetParamId,
                  constituents::bjetTracksParamIds,
                  constituents::bjetSVParamIds);

} // namespace o2::aod

struct BJetTreeCreator {

  Produces<aod::bjetParams> bjetParamsTable;
  Produces<aod::bjetParamsExtra> bjetParamsExtraTable;
  Produces<aod::bjetTracksParams> bjetTracksParamsTable;
  Produces<aod::bjetTracksParamsExtra> bjetTracksExtraTable;
  Produces<aod::bjetSVParams> bjetSVParamsTable;
  Produces<aod::bjetConstituents> bjetConstituentsTable;

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

  Configurable<bool> useQuarkDef{"useQuarkDef", true, "Flag whether to use quarks or hadrons for determining the jet flavor"};

  // track level configurables
  Configurable<float> svPtMin{"svPtMin", 0.5, "minimum SV pT"};

  // jet level configurables
  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT"};
  Configurable<float> jetPtMax{"jetPtMax", 1000.0, "maximum jet pT"};
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};

  Configurable<float> maxConstSV{"maxConstSV", 999.0, "maximum number of SVs to be stored in the table"};

  Configurable<float> svReductionFactor{"svReductionFactor", 1.0, "factor for how many SVs to keep"};
  Configurable<float> eventReductionFactor{"eventReductionFactor", 0.0, "Percentage of events to be removed"};

  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};

  Configurable<bool> produceTree{"produceTree", true, "produce the jet TTree"};

  Configurable<float> vtxRes{"vtxRes", 0.01, "Vertex position resolution (cluster size) for GNN vertex predictions (cm)"};

  Configurable<int> trainingDatasetRatioParam{"trainingDatasetRatioParam", 0, "Parameter for splitting training/evaluation datasets by collisionId"};

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
    registry.add("h2_nSV_jetpT", "Number of secondary vertices;#it{p}_{T,jet} (GeV/#it{c});nSVs", {HistType::kTH2F, {{200, 0., 200.}, {250, 0, 250.0}}});

    registry.add("h2_SIPs2D_jetpT", "2D IP significance;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
    registry.add("h2_SIPs3D_jetpT", "3D IP significance;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
    registry.add("h2_LxyS_jetpT", "Decay length in XY;#it{p}_{T,jet} (GeV/#it{c});S#it{L}_{xy}", {HistType::kTH2F, {{200, 0., 200.}, {100, 0., 100.0}}});
    registry.add("h2_Dispersion_jetpT", "SV dispersion;#it{p}_{T,jet} (GeV/#it{c});Dispersion", {HistType::kTH2F, {{200, 0., 200.}, {100, 0, 0.5}}});
    registry.add("h2_jetMass_jetpT", "Jet mass;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{jet} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 50.0}}});
    registry.add("h2_SVMass_jetpT", "Secondary vertex mass;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{SV} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 10}}});

    if (doprocessMCJets || doprocessMCJetsForGNN) {
      registry.add("h2_SIPs2D_jetpT_bjet", "2D IP significance b-jets;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_SIPs3D_jetpT_bjet", "3D IP significance b-jets;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_LxyS_jetpT_bjet", "Decay length in XY b-jets;#it{p}_{T,jet} (GeV/#it{c});S#it{L}_{xy}", {HistType::kTH2F, {{200, 0., 200.}, {100, 0., 100.0}}});
      registry.add("h2_Dispersion_jetpT_bjet", "SV dispersion b-jets;#it{p}_{T,jet} (GeV/#it{c});Dispersion", {HistType::kTH2F, {{200, 0., 200.}, {100, 0, 0.5}}});
      registry.add("h2_jetMass_jetpT_bjet", "Jet mass b-jets;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{jet} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 50.0}}});
      registry.add("h2_SVMass_jetpT_bjet", "Secondary vertex mass b-jets;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{SV} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 10.0}}});

      registry.add("h2_SIPs2D_jetpT_cjet", "2D IP significance c-jets;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_SIPs3D_jetpT_cjet", "3D IP significance c-jets;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_LxyS_jetpT_cjet", "Decay length in XY c-jets;#it{p}_{T,jet} (GeV/#it{c});S#it{L}_{xy}", {HistType::kTH2F, {{200, 0., 200.}, {100, 0., 100.0}}});
      registry.add("h2_Dispersion_jetpT_cjet", "SV dispersion c-jets;#it{p}_{T,jet} (GeV/#it{c});Dispersion", {HistType::kTH2F, {{200, 0., 200.}, {100, 0, 0.5}}});
      registry.add("h2_jetMass_jetpT_cjet", "Jet mass c-jets;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{jet} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 50.0}}});
      registry.add("h2_SVMass_jetpT_cjet", "Secondary vertex mass c-jets;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{SV} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 10.0}}});

      registry.add("h2_SIPs2D_jetpT_lfjet", "2D IP significance lf-jet;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_SIPs3D_jetpT_lfjet", "3D IP significance lf-jet;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_LxyS_jetpT_lfjet", "Decay length in XY lf-jet;#it{p}_{T,jet} (GeV/#it{c});S#it{L}_{xy}", {HistType::kTH2F, {{200, 0., 200.}, {100, 0., 100.0}}});
      registry.add("h2_Dispersion_jetpT_lfjet", "SV dispersion lf-jet;#it{p}_{T,jet} (GeV/#it{c});Dispersion", {HistType::kTH2F, {{200, 0., 200.}, {100, 0, 0.5}}});
      registry.add("h2_jetMass_jetpT_lfjet", "Jet mass lf-jet;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{jet} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 50.0}}});
      registry.add("h2_SVMass_jetpT_lfjet", "Secondary vertex mass lf-jet;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{SV} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 10.0}}});

      registry.add("h_jetpT_detector_bjet", "Jet transverse momentum b-jets;#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});
      registry.add("h_jetpT_detector_cjet", "Jet transverse momentum c-jets;#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});
      registry.add("h_jetpT_detector_lfjet", "Jet transverse momentum lf-jet;#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});

      registry.add("h_jetpT_particle_bjet", "Jet transverse momentum particle level b-jets;#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});
      registry.add("h_jetpT_particle_cjet", "Jet transverse momentum particle level c-jets;#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});
      registry.add("h_jetpT_particle_lfjet", "Jet transverse momentum particle level lf-jet;#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});

      registry.add("h2_Response_DetjetpT_PartjetpT_bjet", "Response matrix b-jets;#it{p}_{T,jet}^{det} (GeV/#it{c});#it{p}_{T,jet}^{part} (GeV/#it{c})", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h2_Response_DetjetpT_PartjetpT_cjet", "Response matrix c-jets;#it{p}_{T,jet}^{det} (GeV/#it{c});#it{p}_{T,jet}^{part} (GeV/#it{c})", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h2_Response_DetjetpT_PartjetpT_lfjet", "Response matrix lf-jet;#it{p}_{T,jet}^{det} (GeV/#it{c});#it{p}_{T,jet}^{part} (GeV/#it{c})", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}});
    }

    if (doprocessMCJetsForGNN) {
      //+jet
      registry.add("h_jet_pt", "jet_pt;#it{p}_{T}^{ch jet} (GeV/#it{c});Entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_jet_eta", "jet_eta;#it{#eta}_{ch jet};Entries", {HistType::kTH1F, {{200, -2., 2.}}});
      registry.add("h_jet_phi", "jet_phi;#it{#phi}_{ch jet};Entries", {HistType::kTH1F, {{200, 0., o2::constants::math::TwoPI}}});
      registry.add("h_jet_flav", "jet_flav;jet flavor;Entries", {HistType::kTH1F, {{4, 0., 4.}}});
      registry.add("h_n_trks", "n_trks;#it{n}_{tracks};Entries", {HistType::kTH1F, {{50, 0., 50.}}});
      registry.add("h_jet_mass", "jet_mass;#it{m}_{jet} (GeV/#it{c}^2);Entries", {HistType::kTH1F, {{200, 0., 50.}}});
      auto hJetFlavor = registry.get<TH1>(HIST("h_jet_flav"));
      hJetFlavor->GetXaxis()->SetBinLabel(1, "no mcparticle"); // 0
      hJetFlavor->GetXaxis()->SetBinLabel(2, "c-jet");         // 1
      hJetFlavor->GetXaxis()->SetBinLabel(3, "b-jet");         // 2
      hJetFlavor->GetXaxis()->SetBinLabel(4, "lf-jet");        // 3
      registry.add("h_n_vertices", "n_vertices;#it{n}_{vertex};Entries", {HistType::kTH1F, {{50, 0., 50.}}});
      //+trk
      registry.add("h_trk_pt", "trk_pt;#it{p}_{T} (GeV/#it{c});Entries", {HistType::kTH1F, {{200, 0., 100.}}});
      registry.add("h_trk_eta", "trk_eta;#it{#eta};Entries", {HistType::kTH1F, {{200, -2., 2.}}});
      registry.add("h_trk_phi", "trk_phi;#it{#phi};Entries", {HistType::kTH1F, {{200, 0., o2::constants::math::TwoPI}}});
      registry.add("h_trk_charge", "trk_charge;#it{q};Entries", {HistType::kTH1F, {{3, -1.5, 1.5}}});
      registry.add("h_trk_dcaxy", "trk_dcaxy;#it{DCA}_{xy} (cm);Entries", {HistType::kTH1F, {{200, -0.1, 0.1}}});
      registry.add("h_trk_dcaz", "trk_dcaxyz;#it{DCA}_{z} (cm);Entries", {HistType::kTH1F, {{200, -0.1, 0.1}}});
      registry.add("h_trk_sigmadcaxy", "trk_sigmadcaxy;#it{#sigma}_{#it{DCA}_{xy}} (cm);Entries", {HistType::kTH1F, {{200, 0., 0.1}}});
      registry.add("h_trk_sigmadcaz", "trk_sigmadcaxyz;#it{#sigma}_{#it{DCA}_{z}} (cm);Entries", {HistType::kTH1F, {{200, 0., 0.1}}});
      registry.add("h_trk_itsncls", "trk_itsncls;ITS NCls;Entries", {HistType::kTH1F, {{10, 0., 10.}}});
      registry.add("h_trk_tpcncls", "trk_tpcncls;TPC NCls (Found);Entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_trk_tpcncrs", "trk_tpcncrs;TPC NCrossedRows;Entries", {HistType::kTH1F, {{200, 0., 200.}}});
      registry.add("h_trk_itschi2ncl", "trk_itschi2ncl;ITS #it{#chi}^{2}/ndf;Entries", {HistType::kTH1F, {{200, 0., 20.}}});
      registry.add("h_trk_tpcchi2ncl", "trk_tpcchi2ncl;TPC #it{#chi}^{2}/ndf;Entries", {HistType::kTH1F, {{200, 0., 10.}}});
      registry.add("h2_trk_jtrackpt_vs_origtrackpt", "JTracks::pt vs Tracks::pt", {HistType::kTH2F, {{200, 0., 100.}, {200, 0., 100.}}});
      registry.add("h_trk_vtx_index", "trk_vtx_index;Vertex index;Entries", {HistType::kTH1F, {{20, 0., 20.}}});
      registry.add("h_trk_origin", "trk_origin;Track origin;Entries", {HistType::kTH1F, {{5, 0., 5.}}});
      auto hTrackOrigin = registry.get<TH1>(HIST("h_trk_origin"));
      hTrackOrigin->GetXaxis()->SetBinLabel(1, "NotPhysPrim");
      hTrackOrigin->GetXaxis()->SetBinLabel(2, "Charm");
      hTrackOrigin->GetXaxis()->SetBinLabel(3, "Beauty");
      hTrackOrigin->GetXaxis()->SetBinLabel(4, "Primary");
      hTrackOrigin->GetXaxis()->SetBinLabel(5, "OtherSecondary");
    }
  }

  // FIXME filtering only works when you loop directly over the list, but if you loop over it as a constituent they will not be filtered
  Filter collisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter trackCuts = (aod::jtrack::pt > trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter partCuts = (aod::jmcparticle::pt >= trackPtMin && aod::jmcparticle::pt < trackPtMax);
  Filter jetFilter = (aod::jet::pt >= jetPtMin && aod::jet::pt <= jetPtMax && aod::jet::eta < jetEtaMax - aod::jet::r / 100.f && aod::jet::eta > jetEtaMin + aod::jet::r / 100.f);

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

  // Looping over the SV info and writing them to a table
  template <typename AnalysisJet, typename AnyTracks, typename SecondaryVertices>
  void analyzeJetSVInfo(AnalysisJet const& myJet, AnyTracks const& /*allTracks*/, SecondaryVertices const& /*allSVs*/, std::vector<int>& svIndices, int jetFlavor = 0, double eventweight = 1.0)
  {
    using SVType = typename SecondaryVertices::iterator;

    // Min-heap to store the top 30 SVs by decayLengthXY/errorDecayLengthXY
    auto compare = [](SVType& sv1, SVType& sv2) {
      return (sv1.decayLengthXY() / sv1.errorDecayLengthXY()) > (sv2.decayLengthXY() / sv2.errorDecayLengthXY());
    };

    auto svs = myJet.template secondaryVertices_as<SecondaryVertices>();

    // Sort the SVs based on their decay length significance in descending order
    // This is needed in order to select longest SVs since some jets could have thousands of SVs
    std::sort(svs.begin(), svs.end(), compare);

    for (const auto& candSV : svs) {

      if (candSV.pt() < svPtMin) {
        continue;
      }

      double deltaRJetSV = jetutilities::deltaR(myJet, candSV);
      double massSV = candSV.m();
      double energySV = candSV.e();

      if (svIndices.size() < (svReductionFactor * myJet.template tracks_as<AnyTracks>().size()) && svIndices.size() < maxConstSV) {
        if (produceTree) {
          bjetSVParamsTable(bjetParamsTable.lastIndex() + 1, candSV.pt(), deltaRJetSV, massSV, energySV / myJet.energy(), candSV.impactParameterXY(), candSV.cpa(), candSV.chi2PCA(), candSV.dispersion(), candSV.decayLengthXY(), candSV.errorDecayLengthXY(), candSV.decayLength(), candSV.errorDecayLength());
        }
        svIndices.push_back(bjetSVParamsTable.lastIndex());
      }

      registry.fill(HIST("h2_LxyS_jetpT"), myJet.pt(), candSV.decayLengthXY() / candSV.errorDecayLengthXY(), eventweight);
      registry.fill(HIST("h2_Dispersion_jetpT"), myJet.pt(), candSV.chi2PCA(), eventweight);
      registry.fill(HIST("h2_SVMass_jetpT"), myJet.pt(), massSV, eventweight);

      if (doprocessMCJets) {
        if (jetFlavor == JetTaggingSpecies::beauty) {
          registry.fill(HIST("h2_LxyS_jetpT_bjet"), myJet.pt(), candSV.decayLengthXY() / candSV.errorDecayLengthXY(), eventweight);
          registry.fill(HIST("h2_Dispersion_jetpT_bjet"), myJet.pt(), candSV.chi2PCA(), eventweight);
          registry.fill(HIST("h2_SVMass_jetpT_bjet"), myJet.pt(), massSV, eventweight);
        } else if (jetFlavor == JetTaggingSpecies::charm) {
          registry.fill(HIST("h2_LxyS_jetpT_cjet"), myJet.pt(), candSV.decayLengthXY() / candSV.errorDecayLengthXY(), eventweight);
          registry.fill(HIST("h2_Dispersion_jetpT_cjet"), myJet.pt(), candSV.chi2PCA(), eventweight);
          registry.fill(HIST("h2_SVMass_jetpT_cjet"), myJet.pt(), massSV, eventweight);
        } else {
          registry.fill(HIST("h2_LxyS_jetpT_lfjet"), myJet.pt(), candSV.decayLengthXY() / candSV.errorDecayLengthXY(), eventweight);
          registry.fill(HIST("h2_Dispersion_jetpT_lfjet"), myJet.pt(), candSV.chi2PCA(), eventweight);
          registry.fill(HIST("h2_SVMass_jetpT_lfjet"), myJet.pt(), massSV, eventweight);
        }
      }
    }
  }

  using TrackLabelMap = std::unordered_map<std::string, std::vector<int>>;

  template <typename AnyCollision, typename AnalysisJet, typename AnyTracks, typename SecondaryVertices>
  void analyzeJetTrackInfo(AnyCollision const& /*collision*/, AnalysisJet const& analysisJet, AnyTracks const& /*allTracks*/, SecondaryVertices const& /*allSVs*/, std::vector<int>& trackIndices, int jetFlavor = 0, double eventweight = 1.0)
  {

    for (const auto& constituent : analysisJet.template tracks_as<AnyTracks>()) {

      if (constituent.pt() < trackPtMin || !jettaggingutilities::trackAcceptanceWithDca(constituent, maxIPxy, maxIPz)) {
        continue;
      }

      double deltaRJetTrack = jetutilities::deltaR(analysisJet, constituent);
      double dotProduct = RecoDecay::dotProd(std::array<float, 3>{analysisJet.px(), analysisJet.py(), analysisJet.pz()}, std::array<float, 3>{constituent.px(), constituent.py(), constituent.pz()});
      int sign = jettaggingutilities::getGeoSign(analysisJet, constituent);

      float dRClosestSV = 10.;
      for (const auto& candSV : analysisJet.template secondaryVertices_as<SecondaryVertices>()) {
        double deltaRTrackSV = jetutilities::deltaR(constituent, candSV);
        if (deltaRTrackSV < dRClosestSV) {
          dRClosestSV = deltaRTrackSV;
        }
      }

      registry.fill(HIST("h2_SIPs2D_jetpT"), analysisJet.pt(), sign * std::abs(constituent.dcaXY()) / constituent.sigmadcaXY(), eventweight);
      registry.fill(HIST("h2_SIPs3D_jetpT"), analysisJet.pt(), sign * std::abs(constituent.dcaXYZ()) / constituent.sigmadcaXYZ(), eventweight);

      if (doprocessMCJets) {
        if (jetFlavor == JetTaggingSpecies::beauty) {
          registry.fill(HIST("h2_SIPs2D_jetpT_bjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXY()) / constituent.sigmadcaXY(), eventweight);
          registry.fill(HIST("h2_SIPs3D_jetpT_bjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXYZ()) / constituent.sigmadcaXYZ(), eventweight);
        } else if (jetFlavor == JetTaggingSpecies::charm) {
          registry.fill(HIST("h2_SIPs2D_jetpT_cjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXY()) / constituent.sigmadcaXY(), eventweight);
          registry.fill(HIST("h2_SIPs3D_jetpT_cjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXYZ()) / constituent.sigmadcaXYZ(), eventweight);
        } else {
          registry.fill(HIST("h2_SIPs2D_jetpT_lfjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXY()) / constituent.sigmadcaXY(), eventweight);
          registry.fill(HIST("h2_SIPs3D_jetpT_lfjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXYZ()) / constituent.sigmadcaXYZ(), eventweight);
        }
      }

      if (produceTree) {
        bjetTracksParamsTable(bjetParamsTable.lastIndex() + 1, constituent.pt(), constituent.eta(), dotProduct, dotProduct / analysisJet.p(), deltaRJetTrack, std::abs(constituent.dcaXY()) * sign, constituent.sigmadcaXY(), std::abs(constituent.dcaZ()) * sign, constituent.sigmadcaZ(), constituent.sigmadcaXYZ(), constituent.p() / analysisJet.p(), dRClosestSV);
      }
      trackIndices.push_back(bjetTracksParamsTable.lastIndex());
    }
  }

  template <typename AnyCollision, typename AnalysisJet, typename AnyTracks, typename AnyOriginalTracks>
  void analyzeJetTrackInfoForGNN(AnyCollision const& /*collision*/, AnalysisJet const& analysisJet, AnyTracks const& /*allTracks*/, AnyOriginalTracks const&, std::vector<int>& trackIndices, int jetFlavor = 0, double eventweight = 1.0, TrackLabelMap* trkLabels = nullptr)
  {
    int trkIdx = -1;
    for (const auto& constituent : analysisJet.template tracks_as<AnyTracks>()) {

      trkIdx++;

      if (constituent.pt() < trackPtMin || !jettaggingutilities::trackAcceptanceWithDca(constituent, maxIPxy, maxIPz)) {
        continue;
      }

      double deltaRJetTrack = jetutilities::deltaR(analysisJet, constituent);
      double dotProduct = RecoDecay::dotProd(std::array<float, 3>{analysisJet.px(), analysisJet.py(), analysisJet.pz()}, std::array<float, 3>{constituent.px(), constituent.py(), constituent.pz()});
      int sign = jettaggingutilities::getGeoSign(analysisJet, constituent);

      registry.fill(HIST("h2_SIPs2D_jetpT"), analysisJet.pt(), sign * std::abs(constituent.dcaXY()) / constituent.sigmadcaXY(), eventweight);
      registry.fill(HIST("h2_SIPs3D_jetpT"), analysisJet.pt(), sign * std::abs(constituent.dcaXYZ()) / constituent.sigmadcaXYZ(), eventweight);

      if (doprocessMCJetsForGNN) {
        if (jetFlavor == JetTaggingSpecies::beauty) {
          registry.fill(HIST("h2_SIPs2D_jetpT_bjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXY()) / constituent.sigmadcaXY(), eventweight);
          registry.fill(HIST("h2_SIPs3D_jetpT_bjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXYZ()) / constituent.sigmadcaXYZ(), eventweight);
        } else if (jetFlavor == JetTaggingSpecies::charm) {
          registry.fill(HIST("h2_SIPs2D_jetpT_cjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXY()) / constituent.sigmadcaXY(), eventweight);
          registry.fill(HIST("h2_SIPs3D_jetpT_cjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXYZ()) / constituent.sigmadcaXYZ(), eventweight);
        } else {
          registry.fill(HIST("h2_SIPs2D_jetpT_lfjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXY()) / constituent.sigmadcaXY(), eventweight);
          registry.fill(HIST("h2_SIPs3D_jetpT_lfjet"), analysisJet.pt(), sign * std::abs(constituent.dcaXYZ()) / constituent.sigmadcaXYZ(), eventweight);
        }

        auto origConstit = constituent.template track_as<AnyOriginalTracks>();

        //+
        int trkVtxIndex = 0;
        int trkOrigin = 0;
        if (trkLabels != nullptr) {
          trkVtxIndex = (*trkLabels)["trkVtxIndex"][trkIdx];
          trkOrigin = (*trkLabels)["trkOrigin"][trkIdx];
        }

        //+trk
        registry.fill(HIST("h_trk_pt"), constituent.pt(), eventweight);
        registry.fill(HIST("h_trk_eta"), constituent.eta(), eventweight);
        registry.fill(HIST("h_trk_phi"), origConstit.phi(), eventweight);
        registry.fill(HIST("h_trk_charge"), constituent.sign(), eventweight);
        registry.fill(HIST("h_trk_dcaxy"), std::abs(constituent.dcaXY()) * sign, eventweight);
        registry.fill(HIST("h_trk_dcaz"), std::abs(constituent.dcaZ()) * sign, eventweight);
        registry.fill(HIST("h_trk_sigmadcaxy"), constituent.sigmadcaXY(), eventweight);
        registry.fill(HIST("h_trk_sigmadcaz"), constituent.sigmadcaZ(), eventweight);
        registry.fill(HIST("h_trk_itsncls"), origConstit.itsNCls(), eventweight);
        registry.fill(HIST("h_trk_tpcncls"), origConstit.tpcNClsFound(), eventweight);
        registry.fill(HIST("h_trk_tpcncrs"), origConstit.tpcNClsCrossedRows(), eventweight);
        registry.fill(HIST("h_trk_itschi2ncl"), origConstit.itsChi2NCl(), eventweight);
        registry.fill(HIST("h_trk_tpcchi2ncl"), origConstit.tpcChi2NCl(), eventweight);
        registry.fill(HIST("h2_trk_jtrackpt_vs_origtrackpt"), constituent.pt(), origConstit.pt(), eventweight); // jtrack & new extra table are well joined (linear correlation)
        registry.fill(HIST("h_trk_vtx_index"), trkVtxIndex);
        registry.fill(HIST("h_trk_origin"), trkOrigin);

        if (produceTree) {
          bjetTracksExtraTable(/*bjetParamsTable.lastIndex() + 1, */ origConstit.phi(), constituent.sign(), origConstit.itsChi2NCl(), origConstit.tpcChi2NCl(), origConstit.itsNCls(), origConstit.tpcNClsFound(), origConstit.tpcNClsCrossedRows(), trkOrigin, trkVtxIndex); //+
        }
      }

      if (produceTree) {
        bjetTracksParamsTable(bjetParamsTable.lastIndex() + 1, constituent.pt(), constituent.eta(), dotProduct, dotProduct / analysisJet.p(), deltaRJetTrack, std::abs(constituent.dcaXY()) * sign, constituent.sigmadcaXY(), std::abs(constituent.dcaZ()) * sign, constituent.sigmadcaZ(), constituent.sigmadcaXYZ(), constituent.p() / analysisJet.p(), 0.);
      }
      trackIndices.push_back(bjetTracksParamsTable.lastIndex());
    }
  }

  void processDummy(FilteredCollision::iterator const& /*collision*/)
  {
  }
  PROCESS_SWITCH(BJetTreeCreator, processDummy, "Dummy process function turned on by default", true);

  void processDataJets(FilteredCollision::iterator const& collision, DataJets const& alljets, JetTrackswID const& allTracks, aod::DataSecondaryVertex3Prongs const& allSVs)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) || (static_cast<double>(std::rand()) / RAND_MAX < eventReductionFactor)) {
      return;
    }

    registry.fill(HIST("h_vertexZ"), collision.posZ());

    for (const auto& analysisJet : alljets) {

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

      if (static_cast<double>(std::rand()) / RAND_MAX < getReductionFactor(analysisJet.pt())) {
        continue;
      }

      std::vector<int> indicesTracks;
      std::vector<int> indicesSVs;

      analyzeJetSVInfo(analysisJet, allTracks, allSVs, indicesSVs);
      analyzeJetTrackInfo(collision, analysisJet, allTracks, allSVs, indicesTracks);

      registry.fill(HIST("h2_jetMass_jetpT"), analysisJet.pt(), analysisJet.mass());

      int nSVs = analysisJet.template secondaryVertices_as<aod::DataSecondaryVertex3Prongs>().size();

      registry.fill(HIST("h2_nTracks_jetpT"), analysisJet.pt(), indicesTracks.size());
      registry.fill(HIST("h2_nSV_jetpT"), analysisJet.pt(), nSVs < 250 ? nSVs : 249);

      if (produceTree) {
        bjetConstituentsTable(bjetParamsTable.lastIndex() + 1, indicesTracks, indicesSVs);
        bjetParamsTable(analysisJet.pt(), analysisJet.eta(), analysisJet.phi(), indicesTracks.size(), nSVs, analysisJet.mass(), 0, analysisJet.r());
      }
    }
  }
  PROCESS_SWITCH(BJetTreeCreator, processDataJets, "jet information in Data", false);

  using MCDJetTable = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::MCDSecondaryVertex3ProngIndices>>;
  using MCPJetTable = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>;
  using FilteredCollisionMCD = soa::Filtered<soa::Join<aod::JetCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>>;

  Preslice<aod::JMcParticles> mcParticlesPerCollision = aod::jmcparticle::mcCollisionId;
  Preslice<MCPJetTable> mcpJetsPerCollision = aod::jet::mcCollisionId;

  void processMCJets(FilteredCollisionMCD::iterator const& collision, MCDJetTable const& MCDjets, MCPJetTable const& MCPjets, JetTracksMCDwID const& allTracks, aod::JetParticles const& MCParticles, aod::MCDSecondaryVertex3Prongs const& allSVs)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) || (static_cast<double>(std::rand()) / RAND_MAX < eventReductionFactor)) {
      return;
    }

    float eventWeight = collision.weight();

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

      float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
      if (analysisJet.pt() > pTHatMaxMCD * pTHat) {
        continue;
      }

      std::vector<int> indicesTracks;
      std::vector<int> indicesSVs;

      int16_t jetFlavor = 0;

      // JetTracksMCDwID::iterator hftrack;
      // jetFlavor = jettaggingutilities::mcdJetFromHFShower(analysisJet, allTracks, mcParticlesPerColl, (float)(analysisJet.r() / 100.));
      // jetFlavor = jettaggingutilities::jetTrackFromHFShower(analysisJet, nonFilteredTracks, mcParticlesPerColl, hftrack);

      for (const auto& mcpjet : analysisJet.template matchedJetGeo_as<MCPJetTable>()) {
        if (useQuarkDef) {
          jetFlavor = jettaggingutilities::getJetFlavor(mcpjet, mcParticlesPerColl);
        } else {
          jetFlavor = jettaggingutilities::getJetFlavorHadron(mcpjet, mcParticlesPerColl);
          // jetFlavor = jettaggingutilities::mcpJetFromHFShower(mcpjet, mcParticlesPerColl, (float)(mcpjet.r() / 100.));
        }
      }

      if ((jetFlavor != JetTaggingSpecies::charm && jetFlavor != JetTaggingSpecies::beauty) && (static_cast<double>(std::rand()) / RAND_MAX < getReductionFactor(analysisJet.pt()))) {
        continue;
      }

      analyzeJetSVInfo(analysisJet, allTracks, allSVs, indicesSVs, jetFlavor, eventWeight);
      analyzeJetTrackInfo(collision, analysisJet, allTracks, allSVs, indicesTracks, jetFlavor, eventWeight);

      int nSVs = analysisJet.template secondaryVertices_as<aod::MCDSecondaryVertex3Prongs>().size();

      registry.fill(HIST("h2_jetMass_jetpT"), analysisJet.pt(), analysisJet.mass(), eventWeight);

      registry.fill(HIST("h2_nTracks_jetpT"), analysisJet.pt(), indicesTracks.size());
      registry.fill(HIST("h2_nSV_jetpT"), analysisJet.pt(), nSVs < 250 ? nSVs : 249);

      if (jetFlavor == JetTaggingSpecies::beauty) {
        registry.fill(HIST("h2_jetMass_jetpT_bjet"), analysisJet.pt(), analysisJet.mass(), eventWeight);
        registry.fill(HIST("h_jetpT_detector_bjet"), analysisJet.pt(), eventWeight);
      } else if (jetFlavor == JetTaggingSpecies::charm) {
        registry.fill(HIST("h2_jetMass_jetpT_cjet"), analysisJet.pt(), analysisJet.mass(), eventWeight);
        registry.fill(HIST("h_jetpT_detector_cjet"), analysisJet.pt(), eventWeight);
      } else {
        registry.fill(HIST("h2_jetMass_jetpT_lfjet"), analysisJet.pt(), analysisJet.mass(), eventWeight);
        registry.fill(HIST("h_jetpT_detector_lfjet"), analysisJet.pt(), eventWeight);
      }

      for (const auto& mcpjet : analysisJet.template matchedJetGeo_as<MCPJetTable>()) {

        if (mcpjet.pt() > pTHatMaxMCP * pTHat) {
          continue;
        }

        if (jetFlavor == JetTaggingSpecies::beauty) {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_bjet"), analysisJet.pt(), mcpjet.pt(), eventWeight);
        } else if (jetFlavor == JetTaggingSpecies::charm) {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_cjet"), analysisJet.pt(), mcpjet.pt(), eventWeight);
        } else {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_lfjet"), analysisJet.pt(), mcpjet.pt(), eventWeight);
        }
      }

      if (produceTree) {
        bjetConstituentsTable(bjetParamsTable.lastIndex() + 1, indicesTracks, indicesSVs);
        bjetParamsTable(analysisJet.pt(), analysisJet.eta(), analysisJet.phi(), indicesTracks.size(), nSVs, analysisJet.mass(), jetFlavor, analysisJet.r());
        bjetParamsExtraTable(eventWeight);
      }
    }
  }
  PROCESS_SWITCH(BJetTreeCreator, processMCJets, "jet information in MC", false);

  using MCDJetTableNoSV = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCDetectorLevelJetFlavourDef>>;
  using JetParticleswID = soa::Join<aod::JetParticles, aod::JMcParticlePIs>;

  void processMCJetsForGNN(FilteredCollisionMCD::iterator const& collision, aod::JMcCollisions const&, MCDJetTableNoSV const& MCDjets, MCPJetTable const& MCPjets, JetTracksMCDwID const& allTracks, JetParticleswID const& MCParticles, OriginalTracks const& origTracks, aod::McParticles const& origParticles)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits) || (static_cast<double>(std::rand()) / RAND_MAX < eventReductionFactor)) {
      return;
    }

    // Uses only collisionId % trainingDatasetRaioParam == 0 for training dataset
    if (trainingDatasetRatioParam && collision.collisionId() % trainingDatasetRatioParam != 0) {
      return;
    }

    float eventWeight = collision.weight();

    registry.fill(HIST("h_vertexZ"), collision.posZ(), eventWeight);

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

      int16_t jetFlavor = analysisJet.origin();

      if ((jetFlavor != JetTaggingSpecies::charm && jetFlavor != JetTaggingSpecies::beauty) && (static_cast<double>(std::rand()) / RAND_MAX < getReductionFactor(analysisJet.pt()))) {
        continue;
      }

      //+
      TrackLabelMap trkLabels{{"trkVtxIndex", {}}, {"trkOrigin", {}}};
      int nVertices = jettaggingutilities::vertexClustering(collision.template mcCollision_as<aod::JMcCollisions>(), analysisJet, allTracks, MCParticles, origParticles, trkLabels, true, vtxRes, trackPtMin);
      analyzeJetTrackInfoForGNN(collision, analysisJet, allTracks, origTracks, indicesTracks, jetFlavor, eventWeight, &trkLabels);

      registry.fill(HIST("h2_jetMass_jetpT"), analysisJet.pt(), analysisJet.mass(), eventWeight);
      registry.fill(HIST("h2_nTracks_jetpT"), analysisJet.pt(), indicesTracks.size(), eventWeight);

      //+jet
      registry.fill(HIST("h_jet_pt"), analysisJet.pt());
      registry.fill(HIST("h_jet_eta"), analysisJet.eta());
      registry.fill(HIST("h_jet_phi"), analysisJet.phi());
      registry.fill(HIST("h_jet_flav"), jetFlavor);
      registry.fill(HIST("h_n_trks"), indicesTracks.size());
      registry.fill(HIST("h_jet_mass"), analysisJet.mass());
      registry.fill(HIST("h_n_vertices"), nVertices);

      if (jetFlavor == JetTaggingSpecies::beauty) {
        registry.fill(HIST("h2_jetMass_jetpT_bjet"), analysisJet.pt(), analysisJet.mass(), eventWeight);
        registry.fill(HIST("h_jetpT_detector_bjet"), analysisJet.pt(), eventWeight);
      } else if (jetFlavor == JetTaggingSpecies::charm) {
        registry.fill(HIST("h2_jetMass_jetpT_cjet"), analysisJet.pt(), analysisJet.mass(), eventWeight);
        registry.fill(HIST("h_jetpT_detector_cjet"), analysisJet.pt(), eventWeight);
      } else {
        registry.fill(HIST("h2_jetMass_jetpT_lfjet"), analysisJet.pt(), analysisJet.mass(), eventWeight);
        registry.fill(HIST("h_jetpT_detector_lfjet"), analysisJet.pt(), eventWeight);
      }

      if (produceTree) {
        std::vector<int> indicesSVs;
        bjetConstituentsTable(bjetParamsTable.lastIndex() + 1, indicesTracks, indicesSVs);
        bjetParamsTable(analysisJet.pt(), analysisJet.eta(), analysisJet.phi(), indicesTracks.size(), nVertices, analysisJet.mass(), jetFlavor, analysisJet.r());
      }
    }
  }
  PROCESS_SWITCH(BJetTreeCreator, processMCJetsForGNN, "jet information in MC for GNN", false);

  Filter mccollisionFilter = nabs(aod::jmccollision::posZ) < vertexZCut;
  using FilteredCollisionMCP = soa::Filtered<aod::JMcCollisions>;

  void processMCTruthJets(FilteredCollisionMCP::iterator const& collision, MCPJetTable const& MCPjets, aod::JetParticles const& MCParticles)
  {
    float eventWeight = collision.weight();

    for (const auto& mcpjet : MCPjets) {

      bool jetIncluded = false;
      for (const auto& jetR : jetRadiiValues) {
        if (mcpjet.r() == static_cast<int>(jetR * 100)) {
          jetIncluded = true;
          break;
        }
      }

      if (!jetIncluded) {
        continue;
      }

      float pTHat = 10. / (std::pow(eventWeight, 1.0 / pTHatExponent));
      if (mcpjet.pt() > pTHatMaxMCP * pTHat) {
        continue;
      }

      int16_t jetFlavor = 0;
      if (useQuarkDef) {
        jetFlavor = jettaggingutilities::getJetFlavor(mcpjet, MCParticles);
      } else {
        jetFlavor = jettaggingutilities::getJetFlavorHadron(mcpjet, MCParticles);
        // jetFlavor = jettaggingutilities::mcpJetFromHFShower(mcpjet, MCParticles, (float)(mcpjet.r() / 100.));
      }

      if (jetFlavor == JetTaggingSpecies::beauty) {
        registry.fill(HIST("h_jetpT_particle_bjet"), mcpjet.pt(), eventWeight);
      } else if (jetFlavor == JetTaggingSpecies::charm) {
        registry.fill(HIST("h_jetpT_particle_cjet"), mcpjet.pt(), eventWeight);
      } else {
        registry.fill(HIST("h_jetpT_particle_lfjet"), mcpjet.pt(), eventWeight);
      }
    }
  }
  PROCESS_SWITCH(BJetTreeCreator, processMCTruthJets, "truth jet information", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<BJetTreeCreator>(cfgc, TaskName{"bjet-tree-creator"})};
}
