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

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetTaggingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetTagging.h"

#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/RecoDecay.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace jetInfo
{
// DECLARE_SOA_INDEX_COLUMN(JetIndex, jetindex); //! The jet index
DECLARE_SOA_COLUMN(JetpT, jetpt, float);       //! jet pT
DECLARE_SOA_COLUMN(JetEta, jeteta, float);     //! jet eta
DECLARE_SOA_COLUMN(JetPhi, jetphi, float);     //! jet phi
DECLARE_SOA_COLUMN(NTracks, nTracks, int16_t); //! number of charged tracks inside the jet
DECLARE_SOA_COLUMN(NSV, nSV, int16_t);         //! Number of secondary vertices in the jet
DECLARE_SOA_COLUMN(NSD, nSD, int16_t);         //! Number of softdrop splittings in the jet
DECLARE_SOA_COLUMN(Zg, zg, float);             //! The momentum assymetry between the splittings
DECLARE_SOA_COLUMN(Rg, rg, float);             //! The grooming radius
DECLARE_SOA_COLUMN(JetMass, mass, float);      //! The jet mass
DECLARE_SOA_COLUMN(JetFlavor, jetFl, int16_t); //! The jet flavor (b, c, or lf)
DECLARE_SOA_COLUMN(JetR, jetR, int16_t);       //! The jet radius
} // namespace jetInfo

DECLARE_SOA_TABLE(bjetParams, "AOD", "BJETPARAMS",
                  o2::soa::Index<>,
                  jetInfo::JetpT,
                  jetInfo::JetEta,
                  jetInfo::JetPhi,
                  jetInfo::NTracks,
                  jetInfo::NSV,
                  jetInfo::NSD,
                  jetInfo::Zg,
                  jetInfo::Rg,
                  jetInfo::JetMass,
                  jetInfo::JetFlavor,
                  jetInfo::JetR);

using bjetParam = bjetParams::iterator;

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
DECLARE_SOA_COLUMN(SignedIP3D, ip3d, float);                           //! The track signed 3D IP
DECLARE_SOA_COLUMN(SignedIP3DSign, ip3dsigma, float);                  //! The track signed 3D IP significance
DECLARE_SOA_COLUMN(MomFraction, momfraction, float);                   //! The track momentum fraction of the jets
DECLARE_SOA_COLUMN(DeltaRTrackVertex, rtrackvertex, float);            //! DR between the track and the closest SV, to be decided whether to add to or not
// DECLARE_SOA_COLUMN(DCATrackJet, dcatrackjet, float);                              //! The distance between track and jet, unfortunately it cannot be calculated in O2
} // namespace trackInfo

DECLARE_SOA_TABLE(bjetTracksParams, "AOD", "BJETTRACKPARAMS",
                  o2::soa::Index<>,
                  trackInfo::bjetParamId,
                  trackInfo::TrackpT,
                  trackInfo::TrackEta,
                  trackInfo::DotProdTrackJet,
                  trackInfo::DotProdTrackJetOverJet,
                  trackInfo::DeltaRJetTrack,
                  trackInfo::SignedIP2D,
                  trackInfo::SignedIP2DSign,
                  trackInfo::SignedIP3D,
                  trackInfo::SignedIP3DSign,
                  trackInfo::MomFraction,
                  trackInfo::DeltaRTrackVertex);

using bjetTracksParam = bjetTracksParams::iterator;

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
DECLARE_SOA_COLUMN(DecayLength2D, lxy, float);            //! The decay length of the SV in XY
DECLARE_SOA_COLUMN(DecayLength2DError, lxysigma, float);  //! The decay length of the SV in XY significance
DECLARE_SOA_COLUMN(DecayLength3D, lxyz, float);           //! The decay length of the SV in 3D
DECLARE_SOA_COLUMN(DecayLength3DError, lxyzsigma, float); //! The decay length of the SV in 3d significance
// DECLARE_SOA_COLUMN(SVDispersion, svdispersion, float);                              //! The SV dispersion, unfortunately it cannot be calculated in O2
} // namespace SVInfo

DECLARE_SOA_TABLE(bjetSVParams, "AOD", "BJETSVPARAMS",
                  o2::soa::Index<>,
                  SVInfo::bjetParamId,
                  SVInfo::SVpT,
                  SVInfo::DeltaRSVJet,
                  SVInfo::SVMass,
                  SVInfo::SVfE,
                  SVInfo::IPXY,
                  SVInfo::CPA,
                  SVInfo::Chi2PCA,
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
  Produces<aod::bjetTracksParams> bjetTracksParamsTable;
  Produces<aod::bjetSVParams> bjetSVParamsTable;
  Produces<aod::bjetConstituents> bjetConstituentsTable;

  HistogramRegistry registry;

  // event level configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};

  // track level configurables
  Configurable<float> trackPtMin{"trackPtMin", 0.5, "minimum track pT"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track eta"};

  // track level configurables
  Configurable<float> svPtMin{"svPtMin", 0.5, "minimum SV pT"};

  // jet level configurables
  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT"};
  Configurable<float> jetPtMax{"jetPtMax", 1000.0, "maximum jet pT"};
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};

  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};

  Configurable<float> zCut{"zCut", 0.1, "soft drop z cut"};
  Configurable<float> beta{"beta", 0.0, "soft drop beta"};

  Configurable<bool> produceTree{"produceTree", true, "produce the jet TTree"};

  int eventSelection = -1;

  std::vector<fastjet::PseudoJet> jetConstituents;

  std::vector<double> jetRadiiValues;

  void init(InitContext const&)
  {
    jetRadiiValues = (std::vector<double>)jetRadii;

    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));

    registry.add("h_vertexZ", "Vertex Z;#it{Z} (cm)", {HistType::kTH1F, {{40, -20.0, 20.0}}});

    registry.add("h2_nTracks_jetpT", "Number of tracks;#it{p}_{T,jet} (GeV/#it{c});nTracks", {HistType::kTH2F, {{200, 0., 200.}, {100, 0, 100.0}}});
    registry.add("h2_nSV_jetpT", "Number of secondary vertices;#it{p}_{T,jet} (GeV/#it{c});nSVs", {HistType::kTH2F, {{200, 0., 200.}, {250, 0, 250.0}}});

    registry.add("h2_SIPs2D_jetpT", "2D IP significance;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
    registry.add("h2_SIPs3D_jetpT", "3D IP significance;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
    registry.add("h2_LxyS_jetpT", "Decay length in XY;#it{p}_{T,jet} (GeV/#it{c});S#it{L}_{xy}", {HistType::kTH2F, {{200, 0., 200.}, {100, 0., 100.0}}});
    registry.add("h2_Dispersion_jetpT", "SV dispersion;#it{p}_{T,jet} (GeV/#it{c});Dispersion", {HistType::kTH2F, {{200, 0., 200.}, {100, 0, 50.0}}});
    registry.add("h2_jetMass_jetpT", "Jet mass;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{jet} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 50.0}}});
    registry.add("h2_SVMass_jetpT", "Secondary vertex mass;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{SV} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 10}}});

    if (doprocessMCJets) {
      registry.add("h2_SIPs2D_jetpT_bjet", "2D IP significance b-jets;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_SIPs3D_jetpT_bjet", "3D IP significance b-jets;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_LxyS_jetpT_bjet", "Decay length in XY b-jets;#it{p}_{T,jet} (GeV/#it{c});S#it{L}_{xy}", {HistType::kTH2F, {{200, 0., 200.}, {100, 0., 100.0}}});
      registry.add("h2_Dispersion_jetpT_bjet", "SV dispersion b-jets;#it{p}_{T,jet} (GeV/#it{c});Dispersion", {HistType::kTH2F, {{200, 0., 200.}, {100, 0, 50.0}}});
      registry.add("h2_jetMass_jetpT_bjet", "Jet mass b-jets;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{jet} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 50.0}}});
      registry.add("h2_SVMass_jetpT_bjet", "Secondary vertex mass b-jets;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{SV} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 10.0}}});

      registry.add("h2_SIPs2D_jetpT_cjet", "2D IP significance c-jets;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_SIPs3D_jetpT_cjet", "3D IP significance c-jets;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_LxyS_jetpT_cjet", "Decay length in XY c-jets;#it{p}_{T,jet} (GeV/#it{c});S#it{L}_{xy}", {HistType::kTH2F, {{200, 0., 200.}, {100, 0., 100.0}}});
      registry.add("h2_Dispersion_jetpT_cjet", "SV dispersion c-jets;#it{p}_{T,jet} (GeV/#it{c});Dispersion", {HistType::kTH2F, {{200, 0., 200.}, {100, 0, 50.0}}});
      registry.add("h2_jetMass_jetpT_cjet", "Jet mass c-jets;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{jet} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 50.0}}});
      registry.add("h2_SVMass_jetpT_cjet", "Secondary vertex mass c-jets;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{SV} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 10.0}}});

      registry.add("h2_SIPs2D_jetpT_lfjet", "2D IP significance lf-jet;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_SIPs3D_jetpT_lfjet", "3D IP significance lf-jet;#it{p}_{T,jet} (GeV/#it{c});IPs", {HistType::kTH2F, {{200, 0., 200.}, {100, -50.0, 50.0}}});
      registry.add("h2_LxyS_jetpT_lfjet", "Decay length in XY lf-jet;#it{p}_{T,jet} (GeV/#it{c});S#it{L}_{xy}", {HistType::kTH2F, {{200, 0., 200.}, {100, 0., 100.0}}});
      registry.add("h2_Dispersion_jetpT_lfjet", "SV dispersion lf-jet;#it{p}_{T,jet} (GeV/#it{c});Dispersion", {HistType::kTH2F, {{200, 0., 200.}, {100, 0, 50.0}}});
      registry.add("h2_jetMass_jetpT_lfjet", "Jet mass lf-jet;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{jet} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 50.0}}});
      registry.add("h2_SVMass_jetpT_lfjet", "Secondary vertex mass lf-jet;#it{p}_{T,jet} (GeV/#it{c});#it{m}_{SV} (GeV/#it{c}^{2})", {HistType::kTH2F, {{200, 0., 200.}, {50, 0, 10.0}}});

      registry.add("h_jetpT_detector_bjet", "Jet transverse momentum b-jets;#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});
      registry.add("h_jetpT_detector_cjet", "Jet transverse momentum c-jets;#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});
      registry.add("h_jetpT_detector_lfjet", "Jet transverse momentum lf-jet;#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});

      registry.add("h_jetpT_particle_bjet", "Jet transverse momentum particle level b-jets;#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});
      registry.add("h_jetpT_particle_cjet", "Jet transverse momentum particle level c-jets;#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});
      registry.add("h_jetpT_particle_lfjet", "Jet transverse momentum particle level lf-jet;#it{p}_{T,jet} (GeV/#it{c})", {HistType::kTH1F, {{200, 0., 200.0}}});

      registry.add("h2_Zg_jetpT_particle_bjet", "jet #it{Z}_{g} particle level b-jets;#it{p}_{T,jet} (GeV/#it{c});#it{Z}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {10, 0.0, 0.5}}});
      registry.add("h2_Zg_jetpT_particle_cjet", "jet #it{Z}_{g} particle level c-jets;#it{p}_{T,jet} (GeV/#it{c});#it{Z}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {10, 0.0, 0.5}}});
      registry.add("h2_Zg_jetpT_particle_lfjet", "jet #it{Z}_{g} particle level lf-jet;#it{p}_{T,jet} (GeV/#it{c});#it{Z}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {10, 0.0, 0.5}}});

      registry.add("h2_Rg_jetpT_particle_bjet", "jet #it{R}_{g} particle level b-jets;#it{p}_{T,jet} (GeV/#it{c});#it{R}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {10, 0.0, 0.5}}});
      registry.add("h2_Rg_jetpT_particle_cjet", "jet #it{R}_{g} particle level c-jets;#it{p}_{T,jet} (GeV/#it{c});#it{R}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {10, 0.0, 0.5}}});
      registry.add("h2_Rg_jetpT_particle_lfjet", "jet #it{R}_{g} particle level lf-jet;#it{p}_{T,jet} (GeV/#it{c});#it{R}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {10, 0.0, 0.5}}});

      registry.add("h2_nSD_jetpT_particle_bjet", "Jet #it{n}_{SD} particle level b-jets;#it{p}_{T,jet} (GeV/#it{c});#it{n}_{SD}", {HistType::kTH2F, {{200, 0., 200.}, {7, -0.5, 6.5}}});
      registry.add("h2_nSD_jetpT_particle_cjet", "Jet #it{n}_{SD} particle level c-jets;#it{p}_{T,jet} (GeV/#it{c});#it{n}_{SD}", {HistType::kTH2F, {{200, 0., 200.}, {7, -0.5, 6.5}}});
      registry.add("h2_nSD_jetpT_particle_lfjet", "Jet #it{n}_{SD} particle level lf-jet;#it{p}_{T,jet} (GeV/#it{c});#it{n}_{SD}", {HistType::kTH2F, {{200, 0., 200.}, {7, -0.5, 6.5}}});

      registry.add("h2_Zg_jetpT_detector_bjet", "jet #it{Z}_{g} b-jets;#it{p}_{T,jet} (GeV/#it{c});#it{Z}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {10, 0.0, 0.5}}});
      registry.add("h2_Zg_jetpT_detector_cjet", "jet #it{Z}_{g} c-jets;#it{p}_{T,jet} (GeV/#it{c});#it{Z}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {10, 0.0, 0.5}}});
      registry.add("h2_Zg_jetpT_detector_lfjet", "jet #it{Z}_{g} lf-jet;#it{p}_{T,jet} (GeV/#it{c});#it{Z}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {10, 0.0, 0.5}}});

      registry.add("h2_Rg_jetpT_detector_bjet", "jet #it{R}_{g} b-jets;#it{p}_{T,jet} (GeV/#it{c});#it{R}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {10, 0.0, 0.5}}});
      registry.add("h2_Rg_jetpT_detector_cjet", "jet #it{R}_{g} c-jets;#it{p}_{T,jet} (GeV/#it{c});#it{R}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {10, 0.0, 0.5}}});
      registry.add("h2_Rg_jetpT_detector_lfjet", "jet #it{R}_{g} lf-jet;#it{p}_{T,jet} (GeV/#it{c});#it{R}_{g}", {HistType::kTH2F, {{200, 0., 200.}, {10, 0.0, 0.5}}});

      registry.add("h2_nSD_jetpT_detector_bjet", "Jet #it{n}_{SD} b-jets;#it{p}_{T,jet} (GeV/#it{c});#it{n}_{SD}", {HistType::kTH2F, {{200, 0., 200.}, {7, -0.5, 6.5}}});
      registry.add("h2_nSD_jetpT_detector_cjet", "Jet #it{n}_{SD} c-jets;#it{p}_{T,jet} (GeV/#it{c});#it{n}_{SD}", {HistType::kTH2F, {{200, 0., 200.}, {7, -0.5, 6.5}}});
      registry.add("h2_nSD_jetpT_detector_lfjet", "Jet #it{n}_{SD} lf-jet;#it{p}_{T,jet} (GeV/#it{c});#it{n}_{SD}", {HistType::kTH2F, {{200, 0., 200.}, {7, -0.5, 6.5}}});

      registry.add("h2_Response_DetjetpT_PartjetpT_bjet", "Response matrix b-jets;#it{p}_{T,jet}^{det} (GeV/#it{c});#it{p}_{T,jet}^{part} (GeV/#it{c})", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h2_Response_DetjetpT_PartjetpT_cjet", "Response matrix c-jets;#it{p}_{T,jet}^{det} (GeV/#it{c});#it{p}_{T,jet}^{part} (GeV/#it{c})", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}});
      registry.add("h2_Response_DetjetpT_PartjetpT_lfjet", "Response matrix lf-jet;#it{p}_{T,jet}^{det} (GeV/#it{c});#it{p}_{T,jet}^{part} (GeV/#it{c})", {HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}}});
    }
  }

  // FIXME filtering only works when you loop directly over the list, but if you loop over it as a constituent they will not be filtered
  Filter collisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter trackCuts = (aod::jtrack::pt > trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter partCuts = (aod::jmcparticle::pt >= trackPtMin && aod::jmcparticle::pt < trackPtMax);
  Filter jetFilter = (aod::jet::pt >= jetPtMin && aod::jet::pt <= jetPtMax && aod::jet::eta < jetEtaMax - aod::jet::r / 100.f && aod::jet::eta > jetEtaMin + aod::jet::r / 100.f);

  using FilteredCollision = soa::Filtered<soa::Join<aod::JCollisions, aod::JCollisionPIs>>;
  using JetTrackswID = soa::Filtered<soa::Join<JetTracks, aod::JTrackPIs>>;
  using JetTracksMCDwID = soa::Filtered<soa::Join<JetTracksMCD, aod::JTrackPIs>>;
  using OriginalTracks = soa::Join<aod::Tracks, aod::TracksCov, aod::TrackSelection, aod::TracksDCA, aod::TracksDCACov>;
  using DataJets = soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents, aod::DataSecondaryVertex3ProngIndices>>;

  template <typename T>
  std::tuple<double, double, double> jetReclustering(T const& jet)
  {
    std::vector<fastjet::PseudoJet> jetReclustered;
    fastjet::ClusterSequence clusterSeq(jetConstituents, fastjet::JetDefinition(fastjet::JetAlgorithm::cambridge_algorithm, (jet.r() / 100.f) * 5.0));
    jetReclustered = sorted_by_pt(clusterSeq.inclusive_jets());
    fastjet::PseudoJet daughterSubJet = jetReclustered[0];
    fastjet::PseudoJet parentSubJet1;
    fastjet::PseudoJet parentSubJet2;
    bool softDropped = false;
    double nSD = 0.0;
    double Zg = -1.0;
    double Rg = -1.0;
    while (daughterSubJet.has_parents(parentSubJet1, parentSubJet2)) {
      if (parentSubJet1.perp() < parentSubJet2.perp()) {
        std::swap(parentSubJet1, parentSubJet2);
      }

      // delta_R_vec.push_back(delta_R);
      // xkt_vec.push_back(xkt);
      // z_vec.push_back(z);

      auto z = parentSubJet2.perp() / (parentSubJet1.perp() + parentSubJet2.perp());
      auto theta = parentSubJet1.delta_R(parentSubJet2); // TODO add the lund plane
      // double xkt = parentSubJet2.perp() * sin(theta);  // TODO add the lund plane
      if (z >= zCut * TMath::Power(theta / (jet.r() / 100.f), beta)) {
        if (!softDropped) {
          Zg = z;
          Rg = theta;
          softDropped = true;
        }
        nSD++;
      }
      daughterSubJet = parentSubJet1;
    }
    return std::make_tuple(nSD, Zg, Rg);
  }

  // Looping over the SV info and writing them to a table
  template <typename AnalysisJet, typename SecondaryVertices>
  void analyzeJetSVInfo(AnalysisJet const& myJet, SecondaryVertices const& allSVs, std::vector<int>& svIndices, int jetFlavor = 0, double eventweight = 1.0)
  {
    for (const auto& candSV : myJet.template secondaryVertices_as<SecondaryVertices>()) {

      if (candSV.pt() < svPtMin) {
        continue;
      }

      double deltaRJetSV = jetutilities::deltaR(myJet, candSV);
      double massSV = candSV.m();
      double energySV = candSV.e();

      if (produceTree) {
        bjetSVParamsTable(bjetParamsTable.lastIndex(), candSV.pt(), deltaRJetSV, massSV, energySV / myJet.energy(), candSV.impactParameterXY(), candSV.cpa(), candSV.chi2PCA(), candSV.decayLengthXY(), candSV.errorDecayLengthXY(), candSV.decayLength(), candSV.errorDecayLength());
      }
      svIndices.push_back(bjetSVParamsTable.lastIndex());

      registry.fill(HIST("h2_LxyS_jetpT"), myJet.pt(), candSV.decayLengthXY() / candSV.errorDecayLengthXY(), eventweight);
      registry.fill(HIST("h2_Dispersion_jetpT"), myJet.pt(), candSV.chi2PCA(), eventweight);
      registry.fill(HIST("h2_SVMass_jetpT"), myJet.pt(), massSV, eventweight);

      if (doprocessMCJets) {
        if (jetFlavor == 2) {
          registry.fill(HIST("h2_LxyS_jetpT_bjet"), myJet.pt(), candSV.decayLengthXY() / candSV.errorDecayLengthXY(), eventweight);
          registry.fill(HIST("h2_Dispersion_jetpT_bjet"), myJet.pt(), candSV.chi2PCA(), eventweight);
          registry.fill(HIST("h2_SVMass_jetpT_bjet"), myJet.pt(), massSV, eventweight);
        } else if (jetFlavor == 1) {
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

  template <typename AnyCollision, typename AnalysisJet, typename AnyTracks, typename SecondaryVertices>
  void analyzeJetTrackInfo(AnyCollision const& collision, AnalysisJet const& analysisJet, AnyTracks const& allTracks, SecondaryVertices const& allSVs, std::vector<int>& trackIndices, int jetFlavor = 0, double eventweight = 1.0)
  {
    jetConstituents.clear();

    for (auto& jconstituent : analysisJet.template tracks_as<AnyTracks>()) {
      fastjetutilities::fillTracks(jconstituent, jetConstituents, jconstituent.globalIndex());

      if (jconstituent.pt() < trackPtMin) {
        continue;
      }

      auto constituent = jconstituent.template track_as<OriginalTracks>();
      double deltaRJetTrack = jetutilities::deltaR(analysisJet, constituent);
      double dotProduct = RecoDecay::dotProd(std::array<float, 3>{analysisJet.px(), analysisJet.py(), analysisJet.pz()}, std::array<float, 3>{constituent.px(), constituent.py(), constituent.pz()});
      int sign = jettaggingutilities::getGeoSign(collision, analysisJet, constituent);

      float RClosestSV = 10.;
      for (const auto& candSV : analysisJet.template secondaryVertices_as<SecondaryVertices>()) {
        double deltaRTrackSV = jetutilities::deltaR(constituent, candSV);
        if (deltaRTrackSV < RClosestSV) {
          RClosestSV = deltaRTrackSV;
        }
      }

      float dcaXYZ(0.), sigmaDcaXYZ2(0.);
      dcaXYZ = getDcaXYZ(constituent, &sigmaDcaXYZ2);
      // jettaggingutilities::calculateDcaXYZ(dcaXYZ, sigmaDcaXYZ2, constituent.dcaXY(), constituent.dcaZ(), constituent.cYY(), constituent.cZY(), constituent.cZZ(), constituent.sigmaDcaXY2(), constituent.sigmaDcaZ2());

      registry.fill(HIST("h2_SIPs2D_jetpT"), analysisJet.pt(), sign * TMath::Abs(constituent.dcaXY()) / TMath::Sqrt(constituent.sigmaDcaXY2()), eventweight);
      registry.fill(HIST("h2_SIPs3D_jetpT"), analysisJet.pt(), sign * dcaXYZ / TMath::Sqrt(sigmaDcaXYZ2), eventweight);

      if (doprocessMCJets) {
        if (jetFlavor == 2) {
          registry.fill(HIST("h2_SIPs2D_jetpT_bjet"), analysisJet.pt(), sign * TMath::Abs(constituent.dcaXY()) / TMath::Sqrt(constituent.sigmaDcaXY2()), eventweight);
          registry.fill(HIST("h2_SIPs3D_jetpT_bjet"), analysisJet.pt(), sign * dcaXYZ / TMath::Sqrt(sigmaDcaXYZ2), eventweight);
        } else if (jetFlavor == 1) {
          registry.fill(HIST("h2_SIPs2D_jetpT_cjet"), analysisJet.pt(), sign * TMath::Abs(constituent.dcaXY()) / TMath::Sqrt(constituent.sigmaDcaXY2()), eventweight);
          registry.fill(HIST("h2_SIPs3D_jetpT_cjet"), analysisJet.pt(), sign * dcaXYZ / TMath::Sqrt(sigmaDcaXYZ2), eventweight);
        } else {
          registry.fill(HIST("h2_SIPs2D_jetpT_lfjet"), analysisJet.pt(), sign * TMath::Abs(constituent.dcaXY()) / TMath::Sqrt(constituent.sigmaDcaXY2()), eventweight);
          registry.fill(HIST("h2_SIPs3D_jetpT_lfjet"), analysisJet.pt(), sign * dcaXYZ / TMath::Sqrt(sigmaDcaXYZ2), eventweight);
        }
      }

      if (produceTree) {
        bjetTracksParamsTable(bjetParamsTable.lastIndex(), constituent.pt(), constituent.eta(), dotProduct, dotProduct / analysisJet.p(), deltaRJetTrack, TMath::Abs(constituent.dcaXY()) * sign, TMath::Sqrt(constituent.sigmaDcaXY2()), dcaXYZ * sign, TMath::Sqrt(sigmaDcaXYZ2), constituent.p() / analysisJet.p(), RClosestSV);
      }
      trackIndices.push_back(bjetTracksParamsTable.lastIndex());
    }
  }

  void processDummy(FilteredCollision::iterator const& collision)
  {
  }
  PROCESS_SWITCH(BJetTreeCreator, processDummy, "Dummy process function turned on by default", true);

  void processDataJets(FilteredCollision::iterator const& collision, DataJets const& alljets, JetTrackswID const& allTracks, OriginalTracks const& allOrigTracks, aod::DataSecondaryVertex3Prongs const& allSVs)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }

    registry.fill(HIST("h_vertexZ"), collision.posZ());

    for (const auto& analysisJet : alljets) {

      bool jetIncluded = false;
      for (auto jetR : jetRadiiValues) {
        if (analysisJet.r() == static_cast<int>(jetR * 100)) {
          jetIncluded = true;
          break;
        }
      }

      if (!jetIncluded) {
        continue;
      }

      std::vector<int> tracksIndices;
      std::vector<int> SVsIndices;

      analyzeJetSVInfo(analysisJet, allSVs, SVsIndices);
      analyzeJetTrackInfo(collision, analysisJet, allTracks, allSVs, tracksIndices);

      registry.fill(HIST("h2_jetMass_jetpT"), analysisJet.pt(), analysisJet.mass());

      registry.fill(HIST("h2_nTracks_jetpT"), analysisJet.pt(), tracksIndices.size());
      registry.fill(HIST("h2_nSV_jetpT"), analysisJet.pt(), SVsIndices.size() < 250 ? SVsIndices.size() : 249);

      auto [nSD, Zg, Rg] = jetReclustering(analysisJet);

      if (produceTree) {
        bjetConstituentsTable(bjetParamsTable.lastIndex(), tracksIndices, SVsIndices);
        bjetParamsTable(analysisJet.pt(), analysisJet.eta(), analysisJet.phi(), tracksIndices.size(), SVsIndices.size(), nSD, Zg, Rg, analysisJet.mass(), 0, analysisJet.r());
      }
    }
  }
  PROCESS_SWITCH(BJetTreeCreator, processDataJets, "jet information in Data", false);

  using MCDJetTable = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::MCDSecondaryVertex3ProngIndices, aod::ChargedMCDetectorLevelJetEventWeights>>;
  using MCPJetTable = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets, aod::ChargedMCParticleLevelJetEventWeights>>;
  using FilteredCollisionMCD = soa::Filtered<soa::Join<aod::JCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>>;

  Preslice<aod::JMcParticles> McParticlesPerCollision = aod::jmcparticle::mcCollisionId;
  Preslice<MCPJetTable> McPJetsPerCollision = aod::jet::mcCollisionId;

  void processMCJets(FilteredCollisionMCD::iterator const& collision, MCDJetTable const& MCDjets, MCPJetTable const& MCPjets, JetTracksMCDwID const& allTracks, JetParticles const& MCParticles, aod::MCDSecondaryVertex3Prongs const& allSVs, OriginalTracks const& origTracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }

    registry.fill(HIST("h_vertexZ"), collision.posZ());

    auto const mcParticlesPerColl = MCParticles.sliceBy(McParticlesPerCollision, collision.mcCollisionId());
    auto const mcPJetsPerColl = MCPjets.sliceBy(McPJetsPerCollision, collision.mcCollisionId());

    for (const auto& analysisJet : MCDjets) {

      bool jetIncluded = false;
      for (auto jetR : jetRadiiValues) {
        if (analysisJet.r() == static_cast<int>(jetR * 100)) {
          jetIncluded = true;
          break;
        }
      }

      if (!jetIncluded) {
        continue;
      }

      std::vector<int> tracksIndices;
      std::vector<int> SVsIndices;

      float eventWeight = analysisJet.eventWeight();
      int16_t jetFlavor = 0;

      // JetTracksMCDwID::iterator hftrack;
      // jetFlavor = jettaggingutilities::mcdJetFromHFShower(analysisJet, allTracks, mcParticlesPerColl, (float)(analysisJet.r() / 100.));
      // jetFlavor = jettaggingutilities::jetTrackFromHFShower(analysisJet, nonFilteredTracks, mcParticlesPerColl, hftrack);

      for (auto& mcpjet : analysisJet.template matchedJetGeo_as<MCPJetTable>()) {
        jetFlavor = jettaggingutilities::getJetFlavor(mcpjet, mcParticlesPerColl);
        // jetFlavor = jettaggingutilities::mcpJetFromHFShower(mcpjet, mcParticlesPerColl, (float)(mcpjet.r() / 100.));
      }
      analyzeJetSVInfo(analysisJet, allSVs, SVsIndices, jetFlavor, eventWeight);
      analyzeJetTrackInfo(collision, analysisJet, allTracks, allSVs, tracksIndices, jetFlavor, eventWeight);

      auto [nSD, Zg, Rg] = jetReclustering(analysisJet);

      registry.fill(HIST("h2_jetMass_jetpT"), analysisJet.pt(), analysisJet.mass(), eventWeight);

      registry.fill(HIST("h2_nTracks_jetpT"), analysisJet.pt(), tracksIndices.size());
      registry.fill(HIST("h2_nSV_jetpT"), analysisJet.pt(), SVsIndices.size() < 250 ? SVsIndices.size() : 249);

      if (jetFlavor == 2) {
        registry.fill(HIST("h2_jetMass_jetpT_bjet"), analysisJet.pt(), analysisJet.mass(), eventWeight);
        registry.fill(HIST("h_jetpT_detector_bjet"), analysisJet.pt(), eventWeight);
        registry.fill(HIST("h2_Zg_jetpT_detector_bjet"), analysisJet.pt(), Zg, eventWeight);
        registry.fill(HIST("h2_Rg_jetpT_detector_bjet"), analysisJet.pt(), Rg, eventWeight);
        registry.fill(HIST("h2_nSD_jetpT_detector_bjet"), analysisJet.pt(), nSD, eventWeight);
      } else if (jetFlavor == 1) {
        registry.fill(HIST("h2_jetMass_jetpT_cjet"), analysisJet.pt(), analysisJet.mass(), eventWeight);
        registry.fill(HIST("h_jetpT_detector_cjet"), analysisJet.pt(), eventWeight);
        registry.fill(HIST("h2_Zg_jetpT_detector_cjet"), analysisJet.pt(), Zg, eventWeight);
        registry.fill(HIST("h2_Rg_jetpT_detector_cjet"), analysisJet.pt(), Rg, eventWeight);
        registry.fill(HIST("h2_nSD_jetpT_detector_cjet"), analysisJet.pt(), nSD, eventWeight);
      } else {
        registry.fill(HIST("h2_jetMass_jetpT_lfjet"), analysisJet.pt(), analysisJet.mass(), eventWeight);
        registry.fill(HIST("h_jetpT_detector_lfjet"), analysisJet.pt(), eventWeight);
        registry.fill(HIST("h2_Zg_jetpT_detector_lfjet"), analysisJet.pt(), Zg, eventWeight);
        registry.fill(HIST("h2_Rg_jetpT_detector_lfjet"), analysisJet.pt(), Rg, eventWeight);
        registry.fill(HIST("h2_nSD_jetpT_detector_lfjet"), analysisJet.pt(), nSD, eventWeight);
      }

      for (auto& mcpjet : analysisJet.template matchedJetGeo_as<MCPJetTable>()) {
        if (jetFlavor == 2) {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_bjet"), analysisJet.pt(), mcpjet.pt(), eventWeight);
        } else if (jetFlavor == 1) {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_cjet"), analysisJet.pt(), mcpjet.pt(), eventWeight);
        } else {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_lfjet"), analysisJet.pt(), mcpjet.pt(), eventWeight);
        }
      }

      if (produceTree) {
        bjetConstituentsTable(bjetParamsTable.lastIndex(), tracksIndices, SVsIndices);
        bjetParamsTable(analysisJet.pt(), analysisJet.eta(), analysisJet.phi(), tracksIndices.size(), SVsIndices.size(), nSD, Zg, Rg, analysisJet.mass(), jetFlavor, analysisJet.r());
      }
    }

    for (const auto& mcpjet : mcPJetsPerColl) {

      bool jetIncluded = false;
      for (auto jetR : jetRadiiValues) {
        if (mcpjet.r() == static_cast<int>(jetR * 100)) {
          jetIncluded = true;
          break;
        }
      }

      if (!jetIncluded) {
        continue;
      }

      int16_t jetFlavor = 0;
      jetFlavor = jettaggingutilities::getJetFlavor(mcpjet, mcParticlesPerColl);
      // jetFlavor = jettaggingutilities::mcpJetFromHFShower(mcpjet, mcParticlesPerColl, (float)(mcpjet.r() / 100.));

      float eventWeight = mcpjet.eventWeight();
      jetConstituents.clear();
      for (auto& constituent : mcpjet.template tracks_as<JetParticles>()) {
        fastjetutilities::fillTracks(constituent, jetConstituents, constituent.globalIndex());
      }

      auto [nSD, Zg, Rg] = jetReclustering(mcpjet);

      if (jetFlavor == 2) {
        registry.fill(HIST("h_jetpT_particle_bjet"), mcpjet.pt(), eventWeight);
        registry.fill(HIST("h2_Zg_jetpT_particle_bjet"), mcpjet.pt(), Zg, eventWeight);
        registry.fill(HIST("h2_Rg_jetpT_particle_bjet"), mcpjet.pt(), Rg, eventWeight);
        registry.fill(HIST("h2_nSD_jetpT_particle_bjet"), mcpjet.pt(), nSD, eventWeight);
      } else if (jetFlavor == 1) {
        registry.fill(HIST("h_jetpT_particle_cjet"), mcpjet.pt(), eventWeight);
        registry.fill(HIST("h2_Zg_jetpT_particle_cjet"), mcpjet.pt(), Zg, eventWeight);
        registry.fill(HIST("h2_Rg_jetpT_particle_cjet"), mcpjet.pt(), Rg, eventWeight);
        registry.fill(HIST("h2_nSD_jetpT_particle_cjet"), mcpjet.pt(), nSD, eventWeight);
      } else {
        registry.fill(HIST("h_jetpT_particle_lfjet"), mcpjet.pt(), eventWeight);
        registry.fill(HIST("h2_Zg_jetpT_particle_lfjet"), mcpjet.pt(), Zg, eventWeight);
        registry.fill(HIST("h2_Rg_jetpT_particle_lfjet"), mcpjet.pt(), Rg, eventWeight);
        registry.fill(HIST("h2_nSD_jetpT_particle_lfjet"), mcpjet.pt(), nSD, eventWeight);
      }
    }
  }
  PROCESS_SWITCH(BJetTreeCreator, processMCJets, "jet information in MC", false);

  Filter mccollisionFilter = nabs(aod::jmccollision::posZ) < vertexZCut;
  using FilteredCollisionMCP = soa::Filtered<aod::JMcCollisions>;

  void processMCTruthJets(FilteredCollisionMCP::iterator const& collision, MCPJetTable const& MCPjets, JetParticles const& MCParticles)
  {

    for (const auto& mcpjet : MCPjets) {

      bool jetIncluded = false;
      for (auto jetR : jetRadiiValues) {
        if (mcpjet.r() == static_cast<int>(jetR * 100)) {
          jetIncluded = true;
          break;
        }
      }

      if (!jetIncluded) {
        continue;
      }

      int16_t jetFlavor = 0;
      jetFlavor = jettaggingutilities::getJetFlavor(mcpjet, MCParticles);
      // jetFlavor = jettaggingutilities::mcpJetFromHFShower(mcpjet, mcParticlesPerColl, (float)(mcpjet.r() / 100.));

      float eventWeight = mcpjet.eventWeight();
      jetConstituents.clear();
      for (auto& constituent : mcpjet.template tracks_as<JetParticles>()) {
        fastjetutilities::fillTracks(constituent, jetConstituents, constituent.globalIndex());
      }

      auto [nSD, Zg, Rg] = jetReclustering(mcpjet);

      if (jetFlavor == 2) {
        registry.fill(HIST("h_jetpT_particle_bjet"), mcpjet.pt(), eventWeight);
        registry.fill(HIST("h2_Zg_jetpT_particle_bjet"), mcpjet.pt(), Zg, eventWeight);
        registry.fill(HIST("h2_Rg_jetpT_particle_bjet"), mcpjet.pt(), Rg, eventWeight);
        registry.fill(HIST("h2_nSD_jetpT_particle_bjet"), mcpjet.pt(), nSD, eventWeight);
      } else if (jetFlavor == 1) {
        registry.fill(HIST("h_jetpT_particle_cjet"), mcpjet.pt(), eventWeight);
        registry.fill(HIST("h2_Zg_jetpT_particle_cjet"), mcpjet.pt(), Zg, eventWeight);
        registry.fill(HIST("h2_Rg_jetpT_particle_cjet"), mcpjet.pt(), Rg, eventWeight);
        registry.fill(HIST("h2_nSD_jetpT_particle_cjet"), mcpjet.pt(), nSD, eventWeight);
      } else {
        registry.fill(HIST("h_jetpT_particle_lfjet"), mcpjet.pt(), eventWeight);
        registry.fill(HIST("h2_Zg_jetpT_particle_lfjet"), mcpjet.pt(), Zg, eventWeight);
        registry.fill(HIST("h2_Rg_jetpT_particle_lfjet"), mcpjet.pt(), Rg, eventWeight);
        registry.fill(HIST("h2_nSD_jetpT_particle_lfjet"), mcpjet.pt(), nSD, eventWeight);
      }
    }
  }
  PROCESS_SWITCH(BJetTreeCreator, processMCTruthJets, "truth jet information", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<BJetTreeCreator>(cfgc, TaskName{"bjet-tree-creator"})};
}
