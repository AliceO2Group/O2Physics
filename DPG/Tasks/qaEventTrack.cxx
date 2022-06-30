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
/// \file   qaEventTrack.cxx
/// \author Peter Hristov <Peter.Hristov@cern.ch>, CERN
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Henrique J C Zanoli <henrique.zanoli@cern.ch>, Utrecht University
/// \author Mario Krüger <mario.kruger@cern.ch>
/// \author Nicolò Jacazio <nicolo.jacazio@cern.ch>, CERN
/// \brief  Task to produce QA objects for the track and the event properties in the AOD.
///         This task can also be configured to produce a table with reduced information used for correlation studies for track selection
///

#include "qaEventTrack.h"

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisDataModel.h"
#include "ReconstructionDataFormats/DCA.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::dataformats;

// TODO: add PID wagons as dependency + include impact parameter studies (same or separate task in workflow??)

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
// Task declaration
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
struct qaEventTrack {
  // Tables to produce
  Produces<o2::aod::DPGCollisions> tableCollisions;
  Produces<o2::aod::DPGTracks> tableTracks;
  Produces<o2::aod::DPGRecoParticles> tableRecoParticles;
  Produces<o2::aod::DPGNonRecoParticles> tableNonRecoParticles;

  // general steering settings
  Configurable<bool> isRun3{"isRun3", false, "Is Run3 dataset"}; // TODO: derive this from metadata once possible to get rid of the flag

  // options to select specific events
  Configurable<bool> selectGoodEvents{"selectGoodEvents", true, "select good events"};
  // selection specific to the table creation workflow
  Configurable<float> selectMaxVtxZ{"selectMaxVtxZ", 100.f, "Derived data option: select collision in a given Z window"};
  Configurable<int> targetNumberOfEvents{"targetNumberOfEvents", 10000000, "Derived data option: target number of collisions, if the target is met, future collisions will be skipped"};
  Configurable<float> fractionOfSampledEvents{"fractionOfSampledEvents", 1.f, "Derived data option: fraction of events to sample"};

  // options to select only specific tracks
  Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  Configurable<int> selectCharge{"selectCharge", 0, "select charge +1 or -1 (0 means no selection)"};
  Configurable<bool> selectPrim{"selectPrim", false, "select primaries"};
  Configurable<bool> selectSec{"selectSec", false, "select secondaries"};
  Configurable<int> selectPID{"selectPID", 0, "select pid"};
  Configurable<float> minPt{"minPt", -10.f, "Minimum pt of accepted tracks"};
  Configurable<float> maxPt{"maxPt", 1e10f, "Maximum pt of accepted tracks"};
  Configurable<float> minEta{"minEta", -2.f, "Minimum eta of accepted tracks"};
  Configurable<float> maxEta{"maxEta", 2.0f, "Maximum eta of accepted tracks"};
  Configurable<float> minPhi{"minPhi", -1.f, "Minimum phi of accepted tracks"};
  Configurable<float> maxPhi{"maxPhi", 10.f, "Maximum phi of accepted tracks"};

  // configurable binning of histograms
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 5.0, 10.0, 20.0, 50.0}, ""};

  ConfigurableAxis binsVertexPosZ{"binsVertexPosZ", {100, -20., 20.}, ""}; // TODO: do we need this to be configurable?
  ConfigurableAxis binsVertexPosXY{"binsVertexPosXY", {500, -1., 1.}, ""}; // TODO: do we need this to be configurable?
  ConfigurableAxis binsTrackMultiplicity{"binsTrackMultiplcity", {200, 0, 200}, ""};

  // TODO: ask if one can have different filters for both process functions
  Filter trackFilter = (trackSelection.node() == 0) ||
                       ((trackSelection.node() == 1) && requireGlobalTrackInFilter()) ||
                       ((trackSelection.node() == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
                       ((trackSelection.node() == 3) && requireTrackCutInFilter(TrackSelectionFlags::kQualityTracks) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks)) ||
                       ((trackSelection.node() == 4) && requireTrackCutInFilter(TrackSelectionFlags::kQualityTracks)) ||
                       ((trackSelection.node() == 5) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks));

  using TrackIUTable = soa::Join<aod::TracksIU, aod::TrackSelection>;
  Partition<TrackIUTable> tracksIUFiltered = (trackSelection.node() == 0) ||
                                             ((trackSelection.node() == 1) && requireGlobalTrackInFilter()) ||
                                             ((trackSelection.node() == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
                                             ((trackSelection.node() == 3) && requireTrackCutInFilter(TrackSelectionFlags::kQualityTracks) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks)) ||
                                             ((trackSelection.node() == 4) && requireTrackCutInFilter(TrackSelectionFlags::kQualityTracks)) ||
                                             ((trackSelection.node() == 5) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks));

  HistogramRegistry histos;

  void init(InitContext const&);

  // Function to select tracks
  template <bool IS_MC, typename T>
  bool isSelectedTrack(const T& track)
  {
    if (track.pt() < minPt || track.pt() > maxPt) { // Extra pT selection
      return false;
    }
    if (track.eta() < minEta || track.eta() > maxEta) { // Extra Eta selection
      return false;
    }
    if (track.phi() < minPhi || track.phi() > maxPhi) { // Extra Phi selection
      return false;
    }
    if (selectCharge && (selectCharge != track.sign())) {
      return false;
    }
    if constexpr (IS_MC) {
      if (!track.has_mcParticle()) {
        if (selectPrim || selectSec || selectPID) {
          return false;
        } else {
          return true;
        }
      }
      auto particle = track.mcParticle();
      const bool isPrimary = particle.isPhysicalPrimary();
      if (selectPrim && !isPrimary) {
        return false;
      }
      if (selectSec && isPrimary) {
        return false;
      }
      if (selectPID && selectPID != std::abs(particle.pdgCode())) {
        return false;
      }
    }
    return true;
  }

  // Function to select collisions
  template <bool doFill, typename T>
  bool isSelectedCollision(const T& collision)
  {
    if constexpr (doFill) {
      histos.fill(HIST("Events/recoEff"), 1);
    }
    if (selectGoodEvents && !(isRun3 ? collision.sel8() : collision.sel7())) { // currently only sel8 is defined for run3
      return false;
    }
    if constexpr (doFill) {
      histos.fill(HIST("Events/recoEff"), 2);
    }
    return true;
  }

  // General function to fill data and MC histograms
  template <bool IS_MC, typename C, typename T>
  void fillRecoHistograms(const C& collision, const T& tracks, aod::FullTracks const& tracksUnfiltered);

  // Process function for data
  using CollisionTableData = soa::Join<aod::Collisions, aod::EvSels>;
  using TrackTableData = soa::Join<aod::FullTracks, aod::TracksCov, aod::TracksDCA, aod::TrackSelection>;
  void processData(CollisionTableData::iterator const& collision, soa::Filtered<TrackTableData> const& tracks, aod::FullTracks const& tracksUnfiltered)
  {
    fillRecoHistograms<false>(collision, tracks, tracksUnfiltered);
  }
  PROCESS_SWITCH(qaEventTrack, processData, "process data", false);

  // Process function for IU vs DCA track comparison
  void processDataIU(CollisionTableData::iterator const& collision,
                     aod::FullTracks const& tracksUnfiltered, aod::TracksIU const& tracksIU)
  {
    if (!isSelectedCollision<false>(collision)) {
      return;
    }

    int trackIndex = 0;
    for (const auto& trk : tracksUnfiltered) {
      if (!isSelectedTrack<false>(trk)) {
        continue;
      }

      const auto& trkIU = tracksIU.iteratorAt(trackIndex++);
      histos.fill(HIST("Tracks/IU/Pt"), trkIU.pt());
      histos.fill(HIST("Tracks/IU/Eta"), trkIU.eta());
      histos.fill(HIST("Tracks/IU/Phi"), trkIU.phi());

      histos.fill(HIST("Tracks/IU/alpha"), trkIU.alpha());
      histos.fill(HIST("Tracks/IU/x"), trkIU.x());
      histos.fill(HIST("Tracks/IU/y"), trkIU.y());
      histos.fill(HIST("Tracks/IU/z"), trkIU.z());
      histos.fill(HIST("Tracks/IU/signed1Pt"), trkIU.signed1Pt());
      histos.fill(HIST("Tracks/IU/snp"), trkIU.snp());
      histos.fill(HIST("Tracks/IU/tgl"), trkIU.tgl());

      histos.fill(HIST("Tracks/IU/deltaDCA/Pt"), trk.pt(), trkIU.pt() - trk.pt());
      histos.fill(HIST("Tracks/IU/deltaDCA/Eta"), trk.eta(), trkIU.eta() - trk.eta());
      histos.fill(HIST("Tracks/IU/deltaDCA/Phi"), trk.phi(), trkIU.phi() - trk.phi());

      histos.fill(HIST("Tracks/IU/vsDCA/Pt"), trk.pt(), trkIU.pt());
      histos.fill(HIST("Tracks/IU/vsDCA/Eta"), trk.eta(), trkIU.eta());
      histos.fill(HIST("Tracks/IU/vsDCA/Phi"), trk.phi(), trkIU.phi());
    }
  }
  PROCESS_SWITCH(qaEventTrack, processDataIU, "process IU vs DCA comparison", true);

  // Process function for filtered IU
  void processDataIUFiltered(CollisionTableData::iterator const& collision, TrackIUTable const&)
  {
    if (!isSelectedCollision<false>(collision)) {
      return;
    }

    auto tracksIU = tracksIUFiltered->sliceByCached(aod::track::collisionId, collision.globalIndex());

    // int trackIndex = 0;
    for (const auto& trkIU : tracksIU) {
      if (!isSelectedTrack<false>(trkIU)) {
        continue;
      }

      // const auto& trkIU = tracksIU.iteratorAt(trackIndex++);
      histos.fill(HIST("Tracks/IUFiltered/Pt"), trkIU.pt());
      histos.fill(HIST("Tracks/IUFiltered/Eta"), trkIU.eta());
      histos.fill(HIST("Tracks/IUFiltered/Phi"), trkIU.phi());

      histos.fill(HIST("Tracks/IUFiltered/alpha"), trkIU.alpha());
      histos.fill(HIST("Tracks/IUFiltered/x"), trkIU.x());
      histos.fill(HIST("Tracks/IUFiltered/y"), trkIU.y());
      histos.fill(HIST("Tracks/IUFiltered/z"), trkIU.z());
      histos.fill(HIST("Tracks/IUFiltered/signed1Pt"), trkIU.signed1Pt());
      histos.fill(HIST("Tracks/IUFiltered/snp"), trkIU.snp());
      histos.fill(HIST("Tracks/IUFiltered/tgl"), trkIU.tgl());
    }
  }
  PROCESS_SWITCH(qaEventTrack, processDataIUFiltered, "process IU filtered", true);

  // Process function for MC
  using CollisionTableMC = soa::Join<CollisionTableData, aod::McCollisionLabels>;
  using TrackTableMC = soa::Join<TrackTableData, aod::McTrackLabels>;
  void processMC(CollisionTableMC::iterator const& collision, soa::Filtered<TrackTableMC> const& tracks, aod::FullTracks const& tracksUnfiltered,
                 aod::McParticles const& mcParticles, aod::McCollisions const& mcCollisions)
  {
    fillRecoHistograms<true>(collision, tracks, tracksUnfiltered);
  }
  PROCESS_SWITCH(qaEventTrack, processMC, "process mc", true); // FIXME: would like to disable this by default and swich on via --processMC but currently this crashes -> ask experts

  // Process functions for skimming data
  void processTableData(CollisionTableData::iterator const& collision,
                        soa::Filtered<soa::Join<TrackTableData, aod::TOFSignal, aod::TOFEvTime>> const& tracks,
                        aod::BCs const& bcs)
  {
    fillDerivedTable<false>(collision, tracks, 0, bcs);
  }
  PROCESS_SWITCH(qaEventTrack, processTableData, "Process data for table producing", false);

  void processTableMC(CollisionTableMC::iterator const& collision,
                      soa::Filtered<soa::Join<TrackTableMC, aod::TOFSignal, aod::TOFEvTime>> const& tracks,
                      aod::McParticles const& mcParticles,
                      aod::McCollisions const& mcCollisions,
                      aod::BCs const& bcs)
  {
    fillDerivedTable<true>(collision, tracks, mcParticles, bcs);
  }
  PROCESS_SWITCH(qaEventTrack, processTableMC, "Process MC for table producing", false);

  //**************************************************************************************************
  /**
   * Fill reco level tables.
   */
  //**************************************************************************************************
  int nTableEventCounter = 0; // Number of processed events
  template <bool IS_MC, typename C, typename T, typename P>
  void fillDerivedTable(const C& collision, const T& tracks, const P& particles, const aod::BCs&)
  {
    if (!isSelectedCollision<false>(collision)) {
      return;
    }
    if (abs(collision.posZ()) > selectMaxVtxZ) {
      return;
    }
    if (fractionOfSampledEvents < 1.f && (static_cast<float>(rand()) / static_cast<float>(RAND_MAX)) > fractionOfSampledEvents) { // Skip events that are not sampled
      return;
    }
    if (nTableEventCounter > targetNumberOfEvents) { // Skip events if target is reached
      return;
    }
    nTableEventCounter++;

    tableCollisions(collision.posZ(),
                    (isRun3 ? collision.sel8() : collision.sel7()),
                    collision.bc().runNumber());
    int nTracks = 0;
    int particleProduction = 0;

    for (const auto& track : tracks) {
      if (!isSelectedTrack<IS_MC>(track)) {
        continue;
      }
      ++nTracks;
    }
    tableTracks.reserve(nTracks);
    std::vector<int64_t> recoPartIndices(IS_MC ? nTracks : 0);

    if constexpr (IS_MC) { // Running only on MC
      tableRecoParticles.reserve(nTracks);
    }
    int64_t iTrack = 0;
    for (const auto& track : tracks) {
      if (!isSelectedTrack<IS_MC>(track)) {
        continue;
      }
      tableTracks(tableCollisions.lastIndex(),
                  track.pt(), track.eta(), track.phi(), track.pt() * std::sqrt(track.c1Pt21Pt2()),
                  track.flags(), track.sign(),
                  track.dcaXY(), track.dcaZ(), track.length(),
                  track.itsClusterMap(),
                  track.itsChi2NCl(), track.tpcChi2NCl(), track.trdChi2(), track.tofChi2(),
                  track.hasITS(), track.hasTPC(), track.hasTRD(), track.hasTOF(),
                  track.tpcNClsFound(), track.tpcNClsCrossedRows(),
                  track.tpcCrossedRowsOverFindableCls(), track.tpcFoundOverFindableCls(), track.tpcFractionSharedCls(),
                  track.itsNCls(), track.itsNClsInnerBarrel(), track.tpcSignal(), track.tofSignal() - track.tofEvTime());

      if constexpr (IS_MC) { // Running only on MC
        if (track.has_mcParticle()) {
          auto particle = track.mcParticle();
          recoPartIndices[iTrack++] = particle.globalIndex();
          if (particle.isPhysicalPrimary()) {
            particleProduction = 0;
          } else if (particle.getProcess() == 4) {
            particleProduction = 1;
          } else {
            particleProduction = 2;
          }
          tableRecoParticles(particle.pt(), particle.eta(), particle.phi(),
                             particle.pdgCode(), particleProduction);
        } else { // If it does not have the particle we fill with the track values and tag it with -1 in the production
          tableRecoParticles(track.pt(), track.eta(), track.phi(),
                             0, -1);
        }
      }
    }

    // Running only on MC
    if constexpr (IS_MC) {
      if (!collision.has_mcCollision()) {
        return;
      }
      const auto& particlesInCollision = particles.sliceBy(aod::mcparticle::mcCollisionId, collision.mcCollision().globalIndex());
      tableNonRecoParticles.reserve(particlesInCollision.size() - nTracks);
      for (const auto& particle : particlesInCollision) {
        const auto partReconstructed = std::find(recoPartIndices.begin(), recoPartIndices.end(), particle.globalIndex()) != recoPartIndices.end();
        if (partReconstructed) {
          continue;
        }
        if (particle.isPhysicalPrimary()) {
          particleProduction = 0;
        } else if (particle.getProcess() == 4) {
          particleProduction = 1;
        } else {
          particleProduction = 2;
        }
        tableNonRecoParticles(tableCollisions.lastIndex(),
                              particle.pt(), particle.eta(), particle.phi(),
                              particle.pdgCode(), particleProduction,
                              particle.vx(), particle.vy(), particle.vz());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<qaEventTrack>(cfgc)};
}

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
// Task implementation
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------

//**************************************************************************************************
/**
 * Initialize the task.
 */
//**************************************************************************************************
void qaEventTrack::init(InitContext const&)
{
  if (doprocessTableData == true && doprocessTableMC == true) {
    LOGF(fatal, "Cannot enable processTableData and processTableMC at the same time. Please choose one.");
  }
  if (!doprocessData && !doprocessMC) {
    LOGF(info, "No enabled QA, all histograms are disabled");
    return;
  }
  const AxisSpec axisPt{binsPt, "#it{p}_{T} [GeV/c]"};
  const AxisSpec axisEta{180, -0.9, 0.9, "#it{#eta}"};
  const AxisSpec axisPhi{180, 0., 2 * M_PI, "#it{#varphi} [rad]"};
  const AxisSpec axisVertexNumContrib{200, 0, 200, "Number Of contributors to the PV"};
  const AxisSpec axisVertexPosX{binsVertexPosXY, "X [cm]"};
  const AxisSpec axisVertexPosY{binsVertexPosXY, "Y [cm]"};
  const AxisSpec axisVertexPosZ{binsVertexPosZ, "Z [cm]"};
  const AxisSpec axisVertexCov{100, -0.005, 0.005};
  const AxisSpec axisVertexPosReso{100, -0.5, 0.5};
  const AxisSpec axisTrackMultiplicity{binsTrackMultiplicity, "Track Multiplicity"};
  const AxisSpec axisParX{200, -0.36, 0.36, "#it{x} [cm]"};
  const AxisSpec axisParY{200, -0.5, 0.5, "#it{y} [cm]"};
  const AxisSpec axisParZ{200, -11., 11., "#it{z} [cm]"};
  const AxisSpec axisParAlpha{36, -M_PI, M_PI, "#alpha [rad]"};
  const AxisSpec axisParSigned1Pt{200, -8, 8, "#it{q}/#it{p}_{T}"};
  const AxisSpec axisParSnp{11, -0.1, 0.1, "snp"};
  const AxisSpec axisParTgl{200, -1., 1., "tgl"};

  const AxisSpec axisDeltaPt{100, -0.5, 0.5, "#it{p}_{T, rec} - #it{p}_{T, gen}"};
  const AxisSpec axisDeltaEta{100, -0.1, 0.1, "#eta_{rec} - #eta_{gen}"};
  const AxisSpec axisDeltaPhi{100, -0.1, 0.1, "#varphi_{rec} - #varphi_{gen}"};

  // collision
  auto eventRecoEffHist = histos.add<TH1>("Events/recoEff", "", kTH1D, {{2, 0.5, 2.5}});
  eventRecoEffHist->GetXaxis()->SetBinLabel(1, "all");
  eventRecoEffHist->GetXaxis()->SetBinLabel(2, "selected");
  histos.add("Events/posX", "", kTH1D, {axisVertexPosX});
  histos.add("Events/posY", "", kTH1D, {axisVertexPosY});
  histos.add("Events/posZ", "", kTH1D, {axisVertexPosZ});
  histos.add("Events/posXY", "", kTH2D, {axisVertexPosX, axisVertexPosY});
  histos.add("Events/posXvsNContrib", "", kTH2D, {axisVertexPosX, axisVertexNumContrib});
  histos.add("Events/posYvsNContrib", "", kTH2D, {axisVertexPosY, axisVertexNumContrib});
  histos.add("Events/posZvsNContrib", "", kTH2D, {axisVertexPosZ, axisVertexNumContrib});
  histos.add("Events/nContrib", "", kTH1D, {axisVertexNumContrib});
  histos.add("Events/nContribVsMult", "", kTH2D, {axisVertexNumContrib, axisTrackMultiplicity});
  histos.add("Events/vertexChi2", ";#chi^{2}", kTH1D, {{100, 0, 100}});

  histos.add("Events/covXX", ";Cov_{xx} [cm^{2}]", kTH1D, {axisVertexCov});
  histos.add("Events/covXY", ";Cov_{xy} [cm^{2}]", kTH1D, {axisVertexCov});
  histos.add("Events/covXZ", ";Cov_{xz} [cm^{2}]", kTH1D, {axisVertexCov});
  histos.add("Events/covYY", ";Cov_{yy} [cm^{2}]", kTH1D, {axisVertexCov});
  histos.add("Events/covYZ", ";Cov_{yz} [cm^{2}]", kTH1D, {axisVertexCov});
  histos.add("Events/covZZ", ";Cov_{zz} [cm^{2}]", kTH1D, {axisVertexCov});

  histos.add("Events/nTracks", "", kTH1D, {axisTrackMultiplicity});

  if (doprocessMC) {
    histos.add("Events/resoX", ";X_{Rec} - X_{Gen} [cm]", kTH2D, {axisVertexPosReso, axisVertexNumContrib});
    histos.add("Events/resoY", ";Y_{Rec} - Y_{Gen} [cm]", kTH2D, {axisVertexPosReso, axisVertexNumContrib});
    histos.add("Events/resoZ", ";Z_{Rec} - Z_{Gen} [cm]", kTH2D, {axisVertexPosReso, axisVertexNumContrib});
  }

  auto trackRecoEffHist = histos.add<TH1>("Tracks/recoEff", "", kTH1D, {{2, 0.5, 2.5}});
  trackRecoEffHist->GetXaxis()->SetBinLabel(1, "all");
  trackRecoEffHist->GetXaxis()->SetBinLabel(2, "selected");
  trackRecoEffHist->SetBit(TH1::kIsNotW);

  // kine histograms
  histos.add("Tracks/Kine/pt", "#it{p}_{T}", kTH1D, {axisPt});
  histos.add("Tracks/Kine/eta", "#eta", kTH1D, {axisEta});
  histos.add("Tracks/Kine/phi", "#varphi", kTH1D, {axisPhi});
  histos.add("Tracks/Kine/etavspt", "#eta vs #it{p}_{T}", kTH2F, {axisPt, axisEta});
  histos.add("Tracks/Kine/phivspt", "#varphi vs #it{p}_{T}", kTH2F, {axisPt, axisPhi});
  if (doprocessMC) {
    histos.add("Tracks/Kine/resoPt", "", kTH2D, {axisDeltaPt, axisPt});
    histos.add<TH2>("Tracks/Kine/resoEta", "", kTH2D, {axisDeltaEta, axisEta})->GetYaxis()->SetTitle("#eta_{rec}");
    histos.add<TH2>("Tracks/Kine/resoPhi", "", kTH2D, {axisDeltaPhi, axisPhi})->GetYaxis()->SetTitle("#varphi_{rec}");
  }
  histos.add("Tracks/Kine/relativeResoPt", "relative #it{p}_{T} resolution;#sigma{#it{p}}/#it{p}_{T};#it{p}_{T}", kTH2D, {{axisPt, {100, 0., 0.3}}});
  histos.add("Tracks/Kine/relativeResoPtMean", "mean relative #it{p}_{T} resolution;#LT#sigma{#it{p}}/#it{p}_{T}#GT;#it{p}_{T}", kTProfile, {{axisPt}});

  // track histograms
  auto hselAxis = histos.add<TH1>("Tracks/selection", "trackSelection", kTH1F, {{40, 0.5, 40.5}})->GetXaxis();
  hselAxis->SetBinLabel(1, "Tracks read");
  hselAxis->SetBinLabel(2, "Tracks selected");
  hselAxis->SetBinLabel(3, "passedTrackType");
  hselAxis->SetBinLabel(4, "passedPtRange");
  hselAxis->SetBinLabel(5, "passedEtaRange");
  hselAxis->SetBinLabel(6, "passedTPCNCls");
  hselAxis->SetBinLabel(7, "passedTPCCrossedRows");
  hselAxis->SetBinLabel(8, "passedTPCCrossedRowsOverNCls");
  hselAxis->SetBinLabel(9, "passedTPCChi2NDF");
  hselAxis->SetBinLabel(10, "passedTPCRefit");
  hselAxis->SetBinLabel(11, "passedITSNCls");
  hselAxis->SetBinLabel(12, "passedITSChi2NDF");
  hselAxis->SetBinLabel(13, "passedITSRefit");
  hselAxis->SetBinLabel(14, "passedITSHits");
  hselAxis->SetBinLabel(15, "passedGoldenChi2");
  hselAxis->SetBinLabel(16, "passedDCAxy");
  hselAxis->SetBinLabel(17, "passedDCAz");
  hselAxis->SetBinLabel(18, "isGlobalTrack");
  // Now we combine cuts
  hselAxis->SetBinLabel(19, "Summed cuts#rightarrow");
  hselAxis->SetBinLabel(20, "passedTrackType");
  hselAxis->SetBinLabel(21, "passedPtRange");
  hselAxis->SetBinLabel(22, "passedEtaRange");
  hselAxis->SetBinLabel(23, "passedTPCNCls");
  hselAxis->SetBinLabel(24, "passedTPCCrossedRows");
  hselAxis->SetBinLabel(25, "passedTPCCrossedRowsOverNCls");
  hselAxis->SetBinLabel(26, "passedTPCChi2NDF");
  hselAxis->SetBinLabel(27, "passedTPCRefit");
  hselAxis->SetBinLabel(28, "passedITSNCls");
  hselAxis->SetBinLabel(29, "passedITSChi2NDF");
  hselAxis->SetBinLabel(30, "passedITSRefit");
  hselAxis->SetBinLabel(31, "passedITSHits");
  hselAxis->SetBinLabel(32, "passedGoldenChi2");
  hselAxis->SetBinLabel(33, "passedDCAxy");
  hselAxis->SetBinLabel(34, "passedDCAz");
  hselAxis->SetBinLabel(35, "isGlobalTrack");

  histos.add("Tracks/x", "track #it{x} position at dca in local coordinate system", kTH1D, {axisParX});
  histos.add("Tracks/y", "track #it{y} position at dca in local coordinate system", kTH1D, {axisParY});
  histos.add("Tracks/z", "track #it{z} position at dca in local coordinate system", kTH1D, {axisParZ});
  histos.add("Tracks/alpha", "rotation angle of local wrt. global coordinate system", kTH1D, {axisParAlpha});
  histos.add("Tracks/signed1Pt", "track signed 1/#it{p}_{T}", kTH1D, {axisParSigned1Pt});
  histos.add("Tracks/snp", "sinus of track momentum azimuthal angle", kTH1D, {axisParSnp});
  histos.add("Tracks/tgl", "tangent of the track momentum dip angle", kTH1D, {axisParTgl});
  histos.add("Tracks/flags", "track flag;flag bit", kTH1D, {{64, -0.5, 63.5}});
  histos.add("Tracks/dcaXY", "distance of closest approach in #it{xy} plane;#it{dcaXY} [cm];", kTH1D, {{200, -0.15, 0.15}});
  histos.add("Tracks/dcaZ", "distance of closest approach in #it{z};#it{dcaZ} [cm];", kTH1D, {{200, -0.15, 0.15}});

  histos.add("Tracks/dcaXYvsPt", "distance of closest approach in #it{xy} plane;#it{dcaXY} [cm];", kTH2D, {{200, -0.15, 0.15}, axisPt});
  histos.add("Tracks/dcaZvsPt", "distance of closest approach in #it{z};#it{dcaZ} [cm];", kTH2D, {{200, -0.15, 0.15}, axisPt});

  histos.add("Tracks/length", "track length in cm;#it{Length} [cm];", kTH1D, {{400, 0, 1000}});

  // its histograms
  histos.add("Tracks/ITS/itsNCls", "number of found ITS clusters;# clusters ITS", kTH1D, {{8, -0.5, 7.5}});
  histos.add("Tracks/ITS/itsChi2NCl", "chi2 per ITS cluster;chi2 / cluster ITS", kTH1D, {{100, 0, 40}});
  histos.add("Tracks/ITS/itsHits", "No. of hits vs ITS layer;layer ITS", kTH2D, {{8, -1.5, 6.5}, {8, -0.5, 7.5, "No. of hits"}});
  histos.add("Tracks/ITS/itsHitsUnfiltered", "No. of hits vs ITS layer (unfiltered tracks);layer ITS", kTH2D, {{8, -1.5, 6.5}, {8, -0.5, 7.5, "No. of hits"}});
  histos.add("Tracks/ITS/hasITS", "pt distribution of tracks crossing ITS", kTH1D, {axisPt});
  histos.add("Tracks/ITS/hasITSANDhasTPC", "pt distribution of tracks crossing both ITS and TPC", kTH1D, {axisPt});

  // tpc histograms
  histos.add("Tracks/TPC/tpcNClsFindable", "number of findable TPC clusters;# findable clusters TPC", kTH1D, {{165, -0.5, 164.5}});
  histos.add("Tracks/TPC/tpcNClsFound", "number of found TPC clusters;# clusters TPC", kTH1D, {{165, -0.5, 164.5}});
  histos.add("Tracks/TPC/tpcNClsShared", "number of shared TPC clusters;# shared clusters TPC", kTH1D, {{165, -0.5, 164.5}});
  histos.add("Tracks/TPC/tpcCrossedRows", "number of crossed TPC rows;# crossed rows TPC", kTH1D, {{165, -0.5, 164.5}});
  histos.add("Tracks/TPC/tpcFractionSharedCls", "fraction of shared TPC clusters;fraction shared clusters TPC", kTH1D, {{100, 0., 1.}});
  histos.add("Tracks/TPC/tpcCrossedRowsOverFindableCls", "crossed TPC rows over findable clusters;crossed rows / findable clusters TPC", kTH1D, {{60, 0.7, 1.3}});
  histos.add("Tracks/TPC/tpcChi2NCl", "chi2 per cluster in TPC;chi2 / cluster TPC", kTH1D, {{100, 0, 10}});
  histos.add("Tracks/TPC/hasTPC", "pt distribution of tracks crossing TPC", kTH1D, {axisPt});

  // tracks vs tracks @ IU
  if (doprocessDataIU) {
    // Full distributions
    auto h1 = histos.add<TH1>("Tracks/IU/Pt", "IU: Pt", kTH1F, {axisPt});
    h1->GetXaxis()->SetTitle(Form("%s IU", h1->GetXaxis()->GetTitle()));
    h1 = histos.add<TH1>("Tracks/IU/Eta", "IU: Eta", kTH1F, {axisEta});
    h1->GetXaxis()->SetTitle(Form("%s IU", h1->GetXaxis()->GetTitle()));
    h1 = histos.add<TH1>("Tracks/IU/Phi", "IU: Phi", kTH1F, {axisPhi});
    h1->GetXaxis()->SetTitle(Form("%s IU", h1->GetXaxis()->GetTitle()));

    h1 = histos.add<TH1>("Tracks/IU/x", "IU: x", kTH1F, {axisParX});
    h1->GetXaxis()->SetTitle(Form("%s IU", h1->GetXaxis()->GetTitle()));
    h1 = histos.add<TH1>("Tracks/IU/y", "IU: y", kTH1F, {axisParY});
    h1->GetXaxis()->SetTitle(Form("%s IU", h1->GetXaxis()->GetTitle()));
    h1 = histos.add<TH1>("Tracks/IU/z", "IU: z", kTH1F, {axisParZ});
    h1->GetXaxis()->SetTitle(Form("%s IU", h1->GetXaxis()->GetTitle()));
    h1 = histos.add<TH1>("Tracks/IU/alpha", "rotation angle of local wrt. global coordinate system", kTH1F, {axisParAlpha});
    h1->GetXaxis()->SetTitle(Form("%s IU", h1->GetXaxis()->GetTitle()));
    h1 = histos.add<TH1>("Tracks/IU/signed1Pt", "track signed 1/#it{p}_{T}", kTH1F, {axisParSigned1Pt});
    h1->GetXaxis()->SetTitle(Form("%s IU", h1->GetXaxis()->GetTitle()));
    h1 = histos.add<TH1>("Tracks/IU/snp", "sinus of track momentum azimuthal angle", kTH1F, {axisParSnp});
    h1->GetXaxis()->SetTitle(Form("%s IU", h1->GetXaxis()->GetTitle()));
    h1 = histos.add<TH1>("Tracks/IU/tgl", "tangent of the track momentum dip angle", kTH1F, {axisParTgl});
    h1->GetXaxis()->SetTitle(Form("%s IU", h1->GetXaxis()->GetTitle()));

    // Deltas
    histos.add("Tracks/IU/deltaDCA/Pt", "IU - DCA: Pt", kTH2F, {axisPt, {30, -0.15, 0.15, "#it{p}_{T}^{IU} - #it{p}_{T}^{DCA} [GeV/#it{c}]"}});
    histos.add("Tracks/IU/deltaDCA/Eta", "IU - DCA: Eta", kTH2F, {axisEta, {30, -0.15, 0.15, "#it{#eta}^{IU} - #it{#eta}^{DCA}"}});
    histos.add("Tracks/IU/deltaDCA/Phi", "IU - DCA: Phi", kTH2F, {axisPhi, {30, -0.15, 0.15, "#varphi^{IU} - #varphi^{DCA} [rad]"}});
    // Correlations
    auto h2 = histos.add<TH2>("Tracks/IU/vsDCA/Pt", "IU vs DCA: Pt", kTH2F, {axisPt, axisPt});
    h2->GetXaxis()->SetTitle(Form("%s DCA", h2->GetXaxis()->GetTitle()));
    h2->GetYaxis()->SetTitle(Form("%s IU", h2->GetYaxis()->GetTitle()));
    h2 = histos.add<TH2>("Tracks/IU/vsDCA/Eta", "IU vs DCA: Eta", kTH2F, {axisEta, axisEta});
    h2->GetXaxis()->SetTitle(Form("%s DCA", h2->GetXaxis()->GetTitle()));
    h2->GetYaxis()->SetTitle(Form("%s IU", h2->GetYaxis()->GetTitle()));
    h2 = histos.add<TH2>("Tracks/IU/vsDCA/Phi", "IU vs DCA: Phi", kTH2F, {axisPhi, axisPhi});
    h2->GetXaxis()->SetTitle(Form("%s DCA", h2->GetXaxis()->GetTitle()));
    h2->GetYaxis()->SetTitle(Form("%s IU", h2->GetYaxis()->GetTitle()));
  }

  // filtered tracks @ IU
  if (doprocessDataIUFiltered) {
    // Full distributions
    auto h1 = histos.add<TH1>("Tracks/IUFiltered/Pt", "IU: Pt", kTH1F, {axisPt});
    h1->GetXaxis()->SetTitle(Form("%s IU filtered", h1->GetXaxis()->GetTitle()));
    h1 = histos.add<TH1>("Tracks/IUFiltered/Eta", "IU: Eta", kTH1F, {axisEta});
    h1->GetXaxis()->SetTitle(Form("%s IU filtered", h1->GetXaxis()->GetTitle()));
    h1 = histos.add<TH1>("Tracks/IUFiltered/Phi", "IU: Phi", kTH1F, {axisPhi});
    h1->GetXaxis()->SetTitle(Form("%s IU filtered", h1->GetXaxis()->GetTitle()));

    h1 = histos.add<TH1>("Tracks/IUFiltered/x", "IU: x", kTH1F, {axisParX});
    h1->GetXaxis()->SetTitle(Form("%s IU filtered", h1->GetXaxis()->GetTitle()));
    h1 = histos.add<TH1>("Tracks/IUFiltered/y", "IU: y", kTH1F, {axisParY});
    h1->GetXaxis()->SetTitle(Form("%s IU filtered", h1->GetXaxis()->GetTitle()));
    h1 = histos.add<TH1>("Tracks/IUFiltered/z", "IU: z", kTH1F, {axisParZ});
    h1->GetXaxis()->SetTitle(Form("%s IU filtered", h1->GetXaxis()->GetTitle()));
    h1 = histos.add<TH1>("Tracks/IUFiltered/alpha", "rotation angle of local wrt. global coordinate system", kTH1F, {axisParAlpha});
    h1->GetXaxis()->SetTitle(Form("%s IU filtered", h1->GetXaxis()->GetTitle()));
    h1 = histos.add<TH1>("Tracks/IUFiltered/signed1Pt", "track signed 1/#it{p}_{T}", kTH1F, {axisParSigned1Pt});
    h1->GetXaxis()->SetTitle(Form("%s IU filtered", h1->GetXaxis()->GetTitle()));
    h1 = histos.add<TH1>("Tracks/IUFiltered/snp", "sinus of track momentum azimuthal angle", kTH1F, {axisParSnp});
    h1->GetXaxis()->SetTitle(Form("%s IU filtered", h1->GetXaxis()->GetTitle()));
    h1 = histos.add<TH1>("Tracks/IUFiltered/tgl", "tangent of the track momentum dip angle", kTH1F, {axisParTgl});
    h1->GetXaxis()->SetTitle(Form("%s IU filtered", h1->GetXaxis()->GetTitle()));
  }

  // tracks IU filtered
}

//**************************************************************************************************
/**
 * Fill reco level histograms.
 */
//**************************************************************************************************
template <bool IS_MC, typename C, typename T>
void qaEventTrack::fillRecoHistograms(const C& collision, const T& tracks, const aod::FullTracks& tracksUnfiltered)
{
  // fill reco collision related histograms
  if (!isSelectedCollision<true>(collision)) {
    return;
  }

  int nTracks = 0;
  for (const auto& track : tracks) {
    histos.fill(HIST("Tracks/selection"), 1.f);
    if (!isSelectedTrack<IS_MC>(track)) {
      continue;
    }
    histos.fill(HIST("Tracks/selection"), 2.f);
    ++nTracks;
    if (track.passedTrackType()) {
      histos.fill(HIST("Tracks/selection"), 3.f);
    }
    if (track.passedPtRange()) {
      histos.fill(HIST("Tracks/selection"), 4.f);
    }
    if (track.passedEtaRange()) {
      histos.fill(HIST("Tracks/selection"), 5.f);
    }
    if (track.passedTPCNCls()) {
      histos.fill(HIST("Tracks/selection"), 6.f);
    }
    if (track.passedTPCCrossedRows()) {
      histos.fill(HIST("Tracks/selection"), 7.f);
    }
    if (track.passedTPCCrossedRowsOverNCls()) {
      histos.fill(HIST("Tracks/selection"), 8.f);
    }
    if (track.passedTPCChi2NDF()) {
      histos.fill(HIST("Tracks/selection"), 9.f);
    }
    if (track.passedTPCRefit()) {
      histos.fill(HIST("Tracks/selection"), 10.f);
    }
    if (track.passedITSNCls()) {
      histos.fill(HIST("Tracks/selection"), 11.f);
    }
    if (track.passedITSChi2NDF()) {
      histos.fill(HIST("Tracks/selection"), 12.f);
    }
    if (track.passedITSRefit()) {
      histos.fill(HIST("Tracks/selection"), 13.f);
    }
    if (track.passedITSHits()) {
      histos.fill(HIST("Tracks/selection"), 14.f);
    }
    if (track.passedGoldenChi2()) {
      histos.fill(HIST("Tracks/selection"), 15.f);
    }
    if (track.passedDCAxy()) {
      histos.fill(HIST("Tracks/selection"), 16.f);
    }
    if (track.passedDCAz()) {
      histos.fill(HIST("Tracks/selection"), 17.f);
    }
    if (track.isGlobalTrack()) {
      histos.fill(HIST("Tracks/selection"), 18.f);
    }
    // Filling combined cuts
    if (track.passedTrackType()) {
      histos.fill(HIST("Tracks/selection"), 20.f);
    } else {
      continue;
    }
    if (track.passedPtRange()) {
      histos.fill(HIST("Tracks/selection"), 21.f);
    } else {
      continue;
    }
    if (track.passedEtaRange()) {
      histos.fill(HIST("Tracks/selection"), 22.f);
    } else {
      continue;
    }
    if (track.passedTPCNCls()) {
      histos.fill(HIST("Tracks/selection"), 23.f);
    } else {
      continue;
    }
    if (track.passedTPCCrossedRows()) {
      histos.fill(HIST("Tracks/selection"), 24.f);
    } else {
      continue;
    }
    if (track.passedTPCCrossedRowsOverNCls()) {
      histos.fill(HIST("Tracks/selection"), 25.f);
    } else {
      continue;
    }
    if (track.passedTPCChi2NDF()) {
      histos.fill(HIST("Tracks/selection"), 26.f);
    } else {
      continue;
    }
    if (track.passedTPCRefit()) {
      histos.fill(HIST("Tracks/selection"), 27.f);
    } else {
      continue;
    }
    if (track.passedITSNCls()) {
      histos.fill(HIST("Tracks/selection"), 28.f);
    } else {
      continue;
    }
    if (track.passedITSChi2NDF()) {
      histos.fill(HIST("Tracks/selection"), 29.f);
    } else {
      continue;
    }
    if (track.passedITSRefit()) {
      histos.fill(HIST("Tracks/selection"), 30.f);
    } else {
      continue;
    }
    if (track.passedITSHits()) {
      histos.fill(HIST("Tracks/selection"), 31.f);
    } else {
      continue;
    }
    if (track.passedGoldenChi2()) {
      histos.fill(HIST("Tracks/selection"), 32.f);
    } else {
      continue;
    }
    if (track.passedDCAxy()) {
      histos.fill(HIST("Tracks/selection"), 33.f);
    } else {
      continue;
    }
    if (track.passedDCAz()) {
      histos.fill(HIST("Tracks/selection"), 34.f);
    } else {
      continue;
    }
    if (track.isGlobalTrack()) {
      histos.fill(HIST("Tracks/selection"), 35.f);
    }
  }

  histos.fill(HIST("Events/posX"), collision.posX());
  histos.fill(HIST("Events/posY"), collision.posY());
  histos.fill(HIST("Events/posZ"), collision.posZ());
  histos.fill(HIST("Events/posXY"), collision.posX(), collision.posY());

  histos.fill(HIST("Events/posXvsNContrib"), collision.posX(), collision.numContrib());
  histos.fill(HIST("Events/posYvsNContrib"), collision.posY(), collision.numContrib());
  histos.fill(HIST("Events/posZvsNContrib"), collision.posZ(), collision.numContrib());

  histos.fill(HIST("Events/nContrib"), collision.numContrib());
  histos.fill(HIST("Events/nContribVsMult"), collision.numContrib(), nTracks);
  histos.fill(HIST("Events/vertexChi2"), collision.chi2());

  histos.fill(HIST("Events/covXX"), collision.covXX());
  histos.fill(HIST("Events/covXY"), collision.covXY());
  histos.fill(HIST("Events/covXZ"), collision.covXZ());
  histos.fill(HIST("Events/covYY"), collision.covYY());
  histos.fill(HIST("Events/covYZ"), collision.covYZ());
  histos.fill(HIST("Events/covZZ"), collision.covZZ());

  histos.fill(HIST("Events/nTracks"), nTracks);

  // vertex resolution
  if constexpr (IS_MC) {
    if (collision.has_mcCollision()) {
      const auto mcColl = collision.mcCollision();
      histos.fill(HIST("Events/resoX"), collision.posX() - mcColl.posX(), collision.numContrib());
      histos.fill(HIST("Events/resoY"), collision.posY() - mcColl.posY(), collision.numContrib());
      histos.fill(HIST("Events/resoZ"), collision.posZ() - mcColl.posZ(), collision.numContrib());
    }
  }

  histos.fill(HIST("Tracks/recoEff"), 1, tracks.tableSize());
  histos.fill(HIST("Tracks/recoEff"), 2, tracks.size());

  // unfiltered track related histograms
  for (const auto& trackUnfiltered : tracksUnfiltered) {
    // fill ITS variables
    int itsNhits = 0;
    for (unsigned int i = 0; i < 7; i++) {
      if (trackUnfiltered.itsClusterMap() & (1 << i)) {
        itsNhits += 1;
      }
    }
    bool trkHasITS = false;
    for (unsigned int i = 0; i < 7; i++) {
      if (trackUnfiltered.itsClusterMap() & (1 << i)) {
        trkHasITS = true;
        histos.fill(HIST("Tracks/ITS/itsHitsUnfiltered"), i, itsNhits);
      }
    }
    if (!trkHasITS) {
      histos.fill(HIST("Tracks/ITS/itsHitsUnfiltered"), -1, itsNhits);
    }
  }

  // track related histograms
  for (const auto& track : tracks) {
    if (!isSelectedTrack<IS_MC>(track)) {
      continue;
    }
    // fill kinematic variables
    histos.fill(HIST("Tracks/Kine/pt"), track.pt());
    histos.fill(HIST("Tracks/Kine/eta"), track.eta());
    histos.fill(HIST("Tracks/Kine/phi"), track.phi());
    histos.fill(HIST("Tracks/Kine/etavspt"), track.pt(), track.eta());
    histos.fill(HIST("Tracks/Kine/phivspt"), track.pt(), track.phi());
    histos.fill(HIST("Tracks/Kine/relativeResoPt"), track.pt(), track.pt() * std::sqrt(track.c1Pt21Pt2()));
    histos.fill(HIST("Tracks/Kine/relativeResoPtMean"), track.pt(), track.pt() * std::sqrt(track.c1Pt21Pt2()));

    // fill track parameters
    histos.fill(HIST("Tracks/alpha"), track.alpha());
    histos.fill(HIST("Tracks/x"), track.x());
    histos.fill(HIST("Tracks/y"), track.y());
    histos.fill(HIST("Tracks/z"), track.z());
    histos.fill(HIST("Tracks/signed1Pt"), track.signed1Pt());
    histos.fill(HIST("Tracks/snp"), track.snp());
    histos.fill(HIST("Tracks/tgl"), track.tgl());
    for (unsigned int i = 0; i < 64; i++) {
      if (track.flags() & (1 << i)) {
        histos.fill(HIST("Tracks/flags"), i);
      }
    }
    histos.fill(HIST("Tracks/dcaXY"), track.dcaXY());
    histos.fill(HIST("Tracks/dcaZ"), track.dcaZ());
    histos.fill(HIST("Tracks/dcaXYvsPt"), track.dcaXY(), track.pt());
    histos.fill(HIST("Tracks/dcaZvsPt"), track.dcaZ(), track.pt());
    histos.fill(HIST("Tracks/length"), track.length());

    // fill ITS variables
    histos.fill(HIST("Tracks/ITS/itsNCls"), track.itsNCls());
    histos.fill(HIST("Tracks/ITS/itsChi2NCl"), track.itsChi2NCl());
    int itsNhits = 0;
    for (unsigned int i = 0; i < 7; i++) {
      if (track.itsClusterMap() & (1 << i)) {
        itsNhits += 1;
      }
    }
    bool trkHasITS = false;
    for (unsigned int i = 0; i < 7; i++) {
      if (track.itsClusterMap() & (1 << i)) {
        trkHasITS = true;
        histos.fill(HIST("Tracks/ITS/itsHits"), i, itsNhits);
      }
    }
    if (!trkHasITS) {
      histos.fill(HIST("Tracks/ITS/itsHits"), -1, itsNhits);
    }

    // fill TPC variables
    histos.fill(HIST("Tracks/TPC/tpcNClsFindable"), track.tpcNClsFindable());
    histos.fill(HIST("Tracks/TPC/tpcNClsFound"), track.tpcNClsFound());
    histos.fill(HIST("Tracks/TPC/tpcNClsShared"), track.tpcNClsShared());
    histos.fill(HIST("Tracks/TPC/tpcCrossedRows"), track.tpcNClsCrossedRows());
    histos.fill(HIST("Tracks/TPC/tpcCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
    histos.fill(HIST("Tracks/TPC/tpcFractionSharedCls"), track.tpcFractionSharedCls());
    histos.fill(HIST("Tracks/TPC/tpcChi2NCl"), track.tpcChi2NCl());

    if constexpr (IS_MC) {
      if (track.has_mcParticle()) {
        // resolution plots
        auto particle = track.mcParticle();
        histos.fill(HIST("Tracks/Kine/resoPt"), track.pt() - particle.pt(), track.pt());
        histos.fill(HIST("Tracks/Kine/resoEta"), track.eta() - particle.eta(), track.eta());
        histos.fill(HIST("Tracks/Kine/resoPhi"), track.phi() - particle.phi(), track.phi());
      }
    }

    // ITS-TPC matching pt-distributions
    if (track.hasITS()) {
      histos.fill(HIST("Tracks/ITS/hasITS"), track.pt());
    }
    if (track.hasTPC()) {
      histos.fill(HIST("Tracks/TPC/hasTPC"), track.pt());
    }
    if (track.hasITS() && track.hasTPC()) {
      histos.fill(HIST("Tracks/ITS/hasITSANDhasTPC"), track.pt());
    }
  }
}
