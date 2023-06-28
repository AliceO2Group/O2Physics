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
/// \file   qaKFEventTrack.cxx
/// \author Annalena Kalteyer <annalena.sophie.kalteyer@cern.ch>, GSI Darmstadt
/// \brief  Task to test the performance of the KFParticle package
///

#include "Tools/KFparticle/qaKFEventTrack.h"
#include <CCDB/BasicCCDBManager.h>
#include <string>
#include "TableHelper.h"
#include <iostream>
using namespace std;

/// includes O2
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"

/// includes O2Physics
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Tools/KFparticle/KFUtilities.h"

/// includes KFParticle
#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticleBase.h"
#include "KFVertex.h"

#ifndef HomogeneousField

#define HomogeneousField

#endif

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::dataformats;

struct qaKFEventTrack {

  /// general steering settings
  Configurable<bool> isRun3{"isRun3", true, "Is Run3 dataset"};
  Configurable<std::string> ccdbUrl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  int runNumber;
  double magneticField = 0.;

  /// Histogram Configurables
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 24., 36., 50.0}, ""};

  /// option to select good events
  Configurable<bool> eventSelection{"eventSelection", true, "select good events"}; // currently only sel8 is defined for run3
  /// options to select only specific tracks
  Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  /// singe track selections
  Configurable<float> d_pTMin{"d_pTMin", 0.3, "minimum momentum for tracks"};
  Configurable<float> d_etaRange{"d_etaRange", 0.8, "eta Range for tracks"};
  /// Option to write variables in a tree
  Configurable<double> d_DwnSmplFact{"d_DwnSmplFact", 1., "Downsampling factor for tree"};
  Configurable<bool> writeHistograms{"writeHistograms", true, "write histograms"};
  Configurable<bool> writeTree{"writeTree", false, "write daughter variables in a tree"};

  // Define which track selection should be used:
  // 0 -> No track selection is applied
  // 1 kGlobalTrack = kQualityTracks | kPrimaryTracks | kInAcceptanceTracks
  //        kQualityTracks = kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF |
  //                         kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits
  //        kPrimaryTracks = kGoldenChi2 | kDCAxy | kDCAz
  //        kInAcceptanceTracks = kPtRange | kEtaRange
  // 2 kGlobalTrackWoPtEta = kQualityTracks | kPrimaryTracks
  // 3 kGlobalTrackWoDCA = kQualityTracks | kInAcceptanceTracks
  // 4 kQualityTracks
  // 5 kInAcceptanceTracks
  Filter trackFilter = (trackSelection.node() == 0) ||
                       ((trackSelection.node() == 1) && requireGlobalTrackInFilter()) ||
                       ((trackSelection.node() == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
                       ((trackSelection.node() == 3) && requireGlobalTrackWoDCAInFilter()) ||
                       ((trackSelection.node() == 4) && requireQualityTracksInFilter()) ||
                       ((trackSelection.node() == 5) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks));

  Filter eventFilter = (o2::aod::evsel::sel8 == true);

  using CollisionTableData = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
  using TrackTableData = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::TrackSelection>;

  HistogramRegistry histos;
  /// Table to be produced
  Produces<o2::aod::TreeTracks> rowKFTracks;

  int source = 0;
  int pVContrib = 0;

  void initMagneticFieldCCDB(o2::aod::BCsWithTimestamps::iterator const& bc, int& mRunNumber,
                             o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdb, std::string ccdbPathGrp, o2::base::MatLayerCylSet* lut,
                             bool isRun3)
  {

    if (mRunNumber != bc.runNumber()) {

      LOGF(info, "====== initCCDB function called (isRun3==%d)", isRun3);
      if (!isRun3) { // Run 2 GRP object
        o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(ccdbPathGrp, bc.timestamp());
        if (grpo == nullptr) {
          LOGF(fatal, "Run 2 GRP object (type o2::parameters::GRPObject) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
        }
        o2::base::Propagator::initFieldFromGRP(grpo);
        o2::base::Propagator::Instance()->setMatLUT(lut);
        LOGF(info, "Setting magnetic field to %d kG for run %d from its GRP CCDB object (type o2::parameters::GRPObject)", grpo->getNominalL3Field(), bc.runNumber());
      } else { // Run 3 GRP object
        o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbPathGrp, bc.timestamp());
        if (grpo == nullptr) {
          LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
        }
        o2::base::Propagator::initFieldFromGRP(grpo);
        o2::base::Propagator::Instance()->setMatLUT(lut);
        LOGF(info, "Setting magnetic field to current %f A for run %d from its GRP CCDB object (type o2::parameters::GRPMagField)", grpo->getL3Current(), bc.runNumber());
      }
      mRunNumber = bc.runNumber();
    }
  } /// end initMagneticFieldCCDB

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;

    const AxisSpec axisVertexPosX{100, -0.05, 0.05, "X [cm]"};
    const AxisSpec axisVertexPosY{100, -0.05, 0.05, "Y [cm]"};
    const AxisSpec axisVertexPosZ{100, -20., 20., "Z [cm]"};
    const AxisSpec axisVertexNumContrib{160, 0, 160, "Number Of contributors to the PV"};
    const AxisSpec axisVertexCov{100, -0.00005, 0.00005};

    const AxisSpec axisParX{100, -0.1, 0.1, "#it{x} [cm]"};
    const AxisSpec axisParY{100, -0.1, 0.1, "#it{y} [cm]"};
    const AxisSpec axisParZ{100, -20., 20., "#it{z} [cm]"};
    const AxisSpec axisParPX{binsPt, "#it{p}_{x} [GeV/c]"};
    const AxisSpec axisParPY{binsPt, "#it{p}_{y} [GeV/c]"};
    const AxisSpec axisParPZ{binsPt, "#it{p}_{z} [GeV/c]"};

    if (writeHistograms) {
      /// collisions
      histos.add("Events/covXX", ";Cov_{xx} [cm^{2}]", kTH1D, {axisVertexCov});
      histos.add("Events/covXY", ";Cov_{xy} [cm^{2}]", kTH1D, {axisVertexCov});
      histos.add("Events/covYY", ";Cov_{yy} [cm^{2}]", kTH1D, {axisVertexCov});
      histos.add("Events/covXZ", ";Cov_{xz} [cm^{2}]", kTH1D, {axisVertexCov});
      histos.add("Events/covYZ", ";Cov_{yz} [cm^{2}]", kTH1D, {axisVertexCov});
      histos.add("Events/covZZ", ";Cov_{zz} [cm^{2}]", kTH1D, {axisVertexCov});

      histos.add("EventsKF/posX", "", kTH1D, {axisVertexPosX});
      histos.add("EventsKF/posY", "", kTH1D, {axisVertexPosY});
      histos.add("EventsKF/posZ", "", kTH1D, {axisVertexPosZ});
      histos.add("EventsKF/posXY", "", kTH2D, {axisVertexPosX, axisVertexPosY});
      histos.add("EventsKF/nContrib", "", kTH1D, {axisVertexNumContrib});
      histos.add("EventsKF/vertexChi2", ";#chi^{2}", kTH1D, {{100, 0, 100}});
      histos.add("EventsKF/covXX", ";Cov_{xx} [cm^{2}]", kTH1D, {axisVertexCov});
      histos.add("EventsKF/covXY", ";Cov_{xy} [cm^{2}]", kTH1D, {axisVertexCov});
      histos.add("EventsKF/covYY", ";Cov_{yy} [cm^{2}]", kTH1D, {axisVertexCov});
      histos.add("EventsKF/covXZ", ";Cov_{xz} [cm^{2}]", kTH1D, {axisVertexCov});
      histos.add("EventsKF/covYZ", ";Cov_{yz} [cm^{2}]", kTH1D, {axisVertexCov});
      histos.add("EventsKF/covZZ", ";Cov_{zz} [cm^{2}]", kTH1D, {axisVertexCov});

      /// tracks
      histos.add("Tracks/x", "track #it{x} position at dca in local coordinate system", kTH1D, {axisParX});
      histos.add("Tracks/y", "track #it{y} position at dca in local coordinate system", kTH1D, {axisParY});
      histos.add("Tracks/z", "track #it{z} position at dca in local coordinate system", kTH1D, {axisParZ});
      histos.add("Tracks/px", "track #it{p_{x}} momentum at dca in local coordinate system", kTH1D, {axisParPX});
      histos.add("Tracks/py", "track #it{p_{y}} momentum at dca in local coordinate system", kTH1D, {axisParPY});
      histos.add("Tracks/pz", "track #it{p_{z}} momentum at dca in local coordinate system", kTH1D, {axisParPZ});
      histos.add("Tracks/chi2perNDF", "Chi2/NDF of the track;#it{chi2/ndf};", kTH1D, {{100, 0., 25.}});
      histos.add("Tracks/dcaXYToPV", "distance of closest approach in #it{xy} plane;#it{dcaXY} [cm];", kTH1D, {{100, -0.15, 0.15}});
      histos.add("Tracks/dcaToPV", "distance of closest approach ;#it{dca} [cm];", kTH1D, {{100, -0.15, 0.15}});
      histos.add("Tracks/dcaToPVLargeRange", "distance of closest approach ;#it{dca} [cm];", kTH1D, {{200, 0., 20.}});
      histos.add("Tracks/dcaToPVLargeRangeMCBeforeReassignment", "distance of closest approach ;#it{dca} [cm];", kTH1D, {{100, 0., 20.}});
      histos.add("Tracks/dcaToPVLargeRangeMCAfterReassignment", "distance of closest approach ;#it{dca} [cm];", kTH1D, {{100, 0., 20.}});
      histos.add("Tracks/deviationPiToPV", "deviation of Pi to PV", kTH1D, {{100, 0., 10.}});
    }

    auto hSelectionMC = histos.add<TH1>("DZeroCandTopo/SelectionsMC", "Selections MC", kTH1D, {{5, 0.5, 5.5}});
    hSelectionMC->GetXaxis()->SetBinLabel(hSelectionMC->FindBin(1), "All Tracks");
    hSelectionMC->GetXaxis()->SetBinLabel(hSelectionMC->FindBin(2), "No MC particle");
    hSelectionMC->GetXaxis()->SetBinLabel(hSelectionMC->FindBin(3), "Unmatched");
    hSelectionMC->GetXaxis()->SetBinLabel(hSelectionMC->FindBin(4), "Wrong PV");
    hSelectionMC->GetXaxis()->SetBinLabel(hSelectionMC->FindBin(5), "DCA Z > 10cm");
  } /// End init

  /// Function for single track selection
  template <typename T>
  bool isSelectedTracks(const T& track1)
  {
    if (track1.p() < d_pTMin) {
      return false;
    }
    /// Eta range
    if (abs(track1.eta()) > d_etaRange) {
      return false;
    }
    return true;
  }

  template <typename T1, typename T2>
  void fillHistograms(const T1& kfpTrackPi, const T2& KFPion, const T2& KFPV)
  {
    if (writeHistograms) {
      /// fill daughter track parameters
      histos.fill(HIST("Tracks/x"), kfpTrackPi.GetX());
      histos.fill(HIST("Tracks/y"), kfpTrackPi.GetY());
      histos.fill(HIST("Tracks/z"), kfpTrackPi.GetZ());
      histos.fill(HIST("Tracks/px"), kfpTrackPi.GetPx());
      histos.fill(HIST("Tracks/py"), kfpTrackPi.GetPy());
      histos.fill(HIST("Tracks/pz"), kfpTrackPi.GetPz());
      histos.fill(HIST("Tracks/chi2perNDF"), kfpTrackPi.GetChi2perNDF());
      histos.fill(HIST("Tracks/dcaXYToPV"), KFPion.GetDistanceFromVertexXY(KFPV));
      histos.fill(HIST("Tracks/dcaToPV"), KFPion.GetDistanceFromVertex(KFPV));
      histos.fill(HIST("Tracks/dcaToPVLargeRange"), KFPion.GetDistanceFromVertex(KFPV));
      histos.fill(HIST("Tracks/deviationPiToPV"), KFPion.GetDeviationFromVertex(KFPV));
    }
  }

  template <typename T1, typename T2, typename T3, typename T4>
  void writeVarTree(const T1& kfpTrackPi, const T2& KFPion, const T2& KFPV, const T3& track, const int source, const int pVContrib, const int runNumber, const T4& collision)
  {
    const double pseudoRndm = track.pt() * 1000. - (int64_t)(track.pt() * 1000);
    if (pseudoRndm < d_DwnSmplFact) {
      if (writeTree) {
        /// Filling the tree
        rowKFTracks(source,
                    pVContrib,
                    track.dcaXY(),
                    track.dcaZ(),
                    KFPion.GetDistanceFromVertex(KFPV),
                    KFPion.GetDistanceFromVertexXY(KFPV),
                    track.sign(),
                    track.p(),
                    track.eta(),
                    track.phi(),
                    track.tpcSignal(),
                    runNumber,
                    collision.posX(),
                    collision.posY(),
                    collision.posZ(),
                    collision.covXX(),
                    collision.covYY(),
                    collision.covZZ(),
                    collision.numContrib(),
                    collision.chi2());
      }
    }
  }

  /// Process function for data
  void processData(soa::Filtered<CollisionTableData>::iterator const& collision, soa::Filtered<TrackTableData> const& tracks, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (runNumber != bc.runNumber()) {
      initMagneticFieldCCDB(bc, runNumber, ccdb, isRun3 ? ccdbPathGrpMag : ccdbPathGrp, lut, isRun3);
      magneticField = o2::base::Propagator::Instance()->getNominalBz();
/// Set magnetic field for KF vertexing
#ifdef HomogeneousField
      KFParticle::SetField(magneticField);
#endif
    }
    if (writeHistograms) {
      histos.fill(HIST("Events/covXX"), collision.covXX());
      histos.fill(HIST("Events/covXY"), collision.covXY());
      histos.fill(HIST("Events/covXZ"), collision.covXZ());
      histos.fill(HIST("Events/covYY"), collision.covYY());
      histos.fill(HIST("Events/covYZ"), collision.covYZ());
      histos.fill(HIST("Events/covZZ"), collision.covZZ());
    }
    /// set KF primary vertex
    KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
    KFParticle KFPV(kfpVertex);
    if (writeHistograms) {
      /// fill collision parameters
      histos.fill(HIST("EventsKF/posX"), kfpVertex.GetX());
      histos.fill(HIST("EventsKF/posY"), kfpVertex.GetY());
      histos.fill(HIST("EventsKF/posZ"), kfpVertex.GetZ());
      histos.fill(HIST("EventsKF/posXY"), kfpVertex.GetX(), kfpVertex.GetY());
      histos.fill(HIST("EventsKF/nContrib"), kfpVertex.GetNContributors());
      histos.fill(HIST("EventsKF/vertexChi2"), kfpVertex.GetChi2());
      histos.fill(HIST("EventsKF/covXX"), kfpVertex.GetCovariance(0));
      histos.fill(HIST("EventsKF/covXY"), kfpVertex.GetCovariance(1));
      histos.fill(HIST("EventsKF/covYY"), kfpVertex.GetCovariance(2));
      histos.fill(HIST("EventsKF/covXZ"), kfpVertex.GetCovariance(3));
      histos.fill(HIST("EventsKF/covYZ"), kfpVertex.GetCovariance(4));
      histos.fill(HIST("EventsKF/covZZ"), kfpVertex.GetCovariance(5));
    }

    int ntracks = 0;

    for (auto& track : tracks) {
      source = 0;
      pVContrib = 0;

      if (track.isPVContributor()) {
        ntracks = ntracks + 1;
      }

      /// Apply single track selection
      if (!isSelectedTracks(track)) {
        continue;
      }

      KFPTrack kfpTrack;
      kfpTrack = createKFPTrackFromTrack(track);
      KFParticle KFParticleTrack(kfpTrack, 211);

      if (track.hasITS()) {
        source |= kITS;
      }
      if (track.hasTPC()) {
        source |= kTPC;
      }
      if (track.hasTRD()) {
        source |= kTRD;
      }
      if (track.hasTOF()) {
        source |= kTOF;
      }
      if (track.isPVContributor()) {
        pVContrib = 1;
      }

      fillHistograms(kfpTrack, KFParticleTrack, KFPV);
      writeVarTree(kfpTrack, KFParticleTrack, KFPV, track, source, pVContrib, runNumber, collision);
    }
  }
  PROCESS_SWITCH(qaKFEventTrack, processData, "process data", true);

  /// Process function for MC
  using CollisionTableMC = soa::Join<CollisionTableData, aod::McCollisionLabels>;
  using CollisionTableDataMult = soa::Join<aod::Collisions, aod::Mults, aod::McCollisionLabels>;
  using TrackTableMC = soa::Join<TrackTableData, aod::McTrackLabels>;
  // Preslice<aod::McCollisionLabels> perMcCollision = aod::mccollisionlabel::mcCollisionId;
  void processMC(soa::Filtered<CollisionTableMC>::iterator const& collision, soa::Filtered<CollisionTableMC> const& collisions, soa::Filtered<TrackTableMC> const& tracks, aod::McParticles const& mcParticles, aod::McCollisions const& mcCollisions, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (runNumber != bc.runNumber()) {
      initMagneticFieldCCDB(bc, runNumber, ccdb, isRun3 ? ccdbPathGrpMag : ccdbPathGrp, lut, isRun3);
      magneticField = o2::base::Propagator::Instance()->getNominalBz();
/// Set magnetic field for KF vertexing
#ifdef HomogeneousField
      KFParticle::SetField(magneticField);
#endif
    }

    /// Remove Collisions without a MC Collision
    if (!collision.has_mcCollision()) {
      return;
    }
    if (writeHistograms) {
      histos.fill(HIST("Events/covXX"), collision.covXX());
      histos.fill(HIST("Events/covXY"), collision.covXY());
      histos.fill(HIST("Events/covXZ"), collision.covXZ());
      histos.fill(HIST("Events/covYY"), collision.covYY());
      histos.fill(HIST("Events/covYZ"), collision.covYZ());
      histos.fill(HIST("Events/covZZ"), collision.covZZ());
    }
    /// set KF primary vertex
    KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
    KFParticle KFPV(kfpVertex);

    KFPVertex kfpVertexDefault = createKFPVertexFromCollision(collision);
    KFParticle KFPVDefault(kfpVertexDefault);
    if (writeHistograms) {
      /// fill collision parameters
      histos.fill(HIST("EventsKF/posX"), kfpVertex.GetX());
      histos.fill(HIST("EventsKF/posY"), kfpVertex.GetY());
      histos.fill(HIST("EventsKF/posZ"), kfpVertex.GetZ());
      histos.fill(HIST("EventsKF/posXY"), kfpVertex.GetX(), kfpVertex.GetY());
      histos.fill(HIST("EventsKF/nContrib"), kfpVertex.GetNContributors());
      histos.fill(HIST("EventsKF/vertexChi2"), kfpVertex.GetChi2());
      histos.fill(HIST("EventsKF/covXX"), kfpVertex.GetCovariance(0));
      histos.fill(HIST("EventsKF/covXY"), kfpVertex.GetCovariance(1));
      histos.fill(HIST("EventsKF/covYY"), kfpVertex.GetCovariance(2));
      histos.fill(HIST("EventsKF/covXZ"), kfpVertex.GetCovariance(3));
      histos.fill(HIST("EventsKF/covYZ"), kfpVertex.GetCovariance(4));
      histos.fill(HIST("EventsKF/covZZ"), kfpVertex.GetCovariance(5));
    }

    int ntracks = 0;

    for (auto& track : tracks) {

      source = 0;
      pVContrib = 0;

      if (track.isPVContributor()) {
        ntracks = ntracks + 1;
      }

      histos.fill(HIST("DZeroCandTopo/SelectionsMC"), 1.f);

      /// Apply single track selection
      if (!isSelectedTracks(track)) {
        continue;
      }

      // /// Check whether the track was assigned to the true MC PV
      // auto particle = track.mcParticle();
      // auto collMC = particle.mcCollision();
      // auto mcCollID_recoColl = track.collision_as<CollisionTableMC>().mcCollisionId();
      // auto mcCollID_particle = particle.mcCollisionId();
      // bool indexMatchOK = (mcCollID_recoColl == mcCollID_particle);
      // if (!indexMatchOK) {
      //   histos.fill(HIST("DZeroCandTopo/SelectionsMC"), 4.f);
      //   const auto matchedCollisions = collisions.sliceBy(perMcCollision, collMC.globalIndex());
      //   int i = 0;
      //   std::array<float, 5> dcaZ{100, 100, 100, 100, 100};
      //   float min = 100;
      //   for (auto matchedCollision : matchedCollisions) {
      //     dcaZ[i] = abs(matchedCollision.posZ() - collMC.posZ());
      //     if (i == 0) {
      //       min = dcaZ[i];
      //     }
      //     if (i > 0) {
      //       if (dcaZ[i] < dcaZ[i - 1]) {
      //         min = dcaZ[i];
      //       }
      //     }

      //     i = i + 1;
      //   }
      //   if (min > 10.) {
      //     histos.fill(HIST("DZeroCandTopo/SelectionsMC"), 5.f);
      //   }
      //   int j = 0;
      //   for (auto matchedCollision : matchedCollisions) {
      //     if (i == 1) {
      //       kfpVertex = createKFPVertexFromCollision(matchedCollision);
      //       KFParticle KFPVNew(kfpVertex);
      //       KFPV = KFPVNew;
      //     }
      //     if (i > 1) {
      //       if (abs(matchedCollision.posZ() - collMC.posZ()) == min) {
      //         kfpVertex = createKFPVertexFromCollision(matchedCollision);
      //         KFParticle KFPVNew(kfpVertex);
      //         KFPV = KFPVNew;
      //       }
      //     }
      //     j = j + 1;
      //   }
      // }

      /// Remove fake tracks
      if (!track.has_mcParticle()) {
        histos.fill(HIST("DZeroCandTopo/SelectionsMC"), 2.f);
        continue;
      }
      /// Remove unmatched tracks
      if (track.collisionId() <= 0) {
        histos.fill(HIST("DZeroCandTopo/SelectionsMC"), 3.f);
        continue;
      }

      KFPTrack kfpTrack;
      kfpTrack = createKFPTrackFromTrack(track);
      KFParticle KFParticleTrack(kfpTrack, 211);

      if (track.hasITS()) {
        source |= kITS;
      }
      if (track.hasTPC()) {
        source |= kTPC;
      }
      if (track.hasTRD()) {
        source |= kTRD;
      }
      if (track.hasTOF()) {
        source |= kTOF;
      }
      if (track.isPVContributor()) {
        pVContrib = 1;
      }

      fillHistograms(kfpTrack, KFParticleTrack, KFPV);
      writeVarTree(kfpTrack, KFParticleTrack, KFPV, track, source, pVContrib, runNumber, collision);
      if (writeHistograms) {
        histos.fill(HIST("Tracks/dcaToPVLargeRangeMCBeforeReassignment"), KFParticleTrack.GetDistanceFromVertex(KFPVDefault));
        histos.fill(HIST("Tracks/dcaToPVLargeRangeMCAfterReassignment"), KFParticleTrack.GetDistanceFromVertex(KFPV));
      }
    }
  }
  PROCESS_SWITCH(qaKFEventTrack, processMC, "process mc", false);
};

struct qaKFEvent {

  /// general steering settings
  Configurable<bool> isRun3{"isRun3", true, "Is Run3 dataset"};
  Configurable<std::string> ccdbUrl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGeo{"ccdbPathGeo", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  int runNumber;
  double magneticField = 0.;
  int tfID = 0;
  /// option to select good events
  Configurable<bool> eventSelection{"eventSelection", true, "select good events"}; // currently only sel8 is defined for run3
  Configurable<bool> writeTree{"writeTree", true, "write daughter variables in a tree"};

  using CollisionTableData = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;

  Produces<o2::aod::TreeCollisions> rowKFCollisions;

  void initMagneticFieldCCDB(o2::aod::BCsWithTimestamps::iterator const& bc, int& mRunNumber,
                             o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdb, std::string ccdbPathGrp, o2::base::MatLayerCylSet* lut,
                             bool isRun3)
  {

    if (mRunNumber != bc.runNumber()) {

      LOGF(info, "====== initCCDB function called (isRun3==%d)", isRun3);
      if (!isRun3) { // Run 2 GRP object
        o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(ccdbPathGrp, bc.timestamp());
        if (grpo == nullptr) {
          LOGF(fatal, "Run 2 GRP object (type o2::parameters::GRPObject) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
        }
        o2::base::Propagator::initFieldFromGRP(grpo);
        o2::base::Propagator::Instance()->setMatLUT(lut);
        LOGF(info, "Setting magnetic field to %d kG for run %d from its GRP CCDB object (type o2::parameters::GRPObject)", grpo->getNominalL3Field(), bc.runNumber());
      } else { // Run 3 GRP object
        o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbPathGrp, bc.timestamp());
        if (grpo == nullptr) {
          LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
        }
        o2::base::Propagator::initFieldFromGRP(grpo);
        o2::base::Propagator::Instance()->setMatLUT(lut);
        LOGF(info, "Setting magnetic field to current %f A for run %d from its GRP CCDB object (type o2::parameters::GRPMagField)", grpo->getL3Current(), bc.runNumber());
      }
      mRunNumber = bc.runNumber();
    }
  } /// end initMagneticFieldCCDB

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;

  } /// End init

  /// Function to select collisions
  template <typename T>
  bool isSelectedCollision(const T& collision)
  {
    /// Trigger selection
    if (eventSelection && !(isRun3 ? collision.sel8() : collision.sel7())) { // currently only sel8 is defined for run3
      return false;
    }
    return true;
  }

  template <typename T4>
  void writeVarTreeColl(const T4& collision, uint64_t timeColl, double timestamp, double timeDiff, int tfID)
  {
    if (writeTree) {
      /// Filling the tree
      rowKFCollisions(collision.posX(),
                      collision.posY(),
                      collision.posZ(),
                      collision.covXX(),
                      collision.covYY(),
                      collision.covZZ(),
                      collision.numContrib(),
                      collision.multNTracksPV(),
                      collision.chi2(),
                      runNumber,
                      timeColl,
                      timestamp,
                      timeDiff,
                      collision.bcId(),
                      tfID);
    }
  }
  /// Process function for data
  void processCollisions(CollisionTableData const& collisions, aod::BCsWithTimestamps const&)
  {
    uint64_t timeColl = 0;
    tfID = tfID + 1;
    LOGP(info, "processing TF {}", tfID);

    int bc0 = 0;
    for (auto& collisionIndex : collisions) {
      auto bc = collisionIndex.bc_as<aod::BCsWithTimestamps>();
      if (bc0 == 0) {
        bc0 = bc.globalBC();
      }
      // LOGP(info, "{} vs {}", collisionIndex.bcId(), bc.globalIndex());
      if (runNumber != bc.runNumber()) {
        initMagneticFieldCCDB(bc, runNumber, ccdb, isRun3 ? ccdbPathGrpMag : ccdbPathGrp, lut, isRun3);
        magneticField = o2::base::Propagator::Instance()->getNominalBz();
      }
      /// Apply event selection
      if (!isSelectedCollision(collisionIndex)) {
        continue;
      }
      timeColl = uint64_t((bc.globalBC() - bc0) * o2::constants::lhc::LHCBunchSpacingNS) + collisionIndex.collisionTime();
      writeVarTreeColl(collisionIndex, timeColl, bc.timestamp(), collisionIndex.collisionTime(), tfID);
    }
  }
  PROCESS_SWITCH(qaKFEvent, processCollisions, "process collision", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<qaKFEventTrack>(cfgc));
  workflow.push_back(adaptAnalysisTask<qaKFEvent>(cfgc));
  return workflow;
}
