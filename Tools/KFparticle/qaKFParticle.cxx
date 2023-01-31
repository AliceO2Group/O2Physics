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
/// \file   qaKFParticle.cxx
/// \author Annalena Kalteyer <annalena.sophie.kalteyer@cern.ch>, GSI Darmstadt
/// \brief  Task to test the performance of the KFParticle package
///

#include "Tools/KFparticle/qaKFParticle.h"
#include <CCDB/BasicCCDBManager.h>
#include <string>
#include <TDatabasePDG.h>
#include <TPDGCode.h>
#include "TableHelper.h"

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
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/RecoDecay.h"
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

struct qaKFParticle {

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

  /// Histogram Configurables
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 24., 36., 50.0}, ""};

  /// option to select good events
  Configurable<bool> eventSelection{"eventSelection", true, "select good events"}; // currently only sel8 is defined for run3
  /// options to select only specific tracks
  Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  /// singe track selections
  Configurable<float> nSigmaTPCMinPi{"nSigmaTPCMinPi", -3., "min number of sigma in the TPC for pion tracks"};
  Configurable<float> nSigmaTPCMaxPi{"nSigmaTPCMaxPi", 3., "max number of sigma in the TPC for pion tracks"};
  Configurable<float> nSigmaTPCMinKa{"nSigmaTPCMinKa", -3., "min number of sigma in the TPC for kaon tracks"};
  Configurable<float> nSigmaTPCMaxKa{"nSigmaTPCMaxKa", 3., "max number of sigma in the TPC for kaon tracks"};
  Configurable<float> d_pTMin{"d_pTMin", 0.3, "minimum momentum for tracks"};
  Configurable<float> d_etaRange{"d_etaRange", 0.8, "eta Range for tracks"};
  Configurable<float> d_dcaXYTrackPV{"d_dcaXYTrackPV", 2., "DCA XY of the daughter tracks to the PV"};
  Configurable<float> d_dcaZTrackPV{"d_dcaZTrackPV", 10., "DCA Z of the daughter tracks to the PV"};
  /// D0 selections
  Configurable<bool> applySelectionDoWithTopoConst{"applySelectionDoWithTopoConst", true, "Apply selections on the D0 after constraining it to the PV"};
  Configurable<float> d_pTMinD0{"d_pTMinD0", 0., "minimum momentum for D0 candidates"};
  Configurable<float> d_pTMaxD0{"d_pTMaxD0", 36., "maximum momentum for D0 candidates"};
  Configurable<float> d_massMinD0{"d_massMinD0", 1.65, "minimum mass for D0"};
  Configurable<float> d_massMaxD0{"d_massMaxD0", 2.08, "minimum mass for D0"};
  Configurable<float> d_cosPA{"d_cosPA", -1., "minimum cosine Pointing angle for D0"};
  Configurable<float> d_decayLength{"d_decayLength", 0., "minimum decay length for D0"};
  Configurable<float> d_normdecayLength{"d_normdecayLength", 100., "minimum normalised decay length for D0"};
  Configurable<float> d_chi2topoD0{"d_chi2topoD0", 1000., "maximum chi2 topological of D0 to PV"};
  Configurable<float> d_dist3DSVDau{"d_dist3DSVDau", 1000., "maximum geometrical distance 3D daughter tracks at the SV"};
  Configurable<float> d_cosThetaStarPi{"d_cosThetaStarPi", 1000., "maximum cosine theta star of the pion from D0"};
  Configurable<float> d_cosThetaStarKa{"d_cosThetaStarKa", 1000., "maximum cosine theta star of the kaon from D0"};
  Configurable<float> d_distPiToSV{"d_distPiToSV", 1000., "maximum distance Pi to SV"};
  Configurable<float> d_distKaToSV{"d_distKaToSV", 1000., "maximum distance Ka to SV"};
  /// Option to write D0 variables in a tree
  Configurable<double> d_DwnSmplFact{"d_DwnSmplFact", 1., "Downsampling factor for tree"};
  Configurable<bool> writeTree{"writeTree", false, "write daughter variables in a tree"};
  Configurable<bool> writeQAHistograms{"writeQAHistograms", false, "write all QA histograms"};

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

  using CollisionTableData = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
  using TrackTableData = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTOFFullPi, aod::pidTOFFullKa>;

  HistogramRegistry histos;
  /// Table to be produced
  Produces<o2::aod::TreeKF> rowKF;

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
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(ccdbPathGeo);
    }
    runNumber = 0;

    const AxisSpec axisVertexPosX{500, -0.05, 0.05, "X [cm]"};
    const AxisSpec axisVertexPosY{500, -0.05, 0.05, "Y [cm]"};
    const AxisSpec axisVertexPosZ{100, -20., 20., "Z [cm]"};
    const AxisSpec axisVertexNumContrib{160, 0, 160, "Number Of contributors to the PV"};
    const AxisSpec axisVertexCov{100, -0.00005, 0.00005};

    const AxisSpec axisParX{300, -0.1, 0.1, "#it{x} [cm]"};
    const AxisSpec axisParY{200, -0.1, 0.1, "#it{y} [cm]"};
    const AxisSpec axisParZ{200, -20., 20., "#it{z} [cm]"};
    const AxisSpec axisParPX{binsPt, "#it{p}_{x} [GeV/c]"};
    const AxisSpec axisParPY{binsPt, "#it{p}_{y} [GeV/c]"};
    const AxisSpec axisParPZ{binsPt, "#it{p}_{z} [GeV/c]"};

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

    histos.add("Tracks/dcaXYToPV", "distance of closest approach in #it{xy} plane;#it{dcaXY} [cm];", kTH1D, {{200, -0.15, 0.15}});
    histos.add("Tracks/dcaZToPV", "distance of closest approach in #it{z} plane;#it{dcaXY} [cm];", kTH1D, {{200, -10., 10.}});

    /// tracks
    histos.add("TracksKFPi/x", "track #it{x} position at dca in local coordinate system", kTH1D, {axisParX});
    histos.add("TracksKFPi/y", "track #it{y} position at dca in local coordinate system", kTH1D, {axisParY});
    histos.add("TracksKFPi/z", "track #it{z} position at dca in local coordinate system", kTH1D, {axisParZ});
    histos.add("TracksKFPi/px", "track #it{p_{x}} momentum at dca in local coordinate system", kTH1D, {axisParPX});
    histos.add("TracksKFPi/py", "track #it{p_{y}} momentum at dca in local coordinate system", kTH1D, {axisParPY});
    histos.add("TracksKFPi/pz", "track #it{p_{z}} momentum at dca in local coordinate system", kTH1D, {axisParPZ});
    histos.add("TracksKFPi/chi2perNDF", "Chi2/NDF of the track;#it{chi2/ndf};", kTH1D, {{200, 0., 25.}});
    histos.add("TracksKFPi/dcaXYToPV", "distance of closest approach in #it{xy} plane;#it{dcaXY} [cm];", kTH1D, {{200, -0.15, 0.15}});
    histos.add("TracksKFPi/dcaToPV", "distance of closest approach ;#it{dca} [cm];", kTH1D, {{200, -0.15, 0.15}});
    histos.add("TracksKFPi/dcaToPVLargeRange", "distance of closest approach ;#it{dca} [cm];", kTH1D, {{200, 0., 20.}});
    histos.add("TracksKFPi/dcaToPVLargeRangeMCBeforeReassignment", "distance of closest approach ;#it{dca} [cm];", kTH1D, {{200, 0., 20.}});
    histos.add("TracksKFPi/dcaToPVLargeRangeMCAfterReassignment", "distance of closest approach ;#it{dca} [cm];", kTH1D, {{200, 0., 20.}});
    histos.add("TracksKFPi/nSigmaTPC", "nSigmaTPC vs P", kTH2D, {{axisParPX}, {100, -6., 6.}}); // Add array in function
    histos.add("TracksKFPi/cosThetaStarPion", "cosine theta star of the pion from D0", kTH1D, {{100, -1., 1.}});
    histos.add("TracksKFPi/deviationPiToSV", "deviation of Pi to SV", kTH1D, {{200, 0., 10.}});
    histos.add("TracksKFPi/deviationPiToPV", "deviation of Pi to PV", kTH1D, {{200, 0., 10.}});
    histos.add("TracksKFPi/dcaToSV", "distance of Pi to SV", kTH1D, {{100, -0.15, 0.15}});
    histos.add("TracksKFPi/dcaXYToSV", "distance xy of Pi to SV", kTH1D, {{100, -0.15, 0.15}});

    histos.add("TracksKFKa/x", "track #it{x} position at dca in local coordinate system", kTH1D, {axisParX});
    histos.add("TracksKFKa/y", "track #it{y} position at dca in local coordinate system", kTH1D, {axisParY});
    histos.add("TracksKFKa/z", "track #it{z} position at dca in local coordinate system", kTH1D, {axisParZ});
    histos.add("TracksKFKa/px", "track #it{p_{x}} momentum at dca in local coordinate system", kTH1D, {axisParPX});
    histos.add("TracksKFKa/py", "track #it{p_{y}} momentum at dca in local coordinate system", kTH1D, {axisParPY});
    histos.add("TracksKFKa/pz", "track #it{p_{z}} momentum at dca in local coordinate system", kTH1D, {axisParPZ});
    histos.add("TracksKFKa/chi2perNDF", "Chi2/NDF of the track;#it{chi2/ndf};", kTH1D, {{200, 0., 25.}});
    histos.add("TracksKFKa/dcaXYToPV", "distance of closest approach in #it{xy} plane;#it{dcaXY} [cm];", kTH1D, {{200, -0.15, 0.15}});
    histos.add("TracksKFKa/dcaToPV", "distance of closest approach ;#it{dca} [cm];", kTH1D, {{200, -0.15, 0.15}});
    histos.add("TracksKFKa/dcaToPVLargeRange", "distance of closest approach ;#it{dca} [cm];", kTH1D, {{200, 0., 20.}});
    histos.add("TracksKFKa/dcaToPVLargeRangeMCBeforeReassignment", "distance of closest approach ;#it{dca} [cm];", kTH1D, {{200, 0., 20.}});
    histos.add("TracksKFKa/dcaToPVLargeRangeMCAfterReassignment", "distance of closest approach ;#it{dca} [cm];", kTH1D, {{200, 0., 20.}});
    histos.add("TracksKFKa/nSigmaTPC", "nSigmaTPC vs P", kTH2D, {{axisParPX}, {100, -6., 6.}}); // Add Array in function
    histos.add("TracksKFKa/cosThetaStarKaon", "cosine theta star of the kaon from D0", kTH1D, {{100, -1., 1}});
    histos.add("TracksKFKa/deviationKaToSV", "deviation of Ka to SV", kTH1D, {{200, 0., 10.}});
    histos.add("TracksKFKa/deviationKaToPV", "deviation of Ka to PV", kTH1D, {{200, 0., 10.}});
    histos.add("TracksKFKa/dcaToSV", "distance of Ka to SV", kTH1D, {{100, -0.15, 0.15}});
    histos.add("TracksKFKa/dcaXYToSV", "distance xy of Ka to SV", kTH1D, {{100, -0.15, 0.15}});

    histos.add("TracksDaughter/deviationDaugtherTracks", "chi2 in 3D of daughter tracks at the SV", kTH1D, {{200, 0., 0.2}});
    histos.add("TracksDaughter/dcaDaugtherTracks", "distance in 3D of daughter tracks at the SV", kTH1D, {{100, 0., 0.01}});
    histos.add("TracksDaughter/d0pid0ka", "product of impact parameters of daughters to the PV", kTH1D, {{100, -0.003, 0.003}});

    /// D0 candidates

    histos.add("DZeroCandTopo/X", "X [cm]", kTH1D, {axisParX});
    histos.add("DZeroCandTopo/Y", "Y [cm]", kTH1D, {axisParY});
    histos.add("DZeroCandTopo/Z", "Z [cm]", kTH1D, {axisParZ});
    histos.add("DZeroCandTopo/E", "E", kTH1D, {{100, 0., 50.}});
    histos.add("DZeroCandTopo/Chi2", "Chi2", kTH1D, {{100, 0., 100.}});
    histos.add("DZeroCandTopo/NDF", "NDF", kTH1D, {{5, 0., 5.}});
    histos.add("DZeroCandTopo/Chi2OverNDF", "Chi2OverNDF", kTH1D, {{100, 0., 25.}});
    histos.add("DZeroCandTopo/p", "momentum", kTH1D, {axisParPX});
    histos.add("DZeroCandTopo/pt", "transverse momentum", kTH1D, {axisParPX});
    histos.add("DZeroCandTopo/eta", "eta", kTH1D, {{100, -2., 2.}});
    histos.add("DZeroCandTopo/phi", "phi", kTH1D, {{100, 0., 3.6}});
    histos.add("DZeroCandTopo/mass", "mass", kTH1D, {{430, 1.65, 2.08}});
    histos.add("DZeroCandTopo/massvspt", "mass vs pt", kTH2D, {{axisParPX}, {430, 1.65, 2.08}});
    histos.add("DZeroCandTopo/decayLength", "decay length [cm]", kTH1D, {{200, 0., 0.5}});
    histos.add("DZeroCandTopo/decayLengthXY", "decay length in xy plane [cm]", kTH1D, {{200, 0., 0.5}});
    histos.add("DZeroCandTopo/cosPA", "cosine of pointing angle", kTH1D, {{100, -1, 1.}});
    histos.add("DZeroCandTopo/lifetime", "life time", kTH1D, {{100, 0., 0.2}});
    histos.add("DZeroCandTopo/massErr", "error mass", kTH1D, {{100, 0., 0.1}});
    histos.add("DZeroCandTopo/decayLengthErr", "decay length error [cm]", kTH1D, {{200, 0., 0.05}});
    histos.add("DZeroCandTopo/dcaDToPV", "distance to PV", kTH1D, {{100, -0.01, 0.01}});
    histos.add("DZeroCandTopo/deviationDToPV", "deviation to PV", kTH1D, {{200, 0., 5.}});
    histos.add("DZeroCandTopo/dcaDToPVXY", "distance to PV in xy plane", kTH1D, {{100, -0.02, 0.02}});
    histos.add("DZeroCandTopo/deviationDToPVXY", "deviation to PV in xy plane", kTH1D, {{200, 0., 5.}});
    auto hSelection = histos.add<TH1>("DZeroCandTopo/Selections", "Selections", kTH1D, {{25, 0.5, 25.5}});
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(1), "All Collisions");
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(2), "Cov < 0");
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(3), "All Tracks");
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(4), "Tr Pt");
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(5), "Tr Eta");
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(6), "Tr DCAXY");
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(7), "Tr DCAZ");
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(8), "Tr Sign");
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(9), "Tr PID");
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(10), "All Candidates");
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(11), "Di Dau SV");
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(12), "CTS Pi DT");
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(13), "CTS Ka DT");
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(14), "Di Pi SV");
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(15), "Di Ka SV");
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(16), "Prod Imp Dau");
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(17), "Pt D");
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(18), "M D");
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(19), "CPA D");
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(20), "Pt DT");
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(21), "M DT 0");
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(22), "CPA DT");
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(23), "DL DT");
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(24), "NDL DT");
    hSelection->GetXaxis()->SetBinLabel(hSelection->FindBin(25), "Chi2T");
    auto hSelectionMC = histos.add<TH1>("DZeroCandTopo/SelectionsMC", "Selections MC", kTH1D, {{5, 0.5, 5.5}});
    hSelectionMC->GetXaxis()->SetBinLabel(hSelectionMC->FindBin(1), "All Tracks");
    hSelectionMC->GetXaxis()->SetBinLabel(hSelectionMC->FindBin(2), "No MC particle");
    hSelectionMC->GetXaxis()->SetBinLabel(hSelectionMC->FindBin(3), "Unmatched");
    hSelectionMC->GetXaxis()->SetBinLabel(hSelectionMC->FindBin(4), "Wrong PV");
    hSelectionMC->GetXaxis()->SetBinLabel(hSelectionMC->FindBin(5), "DCA Z > 10cm");

    if (writeQAHistograms) {
      histos.add("DZeroCandGeo/X", "X [cm]", kTH1D, {axisParX});
      histos.add("DZeroCandGeo/Y", "Y [cm]", kTH1D, {axisParY});
      histos.add("DZeroCandGeo/Z", "Z [cm]", kTH1D, {axisParZ});
      histos.add("DZeroCandGeo/E", "E", kTH1D, {{100, 0., 50.}});
      histos.add("DZeroCandGeo/Chi2", "Chi2", kTH1D, {{100, 0., 100.}});
      histos.add("DZeroCandGeo/NDF", "NDF", kTH1D, {{5, 0., 5.}});
      histos.add("DZeroCandGeo/Chi2OverNDF", "Chi2OverNDF", kTH1D, {{100, 0., 25.}});
      histos.add("DZeroCandGeo/p", "momentum", kTH1D, {axisParPX});
      histos.add("DZeroCandGeo/pt", "transverse momentum", kTH1D, {axisParPX});
      histos.add("DZeroCandGeo/eta", "eta", kTH1D, {{100, -2., 2.}});
      histos.add("DZeroCandGeo/phi", "phi", kTH1D, {{100, 0., 3.6}});
      histos.add("DZeroCandGeo/mass", "mass", kTH1D, {{430, 1.65, 2.08}});
      histos.add("DZeroCandGeo/massvspt", "mass vs pt", kTH2D, {{axisParPX}, {430, 1.65, 2.08}});
      histos.add("DZeroCandGeo/cosPA", "cosine of pointing angle", kTH1D, {{100, -1, 1.}});
      histos.add("DZeroCandGeo/massErr", "error mass", kTH1D, {{100, 0., 0.1}});
      histos.add("DZeroCandGeo/dcaDToPV", "distance to PV", kTH1D, {{100, -0.01, 0.01}});
      histos.add("DZeroCandGeo/deviationDToPV", "deviation to PV", kTH1D, {{200, 0., 5.}});
      histos.add("DZeroCandGeo/dcaDToPVXY", "distance to PV in xy plane", kTH1D, {{100, -0.02, 0.02}});
      histos.add("DZeroCandGeo/deviationDToPVXY", "deviation to PV in xy plane", kTH1D, {{200, 0., 5.}});

      histos.add("DZeroCandTopoAtSV/X", "X [cm]", kTH1D, {axisParX});
      histos.add("DZeroCandTopoAtSV/Y", "Y [cm]", kTH1D, {axisParY});
      histos.add("DZeroCandTopoAtSV/Z", "Z [cm]", kTH1D, {axisParZ});
      histos.add("DZeroCandTopoAtSV/E", "E", kTH1D, {{100, 0., 50.}});
      histos.add("DZeroCandTopoAtSV/Chi2", "Chi2", kTH1D, {{100, 0., 100.}});
      histos.add("DZeroCandTopoAtSV/NDF", "NDF", kTH1D, {{5, 0., 5.}});
      histos.add("DZeroCandTopoAtSV/Chi2OverNDF", "Chi2OverNDF", kTH1D, {{100, 0., 25.}});
      histos.add("DZeroCandTopoAtSV/p", "momentum", kTH1D, {axisParPX});
      histos.add("DZeroCandTopoAtSV/pt", "transverse momentum", kTH1D, {axisParPX});
      histos.add("DZeroCandTopoAtSV/eta", "eta", kTH1D, {{100, -2., 2.}});
      histos.add("DZeroCandTopoAtSV/phi", "phi", kTH1D, {{100, 0., 3.6}});
      histos.add("DZeroCandTopoAtSV/mass", "mass", kTH1D, {{430, 1.65, 2.08}});
      histos.add("DZeroCandTopoAtSV/massvspt", "mass vs pt", kTH2D, {{axisParPX}, {430, 1.65, 2.08}});
      histos.add("DZeroCandTopoAtSV/cosPA", "cosine of pointing angle", kTH1D, {{100, -1, 1.}});
      histos.add("DZeroCandTopoAtSV/massErr", "error mass", kTH1D, {{100, 0., 0.1}});
      histos.add("DZeroCandTopoAtSV/dcaDToPV", "distance to PV", kTH1D, {{100, -0.01, 0.01}});
      histos.add("DZeroCandTopoAtSV/deviationDToPV", "deviation to PV", kTH1D, {{200, 0., 5.}});
      histos.add("DZeroCandTopoAtSV/dcaDToPVXY", "distance to PV in xy plane", kTH1D, {{100, -0.02, 0.02}});
      histos.add("DZeroCandTopoAtSV/deviationDToPVXY", "deviation to PV in xy plane", kTH1D, {{200, 0., 5.}});
    }
  } /// End init

  /// Function to select collisions
  template <typename T>
  bool isSelectedCollision(const T& collision)
  {
    /// Trigger selection
    if (eventSelection && !(isRun3 ? collision.sel8() : collision.sel7())) { // currently only sel8 is defined for run3
      return false;
    }
    /// Reject collisions with negative covariance matrix elemts on the digonal
    if (collision.covXX() < 0. || collision.covYY() < 0. || collision.covZZ() < 0.) {
      histos.fill(HIST("DZeroCandTopo/Selections"), 2.f);
      return false;
    }
    return true;
  }

  /// Function for single track selection
  template <typename T>
  bool isSelectedTracks(const T& track1, const T& track2)
  {
    if (track1.p() < d_pTMin) {
      histos.fill(HIST("DZeroCandTopo/Selections"), 4.f);
      return false;
    }
    if (track2.p() < d_pTMin) {
      histos.fill(HIST("DZeroCandTopo/Selections"), 4.f);
      return false;
    }
    /// Eta range
    if (abs(track1.eta()) > d_etaRange) {
      histos.fill(HIST("DZeroCandTopo/Selections"), 5.f);
      return false;
    }
    if (abs(track2.eta()) > d_etaRange) {
      histos.fill(HIST("DZeroCandTopo/Selections"), 5.f);
      return false;
    }
    /// DCA XY of the daughter tracks to the primaty vertex
    if (track1.dcaXY() > d_dcaXYTrackPV) {
      histos.fill(HIST("DZeroCandTopo/Selections"), 6.f);
      return false;
    }
    /// DCA XY of the daughter tracks to the primaty vertex
    if (track2.dcaXY() > d_dcaXYTrackPV) {
      histos.fill(HIST("DZeroCandTopo/Selections"), 6.f);
      return false;
    }
    /// DCA Z of the daughter tracks to the primaty vertex
    if (track1.dcaZ() > d_dcaZTrackPV) {
      histos.fill(HIST("DZeroCandTopo/Selections"), 7.f);
      return false;
    }
    /// DCA Z of the daughter tracks to the primaty vertex
    if (track2.dcaZ() > d_dcaZTrackPV) {
      histos.fill(HIST("DZeroCandTopo/Selections"), 7.f);
      return false;
    }
    /// reject if the tracks have the same sign
    if (track1.sign() == track2.sign()) {
      histos.fill(HIST("DZeroCandTopo/Selections"), 8.f);
      return false;
    }
    return true;
  }

  template <typename T>
  bool isSelectedDaughters(const T& KFPion, const T& KFKaon, const T& KFDZero, const T& KFPV)
  {
    /// distance 3D daughter tracks at the secondary vertex
    if (KFPion.GetDistanceFromParticle(KFKaon) > d_dist3DSVDau) {
      histos.fill(HIST("DZeroCandTopo/Selections"), 11.f);
      return false;
    }
    /// distance Pi to SV
    if (KFPion.GetDistanceFromVertex(KFDZero) > d_distPiToSV) {
      histos.fill(HIST("DZeroCandTopo/Selections"), 14.f);
      return false;
    }
    /// distance Ka to SV
    if (KFKaon.GetDistanceFromVertex(KFDZero) > d_distKaToSV) {
      histos.fill(HIST("DZeroCandTopo/Selections"), 15.f);
      return false;
    }
    return true;
  }

  template <typename T>
  bool isSelectedDoGeo(const T& KFDZero, const T& KFPV)
  {
    /// Pt selection
    if (KFDZero.GetPt() < d_pTMinD0 || KFDZero.GetPt() > d_pTMaxD0) {
      histos.fill(HIST("DZeroCandTopo/Selections"), 17.f);
      return false;
    }
    /// Mass window selection
    if (KFDZero.GetMass() < d_massMinD0 || KFDZero.GetMass() > d_massMaxD0) {
      histos.fill(HIST("DZeroCandTopo/Selections"), 18.f);
      return false;
    }
    /// cosine pointing angle selection
    if (cpaFromKF(KFDZero, KFPV) < d_cosPA) {
      histos.fill(HIST("DZeroCandTopo/Selections"), 19.f);
      return false;
    }
    return true;
  }

  template <typename T>
  bool isSelectedDoTopo(const T& KFDZero_PV, const T& KFPion, const T& KFKaon, const T& KFDZero_DecayVtx, const T& KFPV)
  {
    /// cosine theta star of the pion from D0
    if (cosThetaStarFromKF(0, 421, 211, 321, KFDZero_PV, KFPion, KFKaon) > d_cosThetaStarPi) {
      histos.fill(HIST("DZeroCandTopo/Selections"), 12.f);
      return false;
    }
    /// cosine theta star of the kaon from D0
    if (cosThetaStarFromKF(1, 421, 211, 321, KFDZero_PV, KFPion, KFKaon) > d_cosThetaStarKa) {
      histos.fill(HIST("DZeroCandTopo/Selections"), 13.f);
      return false;
    }

    /// Pt selection
    if (KFDZero_PV.GetPt() < d_pTMinD0 || KFDZero_PV.GetPt() > d_pTMaxD0) {
      histos.fill(HIST("DZeroCandTopo/Selections"), 20.f);
      return false;
    }
    /// Mass window selection
    if (KFDZero_PV.GetMass() == 0.) {
      histos.fill(HIST("DZeroCandTopo/Selections"), 21.f);
      return false;
    }
    /// cosine pointing angle selection
    if (cpaFromKF(KFDZero_DecayVtx, KFPV) < d_cosPA) {
      histos.fill(HIST("DZeroCandTopo/Selections"), 22.f);
      return false;
    }
    /// decay length selection
    if (KFDZero_PV.GetDecayLength() < d_decayLength) {
      histos.fill(HIST("DZeroCandTopo/Selections"), 23.f);
      return false;
    }
    /// decay length error selection
    float normdecayLength = KFDZero_PV.GetDecayLength() / KFDZero_PV.GetErrDecayLength();
    if (normdecayLength < d_normdecayLength) {
      histos.fill(HIST("DZeroCandTopo/Selections"), 24.f);
      return false;
    }
    /// chi2 topological of DZero to PV
    float chi2topo = KFDZero_PV.GetChi2() / KFDZero_PV.GetNDF();
    if (chi2topo > d_chi2topoD0) {
      histos.fill(HIST("DZeroCandTopo/Selections"), 25.f);
      return false;
    }
    return true;
  }

  template <typename T1, typename T2>
  void fillHistograms(const T1& kfpTrackPi, const T1& kfpTrackKa, const T2& KFPion, const T2& KFKaon, const T2& KFDZero_PV, const T2& KFDZero, const T2& KFPV, const T2& KFDZero_DecayVtx)
  {
    /// fill daughter track parameters
    histos.fill(HIST("TracksKFPi/x"), kfpTrackPi.GetX());
    histos.fill(HIST("TracksKFPi/y"), kfpTrackPi.GetY());
    histos.fill(HIST("TracksKFPi/z"), kfpTrackPi.GetZ());
    histos.fill(HIST("TracksKFPi/px"), kfpTrackPi.GetPx());
    histos.fill(HIST("TracksKFPi/py"), kfpTrackPi.GetPy());
    histos.fill(HIST("TracksKFPi/pz"), kfpTrackPi.GetPz());
    histos.fill(HIST("TracksKFPi/chi2perNDF"), kfpTrackPi.GetChi2perNDF());
    histos.fill(HIST("TracksKFPi/dcaXYToPV"), KFPion.GetDistanceFromVertexXY(KFPV));
    histos.fill(HIST("TracksKFPi/dcaToPV"), KFPion.GetDistanceFromVertex(KFPV));
    histos.fill(HIST("TracksKFPi/dcaToPVLargeRange"), KFPion.GetDistanceFromVertex(KFPV));
    histos.fill(HIST("TracksKFPi/cosThetaStarPion"), cosThetaStarFromKF(0, 421, 211, 321, KFDZero_PV, KFPion, KFKaon));
    histos.fill(HIST("TracksKFPi/deviationPiToSV"), KFPion.GetDeviationFromVertex(KFDZero));
    histos.fill(HIST("TracksKFPi/deviationPiToPV"), KFPion.GetDeviationFromVertex(KFPV));
    histos.fill(HIST("TracksKFPi/dcaToSV"), KFPion.GetDistanceFromVertex(KFDZero));
    histos.fill(HIST("TracksKFPi/dcaXYToSV"), KFPion.GetDistanceFromVertexXY(KFDZero));

    histos.fill(HIST("TracksKFKa/x"), kfpTrackKa.GetX());
    histos.fill(HIST("TracksKFKa/y"), kfpTrackKa.GetY());
    histos.fill(HIST("TracksKFKa/z"), kfpTrackKa.GetZ());
    histos.fill(HIST("TracksKFKa/px"), kfpTrackKa.GetPx());
    histos.fill(HIST("TracksKFKa/py"), kfpTrackKa.GetPy());
    histos.fill(HIST("TracksKFKa/pz"), kfpTrackKa.GetPz());
    histos.fill(HIST("TracksKFKa/chi2perNDF"), kfpTrackKa.GetChi2perNDF());
    histos.fill(HIST("TracksKFKa/dcaXYToPV"), KFKaon.GetDistanceFromVertexXY(KFPV));
    histos.fill(HIST("TracksKFKa/dcaToPV"), KFKaon.GetDistanceFromVertex(KFPV));
    histos.fill(HIST("TracksKFKa/dcaToPVLargeRange"), KFPion.GetDistanceFromVertex(KFPV));
    histos.fill(HIST("TracksKFKa/cosThetaStarKaon"), cosThetaStarFromKF(1, 421, 211, 321, KFDZero_PV, KFPion, KFKaon));
    histos.fill(HIST("TracksKFKa/deviationKaToSV"), KFKaon.GetDeviationFromVertex(KFDZero));
    histos.fill(HIST("TracksKFKa/deviationKaToPV"), KFKaon.GetDeviationFromVertex(KFPV));
    histos.fill(HIST("TracksKFKa/dcaToSV"), KFKaon.GetDistanceFromVertex(KFDZero));
    histos.fill(HIST("TracksKFKa/dcaXYToSV"), KFKaon.GetDistanceFromVertexXY(KFDZero));

    histos.fill(HIST("TracksDaughter/deviationDaugtherTracks"), KFPion.GetDeviationFromParticle(KFKaon));
    histos.fill(HIST("TracksDaughter/dcaDaugtherTracks"), KFPion.GetDistanceFromParticle(KFKaon));
    float d0pid0ka = KFPion.GetDistanceFromVertexXY(KFPV) * KFKaon.GetDistanceFromVertexXY(KFPV);
    histos.fill(HIST("TracksDaughter/d0pid0ka"), d0pid0ka);

    histos.fill(HIST("DZeroCandTopo/X"), KFDZero_PV.GetX());
    histos.fill(HIST("DZeroCandTopo/Y"), KFDZero_PV.GetY());
    histos.fill(HIST("DZeroCandTopo/Z"), KFDZero_PV.GetZ());
    histos.fill(HIST("DZeroCandTopo/E"), KFDZero_PV.GetE());
    histos.fill(HIST("DZeroCandTopo/Chi2"), KFDZero_PV.GetChi2());
    histos.fill(HIST("DZeroCandTopo/NDF"), KFDZero_PV.GetNDF());
    float chi2topo = KFDZero_PV.GetChi2() / KFDZero_PV.GetNDF();
    histos.fill(HIST("DZeroCandTopo/Chi2OverNDF"), chi2topo);
    histos.fill(HIST("DZeroCandTopo/p"), KFDZero_PV.GetP());
    histos.fill(HIST("DZeroCandTopo/pt"), KFDZero_PV.GetPt());
    histos.fill(HIST("DZeroCandTopo/eta"), KFDZero_PV.GetEta());
    histos.fill(HIST("DZeroCandTopo/phi"), KFDZero_PV.GetPhi());
    histos.fill(HIST("DZeroCandTopo/mass"), KFDZero_PV.GetMass());
    histos.fill(HIST("DZeroCandTopo/massvspt"), KFDZero_PV.GetPt(), KFDZero_PV.GetMass());
    histos.fill(HIST("DZeroCandTopo/decayLength"), KFDZero_PV.GetDecayLength());
    histos.fill(HIST("DZeroCandTopo/decayLengthXY"), KFDZero_PV.GetDecayLengthXY());
    histos.fill(HIST("DZeroCandTopo/cosPA"), cpaFromKF(KFDZero_DecayVtx, KFPV));
    histos.fill(HIST("DZeroCandTopo/lifetime"), KFDZero_PV.GetLifeTime());
    histos.fill(HIST("DZeroCandTopo/massErr"), KFDZero_PV.GetErrMass());
    histos.fill(HIST("DZeroCandTopo/decayLengthErr"), KFDZero_PV.GetErrDecayLength());
    histos.fill(HIST("DZeroCandTopo/dcaDToPV"), KFDZero_PV.GetDistanceFromVertex(KFPV));
    histos.fill(HIST("DZeroCandTopo/deviationDToPV"), KFDZero_PV.GetDeviationFromVertex(KFPV));
    histos.fill(HIST("DZeroCandTopo/dcaDToPVXY"), KFDZero_PV.GetDistanceFromVertexXY(KFPV));
    histos.fill(HIST("DZeroCandTopo/deviationDToPVXY"), KFDZero_PV.GetDeviationFromVertexXY(KFPV));

    if (writeQAHistograms) {
      histos.fill(HIST("DZeroCandGeo/X"), KFDZero.GetX());
      histos.fill(HIST("DZeroCandGeo/Y"), KFDZero.GetY());
      histos.fill(HIST("DZeroCandGeo/Z"), KFDZero.GetZ());
      histos.fill(HIST("DZeroCandGeo/E"), KFDZero.GetE());
      histos.fill(HIST("DZeroCandGeo/Chi2"), KFDZero.GetChi2());
      histos.fill(HIST("DZeroCandGeo/NDF"), KFDZero.GetNDF());
      float chi2geo = KFDZero.GetChi2() / KFDZero.GetNDF();
      histos.fill(HIST("DZeroCandGeo/Chi2OverNDF"), chi2geo);
      histos.fill(HIST("DZeroCandGeo/p"), KFDZero.GetP());
      histos.fill(HIST("DZeroCandGeo/pt"), KFDZero.GetPt());
      histos.fill(HIST("DZeroCandGeo/eta"), KFDZero.GetEta());
      histos.fill(HIST("DZeroCandGeo/phi"), KFDZero.GetPhi());
      histos.fill(HIST("DZeroCandGeo/mass"), KFDZero.GetMass());
      histos.fill(HIST("DZeroCandGeo/massvspt"), KFDZero.GetP(), KFDZero.GetMass());
      histos.fill(HIST("DZeroCandGeo/cosPA"), cpaFromKF(KFDZero, KFPV));
      histos.fill(HIST("DZeroCandGeo/massErr"), KFDZero.GetErrMass());
      histos.fill(HIST("DZeroCandGeo/dcaDToPV"), KFDZero.GetDistanceFromVertex(KFPV));
      histos.fill(HIST("DZeroCandGeo/deviationDToPV"), KFDZero.GetDeviationFromVertex(KFPV));
      histos.fill(HIST("DZeroCandGeo/dcaDToPVXY"), KFDZero.GetDistanceFromVertexXY(KFPV));
      histos.fill(HIST("DZeroCandGeo/deviationDToPVXY"), KFDZero.GetDeviationFromVertexXY(KFPV));

      histos.fill(HIST("DZeroCandTopoAtSV/X"), KFDZero_DecayVtx.GetX());
      histos.fill(HIST("DZeroCandTopoAtSV/Y"), KFDZero_DecayVtx.GetY());
      histos.fill(HIST("DZeroCandTopoAtSV/Z"), KFDZero_DecayVtx.GetZ());
      histos.fill(HIST("DZeroCandTopoAtSV/E"), KFDZero_DecayVtx.GetE());
      histos.fill(HIST("DZeroCandTopoAtSV/Chi2"), KFDZero_DecayVtx.GetChi2());
      histos.fill(HIST("DZeroCandTopoAtSV/NDF"), KFDZero_DecayVtx.GetNDF());
      float chi2Topo = KFDZero_DecayVtx.GetChi2() / KFDZero_DecayVtx.GetNDF();
      histos.fill(HIST("DZeroCandTopoAtSV/Chi2OverNDF"), chi2Topo);
      histos.fill(HIST("DZeroCandTopoAtSV/p"), KFDZero_DecayVtx.GetP());
      histos.fill(HIST("DZeroCandTopoAtSV/pt"), KFDZero_DecayVtx.GetPt());
      histos.fill(HIST("DZeroCandTopoAtSV/eta"), KFDZero_DecayVtx.GetEta());
      histos.fill(HIST("DZeroCandTopoAtSV/phi"), KFDZero_DecayVtx.GetPhi());
      histos.fill(HIST("DZeroCandTopoAtSV/mass"), KFDZero_DecayVtx.GetMass());
      histos.fill(HIST("DZeroCandTopoAtSV/massvspt"), KFDZero_DecayVtx.GetPt(), KFDZero_DecayVtx.GetMass());
      histos.fill(HIST("DZeroCandTopoAtSV/massErr"), KFDZero_DecayVtx.GetErrMass());
      histos.fill(HIST("DZeroCandTopoAtSV/dcaDToPV"), KFDZero_DecayVtx.GetDistanceFromVertex(KFPV));
      histos.fill(HIST("DZeroCandTopoAtSV/deviationDToPV"), KFDZero_DecayVtx.GetDeviationFromVertex(KFPV));
      histos.fill(HIST("DZeroCandTopoAtSV/dcaDToPVXY"), KFDZero_DecayVtx.GetDistanceFromVertexXY(KFPV));
      histos.fill(HIST("DZeroCandTopoAtSV/deviationDToPVXY"), KFDZero_DecayVtx.GetDeviationFromVertexXY(KFPV));
    }
  }
  template <typename T1, typename T2, typename T3>
  void writeVarTree(const T1& kfpTrackPi, const T1& kfpTrackKa, const T2& KFPion, const T2& KFKaon, const T2& KFDZero_PV, const T2& KFDZero, const T2& KFPV, const T2& KFDZero_DecayVtx, float nSigmaPosPi, float nSigmaNegPi, float nSigmaPosKa, float nSigmaNegKa, const T3& track1, const int source)
  {

    float d0pid0ka = KFPion.GetDistanceFromVertexXY(KFPV) * KFKaon.GetDistanceFromVertexXY(KFPV);
    float chi2geo = KFDZero.GetChi2() / KFDZero.GetNDF();
    float normdecayLength = KFDZero_PV.GetDecayLength() / KFDZero_PV.GetErrDecayLength();
    float chi2topo = KFDZero_PV.GetChi2() / KFDZero_PV.GetNDF();
    float chi2Topo = KFDZero_DecayVtx.GetChi2() / KFDZero_DecayVtx.GetNDF();
    const double pseudoRndm = track1.pt() * 1000. - (int64_t)(track1.pt() * 1000);
    if (pseudoRndm < d_DwnSmplFact) {
      if (writeTree) {
        /// Filling the D0 tree
        rowKF(KFPion.GetPt(),
              KFKaon.GetPt(),
              KFPion.GetDistanceFromVertexXY(KFPV),
              KFKaon.GetDistanceFromVertexXY(KFPV),
              nSigmaPosPi,
              nSigmaNegPi,
              nSigmaPosKa,
              nSigmaNegKa,
              KFPion.GetDistanceFromVertexXY(KFDZero),
              KFKaon.GetDistanceFromVertexXY(KFDZero),
              cosThetaStarFromKF(0, 421, 211, 321, KFDZero_PV, KFPion, KFKaon),
              cosThetaStarFromKF(1, 421, 211, 321, KFDZero_PV, KFPion, KFKaon),
              KFPion.GetDistanceFromParticle(KFKaon),
              d0pid0ka,
              KFDZero.GetPt(),
              KFDZero.GetMass(),
              cpaFromKF(KFDZero, KFPV),
              KFDZero.GetDistanceFromVertex(KFPV),
              KFDZero.GetDistanceFromVertexXY(KFPV),
              chi2geo,
              KFDZero_PV.GetPt(),
              KFDZero_PV.GetMass(),
              KFDZero_PV.GetDecayLength(),
              KFDZero_PV.GetDecayLengthXY(),
              cpaFromKF(KFDZero_DecayVtx, KFPV),
              KFDZero_PV.GetLifeTime(),
              normdecayLength,
              KFDZero_PV.GetDistanceFromVertex(KFPV),
              KFDZero_PV.GetDistanceFromVertexXY(KFPV),
              chi2topo,
              source,
              KFDZero_DecayVtx.GetPt(),
              KFDZero_DecayVtx.GetMass(),
              KFDZero_DecayVtx.GetDistanceFromVertex(KFPV),
              KFDZero_DecayVtx.GetDistanceFromVertexXY(KFPV),
              chi2Topo);
      }
    }
  }

  /// Process function for data
  void processData(CollisionTableData::iterator const& collision, soa::Filtered<TrackTableData> const& tracks, aod::BCsWithTimestamps const&)
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
    histos.fill(HIST("DZeroCandTopo/Selections"), 1.f);
    histos.fill(HIST("Events/covXX"), collision.covXX());
    histos.fill(HIST("Events/covXY"), collision.covXY());
    histos.fill(HIST("Events/covXZ"), collision.covXZ());
    histos.fill(HIST("Events/covYY"), collision.covYY());
    histos.fill(HIST("Events/covYZ"), collision.covYZ());
    histos.fill(HIST("Events/covZZ"), collision.covZZ());
    /// Apply event selection
    if (!isSelectedCollision(collision)) {
      return;
    }
    /// set KF primary vertex
    KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
    KFParticle KFPV(kfpVertex);

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

    for (auto& [track1, track2] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(tracks, tracks))) {

      histos.fill(HIST("DZeroCandTopo/Selections"), 3.f);
      histos.fill(HIST("DZeroCandTopo/Selections"), 3.f);

      /// Apply single track selection
      if (!isSelectedTracks(track1, track2)) {
        return;
      }

      histos.fill(HIST("Tracks/dcaXYToPV"), track1.dcaXY());
      histos.fill(HIST("Tracks/dcaXYToPV"), track2.dcaXY());
      histos.fill(HIST("Tracks/dcaZToPV"), track1.dcaZ());
      histos.fill(HIST("Tracks/dcaZToPV"), track2.dcaZ());

      KFPTrack kfpTrackPosPi;
      KFPTrack kfpTrackNegPi;
      KFPTrack kfpTrackPosKa;
      KFPTrack kfpTrackNegKa;

      bool CandD0 = false;
      bool CandD0bar = false;

      float nSigmaPosPi = 0;
      float nSigmaNegPi = 0;
      float nSigmaPosKa = 0;
      float nSigmaNegKa = 0;
      int source = 0;

      /// At the moment pT independent TPC selection. Add a minimum and Maximum momentum
      /// Apply TPC+TOF at higher momenta (TOF still uncalibrated for LHC22f).

      /// Select D0 and D0bar candidates
      if (nSigmaTPCMinPi <= track1.tpcNSigmaPi() && track1.tpcNSigmaPi() <= nSigmaTPCMaxPi && nSigmaTPCMinKa <= track2.tpcNSigmaKa() && track2.tpcNSigmaKa() <= nSigmaTPCMaxKa) {
        if (track1.sign() == 1 && track2.sign() == -1) {
          CandD0 = true;
          source = 1;
          kfpTrackPosPi = createKFPTrackFromTrack(track1);
          kfpTrackNegKa = createKFPTrackFromTrack(track2);
          nSigmaPosPi = track1.tpcNSigmaPi();
          nSigmaNegKa = track2.tpcNSigmaKa();
          histos.fill(HIST("TracksKFPi/nSigmaTPC"), kfpTrackPosPi.GetPt(), nSigmaPosPi);
          histos.fill(HIST("TracksKFKa/nSigmaTPC"), kfpTrackNegKa.GetPt(), nSigmaNegKa);
        } else if (track1.sign() == -1 && track2.sign() == 1) {
          CandD0bar = true;
          source = 2;
          kfpTrackNegPi = createKFPTrackFromTrack(track1);
          kfpTrackPosKa = createKFPTrackFromTrack(track2);
          nSigmaNegPi = track1.tpcNSigmaPi();
          nSigmaPosKa = track2.tpcNSigmaKa();
          histos.fill(HIST("TracksKFPi/nSigmaTPC"), kfpTrackNegPi.GetPt(), nSigmaNegPi);
          histos.fill(HIST("TracksKFKa/nSigmaTPC"), kfpTrackPosKa.GetPt(), nSigmaPosKa);
        } else {
          continue;
        }
      }
      if (nSigmaTPCMinKa <= track1.tpcNSigmaKa() && track1.tpcNSigmaKa() <= nSigmaTPCMaxKa && nSigmaTPCMinPi <= track2.tpcNSigmaPi() && track2.tpcNSigmaPi() <= nSigmaTPCMaxPi) {
        if (track1.sign() == 1 && track2.sign() == -1) {
          CandD0bar = true;
          source = 2;
          kfpTrackNegPi = createKFPTrackFromTrack(track2);
          kfpTrackPosKa = createKFPTrackFromTrack(track1);
          nSigmaNegPi = track2.tpcNSigmaPi();
          nSigmaPosKa = track1.tpcNSigmaKa();
          histos.fill(HIST("TracksKFPi/nSigmaTPC"), kfpTrackNegPi.GetPt(), nSigmaNegPi);
          histos.fill(HIST("TracksKFKa/nSigmaTPC"), kfpTrackPosKa.GetPt(), nSigmaPosKa);
        } else if (track1.sign() == -1 && track2.sign() == 1) {
          CandD0 = true;
          source = 1;
          kfpTrackPosPi = createKFPTrackFromTrack(track2);
          kfpTrackNegKa = createKFPTrackFromTrack(track1);
          nSigmaPosPi = track2.tpcNSigmaPi();
          nSigmaNegKa = track1.tpcNSigmaKa();
          histos.fill(HIST("TracksKFPi/nSigmaTPC"), kfpTrackPosPi.GetPt(), nSigmaPosPi);
          histos.fill(HIST("TracksKFKa/nSigmaTPC"), kfpTrackNegKa.GetPt(), nSigmaNegKa);
        } else {
          continue;
        }
      }
      if (!CandD0 && !CandD0bar) {
        histos.fill(HIST("DZeroCandTopo/Selections"), 9.f);
        continue;
      }
      if (CandD0 && CandD0bar) {
        source = 3;
      }

      KFParticle KFPosPion(kfpTrackPosPi, 211);
      KFParticle KFNegPion(kfpTrackNegPi, 211);
      KFParticle KFPosKaon(kfpTrackPosKa, 321);
      KFParticle KFNegKaon(kfpTrackNegKa, 321);

      int NDaughters = 2;

      if (CandD0) {
        KFParticle KFDZero;
        const KFParticle* D0Daughters[2] = {&KFPosPion, &KFNegKaon};
        KFDZero.SetConstructMethod(2);
        KFDZero.Construct(D0Daughters, NDaughters);
        histos.fill(HIST("DZeroCandTopo/Selections"), 10.f);
        /// Apply daughter selection
        if (!isSelectedDaughters(KFPosPion, KFNegKaon, KFDZero, KFPV)) {
          continue;
        }
        /// Apply selection on geometrically reconstructed D0
        if (!isSelectedDoGeo(KFDZero, KFPV)) {
          continue;
        }
        /// Apply a topological constraint of the D0 to the PV.
        /// Parameters will be given at the primary vertex.
        KFParticle KFDZero_PV = KFDZero;
        KFDZero_PV.SetProductionVertex(KFPV);
        /// Transport the D0 after the topological constraint back to the decay vertex
        KFParticle KFDZero_DecayVtx = KFDZero_PV;
        KFDZero_DecayVtx.TransportToDecayVertex();
        if (applySelectionDoWithTopoConst) {
          /// Apply selection on D0 after constraint to the PV
          if (!isSelectedDoTopo(KFDZero_PV, KFPosPion, KFNegKaon, KFDZero_DecayVtx, KFPV)) {
            continue;
          }
        }

        fillHistograms(kfpTrackPosPi, kfpTrackNegKa, KFPosPion, KFNegKaon, KFDZero_PV, KFDZero, KFPV, KFDZero_DecayVtx);
        writeVarTree(kfpTrackPosPi, kfpTrackNegKa, KFPosPion, KFNegKaon, KFDZero_PV, KFDZero, KFPV, KFDZero_DecayVtx, nSigmaPosPi, nSigmaNegPi, nSigmaPosKa, nSigmaNegKa, track1, source);
      }
      if (CandD0bar) {
        KFParticle KFDZeroBar;
        const KFParticle* D0BarDaughters[2] = {&KFNegPion, &KFPosKaon};
        KFDZeroBar.SetConstructMethod(2);
        KFDZeroBar.Construct(D0BarDaughters, NDaughters);
        histos.fill(HIST("DZeroCandTopo/Selections"), 10.f);
        /// Apply daughter selection
        if (!isSelectedDaughters(KFNegPion, KFPosKaon, KFDZeroBar, KFPV)) {
          continue;
        }
        /// Apply selection on geometrically reconstructed D0
        if (!isSelectedDoGeo(KFDZeroBar, KFPV)) {
          continue;
        }
        /// Apply a topological constraint of the D0 to the PV.
        /// Parameters will be given at the primary vertex.
        KFParticle KFDZeroBar_PV = KFDZeroBar;
        KFDZeroBar_PV.SetProductionVertex(KFPV);
        /// Transport the D0 after the topological constraint back to the decay vertex
        KFParticle KFDZeroBar_DecayVtx = KFDZeroBar_PV;
        KFDZeroBar_DecayVtx.TransportToDecayVertex();
        if (applySelectionDoWithTopoConst) {
          /// Apply selection on D0 after constraint to the PV
          if (!isSelectedDoTopo(KFDZeroBar_PV, KFNegPion, KFPosKaon, KFDZeroBar_DecayVtx, KFPV)) {
            continue;
          }
        }
        fillHistograms(kfpTrackNegPi, kfpTrackPosKa, KFNegPion, KFPosKaon, KFDZeroBar_PV, KFDZeroBar, KFPV, KFDZeroBar_DecayVtx);
        writeVarTree(kfpTrackPosPi, kfpTrackNegKa, KFPosPion, KFNegKaon, KFDZeroBar_PV, KFDZeroBar, KFPV, KFDZeroBar_DecayVtx, nSigmaPosPi, nSigmaNegPi, nSigmaPosKa, nSigmaNegKa, track1, source);
      }
    }
  }
  PROCESS_SWITCH(qaKFParticle, processData, "process data", true);

  /// Process function for MC
  using CollisionTableMC = soa::Join<CollisionTableData, aod::McCollisionLabels>;
  using CollisionTableDataMult = soa::Join<aod::Collisions, aod::Mults, aod::McCollisionLabels>;
  using TrackTableMC = soa::Join<TrackTableData, aod::McTrackLabels>;
  Preslice<aod::McCollisionLabels> perMcCollision = aod::mccollisionlabel::mcCollisionId;
  void processMC(CollisionTableMC::iterator const& collision, CollisionTableMC const& collisions, soa::Filtered<TrackTableMC> const& tracks, aod::McParticles const& mcParticles, aod::McCollisions const& mcCollisions, aod::BCsWithTimestamps const&)
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
    histos.fill(HIST("DZeroCandTopo/Selections"), 1.f);
    histos.fill(HIST("Events/covXX"), collision.covXX());
    histos.fill(HIST("Events/covXY"), collision.covXY());
    histos.fill(HIST("Events/covXZ"), collision.covXZ());
    histos.fill(HIST("Events/covYY"), collision.covYY());
    histos.fill(HIST("Events/covYZ"), collision.covYZ());
    histos.fill(HIST("Events/covZZ"), collision.covZZ());
    /// Apply event selection
    if (!isSelectedCollision(collision)) {
      return;
    }
    /// set KF primary vertex
    KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
    KFParticle KFPV(kfpVertex);

    KFPVertex kfpVertexDefault = createKFPVertexFromCollision(collision);
    KFParticle KFPVDefault(kfpVertexDefault);

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

    for (auto& [track1, track2] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(tracks, tracks))) {

      histos.fill(HIST("DZeroCandTopo/Selections"), 3.f);
      histos.fill(HIST("DZeroCandTopo/Selections"), 3.f);

      histos.fill(HIST("DZeroCandTopo/SelectionsMC"), 1.f);

      /// Apply single track selection
      if (!isSelectedTracks(track1, track2)) {
        return;
      }

      /// Check whether the track was assigned to the true MC PV
      auto particle1 = track1.mcParticle();
      auto particle2 = track2.mcParticle();
      auto collMC = particle1.mcCollision();
      auto mcCollID_recoColl = track1.collision_as<CollisionTableMC>().mcCollisionId();
      auto mcCollID_particle = particle1.mcCollisionId();
      bool indexMatchOK = (mcCollID_recoColl == mcCollID_particle);
      if (!indexMatchOK) {
        histos.fill(HIST("DZeroCandTopo/SelectionsMC"), 4.f);
        const auto matchedCollisions = collisions.sliceBy(perMcCollision, collMC.globalIndex());
        int i = 0;
        std::array<float, 5> dcaZ{100, 100, 100, 100, 100};
        float min = 100;
        for (auto matchedCollision : matchedCollisions) {
          dcaZ[i] = abs(matchedCollision.posZ() - collMC.posZ());
          if (i == 0) {
            min = dcaZ[i];
          }
          if (i > 0) {
            if (dcaZ[i] < dcaZ[i - 1]) {
              min = dcaZ[i];
            }
          }

          i = i + 1;
        }
        if (min > 10.) {
          histos.fill(HIST("DZeroCandTopo/SelectionsMC"), 5.f);
        }
        int j = 0;
        for (auto matchedCollision : matchedCollisions) {
          if (i == 1) {
            kfpVertex = createKFPVertexFromCollision(matchedCollision);
            KFParticle KFPVNew(kfpVertex);
            KFPV = KFPVNew;
          }
          if (i > 1) {
            if (abs(matchedCollision.posZ() - collMC.posZ()) == min) {
              kfpVertex = createKFPVertexFromCollision(matchedCollision);
              KFParticle KFPVNew(kfpVertex);
              KFPV = KFPVNew;
            }
          }
          j = j + 1;
        }
      }

      /// Remove fake tracks
      if (!track1.has_mcParticle() || !track2.has_mcParticle()) {
        histos.fill(HIST("DZeroCandTopo/SelectionsMC"), 2.f);
        continue;
      }
      /// Remove unmatched tracks
      if (track1.collisionId() <= 0 || track1.collisionId() <= 0) {
        histos.fill(HIST("DZeroCandTopo/SelectionsMC"), 3.f);
        continue;
      }

      int8_t sign = 0;
      int8_t flag = RecoDecay::OriginType::None;
      int sourceD0 = 0;
      int sourceD0Bar = 0;

      auto indexRec = RecoDecay::getMatchedMCRec(mcParticles, std::array{track1, track2}, 421, array{211, -321}, true, &sign);
      if (indexRec > -1) {
        auto particle = mcParticles.rawIteratorAt(indexRec);
        flag = RecoDecay::getCharmHadronOrigin(mcParticles, particle);
        if (flag == RecoDecay::OriginType::Prompt) {
          sourceD0 |= kPrompt;
          sourceD0Bar |= kPrompt;
        } else if (flag == RecoDecay::OriginType::NonPrompt) {
          sourceD0 |= kNonPrompt;
          sourceD0Bar |= kNonPrompt;
        }
        if (flag != RecoDecay::OriginType::Prompt && flag != RecoDecay::OriginType::NonPrompt) {
          continue;
        }
      } else {
        continue;
      }

      histos.fill(HIST("Tracks/dcaXYToPV"), track1.dcaXY());
      histos.fill(HIST("Tracks/dcaXYToPV"), track2.dcaXY());
      histos.fill(HIST("Tracks/dcaZToPV"), track1.dcaZ());
      histos.fill(HIST("Tracks/dcaZToPV"), track2.dcaZ());

      KFPTrack kfpTrackPosPi;
      KFPTrack kfpTrackNegPi;
      KFPTrack kfpTrackPosKa;
      KFPTrack kfpTrackNegKa;

      bool CandD0 = false;
      bool CandD0bar = false;
      float nSigmaPosPi = 0.;
      float nSigmaNegPi = 0.;
      float nSigmaPosKa = 0.;
      float nSigmaNegKa = 0.;
      /// At the moment pT independent TPC selection. Add a minimum and Maximum momentum
      /// Apply TPC+TOF at higher momenta.

      int pdgMother = mcParticles.rawIteratorAt(indexRec - mcParticles.offset()).pdgCode();
      /// Select D0 and D0bar candidates
      if (nSigmaTPCMinPi <= track1.tpcNSigmaPi() && track1.tpcNSigmaPi() <= nSigmaTPCMaxPi && nSigmaTPCMinKa <= track2.tpcNSigmaKa() && track2.tpcNSigmaKa() <= nSigmaTPCMaxKa) {
        if (track1.sign() == 1 && track2.sign() == -1) {
          CandD0 = true;
          particle1 = track1.mcParticle();
          particle2 = track2.mcParticle();
          if (pdgMother == -421) {
            sourceD0 |= kReflection;
          }
          kfpTrackPosPi = createKFPTrackFromTrack(track1);
          kfpTrackNegKa = createKFPTrackFromTrack(track2);
          nSigmaPosPi = track1.tpcNSigmaPi();
          nSigmaNegKa = track2.tpcNSigmaKa();
          histos.fill(HIST("TracksKFPi/nSigmaTPC"), kfpTrackPosPi.GetPt(), nSigmaPosPi);
          histos.fill(HIST("TracksKFKa/nSigmaTPC"), kfpTrackNegKa.GetPt(), nSigmaNegKa);
        } else if (track1.sign() == -1 && track2.sign() == 1) {
          CandD0bar = true;
          if (pdgMother == 421) {
            sourceD0Bar |= kReflection;
          }
          kfpTrackNegPi = createKFPTrackFromTrack(track1);
          kfpTrackPosKa = createKFPTrackFromTrack(track2);
          nSigmaNegPi = track1.tpcNSigmaPi();
          nSigmaPosKa = track2.tpcNSigmaKa();
          histos.fill(HIST("TracksKFPi/nSigmaTPC"), kfpTrackNegPi.GetPt(), nSigmaNegPi);
          histos.fill(HIST("TracksKFKa/nSigmaTPC"), kfpTrackPosKa.GetPt(), nSigmaPosKa);
        } else {
          continue;
        }
      }
      if (nSigmaTPCMinKa <= track1.tpcNSigmaKa() && track1.tpcNSigmaKa() <= nSigmaTPCMaxKa && nSigmaTPCMinPi <= track2.tpcNSigmaPi() && track2.tpcNSigmaPi() <= nSigmaTPCMaxPi) {
        if (track1.sign() == 1 && track2.sign() == -1) {
          CandD0bar = true;
          if (pdgMother == 421) {
            sourceD0Bar |= kReflection;
          }
          kfpTrackNegPi = createKFPTrackFromTrack(track2);
          kfpTrackPosKa = createKFPTrackFromTrack(track1);
          nSigmaNegPi = track2.tpcNSigmaPi();
          nSigmaPosKa = track1.tpcNSigmaKa();
          histos.fill(HIST("TracksKFPi/nSigmaTPC"), kfpTrackNegPi.GetPt(), nSigmaNegPi);
          histos.fill(HIST("TracksKFKa/nSigmaTPC"), kfpTrackPosKa.GetPt(), nSigmaPosKa);
        } else if (track1.sign() == -1 && track2.sign() == 1) {
          CandD0 = true;
          if (pdgMother == -421) {
            sourceD0 |= kReflection;
          }
          kfpTrackPosPi = createKFPTrackFromTrack(track2);
          kfpTrackNegKa = createKFPTrackFromTrack(track1);
          nSigmaPosPi = track2.tpcNSigmaPi();
          nSigmaNegKa = track1.tpcNSigmaKa();
          histos.fill(HIST("TracksKFPi/nSigmaTPC"), kfpTrackPosPi.GetPt(), nSigmaPosPi);
          histos.fill(HIST("TracksKFKa/nSigmaTPC"), kfpTrackNegKa.GetPt(), nSigmaNegKa);
        } else {
          continue;
        }
      }
      if (!CandD0 && !CandD0bar) {
        histos.fill(HIST("DZeroCandTopo/Selections"), 9.f);
        continue;
      }

      KFParticle KFPosPion(kfpTrackPosPi, 211);
      KFParticle KFNegPion(kfpTrackNegPi, 211);
      KFParticle KFPosKaon(kfpTrackPosKa, 321);
      KFParticle KFNegKaon(kfpTrackNegKa, 321);

      int NDaughters = 2;

      if (CandD0) {
        KFParticle KFDZero;
        const KFParticle* D0Daughters[2] = {&KFPosPion, &KFNegKaon};
        KFDZero.SetConstructMethod(2);
        KFDZero.Construct(D0Daughters, NDaughters);
        histos.fill(HIST("DZeroCandTopo/Selections"), 10.f);
        /// Apply daughter selection
        if (!isSelectedDaughters(KFPosPion, KFNegKaon, KFDZero, KFPV)) {
          continue;
        }
        /// Apply selection on geometrically reconstructed D0
        if (!isSelectedDoGeo(KFDZero, KFPV)) {
          continue;
        }
        /// Apply a topological constraint of the D0 to the PV.
        /// Parameters will be given at the primary vertex.
        KFParticle KFDZero_PV = KFDZero;
        KFDZero_PV.SetProductionVertex(KFPV);
        /// Transport the D0 after the topological constraint back to the decay vertex
        KFParticle KFDZero_DecayVtx = KFDZero_PV;
        KFDZero_DecayVtx.TransportToDecayVertex();
        if (applySelectionDoWithTopoConst) {
          /// Apply selection on D0 after constraint to the PV
          if (!isSelectedDoTopo(KFDZero_PV, KFPosPion, KFNegKaon, KFDZero_DecayVtx, KFPV)) {
            continue;
          }
        }

        fillHistograms(kfpTrackPosPi, kfpTrackNegKa, KFPosPion, KFNegKaon, KFDZero_PV, KFDZero, KFPV, KFDZero_DecayVtx);
        histos.fill(HIST("TracksKFPi/dcaToPVLargeRangeMCBeforeReassignment"), KFPosPion.GetDistanceFromVertex(KFPVDefault));
        histos.fill(HIST("TracksKFPi/dcaToPVLargeRangeMCAfterReassignment"), KFPosPion.GetDistanceFromVertex(KFPV));
        histos.fill(HIST("TracksKFKa/dcaToPVLargeRangeMCBeforeReassignment"), KFNegKaon.GetDistanceFromVertex(KFPVDefault));
        histos.fill(HIST("TracksKFKa/dcaToPVLargeRangeMCAfterReassignment"), KFNegKaon.GetDistanceFromVertex(KFPV));
        writeVarTree(kfpTrackPosPi, kfpTrackNegKa, KFPosPion, KFNegKaon, KFDZero_PV, KFDZero, KFPV, KFDZero_DecayVtx, nSigmaPosPi, nSigmaNegPi, nSigmaPosKa, nSigmaNegKa, track1, sourceD0);
      }
      if (CandD0bar) {
        KFParticle KFDZeroBar;
        const KFParticle* D0BarDaughters[2] = {&KFNegPion, &KFPosKaon};
        KFDZeroBar.SetConstructMethod(2);
        KFDZeroBar.Construct(D0BarDaughters, NDaughters);
        histos.fill(HIST("DZeroCandTopo/Selections"), 10.f);
        /// Apply daughter selection
        if (!isSelectedDaughters(KFNegPion, KFPosKaon, KFDZeroBar, KFPV)) {
          continue;
        }
        /// Apply selection on geometrically reconstructed D0
        if (!isSelectedDoGeo(KFDZeroBar, KFPV)) {
          continue;
        }
        /// Apply a topological constraint of the D0 to the PV.
        /// Parameters will be given at the primary vertex.
        KFParticle KFDZeroBar_PV = KFDZeroBar;
        KFDZeroBar_PV.SetProductionVertex(KFPV);
        /// Transport the D0 after the topological constraint back to the decay vertex
        KFParticle KFDZeroBar_DecayVtx = KFDZeroBar_PV;
        KFDZeroBar_DecayVtx.TransportToDecayVertex();
        if (applySelectionDoWithTopoConst) {
          /// Apply selection on D0 after constraint to the PV
          if (!isSelectedDoTopo(KFDZeroBar_PV, KFNegPion, KFPosKaon, KFDZeroBar_DecayVtx, KFPV)) {
            continue;
          }
        }
        fillHistograms(kfpTrackNegPi, kfpTrackPosKa, KFNegPion, KFPosKaon, KFDZeroBar_PV, KFDZeroBar, KFPV, KFDZeroBar_DecayVtx);
        histos.fill(HIST("TracksKFPi/dcaToPVLargeRangeMCBeforeReassignment"), KFNegPion.GetDistanceFromVertex(KFPVDefault));
        histos.fill(HIST("TracksKFPi/dcaToPVLargeRangeMCAfterReassignment"), KFNegPion.GetDistanceFromVertex(KFPV));
        histos.fill(HIST("TracksKFKa/dcaToPVLargeRangeMCBeforeReassignment"), KFPosKaon.GetDistanceFromVertex(KFPVDefault));
        histos.fill(HIST("TracksKFKa/dcaToPVLargeRangeMCAfterReassignment"), KFPosKaon.GetDistanceFromVertex(KFPV));
        writeVarTree(kfpTrackPosPi, kfpTrackNegKa, KFPosPion, KFNegKaon, KFDZeroBar_PV, KFDZeroBar, KFPV, KFDZeroBar_DecayVtx, nSigmaPosPi, nSigmaNegPi, nSigmaPosKa, nSigmaNegKa, track1, sourceD0Bar);
      }
    }
  }
  PROCESS_SWITCH(qaKFParticle, processMC, "process mc", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<qaKFParticle>(cfgc)};
}
