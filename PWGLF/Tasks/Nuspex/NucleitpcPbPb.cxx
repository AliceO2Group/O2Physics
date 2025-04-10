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
/// \author jaideep tanwar <jaideep.tanwar@cern.ch>
///
#include <limits>
#include <vector>
#include <string>
#include <cmath>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TF1.h>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "DetectorsBase/Propagator.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "Common/DataModel/PIDResponse.h"
#include "TRandom3.h"
#include "Common/DataModel/CollisionAssociationTables.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

using CollisionsFull = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As>;
/*
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As>;
*/
using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, o2::aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCDe, aod::pidTPCTr, aod::pidTPCHe, aod::pidTPCAl, aod::pidTOFmass, aod::pidTOFbeta>;
namespace
{
static const int number_of_particles = 6;
static const std::vector<std::string> particleNames{"pion", "proton", "deuteron", "triton", "helium3", "alpha"};
static const std::vector<std::string> antiparticleNames{"anti-pion", "anti-proton", "anti-deuteron", "anti-triton", "anti-helium3", "anti-alpha"};
static const std::vector<int> particlePdgCodes{211, 2212, o2::constants::physics::kDeuteron, o2::constants::physics::kTriton, o2::constants::physics::kHelium3, o2::constants::physics::kAlpha};
static const std::vector<double> particleMasses{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron, o2::constants::physics::MassTriton, o2::constants::physics::MassHelium3, o2::constants::physics::MassAlpha};
static const std::vector<int> particleCharge{1, 1, 1, 1, 2, 2};
const int no_BBparam = 6;
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
// default bethbloch parameters
constexpr double betheBlochDefault[number_of_particles][no_BBparam]{
  {13.611469, 3.598765, -0.021138, 2.039562, 0.651040, 0.09},    // pion
  {5.393020, 7.859534, 0.004048, 2.323197, 1.609307, 0.09},      // proton
  {5.393020, 7.859534, 0.004048, 2.323197, 1.609307, 0.09},      // deuteron
  {5.393020, 7.859534, 0.004048, 2.323197, 1.609307, 0.09},      // triton
  {-126.557359, -0.858569, 1.111643, 1.210323, 2.656374, 0.09},  // helium3
  {-126.557359, -0.858569, 1.111643, 1.210323, 2.656374, 0.09}}; // alpha
const int no_trackcuts = 15;
static const std::vector<std::string> trackPIDsettingsNames{"useBBparams", "minITSnCls", "minTPCnCls", "maxTPCchi2", "maxITSchi2", "minRigidity", "maxRigidity", "maxTPCnSigma", "TOFrequiredabove", "minTOFmass", "maxTOFmass", "minDcaToPvXY", "minDcaToPvZ", "minITSclsSize", "maxITSclsSize"};
constexpr double trackPIDsettings[number_of_particles][no_trackcuts]{
  {0, 0, 60, 3.0, 100, 0.15, 1.2, 3.0, 1, 0, 100, 0., 0., 0., 1000},
  {1, 0, 60, 3.0, 100, 0.20, 4.0, 3.0, 1, 0, 100, 0., 0., 0., 1000},
  {1, 0, 60, 3.0, 100, 0.50, 5.0, 3.0, 1, 0, 100, 0., 0., 0., 1000},
  {1, 0, 60, 3.0, 100, 0.50, 5.0, 3.0, 1, 0, 100, 0., 0., 0., 1000},
  {1, 0, 60, 3.0, 100, 0.50, 5.0, 3.0, 1, 0, 100, 0., 0., 0., 1000},
  {1, 0, 60, 3.0, 100, 0.50, 5.0, 3.0, 1, 0, 100, 0., 0., 0., 1000}};
struct Particle {
  TString name;
  int pdgCode, charge;
  double mass, resolution;
  std::vector<double> betheParams;
  bool active;
  Particle(std::string name_, int pdgCode_, double mass_, int charge_, LabeledArray<double> bethe) : name(name_), pdgCode(pdgCode_), charge(charge_), mass(mass_), active(false)
  {
    resolution = bethe.get(name, "resolution");
    betheParams.clear();
    for (unsigned int i = 0; i < 5; i++)
      betheParams.push_back(bethe.get(name, i));
  }
}; // struct Particle
//----------------------------------------------------------------------------------------------------------------
std::vector<std::shared_ptr<TH2>> hDeDx;
std::vector<std::shared_ptr<TH2>> hDeDxanti;
std::vector<std::shared_ptr<TH2>> hnsigma_pt;
std::vector<std::shared_ptr<TH2>> hnsigma_ptanti;
std::vector<std::shared_ptr<TH2>> hdcaXY_pt;
std::vector<std::shared_ptr<TH2>> hdcaXY_ptanti;
std::vector<std::shared_ptr<TH1>> hrapidity;
std::vector<std::shared_ptr<TH2>> hmass_pt;
std::vector<std::shared_ptr<TH2>> hmass_ptanti;
std::vector<std::shared_ptr<TH2>> hdelta_mass;
} // namespace
//----------------------------------------------------------------------------------------------------------------
struct NucleitpcPbPb {
  Preslice<aod::TrackAssoc> perCollision = aod::track_association::collisionId;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> cfgDebug{"cfgDebug", 1, "debug level"};
  Configurable<bool> cfgRigidityCorrection{"cfgRigidityCorrection", false, "apply rigidity correction"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> centcut{"centcut", 80.0f, "centrality cut"};
  Configurable<bool> cfgUsePVcontributors{"cfgUsePVcontributors", true, "use tracks that are PV contibutors"};
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], number_of_particles, no_BBparam, particleNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for light nuclei"};
  Configurable<LabeledArray<double>> cfgTrackPIDsettings{"cfgTrackPIDsettings", {trackPIDsettings[0], number_of_particles, no_trackcuts, particleNames, trackPIDsettingsNames}, "track selection and PID criteria"};
  Configurable<float> maxDcaXYFactor{"maxDcaXYFactor", 2.0f, "DCA xy factor"};
  // CCDB
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<double> bField{"bField", -999, "bz field, -999 is automatic"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> pidPath{"pidPath", "", "Path to the PID response object"};
  //--------------------------------------------------------------------------------------------------------------------
  std::vector<Particle> primaryParticles;
  std::vector<float> primVtx, cents;
  bool collHasCandidate, collPassedEvSel;
  int mRunNumber, occupancy;
  float dBz;
  TRandom3 rand;
  double momn;
  //----------------------------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------------------
  void init(InitContext const&)
  {
    mRunNumber = 0;
    dBz = 0;
    rand.SetSeed(0);
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    for (int i = 0; i < number_of_particles; i++) { // create primaryparticles
      primaryParticles.push_back(Particle(particleNames.at(i), particlePdgCodes.at(i), particleMasses.at(i), particleCharge.at(i), cfgBetheBlochParams));
    }
    std::vector<double> ptBinning = {0.1, 0.5, 1.0, 1.5, 2.0, 2.4, 3.2, 4., 5., 6., 8., 10., 12., 14.};
    std::vector<double> etaBinning = {-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    // define histogram axes
    const AxisSpec axisMagField{10, -10., 10., "magnetic field"};
    const AxisSpec axisNev{3, 0., 5., "Number of events"};
    const AxisSpec axisRigidity{4000, -10., 10., "#it{p}^{TPC}/#it{z}"};
    const AxisSpec axisdEdx{30000, 0, 3000, "d#it{E}/d#it{x}"};
    const AxisSpec axisCent{100, 0, 100, "centrality"};
    const AxisSpec axisVtxZ{100, -20, 20, "z"};
    const AxisSpec axisDCAZ{100, -10, 10, "z"};
    //  const AxisSpec axiseta{100, -1, 1, "eta"};
    const AxisSpec axisrapidity{100, -2, 2, "rapidity"};
    AxisSpec axiseta = {etaBinning, "#eta"};
    AxisSpec axismass = {100, -0.5, 15, "mass^{2}"};
    AxisSpec axisdelta_mass = {100, -6, 6, "#delta mass^{2}"};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec dcaXY = {100, -2, 2, "dcaXY"};
    AxisSpec nsigmaAxis = {160, -20, 20, "n#sigma_{#pi^{+}}"};
    // create histograms
    histos.add("histMagField", "histMagField", kTH1F, {axisMagField});
    histos.add("histNev", "histNev", kTH1F, {axisNev});
    histos.add("histVtxZ", "histVtxZ", kTH1F, {axisVtxZ});
    histos.add("histCentFT0A", "histCentFT0A", kTH1F, {axisCent});
    histos.add("histCentFT0C", "histCentFT0C", kTH1F, {axisCent});
    histos.add("histCentFTOC_cut", "histCentFTOC_cut", kTH1F, {axisCent});
    histos.add("histCentFT0M", "histCentFT0M", kTH1F, {axisCent});
    histos.add("histeta", "histeta", kTH1F, {axiseta});
    histos.add("Tof_signal", "Tof_signal", kTH2F, {axisRigidity, {4000, 0.2, 1.2, "#beta"}});
    histos.add("histDcaZVsPtData_particle", "dcaZ vs Pt (particle)", HistType::kTH2F, {{1000, 0, 20}, {1000, -2.5, 2.5, "dca"}});
    histos.add("histDcaXYVsPtData_particle", "dcaXY vs Pt (particle)", HistType::kTH2F, {{1000, 0, 20}, {1000, -2.0, 2.0, "dca"}});
    histos.add("histDcaZVsPtData_antiparticle", "dcaZ vs Pt (antiparticle)", HistType::kTH2F, {{1000, 0, 20}, {1000, -2.5, 2.5, "dca"}});
    histos.add("histDcaXYVsPtData_antiparticle", "dcaXY vs Pt (antiparticle)", HistType::kTH2F, {{1000, 0, 20}, {1000, -2.0, 2.0, "dca"}});
    hDeDx.resize(2 * number_of_particles + 2);
    hDeDxanti.resize(2 * number_of_particles + 2);
    hnsigma_pt.resize(2 * number_of_particles + 2);
    hnsigma_ptanti.resize(2 * number_of_particles + 2);
    hdcaXY_pt.resize(2 * number_of_particles + 2);
    hdcaXY_ptanti.resize(2 * number_of_particles + 2);
    hrapidity.resize(2 * number_of_particles + 2);
    hmass_pt.resize(2 * number_of_particles + 2);
    hmass_ptanti.resize(2 * number_of_particles + 2);
    hdelta_mass.resize(2 * number_of_particles + 2);
    for (int i = 0; i <= number_of_particles; i++) {
      TString histName = i < number_of_particles ? primaryParticles[i].name : "all";
      hDeDx[2 * i] = histos.add<TH2>(Form("full/histdEdx_%s", histName.Data()), ";p_{TPC}/z (GeV/#it{c}); d#it{E}/d#it{x}", HistType::kTH2F, {axisRigidity, axisdEdx});
      hDeDx[2 * i + 1] = histos.add<TH2>(Form("cuts/histdEdx_%s_Cuts", histName.Data()), ";p_{TPC}/z (GeV/#it{c}); d#it{E}/d#it{x}", HistType::kTH2F, {axisRigidity, axisdEdx});
      hDeDxanti[2 * i] = histos.add<TH2>(Form("antifull/histdEdx_%s", histName.Data()), ";p_{TPC}/z (GeV/#it{c}); d#it{E}/d#it{x}", HistType::kTH2F, {axisRigidity, axisdEdx});
      hDeDxanti[2 * i + 1] = histos.add<TH2>(Form("anticuts/histdEdx_%s_Cuts", histName.Data()), ";p_{TPC}/z (GeV/#it{c}); d#it{E}/d#it{x}", HistType::kTH2F, {axisRigidity, axisdEdx});
    }
    for (int i = 0; i < number_of_particles; i++) {
      TString histName = primaryParticles[i].name;
      hnsigma_pt[2 * i] = histos.add<TH2>(Form("histnsigma_pt/histnsigmaTPC_%s", histName.Data()), ";p_T{TPC} (GeV/#it{c}); TPCnsigma", HistType::kTH2F, {ptAxis, nsigmaAxis});
      hnsigma_ptanti[2 * i] = histos.add<TH2>(Form("histnsigma_ptanti/histnsigmaTPC_%s", histName.Data()), ";p_T{TPC} (GeV/#it{c}); TPCnsigma", HistType::kTH2F, {ptAxis, nsigmaAxis});
      hdcaXY_pt[2 * i] = histos.add<TH2>(Form("histdcaXY_pt/histdcaXY_%s", histName.Data()), ";p_T{TPC} (GeV/#it{c}); dcaXY", HistType::kTH2F, {ptAxis, dcaXY});
      hdcaXY_ptanti[2 * i] = histos.add<TH2>(Form("histdcaXY_ptanti/histdcaXY_%s", histName.Data()), ";p_T{TPC} (GeV/#it{c}); dcaXY", HistType::kTH2F, {ptAxis, dcaXY});
      hmass_pt[2 * i] = histos.add<TH2>(Form("histmass_pt/histmass_%s", histName.Data()), ";p_T{TPC} (GeV/#it{c}); mass^{2}", HistType::kTH2F, {ptAxis, axismass});
      hmass_ptanti[2 * i] = histos.add<TH2>(Form("histmass_ptanti/histmass_%s", histName.Data()), ";p_T{TPC} (GeV/#it{c}); mass^{2}", HistType::kTH2F, {ptAxis, axismass});
      hrapidity[2 * i] = histos.add<TH1>(Form("rapidity/histrapidity_%s", histName.Data()), "; rapidity", HistType::kTH1F, {axisrapidity});
      hdelta_mass[2 * i] = histos.add<TH2>(Form("histdelta/histmass_%s", histName.Data()), ";p_T{TPC} (GeV/#it{c}); #Delta mass", HistType::kTH2F, {ptAxis, axisdelta_mass});
    }
  } // completed void init bracket
  //----------------------------------------------------------------------------------------------------------------
  void findprimaryParticles(aod::TrackAssoc const& tracksByColl, TracksFull const& tracks)
  {
    // track loop, store primary candidates in std::vector
    for (const auto& trackId : tracksByColl) {
      const auto& track = tracks.rawIteratorAt(trackId.trackId());
      /*
        if (!track.isPVContributor())
            continue;
     */
      filldedx(track, number_of_particles);
      if (track.sign() > 0) {
        histos.fill(HIST("histDcaZVsPtData_particle"), track.pt(), track.dcaZ());
        histos.fill(HIST("histDcaXYVsPtData_particle"), track.pt(), track.dcaXY());
      }
      if (track.sign() < 0) {
        histos.fill(HIST("histDcaZVsPtData_antiparticle"), track.pt(), track.dcaZ());
        histos.fill(HIST("histDcaXYVsPtData_antiparticle"), track.pt(), track.dcaXY());
      }
      if (std::abs(track.eta()) > cfgCutEta)
        continue;
      histos.fill(HIST("histeta"), track.eta());
      for (size_t i = 0; i < primaryParticles.size(); i++) {
        if (track.tpcNClsFound() < cfgTrackPIDsettings->get(i, "minTPCnCls"))
          continue;
        if (track.tpcChi2NCl() > cfgTrackPIDsettings->get(i, "maxTPCchi2"))
          continue;
        if (track.itsNCls() < cfgTrackPIDsettings->get(i, "minITSnCls"))
          continue;
        if (track.itsChi2NCl() > cfgTrackPIDsettings->get(i, "maxITSchi2"))
          continue;
        if (getMeanItsClsSize(track) < cfgTrackPIDsettings->get(i, "minITSclsSize"))
          continue;
        if (getMeanItsClsSize(track) > cfgTrackPIDsettings->get(i, "maxITSclsSize"))
          continue;
        if (i == 4 || i == 5) {
          momn = 2 * track.pt();
        } else {
          momn = track.pt();
        }
        bool insideDCAxy = (std::abs(track.dcaXY()) <= (maxDcaXYFactor.value * (0.0105f + 0.0350f / pow(track.pt(), 1.1f))));
        if (!(insideDCAxy) || TMath::Abs(track.dcaZ()) > 2)
          continue;
        if (TMath::Abs(getRapidity(track, i)) > 0.5)
          continue;
        fillhsigma(track, i);
        if (std::abs(getTPCnSigma(track, primaryParticles.at(i))) > cfgTrackPIDsettings->get(i, "maxTPCnSigma"))
          continue;
        filldedx(track, i);
        fillhdcaXY(track, i);
        fillhmass(track, i);
        fillhrapidity(track, i);
        fillhdelta_mass(track, i);
        if (getRigidity(track) < cfgTrackPIDsettings->get(i, "minRigidity") || getRigidity(track) > cfgTrackPIDsettings->get(i, "maxRigidity"))
          continue;
        if (cfgTrackPIDsettings->get(i, "TOFrequiredabove") >= 0 && getRigidity(track) > cfgTrackPIDsettings->get(i, "TOFrequiredabove") && (track.mass() < cfgTrackPIDsettings->get(i, "minTOFmass") || track.mass() > cfgTrackPIDsettings->get(i, "maxTOFmass")))
          continue;
        histos.fill(HIST("Tof_signal"), track.sign() * momn, track.beta());
      }
    } // track loop
  }
  //----------------------------------------------------------------------------------------------------------------
  void processData(CollisionsFull const& collisions, TracksFull const& tracks, aod::BCsWithTimestamps const&, aod::TrackAssoc const& tracksColl)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      initCollision(collision);
      if (!collPassedEvSel)
        continue;
      const uint64_t collIdx = collision.globalIndex();
      auto tracksByColl = tracksColl.sliceBy(perCollision, collIdx);
      findprimaryParticles(tracksByColl, tracks);
      if (!collHasCandidate)
        continue;
      if (collision.centFT0C() > centcut)
        continue;
    }
  }
  PROCESS_SWITCH(NucleitpcPbPb, processData, "data analysis", true);
  //----------------------------------------------------------------------------------------------------------------
  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    auto run3grpTimestamp = bc.timestamp();
    dBz = 0;
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grpTimestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      if (bField < -990) {
        // Fetch magnetic field from ccdb for current collision
        dBz = grpo->getNominalL3Field();
        LOG(info) << "Retrieved GRP for timestamp " << run3grpTimestamp << " with magnetic field of " << dBz << " kZG";
      } else {
        dBz = bField;
      }
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grpTimestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grpTimestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      if (bField < -990) {
        // Fetch magnetic field from ccdb for current collision
        dBz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
        LOG(info) << "Retrieved GRP for timestamp " << run3grpTimestamp << " with magnetic field of " << dBz << " kZG";
      } else {
        dBz = bField;
      }
    }
    mRunNumber = bc.runNumber();
  }
  //----------------------------------------------------------------------------------------------------------------
  template <typename T>
  void initCollision(const T& collision)
  {
    collHasCandidate = false;
    histos.fill(HIST("histMagField"), dBz);
    histos.fill(HIST("histNev"), 0.5);
    collPassedEvSel = collision.sel8() && std::abs(collision.posZ()) < 10;
    if (collision.sel8()) {
      histos.fill(HIST("histNev"), 1.5);
      if (std::abs(collision.posZ()) < 10.0000000000000000) {
        histos.fill(HIST("histNev"), 2.5);
      }
    }
    if (collPassedEvSel) {
      histos.fill(HIST("histVtxZ"), collision.posZ());
      histos.fill(HIST("histCentFT0A"), collision.centFT0A());
      histos.fill(HIST("histCentFT0C"), collision.centFT0C());
      histos.fill(HIST("histCentFT0M"), collision.centFT0M());
      if (collision.centFT0C() < centcut) {
        histos.fill(HIST("histCentFTOC_cut"), collision.centFT0C());
      }
    }
    occupancy = collision.trackOccupancyInTimeRange();
    primVtx.assign({collision.posX(), collision.posY(), collision.posZ()});
    cents.assign({collision.centFT0A(), collision.centFT0C(), collision.centFT0M()});
  }
  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  void filldedx(T const& track, int species)
  {
    const float rigidity = getRigidity(track);
    int idx = 2 * species;
    if (species != 6) {
      auto& hist = (track.sign() > 0) ? hDeDx[idx] : hDeDxanti[idx];
      hist->Fill(track.sign() * rigidity, track.tpcSignal());
    } else {
      hDeDx[idx]->Fill(track.sign() * rigidity, track.tpcSignal());
      hDeDxanti[idx]->Fill(track.sign() * rigidity, track.tpcSignal());
    }
    if (track.tpcNClsFound() < 100 || track.itsNCls() < 2)
      return;

    auto& hist2 = (track.sign() > 0) ? hDeDx[idx + 1] : hDeDxanti[idx + 1];
    hist2->Fill(track.sign() * rigidity, track.tpcSignal());
  }
  template <class T>
  void fillhsigma(T const& track, int species)
  {
    if (track.tpcNClsFound() < 100 || track.itsNCls() < 2)
      return;
    int i = species;
    const float tpcNsigma = getTPCnSigma(track, primaryParticles.at(i));
    double momn;
    if (species == 4 || species == 5) {
      momn = 2 * track.pt();
    } else {
      momn = track.pt();
    }
    if (track.sign() > 0) {
      hnsigma_pt[2 * species]->Fill(momn, tpcNsigma);
    }
    if (track.sign() < 0) {
      hnsigma_ptanti[2 * species]->Fill(momn, tpcNsigma);
    }
  }
  template <class T>
  void fillhdcaXY(T const& track, int species)
  {
    if (track.tpcNClsFound() < 100 || track.itsNCls() < 2)
      return;
    double momn;
    if (species == 4 || species == 5) {
      momn = 2 * track.pt();
    } else {
      momn = track.pt();
    }
    const float dcaXY = track.dcaXY();
    if (track.sign() > 0) {
      hdcaXY_pt[2 * species]->Fill(momn, dcaXY);
    }
    if (track.sign() < 0) {
      hdcaXY_ptanti[2 * species]->Fill(momn, dcaXY);
    }
  }
  template <class T>
  void fillhmass(T const& track, int species)
  {
    if (track.tpcNClsFound() < 100 || track.itsNCls() < 2)
      return;
    double mass;
    if (species == 4 || species == 5) {
      mass = 2 * track.mass();
    } else {
      mass = track.mass();
    }
    double momn;
    if (species == 4 || species == 5) {
      momn = 2 * track.pt();
    } else {
      momn = track.pt();
    }
    if (track.sign() > 0) {
      hmass_pt[2 * species]->Fill(momn, mass * mass);
    }
    if (track.sign() < 0) {
      hmass_ptanti[2 * species]->Fill(momn, mass * mass);
    }
  }
  template <class T>
  void fillhdelta_mass(T const& track, int species)
  {
    if (track.tpcNClsFound() < 100 || track.itsNCls() < 2)
      return;
    double mass;
    if (species == 4 || species == 5) {
      mass = 2 * track.mass();
    } else {
      mass = track.mass();
    }

    double delta_mass = (mass - particleMasses[species]);

    hdelta_mass[2 * species]->Fill(track.pt() * particleCharge[species], delta_mass);
  }
  template <class T>
  void fillhrapidity(T const& track, int species)
  {
    if (track.tpcNClsFound() < 100 || track.itsNCls() < 2)
      return;
    double rap = getRapidity(track, species);
    hrapidity[2 * species]->Fill(rap);
  }
  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  float getTPCnSigma(T const& track, Particle const& particle)
  {
    const float rigidity = getRigidity(track);
    if (!track.hasTPC())
      return -999;
    if (particle.name == "pion" && cfgTrackPIDsettings->get("pion", "useBBparams") == 0)
      return track.tpcNSigmaPi();
    if (particle.name == "proton" && cfgTrackPIDsettings->get("proton", "useBBparams") == 0)
      return track.tpcNSigmaPr();
    if (particle.name == "deuteron" && cfgTrackPIDsettings->get("deuteron", "useBBparams") == 0)
      return track.tpcNSigmaDe();
    if (particle.name == "triton" && cfgTrackPIDsettings->get("triton", "useBBparams") == 0)
      return track.tpcNSigmaTr();
    if (particle.name == "helium3" && cfgTrackPIDsettings->get("helium3", "useBBparams") == 0)
      return track.tpcNSigmaHe();
    if (particle.name == "alpha" && cfgTrackPIDsettings->get("alpha", "useBBparams") == 0)
      return track.tpcNSigmaAl();
    double expBethe{tpc::BetheBlochAleph(static_cast<double>(particle.charge * rigidity / particle.mass), particle.betheParams[0], particle.betheParams[1], particle.betheParams[2], particle.betheParams[3], particle.betheParams[4])};
    double expSigma{expBethe * particle.resolution};
    float sigmaTPC = static_cast<float>((track.tpcSignal() - expBethe) / expSigma);
    return sigmaTPC;
  }
  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  float getMeanItsClsSize(T const& track)
  {
    int sum = 0, n = 0;
    for (int i = 0; i < 8; i++) {
      sum += (track.itsClusterSizes() >> (4 * i) & 15);
      if (track.itsClusterSizes() >> (4 * i) & 15)
        n++;
    }
    return n > 0 ? static_cast<float>(sum) / n : 0.f;
  }
  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  float getRigidity(T const& track)
  {
    if (!cfgRigidityCorrection)
      return track.tpcInnerParam();
    bool hePID = track.pidForTracking() == o2::track::PID::Helium3 || track.pidForTracking() == o2::track::PID::Alpha;
    return hePID ? track.tpcInnerParam() / 2 : track.tpcInnerParam();
  }
  template <class T>
  float getRapidity(T const& track, int species)
  {
    double momn;
    TLorentzVector lorentzVector_particle;
    if (species == 4 || species == 5) {
      momn = 2 * track.pt();
    } else {
      momn = track.pt();
    }
    lorentzVector_particle.SetPtEtaPhiM(momn, track.eta(), track.phi(), particleMasses[species]);
    return lorentzVector_particle.Rapidity();
  }
}; // end of the task here
//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NucleitpcPbPb>(cfgc)};
}
