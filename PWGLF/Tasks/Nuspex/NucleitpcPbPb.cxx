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
/// \file NucleitpcPbPb.cxx
/// \brief Analysis task for light nuclei spectra in Pb–Pb collisions using TPC
/// \author Jaideep Tanwar, <Jaideep.tanwar@cern.ch>
/// \since Jan 2025
///
#include <limits>
#include <string>
#include <vector>
#include <Math/Vector4D.h>
#include <TRandom3.h>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "CCDB/BasicCCDBManager.h"
#include "DCAFitter/DCAFitterN.h"
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using CollisionsFull = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As>;
using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, o2::aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCDe, aod::pidTPCTr, aod::pidTPCHe, aod::pidTPCAl, aod::pidTOFmass>;
//---------------------------------------------------------------------------------------------------------------------------------
namespace
{
static const int nParticles = 6;
static const std::vector<std::string> particleNames{"pion", "proton", "deuteron", "triton", "helion", "alpha"};
static const std::vector<int> particlePdgCodes{211, 2212, o2::constants::physics::kDeuteron, o2::constants::physics::kTriton, o2::constants::physics::kHelium3, o2::constants::physics::kAlpha};
static const std::vector<double> particleMasses{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron, o2::constants::physics::MassTriton, o2::constants::physics::MassHelium3, o2::constants::physics::MassAlpha};
static const std::vector<int> particleCharge{1, 1, 1, 1, 2, 2};
const int nBetheParams = 6;
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
constexpr double betheBlochDefault[nParticles][nBetheParams]{
  {13.611469, 3.598765, -0.021138, 2.039562, 0.651040, 0.09},    // pion
  {5.393020, 7.859534, 0.004048, 2.323197, 1.609307, 0.09},      // proton
  {5.393020, 7.859534, 0.004048, 2.323197, 1.609307, 0.09},      // deuteron
  {5.393020, 7.859534, 0.004048, 2.323197, 1.609307, 0.09},      // triton
  {-126.557359, -0.858569, 1.111643, 1.210323, 2.656374, 0.09},  // helion
  {-126.557359, -0.858569, 1.111643, 1.210323, 2.656374, 0.09}}; // alpha
const int nTrkSettings = 15;
static const std::vector<std::string> trackPIDsettingsNames{"useBBparams", "minITSnCls", "minTPCnCls", "maxTPCchi2", "maxITSchi2", "minRigidity", "maxRigidity", "maxTPCnSigma", "TOFrequiredabove", "minTOFmass", "maxTOFmass", "maxDcaXY", "maxDcaZ", "minITSclsSize", "maxITSclsSize"};
constexpr double trackPIDsettings[nParticles][nTrkSettings]{
  {0, 0, 60, 3.0, 100, 0.15, 1.2, 2.5, -1, 0, 100, 2., 2., 0., 1000},
  {1, 0, 70, 2.5, 100, 0.20, 4.0, 3.0, -1, 0, 100, 2., 2., 0., 1000},
  {1, 0, 70, 5.0, 100, 0.50, 5.0, 3.0, -1, 0, 100, 2., 2., 0., 1000},
  {1, 0, 70, 5.0, 100, 0.50, 5.0, 3.0, -1, 0, 100, 2., 2., 0., 1000},
  {1, 0, 75, 1.5, 100, 0.50, 5.0, 3.0, -1, 0, 100, 2., 2., 0., 1000},
  {1, 0, 70, 1.5, 100, 0.50, 5.0, 3.0, -1, 0, 100, 2., 2., 0., 1000}};
struct PrimParticles {
  TString name;
  int pdgCode, charge;
  double mass, resolution;
  std::vector<double> betheParams;
  bool active;
  PrimParticles(std::string name_, int pdgCode_, double mass_, int charge_, LabeledArray<double> bethe) : name(name_), pdgCode(pdgCode_), charge(charge_), mass(mass_), active(false)
  {
    resolution = bethe.get(name, "resolution");
    betheParams.clear();
    constexpr unsigned int kNSpecies = 5;
    for (unsigned int i = 0; i < kNSpecies; i++)
      betheParams.push_back(bethe.get(name, i));
  }
}; // struct PrimParticles
//----------------------------------------------------------------------------------------------------------------
std::vector<std::shared_ptr<TH2>> hDeDx;
std::vector<std::shared_ptr<TH2>> hNsigmaPt;
std::vector<std::shared_ptr<TH2>> hmass;
} // namespace
//----------------------------------------------------------------------------------------------------------------
struct NucleitpcPbPb {
  Preslice<aod::TrackAssoc> perCollision = aod::track_association::collisionId;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> cfgDebug{"cfgDebug", 1, "debug level"};
  Configurable<bool> cfgRigidityCorrection{"cfgRigidityCorrection", false, "apply rigidity correction"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.9f, "Eta range for tracks"};
  Configurable<bool> cfgUsePVcontributors{"cfgUsePVcontributors", true, "use tracks that are PV contibutors"};
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], nParticles, nBetheParams, particleNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for light nuclei"};
  Configurable<LabeledArray<double>> cfgTrackPIDsettings{"cfgTrackPIDsettings", {trackPIDsettings[0], nParticles, nTrkSettings, particleNames, trackPIDsettingsNames}, "track selection and PID criteria"};
  Configurable<bool> cfgFillDeDxWithoutCut{"cfgFillDeDxWithoutCut", false, "Fill without cut beth bloch"};
  Configurable<bool> cfgFillDeDxWithCut{"cfgFillDeDxWithCut", false, "Fill with cut beth bloch"};
  Configurable<bool> cfgFillnsigma{"cfgFillnsigma", false, "Fill n-sigma histograms"};
  Configurable<bool> cfgFillmass{"cfgFillmass", true, "Fill mass histograms"};
  Configurable<float> centcut{"centcut", 80.0f, "centrality cut"};
  Configurable<float> cfgCutRapidity{"cfgCutRapidity", 0.5f, "Rapidity range"};
  Configurable<float> cfgZvertex{"cfgZvertex", 10, "Min Z Vertex"};
  Configurable<float> cfgtpcNClsFound{"cfgtpcNClsFound", 100.0f, "min. no. of tpcNClsFound"};
  Configurable<float> cfgitsNCls{"cfgitsNCls", 2.0f, "min. no. of itsNCls"};

  // CCDB
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<double> bField{"bField", -999, "bz field, -999 is automatic"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> pidPath{"pidPath", "", "Path to the PID response object"};
  std::vector<PrimParticles> primaryParticles;
  std::vector<float> primVtx, cents;
  bool collHasCandidate, collPassedEvSel;
  int mRunNumber, occupancy;
  float dBz, momn;
  TRandom3 rand;
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
    for (int i = 0; i < nParticles; i++) { // create primaryParticles
      primaryParticles.push_back(PrimParticles(particleNames.at(i), particlePdgCodes.at(i), particleMasses.at(i), particleCharge.at(i), cfgBetheBlochParams));
    }
    // define histogram axes
    const AxisSpec axisMagField{10, -10., 10., "magnetic field"};
    const AxisSpec axisNev{3, 0., 3., "Number of events"};
    const AxisSpec axisRigidity{4000, -10., 10., "#it{p}^{TPC}/#it{z}"};
    const AxisSpec axisdEdx{2000, 0, 2000, "d#it{E}/d#it{x}"};
    const AxisSpec axisCent{100, 0, 100, "centrality"};
    const AxisSpec axisVtxZ{100, -20, 20, "z"};
    const AxisSpec ptAxis{100, 0, 20, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axiseta{100, -1, 1, "eta"};
    const AxisSpec axisrapidity{100, -2, 2, "rapidity"};
    const AxisSpec axismass{100, 0, 20, "mass^{2}"};
    AxisSpec nsigmaAxis = {160, -20, 20, "n#sigma_{#pi^{+}}"};
    // create histograms
    histos.add("histeta", "histeta", kTH1F, {axiseta});
    histos.add("histCentFTOC_cut", "histCentFTOC_cut", kTH1F, {axisCent});
    histos.add("Tof_signal", "Tof_signal", kTH2F, {axisRigidity, {4000, 0.2, 1.2, "#beta"}});
    histos.add("histDcaZVsPtData_particle", "dcaZ vs Pt (particle)", HistType::kTH2F, {{1000, 0, 20}, {1000, -2.5, 2.5, "dca"}});
    histos.add("histDcaXYVsPtData_particle", "dcaXY vs Pt (particle)", HistType::kTH2F, {{1000, 0, 20}, {1000, -2.0, 2.0, "dca"}});
    histos.add("histDcaZVsPtData_antiparticle", "dcaZ vs Pt (antiparticle)", HistType::kTH2F, {{1000, 0, 20}, {1000, -2.5, 2.5, "dca"}});
    histos.add("histDcaXYVsPtData_antiparticle", "dcaXY vs Pt (antiparticle)", HistType::kTH2F, {{1000, 0, 20}, {1000, -2.0, 2.0, "dca"}});
    histos.add("histMagField", "histMagField", kTH1F, {axisMagField});
    histos.add("histNev", "histNev", kTH1F, {axisNev});
    histos.add("histVtxZ", "histVtxZ", kTH1F, {axisVtxZ});
    histos.add("histCentFT0C", "histCentFT0C", kTH1F, {axisCent});
    histos.add("histCentFT0M", "histCentFT0M", kTH1F, {axisCent});
    hDeDx.resize(2 * nParticles + 2);
    hNsigmaPt.resize(2 * nParticles + 2);
    hmass.resize(2 * nParticles + 2);
    for (int i = 0; i < nParticles + 1; i++) {
      TString histName = i < nParticles ? primaryParticles[i].name : "all";
      if (cfgFillDeDxWithoutCut) {
        hDeDx[2 * i] = histos.add<TH2>(Form("dedx/histdEdx_%s", histName.Data()), ";p_{TPC}/z (GeV/#it{c}); d#it{E}/d#it{x}", HistType::kTH2F, {axisRigidity, axisdEdx});
      }
      if (cfgFillDeDxWithCut) {
        hDeDx[2 * i + 1] = histos.add<TH2>(Form("dedx/histdEdx_%s_Cuts", histName.Data()), ";p_{TPC}/z (GeV/#it{c}); d#it{E}/d#it{x}", HistType::kTH2F, {axisRigidity, axisdEdx});
      }
    }
    for (int i = 0; i < nParticles; i++) {
      TString histName = primaryParticles[i].name;
      if (cfgFillnsigma) {
        hNsigmaPt[2 * i] = histos.add<TH2>(Form("histnsigmaTPC_%s", histName.Data()), ";p_T{TPC} (GeV/#it{c}); TPCnsigma", HistType::kTH2F, {ptAxis, nsigmaAxis});
        hNsigmaPt[2 * i + 1] = histos.add<TH2>(Form("histnsigmaTPC_anti_%s", histName.Data()), ";p_T{TPC} (GeV/#it{c}); TPCnsigma", HistType::kTH2F, {ptAxis, nsigmaAxis});
      }
      if (cfgFillmass) {
        hmass[2 * i] = histos.add<TH2>(Form("histmass_pt/histmass_%s", histName.Data()), ";p_T{TPC} (GeV/#it{c}); mass^{2}", HistType::kTH2F, {ptAxis, axismass});
        hmass[2 * i + 1] = histos.add<TH2>(Form("histmass_ptanti/histmass_%s", histName.Data()), ";p_T{TPC} (GeV/#it{c}); mass^{2}", HistType::kTH2F, {ptAxis, axismass});
      }
    }
  }
  //----------------------------------------------------------------------------------------------------------------
  void findprimaryParticles(aod::TrackAssoc const& tracksByColl, TracksFull const& tracks)
  {
    // track loop, store daughter candidates in std::vector
    for (const auto& trackId : tracksByColl) {
      const auto& track = tracks.rawIteratorAt(trackId.trackId());
      filldedx(track, nParticles);
      if (std::abs(track.eta()) > cfgCutEta)
        continue;
      histos.fill(HIST("histeta"), track.eta());
      for (size_t i = 0; i < primaryParticles.size(); i++) {
        if (std::abs(getRapidity(track, i)) > cfgCutRapidity)
          continue;
        bool insideDCAxy = (std::abs(track.dcaXY()) <= (cfgTrackPIDsettings->get(i, "maxDcaXY") * (0.0105f + 0.0350f / std::pow(track.pt(), 1.1f))));
        if (!(insideDCAxy) || std::abs(track.dcaZ()) > cfgTrackPIDsettings->get(i, "maxDcaZ"))
          continue;
        if (track.sign() > 0) {
          histos.fill(HIST("histDcaZVsPtData_particle"), track.pt(), track.dcaZ());
          histos.fill(HIST("histDcaXYVsPtData_particle"), track.pt(), track.dcaXY());
        }
        if (track.sign() < 0) {
          histos.fill(HIST("histDcaZVsPtData_antiparticle"), track.pt(), track.dcaZ());
          histos.fill(HIST("histDcaXYVsPtData_antiparticle"), track.pt(), track.dcaXY());
        }
        if (track.tpcNClsFound() < cfgTrackPIDsettings->get(i, "minTPCnCls"))
          continue;
        if (track.tpcChi2NCl() > cfgTrackPIDsettings->get(i, "maxTPCchi2"))
          continue;
        if (track.itsNCls() < cfgTrackPIDsettings->get(i, "minITSnCls"))
          continue;
        if (track.itsChi2NCl() > cfgTrackPIDsettings->get(i, "maxITSchi2"))
          continue;
        fillnsigma(track, i);
        if (std::abs(getTPCnSigma(track, primaryParticles.at(i))) > cfgTrackPIDsettings->get(i, "maxTPCnSigma"))
          continue;
        filldedx(track, i);
        fillhmass(track, i);
        if (getMeanItsClsSize(track) < cfgTrackPIDsettings->get(i, "minITSclsSize"))
          continue;
        if (getMeanItsClsSize(track) > cfgTrackPIDsettings->get(i, "maxITSclsSize"))
          continue;
        if (getRigidity(track) < cfgTrackPIDsettings->get(i, "minRigidity") || getRigidity(track) > cfgTrackPIDsettings->get(i, "maxRigidity"))
          continue;
        if (cfgTrackPIDsettings->get(i, "TOFrequiredabove") >= 0 && getRigidity(track) > cfgTrackPIDsettings->get(i, "TOFrequiredabove") && (track.mass() < cfgTrackPIDsettings->get(i, "minTOFmass") || track.mass() > cfgTrackPIDsettings->get(i, "maxTOFmass")))
          continue;
      }
    } // track loop
  }
  //----------------------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------------------
  void processData(CollisionsFull const& collisions, TracksFull const& tracks, aod::BCsWithTimestamps const&, aod::TrackAssoc const& tracksColl)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      initCollision(collision);
      if (!collPassedEvSel)
        continue;
      if (collision.centFT0C() > centcut)
        continue;
      histos.fill(HIST("histCentFTOC_cut"), collision.centFT0C());
      const uint64_t collIdx = collision.globalIndex();
      auto tracksByColl = tracksColl.sliceBy(perCollision, collIdx);
      findprimaryParticles(tracksByColl, tracks);
      if (!collHasCandidate)
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
    constexpr float kInvalidBField = -990.f;
    auto run3grpTimestamp = bc.timestamp();
    dBz = 0;
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grpTimestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      if (bField < kInvalidBField) {
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
      if (bField < kInvalidBField) {
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
    collPassedEvSel = collision.sel8() && std::abs(collision.posZ()) < cfgZvertex;
    occupancy = collision.trackOccupancyInTimeRange();
    if (collPassedEvSel) {
      histos.fill(HIST("histNev"), 1.5);
      histos.fill(HIST("histVtxZ"), collision.posZ());
      //  histos.fill(HIST("histCentFT0A"), collision.centFT0A());
      histos.fill(HIST("histCentFT0C"), collision.centFT0C());
      histos.fill(HIST("histCentFT0M"), collision.centFT0M());
      // histos.fill(HIST("histEvents"), collision.centFT0C(), occupancy);
    }
    primVtx.assign({collision.posX(), collision.posY(), collision.posZ()});
    cents.assign({collision.centFT0A(), collision.centFT0C(), collision.centFT0M()});
  }
  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  void filldedx(T const& track, int species)
  {
    const float rigidity = getRigidity(track);
    if (cfgFillDeDxWithoutCut) {
      hDeDx[2 * species]->Fill(track.sign() * rigidity, track.tpcSignal());
    }
    if (track.tpcNClsFound() < cfgtpcNClsFound || track.itsNCls() < cfgitsNCls)
      return;
    if (cfgFillDeDxWithCut) {
      hDeDx[2 * species + 1]->Fill(track.sign() * rigidity, track.tpcSignal());
    }
  }
  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  void fillnsigma(T const& track, int species)
  {
    if (track.tpcNClsFound() < cfgtpcNClsFound || track.itsNCls() < cfgitsNCls)
      return;
    if (cfgFillnsigma) {
      int i = species;
      const float tpcNsigma = getTPCnSigma(track, primaryParticles.at(i));
      double momn;
      int speciesHe3 = 4;
      int speciesHe4 = 5;
      if (species == speciesHe3 || species == speciesHe4) {
        momn = 2 * track.pt();
      } else {
        momn = track.pt();
      }
      if (track.sign() > 0) {
        hNsigmaPt[2 * species]->Fill(momn, tpcNsigma);
      }
      if (track.sign() < 0) {
        hNsigmaPt[2 * species + 1]->Fill(momn, tpcNsigma);
      }
    }
  }
  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  void fillhmass(T const& track, int species)
  {
    if (track.tpcNClsFound() < cfgtpcNClsFound || track.itsNCls() < cfgitsNCls)
      return;
    if (cfgFillmass) {
      double mass;
      int speciesHe3 = 4;
      int speciesHe4 = 5;
      if (species == speciesHe3 || species == speciesHe4) {
        mass = 2 * track.mass();
      } else {
        mass = track.mass();
      }
      double momn;
      if (species == speciesHe3 || species == speciesHe4) {
        momn = 2 * track.pt();
      } else {
        momn = track.pt();
      }
      if (track.sign() > 0) {
        hmass[2 * species]->Fill(momn, mass * mass);
      }
      if (track.sign() < 0) {
        hmass[2 * species + 1]->Fill(momn, mass * mass);
      }
    }
  }
  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  float getTPCnSigma(T const& track, PrimParticles const& particle)
  {
    const float rigidity = getRigidity(track);
    if (!track.hasTPC())
      return -999;
    if (particle.name == "pion" && cfgTrackPIDsettings->get("pion", "useBBparams") < 1)
      return cfgTrackPIDsettings->get("pion", "useBBparams") == 0 ? track.tpcNSigmaPi() : 0;
    if (particle.name == "proton" && cfgTrackPIDsettings->get("proton", "useBBparams") < 1)
      return cfgTrackPIDsettings->get("proton", "useBBparams") == 0 ? track.tpcNSigmaPr() : 0;
    if (particle.name == "deuteron" && cfgTrackPIDsettings->get("deuteron", "useBBparams") < 1)
      return cfgTrackPIDsettings->get("deuteron", "useBBparams") == 0 ? track.tpcNSigmaDe() : 0;
    if (particle.name == "triton" && cfgTrackPIDsettings->get("triton", "useBBparams") < 1)
      return cfgTrackPIDsettings->get("triton", "useBBparams") == 0 ? track.tpcNSigmaTr() : 0;
    if (particle.name == "helion" && cfgTrackPIDsettings->get("helion", "useBBparams") < 1)
      return cfgTrackPIDsettings->get("helion", "useBBparams") == 0 ? track.tpcNSigmaHe() : 0;
    if (particle.name == "alpha" && cfgTrackPIDsettings->get("alpha", "useBBparams") < 1)
      return cfgTrackPIDsettings->get("alpha", "useBBparams") == 0 ? track.tpcNSigmaAl() : 0;

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
    constexpr int kNITSLayers = 8;
    for (int i = 0; i < kNITSLayers; i++) {
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
  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  float getRapidity(T const& track, int species)
  {
    using PtEtaPhiMVector = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>;
    double momn;
    int speciesHe3 = 4;
    int speciesHe4 = 5;
    if (species == speciesHe3 || species == speciesHe4) {
      momn = 2 * track.pt();
    } else {
      momn = track.pt();
    }
    PtEtaPhiMVector lorentzVectorParticle(momn, track.eta(), track.phi(), particleMasses[species]);
    return lorentzVectorParticle.Rapidity();
  }
}; // end of the task here
//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NucleitpcPbPb>(cfgc)};
}
