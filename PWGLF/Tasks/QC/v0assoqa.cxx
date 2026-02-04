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
//
// Strangeness-to-collision association tests
//
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/Utils/strangenessBuilderHelper.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TPCVDriftManager.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <Math/Vector4D.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TProfile.h>

#include <array>
#include <cmath>
#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

// using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCPr>;
using TracksCompleteIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA>;
using TracksCompleteIUMC = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::McTrackLabels>;
using V0DataLabeled = soa::Join<aod::V0Datas, aod::McV0Labels>;
using CascMC = soa::Join<aod::CascDataExt, aod::McCascLabels>;
using TraCascMC = soa::Join<aod::TraCascDatas, aod::McTraCascLabels>;
using RecoedMCCollisions = soa::Join<aod::McCollisions, aod::McCollsExtra>;
using CollisionsWithEvSels = soa::Join<aod::Collisions, aod::EvSels>;

// For MC association in pre-selection
using LabeledTracksExtra = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::McTrackLabels>;

struct v0assoqa {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // helper object
  o2::pwglf::strangenessBuilderHelper straHelper;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  int mRunNumber;
  o2::base::MatLayerCylSet* lut = nullptr;

  // for handling TPC-only tracks (photons)
  o2::aod::common::TPCVDriftManager mVDriftMgr;

  // CCDB options
  struct : ConfigurableGroup {
    std::string prefix = "ccdb";
    Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
    Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  } ccdbConfigurations;

  // V0 building options
  struct : ConfigurableGroup {
    std::string prefix = "v0BuilderOpts";
    Configurable<bool> moveTPCOnlyTracks{"moveTPCOnlyTracks", true, "if dealing with TPC-only tracks, move them according to TPC drift / time info"};

    // baseline conditionals of V0 building
    Configurable<int> minCrossedRows{"minCrossedRows", -1, "minimum TPC crossed rows for daughter tracks"};
    Configurable<float> dcanegtopv{"dcanegtopv", .0, "DCA Neg To PV"};
    Configurable<float> dcapostopv{"dcapostopv", .0, "DCA Pos To PV"};
    Configurable<double> v0cospa{"v0cospa", -2, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0)
    Configurable<float> dcav0dau{"dcav0dau", 10000.0, "DCA V0 Daughters"};
    Configurable<float> v0radius{"v0radius", 0.0, "v0radius"};
    Configurable<float> maxDaughterEta{"maxDaughterEta", 5.0, "Maximum daughter eta (in abs value)"};
  } v0BuilderOpts;

  // cascade building options
  struct : ConfigurableGroup {
    std::string prefix = "cascadeBuilderOpts";
    // conditionals
    Configurable<int> minCrossedRows{"minCrossedRows", 50, "minimum TPC crossed rows for daughter tracks"};
    Configurable<float> dcabachtopv{"dcabachtopv", .05, "DCA Bach To PV"};
    Configurable<float> cascradius{"cascradius", 0.9, "cascradius"};
    Configurable<float> casccospa{"casccospa", 0.95, "casccospa"};
    Configurable<float> dcacascdau{"dcacascdau", 1.0, "DCA cascade Daughters"};
    Configurable<float> lambdaMassWindow{"lambdaMassWindow", .010, "Distance from Lambda mass (does not apply to KF path)"};
    Configurable<float> maxDaughterEta{"maxDaughterEta", 5.0, "Maximum daughter eta (in abs value)"};
  } cascadeBuilderOpts;

  //_______________________________________________________________________
  template <typename TMCParticles>
  int findMotherFromLabels(int const& p1, int const& p2, const int expected_pdg1, const int expected_pdg2, const int expected_mother_pdg, TMCParticles const& mcparticles)
  {
    // encompasses a simple check for labels existing
    if (p1 < 0 || p2 < 0) {
      return -1;
    }
    auto mcParticle1 = mcparticles.rawIteratorAt(p1);
    auto mcParticle2 = mcparticles.rawIteratorAt(p2);
    return (findMother(mcParticle1, mcParticle2, expected_pdg1, expected_pdg2, expected_mother_pdg, mcparticles));
  }

  //_______________________________________________________________________
  template <typename TMCParticle1, typename TMCParticle2, typename TMCParticles>
  int findMother(TMCParticle1 const& p1, TMCParticle2 const& p2, const int expected_pdg1, const int expected_pdg2, const int expected_mother_pdg, TMCParticles const& mcparticles)
  {
    if (p1.globalIndex() == p2.globalIndex())
      return -1;
    if (p1.pdgCode() != expected_pdg1 || p2.pdgCode() != expected_pdg2)
      return -1;
    if (!p1.has_mothers() || !p2.has_mothers())
      return -1;

    int motherid1 = p1.mothersIds()[0];
    auto mother1 = mcparticles.iteratorAt(motherid1);
    int mother1_pdg = mother1.pdgCode();
    int motherid2 = p2.mothersIds()[0];
    auto mother2 = mcparticles.iteratorAt(motherid2);
    int mother2_pdg = mother2.pdgCode();

    if (motherid1 != motherid2 || mother1_pdg != mother2_pdg || mother1_pdg != expected_mother_pdg)
      return -1;

    return motherid1;
  }

  void init(InitContext const&)
  {
    histos.add("hDuplicateCount", "hDuplicateCount", kTH1F, {{50, -0.5f, 49.5f}});
    histos.add("hDuplicateCountType7", "hDuplicateCountType7", kTH1F, {{50, -0.5f, 49.5f}});
    histos.add("hDuplicateCountType7allTPConly", "hDuplicateCountType7allTPConly", kTH1F, {{50, -0.5f, 49.5f}});

    histos.add("hPhotonPt", "hPhotonPt", kTH1F, {{200, 0.0f, 20.0f}});
    histos.add("hPhotonPt_Duplicates", "hPhotonPt_Duplicates", kTH1F, {{200, 0.0f, 20.0f}});
    histos.add("hPhotonPt_withRecoedMcCollision", "hPhotonPt_withRecoedMcCollision", kTH1F, {{200, 0.0f, 20.0f}});
    histos.add("hPhotonPt_withCorrectCollisionCopy", "hPhotonPt_withCorrectCollisionCopy", kTH1F, {{200, 0.0f, 20.0f}});

    histos.add("hPA_All", "hPA_All", kTH1F, {{100, 0.0f, 1.0f}});
    histos.add("hPA_Correct", "hPA_Correct", kTH1F, {{100, 0.0f, 1.0f}});

    // 2D for <de-duplication criteria> vs pT
    histos.add("hPAvsPt_All", "hPAvsPt_All", kTH2F, {{200, 0.0f, 20.0f}, {100, 0.0f, 1.0f}});
    histos.add("hPAvsPt_Correct", "hPAvsPt_Correct", kTH2F, {{200, 0.0f, 20.0f}, {100, 0.0f, 1.0f}});
    histos.add("hDCADaughtersvsPt_All", "hDCADaughtersvsPt_All", kTH2F, {{200, 0.0f, 20.0f}, {100, 0.0f, 5.0f}});
    histos.add("hDCADaughtersvsPt_Correct", "hDCADaughtersvsPt_Correct", kTH2F, {{200, 0.0f, 20.0f}, {100, 0.0f, 5.0f}});
    histos.add("hDCADaughters3DvsPt_All", "hDCADaughters3DvsPt_All", kTH2F, {{200, 0.0f, 20.0f}, {100, 0.0f, 5.0f}});
    histos.add("hDCADaughters3DvsPt_Correct", "hDCADaughters3DvsPt_Correct", kTH2F, {{200, 0.0f, 20.0f}, {100, 0.0f, 5.0f}});
    histos.add("hDCADaughtersXYvsPt_All", "hDCADaughtersXYvsPt_All", kTH2F, {{200, 0.0f, 20.0f}, {100, 0.0f, 5.0f}});
    histos.add("hDCADaughtersXYvsPt_Correct", "hDCADaughtersXYvsPt_Correct", kTH2F, {{200, 0.0f, 20.0f}, {100, 0.0f, 5.0f}});
    histos.add("hDCADaughtersZvsPt_All", "hDCADaughtersZvsPt_All", kTH2F, {{200, 0.0f, 20.0f}, {100, 0.0f, 5.0f}});
    histos.add("hDCADaughtersZvsPt_Correct", "hDCADaughtersZvsPt_Correct", kTH2F, {{200, 0.0f, 20.0f}, {100, 0.0f, 5.0f}});

    // winner-takes-all criteria spectra
    histos.add("hCorrect_BestPA", "hCorrect_BestPA", kTH1F, {{200, 0.0f, 20.0f}});
    histos.add("hCorrect_BestDCADau", "hCorrect_DCADau", kTH1F, {{200, 0.0f, 20.0f}});
    histos.add("hCorrect_BestDCADau3D", "hCorrect_DCADau3D", kTH1F, {{200, 0.0f, 20.0f}});
    histos.add("hCorrect_BestDCADauXY", "hCorrect_DCADauXY", kTH1F, {{200, 0.0f, 20.0f}});
    histos.add("hCorrect_BestDCADauZ", "hCorrect_DCADauZ", kTH1F, {{200, 0.0f, 20.0f}});
    histos.add("hCorrect_BestPAandDCADau3D", "hCorrect_BestPAandDCADau3D", kTH1F, {{200, 0.0f, 20.0f}});

    ccdb->setURL(ccdbConfigurations.ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    // set V0 parameters in the helper
    straHelper.v0selections.minCrossedRows = v0BuilderOpts.minCrossedRows;
    straHelper.v0selections.dcanegtopv = v0BuilderOpts.dcanegtopv;
    straHelper.v0selections.dcapostopv = v0BuilderOpts.dcapostopv;
    straHelper.v0selections.v0cospa = v0BuilderOpts.v0cospa;
    straHelper.v0selections.dcav0dau = v0BuilderOpts.dcav0dau;
    straHelper.v0selections.v0radius = v0BuilderOpts.v0radius;
    straHelper.v0selections.maxDaughterEta = v0BuilderOpts.maxDaughterEta;

    // set cascade parameters in the helper
    straHelper.cascadeselections.minCrossedRows = cascadeBuilderOpts.minCrossedRows;
    straHelper.cascadeselections.dcabachtopv = cascadeBuilderOpts.dcabachtopv;
    straHelper.cascadeselections.cascradius = cascadeBuilderOpts.cascradius;
    straHelper.cascadeselections.casccospa = cascadeBuilderOpts.casccospa;
    straHelper.cascadeselections.dcacascdau = cascadeBuilderOpts.dcacascdau;
    straHelper.cascadeselections.lambdaMassWindow = cascadeBuilderOpts.lambdaMassWindow;
    straHelper.cascadeselections.maxDaughterEta = cascadeBuilderOpts.maxDaughterEta;
  }

  template <typename TCollisions>
  bool initCCDB(aod::BCsWithTimestamps const& bcs, TCollisions const& collisions)
  {
    auto bc = collisions.size() ? collisions.begin().template bc_as<aod::BCsWithTimestamps>() : bcs.begin();
    if (!bcs.size()) {
      LOGF(warn, "No BC found, skipping this DF.");
      return false; // signal to skip this DF
    }

    if (mRunNumber == bc.runNumber()) {
      return true;
    }

    auto timestamp = bc.timestamp();
    o2::parameters::GRPMagField* grpmag = 0x0;

    grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbConfigurations.grpmagPath, timestamp);
    if (!grpmag) {
      LOG(fatal) << "Got nullptr from CCDB for path " << ccdbConfigurations.grpmagPath << " of object GRPMagField for timestamp " << timestamp;
    }
    o2::base::Propagator::initFieldFromGRP(grpmag);

    // Fetch magnetic field from ccdb for current collision
    auto magneticField = o2::base::Propagator::Instance()->getNominalBz();
    LOG(info) << "Retrieved GRP for timestamp " << timestamp << " with magnetic field of " << magneticField << " kG";

    // Set magnetic field value once known
    straHelper.fitter.setBz(magneticField);

    // acquire LUT for this timestamp
    LOG(info) << "Loading material look-up table for timestamp: " << timestamp;
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->getForTimeStamp<o2::base::MatLayerCylSet>(ccdbConfigurations.lutPath, timestamp));
    o2::base::Propagator::Instance()->setMatLUT(lut);
    straHelper.lut = lut;

    LOG(info) << "Fully configured for run: " << bc.runNumber();
    // mmark this run as configured
    mRunNumber = bc.runNumber();

    // initialize only if needed, avoid unnecessary CCDB calls
    mVDriftMgr.init(&ccdb->instance());
    mVDriftMgr.update(timestamp);

    return true;
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels> const& collisions, aod::McCollisions const& mcCollisions, aod::V0s const& V0s, LabeledTracksExtra const& tracks, aod::McParticles const& mcParticles, aod::BCsWithTimestamps const& bcs)
  {
    if (!initCCDB(bcs, collisions))
      return;

    std::vector<o2::pwglf::V0group> v0tableGrouped = o2::pwglf::groupDuplicates(V0s);

    // determine map of McCollisions -> Collisions
    std::vector<std::vector<int>> mcCollToColl(mcCollisions.size());
    for (auto const& collision : collisions) {
      if (collision.mcCollisionId() > -1) {
        // useful to determine if collision has been reconstructed afterwards
        mcCollToColl[collision.mcCollisionId()].push_back(collision.globalIndex());
      }
    }

    // simple inspection of grouped duplicates
    for (size_t iV0 = 0; iV0 < v0tableGrouped.size(); iV0++) {
      // base QA histograms
      histos.fill(HIST("hDuplicateCount"), v0tableGrouped[iV0].collisionIds.size());
      if (v0tableGrouped[iV0].v0Type == 7) {
        histos.fill(HIST("hDuplicateCountType7"), v0tableGrouped[iV0].collisionIds.size());
      }

      // Monte Carlo exclusive: process
      auto pTrack = tracks.rawIteratorAt(v0tableGrouped[iV0].posTrackId);
      auto nTrack = tracks.rawIteratorAt(v0tableGrouped[iV0].negTrackId);
      bool pTrackTPCOnly = (pTrack.hasTPC() && !pTrack.hasITS() && !pTrack.hasTRD() && !pTrack.hasTOF());
      bool nTrackTPCOnly = (nTrack.hasTPC() && !nTrack.hasITS() && !nTrack.hasTRD() && !nTrack.hasTOF());

      if (v0tableGrouped[iV0].v0Type == 7 && pTrackTPCOnly && nTrackTPCOnly) {
        histos.fill(HIST("hDuplicateCountType7allTPConly"), v0tableGrouped[iV0].collisionIds.size());
      }

      int pTrackLabel = pTrack.mcParticleId();
      int nTrackLabel = nTrack.mcParticleId();
      int v0Label = findMotherFromLabels(pTrackLabel, nTrackLabel, -11, 11, 22, mcParticles);
      int correctMcCollision = -1;
      if (v0Label > -1) {
        // this mc particle exists and is a gamma
        auto mcV0 = mcParticles.rawIteratorAt(v0Label);
        correctMcCollision = mcV0.mcCollisionId();

        histos.fill(HIST("hPhotonPt"), mcV0.pt());

        if (mcCollToColl[mcV0.mcCollisionId()].size() > 0) {
          histos.fill(HIST("hPhotonPt_withRecoedMcCollision"), mcV0.pt());
        }

        bool hasCorrectCollisionCopy = false;
        for (size_t ic = 0; ic < v0tableGrouped[iV0].collisionIds.size(); ic++) {
          for (size_t imcc = 0; imcc < mcCollToColl[mcV0.mcCollisionId()].size(); imcc++) {
            if (v0tableGrouped[iV0].collisionIds[ic] == mcCollToColl[mcV0.mcCollisionId()][imcc]) {
              hasCorrectCollisionCopy = true;
            }
          }
        }

        if (hasCorrectCollisionCopy) {
          histos.fill(HIST("hPhotonPt_withCorrectCollisionCopy"), mcV0.pt());
        }

        std::vector<o2::pwglf::v0candidate> v0duplicates; // Vector of v0 candidate duplicates
        std::vector<bool> v0duplicatesCorrectlyAssociated;

        // de-duplication strategy tests start here
        // store best-of index for cross-checking strict de-duplication techniques

        float bestPointingAngle = .99;
        float bestDCADaughters = 1e+6;
        float bestDCADaughters3D = 1e+6;
        float bestDCADaughtersXY = 1e+6;
        float bestDCADaughtersZ = 1e+6;

        bool bestPointingAngleCorrect = false;
        bool bestDCADaughtersCorrect = false;
        bool bestDCADaughters3DCorrect = false;
        bool bestDCADaughtersXYCorrect = false;
        bool bestDCADaughtersZCorrect = false;

        // START OF MAIN DUPLICATE LOOP IS HERE
        for (size_t ic = 0; ic < v0tableGrouped[iV0].collisionIds.size(); ic++) {
          // simple duplicate accounting
          histos.fill(HIST("hPhotonPt_Duplicates"), mcV0.pt());

          // check if candidate is correctly associated
          bool correctlyAssociated = false;
          for (size_t imcc = 0; imcc < mcCollToColl[correctMcCollision].size(); imcc++) {
            if (v0tableGrouped[iV0].collisionIds[ic] == mcCollToColl[correctMcCollision][imcc]) {
              correctlyAssociated = true;
            }
          }
          // store check for correct association
          v0duplicatesCorrectlyAssociated.push_back(correctlyAssociated);

          // actually treat tracks
          auto posTrackPar = getTrackParCov(pTrack);
          auto negTrackPar = getTrackParCov(nTrack);

          auto const& collision = collisions.rawIteratorAt(v0tableGrouped[iV0].collisionIds[ic]);

          // handle TPC-only tracks properly (photon conversions)
          if (v0BuilderOpts.moveTPCOnlyTracks) {
            bool isPosTPCOnly = (pTrack.hasTPC() && !pTrack.hasITS() && !pTrack.hasTRD() && !pTrack.hasTOF());
            if (isPosTPCOnly) {
              // Nota bene: positive is TPC-only -> this entire V0 merits treatment as photon candidate
              posTrackPar.setPID(o2::track::PID::Electron);
              negTrackPar.setPID(o2::track::PID::Electron);

              if (!mVDriftMgr.moveTPCTrack<aod::BCsWithTimestamps, soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>>(collision, pTrack, posTrackPar)) {
                return;
              }
            }

            bool isNegTPCOnly = (nTrack.hasTPC() && !nTrack.hasITS() && !nTrack.hasTRD() && !nTrack.hasTOF());
            if (isNegTPCOnly) {
              // Nota bene: negative is TPC-only -> this entire V0 merits treatment as photon candidate
              posTrackPar.setPID(o2::track::PID::Electron);
              negTrackPar.setPID(o2::track::PID::Electron);

              if (!mVDriftMgr.moveTPCTrack<aod::BCsWithTimestamps, soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>>(collision, nTrack, negTrackPar)) {
                return;
              }
            }
          } // end TPC drift treatment

          // process candidate with helper
          bool buildOK = straHelper.buildV0Candidate(v0tableGrouped[iV0].collisionIds[ic], collision.posX(), collision.posY(), collision.posZ(), pTrack, nTrack, posTrackPar, negTrackPar, true, false);

          v0duplicates.push_back(straHelper.v0);

          float daughterDCA3D = std::hypot(
            straHelper.v0.positivePosition[0] - straHelper.v0.negativePosition[0],
            straHelper.v0.positivePosition[1] - straHelper.v0.negativePosition[1],
            straHelper.v0.positivePosition[2] - straHelper.v0.negativePosition[2]);
          float daughterDCAXY = std::hypot(
            straHelper.v0.positivePosition[0] - straHelper.v0.negativePosition[0],
            straHelper.v0.positivePosition[1] - straHelper.v0.negativePosition[1]);
          float daughterDCAZ = std::abs(
            straHelper.v0.positivePosition[2] - straHelper.v0.negativePosition[2]);

          if (!buildOK) {
            daughterDCA3D = daughterDCAXY = daughterDCAZ = 1e+6;
          }

          histos.fill(HIST("hPA_All"), straHelper.v0.pointingAngle);
          histos.fill(HIST("hPAvsPt_All"), mcV0.pt(), straHelper.v0.pointingAngle);
          histos.fill(HIST("hDCADaughtersvsPt_All"), mcV0.pt(), straHelper.v0.daughterDCA);
          histos.fill(HIST("hDCADaughters3DvsPt_All"), mcV0.pt(), daughterDCA3D);
          histos.fill(HIST("hDCADaughtersXYvsPt_All"), mcV0.pt(), daughterDCAXY);
          histos.fill(HIST("hDCADaughtersZvsPt_All"), mcV0.pt(), daughterDCAZ);

          if (correctlyAssociated) {
            histos.fill(HIST("hPA_Correct"), straHelper.v0.pointingAngle);
            histos.fill(HIST("hPAvsPt_Correct"), mcV0.pt(), straHelper.v0.pointingAngle);
            histos.fill(HIST("hDCADaughtersvsPt_Correct"), mcV0.pt(), straHelper.v0.daughterDCA);
            histos.fill(HIST("hDCADaughters3DvsPt_Correct"), mcV0.pt(), daughterDCA3D);
            histos.fill(HIST("hDCADaughtersXYvsPt_Correct"), mcV0.pt(), daughterDCAXY);
            histos.fill(HIST("hDCADaughtersZvsPt_Correct"), mcV0.pt(), daughterDCAZ);
          }

          // check criteria
          if (straHelper.v0.pointingAngle < bestPointingAngle) {
            bestPointingAngle = straHelper.v0.pointingAngle;
            bestPointingAngleCorrect = correctlyAssociated;
          }
          if (straHelper.v0.daughterDCA < bestDCADaughters) {
            bestDCADaughters = straHelper.v0.daughterDCA;
            bestDCADaughtersCorrect = correctlyAssociated;
          }
          if (daughterDCA3D < bestDCADaughters3D) {
            bestDCADaughters3D = daughterDCA3D;
            bestDCADaughters3DCorrect = correctlyAssociated;
          }
          if (daughterDCAXY < bestDCADaughtersXY) {
            bestDCADaughtersXY = daughterDCAXY;
            bestDCADaughtersXYCorrect = correctlyAssociated;
          }
          if (daughterDCAZ < bestDCADaughtersZ) {
            bestDCADaughtersZ = daughterDCAZ;
            bestDCADaughtersZCorrect = correctlyAssociated;
          }
        } // end duplicate loop

        if (hasCorrectCollisionCopy) {
          // check individual criteria for winner-is-correct
          if (bestPointingAngleCorrect) {
            histos.fill(HIST("hCorrect_BestPA"), mcV0.pt());
          }
          if (bestDCADaughtersCorrect) {
            histos.fill(HIST("hCorrect_BestDCADau"), mcV0.pt());
          }
          if (bestDCADaughters3DCorrect) {
            histos.fill(HIST("hCorrect_BestDCADau3D"), mcV0.pt());
          }
          if (bestDCADaughtersXYCorrect) {
            histos.fill(HIST("hCorrect_BestDCADauXY"), mcV0.pt());
          }
          if (bestDCADaughtersZCorrect) {
            histos.fill(HIST("hCorrect_BestDCADauZ"), mcV0.pt());
          }
          if (bestPointingAngleCorrect && bestDCADaughtersZCorrect) {
            histos.fill(HIST("hCorrect_BestPAandDCADau3D"), mcV0.pt());
          }
        }

        // printout for inspection
        // TString cosPAString = "";
        // for (size_t iCollisionId = 0; iCollisionId < v0tableGrouped[iV0].collisionIds.size(); iCollisionId++) {
        //   cosPAString.Append(Form("%.5f ", v0duplicates[iCollisionId].pointingAngle));
        // }
        // LOGF(info, "#%i (p,n) = (%i,%i), type %i, point. angles: %s", iV0, v0tableGrouped[iV0].posTrackId, v0tableGrouped[iV0].negTrackId, v0tableGrouped[iV0].v0Type, cosPAString.Data());
      } // end this-is-a-mc-gamma check
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<v0assoqa>(cfgc)};
}
