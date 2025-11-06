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
// Test hypertriton task
// =====================
//
// First rough code for hypertriton checks
//
//    Comments, questions, complaints, suggestions?
//    Please write to:
//    david.dobrigkeit.chinellato@cern.ch
//
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
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
using TracksCompleteIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCLfPi, aod::pidTPCLfHe>;
using TracksCompleteIUMC = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCLfPi, aod::pidTPCLfHe, aod::McTrackLabels>;
using V0MC = soa::Join<aod::V0Datas, aod::McV0Labels>;

struct hypertritonAnalysis {
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry registry{
    "registry",
    {
      // Invariant mass
      {"h2dMassHypertriton", "h2dMassHypertriton", {HistType::kTH2F, {{100, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {400, 2.800f, 3.200f, "Inv. Mass (GeV/c^{2})"}}}},
      {"h2dMassAntiHypertriton", "h2dMassAntiHypertriton", {HistType::kTH2F, {{100, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {400, 2.800f, 3.200f, "Inv. Mass (GeV/c^{2})"}}}},
      // Basic QA histogram
      {"h3dPtVsMassHyVsDCAxy", "h3dPtVsMassHyVsDCAxy", {HistType::kTH3F, {{10, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {200, 2.900f, 3.100f, "Inv. Mass (GeV/c^{2})"}, {50, 0.0f, 0.5f, "abs(dcaxy)"}}}},
      {"h3dPtVsMassAHyVsDCAxy", "h3dPtVsMassAHyVsDCAxy", {HistType::kTH3F, {{10, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {200, 2.900f, 3.100f, "Inv. Mass (GeV/c^{2})"}, {50, 0.0f, 0.5f, "abs(dcaxy)"}}}},
      // Very simple QA
      {"h2dHypertritonQAV0Radius", "h2dHypertritonQAV0Radius", {HistType::kTH2F, {{10, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {200, 0, 50}}}},
      {"h2dHypertritonQADCAV0Dau", "h2dHypertritonQADCAV0Dau", {HistType::kTH2F, {{10, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {100, 0, 5}}}},
      {"h2dHypertritonQADCAPosToPV", "h2dHypertritonQADCAPosToPV", {HistType::kTH2F, {{10, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {200, -1, 1}}}},
      {"h2dHypertritonQADCANegToPV", "h2dHypertritonQADCANegToPV", {HistType::kTH2F, {{10, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {200, -10, 10}}}},
      {"h2dHypertritonQADCAToPV", "h2dHypertritonQADCAToPV", {HistType::kTH2F, {{10, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {100, -1, 1}}}},

      {"h2dAntiHypertritonQAV0Radius", "h2dAntiHypertritonQAV0Radius", {HistType::kTH2F, {{10, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {200, 0, 50}}}},
      {"h2dAntiHypertritonQADCAV0Dau", "h2dAntiHypertritonQADCAV0Dau", {HistType::kTH2F, {{10, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {100, 0, 5}}}},
      {"h2dAntiHypertritonQADCAPosToPV", "h2dAntiHypertritonQADCAPosToPV", {HistType::kTH2F, {{10, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {200, -10, 10}}}},
      {"h2dAntiHypertritonQADCANegToPV", "h2dAntiHypertritonQADCANegToPV", {HistType::kTH2F, {{10, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {200, -1, 1}}}},
      {"h2dAntiHypertritonQADCAToPV", "h2dAntiHypertritonQADCAToPV", {HistType::kTH2F, {{10, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {100, -1, 1}}}},

      // Bookkeeping
      {"hEventSelection", "hEventSelection", {HistType::kTH1F, {{3, -0.5f, 2.5f}}}},
      {"V0loopFiltersCounts", "V0loopFiltersCounts", {HistType::kTH1F, {{10, -0.5f, 9.5f}}}},
      // dEdx QA plots
      {"h2dHypertritonQAdEdxHelium", "h2dHypertritonQAdEdxHelium", {HistType::kTH2F, {{50, 0.0f, 10.0f, "p/z (GeV/c)"}, {200, -10.0f, 10.0f, "dE/dx N_{sigma}"}}}},
      {"h2dHypertritonQAdEdxPion", "h2dHypertritonQAdEdxPion", {HistType::kTH2F, {{50, 0.0f, 10.0f, "p/z (GeV/c)"}, {200, -10.0f, 10.0f, "dE/dx N_{sigma}"}}}},
      {"h2dAntiHypertritonQAdEdxHelium", "h2dAntiHypertritonQAdEdxHelium", {HistType::kTH2F, {{50, 0.0f, 10.0f, "p/z (GeV/c)"}, {200, -10.0f, 10.0f, "dE/dx N_{sigma}"}}}},
      {"h2dAntiHypertritonQAdEdxPion", "h2dAntiHypertritonQAdEdxPion", {HistType::kTH2F, {{50, 0.0f, 10.0f, "p/z (GeV/c)"}, {200, -10.0f, 10.0f, "dE/dx N_{sigma}"}}}},
    },
  };

  // Selection criteria
  Configurable<double> v0cospa{"v0cospa", 0.95, "V0 CosPA"}; // very open, please
  Configurable<float> dcav0dau{"dcav0dau", 1.0, "DCA V0 Daughters"};
  Configurable<float> dca3hetopv{"dca3hetopv", -1, "DCA helium To PV"};
  Configurable<float> dcapiontopv{"dcapiontopv", .1, "DCA pion To PV"};
  Configurable<float> dcahyptopv{"dcahyptopv", .2, "DCA Hypertriton To PV"}; // propagated correctly
  Configurable<float> v0radius{"v0radius", 1.0, "v0radius"};
  Configurable<float> rapidity{"rapidity", 0.8, "rapidity"};
  Configurable<float> heliumdEdx{"heliumdEdx", 8, "heliumdEdx"};
  Configurable<float> piondEdx{"piondEdx", 5, "piondEdx"};

  Configurable<bool> event_sel8_selection{"event_sel8_selection", true, "event selection count post sel8 cut"};
  Configurable<bool> event_posZ_selection{"event_posZ_selection", true, "event selection count post poZ cut"};

  // Material correction to use when propagating hypertriton: none
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  o2::track::TrackPar lHyTrack;

  // CCDB options
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  void init(InitContext const&)
  {
    resetHistos();

    registry.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    registry.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(2, "Sel8 cut");
    registry.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(3, "posZ cut");

    mRunNumber = 0;
    d_bz = 0;

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    auto run3grp_timestamp = bc.timestamp();

    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      if (d_bz_input < -990) {
        // Fetch magnetic field from ccdb for current collision
        d_bz = grpo->getNominalL3Field();
        LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
      } else {
        d_bz = d_bz_input;
      }
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      if (d_bz_input < -990) {
        // Fetch magnetic field from ccdb for current collision
        d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
        LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
      } else {
        d_bz = d_bz_input;
      }
    }
    mRunNumber = bc.runNumber();
  }

  int mRunNumber;
  float d_bz;

  enum anastep { kHypAll = 0,
                 kHypRadius,
                 kHypCosPA,
                 kHypDCADaughters,
                 kHypRapidity,
                 kHypTPCdEdx,
                 kHypDauDCAtoPV,
                 kHypDCAtoPV,
                 kHypAllSteps };

  enum evselstep { kEvSelAll = 0,
                   kEvSelBool,
                   kEvSelVtxZ,
                   kEvSelAllSteps };

  // Helper to do bookkeeping and late filling of QA histos
  std::array<long, kHypAllSteps> stats;
  std::array<long, kEvSelAllSteps> evselstats;

  void resetHistos()
  {
    for (Int_t ii = 0; ii < kHypAllSteps; ii++)
      stats[ii] = 0;
    for (Int_t ii = 0; ii < kEvSelAllSteps; ii++)
      evselstats[ii] = 0;
  }

  void fillHistos()
  {
    for (Int_t ii = 0; ii < kHypAllSteps; ii++)
      registry.fill(HIST("V0loopFiltersCounts"), ii, stats[ii]);
    for (Int_t ii = 0; ii < kEvSelAllSteps; ii++)
      registry.fill(HIST("hEventSelection"), ii, evselstats[ii]);
  }

  void processRealData(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::V0Datas const& V0s, TracksCompleteIU const& /*tracks*/, aod::BCsWithTimestamps const&)
  {
    /* check the previous run number */
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    std::array<float, 2> dcaInfo;

    evselstats[kEvSelAll]++;
    if (event_sel8_selection && !collision.sel8()) {
      return;
    }
    evselstats[kEvSelBool]++;
    if (event_posZ_selection && abs(collision.posZ()) > 10.f) { // 10cm
      return;
    }
    evselstats[kEvSelVtxZ]++;

    for (auto& v0 : V0s) {
      stats[kHypAll]++;
      registry.fill(HIST("h2dHypertritonQAV0Radius"), v0.ptHypertriton(), v0.v0radius());
      registry.fill(HIST("h2dHypertritonQADCAV0Dau"), v0.ptHypertriton(), v0.dcaV0daughters());
      registry.fill(HIST("h2dHypertritonQADCAPosToPV"), v0.ptHypertriton(), v0.posTrack_as<TracksCompleteIU>().dcaXY());
      registry.fill(HIST("h2dHypertritonQADCANegToPV"), v0.ptHypertriton(), v0.negTrack_as<TracksCompleteIU>().dcaXY());
      registry.fill(HIST("h2dAntiHypertritonQAV0Radius"), v0.ptAntiHypertriton(), v0.v0radius());
      registry.fill(HIST("h2dAntiHypertritonQADCAV0Dau"), v0.ptAntiHypertriton(), v0.dcaV0daughters());
      registry.fill(HIST("h2dAntiHypertritonQADCAPosToPV"), v0.ptAntiHypertriton(), v0.posTrack_as<TracksCompleteIU>().dcaXY());
      registry.fill(HIST("h2dAntiHypertritonQADCANegToPV"), v0.ptAntiHypertriton(), v0.negTrack_as<TracksCompleteIU>().dcaXY());
      registry.fill(HIST("h2dHypertritonQAdEdxHelium"),
                    TMath::Sqrt(v0.pxpos() * v0.pxpos() + v0.pypos() * v0.pypos() + v0.pzpos() * v0.pzpos()),
                    v0.posTrack_as<TracksCompleteIU>().tpcNSigmaHe());
      registry.fill(HIST("h2dHypertritonQAdEdxPion"),
                    TMath::Sqrt(v0.pxneg() * v0.pxneg() + v0.pyneg() * v0.pyneg() + v0.pzneg() * v0.pzneg()),
                    v0.negTrack_as<TracksCompleteIU>().tpcNSigmaPi());
      registry.fill(HIST("h2dAntiHypertritonQAdEdxHelium"),
                    TMath::Sqrt(v0.pxneg() * v0.pxneg() + v0.pyneg() * v0.pyneg() + v0.pzneg() * v0.pzneg()),
                    v0.negTrack_as<TracksCompleteIU>().tpcNSigmaHe());
      registry.fill(HIST("h2dAntiHypertritonQAdEdxPion"),
                    TMath::Sqrt(v0.pxpos() * v0.pxpos() + v0.pypos() * v0.pypos() + v0.pzpos() * v0.pzpos()),
                    v0.posTrack_as<TracksCompleteIU>().tpcNSigmaPi());
      if (v0.v0radius() > v0radius) {
        stats[kHypRadius]++;
        if (v0.v0cosPA() > v0cospa) {
          stats[kHypCosPA]++;
          if (v0.dcaV0daughters() < dcav0dau) {
            stats[kHypDCADaughters]++;
            //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
            // Hypertriton
            if (TMath::Abs(v0.yHypertriton()) < rapidity) {
              stats[kHypRapidity]++;
              if (v0.posTrack_as<TracksCompleteIU>().tpcNSigmaHe() < heliumdEdx &&
                  TMath::Abs(v0.negTrack_as<TracksCompleteIU>().tpcNSigmaPi()) < piondEdx) {
                stats[kHypTPCdEdx]++;
                if (TMath::Abs(v0.posTrack_as<TracksCompleteIU>().dcaXY()) > dca3hetopv &&
                    TMath::Abs(v0.negTrack_as<TracksCompleteIU>().dcaXY()) > dcapiontopv) {
                  stats[kHypDauDCAtoPV]++;

                  lHyTrack = o2::track::TrackPar(
                    {v0.x(), v0.y(), v0.z()},
                    {2.0f * v0.pxpos() + v0.pxneg(), 2.0f * v0.pypos() + v0.pyneg(), 2.0f * v0.pzpos() + v0.pzneg()},
                    +1, true);

                  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, lHyTrack, 2.f, matCorr, &dcaInfo);
                  registry.fill(HIST("h2dHypertritonQADCAToPV"), v0.ptAntiHypertriton(), dcaInfo[0]);
                  if (TMath::Abs(dcaInfo[0]) < dcahyptopv) {
                    stats[kHypDCAtoPV]++;
                    registry.fill(HIST("h2dMassHypertriton"), v0.ptHypertriton(), v0.mHypertriton());
                    registry.fill(HIST("h3dPtVsMassHyVsDCAxy"), v0.ptHypertriton(), v0.mHypertriton(), TMath::Abs(dcaInfo[0]));
                  }
                }
              }
            }
            //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
            // AntiHypertriton
            if (TMath::Abs(v0.yAntiHypertriton()) < rapidity) {
              if (v0.negTrack_as<TracksCompleteIU>().tpcNSigmaHe() < heliumdEdx &&
                  TMath::Abs(v0.posTrack_as<TracksCompleteIU>().tpcNSigmaPi()) < piondEdx) {
                stats[kHypTPCdEdx]++;
                if (TMath::Abs(v0.posTrack_as<TracksCompleteIU>().dcaXY()) > dcapiontopv &&
                    TMath::Abs(v0.negTrack_as<TracksCompleteIU>().dcaXY()) > dca3hetopv) {
                  stats[kHypDauDCAtoPV]++;

                  lHyTrack = o2::track::TrackPar(
                    {v0.x(), v0.y(), v0.z()},
                    {v0.pxpos() + 2.0f * v0.pxneg(), v0.pypos() + 2.0f * v0.pyneg(), v0.pzpos() + 2.0f * v0.pzneg()},
                    -1, true);

                  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, lHyTrack, 2.f, matCorr, &dcaInfo);
                  registry.fill(HIST("h2dAntiHypertritonQADCAToPV"), v0.ptAntiHypertriton(), dcaInfo[0]);
                  if (TMath::Abs(dcaInfo[0]) < dcahyptopv) {
                    stats[kHypDCAtoPV]++;
                    registry.fill(HIST("h2dMassAntiHypertriton"), v0.ptAntiHypertriton(), v0.mAntiHypertriton());
                    registry.fill(HIST("h3dPtVsMassAHyVsDCAxy"), v0.ptAntiHypertriton(), v0.mAntiHypertriton(), TMath::Abs(dcaInfo[0]));
                  }
                }
              }
            }
          }
        }
      }
    } // end v0 loop
    fillHistos();
    resetHistos();
  }
  PROCESS_SWITCH(hypertritonAnalysis, processRealData, "Regular analysis", false);

  void processMC(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, V0MC const& fullV0s, TracksCompleteIUMC const& /*tracks*/, aod::BCsWithTimestamps const&, aod::McParticles const&)
  {
    /* check the previous run number */
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    std::array<float, 2> dcaInfo;

    evselstats[kEvSelAll]++;
    if (event_sel8_selection && !collision.sel8()) {
      return;
    }
    evselstats[kEvSelBool]++;
    if (event_posZ_selection && abs(collision.posZ()) > 10.f) { // 10cm
      return;
    }
    evselstats[kEvSelVtxZ]++;

    for (auto& v0 : fullV0s) {
      // Get perfect MC particle
      auto posPartTrack = v0.posTrack_as<TracksCompleteIUMC>();
      auto negPartTrack = v0.negTrack_as<TracksCompleteIUMC>();
      if (!v0.has_mcParticle() || !posPartTrack.has_mcParticle() || !negPartTrack.has_mcParticle())
        continue;
      auto v0mc = v0.mcParticle();
      auto posMC = posPartTrack.mcParticle();
      auto negMC = negPartTrack.mcParticle();

      if (TMath::Abs(v0mc.pdgCode()) != 1010010030)
        continue;
      if (posMC.pdgCode() != 1000020030)
        continue;
      if (negMC.pdgCode() != -211)
        continue;
      stats[kHypAll]++;
      if (v0.v0radius() > v0radius) {
        registry.fill(HIST("h2dHypertritonQAV0Radius"), v0.ptHypertriton(), v0.v0radius());
        registry.fill(HIST("h2dHypertritonQADCAV0Dau"), v0.ptHypertriton(), v0.dcaV0daughters());
        registry.fill(HIST("h2dHypertritonQADCAPosToPV"), v0.ptHypertriton(), v0.posTrack_as<TracksCompleteIUMC>().dcaXY());
        registry.fill(HIST("h2dHypertritonQADCANegToPV"), v0.ptHypertriton(), v0.negTrack_as<TracksCompleteIUMC>().dcaXY());
        registry.fill(HIST("h2dAntiHypertritonQAV0Radius"), v0.ptAntiHypertriton(), v0.v0radius());
        registry.fill(HIST("h2dAntiHypertritonQADCAV0Dau"), v0.ptAntiHypertriton(), v0.dcaV0daughters());
        registry.fill(HIST("h2dAntiHypertritonQADCAPosToPV"), v0.ptAntiHypertriton(), v0.posTrack_as<TracksCompleteIUMC>().dcaXY());
        registry.fill(HIST("h2dAntiHypertritonQADCANegToPV"), v0.ptAntiHypertriton(), v0.negTrack_as<TracksCompleteIUMC>().dcaXY());
        registry.fill(HIST("h2dHypertritonQAdEdxHelium"),
                      TMath::Sqrt(v0.pxpos() * v0.pxpos() + v0.pypos() * v0.pypos() + v0.pzpos() * v0.pzpos()),
                      v0.posTrack_as<TracksCompleteIUMC>().tpcNSigmaHe());
        registry.fill(HIST("h2dHypertritonQAdEdxPion"),
                      TMath::Sqrt(v0.pxneg() * v0.pxneg() + v0.pyneg() * v0.pyneg() + v0.pzneg() * v0.pzneg()),
                      v0.negTrack_as<TracksCompleteIUMC>().tpcNSigmaPi());
        registry.fill(HIST("h2dAntiHypertritonQAdEdxHelium"),
                      TMath::Sqrt(v0.pxneg() * v0.pxneg() + v0.pyneg() * v0.pyneg() + v0.pzneg() * v0.pzneg()),
                      v0.negTrack_as<TracksCompleteIUMC>().tpcNSigmaHe());
        registry.fill(HIST("h2dAntiHypertritonQAdEdxPion"),
                      TMath::Sqrt(v0.pxpos() * v0.pxpos() + v0.pypos() * v0.pypos() + v0.pzpos() * v0.pzpos()),
                      v0.posTrack_as<TracksCompleteIUMC>().tpcNSigmaPi());
        stats[kHypRadius]++;
        if (v0.v0cosPA() > v0cospa) {
          stats[kHypCosPA]++;
          if (v0.dcaV0daughters() < dcav0dau) {
            stats[kHypDCADaughters]++;
            //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
            // Hypertriton
            if (TMath::Abs(v0.yHypertriton()) < rapidity) {
              stats[kHypRapidity]++;
              if (v0.posTrack_as<TracksCompleteIUMC>().tpcNSigmaHe() < heliumdEdx &&
                  TMath::Abs(v0.negTrack_as<TracksCompleteIUMC>().tpcNSigmaPi()) < piondEdx) {
                stats[kHypTPCdEdx]++;
                if (TMath::Abs(v0.posTrack_as<TracksCompleteIUMC>().dcaXY()) > dca3hetopv &&
                    TMath::Abs(v0.negTrack_as<TracksCompleteIUMC>().dcaXY()) > dcapiontopv) {
                  stats[kHypDauDCAtoPV]++;

                  lHyTrack = o2::track::TrackPar(
                    {v0.x(), v0.y(), v0.z()},
                    {2.0f * v0.pxpos() + v0.pxneg(), 2.0f * v0.pypos() + v0.pyneg(), 2.0f * v0.pzpos() + v0.pzneg()},
                    +1, true);

                  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, lHyTrack, 2.f, matCorr, &dcaInfo);
                  registry.fill(HIST("h2dHypertritonQADCAToPV"), v0.ptHypertriton(), dcaInfo[0]);
                  if (TMath::Abs(dcaInfo[0]) < dcahyptopv) {
                    stats[kHypDCAtoPV]++;
                    registry.fill(HIST("h2dMassHypertriton"), v0.ptHypertriton(), v0.mHypertriton());
                    registry.fill(HIST("h3dPtVsMassHyVsDCAxy"), v0.ptHypertriton(), v0.mHypertriton(), TMath::Abs(dcaInfo[0]));
                  }
                }
              }
            }
            //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
            // AntiHypertriton
            if (TMath::Abs(v0.yAntiHypertriton()) < rapidity) {
              if (v0.negTrack_as<TracksCompleteIUMC>().tpcNSigmaHe() < heliumdEdx &&
                  TMath::Abs(v0.posTrack_as<TracksCompleteIUMC>().tpcNSigmaPi()) < piondEdx) {
                stats[kHypTPCdEdx]++;
                if (TMath::Abs(v0.posTrack_as<TracksCompleteIUMC>().dcaXY()) > dcapiontopv &&
                    TMath::Abs(v0.negTrack_as<TracksCompleteIUMC>().dcaXY()) > dca3hetopv) {
                  stats[kHypDauDCAtoPV]++;

                  lHyTrack = o2::track::TrackPar(
                    {v0.x(), v0.y(), v0.z()},
                    {v0.pxpos() + 2.0f * v0.pxneg(), v0.pypos() + 2.0f * v0.pyneg(), v0.pzpos() + 2.0f * v0.pzneg()},
                    -1, true);

                  o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, lHyTrack, 2.f, matCorr, &dcaInfo);
                  registry.fill(HIST("h2dAntiHypertritonQADCAToPV"), v0.ptAntiHypertriton(), dcaInfo[0]);
                  if (TMath::Abs(dcaInfo[0]) < dcahyptopv) {
                    stats[kHypDCAtoPV]++;
                    registry.fill(HIST("h2dMassAntiHypertriton"), v0.ptAntiHypertriton(), v0.mAntiHypertriton());
                    registry.fill(HIST("h3dPtVsMassAHyVsDCAxy"), v0.ptAntiHypertriton(), v0.mAntiHypertriton(), TMath::Abs(dcaInfo[0]));
                  }
                }
              }
            }
          }
        }
      }
    } // end v0 loop
    fillHistos();
    resetHistos();
  }
  PROCESS_SWITCH(hypertritonAnalysis, processMC, "MC analysis", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hypertritonAnalysis>(cfgc)};
}
