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
/// \brief create a table applying some basic cuts on the ITS and DCA.
/// \author Sofia Tomassini, Gleb Romanenko, Nicolò Jacazio
/// \since 31 May 2023

#include "PWGCF/Femto3D/DataModel/singletrackselector.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include <Framework/AnalysisDataModel.h>

#include <fairlogger/Logger.h>

#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;
using namespace o2::aod;
//::singletrackselector; // the namespace defined in .h

struct singleTrackSelector {
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<bool> applySkimming{"applySkimming", false, "Skimmed dataset processing"};
  Configurable<std::string> cfgSkimming{"cfgSkimming", "fPD", "Configurable for skimming"};
  Configurable<bool> CBThadronPID{"CBThadronPID", false, "Apply ev. sel. based on RCT flag `hadronPID`"}; // more in Common/CCDB/RCTSelectionFlags.h

  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply rapidity cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  // Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  Configurable<int> centTableToUse{"centTableToUse", 1, "Flag to choose cent./mult.perc. estimator (Run3 only [FTOC for PbPb; FTOM for pp], for Run2 the V0M is used): 0 -> CentFV0As, 1 -> CentFT0Ms, 2 -> CentFT0As, 3 -> CentFT0Cs, 4 -> CentFDDMs, 5 -> CentNTPVs"};
  Configurable<int> multTableToUse{"multTableToUse", 1, "Flag to choose mult. estimator (Run3 only): 0 -> TPCMults, 1 -> MultNTracksPV, 2 -> MultNTracksPVeta1"};
  Configurable<bool> rejectNotPropagatedTrks{"rejectNotPropagatedTrks", true, "rejects tracks that are not propagated to the primary vertex"};
  Configurable<bool> enable_gen_info{"enable_gen_info", false, "Enable MC true info"};
  Configurable<bool> fetchRate{"fetchRate", true, "Fetch the hadronic rate from the CCDB"};

  Configurable<std::vector<int>> _particlesToKeep{"particlesToKeepPDGs", std::vector<int>{2212, 1000010020}, "PDG codes of perticles for which the 'singletrackselector' tables will be created (only proton and deurton are supported now)"};
  Configurable<std::vector<float>> keepWithinNsigmaTPC{"keepWithinNsigmaTPC", std::vector<float>{-4.0f, 4.0f}, "TPC range for preselection of particles specified with PDG"};
  Configurable<std::vector<int>> _particlesToReject{"particlesToRejectPDGs", std::vector<int>{211, 321, 1000020030}, "PDG codes of particles that will be rejected with TOF (only pion, kaon, proton and deurton are supported now)"};
  Configurable<std::vector<float>> rejectWithinNsigmaTOF{"rejectWithinNsigmaTOF", std::vector<float>{-5.0f, 5.0f}, "TOF rejection Nsigma range for particles specified with PDG to be rejected"};

  Configurable<float> _pRemoveTofOutOfRange{"pRemoveTofOutOfRange", 100.f, "momentum starting from which request TOF nSigma to be within the stored range (-10 < Nsigma < 10)"};
  Configurable<std::array<float, 3>> _ptRemoveTofOutOfRange{"ptRemoveTofOutOfRange", {100.f, -10.f, 10.f}, "transverse momentum starting from which request TOF nSigma to be within the stored range (-10 < Nsigma < 10)"};

  Configurable<float> _min_P{"min_P", 0.f, "lower mometum limit"};
  Configurable<float> _max_P{"max_P", 100.f, "upper mometum limit"};
  Configurable<float> _min_Pt{"min_Pt", 0.f, "lower trasnverse mometum limit"};
  Configurable<float> _max_Pt{"max_Pt", 100.f, "upper trasnverse mometum limit"};
  Configurable<float> _eta{"eta", 100.f, "abs eta value limit"};
  Configurable<float> _dcaXY{"dcaXY", 1000.f, "Maximum dca of track in xy"};
  Configurable<float> _dcaZ{"dcaZ", 1000.f, "Maximum dca of track in xy"};
  Configurable<float> _dcaXYmin{"dcaXYmin", -1000.f, "Minimum dca of track in xy"};
  Configurable<float> _dcaZmin{"dcaZmin", -1000.f, "Minimum dca of track in xy"};
  Configurable<float> _maxTofChi2{"maxTofChi2", 10.f, "Maximum TOF Chi2 value -> to remove mismatched tracks"};
  Configurable<float> _vertexZ{"VertexZ", 15.0, "abs vertexZ value limit"};
  Configurable<std::pair<float, float>> _centCut{"centCut", std::pair<float, float>{0.f, 100.f}, "[min., max.] centrality range to keep events within"};

  using Trks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidEvTimeFlags, aod::TracksDCA,
                         aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa,
                         aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullTr, aod::pidTPCFullHe,
                         aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa,
                         aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe,
                         aod::TrackSelection, aod::pidTOFbeta>;

  using CollRun2 = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::CentRun2V0Ms>;
  using CollRun3 = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFDDMs, aod::CentNTPVs>;
  using CollRun3MC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::Mults, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFDDMs, aod::CentNTPVs>;

  Produces<o2::aod::SingleCollSels> tableRowColl;
  Produces<o2::aod::SingleCollExtras> tableRowCollExtra;
  Produces<o2::aod::SingleTrackSels> tableRow;
  Produces<o2::aod::SingleTrkExtras> tableRowExtra;

  Produces<o2::aod::SinglePIDEls> tableRowPIDEl;
  Produces<o2::aod::SinglePIDPis> tableRowPIDPi;
  Produces<o2::aod::SinglePIDKas> tableRowPIDKa;
  Produces<o2::aod::SinglePIDPrs> tableRowPIDPr;
  Produces<o2::aod::SinglePIDDes> tableRowPIDDe;
  Produces<o2::aod::SinglePIDTrs> tableRowPIDTr;
  Produces<o2::aod::SinglePIDHes> tableRowPIDHe;

  Produces<o2::aod::SingleTrkMCs> tableRowMC;
  // Produces<o2::aod::SingleTrkMCExtras> tableRowMCExtra;

  Filter eventFilter = (applyEvSel.node() == 0) ||
                       ((applyEvSel.node() == 1) && (aod::evsel::sel7 == true)) ||
                       ((applyEvSel.node() == 2) && (aod::evsel::sel8 == true));
  Filter vertexFilter = nabs(o2::aod::collision::posZ) < _vertexZ;
  Filter trackFilter = ((o2::aod::track::itsChi2NCl <= 36.f) && (o2::aod::track::itsChi2NCl >= 0.f) && (o2::aod::track::tpcChi2NCl >= 0.f) && (o2::aod::track::tpcChi2NCl <= 4.f));

  Filter pFilter = o2::aod::track::p > _min_P&& o2::aod::track::p < _max_P;
  Filter ptFilter = o2::aod::track::pt > _min_Pt&& o2::aod::track::pt < _max_Pt;
  Filter etaFilter = nabs(o2::aod::track::eta) < _eta;
  Filter dcaFilter = ((nabs(o2::aod::track::dcaXY) <= _dcaXY) && (nabs(o2::aod::track::dcaZ) <= _dcaZ)) &&
                     ((nabs(o2::aod::track::dcaXY) >= _dcaXYmin) && (nabs(o2::aod::track::dcaZ) >= _dcaZmin));
  Filter tofChi2Filter = o2::aod::track::tofChi2 < _maxTofChi2;

  ctpRateFetcher mRateFetcher; // inspired by zdcSP.cxx in PWGLF
  int mRunNumber = 0;
  float d_bz = 0.f;

  std::vector<int> particlesToKeep;
  std::vector<int> particlesToReject;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  SliceCache cache;

  rctsel::RCTFlagsChecker myChecker{"CBT_hadronPID"};

  void init(InitContext&)
  {

    particlesToKeep = _particlesToKeep;
    particlesToReject = _particlesToReject;

    if (applySkimming) {
      zorroSummary.setObject(zorro.getZorroSummary());
    }
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    myChecker.init("CBT_hadronPID", true);

    registry.add("hNEvents", "hNEvents", {HistType::kTH1D, {{2, 0.f, 2.f}}});
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(1, "All");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(2, "Skimmed");

    registry.add("hNTracks", "hNTracks", {HistType::kTH1D, {{2, 0.f, 2.f}}});
    registry.get<TH1>(HIST("hNTracks"))->GetXaxis()->SetBinLabel(1, "All");
    registry.get<TH1>(HIST("hNTracks"))->GetXaxis()->SetBinLabel(2, "Selected");

    if (enable_gen_info) {
      registry.add("hNEvents_MCGen", "hNEvents_MCGen", {HistType::kTH1F, {{1, 0.f, 1.f}}});
      registry.add("hGen_EtaPhiPt_Proton", "Gen (anti)protons in true collisions", {HistType::kTH3F, {{100, -1., 1., "#eta"}, {157, 0., o2::constants::math::TwoPI, "#phi"}, {100, -5.f, 5.f, "p_{T} GeV/c"}}});
      registry.add("hGen_EtaPhiPt_Deuteron", "Gen (anti)deuteron in true collisions", {HistType::kTH3F, {{100, -1., 1., "#eta"}, {157, 0., o2::constants::math::TwoPI, "#phi"}, {100, -5.f, 5.f, "p_{T} GeV/c"}}});
      registry.add("hGen_EtaPhiPt_Helium3", "Gen (anti)Helium3 in true collisions", {HistType::kTH3F, {{100, -1., 1., "#eta"}, {157, 0., o2::constants::math::TwoPI, "#phi"}, {100, -5.f, 5.f, "p_{T} GeV/c"}}});
      registry.add("hReco_EtaPhiPt_Proton", "Gen (anti)protons in reco collisions", {HistType::kTH3F, {{100, -1., 1., "#eta"}, {157, 0., o2::constants::math::TwoPI, "#phi"}, {100, -5.f, 5.f, "p_{T} GeV/c"}}});
      registry.add("hReco_EtaPhiPt_Deuteron", "Gen (anti)deuteron in reco collisions", {HistType::kTH3F, {{100, -1., 1., "#eta"}, {157, 0., o2::constants::math::TwoPI, "#phi"}, {100, -5.f, 5.f, "p_{T} GeV/c"}}});
      registry.add("hReco_EtaPhiPt_Helium3", "Gen (anti)Helium3 in reco collisions", {HistType::kTH3F, {{100, -1., 1., "#eta"}, {157, 0., o2::constants::math::TwoPI, "#phi"}, {100, -5.f, 5.f, "p_{T} GeV/c"}}});
    }
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc) // inspired by PWGLF/TableProducer/lambdakzerobuilder.cxx
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    d_bz = 0.f;

    if (applySkimming) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), cfgSkimming.value);
      zorro.populateHistRegistry(registry, bc.runNumber());
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();
    d_bz = 0.1 * d_bz;
  }

  template <bool isMC, typename TracksType>
  inline void fillTrackTables(TracksType const& tracks)
  {
    bool skip_track = false; // flag used for track rejection

    for (auto& track : tracks) {
      registry.fill(HIST("hNTracks"), 0.5);
      if constexpr (isMC) {
        if (!track.has_mcParticle())
          continue;
      }
      if (rejectNotPropagatedTrks && track.trackType() != aod::track::Track) {
        continue;
      }
      skip_track = false;
      for (auto i : particlesToReject) {
        // if satisfied, want to continue in the upper loop (over tracks) -- skip the current track
        // cannot use simple 'continue' since it will be applied to the current loop, so have to use a flag
        if (o2::aod::singletrackselector::TOFselection(track, std::make_pair(i, rejectWithinNsigmaTOF))) {
          skip_track = true;
          break;
        }
      }

      if (skip_track)
        continue;
      registry.fill(HIST("hNTracks"), 1.5);

      for (auto ii : particlesToKeep)
        if (o2::aod::singletrackselector::TPCselection<false>(track, std::make_pair(ii, keepWithinNsigmaTPC))) {
          if (track.p() > _pRemoveTofOutOfRange && !o2::aod::singletrackselector::TOFselection(track, std::make_pair(ii, std::vector<float>{-10.0, 10.0}), std::vector<float>{-10.0, 10.0}))
            continue;
          if (track.pt() > _ptRemoveTofOutOfRange.value[0] && !o2::aod::singletrackselector::TOFselection(track, std::make_pair(ii, std::vector<float>{_ptRemoveTofOutOfRange.value[1], _ptRemoveTofOutOfRange.value[2]}), std::vector<float>{-10.f, +10.f}))
            continue;

          tableRow(tableRowColl.lastIndex(),
                   track.p(),
                   track.eta(),
                   track.phi(),
                   track.sign(),
                   track.tpcNClsFound(),
                   track.tpcNClsShared(),
                   track.itsClusterMap(),
                   track.itsClusterSizes(),
                   singletrackselector::packSymmetric<singletrackselector::binning::dca>(track.dcaXY()),
                   singletrackselector::packSymmetric<singletrackselector::binning::dca>(track.dcaZ()),
                   singletrackselector::packInTable<singletrackselector::binning::chi2>(track.tpcChi2NCl()),
                   singletrackselector::packInTable<singletrackselector::binning::chi2>(track.itsChi2NCl()),
                   singletrackselector::packInTable<singletrackselector::binning::rowsOverFindable>(track.tpcCrossedRowsOverFindableCls()));

          tableRowExtra(track.tpcInnerParam(),
                        track.tpcSignal(),
                        track.beta());

          tableRowPIDEl(singletrackselector::packSymmetric<singletrackselector::binning::nsigma>(track.tofNSigmaEl()),
                        singletrackselector::packSymmetric<singletrackselector::binning::nsigma>(track.tpcNSigmaEl()));

          tableRowPIDPi(singletrackselector::packSymmetric<singletrackselector::binning::nsigma>(track.tofNSigmaPi()),
                        singletrackselector::packSymmetric<singletrackselector::binning::nsigma>(track.tpcNSigmaPi()));

          tableRowPIDKa(singletrackselector::packSymmetric<singletrackselector::binning::nsigma>(track.tofNSigmaKa()),
                        singletrackselector::packSymmetric<singletrackselector::binning::nsigma>(track.tpcNSigmaKa()));

          tableRowPIDPr(singletrackselector::packSymmetric<singletrackselector::binning::nsigma>(track.tofNSigmaPr()),
                        singletrackselector::packSymmetric<singletrackselector::binning::nsigma>(track.tpcNSigmaPr()));

          tableRowPIDDe(singletrackselector::packSymmetric<singletrackselector::binning::nsigma>(track.tofNSigmaDe()),
                        singletrackselector::packSymmetric<singletrackselector::binning::nsigma>(track.tpcNSigmaDe()));

          tableRowPIDTr(singletrackselector::packSymmetric<singletrackselector::binning::nsigma>(track.tofNSigmaTr()),
                        singletrackselector::packSymmetric<singletrackselector::binning::nsigma>(track.tpcNSigmaTr()));

          tableRowPIDHe(singletrackselector::packSymmetric<singletrackselector::binning::nsigma>(track.tofNSigmaHe()),
                        singletrackselector::packSymmetric<singletrackselector::binning::nsigma>(track.tpcNSigmaHe()));

          if constexpr (isMC) {
            int origin = -1;
            if (track.mcParticle().isPhysicalPrimary()) {
              origin = 0; // primary
            } else {
              if (track.mcParticle().getProcess() == 4)
                origin = 1; // weak
              else
                origin = 2; // material
            }

            if (origin == -1)
              LOGF(fatal, "Could not define the origin (primary/weak decay/material) of the track!!!");

            tableRowMC(track.mcParticle().pdgCode(),
                       origin,
                       track.mcParticle().p(),
                       track.mcParticle().eta(),
                       track.mcParticle().phi());

            // tableRowMCExtra(track.mcParticle().vx(),
            //                 track.mcParticle().vy(),
            //                 track.mcParticle().vz());
          }
          break; // break the loop with particlesToKeep after the 'if' condition is satisfied -- don't want double entries
        }
    }
  }

  void processDataRun2(soa::Filtered<CollRun2>::iterator const& collision,
                       soa::Filtered<Trks> const& tracks,
                       aod::BCsWithTimestamps const&)
  {

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    float centValue = collision.centRun2V0M();
    if (centValue >= _centCut.value.first && centValue <= _centCut.value.second) {

      int multValue = -1;

      switch (multTableToUse) {
        case 0:
          multValue = collision.multTPC();
          break;
        case 1:
          multValue = collision.multNTracksPV();
          break;
        case 2:
          multValue = collision.multNTracksPVeta1();
          break;
        default:
          LOGF(fatal, "Invalid flag for mult. estimator has been choosen. Please check.");
          break;
      }

      tableRowColl(multValue,
                   centValue,
                   collision.posZ(),
                   d_bz);

      fillTrackTables<false>(tracks);
    }
  }
  PROCESS_SWITCH(singleTrackSelector, processDataRun2, "process data Run2", false);

  void processDataRun3(soa::Filtered<CollRun3>::iterator const& collision,
                       soa::Filtered<Trks> const& tracks,
                       aod::BCsWithTimestamps const&)
  {

    const auto& bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    if (!myChecker(*collision) && CBThadronPID)
      return;

    registry.fill(HIST("hNEvents"), 0.5);
    if (applySkimming) {
      if (!zorro.isSelected(bc.globalBC())) {
        return;
      }
    }
    registry.fill(HIST("hNEvents"), 1.5);

    double hadronicRate = 0.;
    if (fetchRate) {
      hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), mRunNumber, "ZNC hadronic") * 1.e-3; // fetch IR
    }
    int occupancy = collision.trackOccupancyInTimeRange();

    float centValue = -100.0f;

    switch (centTableToUse) {
      case 0:
        centValue = collision.centFV0A();
        break;
      case 1:
        centValue = collision.centFT0M();
        break;
      case 2:
        centValue = collision.centFT0A();
        break;
      case 3:
        centValue = collision.centFT0C();
        break;
      case 4:
        centValue = collision.centFDDM();
        break;
      case 5:
        centValue = collision.centNTPV();
        break;
      default:
        LOGF(fatal, "Invalid flag for cent./mult.perc. estimator has been choosen. Please check.");
        break;
    }
    if (centValue >= _centCut.value.first && centValue <= _centCut.value.second) {
      int multValue = -1;

      switch (multTableToUse) {
        case 0:
          multValue = collision.multTPC();
          break;
        case 1:
          multValue = collision.multNTracksPV();
          break;
        case 2:
          multValue = collision.multNTracksPVeta1();
          break;
        default:
          LOGF(fatal, "Invalid flag for mult. estimator has been choosen. Please check.");
          break;
      }

      tableRowColl(multValue,
                   centValue,
                   collision.posZ(),
                   d_bz);

      tableRowCollExtra(collision.selection_raw(),
                        hadronicRate,
                        occupancy);

      fillTrackTables<false>(tracks);
    }
  }
  PROCESS_SWITCH(singleTrackSelector, processDataRun3, "process data Run3", true);

  void processMCRun2(soa::Filtered<CollRun2>::iterator const& collision,
                     soa::Filtered<soa::Join<Trks, aod::McTrackLabels>> const& tracks,
                     aod::McParticles const&, aod::BCsWithTimestamps const&)
  {

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    float centValue = collision.centRun2V0M();
    if (centValue >= _centCut.value.first && centValue <= _centCut.value.second) {
      int multValue = -1;

      switch (multTableToUse) {
        case 0:
          multValue = collision.multTPC();
          break;
        case 1:
          multValue = collision.multNTracksPV();
          break;
        case 2:
          multValue = collision.multNTracksPVeta1();
          break;
        default:
          LOGF(fatal, "Invalid flag for mult. estimator has been choosen. Please check.");
          break;
      }

      tableRowColl(multValue,
                   centValue,
                   collision.posZ(),
                   d_bz);

      fillTrackTables<true>(tracks);
    }
  }
  PROCESS_SWITCH(singleTrackSelector, processMCRun2, "process MC Run2", false);

  void processMCRun3(soa::Filtered<CollRun3MC>::iterator const& collision, aod::McCollisions const&,
                     soa::Filtered<soa::Join<Trks, aod::McTrackLabels>> const& tracks,
                     aod::McParticles const& mcParticles,
                     aod::BCsWithTimestamps const&)
  {

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
    double hadronicRate = 0.;
    if (fetchRate) {
      hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), mRunNumber, "ZNC hadronic") * 1.e-3; // fetch IR
    }
    int occupancy = collision.trackOccupancyInTimeRange();

    float centValue = -100.0f;

    switch (centTableToUse) {
      case 0:
        centValue = collision.centFV0A();
        break;
      case 1:
        centValue = collision.centFT0M();
        break;
      case 2:
        centValue = collision.centFT0A();
        break;
      case 3:
        centValue = collision.centFT0C();
        break;
      case 4:
        centValue = collision.centFDDM();
        break;
      case 5:
        centValue = collision.centNTPV();
        break;
      default:
        LOGF(fatal, "Invalid flag for cent./mult.perc. estimator has been choosen. Please check.");
        break;
    }

    if (centValue >= _centCut.value.first && centValue <= _centCut.value.second) {
      int multValue = -1;

      switch (multTableToUse) {
        case 0:
          multValue = collision.multTPC();
          break;
        case 1:
          multValue = collision.multNTracksPV();
          break;
        case 2:
          multValue = collision.multNTracksPVeta1();
          break;
        default:
          LOGF(fatal, "Invalid flag for mult. estimator has been choosen. Please check.");
          break;
      }

      tableRowColl(multValue,
                   centValue,
                   collision.posZ(),
                   d_bz);

      tableRowCollExtra(collision.selection_raw(),
                        hadronicRate,
                        occupancy);

      fillTrackTables<true>(tracks);

      if (!enable_gen_info) {
        return;
      }

      const auto particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, collision.mcCollision().globalIndex(), cache);

      for (auto& mcParticle : particlesInCollision) { // loop over generated particles

        if (!mcParticle.isPhysicalPrimary()) {
          continue;
        }

        if (mcParticle.pdgCode() == 1000010020) {
          registry.fill(HIST("hReco_EtaPhiPt_Deuteron"), mcParticle.eta(), mcParticle.phi(), mcParticle.pt());
        } else if (mcParticle.pdgCode() == -1000010020) {
          registry.fill(HIST("hReco_EtaPhiPt_Deuteron"), mcParticle.eta(), mcParticle.phi(), mcParticle.pt() * -1);
        }

        if (mcParticle.pdgCode() == 2212) {
          registry.fill(HIST("hReco_EtaPhiPt_Proton"), mcParticle.eta(), mcParticle.phi(), mcParticle.pt());
        } else if (mcParticle.pdgCode() == -2212) {
          registry.fill(HIST("hReco_EtaPhiPt_Proton"), mcParticle.eta(), mcParticle.phi(), mcParticle.pt() * -1);
        }

        if (mcParticle.pdgCode() == 1000020030) {
          registry.fill(HIST("hReco_EtaPhiPt_Helium3"), mcParticle.eta(), mcParticle.phi(), mcParticle.pt());
        } else if (mcParticle.pdgCode() == -1000020030) {
          registry.fill(HIST("hReco_EtaPhiPt_Helium3"), mcParticle.eta(), mcParticle.phi(), mcParticle.pt() * -1);
        }
      }
    }
  }
  PROCESS_SWITCH(singleTrackSelector, processMCRun3, "process MC Run3", false);

  void processGenRun3(aod::McCollisions::iterator const& mcCollision,
                      aod::McParticles const& mcParticles)
  {
    if (!enable_gen_info) {
      return;
    }

    if (std::fabs(mcCollision.posZ()) > _vertexZ) {
      return;
    }

    registry.fill(HIST("hNEvents_MCGen"), 0.5);

    for (auto& mcParticle : mcParticles) { // loop over generated particles

      if (!mcParticle.isPhysicalPrimary()) {
        continue;
      }
      if (mcParticle.pdgCode() == 1000010020) {
        registry.fill(HIST("hGen_EtaPhiPt_Deuteron"), mcParticle.eta(), mcParticle.phi(), mcParticle.pt());
      } else if (mcParticle.pdgCode() == -1000010020) {
        registry.fill(HIST("hGen_EtaPhiPt_Deuteron"), mcParticle.eta(), mcParticle.phi(), mcParticle.pt() * -1);
      }

      if (mcParticle.pdgCode() == 2212) {
        registry.fill(HIST("hGen_EtaPhiPt_Proton"), mcParticle.eta(), mcParticle.phi(), mcParticle.pt());
      } else if (mcParticle.pdgCode() == -2212) {
        registry.fill(HIST("hGen_EtaPhiPt_Proton"), mcParticle.eta(), mcParticle.phi(), mcParticle.pt() * -1);
      }

      if (mcParticle.pdgCode() == 1000020030) {
        registry.fill(HIST("hGen_EtaPhiPt_Helium3"), mcParticle.eta(), mcParticle.phi(), mcParticle.pt());
      } else if (mcParticle.pdgCode() == -1000020030) {
        registry.fill(HIST("hGen_EtaPhiPt_Helium3"), mcParticle.eta(), mcParticle.phi(), mcParticle.pt() * -1);
      }
    }
  }
  PROCESS_SWITCH(singleTrackSelector, processGenRun3, "process MC Gen collisions Run3", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<singleTrackSelector>(cfgc)};
}
