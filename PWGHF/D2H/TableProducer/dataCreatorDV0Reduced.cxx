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

/// \file dataCreatorDV0Reduced.cxx
/// \brief Creation of D-V0 pairs
///
/// \author Luca Aglietta <luca.aglietta@cern.ch>, UniTO Turin

#include <cmath>
#include <map>

#include "DetectorsBase/Propagator.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsBfieldCCDB.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;

// event types
enum Event : uint8_t {
  Processed = 0,
  NoDV0Selected,
  DV0Selected,
  kNEvent
};

enum DecayChannel : uint8_t {
  Ds1ToDstarK0s = 0,
  DsStar2ToDplusK0s,
  Xc3055ToDplusLambda
};

enum V0_type : int8_t {
  K0s = 0,
  Lambda
};

enum D_type : int8_t {
  Dplus = 1,
  Dstar
};

/// Creation of D-V0 pairs
struct HfDataCreatorDV0Reduced {
  // Produces AOD tables to store track information
  Produces<aod::HfRedCollision> hfReducedCollision;  // Defined in PWGLF/DataModel/LFStrangenessTables.h
  Produces<aod::HfOrigColCounts> hfCollisionCounter; // Defined in PWGHF/D2H/DataModel/ReducedDataModel.h

  Produces<aod::HfRedVzeros> hfCandV0;   // Defined in PWGHF/D2H/DataModel/ReducedDataModel.h
  Produces<aod::HfRed3PrNoTrks> hfCandD; // Defined in PWGHF/D2H/DataModel/ReducedDataModel.h

  // CCDB configuration
  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  Configurable<bool> propagateV0toPV{"propagateV0toPV", false, "Enable or disable V0 propagation to V0"};

  int runNumber{0}; // needed to detect if the run changed and trigger update of calibrations etc.
  // selection
  Configurable<int> selectionFlagDplus{"selectionFlagDplus", 7, "Selection Flag for D"};
  Configurable<bool> selectionFlagDstarToD0Pi{"selectionFlagDstarToD0Pi", true, "Selection Flag for D* decay to D0 & Pi"};

  // material correction for track propagation
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  HfHelper hfHelper;

  bool isHfCandResoConfigFilled = false;

  using CandsDplusFiltered = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>>;
  using CandDstarFiltered = soa::Filtered<soa::Join<aod::HfD0FromDstar, aod::HfCandDstar, aod::HfSelDstarToD0Pi>>;

  Filter filterSelectDplus = (aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlagDplus);
  Filter filterSelectedCandDstar = (aod::hf_sel_candidate_dstar::isSelDstarToD0Pi == selectionFlagDstarToD0Pi);

  Preslice<CandsDplusFiltered> candsDplusPerCollision = aod::track_association::collisionId;
  Preslice<CandDstarFiltered> candsDstarPerCollision = aod::track_association::collisionId;
  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  Preslice<aod::V0Datas> candsV0PerCollision = aod::track_association::collisionId;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    // histograms
    constexpr int kNBinsEvents = kNEvent;
    std::string labels[kNBinsEvents];
    labels[Event::Processed] = "processed";
    labels[Event::NoDV0Selected] = "without DV0 pairs";
    labels[Event::DV0Selected] = "with DV0 pairs";
    static const AxisSpec axisEvents = {kNBinsEvents, 0.5, kNBinsEvents + 0.5, ""};
    registry.add("hEvents", "Events;;entries", HistType::kTH1F, {axisEvents});
    for (int iBin = 0; iBin < kNBinsEvents; iBin++) {
      registry.get<TH1>(HIST("hEvents"))->GetXaxis()->SetBinLabel(iBin + 1, labels[iBin].data());
    }
    registry.add("hMassDplus", "Dplus candidates;inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{100, 1.7, 2}}});
    registry.add("hMassDstar", "Dstar candidates;inv. mass (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{100, 0.05, 0.25}}});
    registry.add("hMassK0s", "K0^{s} candidates;inv. mass (#pi^{#plus}#pi^{#minus}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{100, 0.35, 0.65}}});
    registry.add("hMassLambda", "Lambda candidates;inv. mass (p #pi^{#minus}) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{100, 1.05, 1.35}}});
    registry.add("hPtDplus", "D^{#minus} candidates;D^{#minus} candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("hPtDstar", "D^{*} candidates;D^{*} candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("hPtV0", "V0 candidates;V0 candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("hMassDs1", "Ds1 candidates;m_{Ds1} - m_{D^{*}} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{100, 0.45, 0.7}}});
    registry.add("hMassDsStar2", "Ds^{*}2 candidates; Ds^{*}2 - m_{D^{#plus}} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{100, 0.4, 1}}});
    registry.add("hMassXcRes", "XcRes candidates; XcRes - m_{D^{#plus}} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{100, 1., 1.4}}});
    registry.add("hV0_type", "XcRes candidates; XcRes - m_{D^{#plus}} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{100, -3, 3}}});
    registry.add("hD_type", "XcRes candidates; XcRes - m_{D^{#plus}} (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{100, -3, 3}}});

    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    ccdbApi.init(url);
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));
  }

  template <uint8_t DecayChannel, typename CCands>
  void runDataCreation(aod::Collision const& collision,
                       C const& candsD,
                       aod::V0Datas const& V0s,
                       aod::BCsWithTimestamps const& bcs)
  {
    // helpers for ReducedTables filling
    int indexHfReducedCollision = hfReducedCollision.lastIndex() + 1;
    // std::map where the key is the V0.globalIndex() and
    // the value is the V0 index in the table of the selected v0s
    std::map<int64_t, int64_t> selectedV0s;
    bool fillHfReducedCollision = false;
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (runNumber != bc.runNumber()) {
      LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
      o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbPathGrpMag, bc.timestamp());
      if (grpo == nullptr) {
        LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      o2::base::Propagator::initFieldFromGRP(grpo);
      // setMatLUT only after magfield has been initalized
      o2::base::Propagator::Instance()->setMatLUT(lut);
      runNumber = bc.runNumber();
    }

    for (const auto& candD : candsD) {
      // initialize variables depending on decay channel
      bool fillHfCandD = false;
      float invMassD;
      float massD;
      float massV0;
      std::array<float, 3> pVecD;
      std::array<float, 3> secondaryVertexD;
      std::array<int, 3> prongIdsD;
      int8_t v0_type;
      int8_t d_type;

      if constexpr (std::is_same<C, CandDstarFiltered>::value) {
        if (candD.signSoftPi() > 0)
          invMassD = candD.invMassDstar() - candD.invMassD0();
        else
          invMassD = candD.invMassAntiDstar() - candD.invMassD0Bar();
        massD = MassDStar;
        pVecD = candD.pVector();
        secondaryVertexD[0] = candD.xSecondaryVertexD0();
        secondaryVertexD[1] = candD.ySecondaryVertexD0();
        secondaryVertexD[2] = candD.zSecondaryVertexD0();
        prongIdsD[0] = candD.prong0Id();
        prongIdsD[1] = candD.prong1Id();
        prongIdsD[2] = candD.prongPiId();
        v0_type = V0_type::K0s;
        d_type = candD.signSoftPi() * D_type::Dstar;
        massV0 = MassK0Short;
      } else if constexpr (std::is_same<C, CandsDplusFiltered>::value) {
        switch (DecayChannel) {
          case DecayChannel::DsStar2ToDplusK0s:
            invMassD = hfHelper.invMassDplusToPiKPi(candD);
            massD = MassDPlus;
            pVecD = candD.pVector();
            secondaryVertexD[0] = candD.xSecondaryVertex();
            secondaryVertexD[1] = candD.ySecondaryVertex();
            secondaryVertexD[2] = candD.zSecondaryVertex();
            prongIdsD[0] = candD.prong0Id();
            prongIdsD[1] = candD.prong1Id();
            prongIdsD[2] = candD.prong2Id();
            v0_type = V0_type::K0s;
            d_type = candD.sign() * D_type::Dplus;
            massV0 = MassK0Short;
            break;

          case DecayChannel::Xc3055ToDplusLambda:
            invMassD = hfHelper.invMassDplusToPiKPi(candD);
            massD = MassDPlus;
            pVecD = candD.pVector();
            secondaryVertexD[0] = candD.xSecondaryVertex();
            secondaryVertexD[1] = candD.ySecondaryVertex();
            secondaryVertexD[2] = candD.zSecondaryVertex();
            prongIdsD[0] = candD.prong0Id();
            prongIdsD[1] = candD.prong1Id();
            prongIdsD[2] = candD.prong2Id();
            v0_type = candD.sign() * V0_type::Lambda;
            d_type = candD.sign() * D_type::Dplus;
            massV0 = MassLambda0;
            break;

          default:
            LOG(warning) << "Decay channel not valid please choose between Ds1ToDstarK0s, DsStar2ToDplusK0s, Xc3055ToDplusLambda";
            break;
        } // switch
      }   // else if

      // Loop on V0 candidates
      for (const auto& v0 : V0s) {
        // propagate V0 to primary vertex
        std::array<float, 3> pVecV0 = {v0.px(), v0.py(), v0.pz()};
        if (propagateV0toPV) {
          std::array<float, 3> pVecV0Orig = {v0.px(), v0.py(), v0.pz()};
          std::array<float, 3> posVecV0 = {v0.x(), v0.y(), v0.z()};
          gpu::gpustd::array<float, 2> dcaInfo;
          auto trackParK0 = o2::track::TrackPar(posVecV0, pVecV0Orig, 0, true);
          trackParK0.setPID(o2::track::PID::K0);
          trackParK0.setAbsCharge(0);
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParK0, 2.f, matCorr, &dcaInfo);
          getPxPyPz(trackParK0, pVecV0);
        }

        float ptV0 = sqrt(pVecV0[0] * pVecV0[0] + pVecV0[1] * pVecV0[1]);
        auto invMass2DV0 = RecoDecay::m2(std::array{pVecD, pVecV0}, std::array{massD, massV0});
        // LOG(info) << "V0 p before propagation: " << pVecV0Orig[0] << "," << pVecV0Orig[1] << "," << pVecV0Orig[2];
        // LOG(info) << "V0 p after propagation: " << pVecV0[0] << "," << pVecV0[1] << "," << pVecV0[2];

        // fill histos
        registry.fill(HIST("hPtV0"), ptV0);
        registry.fill(HIST("hV0_type"), v0_type);

        switch (DecayChannel) {
          case DecayChannel::Ds1ToDstarK0s:
            registry.fill(HIST("hMassK0s"), v0.mK0Short());
            registry.fill(HIST("hMassDs1"), sqrt(invMass2DV0) - invMassD);
            break;
          case DecayChannel::DsStar2ToDplusK0s:
            registry.fill(HIST("hMassK0s"), v0.mK0Short());
            registry.fill(HIST("hMassDsStar2"), sqrt(invMass2DV0) - invMassD);
            break;
          case DecayChannel::Xc3055ToDplusLambda:
            if (v0_type > 0)
              registry.fill(HIST("hMassLambda"), v0.mLambda());
            else
              registry.fill(HIST("hMassLambda"), v0.mAntiLambda());
            registry.fill(HIST("hMassXcRes"), sqrt(invMass2DV0) - invMassD);
            break;
          default:
            break;
        }

        // fill V0 table
        // if information on V0 already stored, go to next V0
        if (!selectedV0s.count(v0.globalIndex())) {
          hfCandV0(v0.posTrackId(), v0.negTrackId(),
                   indexHfReducedCollision,
                   v0.x(), v0.y(), v0.z(),
                   v0.mK0Short(), v0.mLambda(), v0.mAntiLambda(),
                   pVecV0[0], pVecV0[1], pVecV0[2],
                   v0_type);
          selectedV0s[v0.globalIndex()] = hfCandV0.lastIndex();
        }
        fillHfCandD = true;
      } // V0 loop

      if (fillHfCandD) { // fill candDplus table only once per D candidate, only if at least one V0 is found
        hfCandD(prongIdsD[0], prongIdsD[1], prongIdsD[2],
                indexHfReducedCollision,
                secondaryVertexD[0], secondaryVertexD[1], secondaryVertexD[2],
                invMassD,
                pVecD[0], pVecD[1], pVecD[2],
                d_type);
        fillHfReducedCollision = true;
        switch (DecayChannel) {
          case DecayChannel::Ds1ToDstarK0s:
            registry.fill(HIST("hMassDstar"), invMassD);
            registry.fill(HIST("hPtDstar"), candD.pt());
            break;
          case DecayChannel::DsStar2ToDplusK0s:
            registry.fill(HIST("hMassDplus"), invMassD);
            registry.fill(HIST("hPtDplus"), candD.pt());
            break;
          case DecayChannel::Xc3055ToDplusLambda:
            registry.fill(HIST("hMassDplus"), invMassD);
            registry.fill(HIST("hPtDplus"), candD.pt());
            break;
          default:
            break;
        }
        registry.fill(HIST("hPtDplus"), candD.pt());
        registry.fill(HIST("hD_type"), d_type);
      }
    } // candsD loop
    registry.fill(HIST("hEvents"), 1 + Event::Processed);
    if (!fillHfReducedCollision) {
      registry.fill(HIST("hEvents"), 1 + Event::NoDV0Selected);
      return;
    }
    registry.fill(HIST("hEvents"), 1 + Event::DV0Selected);
    // fill collision table if it contains a DPi pair a minima
    hfReducedCollision(collision.posX(), collision.posY(), collision.posZ(),
                       collision.covXX(), collision.covXY(), collision.covYY(),
                       collision.covXZ(), collision.covYZ(), collision.covZZ(),
                       0);
  }

  void processDsStar2(aod::Collisions const& collisions,
                      CandsDplusFiltered const& candsDplus,
                      aod::TrackAssoc const& trackIndices,
                      aod::V0Datas const& V0s,
                      aod::BCsWithTimestamps const& bcs)
  {
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize());

    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDplus.sliceBy(candsDplusPerCollision, thisCollId);
      auto V0sThisColl = V0s.sliceBy(candsV0PerCollision, thisCollId);
      runDataCreation<DecayChannel::DsStar2ToDplusK0s, CandsDplusFiltered>(collision, candsDThisColl, V0sThisColl, bcs);
    }
  }
  PROCESS_SWITCH(HfDataCreatorDV0Reduced, processDsStar2, "Process DsStar2 to Dplus K0s without MC info and without ML info", true);

  void processDs1(aod::Collisions const& collisions,
                  CandDstarFiltered const& candsDstar,
                  aod::TrackAssoc const& trackIndices,
                  aod::V0Datas const& V0s,
                  aod::BCsWithTimestamps const& bcs)
  {
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize());

    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDstar.sliceBy(candsDstarPerCollision, thisCollId);
      auto V0sThisColl = V0s.sliceBy(candsV0PerCollision, thisCollId);
      runDataCreation<DecayChannel::Ds1ToDstarK0s, CandDstarFiltered>(collision, candsDThisColl, V0sThisColl, bcs);
    }
  }
  PROCESS_SWITCH(HfDataCreatorDV0Reduced, processDs1, "Process Ds1 to DStar K0s without MC info and without ML info", false);

  void processXc(aod::Collisions const& collisions,
                 CandsDplusFiltered const& candsDplus,
                 aod::TrackAssoc const& trackIndices,
                 aod::V0Datas const& V0s,
                 aod::BCsWithTimestamps const& bcs)
  {
    // handle normalization by the right number of collisions
    hfCollisionCounter(collisions.tableSize());

    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto candsDThisColl = candsDplus.sliceBy(candsDplusPerCollision, thisCollId);
      auto V0sThisColl = V0s.sliceBy(candsV0PerCollision, thisCollId);
      runDataCreation<DecayChannel::Xc3055ToDplusLambda, CandsDplusFiltered>(collision, candsDThisColl, V0sThisColl, bcs);
    }
  }
  PROCESS_SWITCH(HfDataCreatorDV0Reduced, processXc, "Process Xc to Dplus Lambda without MC info and without ML info", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfDataCreatorDV0Reduced>(cfgc)};
}
