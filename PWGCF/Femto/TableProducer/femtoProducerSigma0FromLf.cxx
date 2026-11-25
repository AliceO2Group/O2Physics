// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoProducerSigma0FromLf.cxx
/// \brief TEST: converts PWGLF Sigma0 candidates (sigma0builder) into femto tables
///
/// Cross-check of the femto Sigma0 builder against the PWGLF one. To be run with
///   strangenessbuilder + strangederivedbuilder (+ LF V0 ML-score tables) + sigma0builder
/// Writes FCollisions + FTracks + FSigma0s + FSigma0Masks (+ optionally FSigma0Extras), so that the femto pair/triplet
/// tasks can run on output unchanged.
///
/// StraCollisions must be 1:1 with the AO2D collisions (strangederivedbuilder with fillEmptyCollisions=true).
/// Masks are set to all bits, so the mask cuts of downstream Sigma0 partitions must be 0.
/// Only the PCM photon branch of sigma0builder is supported (processRealData*); no MC.
///
/// \author Yakiv Paroviak, CTU in Prague, yakiv.paroviak@cern.ch

#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/sigma0Builder.h"
#include "PWGCF/Femto/Core/trackBuilder.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"
#include "PWGLF/DataModel/LFSigmaTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/VectorUtil.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TH1.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

using namespace o2::analysis::femto;

namespace o2::analysis::femto::lfsigma0inputs
{
using Run3PpCollisions = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels, o2::aod::Mults, o2::aod::CentFT0Ms, o2::aod::CentFT0Cs, o2::aod::CentFT0As>;
using Run3FullPidTracks =
  o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TracksDCA,
                o2::aod::pidTPCFullEl, o2::aod::pidTPCFullPi, o2::aod::pidTPCFullKa, o2::aod::pidTPCFullPr, o2::aod::pidTPCFullDe, o2::aod::pidTPCFullTr, o2::aod::pidTPCFullHe,
                o2::aod::pidTOFFullEl, o2::aod::pidTOFFullPi, o2::aod::pidTOFFullKa, o2::aod::pidTOFFullPr, o2::aod::pidTOFFullDe, o2::aod::pidTOFFullTr, o2::aod::pidTOFFullHe,
                o2::aod::pidTOFbeta, o2::aod::pidTOFmass>;
using LfSigma0s = o2::soa::Join<o2::aod::Sigma0Cores, o2::aod::SigmaCollRef>;
using LfSigma0sWithExtras = o2::soa::Join<o2::aod::Sigma0Cores, o2::aod::SigmaCollRef, o2::aod::Sigma0PhotonExtras, o2::aod::Sigma0LambdaExtras>;
} // namespace o2::analysis::femto::lfsigma0inputs

struct FemtoProducerSigma0FromLf {

  o2::framework::Preslice<lfsigma0inputs::Run3FullPidTracks> perColTracks = o2::aod::track::collisionId;

  // femto collision and track builders, configured like in femto-producer
  collisionbuilder::CollisionBuilder<collisionbuilder::ColSelHistName, collisionbuilder::CollisionFilterHistName> collisionBuilder;
  trackbuilder::TrackBuilder<trackbuilder::TrackSelHistName, trackbuilder::TrackFilterHistName> trackBuilder;

  collisionbuilder::ConfCcdb confCcdb;
  collisionbuilder::CollisionBuilderProducts collisionBuilderProducts;
  collisionbuilder::ConfCollisionTables confCollisionTables;
  collisionbuilder::ConfCollisionFilters confCollisionFilters;
  collisionbuilder::ConfCollisionBits confCollisionBits;
  collisionbuilder::ConfCollisionRctFlags confCollisionRctFlags;

  trackbuilder::TrackBuilderProducts trackBuilderProducts;
  trackbuilder::ConfTrackTables confTrackTables;
  trackbuilder::ConfTrackBits confTrackBits;
  trackbuilder::ConfTrackFilters confTrackFilters;

  // femto Sigma0 tables and filters (same configurables as femto-producer-sigma0)
  sigma0builder::Sigma0BuilderProducts sigma0Products;
  sigma0builder::ConfSigma0Tables confSigma0Tables;
  sigma0builder::ConfSigma0Filters confSigma0Filters;

  struct : o2::framework::ConfigurableGroup {
    std::string prefix = std::string("LfSigma0Conversion");
    o2::framework::Configurable<float> lambdaMassWindow{"lambdaMassWindow", 0.015f, "A candidate is written for each Lambda/AntiLambda hypothesis with |m - m_Lambda(PDG)| below this (GeV/c^2). Set <= 0 to use only the hypothesis closer to the PDG mass"};
  } confConversion;

  o2::framework::HistogramRegistry hRegistry{"FemtoProducerSigma0FromLf", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb = {};

  bool mProduceSigma0s = false;
  bool mProduceSigma0Masks = false;
  bool mProduceSigma0Extras = false;

  // AO2D collision index -> rows of Sigma0Cores, rebuilt every dataframe
  std::unordered_map<int64_t, std::vector<int64_t>> mSigma0sPerCollision;

  void init(o2::framework::InitContext& context)
  {
    ccdb->setURL(confCcdb.ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    collisionBuilder.init(&hRegistry, confCollisionFilters, confCollisionBits, confCollisionRctFlags, confCcdb, confCollisionTables, context);
    trackBuilder.init(&hRegistry, confTrackBits, confTrackFilters, confTrackTables, context, collisionBuilder);

    mProduceSigma0s = utils::enableTable("FSigma0s_001", confSigma0Tables.produceSigma0s.value, context);
    mProduceSigma0Masks = utils::enableTable("FSigma0Masks_001", confSigma0Tables.produceSigma0Masks.value, context);
    mProduceSigma0Extras = utils::enableTable("FSigma0Extras_001", confSigma0Tables.produceSigma0Extras.value, context);
    if (utils::enableTable("FLiteSigma0s_001", confSigma0Tables.produceLiteSigma0s.value, context) ||
        utils::enableTable("FSigma0Labels_001", confSigma0Tables.produceSigma0Labels.value, context)) {
      LOG(fatal) << "FLiteSigma0s and FSigma0Labels are not supported by the LF Sigma0 converter";
    }
    if (static_cast<int>(doprocess) + static_cast<int>(doprocessWithExtras) != 1) {
      LOG(fatal) << "Enable exactly one of process and processWithExtras";
    }
    if (mProduceSigma0Extras && !doprocessWithExtras) {
      LOG(fatal) << "FSigma0Extras requested, enable processWithExtras";
    }

    hRegistry.add("hLfSigma0s", "LF Sigma0 conversion;;Entries", o2::framework::HistType::kTH1D, {{6, 0.5, 6.5}});
    auto hist = hRegistry.get<TH1>(HIST("hLfSigma0s"));
    hist->GetXaxis()->SetBinLabel(1, "LF candidates");
    hist->GetXaxis()->SetBinLabel(2, "In accepted collision");
    hist->GetXaxis()->SetBinLabel(3, "Lambda hypotheses");
    hist->GetXaxis()->SetBinLabel(4, "AntiLambda hypotheses");
    hist->GetXaxis()->SetBinLabel(5, "Written (pass filters)");
    hist->GetXaxis()->SetBinLabel(6, "Lambda V0 in other collision");
  }

  template <typename T1, typename T2>
  void mapSigma0sToCollisions(T1 const& lfSigma0s, T2 const& v0Indices)
  {
    mSigma0sPerCollision.clear();
    for (const auto& sigma0 : lfSigma0s) {
      hRegistry.fill(HIST("hLfSigma0s"), 1.);
      // StraCollisions are 1:1 with the AO2D collisions (checked in processDataframe)
      mSigma0sPerCollision[sigma0.straCollisionId()].push_back(sigma0.globalIndex());
      // the Lambda daughters are looked up via V0Indices; it has to point to the same collision
      if (v0Indices.rawIteratorAt(sigma0.lambdaV0ID()).collisionId() != sigma0.straCollisionId()) {
        hRegistry.fill(HIST("hLfSigma0s"), 6.);
      }
    }
  }

  template <typename T1, typename T2, typename T3, typename T4>
  void fillSigma0s(T1 const& col, T2 const& lfSigma0s, T3 const& v0Indices, T4 const& tracks)
  {
    auto candidates = mSigma0sPerCollision.find(col.globalIndex());
    if (candidates == mSigma0sPerCollision.end()) {
      return;
    }
    for (const auto row : candidates->second) {
      auto sigma0 = lfSigma0s.rawIteratorAt(row);
      hRegistry.fill(HIST("hLfSigma0s"), 2.);

      // Lambda hypotheses accepted for this candidate
      const float dLambda = std::abs(sigma0.lambdaMass() - o2::constants::physics::MassLambda0);
      const float dAntiLambda = std::abs(sigma0.antilambdaMass() - o2::constants::physics::MassLambda0);
      std::array<bool, 2> hypotheses = {false, false}; // {Lambda, AntiLambda}
      if (confConversion.lambdaMassWindow.value > 0.f) {
        hypotheses = {dLambda < confConversion.lambdaMassWindow.value, dAntiLambda < confConversion.lambdaMassWindow.value};
      } else {
        hypotheses = {dLambda <= dAntiLambda, dAntiLambda < dLambda};
      }

      // Sigma0 kinematics with the femto convention (PDG Lambda mass, massless photon)
      const ROOT::Math::PxPyPzMVector vecLambda{sigma0.lambdaPx(), sigma0.lambdaPy(), sigma0.lambdaPz(), o2::constants::physics::MassLambda0};
      const ROOT::Math::PxPyPzMVector vecPhoton{sigma0.photonPx(), sigma0.photonPy(), sigma0.photonPz(), o2::constants::physics::MassPhoton};
      const auto vecSigma0 = vecLambda + vecPhoton;
      const float pt = vecSigma0.Pt();
      const float eta = vecSigma0.Eta();
      const float phi = RecoDecay::constrainAngle(vecSigma0.Phi());
      const float mass = vecSigma0.M();
      if (!(mass > confSigma0Filters.massMin.value && mass < confSigma0Filters.massMax.value &&
            pt > confSigma0Filters.ptMin.value && pt < confSigma0Filters.ptMax.value &&
            eta > confSigma0Filters.etaMin.value && eta < confSigma0Filters.etaMax.value &&
            phi >= confSigma0Filters.phiMin.value && phi < confSigma0Filters.phiMax.value)) {
        continue;
      }

      for (int iHyp = 0; iHyp < 2; iHyp++) {
        if (!hypotheses.at(iHyp)) {
          continue;
        }
        hRegistry.fill(HIST("hLfSigma0s"), 3. + iHyp);
        const float sign = iHyp == 0 ? 1.f : -1.f;

        // register the Lambda daughters in FTracks (writes the collision first if not done yet)
        collisionBuilder.template fillCollision<modes::System::kPP_Run3>(collisionBuilderProducts, col);
        auto lambdaIndices = v0Indices.rawIteratorAt(sigma0.lambdaV0ID());
        auto posDaughter = tracks.rawIteratorAt(lambdaIndices.posTrackId());
        auto negDaughter = tracks.rawIteratorAt(lambdaIndices.negTrackId());
        const int64_t posDauIndex = trackBuilder.template getDaughterIndex<modes::Track::kV0Daughter>(posDaughter, trackBuilderProducts, collisionBuilder);
        const int64_t negDauIndex = trackBuilder.template getDaughterIndex<modes::Track::kV0Daughter>(negDaughter, trackBuilderProducts, collisionBuilder);

        if (mProduceSigma0s) {
          sigma0Products.producedSigma0s(collisionBuilder.collisionIndex(), sign * pt, eta, phi, mass, posDauIndex, negDauIndex, -1);
        }
        if (mProduceSigma0Masks) {
          sigma0Products.producedSigma0Masks(std::numeric_limits<datatypes::V0MaskType>::max(), std::numeric_limits<datatypes::PhotonMaskType>::max());
        }
        if (mProduceSigma0Extras) {
          if constexpr (requires { sigma0.lambdaCosPA(); sigma0.photonCosPA(); }) {
            sigma0Products.producedSigma0Extras(iHyp == 0 ? sigma0.lambdaMass() : sigma0.antilambdaMass(),
                                                vecLambda.Pt(),
                                                sigma0.lambdaCosPA(),
                                                sigma0.lambdaDCADau(),
                                                sigma0.lambdaRadius(),
                                                vecPhoton.Pt(),
                                                sigma0.photonCosPA(),
                                                sigma0.photonRadius(),
                                                -1.f, // photon chi2/NDF is not available for LF photons
                                                ROOT::Math::VectorUtil::Angle(vecLambda, vecPhoton));
          }
        }
        hRegistry.fill(HIST("hLfSigma0s"), 5.);
      }
    }
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
  void processDataframe(T1 const& cols, T2 const& bcs, T3 const& tracks, T4 const& lfSigma0s, T5 const& v0Indices, T6 const& v0Cores, T7 const& straCols)
  {
    if (straCols.size() != cols.size()) {
      LOG(fatal) << "StraCollisions (" << straCols.size() << ") are not 1:1 with the AO2D collisions (" << cols.size() << "); "
                 << "run strangederivedbuilder with fillEmptyCollisions=true and fillOnlySelectedCollisions=false";
    }
    if (v0Indices.size() != v0Cores.size()) {
      LOG(fatal) << "V0Indices (" << v0Indices.size() << ") and V0Cores (" << v0Cores.size() << ") are not row-aligned; "
                 << "both have to be produced by strangenessbuilder in this workflow";
    }
    trackBuilder.reset(tracks);
    mapSigma0sToCollisions(lfSigma0s, v0Indices);
    for (const auto& col : cols) {
      auto tracksThisCol = tracks.sliceBy(perColTracks, col.globalIndex());
      collisionBuilder.reset();
      auto bc = col.template bc_as<T2>();
      collisionBuilder.template initCollision<modes::System::kPP_Run3>(bc, col, tracksThisCol, ccdb, hRegistry);
      if (!collisionBuilder.checkCollision(col)) {
        continue;
      }
      if (collisionBuilder.isPassThrough()) {
        collisionBuilder.template fillCollision<modes::System::kPP_Run3>(collisionBuilderProducts, col);
      }
      auto tracksWithItsPid = o2::soa::Attach<lfsigma0inputs::Run3FullPidTracks, o2::aod::pidits::ITSNSigmaEl, o2::aod::pidits::ITSNSigmaPi, o2::aod::pidits::ITSNSigmaKa,
                                              o2::aod::pidits::ITSNSigmaPr, o2::aod::pidits::ITSNSigmaDe, o2::aod::pidits::ITSNSigmaTr, o2::aod::pidits::ITSNSigmaHe>(tracksThisCol);
      trackBuilder.template fillTracks<modes::System::kPP_Run3>(col, collisionBuilder, collisionBuilderProducts, tracksWithItsPid, trackBuilderProducts);
      fillSigma0s(col, lfSigma0s, v0Indices, tracks);
    }
  }

  void process(lfsigma0inputs::Run3PpCollisions const& cols, o2::aod::BCsWithTimestamps const& bcs, lfsigma0inputs::Run3FullPidTracks const& tracks,
               lfsigma0inputs::LfSigma0s const& lfSigma0s, o2::aod::V0Indices const& v0Indices, o2::aod::V0Cores const& v0Cores, o2::aod::StraCollisions const& straCols)
  {
    processDataframe(cols, bcs, tracks, lfSigma0s, v0Indices, v0Cores, straCols);
  }
  PROCESS_SWITCH(FemtoProducerSigma0FromLf, process, "Convert LF Sigma0s (pp, Run 3)", true);

  void processWithExtras(lfsigma0inputs::Run3PpCollisions const& cols, o2::aod::BCsWithTimestamps const& bcs, lfsigma0inputs::Run3FullPidTracks const& tracks,
                         lfsigma0inputs::LfSigma0sWithExtras const& lfSigma0s, o2::aod::V0Indices const& v0Indices, o2::aod::V0Cores const& v0Cores, o2::aod::StraCollisions const& straCols)
  {
    processDataframe(cols, bcs, tracks, lfSigma0s, v0Indices, v0Cores, straCols);
  }
  PROCESS_SWITCH(FemtoProducerSigma0FromLf, processWithExtras, "Convert LF Sigma0s (pp, Run 3) and fill FSigma0Extras from the LF photon/lambda extras", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  return o2::framework::WorkflowSpec{adaptAnalysisTask<FemtoProducerSigma0FromLf>(context)};
}
