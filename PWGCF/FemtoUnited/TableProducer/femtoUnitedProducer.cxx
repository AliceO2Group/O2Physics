// Copyright 2019-2024 CERN andhcopyright holders of ALICE O2.
//
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoUnitedProducer.cxx
/// \brief Tasks that produces the all femto tables
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include "PWGCF/FemtoUnited/Core/cascadeSelection.h"
#include "PWGCF/FemtoUnited/Core/collisionSelection.h"
#include "PWGCF/FemtoUnited/Core/dataTypes.h"
#include "PWGCF/FemtoUnited/Core/femtoUtils.h"
#include "PWGCF/FemtoUnited/Core/trackSelection.h"
#include "PWGCF/FemtoUnited/Core/twoTrackResonanceSelection.h"
#include "PWGCF/FemtoUnited/Core/v0Selection.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoCascadesDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoCollisionsDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTwoTrackResonancesDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoV0sDerived.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/Expressions.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "fairlogger/Logger.h"

#include <string>
#include <unordered_map>
#include <vector>

using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::analysis::femtounited;

namespace o2::analysis::femtounited
{
namespace consumeddata
{
using Run3PpCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms>;
using Run3PpWithoutCentCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;

using Run3FullPidTracks =
  soa::Join<Tracks, TracksExtra, TracksDCA,
            pidTPCFullEl, pidTPCFullPi, pidTPCFullKa, pidTPCFullPr, pidTPCFullDe, pidTPCFullTr, pidTPCFullHe,
            pidTOFFullEl, pidTOFFullPi, pidTOFFullKa, pidTOFFullPr, pidTOFFullDe, pidTOFFullTr, pidTOFFullHe,
            pidTOFbeta>;

using Run3PpVzeros = V0Datas;

using Run3PpCascades = CascDatas;

} // namespace consumeddata
} // namespace o2::analysis::femtounited

struct FemtoUnitedProducer {
  SliceCache cache;
  // preslicing
  Preslice<Tracks> perColTracks = track::collisionId;

  // produced objectes
  struct : ProducesGroup {
    Produces<FUCols> producedCollision;

    Produces<FUTracks> producedTracks;
    Produces<FUTrackMasks> producedTrackMasks;
    Produces<FUTrackDCAs> producedTrackDCAs;
    Produces<FUTrackExtras> producedTrackExtras;
    Produces<FUTrackPids> producedTrackPids;

    Produces<FULambdas> producedLambdas;
    Produces<FULambdaMasks> producedLambdaMasks;
    Produces<FULambdaExtras> producedLambdaExtras;

    Produces<FUK0shorts> producedK0shorts;
    Produces<FUK0shortMasks> producedK0shortMasks;
    Produces<FUK0shortExtras> producedK0shortExtras;

    Produces<FUXis> producedXis;
    Produces<FUXiMasks> producedXiMasks;
    Produces<FUXiExtras> producedXiExtras;

    Produces<FUOmegas> producedOmegas;
    Produces<FUOmegaMasks> producedOmegaMasks;
    Produces<FUOmegaExtras> producedOmegaExtras;

    Produces<FUPhis> producedPhis;
    Produces<FUPhiMasks> producedPhiMasks;

    Produces<FUKstars> producedKstars;
    Produces<FUKstarMasks> producedKstarMasks;

    Produces<FURhos> producedRhos;
    Produces<FURhoMasks> producedRhoMasks;
  } products;

  // configurables
  struct : ConfigurableGroup {
    std::string prefix = std::string("General");
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "URL to ccdb"};
    Configurable<std::string> grpPath{"grpPath", "GLO/Config/GRPMagField", "Path to GRP object (Run3 -> GLO/Config/GRPMagField/Run2 -> GLO/GRP/GRP"};
    Configurable<bool> produceExtraTables{"produceExtraTables", false, "Flag to produce extra tables (for all actived tables)"};
    Configurable<bool> produceK0short{"produceK0short", false, "Flag to produce K0shorts"};
    Configurable<bool> produceLambda{"produceLambda", false, "Flag to produce Lambda"};
    Configurable<bool> produceXi{"produceXi", false, "Flag to produce Xi"};
    Configurable<bool> produceOmega{"produceOmega", false, "Flag to produce Omega"};
    Configurable<bool> producePhi{"producePhi", false, "Flag to produce Phi"};
    Configurable<bool> produceKstar0{"produceKstar0", false, "Flag to produce Kstar"};
    Configurable<bool> produceRho0{"produceRho0", false, "Flag to produce Rho"};
  } ConfOptions;

  // Event selections
  collisionselection::ConfCollisionSelection confCollisionFilter;
  Filter collisionFilter = o2::aod::collision::posZ >= confCollisionFilter.vtxZMin &&
                           o2::aod::collision::posZ <= confCollisionFilter.vtxZMax;
  collisionselection::CollisionSelection collisionSel;

  // filters for tracks
  trackselection::ConfTrackFilters confTrackFilters;
  Filter trackFilter = track::pt >= confTrackFilters.ptMin && track::pt <= confTrackFilters.ptMax &&
                       track::eta >= confTrackFilters.etaMin && track::eta <= confTrackFilters.etaMax &&
                       track::phi >= confTrackFilters.phiMin && track::phi <= confTrackFilters.phiMax;
  // track bits
  trackselection::ConfTrackBits confTrackBits;
  trackselection::TrackSelection trackSel;

  // lambda filters
  // most v0 columns are now dynamic columns, so we cannot prefilter anymore
  v0selection::ConfV0Filters confV0Filters;

  // K0short bits
  v0selection::ConfK0shortBits confK0shortBits;
  v0selection::V0Selection<modes::V0::kK0short> k0shortSel;

  // lambda bits
  v0selection::ConfLambdaBits confLambdaBits;
  v0selection::V0Selection<modes::V0::kLambda> lambdaSel;
  v0selection::V0Selection<modes::V0::kAntiLambda> antiLambdaSel;

  // cascade filters
  cascadeselection::ConfCascadeFilters confCascadeFilters;

  // xi bits
  cascadeselection::ConfXiBits confXiBits;
  cascadeselection::CascadeSelection<modes::Cascade::kXi> xiSel;

  // omega bits
  cascadeselection::ConfOmegaBits confOmegaBits;
  cascadeselection::CascadeSelection<modes::Cascade::kOmega> omegaSel;

  // resonance filters
  twotrackresonanceselection::ConfTwoTrackResonanceDaughterFilters confResonanceDaughterFilters;
  twotrackresonanceselection::ConfRhoFilters confRhoFilters;
  twotrackresonanceselection::ConfPhiFilters confPhiFilters;
  twotrackresonanceselection::ConfKstarFilters confKstarFilters;

  // rho bits
  twotrackresonanceselection::ConfRho0Bits confRho0Bits;
  twotrackresonanceselection::TwoTrackResonanceSelection<modes::TwoTrackResonance::kRho0> rho0Sels;

  // phi bits
  twotrackresonanceselection::ConfPhiBits confPhiBits;
  twotrackresonanceselection::TwoTrackResonanceSelection<modes::TwoTrackResonance::kPhi> phiSels;

  // kstar bits
  twotrackresonanceselection::ConfKstar0Bits confKstar0Bits;
  twotrackresonanceselection::TwoTrackResonanceSelection<modes::TwoTrackResonance::kKstar0> kstar0Sels;
  twotrackresonanceselection::TwoTrackResonanceSelection<modes::TwoTrackResonance::kKstarBar0> kstarBar0Sels;

  Partition<Filtered<consumeddata::Run3FullPidTracks>> partitionPositiveDaughters =
    (track::signed1Pt > 0.f) &&
    (track::pt > confResonanceDaughterFilters.posDauPtMin && track::pt < confResonanceDaughterFilters.posDauPtMax) &&
    (track::eta > confResonanceDaughterFilters.posDauEtaMin && track::eta < confResonanceDaughterFilters.posDauEtaMax) &&
    (track::phi > confResonanceDaughterFilters.posDauPhiMin && track::phi < confResonanceDaughterFilters.posDauPhiMax);
  Partition<Filtered<consumeddata::Run3FullPidTracks>> partitionNegativeDaughters =
    (track::signed1Pt < 0.f) &&
    (track::pt > confResonanceDaughterFilters.negDauPtMin && track::pt < confResonanceDaughterFilters.negDauPtMax) &&
    (track::eta > confResonanceDaughterFilters.negDauEtaMin && track::eta < confResonanceDaughterFilters.negDauEtaMax) &&
    (track::phi > confResonanceDaughterFilters.negDauPhiMin && track::phi < confResonanceDaughterFilters.negDauPhiMax);

  // histogramming
  HistogramRegistry hRegistry{"FemtoProducer", {}, OutputObjHandlingPolicy::AnalysisObject};

  // data members
  int runNumber = -1;
  float magField = 0.f;
  Service<o2::ccdb::BasicCCDBManager> ccdb;            /// Accessing the CCDB
  std::unordered_map<int64_t, int64_t> indexMapTracks; // for mapping tracks to lambdas, cascades and resonances

  bool produceV0s = false;
  bool produceResonances = false;
  bool produceCascades = false;

  void initFromCcdb(o2::aod::BCsWithTimestamps::iterator const& bc)
  {
    if (runNumber == bc.runNumber())
      return;
    auto timestamp = bc.timestamp();

    static o2::parameters::GRPMagField* grpo = nullptr;
    grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ConfOptions.grpPath.value, timestamp);
    if (grpo == nullptr) {
      LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
      return;
    }
    magField = 0.1 * grpo->getNominalL3Field(); // get magnetic field in tesla
    runNumber = bc.runNumber();
  };

  void init(InitContext& /*contex*/)
  {
    // init ccdb
    ccdb->setURL(ConfOptions.ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    // collision selection
    collisionSel.configure(confCollisionFilter);

    // init track selection objects
    trackSel.configure(confTrackBits, confTrackFilters);

    // init v0 selection ojects
    lambdaSel.configure(confLambdaBits, confV0Filters);
    antiLambdaSel.configure(confLambdaBits, confV0Filters);
    k0shortSel.configure(confK0shortBits, confV0Filters);

    // cascade selections
    xiSel.configure(confXiBits, confCascadeFilters);
    omegaSel.configure(confOmegaBits, confCascadeFilters);

    // resonance selections
    rho0Sels.configure(confRho0Bits, confRhoFilters, confResonanceDaughterFilters);
    phiSels.configure(confPhiBits, confPhiFilters, confResonanceDaughterFilters);
    kstar0Sels.configure(confKstar0Bits, confKstarFilters, confResonanceDaughterFilters);
    kstarBar0Sels.configure(confKstar0Bits, confKstarFilters, confResonanceDaughterFilters);

    trackSel.printSelections(trackselection::TrackSelsName, trackselection::TrackSelsToString);

    // tracks are produced by default
    if (ConfOptions.produceK0short.value || ConfOptions.produceLambda.value) {
      produceV0s = true;
      k0shortSel.printSelections(v0selection::k0shortSelsName, v0selection::V0SelesNames);
      lambdaSel.printSelections(v0selection::lambdaSelsName, v0selection::V0SelesNames);
      antiLambdaSel.printSelections(v0selection::antiLambdaSelsName, v0selection::V0SelesNames);
    }
    if (ConfOptions.producePhi.value || ConfOptions.produceRho0.value || ConfOptions.produceKstar0.value) {
      produceResonances = true;
    }
    if (ConfOptions.produceXi.value || ConfOptions.produceOmega.value) {
      produceCascades = true;
      xiSel.printSelections(cascadeselection::xiSelsName, cascadeselection::CascadeSelsNames);
      omegaSel.printSelections(cascadeselection::omegaSelsName, cascadeselection::CascadeSelsNames);
    }
    if ((doprocessTracksRun3pp + doprocessTracksV0sRun3pp + doprocessTracksV0sCascadesRun3pp) > 1) {
      LOG(fatal) << "Only one process function can be activated.";
    }
  }

  template <modes::Mode mode, modes::Track type, typename T>
  int64_t getDaughterIndex(const T& daughter)
  {
    auto result = utils::getTrackIndex(daughter.globalIndex(), indexMapTracks);
    if (result) {
      return result.value();
    } else {
      fillTrack<mode, type>(daughter);
      int64_t idx = products.producedTracks.lastIndex();
      indexMapTracks.emplace(daughter.globalIndex(), idx);
      return idx;
    }
  }

  template <modes::System sys, typename T>
  void fillCollision(T const& col)
  {
    if constexpr (!modes::isFlagSet(sys, modes::System::kNoCentCal)) {
      products.producedCollision(col.posZ(),
                                 col.multNTracksPV(),
                                 col.centFT0M(),
                                 collisionSel.getSphericity(),
                                 collisionSel.getMagneticField());
    }

    if constexpr (modes::isFlagSet(sys, modes::System::kNoCentCal)) {
      products.producedCollision(col.posZ(),
                                 col.multNTracksPV(),
                                 0,
                                 collisionSel.getSphericity(),
                                 collisionSel.getMagneticField());
    }
  }

  template <modes::Mode mode, modes::Track type, typename T1>
  void fillTrack(T1 const& track)
  {
    if constexpr (modes::isFlagSet(mode, modes::Mode::kANALYSIS)) {
      products.producedTracks(products.producedCollision.lastIndex(),
                              track.pt() * track.sign(),
                              track.eta(),
                              track.phi());
      if constexpr (type == modes::Track::kPrimaryTrack) {
        products.producedTrackMasks(trackSel.getBitmask());
      } else {
        products.producedTrackMasks(static_cast<femtodatatypes::TrackMaskType>(0u));
      }
    }

    if constexpr (modes::isFlagSet(mode, modes::Mode::kQA)) {
      products.producedTrackDCAs(track.dcaXY(), track.dcaZ());
      products.producedTrackExtras(track.isPVContributor(),
                                   track.itsNCls(),
                                   track.itsNClsInnerBarrel(),
                                   track.itsChi2NCl(),
                                   track.itsClusterSizes(),
                                   track.tpcSignal(),
                                   track.tpcInnerParam(),
                                   track.tpcNClsFound(),
                                   track.tpcNClsCrossedRows(),
                                   track.tpcNClsShared(),
                                   track.beta());

      if constexpr (type == modes::Track::kPrimaryTrack) {
        products.producedTrackPids(track.itsNSigmaEl(),
                                   track.itsNSigmaPi(),
                                   track.itsNSigmaKa(),
                                   track.itsNSigmaPr(),
                                   track.itsNSigmaDe(),
                                   track.itsNSigmaTr(),
                                   track.itsNSigmaHe(),
                                   track.tpcNSigmaEl(),
                                   track.tpcNSigmaPi(),
                                   track.tpcNSigmaKa(),
                                   track.tpcNSigmaPr(),
                                   track.tpcNSigmaDe(),
                                   track.tpcNSigmaTr(),
                                   track.tpcNSigmaHe(),
                                   track.tofNSigmaEl(),
                                   track.tofNSigmaPi(),
                                   track.tofNSigmaKa(),
                                   track.tofNSigmaPr(),
                                   track.tofNSigmaDe(),
                                   track.tofNSigmaTr(),
                                   track.tofNSigmaHe());
      } else {
        products.producedTrackPids(0,
                                   0,
                                   0,
                                   0,
                                   0,
                                   0,
                                   0,
                                   track.tpcNSigmaEl(),
                                   track.tpcNSigmaPi(),
                                   track.tpcNSigmaKa(),
                                   track.tpcNSigmaPr(),
                                   track.tpcNSigmaDe(),
                                   track.tpcNSigmaTr(),
                                   track.tpcNSigmaHe(),
                                   track.tofNSigmaEl(),
                                   track.tofNSigmaPi(),
                                   track.tofNSigmaKa(),
                                   track.tofNSigmaPr(),
                                   track.tofNSigmaDe(),
                                   track.tofNSigmaTr(),
                                   track.tofNSigmaHe());
      }
    }
    indexMapTracks.emplace(track.globalIndex(), products.producedTracks.lastIndex());
  }

  template <modes::Mode mode, typename T>
  void fillTracks(T const& tracks)
  {
    for (const auto& track : tracks) {
      if (!trackSel.hasTofAboveThreshold(track)) {
        continue;
      }
      trackSel.applySelections(track);
      if (!trackSel.passesAllRequiredSelections()) {
        continue;
      }
      fillTrack<mode, modes::Track::kPrimaryTrack>(track);
    }
  }

  template <modes::Mode mode, typename T>
  void fillLambda(T const& v0, float sign, int posDaughterIndex, int negDaughterIndex)
  {
    float mass, massAnti;
    o2::aod::femtodatatypes::V0MaskType mask;
    if (sign > 0.f) {
      mass = v0.mLambda();
      massAnti = v0.mAntiLambda();
      mask = lambdaSel.getBitmask();
    } else {
      mass = v0.mAntiLambda();
      massAnti = v0.mLambda();
      mask = antiLambdaSel.getBitmask();
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kANALYSIS)) {
      products.producedLambdas(products.producedCollision.lastIndex(),
                               sign * v0.pt(),
                               v0.eta(),
                               v0.phi(),
                               mass,
                               posDaughterIndex,
                               negDaughterIndex);
      products.producedLambdaMasks(mask);
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQA)) {
      products.producedLambdaExtras(
        massAnti,
        v0.mK0Short(),
        v0.v0cosPA(),
        v0.dcaV0daughters(),
        v0.v0radius(),
        v0.x(),
        v0.y(),
        v0.z());
    }
  }

  template <modes::Mode mode, typename T>
  void fillK0short(T const& v0, int posDaughterIndex, int negDaughterIndex)
  {
    if constexpr (modes::isFlagSet(mode, modes::Mode::kANALYSIS)) {
      products.producedK0shorts(products.producedCollision.lastIndex(),
                                v0.pt(),
                                v0.eta(),
                                v0.phi(),
                                v0.mK0Short(),
                                posDaughterIndex,
                                negDaughterIndex);
      products.producedK0shortMasks(k0shortSel.getBitmask());
    }
    if constexpr (modes::isFlagSet(mode, modes::Mode::kQA)) {
      products.producedK0shortExtras(
        v0.mLambda(),
        v0.mAntiLambda(),
        v0.v0cosPA(),
        v0.dcaV0daughters(),
        v0.v0radius(),
        v0.x(),
        v0.y(),
        v0.z());
    }
  }

  template <modes::Mode mode, typename T1, typename T2>
  void fillV0s(T1 const& v0s, T2 const& fullTracks)
  {
    int64_t posDaughterIndex = 0;
    int64_t negDaughterIndex = 0;
    for (const auto& v0 : v0s) {
      if (v0.pt() < confV0Filters.ptMin.value || v0.pt() > confV0Filters.ptMax.value ||
          v0.eta() < confV0Filters.etaMin.value || v0.eta() > confV0Filters.etaMax.value ||
          v0.phi() < confV0Filters.phiMin.value || v0.phi() > confV0Filters.phiMax.value) {
        continue;
      }
      auto posDaughter = v0.template posTrack_as<T2>();
      auto negDaughter = v0.template negTrack_as<T2>();
      if (ConfOptions.produceLambda.value) {
        lambdaSel.applySelections(v0, fullTracks);
        if (lambdaSel.passesAllRequiredSelections() && lambdaSel.checkHypothesis(v0)) {
          posDaughterIndex = getDaughterIndex<mode, modes::Track::kV0Daughter>(posDaughter);
          negDaughterIndex = getDaughterIndex<mode, modes::Track::kV0Daughter>(negDaughter);
          fillLambda<mode>(v0, 1.f, posDaughterIndex, negDaughterIndex);
        }
        antiLambdaSel.applySelections(v0, fullTracks);
        if (antiLambdaSel.passesAllRequiredSelections() && antiLambdaSel.checkHypothesis(v0)) {
          posDaughterIndex = getDaughterIndex<mode, modes::Track::kV0Daughter>(posDaughter);
          negDaughterIndex = getDaughterIndex<mode, modes::Track::kV0Daughter>(negDaughter);
          fillLambda<mode>(v0, -1.f, posDaughterIndex, negDaughterIndex);
        }
      }
      if (ConfOptions.produceK0short.value) {
        k0shortSel.applySelections(v0, fullTracks);
        if (k0shortSel.passesAllRequiredSelections() && k0shortSel.checkHypothesis(v0)) {
          posDaughterIndex = getDaughterIndex<mode, modes::Track::kV0Daughter>(posDaughter);
          negDaughterIndex = getDaughterIndex<mode, modes::Track::kV0Daughter>(negDaughter);
          fillK0short<mode>(v0, posDaughterIndex, negDaughterIndex);
        }
      }
    }
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3>
  void fillCascades(T1 const& fullCascades, T2 const& fullTracks, T3 const& col)
  {
    int64_t bachelorIndex = 0;
    int64_t posDaughterIndex = 0;
    int64_t negDaughterIndex = 0;
    for (const auto& cascade : fullCascades) {
      if (cascade.pt() < confCascadeFilters.ptMin.value || cascade.pt() > confCascadeFilters.ptMax.value ||
          cascade.eta() < confCascadeFilters.etaMin.value || cascade.eta() > confCascadeFilters.etaMax.value ||
          cascade.phi() < confCascadeFilters.phiMin.value || cascade.phi() > confCascadeFilters.phiMax.value ||
          cascade.mLambda() < confCascadeFilters.massLambdaMin.value || cascade.mLambda() > confCascadeFilters.massLambdaMax.value) {
        continue;
      }
      auto bachelor = cascade.template bachelor_as<T2>();
      auto posDaughter = cascade.template posTrack_as<T2>();
      auto negDaughter = cascade.template negTrack_as<T2>();

      if (ConfOptions.produceXi.value) {
        xiSel.applySelections(cascade, fullTracks, col);
        if (xiSel.passesAllRequiredSelections() && xiSel.checkHypothesis(cascade)) {
          bachelorIndex = getDaughterIndex<mode, modes::Track::kCascadeBachelor>(bachelor);
          posDaughterIndex = getDaughterIndex<mode, modes::Track::kV0Daughter>(posDaughter);
          negDaughterIndex = getDaughterIndex<mode, modes::Track::kV0Daughter>(negDaughter);
          fillCascade<mode, o2::analysis::femtounited::modes::Cascade::kXi>(
            cascade, bachelorIndex, posDaughterIndex, negDaughterIndex, col);
        }
      }

      if (ConfOptions.produceOmega.value) {
        omegaSel.applySelections(cascade, fullTracks, col);
        if (omegaSel.passesAllRequiredSelections() && omegaSel.checkHypothesis(cascade)) {
          bachelorIndex = getDaughterIndex<mode, modes::Track::kCascadeBachelor>(bachelor);
          posDaughterIndex = getDaughterIndex<mode, modes::Track::kV0Daughter>(posDaughter);
          negDaughterIndex = getDaughterIndex<mode, modes::Track::kV0Daughter>(negDaughter);
          fillCascade<mode, o2::analysis::femtounited::modes::Cascade::kOmega>(
            cascade, bachelorIndex, posDaughterIndex, negDaughterIndex, col);
        }
      }
    }
  }

  template <modes::Mode mode, modes::Cascade C, typename T1, typename T2>
  void fillCascade(T1 const& cascade, int bachelorIndex, int posDaughterIndex, int negDaughterIndex, T2 const& col)
  {
    if constexpr (modes::isFlagSet(mode, modes::Mode::kANALYSIS)) {
      if constexpr (modes::isFlagSet(C, modes::Cascade::kXi)) {
        products.producedXis(products.producedCollision.lastIndex(),
                             cascade.sign() * cascade.pt(),
                             cascade.eta(),
                             cascade.phi(),
                             cascade.mXi(),
                             bachelorIndex,
                             posDaughterIndex,
                             negDaughterIndex);
        products.producedXiMasks(xiSel.getBitmask());
      }
      if constexpr (modes::isFlagSet(C, modes::Cascade::kOmega)) {
        products.producedOmegas(products.producedCollision.lastIndex(),
                                cascade.sign() * cascade.pt(),
                                cascade.eta(),
                                cascade.phi(),
                                cascade.mOmega(),
                                bachelorIndex,
                                posDaughterIndex,
                                negDaughterIndex);
        products.producedOmegaMasks(omegaSel.getBitmask());
      }
    }

    if constexpr (modes::isFlagSet(mode, modes::Mode::kQA)) {
      if constexpr (modes::isFlagSet(C, modes::Cascade::kXi)) {
        products.producedXiExtras(
          cascade.mOmega(),
          cascade.casccosPA(col.posX(), col.posY(), col.posZ()),
          cascade.dcacascdaughters(),
          cascade.cascradius(),
          cascade.v0cosPA(col.posX(), col.posY(), col.posZ()),
          cascade.dcaV0daughters(),
          cascade.v0radius(),
          cascade.dcav0topv(col.posY(), col.posY(), col.posZ()));
      }
      if constexpr (modes::isFlagSet(C, modes::Cascade::kOmega)) {
        products.producedOmegaExtras(
          cascade.mXi(),
          cascade.casccosPA(col.posX(), col.posY(), col.posZ()),
          cascade.dcacascdaughters(),
          cascade.cascradius(),
          cascade.v0cosPA(col.posX(), col.posY(), col.posZ()),
          cascade.dcaV0daughters(),
          cascade.v0radius(),
          cascade.dcav0topv(col.posY(), col.posY(), col.posZ()));
      }
    }
  }

  template <modes::Mode mode, modes::TwoTrackResonance reso, typename T>
  void fillResonance(T const& posDaughter, T const& negDaughter)
  {
    int64_t posDaughterIndex = 0;
    int64_t negDaughterIndex = 0;
    if constexpr (modes::isFlagSet(reso, modes::TwoTrackResonance::kRho0)) {
      if (!rho0Sels.hasTofAboveThreshold(posDaughter, negDaughter)) {
        return;
      }
      rho0Sels.applySelections(posDaughter, negDaughter);
      if (!rho0Sels.passesAllRequiredSelections()) {
        return;
      }
      rho0Sels.reconstructResonance(posDaughter, negDaughter);
      if (!rho0Sels.checkFilters() || !rho0Sels.checkHypothesis()) {
        return;
      }
      posDaughterIndex = getDaughterIndex<mode, modes::Track::kResonanceDaughter>(posDaughter);
      negDaughterIndex = getDaughterIndex<mode, modes::Track::kResonanceDaughter>(negDaughter);
      products.producedRhos(
        products.producedCollision.lastIndex(),
        rho0Sels.getPt(),
        rho0Sels.getEta(),
        rho0Sels.getPhi(),
        rho0Sels.getMass(),
        posDaughterIndex,
        negDaughterIndex,
        rho0Sels.getPosDauMomAboveThres(),
        rho0Sels.getNegDauMomAboveThres());
      products.producedRhoMasks(rho0Sels.getBitmask());
    }
    if constexpr (modes::isFlagSet(reso, modes::TwoTrackResonance::kPhi)) {
      if (!phiSels.hasTofAboveThreshold(posDaughter, negDaughter)) {
        return;
      }
      phiSels.applySelections(posDaughter, negDaughter);
      if (!phiSels.passesAllRequiredSelections()) {
        return;
      }
      phiSels.reconstructResonance(posDaughter, negDaughter);
      if (!phiSels.checkFilters() || !phiSels.checkHypothesis()) {
        return;
      }
      posDaughterIndex = getDaughterIndex<mode, modes::Track::kResonanceDaughter>(posDaughter);
      negDaughterIndex = getDaughterIndex<mode, modes::Track::kResonanceDaughter>(negDaughter);
      products.producedPhis(
        products.producedCollision.lastIndex(),
        phiSels.getPt(),
        phiSels.getEta(),
        phiSels.getPhi(),
        phiSels.getMass(),
        posDaughterIndex,
        negDaughterIndex,
        phiSels.getPosDauMomAboveThres(),
        phiSels.getNegDauMomAboveThres());
      products.producedPhiMasks(phiSels.getBitmask());
    }
    if constexpr (modes::isFlagSet(reso, modes::TwoTrackResonance::kKstar0)) {
      if (!kstar0Sels.hasTofAboveThreshold(posDaughter, negDaughter)) {
        return;
      }
      kstar0Sels.applySelections(posDaughter, negDaughter);
      if (!kstar0Sels.passesAllRequiredSelections()) {
        return;
      }
      kstar0Sels.reconstructResonance(posDaughter, negDaughter);
      if (kstar0Sels.checkFilters() && kstar0Sels.checkHypothesis()) {
        posDaughterIndex = getDaughterIndex<mode, modes::Track::kResonanceDaughter>(posDaughter);
        negDaughterIndex = getDaughterIndex<mode, modes::Track::kResonanceDaughter>(negDaughter);
        products.producedKstars(
          products.producedCollision.lastIndex(),
          kstar0Sels.getPt(),
          kstar0Sels.getEta(),
          kstar0Sels.getPhi(),
          kstar0Sels.getMass(),
          posDaughterIndex,
          negDaughterIndex,
          kstar0Sels.getPosDauMomAboveThres(),
          kstar0Sels.getNegDauMomAboveThres());
        products.producedKstarMasks(kstar0Sels.getBitmask());
      }
    }
    if constexpr (modes::isFlagSet(reso, modes::TwoTrackResonance::kKstarBar0)) {
      // try kstar0
      if (!kstarBar0Sels.hasTofAboveThreshold(posDaughter, negDaughter)) {
        return;
      }
      kstarBar0Sels.applySelections(posDaughter, negDaughter);
      if (!kstarBar0Sels.passesAllRequiredSelections()) {
        return;
      }
      kstarBar0Sels.reconstructResonance(posDaughter, negDaughter);
      if (kstarBar0Sels.checkFilters() && kstarBar0Sels.checkHypothesis()) {
        posDaughterIndex = getDaughterIndex<mode, modes::Track::kResonanceDaughter>(posDaughter);
        negDaughterIndex = getDaughterIndex<mode, modes::Track::kResonanceDaughter>(negDaughter);
        products.producedKstars(
          products.producedCollision.lastIndex(),
          -1.f * kstarBar0Sels.getPt(),
          kstarBar0Sels.getEta(),
          kstarBar0Sels.getPhi(),
          kstarBar0Sels.getMass(),
          posDaughterIndex,
          negDaughterIndex,
          kstarBar0Sels.getPosDauMomAboveThres(),
          kstarBar0Sels.getNegDauMomAboveThres());
        products.producedKstarMasks(kstarBar0Sels.getBitmask());
      }
    }
  }

  template <modes::Mode mode, typename T1, typename T2>
  void fillResonances(T1 const& col, T2 const& /*tracks*/)
  {
    auto groupPositiveTracks = partitionPositiveDaughters->sliceByCached(track::collisionId, col.globalIndex(), cache);
    auto groupNegativeTracks = partitionNegativeDaughters->sliceByCached(track::collisionId, col.globalIndex(), cache);
    for (auto const& positiveTrack : groupPositiveTracks) {
      for (auto const& negativeTrack : groupNegativeTracks) {
        if (ConfOptions.produceRho0.value) {
          fillResonance<mode, modes::TwoTrackResonance::kRho0>(positiveTrack, negativeTrack);
        }
        if (ConfOptions.producePhi.value) {
          fillResonance<mode, modes::TwoTrackResonance::kPhi>(positiveTrack, negativeTrack);
        }
        if (ConfOptions.produceKstar0.value) {
          fillResonance<mode, modes::TwoTrackResonance::kKstar0>(positiveTrack, negativeTrack);
        }
      }
    }
  }

  template <modes::System system, modes::Mode mode, typename T1, typename T2, typename T3>
  void processTracks(T1 const& col, T2 const& /* fullBcs*/, T3 const& fullTracks)
  {
    initFromCcdb(col.template bc_as<T2>());
    collisionSel.setMagneticField(magField);
    collisionSel.setSphericity(fullTracks);
    if (!collisionSel.checkCuts<modes::System::kPP_Run3>(col)) {
      return;
    }
    fillCollision<system>(col);
    indexMapTracks.clear();
    fillTracks<mode>(fullTracks);
    if (produceResonances) {
      fillResonances<mode>(col, fullTracks);
    }
  }

  template <modes::System system, modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5>
  void processTracksV0s(T1 const& col, T2 const& /* fullBcs*/, T3 const& fullTracks, T4 const& fullTracksWithItsPid, T5 const& fullV0s)
  {
    initFromCcdb(col.template bc_as<T2>());
    collisionSel.setMagneticField(magField);
    collisionSel.setSphericity(fullTracks);
    if (!collisionSel.checkCuts<modes::System::kPP_Run3>(col)) {
      return;
    }
    fillCollision<system>(col);
    indexMapTracks.clear();
    fillTracks<mode>(fullTracksWithItsPid);
    if (produceResonances) {
      fillResonances<mode>(col, fullTracksWithItsPid);
    }
    if (produceV0s) {
      fillV0s<mode>(fullV0s, fullTracks);
    }
  }

  template <modes::System system, modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void processTracksV0sCascades(T1 const& col, T2 const& /* fullBcs*/, T3 const& fullTracks, T4 const& fullTracksWithItsPid, T5 const& fullV0s, T6 const& fullCascades)
  {
    initFromCcdb(col.template bc_as<T2>());
    collisionSel.setMagneticField(magField);
    collisionSel.setSphericity(fullTracksWithItsPid);
    if (!collisionSel.checkCuts<modes::System::kPP_Run3>(col)) {
      return;
    }
    fillCollision<system>(col);
    indexMapTracks.clear();
    fillTracks<mode>(fullTracksWithItsPid);
    if (produceResonances) {
      fillResonances<mode>(col, fullTracksWithItsPid);
    }
    if (produceV0s) {
      fillV0s<mode>(fullV0s, fullTracks);
    }
    if (produceCascades) {
      fillCascades<mode>(fullCascades, fullTracks, col);
    }
  }

  // proccess tracks
  void processTracksRun3pp(Filtered<consumeddata::Run3PpCollisions>::iterator const& col,
                           BCsWithTimestamps const& bcs,
                           Filtered<consumeddata::Run3FullPidTracks> const& tracks)
  {
    // its pid information is generated dynamically, so we need to add it here
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3FullPidTracks, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                            pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);
    if (ConfOptions.produceExtraTables.value) {
      processTracks<modes::System::kPP_Run3, modes::Mode::kANALYSIS_QA>(col, bcs, tracksWithItsPid);
    } else {
      processTracks<modes::System::kPP_Run3, modes::Mode::kANALYSIS>(col, bcs, tracksWithItsPid);
    }
  }
  PROCESS_SWITCH(FemtoUnitedProducer, processTracksRun3pp, "Process tracks", true);

  // process tracks and v0s
  void processTracksV0sRun3pp(Filtered<consumeddata::Run3PpCollisions>::iterator const& col,
                              BCsWithTimestamps const& bcs,
                              Filtered<consumeddata::Run3FullPidTracks> const& tracks,
                              consumeddata::Run3PpVzeros const& v0s)
  {
    // its pid information is generated dynamically, so we need to add it here
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3FullPidTracks, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                            pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);
    if (ConfOptions.produceExtraTables.value) {
      processTracksV0s<modes::System::kPP_Run3, modes::Mode::kANALYSIS_QA>(col, bcs, tracks, tracksWithItsPid, v0s);
    } else {
      processTracksV0s<modes::System::kPP_Run3, modes::Mode::kANALYSIS>(col, bcs, tracks, tracksWithItsPid, v0s);
    }
  };
  PROCESS_SWITCH(FemtoUnitedProducer, processTracksV0sRun3pp, "Process tracks and v0s", false);

  // process tracks, v0s and casacades
  void processTracksV0sCascadesRun3pp(Filtered<consumeddata::Run3PpCollisions>::iterator const& col,
                                      BCsWithTimestamps const& bcs,
                                      Filtered<consumeddata::Run3FullPidTracks> const& tracks,
                                      consumeddata::Run3PpVzeros const& v0s,
                                      consumeddata::Run3PpCascades const& cascades)
  {
    // its pid information is generated dynamically, so we need to add it here
    auto tracksWithItsPid = o2::soa::Attach<consumeddata::Run3FullPidTracks, pidits::ITSNSigmaEl, pidits::ITSNSigmaPi,
                                            pidits::ITSNSigmaKa, pidits::ITSNSigmaPr, pidits::ITSNSigmaDe, pidits::ITSNSigmaTr, pidits::ITSNSigmaHe>(tracks);
    if (ConfOptions.produceExtraTables.value) {
      processTracksV0sCascades<modes::System::kPP_Run3, modes::Mode::kANALYSIS_QA>(col, bcs, tracks, tracksWithItsPid, v0s, cascades);
    } else {
      processTracksV0sCascades<modes::System::kPP_Run3, modes::Mode::kANALYSIS>(col, bcs, tracks, tracksWithItsPid, v0s, cascades);
    }
  }
  PROCESS_SWITCH(FemtoUnitedProducer, processTracksV0sCascadesRun3pp, "Provide Tracks, V0s and Cascades for Run3", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<FemtoUnitedProducer>(cfgc)};
  return workflow;
}
