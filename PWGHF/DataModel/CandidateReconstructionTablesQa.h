// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CandidateReconstructionTablesQa.h
/// \brief QA candidate tables for Xic0/Omegac0 cascade decays
///
/// \author Fabio Catalano <fabio.catalano@cern.ch>, University of Houston

#ifndef PWGHF_DATAMODEL_CANDIDATERECONSTRUCTIONTABLESQA_H_
#define PWGHF_DATAMODEL_CANDIDATERECONSTRUCTIONTABLESQA_H_

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>

#include <array>
#include <cstdint>

namespace o2::aod
{
namespace hf_cand_casc_bach_qa
{
DECLARE_SOA_COLUMN(XCascadeVertex, xCascadeVertex, float);
DECLARE_SOA_COLUMN(YCascadeVertex, yCascadeVertex, float);
DECLARE_SOA_COLUMN(ZCascadeVertex, zCascadeVertex, float);
DECLARE_SOA_COLUMN(XV0Vertex, xV0Vertex, float);
DECLARE_SOA_COLUMN(YV0Vertex, yV0Vertex, float);
DECLARE_SOA_COLUMN(ZV0Vertex, zV0Vertex, float);
DECLARE_SOA_COLUMN(MassCascade, massCascade, float);
DECLARE_SOA_COLUMN(MassLambda, massLambda, float);
DECLARE_SOA_COLUMN(EtaCascade, etaCascade, float);
DECLARE_SOA_COLUMN(EtaV0, etaV0, float);
DECLARE_SOA_COLUMN(CpaCascade, cpaCascade, float);
DECLARE_SOA_COLUMN(CpaXYCascade, cpaXYCascade, float);
DECLARE_SOA_COLUMN(CpaV0, cpaV0, float);
DECLARE_SOA_COLUMN(CpaXYV0, cpaXYV0, float);
DECLARE_SOA_COLUMN(CpaCascadeToCharmHadron, cpaCascadeToCharmHadron, float);
DECLARE_SOA_COLUMN(CpaXYCascadeToCharmHadron, cpaXYCascadeToCharmHadron, float);
DECLARE_SOA_COLUMN(CpaV0ToCascade, cpaV0ToCascade, float);
DECLARE_SOA_COLUMN(CpaXYV0ToCascade, cpaXYV0ToCascade, float);
DECLARE_SOA_COLUMN(ImpactParameterXYCascadeBach, impactParameterXYCascadeBach, float);
DECLARE_SOA_COLUMN(ImpactParameterZCascadeBach, impactParameterZCascadeBach, float);
DECLARE_SOA_COLUMN(ImpactParameterXYV0PosDau, impactParameterXYV0PosDau, float);
DECLARE_SOA_COLUMN(ImpactParameterZV0PosDau, impactParameterZV0PosDau, float);
DECLARE_SOA_COLUMN(ImpactParameterXYV0NegDau, impactParameterXYV0NegDau, float);
DECLARE_SOA_COLUMN(ImpactParameterZV0NegDau, impactParameterZV0NegDau, float);
DECLARE_SOA_COLUMN(DcaV0Daughters, dcaV0Daughters, float);
DECLARE_SOA_COLUMN(DcaCascadeDaughters, dcaCascadeDaughters, float);
DECLARE_SOA_COLUMN(DcaDaughters, dcaDaughters, float);

DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthCascade, decayLengthCascade, //! Cascade decay length from the parent vertex
                           [](float pvx, float pvy, float pvz, float svx, float svy, float svz) -> float {
                             return RecoDecay::distance(std::array{pvx, pvy, pvz}, std::array{svx, svy, svz});
                           });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthXYCascade, decayLengthXYCascade, //! Cascade transverse decay length from the parent vertex
                           [](float pvx, float pvy, float svx, float svy) -> float {
                             return RecoDecay::distanceXY(std::array{pvx, pvy}, std::array{svx, svy});
                           });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthV0, decayLengthV0, //! Lambda decay length from the cascade vertex
                           [](float pvx, float pvy, float pvz, float svx, float svy, float svz) -> float {
                             return RecoDecay::distance(std::array{pvx, pvy, pvz}, std::array{svx, svy, svz});
                           });
DECLARE_SOA_DYNAMIC_COLUMN(DecayLengthXYV0, decayLengthXYV0, //! Lambda transverse decay length from the cascade vertex
                           [](float pvx, float pvy, float svx, float svy) -> float {
                             return RecoDecay::distanceXY(std::array{pvx, pvy}, std::array{svx, svy});
                           });

// KFParticle specific columns
DECLARE_SOA_COLUMN(Chi2Geo, chi2Geo, float);
DECLARE_SOA_COLUMN(Chi2GeoCascade, chi2GeoCascade, float);
DECLARE_SOA_COLUMN(Chi2GeoV0, chi2GeoV0, float);
DECLARE_SOA_COLUMN(NdfGeo, ndfGeo, int32_t);
DECLARE_SOA_COLUMN(NdfGeoCascade, ndfGeoCascade, int32_t);
DECLARE_SOA_COLUMN(NdfGeoV0, ndfGeoV0, int32_t);
DECLARE_SOA_COLUMN(DecayLengthNormKf, decayLengthNormKf, float);
DECLARE_SOA_COLUMN(DecayLengthNormKfCascade, decayLengthNormKfCascade, float);
DECLARE_SOA_COLUMN(DecayLengthNormKfV0, decayLengthNormKfV0, float);
DECLARE_SOA_COLUMN(Chi2TopoToPv, chi2TopoToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoCascadeToPv, chi2TopoCascadeToPv, float);
DECLARE_SOA_COLUMN(Chi2TopoV0ToPv, chi2TopoV0ToPv, float);
DECLARE_SOA_COLUMN(NdfTopoToPv, ndfTopoToPv, int32_t);
DECLARE_SOA_COLUMN(NdfTopoCascadeToPv, ndfTopoCascadeToPv, int32_t);
DECLARE_SOA_COLUMN(NdfTopoV0ToPv, ndfTopoV0ToPv, int32_t);
DECLARE_SOA_COLUMN(Chi2TopoCascadeToCharmHadron, chi2TopoCascadeToCharmHadron, float);
DECLARE_SOA_COLUMN(Chi2TopoV0ToCascade, chi2TopoV0ToCascade, float);
DECLARE_SOA_COLUMN(NdfTopoCascadeToCharmHadron, ndfTopoCascadeToCharmHadron, int32_t);
DECLARE_SOA_COLUMN(NdfTopoV0ToCascade, ndfTopoV0ToCascade, int32_t);

DECLARE_SOA_DYNAMIC_COLUMN(Chi2GeoOverNdf, chi2GeoOverNdf,
                           [](float chi2, int32_t ndf) -> float { return chi2 / ndf; });
DECLARE_SOA_DYNAMIC_COLUMN(Chi2GeoOverNdfCascade, chi2GeoOverNdfCascade,
                           [](float chi2, int32_t ndf) -> float { return chi2 / ndf; });
DECLARE_SOA_DYNAMIC_COLUMN(Chi2GeoOverNdfV0, chi2GeoOverNdfV0,
                           [](float chi2, int32_t ndf) -> float { return chi2 / ndf; });
DECLARE_SOA_DYNAMIC_COLUMN(Chi2TopoOverNdfToPv, chi2TopoOverNdfToPv,
                           [](float chi2, int32_t ndf) -> float { return chi2 / ndf; });
DECLARE_SOA_DYNAMIC_COLUMN(Chi2TopoOverNdfCascadeToPv, chi2TopoOverNdfCascadeToPv,
                           [](float chi2, int32_t ndf) -> float { return chi2 / ndf; });
DECLARE_SOA_DYNAMIC_COLUMN(Chi2TopoOverNdfV0ToPv, chi2TopoOverNdfV0ToPv,
                           [](float chi2, int32_t ndf) -> float { return chi2 / ndf; });
DECLARE_SOA_DYNAMIC_COLUMN(Chi2TopoOverNdfCascadeToCharmHadron, chi2TopoOverNdfCascadeToCharmHadron,
                           [](float chi2, int32_t ndf) -> float { return chi2 / ndf; });
DECLARE_SOA_DYNAMIC_COLUMN(Chi2TopoOverNdfV0ToCascade, chi2TopoOverNdfV0ToCascade,
                           [](float chi2, int32_t ndf) -> float { return chi2 / ndf; });

} // namespace hf_cand_casc_bach_qa

// Common QA candidate table for XiPi, OmegaPi and OmegaK
// Convention: Prong0 is the Cascade, Prong1 is the bachelor track
DECLARE_SOA_TABLE(HfCandCascBachBaseQa, "AOD", "HFCASCBACHBASE",
                  o2::soa::Index<>,
                  hf_cand::CollisionId,
                  collision::PosX, collision::PosY, collision::PosZ,
                  hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex,
                  hf_cand_casc_bach_qa::XCascadeVertex, hf_cand_casc_bach_qa::YCascadeVertex, hf_cand_casc_bach_qa::ZCascadeVertex,
                  hf_cand_casc_bach_qa::XV0Vertex, hf_cand_casc_bach_qa::YV0Vertex, hf_cand_casc_bach_qa::ZV0Vertex,
                  hf_cand::ErrorDecayLength, hf_cand::ErrorDecayLengthXY,
                  hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0,
                  hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1,
                  hf_cand::ImpactParameter0, hf_cand::ImpactParameter1,
                  hf_cand::ErrorImpactParameter0, hf_cand::ErrorImpactParameter1,
                  hf_cand::ImpactParameterZ0, hf_cand::ImpactParameterZ1,
                  hf_cand::ErrorImpactParameterZ0, hf_cand::ErrorImpactParameterZ1,
                  hf_cand_casc_bach_qa::ImpactParameterXYV0PosDau, hf_cand_casc_bach_qa::ImpactParameterXYV0NegDau, hf_cand_casc_bach_qa::ImpactParameterXYCascadeBach,
                  hf_cand_casc_bach_qa::ImpactParameterZV0PosDau, hf_cand_casc_bach_qa::ImpactParameterZV0NegDau, hf_cand_casc_bach_qa::ImpactParameterZCascadeBach,
                  hf_cand_casc_bach_qa::DcaDaughters, hf_cand_casc_bach_qa::DcaCascadeDaughters, hf_cand_casc_bach_qa::DcaV0Daughters,
                  hf_cand_casc_bach_qa::CpaCascade, hf_cand_casc_bach_qa::CpaXYCascade,
                  hf_cand_casc_bach_qa::CpaV0, hf_cand_casc_bach_qa::CpaXYV0,
                  hf_cand_casc_bach_qa::CpaCascadeToCharmHadron, hf_cand_casc_bach_qa::CpaXYCascadeToCharmHadron,
                  hf_cand_casc_bach_qa::CpaV0ToCascade, hf_cand_casc_bach_qa::CpaXYV0ToCascade,
                  hf_cand_casc_bach_qa::MassCascade, hf_cand_casc_bach_qa::MassLambda,
                  hf_cand_casc_bach_qa::EtaCascade, hf_cand_casc_bach_qa::EtaV0,
                  hf_track_index::Prong1Id, cascdata::BachelorId, cascdata::PosTrackId, cascdata::NegTrackId,
                  hf_track_index::HFflag,
                  /* dynamic columns */
                  hf_cand_2prong::M<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand::Pt<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::P<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::PVector<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Eta<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Phi<hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  hf_cand::Y<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::E<hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Ct<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand_2prong::ImpactParameterProduct<hf_cand::ImpactParameter0, hf_cand::ImpactParameter1>,
                  hf_cand_2prong::CosThetaStar<hf_cand::PxProng0, hf_cand::PyProng0, hf_cand::PzProng0, hf_cand::PxProng1, hf_cand::PyProng1, hf_cand::PzProng1>,
                  hf_cand::DecayLength<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex>,
                  hf_cand::DecayLengthXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex>,
                  hf_cand::DecayLengthNormalised<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand::ErrorDecayLength>,
                  hf_cand::DecayLengthXYNormalised<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ErrorDecayLengthXY>,
                  hf_cand::ImpactParameterXY<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::Cpa<collision::PosX, collision::PosY, collision::PosZ, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz>,
                  hf_cand::CpaXY<collision::PosX, collision::PosY, hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand_2prong::Px, hf_cand_2prong::Py>,
                  // TODO: investigate if these decay lenghts should not be computed w.r.t. the primary vertex
                  hf_cand_casc_bach_qa::DecayLengthCascade<collision::PosX, collision::PosY, collision::PosZ, hf_cand_casc_bach_qa::XCascadeVertex, hf_cand_casc_bach_qa::YCascadeVertex, hf_cand_casc_bach_qa::ZCascadeVertex>,
                  hf_cand_casc_bach_qa::DecayLengthXYCascade<collision::PosX, collision::PosY, hf_cand_casc_bach_qa::XCascadeVertex, hf_cand_casc_bach_qa::YCascadeVertex>,
                  hf_cand_casc_bach_qa::DecayLengthV0<collision::PosX, collision::PosY, collision::PosZ, hf_cand_casc_bach_qa::XV0Vertex, hf_cand_casc_bach_qa::YV0Vertex, hf_cand_casc_bach_qa::ZV0Vertex>,
                  hf_cand_casc_bach_qa::DecayLengthXYV0<collision::PosX, collision::PosY, hf_cand_casc_bach_qa::XV0Vertex, hf_cand_casc_bach_qa::YV0Vertex>,
                  hf_cand::ImpactParameterNormalised0<hf_cand::ImpactParameter0, hf_cand::ErrorImpactParameter0>,
                  hf_cand::ImpactParameterNormalised1<hf_cand::ImpactParameter1, hf_cand::ErrorImpactParameter1>,
                  hf_cand::ImpactParameterZNormalised0<hf_cand::ImpactParameterZ0, hf_cand::ErrorImpactParameterZ0>,
                  hf_cand::ImpactParameterZNormalised1<hf_cand::ImpactParameterZ1, hf_cand::ErrorImpactParameterZ1>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandCascBachExtQa, HfCandCascBachBaseQa, "HFCASCBACHEXT",
                                hf_cand_2prong::Px, hf_cand_2prong::Py, hf_cand_2prong::Pz);

using HfCandCascBachQa = HfCandCascBachExtQa;

// KF specific information for HfCandCascBach
DECLARE_SOA_TABLE(HfCandCascBachKfQa, "AOD", "HFCASCBACHKF",
                  hf_cand_casc_bach_qa::Chi2Geo, hf_cand_casc_bach_qa::Chi2GeoCascade, hf_cand_casc_bach_qa::Chi2GeoV0,
                  hf_cand_casc_bach_qa::NdfGeo, hf_cand_casc_bach_qa::NdfGeoCascade, hf_cand_casc_bach_qa::NdfGeoV0,
                  hf_cand_casc_bach_qa::DecayLengthNormKf, hf_cand_casc_bach_qa::DecayLengthNormKfCascade, hf_cand_casc_bach_qa::DecayLengthNormKfV0,
                  hf_cand_casc_bach_qa::Chi2TopoToPv, hf_cand_casc_bach_qa::Chi2TopoCascadeToPv, hf_cand_casc_bach_qa::Chi2TopoV0ToPv,
                  hf_cand_casc_bach_qa::NdfTopoToPv, hf_cand_casc_bach_qa::NdfTopoCascadeToPv, hf_cand_casc_bach_qa::NdfTopoV0ToPv,
                  hf_cand_casc_bach_qa::Chi2TopoCascadeToCharmHadron, hf_cand_casc_bach_qa::Chi2TopoV0ToCascade,
                  hf_cand_casc_bach_qa::NdfTopoCascadeToCharmHadron, hf_cand_casc_bach_qa::NdfTopoV0ToCascade,
                  /* dynamic columns */
                  hf_cand_casc_bach_qa::Chi2GeoOverNdf<hf_cand_casc_bach_qa::Chi2Geo, hf_cand_casc_bach_qa::NdfGeo>,
                  hf_cand_casc_bach_qa::Chi2GeoOverNdfCascade<hf_cand_casc_bach_qa::Chi2GeoCascade, hf_cand_casc_bach_qa::NdfGeoCascade>,
                  hf_cand_casc_bach_qa::Chi2GeoOverNdfV0<hf_cand_casc_bach_qa::Chi2GeoV0, hf_cand_casc_bach_qa::NdfGeoV0>,
                  hf_cand_casc_bach_qa::Chi2TopoOverNdfToPv<hf_cand_casc_bach_qa::Chi2TopoToPv, hf_cand_casc_bach_qa::NdfTopoToPv>,
                  hf_cand_casc_bach_qa::Chi2TopoOverNdfCascadeToPv<hf_cand_casc_bach_qa::Chi2TopoCascadeToPv, hf_cand_casc_bach_qa::NdfTopoCascadeToPv>,
                  hf_cand_casc_bach_qa::Chi2TopoOverNdfV0ToPv<hf_cand_casc_bach_qa::Chi2TopoV0ToPv, hf_cand_casc_bach_qa::NdfTopoV0ToPv>,
                  hf_cand_casc_bach_qa::Chi2TopoOverNdfCascadeToCharmHadron<hf_cand_casc_bach_qa::Chi2TopoCascadeToCharmHadron, hf_cand_casc_bach_qa::NdfTopoCascadeToCharmHadron>,
                  hf_cand_casc_bach_qa::Chi2TopoOverNdfV0ToCascade<hf_cand_casc_bach_qa::Chi2TopoV0ToCascade, hf_cand_casc_bach_qa::NdfTopoV0ToCascade>);

// table with results of reconstruction level MC matching
DECLARE_SOA_TABLE(HfCandCascBachMCRecQa, "AOD", "HFCASCBACHMCREC", //!
                  hf_cand_mc_flag::FlagMcMatchRec,
                  hf_cand_mc_flag::OriginMcRec,
                  hf_cand_mc_flag::FlagMcDecayChanRec,
                  hf_cand_mc_flag::PtBhadMotherPart,
                  hf_cand_mc_flag::PdgBhadMotherPart,
                  hf_cand_mc_flag::NTracksDecayed,
                  hf_cand_mc_flag::NInteractionsWithMaterial);

// table with results of generator level MC matching
DECLARE_SOA_TABLE(HfCandCascBachMCGenQa, "AOD", "HFCASCBACHMCGEN", //!
                  hf_cand_mc_flag::FlagMcMatchGen,
                  hf_cand_mc_flag::OriginMcGen,
                  hf_cand_mc_flag::FlagMcDecayChanGen,
                  hf_cand_mc_flag::IdxBhadMotherPart);

} // namespace o2::aod

struct HfHelperQa {
  template <typename T>
  static auto invMassXiczeroToXiPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassXiMinus, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  static auto invMassOmegaczeroToXiPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassXiMinus, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  static auto invMassOmegaczeroToOmegaPi(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassOmegaMinus, o2::constants::physics::MassPiPlus});
  }

  template <typename T>
  static auto invMassOmegaczeroToOmegaK(const T& candidate)
  {
    return candidate.m(std::array{o2::constants::physics::MassOmegaMinus, o2::constants::physics::MassKPlus});
  }

  template <typename T>
  static auto yXiczero(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassXiC0);
  }

  template <typename T>
  static auto eXiczero(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassXiC0);
  }

  template <typename T>
  static auto ctXiczero(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassXiC0);
  }

  template <typename T>
  static auto cosThetaStarXiczeroToXiPi(const T& candidate)
  {
    return candidate.cosThetaStar(std::array{o2::constants::physics::MassXiMinus, o2::constants::physics::MassPiPlus}, o2::constants::physics::MassXiC0, 1);
  }

  template <typename T>
  static auto yOmegaczero(const T& candidate)
  {
    return candidate.y(o2::constants::physics::MassOmegaC0);
  }

  template <typename T>
  static auto eOmegaczero(const T& candidate)
  {
    return candidate.e(o2::constants::physics::MassOmegaC0);
  }

  template <typename T>
  static auto ctOmegaczero(const T& candidate)
  {
    return candidate.ct(o2::constants::physics::MassOmegaC0);
  }

  template <typename T>
  static auto cosThetaStarOmegaczeroToXiPi(const T& candidate)
  {
    return candidate.cosThetaStar(std::array{o2::constants::physics::MassXiMinus, o2::constants::physics::MassPiPlus}, o2::constants::physics::MassOmegaC0, 1);
  }

  template <typename T>
  static auto cosThetaStarOmegaczeroToOmegaPi(const T& candidate)
  {
    return candidate.cosThetaStar(std::array{o2::constants::physics::MassOmegaMinus, o2::constants::physics::MassPiPlus}, o2::constants::physics::MassOmegaC0, 1);
  }

  template <typename T>
  static auto cosThetaStarOmegaczeroToOmegaK(const T& candidate)
  {
    return candidate.cosThetaStar(std::array{o2::constants::physics::MassOmegaMinus, o2::constants::physics::MassKPlus}, o2::constants::physics::MassOmegaC0, 1);
  }
};

namespace o2::hf_decay
{
using HfDecayChannel = int8_t;

/// @brief XiC0 and OmegaC0 candidates: main channels
namespace hf_cand_casc_bach_qa
{
enum DecayChannelMain : HfDecayChannel {
  // XiC0
  Xic0ToXiPi = 1,   // Xi- π+ → Λ π- π+ → p π- π- π+
  Xic0ToOmegaK = 2, // Ω- K+ → Λ K- K+ → p π- K- K+
  // OmegaC0
  Omegac0ToOmegaPi = 3,         // Ω- π+ → Λ π- π+ → p π- π- π+
  Omegac0ToOmegaK = 4,          // Ω- K+ → Λ K- K+ → p π- K- K+
  Omegac0ToXiPi = 5,            // Xi- π+ → Λ π- π+ → p π- π- π+
  NChannelsMain = Omegac0ToXiPi // last channel
};
} // namespace hf_cand_casc_bach_qa

} // namespace o2::hf_decay

// TODO: remove once the candidate creator is updated
namespace o2::aod
{
namespace hf_cand_xic0_omegac0
{
enum DecayType {
  XiczeroToXiPi = 0,
  OmegaczeroToXiPi,
  OmegaczeroToOmegaPi,
  OmegaczeroToOmegaK,
  OmegaczeroToOmegaPiOneMu,
  XiczeroToOmegaK
};
} // namespace hf_cand_xic0_omegac0
} // namespace o2::aod

#endif // PWGHF_DATAMODEL_CANDIDATERECONSTRUCTIONTABLESQA_H_
