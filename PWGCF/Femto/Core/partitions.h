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

/// \file partitions.h
/// \brief common partition definitons
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_PARTITIONS_H_
#define PWGCF_FEMTO_CORE_PARTITIONS_H_

#include <Framework/Expressions.h>

// collsion selection
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MAKE_COLLISION_FILTER(selection)                                                                                \
  (o2::aod::femtocollisions::posZ >= (selection).vtxZMin && o2::aod::femtocollisions::posZ <= (selection).vtxZMax) &&   \
    (o2::aod::femtocollisions::mult >= (selection).multMin && o2::aod::femtocollisions::mult <= (selection).multMax) && \
    (o2::aod::femtocollisions::cent >= (selection).centMin && o2::aod::femtocollisions::cent <= (selection).centMax) && \
    (o2::aod::femtocollisions::magField >= o2::framework::expressions::as<int8_t>((selection).magFieldMin) &&           \
     o2::aod::femtocollisions::magField <= o2::framework::expressions::as<int8_t>((selection).magFieldMax)) &&          \
    ncheckbit(o2::aod::femtocollisions::mask, (selection).collisionMask)

// macro for track momentum, i.e. ||q|*pT/q| * cosh(eta)
// there is no ncosh function, so we have to make our own, i.e. cosh(x) = (exp(x)+exp(-x))/2
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define TRACK_MOMENTUM(chargeAbs, signedPt, eta) (nabs((chargeAbs) * (signedPt)) * (nexp(eta) + nexp(-1.f * (eta))) / 2.f)

// standard track partition
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MAKE_TRACK_PARTITION(selection)                                                                                                                                                  \
  ifnode((selection).chargeSign.node() != 0, ifnode((selection).chargeSign.node() > 0, o2::aod::femtobase::stored::signedPt > 0.f, o2::aod::femtobase::stored::signedPt < 0.f), true) && \
    (nabs((selection).chargeAbs * o2::aod::femtobase::stored::signedPt) > (selection).ptMin) &&                                                                                          \
    (nabs((selection).chargeAbs * o2::aod::femtobase::stored::signedPt) < (selection).ptMax) &&                                                                                          \
    (o2::aod::femtobase::stored::eta > (selection).etaMin) &&                                                                                                                            \
    (o2::aod::femtobase::stored::eta < (selection).etaMax) &&                                                                                                                            \
    (o2::aod::femtobase::stored::phi > (selection).phiMin) &&                                                                                                                            \
    (o2::aod::femtobase::stored::phi < (selection).phiMax) &&                                                                                                                            \
    ifnode(TRACK_MOMENTUM((selection).chargeAbs, o2::aod::femtobase::stored::signedPt, o2::aod::femtobase::stored::eta) <= (selection).pidThres,                                         \
           ncheckbit(o2::aod::femtotracks::mask, (selection).maskLowMomentum) &&                                                                                                         \
             (o2::aod::femtotracks::mask & (selection).rejectionMaskLowMomentum) == static_cast<o2::analysis::femto::datatypes::TrackMaskType>(0),                                       \
           ncheckbit(o2::aod::femtotracks::mask, (selection).maskHighMomentum) &&                                                                                                        \
             (o2::aod::femtotracks::mask & (selection).rejectionMaskHighMomentum) == static_cast<o2::analysis::femto::datatypes::TrackMaskType>(0))

// track partition with optional mass cut
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MAKE_TRACK_PARTITION_WITH_MASS(selection)               \
  MAKE_TRACK_PARTITION((selection)) &&                          \
    (o2::aod::femtobase::stored::mass > (selection).massMin) && \
    (o2::aod::femtobase::stored::mass < (selection).massMax)

// partition for resonance daughters (pion+/-, kaon+/- legs): unlike MAKE_TRACK_PARTITION, chargeSign
// and chargeAbs are compile-time constants on the selection struct (ConfResonanceDaughterSelection in
// trackBuilder.h), not runtime Configurables, so the sign comparison is just a comparison to a constant
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MAKE_RESONANCE_DAUGHTER_PARTITION(selection)                                                                                               \
  ((selection).chargeSign > 0 ? o2::aod::femtobase::stored::signedPt > 0.f : o2::aod::femtobase::stored::signedPt < 0.f) &&                        \
    (nabs((selection).chargeAbs * o2::aod::femtobase::stored::signedPt) > (selection).ptMin) &&                                                    \
    (nabs((selection).chargeAbs * o2::aod::femtobase::stored::signedPt) < (selection).ptMax) &&                                                    \
    (o2::aod::femtobase::stored::eta > (selection).etaMin) &&                                                                                      \
    (o2::aod::femtobase::stored::eta < (selection).etaMax) &&                                                                                      \
    (o2::aod::femtobase::stored::phi > (selection).phiMin) &&                                                                                      \
    (o2::aod::femtobase::stored::phi < (selection).phiMax) &&                                                                                      \
    ifnode(TRACK_MOMENTUM((selection).chargeAbs, o2::aod::femtobase::stored::signedPt, o2::aod::femtobase::stored::eta) <= (selection).pidThres,   \
           ncheckbit(o2::aod::femtotracks::mask, (selection).maskLowMomentum) &&                                                                   \
             (o2::aod::femtotracks::mask & (selection).rejectionMaskLowMomentum) == static_cast<o2::analysis::femto::datatypes::TrackMaskType>(0), \
           ncheckbit(o2::aod::femtotracks::mask, (selection).maskHighMomentum) &&                                                                  \
             (o2::aod::femtotracks::mask & (selection).rejectionMaskHighMomentum) == static_cast<o2::analysis::femto::datatypes::TrackMaskType>(0))

// partition for phis and rhos, i.e. resonance that are their own antiparticle -- both daughters are
// the same species, so the selection carries a single dauMask pair applied to both stored columns
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MAKE_RESONANCE_0_PARTITON(selection)                                                          \
  (o2::aod::femtobase::stored::pt > (selection).ptMin) &&                                             \
    (o2::aod::femtobase::stored::pt < (selection).ptMax) &&                                           \
    (o2::aod::femtobase::stored::eta > (selection).etaMin) &&                                         \
    (o2::aod::femtobase::stored::eta < (selection).etaMax) &&                                         \
    (o2::aod::femtobase::stored::phi > (selection).phiMin) &&                                         \
    (o2::aod::femtobase::stored::phi < (selection).phiMax) &&                                         \
    (o2::aod::femtobase::stored::mass > (selection).massMin) &&                                       \
    (o2::aod::femtobase::stored::mass < (selection).massMax) &&                                       \
    ifnode(o2::aod::femtotwotrackresonances::posDauHasHighMomentum,                                   \
           ncheckbit(o2::aod::femtotwotrackresonances::maskPosDau, (selection).dauMaskAboveThres),    \
           ncheckbit(o2::aod::femtotwotrackresonances::maskPosDau, (selection).dauMaskBelowThres)) && \
    ifnode(o2::aod::femtotwotrackresonances::negDauHasHighMomentum,                                   \
           ncheckbit(o2::aod::femtotwotrackresonances::maskNegDau, (selection).dauMaskAboveThres),    \
           ncheckbit(o2::aod::femtotwotrackresonances::maskNegDau, (selection).dauMaskBelowThres))

// partition for kstars, they have distinct antiparticle -- K*0 (K+ pi-) and K*0bar (K- pi+) swap
// which physical daughter (pos/neg) is the kaon vs the pion, so the selection carries masks named by
// species (kaon/pion) and this macro resolves pos/neg <-> kaon/pion from the candidate's own
// signedPt, the same way it's reconstructed (no need for the user to reason about this downstream)
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MAKE_RESONANCE_1_PARTITON(selection)                                                                                                   \
  ifnode((selection).sign.node() != 0,                                                                                                         \
         ifnode((selection).sign.node() > 0, o2::aod::femtobase::stored::signedPt > 0.f, o2::aod::femtobase::stored::signedPt < 0.f), true) && \
    (nabs(o2::aod::femtobase::stored::signedPt) > (selection).ptMin) &&                                                                        \
    (nabs(o2::aod::femtobase::stored::signedPt) < (selection).ptMax) &&                                                                        \
    (o2::aod::femtobase::stored::eta > (selection).etaMin) &&                                                                                  \
    (o2::aod::femtobase::stored::eta < (selection).etaMax) &&                                                                                  \
    (o2::aod::femtobase::stored::phi > (selection).phiMin) &&                                                                                  \
    (o2::aod::femtobase::stored::phi < (selection).phiMax) &&                                                                                  \
    (o2::aod::femtobase::stored::mass > (selection).massMin) &&                                                                                \
    (o2::aod::femtobase::stored::mass < (selection).massMax) &&                                                                                \
    ifnode(o2::aod::femtobase::stored::signedPt > 0.f,                                                                                         \
           ifnode(o2::aod::femtotwotrackresonances::posDauHasHighMomentum,                                                                     \
                  ncheckbit(o2::aod::femtotwotrackresonances::maskPosDau, (selection).kaonDauMaskAboveThres),                                  \
                  ncheckbit(o2::aod::femtotwotrackresonances::maskPosDau, (selection).kaonDauMaskBelowThres)) &&                               \
             ifnode(o2::aod::femtotwotrackresonances::negDauHasHighMomentum,                                                                   \
                    ncheckbit(o2::aod::femtotwotrackresonances::maskNegDau, (selection).pionDauMaskAboveThres),                                \
                    ncheckbit(o2::aod::femtotwotrackresonances::maskNegDau, (selection).pionDauMaskBelowThres)),                               \
           ifnode(o2::aod::femtotwotrackresonances::posDauHasHighMomentum,                                                                     \
                  ncheckbit(o2::aod::femtotwotrackresonances::maskPosDau, (selection).pionDauMaskAboveThres),                                  \
                  ncheckbit(o2::aod::femtotwotrackresonances::maskPosDau, (selection).pionDauMaskBelowThres)) &&                               \
             ifnode(o2::aod::femtotwotrackresonances::negDauHasHighMomentum,                                                                   \
                    ncheckbit(o2::aod::femtotwotrackresonances::maskNegDau, (selection).kaonDauMaskAboveThres),                                \
                    ncheckbit(o2::aod::femtotwotrackresonances::maskNegDau, (selection).kaonDauMaskBelowThres)))

// partition for resonances built from two photons (pi0, eta, ...): unlike MAKE_RESONANCE_0/1_PARTITON
// there is no sign and no momentum-threshold PID switch, and just a single mask check -- both
// daughters are plain photons with the same selection requirement, so the candidate only stores the
// AND of both legs' own masks (see femtotwophotonresonances::MaskDau)
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MAKE_TWOPHOTONRESONANCE_PARTITION(selection)            \
  (o2::aod::femtobase::stored::pt > (selection).ptMin) &&       \
    (o2::aod::femtobase::stored::pt < (selection).ptMax) &&     \
    (o2::aod::femtobase::stored::eta > (selection).etaMin) &&   \
    (o2::aod::femtobase::stored::eta < (selection).etaMax) &&   \
    (o2::aod::femtobase::stored::phi > (selection).phiMin) &&   \
    (o2::aod::femtobase::stored::phi < (selection).phiMax) &&   \
    (o2::aod::femtobase::stored::mass > (selection).massMin) && \
    (o2::aod::femtobase::stored::mass < (selection).massMax) && \
    ncheckbit(o2::aod::femtotwophotonresonances::maskDau, (selection).dauMask)

// partition for resonances built from two tracks and a two-photon resonance (e.g. omega(782) ->
// pi+ pi- pi0): same pattern as MAKE_RESONANCE_0_PARTITON (both charged daughters are the same
// species -- pions, just opposite charge, no kstar-style species swap -- so a single dauMask pair
// covers both), plus a 3rd mask for the neutral daughter (no momentum-threshold split there)
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MAKE_TRACKTRACKTWOPHOTONRESONANCE_PARTITION(selection)                                                   \
  (o2::aod::femtobase::stored::pt > (selection).ptMin) &&                                                        \
    (o2::aod::femtobase::stored::pt < (selection).ptMax) &&                                                      \
    (o2::aod::femtobase::stored::eta > (selection).etaMin) &&                                                    \
    (o2::aod::femtobase::stored::eta < (selection).etaMax) &&                                                    \
    (o2::aod::femtobase::stored::phi > (selection).phiMin) &&                                                    \
    (o2::aod::femtobase::stored::phi < (selection).phiMax) &&                                                    \
    (o2::aod::femtobase::stored::mass > (selection).massMin) &&                                                  \
    (o2::aod::femtobase::stored::mass < (selection).massMax) &&                                                  \
    ifnode(o2::aod::femtotrackstwophotonresonances::posDauHasHighMomentum,                                       \
           ncheckbit(o2::aod::femtotrackstwophotonresonances::maskPosDau, (selection).pionDauMaskAboveThres),    \
           ncheckbit(o2::aod::femtotrackstwophotonresonances::maskPosDau, (selection).pionDauMaskBelowThres)) && \
    ifnode(o2::aod::femtotrackstwophotonresonances::negDauHasHighMomentum,                                       \
           ncheckbit(o2::aod::femtotrackstwophotonresonances::maskNegDau, (selection).pionDauMaskAboveThres),    \
           ncheckbit(o2::aod::femtotrackstwophotonresonances::maskNegDau, (selection).pionDauMaskBelowThres)) && \
    ncheckbit(o2::aod::femtotrackstwophotonresonances::maskNeutralDau, (selection).neutralDauMask)

// partition for lambdas
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MAKE_LAMBDA_PARTITION(selection)                                                                                                       \
  ifnode((selection).sign.node() != 0,                                                                                                         \
         ifnode((selection).sign.node() > 0, o2::aod::femtobase::stored::signedPt > 0.f, o2::aod::femtobase::stored::signedPt < 0.f), true) && \
    (nabs(o2::aod::femtobase::stored::signedPt) > (selection).ptMin) &&                                                                        \
    (nabs(o2::aod::femtobase::stored::signedPt) < (selection).ptMax) &&                                                                        \
    (o2::aod::femtobase::stored::eta > (selection).etaMin) &&                                                                                  \
    (o2::aod::femtobase::stored::eta < (selection).etaMax) &&                                                                                  \
    (o2::aod::femtobase::stored::phi > (selection).phiMin) &&                                                                                  \
    (o2::aod::femtobase::stored::phi < (selection).phiMax) &&                                                                                  \
    (o2::aod::femtobase::stored::mass > (selection).massMin) &&                                                                                \
    (o2::aod::femtobase::stored::mass < (selection).massMax) &&                                                                                \
    ncheckbit(o2::aod::femtov0s::mask, (selection).mask)

// partition for k0shorts
// need special partition since k0shorts have no antiparticle
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MAKE_K0SHORT_PARTITION(selection)                       \
  (o2::aod::femtobase::stored::pt > (selection).ptMin) &&       \
    (o2::aod::femtobase::stored::pt < (selection).ptMax) &&     \
    (o2::aod::femtobase::stored::eta > (selection).etaMin) &&   \
    (o2::aod::femtobase::stored::eta < (selection).etaMax) &&   \
    (o2::aod::femtobase::stored::phi > (selection).phiMin) &&   \
    (o2::aod::femtobase::stored::phi < (selection).phiMax) &&   \
    (o2::aod::femtobase::stored::mass > (selection).massMin) && \
    (o2::aod::femtobase::stored::mass < (selection).massMax) && \
    ncheckbit(o2::aod::femtov0s::mask, (selection).mask)

// partition for photons (PCM)
// no sign/mass handling: photons are their own antiparticle and have no mass window
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MAKE_PHOTON_PARTITION(selection)                      \
  (o2::aod::femtobase::stored::pt > (selection).ptMin) &&     \
    (o2::aod::femtobase::stored::pt < (selection).ptMax) &&   \
    (o2::aod::femtobase::stored::eta > (selection).etaMin) && \
    (o2::aod::femtobase::stored::eta < (selection).etaMax) && \
    (o2::aod::femtobase::stored::phi > (selection).phiMin) && \
    (o2::aod::femtobase::stored::phi < (selection).phiMax) && \
    ncheckbit(o2::aod::femtophotons::mask, (selection).mask)

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MAKE_CASCADE_PARTITION(selection)                                                                                                      \
  ifnode((selection).sign.node() != 0,                                                                                                         \
         ifnode((selection).sign.node() > 0, o2::aod::femtobase::stored::signedPt > 0.f, o2::aod::femtobase::stored::signedPt < 0.f), true) && \
    (nabs(o2::aod::femtobase::stored::signedPt) > (selection).ptMin) &&                                                                        \
    (nabs(o2::aod::femtobase::stored::signedPt) < (selection).ptMax) &&                                                                        \
    (o2::aod::femtobase::stored::eta > (selection).etaMin) &&                                                                                  \
    (o2::aod::femtobase::stored::eta < (selection).etaMax) &&                                                                                  \
    (o2::aod::femtobase::stored::phi > (selection).phiMin) &&                                                                                  \
    (o2::aod::femtobase::stored::phi < (selection).phiMax) &&                                                                                  \
    (o2::aod::femtobase::stored::mass > (selection).massMin) &&                                                                                \
    (o2::aod::femtobase::stored::mass < (selection).massMax) &&                                                                                \
    ncheckbit(o2::aod::femtocascades::mask, (selection).mask)

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MAKE_SIGMA_PARTITION(selection)                                                                                                        \
  ifnode((selection).sign.node() != 0,                                                                                                         \
         ifnode((selection).sign.node() > 0, o2::aod::femtobase::stored::signedPt > 0.f, o2::aod::femtobase::stored::signedPt < 0.f), true) && \
    (nabs(o2::aod::femtobase::stored::signedPt) > (selection).ptMin) &&                                                                        \
    (nabs(o2::aod::femtobase::stored::signedPt) < (selection).ptMax) &&                                                                        \
    (o2::aod::femtobase::stored::eta > (selection).etaMin) &&                                                                                  \
    (o2::aod::femtobase::stored::eta < (selection).etaMax) &&                                                                                  \
    (o2::aod::femtobase::stored::phi > (selection).phiMin) &&                                                                                  \
    (o2::aod::femtobase::stored::phi < (selection).phiMax) &&                                                                                  \
    (o2::aod::femtobase::stored::mass > (selection).massMin) &&                                                                                \
    (o2::aod::femtobase::stored::mass < (selection).massMax) &&                                                                                \
    ncheckbit(o2::aod::femtokinks::mask, (selection).mask)

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MAKE_SIGMAPLUS_PARTITION(selection)                                                                                                    \
  ifnode((selection).sign.node() != 0,                                                                                                         \
         ifnode((selection).sign.node() > 0, o2::aod::femtobase::stored::signedPt > 0.f, o2::aod::femtobase::stored::signedPt < 0.f), true) && \
    (nabs(o2::aod::femtobase::stored::signedPt) > (selection).ptMin) &&                                                                        \
    (nabs(o2::aod::femtobase::stored::signedPt) < (selection).ptMax) &&                                                                        \
    (o2::aod::femtobase::stored::eta > (selection).etaMin) &&                                                                                  \
    (o2::aod::femtobase::stored::eta < (selection).etaMax) &&                                                                                  \
    (o2::aod::femtobase::stored::phi > (selection).phiMin) &&                                                                                  \
    (o2::aod::femtobase::stored::phi < (selection).phiMax) &&                                                                                  \
    (o2::aod::femtobase::stored::mass > (selection).massMin) &&                                                                                \
    (o2::aod::femtobase::stored::mass < (selection).massMax) &&                                                                                \
    ncheckbit(o2::aod::femtokinks::mask, (selection).mask)

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MAKE_SIGMAPLUSPCM_PARTITION(selection)                                                                                                 \
  ifnode((selection).sign.node() != 0,                                                                                                         \
         ifnode((selection).sign.node() > 0, o2::aod::femtobase::stored::signedPt > 0.f, o2::aod::femtobase::stored::signedPt < 0.f), true) && \
    (nabs(o2::aod::femtobase::stored::signedPt) > (selection).ptMin) &&                                                                        \
    (nabs(o2::aod::femtobase::stored::signedPt) < (selection).ptMax) &&                                                                        \
    (o2::aod::femtobase::stored::eta > (selection).etaMin) &&                                                                                  \
    (o2::aod::femtobase::stored::eta < (selection).etaMax) &&                                                                                  \
    (o2::aod::femtobase::stored::phi > (selection).phiMin) &&                                                                                  \
    (o2::aod::femtobase::stored::phi < (selection).phiMax) &&                                                                                  \
    (o2::aod::femtobase::stored::mass > (selection).massMin) &&                                                                                \
    (o2::aod::femtobase::stored::mass < (selection).massMax) &&                                                                                \
    ncheckbit(o2::aod::femtokinks::mask, (selection).mask)

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MAKE_SIGMAPLUSPCM_ML_PARTITION(selection)                                                                                              \
  ifnode((selection).sign.node() != 0,                                                                                                         \
         ifnode((selection).sign.node() > 0, o2::aod::femtobase::stored::signedPt > 0.f, o2::aod::femtobase::stored::signedPt < 0.f), true) && \
    (nabs(o2::aod::femtobase::stored::signedPt) > (selection).ptMin) &&                                                                        \
    (nabs(o2::aod::femtobase::stored::signedPt) < (selection).ptMax) &&                                                                        \
    (o2::aod::femtobase::stored::eta > (selection).etaMin) &&                                                                                  \
    (o2::aod::femtobase::stored::eta < (selection).etaMax) &&                                                                                  \
    (o2::aod::femtobase::stored::phi > (selection).phiMin) &&                                                                                  \
    (o2::aod::femtobase::stored::phi < (selection).phiMax) &&                                                                                  \
    (o2::aod::femtobase::stored::mass > (selection).massMin) &&                                                                                \
    (o2::aod::femtobase::stored::mass < (selection).massMax) &&                                                                                \
    (o2::aod::femtokinks::mlScore > (selection).mlScoreMin) &&                                                                                 \
    ncheckbit(o2::aod::femtokinks::mask, (selection).mask)

// partition for D0
// sign of signedPt encodes the hypothesis: D0 (+) / D0bar (-)
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MAKE_D0_PARTITION(selection)                                                                                                           \
  ifnode((selection).sign.node() != 0,                                                                                                         \
         ifnode((selection).sign.node() > 0, o2::aod::femtobase::stored::signedPt > 0.f, o2::aod::femtobase::stored::signedPt < 0.f), true) && \
    (nabs(o2::aod::femtobase::stored::signedPt) > (selection).ptMin) &&                                                                        \
    (nabs(o2::aod::femtobase::stored::signedPt) < (selection).ptMax) &&                                                                        \
    (o2::aod::femtobase::stored::eta > (selection).etaMin) &&                                                                                  \
    (o2::aod::femtobase::stored::eta < (selection).etaMax) &&                                                                                  \
    (o2::aod::femtobase::stored::phi > (selection).phiMin) &&                                                                                  \
    (o2::aod::femtobase::stored::phi < (selection).phiMax) &&                                                                                  \
    (o2::aod::femtobase::stored::mass > (selection).massMin) &&                                                                                \
    (o2::aod::femtobase::stored::mass < (selection).massMax) &&                                                                                \
    ncheckbit(o2::aod::femtocharmhadrons::mask, (selection).mask)

// partition for 3-prong charm hadrons
// sign of signedPt encodes particle (+) / antiparticle (-)
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MAKE_CHARM3PRONG_PARTITION(selection)                                                                                                  \
  ifnode((selection).sign.node() != 0,                                                                                                         \
         ifnode((selection).sign.node() > 0, o2::aod::femtobase::stored::signedPt > 0.f, o2::aod::femtobase::stored::signedPt < 0.f), true) && \
    (nabs(o2::aod::femtobase::stored::signedPt) > (selection).ptMin) &&                                                                        \
    (nabs(o2::aod::femtobase::stored::signedPt) < (selection).ptMax) &&                                                                        \
    (o2::aod::femtobase::stored::eta > (selection).etaMin) &&                                                                                  \
    (o2::aod::femtobase::stored::eta < (selection).etaMax) &&                                                                                  \
    (o2::aod::femtobase::stored::phi > (selection).phiMin) &&                                                                                  \
    (o2::aod::femtobase::stored::phi < (selection).phiMax) &&                                                                                  \
    (o2::aod::femtobase::stored::mass > (selection).massMin) &&                                                                                \
    (o2::aod::femtobase::stored::mass < (selection).massMax) &&                                                                                \
    ncheckbit(o2::aod::femtocharmhadrons::mask, (selection).mask)

// macros for mc collisions (mc only)
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MAKE_MC_COLLISION_FILTER(selection)                                                                             \
  (o2::aod::femtocollisions::posZ >= (selection).vtxZMin && o2::aod::femtocollisions::posZ <= (selection).vtxZMax) &&   \
    (o2::aod::femtocollisions::mult >= (selection).multMin && o2::aod::femtocollisions::mult <= (selection).multMax) && \
    (o2::aod::femtocollisions::cent >= (selection).centMin && o2::aod::femtocollisions::cent <= (selection).centMax)

// macros for mc particle (mc only)
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MAKE_MC_PARTICLE_PARTITION(selection)                                                                                                                                              \
  ifnode((selection).requireOrigin.node(), o2::aod::femtomcparticle::origin == o2::framework::expressions::as<uint8_t>((selection).origin), true) &&                                       \
    ifnode((selection).pdgCodeAbs.node() == 0, true, (selection).pdgCodeAbs == nabs(o2::aod::femtomcparticle::pdgCode)) &&                                                                 \
    ifnode((selection).chargeSign.node() != 0, ifnode((selection).chargeSign.node() > 0, o2::aod::femtobase::stored::signedPt > 0.f, o2::aod::femtobase::stored::signedPt < 0.f), true) && \
    (nabs(o2::aod::femtobase::stored::signedPt) > (selection).ptMin) &&                                                                                                                    \
    (nabs(o2::aod::femtobase::stored::signedPt) < (selection).ptMax) &&                                                                                                                    \
    (o2::aod::femtobase::stored::eta > (selection).etaMin) &&                                                                                                                              \
    (o2::aod::femtobase::stored::eta < (selection).etaMax) &&                                                                                                                              \
    (o2::aod::femtobase::stored::phi > (selection).phiMin) &&                                                                                                                              \
    (o2::aod::femtobase::stored::phi < (selection).phiMax)

#endif // PWGCF_FEMTO_CORE_PARTITIONS_H_
