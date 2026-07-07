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

// partition for phis and rhos, i.e. resonance that are their own antiparticle
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MAKE_RESONANCE_0_PARTITON(selection)                                                             \
  (o2::aod::femtobase::stored::pt > (selection).ptMin) &&                                                \
    (o2::aod::femtobase::stored::pt < (selection).ptMax) &&                                              \
    (o2::aod::femtobase::stored::eta > (selection).etaMin) &&                                            \
    (o2::aod::femtobase::stored::eta < (selection).etaMax) &&                                            \
    (o2::aod::femtobase::stored::phi > (selection).phiMin) &&                                            \
    (o2::aod::femtobase::stored::phi < (selection).phiMax) &&                                            \
    (o2::aod::femtobase::stored::mass > (selection).massMin) &&                                          \
    (o2::aod::femtobase::stored::mass < (selection).massMax) &&                                          \
    ifnode(o2::aod::femtotwotrackresonances::posDauHasHighMomentum,                                      \
           ncheckbit(o2::aod::femtotwotrackresonances::maskPosDau, (selection).posDauMaskAboveThres),    \
           ncheckbit(o2::aod::femtotwotrackresonances::maskPosDau, (selection).posDauMaskBelowThres)) && \
    ifnode(o2::aod::femtotwotrackresonances::negDauHasHighMomentum,                                      \
           ncheckbit(o2::aod::femtotwotrackresonances::maskNegDau, (selection).negDauMaskAboveThres),    \
           ncheckbit(o2::aod::femtotwotrackresonances::maskNegDau, (selection).negDauMaskBelowThres))

// partition for kstars, they have distinct antiparticle
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
    ifnode(o2::aod::femtotwotrackresonances::posDauHasHighMomentum,                                                                            \
           ncheckbit(o2::aod::femtotwotrackresonances::maskPosDau, (selection).posDauMaskAboveThres),                                          \
           ncheckbit(o2::aod::femtotwotrackresonances::maskPosDau, (selection).posDauMaskBelowThres)) &&                                       \
    ifnode(o2::aod::femtotwotrackresonances::negDauHasHighMomentum,                                                                            \
           ncheckbit(o2::aod::femtotwotrackresonances::maskNegDau, (selection).negDauMaskAboveThres),                                          \
           ncheckbit(o2::aod::femtotwotrackresonances::maskNegDau, (selection).negDauMaskBelowThres))

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

// macros for mc collisions (mc only)
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MAKE_MC_COLLISION_FILTER(selection)                                                                             \
  (o2::aod::femtocollisions::posZ >= (selection).vtxZMin && o2::aod::femtocollisions::posZ <= (selection).vtxZMax) &&   \
    (o2::aod::femtocollisions::mult >= (selection).multMin && o2::aod::femtocollisions::mult <= (selection).multMax) && \
    (o2::aod::femtocollisions::cent >= (selection).centMin && o2::aod::femtocollisions::cent <= (selection).centMax)

// macros for mc particle (mc only)
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define MAKE_MC_PARTICLE_PARTITION(selection)                                                                                                                                              \
  ifnode((selection).pdgCodeAbs.node() == 0, true, (selection).pdgCodeAbs == nabs(o2::aod::femtomcparticle::pdgCode)) &&                                                                   \
    ifnode((selection).chargeSign.node() != 0, ifnode((selection).chargeSign.node() > 0, o2::aod::femtobase::stored::signedPt > 0.f, o2::aod::femtobase::stored::signedPt < 0.f), true) && \
    (nabs(o2::aod::femtobase::stored::signedPt) > (selection).ptMin) &&                                                                                                                    \
    (nabs(o2::aod::femtobase::stored::signedPt) < (selection).ptMax) &&                                                                                                                    \
    (o2::aod::femtobase::stored::eta > (selection).etaMin) &&                                                                                                                              \
    (o2::aod::femtobase::stored::eta < (selection).etaMax) &&                                                                                                                              \
    (o2::aod::femtobase::stored::phi > (selection).phiMin) &&                                                                                                                              \
    (o2::aod::femtobase::stored::phi < (selection).phiMax)

#endif // PWGCF_FEMTO_CORE_PARTITIONS_H_
