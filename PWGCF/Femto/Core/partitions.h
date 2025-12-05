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
#define MAKE_COLLISION_FILTER(selection)                                                                                                                    \
  (femtocollisions::posZ >= selection.vtxZMin && femtocollisions::posZ <= selection.vtxZMax) &&                                                             \
    (femtocollisions::mult >= selection.multMin && femtocollisions::mult <= selection.multMax) &&                                                           \
    (femtocollisions::cent >= selection.centMin && femtocollisions::cent <= selection.centMax) &&                                                           \
    (femtocollisions::magField >= static_cast<int8_t>(selection.magFieldMin) && femtocollisions::magField <= static_cast<int8_t>(selection.magFieldMax)) && \
    ncheckbit(femtocollisions::mask, selection.collisionMask)

// standard track partition
#define MAKE_TRACK_PARTITION(selection)                                                                                                                                         \
  ifnode(selection.chargeSign.node() != 0, ifnode(selection.chargeSign.node() > 0, femtobase::stored::signedPt > 0.f, femtobase::stored::signedPt < 0.f), true) &&              \
    (nabs(selection.chargeAbs.node() * femtobase::stored::signedPt) > selection.ptMin) &&                                                                                       \
    (nabs(selection.chargeAbs.node() * femtobase::stored::signedPt) < selection.ptMax) &&                                                                                       \
    (femtobase::stored::eta > selection.etaMin) &&                                                                                                                              \
    (femtobase::stored::eta < selection.etaMax) &&                                                                                                                              \
    (femtobase::stored::phi > selection.phiMin) &&                                                                                                                              \
    (femtobase::stored::phi < selection.phiMax) &&                                                                                                                              \
    ifnode(nabs(selection.chargeAbs.node() * femtobase::stored::signedPt) * (nexp(femtobase::stored::eta) + nexp(-1.f * femtobase::stored::eta)) / (2.f) <= selection.pidThres, \
           ncheckbit(femtotracks::mask, selection.maskLowMomentum),                                                                                                             \
           ncheckbit(femtotracks::mask, selection.maskHighMomentum))

// partition for phis and rhos, i.e. resonance that are their own antiparticle
#define MAKE_RESONANCE_0_PARTITON(selection)                                            \
  (femtobase::stored::pt > selection.ptMin) &&                                          \
    (femtobase::stored::pt < selection.ptMax) &&                                        \
    (femtobase::stored::eta > selection.etaMin) &&                                      \
    (femtobase::stored::eta < selection.etaMax) &&                                      \
    (femtobase::stored::phi > selection.phiMin) &&                                      \
    (femtobase::stored::phi < selection.phiMax) &&                                      \
    (femtobase::stored::mass > selection.massMin) &&                                    \
    (femtobase::stored::mass < selection.massMax) &&                                    \
    ifnode(ncheckbit(femtotwotrackresonances::mask, selection.posDauBitForThres),       \
           ncheckbit(femtotwotrackresonances::mask, selection.posDauMaskAboveThres),    \
           ncheckbit(femtotwotrackresonances::mask, selection.posDauMaskBelowThres)) && \
    ifnode(ncheckbit(femtotwotrackresonances::mask, selection.negDauBitForThres),       \
           ncheckbit(femtotwotrackresonances::mask, selection.negDauMaskAboveThres),    \
           ncheckbit(femtotwotrackresonances::mask, selection.negDauMaskBelowThres))

// partition for kstars, they have distinct antiparticle
#define MAKE_RESONANCE_1_PARTITON(selection)                                                                               \
  ifnode(selection.sign.node() != 0,                                                                                       \
         ifnode(selection.sign.node() > 0, femtobase::stored::signedPt > 0.f, femtobase::stored::signedPt < 0.f), true) && \
    (nabs(femtobase::stored::signedPt) > selection.ptMin) &&                                                               \
    (nabs(femtobase::stored::signedPt) < selection.ptMax) &&                                                               \
    (femtobase::stored::eta > selection.etaMin) &&                                                                         \
    (femtobase::stored::eta < selection.etaMax) &&                                                                         \
    (femtobase::stored::phi > selection.phiMin) &&                                                                         \
    (femtobase::stored::phi < selection.phiMax) &&                                                                         \
    (femtobase::stored::mass > selection.massMin) &&                                                                       \
    (femtobase::stored::mass < selection.massMax) &&                                                                       \
    ifnode(ncheckbit(femtotwotrackresonances::mask, selection.posDauBitForThres),                                          \
           ncheckbit(femtotwotrackresonances::mask, selection.posDauMaskAboveThres),                                       \
           ncheckbit(femtotwotrackresonances::mask, selection.posDauMaskBelowThres)) &&                                    \
    ifnode(ncheckbit(femtotwotrackresonances::mask, selection.negDauBitForThres),                                          \
           ncheckbit(femtotwotrackresonances::mask, selection.negDauMaskAboveThres),                                       \
           ncheckbit(femtotwotrackresonances::mask, selection.negDauMaskBelowThres))

// partition for lambdas
#define MAKE_LAMBDA_PARTITION(selection)                                                                                   \
  ifnode(selection.sign.node() != 0,                                                                                       \
         ifnode(selection.sign.node() > 0, femtobase::stored::signedPt > 0.f, femtobase::stored::signedPt < 0.f), true) && \
    (nabs(femtobase::stored::signedPt) > selection.ptMin) &&                                                               \
    (nabs(femtobase::stored::signedPt) < selection.ptMax) &&                                                               \
    (femtobase::stored::eta > selection.etaMin) &&                                                                         \
    (femtobase::stored::eta < selection.etaMax) &&                                                                         \
    (femtobase::stored::phi > selection.phiMin) &&                                                                         \
    (femtobase::stored::phi < selection.phiMax) &&                                                                         \
    (femtobase::stored::mass > selection.massMin) &&                                                                       \
    (femtobase::stored::mass < selection.massMax) &&                                                                       \
    ncheckbit(femtov0s::mask, selection.mask)

// partition for k0shorts
// need special partition since k0shorts have no antiparticle
#define MAKE_K0SHORT_PARTITION(selection)            \
  (femtobase::stored::pt > selection.ptMin) &&       \
    (femtobase::stored::pt < selection.ptMax) &&     \
    (femtobase::stored::eta > selection.etaMin) &&   \
    (femtobase::stored::eta < selection.etaMax) &&   \
    (femtobase::stored::phi > selection.phiMin) &&   \
    (femtobase::stored::phi < selection.phiMax) &&   \
    (femtobase::stored::mass > selection.massMin) && \
    (femtobase::stored::mass < selection.massMax) && \
    ncheckbit(femtov0s::mask, selection.mask)

#define MAKE_CASCADE_PARTITION(selection)                                                                                  \
  ifnode(selection.sign.node() != 0,                                                                                       \
         ifnode(selection.sign.node() > 0, femtobase::stored::signedPt > 0.f, femtobase::stored::signedPt < 0.f), true) && \
    (nabs(femtobase::stored::signedPt) > selection.ptMin) &&                                                               \
    (nabs(femtobase::stored::signedPt) < selection.ptMax) &&                                                               \
    (femtobase::stored::eta > selection.etaMin) &&                                                                         \
    (femtobase::stored::eta < selection.etaMax) &&                                                                         \
    (femtobase::stored::phi > selection.phiMin) &&                                                                         \
    (femtobase::stored::phi < selection.phiMax) &&                                                                         \
    (femtobase::stored::mass > selection.massMin) &&                                                                       \
    (femtobase::stored::mass < selection.massMax) &&                                                                       \
    ncheckbit(femtocascades::mask, selection.mask)

#define MAKE_SIGMA_PARTITION(selection)                                                                                    \
  ifnode(selection.sign.node() != 0,                                                                                       \
         ifnode(selection.sign.node() > 0, femtobase::stored::signedPt > 0.f, femtobase::stored::signedPt < 0.f), true) && \
    (nabs(femtobase::stored::signedPt) > selection.ptMin) &&                                                               \
    (nabs(femtobase::stored::signedPt) < selection.ptMax) &&                                                               \
    (femtobase::stored::eta > selection.etaMin) &&                                                                         \
    (femtobase::stored::eta < selection.etaMax) &&                                                                         \
    (femtobase::stored::phi > selection.phiMin) &&                                                                         \
    (femtobase::stored::phi < selection.phiMax) &&                                                                         \
    (femtobase::stored::mass > selection.massMin) &&                                                                       \
    (femtobase::stored::mass < selection.massMax) &&                                                                       \
    ncheckbit(femtokinks::mask, selection.mask)

#define MAKE_SIGMAPLUS_PARTITION(selection)                                                                                \
  ifnode(selection.sign.node() != 0,                                                                                       \
         ifnode(selection.sign.node() > 0, femtobase::stored::signedPt > 0.f, femtobase::stored::signedPt < 0.f), true) && \
    (nabs(femtobase::stored::signedPt) > selection.ptMin) &&                                                               \
    (nabs(femtobase::stored::signedPt) < selection.ptMax) &&                                                               \
    (femtobase::stored::eta > selection.etaMin) &&                                                                         \
    (femtobase::stored::eta < selection.etaMax) &&                                                                         \
    (femtobase::stored::phi > selection.phiMin) &&                                                                         \
    (femtobase::stored::phi < selection.phiMax) &&                                                                         \
    (femtobase::stored::mass > selection.massMin) &&                                                                       \
    (femtobase::stored::mass < selection.massMax) &&                                                                       \
    ncheckbit(femtokinks::mask, selection.mask)

#endif // PWGCF_FEMTO_CORE_PARTITIONS_H_
