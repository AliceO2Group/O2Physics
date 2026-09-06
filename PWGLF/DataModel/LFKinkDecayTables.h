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

///
/// \file LFKinkDecayTables.h
/// \brief Slim tables for kinks
/// \author Francesco Mazzaschi <francesco.mazzaschi@cern.ch>
///

#ifndef PWGLF_DATAMODEL_LFKINKDECAYTABLES_H_
#define PWGLF_DATAMODEL_LFKINKDECAYTABLES_H_

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>

#include <array>
#include <cmath>

namespace o2::aod
{

namespace kinkcand
{

DECLARE_SOA_INDEX_COLUMN_FULL(TrackMoth, trackMoth, int, TracksIU, "_Moth"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(TrackDaug, trackDaug, int, TracksIU, "_Daug"); //!
DECLARE_SOA_INDEX_COLUMN(Collision, collision);                              //!

DECLARE_SOA_COLUMN(XDecVtx, xDecVtx, float);         //! Decay vertex of the candidate (x direction)
DECLARE_SOA_COLUMN(YDecVtx, yDecVtx, float);         //! Decay vertex of the candidate (y direction)
DECLARE_SOA_COLUMN(ZDecVtx, zDecVtx, float);         //! Decay vertex of the candidate (z direction)
DECLARE_SOA_COLUMN(PxMoth, pxMoth, float);           //! Px of the mother kink
DECLARE_SOA_COLUMN(PyMoth, pyMoth, float);           //! Py of the mother kink
DECLARE_SOA_COLUMN(PzMoth, pzMoth, float);           //! Pz of the mother kink
DECLARE_SOA_COLUMN(PxDaug, pxDaug, float);           //! Px of the daughter kink
DECLARE_SOA_COLUMN(PyDaug, pyDaug, float);           //! Py of the daughter kink
DECLARE_SOA_COLUMN(PzDaug, pzDaug, float);           //! Pz of the daughter kink
DECLARE_SOA_COLUMN(MothSign, mothSign, int);         //! Sign of the mother kink
DECLARE_SOA_COLUMN(DcaMothPv, dcaMothPv, float);     //! DCA of the mother to the primary vertex
DECLARE_SOA_COLUMN(DcaDaugPv, dcaDaugPv, float);     //! DCA of the daughter kink to the primary vertex
DECLARE_SOA_COLUMN(DcaKinkTopo, dcaKinkTopo, float); //! DCA of the kink topology

DECLARE_SOA_COLUMN(NSigmaTPCPi, nSigmaTPCPi, float); //! Number of sigmas for the pion candidate from Sigma kink in TPC
DECLARE_SOA_COLUMN(NSigmaTPCPr, nSigmaTPCPr, float); //! Number of sigmas for the proton candidate from Sigma kink in TPC
DECLARE_SOA_COLUMN(NSigmaTPCKa, nSigmaTPCKa, float); //! Number of sigmas for the kaon candidate from Sigma kink in TPC
DECLARE_SOA_COLUMN(NSigmaTOFPi, nSigmaTOFPi, float); //! Number of sigmas for the pion candidate from Sigma kink in TOF
DECLARE_SOA_COLUMN(NSigmaTOFPr, nSigmaTOFPr, float); //! Number of sigmas for the proton candidate from Sigma kink in TOF
DECLARE_SOA_COLUMN(NSigmaTOFKa, nSigmaTOFKa, float); //! Number of sigmas for the kaon candidate from Sigma kink in TOF

// MC Columns
DECLARE_SOA_COLUMN(MothPdgCode, mothPdgCode, int);            //! PDG code of the Sigma daughter
DECLARE_SOA_COLUMN(DaugPdgCode, daugPdgCode, int);            //! PDG code of the kink daughter
DECLARE_SOA_COLUMN(PtMC, ptMC, float);                        //! pT of the candidate in MC
DECLARE_SOA_COLUMN(PzMC, pzMC, float);                        //! pZ of the candidate in MC
DECLARE_SOA_COLUMN(MassMC, massMC, float);                    //! Invariant mass of the candidate in MC
DECLARE_SOA_COLUMN(DecayRadiusMC, decayRadiusMC, float);      //! Decay radius of the candidate in MC
DECLARE_SOA_COLUMN(CollisionIdCheck, collisionIdCheck, bool); //! Check if mcDaughter collision ID matches the reconstructed collision ID

// DYNAMIC COLUMNS

DECLARE_SOA_DYNAMIC_COLUMN(PxDaugNeut, pxDaugNeut, //! Px of the daughter neutral particle
                           [](float pxmoth, float pxdau) -> float { return pxmoth - pxdau; });

DECLARE_SOA_DYNAMIC_COLUMN(PyDaugNeut, pyDaugNeut, //! Py of the daughter neutral particle
                           [](float pymoth, float pydau) -> float { return pymoth - pydau; });

DECLARE_SOA_DYNAMIC_COLUMN(PzDaugNeut, pzDaugNeut, //! Pz of the daughter neutral particle
                           [](float pzmoth, float pzdau) -> float { return pzmoth - pzdau; });

DECLARE_SOA_DYNAMIC_COLUMN(PtMoth, ptMoth, //! pT of the mother kink
                           [](float pxmoth, float pymoth) -> float { return std::hypot(pxmoth, pymoth); });

DECLARE_SOA_DYNAMIC_COLUMN(PtDaug, ptDaug, //!
                           [](float pxdaug, float pydaug) -> float { return std::hypot(pxdaug, pydaug); });

DECLARE_SOA_DYNAMIC_COLUMN(MSigmaMinus, mSigmaMinus, //! mass under sigma minus hypothesis
                           [](float pxmoth, float pymoth, float pzmoth, float pxch, float pych, float pzch) -> float {
                            float pxneut = pxmoth - pxch;
                            float pyneut = pymoth - pych;
                            float pzneut = pzmoth - pzch;
                            return RecoDecay::m(std::array{std::array{pxch, pych, pzch}, std::array{pxneut, pyneut, pzneut}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassNeutron}); });

DECLARE_SOA_DYNAMIC_COLUMN(MSigmaPlus, mSigmaPlus, //! mass under sigma plus hypothesis
                           [](float pxmoth, float pymoth, float pzmoth, float pxch, float pych, float pzch) -> float {
                            float pxneut = pxmoth - pxch;
                            float pyneut = pymoth - pych;
                            float pzneut = pzmoth - pzch;
                            return RecoDecay::m(std::array{std::array{pxch, pych, pzch}, std::array{pxneut, pyneut, pzneut}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionNeutral}); });

DECLARE_SOA_DYNAMIC_COLUMN(MXiMinus, mXiMinus, //! mass under Xi minus hypothesis
                           [](float pxmoth, float pymoth, float pzmoth, float pxch, float pych, float pzch) -> float {
                            float pxneut = pxmoth - pxch;
                            float pyneut = pymoth - pych;
                            float pzneut = pzmoth - pzch;
                            return RecoDecay::m(std::array{std::array{pxch, pych, pzch}, std::array{pxneut, pyneut, pzneut}}, std::array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassLambda}); });

} // namespace kinkcand

DECLARE_SOA_TABLE(KinkCands, "AOD", "KINKCANDS",
                  o2::soa::Index<>, kinkcand::CollisionId, kinkcand::TrackMothId, kinkcand::TrackDaugId,
                  kinkcand::XDecVtx, kinkcand::YDecVtx, kinkcand::ZDecVtx,
                  kinkcand::MothSign, kinkcand::PxMoth, kinkcand::PyMoth, kinkcand::PzMoth,
                  kinkcand::PxDaug, kinkcand::PyDaug, kinkcand::PzDaug,
                  kinkcand::DcaMothPv, kinkcand::DcaDaugPv, kinkcand::DcaKinkTopo,

                  // dynamic columns
                  kinkcand::PxDaugNeut<kinkcand::PxMoth, kinkcand::PxDaug>,
                  kinkcand::PyDaugNeut<kinkcand::PyMoth, kinkcand::PyDaug>,
                  kinkcand::PzDaugNeut<kinkcand::PzMoth, kinkcand::PzDaug>,
                  kinkcand::PtMoth<kinkcand::PxMoth, kinkcand::PyMoth>,
                  kinkcand::PtDaug<kinkcand::PxDaug, kinkcand::PyDaug>,
                  kinkcand::MSigmaMinus<kinkcand::PxMoth, kinkcand::PyMoth, kinkcand::PzMoth, kinkcand::PxDaug, kinkcand::PyDaug, kinkcand::PzDaug>,
                  kinkcand::MSigmaPlus<kinkcand::PxMoth, kinkcand::PyMoth, kinkcand::PzMoth, kinkcand::PxDaug, kinkcand::PyDaug, kinkcand::PzDaug>,
                  kinkcand::MXiMinus<kinkcand::PxMoth, kinkcand::PyMoth, kinkcand::PzMoth, kinkcand::PxDaug, kinkcand::PyDaug, kinkcand::PzDaug>);

DECLARE_SOA_TABLE(KinkCandsUnbound, "AOD", "UBKINKCANDS",
                  o2::soa::Index<>, kinkcand::XDecVtx, kinkcand::YDecVtx, kinkcand::ZDecVtx,
                  kinkcand::MothSign, kinkcand::PxMoth, kinkcand::PyMoth, kinkcand::PzMoth,
                  kinkcand::PxDaug, kinkcand::PyDaug, kinkcand::PzDaug,
                  kinkcand::DcaMothPv, kinkcand::DcaDaugPv, kinkcand::DcaKinkTopo,

                  // dynamic columns
                  kinkcand::PxDaugNeut<kinkcand::PxMoth, kinkcand::PxDaug>,
                  kinkcand::PyDaugNeut<kinkcand::PyMoth, kinkcand::PyDaug>,
                  kinkcand::PzDaugNeut<kinkcand::PzMoth, kinkcand::PzDaug>,
                  kinkcand::PtMoth<kinkcand::PxMoth, kinkcand::PyMoth>,
                  kinkcand::PtDaug<kinkcand::PxDaug, kinkcand::PyDaug>,
                  kinkcand::MSigmaMinus<kinkcand::PxMoth, kinkcand::PyMoth, kinkcand::PzMoth, kinkcand::PxDaug, kinkcand::PyDaug, kinkcand::PzDaug>,
                  kinkcand::MSigmaPlus<kinkcand::PxMoth, kinkcand::PyMoth, kinkcand::PzMoth, kinkcand::PxDaug, kinkcand::PyDaug, kinkcand::PzDaug>,
                  kinkcand::MXiMinus<kinkcand::PxMoth, kinkcand::PyMoth, kinkcand::PzMoth, kinkcand::PxDaug, kinkcand::PyDaug, kinkcand::PzDaug>);

DECLARE_SOA_TABLE(SlimKinkCands, "AOD", "SLIMKINKCANDS",
                  kinkcand::XDecVtx, kinkcand::YDecVtx, kinkcand::ZDecVtx,
                  kinkcand::PxMoth, kinkcand::PyMoth, kinkcand::PzMoth,
                  kinkcand::PxDaug, kinkcand::PyDaug, kinkcand::PzDaug,
                  kinkcand::DcaMothPv, kinkcand::DcaDaugPv, kinkcand::DcaKinkTopo,
                  kinkcand::MothSign,
                  kinkcand::NSigmaTPCPi, kinkcand::NSigmaTPCPr, kinkcand::NSigmaTPCKa,
                  kinkcand::NSigmaTOFPi, kinkcand::NSigmaTOFPr, kinkcand::NSigmaTOFKa);

DECLARE_SOA_TABLE(SlimKinkCandsMC, "AOD", "SLIMKINKCANDSMC",
                  kinkcand::XDecVtx, kinkcand::YDecVtx, kinkcand::ZDecVtx,
                  kinkcand::PxMoth, kinkcand::PyMoth, kinkcand::PzMoth,
                  kinkcand::PxDaug, kinkcand::PyDaug, kinkcand::PzDaug,
                  kinkcand::DcaMothPv, kinkcand::DcaDaugPv, kinkcand::DcaKinkTopo,
                  kinkcand::MothSign,
                  kinkcand::NSigmaTPCPi, kinkcand::NSigmaTPCPr, kinkcand::NSigmaTPCKa,
                  kinkcand::NSigmaTOFPi, kinkcand::NSigmaTOFPr, kinkcand::NSigmaTOFKa,
                  kinkcand::MothPdgCode, kinkcand::DaugPdgCode,
                  kinkcand::PtMC, kinkcand::PzMC, kinkcand::MassMC, kinkcand::DecayRadiusMC, kinkcand::CollisionIdCheck);

namespace sigmapluscand
{

DECLARE_SOA_COLUMN(XDecVtx, xDecVtx, float);               //! Decay vertex of the candidate (x direction)
DECLARE_SOA_COLUMN(YDecVtx, yDecVtx, float);               //! Decay vertex of the candidate (y direction)
DECLARE_SOA_COLUMN(ZDecVtx, zDecVtx, float);               //! Decay vertex of the candidate (z direction)
DECLARE_SOA_COLUMN(Radius, radius, float);                 //! Decay radius of the candidate (cm)
DECLARE_SOA_COLUMN(FlightDistance, flightDistance, float); //! Flight distance of the candidate (PV to decay vertex, cm)
DECLARE_SOA_COLUMN(DcaProtonGamma, dcaProtonGamma, float); //! DCA between proton and photon at the fitted vertex (cm)
DECLARE_SOA_COLUMN(Chi2, chi2, float);                     //! chi2 of the proton-photon vertex fit

DECLARE_SOA_COLUMN(PxProton, pxProton, float); //! Px of the proton
DECLARE_SOA_COLUMN(PyProton, pyProton, float); //! Py of the proton
DECLARE_SOA_COLUMN(PzProton, pzProton, float); //! Pz of the proton
DECLARE_SOA_COLUMN(PxGamma1, pxGamma1, float); //! Px of the measured (PCM) photon
DECLARE_SOA_COLUMN(PyGamma1, pyGamma1, float); //! Py of the measured (PCM) photon
DECLARE_SOA_COLUMN(PzGamma1, pzGamma1, float); //! Pz of the measured (PCM) photon
DECLARE_SOA_COLUMN(PxGamma2, pxGamma2, float); //! Px of the reconstructed (missing) photon
DECLARE_SOA_COLUMN(PyGamma2, pyGamma2, float); //! Py of the reconstructed (missing) photon
DECLARE_SOA_COLUMN(PzGamma2, pzGamma2, float); //! Pz of the reconstructed (missing) photon

DECLARE_SOA_COLUMN(NSigmaTPCProton, nSigmaTPCProton, float); //! TPC nSigma of the proton track
DECLARE_SOA_COLUMN(NSigmaTOFProton, nSigmaTOFProton, float); //! TOF nSigma of the proton track
DECLARE_SOA_COLUMN(NSigmaTPCElPos, nSigmaTPCElPos, float);   //! TPC nSigma_el of the photon's positive daughter
DECLARE_SOA_COLUMN(NSigmaTPCElNeg, nSigmaTPCElNeg, float);   //! TPC nSigma_el of the photon's negative daughter

DECLARE_SOA_COLUMN(PhotonMass, photonMass, float);                   //! Invariant mass of the measured photon (V0) candidate (GeV/c^2)
DECLARE_SOA_COLUMN(PhotonAlpha, photonAlpha, float);                 //! Armenteros-Podolanski alpha of the measured photon
DECLARE_SOA_COLUMN(PhotonQt, photonQt, float);                       //! Armenteros-Podolanski qT of the measured photon (GeV/c)
DECLARE_SOA_COLUMN(PhotonConvRadius, photonConvRadius, float);       //! Conversion radius of the measured photon (cm)
DECLARE_SOA_COLUMN(PhotonOpeningAngle, photonOpeningAngle, float);   //! Opening angle between the photon's e+e- daughters (rad)
DECLARE_SOA_COLUMN(PhotonPointingAngle, photonPointingAngle, float); //! Angle between the photon momentum and the line from its conversion point to the candidate decay vertex (rad)
DECLARE_SOA_COLUMN(PhotonDcaToPV, photonDcaToPV, float);             //! DCA of the photon's flight line to the primary vertex (cm)

DECLARE_SOA_COLUMN(ProtonItsNCls, protonItsNCls, uint8_t); //! Number of ITS clusters of the proton track
DECLARE_SOA_COLUMN(ProtonTpcNCls, protonTpcNCls, int16_t); //! Number of found TPC clusters of the proton track
DECLARE_SOA_COLUMN(ProtonDcaXY, protonDcaXY, float);       //! DCA of the proton track to the primary vertex, xy (cm)
DECLARE_SOA_COLUMN(ProtonDcaZ, protonDcaZ, float);         //! DCA of the proton track to the primary vertex, z (cm)

DECLARE_SOA_COLUMN(PhotonPosItsNCls, photonPosItsNCls, uint8_t); //! Number of ITS clusters of the photon's positive daughter
DECLARE_SOA_COLUMN(PhotonPosTpcNCls, photonPosTpcNCls, int16_t); //! Number of found TPC clusters of the photon's positive daughter
DECLARE_SOA_COLUMN(PhotonNegItsNCls, photonNegItsNCls, uint8_t); //! Number of ITS clusters of the photon's negative daughter
DECLARE_SOA_COLUMN(PhotonNegTpcNCls, photonNegTpcNCls, int16_t); //! Number of found TPC clusters of the photon's negative daughter

// MC columns
DECLARE_SOA_COLUMN(IsSignal, isSignal, bool); //! True if proton and photon are MC-truth matched to the same Sigma+

DECLARE_SOA_COLUMN(ProtonPdgCode, protonPdgCode, int);             //! PDG code of the proton's MC particle
DECLARE_SOA_COLUMN(ProtonMotherPdgCode, protonMotherPdgCode, int); //! PDG code of the proton's MC mother
DECLARE_SOA_COLUMN(GammaPdgCode, gammaPdgCode, int);               //! PDG code of the measured photon's MC particle
DECLARE_SOA_COLUMN(GammaMotherPdgCode, gammaMotherPdgCode, int);   //! PDG code of the photon's MC mother (expected: pi0)
DECLARE_SOA_COLUMN(GammaGMotherPdgCode, gammaGMotherPdgCode, int); //! PDG code of the photon's MC grandmother (expected: Sigma+)

DECLARE_SOA_COLUMN(XDecVtxMC, xDecVtxMC, float);   //! MC-truth Sigma+ decay vertex (x direction)
DECLARE_SOA_COLUMN(YDecVtxMC, yDecVtxMC, float);   //! MC-truth Sigma+ decay vertex (y direction)
DECLARE_SOA_COLUMN(ZDecVtxMC, zDecVtxMC, float);   //! MC-truth Sigma+ decay vertex (z direction)
DECLARE_SOA_COLUMN(PxProtonMC, pxProtonMC, float); //! MC-truth proton Px
DECLARE_SOA_COLUMN(PyProtonMC, pyProtonMC, float); //! MC-truth proton Py
DECLARE_SOA_COLUMN(PzProtonMC, pzProtonMC, float); //! MC-truth proton Pz
DECLARE_SOA_COLUMN(PxGammaMC, pxGammaMC, float);   //! MC-truth momentum of the measured photon (Px)
DECLARE_SOA_COLUMN(PyGammaMC, pyGammaMC, float);   //! MC-truth momentum of the measured photon (Py)
DECLARE_SOA_COLUMN(PzGammaMC, pzGammaMC, float);   //! MC-truth momentum of the measured photon (Pz)

// DYNAMIC COLUMNS

DECLARE_SOA_DYNAMIC_COLUMN(PxSigmaPlus, pxSigmaPlus, //! Px of the Sigma+ candidate
                           [](float pxProton, float pxGamma1, float pxGamma2) -> float { return pxProton + pxGamma1 + pxGamma2; });

DECLARE_SOA_DYNAMIC_COLUMN(PySigmaPlus, pySigmaPlus, //! Py of the Sigma+ candidate
                           [](float pyProton, float pyGamma1, float pyGamma2) -> float { return pyProton + pyGamma1 + pyGamma2; });

DECLARE_SOA_DYNAMIC_COLUMN(PzSigmaPlus, pzSigmaPlus, //! Pz of the Sigma+ candidate
                           [](float pzProton, float pzGamma1, float pzGamma2) -> float { return pzProton + pzGamma1 + pzGamma2; });

DECLARE_SOA_DYNAMIC_COLUMN(PtSigmaPlus, ptSigmaPlus, //! pT of the Sigma+ candidate
                           [](float pxProton, float pxGamma1, float pxGamma2, float pyProton, float pyGamma1, float pyGamma2) -> float { return std::hypot(pxProton + pxGamma1 + pxGamma2, pyProton + pyGamma1 + pyGamma2); });

DECLARE_SOA_DYNAMIC_COLUMN(MassSigmaPlus, massSigmaPlus, //! Invariant mass of the Sigma+ candidate (p + gamma1 + gamma2)
                           [](float pxProton, float pyProton, float pzProton,
                              float pxGamma1, float pyGamma1, float pzGamma1,
                              float pxGamma2, float pyGamma2, float pzGamma2) -> float { return RecoDecay::m(std::array{std::array{pxProton, pyProton, pzProton},
                                                                                                                        std::array{pxGamma1, pyGamma1, pzGamma1},
                                                                                                                        std::array{pxGamma2, pyGamma2, pzGamma2}},
                                                                                                             std::array{o2::constants::physics::MassProton, o2::constants::physics::MassGamma, o2::constants::physics::MassGamma}); });

} // namespace sigmapluscand

DECLARE_SOA_TABLE(SigmaPlusCands, "AOD", "SIGMAPLUSCANDS",
                  sigmapluscand::XDecVtx, sigmapluscand::YDecVtx, sigmapluscand::ZDecVtx,
                  sigmapluscand::Radius, sigmapluscand::FlightDistance,
                  sigmapluscand::DcaProtonGamma, sigmapluscand::Chi2,
                  sigmapluscand::PxProton, sigmapluscand::PyProton, sigmapluscand::PzProton,
                  sigmapluscand::PxGamma1, sigmapluscand::PyGamma1, sigmapluscand::PzGamma1,
                  sigmapluscand::PxGamma2, sigmapluscand::PyGamma2, sigmapluscand::PzGamma2,
                  sigmapluscand::NSigmaTPCProton, sigmapluscand::NSigmaTOFProton,
                  sigmapluscand::NSigmaTPCElPos, sigmapluscand::NSigmaTPCElNeg,
                  sigmapluscand::PhotonMass, sigmapluscand::PhotonAlpha, sigmapluscand::PhotonQt, sigmapluscand::PhotonConvRadius,
                  sigmapluscand::PhotonOpeningAngle, sigmapluscand::PhotonPointingAngle, sigmapluscand::PhotonDcaToPV,
                  sigmapluscand::ProtonItsNCls, sigmapluscand::ProtonTpcNCls, sigmapluscand::ProtonDcaXY, sigmapluscand::ProtonDcaZ,
                  sigmapluscand::PhotonPosItsNCls, sigmapluscand::PhotonPosTpcNCls, sigmapluscand::PhotonNegItsNCls, sigmapluscand::PhotonNegTpcNCls,

                  // dynamic columns
                  sigmapluscand::PxSigmaPlus<sigmapluscand::PxProton, sigmapluscand::PxGamma1, sigmapluscand::PxGamma2>,
                  sigmapluscand::PySigmaPlus<sigmapluscand::PyProton, sigmapluscand::PyGamma1, sigmapluscand::PyGamma2>,
                  sigmapluscand::PzSigmaPlus<sigmapluscand::PzProton, sigmapluscand::PzGamma1, sigmapluscand::PzGamma2>,
                  sigmapluscand::PtSigmaPlus<sigmapluscand::PxProton, sigmapluscand::PxGamma1, sigmapluscand::PxGamma2, sigmapluscand::PyProton, sigmapluscand::PyGamma1, sigmapluscand::PyGamma2>,
                  sigmapluscand::MassSigmaPlus<sigmapluscand::PxProton, sigmapluscand::PyProton, sigmapluscand::PzProton, sigmapluscand::PxGamma1, sigmapluscand::PyGamma1, sigmapluscand::PzGamma1, sigmapluscand::PxGamma2, sigmapluscand::PyGamma2, sigmapluscand::PzGamma2>);

DECLARE_SOA_TABLE(SigmaPlusCandsMC, "AOD", "SIGMAPLUSMC",
                  sigmapluscand::XDecVtx, sigmapluscand::YDecVtx, sigmapluscand::ZDecVtx,
                  sigmapluscand::Radius, sigmapluscand::FlightDistance,
                  sigmapluscand::DcaProtonGamma, sigmapluscand::Chi2,
                  sigmapluscand::PxProton, sigmapluscand::PyProton, sigmapluscand::PzProton,
                  sigmapluscand::PxGamma1, sigmapluscand::PyGamma1, sigmapluscand::PzGamma1,
                  sigmapluscand::PxGamma2, sigmapluscand::PyGamma2, sigmapluscand::PzGamma2,
                  sigmapluscand::NSigmaTPCProton, sigmapluscand::NSigmaTOFProton,
                  sigmapluscand::NSigmaTPCElPos, sigmapluscand::NSigmaTPCElNeg,
                  sigmapluscand::PhotonMass, sigmapluscand::PhotonAlpha, sigmapluscand::PhotonQt, sigmapluscand::PhotonConvRadius,
                  sigmapluscand::PhotonOpeningAngle, sigmapluscand::PhotonPointingAngle, sigmapluscand::PhotonDcaToPV,
                  sigmapluscand::ProtonItsNCls, sigmapluscand::ProtonTpcNCls, sigmapluscand::ProtonDcaXY, sigmapluscand::ProtonDcaZ,
                  sigmapluscand::PhotonPosItsNCls, sigmapluscand::PhotonPosTpcNCls, sigmapluscand::PhotonNegItsNCls, sigmapluscand::PhotonNegTpcNCls,
                  sigmapluscand::IsSignal,
                  sigmapluscand::ProtonPdgCode, sigmapluscand::ProtonMotherPdgCode,
                  sigmapluscand::GammaPdgCode, sigmapluscand::GammaMotherPdgCode, sigmapluscand::GammaGMotherPdgCode,
                  sigmapluscand::XDecVtxMC, sigmapluscand::YDecVtxMC, sigmapluscand::ZDecVtxMC,
                  sigmapluscand::PxProtonMC, sigmapluscand::PyProtonMC, sigmapluscand::PzProtonMC,
                  sigmapluscand::PxGammaMC, sigmapluscand::PyGammaMC, sigmapluscand::PzGammaMC,

                  // dynamic columns
                  sigmapluscand::PxSigmaPlus<sigmapluscand::PxProton, sigmapluscand::PxGamma1, sigmapluscand::PxGamma2>,
                  sigmapluscand::PySigmaPlus<sigmapluscand::PyProton, sigmapluscand::PyGamma1, sigmapluscand::PyGamma2>,
                  sigmapluscand::PzSigmaPlus<sigmapluscand::PzProton, sigmapluscand::PzGamma1, sigmapluscand::PzGamma2>,
                  sigmapluscand::PtSigmaPlus<sigmapluscand::PxProton, sigmapluscand::PxGamma1, sigmapluscand::PxGamma2, sigmapluscand::PyProton, sigmapluscand::PyGamma1, sigmapluscand::PyGamma2>,
                  sigmapluscand::MassSigmaPlus<sigmapluscand::PxProton, sigmapluscand::PyProton, sigmapluscand::PzProton, sigmapluscand::PxGamma1, sigmapluscand::PyGamma1, sigmapluscand::PzGamma1, sigmapluscand::PxGamma2, sigmapluscand::PyGamma2, sigmapluscand::PzGamma2>);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFKINKDECAYTABLES_H_
