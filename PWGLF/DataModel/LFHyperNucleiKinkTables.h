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
/// \file LFHyperNucleiKinkTables.h
/// \brief Slim hypernuclei kink tables
/// \author Yuanzhe Wang <yuanzhe.wang@cern.ch>

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"

#ifndef PWGLF_DATAMODEL_LFHYPERNUCLEIKINKTABLES_H_
#define PWGLF_DATAMODEL_LFHYPERNUCLEIKINKTABLES_H_

namespace o2::aod
{

namespace hyperkink
{
DECLARE_SOA_COLUMN(MagPolarity, magPolarity, int8_t);                   //! Magnetic field polarity
DECLARE_SOA_COLUMN(XPV, xPV, float);                                    //! Primary vertex of the candidate (x direction)
DECLARE_SOA_COLUMN(YPV, yPV, float);                                    //! Primary vertex of the candidate (y direction)
DECLARE_SOA_COLUMN(ZPV, zPV, float);                                    //! Primary vertex of the candidate (z direction)
DECLARE_SOA_COLUMN(XSV, xSV, float);                                    //! Decay vertex of the candidate (x direction)
DECLARE_SOA_COLUMN(YSV, ySV, float);                                    //! Decay vertex of the candidate (y direction)
DECLARE_SOA_COLUMN(ZSV, zSV, float);                                    //! Decay vertex of the candidate (z direction)
DECLARE_SOA_COLUMN(XMothIU, xMothIU, float);                            //! X of the mother track at the radii of ITS layer which has the outermost update
DECLARE_SOA_COLUMN(YMothIU, yMothIU, float);                            //! Y of the mother track at the radii of ITS layer which has the outermost update
DECLARE_SOA_COLUMN(ZMothIU, zMothIU, float);                            //! Z of the mother track at the radii of ITS layer which has the outermost update
DECLARE_SOA_COLUMN(PxMothSV, pxMothSV, float);                          //! Px of the mother track at the decay vertex
DECLARE_SOA_COLUMN(PyMothSV, pyMothSV, float);                          //! Py of the mother track at the decay vertex
DECLARE_SOA_COLUMN(PzMothSV, pzMothSV, float);                          //! Pz of the mother track at the decay vertex
DECLARE_SOA_COLUMN(PxDaugSV, pxDaugSV, float);                          //! Px of the daughter track at the decay vertex
DECLARE_SOA_COLUMN(PyDaugSV, pyDaugSV, float);                          //! Py of the daughter track at the decay vertex
DECLARE_SOA_COLUMN(PzDaugSV, pzDaugSV, float);                          //! Pz of the daughter track at the decay vertex
DECLARE_SOA_COLUMN(IsMatter, isMatter, bool);                           //! bool: true for matter
DECLARE_SOA_COLUMN(DcaMothPv, dcaMothPv, float);                        //! DCA of the mother to the primary vertex
DECLARE_SOA_COLUMN(DcaDaugPv, dcaDaugPv, float);                        //! DCA of the daughter kink to the primary vertex
DECLARE_SOA_COLUMN(DcaKinkTopo, dcaKinkTopo, float);                    //! DCA of the kink topology
DECLARE_SOA_COLUMN(ItsChi2Moth, itsChi2Moth, float);                    //! ITS chi2 of the mother track
DECLARE_SOA_COLUMN(ItsClusterSizesMoth, itsClusterSizesMoth, uint32_t); //! ITS cluster size of the mother track
DECLARE_SOA_COLUMN(ItsClusterSizesDaug, itsClusterSizesDaug, uint32_t); //! ITS cluster size of the daughter track
DECLARE_SOA_COLUMN(TpcMomDaug, tpcMomDaug, float);                      //! TPC momentum of the daughter track
DECLARE_SOA_COLUMN(TpcSignalDaug, tpcSignalDaug, float);                //! TPC signal of the daughter track
DECLARE_SOA_COLUMN(TpcNClsPIDDaug, tpcNClsPIDDaug, int16_t);            //! Number of TPC clusters used for PID of the daughter track
DECLARE_SOA_COLUMN(NSigmaTPCDaug, nSigmaTPCDaug, float);                //! Number of tpc sigmas of the daughter track
DECLARE_SOA_COLUMN(NSigmaITSDaug, nSigmaITSDaug, float);                //! Number of ITS sigmas of the daughter track
DECLARE_SOA_COLUMN(NSigmaTOFDaug, nSigmaTOFDaug, float);                //! Number of TOF sigmas of the daughter track

DECLARE_SOA_COLUMN(IsSignal, isSignal, bool);                   //! bool: true for hyperhelium4signal
DECLARE_SOA_COLUMN(IsSignalReco, isSignalReco, bool);           //! bool: true if the signal is reconstructed
DECLARE_SOA_COLUMN(IsCollReco, isCollReco, bool);               //! bool: true if the collision is reconstructed
DECLARE_SOA_COLUMN(IsSurvEvSelection, isSurvEvSelection, bool); //! bool: true for the collision passed the event selection
DECLARE_SOA_COLUMN(TrueXSV, trueXSV, float);                    //! true x decay vertex
DECLARE_SOA_COLUMN(TrueYSV, trueYSV, float);                    //! true y decay vertex
DECLARE_SOA_COLUMN(TrueZSV, trueZSV, float);                    //! true z decay vertex
DECLARE_SOA_COLUMN(TruePxMothPV, truePxMothPV, float);          //! Generated px of the mother track
DECLARE_SOA_COLUMN(TruePyMothPV, truePyMothPV, float);          //! Generated py of the mother track
DECLARE_SOA_COLUMN(TruePzMothPV, truePzMothPV, float);          //! Generated pz of the mother track
DECLARE_SOA_COLUMN(TruePxMothSV, truePxMothSV, float);          //! true px of the mother track at the decay vertex
DECLARE_SOA_COLUMN(TruePyMothSV, truePyMothSV, float);          //! true py of the mother track at the decay vertex
DECLARE_SOA_COLUMN(TruePzMothSV, truePzMothSV, float);          //! true pz of the mother track at the decay vertex
DECLARE_SOA_COLUMN(TruePxDaugSV, truePxDaugSV, float);          //! true px of the daughter track at the decay vertex
DECLARE_SOA_COLUMN(TruePyDaugSV, truePyDaugSV, float);          //! true py of the daughter track at the decay vertex
DECLARE_SOA_COLUMN(TruePzDaugSV, truePzDaugSV, float);          //! true pz of the daughter track at the decay vertex
DECLARE_SOA_COLUMN(IsMothReco, isMothReco, bool);               //! bool: true if the mother track is reconstructed
DECLARE_SOA_COLUMN(PxMothPV, pxMothPV, float);                  //! reconstructed px of the mother track at the primary vertex
DECLARE_SOA_COLUMN(PyMothPV, pyMothPV, float);                  //! reconstructed py of the mother track at the primary vertex
DECLARE_SOA_COLUMN(PzMothPV, pzMothPV, float);                  //! reconstructed pz of the mother track at the primary vertex
DECLARE_SOA_COLUMN(UpdatePxMothPV, updatePxMothPV, float);      //! updated px of the mother track at the primary vertex after update using PV
DECLARE_SOA_COLUMN(UpdatePyMothPV, updatePyMothPV, float);      //! updated py of the mother track at the primary vertex after update using PV
DECLARE_SOA_COLUMN(UpdatePzMothPV, updatePzMothPV, float);      //! updated pz of the mother track at the primary vertex after update using PV
} // namespace hyperkink

DECLARE_SOA_TABLE(HypKinkCand, "AOD", "HYPKINKCANDS",
                  o2::soa::Index<>,
                  hyperkink::MagPolarity,
                  hyperkink::XPV, hyperkink::YPV, hyperkink::ZPV,
                  hyperkink::XSV, hyperkink::YSV, hyperkink::ZSV,
                  hyperkink::IsMatter,
                  hyperkink::XMothIU, hyperkink::YMothIU, hyperkink::ZMothIU,
                  hyperkink::PxMothSV, hyperkink::PyMothSV, hyperkink::PzMothSV,
                  hyperkink::PxDaugSV, hyperkink::PyDaugSV, hyperkink::PzDaugSV,
                  hyperkink::DcaMothPv, hyperkink::DcaDaugPv, hyperkink::DcaKinkTopo,
                  hyperkink::ItsChi2Moth, hyperkink::ItsClusterSizesMoth, hyperkink::ItsClusterSizesDaug,
                  hyperkink::TpcMomDaug, hyperkink::TpcSignalDaug, hyperkink::TpcNClsPIDDaug,
                  hyperkink::NSigmaTPCDaug, hyperkink::NSigmaITSDaug, hyperkink::NSigmaTOFDaug,
                  hyperkink::PxMothPV, hyperkink::PyMothPV, hyperkink::PzMothPV,
                  hyperkink::UpdatePxMothPV, hyperkink::UpdatePyMothPV, hyperkink::UpdatePzMothPV);

DECLARE_SOA_TABLE(MCHypKinkCand, "AOD", "MCHYPKINKCANDS",
                  o2::soa::Index<>,
                  hyperkink::MagPolarity,
                  hyperkink::XPV, hyperkink::YPV, hyperkink::ZPV,
                  hyperkink::XSV, hyperkink::YSV, hyperkink::ZSV,
                  hyperkink::IsMatter,
                  hyperkink::XMothIU, hyperkink::YMothIU, hyperkink::ZMothIU,
                  hyperkink::PxMothSV, hyperkink::PyMothSV, hyperkink::PzMothSV,
                  hyperkink::PxDaugSV, hyperkink::PyDaugSV, hyperkink::PzDaugSV,
                  hyperkink::DcaMothPv, hyperkink::DcaDaugPv, hyperkink::DcaKinkTopo,
                  hyperkink::ItsChi2Moth, hyperkink::ItsClusterSizesMoth, hyperkink::ItsClusterSizesDaug,
                  hyperkink::TpcMomDaug, hyperkink::TpcSignalDaug, hyperkink::TpcNClsPIDDaug,
                  hyperkink::NSigmaTPCDaug, hyperkink::NSigmaITSDaug, hyperkink::NSigmaTOFDaug,
                  hyperkink::IsSignal, hyperkink::IsSignalReco, hyperkink::IsCollReco, hyperkink::IsSurvEvSelection,
                  hyperkink::TrueXSV, hyperkink::TrueYSV, hyperkink::TrueZSV,
                  hyperkink::TruePxMothPV, hyperkink::TruePyMothPV, hyperkink::TruePzMothPV,
                  hyperkink::TruePxMothSV, hyperkink::TruePyMothSV, hyperkink::TruePzMothSV,
                  hyperkink::TruePxDaugSV, hyperkink::TruePyDaugSV, hyperkink::TruePzDaugSV,
                  hyperkink::IsMothReco, hyperkink::PxMothPV, hyperkink::PyMothPV, hyperkink::PzMothPV,
                  hyperkink::UpdatePxMothPV, hyperkink::UpdatePyMothPV, hyperkink::UpdatePzMothPV);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFHYPERNUCLEIKINKTABLES_H_
