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
/// \file LFHyperhelium4sigmaTables.h
/// \brief Slim hyperhelium4sigma tables
/// \author Yuanzhe Wang <yuanzhe.wang@cern.ch>

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"

#ifndef PWGLF_DATAMODEL_LFHYPERHELIUM4SIGMATABLES_H_
#define PWGLF_DATAMODEL_LFHYPERHELIUM4SIGMATABLES_H_

namespace o2::aod
{

namespace he4scand
{
DECLARE_SOA_COLUMN(XPrimVtx, xPrimVtx, float);                            // Primary vertex of the candidate (x direction)
DECLARE_SOA_COLUMN(YPrimVtx, yPrimVtx, float);                            // Primary vertex of the candidate (y direction)
DECLARE_SOA_COLUMN(ZPrimVtx, zPrimVtx, float);                            // Primary vertex of the candidate (z direction)
DECLARE_SOA_COLUMN(XDecVtx, xDecVtx, float);                              // Decay vertex of the candidate (x direction)
DECLARE_SOA_COLUMN(YDecVtx, yDecVtx, float);                              // Decay vertex of the candidate (y direction)
DECLARE_SOA_COLUMN(ZDecVtx, zDecVtx, float);                              // Decay vertex of the candidate (z direction)
DECLARE_SOA_COLUMN(XMoth, xMoth, float);                                  // X of the mother track at the radii of ITS layer which has the outermost update
DECLARE_SOA_COLUMN(YMoth, yMoth, float);                                  // Y of the mother track at the radii of ITS layer which has the outermost update
DECLARE_SOA_COLUMN(ZMoth, zMoth, float);                                  // Z of the mother track at the radii of ITS layer which has the outermost update
DECLARE_SOA_COLUMN(PxMoth, pxMoth, float);                                //! Px of the mother track at the decay vertex
DECLARE_SOA_COLUMN(PyMoth, pyMoth, float);                                //! Py of the mother track at the decay vertex
DECLARE_SOA_COLUMN(PzMoth, pzMoth, float);                                //! Pz of the mother track at the decay vertex
DECLARE_SOA_COLUMN(PxAlpha, pxAlpha, float);                              //! Px of the daughter alpha track at the decay vertex
DECLARE_SOA_COLUMN(PyAlpha, pyAlpha, float);                              //! Py of the daughter alpha track at the decay vertex
DECLARE_SOA_COLUMN(PzAlpha, pzAlpha, float);                              //! Pz of the daughter alpha track at the decay vertex
DECLARE_SOA_COLUMN(IsMatter, isMatter, bool);                             // bool: true for matter
DECLARE_SOA_COLUMN(DcaMothPv, dcaMothPv, float);                          //! DCA of the mother to the primary vertex
DECLARE_SOA_COLUMN(DcaAlphaPv, dcaAlphaPv, float);                        //! DCA of the daughter kink to the primary vertex
DECLARE_SOA_COLUMN(DcaKinkTopo, dcaKinkTopo, float);                      //! DCA of the kink topology
DECLARE_SOA_COLUMN(ItsChi2Moth, itsChi2Moth, float);                      // ITS chi2 of the mother track
DECLARE_SOA_COLUMN(ItsClusterSizesMoth, itsClusterSizesMoth, uint32_t);   // ITS cluster size of the mother track
DECLARE_SOA_COLUMN(ItsClusterSizesAlpha, itsClusterSizesAlpha, uint32_t); // ITS cluster size of the daughter alpha track
DECLARE_SOA_COLUMN(NSigmaTPCAlpha, nSigmaTPCAlpha, float);                // Number of tpc sigmas of the daughter alpha track
DECLARE_SOA_COLUMN(NSigmaITSAlpha, nSigmaITSAlpha, float);                // Number of ITS sigmas of the daughter alpha track

DECLARE_SOA_COLUMN(IsSignal, isSignal, bool);                   // bool: true for hyperhelium4signal
DECLARE_SOA_COLUMN(IsSignalReco, isSignalReco, bool);           // bool: true if the signal is reconstructed
DECLARE_SOA_COLUMN(IsCollReco, isCollReco, bool);               // bool: true if the collision is reconstructed
DECLARE_SOA_COLUMN(IsSurvEvSelection, isSurvEvSelection, bool); // bool: true for the collision passed the event selection
DECLARE_SOA_COLUMN(TrueXDecVtx, trueXDecVtx, float);            // true x decay vertex
DECLARE_SOA_COLUMN(TrueYDecVtx, trueYDecVtx, float);            // true y decay vertex
DECLARE_SOA_COLUMN(TrueZDecVtx, trueZDecVtx, float);            // true z decay vertex
DECLARE_SOA_COLUMN(GenPxMoth, genPxMoth, float);                // Generated px of the mother track
DECLARE_SOA_COLUMN(GenPyMoth, genPyMoth, float);                // Generated py of the mother track
DECLARE_SOA_COLUMN(GenPzMoth, genPzMoth, float);                // Generated pz of the mother track
DECLARE_SOA_COLUMN(TruePxMoth, truePxMoth, float);              // true px of the mother track at the decay vertex
DECLARE_SOA_COLUMN(TruePyMoth, truePyMoth, float);              // true py of the mother track at the decay vertex
DECLARE_SOA_COLUMN(TruePzMoth, truePzMoth, float);              // true pz of the mother track at the decay vertex
DECLARE_SOA_COLUMN(GenPxAlpha, genPxAlpha, float);              // true px of the daughter alpha track
DECLARE_SOA_COLUMN(GenPyAlpha, genPyAlpha, float);              // true py of the daughter alpha track
DECLARE_SOA_COLUMN(GenPzAlpha, genPzAlpha, float);              // true pz of the daughter alpha track
DECLARE_SOA_COLUMN(IsMothReco, isMothReco, bool);               // bool: true if the mother track is reconstructed
DECLARE_SOA_COLUMN(RecoPtMoth, recoPtMoth, float);              // reconstructed pt of the mother track
DECLARE_SOA_COLUMN(RecoPzMoth, recoPzMoth, float);              // reconstructed pz of the mother track
} // namespace he4scand

DECLARE_SOA_TABLE(He4S2BCands, "AOD", "HE4S2BCANDS",
                  o2::soa::Index<>,
                  he4scand::XPrimVtx, he4scand::YPrimVtx, he4scand::ZPrimVtx,
                  he4scand::XDecVtx, he4scand::YDecVtx, he4scand::ZDecVtx,
                  he4scand::IsMatter,
                  he4scand::XMoth, he4scand::YMoth, he4scand::ZMoth,
                  he4scand::PxMoth, he4scand::PyMoth, he4scand::PzMoth,
                  he4scand::PxAlpha, he4scand::PyAlpha, he4scand::PzAlpha,
                  he4scand::DcaMothPv, he4scand::DcaAlphaPv, he4scand::DcaKinkTopo,
                  he4scand::ItsChi2Moth, he4scand::ItsClusterSizesMoth, he4scand::ItsClusterSizesAlpha,
                  he4scand::NSigmaTPCAlpha, he4scand::NSigmaITSAlpha);

DECLARE_SOA_TABLE(MCHe4S2BCands, "AOD", "MCHE4S2BCANDS",
                  o2::soa::Index<>,
                  he4scand::XPrimVtx, he4scand::YPrimVtx, he4scand::ZPrimVtx,
                  he4scand::XDecVtx, he4scand::YDecVtx, he4scand::ZDecVtx,
                  he4scand::IsMatter,
                  he4scand::XMoth, he4scand::YMoth, he4scand::ZMoth,
                  he4scand::PxMoth, he4scand::PyMoth, he4scand::PzMoth,
                  he4scand::PxAlpha, he4scand::PyAlpha, he4scand::PzAlpha,
                  he4scand::DcaMothPv, he4scand::DcaAlphaPv, he4scand::DcaKinkTopo,
                  he4scand::ItsChi2Moth, he4scand::ItsClusterSizesMoth, he4scand::ItsClusterSizesAlpha,
                  he4scand::NSigmaTPCAlpha, he4scand::NSigmaITSAlpha,
                  he4scand::IsSignal, he4scand::IsSignalReco, he4scand::IsCollReco, he4scand::IsSurvEvSelection,
                  he4scand::TrueXDecVtx, he4scand::TrueYDecVtx, he4scand::TrueZDecVtx,
                  he4scand::GenPxMoth, he4scand::GenPyMoth, he4scand::GenPzMoth,
                  he4scand::TruePxMoth, he4scand::TruePyMoth, he4scand::TruePzMoth,
                  he4scand::GenPxAlpha, he4scand::GenPyAlpha, he4scand::GenPzAlpha,
                  he4scand::IsMothReco, he4scand::RecoPtMoth, he4scand::RecoPzMoth);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFHYPERHELIUM4SIGMATABLES_H_
