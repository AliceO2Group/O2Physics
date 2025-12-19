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
/// \file   OTFTOF.h
/// \author David Dobrigkeit Chinellato
/// \author Nicolo Jacazio
/// \since  11/05/2023
/// \brief  Set of tables for the ALICE3 OTFTOF information
///

#ifndef ALICE3_DATAMODEL_OTFTOF_H_
#define ALICE3_DATAMODEL_OTFTOF_H_

// O2 includes
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace upgrade_tof
{
DECLARE_SOA_COLUMN(InnerTOFTrackTime, innerTOFTrackTime, float);     //! Track time generated at the InnerTOF
DECLARE_SOA_COLUMN(InnerTOFTrackLength, innerTOFTrackLength, float); //! track length for calculation of InnerTOF (generated)
DECLARE_SOA_COLUMN(OuterTOFTrackTime, outerTOFTrackTime, float);     //! Track time generated at the OuterTOF
DECLARE_SOA_COLUMN(OuterTOFTrackLength, outerTOFTrackLength, float); //! track length for calculation of OuterTOF (generated)

DECLARE_SOA_COLUMN(TOFEventTime, tofEventTime, float);                       //! Event time reconstructed with the TOF
DECLARE_SOA_COLUMN(TOFEventTimeErr, tofEventTimeErr, float);                 //! Uncertainty on the event time reconstructed with the TOF
DECLARE_SOA_COLUMN(NSigmaElectronInnerTOF, nSigmaElectronInnerTOF, float);   //! NSigma electron InnerTOF
DECLARE_SOA_COLUMN(NSigmaMuonInnerTOF, nSigmaMuonInnerTOF, float);           //! NSigma muon InnerTOF
DECLARE_SOA_COLUMN(NSigmaPionInnerTOF, nSigmaPionInnerTOF, float);           //! NSigma pion InnerTOF
DECLARE_SOA_COLUMN(NSigmaKaonInnerTOF, nSigmaKaonInnerTOF, float);           //! NSigma kaon InnerTOF
DECLARE_SOA_COLUMN(NSigmaProtonInnerTOF, nSigmaProtonInnerTOF, float);       //! NSigma proton InnerTOF
DECLARE_SOA_COLUMN(NSigmaDeuteronInnerTOF, nSigmaDeuteronInnerTOF, float);   //! NSigma deuteron InnerTOF
DECLARE_SOA_COLUMN(NSigmaTritonInnerTOF, nSigmaTritonInnerTOF, float);       //! NSigma triton InnerTOF
DECLARE_SOA_COLUMN(NSigmaHelium3InnerTOF, nSigmaHelium3InnerTOF, float);     //! NSigma helium3 InnerTOF
DECLARE_SOA_COLUMN(NSigmaAlphaInnerTOF, nSigmaAlphaInnerTOF, float);         //! NSigma alpha InnerTOF
DECLARE_SOA_COLUMN(InnerTOFTrackTimeReco, innerTOFTrackTimeReco, float);     //! Track time measured at the InnerTOF
DECLARE_SOA_COLUMN(InnerTOFTrackLengthReco, innerTOFTrackLengthReco, float); //! track length for calculation of InnerTOF (reconstructed)

DECLARE_SOA_COLUMN(InnerTOFExpectedTimeEl, innerTOFExpectedTimeEl, float);   //! Reconstructed expected time at the InnerTOF for the Electron mass hypotheses
DECLARE_SOA_COLUMN(InnerTOFExpectedTimeMu, innerTOFExpectedTimeMu, float);   //! Reconstructed expected time at the InnerTOF for the Muon mass hypotheses
DECLARE_SOA_COLUMN(InnerTOFExpectedTimePi, innerTOFExpectedTimePi, float);   //! Reconstructed expected time at the InnerTOF for the Pion mass hypotheses
DECLARE_SOA_COLUMN(InnerTOFExpectedTimeKa, innerTOFExpectedTimeKa, float);   //! Reconstructed expected time at the InnerTOF for the Kaon mass hypotheses
DECLARE_SOA_COLUMN(InnerTOFExpectedTimePr, innerTOFExpectedTimePr, float);   //! Reconstructed expected time at the InnerTOF for the Proton mass hypotheses
DECLARE_SOA_COLUMN(InnerTOFExpectedTimeDe, innerTOFExpectedTimeDe, float);   //! Reconstructed expected time at the InnerTOF for the Deuteron mass hypotheses
DECLARE_SOA_COLUMN(InnerTOFExpectedTimeTr, innerTOFExpectedTimeTr, float);   //! Reconstructed expected time at the InnerTOF for the Triton mass hypotheses
DECLARE_SOA_COLUMN(InnerTOFExpectedTimeHe3, innerTOFExpectedTimeHe3, float); //! Reconstructed expected time at the InnerTOF for the Helium3 mass hypotheses
DECLARE_SOA_COLUMN(InnerTOFExpectedTimeAl, innerTOFExpectedTimeAl, float);   //! Reconstructed expected time at the InnerTOF for the Alpha mass hypotheses

DECLARE_SOA_COLUMN(NSigmaElectronOuterTOF, nSigmaElectronOuterTOF, float);   //! NSigma electron OuterTOF
DECLARE_SOA_COLUMN(NSigmaMuonOuterTOF, nSigmaMuonOuterTOF, float);           //! NSigma muon OuterTOF
DECLARE_SOA_COLUMN(NSigmaPionOuterTOF, nSigmaPionOuterTOF, float);           //! NSigma pion OuterTOF
DECLARE_SOA_COLUMN(NSigmaKaonOuterTOF, nSigmaKaonOuterTOF, float);           //! NSigma kaon OuterTOF
DECLARE_SOA_COLUMN(NSigmaProtonOuterTOF, nSigmaProtonOuterTOF, float);       //! NSigma proton OuterTOF
DECLARE_SOA_COLUMN(NSigmaDeuteronOuterTOF, nSigmaDeuteronOuterTOF, float);   //! NSigma deuteron OuterTOF
DECLARE_SOA_COLUMN(NSigmaTritonOuterTOF, nSigmaTritonOuterTOF, float);       //! NSigma triton OuterTOF
DECLARE_SOA_COLUMN(NSigmaHelium3OuterTOF, nSigmaHelium3OuterTOF, float);     //! NSigma helium3 OuterTOF
DECLARE_SOA_COLUMN(NSigmaAlphaOuterTOF, nSigmaAlphaOuterTOF, float);         //! NSigma alpha OuterTOF
DECLARE_SOA_COLUMN(OuterTOFTrackTimeReco, outerTOFTrackTimeReco, float);     //! Track time measured at the OuterTOF
DECLARE_SOA_COLUMN(OuterTOFTrackLengthReco, outerTOFTrackLengthReco, float); //! track length for calculation of OuterTOF (reconstructed)

DECLARE_SOA_COLUMN(OuterTOFExpectedTimeEl, outerTOFExpectedTimeEl, float);   //! Reconstructed expected time at the OuterTOF for the Electron mass hypotheses
DECLARE_SOA_COLUMN(OuterTOFExpectedTimeMu, outerTOFExpectedTimeMu, float);   //! Reconstructed expected time at the OuterTOF for the Muon mass hypotheses
DECLARE_SOA_COLUMN(OuterTOFExpectedTimePi, outerTOFExpectedTimePi, float);   //! Reconstructed expected time at the OuterTOF for the Pion mass hypotheses
DECLARE_SOA_COLUMN(OuterTOFExpectedTimeKa, outerTOFExpectedTimeKa, float);   //! Reconstructed expected time at the OuterTOF for the Kaon mass hypotheses
DECLARE_SOA_COLUMN(OuterTOFExpectedTimePr, outerTOFExpectedTimePr, float);   //! Reconstructed expected time at the OuterTOF for the Proton mass hypotheses
DECLARE_SOA_COLUMN(OuterTOFExpectedTimeDe, outerTOFExpectedTimeDe, float);   //! Reconstructed expected time at the OuterTOF for the Deuteron mass hypotheses
DECLARE_SOA_COLUMN(OuterTOFExpectedTimeTr, outerTOFExpectedTimeTr, float);   //! Reconstructed expected time at the OuterTOF for the Triton mass hypotheses
DECLARE_SOA_COLUMN(OuterTOFExpectedTimeHe3, outerTOFExpectedTimeHe3, float); //! Reconstructed expected time at the OuterTOF for the Helium3 mass hypotheses
DECLARE_SOA_COLUMN(OuterTOFExpectedTimeAl, outerTOFExpectedTimeAl, float);   //! Reconstructed expected time at the OuterTOF for the Alpha mass hypotheses
DECLARE_SOA_DYNAMIC_COLUMN(NSigmaInnerTOF, nSigmaInnerTOF,                   //! General function to get the nSigma for the InnerTOF
                           [](const float el,
                              const float mu,
                              const float pi,
                              const float ka,
                              const float pr,
                              const float de,
                              const float tr,
                              const float he3,
                              const float al,
                              const int id) -> float {
                             switch (std::abs(id)) {
                               case 0:
                                 return el;
                               case 1:
                                 return mu;
                               case 2:
                                 return pi;
                               case 3:
                                 return ka;
                               case 4:
                                 return pr;
                               case 5:
                                 return de;
                               case 6:
                                 return tr;
                               case 7:
                                 return he3;
                               case 8:
                                 return al;
                               default:
                                 LOG(fatal) << "Unrecognized PDG code for InnerTOF";
                                 return 999.f;
                             }
                           });
DECLARE_SOA_DYNAMIC_COLUMN(NSigmaOuterTOF, nSigmaOuterTOF, //! General function to get the nSigma for the OuterTOF
                           [](const float el,
                              const float mu,
                              const float pi,
                              const float ka,
                              const float pr,
                              const float de,
                              const float tr,
                              const float he3,
                              const float al,
                              const int id) -> float {
                             switch (std::abs(id)) {
                               case 0:
                                 return el;
                               case 1:
                                 return mu;
                               case 2:
                                 return pi;
                               case 3:
                                 return ka;
                               case 4:
                                 return pr;
                               case 5:
                                 return de;
                               case 6:
                                 return tr;
                               case 7:
                                 return he3;
                               case 8:
                                 return al;
                               default:
                                 LOG(fatal) << "Unrecognized PDG code for InnerTOF";
                                 return 999.f;
                             }
                           });

} // namespace upgrade_tof

DECLARE_SOA_TABLE(UpgradeTofMCs, "AOD", "UPGRADETOFMC",
                  upgrade_tof::InnerTOFTrackTime,
                  upgrade_tof::InnerTOFTrackLength,
                  upgrade_tof::OuterTOFTrackTime,
                  upgrade_tof::OuterTOFTrackLength);

DECLARE_SOA_TABLE(UpgradeTofs, "AOD", "UPGRADETOF",
                  upgrade_tof::TOFEventTime,
                  upgrade_tof::TOFEventTimeErr,
                  upgrade_tof::NSigmaElectronInnerTOF,
                  upgrade_tof::NSigmaMuonInnerTOF,
                  upgrade_tof::NSigmaPionInnerTOF,
                  upgrade_tof::NSigmaKaonInnerTOF,
                  upgrade_tof::NSigmaProtonInnerTOF,
                  upgrade_tof::NSigmaDeuteronInnerTOF,
                  upgrade_tof::NSigmaTritonInnerTOF,
                  upgrade_tof::NSigmaHelium3InnerTOF,
                  upgrade_tof::NSigmaAlphaInnerTOF,
                  upgrade_tof::InnerTOFTrackTimeReco,
                  upgrade_tof::InnerTOFTrackLengthReco,
                  upgrade_tof::NSigmaElectronOuterTOF,
                  upgrade_tof::NSigmaMuonOuterTOF,
                  upgrade_tof::NSigmaPionOuterTOF,
                  upgrade_tof::NSigmaKaonOuterTOF,
                  upgrade_tof::NSigmaProtonOuterTOF,
                  upgrade_tof::NSigmaDeuteronOuterTOF,
                  upgrade_tof::NSigmaTritonOuterTOF,
                  upgrade_tof::NSigmaHelium3OuterTOF,
                  upgrade_tof::NSigmaAlphaOuterTOF,
                  upgrade_tof::OuterTOFTrackTimeReco,
                  upgrade_tof::OuterTOFTrackLengthReco,
                  upgrade_tof::NSigmaInnerTOF<upgrade_tof::NSigmaElectronInnerTOF,
                                              upgrade_tof::NSigmaMuonInnerTOF,
                                              upgrade_tof::NSigmaPionInnerTOF,
                                              upgrade_tof::NSigmaKaonInnerTOF,
                                              upgrade_tof::NSigmaProtonInnerTOF,
                                              upgrade_tof::NSigmaDeuteronInnerTOF,
                                              upgrade_tof::NSigmaTritonInnerTOF,
                                              upgrade_tof::NSigmaHelium3InnerTOF,
                                              upgrade_tof::NSigmaAlphaInnerTOF>,
                  upgrade_tof::NSigmaOuterTOF<upgrade_tof::NSigmaElectronOuterTOF,
                                              upgrade_tof::NSigmaMuonOuterTOF,
                                              upgrade_tof::NSigmaPionOuterTOF,
                                              upgrade_tof::NSigmaKaonOuterTOF,
                                              upgrade_tof::NSigmaProtonOuterTOF,
                                              upgrade_tof::NSigmaDeuteronOuterTOF,
                                              upgrade_tof::NSigmaTritonOuterTOF,
                                              upgrade_tof::NSigmaHelium3OuterTOF,
                                              upgrade_tof::NSigmaAlphaOuterTOF>);

DECLARE_SOA_TABLE(UpgradeTofExpectedTimes, "AOD", "UPGRADETOFEXPT",
                  upgrade_tof::InnerTOFExpectedTimeEl,
                  upgrade_tof::InnerTOFExpectedTimeMu,
                  upgrade_tof::InnerTOFExpectedTimePi,
                  upgrade_tof::InnerTOFExpectedTimeKa,
                  upgrade_tof::InnerTOFExpectedTimePr,
                  upgrade_tof::InnerTOFExpectedTimeDe,
                  upgrade_tof::InnerTOFExpectedTimeTr,
                  upgrade_tof::InnerTOFExpectedTimeHe3,
                  upgrade_tof::InnerTOFExpectedTimeAl,
                  upgrade_tof::OuterTOFExpectedTimeEl,
                  upgrade_tof::OuterTOFExpectedTimeMu,
                  upgrade_tof::OuterTOFExpectedTimePi,
                  upgrade_tof::OuterTOFExpectedTimeKa,
                  upgrade_tof::OuterTOFExpectedTimePr,
                  upgrade_tof::OuterTOFExpectedTimeDe,
                  upgrade_tof::OuterTOFExpectedTimeTr,
                  upgrade_tof::OuterTOFExpectedTimeHe3,
                  upgrade_tof::OuterTOFExpectedTimeAl);

using UpgradeTofMC = UpgradeTofMCs::iterator;
using UpgradeTof = UpgradeTofs::iterator;
using UpgradeTofExpectedTime = UpgradeTofExpectedTimes::iterator;

} // namespace o2::aod

#endif // ALICE3_DATAMODEL_OTFTOF_H_
