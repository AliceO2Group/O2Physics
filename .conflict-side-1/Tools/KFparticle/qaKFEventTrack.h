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

/// \file qaKFEventTrack.h
/// \author Annalena Kalteyer <annalena.sophie.kalteyer@cern.ch>

#ifndef TOOLS_KFPARTICLE_QAKFEVENTTRACK_H_
#define TOOLS_KFPARTICLE_QAKFEVENTTRACK_H_

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/trackUtilities.h"

enum FlagsTracks {
  kITS = BIT(0),
  kTPC = BIT(1),
  kTRD = BIT(2),
  kTOF = BIT(3),
};

namespace o2::aod
{
namespace kfeventtrack
{
DECLARE_SOA_COLUMN(FLAGDetector, FlagsDetector, float);
DECLARE_SOA_COLUMN(FLAGPVContrib, FlagPVContrib, float);
DECLARE_SOA_COLUMN(DCAXYAOD, DcaXYAod, float);
DECLARE_SOA_COLUMN(DCAZAOD, DcaZAod, float);
DECLARE_SOA_COLUMN(DCA3DKF, Dca3dkf, float);
DECLARE_SOA_COLUMN(DCAXYKF, DcaXYkf, float);
DECLARE_SOA_COLUMN(CHARGE, Charge, float);
DECLARE_SOA_COLUMN(MOMENTUM, Momentum, float);
DECLARE_SOA_COLUMN(ETA, Eta, float);
DECLARE_SOA_COLUMN(PHI, Phi, float);
DECLARE_SOA_COLUMN(TPCSIGNAL, Tpcsignal, float);
DECLARE_SOA_COLUMN(RUNNUMBER, Runnumber, float);
DECLARE_SOA_COLUMN(TIMECOLL, TimeColl, uint64_t);
DECLARE_SOA_COLUMN(TIMESTAMP, TimeStamp, double);
DECLARE_SOA_COLUMN(TIMEDIFF, TimeDiff, double);
DECLARE_SOA_COLUMN(BCID, BCid, int);
DECLARE_SOA_COLUMN(TFID, Tfid, int);
DECLARE_SOA_COLUMN(XPV, Xpv, float);
DECLARE_SOA_COLUMN(YPV, Ypv, float);
DECLARE_SOA_COLUMN(ZPV, Zpv, float);
DECLARE_SOA_COLUMN(COVXX, CovXX, float);
DECLARE_SOA_COLUMN(COVYY, CovYY, float);
DECLARE_SOA_COLUMN(COVZZ, CovZZ, float);
DECLARE_SOA_COLUMN(NCONTRIB, Ncontrib, float);
DECLARE_SOA_COLUMN(NTRACKS, Ntracks, int);
DECLARE_SOA_COLUMN(CHI2, Chi2, float);

} // namespace kfeventtrack

DECLARE_SOA_TABLE(TreeTracks, "AOD", "TREETRACKS",
                  kfeventtrack::FLAGDetector,
                  kfeventtrack::FLAGPVContrib,
                  kfeventtrack::DCAXYAOD,
                  kfeventtrack::DCAZAOD,
                  kfeventtrack::DCA3DKF,
                  kfeventtrack::DCAXYKF,
                  kfeventtrack::CHARGE,
                  kfeventtrack::MOMENTUM,
                  kfeventtrack::ETA,
                  kfeventtrack::PHI,
                  kfeventtrack::TPCSIGNAL,
                  kfeventtrack::RUNNUMBER,
                  kfeventtrack::XPV,
                  kfeventtrack::YPV,
                  kfeventtrack::ZPV,
                  kfeventtrack::COVXX,
                  kfeventtrack::COVYY,
                  kfeventtrack::COVZZ,
                  kfeventtrack::NCONTRIB,
                  kfeventtrack::CHI2);

DECLARE_SOA_TABLE(TreeCollisions, "AOD", "TREECOLLISIONS",
                  kfeventtrack::XPV,
                  kfeventtrack::YPV,
                  kfeventtrack::ZPV,
                  kfeventtrack::COVXX,
                  kfeventtrack::COVYY,
                  kfeventtrack::COVZZ,
                  kfeventtrack::NCONTRIB,
                  kfeventtrack::NTRACKS,
                  kfeventtrack::CHI2,
                  kfeventtrack::RUNNUMBER,
                  kfeventtrack::TIMECOLL,
                  kfeventtrack::TIMESTAMP,
                  kfeventtrack::TIMEDIFF,
                  kfeventtrack::BCID,
                  kfeventtrack::TFID);
} // namespace o2::aod

#endif // TOOLS_KFPARTICLE_QAKFEVENTTRACK_H_
