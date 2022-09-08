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
/// \file   GlobalRunningParameters.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Set of tables to carry the global information used when running e.g. magnetic field
///

#ifndef O2_ANALYSIS_GRP_H_
#define O2_ANALYSIS_GRP_H_

#include "Framework/AnalysisDataModel.h"
#include "DataFormatsParameters/GRPMagField.h"

namespace o2::aod
{
namespace grp
{
DECLARE_SOA_COLUMN(L3Current, l3Current, o2::units::Current_t);         //! Current in the L3 magnet -> Magnetic field in the barrel
DECLARE_SOA_COLUMN(DipoleCurrent, dipoleCurrent, o2::units::Current_t); //! Current in the dipole magnet -> Magnetic field in the dipole
DECLARE_SOA_COLUMN(FieldUniformity, fieldUniformity, bool);             //! Flag for uniformity of the field

template <typename BcType>
o2::parameters::GRPMagField getGrpMagField(const BcType& bc)
{
  o2::parameters::GRPMagField grpmag;
  grpmag.setL3Current(bc.l3Current());
  grpmag.setDipoleCurrent(bc.dipoleCurrent());
  grpmag.setFieldUniformity(bc.fieldUniformity());
  return grpmag;
}

template <typename BcType>
void initPropagator(const BcType& bc)
{
  o2::base::Propagator::initFieldFromGRP(&getGrpMagField(bc));
}

} // namespace grp

DECLARE_SOA_TABLE(GrpInfos, "AOD", "GRPINFOS", //! List of global parameters information that can be joined to BCs
                  grp::L3Current,
                  grp::DipoleCurrent,
                  grp::FieldUniformity);

} // namespace o2::aod

#endif