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
/// \file   PIDResponse.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Set of tables, tasks and utilities to provide the interface between
///         the analysis data model and the PID response. This is interim. To be removed
///

#ifndef COMMON_DATAMODEL_PIDRESPONSE_H_
#define COMMON_DATAMODEL_PIDRESPONSE_H_

#include "Common/DataModel/PIDResponseTOF.h" // FIXME: remove
#include "Common/DataModel/PIDResponseTPC.h" // FIXME: remove

namespace o2::aod
{
namespace pidutils
{
// Function to pack a float into a binned value in table (interim solution)
template <typename binningType, typename T>
void packInTable(const float& valueToBin, T& table)
{
  binningType::packInTable(valueToBin, table);
}

// Function to unpack a binned value into a float
template <typename binningType>
float unPackInTable(const typename binningType::binned_t& valueToUnpack)
{
  return binningType::unPackInTable(valueToUnpack);
}
} // namespace pidutils
} // namespace o2::aod

#endif // COMMON_DATAMODEL_PIDRESPONSE_H_
