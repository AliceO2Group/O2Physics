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

/// \file ZDCExtra.h
/// \brief ZDC extra table producer
/// \author Chiara Oppedisano <chiara.oppedisano@cern.ch>, INFN Torino
/// \author Uliana Dmitrieva <uliana.dmitrieva@cern.ch>, INFN Torino

#ifndef COMMON_DATAMODEL_ZDCEXTRA_H_
#define COMMON_DATAMODEL_ZDCEXTRA_H_

#include <Framework/AnalysisDataModel.h>

#include <cstdint>

namespace o2::aod
{
namespace zdcextra
{
DECLARE_SOA_COLUMN(ZnaTowC, znaTowC, float);               //! Common tower ZNA
DECLARE_SOA_COLUMN(ZnaTow1, znaTow1, float);               //! Tower 1 ZNA
DECLARE_SOA_COLUMN(ZnaTow2, znaTow2, float);               //! Tower 2 ZNA
DECLARE_SOA_COLUMN(ZnaTow3, znaTow3, float);               //! Tower 3 ZNA
DECLARE_SOA_COLUMN(ZnaTow4, znaTow4, float);               //! Tower 4 ZNA
DECLARE_SOA_COLUMN(ZnaTdc, znaTdc, float);                 //! TDC ZNA
DECLARE_SOA_COLUMN(ZnaQx, znaQx, float);                   //! Q-vector X ZNA
DECLARE_SOA_COLUMN(ZnaQy, znaQy, float);                   //! Q-vector Y ZNA
DECLARE_SOA_COLUMN(ZncTowC, zncTowC, float);               //! Common tower ZNC
DECLARE_SOA_COLUMN(ZncTow1, zncTow1, float);               //! Tower 1 ZNC
DECLARE_SOA_COLUMN(ZncTow2, zncTow2, float);               //! Tower 2 ZNC
DECLARE_SOA_COLUMN(ZncTow3, zncTow3, float);               //! Tower 3 ZNC
DECLARE_SOA_COLUMN(ZncTow4, zncTow4, float);               //! Tower 4 ZNC
DECLARE_SOA_COLUMN(ZncTdc, zncTdc, float);                 //! TDC ZNC
DECLARE_SOA_COLUMN(ZncQx, zncQx, float);                   //! Q-vector X ZNC
DECLARE_SOA_COLUMN(ZncQy, zncQy, float);                   //! Q-vector Y ZNC
DECLARE_SOA_COLUMN(Centrality, centrality, float);         //! Centrality
DECLARE_SOA_COLUMN(Vx, vx, float);                         //! Vertex X
DECLARE_SOA_COLUMN(Vy, vy, float);                         //! Vertex Y
DECLARE_SOA_COLUMN(Vz, vz, float);                         //! Vertex Z
DECLARE_SOA_COLUMN(Timestamp, timestamp, uint64_t);        //! Timestamp
DECLARE_SOA_COLUMN(RunNumber, runNumber, uint32_t);        //! Run Number
DECLARE_SOA_COLUMN(SelectionBits, selectionBits, uint8_t); //! Selection Flags
} // namespace zdcextra

DECLARE_SOA_TABLE(ZdcExtras, "AOD", "ZDCEXTRA", o2::soa::Index<>,
                  zdcextra::ZnaTowC,
                  zdcextra::ZnaTow1,
                  zdcextra::ZnaTow2,
                  zdcextra::ZnaTow3,
                  zdcextra::ZnaTow4,
                  zdcextra::ZnaTdc,
                  zdcextra::ZnaQx,
                  zdcextra::ZnaQy,
                  zdcextra::ZncTowC,
                  zdcextra::ZncTow1,
                  zdcextra::ZncTow2,
                  zdcextra::ZncTow3,
                  zdcextra::ZncTow4,
                  zdcextra::ZncTdc,
                  zdcextra::ZncQx,
                  zdcextra::ZncQy,
                  zdcextra::Centrality,
                  zdcextra::Vx,
                  zdcextra::Vy,
                  zdcextra::Vz,
                  zdcextra::Timestamp,
                  zdcextra::RunNumber,
                  zdcextra::SelectionBits);
} // namespace o2::aod

#endif // COMMON_DATAMODEL_ZDCEXTRA_H_
