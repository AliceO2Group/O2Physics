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

#ifndef PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_GLOBALCONSTANTS_H_
#define PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_GLOBALCONSTANTS_H_

const int gMaxCorrelator = 12;
const int gMaxHarmonic = 9;
const int gMaxIndex = 300;              // per order, used only in Test0
const int gMaxNoBinsKine = 1000;        // max number of bins for differential q-vector
const int gMaxBinsDiffWeights = 100;    // max number of bins for differential weights, see MakeWeights.C
const int gMaxNumberEtaSeparations = 9; // max number of different eta separations used to calculated 2p corr. with eta separations

#endif // PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_GLOBALCONSTANTS_H_
