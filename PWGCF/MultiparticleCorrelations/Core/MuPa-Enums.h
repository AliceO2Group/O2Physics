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

enum eRecoSim { eRec = 0,
                eSim = 1 };

enum eBeforeAfter { eBefore = 0,
                    eAfter = 1 };

enum eDefaultColors { eColor = kBlack,
                      eFillColor = kGray };

enum eWeights { wPHI = 0,
                wPT = 1,
                wETA = 2,
                eWeights_N };

enum eEventHistograms { eNumberOfEvents,
                        eTotalMultiplicity,
                        eSelectedParticles,
                        eCentrality,
                        eVertex_x,
                        eVertex_y,
                        eVertex_z,
                        eEventHistograms_N };

enum eParticleHistograms { ePhi,
                           ePt,
                           eEta,
                           eParticleHistograms_N };
