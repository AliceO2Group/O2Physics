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

Configurable<float> Vz_min{"Vz_min", -10.0, "minimum vertex z range [cm]"};
Configurable<float> Vz_max{"Vz_max", 10.0, "maximum vertex z range [cm]"};
Configurable<float> pt_min{"pt_min", 0.2, "minimum track pt value [GeV/c]"};
Configurable<float> pt_max{"pt_max", 5.0, "maximum track pt value [GeV/c]"};
