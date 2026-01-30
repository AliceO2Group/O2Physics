// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file jetSubstructureHFPCH.h
/// \brief Precompiled header for jet substructure HF tasks. Reduces compilation time significantly.
///
/// \author Sergio Garcia <sergio.garcia@cern.ch>

#ifndef PWGJE_TASKS_JETSUBSTRUCTUREHFPCH_H_
#define PWGJE_TASKS_JETSUBSTRUCTUREHFPCH_H_

// Heavy PWGJE DataModel headers
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/DataModel/JetSubtraction.h"

// Main template header
#include "PWGJE/Tasks/jetSubstructureHF.h"

// Framework headers
#include <Framework/ASoA.h>
#include <Framework/AnalysisTask.h>
#include <Framework/ConfigContext.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/runDataProcessing.h>

// Standard library
#include <vector>

#endif // PWGJE_TASKS_JETSUBSTRUCTUREHFPCH_H_
