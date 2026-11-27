// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoProducerSigma0.cxx
/// \brief Derived-to-derived task: builds femto Sigma0 -> Lambda + gamma candidates from FLambdas and FPhotons
///
/// Needs the full-precision FLambdas, FLambdaMasks, FPhotons and FPhotonMasks as input (i.e. the
/// producer ran with LambdaTables.produceLambdas and PhotonTables.producePhotons, not the Lite variants,
/// or the Lite tables were converted back before running this task).
/// FSigma0Extras additionally needs FLambdaExtras and FPhotonExtras (processDataWithExtras).
/// TODO: FSigma0Labels: the Sigma0 currently gets the label of its Lambda daughter since photons are not implemented for femto MC.
/// \author Yakiv Paroviak, CTU in Prague, yakiv.paroviak@cern.ch

#include "PWGCF/Femto/Core/collisionBuilder.h"
#include "PWGCF/Femto/Core/partitions.h"
#include "PWGCF/Femto/Core/sigma0Builder.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

using namespace o2::analysis::femto;

struct FemtoProducerSigma0 {

  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;
  using FilteredFemtoCollisions = o2::soa::Filtered<FemtoCollisions>;
  using FilteredFemtoCollision = FilteredFemtoCollisions::iterator;

  using FemtoLambdas = o2::soa::Join<o2::aod::FLambdas, o2::aod::FLambdaMasks>;
  using FemtoPhotons = o2::soa::Join<o2::aod::FPhotons, o2::aod::FPhotonMasks>;
  using FemtoLambdasWithExtras = o2::soa::Join<FemtoLambdas, o2::aod::FLambdaExtras>;
  using FemtoPhotonsWithExtras = o2::soa::Join<FemtoPhotons, o2::aod::FPhotonExtras>;
  using FemtoLambdasWithLabel = o2::soa::Join<FemtoLambdas, o2::aod::FLambdaLabels>;
  using FemtoLambdasWithExtrasWithLabel = o2::soa::Join<FemtoLambdasWithExtras, o2::aod::FLambdaLabels>;

  o2::framework::SliceCache cache;

  // setup collisions
  collisionbuilder::ConfCollisionSelection collisionSelection;
  o2::framework::expressions::Filter collisionFilter = MAKE_COLLISION_FILTER(collisionSelection);

  // setup for Lambda daughters
  sigma0builder::ConfSigma0LambdaSelection confLambdaSelection;
  o2::framework::Partition<FemtoLambdas> lambdaPartition = MAKE_LAMBDA_PARTITION(confLambdaSelection);
  o2::framework::Partition<FemtoLambdasWithExtras> lambdaWithExtrasPartition = MAKE_LAMBDA_PARTITION(confLambdaSelection);
  o2::framework::Preslice<FemtoLambdas> perColLambdas = o2::aod::femtobase::stored::fColId;
  o2::framework::Preslice<FemtoLambdasWithExtras> perColLambdasWithExtras = o2::aod::femtobase::stored::fColId;
  o2::framework::Partition<FemtoLambdasWithLabel> lambdaWithLabelPartition = MAKE_LAMBDA_PARTITION(confLambdaSelection);
  o2::framework::Partition<FemtoLambdasWithExtrasWithLabel> lambdaWithExtrasWithLabelPartition = MAKE_LAMBDA_PARTITION(confLambdaSelection);
  o2::framework::Preslice<FemtoLambdasWithLabel> perColLambdasWithLabel = o2::aod::femtobase::stored::fColId;
  o2::framework::Preslice<FemtoLambdasWithExtrasWithLabel> perColLambdasWithExtrasWithLabel = o2::aod::femtobase::stored::fColId;

  // setup for photon daughters
  sigma0builder::ConfSigma0PhotonSelection confPhotonSelection;
  o2::framework::Partition<FemtoPhotons> photonPartition = MAKE_PHOTON_PARTITION(confPhotonSelection);
  o2::framework::Partition<FemtoPhotonsWithExtras> photonWithExtrasPartition = MAKE_PHOTON_PARTITION(confPhotonSelection);
  o2::framework::Preslice<FemtoPhotons> perColPhotons = o2::aod::femtobase::stored::fColId;
  o2::framework::Preslice<FemtoPhotonsWithExtras> perColPhotonsWithExtras = o2::aod::femtobase::stored::fColId;

  // Sigma0 builder
  sigma0builder::ConfSigma0Filters confSigma0Filter;
  sigma0builder::ConfSigma0Bits confSigma0Bits;
  sigma0builder::ConfSigma0Qa confSigma0Qa;
  sigma0builder::ConfSigma0Tables confSigma0Tables;
  sigma0builder::Sigma0BuilderProducts sigma0BuilderProducts;
  sigma0builder::Sigma0Builder sigma0Builder;

  // filter histograms of the input collisions, Lambdas and photons
  sigma0builder::Sigma0InputFilters sigma0InputFilters;

  o2::framework::HistogramRegistry hRegistry;

  void init(o2::framework::InitContext& context)
  {
    if (static_cast<int>(doprocessData) + static_cast<int>(doprocessDataWithExtras) + static_cast<int>(doprocessMc) + static_cast<int>(doprocessMcWithExtras) > 1) {
      LOG(fatal) << "Enable only one of processData, processDataWithExtras, processMc and processMcWithExtras";
    }
    sigma0InputFilters.init(&hRegistry, collisionSelection, confLambdaSelection, confPhotonSelection);
    sigma0Builder.init(&hRegistry, confLambdaSelection, confPhotonSelection, confSigma0Filter, confSigma0Bits, confSigma0Qa, confSigma0Tables, context);
  }

  // documents the filters applied to the input tables, independent of the process function used to build the Sigma0s
  void processInputFilters(FemtoCollisions const& cols, FemtoLambdas const& lambdas, FemtoPhotons const& photons)
  {
    sigma0InputFilters.fill(cols, lambdas, photons);
  }
  PROCESS_SWITCH(FemtoProducerSigma0, processInputFilters, "Fill the filter histograms of the input collisions, Lambdas and photons", true);

  void processData(FilteredFemtoCollision const& col, FemtoLambdas const& /*lambdas*/, FemtoPhotons const& /*photons*/)
  {
    sigma0Builder.fillSigma0s(col, sigma0BuilderProducts, lambdaPartition, photonPartition, cache);
  }
  PROCESS_SWITCH(FemtoProducerSigma0, processData, "Build Sigma0 candidates", true);

  // copies QA information of the Lambda and photon daughters
  void processDataWithExtras(FilteredFemtoCollision const& col, FemtoLambdasWithExtras const& /*lambdas*/, FemtoPhotonsWithExtras const& /*photons*/)
  {
    sigma0Builder.fillSigma0s(col, sigma0BuilderProducts, lambdaWithExtrasPartition, photonWithExtrasPartition, cache);
  }
  PROCESS_SWITCH(FemtoProducerSigma0, processDataWithExtras, "Build Sigma0 candidates and fill FSigma0Extras (needs FLambdaExtras and FPhotonExtras)", false);

  // the Sigma0 is labelled with the mc particle of its Lambda daughter
  void processMc(FilteredFemtoCollision const& col, FemtoLambdasWithLabel const& /*lambdas*/, FemtoPhotons const& /*photons*/)
  {
    sigma0Builder.fillSigma0s(col, sigma0BuilderProducts, lambdaWithLabelPartition, photonPartition, cache);
  }
  PROCESS_SWITCH(FemtoProducerSigma0, processMc, "Build Sigma0 candidates and fill FSigma0Labels (needs FLambdaLabels)", false);

  void processMcWithExtras(FilteredFemtoCollision const& col, FemtoLambdasWithExtrasWithLabel const& /*lambdas*/, FemtoPhotonsWithExtras const& /*photons*/)
  {
    sigma0Builder.fillSigma0s(col, sigma0BuilderProducts, lambdaWithExtrasWithLabelPartition, photonWithExtrasPartition, cache);
  }
  PROCESS_SWITCH(FemtoProducerSigma0, processMcWithExtras, "Build Sigma0 candidates and fill FSigma0Extras and FSigma0Labels (needs FLambdaExtras, FPhotonExtras and FLambdaLabels)", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{adaptAnalysisTask<FemtoProducerSigma0>(context)};
  return workflow;
}
