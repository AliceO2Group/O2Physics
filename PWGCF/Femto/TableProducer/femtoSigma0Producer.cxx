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

/// \file femtoSigma0Producer.cxx
/// \brief Derived-to-derived task: builds femto Sigma0 -> Lambda + gamma candidates from FLambdas and FPhotons
///
/// Needs the full-precision FLambdas, FLambdaMasks, FPhotons and FPhotonMasks as input (i.e. the
/// producer ran with LambdaTables.produceLambdas and PhotonTables.producePhotons, not the Lite variants,
/// or the Lite tables were converted back before running this task).
/// \author Yakiv Paroviak, CTU in Prague, yakiv.paroviak@cern.ch

#include "PWGCF/Femto/Core/sigma0Builder.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::analysis::femto;

struct FemtoSigma0Producer {
  using Lambdas = o2::soa::Join<o2::aod::FLambdas, o2::aod::FLambdaMasks>;
  using Photons = o2::soa::Join<o2::aod::FPhotons, o2::aod::FPhotonMasks>;

  sigma0builder::Sigma0BuilderProducts products;
  sigma0builder::ConfSigma0Filters filters;
  sigma0builder::ConfSigma0Bits bits;
  sigma0builder::ConfSigma0Tables tables;
  sigma0builder::Sigma0Builder sigma0Builder;

  HistogramRegistry registry{"Sigma0", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext& initContext)
  {
    sigma0Builder.init(&registry, filters, bits, tables, initContext);
  }

  // void process(Lambdas const& lambdas, Photons const& photons, o2::aod::FTracks const& tracks, o2::aod::FCols const& /**/, o2::aod::FColMasks const& /**/, o2::aod::FTrackMasks const& /**/)
  // {
  //   sigma0Builder.fillSigma0s(products, lambdas, photons, tracks);
  // }
  void process(Lambdas const& lambdas, Photons const& photons, o2::aod::FTracks const& tracks, o2::aod::FCols const& /**/, o2::aod::FColMasks const& /**/, o2::aod::FTrackMasks const& /**/)
  {
    sigma0Builder.fillSigma0s(products, lambdas, photons, tracks);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<FemtoSigma0Producer>(cfgc)};
}
