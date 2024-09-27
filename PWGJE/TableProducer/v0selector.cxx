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

/// \brief Task to select V0s based on cuts
///
/// \author Gijs van Weelden <g.van.weelden@cern.ch>

#include <vector>
#include <string>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/RecoDecay.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct V0SelectorTask {
  // FIXME: Make the correct table
  // * Should contain flag that shows if V0 is selected or perhaps a value indicating K0S, Lambda, LambdaBar, Background
  Produces<V0SelectionTableType> v0SelectionTable{"v0SelectionTable", "V0 selection table"};

  Configurable<std::vector<double>> K0SPtBins{"K0SPtBins", {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0}, "K0S pt Vals"};
  Configurable<std::vector<double>> K0SRVals{"K0SRVals", {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, "K0S min R values"};
  Configurable<std::vector<double>> K0SCtauVals{"K0SCtauVals", {20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0}, "K0S max ctau values"};
  Configurable<std::vector<double>> K0SCosPAVals{"K0SCosPAVals", {0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997}, "K0S min cosPA values"};
  Configurable<std::vector<double>> K0SDCAVals{"K0SDCAVals", {0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.15, 0.15, 0.10, 0.10}, "K0S min DCA +- values"};
  Configurable<std::vector<double>> K0SDCAdVals{"K0SDCAdVals", {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, "K0S max DCAd values"};

  Configurable<std::vector<double>> LambdaPtBins{"LambdaPtBins", {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0}, "Lambda pt Vals"};
  Configurable<std::vector<double>> LambdaRVals{"LambdaRVals", {1.0, 10.0, 10.0, 10.0, 10.0, 10.0, 20.0, 20.0}, "Lambda min R values"};
  Configurable<std::vector<double>> LambdaCtauVals{"LambdaCtauVals", {22.5, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0}, "Lambda max ctau values"};
  Configurable<std::vector<double>> LambdaCosPAVals{"LambdaCosPAVals", {0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997}, "Lambda min cosPA values"};
  Configurable<std::vector<double>> LambdaDCAVals{"LambdaDCApVals", {0.20, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10}, "Lambda min DCA+ values"};
  Configurable<std::vector<double>> LambdaDCAVals{"LambdaDCAnVals", {0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.15, 0.15}, "Lambda min DCA- values"};
  Configurable<std::vector<double>> LambdaDCAdVals{"LambdaDCAdVals", {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, "Lambda max DCAd values"};

  Configurable<std::vector<double>> AntiLambdaPtBins{"AntiLambdaPtBins", {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0}, "AntiLambda pt Vals"};
  Configurable<std::vector<double>> AntiLambdaRVals{"AntiLambdaRVals", {10.0, 10.0, 10.0, 10.0, 20.0, 20.0, 20.0, 20.0}, "AntiLambda min R values"};
  Configurable<std::vector<double>> AntiLambdaCtauVals{"AntiLambdaCtauVals", {22.5, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0}, "AntiLambda max ctau values"};
  Configurable<std::vector<double>> AntiLambdaCosPAVals{"AntiLambdaCosPAVals", {0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997}, "AntiLambda min cosPA values"};
  Configurable<std::vector<double>> AntiLambdaDCAVals{"AntiLambdaDCApVals", {0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20}, "AntiLambda min DCA+ values"};
  Configurable<std::vector<double>> AntiLambdaDCAVals{"AntiLambdaDCAnVals", {0.20, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10}, "AntiLambda min DCA- values"};
  Configurable<std::vector<double>> AntiLambdaDCAdVals{"AntiLambdaDCAdVals", {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, "AntiLambda max DCAd values"};

  // FIXME: Use reasonable values
  // Can use these for defining mass signal region as well. Should this take a peak widdth or an upper and lower bound? Latter accounts for shifts in peak position
  Configurable<std::vector<double>> K0SCMCVals{"K0SCMCVals", {10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0}, "K0S competing mass cut values (MeV)"};
  Configurable<std::vector<double>> LambdaCMCVals{"LambdaCMCVals", {10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0}, "Lambda competing mass cut values (MeV)"};
  Configurable<std::vector<double>> AntiLambdaCMCVals{"AntiLambdaCMCVals", {10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0}, "AntiLambda competing mass cut values (MeV)"};

  void init(InitContext const&)
  {}

  template <typename T, typename U>
  bool K0SCuts(T const& collision, U const& v0)
  {
    if (v0.pt() < K0SPtBins[0] || v0.pt() > K0SPtBins[K0SPtBins.size() - 1]) {
      return false;
    }
    int ptBin = std::distance(K0SPtBins.begin(), std::upper_bound(K0SPtBins.begin(), K0SPtBins.end(), v0.pt())) - 1;
    if (v0.v0radius() < K0SRVals[ptBin]) {
      return false;
    }
    if (v0.v0cosPA() < K0SCosPAVals[ptBin]) {
      return false;
    }
    if (v0.dcaV0daughters() > K0SDCAdVals[ptBin]) {
      return false;
    }
    if (TMath::Abs(v0.dcapostopv()) < K0SDCAVals[ptBin]) {
      return false;
    }
    if (TMath::Abs(v0.dcanegtopv()) < K0SDCAVals[ptBin]) {
      return false;
    }
    double ctau = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short;
    if (ctau < K0SCtauVals[ptBin]) {
      return false;
    }
    // TODO: Add mass cut and competing mass cut (if requested)
    return true;
  }
  template <typename T, typename U>
  bool LambdaCuts(T const& collision, U const& v0)
  {
    if (v0.pt() < LambdaPtBins[0] || v0.pt() > LambdaPtBins[LambdaPtBins.size() - 1]) {
      return false;
    }
    int ptBin = std::distance(LambdaPtBins.begin(), std::upper_bound(LambdaPtBins.begin(), LambdaPtBins.end(), v0.pt())) - 1;
    if (v0.v0radius() < LambdaRVals[ptBin]) {
      return false;
    }
    if (v0.v0cosPA() < LambdaCosPAVals[ptBin]) {
      return false;
    }
    if (v0.dcaV0daughters() > LambdaDCAdVals[ptBin]) {
      return false;
    }
    if (TMath::Abs(v0.dcapostopv()) < LambdaDCApVals[ptBin]) {
      return false;
    }
    if (TMath::Abs(v0.dcanegtopv()) < LambdaDCAnVals[ptBin]) {
      return false;
    }
    double ctau = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
    if (ctau < LambdaCtauVals[ptBin]) {
      return false;
    }
    // TODO: Add mass cut and competing mass cut (if requested)
    return true;
  }
  template <typename T, typename U>
  bool AntiLambdaCuts(T const& collision, U const& v0)
  {
    if (v0.pt() < AntiLambdaPtBins[0] || v0.pt() > AntiLambdaPtBins[AntiLambdaPtBins.size() - 1]) {
      return false;
    }
    int ptBin = std::distance(AntiLambdaPtBins.begin(), std::upper_bound(AntiLambdaPtBins.begin(), AntiLambdaPtBins.end(), v0.pt())) - 1;
    if (v0.v0radius() < AntiLambdaRVals[ptBin]) {
      return false;
    }
    if (v0.v0cosPA() < AntiLambdaCosPAVals[ptBin]) {
      return false;
    }
    if (v0.dcaV0daughters() > AntiLambdaDCAdVals[ptBin]) {
      return false;
    }
    if (TMath::Abs(v0.dcapostopv()) < AntiLambdaDCApVals[ptBin]) {
      return false;
    }
    if (TMath::Abs(v0.dcanegtopv()) < AntiLambdaDCAnVals[ptBin]) {
      return false;
    }
    double ctau = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassAntiLambda0;
    if (ctau < AntiLambdaCtauVals[ptBin]) {
      return false;
    }
    // TODO: Add mass cut and competing mass cut (if requested)
    return true;
  }

  template <typename T, typename U>
  void processV0(T const& collision, U const& v0s)
  {
    bool flag = false;
    for (const auto& v0 : v0s) {
      bool candidateK0S = K0SCuts(collision, v0);
      bool candidateLambda = LambdaCuts(collision, v0);
      bool candidateAntiLambda = AntiLambdaCuts(collision, v0);
      bool flag = candidateK0S + candidateLambda + candidateAntiLambda;
      // Fill table
      v0SelectionTable(flag);
    }
  }
  PROCESS_SWITCH(V0SelectorTask, processV0, "flags V0 candidates potentially signal", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<V0SelectorTask>(cfgc, TaskName{"v0-selector"})};
}
