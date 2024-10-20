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
#include <TRandom3.h>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/RecoDecay.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/V0SelectorTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct V0SelectorTask {
  Produces<aod::V0SignalFlags> v0FlagTable;

  Configurable<std::vector<float>> K0SPtBins{"K0SPtBins", {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0}, "K0S pt Vals"};
  Configurable<std::vector<float>> K0SRVals{"K0SRVals", {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, "K0S min R values"};
  Configurable<std::vector<float>> K0SCtauVals{"K0SCtauVals", {20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0}, "K0S max ctau values"};
  Configurable<std::vector<float>> K0SCosPAVals{"K0SCosPAVals", {0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997}, "K0S min cosPA values"};
  Configurable<std::vector<float>> K0SDCAVals{"K0SDCAVals", {0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.15, 0.15, 0.10, 0.10}, "K0S min DCA +- values"};
  Configurable<std::vector<float>> K0SDCAdVals{"K0SDCAdVals", {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, "K0S max DCAd values"};

  Configurable<std::vector<float>> LambdaPtBins{"LambdaPtBins", {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0}, "Lambda pt Vals"};
  Configurable<std::vector<float>> LambdaRVals{"LambdaRVals", {1.0, 10.0, 10.0, 10.0, 10.0, 10.0, 20.0, 20.0}, "Lambda min R values"};
  Configurable<std::vector<float>> LambdaCtauVals{"LambdaCtauVals", {22.5, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0}, "Lambda max ctau values"};
  Configurable<std::vector<float>> LambdaCosPAVals{"LambdaCosPAVals", {0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997}, "Lambda min cosPA values"};
  Configurable<std::vector<float>> LambdaDCApVals{"LambdaDCApVals", {0.20, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10}, "Lambda min DCA+ values"};
  Configurable<std::vector<float>> LambdaDCAnVals{"LambdaDCAnVals", {0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.15, 0.15}, "Lambda min DCA- values"};
  Configurable<std::vector<float>> LambdaDCAdVals{"LambdaDCAdVals", {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, "Lambda max DCAd values"};

  Configurable<std::vector<float>> AntiLambdaPtBins{"AntiLambdaPtBins", {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0}, "AntiLambda pt Vals"};
  Configurable<std::vector<float>> AntiLambdaRVals{"AntiLambdaRVals", {10.0, 10.0, 10.0, 10.0, 20.0, 20.0, 20.0, 20.0}, "AntiLambda min R values"};
  Configurable<std::vector<float>> AntiLambdaCtauVals{"AntiLambdaCtauVals", {22.5, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0}, "AntiLambda max ctau values"};
  Configurable<std::vector<float>> AntiLambdaCosPAVals{"AntiLambdaCosPAVals", {0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997}, "AntiLambda min cosPA values"};
  Configurable<std::vector<float>> AntiLambdaDCApVals{"AntiLambdaDCApVals", {0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20}, "AntiLambda min DCA+ values"};
  Configurable<std::vector<float>> AntiLambdaDCAnVals{"AntiLambdaDCAnVals", {0.20, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10}, "AntiLambda min DCA- values"};
  Configurable<std::vector<float>> AntiLambdaDCAdVals{"AntiLambdaDCAdVals", {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, "AntiLambda max DCAd values"};

  Configurable<bool> massCuts{"massCuts", true, "Apply mass cuts"};
  Configurable<bool> competingMassCuts{"competingMassCuts", true, "Apply competing mass cuts"};
  Configurable<std::vector<float>> K0SMassLowVals{"K0SMassLowVals", {0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4}, "K0S mass cut lower values (MeV)"};
  Configurable<std::vector<float>> K0SMassHighVals{"K0SMassHighVals", {0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6}, "K0S mass cut upper values (MeV)"};
  Configurable<std::vector<float>> LambdaMassLowVals{"LambdaMassLowVals", {1.08, 1.08, 1.08, 1.08, 1.08, 1.08, 1.08, 1.08}, "Lambda mass cut lower values (MeV)"};
  Configurable<std::vector<float>> LambdaMassHighVals{"LambdaMassHighVals", {1.125, 1.125, 1.125, 1.125, 1.125, 1.125, 1.125, 1.125}, "Lambda mass cut upper values (MeV)"};
  Configurable<std::vector<float>> AntiLambdaMassLowVals{"AntiLambdaMassLowVals", {1.08, 1.08, 1.08, 1.08, 1.08, 1.08, 1.08, 1.08}, "AntiLambda mass cut lower values (MeV)"};
  Configurable<std::vector<float>> AntiLambdaMassHighVals{"AntiLambdaMassHighVals", {1.125, 1.125, 1.125, 1.125, 1.125, 1.125, 1.125, 1.125}, "AntiLambda mass cut upper values (MeV)"};

  Configurable<bool> randomSelection{"randomSelection", true, "Randomly select V0s"};
  Configurable<std::vector<float>> K0SFraction{"randomSelectionFraction", {2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0}, "Fraction of K0S to randomly select"};
  Configurable<std::vector<float>> LambdaFraction{"LambdaFraction", {2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0}, "Fraction of Lambda to randomly select"};
  Configurable<std::vector<float>> AntiLambdaFraction{"AntiLambdaFraction", {2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0}, "Fraction of AntiLambda to randomly select"};

  void init(InitContext const&)
  {
  }

  template <typename T, typename U>
  bool K0SCuts(T const& collision, U const& v0)
  {
    if (v0.pt() < K0SPtBins->at(0) || v0.pt() > K0SPtBins->at(K0SPtBins->size() - 1)) {
      return false;
    }
    int ptBin = std::distance(K0SPtBins->begin(), std::upper_bound(K0SPtBins->begin(), K0SPtBins->end(), v0.pt())) - 1;
    if (v0.v0radius() < K0SRVals->at(ptBin)) {
      return false;
    }
    if (v0.v0cosPA() < K0SCosPAVals->at(ptBin)) {
      return false;
    }
    if (v0.dcaV0daughters() > K0SDCAdVals->at(ptBin)) {
      return false;
    }
    if (TMath::Abs(v0.dcapostopv()) < K0SDCAVals->at(ptBin)) {
      return false;
    }
    if (TMath::Abs(v0.dcanegtopv()) < K0SDCAVals->at(ptBin)) {
      return false;
    }
    float ctau = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short;
    if (ctau < K0SCtauVals->at(ptBin)) {
      return false;
    }
    // Apply mass cuts only if requested
    if (!massCuts) {
      return true;
    }
    if (v0.mK0Short() < K0SMassLowVals->at(ptBin) || v0.mK0Short() > K0SMassHighVals->at(ptBin)) {
      return false;
    }
    if (!competingMassCuts) {
      return true;
    }
    if (v0.mLambda() > LambdaMassLowVals->at(ptBin) && v0.mLambda() < LambdaMassHighVals->at(ptBin)) {
      return false;
    }
    if (v0.mAntiLambda() > AntiLambdaMassLowVals->at(ptBin) && v0.mAntiLambda() < AntiLambdaMassHighVals->at(ptBin)) {
      return false;
    }
    return true;
  }
  template <typename T, typename U>
  bool LambdaCuts(T const& collision, U const& v0)
  {
    if (v0.pt() < LambdaPtBins->at(0) || v0.pt() > LambdaPtBins->at(LambdaPtBins->size() - 1)) {
      return false;
    }
    int ptBin = std::distance(LambdaPtBins->begin(), std::upper_bound(LambdaPtBins->begin(), LambdaPtBins->end(), v0.pt())) - 1;
    if (v0.v0radius() < LambdaRVals->at(ptBin)) {
      return false;
    }
    if (v0.v0cosPA() < LambdaCosPAVals->at(ptBin)) {
      return false;
    }
    if (v0.dcaV0daughters() > LambdaDCAdVals->at(ptBin)) {
      return false;
    }
    if (TMath::Abs(v0.dcapostopv()) < LambdaDCApVals->at(ptBin)) {
      return false;
    }
    if (TMath::Abs(v0.dcanegtopv()) < LambdaDCAnVals->at(ptBin)) {
      return false;
    }
    float ctau = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
    if (ctau < LambdaCtauVals->at(ptBin)) {
      return false;
    }
    // Apply mass cuts only if requested
    if (!massCuts) {
      return true;
    }
    if (v0.mLambda() < LambdaMassLowVals->at(ptBin) || v0.mLambda() > LambdaMassHighVals->at(ptBin)) {
      return false;
    }
    if (!competingMassCuts) {
      return true;
    }
    if (v0.mK0Short() > K0SMassLowVals->at(ptBin) && v0.mK0Short() < K0SMassHighVals->at(ptBin)) {
      return false;
    }
    return true;
  }
  template <typename T, typename U>
  bool AntiLambdaCuts(T const& collision, U const& v0)
  {
    if (v0.pt() < AntiLambdaPtBins->at(0) || v0.pt() > AntiLambdaPtBins->at(AntiLambdaPtBins->size() - 1)) {
      return false;
    }
    int ptBin = std::distance(AntiLambdaPtBins->begin(), std::upper_bound(AntiLambdaPtBins->begin(), AntiLambdaPtBins->end(), v0.pt())) - 1;
    if (v0.v0radius() < AntiLambdaRVals->at(ptBin)) {
      return false;
    }
    if (v0.v0cosPA() < AntiLambdaCosPAVals->at(ptBin)) {
      return false;
    }
    if (v0.dcaV0daughters() > AntiLambdaDCAdVals->at(ptBin)) {
      return false;
    }
    if (TMath::Abs(v0.dcapostopv()) < AntiLambdaDCApVals->at(ptBin)) {
      return false;
    }
    if (TMath::Abs(v0.dcanegtopv()) < AntiLambdaDCAnVals->at(ptBin)) {
      return false;
    }
    float ctau = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
    if (ctau < AntiLambdaCtauVals->at(ptBin)) {
      return false;
    }
    // Apply mass cuts only if requested
    if (!massCuts) {
      return true;
    }
    if (v0.mAntiLambda() < AntiLambdaMassLowVals->at(ptBin) || v0.mAntiLambda() > AntiLambdaMassHighVals->at(ptBin)) {
      return false;
    }
    if (!competingMassCuts) {
      return true;
    }
    if (v0.mK0Short() > K0SMassLowVals->at(ptBin) && v0.mK0Short() < K0SMassHighVals->at(ptBin)) {
      return false;
    }
    return true;
  }
  template <typename T>
  bool RandomlyReject(T const& v0, uint8_t flag)
  {
    if (!(flag & aod::v0flags::FK0S || flag & aod::v0flags::FLAMBDA || flag & aod::v0flags::FANTILAMBDA)) {
      return true;
    }

    // In case of multiple candidate types, only check the lowest threshold value
    float threshold = 2.;
    if (flag & aod::v0flags::FK0S) {
      int ptBin = std::distance(K0SPtBins->begin(), std::upper_bound(K0SPtBins->begin(), K0SPtBins->end(), v0.pt())) - 1;
      if (threshold > K0SFraction->at(ptBin)) {
        threshold = K0SFraction->at(ptBin);
      }
    }
    if (flag & aod::v0flags::FLAMBDA) {
      int ptBin = std::distance(LambdaPtBins->begin(), std::upper_bound(LambdaPtBins->begin(), LambdaPtBins->end(), v0.pt())) - 1;
      if (threshold > LambdaFraction->at(ptBin)) {
        threshold = LambdaFraction->at(ptBin);
      }
    }
    if (flag & aod::v0flags::FANTILAMBDA) {
      int ptBin = std::distance(AntiLambdaPtBins->begin(), std::upper_bound(AntiLambdaPtBins->begin(), AntiLambdaPtBins->end(), v0.pt())) - 1;
      if (threshold > AntiLambdaFraction->at(ptBin)) {
        threshold = AntiLambdaFraction->at(ptBin);
      }
    }
    return (gRandom->Uniform() > threshold);
  }

  void processV0(aod::Collision const& collision, aod::V0Datas const& v0s)
  {
    for (const auto& v0 : v0s) {
      bool candidateK0S = K0SCuts(collision, v0);
      bool candidateLambda = LambdaCuts(collision, v0);
      bool candidateAntiLambda = AntiLambdaCuts(collision, v0);
      uint8_t flag = 0;
      flag += candidateK0S * aod::v0flags::FK0S;
      flag += candidateLambda * aod::v0flags::FLAMBDA;
      flag += candidateAntiLambda * aod::v0flags::FANTILAMBDA;

      if (candidateK0S + candidateLambda + candidateAntiLambda == 0) {
        flag += aod::v0flags::FREJECTED;
      } else if (randomSelection) {
        flag += RandomlyReject(v0, flag) * aod::v0flags::FREJECTED;
      }
      v0FlagTable(flag);
    }
  }
  PROCESS_SWITCH(V0SelectorTask, processV0, "flags V0 candidates as potential signal", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<V0SelectorTask>(cfgc, TaskName{"v0-selector"})};
}
