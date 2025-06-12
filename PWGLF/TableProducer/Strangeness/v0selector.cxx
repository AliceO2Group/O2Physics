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
#include <algorithm>
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
  Configurable<std::vector<float>> K0SRminVals{"K0SRminVals", {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, "K0S min R values"};
  Configurable<std::vector<float>> K0SRmaxVals{"K0SRmaxVals", {40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0}, "K0S max R values"};
  Configurable<std::vector<float>> K0SCtauminVals{"K0SCtauminVals", {-1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3}, "K0S min ctau values"};
  Configurable<std::vector<float>> K0SCtaumaxVals{"K0SCtaumaxVals", {20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0}, "K0S max ctau values"};
  Configurable<std::vector<float>> K0SCosPAminVals{"K0SCosPAminVals", {0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997}, "K0S min cosPA values"};
  Configurable<std::vector<float>> K0SCosPAmaxVals{"K0SCosPAmaxVals", {-1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3}, "K0S max cosPA values"};
  Configurable<std::vector<float>> K0SDCAminVals{"K0SDCAminVals", {0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.15, 0.15, 0.10, 0.10}, "K0S min DCA +- values"};
  Configurable<std::vector<float>> K0SDCAmaxVals{"K0SDCAmaxVals", {-1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3}, "K0S max DCA +- values"};
  Configurable<std::vector<float>> K0SDCAdminVals{"K0SDCAdminVals", {-1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3}, "K0S min DCAd values"};
  Configurable<std::vector<float>> K0SDCAdmaxVals{"K0SDCAdmaxVals", {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, "K0S max DCAd values"};

  Configurable<std::vector<float>> LambdaPtBins{"LambdaPtBins", {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0}, "Lambda pt Vals"};
  Configurable<std::vector<float>> LambdaRminVals{"LambdaRminVals", {1.0, 10.0, 10.0, 10.0, 10.0, 10.0, 20.0, 20.0}, "Lambda min R values"};
  Configurable<std::vector<float>> LambdaRmaxVals{"LambdaRmaxVals", {40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0}, "Lambda max R values"};
  Configurable<std::vector<float>> LambdaCtauminVals{"LambdaCtauminVals", {-1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3}, "Lambda min ctau values"};
  Configurable<std::vector<float>> LambdaCtaumaxVals{"LambdaCtaumaxVals", {22.5, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0}, "Lambda max ctau values"};
  Configurable<std::vector<float>> LambdaCosPAminVals{"LambdaCosPAminVals", {0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997}, "Lambda min cosPA values"};
  Configurable<std::vector<float>> LambdaCosPAmaxVals{"LambdaCosPAmaxVals", {-1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3}, "Lambda max cosPA values"};
  Configurable<std::vector<float>> LambdaDCApminVals{"LambdaDCApminVals", {0.20, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10}, "Lambda min DCA+ values"};
  Configurable<std::vector<float>> LambdaDCApmaxVals{"LambdaDCApmaxVals", {-1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3}, "Lambda max DCA+ values"};
  Configurable<std::vector<float>> LambdaDCAnminVals{"LambdaDCAnminVals", {0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.15, 0.15}, "Lambda min DCA- values"};
  Configurable<std::vector<float>> LambdaDCAnmaxVals{"LambdaDCAnmaxVals", {-1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3}, "Lambda max DCA- values"};
  Configurable<std::vector<float>> LambdaDCAdminVals{"LambdaDCAdminVals", {-1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3}, "Lambda max DCAd values"};
  Configurable<std::vector<float>> LambdaDCAdmaxVals{"LambdaDCAdmaxVals", {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, "Lambda max DCAd values"};

  Configurable<std::vector<float>> AntiLambdaPtBins{"AntiLambdaPtBins", {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 15.0, 20.0}, "AntiLambda pt Vals"};
  Configurable<std::vector<float>> AntiLambdaRminVals{"AntiLambdaRminVals", {10.0, 10.0, 10.0, 10.0, 20.0, 20.0, 20.0, 20.0}, "AntiLambda min R values"};
  Configurable<std::vector<float>> AntiLambdaRmaxVals{"AntiLambdaRmaxVals", {40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0}, "AntiLambda max R values"};
  Configurable<std::vector<float>> AntiLambdaCtauminVals{"AntiLambdaCtauminVals", {-1e-3, -1e-3, -1e-3, -1e-3, -1e-3, -1e-3, -1e-3, -1e-3}, "AntiLambda min ctau values"};
  Configurable<std::vector<float>> AntiLambdaCtaumaxVals{"AntiLambdaCtaumaxVals", {22.5, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0, 25.0}, "AntiLambda max ctau values"};
  Configurable<std::vector<float>> AntiLambdaCosPAminVals{"AntiLambdaCosPAminVals", {0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997, 0.997}, "AntiLambda min cosPA values"};
  Configurable<std::vector<float>> AntiLambdaCosPAmaxVals{"AntiLambdaCosPAmaxVals", {-1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3}, "AntiLambda max cosPA values"};
  Configurable<std::vector<float>> AntiLambdaDCApminVals{"AntiLambdaDCApminVals", {0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20}, "AntiLambda min DCA+ values"};
  Configurable<std::vector<float>> AntiLambdaDCApmaxVals{"AntiLambdaDCApmaxVals", {-1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3}, "AntiLambda max DCA+ values"};
  Configurable<std::vector<float>> AntiLambdaDCAnminVals{"AntiLambdaDCAnminVals", {0.20, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10}, "AntiLambda min DCA- values"};
  Configurable<std::vector<float>> AntiLambdaDCAnmaxVals{"AntiLambdaDCAnmaxVals", {-1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3}, "AntiLambda max DCA- values"};
  Configurable<std::vector<float>> AntiLambdaDCAdminVals{"AntiLambdaDCAdminVals", {-1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3, -1e3}, "AntiLambda min DCAd values"};
  Configurable<std::vector<float>> AntiLambdaDCAdmaxVals{"AntiLambdaDCAdmaxVals", {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, "AntiLambda max DCAd values"};

  Configurable<bool> massCuts{"massCuts", true, "Apply mass cuts"};
  Configurable<bool> competingMassCuts{"competingMassCuts", true, "Apply competing mass cuts"};
  Configurable<std::vector<float>> K0SMassLowVals{"K0SMassLowVals", {0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4}, "K0S mass cut lower values (MeV)"};
  Configurable<std::vector<float>> K0SMassHighVals{"K0SMassHighVals", {0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6}, "K0S mass cut upper values (MeV)"};
  Configurable<std::vector<float>> LambdaMassLowVals{"LambdaMassLowVals", {1.08, 1.08, 1.08, 1.08, 1.08, 1.08, 1.08, 1.08}, "Lambda mass cut lower values (MeV)"};
  Configurable<std::vector<float>> LambdaMassHighVals{"LambdaMassHighVals", {1.125, 1.125, 1.125, 1.125, 1.125, 1.125, 1.125, 1.125}, "Lambda mass cut upper values (MeV)"};
  Configurable<std::vector<float>> AntiLambdaMassLowVals{"AntiLambdaMassLowVals", {1.08, 1.08, 1.08, 1.08, 1.08, 1.08, 1.08, 1.08}, "AntiLambda mass cut lower values (MeV)"};
  Configurable<std::vector<float>> AntiLambdaMassHighVals{"AntiLambdaMassHighVals", {1.125, 1.125, 1.125, 1.125, 1.125, 1.125, 1.125, 1.125}, "AntiLambda mass cut upper values (MeV)"};

  Configurable<bool> randomSelection{"randomSelection", true, "Randomly select V0s"};
  Configurable<std::vector<float>> K0SFraction{"K0SFraction", {2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0}, "Fraction of K0S to randomly select"};
  Configurable<std::vector<float>> LambdaFraction{"LambdaFraction", {2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0}, "Fraction of Lambda to randomly select"};
  Configurable<std::vector<float>> AntiLambdaFraction{"AntiLambdaFraction", {2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0}, "Fraction of AntiLambda to randomly select"};

  void init(InitContext const&)
  {
  }

  int getPtBin(float pt, std::vector<float> ptBins)
  {
    if (pt < ptBins.at(0))
      return -1;
    if (pt > ptBins.at(ptBins.size() - 1))
      return -2;

    for (unsigned int i = 0; i < ptBins.size() - 1; i++) {
      if (pt >= ptBins.at(i) && pt < ptBins.at(i + 1)) {
        return i;
      }
    }
    return -3;
  }
  bool cuts(std::vector<float> values, std::vector<float> mincuts, std::vector<float> maxcuts)
  {
    for (unsigned int i = 0; i < values.size(); i++) {
      float val = values.at(i);
      float min = mincuts.at(i);
      float max = maxcuts.at(i);

      if (val < min && min > -1e2)
        return false;
      if (val > max && max > -1e2)
        return false;
    }
    return true;
  }
  template <typename T, typename U>
  bool K0SCuts(T const& collision, U const& v0)
  {
    int ptBin = getPtBin(v0.pt(), K0SPtBins);
    if (ptBin < 0)
      return false;

    // This is the only time we need to check min and max simultaneously
    // K0S and Lambda(bar) do not share pt binning, so check the ptBin for Lambda(bar) separately
    if (competingMassCuts) {
      int ptBinCMC = getPtBin(v0.pt(), LambdaPtBins);
      if (ptBinCMC >= 0) { // TODO: Should we still do CMC when v0 is out of pT range for Lambda(bar)?
        if (v0.mLambda() > LambdaMassLowVals->at(ptBinCMC) && v0.mLambda() < LambdaMassHighVals->at(ptBinCMC)) {
          return false;
        }
      }
      ptBinCMC = getPtBin(v0.pt(), AntiLambdaPtBins);
      if (ptBinCMC >= 0) {
        if (v0.mAntiLambda() > AntiLambdaMassLowVals->at(ptBinCMC) && v0.mAntiLambda() < AntiLambdaMassHighVals->at(ptBinCMC)) {
          return false;
        }
      }
    }

    float rmin = K0SRminVals->at(ptBin);
    float rmax = K0SRmaxVals->at(ptBin);
    float ctaumin = K0SCtauminVals->at(ptBin);
    float ctaumax = K0SCtaumaxVals->at(ptBin);
    float cospamin = K0SCosPAminVals->at(ptBin);
    float cospamax = K0SCosPAmaxVals->at(ptBin);
    float dcapmin = K0SDCAminVals->at(ptBin);
    float dcapmax = K0SDCAmaxVals->at(ptBin);
    float dcanmin = K0SDCAminVals->at(ptBin);
    float dcanmax = K0SDCAmaxVals->at(ptBin);
    float dcadmin = K0SDCAdminVals->at(ptBin);
    float dcadmax = K0SDCAdmaxVals->at(ptBin);

    float ctau = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short;
    std::vector<float> vals = {v0.v0radius(), ctau, v0.v0cosPA(), v0.dcaV0daughters(), TMath::Abs(v0.dcapostopv()), TMath::Abs(v0.dcanegtopv())};
    std::vector<float> mincuts = {rmin, ctaumin, cospamin, dcadmin, dcapmin, dcanmin};
    std::vector<float> maxcuts = {rmax, ctaumax, cospamax, dcadmax, dcapmax, dcanmax};

    if (massCuts) {
      vals.push_back(v0.mK0Short());
      mincuts.push_back(K0SMassLowVals->at(ptBin));
      maxcuts.push_back(K0SMassHighVals->at(ptBin));
    }

    return cuts(vals, mincuts, maxcuts);
  }
  template <typename T, typename U>
  bool LambdaCuts(T const& collision, U const& v0)
  {
    int ptBin = getPtBin(v0.pt(), LambdaPtBins);
    if (ptBin < 0)
      return false;

    if (competingMassCuts) {
      int ptBinCMC = getPtBin(v0.pt(), K0SPtBins);
      if (ptBinCMC >= 0) {
        if (v0.mK0Short() > K0SMassLowVals->at(ptBinCMC) && v0.mK0Short() < K0SMassHighVals->at(ptBinCMC))
          return false;
      }
    }

    float rmin = LambdaRminVals->at(ptBin);
    float rmax = LambdaRmaxVals->at(ptBin);
    float ctaumin = LambdaCtauminVals->at(ptBin);
    float ctaumax = LambdaCtaumaxVals->at(ptBin);
    float cospamin = LambdaCosPAminVals->at(ptBin);
    float cospamax = LambdaCosPAmaxVals->at(ptBin);
    float dcapmin = LambdaDCApminVals->at(ptBin);
    float dcapmax = LambdaDCApmaxVals->at(ptBin);
    float dcanmin = LambdaDCAnminVals->at(ptBin);
    float dcanmax = LambdaDCAnmaxVals->at(ptBin);
    float dcadmin = LambdaDCAdminVals->at(ptBin);
    float dcadmax = LambdaDCAdmaxVals->at(ptBin);

    float ctau = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
    std::vector<float> vals = {v0.v0radius(), ctau, v0.v0cosPA(), v0.dcaV0daughters(), TMath::Abs(v0.dcapostopv()), TMath::Abs(v0.dcanegtopv())};
    std::vector<float> mincuts = {rmin, ctaumin, cospamin, dcadmin, dcapmin, dcanmin};
    std::vector<float> maxcuts = {rmax, ctaumax, cospamax, dcadmax, dcapmax, dcanmax};

    if (massCuts) {
      vals.push_back(v0.mLambda());
      mincuts.push_back(LambdaMassLowVals->at(ptBin));
      maxcuts.push_back(LambdaMassHighVals->at(ptBin));
    }
    return cuts(vals, mincuts, maxcuts);
  }
  template <typename T, typename U>
  bool AntiLambdaCuts(T const& collision, U const& v0)
  {
    int ptBin = getPtBin(v0.pt(), AntiLambdaPtBins);
    if (ptBin < 0)
      return false;

    if (competingMassCuts) {
      int ptBinCMC = getPtBin(v0.pt(), K0SPtBins);
      if (ptBinCMC >= 0) {
        if (v0.mK0Short() > K0SMassLowVals->at(ptBinCMC) && v0.mK0Short() < K0SMassHighVals->at(ptBinCMC))
          return false;
      }
    }

    float rmin = AntiLambdaRminVals->at(ptBin);
    float rmax = AntiLambdaRmaxVals->at(ptBin);
    float ctaumin = AntiLambdaCtauminVals->at(ptBin);
    float ctaumax = AntiLambdaCtaumaxVals->at(ptBin);
    float cospamin = AntiLambdaCosPAminVals->at(ptBin);
    float cospamax = AntiLambdaCosPAmaxVals->at(ptBin);
    float dcapmin = AntiLambdaDCApminVals->at(ptBin);
    float dcapmax = AntiLambdaDCApmaxVals->at(ptBin);
    float dcanmin = AntiLambdaDCAnminVals->at(ptBin);
    float dcanmax = AntiLambdaDCAnmaxVals->at(ptBin);
    float dcadmin = AntiLambdaDCAdminVals->at(ptBin);
    float dcadmax = AntiLambdaDCAdmaxVals->at(ptBin);

    float ctau = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
    std::vector<float> vals = {v0.v0radius(), ctau, v0.v0cosPA(), v0.dcaV0daughters(), TMath::Abs(v0.dcapostopv()), TMath::Abs(v0.dcanegtopv())};
    std::vector<float> mincuts = {rmin, ctaumin, cospamin, dcadmin, dcapmin, dcanmin};
    std::vector<float> maxcuts = {rmax, ctaumax, cospamax, dcadmax, dcapmax, dcanmax};

    if (competingMassCuts) {
      vals.push_back(v0.mAntiLambda());
      mincuts.push_back(AntiLambdaMassLowVals->at(ptBin));
      maxcuts.push_back(AntiLambdaMassHighVals->at(ptBin));
    }
    return cuts(vals, mincuts, maxcuts);
  }
  template <typename T>
  bool RandomlyReject(T const& v0, uint8_t flag)
  {
    // In case of multiple candidate types, only check the lowest threshold value
    float threshold = 2.;
    if (flag & aod::v0flags::FK0S) {
      int ptBin = getPtBin(v0.pt(), K0SPtBins);
      threshold = std::min(threshold, K0SFraction->at(ptBin));
    }
    if (flag & aod::v0flags::FLAMBDA) {
      int ptBin = getPtBin(v0.pt(), LambdaPtBins);
      threshold = std::min(threshold, LambdaFraction->at(ptBin));
    }
    if (flag & aod::v0flags::FANTILAMBDA) {
      int ptBin = getPtBin(v0.pt(), AntiLambdaPtBins);
      threshold = std::min(threshold, AntiLambdaFraction->at(ptBin));
    }
    // If gRandom > threshold, reject the candidate
    return (gRandom->Uniform() > threshold);
  }

  void processV0(aod::Collision const& collision, aod::V0Datas const& v0s)
  {
    for (const auto& v0 : v0s) {
      uint8_t flag = 0;
      flag += K0SCuts(collision, v0) * aod::v0flags::FK0S;
      flag += LambdaCuts(collision, v0) * aod::v0flags::FLAMBDA;
      flag += AntiLambdaCuts(collision, v0) * aod::v0flags::FANTILAMBDA;

      if (flag == 0)
        flag += aod::v0flags::FREJECTED;
      else
        flag += RandomlyReject(v0, flag) * aod::v0flags::FREJECTED;

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
