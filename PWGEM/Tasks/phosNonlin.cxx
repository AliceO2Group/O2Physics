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

/// \file phosNonlin.cxx
/// \brief task to calculate PHOS non-lienarity based on pi0 peak position
/// \author Dmitri Peresunko <Dmitri.Peresunko@cern.ch>
///

#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/DataModel/CaloClusters.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "PHOSBase/Geometry.h"

#include <TLorentzVector.h>

#include <array>
#include <climits>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod::evsel;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct PhosNonlin {
  Configurable<bool> skimmedProcessing{"skimmedProcessing", true, "Skimmed dataset processing"};
  Configurable<std::string> trigName{"trigName", "fPHOSPhoton", "name of offline trigger"};
  Configurable<std::string> zorroCCDBpath{"zorroCCDBpath", "/Users/m/mpuccio/EventFiltering/OTS/", "path to the zorro ccdb objects"};
  Configurable<int> evSelTrig{"evSelTrig", kTVXinPHOS, "Select events with this trigger"};
  Configurable<int> paramType{"paramType", 0, "Functional form 0: a-la data, 1: a-la MC"};
  Configurable<float> minCluE{"minCluE", 0.1, "Minimum cluster energy for analysis"};
  Configurable<float> minCluTime{"minCluTime", -25.e-9, "Min. cluster time"};
  Configurable<float> maxCluTime{"maxCluTime", 25.e-9, "Max. cluster time"};
  Configurable<int> minCluNcell{"minCluNcell", 1, "min cells in cluster"};
  Configurable<float> minM02{"minM02", 0.2, "Min disp M02 cut"};
  Configurable<int> nMixedEvents{"nMixedEvents", 2, "number of events to mix"};
  Configurable<float> mA{"mA", 9.34913e-01, "A"};
  Configurable<float> mdAi{"mdAi", 0., "A var. vs i"};
  Configurable<float> mdAj{"mdAj", 0., "A var. vs j"};
  Configurable<float> mB{"mB", 8.3e-02, "B"};
  Configurable<float> mdBi{"mdBi", 0., "B var. vs i"};
  Configurable<float> mdBj{"mdBj", 0., "B var. vs j"};
  Configurable<float> mC{"mC", -10.5e-03, "C"};
  Configurable<float> mdCi{"mdCi", 0.0003, "C var. vs i"};
  Configurable<float> mdCj{"mdCj", 0., "C var. vs j"};
  Configurable<float> mD{"mD", 4.5e-03, "D"};
  Configurable<float> mdDk{"mdDk", 3.e-04, "D var. vs k"};
  Configurable<float> mdDl{"mdDl", 0., "D var. vs l"};
  Configurable<float> mE{"mE", 1.4e+00, "E"};
  Configurable<float> mdEk{"mdEk", 0., "E var. vs k"};
  Configurable<float> mdEl{"mdEl", 0.1, "E var. vs l"};
  Configurable<float> mF{"mF", 7.09014e-01, "F"};
  Configurable<float> mdFk{"mdFk", 0., "F var. vs k"};
  Configurable<float> mdFl{"mdFl", 0., "F var. vs l"};
  Configurable<float> mG{"mG", 10.05600e-03, "G"};
  Configurable<float> mdGk{"mdGk", 0., "G var. vs k"};
  Configurable<float> mdGl{"mdGl", 0., "G var. vs l"};
  Configurable<float> mS{"mS", 0.8e-4, "S"};
  Configurable<float> mdSi{"mdSi", 0., "S var. vs i"};
  Configurable<float> mdSj{"mdSj", 2.e-5, "S var. vs j"};

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using BCsWithBcSels = soa::Join<aod::BCs, aod::BcSels>;
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry mHistManager1{"phosNonlinHistograms"};

  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  // class to keep photon candidate parameters
  class Photon : public TLorentzVector
  {
   public:
    Photon() = default;
    Photon(const Photon& p) = default;
    Photon(double px, double py, double pz, double e, int m) : TLorentzVector(px, py, pz, e), mod(m) {}

   public:
    int mod;
  };

  int mRunNumber = -1;   // Current run number
  int mixedEventBin = 0; // Which list of Mixed use for mixing
  std::vector<Photon> mCurEvent;
  static constexpr int kMaxMixBins = 10; // maximal number of kinds of events for mixing
  std::array<std::deque<std::vector<Photon>>, kMaxMixBins> mMixedEvents;

  // fast access to histos
  static constexpr int kNp = 10;
  std::array<TH2*, kNp * kNp> hReIJ, hReKL, hReMIJ, hReMKL;
  TH2* hMi;
  std::vector<double> pt = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2,
                            1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0,
                            6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 10., 11., 12., 13., 14., 15., 16., 18., 20., 22., 24., 26., 28.,
                            30., 34., 38., 42., 46., 50.};

  /// \brief Create output histograms
  void init(InitContext const&)
  {
    LOG(info) << "Initializing PHOS nonlin analysis task ...";

    zorroSummary.setObject(zorro.getZorroSummary());
    zorro.setBaseCCDBPath(zorroCCDBpath.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    const AxisSpec ptAxis{pt, "p_{T} (GeV/c)"},
      mggAxis{150, 0., 0.3, "m_{#gamma#gamma} (GeV/c^{2})"};

    for (int i = 0; i < kNp; i++) {
      for (int j = 0; j < kNp; j++) {
        hReIJ[i * kNp + j] = std::get<std::shared_ptr<TH2>>(mHistManager1.add(Form("hRe_a%d_b%d", i, j), "inv mass",
                                                                              HistType::kTH2F, {mggAxis, ptAxis}))
                               .get();
        hReKL[i * kNp + j] = std::get<std::shared_ptr<TH2>>(mHistManager1.add(Form("hRe_c%d_d%d", i, j), "inv mass",
                                                                              HistType::kTH2F, {mggAxis, ptAxis}))
                               .get();
        // hReMIJ[i*kNp+j] = std::get<std::shared_ptr<TH2>>(mHistManager2.add(Form("hReM_a%d_b%d",i,j), "inv mass",
        //           HistType::kTH2F, {mggAxis, ptAxis})).get();
        // hReMKL[i*kNp+j] = std::get<std::shared_ptr<TH2>>(mHistManager2.add(Form("hReM_c%d_d%d",i,j), "inv mass",
        //           HistType::kTH2F, {mggAxis, ptAxis})).get();
      }
    }
    hMi = std::get<std::shared_ptr<TH2>>(mHistManager1.add("hMi", "inv mass",
                                                           HistType::kTH2F, {mggAxis, ptAxis}))
            .get();
  }

  /// \brief Process PHOS data
  void process(SelCollisions::iterator const& col,
               aod::CaloClusters const& clusters,
               aod::BCsWithTimestamps const&)
  {
    // Fill number of events of different kind
    if (skimmedProcessing) {
      auto bc = col.template bc_as<aod::BCsWithTimestamps>();
      if (mRunNumber != bc.runNumber()) {
        zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), trigName);
        zorro.populateHistRegistry(mHistManager1, bc.runNumber());
        mRunNumber = bc.runNumber();
      }

      if (!zorro.isSelected(bc.globalBC())) {
        return; ///
      }
    } else {
      if (!col.selection_bit(evSelTrig)) {
        return;
      }
    }

    mixedEventBin = findMixedEventBin(col.posZ());

    mCurEvent.clear();
    int i, j, k, l;
    for (auto const& clu : clusters) {
      if (clu.e() < minCluE ||
          clu.ncell() < minCluNcell ||
          clu.time() > maxCluTime || clu.time() < minCluTime ||
          clu.m02() < minM02) {
        continue;
      }
      Photon ph1(clu.px(), clu.py(), clu.pz(), clu.e(), clu.mod());
      // Mix with other photons added to stack
      for (auto const& ph2 : mCurEvent) {
        double m = (ph1 + ph2).M();
        double pt1 = ph1.Pt();
        double pt2 = ph2.Pt();

        k = kNp / 2;
        l = kNp / 2;
        for (i = 0; i < kNp; i++) {
          for (j = 0; j < kNp; j++) {
            if (ph1.E() * nonLin(ph1.E(), i, j, k, l) > minCluE && ph2.E() * nonLin(ph2.E(), i, j, k, l) > minCluE) {
              double m12 = m * std::sqrt(nonLin(ph1.E(), i, j, k, l) * nonLin(ph2.E(), i, j, k, l));
              hReIJ[i * kNp + j]->Fill(m12, pt1);
              hReIJ[i * kNp + j]->Fill(m12, pt2);
              // if(ph1.mod==ph2.mod){
              //   hReMIJ[i*kNp + j]->Fill(m12,pt1);
              //   hReMIJ[i*kNp + j]->Fill(m12,pt2);
              // }
            }
          }
        }
        i = kNp / 2;
        j = kNp / 2;
        for (k = 0; k < kNp; k++) {
          for (l = 0; l < kNp; l++) {
            if (ph1.E() * nonLin(ph1.E(), i, j, k, l) > minCluE && ph2.E() * nonLin(ph2.E(), i, j, k, l) > minCluE) {
              double m12 = m * std::sqrt(nonLin(ph1.E(), i, j, k, l) * nonLin(ph2.E(), i, j, k, l));
              hReKL[k * kNp + l]->Fill(m12, pt1);
              hReKL[k * kNp + l]->Fill(m12, pt2);
              // if(ph1.mod==ph2.mod){
              //   hReMKL[k*kNp + l]->Fill(m12,pt1);
              //   hReMKL[k*kNp + l]->Fill(m12,pt2);
              // }
            }
          }
        }
      }
      // Add photon to stack
      mCurEvent.emplace_back(ph1);
    }

    // Mixed
    for (const auto& ph1 : mCurEvent) {
      for (const auto& mixEvent : mMixedEvents[mixedEventBin]) {
        for (const auto& ph2 : mixEvent) {
          double m = (ph1 + ph2).M();
          double pt1 = ph1.Pt();
          double pt2 = ph2.Pt();

          i = kNp / 2;
          j = kNp / 2;
          k = kNp / 2;
          l = kNp / 2;
          if (ph1.E() * nonLin(ph1.E(), i, j, k, l) > minCluE && ph2.E() * nonLin(ph2.E(), i, j, k, l) > minCluE) {
            double m12 = m * std::sqrt(nonLin(ph1.E(), i, j, k, l) * nonLin(ph2.E(), i, j, k, l));
            hMi->Fill(m12, pt1);
            hMi->Fill(m12, pt2);
          }
        }
      }
    }

    // Fill events to store and remove oldest to keep buffer size
    if (mCurEvent.size() > 0) {
      mMixedEvents[mixedEventBin].emplace_back(mCurEvent);
      if (mMixedEvents[mixedEventBin].size() > static_cast<size_t>(nMixedEvents)) {
        mMixedEvents[mixedEventBin].pop_front();
      }
    }
  }

  //_____________________________________________________________________________
  int findMixedEventBin(double vtxZ)
  {
    // calculate index for event mixing
    const double zwidth = 1.; // Width of zvtx bin
    int res = static_cast<int>((vtxZ + 10.) / zwidth);

    if (res < 0)
      return 0;
    if (res >= kMaxMixBins)
      return kMaxMixBins - 1;
    return res;
  }
  //_____________________________________________________________________________
  double nonLin(double en, int i, int j, int k, int l)
  {
    if (en <= 0.)
      return 0.;
    if (paramType == 0) {
      const double a = mA + mdAi * (i - kNp / 2) + mdAj * (j - kNp / 2);
      const double b = mB + mdBi * (i - kNp / 2) + mdBj * (j - kNp / 2);
      const double c = mC + mdCi * (i - kNp / 2) + mdCj * (j - kNp / 2);
      const double d = mD + mdDk * (k - kNp / 2) + mdDl * (l - kNp / 2);
      const double e = mE + mdEk * (k - kNp / 2) + mdEl * (l - kNp / 2);
      const double f = mF + mdFk * (k - kNp / 2) + mdFl * (l - kNp / 2);
      const double g = mG + mdGk * (k - kNp / 2) + mdGl * (l - kNp / 2);
      const double s = mS + mdSi * (i - kNp / 2) + mdSj * (j - kNp / 2);
      return a + b / en + c / (en * en) + d / ((en - e) * (en - e) + f * f) + g * en + s / (en * en * en * en);
    }
    if (paramType == 1) {
      const double a = mA + mdAi * (i - kNp / 2) + mdAj * (j - kNp / 2);
      const double b = mB + mdBi * (i - kNp / 2) + mdBj * (j - kNp / 2);
      const double c = mC + mdCi * (i - kNp / 2) + mdCj * (j - kNp / 2);
      const double d = mD + mdDk * (k - kNp / 2) + mdDl * (l - kNp / 2);
      const double e = mE + mdEk * (k - kNp / 2) + mdEl * (l - kNp / 2);
      return a + b / std::sqrt(en) + c / en + d / (en * std::sqrt(en)) + e / (en * en);
    }
    return 0.;
  }
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return o2::framework::WorkflowSpec{
    o2::framework::adaptAnalysisTask<PhosNonlin>(cfgc)};
}
