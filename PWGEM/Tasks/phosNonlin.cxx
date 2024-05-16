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

#include <climits>
#include <cstdlib>
#include <map>
#include <memory>
#include <vector>
#include <array>
#include <TLorentzVector.h>
#include "Common/DataModel/CaloClusters.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include "Framework/ConfigParamSpec.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"

#include "PHOSBase/Geometry.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPLHCIFData.h"

/// \struct PHOS pi0 analysis
/// \brief Monitoring task for PHOS related quantities
/// \author Dmitri Peresunko, NRC "Kurchatov institute"
/// \since Nov, 2022
///

using namespace o2;
using namespace o2::aod::evsel;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct phosNonlin {
  Configurable<int> mEvSelTrig{"mEvSelTrig", kTVXinPHOS, "Select events with this trigger"};
  Configurable<int> mParamType{"mParamType", 0, "Functional form 0: a-la data, 1: a-la MC"};
  Configurable<float> mMinCluE{"mMinCluE", 0.1, "Minimum cluster energy for analysis"};
  Configurable<float> mMinCluTime{"minCluTime", -25.e-9, "Min. cluster time"};
  Configurable<float> mMaxCluTime{"maxCluTime", 25.e-9, "Max. cluster time"};
  Configurable<int> mMinCluNcell{"minCluNcell", 1, "min cells in cluster"};
  Configurable<float> mMinM02{"minM02", 0.2, "Min disp M02 cut"};
  Configurable<int> mNMixedEvents{"nMixedEvents", 2, "number of events to mix"};
  Configurable<bool> mSelectOneCollPerBC{"selectOneColPerBC", true, "skip multiple coll. per bc"};
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

  HistogramRegistry mHistManager1{"phosNonlinHistograms"};

  // class to keep photon candidate parameters
  class photon : public TLorentzVector
  {
   public:
    photon() = default;
    photon(const photon& p) = default;
    photon(double px, double py, double pz, double e, int m) : TLorentzVector(px, py, pz, e), mod(m) {}

   public:
    int mod;
  };

  int mixedEventBin = 0; // Which list of Mixed use for mixing
  std::vector<photon> mCurEvent;
  static constexpr int nMaxMixBins = 10; // maximal number of kinds of events for mixing
  std::array<std::deque<std::vector<photon>>, nMaxMixBins> mMixedEvents;

  // fast access to histos
  static constexpr int mNp = 10;
  std::array<TH2*, mNp * mNp> hReIJ, hReKL, hReMIJ, hReMKL;
  TH2* hMi;
  std::vector<double> pt = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2,
                            1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0,
                            6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 10., 11., 12., 13., 14., 15., 16., 18., 20., 22., 24., 26., 28.,
                            30., 34., 38., 42., 46., 50.};

  /// \brief Create output histograms
  void init(InitContext const&)
  {
    LOG(info) << "Initializing PHOS nonlin analysis task ...";

    const AxisSpec
      ptAxis{pt, "p_{T} (GeV/c)"},
      mggAxis{150, 0., 0.3, "m_{#gamma#gamma} (GeV/c^{2})"};

    for (int i = 0; i < mNp; i++) {
      for (int j = 0; j < mNp; j++) {
        hReIJ[i * mNp + j] = std::get<std::shared_ptr<TH2>>(mHistManager1.add(Form("hRe_a%d_b%d", i, j), "inv mass",
                                                                              HistType::kTH2F, {mggAxis, ptAxis}))
                               .get();
        hReKL[i * mNp + j] = std::get<std::shared_ptr<TH2>>(mHistManager1.add(Form("hRe_c%d_d%d", i, j), "inv mass",
                                                                              HistType::kTH2F, {mggAxis, ptAxis}))
                               .get();
        // hReMIJ[i*mNp+j] = std::get<std::shared_ptr<TH2>>(mHistManager2.add(Form("hReM_a%d_b%d",i,j), "inv mass",
        //           HistType::kTH2F, {mggAxis, ptAxis})).get();
        // hReMKL[i*mNp+j] = std::get<std::shared_ptr<TH2>>(mHistManager2.add(Form("hReM_c%d_d%d",i,j), "inv mass",
        //           HistType::kTH2F, {mggAxis, ptAxis})).get();
      }
    }
    hMi = std::get<std::shared_ptr<TH2>>(mHistManager1.add("hMi", "inv mass",
                                                           HistType::kTH2F, {mggAxis, ptAxis}))
            .get();
  }

  /// \brief Process PHOS data
  void process(SelCollisions::iterator const& col,
               aod::CaloClusters const& clusters)
  {
    // Fill number of events of different kind
    if (!col.alias_bit(mEvSelTrig)) {
      return;
    }

    mixedEventBin = findMixedEventBin(col.posZ());

    mCurEvent.clear();
    int i, j, k, l;
    for (const auto& clu : clusters) {
      if (clu.e() < mMinCluE ||
          clu.ncell() < mMinCluNcell ||
          clu.time() > mMaxCluTime || clu.time() < mMinCluTime ||
          clu.m02() < mMinM02) {
        continue;
      }
      photon ph1(clu.px(), clu.py(), clu.pz(), clu.e(), clu.mod());
      // Mix with other photons added to stack
      for (auto ph2 : mCurEvent) {
        double m = (ph1 + ph2).M();
        double pt1 = ph1.Pt();
        double pt2 = ph2.Pt();

        k = mNp / 2;
        l = mNp / 2;
        for (i = 0; i < mNp; i++) {
          for (j = 0; j < mNp; j++) {
            if (ph1.E() * NonLin(ph1.E(), i, j, k, l) > mMinCluE && ph2.E() * NonLin(ph2.E(), i, j, k, l) > mMinCluE) {
              Double_t m12 = m * TMath::Sqrt(NonLin(ph1.E(), i, j, k, l) * NonLin(ph2.E(), i, j, k, l));
              hReIJ[i * mNp + j]->Fill(m12, pt1);
              hReIJ[i * mNp + j]->Fill(m12, pt2);
              // if(ph1.mod==ph2.mod){
              //   hReMIJ[i*mNp + j]->Fill(m12,pt1);
              //   hReMIJ[i*mNp + j]->Fill(m12,pt2);
              // }
            }
          }
        }
        i = mNp / 2;
        j = mNp / 2;
        for (k = 0; k < mNp; k++) {
          for (l = 0; l < mNp; l++) {
            if (ph1.E() * NonLin(ph1.E(), i, j, k, l) > mMinCluE && ph2.E() * NonLin(ph2.E(), i, j, k, l) > mMinCluE) {
              Double_t m12 = m * TMath::Sqrt(NonLin(ph1.E(), i, j, k, l) * NonLin(ph2.E(), i, j, k, l));
              hReKL[k * mNp + l]->Fill(m12, pt1);
              hReKL[k * mNp + l]->Fill(m12, pt2);
              // if(ph1.mod==ph2.mod){
              //   hReMKL[k*mNp + l]->Fill(m12,pt1);
              //   hReMKL[k*mNp + l]->Fill(m12,pt2);
              // }
            }
          }
        }
      }
      // Add photon to stack
      mCurEvent.emplace_back(ph1);
    }

    // Mixed
    for (auto ph1 : mCurEvent) {
      for (auto mixEvent : mMixedEvents[mixedEventBin]) {
        for (auto ph2 : mixEvent) {
          double m = (ph1 + ph2).M();
          double pt1 = ph1.Pt();
          double pt2 = ph2.Pt();

          i = mNp / 2;
          j = mNp / 2;
          k = mNp / 2;
          l = mNp / 2;
          if (ph1.E() * NonLin(ph1.E(), i, j, k, l) > mMinCluE && ph2.E() * NonLin(ph2.E(), i, j, k, l) > mMinCluE) {
            Double_t m12 = m * TMath::Sqrt(NonLin(ph1.E(), i, j, k, l) * NonLin(ph2.E(), i, j, k, l));
            hMi->Fill(m12, pt1);
            hMi->Fill(m12, pt2);
          }
        }
      }
    }

    // Fill events to store and remove oldest to keep buffer size
    if (mCurEvent.size() > 0) {
      mMixedEvents[mixedEventBin].emplace_back(mCurEvent);
      if (mMixedEvents[mixedEventBin].size() > static_cast<size_t>(mNMixedEvents)) {
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
    if (res >= nMaxMixBins)
      return nMaxMixBins - 1;
    return res;
  }
  //_____________________________________________________________________________
  double NonLin(double en, int i, int j, int k, int l)
  {
    if (en <= 0.)
      return 0.;
    if (mParamType == 0) {
      const Double_t a = mA + mdAi * (i - mNp / 2) + mdAj * (j - mNp / 2);
      const Double_t b = mB + mdBi * (i - mNp / 2) + mdBj * (j - mNp / 2);
      const Double_t c = mC + mdCi * (i - mNp / 2) + mdCj * (j - mNp / 2);
      const Double_t d = mD + mdDk * (k - mNp / 2) + mdDl * (l - mNp / 2);
      const Double_t e = mE + mdEk * (k - mNp / 2) + mdEl * (l - mNp / 2);
      const Double_t f = mF + mdFk * (k - mNp / 2) + mdFl * (l - mNp / 2);
      const Double_t g = mG + mdGk * (k - mNp / 2) + mdGl * (l - mNp / 2);
      const Double_t s = mS + mdSi * (i - mNp / 2) + mdSj * (j - mNp / 2);
      return a + b / en + c / (en * en) + d / ((en - e) * (en - e) + f * f) + g * en + s / (en * en * en * en);
    }
    if (mParamType == 1) {
      const Double_t a = mA + mdAi * (i - mNp / 2) + mdAj * (j - mNp / 2);
      const Double_t b = mB + mdBi * (i - mNp / 2) + mdBj * (j - mNp / 2);
      const Double_t c = mC + mdCi * (i - mNp / 2) + mdCj * (j - mNp / 2);
      const Double_t d = mD + mdDk * (k - mNp / 2) + mdDl * (l - mNp / 2);
      const Double_t e = mE + mdEk * (k - mNp / 2) + mdEl * (l - mNp / 2);
      return a + b / TMath::Sqrt(en) + c / en + d / (en * TMath::Sqrt(en)) + e / (en * en);
    }
    return 0.;
  }
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return o2::framework::WorkflowSpec{
    o2::framework::adaptAnalysisTask<phosNonlin>(cfgc)};
}
