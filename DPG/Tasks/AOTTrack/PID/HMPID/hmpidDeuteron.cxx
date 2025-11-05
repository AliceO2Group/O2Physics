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

#include "tableHMPID.h"

#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

#include <TMath.h>
#include <TString.h>

using namespace o2;
using namespace o2::framework;

struct HmpidDeuteron {
  HistogramRegistry registryDA{"registryDA", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> nsigmaTPCMin{"nsigmaTPCMin", -3.0, "nsigmaTPCMin"};
  Configurable<float> nsigmaTPCMax{"nsigmaTPCMax", +3.0, "nsigmaTPCMax"};
  Configurable<float> nsigmaTOFMin{"nsigmaTOFMin", -3.0, "nsigmaTOFMin"};
  Configurable<float> nsigmaTOFMax{"nsigmaTOFMax", +3.5, "nsigmaTOFMax"};
  Configurable<float> minReqClusterITS{"minReqClusterITS", 4.0, "min number of clusters required in ITS"};
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 50.0f, "minTPCnClsFound"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.0f, "min number of crossed rows TPC"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  Configurable<float> maxDCAxy{"maxDCAxy", 0.5f, "maxDCAxy"};
  Configurable<float> maxDCAz{"maxDCAz", 0.5f, "maxDCAz"};

  void init(InitContext const&)
  {
    // Deuteron Pos
    registryDA.add("incomingDe_Pos_8cm", "incomingDe_Pos_8cm", HistType::kTH1F, {{490, 0.1, 5.0, "#it{p} (GeV/#it{c})"}});
    registryDA.add("incomingDe_Pos_4cm", "incomingDe_Pos_4cm", HistType::kTH1F, {{490, 0.1, 5.0, "#it{p} (GeV/#it{c})"}});

    registryDA.add("De_Pos_deltaR_8cm", "De_Pos_deltaR_8cm", HistType::kTH1F, {{300, 0.0, 30.0, "#Delta R (cm)"}});
    registryDA.add("De_Pos_deltaR_4cm", "De_Pos_deltaR_4cm", HistType::kTH1F, {{300, 0.0, 30.0, "#Delta R (cm)"}});

    registryDA.add("survivingDe_Pos_8cm", "survivingDe_Pos_8cm", HistType::kTH2F, {{490, 0.1, 5.0, "#it{p} (GeV/#it{c})"}, {300, 0.0, 30.0, "#Delta R (cm)"}});
    registryDA.add("survivingDe_Pos_4cm", "survivingDe_Pos_4cm", HistType::kTH2F, {{490, 0.1, 5.0, "#it{p} (GeV/#it{c})"}, {300, 0.0, 30.0, "#Delta R (cm)"}});
    registryDA.add("De_Pos_Q_8cm", "De_Pos_Q_8cm", HistType::kTH2F, {{490, 0.1, 5.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Q (ADC)"}});
    registryDA.add("De_Pos_Q_4cm", "De_Pos_Q_4cm", HistType::kTH2F, {{490, 0.1, 5.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Q (ADC)"}});
    registryDA.add("De_Pos_ClsSize_8cm", "De_Pos_ClsSize_8cm", HistType::kTH2F, {{490, 0.1, 5.0, "#it{p} (GeV/#it{c})"}, {20, 0.0, 20., "Cls size"}});
    registryDA.add("De_Pos_ClsSize_4cm", "De_Pos_ClsSize_4cm", HistType::kTH2F, {{490, 0.1, 5.0, "#it{p} (GeV/#it{c})"}, {20, 0.0, 20.0, "Cls size"}});
    registryDA.add("De_Pos_momentum", "De_Pos_momentum", HistType::kTH2F, {{100, 0.0, 5.0, "#it{p}_{vtx} (GeV/#it{c})"}, {100, 0.0, 5.0, "#it{p}_{hmpid} (GeV/#it{c})"}});

    registryDA.add("nSigmaTPC_vs_momHMPID_noCut_DePos", "nSigmaTPC_vs_momHMPID_noCut_DePos", HistType::kTH2F, {{490, 0.1, 5.0, "#it{p} (GeV/#it{c})"}, {20, -5., 5.0, "n#sigma_TPC"}});
    registryDA.add("nSigmaTOF_vs_momHMPID_noCut_DePos", "nSigmaTOF_vs_momHMPID_noCut_DePos", HistType::kTH2F, {{490, 0.1, 5.0, "#it{p} (GeV/#it{c})"}, {20, -5., 5.0, "n#sigma_TOF"}});
    registryDA.add("nSigmaTPC_vs_momHMPID_Cut_DePos", "nSigmaTPC_vs_momHMPID_Cut_DePos", HistType::kTH2F, {{490, 0.1, 5.0, "#it{p} (GeV/#it{c})"}, {20, -5., 5.0, "n#sigma_TPC"}});
    registryDA.add("nSigmaTOF_vs_momHMPID_Cut_DePos", "nSigmaTOF_vs_momHMPID_Cut_DePos", HistType::kTH2F, {{490, 0.1, 5.0, "#it{p} (GeV/#it{c})"}, {20, -5., 5.0, "n#sigma_TOF"}});

    // Deuteron Neg
    registryDA.add("incomingDe_Neg_8cm", "incomingDe_Neg_8cm", HistType::kTH1F, {{490, 0.1, 5.0, "#it{p} (GeV/#it{c})"}});
    registryDA.add("incomingDe_Neg_4cm", "incomingDe_Neg_4cm", HistType::kTH1F, {{490, 0.1, 5.0, "#it{p} (GeV/#it{c})"}});

    // plot aggiunti
    registryDA.add("De_Neg_deltaR_8cm", "De_Neg_deltaR_8cm", HistType::kTH1F, {{300, 0.0, 30.0, "#Delta R (cm)"}});
    registryDA.add("De_Neg_deltaR_4cm", "De_Neg_deltaR_4cm", HistType::kTH1F, {{300, 0.0, 30.0, "#Delta R (cm)"}});

    registryDA.add("survivingDe_Neg_8cm", "survivingDe_Neg_8cm", HistType::kTH2F, {{490, 0.1, 5.0, "#it{p} (GeV/#it{c})"}, {300, 0.0, 30.0, "#Delta R (cm)"}});
    registryDA.add("survivingDe_Neg_4cm", "survivingDe_Neg_4cm", HistType::kTH2F, {{490, 0.1, 5.0, "#it{p} (GeV/#it{c})"}, {300, 0.0, 30.0, "#Delta R (cm)"}});
    registryDA.add("De_Neg_Q_8cm", "De_Neg_Q_8cm", HistType::kTH2F, {{490, 0.1, 5.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Q (ADC)"}});
    registryDA.add("De_Neg_Q_4cm", "De_Neg_Q_4cm", HistType::kTH2F, {{490, 0.1, 5.0, "#it{p} (GeV/#it{c})"}, {200, 0.0, 2000.0, "Q (ADC)"}});
    registryDA.add("De_Neg_ClsSize_8cm", "De_Neg_ClsSize_8cm", HistType::kTH2F, {{490, 0.1, 5.0, "#it{p} (GeV/#it{c})"}, {20, 0.0, 20.0, "Cls size"}});
    registryDA.add("De_Neg_ClsSize_4cm", "De_Neg_ClsSize_4cm", HistType::kTH2F, {{490, 0.1, 5.0, "#it{p} (GeV/#it{c})"}, {20, 0.0, 20.0, "Cls size"}});
    registryDA.add("De_Neg_momentum", "De_Neg_momentum", HistType::kTH2F, {{100, 0.0, 5.0, "#it{p}_{vtx} (GeV/#it{c})"}, {100, 0.0, 5.0, "#it{p}_{hmpid} (GeV/#it{c})"}});

    registryDA.add("nSigmaTPC_vs_momHMPID_noCut_DeNeg", "nSigmaTPC_vs_momHMPID_noCut_DeNeg", HistType::kTH2F, {{490, 0.1, 5.0, "#it{p} (GeV/#it{c})"}, {20, -5.0, 5.0, "n#sigma_TPC"}});
    registryDA.add("nSigmaTOF_vs_momHMPID_noCut_DeNeg", "nSigmaTOF_vs_momHMPID_noCut_DeNeg", HistType::kTH2F, {{490, 0.1, 5.0, "#it{p} (GeV/#it{c})"}, {20, -5., 5.0, "n#sigma_TOF"}});
    registryDA.add("nSigmaTPC_vs_momHMPID_Cut_DeNeg", "nSigmaTPC_vs_momHMPID_Cut_DeNeg", HistType::kTH2F, {{490, 0.1, 5.0, "#it{p} (GeV/#it{c})"}, {20, -5., 5.0, "n#sigma_TPC"}});
    registryDA.add("nSigmaTOF_vs_momHMPID_Cut_DeNeg", "nSigmaTOF_vs_momHMPID_Cut_DeNeg", HistType::kTH2F, {{490, 0.1, 5.0, "#it{p} (GeV/#it{c})"}, {20, -5., 5.0, "n#sigma_TOF"}});

    registryDA.add("hmpidCkovvsMom", "hmpidCkovvsMom", kTH2F, {{500, 0, 10., "#it{p} (GeV/#it{c})"}, {800, 0., 0.8, "#theta_{Ch} (rad)"}});
  }

  void process(aod::HmpidAnalysis const& hmpidTable)
  {
    for (const auto& hmpid : hmpidTable) {

      // filters on primary tracks
      if (hmpid.itsNCluster() < minReqClusterITS)
        continue;
      if (hmpid.tpcNCluster() < minTPCnClsFound)
        continue;
      if (hmpid.tpcNClsCrossedRows() < minNCrossedRowsTPC)
        continue;
      if (hmpid.tpcChi2() > maxChi2TPC)
        continue;
      if (hmpid.itsChi2() > maxChi2ITS)
        continue;
      if (TMath::Abs(hmpid.dcaXY()) > maxDCAxy)
        continue;
      if (TMath::Abs(hmpid.dcaZ()) > maxDCAz)
        continue;

      // plots nsigma before cuts
      if (hmpid.momentumHmpid() > 0) {
        registryDA.fill(HIST("nSigmaTPC_vs_momHMPID_noCut_DePos"), std::fabs(hmpid.momentumHmpid()), hmpid.tpcNSigmaDe());
        registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_noCut_DePos"), std::fabs(hmpid.momentumHmpid()), hmpid.tofNSigmaDe());
      }

      if (hmpid.momentumHmpid() < 0) {
        registryDA.fill(HIST("nSigmaTPC_vs_momHMPID_noCut_DeNeg"), std::fabs(hmpid.momentumHmpid()), hmpid.tpcNSigmaDe());
        registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_noCut_DeNeg"), std::fabs(hmpid.momentumHmpid()), hmpid.tofNSigmaDe());
      }

      // deuteron candidate cuts
      if (hmpid.tpcNSigmaDe() < nsigmaTPCMin || hmpid.tpcNSigmaDe() > nsigmaTPCMax)
        continue;

      if (hmpid.momentumHmpid() > 0)
        registryDA.fill(HIST("nSigmaTPC_vs_momHMPID_Cut_DePos"), std::fabs(hmpid.momentumHmpid()), hmpid.tpcNSigmaDe());
      else
        registryDA.fill(HIST("nSigmaTPC_vs_momHMPID_Cut_DeNeg"), std::fabs(hmpid.momentumHmpid()), hmpid.tpcNSigmaDe());

      if (hmpid.tofNSigmaDe() < nsigmaTOFMin || hmpid.tofNSigmaDe() > nsigmaTOFMax)
        continue;

      if (hmpid.momentumHmpid() > 0)
        registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_Cut_DePos"), std::fabs(hmpid.momentumHmpid()), hmpid.tofNSigmaDe());
      else
        registryDA.fill(HIST("nSigmaTOF_vs_momHMPID_Cut_DeNeg"), std::fabs(hmpid.momentumHmpid()), hmpid.tofNSigmaDe());

      // plot changle vs p
      registryDA.fill(HIST("hmpidCkovvsMom"), hmpid.momentumHmpid(), hmpid.chAngle());

      // absorbers
      const float dx = hmpid.xTrack() - hmpid.xMip();
      const float dy = hmpid.yTrack() - hmpid.yMip();
      const float dr = TMath::Sqrt(dx * dx + dy * dy);

      int abs4cm = 2, abs8cm = 4;
      bool hmpidAbs8cm = true, hmpidAbs4cm = true;

      if (hmpid.momentumHmpid() > 0) {
        registryDA.fill(HIST("De_Pos_momentum"), hmpid.momentumTrack(), std::fabs(hmpid.momentumHmpid()));

        if (hmpidAbs8cm && hmpid.chamber() == abs8cm) {
          registryDA.fill(HIST("incomingDe_Pos_8cm"), std::fabs(hmpid.momentumHmpid()));
          registryDA.fill(HIST("survivingDe_Pos_8cm"), std::fabs(hmpid.momentumHmpid()), dr);
          registryDA.fill(HIST("De_Pos_Q_8cm"), std::fabs(hmpid.momentumHmpid()), hmpid.chargeMip());
          registryDA.fill(HIST("De_Pos_ClsSize_8cm"), std::fabs(hmpid.momentumHmpid()), hmpid.clusterSize());
          registryDA.fill(HIST("De_Pos_deltaR_8cm"), dr);
        }
        if (hmpidAbs4cm && hmpid.chamber() == abs4cm) {
          registryDA.fill(HIST("incomingDe_Pos_4cm"), std::fabs(hmpid.momentumHmpid()));
          registryDA.fill(HIST("survivingDe_Pos_4cm"), std::fabs(hmpid.momentumHmpid()), dr);
          registryDA.fill(HIST("De_Pos_Q_4cm"), std::fabs(hmpid.momentumHmpid()), hmpid.chargeMip());
          registryDA.fill(HIST("De_Pos_ClsSize_4cm"), std::fabs(hmpid.momentumHmpid()), hmpid.clusterSize());
          registryDA.fill(HIST("De_Pos_deltaR_4cm"), dr);
        }
      }

      if (hmpid.momentumHmpid() < 0) {
        registryDA.fill(HIST("De_Neg_momentum"), hmpid.momentumTrack(), std::fabs(hmpid.momentumHmpid()));

        if (hmpidAbs8cm && hmpid.chamber() == abs8cm) {
          registryDA.fill(HIST("incomingDe_Neg_8cm"), std::fabs(hmpid.momentumHmpid()));
          registryDA.fill(HIST("survivingDe_Neg_8cm"), std::fabs(hmpid.momentumHmpid()), dr);
          registryDA.fill(HIST("De_Neg_Q_8cm"), std::fabs(hmpid.momentumHmpid()), hmpid.chargeMip());
          registryDA.fill(HIST("De_Neg_ClsSize_8cm"), std::fabs(hmpid.momentumHmpid()), hmpid.clusterSize());
          registryDA.fill(HIST("De_Neg_deltaR_8cm"), dr);
        }
        if (hmpidAbs4cm && hmpid.chamber() == abs4cm) {
          registryDA.fill(HIST("incomingDe_Neg_4cm"), std::fabs(hmpid.momentumHmpid()));
          registryDA.fill(HIST("survivingDe_Neg_4cm"), std::fabs(hmpid.momentumHmpid()), dr);
          registryDA.fill(HIST("De_Neg_Q_4cm"), std::fabs(hmpid.momentumHmpid()), hmpid.chargeMip());
          registryDA.fill(HIST("De_Neg_ClsSize_4cm"), std::fabs(hmpid.momentumHmpid()), hmpid.clusterSize());
          registryDA.fill(HIST("De_Neg_deltaR_4cm"), dr);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HmpidDeuteron>(cfgc),
  };
}
