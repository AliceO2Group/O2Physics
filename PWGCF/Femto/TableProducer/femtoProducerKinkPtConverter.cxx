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

/// \file femtoProducerKinkPtConverter.cxx
/// \brief Task that converts FSigmas_001 to FSigmas_002 with recalculated pT
/// \author Henrik Fribert, TU München, henrik.fribert@tum.de

#include "PWGCF/Femto/Core/femtoUtils.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <Math/Vector4D.h>

using namespace o2::analysis::femto;

struct FemtoProducerKinkPtConverter {

  o2::framework::Produces<o2::aod::FSigmas> producedSigmas;

  o2::framework::Configurable<bool> confUseRecalculatedPt{"confUseRecalculatedPt", true, "Use recalculated pT from kinematic constraints"};

  o2::framework::HistogramRegistry mHistogramRegistry{"FemtoSigmaConverter", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    mHistogramRegistry.add("hPtOriginal", "Original pT;p_{T} (GeV/c);Counts", o2::framework::kTH1F, {{100, 0, 10}});
    mHistogramRegistry.add("hPtRecalculated", "Recalculated pT;p_{T} (GeV/c);Counts", o2::framework::kTH1F, {{100, 0, 10}});
    mHistogramRegistry.add("hPtRatio", "pT Ratio (recalc/orig);p_{T,recalc} / p_{T,orig};Counts", o2::framework::kTH1F, {{200, 0, 2}});
    mHistogramRegistry.add("hRecalcSuccess", "Recalculation Success;Success (0=fail, 1=success);Counts", o2::framework::kTH1I, {{2, -0.5, 1.5}});
  }

  void process(o2::aod::FSigmas_001 const& sigmasV1,
               o2::aod::FTracks const& tracks)
  {
    for (const auto& sigma : sigmasV1) {

      float signedPtToUse = sigma.signedPt();

      if (confUseRecalculatedPt) {
        auto chaDaughter = tracks.rawIteratorAt(sigma.chaDauId() - tracks.offset());

        float pxDaug = chaDaughter.pt() * std::cos(chaDaughter.phi());
        float pyDaug = chaDaughter.pt() * std::sin(chaDaughter.phi());
        float pzDaug = chaDaughter.pt() * std::sinh(chaDaughter.eta());

        float pxMoth = sigma.pt() * std::cos(sigma.phi());
        float pyMoth = sigma.pt() * std::sin(sigma.phi());
        float pzMoth = sigma.pt() * std::sinh(sigma.eta());

        float ptRecalc = utils::calcPtnew(pxMoth, pyMoth, pzMoth, pxDaug, pyDaug, pzDaug);

        if (ptRecalc > 0) {
          signedPtToUse = ptRecalc * utils::signum(sigma.signedPt());

          mHistogramRegistry.fill(HIST("hPtOriginal"), sigma.pt());
          mHistogramRegistry.fill(HIST("hPtRecalculated"), ptRecalc);
          mHistogramRegistry.fill(HIST("hPtRatio"), ptRecalc / sigma.pt());
          mHistogramRegistry.fill(HIST("hRecalcSuccess"), 1);
        } else {
          mHistogramRegistry.fill(HIST("hRecalcSuccess"), 0);
        }
      }

      producedSigmas(sigma.fColId(),
                     signedPtToUse,
                     sigma.eta(),
                     sigma.phi(),
                     sigma.mass(),
                     sigma.chaDauId());
    }
  }
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  o2::framework::WorkflowSpec workflow{adaptAnalysisTask<FemtoProducerKinkPtConverter>(cfgc)};
  return workflow;
}
