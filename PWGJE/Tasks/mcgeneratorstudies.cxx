// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// Task that produces the generated pT spectrum of a given particle for MC studies
//
/// \author Nicolas Strangmann <nicolas.strangmann@cern.ch>, Goethe University Frankfurt / Oak Ridge National Laoratory

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"

#include "PWGJE/DataModel/EMCALMatchedCollisions.h"

#include "DetectorsBase/GeometryManager.h"
#include "EMCALBase/Geometry.h"

#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyMCCollisions = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::EMCALMatchedCollisions>;

struct MCGeneratorStudies {
  HistogramRegistry mHistManager{"MCGeneratorStudyHistograms"};

  Configurable<float> mVertexCut{"vertexCut", 10.f, "apply z-vertex cut with value in cm"};
  Configurable<float> mRapidityCut{"rapidityCut", 0.9f, "Maximum absolute rapidity of counted generated particles"};
  Configurable<int> mSelectedParticleCode{"particlePDGCode", 111, "PDG code of the particle to be investigated"};
  Configurable<bool> mRequireGammaGammaDecay{"requireGammaGammaDecay", false, "Only count generated particles that decayed into two photons"};
  Configurable<bool> mRequireEMCCellContent{"requireEMCCellContent", false, "Ask forEMCal cell content instead of the kTVXinEMC trigger"};

  void init(InitContext const&)
  {
    AxisSpec pTAxis{250, 0., 25., "#it{p}_{T} (GeV/#it{c})"};

    auto hCollisionCounter = mHistManager.add<TH1>("hCollisionCounter", "Number of collisions after event cuts", HistType::kTH1F, {{7, 0.5, 7.5}});
    hCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
    hCollisionCounter->GetXaxis()->SetBinLabel(2, "TVX");
    hCollisionCounter->GetXaxis()->SetBinLabel(3, "T zSmall");
    hCollisionCounter->GetXaxis()->SetBinLabel(4, "Tz zGood");
    hCollisionCounter->GetXaxis()->SetBinLabel(5, "Tzz EMCal");
    hCollisionCounter->GetXaxis()->SetBinLabel(6, "TzzE Sel8");
    hCollisionCounter->GetXaxis()->SetBinLabel(7, "TzzES Unique");
    TString mesonLatexString = (TString)mSelectedParticleCode;
    switch (mSelectedParticleCode) {
      case 111:
        mesonLatexString = "#pi^{0}";
        break;
      case 221:
        mesonLatexString = "#eta";
        break;
    }
    mHistManager.add("hpT_all", Form("Generated %s in all collisions", mesonLatexString.Data()), HistType::kTH1F, {pTAxis});
    mHistManager.add("hpT_TVX", Form("Generated %s in TVX triggered collisions", mesonLatexString.Data()), HistType::kTH1F, {pTAxis});
    mHistManager.add("hpT_T_zsmall", Form("Generated %s in TVX collisions with z < 10cm", mesonLatexString.Data()), HistType::kTH1F, {pTAxis});
    mHistManager.add("hpT_T_z_zGood", Form("Generated %s in TVX collisions with good z < 10cm", mesonLatexString.Data()), HistType::kTH1F, {pTAxis});
    mHistManager.add("hpTAccepted_T_z_z", Form("Accepted (EMCal) %s in TVX collisions with good z < 10cm", mesonLatexString.Data()), HistType::kTH1F, {pTAxis});
    mHistManager.add("hpT_T_z_z_EMCal", Form("Generated %s in TVXinEMC collisions with good z < 10cm", mesonLatexString.Data()), HistType::kTH1F, {pTAxis});
    mHistManager.add("hpT_T_z_z_E_Sel8", Form("Generated %s in TVXinEMC collisions with good z < 10cm and Sel8", mesonLatexString.Data()), HistType::kTH1F, {pTAxis});
    mHistManager.add("hpT_T_z_z_E_S_Unique", Form("Generated %s in unique TVXinEMC collisions with good z < 10cm and Sel8", mesonLatexString.Data()), HistType::kTH1F, {pTAxis});
    mHistManager.add("hpTAccepted_T_z_z_E_S_U", Form("Accepted %s in unique TVXinEMC collisions with good z < 10cm and Sel8", mesonLatexString.Data()), HistType::kTH1F, {pTAxis});

    auto hEMCollisionCounter = mHistManager.add<TH1>("hEMCollisionCounter", "collision counter;;Number of events", kTH1F, {{13, 0.5, 13.5}}, false);
    hEMCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
    hEMCollisionCounter->GetXaxis()->SetBinLabel(2, "No TF border");
    hEMCollisionCounter->GetXaxis()->SetBinLabel(3, "No ITS ROF border");
    hEMCollisionCounter->GetXaxis()->SetBinLabel(4, "No Same Bunch Pileup");
    hEMCollisionCounter->GetXaxis()->SetBinLabel(5, "Is Vertex ITSTPC");
    hEMCollisionCounter->GetXaxis()->SetBinLabel(6, "Is Good Zvtx FT0vsPV");
    hEMCollisionCounter->GetXaxis()->SetBinLabel(7, "FT0AND");
    hEMCollisionCounter->GetXaxis()->SetBinLabel(8, "sel8");
    hEMCollisionCounter->GetXaxis()->SetBinLabel(9, "|Z_{vtx}| < 10 cm");
    hEMCollisionCounter->GetXaxis()->SetBinLabel(10, "EMC MB Readout");
    hEMCollisionCounter->GetXaxis()->SetBinLabel(11, "EMC L0 Triggered");
    hEMCollisionCounter->GetXaxis()->SetBinLabel(12, "EMC Cell Content");
    hEMCollisionCounter->GetXaxis()->SetBinLabel(13, "accepted");

    o2::emcal::Geometry::GetInstanceFromRunNumber(300000);
  }

  PresliceUnsorted<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;

  void process(MyMCCollisions::iterator const& collision, aod::McCollisions const&, aod::McParticles const& mcParticles)
  {
    fillEventHistogram(&mHistManager, collision);

    auto mcCollision = collision.mcCollision();
    auto mcParticles_inColl = mcParticles.sliceBy(perMcCollision, mcCollision.globalIndex());

    for (auto& mcParticle : mcParticles_inColl) {
      if (mcParticle.pdgCode() != mSelectedParticleCode || fabs(mcParticle.y()) > mRapidityCut)
        continue;
      if (!mcParticle.isPhysicalPrimary() && !mcParticle.producedByGenerator())
        continue;
      if (mRequireGammaGammaDecay && !isGammaGammaDecay(mcParticle, mcParticles))
        continue;

      mHistManager.fill(HIST("hpT_all"), mcParticle.pt());
      if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
        mHistManager.fill(HIST("hpT_TVX"), mcParticle.pt());
        if (abs(collision.posZ()) < mVertexCut) {
          mHistManager.fill(HIST("hpT_T_zsmall"), mcParticle.pt());
          if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
            mHistManager.fill(HIST("hpT_T_z_zGood"), mcParticle.pt());
            if (isAccepted(mcParticle, mcParticles))
              mHistManager.fill(HIST("hpTAccepted_T_z_z"), mcParticle.pt());
            if (mRequireEMCCellContent ? collision.isemcreadout() : collision.alias_bit(kTVXinEMC)) {
              mHistManager.fill(HIST("hpT_T_z_z_EMCal"), mcParticle.pt());
              if (collision.sel8()) {
                mHistManager.fill(HIST("hpT_T_z_z_E_Sel8"), mcParticle.pt());
                if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
                  mHistManager.fill(HIST("hpT_T_z_z_E_S_Unique"), mcParticle.pt());
                  if (isAccepted(mcParticle, mcParticles))
                    mHistManager.fill(HIST("hpTAccepted_T_z_z_E_S_U"), mcParticle.pt());
                }
              }
            }
          }
        }
      }
    }
  }

  template <typename TMCParticle, typename TMCParticles>
  bool isGammaGammaDecay(TMCParticle mcParticle, TMCParticles mcParticles)
  {
    auto daughtersIds = mcParticle.daughtersIds();
    if (daughtersIds.size() != 2)
      return false;
    for (auto& daughterId : daughtersIds) {
      if (mcParticles.iteratorAt(daughterId).pdgCode() != 22)
        return false;
    }
    return true;
  }

  template <typename TMCParticle, typename TMCParticles>
  bool isAccepted(TMCParticle mcParticle, TMCParticles mcParticles)
  {
    auto daughtersIds = mcParticle.daughtersIds();
    if (daughtersIds.size() != 2)
      return false;
    for (auto& daughterId : daughtersIds) {
      if (mcParticles.iteratorAt(daughterId).pdgCode() != 22)
        return false;
      int iCellID = -1;
      try {
        iCellID = emcal::Geometry::GetInstance()->GetAbsCellIdFromEtaPhi(mcParticles.iteratorAt(daughterId).eta(), mcParticles.iteratorAt(daughterId).phi());
      } catch (emcal::InvalidPositionException& e) {
        iCellID = -1;
      }
      if (iCellID == -1)
        return false;
    }
    return true;
  }

  void fillEventHistogram(HistogramRegistry* fRegistry, MyMCCollisions::iterator const& collision)
  {
    fRegistry->fill(HIST("hEMCollisionCounter"), 1.0);
    if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))
      fRegistry->fill(HIST("hEMCollisionCounter"), 2.0);
    if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))
      fRegistry->fill(HIST("hEMCollisionCounter"), 3.0);
    if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
      fRegistry->fill(HIST("hEMCollisionCounter"), 4.0);
    if (collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
      fRegistry->fill(HIST("hEMCollisionCounter"), 5.0);
    if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
      fRegistry->fill(HIST("hEMCollisionCounter"), 6.0);
    if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX))
      fRegistry->fill(HIST("hEMCollisionCounter"), 7.0);
    if (collision.sel8())
      fRegistry->fill(HIST("hEMCollisionCounter"), 8.0);
    if (abs(collision.posZ()) < 10.0)
      fRegistry->fill(HIST("hEMCollisionCounter"), 9.0);
    if (collision.alias_bit(kTVXinEMC))
      fRegistry->fill(HIST("hEMCollisionCounter"), 10.0);
    if (collision.alias_bit(kEMC7) || collision.alias_bit(kDMC7))
      fRegistry->fill(HIST("hEMCollisionCounter"), 11.0);
    if (collision.isemcreadout())
      fRegistry->fill(HIST("hEMCollisionCounter"), 12.0);
    fRegistry->fill(HIST("hEMCollisionCounter"), 13.0);

    fRegistry->fill(HIST("hCollisionCounter"), 1);
    if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      fRegistry->fill(HIST("hCollisionCounter"), 2);
      if (abs(collision.posZ()) < mVertexCut) {
        fRegistry->fill(HIST("hCollisionCounter"), 3);
        if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
          fRegistry->fill(HIST("hCollisionCounter"), 4);
        if (mRequireEMCCellContent ? collision.isemcreadout() : collision.alias_bit(kTVXinEMC)) {
          fRegistry->fill(HIST("hCollisionCounter"), 5);
          if (collision.sel8()) {
            fRegistry->fill(HIST("hCollisionCounter"), 6);
            if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
              fRegistry->fill(HIST("hCollisionCounter"), 7);
            }
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MCGeneratorStudies>(cfgc, TaskName{"mc-generator-studies"})};
}
