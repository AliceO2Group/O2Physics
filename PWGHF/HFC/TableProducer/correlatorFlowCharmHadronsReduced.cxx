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

/// \file correlatorFlowCharmHadronsReduced.cxx
/// \brief CharmHadrons-Hadrons correlator tree creator for data analyses
/// \author Marcello Di Costanzo <marcello.di.costanzo@cern.ch>, Politecnico and INFN Torino

#include "PWGHF/HFC/DataModel/DerivedDataCorrelationTables.h"

#include "Common/Core/RecoDecay.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TObject.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using BinningCentPosZ = ColumnBinningPolicy<aod::hf_collisions_reduced::PosZ, aod::hf_collisions_reduced::Centrality>;
using BinningMultPosZ = ColumnBinningPolicy<aod::hf_collisions_reduced::PosZ, aod::hf_collisions_reduced::Multiplicity>;

/// Get charm candidate or hadron track pT
/// \param track is the candidate
template <typename TTrack>
double getPt(const TTrack& track)
{
  if constexpr (requires { track.ptAssoc(); }) {
    return track.ptAssoc();
  } else {
    return track.ptTrig();
  }
}

/// Get charm candidate or hadron track eta
/// \param track is the candidate
template <typename TTrack>
double getEta(const TTrack& track)
{
  if constexpr (requires { track.etaAssoc(); }) {
    return track.etaAssoc();
  } else {
    return track.etaTrig();
  }
}

/// Get charm candidate or hadron track phi
/// \param track is the candidate
template <typename TTrack>
double getPhi(const TTrack& track)
{
  if constexpr (requires { track.phiAssoc(); }) {
    return track.phiAssoc();
  } else {
    return track.phiTrig();
  }
}

enum HistoEffType {
  kTH1 = 0,
  kTH2,
  kTH3,
  kTHnSparse
};

struct HfCorrelatorFlowCharmHadronsReduced {

  Configurable<int> numberEventsMixed{"numberEventsMixed", 5, "Number of events mixed in ME process"};
  Configurable<std::vector<double>> bdtScore0PtMaxs{"bdtScore0PtMaxs", std::vector<double>{0.1, 0.1, 0.1, 0.1, 0.1}, "pT-differential maximum score 0 for charm candidates"};
  Configurable<std::vector<double>> bdtScore1PtMins{"bdtScore1PtMins", std::vector<double>{0.1, 0.1, 0.1, 0.1, 0.1}, "pT-differential minimum score 1 for charm candidates"};
  Configurable<bool> selCentrality{"selCentrality", false, "Use centrality to select events"};
  Configurable<bool> selMultiplicity{"selMultiplicity", false, "Use multiplicity to select events"};
  Configurable<float> centMin{"centMin", 0, "min. centrality"};
  Configurable<float> centMax{"centMax", 10., "max. centrality"};
  Configurable<float> multMin{"multMin", 0, "min. multiplicity"};
  Configurable<float> multMax{"multMax", 10., "max. multiplicity"};
  Configurable<float> deltaEtaAbsMin{"deltaEtaAbsMin", 0.5, "min. pair delta eta"};
  Configurable<float> deltaEtaAbsMax{"deltaEtaAbsMax", 2., "max. pair delta eta"};
  Configurable<float> dcaXYTrackMax{"dcaXYTrackMax", 1., "max. track DCA XY"};
  Configurable<float> dcaZTrackMax{"dcaZTrackMax", 1., "max. track DCA Z"};
  Configurable<int> tpcCrossedRowsMin{"tpcCrossedRowsMin", 1, "min. TPC crossed rows"};
  Configurable<int> itsNClsMin{"itsNClsMin", 1, "min. ITS clusters"};
  Configurable<float> downSamplePairs{"downSamplePairs", 1., "Fraction of pairs to keep"};
  Configurable<float> ptMaxForDownSample{"ptMaxForDownSample", 10., "Maximum pt for the application of the downsampling factor"};
  Configurable<float> centMultMaxForDownSample{"centMultMaxForDownSample", 10., "Maximum centrality for the application of the downsampling factor"};

  // Efficiency correction
  Configurable<bool> doEfficiencyCorrection{"doEfficiencyCorrection", false, "Flag for applying efficiency corrections"};
  Configurable<int> effHistTypeCharm{"effHistTypeCharm", 0, "Histogram type for charm efficiency (0=TH1, 1=TH2, 2=TH3, 3=THnSparse)"};
  Configurable<int> effHistTypeAssoc{"effHistTypeAssoc", 0, "Histogram type for associated efficiency (0=TH1, 1=TH2, 2=TH3, 3=THnSparse)"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> assocEffCcdbPath{"assocEffCcdbPath", "", "CCDB path for associated efficiency"};
  Configurable<std::string> promptEffCcdbPath{"promptEffCcdbPath", "", "CCDB path for trigger efficiency"};
  Configurable<int64_t> timestampCcdb{"timestampCcdb", -1, "timestamp of the efficiency files used to query in CCDB"};
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  SliceCache cache;

  Service<ccdb::BasicCCDBManager> ccdb{};
  std::shared_ptr<TObject> hEffPrompt = nullptr;
  std::shared_ptr<TObject> hEffAssoc = nullptr;

  int poolBins{0};
  std::vector<double> mPtEdgesTrig;

  using SameEvtPairsChHadWEffInfo = soa::Filtered<soa::Join<aod::HfcRedSEChBases, aod::HfcRedAssBases, aod::HfcRedAssTracks>>;
  using SameEvtPairsChHad = soa::Filtered<soa::Join<aod::HfcRedSEChBases, aod::HfcRedAssTracks>>;
  using SameEvtPairsHadHad = soa::Filtered<soa::Join<aod::HfcRedSEHadBases, aod::HfcRedAssTracks>>;
  using AssocTracks = soa::Filtered<soa::Join<aod::HfcRedAssBases, aod::HfcRedAssTracks>>;
  using TrigCharmCands = soa::Join<aod::HfcRedTrigBases, aod::HfcRedTrigCharms>;

  Filter filterAssocTracks = (nabs(aod::hf_correl_charm_had_reduced::dcaXYAssoc) < dcaXYTrackMax) && (nabs(aod::hf_correl_charm_had_reduced::dcaZAssoc) < dcaZTrackMax) && (aod::hf_correl_charm_had_reduced::nTpcCrossedRowsAssoc > tpcCrossedRowsMin) && (aod::hf_correl_charm_had_reduced::itsNClsAssoc > itsNClsMin);
  Filter filterTrigTracks = (nabs(aod::hf_correl_charm_had_reduced::dcaXYTrig) < dcaXYTrackMax) && (nabs(aod::hf_correl_charm_had_reduced::dcaZTrig) < dcaZTrackMax) && (aod::hf_correl_charm_had_reduced::nTpcCrossedRowsTrig > tpcCrossedRowsMin) && (aod::hf_correl_charm_had_reduced::itsNClsTrig > itsNClsMin);
  Filter filterSameEvtPairs = (nabs(aod::hf_correl_charm_had_reduced::deltaEta) > deltaEtaAbsMin) && (nabs(aod::hf_correl_charm_had_reduced::deltaEta) < deltaEtaAbsMax);

  Preslice<AssocTracks> assocTracksPerCol = aod::hf_correl_charm_had_reduced::hfcRedCorrCollId;
  Preslice<TrigCharmCands> trigCharmCandsPerCol = aod::hf_correl_charm_had_reduced::hfcRedCorrCollId;

  ConfigurableAxis binsPtTrig{"binsPtTrig", {VARIABLE_WIDTH, 0., 1., 2., 3., 4., 5., 6., 8., 10., 12., 16., 24.}, "Trigger candidate pT bins"};
  ConfigurableAxis binsPtAssoc{"binsPtAssoc", {VARIABLE_WIDTH, 0.3, 1., 2., 50.}, "Associated particle pT bins"};
  ConfigurableAxis binsZPool{"binsZPool", {VARIABLE_WIDTH, -10.0, -2.5, 2.5, 10.0}, "Z vertex position pools"};
  ConfigurableAxis binsMultPool{"binsMultPool", {VARIABLE_WIDTH, 0., 900., 1800., 6000.}, "Event multiplicity pools (FT0M)"};
  ConfigurableAxis binsCentPool{"binsCentPool", {VARIABLE_WIDTH, 0., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100}, "Event centrality pools"};
  ConfigurableAxis binsInvMass{"binsInvMass", {300, 1.6, 2.2}, "Invariant mass bins"};
  ConfigurableAxis binsMultFT0M{"binsMultFT0M", {100, 0., 10000.}, "Multiplicity as FT0M signal amplitude"};
  ConfigurableAxis binsCent{"binsCent", {100, 0., 100.}, "Centrality bins"};
  ConfigurableAxis binsPosZ{"binsPosZ", {100, -10., 10.}, "Primary vertex z coordinate"};
  ConfigurableAxis binsDeltaEta{"binsDeltaEta", {100, -2., 2.}, "Delta Eta bins"};
  ConfigurableAxis binsDeltaPhi{"binsDeltaPhi", {64, -o2::constants::math::PIHalf, 3. * o2::constants::math::PIHalf}, "Delta Phi bins"};
  ConfigurableAxis binsMlOne{"binsMlOne", {100, 0., 1.}, "ML score index 1 bins"};
  ConfigurableAxis binsMlTwo{"binsMlTwo", {100, 0., 1.}, "ML score index 2 bins"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  std::shared_ptr<TObject> loadEfficiencyHistogram(const std::string& path, int histoType)
  {
    TObject* object = nullptr;

    switch (histoType) {
      case kTH1:
        object = ccdb->getForTimeStamp<TH1F>(Form("%s/TH1/", path.c_str()), timestampCcdb);
        break;
      case kTH2:
        object = ccdb->getForTimeStamp<TH2F>(Form("%s/TH2/", path.c_str()), timestampCcdb);
        break;
      case kTH3:
        object = ccdb->getForTimeStamp<TH3F>(Form("%s/TH3/", path.c_str()), timestampCcdb);
        break;
      case kTHnSparse:
        object = ccdb->getForTimeStamp<THnSparseF>(Form("%s/THnSparse/", path.c_str()), timestampCcdb);
        break;
      default:
        LOGF(fatal, "Invalid efficiency histogram type: %d", histoType);
    }

    if (!object) {
      LOGF(fatal, "Could not load efficiency histogram from %s", path.c_str());
    }

    return std::shared_ptr<TObject>(object);
  }

  void init(InitContext&)
  {
    if (selCentrality && selMultiplicity) {
      LOGP(fatal, "You cannot select events based on both centrality and multiplicity! Please check your configuration!");
    }
    if (!selCentrality && !selMultiplicity) {
      LOGP(fatal, "You must select events based on either centrality or multiplicity! Please check your configuration!");
    }

    if ((doprocessSECharmHadWCentMix && doprocessMECharmHadWMultMix) ||
        (doprocessSECharmHadWMultMix && doprocessMECharmHadWCentMix) ||
        (doprocessSEHadHadWCentMix && doprocessMEHadHadWMultMix) ||
        (doprocessSEHadHadWMultMix && doprocessMEHadHadWCentMix)) {
      LOGP(fatal, "You cannot mix centrality and multiplicity mixing in the same processing! Please check your configuration!");
    }

    mPtEdgesTrig.assign(binsPtTrig.value.begin() + 1, binsPtTrig.value.end());
    if ((mPtEdgesTrig.size() != (bdtScore0PtMaxs.value.size() + 1) ||
         mPtEdgesTrig.size() != (bdtScore1PtMins.value.size() + 1))) {
      LOG(info) << "mPtEdgesTrig.size(): " << mPtEdgesTrig.size();
      LOG(info) << "bdtScore0PtMaxs.value.size(): " << bdtScore0PtMaxs.value.size();
      LOG(info) << "bdtScore1PtMins.value.size(): " << bdtScore1PtMins.value.size();
      LOGP(fatal, "The size of bdtScore0PtMaxs and bdtScore1PtMins must be the one of binsPtTrig minus one!");
    }

    if (doprocessSECharmHadWCentMix || doprocessSEHadHadWCentMix || doprocessMECharmHadWCentMix || doprocessMEHadHadWCentMix) {
      poolBins = (binsCentPool->size() - 2) * (binsZPool->size() - 2);
    } else {
      poolBins = (binsMultPool->size() - 2) * (binsZPool->size() - 2);
    }

    const AxisSpec axisInvMass{binsInvMass, "Inv. mass (GeV/#it{c}^{2})"};
    const AxisSpec axisCent = {binsCent, "Centrality"};
    const AxisSpec axisMultFT0M = {binsMultFT0M, "MultiplicityFT0M"};
    const AxisSpec axisPosZ = {binsPosZ, "PosZ"};
    const AxisSpec axisPoolBin = {poolBins, 0., static_cast<float>(poolBins), "PoolBin"};
    const AxisSpec axisDeltaEta = {binsDeltaEta, "#Delta#it{#eta}"};
    const AxisSpec axisDeltaPhi = {binsDeltaPhi, "#Delta#it{#varphi}"};
    const AxisSpec axisPtTrig = {binsPtTrig, "#it{p}_{T} Trig (GeV/#it{c})"};
    const AxisSpec axisPtAssoc = {binsPtAssoc, "#it{p}_{T} Assoc (GeV/#it{c})"};
    const AxisSpec axisMlOne{binsMlOne, "bdtScore0"};
    const AxisSpec axisMlTwo{binsMlTwo, "bdtScore1"};

    // Histograms for data analysis
    if (doprocessSECharmHadWCentMix || doprocessSECharmHadWCentMixWEffInfo || doprocessSEHadHadWCentMix) {
      registry.add("hCentPoolBinSE", "Centrality SE", {HistType::kTH2F, {{axisCent}, {axisPoolBin}}});
    } else if (doprocessSECharmHadWMultMix || doprocessSECharmHadWMultMixWEffInfo || doprocessSEHadHadWMultMix) {
      registry.add("hMultFT0MPoolBinSE", "Multiplicity FT0M SE", {HistType::kTH2F, {{axisMultFT0M}, {axisPoolBin}}});
    } else if (doprocessMECharmHadWCentMix || doprocessMECharmHadWCentMixWEffInfo || doprocessMEHadHadWCentMix) {
      registry.add("hCentPoolBinME", "Centrality ME", {HistType::kTH2F, {{axisCent}, {axisPoolBin}}});
    } else if (doprocessMECharmHadWMultMix || doprocessMECharmHadWMultMixWEffInfo || doprocessMEHadHadWMultMix) {
      registry.add("hMultFT0MPoolBinME", "Multiplicity FT0M ME", {HistType::kTH2F, {{axisMultFT0M}, {axisPoolBin}}});
    }
    registry.add("hZVtxPoolBinSE", "z vertex SE", {HistType::kTH2F, {{axisPosZ}, {axisPoolBin}}});
    registry.add("hZVtxPoolBinME", "z vertex ME", {HistType::kTH2F, {{axisPosZ}, {axisPoolBin}}});
    registry.add("hCollisionsCentMultAll", "Collision centrality vs multiplicity", {HistType::kTH2F, {{axisCent}, {axisMultFT0M}}});
    registry.add("hCollisionsCentMultSel", "Collision centrality vs multiplicity", {HistType::kTH2F, {{axisCent}, {axisMultFT0M}}});
    registry.add("hPoolBinTrigSE", "Trigger candidates pool bin SE", {HistType::kTH1F, {axisPoolBin}});
    registry.add("hPoolBinTrigME", "Trigger candidates pool bin ME", {HistType::kTH1F, {axisPoolBin}});
    registry.add("hPoolBinAssocSE", "Associated particles pool bin SE", {HistType::kTH1F, {axisPoolBin}});
    registry.add("hPoolBinAssocME", "Associated particles pool bin ME", {HistType::kTH1F, {axisPoolBin}});

    std::vector<AxisSpec> axesTrigger = {axisInvMass, axisPtTrig, axisMlOne, axisMlTwo};
    std::vector<AxisSpec> axes = {axisPoolBin, axisPtTrig, axisPtAssoc, axisDeltaEta, axisDeltaPhi};
    if (doprocessSEHadHadWCentMix || doprocessSEHadHadWMultMix) {
      registry.add("hSparseCorrelationsSEHadHad", "THn for SE Had-Had correlations", HistType::kTHnSparseF, axes);
      if (doEfficiencyCorrection) {
        registry.get<THnSparse>(HIST("hSparseCorrelationsSEHadHad"))->Sumw2();
      }
    } else if (doprocessMEHadHadWCentMix || doprocessMEHadHadWMultMix) {
      registry.add("hSparseCorrelationsMEHadHad", "THn for ME Had-Had correlations", HistType::kTHnSparseF, axes);
      if (doEfficiencyCorrection) {
        registry.get<THnSparse>(HIST("hSparseCorrelationsMEHadHad"))->Sumw2();
      }
    } else {
      axes.insert(axes.end(), {axisInvMass, axisMlOne, axisMlTwo});
      if (doprocessSECharmHadWCentMix || doprocessSECharmHadWCentMixWEffInfo || doprocessSECharmHadWMultMix || doprocessSECharmHadWMultMixWEffInfo) {
        registry.add("hSparseCorrelationsSECharmHad", "THn for SE Charm-Had correlations", HistType::kTHnSparseF, axes);
        if (doEfficiencyCorrection) {
          registry.get<THnSparse>(HIST("hSparseCorrelationsSECharmHad"))->Sumw2();
        }
      } else if (doprocessMECharmHadWCentMix || doprocessMECharmHadWMultMixWEffInfo || doprocessMECharmHadWMultMix || doprocessMECharmHadWCentMixWEffInfo) {
        registry.add("hSparseCorrelationsMECharmHad", "THn for ME Charm-Had correlations", HistType::kTHnSparseF, axes);
        if (doEfficiencyCorrection) {
          registry.get<THnSparse>(HIST("hSparseCorrelationsMECharmHad"))->Sumw2();
        }
      }
      if (doprocessCharmTriggers || doprocessCharmTriggersWEffInfo) {
        registry.add("hSparseTrigCandsCharm", "THn for Charm trigger candidates", HistType::kTHnSparseF, axesTrigger);
        if (doEfficiencyCorrection) {
          registry.get<THnSparse>(HIST("hSparseTrigCandsCharm"))->Sumw2();
        }
      }
    }

    if (doEfficiencyCorrection) {
      hEffPrompt = loadEfficiencyHistogram(promptEffCcdbPath, effHistTypeCharm);
      hEffAssoc = loadEfficiencyHistogram(assocEffCcdbPath, effHistTypeAssoc);
    }
  }

  template <bool IsCharm, typename TCand>
  double getEfficiency(TObject* hist, int histoType, const TCand& cand) //, const std::vector<double>& values)
  {
    if (!hist) {
      LOG(debug) << "Candidate does not have required efficiency histogram.";
      return 0.;
    }

    std::vector<double> kineFeatures{};

    // Charm triggers
    if constexpr (IsCharm) {
      kineFeatures = {cand.ptTrig(), cand.bdtScore0Trig(), cand.bdtScore1Trig(), cand.etaTrig(), cand.phiTrig()};
    } else {
      auto collision = cand.template hfcRedCorrColl_as<o2::aod::HfcRedCorrColls>();
      kineFeatures = {cand.ptAssoc(), cand.etaAssoc(), cand.phiAssoc(), collision.posZ()};
    }

    switch (histoType) {

      case kTH1: {
        auto* h = dynamic_cast<TH1*>(hist);
        if (!h || kineFeatures.size() < 1) {
          return 0.;
        }
        return h->GetBinContent(h->FindBin(kineFeatures[0]));
      }

      case kTH2: {
        auto* h = dynamic_cast<TH2*>(hist);
        if (!h || kineFeatures.size() < 2) {
          return 0.;
        }
        return h->GetBinContent(h->FindBin(kineFeatures[0], kineFeatures[1]));
      }

      case kTH3: {
        auto* h = dynamic_cast<TH3*>(hist);
        if (!h || kineFeatures.size() < 3) {
          return 0.;
        }
        return h->GetBinContent(
          h->FindBin(kineFeatures[0], kineFeatures[1], kineFeatures[2]));
      }

      case kTHnSparse: {
        auto* h = dynamic_cast<THnSparse*>(hist);
        if (!h || kineFeatures.size() < static_cast<size_t>(h->GetNdimensions())) {
          return 0.;
        }

        std::vector<Int_t> bins(h->GetNdimensions());

        for (int i = 0; i < h->GetNdimensions(); ++i) {
          bins[i] = h->GetAxis(i)->FindBin(kineFeatures[i]);
        }

        return h->GetBinContent(bins.data());
      }

      default:
        return 0.;
    }
  }

  bool isValidEfficiency(double efficiency)
  {
    return std::isfinite(efficiency) && efficiency > 0.;
  }

  template <typename TColl>
  bool isSelCollision(const TColl& collision)
  {
    bool isSel = false;
    if (selCentrality) {
      isSel = (collision.centrality() >= centMin && collision.centrality() <= centMax);
    } else if (selMultiplicity) {
      isSel = (collision.multiplicity() >= multMin && collision.multiplicity() <= multMax);
    }
    return isSel;
  }

  /// Get the binning pool associated to the collision
  /// \param collision is the collision
  /// \param binPolicy is the binning policy for the correlation
  template <bool IsMixedEvent, typename TColl, typename TBinningType>
  int getPoolBin(const TColl& collision, const TBinningType& binPolicy)
  {
    int poolBin{0};
    if constexpr (std::is_same_v<TBinningType, BinningCentPosZ>) {
      poolBin = binPolicy.getBin(std::make_tuple(collision.posZ(), collision.centrality()));
      if constexpr (IsMixedEvent) {
        registry.fill(HIST("hCentPoolBinME"), collision.centrality(), poolBin);
        registry.fill(HIST("hZVtxPoolBinME"), collision.posZ(), poolBin);
      } else {
        registry.fill(HIST("hCentPoolBinSE"), collision.centrality(), poolBin);
        registry.fill(HIST("hZVtxPoolBinSE"), collision.posZ(), poolBin);
      }
    } else if constexpr (std::is_same_v<TBinningType, BinningMultPosZ>) {
      poolBin = binPolicy.getBin(std::make_tuple(collision.posZ(), collision.multiplicity()));
      if constexpr (IsMixedEvent) {
        registry.fill(HIST("hMultFT0MPoolBinME"), collision.multiplicity(), poolBin);
        registry.fill(HIST("hZVtxPoolBinME"), collision.posZ(), poolBin);
      } else {
        registry.fill(HIST("hMultFT0MPoolBinSE"), collision.multiplicity(), poolBin);
        registry.fill(HIST("hZVtxPoolBinSE"), collision.posZ(), poolBin);
      }
    }
    return poolBin;
  }

  /// Apply pT-differential ML BDT bkg score cut
  /// \param ptTrig is the pT of the charm candidate
  template <typename TCand>
  bool isSelBdtScoreCut(TCand const& cand,
                        double ptTrig)
  {
    for (size_t iPt = 0; iPt < mPtEdgesTrig.size() - 1; iPt++) {
      if (ptTrig >= mPtEdgesTrig[iPt] && ptTrig < mPtEdgesTrig[iPt + 1]) {
        return (cand.bdtScore0Trig() < bdtScore0PtMaxs.value[iPt]) && (cand.bdtScore1Trig() > bdtScore1PtMins.value[iPt]);
      }
    }
    return false;
  }

  /// Save info for Same Event pairs
  /// \param collisions are the selected collisions
  /// \param trigCands are the selected trigger candidates
  /// \param assocTracks are the selected associated tracks
  /// \param binPolicy is the binning policy for the correlation
  template <bool CorrectEff, typename TPair, typename TTrigCand, typename TBinningType>
  void fillSameEvent(TPair const& pair,
                     TTrigCand const& trigCand,
                     const TBinningType& binPolicy)
  {
    auto collision = pair.template hfcRedCorrColl_as<o2::aod::HfcRedCorrColls>();
    registry.fill(HIST("hCollisionsCentMultAll"), collision.centrality(), collision.multiplicity());
    if (!isSelCollision(collision)) {
      return;
    }
    registry.fill(HIST("hCollisionsCentMultSel"), collision.centrality(), collision.multiplicity());

    double const ptTrig = trigCand.ptTrig();
    if constexpr (requires { trigCand.bdtScore0Trig(); }) { // ML selection on bkg score for Charm-Had case
      if (!isSelBdtScoreCut(trigCand, ptTrig)) {
        return;
      }
    }
    if (downSamplePairs < 1.) {
      float const pseudoRndm = ptTrig * 1000. - static_cast<int64_t>(ptTrig * 1000);
      float collCentMult = selCentrality ? collision.centrality() : collision.multiplicity();
      if (ptTrig < ptMaxForDownSample && collCentMult < centMultMaxForDownSample && pseudoRndm >= downSamplePairs) {
        return;
      }
    }
    int const poolBin = getPoolBin<false>(collision, binPolicy);
    registry.fill(HIST("hPoolBinTrigSE"), poolBin);
    registry.fill(HIST("hPoolBinAssocSE"), poolBin);

    double effAssoc{1.f};
    if constexpr (CorrectEff) {
      if (doEfficiencyCorrection) {
        effAssoc = getEfficiency<false>(hEffAssoc.get(), effHistTypeAssoc.value, pair);
        if (!isValidEfficiency(effAssoc)) {
          LOGF(warn, "Invalid efficiency value for associated track: effAssoc = %f", effAssoc);
          return;
        }
      }
    }

    if constexpr (requires { trigCand.bdtScore0Trig(); }) { // Separate Charm-Had and Had-Had cases

      double effTrig{1.f};
      if constexpr (CorrectEff) {
        if (doEfficiencyCorrection) {
          effTrig = getEfficiency<true>(hEffPrompt.get(), effHistTypeCharm.value, trigCand);
          if (!isValidEfficiency(effTrig)) {
            LOGF(warn, "Invalid efficiency value for trigger candidate: effTrig = %f", effTrig);
            return;
          }
        }
      }

      double weight = doEfficiencyCorrection ? 1 / (effTrig * effAssoc) : 1.f;

      std::vector<double> values = {static_cast<double>(poolBin), ptTrig, pair.ptAssoc(), pair.deltaEta(), pair.deltaPhi(), trigCand.invMassTrig(), trigCand.bdtScore0Trig(), trigCand.bdtScore1Trig()};

      registry.get<THnSparse>(HIST("hSparseCorrelationsSECharmHad"))->Fill(values.data(), weight);
    } else {

      double effTrig{1.f};
      if constexpr (CorrectEff) {
        if (doEfficiencyCorrection) {
          effTrig = getEfficiency<false>(hEffAssoc.get(), effHistTypeAssoc.value, trigCand);
          if (!isValidEfficiency(effTrig)) {
            LOGF(warn, "Invalid efficiency value for trigger candidate: effTrig = %f", effTrig);
            return;
          }
        }
      }

      double weight = doEfficiencyCorrection ? 1 / (effTrig * effAssoc) : 1.f;

      std::vector<double> valuesHadHad = {static_cast<double>(poolBin), ptTrig, pair.ptAssoc(), pair.deltaEta(), pair.deltaPhi()};
      registry.get<THnSparse>(HIST("hSparseCorrelationsSEHadHad"))->Fill(valuesHadHad.data(), weight);
    }
  }

  /// Save info for Mixed Event pairs
  /// \param collisions are the selected collisions
  /// \param pairs are the mixed event pairs of trigger candidates and associated tracks
  /// \param binPolicy is the binning policy for the correlation
  template <bool CorrectEff, typename TPairs, typename TBinningType>
  void fillMixedEvent(TPairs const& pairs,
                      const TBinningType& binPolicy)
  {
    for (const auto& [trigColl, trigCands, assocColl, assocTracks] : pairs) {
      if (trigCands.size() == 0 || assocTracks.size() == 0) {
        continue;
      }
      registry.fill(HIST("hCollisionsCentMultAll"), trigColl.centrality(), trigColl.multiplicity());
      registry.fill(HIST("hCollisionsCentMultAll"), assocColl.centrality(), assocColl.multiplicity());
      if (!isSelCollision(trigColl) || !isSelCollision(assocColl)) {
        continue;
      }
      registry.fill(HIST("hCollisionsCentMultSel"), trigColl.centrality(), trigColl.multiplicity());
      registry.fill(HIST("hCollisionsCentMultSel"), assocColl.centrality(), assocColl.multiplicity());

      int const poolBinTrig = getPoolBin<true>(trigColl, binPolicy);
      int const poolBinAssoc = getPoolBin<true>(assocColl, binPolicy);
      if (poolBinAssoc != poolBinTrig) {
        LOGF(info, "Error, poolBins are different");
        continue;
      }
      registry.fill(HIST("hPoolBinTrigME"), poolBinTrig);
      registry.fill(HIST("hPoolBinAssocME"), poolBinAssoc);

      for (const auto& [trigCand, assocTrack] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(trigCands, assocTracks))) {
        // LOGF(info, "Mixed event tracks pair: (%d, %d) from events (%d, %d), track event: (%d, %d)", trigCand.index(), assocTrack.index(), trigColl.index(), assocColl.index(), trigCand.hfcRedFlowCollId(), assocTrack.hfcRedFlowCollId());
        double const deltaEta = getEta(assocTrack) - getEta(trigCand);
        if (std::abs(deltaEta) < deltaEtaAbsMin || std::abs(deltaEta) > deltaEtaAbsMax) {
          continue;
        }
        double const ptTrig = getPt(trigCand);
        if constexpr (requires { trigCand.bdtScore0Trig(); }) { // ML selection on bkg score for Charm-Had case
          if (!isSelBdtScoreCut(trigCand, ptTrig)) {
            continue;
          }
        }
        double const ptAssoc = getPt(assocTrack);
        if (downSamplePairs < 1.) {
          float const pseudoRndm = ptAssoc * 1000. - static_cast<int64_t>(ptAssoc * 1000);
          float trigCollCentMult = selCentrality ? trigColl.centrality() : trigColl.multiplicity();
          float assocCollCentMult = selCentrality ? assocColl.centrality() : assocColl.multiplicity();
          if (ptTrig < ptMaxForDownSample && trigCollCentMult < centMultMaxForDownSample &&
              assocCollCentMult < centMultMaxForDownSample && pseudoRndm >= downSamplePairs) {
            continue;
          }
        }
        double const deltaPhi = RecoDecay::constrainAngle(getPhi(assocTrack) - getPhi(trigCand), -o2::constants::math::PIHalf);

        double effAssoc{1.f};
        if (doEfficiencyCorrection) {
          effAssoc = getEfficiency<false>(hEffAssoc.get(), effHistTypeAssoc.value, assocTrack);
          if (!isValidEfficiency(effAssoc)) {
            LOGF(warn, "Invalid efficiency value for associated track: effAssoc = %f", effAssoc);
            continue;
          }
        }

        if constexpr (requires { trigCand.bdtScore0Trig(); }) { // Separate Charm-Had and Had-Had cases

          double effTrig{1.f};
          if (doEfficiencyCorrection) {
            effTrig = getEfficiency<true>(hEffPrompt.get(), effHistTypeCharm.value, trigCand);
            if (!isValidEfficiency(effTrig)) {
              LOGF(warn, "Invalid efficiency value for trigger candidate: effTrig = %f", effTrig);
              continue;
            }
          }

          double weight = doEfficiencyCorrection ? 1 / (effTrig * effAssoc) : 1.f;

          std::vector<double> values = {static_cast<double>(poolBinTrig), ptTrig, ptAssoc, deltaEta, deltaPhi, trigCand.invMassTrig(), trigCand.bdtScore0Trig(), trigCand.bdtScore1Trig()};
          registry.get<THnSparse>(HIST("hSparseCorrelationsMECharmHad"))->Fill(values.data(), weight);
        } else {

          double effTrig{1.f};
          if (doEfficiencyCorrection) {
            effTrig = getEfficiency<false>(hEffAssoc.get(), effHistTypeAssoc.value, trigCand);
            if (!isValidEfficiency(effTrig)) {
              LOGF(warn, "Invalid efficiency value for trigger candidate: effTrig = %f", effTrig);
              continue;
            }
          }

          double weight = doEfficiencyCorrection ? 1 / (effTrig * effAssoc) : 1.f;

          std::vector<double> valuesHadHad = {static_cast<double>(poolBinTrig), ptTrig, ptAssoc, deltaEta, deltaPhi};
          registry.get<THnSparse>(HIST("hSparseCorrelationsMEHadHad"))->Fill(valuesHadHad.data(), weight);
        }
      }
    }
  }

  void processSECharmHadWMultMix(SameEvtPairsChHad::iterator const& pair,
                                 aod::HfcRedTrigCharms const&,
                                 aod::HfcRedCorrColls const&)
  {
    BinningMultPosZ binPolicyPosZMult{{binsZPool, binsMultPool}};
    auto trigCand = pair.template hfcRedTrigCharm_as<aod::HfcRedTrigCharms>();
    fillSameEvent<false>(pair, trigCand, binPolicyPosZMult);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processSECharmHadWMultMix, "Process Same Event for Charm-Had with multiplicity pools", false);

  void processSECharmHadWMultMixWEffInfo(SameEvtPairsChHadWEffInfo::iterator const& pair,
                                         TrigCharmCands const&,
                                         aod::HfcRedCorrColls const&)
  {
    BinningMultPosZ binPolicyPosZMult{{binsZPool, binsMultPool}};
    auto trigCand = pair.template hfcRedTrigCharm_as<TrigCharmCands>();
    fillSameEvent<true>(pair, trigCand, binPolicyPosZMult);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processSECharmHadWMultMixWEffInfo, "Process Same Event for Charm-Had with multiplicity pools", false);

  void processSEHadHadWMultMix(SameEvtPairsHadHad::iterator const& pair,
                               aod::HfcRedTrigTracks const&,
                               aod::HfcRedCorrColls const&)
  {
    BinningMultPosZ binPolicyPosZMult{{binsZPool, binsMultPool}};
    auto trigCand = pair.template hfcRedTrigTrack_as<aod::HfcRedTrigTracks>();
    fillSameEvent<false>(pair, trigCand, binPolicyPosZMult);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processSEHadHadWMultMix, "Process Same Event for Had-Had with multiplicity pools", false);

  void processSECharmHadWCentMix(SameEvtPairsChHad::iterator const& pair,
                                 aod::HfcRedTrigCharms const&,
                                 aod::HfcRedCorrColls const&)
  {
    BinningCentPosZ binPolicyPosZCent{{binsZPool, binsCentPool}};
    auto trigCand = pair.template hfcRedTrigCharm_as<aod::HfcRedTrigCharms>();
    fillSameEvent<false>(pair, trigCand, binPolicyPosZCent);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processSECharmHadWCentMix, "Process Same Event for Charm-Had with centrality pools", true);

  void processSECharmHadWCentMixWEffInfo(SameEvtPairsChHadWEffInfo::iterator const& pair,
                                         TrigCharmCands const&,
                                         aod::HfcRedCorrColls const&)
  {
    BinningCentPosZ binPolicyPosZCent{{binsZPool, binsCentPool}};
    auto trigCand = pair.template hfcRedTrigCharm_as<TrigCharmCands>();
    fillSameEvent<true>(pair, trigCand, binPolicyPosZCent);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processSECharmHadWCentMixWEffInfo, "Process Same Event for Charm-Had with centrality pools", false);

  void processSEHadHadWCentMix(SameEvtPairsHadHad::iterator const& pair,
                               aod::HfcRedTrigTracks const&,
                               aod::HfcRedCorrColls const&)
  {
    BinningCentPosZ binPolicyPosZCent{{binsZPool, binsCentPool}};
    auto trigCand = pair.template hfcRedTrigTrack_as<aod::HfcRedTrigTracks>();
    fillSameEvent<false>(pair, trigCand, binPolicyPosZCent);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processSEHadHadWCentMix, "Process Same Event for Had-Had with centrality pools", false);

  void processMECharmHadWCentMix(aod::HfcRedCorrColls const& collisions,
                                 TrigCharmCands const& candidates,
                                 AssocTracks const& tracks)
  {
    BinningCentPosZ binPolicyPosZCent{{binsZPool, binsCentPool}};
    auto pairsTuple = std::make_tuple(candidates, tracks);
    Pair<aod::HfcRedCorrColls, TrigCharmCands, AssocTracks, BinningCentPosZ> const pairs{binPolicyPosZCent, numberEventsMixed, -1, collisions, pairsTuple, &cache};
    fillMixedEvent<false>(pairs, binPolicyPosZCent);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processMECharmHadWCentMix, "Process Mixed Event for Charm-Had with centrality pools", false);

  void processMECharmHadWCentMixWEffInfo(aod::HfcRedCorrColls const& collisions,
                                         TrigCharmCands const& candidates,
                                         AssocTracks const& tracks)
  {
    BinningCentPosZ binPolicyPosZCent{{binsZPool, binsCentPool}};
    auto pairsTuple = std::make_tuple(candidates, tracks);
    Pair<aod::HfcRedCorrColls, TrigCharmCands, AssocTracks, BinningCentPosZ> const pairs{binPolicyPosZCent, numberEventsMixed, -1, collisions, pairsTuple, &cache};
    fillMixedEvent<true>(pairs, binPolicyPosZCent);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processMECharmHadWCentMixWEffInfo, "Process Mixed Event for Charm-Had with centrality pools and efficiency information", false);

  void processMECharmHadWMultMix(aod::HfcRedCorrColls const& collisions,
                                 TrigCharmCands const& candidates,
                                 AssocTracks const& tracks)
  {
    BinningMultPosZ binPolicyPosZMult{{binsZPool, binsMultPool}};
    auto pairsTuple = std::make_tuple(candidates, tracks);
    Pair<aod::HfcRedCorrColls, TrigCharmCands, AssocTracks, BinningMultPosZ> const pairs{binPolicyPosZMult, numberEventsMixed, -1, collisions, pairsTuple, &cache};
    fillMixedEvent<false>(pairs, binPolicyPosZMult);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processMECharmHadWMultMix, "Process Mixed Event for Charm-Had with multiplicity pools", false);

  void processMECharmHadWMultMixWEffInfo(aod::HfcRedCorrColls const& collisions,
                                         TrigCharmCands const& candidates,
                                         AssocTracks const& tracks)
  {
    BinningMultPosZ binPolicyPosZMult{{binsZPool, binsMultPool}};
    auto pairsTuple = std::make_tuple(candidates, tracks);
    Pair<aod::HfcRedCorrColls, TrigCharmCands, AssocTracks, BinningMultPosZ> const pairs{binPolicyPosZMult, numberEventsMixed, -1, collisions, pairsTuple, &cache};
    fillMixedEvent<true>(pairs, binPolicyPosZMult);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processMECharmHadWMultMixWEffInfo, "Process Mixed Event for Charm-Had with multiplicity pools and efficiency information", false);

  void processMEHadHadWCentMix(aod::HfcRedCorrColls const& collisions,
                               AssocTracks const& tracks)
  {
    BinningCentPosZ binPolicyPosZCent{{binsZPool, binsCentPool}};
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<aod::HfcRedCorrColls, AssocTracks, BinningCentPosZ> const pairs{binPolicyPosZCent, numberEventsMixed, -1, collisions, tracksTuple, &cache};
    fillMixedEvent<false>(pairs, binPolicyPosZCent);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processMEHadHadWCentMix, "Process Mixed Event for Had-Had with centrality pools", false);

  void processMEHadHadWMultMix(aod::HfcRedCorrColls const& collisions,
                               AssocTracks const& tracks)
  {
    BinningMultPosZ binPolicyPosZMult{{binsZPool, binsMultPool}};
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<aod::HfcRedCorrColls, AssocTracks, BinningMultPosZ> const pairs{binPolicyPosZMult, numberEventsMixed, -1, collisions, tracksTuple, &cache};
    fillMixedEvent<false>(pairs, binPolicyPosZMult);
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processMEHadHadWMultMix, "Process Mixed Event for Had-Had with multiplicity pools", false);

  void processCharmTriggers(aod::HfcRedTrigCharms const& trigCands,
                            aod::HfcRedCorrColls const&)
  {

    for (const auto& trigCand : trigCands) {
      auto collision = trigCand.template hfcRedCorrColl_as<o2::aod::HfcRedCorrColls>();
      if (!isSelCollision(collision)) {
        continue;
      }
      if (!isSelBdtScoreCut(trigCand, trigCand.ptTrig())) {
        continue;
      }

      std::vector<double> values = {trigCand.invMassTrig(), trigCand.ptTrig(), trigCand.bdtScore0Trig(), trigCand.bdtScore1Trig()};
      registry.get<THnSparse>(HIST("hSparseTrigCandsCharm"))->Fill(values.data());
    }
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processCharmTriggers, "Process charm trigger info", false);

  void processCharmTriggersWEffInfo(TrigCharmCands const& trigCands,
                                    aod::HfcRedCorrColls const&)
  {

    for (const auto& trigCand : trigCands) {
      auto collision = trigCand.template hfcRedCorrColl_as<o2::aod::HfcRedCorrColls>();
      if (!isSelCollision(collision)) {
        continue;
      }
      if (!isSelBdtScoreCut(trigCand, trigCand.ptTrig())) {
        continue;
      }

      double effTrig = doEfficiencyCorrection ? getEfficiency<true>(hEffPrompt.get(), effHistTypeCharm.value, trigCand) : 1.;
      if (doEfficiencyCorrection && !isValidEfficiency(effTrig)) {
        LOGF(warn, "Invalid efficiency value for trigger candidate: effTrig = %f", effTrig);
        continue;
      }
      std::vector<double> values = {trigCand.invMassTrig(), trigCand.ptTrig(), trigCand.bdtScore0Trig(), trigCand.bdtScore1Trig()};
      double weight = doEfficiencyCorrection ? 1. / effTrig : 1.;
      registry.get<THnSparse>(HIST("hSparseTrigCandsCharm"))->Fill(values.data(), weight);
    }
  }
  PROCESS_SWITCH(HfCorrelatorFlowCharmHadronsReduced, processCharmTriggersWEffInfo, "Process charm trigger info with efficiency", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfCorrelatorFlowCharmHadronsReduced>(cfgc)};
}
