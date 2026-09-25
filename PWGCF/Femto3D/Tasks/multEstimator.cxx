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
//
/// \brief create a table applying some basic cuts on the ITS and DCA.
/// \author Sofia Tomassini, Gleb Romanenko, Nicolò Jacazio
/// \since 31 May 2023

#include "PWGCF/Femto3D/DataModel/PIDutils.h"
#include "PWGCF/Femto3D/DataModel/singletrackselector.h"

#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>

#include <array>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;
using namespace o2::aod;
//::singletrackselector; // the namespace defined in .h

struct multEstimator {
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<bool> CBThadronPID{"CBThadronPID", false, "Apply ev. sel. based on RCT flag `hadronPID`"}; // more in Common/CCDB/RCTSelectionFlags.h
  Configurable<int> centTableToUse{"centTableToUse", 1, "Flag to choose cent./mult.perc. estimator (Run3 only [FTOC for PbPb; FTOM for pp], for Run2 the V0M is used): 0 -> CentFV0As, 1 -> CentFT0Ms, 2 -> CentFT0As, 3 -> CentFT0Cs, 4 -> CentFDDMs, 5 -> CentNTPVs"};

  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply rapidity cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<float> _eta{"eta", 0.5, "abs eta value limit"};
  Configurable<float> _vertexZ{"VertexZ", 10.0, "abs vertexZ value limit"};

  // the cuts below are used to count and for the conditional mult ...

  Configurable<float> _maxTofChi2{"maxTofChi2", 10.f, "Maximum TOF Chi2 value -> to remove mismatched tracks"};

  Configurable<bool> rejectNotPropagatedTrks{"rejectNotPropagatedTrks", true, "rejects tracks that are not propagated to the primary vertex"};
  Configurable<bool> _removeSameBunchPileup{"removeSameBunchPileup", false, ""};
  Configurable<bool> _requestGoodZvtxFT0vsPV{"requestGoodZvtxFT0vsPV", false, ""};
  Configurable<bool> _requestVertexITSTPC{"requestVertexITSTPC", false, ""};
  Configurable<int> _requestVertexTOForTRDmatched{"requestVertexTOFmatched", 0, "0 -> no selectio; 1 -> vertex is matched to TOF or TRD; 2 -> matched to both;"};
  Configurable<bool> _requestNoCollInTimeRangeStandard{"requestNoCollInTimeRangeStandard", false, ""};
  Configurable<bool> _requestIsGoodITSLayersAll{"requestIsGoodITSLayersAll", false, "cut time intervals with dead ITS staves"};
  // Configurable<bool> fetchRate{"fetchRate", true, "Fetch the hadronic rate from the CCDB"};
  Configurable<std::pair<float, float>> _IRcut{"IRcut", std::pair<float, float>{0.f, 100.f}, "[min., max.] IR range to keep events within"};
  Configurable<std::pair<int, int>> _OccupancyCut{"OccupancyCut", std::pair<int, int>{0, 10000}, "[min., max.] occupancy range to keep events within"};

  Configurable<int> _sign{"sign", 1, "sign of a track"};
  Configurable<float> _min_P{"min_P", 0.0, "lower mometum limit"};
  Configurable<float> _max_P{"max_P", 100.0, "upper mometum limit"};
  Configurable<std::vector<float>> _dcaXY{"dcaXY", std::vector<float>{0.3f, 0.0f, 0.0f}, "abs dcaXY value limit; formula: [0] + [1]*pT^[2]"};
  Configurable<std::vector<float>> _dcaZ{"dcaZ", std::vector<float>{0.3f, 0.0f, 0.0f}, "abs dcaZ value limit; formula: [0] + [1]*pT^[2]"};
  Configurable<int16_t> _tpcNClsFound{"minTpcNClsFound", 0, "minimum allowed number of TPC clasters"};
  Configurable<float> _tpcChi2NCl{"tpcChi2NCl", 100.0, "upper limit for chi2 value of a fit over TPC clasters"};
  Configurable<float> _tpcCrossedRowsOverFindableCls{"tpcCrossedRowsOverFindableCls", 0, "lower limit of TPC CrossedRows/FindableCls value"};
  Configurable<float> _tpcFractionSharedCls{"maxTpcFractionSharedCls", 0.4, "maximum fraction of TPC shared clasters"};
  Configurable<int> _itsNCls{"minItsNCls", 0, "minimum allowed number of ITS clasters for a track"};
  Configurable<float> _itsChi2NCl{"itsChi2NCl", 100.0, "upper limit for chi2 value of a fit over ITS clasters for a track"};

  Configurable<std::vector<float>> _dcaXYCond{"dcaXYCond", std::vector<float>{0.3f, 0.0f, 0.0f}, "abs dcaXY value limit for conditioned multiplicity; formula: [0] + [1]*pT^[2]"};
  Configurable<std::vector<float>> _dcaZCond{"dcaZCond", std::vector<float>{0.3f, 0.0f, 0.0f}, "abs dcaZ value limit for conditioned multiplicity; formula: [0] + [1]*pT^[2]"};
  Configurable<int16_t> _tpcNClsFoundCond{"minTpcNClsFoundCond", 0, "minimum allowed number of TPC clasters for conditioned multiplicity"};
  Configurable<float> _tpcChi2NClCond{"tpcChi2NClCond", 100.0, "upper limit for chi2 value of a fit over TPC clasters for conditioned multiplicity"};
  Configurable<float> _tpcCrossedRowsOverFindableClsCond{"tpcCrossedRowsOverFindableClsCond", 0, "lower limit of TPC CrossedRows/FindableCls value for conditioned multiplicity"};
  Configurable<int> _itsNClsCond{"minItsNClsCond", 0, "minimum allowed number of ITS clasters for a track for conditioned multiplicity"};
  Configurable<float> _itsChi2NClCond{"itsChi2NClCond", 100.0, "upper limit for chi2 value of a fit over ITS clasters for a track for conditioned multiplicity"};

  Configurable<int> _particlePDG{"particlePDG", 2212, "PDG code of a particle to perform PID for (only pion, kaon, proton and deurton are supported now)"};
  Configurable<std::vector<float>> _tpcNSigma{"tpcNSigma", std::vector<float>{-4.0f, 4.0f}, "Nsigma range in TPC before the TOF is used"};
  Configurable<std::vector<float>> _itsNSigma{"itsNSigma", std::vector<float>{-10.0f, 10.0f}, "Nsigma range in ITS to use along with TPC"};
  Configurable<float> _PIDtrshld{"PIDtrshld", 10.0, "value of momentum from which the PID is done with TOF (before that only TPC is used)"};
  Configurable<std::vector<float>> _tofNSigma{"tofNSigma", std::vector<float>{-4.0f, 4.0f}, "Nsigma range in TOF"};
  Configurable<std::vector<float>> _tpcNSigmaResidual{"tpcNSigmaResidual", std::vector<float>{-5.0f, 5.0f}, "residual TPC Nsigma cut to use with the TOF"};

  Configurable<int> _particlePDGtoReject{"particlePDGtoReject", 211, "PDG codes of perticles that will be rejected with TOF (only pion, kaon, proton and deurton are supported now)"};
  Configurable<std::vector<float>> _rejectWithinNsigmaTOF{"rejectWithinNsigmaTOF", std::vector<float>{-0.0f, 0.0f}, "TOF rejection Nsigma range for particles specified with PDG to be rejected"};

  Configurable<int> _minNcandidates{"minNcandidates", 1, "minimum namber of candidates for the conditional multiplicity"};

  std::pair<int, std::vector<float>> TPCcuts;
  std::pair<int, std::vector<float>> TOFcuts;

  std::shared_ptr<TH2> ITShisto;
  std::shared_ptr<TH2> TPChisto;
  std::shared_ptr<TH2> TOFhisto;

  Configurable<std::pair<float, float>> _centCut{"centCut", std::pair<float, float>{0.f, 100.f}, "[min., max.] centrality range to keep events within"};

  using Trks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidEvTimeFlags, aod::TracksDCA,
                         aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa,
                         aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullTr, aod::pidTPCFullHe,
                         aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa,
                         aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe,
                         aod::TrackSelection, aod::pidTOFbeta>;

  using CollRun3 = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFDDMs, aod::CentNTPVs>;

  Filter eventFilter = (applyEvSel.node() == 0) ||
                       ((applyEvSel.node() == 1) && (aod::evsel::sel7 == true)) ||
                       ((applyEvSel.node() == 2) && (aod::evsel::sel8 == true));
  Filter vertexFilter = nabs(o2::aod::collision::posZ) < _vertexZ;

  Filter etaFilter = nabs(o2::aod::track::eta) < _eta;

  ctpRateFetcher mRateFetcher; // inspired by zdcSP.cxx in PWGLF
  int mRunNumber = 0;

  int tot_counter = 0;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};
  SliceCache cache;

  rctsel::RCTFlagsChecker myChecker{"CBT_hadronPID"};

  std::shared_ptr<TH2> Nch_vs_cent_vs_eta;
  std::shared_ptr<TH2> Nch_vs_cent_vs_eta_conditional;
  std::shared_ptr<TH2> tmp_histo_per_event;
  std::shared_ptr<TH1> Events_vs_cent;
  std::shared_ptr<TH1> Events_vs_cent_conditional;

  void init(InitContext&)
  {

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    myChecker.init("CBT_hadronPID", true);

    TPCcuts = std::make_pair(_particlePDG.value, _tpcNSigma.value);
    TOFcuts = std::make_pair(_particlePDG.value, _tofNSigma.value);

    Nch_vs_cent_vs_eta = registry.add<TH2>("Nch_vs_cent_vs_eta", "Nch_vs_cent_vs_eta", kTH2F, {{100, 0.0, 100.0, "cent"}, {200, -1.0, 1.0, "deta"}});
    Nch_vs_cent_vs_eta_conditional = registry.add<TH2>("Nch_vs_cent_vs_eta_conditional", "Nch_vs_cent_vs_eta_conditional", kTH2F, {{100, 0.0, 100.0, "cent"}, {200, -1.0, 1.0, "deta"}});
    tmp_histo_per_event = std::make_shared<TH2F>(TH2F("tmp_histo_per_event", "tmp_histo_per_event", 100, 0.0, 100.0, 200, -1.0, 1.0));

    Events_vs_cent = std::make_shared<TH1F>(TH1F("Events_vs_cent", "Events_vs_cent", 100, 0.0, 100.));
    Events_vs_cent_conditional = std::make_shared<TH1F>(TH1F("Events_vs_cent_conditional", "Events_vs_cent_conditional", 100, 0.0, 100.));

    ITShisto = registry.add<TH2>(Form("nsigmaITS_PDG%i", _particlePDG.value), Form("nsigmaITS_PDG%i", _particlePDG.value), kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
    TPChisto = registry.add<TH2>(Form("nsigmaTPC_PDG%i", _particlePDG.value), Form("nsigmaTPC_PDG%i", _particlePDG.value), kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
    TOFhisto = registry.add<TH2>(Form("nsigmaTOF_PDG%i", _particlePDG.value), Form("nsigmaTOF_PDG%i", _particlePDG.value), kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc) // inspired by PWGLF/TableProducer/lambdakzerobuilder.cxx
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
  }

  void process(soa::Filtered<CollRun3>::iterator const& collision,
               soa::Filtered<Trks> const& tracks,
               aod::BCsWithTimestamps const&)

  {

    const auto& bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    if (!myChecker(*collision) && CBThadronPID)
      return;

    if (_removeSameBunchPileup && !collision.selection_bit(evsel::kNoSameBunchPileup))
      return;
    if (_requestGoodZvtxFT0vsPV && !collision.selection_bit(evsel::kIsGoodZvtxFT0vsPV))
      return;
    if (std::fabs(collision.posZ()) > _vertexZ)
      return;

    int occupancy = collision.trackOccupancyInTimeRange();

    float centValue = -100.0f;

    switch (centTableToUse) {
      case 0:
        centValue = collision.centFV0A();
        break;
      case 1:
        centValue = collision.centFT0M();
        break;
      case 2:
        centValue = collision.centFT0A();
        break;
      case 3:
        centValue = collision.centFT0C();
        break;
      case 4:
        centValue = collision.centFDDM();
        break;
      case 5:
        centValue = collision.centNTPV();
        break;
      default:
        LOGF(fatal, "Invalid flag for cent./mult.perc. estimator has been choosen. Please check.");
        break;
    }

    if (centValue < _centCut.value.first || centValue >= _centCut.value.second)
      return;

    Events_vs_cent->Fill(centValue);
    // ============================ dNch/deta no cuts ================================

    int cadidates_counter = 0;

    // auto tracksWithPid = soa::Attach<soa::Filtered<Trks>, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaKa, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe, aod::pidits::ITSNSigmaTr, aod::pidits::ITSNSigmaHe>(tracks);
    auto tracksWithPid = soa::Attach<soa::Filtered<Trks>, aod::pidits::ITSNSigmaPr>(tracks);

    for (const auto& track : tracksWithPid) {

      // ============================ dNch/deta normal ================================

      if (rejectNotPropagatedTrks && track.trackType() != aod::track::Track) {
        continue;
      }

      if ((track.tpcNClsFound()) < _tpcNClsFound || (track.itsNCls()) < _itsNCls)
        continue;
      if ((track.itsChi2NCl()) > _itsChi2NCl || (track.tpcChi2NCl()) > _tpcChi2NCl)
        continue;
      if ((track.tpcCrossedRowsOverFindableCls()) < _tpcCrossedRowsOverFindableCls)
        continue;
      if (std::fabs(track.dcaXY()) > _dcaXY.value[0] + _dcaXY.value[1] * std::pow(track.pt(), _dcaXY.value[2]) || std::fabs(track.dcaZ()) > _dcaZ.value[0] + _dcaZ.value[1] * std::pow(track.pt(), _dcaZ.value[2]))
        continue;

      if (track.sign() != 0) {
        tmp_histo_per_event->Fill(centValue, track.eta());
        tot_counter++;
      }

      // ============================ dNch/deta conditioned (additional cuts) ================================

      // analisiys track cuts
      if (std::fabs(track.eta()) > _eta) // because soa::Attach operation breaks the filtering
        continue;
      if (track.sign() != _sign)
        continue;
      if (track.p() < _min_P || track.p() > _max_P)
        continue;
      if (_requestVertexITSTPC && !collision.selection_bit(evsel::kIsVertexITSTPC))
        continue;
      if (_requestVertexTOForTRDmatched > static_cast<int>(collision.selection_bit(evsel::kIsVertexTOFmatched)) + static_cast<int>(collision.selection_bit(evsel::kIsVertexTRDmatched)))
        continue;
      if (_requestNoCollInTimeRangeStandard && !collision.selection_bit(evsel::kNoCollInTimeRangeStandard))
        continue;
      if (_requestIsGoodITSLayersAll && !collision.selection_bit(evsel::kIsGoodITSLayersAll))
        continue;

      if (occupancy < _OccupancyCut.value.first || occupancy >= _OccupancyCut.value.second)
        continue;
      if (track.tofChi2() > _maxTofChi2)
        continue;
      if ((track.tpcNClsFound()) < _tpcNClsFoundCond || (track.itsNCls()) < _itsNClsCond)
        continue;
      if ((track.itsChi2NCl()) > _itsChi2NClCond || (track.tpcChi2NCl()) > _tpcChi2NClCond)
        continue;
      if ((track.tpcCrossedRowsOverFindableCls()) < _tpcCrossedRowsOverFindableClsCond)
        continue;
      if (std::fabs(track.dcaXY()) > _dcaXYCond.value[0] + _dcaXYCond.value[1] * std::pow(track.pt(), _dcaXYCond.value[2]) || std::fabs(track.dcaZ()) > _dcaZCond.value[0] + _dcaZCond.value[1] * std::pow(track.pt(), _dcaZCond.value[2]))
        continue;
      // ========================== PID cuts ==========================

      if (o2::aod::singletrackselector::TOFselection(track, std::make_pair(_particlePDGtoReject, _rejectWithinNsigmaTOF)))
        continue;

      if (track.p() < _PIDtrshld ? o2::aod::singletrackselector::TPCselection<true>(track, TPCcuts, _itsNSigma.value) : o2::aod::singletrackselector::TOFselection(track, TOFcuts, _tpcNSigmaResidual.value)) {
        ITShisto->Fill(track.p(), o2::aod::singletrackselector::getITSNsigma(track, _particlePDG));
        TPChisto->Fill(track.p(), o2::aod::singletrackselector::getTPCNsigma(track, _particlePDG));
        TOFhisto->Fill(track.p(), o2::aod::singletrackselector::getTOFNsigma(track, _particlePDG));
        cadidates_counter++;
      }
    }
    Nch_vs_cent_vs_eta->Add(tmp_histo_per_event.get());
    if (Nch_vs_cent_vs_eta->GetEntries() != tot_counter)
      LOGF(fatal, "tot counter != entries");
    if (cadidates_counter >= _minNcandidates) {
      Nch_vs_cent_vs_eta_conditional->Add(tmp_histo_per_event.get());
      Events_vs_cent_conditional->Fill(centValue);
    }

    tmp_histo_per_event->Reset();
    if (tmp_histo_per_event->GetEntries())
      LOGF(fatal, "not cleared");
    cadidates_counter = 0;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<multEstimator>(cfgc)};
}
