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
/// \brief Femto3D pair mixing task
/// \author Sofia Tomassini, Gleb Romanenko, Nicolò Jacazio
/// \since 31 May 2023

#include <ctime>
#include <algorithm> // std::random_shuffle
#include <random>
#include <chrono>
#include <vector>
#include <map>
#include <memory>
#include <utility>
#include <TParameter.h>
#include <TH1F.h>

#include "PWGCF/Femto3D/Core/femto3dPairTask.h"
#include "PWGCF/Femto3D/DataModel/singletrackselector.h"
#include "TLorentzVector.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoA.h"
#include "Framework/DataTypes.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/Expressions.h"
#include "Framework/StaticFor.h"
#include "MathUtils/Utils.h"
#include "Common/DataModel/Multiplicity.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct FemtoCorrelations {
  // using allinfo = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::pidTPCFullPr, aod::TOFSignal, aod::TracksDCA, aod::pidTOFFullPr, aod::pidTOFbeta, aod::pidTOFFullKa, aod::pidTPCFullKa, aod::pidTOFFullDe, aod::pidTPCFullDe>; // aod::pidTPCPr
  /// Construct a registry object with direct declaration
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<bool> _removeSameBunchPileup{"removeSameBunchPileup", false, ""};
  Configurable<bool> _requestGoodZvtxFT0vsPV{"requestGoodZvtxFT0vsPV", false, ""};
  Configurable<bool> _requestVertexITSTPC{"requestVertexITSTPC", false, ""};
  Configurable<int> _requestVertexTOForTRDmatched{"requestVertexTOFmatched", 0, "0 -> no selectio; 1 -> vertex is matched to TOF or TRD; 2 -> matched to both;"};
  Configurable<bool> _requestNoCollInTimeRangeStandard{"requestNoCollInTimeRangeStandard", false, ""};
  Configurable<std::pair<float, float>> _IRcut{"IRcut", std::pair<float, float>{0.f, 100.f}, "[min., max.] IR range to keep events within"};
  Configurable<std::pair<int, int>> _OccupancyCut{"OccupancyCut", std::pair<int, int>{0, 10000}, "[min., max.] occupancy range to keep events within"};

  Configurable<float> _min_P{"min_P", 0.0, "lower mometum limit"};
  Configurable<float> _max_P{"max_P", 100.0, "upper mometum limit"};
  Configurable<float> _eta{"eta", 100.0, "abs eta value limit"};
  Configurable<std::vector<float>> _dcaXY{"dcaXY", std::vector<float>{0.3f, 0.0f, 0.0f}, "abs dcaXY value limit; formula: [0] + [1]*pT^[2]"};
  Configurable<std::vector<float>> _dcaZ{"dcaZ", std::vector<float>{0.3f, 0.0f, 0.0f}, "abs dcaZ value limit; formula: [0] + [1]*pT^[2]"};
  Configurable<int16_t> _tpcNClsFound{"minTpcNClsFound", 0, "minimum allowed number of TPC clasters"};
  Configurable<float> _tpcChi2NCl{"tpcChi2NCl", 100.0, "upper limit for chi2 value of a fit over TPC clasters"};
  Configurable<float> _tpcCrossedRowsOverFindableCls{"tpcCrossedRowsOverFindableCls", 0, "lower limit of TPC CrossedRows/FindableCls value"};
  Configurable<float> _tpcFractionSharedCls{"maxTpcFractionSharedCls", 0.4, "maximum fraction of TPC shared clasters"};
  Configurable<int> _itsNCls{"minItsNCls", 0, "minimum allowed number of ITS clasters"};
  Configurable<float> _itsChi2NCl{"itsChi2NCl", 100.0, "upper limit for chi2 value of a fit over ITS clasters"};
  Configurable<float> _vertexZ{"VertexZ", 10.0, "abs vertexZ value limit"};

  Configurable<int> _sign_1{"sign_1", 1, "sign of the first particle in a pair"};
  Configurable<int> _particlePDG_1{"particlePDG_1", 2212, "PDG code of the first particle in a pair to perform PID for (only pion, kaon, proton and deurton are supported now)"};
  Configurable<std::vector<float>> _tpcNSigma_1{"tpcNSigma_1", std::vector<float>{-3.0f, 3.0f}, "first particle PID: Nsigma range in TPC before the TOF is used"};
  Configurable<std::vector<float>> _itsNSigma_1{"itsNSigma_1", std::vector<float>{-10.0f, 10.0f}, "first particle PID: Nsigma range in ITS with TPC is used"};
  Configurable<float> _PIDtrshld_1{"PIDtrshld_1", 10.0, "first particle PID: value of momentum from which the PID is done with TOF (before that only TPC is used)"};
  Configurable<std::vector<float>> _tofNSigma_1{"tofNSigma_1", std::vector<float>{-3.0f, 3.0f}, "first particle PID: Nsigma range in TOF"};
  Configurable<std::vector<float>> _tpcNSigmaResidual_1{"tpcNSigmaResidual_1", std::vector<float>{-5.0f, 5.0f}, "first particle PID: residual TPC Nsigma cut to use with the TOF"};

  Configurable<int> _sign_2{"sign_2", 1, "sign of the second particle in a pair"};
  Configurable<int> _particlePDG_2{"particlePDG_2", 2212, "PDG code of the second particle in a pair to perform PID for (only pion, kaon, proton and deurton are supported now)"};
  Configurable<std::vector<float>> _tpcNSigma_2{"tpcNSigma_2", std::vector<float>{-3.0f, 3.0f}, "second particle PID: Nsigma range in TPC before the TOF is used"};
  Configurable<std::vector<float>> _itsNSigma_2{"itsNSigma_2", std::vector<float>{-10.0f, 10.0f}, "first particle PID: Nsigma range in ITS with TPC is used"};
  Configurable<float> _PIDtrshld_2{"PIDtrshld_2", 10.0, "second particle PID: value of momentum from which the PID is done with TOF (before that only TPC is used)"};
  Configurable<std::vector<float>> _tofNSigma_2{"tofNSigma_2", std::vector<float>{-3.0f, 3.0f}, "second particle PID: Nsigma range in TOF"};
  Configurable<std::vector<float>> _tpcNSigmaResidual_2{"tpcNSigmaResidual_2", std::vector<float>{-5.0f, 5.0f}, "second particle PID: residual TPC Nsigma cut to use with the TOF"};

  Configurable<int> _particlePDGtoReject{"particlePDGtoRejectFromSecond", 0, "applied only if the particles are non-identical and only to the second particle in the pair!!!"};
  Configurable<std::vector<float>> _rejectWithinNsigmaTOF{"rejectWithinNsigmaTOF", std::vector<float>{-0.0f, 0.0f}, "TOF rejection Nsigma range for the particle specified with PDG to be rejected"};

  Configurable<float> _deta{"deta", 0.01, "minimum allowed defference in eta between two tracks in a pair"};
  Configurable<float> _dphi{"dphi", 0.01, "minimum allowed defference in phi_star between two tracks in a pair"};
  // Configurable<float> _radiusTPC{"radiusTPC", 1.2, "TPC radius to calculate phi_star for"};
  Configurable<float> _avgSepTPC{"avgSepTPC", 10, "average sep. (cm) in TPC"};

  Configurable<int> _vertexNbinsToMix{"vertexNbinsToMix", 10, "Number of vertexZ bins for the mixing"};
  Configurable<std::vector<float>> _centBins{"multBins", std::vector<float>{0.0f, 100.0f}, "multiplicity percentile/centrality binning (min:0, max:100)"};
  Configurable<int> _multNsubBins{"multSubBins", 1, "number of sub-bins to perform the mixing within"};
  Configurable<std::vector<float>> _kTbins{"kTbins", std::vector<float>{0.0f, 100.0f}, "pair transverse momentum kT binning"};
  ConfigurableAxis CFkStarBinning{"CFkStarBinning", {500, 0.005, 5.005}, "k* binning of the CF (Nbins, lowlimit, uplimit)"};

  Configurable<bool> _fill3dCF{"fill3dCF", false, "flag for filling 3D LCMS histos: true -- fill; false -- not"};
  Configurable<int> _fill3dAddHistos{"fill3dAddHistos", 1, "flag for filling additional 3D histos: 0 -- nothing; 1 -- Q_LCMS vs. k*; 2 -- Q_LCMS vs. Gamma_out (currently testing)"};
  Configurable<int> _fillDetaDphi{"fillDetaDphi", -1, "flag for filling dEta(dPhi*) histos: '-1' -- don't fill; '0' -- fill before the cut; '1' -- fill after the cut; '2' -- fill before & after the cut"};
  ConfigurableAxis CF3DqLCMSBinning{"CF3DqLCMSBinning", {60, -0.3, 0.3}, "q_out/side/long binning of the CF 3D in LCMS (Nbins, lowlimit, uplimit)"};
  // the next configarable is responsible for skipping (pseudo)randomly chosen ($value -1) pairs of events in the mixing process
  // migth be useful (for the sake of execution time and the output file size ...) in case of too many events per DF since in the SE thacks are mixed in N_ev and in the ME 0.5*N_ev*(N_ev - 1)
  // note that in the SE number of pairs per event = 0.5*N_trks*(N_trks - 1) BUT in the ME its = N_trks^2 so final differense in the SE & ME pairs histos will be more
  // P.S.: the explanation might be not very clear (sorry for that) as well as the name of the variable so, feel free to change it!
  // P.P.S. the chosen way of optimizing the mixing midgt not be the correct one -- feel free to propose the right one!
  // P.P.P.S. choose wisely....
  // P.P.P.P.S this way is still being testing i might be reconsidered; might change in the future, keep looking at the source code
  Configurable<unsigned int> _MEreductionFactor{"MEreductionFactor", 1, "only one (pseudo)randomly choosen event out per pair $value events will be processed and contribute to the final mixing (if < 2 -> all the possible event pairs (per vertex&cent bin) will be processed); implemented for the sake of efficiency; look at the source code;"};

  bool IsIdentical;

  std::pair<int, std::vector<float>> TPCcuts_1;
  std::pair<int, std::vector<float>> TOFcuts_1;

  std::pair<int, std::vector<float>> TPCcuts_2;
  std::pair<int, std::vector<float>> TOFcuts_2;

  using FilteredCollisions = soa::Join<aod::SingleCollSels, aod::SingleCollExtras>;
  using FilteredTracks = soa::Join<aod::SingleTrackSels, aod::SinglePIDPis, aod::SinglePIDKas, aod::SinglePIDPrs, aod::SinglePIDDes, aod::SinglePIDTrs, aod::SinglePIDHes>;

  typedef std::shared_ptr<soa::Filtered<FilteredTracks>::iterator> trkType;
  typedef std::shared_ptr<soa::Filtered<FilteredCollisions>::iterator> colType;

  std::map<int64_t, std::vector<trkType>> selectedtracks_1;
  std::map<int64_t, std::vector<trkType>> selectedtracks_2;
  std::map<std::pair<int, float>, std::vector<colType>> mixbins;

  std::unique_ptr<o2::aod::singletrackselector::FemtoPair<trkType>> Pair = std::make_unique<o2::aod::singletrackselector::FemtoPair<trkType>>();

  Filter pFilter = o2::aod::singletrackselector::p > _min_P&& o2::aod::singletrackselector::p < _max_P;
  Filter etaFilter = nabs(o2::aod::singletrackselector::eta) < _eta;

  Filter tpcTrkFilter = o2::aod::singletrackselector::tpcNClsFound >= _tpcNClsFound &&
                        o2::aod::singletrackselector::unPack<singletrackselector::binning::chi2>(o2::aod::singletrackselector::storedTpcChi2NCl) < _tpcChi2NCl &&
                        o2::aod::singletrackselector::unPack<singletrackselector::binning::rowsOverFindable>(o2::aod::singletrackselector::storedTpcCrossedRowsOverFindableCls) > _tpcCrossedRowsOverFindableCls;

  Filter itsTrkFilter = o2::aod::singletrackselector::unPack<singletrackselector::binning::chi2>(o2::aod::singletrackselector::storedItsChi2NCl) < _itsChi2NCl;

  Filter vertexFilter = nabs(o2::aod::singletrackselector::posZ) < _vertexZ;

  std::vector<std::shared_ptr<TH1>> MultHistos;
  std::vector<std::vector<std::shared_ptr<TH1>>> kThistos;
  std::vector<std::vector<std::shared_ptr<TH1>>> mThistos; // test
  std::vector<std::vector<std::shared_ptr<TH1>>> SEhistos_1D;
  std::vector<std::vector<std::shared_ptr<TH1>>> MEhistos_1D;

  std::vector<std::vector<std::shared_ptr<TH3>>> SEhistos_3D;
  std::vector<std::vector<std::shared_ptr<TH3>>> MEhistos_3D;
  std::vector<std::vector<std::shared_ptr<TH3>>> Add3dHistos;

  std::vector<std::vector<std::shared_ptr<TH2>>> DoubleTrack_SE_histos_BC; // BC -- before cutting
  std::vector<std::vector<std::shared_ptr<TH2>>> DoubleTrack_ME_histos_BC; // BC -- before cutting

  std::vector<std::vector<std::shared_ptr<TH2>>> DoubleTrack_SE_histos_AC; // AC -- after cutting
  std::vector<std::vector<std::shared_ptr<TH2>>> DoubleTrack_ME_histos_AC; // AC -- after cutting

  void init(o2::framework::InitContext&)
  {

    if (_centBins.value.size() < 2)
      LOGF(fatal, "The configured number of multiplicity/centrality bins in the array is less than 2 !!!");
    if (_kTbins.value.size() < 2)
      LOGF(fatal, "The configured number of kT bins in the array is less than 2 !!!");
    if (_vertexNbinsToMix.value < 1)
      LOGF(fatal, "The configured number of VertexZ bins is less than 1 !!!");

    IsIdentical = (_sign_1 * _particlePDG_1 == _sign_2 * _particlePDG_2);

    Pair->SetIdentical(IsIdentical);
    Pair->SetPDG1(_particlePDG_1);
    Pair->SetPDG2(_particlePDG_2);

    TPCcuts_1 = std::make_pair(_particlePDG_1, _tpcNSigma_1);
    TOFcuts_1 = std::make_pair(_particlePDG_1, _tofNSigma_1);
    TPCcuts_2 = std::make_pair(_particlePDG_2, _tpcNSigma_2);
    TOFcuts_2 = std::make_pair(_particlePDG_2, _tofNSigma_2);

    for (unsigned int i = 0; i < _centBins.value.size() - 1; i++) {
      std::vector<std::shared_ptr<TH1>> SEperMult_1D;
      std::vector<std::shared_ptr<TH1>> MEperMult_1D;
      std::vector<std::shared_ptr<TH1>> kTperMult;
      std::vector<std::shared_ptr<TH1>> mTperMult; // test

      auto hMult = registry.add<TH1>(Form("Cent%i/TPCMult_cent%i", i, i), Form("TPCMult_cent%i", i), kTH1F, {{5001, -0.5, 5000.5, "Mult."}});
      MultHistos.push_back(std::move(hMult));

      for (unsigned int j = 0; j < _kTbins.value.size() - 1; j++) {
        auto hSE_1D = registry.add<TH1>(Form("Cent%i/SE_1D_cent%i_kT%i", i, i, j), Form("SE_1D_cent%i_kT%i", i, j), kTH1F, {{CFkStarBinning, "k* (GeV/c)"}});
        auto hME_1D = registry.add<TH1>(Form("Cent%i/ME_1D_cent%i_kT%i", i, i, j), Form("ME_1D_cent%i_kT%i", i, j), kTH1F, {{CFkStarBinning, "k* (GeV/c)"}});
        auto hkT = registry.add<TH1>(Form("Cent%i/kT_cent%i_kT%i", i, i, j), Form("kT_cent%i_kT%i", i, j), kTH1F, {{500, 0., 5., "kT"}});
        auto hmT = registry.add<TH1>(Form("Cent%i/mT_test_cent%i_kT%i", i, i, j), Form("mT_cent%i_kT%i", i, j), kTH1F, {{500, 0., 5., "mT"}}); // test
        SEperMult_1D.push_back(std::move(hSE_1D));
        MEperMult_1D.push_back(std::move(hME_1D));
        kTperMult.push_back(std::move(hkT));
        mTperMult.push_back(std::move(hmT)); // test
      }

      SEhistos_1D.push_back(std::move(SEperMult_1D));
      MEhistos_1D.push_back(std::move(MEperMult_1D));
      kThistos.push_back(std::move(kTperMult));
      mThistos.push_back(std::move(mTperMult)); // test

      if (_fill3dCF) {
        std::vector<std::shared_ptr<TH3>> SEperMult_3D;
        std::vector<std::shared_ptr<TH3>> MEperMult_3D;
        std::vector<std::shared_ptr<TH3>> Add3dHistosperMult;

        for (unsigned int j = 0; j < _kTbins.value.size() - 1; j++) {
          auto hSE_3D = registry.add<TH3>(Form("Cent%i/SE_3D_cent%i_kT%i", i, i, j), Form("SE_3D_cent%i_kT%i", i, j), kTH3F, {{CF3DqLCMSBinning, "q_out (GeV/c)"}, {CF3DqLCMSBinning, "q_side (GeV/c)"}, {CF3DqLCMSBinning, "q_long (GeV/c)"}});
          auto hME_3D = registry.add<TH3>(Form("Cent%i/ME_3D_cent%i_kT%i", i, i, j), Form("ME_3D_cent%i_kT%i", i, j), kTH3F, {{CF3DqLCMSBinning, "q_out (GeV/c)"}, {CF3DqLCMSBinning, "q_side (GeV/c)"}, {CF3DqLCMSBinning, "q_long (GeV/c)"}});
          SEperMult_3D.push_back(std::move(hSE_3D));
          MEperMult_3D.push_back(std::move(hME_3D));

          if (_fill3dAddHistos == 1) {
            auto hAdd3dHistos = registry.add<TH3>(Form("Cent%i/qLCMSvskStar_cent%i_kT%i", i, i, j), Form("qLCMSvskStar_cent%i_kT%i", i, j), kTH3F, {{CF3DqLCMSBinning, "q_out (GeV/c)"}, {CF3DqLCMSBinning, "q_side (GeV/c)"}, {CF3DqLCMSBinning, "q_long (GeV/c)"}});
            Add3dHistosperMult.push_back(std::move(hAdd3dHistos));
          } else if (_fill3dAddHistos == 2) {
            auto hAdd3dHistos = registry.add<TH3>(Form("Cent%i/qLCMSvsGout_cent%i_kT%i", i, i, j), Form("qLCMSvsGout_cent%i_kT%i", i, j), kTH3F, {{CF3DqLCMSBinning, "q_out (GeV/c)"}, {CF3DqLCMSBinning, "q_side (GeV/c)"}, {CF3DqLCMSBinning, "q_long (GeV/c)"}});
            Add3dHistosperMult.push_back(std::move(hAdd3dHistos));
          }
        }
        SEhistos_3D.push_back(std::move(SEperMult_3D));
        MEhistos_3D.push_back(std::move(MEperMult_3D));
        if (_fill3dAddHistos != 0)
          Add3dHistos.push_back(std::move(Add3dHistosperMult));
      }

      if (_fillDetaDphi > -1) {
        TString suffix[] = {"_nocut", "_cut"};

        for (int f = 0; f < 2; f++) {
          std::vector<std::shared_ptr<TH2>> DoubleTrack_SE_histos_perMult;
          std::vector<std::shared_ptr<TH2>> DoubleTrack_ME_histos_perMult;

          if (_fillDetaDphi == 0 && f == 1)
            continue;
          if (_fillDetaDphi == 1 && f == 0)
            continue;

          for (unsigned int j = 0; j < _kTbins.value.size() - 1; j++) {
            auto hDblTrk_SE = registry.add<TH2>(Form("Cent%i/DoubleTrackEffects/SE_cent%i_kT%i", i, i, j) + suffix[f], Form("DoubleTrackEffects_deta(dphi*)_SE_cent%i_kT%i", i, j) + suffix[f], kTH2F, {{101, -0.2, 0.2, "dphi*"}, {101, -0.2, 0.2, "deta"}});
            auto hDblTrk_ME = registry.add<TH2>(Form("Cent%i/DoubleTrackEffects/ME_cent%i_kT%i", i, i, j) + suffix[f], Form("DoubleTrackEffects_deta(dphi*)_ME_cent%i_kT%i", i, j) + suffix[f], kTH2F, {{101, -0.2, 0.2, "dphi*"}, {101, -0.2, 0.2, "deta"}});

            DoubleTrack_SE_histos_perMult.push_back(std::move(hDblTrk_SE));
            DoubleTrack_ME_histos_perMult.push_back(std::move(hDblTrk_ME));
          }
          if (f == 0) {
            DoubleTrack_SE_histos_BC.push_back(std::move(DoubleTrack_SE_histos_perMult)); // BC -- before cutting
            DoubleTrack_ME_histos_BC.push_back(std::move(DoubleTrack_ME_histos_perMult));
          } else {
            DoubleTrack_SE_histos_AC.push_back(std::move(DoubleTrack_SE_histos_perMult)); // AC -- after cutting
            DoubleTrack_ME_histos_AC.push_back(std::move(DoubleTrack_ME_histos_perMult));
          }
        }
      }
    }

    registry.add("p_first", Form("p_%i", static_cast<int>(_particlePDG_1)), kTH1F, {{100, 0., 5., "p"}});
    registry.add("nsigmaTOF_first", Form("nsigmaTOF_%i", static_cast<int>(_particlePDG_1)), kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
    registry.add("nsigmaTPC_first", Form("nsigmaTPC_%i", static_cast<int>(_particlePDG_1)), kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
    if (!IsIdentical) {
      registry.add("p_second", Form("p_%i", static_cast<int>(_particlePDG_2)), kTH1F, {{100, 0., 5., "p"}});
      registry.add("nsigmaTOF_second", Form("nsigmaTOF_%i", static_cast<int>(_particlePDG_2)), kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
      registry.add("nsigmaTPC_second", Form("nsigmaTPC_%i", static_cast<int>(_particlePDG_2)), kTH2F, {{100, 0., 5.}, {100, -10., 10.}});
    }
  }

  template <typename Type>
  void mixTracks(Type const& tracks, unsigned int multBin)
  { // template for identical particles from the same collision
    if (multBin > SEhistos_1D.size())
      LOGF(fatal, "multBin value passed to the mixTracks function exceeds the configured number of Cent. bins (1D)");
    if (_fill3dCF && multBin > SEhistos_3D.size())
      LOGF(fatal, "multBin value passed to the mixTracks function exceeds the configured number of Cent. bins (3D)");

    for (unsigned int ii = 0; ii < tracks.size(); ii++) { // nested loop for all the combinations
      for (unsigned int iii = ii + 1; iii < tracks.size(); iii++) {

        Pair->SetPair(tracks[ii], tracks[iii]);
        float pair_kT = Pair->GetKt();

        if (pair_kT < *_kTbins.value.begin() || pair_kT >= *(_kTbins.value.end() - 1))
          continue;

        unsigned int kTbin = o2::aod::singletrackselector::getBinIndex<unsigned int>(pair_kT, _kTbins);
        if (kTbin > SEhistos_1D[multBin].size())
          LOGF(fatal, "kTbin value obtained for a pair exceeds the configured number of kT bins (1D)");
        if (_fill3dCF && kTbin > SEhistos_3D[multBin].size())
          LOGF(fatal, "kTbin value obtained for a pair exceeds the configured number of kT bins (3D)");

        if (_fillDetaDphi % 2 == 0)
          DoubleTrack_SE_histos_BC[multBin][kTbin]->Fill(Pair->GetAvgPhiStarDiff(), Pair->GetEtaDiff());

        if (_deta > 0 && _dphi > 0 && Pair->IsClosePair(_deta, _dphi))
          continue;
        if (_avgSepTPC > 0 && Pair->IsClosePair(_avgSepTPC))
          continue;

        if (_fillDetaDphi > 0)
          DoubleTrack_SE_histos_AC[multBin][kTbin]->Fill(Pair->GetAvgPhiStarDiff(), Pair->GetEtaDiff());

        kThistos[multBin][kTbin]->Fill(pair_kT);
        mThistos[multBin][kTbin]->Fill(Pair->GetMt());       // test
        SEhistos_1D[multBin][kTbin]->Fill(Pair->GetKstar()); // close pair rejection and fillig the SE histo

        if (_fill3dCF) {
          std::mt19937 mt(std::chrono::steady_clock::now().time_since_epoch().count());
          TVector3 qLCMS = std::pow(-1, (mt() % 2)) * Pair->GetQLCMS(); // introducing randomness to the pair order ([first, second]); important only for 3D because if there are any sudden order/correlation in the tables, it could couse unwanted asymmetries in the final 3d rel. momentum distributions; irrelevant in 1D case because the absolute value of the rel.momentum is taken
          SEhistos_3D[multBin][kTbin]->Fill(qLCMS.X(), qLCMS.Y(), qLCMS.Z());
        }
        Pair->ResetPair();
      }
    }
  }

  template <int SE_or_ME, typename Type>
  void mixTracks(Type const& tracks1, Type const& tracks2, unsigned int multBin)
  { // last value: 0 -- SE; 1 -- ME
    if (multBin > SEhistos_1D.size())
      LOGF(fatal, "multBin value passed to the mixTracks function exceeds the configured number of Cent. bins (1D)");
    if (_fill3dCF && multBin > SEhistos_3D.size())
      LOGF(fatal, "multBin value passed to the mixTracks function exceeds the configured number of Cent. bins (3D)");

    for (auto ii : tracks1) {
      for (auto iii : tracks2) {

        Pair->SetPair(ii, iii);
        float pair_kT = Pair->GetKt();

        if (pair_kT < *_kTbins.value.begin() || pair_kT >= *(_kTbins.value.end() - 1))
          continue;

        unsigned int kTbin = o2::aod::singletrackselector::getBinIndex<unsigned int>(pair_kT, _kTbins);
        if (kTbin > SEhistos_1D[multBin].size())
          LOGF(fatal, "kTbin value obtained for a pair exceeds the configured number of kT bins (1D)");
        if (_fill3dCF && kTbin > SEhistos_3D[multBin].size())
          LOGF(fatal, "kTbin value obtained for a pair exceeds the configured number of kT bins (3D)");

        if (_fillDetaDphi % 2 == 0) {
          if (!SE_or_ME)
            DoubleTrack_SE_histos_BC[multBin][kTbin]->Fill(Pair->GetAvgPhiStarDiff(), Pair->GetEtaDiff());
          else
            DoubleTrack_ME_histos_BC[multBin][kTbin]->Fill(Pair->GetAvgPhiStarDiff(), Pair->GetEtaDiff());
        }

        if (_deta > 0 && _dphi > 0 && Pair->IsClosePair(_deta, _dphi))
          continue;
        if (_avgSepTPC > 0 && Pair->IsClosePair(_avgSepTPC))
          continue;

        if (_fillDetaDphi > 0) {
          if (!SE_or_ME)
            DoubleTrack_SE_histos_AC[multBin][kTbin]->Fill(Pair->GetAvgPhiStarDiff(), Pair->GetEtaDiff());
          else
            DoubleTrack_ME_histos_AC[multBin][kTbin]->Fill(Pair->GetAvgPhiStarDiff(), Pair->GetEtaDiff());
        }

        if (!SE_or_ME) {
          SEhistos_1D[multBin][kTbin]->Fill(Pair->GetKstar());
          kThistos[multBin][kTbin]->Fill(pair_kT);
          mThistos[multBin][kTbin]->Fill(Pair->GetMt()); // test

          if (_fill3dCF) {
            std::mt19937 mt(std::chrono::steady_clock::now().time_since_epoch().count());
            TVector3 qLCMS = std::pow(-1, (mt() % 2)) * Pair->GetQLCMS(); // introducing randomness to the pair order ([first, second]); important only for 3D because if there are any sudden order/correlation in the tables, it could couse unwanted asymmetries in the final 3d rel. momentum distributions; irrelevant in 1D case because the absolute value of the rel.momentum is taken
            SEhistos_3D[multBin][kTbin]->Fill(qLCMS.X(), qLCMS.Y(), qLCMS.Z());
          }
        } else {
          MEhistos_1D[multBin][kTbin]->Fill(Pair->GetKstar());

          if (_fill3dCF) {
            std::mt19937 mt(std::chrono::steady_clock::now().time_since_epoch().count());
            TVector3 qLCMS = std::pow(-1, (mt() % 2)) * Pair->GetQLCMS(); // introducing randomness to the pair order ([first, second]); important only for 3D because if there are any sudden order/correlation in the tables, it could couse unwanted asymmetries in the final 3d rel. momentum distributions; irrelevant in 1D case because the absolute value of the rel.momentum is taken
            MEhistos_3D[multBin][kTbin]->Fill(qLCMS.X(), qLCMS.Y(), qLCMS.Z());
            if (_fill3dAddHistos == 1)
              Add3dHistos[multBin][kTbin]->Fill(qLCMS.X(), qLCMS.Y(), qLCMS.Z(), Pair->GetKstar());
            else if (_fill3dAddHistos == 2)
              Add3dHistos[multBin][kTbin]->Fill(qLCMS.X(), qLCMS.Y(), qLCMS.Z(), Pair->GetGammaOut());
          }
        }
        Pair->ResetPair();
      }
    }
  }

  void process(soa::Filtered<FilteredCollisions> const& collisions, soa::Filtered<FilteredTracks> const& tracks)
  {
    if (_particlePDG_1 == 0 || _particlePDG_2 == 0)
      LOGF(fatal, "One of passed PDG is 0!!!");

    for (auto track : tracks) {
      if (std::fabs(track.template singleCollSel_as<soa::Filtered<FilteredCollisions>>().posZ()) > _vertexZ)
        continue;
      if (_removeSameBunchPileup && !track.template singleCollSel_as<soa::Filtered<FilteredCollisions>>().isNoSameBunchPileup())
        continue;
      if (_requestGoodZvtxFT0vsPV && !track.template singleCollSel_as<soa::Filtered<FilteredCollisions>>().isGoodZvtxFT0vsPV())
        continue;
      if (_requestVertexITSTPC && !track.template singleCollSel_as<soa::Filtered<FilteredCollisions>>().isVertexITSTPC())
        continue;
      if (_requestVertexTOForTRDmatched > track.template singleCollSel_as<soa::Filtered<FilteredCollisions>>().isVertexTOForTRDmatched())
        continue;
      if (_requestNoCollInTimeRangeStandard && !track.template singleCollSel_as<soa::Filtered<FilteredCollisions>>().noCollInTimeRangeStandard())
        continue;
      if (track.tpcFractionSharedCls() > _tpcFractionSharedCls || track.itsNCls() < _itsNCls)
        continue;
      if (track.template singleCollSel_as<soa::Filtered<FilteredCollisions>>().multPerc() < *_centBins.value.begin() || track.template singleCollSel_as<soa::Filtered<FilteredCollisions>>().multPerc() >= *(_centBins.value.end() - 1))
        continue;
      if (track.template singleCollSel_as<soa::Filtered<FilteredCollisions>>().hadronicRate() < _IRcut.value.first || track.template singleCollSel_as<soa::Filtered<FilteredCollisions>>().hadronicRate() >= _IRcut.value.second)
        continue;
      if (track.template singleCollSel_as<soa::Filtered<FilteredCollisions>>().occupancy() < _OccupancyCut.value.first || track.template singleCollSel_as<soa::Filtered<FilteredCollisions>>().occupancy() >= _OccupancyCut.value.second)
        continue;
      if (std::fabs(track.dcaXY()) > _dcaXY.value[0] + _dcaXY.value[1] * std::pow(track.pt(), _dcaXY.value[2]) || std::fabs(track.dcaZ()) > _dcaZ.value[0] + _dcaZ.value[1] * std::pow(track.pt(), _dcaZ.value[2]))
        continue;

      if (track.sign() == _sign_1 && (track.p() < _PIDtrshld_1 ? o2::aod::singletrackselector::TPCselection<true>(track, TPCcuts_1, _itsNSigma_1.value) : o2::aod::singletrackselector::TOFselection(track, TOFcuts_1, _tpcNSigmaResidual_1.value))) { // filling the map: eventID <-> selected particles1
        selectedtracks_1[track.singleCollSelId()].push_back(std::make_shared<decltype(track)>(track));

        registry.fill(HIST("p_first"), track.p());
        if (_particlePDG_1 == 211) {
          registry.fill(HIST("nsigmaTOF_first"), track.p(), track.tofNSigmaPi());
          registry.fill(HIST("nsigmaTPC_first"), track.p(), track.tpcNSigmaPi());
        }
        if (_particlePDG_1 == 321) {
          registry.fill(HIST("nsigmaTOF_first"), track.p(), track.tofNSigmaKa());
          registry.fill(HIST("nsigmaTPC_first"), track.p(), track.tpcNSigmaKa());
        }
        if (_particlePDG_1 == 2212) {
          registry.fill(HIST("nsigmaTOF_first"), track.p(), track.tofNSigmaPr());
          registry.fill(HIST("nsigmaTPC_first"), track.p(), track.tpcNSigmaPr());
        }
        if (_particlePDG_1 == 1000010020) {
          registry.fill(HIST("nsigmaTOF_first"), track.p(), track.tofNSigmaDe());
          registry.fill(HIST("nsigmaTPC_first"), track.p(), track.tpcNSigmaDe());
        }
      }

      if (IsIdentical) {
        continue;
      } else if (track.sign() != _sign_2 && !TOFselection(track, std::make_pair(_particlePDGtoReject, _rejectWithinNsigmaTOF)) && (track.p() < _PIDtrshld_2 ? o2::aod::singletrackselector::TPCselection<true>(track, TPCcuts_2, _itsNSigma_2.value) : o2::aod::singletrackselector::TOFselection(track, TOFcuts_2, _tpcNSigmaResidual_2.value))) { // filling the map: eventID <-> selected particles2 if (see condition above ^)
        selectedtracks_2[track.singleCollSelId()].push_back(std::make_shared<decltype(track)>(track));

        registry.fill(HIST("p_second"), track.p());
        if (_particlePDG_2 == 211) {
          registry.fill(HIST("nsigmaTOF_second"), track.p(), track.tofNSigmaPi());
          registry.fill(HIST("nsigmaTPC_second"), track.p(), track.tpcNSigmaPi());
        }
        if (_particlePDG_2 == 321) {
          registry.fill(HIST("nsigmaTOF_second"), track.p(), track.tofNSigmaKa());
          registry.fill(HIST("nsigmaTPC_second"), track.p(), track.tpcNSigmaKa());
        }
        if (_particlePDG_2 == 2212) {
          registry.fill(HIST("nsigmaTOF_second"), track.p(), track.tofNSigmaPr());
          registry.fill(HIST("nsigmaTPC_second"), track.p(), track.tpcNSigmaPr());
        }
        if (_particlePDG_2 == 1000010020) {
          registry.fill(HIST("nsigmaTOF_second"), track.p(), track.tofNSigmaDe());
          registry.fill(HIST("nsigmaTPC_second"), track.p(), track.tpcNSigmaDe());
        }
      }
    }

    for (auto collision : collisions) {
      if (collision.multPerc() < *_centBins.value.begin() || collision.multPerc() >= *(_centBins.value.end() - 1))
        continue;
      if (collision.hadronicRate() < _IRcut.value.first || collision.hadronicRate() >= _IRcut.value.second)
        continue;
      if (collision.occupancy() < _OccupancyCut.value.first || collision.occupancy() >= _OccupancyCut.value.second)
        continue;

      if (_removeSameBunchPileup && !collision.isNoSameBunchPileup())
        continue;
      if (_requestGoodZvtxFT0vsPV && !collision.isGoodZvtxFT0vsPV())
        continue;
      if (_requestVertexITSTPC && !collision.isVertexITSTPC())
        continue;
      if (_requestVertexTOForTRDmatched > collision.isVertexTOForTRDmatched())
        continue;
      if (_requestNoCollInTimeRangeStandard && !collision.noCollInTimeRangeStandard())
        continue;

      if (selectedtracks_1.find(collision.globalIndex()) == selectedtracks_1.end()) {
        if (IsIdentical)
          continue;
        else if (selectedtracks_2.find(collision.globalIndex()) == selectedtracks_2.end())
          continue;
      }
      int vertexBinToMix = std::floor((collision.posZ() + _vertexZ) / (2 * _vertexZ / _vertexNbinsToMix));
      float centBinToMix = o2::aod::singletrackselector::getBinIndex<float>(collision.multPerc(), _centBins, _multNsubBins);

      mixbins[std::pair<int, float>{vertexBinToMix, centBinToMix}].push_back(std::make_shared<decltype(collision)>(collision));
    }

    //====================================== mixing starts here ======================================

    if (IsIdentical) { //====================================== mixing identical ======================================

      for (auto i = mixbins.begin(); i != mixbins.end(); i++) { // iterating over all vertex&mult bins
        unsigned int EvPerBin = (i->second).size();

        for (unsigned int indx1 = 0; indx1 < EvPerBin; indx1++) { // loop over all the events in each vertex&mult bin

          auto col1 = (i->second)[indx1];

          Pair->SetMagField1(col1->magField());
          Pair->SetMagField2(col1->magField());

          unsigned int centBin = std::floor((i->first).second);
          MultHistos[centBin]->Fill(col1->mult());

          mixTracks(selectedtracks_1[col1->index()], centBin); // mixing SE identical

          for (unsigned int indx2 = indx1 + 1; indx2 < EvPerBin; indx2++) { // nested loop for all the combinations of collisions in a chosen mult/vertex bin
            if (_MEreductionFactor.value > 1) {
              std::mt19937 mt(std::chrono::steady_clock::now().time_since_epoch().count());
              if ((mt() % (_MEreductionFactor.value + 1)) < _MEreductionFactor.value)
                continue;
            }

            auto col2 = (i->second)[indx2];

            Pair->SetMagField2(col2->magField());
            mixTracks<1>(selectedtracks_1[col1->index()], selectedtracks_1[col2->index()], centBin); // mixing ME identical, in <> brackets: 0 -- SE; 1 -- ME
          }
        }
      }

    } else { //====================================== mixing non-identical ======================================

      for (auto i = mixbins.begin(); i != mixbins.end(); i++) { // iterating over all vertex&mult bins
        unsigned int EvPerBin = (i->second).size();

        for (unsigned int indx1 = 0; indx1 < EvPerBin; indx1++) { // loop over all the events in each vertex&mult bin

          auto col1 = (i->second)[indx1];

          Pair->SetMagField1(col1->magField());
          Pair->SetMagField2(col1->magField());

          unsigned int centBin = std::floor((i->first).second);
          MultHistos[centBin]->Fill(col1->mult());

          mixTracks<0>(selectedtracks_1[col1->index()], selectedtracks_2[col1->index()], centBin); // mixing SE non-identical, in <> brackets: 0 -- SE; 1 -- ME

          for (unsigned int indx2 = indx1 + 1; indx2 < EvPerBin; indx2++) { // nested loop for all the combinations of collisions in a chosen mult/vertex bin
            if (_MEreductionFactor.value > 1) {
              std::mt19937 mt(std::chrono::steady_clock::now().time_since_epoch().count());
              if (mt() % (_MEreductionFactor.value + 1) < _MEreductionFactor.value)
                continue;
            }

            auto col2 = (i->second)[indx2];

            Pair->SetMagField2(col2->magField());
            mixTracks<1>(selectedtracks_1[col1->index()], selectedtracks_2[col2->index()], centBin); // mixing ME non-identical, in <> brackets: 0 -- SE; 1 -- ME
          }
        }
      }

    } //====================================== end of mixing non-identical ======================================

    // clearing up
    for (auto i = selectedtracks_1.begin(); i != selectedtracks_1.end(); i++)
      (i->second).clear();
    selectedtracks_1.clear();

    if (!IsIdentical) {
      for (auto i = selectedtracks_2.begin(); i != selectedtracks_2.end(); i++)
        (i->second).clear();
      selectedtracks_2.clear();
    }

    for (auto i = mixbins.begin(); i != mixbins.end(); i++)
      (i->second).clear();
    mixbins.clear();
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<FemtoCorrelations>(cfgc)};
}
