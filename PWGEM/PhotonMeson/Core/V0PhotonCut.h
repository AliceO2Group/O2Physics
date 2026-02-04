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

/// \file V0PhotonCut.h
/// \brief Header of class for V0 photon selection.
/// \author D. Sekihata, daiki.sekihata@cern.ch

#ifndef PWGEM_PHOTONMESON_CORE_V0PHOTONCUT_H_
#define PWGEM_PHOTONMESON_CORE_V0PHOTONCUT_H_

#include "PWGEM/PhotonMeson/Core/EMBitFlags.h"
#include "PWGEM/PhotonMeson/Core/EmMlResponsePCM.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCandidate.h"
#include "PWGEM/PhotonMeson/Utils/TrackSelection.h"

#include <CCDB/CcdbApi.h>
#include <Framework/ASoA.h>
#include <Framework/Array2D.h>

#include <TMath.h>
#include <TNamed.h>

#include <fairlogger/Logger.h>

#include <Rtypes.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace o2::pwgem::photonmeson;

namespace o2::analysis
{

// namespace per channel
namespace em_cuts_ml
{
// direction of the cut
enum CutDirection {
  CutGreater = 0, // require score < cut value
  CutSmaller,     // require score > cut value
  CutNot          // do not cut on score
};

static constexpr int NBins = 12;

static constexpr int NBinsPt = 12;
static constexpr int NCutScores = 2;
// default values for the pT bin edges, offset by 1 from the bin numbers in cuts array
constexpr double BinsPt[NBinsPt + 1] = {
  0.,
  0.25,
  0.5,
  0.75,
  1.,
  1.5,
  2.,
  4.,
  6.,
  10.,
  20.,
  50.,
  100.};
const auto vecBinsPt = std::vector<double>{BinsPt, BinsPt + NBinsPt + 1};
static constexpr int NBinsCent = 11;
constexpr double BinsCent[NBinsCent + 1] = {
  0.,
  5,
  10,
  20,
  30,
  40,
  50,
  60,
  70,
  80,
  90,
  100.};
const auto vecBinsCent = std::vector<double>{BinsCent, BinsCent + NBinsCent + 1};

// default values for the ML model paths, one model per pT bin
static const std::vector<std::string> modelPaths = {
  ""};

// default values for the cut directions
constexpr int CutDir[NCutScores] = {CutGreater, CutSmaller};
const auto vecCutDir = std::vector<int>{CutDir, CutDir + NCutScores};

// default values for the cuts
constexpr double Cuts[NBins][NCutScores] = {
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5},
  {0.5, 0.5}};

// row labels
static const std::vector<std::string> labelsPt = {
  "pT bin 0",
  "pT bin 1",
  "pT bin 2",
  "pT bin 3",
  "pT bin 4",
  "pT bin 5",
  "pT bin 6",
  "pT bin 7",
  "pT bin 8",
  "pT bin 9",
  "pT bin 10",
  "pT bin 11"};
// labels
static const std::vector<std::string> labelsCent = {
  "Cent bin 0",
  "Cent bin 1",
  "Cent bin 2",
  "Cent bin 3",
  "Cent bin 4",
  "Cent bin 5",
  "Cent bin 6",
  "Cent bin 7",
  "Cent bin 8",
  "Cent bin 9",
  "Cent bin 10"};

// column labels
static const std::vector<std::string> labelsCutScore = {"score primary photons", "score background"};
} // namespace em_cuts_ml

} // namespace o2::analysis

class V0PhotonCut : public TNamed
{
 public:
  V0PhotonCut() = default;
  V0PhotonCut(const char* name, const char* title) : TNamed(name, title) {}
  ~V0PhotonCut() override { delete mEmMlResponse; };

  enum class V0PhotonCuts : int {
    // v0 cut
    kMee = 0,
    kV0PtRange,
    kV0EtaRange,
    kAP,
    kPsiPair,
    kPhiV,
    kRxy,
    kCosPA,
    kPCA,
    kChi2KF,
    kRZLine,
    kOnWwireIB,
    kOnWwireOB,
    // leg cut
    kTrackPtRange,
    kTrackEtaRange,
    kTPCNCls,
    kTPCCrossedRows,
    kTPCCrossedRowsOverNCls,
    kTPCFracSharedClusters,
    kTPCChi2NDF,
    kTPCNsigmaEl,
    kTPCNsigmaPi,
    kDCAxy,
    kDCAz,
    kITSNCls,
    kITSChi2NDF,
    kITSClusterSize,
    kRequireITSTPC,
    kRequireITSonly,
    kRequireTPConly,
    kRequireTPCTRD,
    kRequireTPCTOF,
    kNCuts
  };

  /// \brief check if given v0 photon survives all cuts
  /// \param flags EMBitFlags where results will be stored
  /// \param v0s v0 photon table to check
  template <o2::soa::is_table TV0, typename TLeg>
  void AreSelectedRunning(EMBitFlags& flags, TV0 const& v0s) const
  {
    size_t iV0 = 0;
    for (const auto& v0 : v0s) {
      if (!IsSelected<decltype(v0), TLeg>(v0)) {
        flags.set(iV0);
      }
      ++iV0;
    }
  }

  template <o2::soa::is_iterator TV0, typename TLeg>
  bool IsSelected(TV0 const& v0) const
  {
    if (!IsSelectedV0(v0, V0PhotonCuts::kV0PtRange)) {
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kV0EtaRange)) {
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kMee)) {
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kAP)) {
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kPsiPair)) {
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kPhiV)) {
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kRxy)) {
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kCosPA)) {
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kPCA)) {
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kChi2KF)) {
      return false;
    }
    if (!IsSelectedV0(v0, V0PhotonCuts::kRZLine)) {
      return false;
    }

    if (mIsOnWwireIB && mIsOnWwireOB) {
      if (!IsSelectedV0(v0, V0PhotonCuts::kOnWwireIB) && !IsSelectedV0(v0, V0PhotonCuts::kOnWwireOB)) {
        return false;
      }
    } else if (mIsOnWwireIB) {
      if (!IsSelectedV0(v0, V0PhotonCuts::kOnWwireIB)) {
        return false;
      }
    } else if (mIsOnWwireOB) {
      if (!IsSelectedV0(v0, V0PhotonCuts::kOnWwireOB)) {
        return false;
      }
    }

    auto pos = v0.template posTrack_as<TLeg>();
    auto ele = v0.template negTrack_as<TLeg>();

    for (auto& track : {pos, ele}) {
      if (!IsSelectedTrack(track, V0PhotonCuts::kTrackPtRange)) {
        return false;
      }
      if (!IsSelectedTrack(track, V0PhotonCuts::kTrackEtaRange)) {
        return false;
      }
      if (!IsSelectedTrack(track, V0PhotonCuts::kDCAxy)) {
        return false;
      }
      if (!IsSelectedTrack(track, V0PhotonCuts::kDCAz)) {
        return false;
      }
      if (!track.hasITS() && !track.hasTPC()) { // track has to be ITSonly or TPConly or ITS-TPC
        return false;
      }
      if (mDisableITSonly && isITSonlyTrack(track)) {
        return false;
      }
      if (mDisableTPConly && isTPConlyTrack(track)) {
        return false;
      }

      if (mRejectITSib) {
        auto hits_ib = std::count_if(its_ib_Requirement.second.begin(), its_ib_Requirement.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
        auto hits_ob = std::count_if(its_ob_Requirement.second.begin(), its_ob_Requirement.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
        bool its_ob_only = (hits_ib <= its_ib_Requirement.first) && (hits_ob >= its_ob_Requirement.first);
        if (isITSonlyTrack(track) && !its_ob_only) { // ITSonly tracks should not have any ITSib hits.
          return false;
        }

        auto hits_ob_itstpc = std::count_if(its_ob_Requirement_ITSTPC.second.begin(), its_ob_Requirement_ITSTPC.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
        bool its_ob_only_itstpc = (hits_ib <= its_ib_Requirement.first) && (hits_ob_itstpc >= its_ob_Requirement_ITSTPC.first);
        if (isITSTPCTrack(track) && !its_ob_only_itstpc) { // ITSTPC tracks should not have any ITSib hits.
          return false;
        }
      }

      if (track.hasITS() && !CheckITSCuts(track)) {
        return false;
      }
      if (track.hasTPC() && !CheckTPCCuts(track)) {
        return false;
      }
      if (track.hasITS() && !track.hasTPC() && (track.hasTRD() || track.hasTOF())) { // remove ITS-TRD, ITS-TOF, ITS-TRD-TOF that are unrealistic tracks.
        return false;
      }

      if (mIsOnWwireIB && !CheckITSCuts(track)) { // photon conversion on ibw requires ITS hits.
        return false;
      }

      if (mRequireITSonly && !IsSelectedTrack(track, V0PhotonCuts::kRequireITSonly)) {
        return false;
      }
      if (mRequireITSTPC && !IsSelectedTrack(track, V0PhotonCuts::kRequireITSTPC)) {
        return false;
      }
      if (mRequireTPConly && !IsSelectedTrack(track, V0PhotonCuts::kRequireTPConly)) {
        return false;
      }
      if (mRequireTPCTRD && !IsSelectedTrack(track, V0PhotonCuts::kRequireTPCTRD)) {
        return false;
      }
      if (mRequireTPCTOF && !IsSelectedTrack(track, V0PhotonCuts::kRequireTPCTOF)) {
        return false;
      }
    }
    if (mApplyMlCuts) {
      if (!mEmMlResponse) {
        LOG(error) << "EM ML Response is not initialized!";
        return false;
      }
      bool mIsSelectedMl = false;
      std::vector<float> mOutputML;
      V0PhotonCandidate v0photoncandidate(v0, pos, ele, mCentFT0A, mCentFT0C, mCentFT0M, mD_Bz);
      std::vector<float> mlInputFeatures = mEmMlResponse->getInputFeatures(v0photoncandidate, pos, ele);
      if (mUse2DBinning) {
        if (mCentralityTypeMl == "CentFT0C") {
          mIsSelectedMl = mEmMlResponse->isSelectedMl(mlInputFeatures, v0photoncandidate.getPt(), v0photoncandidate.getCentFT0C(), mOutputML);
        } else if (mCentralityTypeMl == "CentFT0A") {
          mIsSelectedMl = mEmMlResponse->isSelectedMl(mlInputFeatures, v0photoncandidate.getPt(), v0photoncandidate.getCentFT0A(), mOutputML);
        } else if (mCentralityTypeMl == "CentFT0M") {
          mIsSelectedMl = mEmMlResponse->isSelectedMl(mlInputFeatures, v0photoncandidate.getPt(), v0photoncandidate.getCentFT0M(), mOutputML);
        } else {
          LOG(fatal) << "Unsupported centTypePCMMl: " << mCentralityTypeMl << " , please choose from CentFT0C, CentFT0A, CentFT0M.";
        }
      } else {
        mIsSelectedMl = mEmMlResponse->isSelectedMl(mlInputFeatures, v0photoncandidate.getPt(), mOutputML);
      }
      if (!mIsSelectedMl) {
        return false;
      }
    }
    return true;
  }

  template <typename T>
  bool CheckITSCuts(T const& track) const
  {
    if (!IsSelectedTrack(track, V0PhotonCuts::kITSNCls)) {
      return false;
    }
    if (!IsSelectedTrack(track, V0PhotonCuts::kITSChi2NDF)) {
      return false;
    }
    if (!IsSelectedTrack(track, V0PhotonCuts::kITSClusterSize)) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool CheckTPCCuts(T const& track) const
  {
    if (!IsSelectedTrack(track, V0PhotonCuts::kTPCNCls)) {
      return false;
    }
    if (!IsSelectedTrack(track, V0PhotonCuts::kTPCCrossedRows)) {
      return false;
    }
    if (!IsSelectedTrack(track, V0PhotonCuts::kTPCCrossedRowsOverNCls)) {
      return false;
    }
    if (!IsSelectedTrack(track, V0PhotonCuts::kTPCFracSharedClusters)) {
      return false;
    }
    if (!IsSelectedTrack(track, V0PhotonCuts::kTPCChi2NDF)) {
      return false;
    }
    if (!IsSelectedTrack(track, V0PhotonCuts::kTPCNsigmaEl)) {
      return false;
    }
    if (!IsSelectedTrack(track, V0PhotonCuts::kTPCNsigmaPi)) {
      return false;
    }
    return true;
  }

  template <o2::soa::is_iterator T>
  bool IsSelectedV0(T const& v0, const V0PhotonCuts& cut) const
  {
    switch (cut) {
      case V0PhotonCuts::kV0PtRange:
        return v0.pt() >= mMinV0Pt && v0.pt() <= mMaxV0Pt;

      case V0PhotonCuts::kV0EtaRange:
        return v0.eta() >= mMinV0Eta && v0.eta() <= mMaxV0Eta;

      case V0PhotonCuts::kMee:
        return v0.mGamma() < mMaxQt * 2.f;

      case V0PhotonCuts::kAP:
        return pow(v0.alpha() / mMaxAlpha, 2) + pow(v0.qtarm() / mMaxQt, 2) < 1.0;

      case V0PhotonCuts::kPsiPair:
        return true;

      case V0PhotonCuts::kPhiV:
        return true;

      case V0PhotonCuts::kRxy: {
        if (v0.v0radius() < mMinRxy || mMaxRxy < v0.v0radius()) {
          return false;
        }
        return true;
      }

      case V0PhotonCuts::kCosPA:
        return v0.cospa() >= mMinCosPA;

      case V0PhotonCuts::kPCA:
        return v0.pca() <= mMaxPCA;

      case V0PhotonCuts::kChi2KF:
        return v0.chiSquareNDF() <= mMaxChi2KF;

      case V0PhotonCuts::kRZLine:
        return v0.v0radius() > std::fabs(v0.vz()) * std::tan(2 * std::atan(std::exp(-mMaxV0Eta))) - mMaxMarginZ;

      case V0PhotonCuts::kOnWwireIB: {
        const float margin_xy = 1.0; // cm
        // const float margin_z = 20.0; // cm
        // const float rxy_min = 5.506; // cm
        // const float rxy_max = 14.846;         // cm
        // const float z_min = -17.56; // cm
        // const float z_max = +31.15;           // cm
        float x = std::fabs(v0.vx()); // cm, measured secondary vertex of gamma->ee
        float y = v0.vy();            // cm, measured secondary vertex of gamma->ee
        float z = v0.vz();            // cm, measured secondary vertex of gamma->ee

        float rxy = sqrt(x * x + y * y);
        if (rxy < 7.0 || 14.0 < rxy) {
          return false;
        }

        // r = 0.192 * z + 8.88 (cm) expected wire position in RZ plane.TMath::Tan(10.86 * TMath::DegToRad()) = 0.192
        if (rxy > 0.192 * z + 14.0) { // upper limit
          return false;
        }

        float dxy = std::fabs(1.0 * y - x * std::tan(-8.52 * TMath::DegToRad())) / sqrt(pow(1.0, 2) + pow(std::tan(-8.52 * TMath::DegToRad()), 2));
        return !(dxy > margin_xy);
      }
      case V0PhotonCuts::kOnWwireOB: {
        const float margin_xy = 1.0;                                      // cm
        const float rxy_exp = 30.8;                                       // cm
        const float x_exp = rxy_exp * std::cos(-1.3 * TMath::DegToRad()); // cm, expected position x of W wire
        const float y_exp = rxy_exp * std::sin(-1.3 * TMath::DegToRad()); // cm, expected position y of W wire
        // const float z_min = -47.0;                                          // cm
        // const float z_max = +47.0;                                          // cm
        float x = v0.vx(); // cm, measured secondary vertex of gamma->ee
        float y = v0.vy(); // cm, measured secondary vertex of gamma->ee
        // float z = v0.vz();                                                  // cm, measured secondary vertex of gamma->ee

        // float rxy = sqrt(x * x + y * y);
        // if (rxy < 28.0 || 33.0 < rxy) {
        //   return false;
        // }

        float dxy = std::sqrt(std::pow(x - x_exp, 2) + std::pow(y - y_exp, 2));
        return !(dxy > margin_xy);
      }
      default:
        return false;
    }
  }

  // Temporary function to check if track passes a given selection criteria. To be replaced by framework filters.
  template <typename T>
  bool IsSelectedTrack(T const& track, const V0PhotonCuts& cut) const
  {
    switch (cut) {
      case V0PhotonCuts::kTrackPtRange:
        return track.pt() > mMinTrackPt && track.pt() < mMaxTrackPt;

      case V0PhotonCuts::kTrackEtaRange:
        return track.eta() > mMinTrackEta && track.eta() < mMaxTrackEta;

      case V0PhotonCuts::kTPCNCls:
        return track.tpcNClsFound() >= mMinNClustersTPC;

      case V0PhotonCuts::kTPCCrossedRows:
        return track.tpcNClsCrossedRows() >= mMinNCrossedRowsTPC;

      case V0PhotonCuts::kTPCCrossedRowsOverNCls:
        return track.tpcCrossedRowsOverFindableCls() >= mMinNCrossedRowsOverFindableClustersTPC;

      case V0PhotonCuts::kTPCFracSharedClusters:
        return track.tpcFractionSharedCls() < mMaxFracSharedClustersTPC;

      case V0PhotonCuts::kTPCChi2NDF:
        return mMinChi2PerClusterTPC < track.tpcChi2NCl() && track.tpcChi2NCl() < mMaxChi2PerClusterTPC;

      case V0PhotonCuts::kTPCNsigmaEl:
        return track.tpcNSigmaEl() > mMinTPCNsigmaEl && track.tpcNSigmaEl() < mMaxTPCNsigmaEl;

      case V0PhotonCuts::kTPCNsigmaPi:
        return track.tpcNSigmaPi() > mMinTPCNsigmaPi && track.tpcNSigmaPi() < mMaxTPCNsigmaPi;

      case V0PhotonCuts::kDCAxy:
        return std::fabs(track.dcaXY()) < ((mMaxDcaXYPtDep) ? mMaxDcaXYPtDep(track.pt()) : mMaxDcaXY);

      case V0PhotonCuts::kDCAz:
        return std::fabs(track.dcaZ()) < mMaxDcaZ;

      case V0PhotonCuts::kITSNCls:
        return mMinNClustersITS <= track.itsNCls() && track.itsNCls() <= mMaxNClustersITS;

      case V0PhotonCuts::kITSChi2NDF:
        return mMinChi2PerClusterITS < track.itsChi2NCl() && track.itsChi2NCl() < mMaxChi2PerClusterITS;

      case V0PhotonCuts::kITSClusterSize: {
        if (!isITSonlyTrack(track)) {
          return true;
        }
        return mMinMeanClusterSizeITS < track.meanClusterSizeITSob() * std::cos(std::atan(track.tgl())) && track.meanClusterSizeITSob() * std::cos(std::atan(track.tgl())) < mMaxMeanClusterSizeITS;
      }

      case V0PhotonCuts::kRequireITSTPC:
        return isITSTPCTrack(track);

      case V0PhotonCuts::kRequireITSonly:
        return isITSonlyTrack(track);

      case V0PhotonCuts::kRequireTPConly:
        return isTPConlyTrack(track);

      case V0PhotonCuts::kRequireTPCTRD:
        return isTPCTRDTrack(track);

      case V0PhotonCuts::kRequireTPCTOF:
        return isTPCTOFTrack(track);

      default:
        return false;
    }
  }

  void initV0MlModels(o2::ccdb::CcdbApi& ccdbApi)
  {
    if (!mEmMlResponse) {
      mEmMlResponse = new o2::analysis::EmMlResponsePCM<float>();
    }
    if (mUse2DBinning) {
      int binsNPt = static_cast<int>(mBinsPtMl.size()) - 1;
      int binsNCent = static_cast<int>(mBinsCentMl.size()) - 1;
      int binsN = binsNPt * binsNCent;
      if (binsN * static_cast<int>(mCutDirMl.size()) != static_cast<int>(mCutsMlFlat.size())) {
        LOG(fatal) << "Mismatch in number of bins and cuts provided for 2D ML application: binsN * mCutDirMl: " << int(binsN) * int(mCutDirMl.size()) << " bins vs. mCutsMlFlat: " << mCutsMlFlat.size() << " cuts";
      }
      if (binsN != static_cast<int>(mOnnxFileNames.size())) {
        LOG(fatal) << "Mismatch in number of bins and ONNX files provided for 2D ML application: binsN " << binsN << " bins vs. mOnnxFileNames: " << mOnnxFileNames.size() << " ONNX files";
      }
      if (binsN != static_cast<int>(mLabelsBinsMl.size())) {
        LOG(fatal) << "Mismatch in number of bins and labels provided for 2D ML application: binsN:" << binsN << " bins vs. mLabelsBinsMl: " << mLabelsBinsMl.size() << " labels";
      }
      if (static_cast<int>(mCutDirMl.size()) != mNClassesMl) {
        LOG(fatal) << "Mismatch in number of classes and cut directions provided for 2D ML application: mNClassesMl: " << mNClassesMl << " classes vs. mCutDirMl: " << mCutDirMl.size() << " cut directions";
      }
      if (static_cast<int>(mLabelsCutScoresMl.size()) != mNClassesMl) {
        LOG(fatal) << "Mismatch in number of labels for cut scores and number of classes provided for 2D ML application: mNClassesMl: " << mNClassesMl << " classes vs. mLabelsCutScoresMl: " << mLabelsCutScoresMl.size() << " labels";
      }
      o2::framework::LabeledArray<double> mCutsMl(mCutsMlFlat.data(), binsN, mNClassesMl, mLabelsBinsMl, mLabelsCutScoresMl);
      mEmMlResponse->configure2D(mBinsPtMl, mBinsCentMl, mCutsMl, mCutDirMl, mNClassesMl);
    } else {
      int binsNPt = static_cast<int>(mBinsPtMl.size()) - 1;
      if (binsNPt * static_cast<int>(mCutDirMl.size()) != static_cast<int>(mCutsMlFlat.size())) {
        LOG(fatal) << "Mismatch in number of pT bins and cuts provided for ML application: binsNPt * mCutDirMl:" << binsNPt * mCutDirMl.size() << " bins vs. mCutsMlFlat: " << mCutsMlFlat.size() << " cuts";
      }
      if (binsNPt != static_cast<int>(mOnnxFileNames.size())) {
        LOG(fatal) << "Mismatch in number of pT bins and ONNX files provided for ML application: binsNPt " << binsNPt << " bins vs. mOnnxFileNames: " << mOnnxFileNames.size() << " ONNX files";
      }
      if (binsNPt != static_cast<int>(mLabelsBinsMl.size())) {
        LOG(fatal) << "Mismatch in number of pT bins and labels provided for ML application: binsNPt:" << binsNPt << " bins vs. mLabelsBinsMl: " << mLabelsBinsMl.size() << " labels";
      }
      if (mNClassesMl != static_cast<int>(mCutDirMl.size())) {
        LOG(fatal) << "Mismatch in number of classes and cut directions provided for ML application: mNClassesMl: " << mNClassesMl << " classes vs. mCutDirMl: " << mCutDirMl.size() << " cut directions";
      }
      if (static_cast<int>(mLabelsCutScoresMl.size()) != mNClassesMl) {
        LOG(fatal) << "Mismatch in number of labels for cut scores and number of classes provided for ML application: mNClassesMl:" << mNClassesMl << " classes vs. mLabelsCutScoresMl: " << mLabelsCutScoresMl.size() << " labels";
      }
      o2::framework::LabeledArray<double> mCutsMl(mCutsMlFlat.data(), binsNPt, mNClassesMl, mLabelsBinsMl, mLabelsCutScoresMl);
      mEmMlResponse->configure(mBinsPtMl, mCutsMl, mCutDirMl, mNClassesMl);
    }
    if (mLoadMlModelsFromCCDB) {
      ccdbApi.init(mCcdbUrl);
      mEmMlResponse->setModelPathsCCDB(mOnnxFileNames, ccdbApi, mModelPathsCCDB, mTimestampCCDB);
    } else {
      mEmMlResponse->setModelPathsLocal(mOnnxFileNames);
    }
    mEmMlResponse->cacheInputFeaturesIndices(mNamesInputFeatures);
    mEmMlResponse->init();
  }

  // Setters
  void SetV0PtRange(float minPt = 0.f, float maxPt = 1e10f);
  void SetV0EtaRange(float minEta = -1e10f, float maxEta = 1e10f);
  void SetMeeRange(float min = 0.f, float max = 0.1);
  void SetPsiPairRange(float min = -3.15, float max = +3.15);
  void SetPhivPairRange(float min = 0.f, float max = +3.15);
  void SetAPRange(float max_alpha = 0.95, float max_qt = 0.05); // Armenteros Podolanski
  void SetRxyRange(float min = 0.f, float max = 180.f);
  void SetMinCosPA(float min = 0.95);
  void SetMaxPCA(float max = 2.f);
  void SetMaxChi2KF(float max = 1e+10);
  void SetMaxMarginZ(float max = 7.f);
  void SetMaxMeePsiPairDep(std::function<float(float)> psiDepCut);
  void SetOnWwireIB(bool flag = false);
  void SetOnWwireOB(bool flag = false);
  void RejectITSib(bool flag = false);

  void SetTrackPtRange(float minPt = 0.f, float maxPt = 1e10f);
  void SetTrackEtaRange(float minEta = -1e10f, float maxEta = 1e10f);
  void SetMinNClustersTPC(int minNClustersTPC);
  void SetMinNCrossedRowsTPC(int minNCrossedRowsTPC);
  void SetMinNCrossedRowsOverFindableClustersTPC(float minNCrossedRowsOverFindableClustersTPC);
  void SetMaxFracSharedClustersTPC(float max);
  void SetChi2PerClusterTPC(float min, float max);
  void SetNClustersITS(int min, int max);
  void SetChi2PerClusterITS(float min, float max);
  void SetMeanClusterSizeITSob(float min, float max);

  void SetTPCNsigmaElRange(float min = -3, float max = +3);
  void SetTPCNsigmaPiRange(float min = -1e+10, float max = 1e+10);

  void SetMaxDcaXY(float maxDcaXY);
  void SetMaxDcaZ(float maxDcaZ);
  void SetMaxDcaXYPtDep(std::function<float(float)> ptDepCut);
  void SetIsWithinBeamPipe(bool flag);
  void SetRequireITSTPC(bool flag);
  void SetRequireITSonly(bool flag);
  void SetRequireTPConly(bool flag);
  void SetRequireTPCTRD(bool flag);
  void SetRequireTPCTOF(bool flag);
  void SetDisableITSonly(bool flag);
  void SetDisableTPConly(bool flag);

  void SetApplyMlCuts(bool flag = false);
  void SetUse2DBinning(bool flag = true);
  void SetLoadMlModelsFromCCDB(bool flag = true);
  void SetNClassesMl(int nClasses);
  void SetMlTimestampCCDB(int timestamp);
  void SetCentrality(float centFT0A, float centFT0C, float centFT0M);
  void SetD_Bz(float d_bz);
  void SetCcdbUrl(const std::string& url = "http://alice-ccdb.cern.ch");
  void SetCentralityTypeMl(const std::string& centType);
  void SetCutDirMl(const std::vector<int>& cutDirMl);
  void SetMlModelPathsCCDB(const std::vector<std::string>& modelPaths);
  void SetMlOnnxFileNames(const std::vector<std::string>& onnxFileNamesVec);
  void SetLabelsBinsMl(const std::vector<std::string>& labelsBins);
  void SetLabelsCutScoresMl(const std::vector<std::string>& labelsCutScores);
  void SetBinsPtMl(const std::vector<double>& binsPt);
  void SetBinsCentMl(const std::vector<double>& binsCent);
  void SetCutsMl(const std::vector<double>& cutsMlFlat);
  void SetNamesInputFeatures(const std::vector<std::string>& namesInputFeaturesVec);

 private:
  static const std::pair<int8_t, std::set<uint8_t>> its_ib_Requirement;
  static const std::pair<int8_t, std::set<uint8_t>> its_ob_Requirement;
  static const std::pair<int8_t, std::set<uint8_t>> its_ob_Requirement_ITSTPC;
  // v0 cuts
  float mMinMee{0.f}, mMaxMee{0.1f};
  float mMinV0Pt{0.f}, mMaxV0Pt{1e10f};      // range in pT
  float mMinV0Eta{-1e10f}, mMaxV0Eta{1e10f}; // range in eta
  float mMaxAlpha{0.95}, mMaxQt{0.05};
  float mMinPsiPair{-3.15}, mMaxPsiPair{+3.15};
  float mMinPhivPair{0.f}, mMaxPhivPair{+3.15};
  float mMinRxy{0.f}, mMaxRxy{180.f};
  float mMinCosPA{0.95};
  float mMaxPCA{2.f};
  float mMaxChi2KF{1e+10};
  float mMaxMarginZ{7.f};
  std::function<float(float)> mMaxMeePsiPairDep{}; // max mee as a function of psipair
  bool mIsOnWwireIB{false};
  bool mIsOnWwireOB{false};
  bool mRejectITSib{false};

  // ML cuts
  bool mApplyMlCuts{false};
  bool mUse2DBinning{true};
  bool mLoadMlModelsFromCCDB{true};
  int mTimestampCCDB{-1};
  int mNClassesMl{static_cast<int>(o2::analysis::em_cuts_ml::NCutScores)};
  float mCentFT0A{0.f};
  float mCentFT0C{0.f};
  float mCentFT0M{0.f};
  float mD_Bz{0.f};
  std::string mCcdbUrl{"http://alice-ccdb.cern.ch"};
  std::string mCentralityTypeMl{"CentFT0C"};
  std::vector<int> mCutDirMl{std::vector<int>{o2::analysis::em_cuts_ml::vecCutDir}};
  std::vector<std::string> mModelPathsCCDB{std::vector<std::string>{"path_ccdb/BDT_PCM/"}};
  std::vector<std::string> mOnnxFileNames{std::vector<std::string>{"ModelHandler_onnx_PCM.onnx"}};
  std::vector<std::string> mNamesInputFeatures{std::vector<std::string>{"feature1", "feature2"}};
  std::vector<std::string> mLabelsBinsMl{std::vector<std::string>{"bin 0", "bin 1"}};
  std::vector<std::string> mLabelsCutScoresMl{std::vector<std::string>{"score primary photons", "score background"}};
  std::vector<double> mBinsPtMl{std::vector<double>{o2::analysis::em_cuts_ml::vecBinsPt}};
  std::vector<double> mBinsCentMl{std::vector<double>{o2::analysis::em_cuts_ml::vecBinsCent}};
  std::vector<double> mCutsMlFlat{std::vector<double>{0.5}};
  o2::analysis::EmMlResponsePCM<float>* mEmMlResponse{nullptr};

  // pid cuts
  float mMinTPCNsigmaEl{-5}, mMaxTPCNsigmaEl{+5};
  float mMinTPCNsigmaPi{-1e+10}, mMaxTPCNsigmaPi{+1e+10};

  // kinematic cuts
  float mMinTrackPt{0.f}, mMaxTrackPt{1e10f};      // range in pT
  float mMinTrackEta{-1e10f}, mMaxTrackEta{1e10f}; // range in eta

  // track quality cuts
  int mMinNClustersTPC{0};                                             // min number of TPC clusters
  int mMinNCrossedRowsTPC{0};                                          // min number of crossed rows in TPC
  float mMinChi2PerClusterTPC{-1e10f}, mMaxChi2PerClusterTPC{1e10f};   // max tpc fit chi2 per TPC cluster
  float mMinNCrossedRowsOverFindableClustersTPC{0.f};                  // min ratio crossed rows / findable clusters
  float mMaxFracSharedClustersTPC{999.f};                              // max ratio shared clusters / clusters in TPC
  int mMinNClustersITS{0}, mMaxNClustersITS{7};                        // range in number of ITS clusters
  float mMinChi2PerClusterITS{-1e10f}, mMaxChi2PerClusterITS{1e10f};   // max its fit chi2 per ITS cluster
  float mMinMeanClusterSizeITS{-1e10f}, mMaxMeanClusterSizeITS{1e10f}; // max <its cluster size> x cos(Lmabda)

  float mMaxDcaXY{1e10f};                       // max dca in xy plane
  float mMaxDcaZ{1e10f};                        // max dca in z direction
  std::function<float(float)> mMaxDcaXYPtDep{}; // max dca in xy plane as function of pT
  bool mRequireITSTPC{false};
  bool mRequireITSonly{false};
  bool mRequireTPConly{false};
  bool mRequireTPCTRD{false};
  bool mRequireTPCTOF{false};
  bool mDisableITSonly{false};
  bool mDisableTPConly{false};

  ClassDef(V0PhotonCut, 4);
};

#endif // PWGEM_PHOTONMESON_CORE_V0PHOTONCUT_H_
