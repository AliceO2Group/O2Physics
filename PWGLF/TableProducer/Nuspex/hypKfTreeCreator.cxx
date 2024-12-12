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
///
/// \brief Creates flat tree for ML analysis
/// \authors Janik Ditzel <jditzel@cern.ch> and Michael Hartung <mhartung@cern.ch>

#include <vector>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DCAFitter/DCAFitterN.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGLF/DataModel/LFHypernucleiKfTables.h"

using namespace o2;
using namespace o2::framework;
typedef std::array<float, 3> arr3;

namespace
{
std::vector<std::shared_ptr<TH1>> hPt;

struct trackProperties {
  trackProperties() : X(0), Y(0), Z(0), Px(0), Py(0), Pz(0), TPCnCls(0), ITSnCls(0), TPCchi2(0), ITSchi2(0), ITSmeanClsSize(0), ITSmeanClsSizeL(0), Rigidity(0), TPCsignal(0), TPCnSigma(0), TPCnSigmaNhp(0), TPCnSigmaNlp(0), TOFmass(0), DcaXY(0), DcaZ(0), IsPvContributor(0), SubMass(0) {}
  float X, Y, Z, Px, Py, Pz;
  uint8_t TPCnCls, ITSnCls;
  float TPCchi2, ITSchi2, ITSmeanClsSize, ITSmeanClsSizeL;
  float Rigidity, TPCsignal, TPCnSigma, TPCnSigmaNhp, TPCnSigmaNlp;
  float TOFmass, DcaXY, DcaZ;
  bool IsPvContributor;
  float SubMass;
};

struct hyperNucleus {
  hyperNucleus() : PdgCode(0), IsReconstructed(0), GlobalIndex(0), Species(0), IsMatter(0), PassedEvSel(0), IsMatterMC(0), PassedEvSelMC(0), IsPhysicalPrimary(0), CollisionMcTrue(0), Mass(0), Y(0), Pt(0), Ct(0), YGen(0), PtGen(0), CtGen(0), CpaPvGen(0), CpaPv(0), CpaSv(0), MaxDcaTracks(0), MaxDcaTracksSV(0), DcaToPvXY(0), DcaToPvZ(0), DcaToVtxXY(0), DcaToVtxZ(0), DevToPvXY(0), Chi2(0), Pvx(0), Pvy(0), Pvz(0), Svx(0), Svy(0), Svz(0), Px(0), Py(0), Pz(0), PvxGen(0), PvyGen(0), PvzGen(0), SvxGen(0), SvyGen(0), SvzGen(0), PxGen(0), PyGen(0), PzGen(0), NsingleDaughters(0), NcascadeDaughters(0), McTrue(0), MCTrueVtx(0), McPhysicalPrimary(0), HypNucDaughter(0) {}
  int PdgCode, IsReconstructed, GlobalIndex;
  uint8_t Species;
  bool IsMatter, PassedEvSel, IsMatterMC, PassedEvSelMC, IsPhysicalPrimary, CollisionMcTrue;
  float Mass, Y, Pt, Ct, YGen, PtGen, CtGen, CpaPvGen, CpaPv, CpaSv, MaxDcaTracks, MaxDcaTracksSV;
  float DcaToPvXY, DcaToPvZ, DcaToVtxXY, DcaToVtxZ, DevToPvXY, Chi2;
  float Pvx, Pvy, Pvz, Svx, Svy, Svz, Px, Py, Pz;
  float PvxGen, PvyGen, PvzGen, SvxGen, SvyGen, SvzGen, PxGen, PyGen, PzGen;
  int NsingleDaughters, NcascadeDaughters;
  bool McTrue, MCTrueVtx, McPhysicalPrimary;
  std::vector<trackProperties> daughterTracks;
  std::vector<float> subDaughterMassVec;
  hyperNucleus* HypNucDaughter;
  ~hyperNucleus()
  {
    if (HypNucDaughter)
      delete HypNucDaughter;
  }
};
} // namespace
namespace o2::aod
{
namespace hypkftree
{
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(Ct, ct, float);
DECLARE_SOA_COLUMN(YGen, yGen, float);
DECLARE_SOA_COLUMN(PtGen, ptGen, float);
DECLARE_SOA_COLUMN(CtGen, ctGen, float);
DECLARE_SOA_COLUMN(CpaPvGen, cpaPvGen, float);
DECLARE_SOA_COLUMN(DcaTracks, dcaTracks, float);
DECLARE_SOA_COLUMN(DcaTrackSv, dcaTrackSv, float);
DECLARE_SOA_COLUMN(CosPa, cosPa, double);
DECLARE_SOA_COLUMN(McTrue, mcTrue, bool);
DECLARE_SOA_COLUMN(Pvx, pvx, float);
DECLARE_SOA_COLUMN(Pvy, pvy, float);
DECLARE_SOA_COLUMN(Pvz, pvz, float);
DECLARE_SOA_COLUMN(Tvx, tvx, float);
DECLARE_SOA_COLUMN(Tvy, tvy, float);
DECLARE_SOA_COLUMN(Tvz, tvz, float);
DECLARE_SOA_COLUMN(PxGen, pxGen, float);
DECLARE_SOA_COLUMN(PyGen, pyGen, float);
DECLARE_SOA_COLUMN(PzGen, pzGen, float);
DECLARE_SOA_COLUMN(PvxGen, pvxGen, float);
DECLARE_SOA_COLUMN(PvyGen, pvyGen, float);
DECLARE_SOA_COLUMN(PvzGen, pvzGen, float);
DECLARE_SOA_COLUMN(SvxGen, svxGen, float);
DECLARE_SOA_COLUMN(SvyGen, svyGen, float);
DECLARE_SOA_COLUMN(SvzGen, svzGen, float);
DECLARE_SOA_COLUMN(TvxGen, tvxGen, float);
DECLARE_SOA_COLUMN(TvyGen, tvyGen, float);
DECLARE_SOA_COLUMN(TvzGen, tvzGen, float);
DECLARE_SOA_COLUMN(PassedEvSelMC, passedEvSelMC, bool);
DECLARE_SOA_COLUMN(IsMatter, isMatter, bool);
DECLARE_SOA_COLUMN(IsMatterGen, isMatterGen, bool);
DECLARE_SOA_COLUMN(IsReconstructed, isReconstructed, bool);
DECLARE_SOA_COLUMN(CollMcTrue, collMcTrue, bool);
DECLARE_SOA_COLUMN(D1X, d1X, float);
DECLARE_SOA_COLUMN(D1Y, d1Y, float);
DECLARE_SOA_COLUMN(D1Z, d1Z, float);
DECLARE_SOA_COLUMN(D1Px, d1Px, float);
DECLARE_SOA_COLUMN(D1Py, d1Py, float);
DECLARE_SOA_COLUMN(D1Pz, d1Pz, float);
DECLARE_SOA_COLUMN(D1TPCnCls, d1TPCnCls, uint8_t);
DECLARE_SOA_COLUMN(D1TPCchi2, d1TPCchi2, float);
DECLARE_SOA_COLUMN(D1ITSnCls, d1ITSnCls, uint8_t);
DECLARE_SOA_COLUMN(D1ITSchi2, d1ITSchi2, float);
DECLARE_SOA_COLUMN(D1ITSmeanClsSize, d1ITSmeanClsSize, float);
DECLARE_SOA_COLUMN(D1ITSmeanClsSizeL, d1ITSmeanClsSizeL, float);
DECLARE_SOA_COLUMN(D1Rigidity, d1Rigidity, float);
DECLARE_SOA_COLUMN(D1TPCsignal, d1TPCsignal, float);
DECLARE_SOA_COLUMN(D1TPCnSigma, d1TPCnSigma, float);
DECLARE_SOA_COLUMN(D1TPCnSigmaNhp, d1TPCnSigmaNhp, float);
DECLARE_SOA_COLUMN(D1TPCnSigmaNlp, d1TPCnSigmaNlp, float);
DECLARE_SOA_COLUMN(D1TOFmass, d1TOFmass, float);
DECLARE_SOA_COLUMN(D1DcaXY, d1DcaXY, float);
DECLARE_SOA_COLUMN(D1DcaZ, d1DcaZ, float);
DECLARE_SOA_COLUMN(D1IsPvContributor, d1IsPvContributor, bool);
DECLARE_SOA_COLUMN(D2X, d2X, float);
DECLARE_SOA_COLUMN(D2Y, d2Y, float);
DECLARE_SOA_COLUMN(D2Z, d2Z, float);
DECLARE_SOA_COLUMN(D2Px, d2Px, float);
DECLARE_SOA_COLUMN(D2Py, d2Py, float);
DECLARE_SOA_COLUMN(D2Pz, d2Pz, float);
DECLARE_SOA_COLUMN(D2TPCnCls, d2TPCnCls, uint8_t);
DECLARE_SOA_COLUMN(D2TPCchi2, d2TPCchi2, float);
DECLARE_SOA_COLUMN(D2ITSnCls, d2ITSnCls, uint8_t);
DECLARE_SOA_COLUMN(D2ITSchi2, d2ITSchi2, float);
DECLARE_SOA_COLUMN(D2ITSmeanClsSize, d2ITSmeanClsSize, float);
DECLARE_SOA_COLUMN(D2ITSmeanClsSizeL, d2ITSmeanClsSizeL, float);
DECLARE_SOA_COLUMN(D2Rigidity, d2Rigidity, float);
DECLARE_SOA_COLUMN(D2TPCsignal, d2TPCsignal, float);
DECLARE_SOA_COLUMN(D2TPCnSigma, d2TPCnSigma, float);
DECLARE_SOA_COLUMN(D2TPCnSigmaNhp, d2TPCnSigmaNhp, float);
DECLARE_SOA_COLUMN(D2TPCnSigmaNlp, d2TPCnSigmaNlp, float);
DECLARE_SOA_COLUMN(D2TOFmass, d2TOFmass, float);
DECLARE_SOA_COLUMN(D2DcaXY, d2DcaXY, float);
DECLARE_SOA_COLUMN(D2DcaZ, d2DcaZ, float);
DECLARE_SOA_COLUMN(D2IsPvContributor, d2IsPvContributor, bool);
DECLARE_SOA_COLUMN(D3X, d3X, float);
DECLARE_SOA_COLUMN(D3Y, d3Y, float);
DECLARE_SOA_COLUMN(D3Z, d3Z, float);
DECLARE_SOA_COLUMN(D3Px, d3Px, float);
DECLARE_SOA_COLUMN(D3Py, d3Py, float);
DECLARE_SOA_COLUMN(D3Pz, d3Pz, float);
DECLARE_SOA_COLUMN(D3TPCnCls, d3TPCnCls, uint8_t);
DECLARE_SOA_COLUMN(D3TPCchi2, d3TPCchi2, float);
DECLARE_SOA_COLUMN(D3ITSnCls, d3ITSnCls, uint8_t);
DECLARE_SOA_COLUMN(D3ITSchi2, d3ITSchi2, float);
DECLARE_SOA_COLUMN(D3ITSmeanClsSize, d3ITSmeanClsSize, float);
DECLARE_SOA_COLUMN(D3ITSmeanClsSizeL, d3ITSmeanClsSizeL, float);
DECLARE_SOA_COLUMN(D3Rigidity, d3Rigidity, float);
DECLARE_SOA_COLUMN(D3TPCsignal, d3TPCsignal, float);
DECLARE_SOA_COLUMN(D3TPCnSigma, d3TPCnSigma, float);
DECLARE_SOA_COLUMN(D3TPCnSigmaNhp, d3TPCnSigmaNhp, float);
DECLARE_SOA_COLUMN(D3TPCnSigmaNlp, d3TPCnSigmaNlp, float);
DECLARE_SOA_COLUMN(D3TOFmass, d3TOFmass, float);
DECLARE_SOA_COLUMN(D3DcaXY, d3DcaXY, float);
DECLARE_SOA_COLUMN(D3DcaZ, d3DcaZ, float);
DECLARE_SOA_COLUMN(D1D2Mass, d1d2Mass, float);
DECLARE_SOA_COLUMN(D1D3Mass, d1d3Mass, float);
DECLARE_SOA_COLUMN(D2D3Mass, d2d3Mass, float);
DECLARE_SOA_COLUMN(D3IsPvContributor, d3IsPvContributor, bool);
DECLARE_SOA_COLUMN(D0X, d0X, float);
DECLARE_SOA_COLUMN(D0Y, d0Y, float);
DECLARE_SOA_COLUMN(D0Z, d0Z, float);
DECLARE_SOA_COLUMN(D0Px, d0Px, float);
DECLARE_SOA_COLUMN(D0Py, d0Py, float);
DECLARE_SOA_COLUMN(D0Pz, d0Pz, float);
DECLARE_SOA_COLUMN(D0Mass, d0Mass, float);
DECLARE_SOA_COLUMN(D0Ct, d0ct, float);
DECLARE_SOA_COLUMN(D0CosPA, d0cosPa, float);
DECLARE_SOA_COLUMN(D0DcaTracks, d0dcaTracks, float);
DECLARE_SOA_COLUMN(D0DcaTracksTv, d0dcaTracksTv, float);
DECLARE_SOA_COLUMN(D0DcaToPvXY, d0dcaToPvXY, float);
DECLARE_SOA_COLUMN(D0DcaToPvZ, d0dcaToPvZ, float);
DECLARE_SOA_COLUMN(D0DcaToSvXY, d0dcaToSvXY, float);
DECLARE_SOA_COLUMN(D0DcaToSvZ, d0dcaToSvZ, float);
DECLARE_SOA_COLUMN(D0Chi2, d0chi2, float);
DECLARE_SOA_COLUMN(SD1X, sd1X, float);
DECLARE_SOA_COLUMN(SD1Y, sd1Y, float);
DECLARE_SOA_COLUMN(SD1Z, sd1Z, float);
DECLARE_SOA_COLUMN(SD1Px, sd1Px, float);
DECLARE_SOA_COLUMN(SD1Py, sd1Py, float);
DECLARE_SOA_COLUMN(SD1Pz, sd1Pz, float);
DECLARE_SOA_COLUMN(SD1TPCnCls, sd1TPCnCls, uint8_t);
DECLARE_SOA_COLUMN(SD1TPCchi2, sd1TPCchi2, float);
DECLARE_SOA_COLUMN(SD1ITSnCls, sd1ITSnCls, uint8_t);
DECLARE_SOA_COLUMN(SD1ITSchi2, sd1ITSchi2, float);
DECLARE_SOA_COLUMN(SD1ITSmeanClsSize, sd1ITSmeanClsSize, float);
DECLARE_SOA_COLUMN(SD1ITSmeanClsSizeL, sd1ITSmeanClsSizeL, float);
DECLARE_SOA_COLUMN(SD1Rigidity, sd1Rigidity, float);
DECLARE_SOA_COLUMN(SD1TPCsignal, sd1TPCsignal, float);
DECLARE_SOA_COLUMN(SD1TPCnSigma, sd1TPCnSigma, float);
DECLARE_SOA_COLUMN(SD1TPCnSigmaNhp, sd1TPCnSigmaNhp, float);
DECLARE_SOA_COLUMN(SD1TPCnSigmaNlp, sd1TPCnSigmaNlp, float);
DECLARE_SOA_COLUMN(SD1TOFmass, sd1TOFmass, float);
DECLARE_SOA_COLUMN(SD1DcaXY, sd1DcaXY, float);
DECLARE_SOA_COLUMN(SD1DcaZ, sd1DcaZ, float);
DECLARE_SOA_COLUMN(SD1IsPvContributor, sd1IsPvContributor, bool);
DECLARE_SOA_COLUMN(SD2X, sd2X, float);
DECLARE_SOA_COLUMN(SD2Y, sd2Y, float);
DECLARE_SOA_COLUMN(SD2Z, sd2Z, float);
DECLARE_SOA_COLUMN(SD2Px, sd2Px, float);
DECLARE_SOA_COLUMN(SD2Py, sd2Py, float);
DECLARE_SOA_COLUMN(SD2Pz, sd2Pz, float);
DECLARE_SOA_COLUMN(SD2TPCnCls, sd2TPCnCls, uint8_t);
DECLARE_SOA_COLUMN(SD2TPCchi2, sd2TPCchi2, float);
DECLARE_SOA_COLUMN(SD2ITSnCls, sd2ITSnCls, uint8_t);
DECLARE_SOA_COLUMN(SD2ITSchi2, sd2ITSchi2, float);
DECLARE_SOA_COLUMN(SD2ITSmeanClsSize, sd2ITSmeanClsSize, float);
DECLARE_SOA_COLUMN(SD2ITSmeanClsSizeL, sd2ITSmeanClsSizeL, float);
DECLARE_SOA_COLUMN(SD2Rigidity, sd2Rigidity, float);
DECLARE_SOA_COLUMN(SD2TPCsignal, sd2TPCsignal, float);
DECLARE_SOA_COLUMN(SD2TPCnSigma, sd2TPCnSigma, float);
DECLARE_SOA_COLUMN(SD2TPCnSigmaNhp, sd2TPCnSigmaNhp, float);
DECLARE_SOA_COLUMN(SD2TPCnSigmaNlp, sd2TPCnSigmaNlp, float);
DECLARE_SOA_COLUMN(SD2TOFmass, sd2TOFmass, float);
DECLARE_SOA_COLUMN(SD2DcaXY, sd2DcaXY, float);
DECLARE_SOA_COLUMN(SD2DcaZ, sd2DcaZ, float);
DECLARE_SOA_COLUMN(SD2IsPvContributor, sd2IsPvContributor, bool);
DECLARE_SOA_COLUMN(SD3X, sd3X, float);
DECLARE_SOA_COLUMN(SD3Y, sd3Y, float);
DECLARE_SOA_COLUMN(SD3Z, sd3Z, float);
DECLARE_SOA_COLUMN(SD3Px, sd3Px, float);
DECLARE_SOA_COLUMN(SD3Py, sd3Py, float);
DECLARE_SOA_COLUMN(SD3Pz, sd3Pz, float);
DECLARE_SOA_COLUMN(SD3TPCnCls, sd3TPCnCls, uint8_t);
DECLARE_SOA_COLUMN(SD3TPCchi2, sd3TPCchi2, float);
DECLARE_SOA_COLUMN(SD3ITSnCls, sd3ITSnCls, uint8_t);
DECLARE_SOA_COLUMN(SD3ITSchi2, sd3ITSchi2, float);
DECLARE_SOA_COLUMN(SD3ITSmeanClsSize, sd3ITSmeanClsSize, float);
DECLARE_SOA_COLUMN(SD3ITSmeanClsSizeL, sd3ITSmeanClsSizeL, float);
DECLARE_SOA_COLUMN(SD3Rigidity, sd3Rigidity, float);
DECLARE_SOA_COLUMN(SD3TPCsignal, sd3TPCsignal, float);
DECLARE_SOA_COLUMN(SD3TPCnSigma, sd3TPCnSigma, float);
DECLARE_SOA_COLUMN(SD3TPCnSigmaNhp, sd3TPCnSigmaNhp, float);
DECLARE_SOA_COLUMN(SD3TPCnSigmaNlp, sd3TPCnSigmaNlp, float);
DECLARE_SOA_COLUMN(SD3TOFmass, sd3TOFmass, float);
DECLARE_SOA_COLUMN(SD3DcaXY, sd3DcaXY, float);
DECLARE_SOA_COLUMN(SD3DcaZ, sd3DcaZ, float);
DECLARE_SOA_COLUMN(SD3IsPvContributor, sd3IsPvContributor, bool);
DECLARE_SOA_COLUMN(SD1SD2Mass, sd1sd2Mass, float);
DECLARE_SOA_COLUMN(SD1SD3Mass, sd1sd3Mass, float);
DECLARE_SOA_COLUMN(SD2SD3Mass, sd2sd3Mass, float);
} // namespace hypkftree

#define HYPKFGENBASE mcparticle::PdgCode, hypkftree::IsMatterGen, hypkftree::IsReconstructed, hykfmc::IsPhysicalPrimary, hypkftree::PassedEvSelMC, hypkftree::YGen, hypkftree::PtGen, hypkftree::CtGen

#define HYPKFGENEXT hypkftree::CpaPvGen, hypkftree::PxGen, hypkftree::PyGen, hypkftree::PzGen, hypkftree::PvxGen, hypkftree::PvyGen, hypkftree::PvzGen, hypkftree::SvxGen, hypkftree::SvyGen, hypkftree::SvzGen

#define HYPKFGENCAS hypkftree::TvxGen, hypkftree::TvyGen, hypkftree::TvzGen

#define HYPKFHYPNUC hykfmc::Species, hypkftree::IsMatter, hykfmcColl::PassedEvSel, hykfhyp::Mass, hypkftree::Y, track::Pt, hypkftree::Ct, hypkftree::CosPa, hypkftree::DcaTracks, hykfhyp::DcaToPvXY, hykfhyp::DcaToPvZ, hykfhyp::DevToPvXY, hykfhyp::Chi2, hypkftree::Pvx, hypkftree::Pvy, hypkftree::Pvz, hykfmc::Svx, hykfmc::Svy, hykfmc::Svz, hykfhyp::Px, hykfhyp::Py, hykfhyp::Pz, hypkftree::CollMcTrue

#define HYPKFHYPNUCMC hypkftree::McTrue, hykfmc::IsPhysicalPrimary

#define HYPKFD0 hypkftree::Tvx, hypkftree::Tvy, hypkftree::Tvz, hypkftree::D0X, hypkftree::D0Y, hypkftree::D0Z, hypkftree::D0Px, hypkftree::D0Py, hypkftree::D0Pz, hypkftree::D0Mass, hypkftree::D0Ct, hypkftree::D0CosPA, hypkftree::D0DcaTracks, hypkftree::D0DcaToPvXY, hypkftree::D0DcaToPvZ, hypkftree::D0DcaToSvXY, hypkftree::D0DcaToSvZ, hypkftree::D0Chi2

#define HYPKFD1 hypkftree::D1X, hypkftree::D1Y, hypkftree::D1Z, hypkftree::D1Px, hypkftree::D1Py, hypkftree::D1Pz, hypkftree::D1TPCnCls, hypkftree::D1TPCchi2, hypkftree::D1ITSnCls, hypkftree::D1ITSchi2, hypkftree::D1ITSmeanClsSize, hypkftree::D1ITSmeanClsSizeL, hypkftree::D1Rigidity, hypkftree::D1TPCsignal, hypkftree::D1TPCnSigma, hypkftree::D1TPCnSigmaNhp, hypkftree::D1TPCnSigmaNlp, hypkftree::D1TOFmass, hypkftree::D1DcaXY, hypkftree::D1DcaZ, hypkftree::D1IsPvContributor

#define HYPKFD2 hypkftree::D2X, hypkftree::D2Y, hypkftree::D2Z, hypkftree::D2Px, hypkftree::D2Py, hypkftree::D2Pz, hypkftree::D2TPCnCls, hypkftree::D2TPCchi2, hypkftree::D2ITSnCls, hypkftree::D2ITSchi2, hypkftree::D2ITSmeanClsSize, hypkftree::D2ITSmeanClsSizeL, hypkftree::D2Rigidity, hypkftree::D2TPCsignal, hypkftree::D2TPCnSigma, hypkftree::D2TPCnSigmaNhp, hypkftree::D2TPCnSigmaNlp, hypkftree::D2TOFmass, hypkftree::D2DcaXY, hypkftree::D2DcaZ, hypkftree::D2IsPvContributor

#define HYPKFD3 hypkftree::D3X, hypkftree::D3Y, hypkftree::D3Z, hypkftree::D3Px, hypkftree::D3Py, hypkftree::D3Pz, hypkftree::D3TPCnCls, hypkftree::D3TPCchi2, hypkftree::D3ITSnCls, hypkftree::D3ITSchi2, hypkftree::D3ITSmeanClsSize, hypkftree::D3ITSmeanClsSizeL, hypkftree::D3Rigidity, hypkftree::D3TPCsignal, hypkftree::D3TPCnSigma, hypkftree::D3TPCnSigmaNhp, hypkftree::D3TPCnSigmaNlp, hypkftree::D3TOFmass, hypkftree::D3DcaXY, hypkftree::D3DcaZ, hypkftree::D3IsPvContributor

#define HYPKFSD1 hypkftree::SD1X, hypkftree::SD1Y, hypkftree::SD1Z, hypkftree::SD1Px, hypkftree::SD1Py, hypkftree::SD1Pz, hypkftree::SD1TPCnCls, hypkftree::SD1TPCchi2, hypkftree::SD1ITSnCls, hypkftree::SD1ITSchi2, hypkftree::SD1ITSmeanClsSize, hypkftree::SD1ITSmeanClsSizeL, hypkftree::SD1Rigidity, hypkftree::SD1TPCsignal, hypkftree::SD1TPCnSigma, hypkftree::SD1TPCnSigmaNhp, hypkftree::SD1TPCnSigmaNlp, hypkftree::SD1TOFmass, hypkftree::SD1DcaXY, hypkftree::SD1DcaZ, hypkftree::SD1IsPvContributor

#define HYPKFSD2 hypkftree::SD2X, hypkftree::SD2Y, hypkftree::SD2Z, hypkftree::SD2Px, hypkftree::SD2Py, hypkftree::SD2Pz, hypkftree::SD2TPCnCls, hypkftree::SD2TPCchi2, hypkftree::SD2ITSnCls, hypkftree::SD2ITSchi2, hypkftree::SD2ITSmeanClsSize, hypkftree::SD2ITSmeanClsSizeL, hypkftree::SD2Rigidity, hypkftree::SD2TPCsignal, hypkftree::SD2TPCnSigma, hypkftree::SD2TPCnSigmaNhp, hypkftree::SD2TPCnSigmaNlp, hypkftree::SD2TOFmass, hypkftree::SD2DcaXY, hypkftree::SD2DcaZ, hypkftree::SD2IsPvContributor

#define HYPKFSD3 hypkftree::SD3X, hypkftree::SD3Y, hypkftree::SD3Z, hypkftree::SD3Px, hypkftree::SD3Py, hypkftree::SD3Pz, hypkftree::SD3TPCnCls, hypkftree::SD3TPCchi2, hypkftree::SD3ITSnCls, hypkftree::SD3ITSchi2, hypkftree::SD3ITSmeanClsSize, hypkftree::SD3ITSmeanClsSizeL, hypkftree::SD3Rigidity, hypkftree::SD3TPCsignal, hypkftree::SD3TPCnSigma, hypkftree::SD3TPCnSigmaNhp, hypkftree::SD3TPCnSigmaNlp, hypkftree::SD3TOFmass, hypkftree::SD3DcaXY, hypkftree::SD3DcaZ, hypkftree::SD3IsPvContributor

#define HYPKFSDMASS hypkftree::D1D2Mass, hypkftree::D1D3Mass, hypkftree::D2D3Mass
#define HYPKFSSDMASS hypkftree::SD1SD2Mass, hypkftree::SD1SD3Mass, hypkftree::SD2SD3Mass

DECLARE_SOA_TABLE(HypKfGens, "AOD", "HYPKFGEN", HYPKFGENBASE);
using HypKfGen = HypKfGens::iterator;

DECLARE_SOA_TABLE(HypKfSingleTwoBodyCandidates, "AOD", "HYPKFCAND2", HYPKFHYPNUC, HYPKFHYPNUCMC, HYPKFD1, HYPKFD2);
using HypKfSingleTwoBodyCandidate = HypKfSingleTwoBodyCandidates::iterator;

DECLARE_SOA_TABLE(HypKfMcSingleTwoBodyCandidates, "AOD", "HYPKFMCCAND2", HYPKFGENBASE, HYPKFGENEXT, HYPKFHYPNUC, HYPKFD1, HYPKFD2);
using HypKfMcSingleTwoBodyCandidate = HypKfMcSingleTwoBodyCandidates::iterator;

DECLARE_SOA_TABLE(HypKfSingleThreeBodyCandidates, "AOD", "HYPKFCAND3", HYPKFHYPNUC, HYPKFHYPNUCMC, HYPKFD1, HYPKFD2, HYPKFD3, HYPKFSDMASS);
using HypKfSingleThreeBodyCandidate = HypKfSingleThreeBodyCandidates::iterator;

DECLARE_SOA_TABLE(HypKfMcSingleThreeBodyCandidates, "AOD", "HYPKFMCCAND3", HYPKFGENBASE, HYPKFGENEXT, HYPKFHYPNUC, HYPKFD1, HYPKFD2, HYPKFD3, HYPKFSDMASS);
using HypKfMcSingleThreeBodyCandidate = HypKfMcSingleThreeBodyCandidates::iterator;

DECLARE_SOA_TABLE(HypKfCascadeTwoThreeCandidates, "AOD", "HYPKFCAND23", HYPKFHYPNUC, HYPKFHYPNUCMC, HYPKFD0, HYPKFD1, HYPKFSD1, HYPKFSD2, HYPKFSD3, HYPKFSSDMASS);
using HypKfCascadeTwoThreeCandidate = HypKfCascadeTwoThreeCandidates::iterator;

DECLARE_SOA_TABLE(HypKfMcCascadeTwoThreeCandidates, "AOD", "HYPKFMCCAND23", HYPKFGENBASE, HYPKFGENEXT, HYPKFHYPNUC, HYPKFD0, HYPKFD1, HYPKFSD1, HYPKFSD2, HYPKFSD3, HYPKFSSDMASS);
using HypKfMcCascadeTwoThreeCandidate = HypKfMcCascadeTwoThreeCandidates::iterator;

DECLARE_SOA_TABLE(HypKfCascadeThreeTwoCandidates, "AOD", "HYPKFCAND32", HYPKFHYPNUC, HYPKFHYPNUCMC, HYPKFD0, HYPKFD1, HYPKFD2, HYPKFSDMASS, HYPKFSD1, HYPKFSD2);
using HypKfCascadeThreeTwoCandidate = HypKfCascadeThreeTwoCandidates::iterator;

DECLARE_SOA_TABLE(HypKfMcCascadeThreeTwoCandidates, "AOD", "HYPKFMCCAND32", HYPKFGENBASE, HYPKFGENEXT, HYPKFHYPNUC, HYPKFD0, HYPKFD1, HYPKFD2, HYPKFSDMASS, HYPKFSD1, HYPKFSD2);
using HypKfMcCascadeThreeTwoCandidate = HypKfMcCascadeThreeTwoCandidates::iterator;
} // namespace o2::aod

struct hypKfTreeCreator {

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Produces<aod::HypKfGens> outputMcGenTable;
  Produces<aod::HypKfSingleTwoBodyCandidates> outputTableTwo;
  Produces<aod::HypKfMcSingleTwoBodyCandidates> outputTableMcTwo;
  Produces<aod::HypKfSingleThreeBodyCandidates> outputTableThree;
  Produces<aod::HypKfMcSingleThreeBodyCandidates> outputTableMcThree;
  Produces<aod::HypKfCascadeTwoThreeCandidates> outputTableTwoThree;
  Produces<aod::HypKfMcCascadeTwoThreeCandidates> outputTableMcTwoThree;
  Produces<aod::HypKfCascadeThreeTwoCandidates> outputTableThreeTwo;
  Produces<aod::HypKfMcCascadeThreeTwoCandidates> outputTableMcThreeTwo;
  PresliceUnsorted<aod::HypKfHypNucs> perMcParticle = aod::hykfhyp::hypKfMcPartId;

  Configurable<int> cfgSpecies{"cfgSpecies", 0, "Select species"};
  Configurable<int> cfgNprimDaughters{"cfgNprimDaughters", 0, "Number of primary daughters"};
  Configurable<int> cfgNsecDaughters{"cfgNsecDaughters", 0, "Number of secondary daughters (cascades only)"};
  Configurable<bool> cfgMCGenerated{"cfgMCGenerated", false, "create MC generated tree"};
  Configurable<bool> cfgMCReconstructed{"cfgMCReconstructed", false, "create MC reconstructed tree"};
  Configurable<bool> cfgMCCombined{"cfgMCCombined", false, "create MC tree containig generated and reconstructed"};

  bool isMC;
  //___________________________________________________________________________________________________________________________________________________________

  void init(InitContext const&)
  {
    const AxisSpec axisPt{10, 0., 10., "#it{p}_{T} (GeV/#it{c})"};
    hPt.resize(3);
    hPt[0] = histos.add<TH1>("hGen", "", HistType::kTH1F, {axisPt});
    hPt[0]->Sumw2();
    hPt[1] = histos.add<TH1>("hRec", "", HistType::kTH1F, {axisPt});
    hPt[1]->Sumw2();
    hPt[2] = histos.add<TH1>("hEff", "", HistType::kTH1F, {axisPt});
    isMC = false;
  }
  //___________________________________________________________________________________________________________________________________________________________

  void processData(aod::HypKfHypNucs const& hypNucs, aod::HypKfColls const& hypKfColls, aod::HypKfTracks const& hypKfTrks, aod::HypKfDaughtAdds const& hypKfDAdd, aod::HypKfSubDs const& hypKfDSub)
  {
    for (auto& hypNuc : hypNucs) {
      if (std::abs(hypNuc.species()) != cfgSpecies)
        continue;
      hyperNucleus candidate, hypDaughter, dummy;
      fillCandidate(candidate, hypDaughter, hypNuc, hypNucs, hypKfColls, hypKfTrks, hypKfDAdd, hypKfDSub);
      if (cfgNsecDaughters) {
        fillCandidate(hypDaughter, dummy, hypNucs.rawIteratorAt(hypNuc.hypDaughterId()), hypNucs, hypKfColls, hypKfTrks, hypKfDAdd, hypKfDSub);
      }
      fillTable(candidate, hypDaughter);
    }
  }
  PROCESS_SWITCH(hypKfTreeCreator, processData, "single tree", false);
  //___________________________________________________________________________________________________________________________________________________________
  void fillTable(hyperNucleus& cand, hyperNucleus& hypDaughter)
  {
    if (isMC && cfgMCGenerated)
      outputMcGenTable(
        cand.PdgCode, cand.IsMatterMC, cand.IsReconstructed, cand.IsPhysicalPrimary, cand.PassedEvSelMC, cand.YGen, cand.PtGen, cand.CtGen);

    if (!cand.IsReconstructed) {
      cand.daughterTracks.resize(4);
      cand.subDaughterMassVec.resize(4);
      hypDaughter.daughterTracks.resize(4);
      hypDaughter.subDaughterMassVec.resize(4);
    }

    if (cfgNprimDaughters == 2 && cfgNsecDaughters == 0) {
      const auto& d1 = cand.daughterTracks.at(0);
      const auto& d2 = cand.daughterTracks.at(1);
      if (!isMC || (isMC && cfgMCReconstructed && cand.IsReconstructed))
        outputTableTwo(
          cand.Species, cand.IsMatter, cand.PassedEvSel, cand.Mass, cand.Y, cand.Pt, cand.Ct, cand.CpaPv, cand.MaxDcaTracks, cand.DcaToPvXY, cand.DcaToPvZ, cand.DevToPvXY,
          cand.Chi2, cand.Pvx, cand.Pvy, cand.Pvz, cand.Svx, cand.Svy, cand.Svz, cand.Px, cand.Py, cand.Pz, cand.CollisionMcTrue, cand.McTrue, cand.McPhysicalPrimary,
          d1.X, d1.Y, d1.Z, d1.Px, d1.Py, d1.Pz, d1.TPCnCls, d1.TPCchi2, d1.ITSnCls, d1.ITSchi2, d1.ITSmeanClsSize, d1.ITSmeanClsSizeL,
          d1.Rigidity, d1.TPCsignal, d1.TPCnSigma, d1.TPCnSigmaNhp, d1.TPCnSigmaNlp, d1.TOFmass, d1.DcaXY, d1.DcaZ, d1.IsPvContributor,
          d2.X, d2.Y, d2.Z, d2.Px, d2.Py, d2.Pz, d2.TPCnCls, d2.TPCchi2, d2.ITSnCls, d2.ITSchi2, d2.ITSmeanClsSize, d2.ITSmeanClsSizeL,
          d2.Rigidity, d2.TPCsignal, d2.TPCnSigma, d2.TPCnSigmaNhp, d2.TPCnSigmaNlp, d2.TOFmass, d2.DcaXY, d2.DcaZ, d2.IsPvContributor);
      if (isMC && cfgMCCombined)
        outputTableMcTwo(
          cand.PdgCode, cand.IsMatterMC, cand.IsReconstructed, cand.IsPhysicalPrimary, cand.PassedEvSelMC, cand.YGen, cand.PtGen, cand.CtGen,
          cand.CpaPvGen, cand.PxGen, cand.PyGen, cand.PzGen, cand.PvxGen, cand.PvyGen, cand.PvzGen, cand.SvxGen, cand.SvyGen, cand.SvzGen,
          cand.Species, cand.IsMatter, cand.PassedEvSel, cand.Mass, cand.Y, cand.Pt, cand.Ct, cand.CpaPv, cand.MaxDcaTracks, cand.DcaToPvXY, cand.DcaToPvZ, cand.DevToPvXY,
          cand.Chi2, cand.Pvx, cand.Pvy, cand.Pvz, cand.Svx, cand.Svy, cand.Svz, cand.Px, cand.Py, cand.Pz, cand.CollisionMcTrue,
          d1.X, d1.Y, d1.Z, d1.Px, d1.Py, d1.Pz, d1.TPCnCls, d1.TPCchi2, d1.ITSnCls, d1.ITSchi2, d1.ITSmeanClsSize, d1.ITSmeanClsSizeL,
          d1.Rigidity, d1.TPCsignal, d1.TPCnSigma, d1.TPCnSigmaNhp, d1.TPCnSigmaNlp, d1.TOFmass, d1.DcaXY, d1.DcaZ, d1.IsPvContributor,
          d2.X, d2.Y, d2.Z, d2.Px, d2.Py, d2.Pz, d2.TPCnCls, d2.TPCchi2, d2.ITSnCls, d2.ITSchi2, d2.ITSmeanClsSize, d2.ITSmeanClsSizeL,
          d2.Rigidity, d2.TPCsignal, d2.TPCnSigma, d2.TPCnSigmaNhp, d2.TPCnSigmaNlp, d2.TOFmass, d2.DcaXY, d2.DcaZ, d2.IsPvContributor);
    }
    if (cfgNprimDaughters == 3 && cfgNsecDaughters == 0) {
      const auto& d1 = cand.daughterTracks.at(0);
      const auto& d2 = cand.daughterTracks.at(1);
      const auto& d3 = cand.daughterTracks.at(2);
      if (!isMC || (isMC && cfgMCReconstructed && cand.IsReconstructed))
        outputTableThree(
          cand.Species, cand.IsMatter, cand.PassedEvSel, cand.Mass, cand.Y, cand.Pt, cand.Ct, cand.CpaPv, cand.MaxDcaTracks, cand.DcaToPvXY, cand.DcaToPvZ, cand.DevToPvXY,
          cand.Chi2, cand.Pvx, cand.Pvy, cand.Pvz, cand.Svx, cand.Svy, cand.Svz, cand.Px, cand.Py, cand.Pz, cand.CollisionMcTrue, cand.McTrue, cand.McPhysicalPrimary,
          d1.X, d1.Y, d1.Z, d1.Px, d1.Py, d1.Pz, d1.TPCnCls, d1.TPCchi2, d1.ITSnCls, d1.ITSchi2, d1.ITSmeanClsSize, d1.ITSmeanClsSizeL,
          d1.Rigidity, d1.TPCsignal, d1.TPCnSigma, d1.TPCnSigmaNhp, d1.TPCnSigmaNlp, d1.TOFmass, d1.DcaXY, d1.DcaZ, d1.IsPvContributor,
          d2.X, d2.Y, d2.Z, d2.Px, d2.Py, d2.Pz, d2.TPCnCls, d2.TPCchi2, d2.ITSnCls, d2.ITSchi2, d2.ITSmeanClsSize, d2.ITSmeanClsSizeL,
          d2.Rigidity, d2.TPCsignal, d2.TPCnSigma, d2.TPCnSigmaNhp, d2.TPCnSigmaNlp, d2.TOFmass, d2.DcaXY, d2.DcaZ, d2.IsPvContributor,
          d3.X, d3.Y, d3.Z, d3.Px, d3.Py, d3.Pz, d3.TPCnCls, d3.TPCchi2, d3.ITSnCls, d3.ITSchi2, d3.ITSmeanClsSize, d3.ITSmeanClsSizeL,
          d3.Rigidity, d3.TPCsignal, d3.TPCnSigma, d3.TPCnSigmaNhp, d3.TPCnSigmaNlp, d3.TOFmass, d3.DcaXY, d3.DcaZ, d3.IsPvContributor,
          d1.SubMass, d2.SubMass, d3.SubMass);
      if (isMC && cfgMCCombined)
        outputTableMcThree(
          cand.PdgCode, cand.IsMatterMC, cand.IsReconstructed, cand.IsPhysicalPrimary, cand.PassedEvSelMC, cand.YGen, cand.PtGen, cand.CtGen,
          cand.CpaPvGen, cand.PxGen, cand.PyGen, cand.PzGen, cand.PvxGen, cand.PvyGen, cand.PvzGen, cand.SvxGen, cand.SvyGen, cand.SvzGen,
          cand.Species, cand.IsMatter, cand.PassedEvSel, cand.Mass, cand.Y, cand.Pt, cand.Ct, cand.CpaPv, cand.MaxDcaTracks, cand.DcaToPvXY, cand.DcaToPvZ, cand.DevToPvXY,
          cand.Chi2, cand.Pvx, cand.Pvy, cand.Pvz, cand.Svx, cand.Svy, cand.Svz, cand.Px, cand.Py, cand.Pz, cand.CollisionMcTrue,
          d1.X, d1.Y, d1.Z, d1.Px, d1.Py, d1.Pz, d1.TPCnCls, d1.TPCchi2, d1.ITSnCls, d1.ITSchi2, d1.ITSmeanClsSize, d1.ITSmeanClsSizeL,
          d1.Rigidity, d1.TPCsignal, d1.TPCnSigma, d1.TPCnSigmaNhp, d1.TPCnSigmaNlp, d1.TOFmass, d1.DcaXY, d1.DcaZ, d1.IsPvContributor,
          d2.X, d2.Y, d2.Z, d2.Px, d2.Py, d2.Pz, d2.TPCnCls, d2.TPCchi2, d2.ITSnCls, d2.ITSchi2, d2.ITSmeanClsSize, d2.ITSmeanClsSizeL,
          d2.Rigidity, d2.TPCsignal, d2.TPCnSigma, d2.TPCnSigmaNhp, d2.TPCnSigmaNlp, d2.TOFmass, d2.DcaXY, d2.DcaZ, d2.IsPvContributor,
          d3.X, d3.Y, d3.Z, d3.Px, d3.Py, d3.Pz, d3.TPCnCls, d3.TPCchi2, d3.ITSnCls, d3.ITSchi2, d3.ITSmeanClsSize, d3.ITSmeanClsSizeL,
          d3.Rigidity, d3.TPCsignal, d3.TPCnSigma, d3.TPCnSigmaNhp, d3.TPCnSigmaNlp, d3.TOFmass, d3.DcaXY, d3.DcaZ, d3.IsPvContributor,
          d1.SubMass, d2.SubMass, d3.SubMass);
    }
    if (cfgNprimDaughters == 2 && cfgNsecDaughters == 3) {
      const auto& d0 = cand.daughterTracks.at(0);
      const auto& d1 = cand.daughterTracks.at(1);
      const auto& sd1 = hypDaughter.daughterTracks.at(0);
      const auto& sd2 = hypDaughter.daughterTracks.at(1);
      const auto& sd3 = hypDaughter.daughterTracks.at(2);
      if (!isMC || (isMC && cfgMCReconstructed && cand.IsReconstructed))
        outputTableTwoThree(
          cand.Species, cand.IsMatter, cand.PassedEvSel, cand.Mass, cand.Y, cand.Pt, cand.Ct, cand.CpaPv, cand.MaxDcaTracks, cand.DcaToPvXY, cand.DcaToPvZ, cand.DevToPvXY,
          cand.Chi2, cand.Pvx, cand.Pvy, cand.Pvz, cand.Svx, cand.Svy, cand.Svz, cand.Px, cand.Py, cand.Pz, cand.CollisionMcTrue, cand.McTrue, cand.McPhysicalPrimary,
          hypDaughter.Svx, hypDaughter.Svy, hypDaughter.Svz, d0.X, d0.Y, d0.Z, d0.Px, d0.Py, d0.Pz, hypDaughter.Mass, hypDaughter.Ct, hypDaughter.CpaPv,
          hypDaughter.MaxDcaTracks, hypDaughter.DcaToPvXY, hypDaughter.DcaToPvZ, hypDaughter.DcaToVtxXY, hypDaughter.DcaToVtxZ, hypDaughter.Chi2,
          d1.X, d1.Y, d1.Z, d1.Px, d1.Py, d1.Pz, d1.TPCnCls, d1.TPCchi2, d1.ITSnCls, d1.ITSchi2, d1.ITSmeanClsSize, d1.ITSmeanClsSizeL,
          d1.Rigidity, d1.TPCsignal, d1.TPCnSigma, d1.TPCnSigmaNhp, d1.TPCnSigmaNlp, d1.TOFmass, d1.DcaXY, d1.DcaZ, d1.IsPvContributor,
          sd1.X, sd1.Y, sd1.Z, sd1.Px, sd1.Py, sd1.Pz, sd1.TPCnCls, sd1.TPCchi2, sd1.ITSnCls, sd1.ITSchi2, sd1.ITSmeanClsSize, sd1.ITSmeanClsSizeL,
          sd1.Rigidity, sd1.TPCsignal, sd1.TPCnSigma, sd1.TPCnSigmaNhp, sd1.TPCnSigmaNlp, sd1.TOFmass, sd1.DcaXY, sd1.DcaZ, sd1.IsPvContributor,
          sd2.X, sd2.Y, sd2.Z, sd2.Px, sd2.Py, sd2.Pz, sd2.TPCnCls, sd2.TPCchi2, sd2.ITSnCls, sd2.ITSchi2, sd2.ITSmeanClsSize, sd2.ITSmeanClsSizeL,
          sd2.Rigidity, sd2.TPCsignal, sd2.TPCnSigma, sd2.TPCnSigmaNhp, sd2.TPCnSigmaNlp, sd2.TOFmass, sd2.DcaXY, sd2.DcaZ, sd2.IsPvContributor,
          sd3.X, sd3.Y, sd3.Z, sd3.Px, sd3.Py, sd3.Pz, sd3.TPCnCls, sd3.TPCchi2, sd3.ITSnCls, sd3.ITSchi2, sd3.ITSmeanClsSize, sd3.ITSmeanClsSizeL,
          sd3.Rigidity, sd3.TPCsignal, sd3.TPCnSigma, sd3.TPCnSigmaNhp, sd3.TPCnSigmaNlp, sd3.TOFmass, sd3.DcaXY, sd3.DcaZ, sd3.IsPvContributor,
          sd1.SubMass, sd2.SubMass, sd3.SubMass);
      if (isMC && cfgMCCombined)
        outputTableMcTwoThree(
          cand.PdgCode, cand.IsMatterMC, cand.IsReconstructed, cand.IsPhysicalPrimary, cand.PassedEvSelMC, cand.YGen, cand.PtGen, cand.CtGen,
          cand.CpaPvGen, cand.PxGen, cand.PyGen, cand.PzGen, cand.PvxGen, cand.PvyGen, cand.PvzGen, cand.SvxGen, cand.SvyGen, cand.SvzGen,
          cand.Species, cand.IsMatter, cand.PassedEvSel, cand.Mass, cand.Y, cand.Pt, cand.Ct, cand.CpaPv, cand.MaxDcaTracks, cand.DcaToPvXY, cand.DcaToPvZ, cand.DevToPvXY,
          cand.Chi2, cand.Pvx, cand.Pvy, cand.Pvz, cand.Svx, cand.Svy, cand.Svz, cand.Px, cand.Py, cand.Pz, cand.CollisionMcTrue,
          hypDaughter.Svx, hypDaughter.Svy, hypDaughter.Svz, d0.X, d0.Y, d0.Z, d0.Px, d0.Py, d0.Pz, hypDaughter.Mass, hypDaughter.Ct, hypDaughter.CpaPv,
          hypDaughter.MaxDcaTracks, hypDaughter.DcaToPvXY, hypDaughter.DcaToPvZ, hypDaughter.DcaToVtxXY, hypDaughter.DcaToVtxZ, hypDaughter.Chi2,
          d1.X, d1.Y, d1.Z, d1.Px, d1.Py, d1.Pz, d1.TPCnCls, d1.TPCchi2, d1.ITSnCls, d1.ITSchi2, d1.ITSmeanClsSize, d1.ITSmeanClsSizeL,
          d1.Rigidity, d1.TPCsignal, d1.TPCnSigma, d1.TPCnSigmaNhp, d1.TPCnSigmaNlp, d1.TOFmass, d1.DcaXY, d1.DcaZ, d1.IsPvContributor,
          sd1.X, sd1.Y, sd1.Z, sd1.Px, sd1.Py, sd1.Pz, sd1.TPCnCls, sd1.TPCchi2, sd1.ITSnCls, sd1.ITSchi2, sd1.ITSmeanClsSize, sd1.ITSmeanClsSizeL,
          sd1.Rigidity, sd1.TPCsignal, sd1.TPCnSigma, sd1.TPCnSigmaNhp, sd1.TPCnSigmaNlp, sd1.TOFmass, sd1.DcaXY, sd1.DcaZ, sd1.IsPvContributor,
          sd2.X, sd2.Y, sd2.Z, sd2.Px, sd2.Py, sd2.Pz, sd2.TPCnCls, sd2.TPCchi2, sd2.ITSnCls, sd2.ITSchi2, sd2.ITSmeanClsSize, sd2.ITSmeanClsSizeL,
          sd2.Rigidity, sd2.TPCsignal, sd2.TPCnSigma, sd2.TPCnSigmaNhp, sd2.TPCnSigmaNlp, sd2.TOFmass, sd2.DcaXY, sd2.DcaZ, sd2.IsPvContributor,
          sd3.X, sd3.Y, sd3.Z, sd3.Px, sd3.Py, sd3.Pz, sd3.TPCnCls, sd3.TPCchi2, sd3.ITSnCls, sd3.ITSchi2, sd3.ITSmeanClsSize, sd3.ITSmeanClsSizeL,
          sd3.Rigidity, sd3.TPCsignal, sd3.TPCnSigma, sd3.TPCnSigmaNhp, sd3.TPCnSigmaNlp, sd3.TOFmass, sd3.DcaXY, sd3.DcaZ, sd3.IsPvContributor,
          sd1.SubMass, sd2.SubMass, sd3.SubMass);
    }
    if (cfgNprimDaughters == 3 && cfgNsecDaughters == 1) {
      const auto& d0 = cand.daughterTracks.at(0);
      const auto& d1 = cand.daughterTracks.at(1);
      const auto& d2 = cand.daughterTracks.at(2);
      const auto& sd1 = hypDaughter.daughterTracks.at(0);
      const auto& sd2 = hypDaughter.daughterTracks.at(1);
      if (!isMC || (isMC && cfgMCReconstructed && cand.IsReconstructed))
        outputTableTwoThree(
          cand.Species, cand.IsMatter, cand.PassedEvSel, cand.Mass, cand.Y, cand.Pt, cand.Ct, cand.CpaPv, cand.MaxDcaTracks, cand.DcaToPvXY, cand.DcaToPvZ, cand.DevToPvXY,
          cand.Chi2, cand.Pvx, cand.Pvy, cand.Pvz, cand.Svx, cand.Svy, cand.Svz, cand.Px, cand.Py, cand.Pz, cand.CollisionMcTrue, cand.McTrue, cand.McPhysicalPrimary,
          hypDaughter.Svx, hypDaughter.Svy, hypDaughter.Svz, d0.X, d0.Y, d0.Z, d0.Px, d0.Py, d0.Pz, hypDaughter.Mass, hypDaughter.Ct, hypDaughter.CpaPv,
          hypDaughter.MaxDcaTracks, hypDaughter.DcaToPvXY, hypDaughter.DcaToPvZ, hypDaughter.DcaToVtxXY, hypDaughter.DcaToVtxZ, hypDaughter.Chi2,
          d1.X, d1.Y, d1.Z, d1.Px, d1.Py, d1.Pz, d1.TPCnCls, d1.TPCchi2, d1.ITSnCls, d1.ITSchi2, d1.ITSmeanClsSize, d1.ITSmeanClsSizeL,
          d1.Rigidity, d1.TPCsignal, d1.TPCnSigma, d1.TPCnSigmaNhp, d1.TPCnSigmaNlp, d1.TOFmass, d1.DcaXY, d1.DcaZ, d1.IsPvContributor,
          d2.X, d2.Y, d2.Z, d2.Px, d2.Py, d2.Pz, d2.TPCnCls, d2.TPCchi2, d2.ITSnCls, d2.ITSchi2, d2.ITSmeanClsSize, d2.ITSmeanClsSizeL,
          d2.Rigidity, d2.TPCsignal, d2.TPCnSigma, d2.TPCnSigmaNhp, d2.TPCnSigmaNlp, d2.TOFmass, d2.DcaXY, d2.DcaZ, d2.IsPvContributor,
          d0.SubMass, d1.SubMass, d2.SubMass,
          sd1.X, sd1.Y, sd1.Z, sd1.Px, sd1.Py, sd1.Pz, sd1.TPCnCls, sd1.TPCchi2, sd1.ITSnCls, sd1.ITSchi2, sd1.ITSmeanClsSize, sd1.ITSmeanClsSizeL,
          sd1.Rigidity, sd1.TPCsignal, sd1.TPCnSigma, sd1.TPCnSigmaNhp, sd1.TPCnSigmaNlp, sd1.TOFmass, sd1.DcaXY, sd1.DcaZ, sd1.IsPvContributor,
          sd2.X, sd2.Y, sd2.Z, sd2.Px, sd2.Py, sd2.Pz, sd2.TPCnCls, sd2.TPCchi2, sd2.ITSnCls, sd2.ITSchi2, sd2.ITSmeanClsSize, sd2.ITSmeanClsSizeL,
          sd2.Rigidity, sd2.TPCsignal, sd2.TPCnSigma, sd2.TPCnSigmaNhp, sd2.TPCnSigmaNlp, sd2.TOFmass, sd2.DcaXY, sd2.DcaZ, sd2.IsPvContributor);
      if (isMC && cfgMCCombined)
        outputTableMcTwoThree(
          cand.PdgCode, cand.IsMatterMC, cand.IsReconstructed, cand.IsPhysicalPrimary, cand.PassedEvSelMC, cand.YGen, cand.PtGen, cand.CtGen,
          cand.CpaPvGen, cand.PxGen, cand.PyGen, cand.PzGen, cand.PvxGen, cand.PvyGen, cand.PvzGen, cand.SvxGen, cand.SvyGen, cand.SvzGen,
          cand.Species, cand.IsMatter, cand.PassedEvSel, cand.Mass, cand.Y, cand.Pt, cand.Ct, cand.CpaPv, cand.MaxDcaTracks, cand.DcaToPvXY, cand.DcaToPvZ, cand.DevToPvXY,
          cand.Chi2, cand.Pvx, cand.Pvy, cand.Pvz, cand.Svx, cand.Svy, cand.Svz, cand.Px, cand.Py, cand.Pz, cand.CollisionMcTrue,
          hypDaughter.Svx, hypDaughter.Svy, hypDaughter.Svz, d0.X, d0.Y, d0.Z, d0.Px, d0.Py, d0.Pz, hypDaughter.Mass, hypDaughter.Ct, hypDaughter.CpaPv,
          hypDaughter.MaxDcaTracks, hypDaughter.DcaToPvXY, hypDaughter.DcaToPvZ, hypDaughter.DcaToVtxXY, hypDaughter.DcaToVtxZ, hypDaughter.Chi2,
          d1.X, d1.Y, d1.Z, d1.Px, d1.Py, d1.Pz, d1.TPCnCls, d1.TPCchi2, d1.ITSnCls, d1.ITSchi2, d1.ITSmeanClsSize, d1.ITSmeanClsSizeL,
          d1.Rigidity, d1.TPCsignal, d1.TPCnSigma, d1.TPCnSigmaNhp, d1.TPCnSigmaNlp, d1.TOFmass, d1.DcaXY, d1.DcaZ, d1.IsPvContributor,
          d2.X, d2.Y, d2.Z, d2.Px, d2.Py, d2.Pz, d2.TPCnCls, d2.TPCchi2, d2.ITSnCls, d2.ITSchi2, d2.ITSmeanClsSize, d2.ITSmeanClsSizeL,
          d2.Rigidity, d2.TPCsignal, d2.TPCnSigma, d2.TPCnSigmaNhp, d2.TPCnSigmaNlp, d2.TOFmass, d2.DcaXY, d2.DcaZ, d2.IsPvContributor,
          d0.SubMass, d1.SubMass, d2.SubMass,
          sd1.X, sd1.Y, sd1.Z, sd1.Px, sd1.Py, sd1.Pz, sd1.TPCnCls, sd1.TPCchi2, sd1.ITSnCls, sd1.ITSchi2, sd1.ITSmeanClsSize, sd1.ITSmeanClsSizeL,
          sd1.Rigidity, sd1.TPCsignal, sd1.TPCnSigma, sd1.TPCnSigmaNhp, sd1.TPCnSigmaNlp, sd1.TOFmass, sd1.DcaXY, sd1.DcaZ, sd1.IsPvContributor,
          sd2.X, sd2.Y, sd2.Z, sd2.Px, sd2.Py, sd2.Pz, sd2.TPCnCls, sd2.TPCchi2, sd2.ITSnCls, sd2.ITSchi2, sd2.ITSmeanClsSize, sd2.ITSmeanClsSizeL,
          sd2.Rigidity, sd2.TPCsignal, sd2.TPCnSigma, sd2.TPCnSigmaNhp, sd2.TPCnSigmaNlp, sd2.TOFmass, sd2.DcaXY, sd2.DcaZ, sd2.IsPvContributor);
    }
  }
  //___________________________________________________________________________________________________________________________________________________________

  void fillCandidate(hyperNucleus& cand, hyperNucleus& /*hypDaughter*/, aod::HypKfHypNuc const& hypNuc, aod::HypKfHypNucs const&, aod::HypKfColls const&, aod::HypKfTracks const&, aod::HypKfDaughtAdds const&, aod::HypKfSubDs const&)
  {
    cand.daughterTracks.clear();
    cand.subDaughterMassVec.clear();
    auto coll = hypNuc.hypKfColl();
    auto addOns = hypNuc.addons_as<aod::HypKfDaughtAdds>();
    auto posVec = posVector(addOns);
    cand.Species = std::abs(hypNuc.species());
    cand.IsMatter = hypNuc.isMatter();
    cand.PassedEvSel = coll.passedEvSel();
    cand.Mass = hypNuc.mass();
    cand.Y = hypNuc.y();
    cand.Pt = hypNuc.pt();
    cand.Ct = ct(coll, hypNuc);
    cand.CpaPv = cpa(coll, hypNuc);
    cand.MaxDcaTracks = maxValue(dcaTrackSvAll(posVec, hypNuc, "XY"));
    cand.DcaToPvXY = hypNuc.dcaToPvXY();
    cand.DcaToPvZ = hypNuc.dcaToPvZ();
    cand.DcaToVtxXY = hypNuc.dcaToVtxXY();
    cand.DcaToVtxZ = hypNuc.dcaToVtxZ();
    cand.DevToPvXY = hypNuc.devToPvXY();
    cand.Chi2 = hypNuc.chi2();
    cand.Pvx = coll.posX();
    cand.Pvy = coll.posY();
    cand.Pvz = coll.posZ();
    cand.Svx = hypNuc.svx();
    cand.Svy = hypNuc.svy();
    cand.Svz = hypNuc.svz();
    cand.Px = hypNuc.px();
    cand.Py = hypNuc.py();
    cand.Pz = hypNuc.pz();
    if (cfgNsecDaughters) {
      trackProperties hypDaughter;
      cand.daughterTracks.push_back(hypDaughter);
    }
    auto daughterTracks = hypNuc.daughterTracks_as<aod::HypKfTracks>();
    for (auto& track : daughterTracks) {
      trackProperties daughter;
      daughter.TPCnCls = track.tpcNcluster();
      daughter.ITSnCls = track.itsNcluster();
      daughter.TPCchi2 = track.tpcChi2NCl();
      daughter.ITSchi2 = track.itsChi2NCl();
      daughter.ITSmeanClsSize = track.itsMeanClsSize();
      daughter.ITSmeanClsSizeL = track.itsMeanClsSize() * track.lambda();
      daughter.Rigidity = track.rigidity();
      daughter.TPCsignal = track.tpcSignal();
      daughter.TPCnSigma = track.tpcNsigma();
      daughter.TPCnSigmaNhp = track.tpcNsigmaNhp();
      daughter.TPCnSigmaNlp = track.tpcNsigmaNlp();
      daughter.TOFmass = track.tofMass();
      daughter.DcaXY = track.dcaXY();
      daughter.DcaZ = track.dcaZ();
      daughter.IsPvContributor = track.isPVContributor();
      cand.daughterTracks.push_back(daughter);
    }
    int trackCount = 0;
    for (auto& addOn : addOns) {
      cand.daughterTracks.at(trackCount).X = addOn.x();
      cand.daughterTracks.at(trackCount).Y = addOn.y();
      cand.daughterTracks.at(trackCount).Z = addOn.z();
      cand.daughterTracks.at(trackCount).Px = addOn.px();
      cand.daughterTracks.at(trackCount).Py = addOn.py();
      cand.daughterTracks.at(trackCount).Pz = addOn.py();
      trackCount++;
    }
    cand.NsingleDaughters = trackCount;
    if (cand.NsingleDaughters < 3)
      return;

    trackCount = 0;
    auto subDaughters = hypNuc.subDaughters_as<aod::HypKfSubDs>();
    for (auto& subDaughter : subDaughters) {
      cand.daughterTracks.at(trackCount++).SubMass = subDaughter.subMass();
    }
  }
  //___________________________________________________________________________________________________________________________________________________________

  void processMC(aod::HypKfMcParts const& mcHypNucs, aod::HypKfHypNucs const& hypNucs, aod::HypKfMcColls const&, aod::HypKfColls const& hypKfColls, aod::HypKfTracks const& hypKfTrks, aod::HypKfDaughtAdds const& hypKfDAdd, aod::HypKfSubDs const& hypKfDSub)
  {
    isMC = true;
    for (auto& mcHypNuc : mcHypNucs) {
      if (std::abs(mcHypNuc.species()) != cfgSpecies)
        continue;
      auto mcColl = mcHypNuc.hypKfMcColl();
      const auto mcParticleIdx = mcHypNuc.globalIndex();
      auto hypNucsByMc = hypNucs.sliceBy(perMcParticle, mcParticleIdx);
      hyperNucleus candidate, hypDaughter, dummy;
      candidate.PdgCode = mcHypNuc.pdgCode();
      candidate.IsMatterMC = mcHypNuc.isMatter();
      candidate.IsPhysicalPrimary = mcHypNuc.isPhysicalPrimary();
      candidate.PassedEvSelMC = mcColl.passedEvSel();
      candidate.YGen = mcHypNuc.y();
      candidate.PtGen = mcHypNuc.pt();
      candidate.CtGen = ct(mcColl, mcHypNuc);
      candidate.IsReconstructed = 0;
      candidate.CpaPvGen = cpa(mcColl, mcHypNuc);
      candidate.PxGen = mcHypNuc.px();
      candidate.PyGen = mcHypNuc.py();
      candidate.PzGen = mcHypNuc.pz();
      candidate.PvxGen = mcColl.posX();
      candidate.PvyGen = mcColl.posY();
      candidate.PvzGen = mcColl.posZ();
      candidate.SvxGen = mcHypNuc.svx();
      candidate.SvyGen = mcHypNuc.svy();
      candidate.SvzGen = mcHypNuc.svz();
      for (auto& hypNuc : hypNucsByMc) {
        auto coll = hypNuc.hypKfColl();
        if (coll.hypKfMcCollId() == mcHypNuc.hypKfMcCollId()) {
          candidate.CollisionMcTrue = true;
        }
        candidate.IsReconstructed++;
        fillCandidate(candidate, hypDaughter, hypNucs.rawIteratorAt(hypNuc.globalIndex()), hypNucs, hypKfColls, hypKfTrks, hypKfDAdd, hypKfDSub);
        if (cfgNsecDaughters) {
          fillCandidate(hypDaughter, dummy, hypNucs.rawIteratorAt(hypNuc.hypDaughterId()), hypNucs, hypKfColls, hypKfTrks, hypKfDAdd, hypKfDSub);
        }
      }
      fillTable(candidate, hypDaughter);
      hPt[0]->Fill(mcHypNuc.pt());
      if (candidate.IsReconstructed)
        hPt[1]->Fill(candidate.Pt);
    }
    hPt[2]->Divide(hPt[1].get(), hPt[0].get());
  }
  PROCESS_SWITCH(hypKfTreeCreator, processMC, "MC Gen tree", false);

  //___________________________________________________________________________________________________________________________________________________________
  std::vector<float> dcaTracksAll(std::vector<arr3>& posVec, TString opt = "")
  {
    std::vector<float> vec;
    int n = posVec.size();
    for (int i = 0; i < (n - 1); i++) {
      for (int j = (i + 1); j < n; j++) {
        vec.push_back(dcaTracks(posVec, i, j, opt));
      }
    }
    return vec;
  }
  template <class T>
  std::vector<float> dcaTrackSvAll(std::vector<arr3>& posVec, T const& hypNuc, TString opt = "")
  {
    std::vector<float> vec;
    for (size_t i = 0; i < posVec.size(); i++) {
      vec.push_back(dcaTrackSv(posVec, i, hypNuc, opt));
    }
    return vec;
  }

  float maxValue(std::vector<float> vec)
  {
    return *max_element(vec.begin(), vec.end());
  }
  float meanValue(std::vector<float> vec)
  {
    float sum = 0;
    for (auto value : vec)
      sum += value;
    return sum / vec.size();
  }
  float mean2Value(std::vector<float> vec)
  {
    float sum = 0;
    for (auto value : vec)
      sum += (value * value);
    return TMath::Sqrt(sum / vec.size());
  }

  float dcaTracks(std::vector<arr3> v, int track1, int track2, TString opt = "XY")
  {
    if (opt == "XY")
      return RecoDecay::distanceXY(v.at(track1), v.at(track2));
    else if (opt == "Z")
      return std::abs(v.at(track1).at(2) - v.at(track2).at(2));
    else
      return RecoDecay::distance(v.at(track1), v.at(track2));
  }
  template <class T>
  float dcaTrackSv(std::vector<arr3>& v, int track, T const& hypNuc, TString opt = "")
  {
    if (opt == "XY")
      return RecoDecay::distanceXY(v.at(track), decayVtx(hypNuc));
    else if (opt == "Z")
      return std::abs(v.at(track).at(2) - decayVtx(hypNuc).at(2));
    else
      return RecoDecay::distance(v.at(track), decayVtx(hypNuc));
  }
  template <class T>
  std::vector<arr3> posVector(T const& addons)
  {
    std::vector<arr3> v;
    for (auto& pos : addons) {
      v.push_back(std::array{pos.x(), pos.y(), pos.z()});
    }
    return v;
  }
  template <class T>
  arr3 primVtx(T const& coll)
  {
    return std::array{coll.posX(), coll.posY(), coll.posZ()};
  }
  template <class T>
  arr3 decayVtx(T const& hypNuc)
  {
    return std::array{hypNuc.svx(), hypNuc.svy(), hypNuc.svz()};
  }
  template <class T>
  arr3 momenta(T const& hypNuc)
  {
    return std::array{hypNuc.px(), hypNuc.py(), hypNuc.pz()};
  }
  template <class TColl, class TPart>
  float decayLength(TColl const& coll, TPart const& hypNuc)
  {
    return RecoDecay::distance(primVtx(coll), decayVtx(hypNuc));
  }
  template <class TColl, class TPart>
  float ct(TColl const& coll, TPart const& hypNuc)
  {
    return RecoDecay::ct(momenta(hypNuc), decayLength(coll, hypNuc), hypNuc.mass());
  }
  template <class TColl, class TPart>
  double cpa(TColl const& coll, TPart const& hypNuc)
  {
    return RecoDecay::cpa(primVtx(coll), decayVtx(hypNuc), momenta(hypNuc));
  }
  double paArcMin(double cosPa)
  {
    // returns the pointing angle (in arc min) for a given cosPa
    return TMath::ACos(cosPa) * 60 * 180 / TMath::Pi();
  }
  // only for Cascades
  template <class TPart>
  float decayLength(TPart const& mother, TPart const& daughter)
  {
    return RecoDecay::distance(decayVtx(mother), decayVtx(daughter));
  }
  template <class TPart>
  float ct(TPart const& mother, TPart const& daughter)
  {
    return RecoDecay::ct(momenta(daughter), decayLength(mother, daughter), daughter.mass());
  }
  template <class TPart>
  double cpa(TPart const& mother, TPart const& daughter)
  {
    return RecoDecay::cpa(decayVtx(mother), decayVtx(daughter), momenta(daughter));
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hypKfTreeCreator>(cfgc)};
}
