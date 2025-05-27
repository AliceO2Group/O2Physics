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
/// \file hypKfTreeCreator.cxx
/// \brief Creates flat tree for ML analysis
/// \author Janik Ditzel <jditzel@cern.ch> and Michael Hartung <mhartung@cern.ch>

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

struct TrackProperties {
  TrackProperties() : x(0), y(0), z(0), px(0), py(0), pz(0), tpcNcls(0), itsNcls(0), tpcChi2(0), itsChi2(0), itsMeanClsSize(0), itsMeanClsSizeL(0), rigidity(0), tpcSignal(0), tpcNsigma(0), tpcNsigmaNhp(0), tpcNsigmaNlp(0), tofMass(0), dcaXY(0), dcaZ(0), isPvContributor(0), subMass(0) {}
  float x, y, z, px, py, pz;
  uint8_t tpcNcls, itsNcls;
  float tpcChi2, itsChi2, itsMeanClsSize, itsMeanClsSizeL;
  float rigidity, tpcSignal, tpcNsigma, tpcNsigmaNhp, tpcNsigmaNlp;
  float tofMass, dcaXY, dcaZ;
  bool isPvContributor;
  float subMass;
};

struct HyperNucleus {
  HyperNucleus() : pdgCode(0), isReconstructed(0), globalIndex(0), species(0), isPrimaryCandidate(0), isMatter(0), passedEvSel(0), isMatterMC(0), passedEvSelMC(0), isPhysicalPrimary(0), collisionMcTrue(0), mass(0), y(0), pt(0), ct(0), yGen(0), ptGen(0), ctGen(0), cpaPvGen(0), cpaPv(0), cpaSv(0), maxDcaTracks(0), maxDcaTracksSv(0), dcaToPvXY(0), dcaToPvZ(0), dcaToVtxXY(0), dcaToVtxZ(0), devToPvXY(0), chi2(0), pvx(0), pvy(0), pvz(0), svx(0), svy(0), svz(0), px(0), py(0), pz(0), pvxGen(0), pvyGen(0), pvzGen(0), svxGen(0), svyGen(0), svzGen(0), pxGen(0), pyGen(0), pzGen(0), nSingleDaughters(0), nCascadeDaughters(0), mcTrue(0), mcTrueVtx(0), mcPhysicalPrimary(0), hypNucDaughter(0) {}
  int pdgCode, isReconstructed, globalIndex;
  uint8_t species;
  bool isPrimaryCandidate, isMatter, passedEvSel, isMatterMC, passedEvSelMC, isPhysicalPrimary, collisionMcTrue;
  float mass, y, pt, ct, yGen, ptGen, ctGen, cpaPvGen, cpaPv, cpaSv, maxDcaTracks, maxDcaTracksSv;
  float dcaToPvXY, dcaToPvZ, dcaToVtxXY, dcaToVtxZ, devToPvXY, chi2;
  float pvx, pvy, pvz, svx, svy, svz, px, py, pz;
  float pvxGen, pvyGen, pvzGen, svxGen, svyGen, svzGen, pxGen, pyGen, pzGen;
  int nSingleDaughters, nCascadeDaughters, cent, occu;
  bool mcTrue, mcTrueVtx, mcPhysicalPrimary;
  std::vector<TrackProperties> daughterTracks;
  std::vector<float> subDaughterMassVec;
  HyperNucleus* hypNucDaughter;
  ~HyperNucleus()
  {
    if (hypNucDaughter)
      delete hypNucDaughter;
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
DECLARE_SOA_COLUMN(Centrality, centrality, int);
DECLARE_SOA_COLUMN(Occupancy, occupancy, int);
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
DECLARE_SOA_COLUMN(D1d2Mass, d1d2Mass, float);
DECLARE_SOA_COLUMN(D1d3Mass, d1d3Mass, float);
DECLARE_SOA_COLUMN(D2d3Mass, d2d3Mass, float);
DECLARE_SOA_COLUMN(D3IsPvContributor, d3IsPvContributor, bool);
DECLARE_SOA_COLUMN(D0X, d0X, float);
DECLARE_SOA_COLUMN(D0Y, d0Y, float);
DECLARE_SOA_COLUMN(D0Z, d0Z, float);
DECLARE_SOA_COLUMN(D0Px, d0Px, float);
DECLARE_SOA_COLUMN(D0Py, d0Py, float);
DECLARE_SOA_COLUMN(D0Pz, d0Pz, float);
DECLARE_SOA_COLUMN(D0Mass, d0Mass, float);
DECLARE_SOA_COLUMN(D0ct, d0ct, float);
DECLARE_SOA_COLUMN(D0cosPa, d0cosPa, float);
DECLARE_SOA_COLUMN(D0dcaTracks, d0dcaTracks, float);
DECLARE_SOA_COLUMN(D0dcaTracksTv, d0dcaTracksTv, float);
DECLARE_SOA_COLUMN(D0dcaToPvXY, d0dcaToPvXY, float);
DECLARE_SOA_COLUMN(D0dcaToPvZ, d0dcaToPvZ, float);
DECLARE_SOA_COLUMN(D0dcaToSvXY, d0dcaToSvXY, float);
DECLARE_SOA_COLUMN(D0dcaToSvZ, d0dcaToSvZ, float);
DECLARE_SOA_COLUMN(D0chi2, d0chi2, float);
DECLARE_SOA_COLUMN(Sd1X, sd1X, float);
DECLARE_SOA_COLUMN(Sd1Y, sd1Y, float);
DECLARE_SOA_COLUMN(Sd1Z, sd1Z, float);
DECLARE_SOA_COLUMN(Sd1Px, sd1Px, float);
DECLARE_SOA_COLUMN(Sd1Py, sd1Py, float);
DECLARE_SOA_COLUMN(Sd1Pz, sd1Pz, float);
DECLARE_SOA_COLUMN(Sd1TPCnCls, sd1TPCnCls, uint8_t);
DECLARE_SOA_COLUMN(Sd1TPCchi2, sd1TPCchi2, float);
DECLARE_SOA_COLUMN(Sd1ITSnCls, sd1ITSnCls, uint8_t);
DECLARE_SOA_COLUMN(Sd1ITSchi2, sd1ITSchi2, float);
DECLARE_SOA_COLUMN(Sd1ITSmeanClsSize, sd1ITSmeanClsSize, float);
DECLARE_SOA_COLUMN(Sd1ITSmeanClsSizeL, sd1ITSmeanClsSizeL, float);
DECLARE_SOA_COLUMN(Sd1Rigidity, sd1Rigidity, float);
DECLARE_SOA_COLUMN(Sd1TPCsignal, sd1TPCsignal, float);
DECLARE_SOA_COLUMN(Sd1TPCnSigma, sd1TPCnSigma, float);
DECLARE_SOA_COLUMN(Sd1TPCnSigmaNhp, sd1TPCnSigmaNhp, float);
DECLARE_SOA_COLUMN(Sd1TPCnSigmaNlp, sd1TPCnSigmaNlp, float);
DECLARE_SOA_COLUMN(Sd1TOFmass, sd1TOFmass, float);
DECLARE_SOA_COLUMN(Sd1DcaXY, sd1DcaXY, float);
DECLARE_SOA_COLUMN(Sd1DcaZ, sd1DcaZ, float);
DECLARE_SOA_COLUMN(Sd1IsPvContributor, sd1IsPvContributor, bool);
DECLARE_SOA_COLUMN(Sd2X, sd2X, float);
DECLARE_SOA_COLUMN(Sd2Y, sd2Y, float);
DECLARE_SOA_COLUMN(Sd2Z, sd2Z, float);
DECLARE_SOA_COLUMN(Sd2Px, sd2Px, float);
DECLARE_SOA_COLUMN(Sd2Py, sd2Py, float);
DECLARE_SOA_COLUMN(Sd2Pz, sd2Pz, float);
DECLARE_SOA_COLUMN(Sd2TPCnCls, sd2TPCnCls, uint8_t);
DECLARE_SOA_COLUMN(Sd2TPCchi2, sd2TPCchi2, float);
DECLARE_SOA_COLUMN(Sd2ITSnCls, sd2ITSnCls, uint8_t);
DECLARE_SOA_COLUMN(Sd2ITSchi2, sd2ITSchi2, float);
DECLARE_SOA_COLUMN(Sd2ITSmeanClsSize, sd2ITSmeanClsSize, float);
DECLARE_SOA_COLUMN(Sd2ITSmeanClsSizeL, sd2ITSmeanClsSizeL, float);
DECLARE_SOA_COLUMN(Sd2Rigidity, sd2Rigidity, float);
DECLARE_SOA_COLUMN(Sd2TPCsignal, sd2TPCsignal, float);
DECLARE_SOA_COLUMN(Sd2TPCnSigma, sd2TPCnSigma, float);
DECLARE_SOA_COLUMN(Sd2TPCnSigmaNhp, sd2TPCnSigmaNhp, float);
DECLARE_SOA_COLUMN(Sd2TPCnSigmaNlp, sd2TPCnSigmaNlp, float);
DECLARE_SOA_COLUMN(Sd2TOFmass, sd2TOFmass, float);
DECLARE_SOA_COLUMN(Sd2DcaXY, sd2DcaXY, float);
DECLARE_SOA_COLUMN(Sd2DcaZ, sd2DcaZ, float);
DECLARE_SOA_COLUMN(Sd2IsPvContributor, sd2IsPvContributor, bool);
DECLARE_SOA_COLUMN(Sd3X, sd3X, float);
DECLARE_SOA_COLUMN(Sd3Y, sd3Y, float);
DECLARE_SOA_COLUMN(Sd3Z, sd3Z, float);
DECLARE_SOA_COLUMN(Sd3Px, sd3Px, float);
DECLARE_SOA_COLUMN(Sd3Py, sd3Py, float);
DECLARE_SOA_COLUMN(Sd3Pz, sd3Pz, float);
DECLARE_SOA_COLUMN(Sd3TPCnCls, sd3TPCnCls, uint8_t);
DECLARE_SOA_COLUMN(Sd3TPCchi2, sd3TPCchi2, float);
DECLARE_SOA_COLUMN(Sd3ITSnCls, sd3ITSnCls, uint8_t);
DECLARE_SOA_COLUMN(Sd3ITSchi2, sd3ITSchi2, float);
DECLARE_SOA_COLUMN(Sd3ITSmeanClsSize, sd3ITSmeanClsSize, float);
DECLARE_SOA_COLUMN(Sd3ITSmeanClsSizeL, sd3ITSmeanClsSizeL, float);
DECLARE_SOA_COLUMN(Sd3Rigidity, sd3Rigidity, float);
DECLARE_SOA_COLUMN(Sd3TPCsignal, sd3TPCsignal, float);
DECLARE_SOA_COLUMN(Sd3TPCnSigma, sd3TPCnSigma, float);
DECLARE_SOA_COLUMN(Sd3TPCnSigmaNhp, sd3TPCnSigmaNhp, float);
DECLARE_SOA_COLUMN(Sd3TPCnSigmaNlp, sd3TPCnSigmaNlp, float);
DECLARE_SOA_COLUMN(Sd3TOFmass, sd3TOFmass, float);
DECLARE_SOA_COLUMN(Sd3DcaXY, sd3DcaXY, float);
DECLARE_SOA_COLUMN(Sd3DcaZ, sd3DcaZ, float);
DECLARE_SOA_COLUMN(Sd3IsPvContributor, sd3IsPvContributor, bool);
DECLARE_SOA_COLUMN(Sd1sd2Mass, sd1sd2Mass, float);
DECLARE_SOA_COLUMN(Sd1sd3Mass, sd1sd3Mass, float);
DECLARE_SOA_COLUMN(Sd2sd3Mass, sd2sd3Mass, float);
} // namespace hypkftree

#define HYPKFGENBASE mcparticle::PdgCode, hypkftree::IsMatterGen, hypkftree::IsReconstructed, hykfmc::IsPhysicalPrimary, hypkftree::PassedEvSelMC, hypkftree::YGen, hypkftree::PtGen, hypkftree::CtGen

#define HYPKFGENEXT hypkftree::CpaPvGen, hypkftree::PxGen, hypkftree::PyGen, hypkftree::PzGen, hypkftree::PvxGen, hypkftree::PvyGen, hypkftree::PvzGen, hypkftree::SvxGen, hypkftree::SvyGen, hypkftree::SvzGen

#define HYPKFGENCAS hypkftree::TvxGen, hypkftree::TvyGen, hypkftree::TvzGen

#define HYPKFHYPNUC hykfmc::Species, hypkftree::IsMatter, hypkftree::Centrality, hypkftree::Occupancy, hykfmccoll::PassedEvSel, hykfhyp::Mass, hypkftree::Y, track::Pt, hypkftree::Ct, hypkftree::CosPa, hypkftree::DcaTracks, hypkftree::DcaTrackSv, hykfhyp::DcaToPvXY, hykfhyp::DcaToPvZ, hykfhyp::DevToPvXY, hykfhyp::Chi2, hypkftree::Pvx, hypkftree::Pvy, hypkftree::Pvz, hykfmc::Svx, hykfmc::Svy, hykfmc::Svz, hykfhyp::Px, hykfhyp::Py, hykfhyp::Pz, hypkftree::CollMcTrue

#define HYPKFHYPNUCMC hypkftree::McTrue, hykfmc::IsPhysicalPrimary

#define HYPKFD0 hypkftree::Tvx, hypkftree::Tvy, hypkftree::Tvz, hypkftree::D0X, hypkftree::D0Y, hypkftree::D0Z, hypkftree::D0Px, hypkftree::D0Py, hypkftree::D0Pz, hypkftree::D0Mass, hypkftree::D0ct, hypkftree::D0cosPa, hypkftree::D0dcaTracks, hypkftree::D0dcaToPvXY, hypkftree::D0dcaToPvZ, hypkftree::D0dcaToSvXY, hypkftree::D0dcaToSvZ, hypkftree::D0chi2

#define HYPKFD1 hypkftree::D1X, hypkftree::D1Y, hypkftree::D1Z, hypkftree::D1Px, hypkftree::D1Py, hypkftree::D1Pz, hypkftree::D1TPCnCls, hypkftree::D1TPCchi2, hypkftree::D1ITSnCls, hypkftree::D1ITSchi2, hypkftree::D1ITSmeanClsSizeL, hypkftree::D1Rigidity, hypkftree::D1TPCsignal, hypkftree::D1TPCnSigma, hypkftree::D1TPCnSigmaNhp, hypkftree::D1TPCnSigmaNlp, hypkftree::D1TOFmass, hypkftree::D1DcaXY, hypkftree::D1DcaZ, hypkftree::D1IsPvContributor

#define HYPKFD2 hypkftree::D2X, hypkftree::D2Y, hypkftree::D2Z, hypkftree::D2Px, hypkftree::D2Py, hypkftree::D2Pz, hypkftree::D2TPCnCls, hypkftree::D2TPCchi2, hypkftree::D2ITSnCls, hypkftree::D2ITSchi2, hypkftree::D2ITSmeanClsSizeL, hypkftree::D2Rigidity, hypkftree::D2TPCsignal, hypkftree::D2TPCnSigma, hypkftree::D2TPCnSigmaNhp, hypkftree::D2TPCnSigmaNlp, hypkftree::D2TOFmass, hypkftree::D2DcaXY, hypkftree::D2DcaZ, hypkftree::D2IsPvContributor

#define HYPKFD3 hypkftree::D3X, hypkftree::D3Y, hypkftree::D3Z, hypkftree::D3Px, hypkftree::D3Py, hypkftree::D3Pz, hypkftree::D3TPCnCls, hypkftree::D3TPCchi2, hypkftree::D3ITSnCls, hypkftree::D3ITSchi2, hypkftree::D3ITSmeanClsSizeL, hypkftree::D3Rigidity, hypkftree::D3TPCsignal, hypkftree::D3TPCnSigma, hypkftree::D3TPCnSigmaNhp, hypkftree::D3TPCnSigmaNlp, hypkftree::D3TOFmass, hypkftree::D3DcaXY, hypkftree::D3DcaZ, hypkftree::D3IsPvContributor

#define HYPKFSD1 hypkftree::Sd1X, hypkftree::Sd1Y, hypkftree::Sd1Z, hypkftree::Sd1Px, hypkftree::Sd1Py, hypkftree::Sd1Pz, hypkftree::Sd1TPCnCls, hypkftree::Sd1TPCchi2, hypkftree::Sd1ITSnCls, hypkftree::Sd1ITSchi2, hypkftree::Sd1ITSmeanClsSizeL, hypkftree::Sd1Rigidity, hypkftree::Sd1TPCsignal, hypkftree::Sd1TPCnSigma, hypkftree::Sd1TPCnSigmaNhp, hypkftree::Sd1TPCnSigmaNlp, hypkftree::Sd1TOFmass, hypkftree::Sd1DcaXY, hypkftree::Sd1DcaZ, hypkftree::Sd1IsPvContributor

#define HYPKFSD2 hypkftree::Sd2X, hypkftree::Sd2Y, hypkftree::Sd2Z, hypkftree::Sd2Px, hypkftree::Sd2Py, hypkftree::Sd2Pz, hypkftree::Sd2TPCnCls, hypkftree::Sd2TPCchi2, hypkftree::Sd2ITSnCls, hypkftree::Sd2ITSchi2, hypkftree::Sd2ITSmeanClsSizeL, hypkftree::Sd2Rigidity, hypkftree::Sd2TPCsignal, hypkftree::Sd2TPCnSigma, hypkftree::Sd2TPCnSigmaNhp, hypkftree::Sd2TPCnSigmaNlp, hypkftree::Sd2TOFmass, hypkftree::Sd2DcaXY, hypkftree::Sd2DcaZ, hypkftree::Sd2IsPvContributor

#define HYPKFSD3 hypkftree::Sd3X, hypkftree::Sd3Y, hypkftree::Sd3Z, hypkftree::Sd3Px, hypkftree::Sd3Py, hypkftree::Sd3Pz, hypkftree::Sd3TPCnCls, hypkftree::Sd3TPCchi2, hypkftree::Sd3ITSnCls, hypkftree::Sd3ITSchi2, hypkftree::Sd3ITSmeanClsSizeL, hypkftree::Sd3Rigidity, hypkftree::Sd3TPCsignal, hypkftree::Sd3TPCnSigma, hypkftree::Sd3TPCnSigmaNhp, hypkftree::Sd3TPCnSigmaNlp, hypkftree::Sd3TOFmass, hypkftree::Sd3DcaXY, hypkftree::Sd3DcaZ, hypkftree::Sd3IsPvContributor

#define HYPKFSDMASS hypkftree::D1d2Mass, hypkftree::D1d3Mass, hypkftree::D2d3Mass
#define HYPKFSSDMASS hypkftree::Sd1sd2Mass, hypkftree::Sd1sd3Mass, hypkftree::Sd2sd3Mass

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

struct HypKfTreeCreator {

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
    for (const auto& hypNuc : hypNucs) {
      if (cfgSpecies && std::abs(hypNuc.species()) != cfgSpecies)
        continue;
      HyperNucleus candidate, hypNucDaughter;
      fillCandidatePrim(candidate, hypNuc, hypNucs, hypKfColls, hypKfTrks, hypKfDAdd, hypKfDSub);
      if (hypNuc.hypDaughterId() >= 0) {
        fillCandidateSec(hypNucDaughter, hypNucs.rawIteratorAt(hypNuc.hypDaughterId()), hypNuc, hypNucs, hypKfColls, hypKfTrks, hypKfDAdd, hypKfDSub);
      }
      fillTable(candidate, hypNucDaughter);
    }
  }
  PROCESS_SWITCH(HypKfTreeCreator, processData, "single tree", false);
  //___________________________________________________________________________________________________________________________________________________________
  void fillTable(HyperNucleus& cand, HyperNucleus& hypDaughter)
  {
    if (isMC && cfgMCGenerated)
      outputMcGenTable(
        cand.pdgCode, cand.isMatterMC, cand.isReconstructed, cand.isPhysicalPrimary, cand.passedEvSelMC, cand.yGen, cand.ptGen, cand.ctGen);

    if (!cand.isReconstructed) {
      cand.daughterTracks.resize(4);
      cand.subDaughterMassVec.resize(4);
      hypDaughter.daughterTracks.resize(4);
      hypDaughter.subDaughterMassVec.resize(4);
    }

    if (cfgNprimDaughters == 2 && cfgNsecDaughters == 0) {
      const auto& d1 = cand.daughterTracks.at(0);
      const auto& d2 = cand.daughterTracks.at(1);
      if (!isMC || (isMC && cfgMCReconstructed && cand.isReconstructed))
        outputTableTwo(
          cand.species, cand.isMatter, cand.cent, cand.occu, cand.passedEvSel, cand.mass, cand.y, cand.pt, cand.ct, cand.cpaPv, cand.maxDcaTracks, cand.maxDcaTracksSv, cand.dcaToPvXY,
          cand.dcaToPvZ, cand.devToPvXY, cand.chi2, cand.pvx, cand.pvy, cand.pvz, cand.svx, cand.svy, cand.svz, cand.px, cand.py, cand.pz, cand.collisionMcTrue,
          cand.mcTrue, cand.mcPhysicalPrimary,
          d1.x, d1.y, d1.z, d1.px, d1.py, d1.pz, d1.tpcNcls, d1.tpcChi2, d1.itsNcls, d1.itsChi2, d1.itsMeanClsSizeL,
          d1.rigidity, d1.tpcSignal, d1.tpcNsigma, d1.tpcNsigmaNhp, d1.tpcNsigmaNlp, d1.tofMass, d1.dcaXY, d1.dcaZ, d1.isPvContributor,
          d2.x, d2.y, d2.z, d2.px, d2.py, d2.pz, d2.tpcNcls, d2.tpcChi2, d2.itsNcls, d2.itsChi2, d2.itsMeanClsSizeL,
          d2.rigidity, d2.tpcSignal, d2.tpcNsigma, d2.tpcNsigmaNhp, d2.tpcNsigmaNlp, d2.tofMass, d2.dcaXY, d2.dcaZ, d2.isPvContributor);
      if (isMC && cfgMCCombined)
        outputTableMcTwo(
          cand.pdgCode, cand.isMatterMC, cand.isReconstructed, cand.isPhysicalPrimary, cand.passedEvSelMC, cand.yGen, cand.ptGen, cand.ctGen,
          cand.cpaPvGen, cand.pxGen, cand.pyGen, cand.pzGen, cand.pvxGen, cand.pvyGen, cand.pvzGen, cand.svxGen, cand.svyGen, cand.svzGen,
          cand.species, cand.isMatter, cand.cent, cand.occu, cand.passedEvSel, cand.mass, cand.y, cand.pt, cand.ct, cand.cpaPv, cand.maxDcaTracks, cand.maxDcaTracksSv,
          cand.dcaToPvXY, cand.dcaToPvZ, cand.devToPvXY,
          cand.chi2, cand.pvx, cand.pvy, cand.pvz, cand.svx, cand.svy, cand.svz, cand.px, cand.py, cand.pz, cand.collisionMcTrue,
          d1.x, d1.y, d1.z, d1.px, d1.py, d1.pz, d1.tpcNcls, d1.tpcChi2, d1.itsNcls, d1.itsChi2, d1.itsMeanClsSizeL,
          d1.rigidity, d1.tpcSignal, d1.tpcNsigma, d1.tpcNsigmaNhp, d1.tpcNsigmaNlp, d1.tofMass, d1.dcaXY, d1.dcaZ, d1.isPvContributor,
          d2.x, d2.y, d2.z, d2.px, d2.py, d2.pz, d2.tpcNcls, d2.tpcChi2, d2.itsNcls, d2.itsChi2, d2.itsMeanClsSizeL,
          d2.rigidity, d2.tpcSignal, d2.tpcNsigma, d2.tpcNsigmaNhp, d2.tpcNsigmaNlp, d2.tofMass, d2.dcaXY, d2.dcaZ, d2.isPvContributor);
    }
    if (cand.isPrimaryCandidate && ((cfgNprimDaughters == 3 && cfgNsecDaughters == 0) || (cfgNsecDaughters == 3 && cfgSpecies == 0))) {
      const auto& d1 = cand.daughterTracks.at(0);
      const auto& d2 = cand.daughterTracks.at(1);
      const auto& d3 = cand.daughterTracks.at(2);
      if (!isMC || (isMC && cfgMCReconstructed && cand.isReconstructed))
        outputTableThree(
          cand.species, cand.isMatter, cand.cent, cand.occu, cand.passedEvSel, cand.mass, cand.y, cand.pt, cand.ct, cand.cpaPv, cand.maxDcaTracks, cand.maxDcaTracksSv, cand.dcaToPvXY,
          cand.dcaToPvZ, cand.devToPvXY, cand.chi2, cand.pvx, cand.pvy, cand.pvz, cand.svx, cand.svy, cand.svz, cand.px, cand.py, cand.pz, cand.collisionMcTrue,
          cand.mcTrue, cand.mcPhysicalPrimary,
          d1.x, d1.y, d1.z, d1.px, d1.py, d1.pz, d1.tpcNcls, d1.tpcChi2, d1.itsNcls, d1.itsChi2, d1.itsMeanClsSizeL,
          d1.rigidity, d1.tpcSignal, d1.tpcNsigma, d1.tpcNsigmaNhp, d1.tpcNsigmaNlp, d1.tofMass, d1.dcaXY, d1.dcaZ, d1.isPvContributor,
          d2.x, d2.y, d2.z, d2.px, d2.py, d2.pz, d2.tpcNcls, d2.tpcChi2, d2.itsNcls, d2.itsChi2, d2.itsMeanClsSizeL,
          d2.rigidity, d2.tpcSignal, d2.tpcNsigma, d2.tpcNsigmaNhp, d2.tpcNsigmaNlp, d2.tofMass, d2.dcaXY, d2.dcaZ, d2.isPvContributor,
          d3.x, d3.y, d3.z, d3.px, d3.py, d3.pz, d3.tpcNcls, d3.tpcChi2, d3.itsNcls, d3.itsChi2, d3.itsMeanClsSizeL,
          d3.rigidity, d3.tpcSignal, d3.tpcNsigma, d3.tpcNsigmaNhp, d3.tpcNsigmaNlp, d3.tofMass, d3.dcaXY, d3.dcaZ, d3.isPvContributor,
          d1.subMass, d2.subMass, d3.subMass);
      if (isMC && cfgMCCombined)
        outputTableMcThree(
          cand.pdgCode, cand.isMatterMC, cand.isReconstructed, cand.isPhysicalPrimary, cand.passedEvSelMC, cand.yGen, cand.ptGen, cand.ctGen,
          cand.cpaPvGen, cand.pxGen, cand.pyGen, cand.pzGen, cand.pvxGen, cand.pvyGen, cand.pvzGen, cand.svxGen, cand.svyGen, cand.svzGen,
          cand.species, cand.isMatter, cand.cent, cand.occu, cand.passedEvSel, cand.mass, cand.y, cand.pt, cand.ct, cand.cpaPv, cand.maxDcaTracks, cand.maxDcaTracksSv,
          cand.dcaToPvXY, cand.dcaToPvZ, cand.devToPvXY, cand.chi2, cand.pvx, cand.pvy, cand.pvz, cand.svx, cand.svy, cand.svz, cand.px, cand.py,
          cand.pz, cand.collisionMcTrue,
          d1.x, d1.y, d1.z, d1.px, d1.py, d1.pz, d1.tpcNcls, d1.tpcChi2, d1.itsNcls, d1.itsChi2, d1.itsMeanClsSizeL,
          d1.rigidity, d1.tpcSignal, d1.tpcNsigma, d1.tpcNsigmaNhp, d1.tpcNsigmaNlp, d1.tofMass, d1.dcaXY, d1.dcaZ, d1.isPvContributor,
          d2.x, d2.y, d2.z, d2.px, d2.py, d2.pz, d2.tpcNcls, d2.tpcChi2, d2.itsNcls, d2.itsChi2, d2.itsMeanClsSizeL,
          d2.rigidity, d2.tpcSignal, d2.tpcNsigma, d2.tpcNsigmaNhp, d2.tpcNsigmaNlp, d2.tofMass, d2.dcaXY, d2.dcaZ, d2.isPvContributor,
          d3.x, d3.y, d3.z, d3.px, d3.py, d3.pz, d3.tpcNcls, d3.tpcChi2, d3.itsNcls, d3.itsChi2, d3.itsMeanClsSizeL,
          d3.rigidity, d3.tpcSignal, d3.tpcNsigma, d3.tpcNsigmaNhp, d3.tpcNsigmaNlp, d3.tofMass, d3.dcaXY, d3.dcaZ, d3.isPvContributor,
          d1.subMass, d2.subMass, d3.subMass);
    }
    if (cfgNprimDaughters == 2 && cfgNsecDaughters == 3) {
      const auto& d0 = cand.daughterTracks.at(0);
      const auto& d1 = cand.daughterTracks.at(1);
      const auto& sd1 = hypDaughter.daughterTracks.at(0);
      const auto& sd2 = hypDaughter.daughterTracks.at(1);
      const auto& sd3 = hypDaughter.daughterTracks.at(2);
      if (!isMC || (isMC && cfgMCReconstructed && cand.isReconstructed))
        outputTableTwoThree(
          cand.species, cand.isMatter, cand.cent, cand.occu, cand.passedEvSel, cand.mass, cand.y, cand.pt, cand.ct, cand.cpaPv, cand.maxDcaTracks, cand.maxDcaTracksSv, cand.dcaToPvXY,
          cand.dcaToPvZ, cand.devToPvXY, cand.chi2, cand.pvx, cand.pvy, cand.pvz, cand.svx, cand.svy, cand.svz, cand.px, cand.py, cand.pz, cand.collisionMcTrue,
          cand.mcTrue, cand.mcPhysicalPrimary,
          hypDaughter.svx, hypDaughter.svy, hypDaughter.svz, d0.x, d0.y, d0.z, d0.px, d0.py, d0.pz, hypDaughter.mass, hypDaughter.ct, hypDaughter.cpaPv,
          hypDaughter.maxDcaTracks, hypDaughter.dcaToPvXY, hypDaughter.dcaToPvZ, hypDaughter.dcaToVtxXY, hypDaughter.dcaToVtxZ, hypDaughter.chi2,
          d1.x, d1.y, d1.z, d1.px, d1.py, d1.pz, d1.tpcNcls, d1.tpcChi2, d1.itsNcls, d1.itsChi2, d1.itsMeanClsSizeL,
          d1.rigidity, d1.tpcSignal, d1.tpcNsigma, d1.tpcNsigmaNhp, d1.tpcNsigmaNlp, d1.tofMass, d1.dcaXY, d1.dcaZ, d1.isPvContributor,
          sd1.x, sd1.y, sd1.z, sd1.px, sd1.py, sd1.pz, sd1.tpcNcls, sd1.tpcChi2, sd1.itsNcls, sd1.itsChi2, sd1.itsMeanClsSizeL,
          sd1.rigidity, sd1.tpcSignal, sd1.tpcNsigma, sd1.tpcNsigmaNhp, sd1.tpcNsigmaNlp, sd1.tofMass, sd1.dcaXY, sd1.dcaZ, sd1.isPvContributor,
          sd2.x, sd2.y, sd2.z, sd2.px, sd2.py, sd2.pz, sd2.tpcNcls, sd2.tpcChi2, sd2.itsNcls, sd2.itsChi2, sd2.itsMeanClsSizeL,
          sd2.rigidity, sd2.tpcSignal, sd2.tpcNsigma, sd2.tpcNsigmaNhp, sd2.tpcNsigmaNlp, sd2.tofMass, sd2.dcaXY, sd2.dcaZ, sd2.isPvContributor,
          sd3.x, sd3.y, sd3.z, sd3.px, sd3.py, sd3.pz, sd3.tpcNcls, sd3.tpcChi2, sd3.itsNcls, sd3.itsChi2, sd3.itsMeanClsSizeL,
          sd3.rigidity, sd3.tpcSignal, sd3.tpcNsigma, sd3.tpcNsigmaNhp, sd3.tpcNsigmaNlp, sd3.tofMass, sd3.dcaXY, sd3.dcaZ, sd3.isPvContributor,
          sd1.subMass, sd2.subMass, sd3.subMass);
      if (isMC && cfgMCCombined)
        outputTableMcTwoThree(
          cand.pdgCode, cand.isMatterMC, cand.isReconstructed, cand.isPhysicalPrimary, cand.passedEvSelMC, cand.yGen, cand.ptGen, cand.ctGen,
          cand.cpaPvGen, cand.pxGen, cand.pyGen, cand.pzGen, cand.pvxGen, cand.pvyGen, cand.pvzGen, cand.svxGen, cand.svyGen, cand.svzGen,
          cand.species, cand.isMatter, cand.cent, cand.occu, cand.passedEvSel, cand.mass, cand.y, cand.pt, cand.ct, cand.cpaPv, cand.maxDcaTracks, cand.maxDcaTracksSv, cand.dcaToPvXY,
          cand.dcaToPvZ, cand.devToPvXY,
          cand.chi2, cand.pvx, cand.pvy, cand.pvz, cand.svx, cand.svy, cand.svz, cand.px, cand.py, cand.pz, cand.collisionMcTrue,
          hypDaughter.svx, hypDaughter.svy, hypDaughter.svz, d0.x, d0.y, d0.z, d0.px, d0.py, d0.pz, hypDaughter.mass, hypDaughter.ct, hypDaughter.cpaPv,
          hypDaughter.maxDcaTracks, hypDaughter.dcaToPvXY, hypDaughter.dcaToPvZ, hypDaughter.dcaToVtxXY, hypDaughter.dcaToVtxZ, hypDaughter.chi2,
          d1.x, d1.y, d1.z, d1.px, d1.py, d1.pz, d1.tpcNcls, d1.tpcChi2, d1.itsNcls, d1.itsChi2, d1.itsMeanClsSizeL,
          d1.rigidity, d1.tpcSignal, d1.tpcNsigma, d1.tpcNsigmaNhp, d1.tpcNsigmaNlp, d1.tofMass, d1.dcaXY, d1.dcaZ, d1.isPvContributor,
          sd1.x, sd1.y, sd1.z, sd1.px, sd1.py, sd1.pz, sd1.tpcNcls, sd1.tpcChi2, sd1.itsNcls, sd1.itsChi2, sd1.itsMeanClsSizeL,
          sd1.rigidity, sd1.tpcSignal, sd1.tpcNsigma, sd1.tpcNsigmaNhp, sd1.tpcNsigmaNlp, sd1.tofMass, sd1.dcaXY, sd1.dcaZ, sd1.isPvContributor,
          sd2.x, sd2.y, sd2.z, sd2.px, sd2.py, sd2.pz, sd2.tpcNcls, sd2.tpcChi2, sd2.itsNcls, sd2.itsChi2, sd2.itsMeanClsSizeL,
          sd2.rigidity, sd2.tpcSignal, sd2.tpcNsigma, sd2.tpcNsigmaNhp, sd2.tpcNsigmaNlp, sd2.tofMass, sd2.dcaXY, sd2.dcaZ, sd2.isPvContributor,
          sd3.x, sd3.y, sd3.z, sd3.px, sd3.py, sd3.pz, sd3.tpcNcls, sd3.tpcChi2, sd3.itsNcls, sd3.itsChi2, sd3.itsMeanClsSizeL,
          sd3.rigidity, sd3.tpcSignal, sd3.tpcNsigma, sd3.tpcNsigmaNhp, sd3.tpcNsigmaNlp, sd3.tofMass, sd3.dcaXY, sd3.dcaZ, sd3.isPvContributor,
          sd1.subMass, sd2.subMass, sd3.subMass);
    }
    if (cfgNprimDaughters == 3 && cfgNsecDaughters == 1) {
      const auto& d0 = cand.daughterTracks.at(0);
      const auto& d1 = cand.daughterTracks.at(1);
      const auto& d2 = cand.daughterTracks.at(2);
      const auto& sd1 = hypDaughter.daughterTracks.at(0);
      const auto& sd2 = hypDaughter.daughterTracks.at(1);
      if (!isMC || (isMC && cfgMCReconstructed && cand.isReconstructed))
        outputTableThreeTwo(
          cand.species, cand.isMatter, cand.cent, cand.occu, cand.passedEvSel, cand.mass, cand.y, cand.pt, cand.ct, cand.cpaPv, cand.maxDcaTracks, cand.maxDcaTracksSv, cand.dcaToPvXY,
          cand.dcaToPvZ, cand.devToPvXY, cand.chi2, cand.pvx, cand.pvy, cand.pvz, cand.svx, cand.svy, cand.svz, cand.px, cand.py, cand.pz, cand.collisionMcTrue,
          cand.mcTrue, cand.mcPhysicalPrimary, hypDaughter.svx, hypDaughter.svy, hypDaughter.svz, d0.x, d0.y, d0.z, d0.px, d0.py, d0.pz, hypDaughter.mass, hypDaughter.ct,
          hypDaughter.cpaPv, hypDaughter.maxDcaTracks, hypDaughter.dcaToPvXY, hypDaughter.dcaToPvZ, hypDaughter.dcaToVtxXY, hypDaughter.dcaToVtxZ, hypDaughter.chi2,
          d1.x, d1.y, d1.z, d1.px, d1.py, d1.pz, d1.tpcNcls, d1.tpcChi2, d1.itsNcls, d1.itsChi2, d1.itsMeanClsSizeL,
          d1.rigidity, d1.tpcSignal, d1.tpcNsigma, d1.tpcNsigmaNhp, d1.tpcNsigmaNlp, d1.tofMass, d1.dcaXY, d1.dcaZ, d1.isPvContributor,
          d2.x, d2.y, d2.z, d2.px, d2.py, d2.pz, d2.tpcNcls, d2.tpcChi2, d2.itsNcls, d2.itsChi2, d2.itsMeanClsSizeL,
          d2.rigidity, d2.tpcSignal, d2.tpcNsigma, d2.tpcNsigmaNhp, d2.tpcNsigmaNlp, d2.tofMass, d2.dcaXY, d2.dcaZ, d2.isPvContributor,
          d0.subMass, d1.subMass, d2.subMass,
          sd1.x, sd1.y, sd1.z, sd1.px, sd1.py, sd1.pz, sd1.tpcNcls, sd1.tpcChi2, sd1.itsNcls, sd1.itsChi2, sd1.itsMeanClsSizeL,
          sd1.rigidity, sd1.tpcSignal, sd1.tpcNsigma, sd1.tpcNsigmaNhp, sd1.tpcNsigmaNlp, sd1.tofMass, sd1.dcaXY, sd1.dcaZ, sd1.isPvContributor,
          sd2.x, sd2.y, sd2.z, sd2.px, sd2.py, sd2.pz, sd2.tpcNcls, sd2.tpcChi2, sd2.itsNcls, sd2.itsChi2, sd2.itsMeanClsSizeL,
          sd2.rigidity, sd2.tpcSignal, sd2.tpcNsigma, sd2.tpcNsigmaNhp, sd2.tpcNsigmaNlp, sd2.tofMass, sd2.dcaXY, sd2.dcaZ, sd2.isPvContributor);
      if (isMC && cfgMCCombined)
        outputTableMcThreeTwo(
          cand.pdgCode, cand.isMatterMC, cand.isReconstructed, cand.isPhysicalPrimary, cand.passedEvSelMC, cand.yGen, cand.ptGen, cand.ctGen,
          cand.cpaPvGen, cand.pxGen, cand.pyGen, cand.pzGen, cand.pvxGen, cand.pvyGen, cand.pvzGen, cand.svxGen, cand.svyGen, cand.svzGen,
          cand.species, cand.isMatter, cand.cent, cand.occu, cand.passedEvSel, cand.mass, cand.y, cand.pt, cand.ct, cand.cpaPv, cand.maxDcaTracks, cand.maxDcaTracksSv, cand.dcaToPvXY,
          cand.dcaToPvZ, cand.devToPvXY, cand.chi2, cand.pvx, cand.pvy, cand.pvz, cand.svx, cand.svy, cand.svz, cand.px, cand.py, cand.pz, cand.collisionMcTrue,
          hypDaughter.svx, hypDaughter.svy, hypDaughter.svz, d0.x, d0.y, d0.z, d0.px, d0.py, d0.pz, hypDaughter.mass, hypDaughter.ct, hypDaughter.cpaPv,
          hypDaughter.maxDcaTracks, hypDaughter.dcaToPvXY, hypDaughter.dcaToPvZ, hypDaughter.dcaToVtxXY, hypDaughter.dcaToVtxZ, hypDaughter.chi2,
          d1.x, d1.y, d1.z, d1.px, d1.py, d1.pz, d1.tpcNcls, d1.tpcChi2, d1.itsNcls, d1.itsChi2, d1.itsMeanClsSizeL,
          d1.rigidity, d1.tpcSignal, d1.tpcNsigma, d1.tpcNsigmaNhp, d1.tpcNsigmaNlp, d1.tofMass, d1.dcaXY, d1.dcaZ, d1.isPvContributor,
          d2.x, d2.y, d2.z, d2.px, d2.py, d2.pz, d2.tpcNcls, d2.tpcChi2, d2.itsNcls, d2.itsChi2, d2.itsMeanClsSizeL,
          d2.rigidity, d2.tpcSignal, d2.tpcNsigma, d2.tpcNsigmaNhp, d2.tpcNsigmaNlp, d2.tofMass, d2.dcaXY, d2.dcaZ, d2.isPvContributor,
          d0.subMass, d1.subMass, d2.subMass,
          sd1.x, sd1.y, sd1.z, sd1.px, sd1.py, sd1.pz, sd1.tpcNcls, sd1.tpcChi2, sd1.itsNcls, sd1.itsChi2, sd1.itsMeanClsSizeL,
          sd1.rigidity, sd1.tpcSignal, sd1.tpcNsigma, sd1.tpcNsigmaNhp, sd1.tpcNsigmaNlp, sd1.tofMass, sd1.dcaXY, sd1.dcaZ, sd1.isPvContributor,
          sd2.x, sd2.y, sd2.z, sd2.px, sd2.py, sd2.pz, sd2.tpcNcls, sd2.tpcChi2, sd2.itsNcls, sd2.itsChi2, sd2.itsMeanClsSizeL,
          sd2.rigidity, sd2.tpcSignal, sd2.tpcNsigma, sd2.tpcNsigmaNhp, sd2.tpcNsigmaNlp, sd2.tofMass, sd2.dcaXY, sd2.dcaZ, sd2.isPvContributor);
    }
  }
  //___________________________________________________________________________________________________________________________________________________________
  void fillCandidatePrim(HyperNucleus& cand, aod::HypKfHypNuc const& hypNuc, aod::HypKfHypNucs const& hypNucs, aod::HypKfColls const& colls, aod::HypKfTracks const& tracks, aod::HypKfDaughtAdds const& daughterAdds, aod::HypKfSubDs const& subDs)
  {
    auto coll = hypNuc.hypKfColl();
    cand.ct = ct(coll, hypNuc);
    cand.cpaPv = cpa(coll, hypNuc);
    fillCandidate(cand, hypNuc, hypNucs, colls, tracks, daughterAdds, subDs);
  }
  //___________________________________________________________________________________________________________________________________________________________
  void fillCandidateSec(HyperNucleus& cand, aod::HypKfHypNuc const& hypNuc, aod::HypKfHypNuc const& mother, aod::HypKfHypNucs const& hypNucs, aod::HypKfColls const& colls, aod::HypKfTracks const& tracks, aod::HypKfDaughtAdds const& daughterAdds, aod::HypKfSubDs const& subDs)
  {
    cand.ct = ct(mother, hypNuc);
    cand.cpaPv = cpa(mother, hypNuc);
    fillCandidate(cand, hypNuc, hypNucs, colls, tracks, daughterAdds, subDs);
  }
  //___________________________________________________________________________________________________________________________________________________________
  void fillCandidate(HyperNucleus& cand, aod::HypKfHypNuc const& hypNuc, aod::HypKfHypNucs const&, aod::HypKfColls const&, aod::HypKfTracks const&, aod::HypKfDaughtAdds const&, aod::HypKfSubDs const&)
  {
    cand.daughterTracks.clear();
    cand.subDaughterMassVec.clear();
    auto coll = hypNuc.hypKfColl();
    auto addOns = hypNuc.hypKfDaughtAdd_as<aod::HypKfDaughtAdds>();
    auto posVec = posVector(addOns);
    cand.species = std::abs(hypNuc.species());
    cand.isPrimaryCandidate = hypNuc.primary();
    cand.isMatter = hypNuc.isMatter();
    cand.cent = coll.centFT0C();
    cand.occu = coll.occupancy();
    cand.passedEvSel = coll.passedEvSel();
    cand.mass = hypNuc.mass();
    cand.y = hypNuc.y();
    cand.pt = hypNuc.pt();
    cand.maxDcaTracks = maxValue(dcaTracksAll(posVec, "XY"));
    cand.maxDcaTracksSv = maxValue(dcaTrackSvAll(posVec, hypNuc, "XY"));
    cand.dcaToPvXY = hypNuc.dcaToPvXY();
    cand.dcaToPvZ = hypNuc.dcaToPvZ();
    cand.dcaToVtxXY = hypNuc.dcaToVtxXY();
    cand.dcaToVtxZ = hypNuc.dcaToVtxZ();
    cand.devToPvXY = hypNuc.devToPvXY();
    cand.chi2 = hypNuc.chi2();
    cand.pvx = coll.posX();
    cand.pvy = coll.posY();
    cand.pvz = coll.posZ();
    cand.svx = hypNuc.svx();
    cand.svy = hypNuc.svy();
    cand.svz = hypNuc.svz();
    cand.px = hypNuc.px();
    cand.py = hypNuc.py();
    cand.pz = hypNuc.pz();
    if (hypNuc.hypDaughterId() >= 0) {
      TrackProperties hypDaughter;
      cand.daughterTracks.push_back(hypDaughter);
    }
    auto daughterTracks = hypNuc.hypKfTrack_as<aod::HypKfTracks>();
    for (const auto& track : daughterTracks) {
      TrackProperties daughter;
      daughter.tpcNcls = track.tpcNcluster();
      daughter.itsNcls = track.itsNcluster();
      daughter.tpcChi2 = track.tpcChi2NCl();
      daughter.itsChi2 = track.itsChi2NCl();
      daughter.itsMeanClsSize = track.itsMeanClsSize();
      daughter.itsMeanClsSizeL = track.itsMeanClsSize() * track.lambda();
      daughter.rigidity = track.rigidity();
      daughter.tpcSignal = track.tpcSignal();
      daughter.tpcNsigma = track.tpcNsigma();
      daughter.tpcNsigmaNhp = track.tpcNsigmaNhp();
      daughter.tpcNsigmaNlp = track.tpcNsigmaNlp();
      daughter.tofMass = track.tofMass();
      daughter.dcaXY = track.dcaXY();
      daughter.dcaZ = track.dcaZ();
      daughter.isPvContributor = track.isPVContributor();
      cand.daughterTracks.push_back(daughter);
    }
    int trackCount = 0;
    for (const auto& addOn : addOns) {
      cand.daughterTracks.at(trackCount).x = addOn.x();
      cand.daughterTracks.at(trackCount).y = addOn.y();
      cand.daughterTracks.at(trackCount).z = addOn.z();
      cand.daughterTracks.at(trackCount).px = addOn.px();
      cand.daughterTracks.at(trackCount).py = addOn.py();
      cand.daughterTracks.at(trackCount).pz = addOn.py();
      trackCount++;
    }
    cand.nSingleDaughters = trackCount;
    if (cand.nSingleDaughters < 3)
      return;

    trackCount = 0;
    auto subDaughters = hypNuc.hypKfSubD_as<aod::HypKfSubDs>();
    for (const auto& subDaughter : subDaughters) {
      cand.daughterTracks.at(trackCount++).subMass = subDaughter.subMass();
    }
  }
  //___________________________________________________________________________________________________________________________________________________________

  void processMC(aod::HypKfMcParts const& mcHypNucs, aod::HypKfHypNucs const& hypNucs, aod::HypKfMcColls const&, aod::HypKfColls const& hypKfColls, aod::HypKfTracks const& hypKfTrks, aod::HypKfDaughtAdds const& hypKfDAdd, aod::HypKfSubDs const& hypKfDSub)
  {
    isMC = true;
    for (const auto& mcHypNuc : mcHypNucs) {
      if (cfgSpecies && std::abs(mcHypNuc.species()) != cfgSpecies)
        continue;
      auto mcColl = mcHypNuc.hypKfMcColl();
      const auto mcParticleIdx = mcHypNuc.globalIndex();
      auto hypNucsByMc = hypNucs.sliceBy(perMcParticle, mcParticleIdx);
      HyperNucleus candidate, hypNucDaughter;
      candidate.pdgCode = mcHypNuc.pdgCode();
      candidate.isMatterMC = mcHypNuc.isMatter();
      candidate.isPhysicalPrimary = mcHypNuc.isPhysicalPrimary();
      candidate.passedEvSelMC = mcColl.passedEvSel();
      candidate.yGen = mcHypNuc.y();
      candidate.ptGen = mcHypNuc.pt();
      candidate.ctGen = ct(mcColl, mcHypNuc);
      candidate.isReconstructed = 0;
      candidate.cpaPvGen = cpa(mcColl, mcHypNuc);
      candidate.pxGen = mcHypNuc.px();
      candidate.pyGen = mcHypNuc.py();
      candidate.pzGen = mcHypNuc.pz();
      candidate.pvxGen = mcColl.posX();
      candidate.pvyGen = mcColl.posY();
      candidate.pvzGen = mcColl.posZ();
      candidate.svxGen = mcHypNuc.svx();
      candidate.svyGen = mcHypNuc.svy();
      candidate.svzGen = mcHypNuc.svz();
      for (const auto& hypNuc : hypNucsByMc) {
        auto coll = hypNuc.hypKfColl();
        if (coll.hypKfMcCollId() == mcHypNuc.hypKfMcCollId()) {
          candidate.collisionMcTrue = true;
        }
        candidate.isReconstructed++;
        fillCandidatePrim(candidate, hypNucs.rawIteratorAt(hypNuc.globalIndex()), hypNucs, hypKfColls, hypKfTrks, hypKfDAdd, hypKfDSub);
        if (hypNuc.hypDaughterId() >= 0) {
          fillCandidateSec(hypNucDaughter, hypNucs.rawIteratorAt(hypNuc.hypDaughterId()), hypNucs.rawIteratorAt(hypNuc.globalIndex()), hypNucs, hypKfColls, hypKfTrks, hypKfDAdd, hypKfDSub);
        }
      }
      fillTable(candidate, hypNucDaughter);
      hPt[0]->Fill(mcHypNuc.pt());
      if (candidate.isReconstructed)
        hPt[1]->Fill(candidate.pt);
    }
    hPt[2]->Divide(hPt[1].get(), hPt[0].get());
  }
  PROCESS_SWITCH(HypKfTreeCreator, processMC, "MC Gen tree", false);

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
    for (const auto& value : vec)
      sum += value;
    return sum / vec.size();
  }
  float mean2Value(std::vector<float> vec)
  {
    float sum = 0;
    for (const auto& value : vec)
      sum += (value * value);
    return std::sqrt(sum / vec.size());
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
    for (const auto& pos : addons) {
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
    adaptAnalysisTask<HypKfTreeCreator>(cfgc)};
}
