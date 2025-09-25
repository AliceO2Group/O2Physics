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

/// \file D0CalibTables.h
/// \brief Definitions of derived tables produced by data creator for D0 calibration studies
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#ifndef DPG_TASKS_AOTTRACK_D0CALIBTABLES_H_
#define DPG_TASKS_AOTTRACK_D0CALIBTABLES_H_

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>

#include <sys/types.h>

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace o2
{
namespace hf_calib
{
enum D0MassHypo : uint8_t {
  D0 = 1,
  D0Bar,
  D0AndD0Bar,
  ND0MassHypos
};

const float toMicrometers = 10000.; // from cm to µm

/// It compresses a value to a int8_t with a given precision
///\param origValue is the original values
///\param precision is the desired precision
///\return The value compressed to a int8_t
template <typename T>
int8_t getCompressedInt8(T origValue, double precision)
{
  int roundValue = static_cast<int>(std::round(origValue / precision));
  return static_cast<int8_t>(std::clamp(roundValue, -128, 128));
}

/// It compresses a value to a uint8_t with a given precision
///\param origValue is the original values
///\param precision is the desired precision
///\return The value compressed to a uint8_t
template <typename T>
uint8_t getCompressedUint8(T origValue, double precision)
{
  int roundValue = static_cast<int>(std::round(origValue / precision));
  return static_cast<uint8_t>(std::clamp(roundValue, 0, 255));
}

/// It compresses a value to a uint16_t with a given precision
///\param origValue is the original values
///\param precision is the desired precision
///\return The value compressed to a uint16_t
template <typename T>
uint16_t getCompressedUint16(T origValue, double precision)
{
  int roundValue = static_cast<int>(std::round(origValue / precision));
  return static_cast<uint16_t>(std::clamp(roundValue, 0, 65535));
}

/// It uses a sinh-based scaling function, which provides a compromise between fixed-step and relative quantization.
// This approach reflects typical resolution formulas and is well-suited for detector calibration data.
///\param origValue is the original value
///\param sigma0 is a asinh parameter
///\param sigma1 is a asinh parameter
///\param clampMin is the maximum value
///\param clampMax is the minimum value
///\return The value compressed
int codeSqrtScaling(float origValue, float sigma0, float sigma1, int clampMin, int clampMax)
{
  float codeF = std::asinh((sigma1 * origValue) / sigma0) / sigma0;
  return std::clamp(static_cast<int>(std::round(codeF)), clampMin, clampMax);
}

/// It compresses the decay length (10 micron precision)
///\param decLen is the decay length in cm
///\return The decay length compressed to a uint8_t with 10 micron precision
template <typename T>
uint8_t getCompressedDecayLength(T decLen)
{
  return getCompressedUint8(decLen * hf_calib::toMicrometers, 10);
}

/// It compresses the normalised decay length (0.5 precision)
///\param normDecLen is the normalised decay length
///\return The normalised decay length compressed to a uint8_t with 0.5 precision
template <typename T>
uint8_t getCompressedNormDecayLength(T normDecLen)
{
  return getCompressedUint8(normDecLen, 0.5);
}

/// It compresses the pointing angle (0.005 precision)
///\param pointAngle is the pointing angle
///\return The pointing angle compressed to a uint8_t with 0.005 precision
template <typename T>
uint8_t getCompressedPointingAngle(T pointAngle)
{
  return getCompressedUint8(pointAngle, 0.005);
}

/// It compresses the cosine of pointing angle (0.001 precision)
///\param cosPa is the cosine of pointing angle
///\return The cosine of pointing angle compressed to a uint8_t with 0.001 precision
template <typename T>
int8_t getCompressedCosPa(T cosPa)
{
  return getCompressedUint8(cosPa - 0.75, 0.001); // in the range from 0.75 to 1
}

/// It compresses the chi2
///\param chi2 is the chi2
///\return The chi2 compressed to a uint8_t
template <typename T>
int8_t getCompressedChi2(T chi2)
{
  uint8_t compressedChi2 = static_cast<uint8_t>(codeSqrtScaling(chi2, 0.015, 0.015, 0, 255));
  return compressedChi2;
}

/// It compresses the number of sigma
///\param numSigma is the number of sigma
///\return The number of sigma compressed to a int8_t
template <typename T>
int8_t getCompressedNumSigmaPid(T numSigma)
{
  int8_t compressedNumSigma = static_cast<int8_t>(codeSqrtScaling(numSigma, 0.05, 0.05, -128, 128));
  return compressedNumSigma;
}

/// It compresses the bdt score (1./65535 precision)
///\param bdtScore is the bdt score
///\return The bdt score compressed to a uint16_t with 1./65535 precision
template <typename T>
uint16_t getCompressedBdtScoreBkg(T bdtScore)
{
  return getCompressedUint16(bdtScore, 1. / 65535);
}

/// It compresses the bdt score (1./255 precision)
///\param bdtScore is the bdt score
///\return The bdt score compressed to a uint8_t with 1./255 precision
template <typename T>
uint8_t getCompressedBdtScoreSgn(T bdtScore)
{
  return getCompressedUint8(bdtScore, 1. / 255);
}

/// It compresses the occupancy value
///\param occupancy is the occupancy value
///\return The number of occupancy compressed to a uint8_t
template <typename T>
uint8_t getCompressedOccupancy(T occupancy)
{
  uint8_t compressedOcc = static_cast<uint8_t>(codeSqrtScaling(occupancy, 0.04, 0.04, 0, 255));
  return compressedOcc;
}

static constexpr int NBinsPtTrack = 6;
static constexpr int NCutVarsTrack = 4;
constexpr float BinsPtTrack[NBinsPtTrack + 1] = {
  0,
  0.5,
  1.0,
  1.5,
  2.0,
  3.0,
  1000.0};
auto vecBinsPtTrack = std::vector<float>{BinsPtTrack, BinsPtTrack + NBinsPtTrack + 1};

// default values for the dca_xy and dca_z cuts of displaced tracks
constexpr float CutsTrack[NBinsPtTrack][NCutVarsTrack] = {{0.0015, 2., 0.0000, 2.},  /* 0   < pt < 0.5 */
                                                          {0.0015, 2., 0.0000, 2.},  /* 0.5 < pt < 1 */
                                                          {0.0015, 2., 0.0000, 2.},  /* 1   < pt < 1.5 */
                                                          {0.0015, 2., 0.0000, 2.},  /* 1.5 < pt < 2 */
                                                          {0.0000, 2., 0.0000, 2.},  /* 2   < pt < 3 */
                                                          {0.0000, 2., 0.0000, 2.}}; /* 3   < pt < 1000 */
// row labels
static const std::vector<std::string> labelsPtTrack{};

// column labels
static const std::vector<std::string> labelsCutVarTrack = {"min_dcaxytoprimary", "max_dcaxytoprimary", "min_dcaztoprimary", "max_dcaztoprimary"};

static constexpr int NBinsPtCand = 10;
static constexpr int NCutVarsCand = 10;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr float BinsPtCand[NBinsPtCand + 1] = {
  0,
  1.0,
  2.0,
  3.0,
  4.0,
  6.0,
  8.0,
  12.0,
  24.0,
  50.0,
  1000.0};
auto vecBinsPtCand = std::vector<float>{BinsPtCand, BinsPtCand + NBinsPtCand + 1};

// default values for the cuts
constexpr float CutsCand[NBinsPtCand][NCutVarsCand] = {{0.400, 0., 10., 10., 0.97, 0.97, 0, 2, 0.01, 0.01},  /* 0  < pT < 1    */
                                                       {0.400, 0., 10., 10., 0.97, 0.97, 0, 2, 0.01, 0.01},  /* 1  < pT < 2    */
                                                       {0.400, 0., 10., 10., 0.95, 0.95, 0, 2, 0.01, 0.01},  /* 2  < pT < 3    */
                                                       {0.400, 0., 10., 10., 0.95, 0.95, 0, 2, 0.01, 0.01},  /* 3  < pT < 4    */
                                                       {0.400, 0., 10., 10., 0.95, 0.95, 0, 2, 0.01, 0.01},  /* 4  < pT < 6    */
                                                       {0.400, 0., 10., 10., 0.95, 0.95, 0, 2, 0.01, 0.01},  /* 6  < pT < 8    */
                                                       {0.400, 0., 10., 10., 0.95, 0.95, 0, 2, 0.01, 0.01},  /* 8  < pT < 12   */
                                                       {0.400, 0., 10., 10., 0.95, 0.95, 0, 2, 0.01, 0.01},  /* 12 < pT < 24   */
                                                       {0.400, 0., 10., 10., 0.95, 0.95, 0, 2, 0.01, 0.01},  /* 24 < pT < 50   */
                                                       {0.400, 0., 10., 10., 0.95, 0.95, 0, 2, 0.01, 0.01}}; /* 50 < pT < 1000 */

// row labels
static const std::vector<std::string> labelsPtCand = {
  "pT bin 0",
  "pT bin 1",
  "pT bin 2",
  "pT bin 3",
  "pT bin 4",
  "pT bin 5",
  "pT bin 6",
  "pT bin 7",
  "pT bin 8",
  "pT bin 9"};

// column labels
static const std::vector<std::string> labelsCutVarCand = {"delta inv. mass", "max d0d0", "max pointing angle", "max pointing angle XY", "min cos pointing angle", "min cos pointing angle XY", "min norm decay length", "min norm decay length XY", "min decay length", "min decay length XY"};

static constexpr int NBinsPtMl = 10;
// default values for the pT bin edges (can be used to configure histogram axis)
// offset by 1 from the bin numbers in cuts array
constexpr double BinsPtMl[NBinsPtMl + 1] = {
  0,
  1.0,
  2.0,
  3.0,
  4.0,
  6.0,
  8.0,
  12.0,
  24.0,
  50.0,
  1000.0};
auto vecBinsPtMl = std::vector<double>{BinsPtMl, BinsPtMl + NBinsPtMl + 1};

// default values for the cuts
constexpr double CutsMl[NBinsPtMl][3] = {{1., 0., 0.},  /* 0  < pT < 1    */
                                         {1., 0., 0.},  /* 1  < pT < 2    */
                                         {1., 0., 0.},  /* 2  < pT < 3    */
                                         {1., 0., 0.},  /* 3  < pT < 4    */
                                         {1., 0., 0.},  /* 4  < pT < 6    */
                                         {1., 0., 0.},  /* 6  < pT < 8    */
                                         {1., 0., 0.},  /* 8  < pT < 12   */
                                         {1., 0., 0.},  /* 12 < pT < 24   */
                                         {1., 0., 0.},  /* 24 < pT < 50   */
                                         {1., 0., 0.}}; /* 50 < pT < 1000 */

// row labels
static const std::vector<std::string> labelsPtMl = {
  "pT bin 0",
  "pT bin 1",
  "pT bin 2",
  "pT bin 3",
  "pT bin 4",
  "pT bin 5",
  "pT bin 6",
  "pT bin 7",
  "pT bin 8",
  "pT bin 9"};

// column labels
static const std::vector<std::string> labelsCutMl = {"max BDT score bkg", "min BDT score prompt", "min BDT score nonprompt"};
} // namespace hf_calib

namespace aod
{
namespace hf_calib
{
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);                 //! Run number
DECLARE_SOA_COLUMN(Orbit, orbit, uint32_t);                    //! orbit ID
DECLARE_SOA_COLUMN(CentFT0C, centFT0C, uint8_t);               //! FTOC centrality
DECLARE_SOA_COLUMN(OccupancyTracks, occupancyTracks, uint8_t); //! FT0 occupancy
DECLARE_SOA_COLUMN(OccupancyFT0C, occupancyFT0C, uint8_t);     //! FT0 occupancy
} // namespace hf_calib

DECLARE_SOA_TABLE(D0CalibColls, "AOD", "D0CALIBCOLL",
                  o2::soa::Index<>,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  collision::CovXX,
                  collision::CovXY,
                  collision::CovXZ,
                  collision::CovYY,
                  collision::CovYZ,
                  collision::CovZZ,
                  collision::NumContrib,
                  hf_calib::CentFT0C,
                  hf_calib::OccupancyTracks,
                  hf_calib::OccupancyFT0C,
                  hf_calib::Orbit,
                  hf_calib::RunNumber);

namespace hf_calib
{
DECLARE_SOA_INDEX_COLUMN_FULL(Collision, collision, int, D0CalibColls, ""); //! Index of collision
DECLARE_SOA_COLUMN(TpcNumSigmaPi, tpcNumSigmaPi, int8_t);                   //! compressed NsigmaTPC for pions
DECLARE_SOA_COLUMN(TpcNumSigmaKa, tpcNumSigmaKa, int8_t);                   //! compressed NsigmaTPC for kaons
DECLARE_SOA_COLUMN(TofNumSigmaPi, tofNumSigmaPi, int8_t);                   //! compressed NsigmaTOF for pions
DECLARE_SOA_COLUMN(TofNumSigmaKa, tofNumSigmaKa, int8_t);                   //! compressed NsigmaTOF for kaons
DECLARE_SOA_COLUMN(ITSChi2NCl, itsChi2NCl, uint8_t);                        //! compressed NsigmaTOF for kaons // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(TPCChi2NCl, tpcChi2NCl, uint8_t);                        //! compressed NsigmaTOF for kaons // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(TRDChi2, trdChi2, uint8_t);                              //! compressed NsigmaTOF for kaons // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(TOFChi2, tofChi2, uint8_t);                              //! compressed NsigmaTOF for kaons // o2-linter: disable=name/o2-column
DECLARE_SOA_COLUMN(CmoPrimUnfm80, cmoPrimUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoFV0AUnfm80, cmoFV0AUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoFT0AUnfm80, cmoFT0AUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoFT0CUnfm80, cmoFT0CUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoPrimUnfm80, cwmoPrimUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoFV0AUnfm80, cwmoFV0AUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoFT0AUnfm80, cwmoFT0AUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoFT0CUnfm80, cwmoFT0CUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CmoRobustT0V0PrimUnfm80, cmoRobustT0V0PrimUnfm80, uint8_t);
DECLARE_SOA_COLUMN(CwmoRobustT0V0PrimUnfm80, cwmoRobustT0V0PrimUnfm80, uint8_t);
} // namespace hf_calib

DECLARE_SOA_TABLE(D0CalibTracks, "AOD", "D0CALIBTRACK",
                  o2::soa::Index<>,
                  /// *** collision index
                  hf_calib::CollisionId,
                  /// *** track pars
                  track::X,
                  track::Alpha,
                  track::Y,
                  track::Z,
                  track::Snp,
                  track::Tgl,
                  track::Signed1Pt,
                  /// *** track covs
                  track::CYY,
                  track::CZY,
                  track::CZZ,
                  track::CSnpY,
                  track::CSnpZ,
                  track::CSnpSnp,
                  track::CTglY,
                  track::CTglZ,
                  track::CTglSnp,
                  track::CTglTgl,
                  track::C1PtY,
                  track::C1PtZ,
                  track::C1PtSnp,
                  track::C1PtTgl,
                  track::C1Pt21Pt2,
                  /// *** track extra (static)
                  track::TPCInnerParam,
                  track::Flags,
                  track::ITSClusterSizes,
                  track::TPCNClsFindable,
                  track::TPCNClsFindableMinusFound,
                  track::TPCNClsFindableMinusCrossedRows,
                  track::TPCNClsShared,
                  track::TRDPattern,
                  hf_calib::ITSChi2NCl,
                  hf_calib::TPCChi2NCl,
                  hf_calib::TRDChi2,
                  hf_calib::TOFChi2,
                  track::TPCSignal,
                  track::TRDSignal,
                  track::Length,
                  track::TOFExpMom,
                  track::TrackTime,
                  track::TrackTimeRes,
                  /// *** track QA
                  trackqa::TPCTime0,
                  trackqa::TPCdEdxNorm,
                  trackqa::TPCDCAR,
                  trackqa::TPCDCAZ,
                  trackqa::TPCClusterByteMask,
                  trackqa::TPCdEdxMax0R,
                  trackqa::TPCdEdxMax1R,
                  trackqa::TPCdEdxMax2R,
                  trackqa::TPCdEdxMax3R,
                  trackqa::TPCdEdxTot0R,
                  trackqa::TPCdEdxTot1R,
                  trackqa::TPCdEdxTot2R,
                  trackqa::TPCdEdxTot3R,
                  trackqa::DeltaRefContParamY,
                  trackqa::DeltaRefContParamZ,
                  trackqa::DeltaRefContParamSnp,
                  trackqa::DeltaRefContParamTgl,
                  trackqa::DeltaRefContParamQ2Pt,
                  trackqa::DeltaRefGloParamY,
                  trackqa::DeltaRefGloParamZ,
                  trackqa::DeltaRefGloParamSnp,
                  trackqa::DeltaRefGloParamTgl,
                  trackqa::DeltaRefGloParamQ2Pt,
                  trackqa::DeltaTOFdX,
                  trackqa::DeltaTOFdZ,
                  /// *** DCA, Nsigma
                  track::DcaXY,
                  track::DcaZ,
                  hf_calib::TpcNumSigmaPi,
                  hf_calib::TpcNumSigmaKa,
                  hf_calib::TofNumSigmaPi,
                  hf_calib::TofNumSigmaKa,
                  /// *** Occupancy variables
                  hf_calib::CmoPrimUnfm80,
                  hf_calib::CmoFV0AUnfm80,
                  hf_calib::CmoFT0AUnfm80,
                  hf_calib::CmoFT0CUnfm80,
                  hf_calib::CwmoPrimUnfm80,
                  hf_calib::CwmoFV0AUnfm80,
                  hf_calib::CwmoFT0AUnfm80,
                  hf_calib::CwmoFT0CUnfm80,
                  hf_calib::CmoRobustT0V0PrimUnfm80,
                  hf_calib::CwmoRobustT0V0PrimUnfm80);

namespace hf_calib
{
DECLARE_SOA_INDEX_COLUMN_FULL(TrackPos, trackPos, int, D0CalibTracks, "_Pos"); //! Index of positive track
DECLARE_SOA_INDEX_COLUMN_FULL(TrackNeg, trackNeg, int, D0CalibTracks, "_Neg"); //! Index of negative track
DECLARE_SOA_COLUMN(MassHypo, massHypo, uint8_t);                               //! mass hypothesis for D0 (D0, D0bar, or both)
DECLARE_SOA_COLUMN(Pt, pt, float);                                             //! D0-candidate pT
DECLARE_SOA_COLUMN(Eta, eta, float);                                           //! D0-candidate eta
DECLARE_SOA_COLUMN(Phi, phi, float);                                           //! D0-candidate phi
DECLARE_SOA_COLUMN(InvMassD0, invMassD0, float);                               //! invariant mass (D0 hypothesis)
DECLARE_SOA_COLUMN(InvMassD0bar, invMassD0bar, float);                         //! invariant mass (D0bar hypothesis)
DECLARE_SOA_COLUMN(DecLength, decLength, uint8_t);                             //! compressed decay length
DECLARE_SOA_COLUMN(DecLengthXY, decLengthXY, uint8_t);                         //! compressed decay length XY
DECLARE_SOA_COLUMN(NormDecLength, normDecLength, uint8_t);                     //! compressed normalised decay length
DECLARE_SOA_COLUMN(NormDecLengthXY, normDecLengthXY, uint8_t);                 //! compressed normalised decay length XY
DECLARE_SOA_COLUMN(CosPa, cosPa, uint8_t);                                     //! compressed cosine of pointing angle
DECLARE_SOA_COLUMN(CosPaXY, cosPaXY, uint8_t);                                 //! compressed cosine of pointing angle XY
DECLARE_SOA_COLUMN(PointingAngle, pointingAngle, uint8_t);                     //! compressed pointing angle
DECLARE_SOA_COLUMN(PointingAngleXY, pointingAngleXY, uint8_t);                 //! compressed pointing angle XY
DECLARE_SOA_COLUMN(DecVtxChi2, decVtxChi2, uint8_t);                           //! compressed decay vertex chi2
DECLARE_SOA_COLUMN(BdtScoreBkgD0, bdtScoreBkgD0, uint16_t);                    //! compressed BDT score (bkg, D0 mass hypo)
DECLARE_SOA_COLUMN(BdtScorePromptD0, bdtScorePromptD0, uint8_t);               //! compressed BDT score (prompt, D0 mass hypo)
DECLARE_SOA_COLUMN(BdtScoreNonpromptD0, bdtScoreNonpromptD0, uint8_t);         //! compressed BDT score (non-prompt, D0 mass hypo)
DECLARE_SOA_COLUMN(BdtScoreBkgD0bar, bdtScoreBkgD0bar, uint16_t);              //! compressed BDT score (bkg, D0bar mass hypo)
DECLARE_SOA_COLUMN(BdtScorePromptD0bar, bdtScorePromptD0bar, uint8_t);         //! compressed BDT score (prompt, D0bar mass hypo)
DECLARE_SOA_COLUMN(BdtScoreNonpromptD0bar, bdtScoreNonpromptD0bar, uint8_t);   //! compressed BDT score (non-prompt, D0bar mass hypo)
} // namespace hf_calib

DECLARE_SOA_TABLE(D0CalibCands, "AOD", "D0CALIBCAND",
                  o2::soa::Index<>,
                  hf_calib::CollisionId,
                  hf_calib::TrackPosId,
                  hf_calib::TrackNegId,
                  hf_calib::MassHypo,
                  hf_calib::Pt,
                  hf_calib::Eta,
                  hf_calib::Phi,
                  hf_calib::InvMassD0,
                  hf_calib::InvMassD0bar,
                  hf_calib::DecLength,
                  hf_calib::DecLengthXY,
                  hf_calib::NormDecLength,
                  hf_calib::NormDecLengthXY,
                  hf_calib::CosPa,
                  hf_calib::CosPaXY,
                  hf_calib::PointingAngle,
                  hf_calib::PointingAngleXY,
                  hf_calib::DecVtxChi2,
                  hf_calib::BdtScoreBkgD0,
                  hf_calib::BdtScorePromptD0,
                  hf_calib::BdtScoreNonpromptD0,
                  hf_calib::BdtScoreBkgD0bar,
                  hf_calib::BdtScorePromptD0bar,
                  hf_calib::BdtScoreNonpromptD0bar);
} // namespace aod
} // namespace o2
#endif // DPG_TASKS_AOTTRACK_D0CALIBTABLES_H_
