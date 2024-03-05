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
/// \file nuclei_flow.cxx
/// \brief Measure flow of nuclei using LFSlimNucleiTable data format
///

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/Configurable.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "PWGLF/DataModel/LFSlimNucleiTables.h"
#include "Common/Core/EventPlaneHelper.h"
#include "TMath.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

namespace nuclei_spectra
{

constexpr float charges[5]{1.f, 1.f, 1.f, 2.f, 2.f};

constexpr float masses[5]{MassProton, MassDeuteron, MassTriton, MassHelium3, MassAlpha};

constexpr float bbMomScalingDefault[5][2]{
  {1., 1.},
  {1., 1.},
  {1., 1.},
  {1., 1.},
  {1., 1.}};

constexpr float betheBlochDefault[5][6]{
  {-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32},
  {-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32},
  {-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32},
  {-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32},
  {-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};

static const std::vector<std::string> names{"proton", "deuteron", "triton", "He3", "alpha"};
static const std::vector<std::string> chargeLabelNames{"Positive", "Negative"};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};

enum CandBits {
  kProton = BIT(0),
  kDeuteron = BIT(1),
  kTriton = BIT(2),
  kHe3 = BIT(3),
  kHe4 = BIT(4),
  kHasTOF = BIT(5),
  kIsReconstructed = BIT(6),
  kIsAmbiguous = BIT(7), /// just a placeholder now
  kPositive = BIT(8),
  kIsPhysicalPrimary = BIT(9), /// MC flags starting from the second half of the short
  kIsSecondaryFromMaterial = BIT(10),
  kIsSecondaryFromWeakDecay = BIT(11) /// the last 4 bits are reserved for the PID in tracking
};

} // namespace nuclei_spectra

namespace flow
{
enum kDetector {
  kFV0A = 0,
  kFT0M = 1,
  kFT0A = 2,
  kFT0C = 3,
  kTPCpos = 4,
  kTPCneg = 5
};
} // namespace flow

struct nucleiFlow {

  Configurable<int> cfgCentDetector{"cfgCentDetector", 0, "Detector for centrality estimation (FV0A: 0, FT0M: 1, FT0A: 2, FT0C: 3)"};
  Configurable<int> cfgHarmonic{"cfgHarmonic", 2, "cfgHarmonic number"};
  Configurable<int> cfgQvecDetector{"cfgQvecDetector", 0, "Detector for Q vector estimation (FV0A: 0, FT0M: 1, FT0A: 2, FT0C: 3, TPC Pos: 4, TPC Neg: 5)"};
  Configurable<int> cfgSpecies{"cfgSpecies", 3, "Species under study (proton: 0, deuteron: 1, triton: 2, helion: 3, alpha: 4)"};

  Configurable<float> cfgNclusTPCcut{"cfgNclusTPCcut", 70, "Minimum number of TPC clusters"};
  Configurable<float> cfgDCAxyCut{"cfgDCAxyCut", 0.1, "Cut on DCAxy (cm)"};
  Configurable<float> cfgDCAzCut{"cfgDCAzCut", 1, "Cut on DCAz (cm)"};
  Configurable<int> cfgItsClusSizeCut{"cfgItsClusSizeCut", 4, "Cut on the average ITS cluster size."};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  ConfigurableAxis nSigmaBins{"nSigmaBins", {200, -5.f, 5.f}, "Binning for n sigma"};
  ConfigurableAxis centBins{"centBins", {111, -0.5f, 110.5f}, "Binning for centrality"};
  ConfigurableAxis ptBins{"ptBins", {50, 0.f, 5.f}, "Binning for pt"};
  ConfigurableAxis spBins{"spBins", {100, -1.f, 1.f}, "Binning for scalar product"};

  /// \brief momentum scaling-factor for TPC Bethe-Bloch
  Configurable<LabeledArray<float>> cfgMomentumScalingBetheBloch{"cfgMomentumScalingBetheBloch", {nuclei_spectra::bbMomScalingDefault[0], 5, 2, nuclei_spectra::names, nuclei_spectra::chargeLabelNames}, "TPC Bethe-Bloch momentum scaling for light nuclei"};

  /// \brief bethe-bloch parameters
  Configurable<LabeledArray<float>> cfgBetheBlochParams{"cfgBetheBlochParams", {nuclei_spectra::betheBlochDefault[0], 5, 6, nuclei_spectra::names, nuclei_spectra::betheBlochParNames}, "TPC Bethe-Bloch parameterisation for light nuclei"};

  // Selected nuclei tracks with flow info
  using TracksWithFlow = soa::Join<aod::NucleiTable, aod::NucleiTableFlow>;
  using TrackWithFlow = TracksWithFlow::iterator;

  EventPlaneHelper epHelper;

  /// \brief Get n-sigma TPC
  /// \tparam T type for the track
  /// \param candidate track candidate
  /// \param iSpecies 0: proton, 1: deuteron, 2: triton, 3: He3, 4: Alpha
  /// \param iCharge 0: positive, 1: negative
  /// \return n-sigma TPC for candidates with iSpecies hypothesis
  template <typename T>
  float getNSigmaTPC(const T& candidate, int iSpecies, int iCharge)
  {
    float scaling_factor = nuclei_spectra::charges[iSpecies] * cfgMomentumScalingBetheBloch->get(iSpecies, iCharge) / nuclei_spectra::masses[iSpecies];

    float expTPCSignal = o2::tpc::BetheBlochAleph(static_cast<float>(candidate.tpcInnerParam() * scaling_factor),
                                                  cfgBetheBlochParams->get(iSpecies, 0u),
                                                  cfgBetheBlochParams->get(iSpecies, 1u),
                                                  cfgBetheBlochParams->get(iSpecies, 2u),
                                                  cfgBetheBlochParams->get(iSpecies, 3u),
                                                  cfgBetheBlochParams->get(iSpecies, 4u));
    float resolutionTPC{expTPCSignal * cfgBetheBlochParams->get(iSpecies, 5u)};
    return static_cast<float>((candidate.tpcSignal() - expTPCSignal) / resolutionTPC);
  }

  /// @brief Get average ITS cluster size
  /// @tparam T type for the track
  /// @param track
  /// @return average cluster size in ITS
  template <class T>
  float getITSClSize(T const& track)
  {
    float sum{0.f};
    int nClus = 0;
    for (int iL{0}; iL < 6; ++iL) {
      auto size = (track.itsClusterSizes() >> (iL * 4)) & 0xf;
      if (size > 0) {
        nClus++;
        sum += size;
      }
    }
    return sum / nClus;
  }

  /// @brief Get average ITS cluster size
  /// @tparam T type for the track
  /// @param track
  /// @return true if the candidates passes all the selections
  template <class T>
  bool selectTrack(T const& track)
  {
    if (track.tpcNCls() < cfgNclusTPCcut)
      return false;
    if (track.dcaxy() > cfgDCAxyCut)
      return false;
    if (track.dcaz() > cfgDCAzCut)
      return false;
    if (getITSClSize(track) < cfgItsClusSizeCut)
      return false;
    return true;
  }

  /// \brief Get the centrality with the selected detector
  /// \param track track with the centrality information
  /// \return centrality expressed in percentage
  float getRefCentrality(TrackWithFlow const& track)
  {
    float cent = -999.;
    switch (cfgCentDetector) {
      case flow::kDetector::kFV0A:
        cent = track.centFV0A();
        break;
      case flow::kDetector::kFT0M:
        cent = track.centFT0M();
        break;
      case flow::kDetector::kFT0A:
        cent = track.centFT0A();
        break;
      case flow::kDetector::kFT0C:
        cent = track.centFT0C();
        break;
      default:
        LOG(warning) << "Centrality estimator not valid. Possible values are V0A, T0M, T0A, T0C. Fallback to V0A";
        cent = track.centFV0A();
        break;
    }
    return cent;
  }

  /// \brief Get the Q vector etimated with a particular detector
  /// \param track track with the Q vector information
  /// \return Q vector in format {x, y}
  std::vector<float> getQvec(TrackWithFlow const& track, int detector)
  {
    float xQvec = -999.;
    float yQvec = -999.;
    float amplQvec = -999.;
    switch (detector) {
      case flow::kDetector::kFV0A:
        xQvec = track.xQvecFV0A();
        yQvec = track.yQvecFV0A();
        amplQvec = track.amplQvecFV0A();
        break;
      case flow::kDetector::kFT0M:
        xQvec = track.xQvecFT0M();
        yQvec = track.yQvecFT0M();
        amplQvec = track.amplQvecFT0M();
        break;
      case flow::kDetector::kFT0A:
        xQvec = track.xQvecFT0A();
        yQvec = track.yQvecFT0A();
        amplQvec = track.amplQvecFT0A();
        break;
      case flow::kDetector::kFT0C:
        xQvec = track.xQvecFT0M();
        yQvec = track.yQvecFT0M();
        amplQvec = track.amplQvecFT0M();
        break;
      case flow::kDetector::kTPCpos:
        xQvec = track.xQvecTPCpos();
        yQvec = track.yQvecTPCpos();
        amplQvec = track.amplQvecTPCpos();
        break;
      case flow::kDetector::kTPCneg:
        xQvec = track.xQvecTPCneg();
        yQvec = track.yQvecTPCneg();
        amplQvec = track.amplQvecTPCneg();
        break;
      default:
        LOG(warning) << "Q vector estimator not valid. Please choose between FV0A, FT0M, FT0A, FT0C, TPC Pos, TPC Neg. Fallback to FV0A";
        xQvec = track.xQvecFV0A();
        yQvec = track.yQvecFV0A();
        amplQvec = track.amplQvecFV0A();
        break;
    }
    return {xQvec, yQvec, amplQvec};
  }

  void init(o2::framework::InitContext&)
  {
    const AxisSpec nSigmaTPCHe3Axis{nSigmaBins, "n#sigma_{TPC}({}^{3}He)"};
    const AxisSpec ptAxis{ptBins, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec centAxis{centBins, "centrality(%)"};
    const AxisSpec FT0AspAxis{spBins, "#hat{u}_{2} #upoint #vec{Q}_{2}^{FT0A}"};
    const AxisSpec FT0CspAxis{spBins, "#hat{u}_{2} #upoint #vec{Q}_{2}^{FT0C}"};
    const AxisSpec FV0AspAxis{spBins, "#hat{u}_{2} #upoint #vec{Q}_{2}^{FV0A}"};

    histos.add("hSpFT0AvsNsigmaHe3VsPtvsCent", "", HistType::kTHnSparseF, {FT0AspAxis, nSigmaTPCHe3Axis, ptAxis, centAxis});
    histos.add("hSpFT0CvsNsigmaHe3VsPtvsCent", "", HistType::kTHnSparseF, {FT0CspAxis, nSigmaTPCHe3Axis, ptAxis, centAxis});
    histos.add("hSpFV0AvsNsigmaHe3VsPtvsCent", "", HistType::kTHnSparseF, {FV0AspAxis, nSigmaTPCHe3Axis, ptAxis, centAxis});
  }

  void process(TracksWithFlow const& tracks)
  {
    for (auto& track : tracks) {

      if (!selectTrack(track))
        return;

      // Select the centrality value to be stored in the output histograms
      float ref_cent = getRefCentrality(track);

      // Get event plane with T0A
      std::vector<float> qVecFT0A = getQvec(track, flow::kDetector::kFT0A);
      float xQvecFT0A = qVecFT0A[0];
      float yQvecFT0A = qVecFT0A[1];

      // Get event plane with T0C
      std::vector<float> qVecFT0C = getQvec(track, flow::kDetector::kFT0C);
      float xQvecFT0C = qVecFT0C[0];
      float yQvecFT0C = qVecFT0C[1];
      // Get event plane with V0A
      std::vector<float> qVecFV0A = getQvec(track, flow::kDetector::kFV0A);
      float xQvecFV0A = qVecFV0A[0];
      float yQvecFV0A = qVecFV0A[1];

      // Get candidate vector
      float xCandVec = TMath::Cos(cfgHarmonic * track.phi());
      float yCandVec = TMath::Sin(cfgHarmonic * track.phi());

      // Get scalar products for different event planes
      float spFT0A = xCandVec * xQvecFT0A + yCandVec * yQvecFT0A;
      float spFT0C = xCandVec * xQvecFT0C + yCandVec * yQvecFT0C;
      float spFV0A = xCandVec * xQvecFV0A + yCandVec * yQvecFV0A;

      // Get candidate info
      int iCharge = (track.flags() & nuclei_spectra::CandBits::kPositive) ? 0 : 1;
      float nSigmaTPC = getNSigmaTPC(track, cfgSpecies, iCharge);
      float pt = track.pt() * nuclei_spectra::charges[cfgSpecies];

      // Fill relevant histograms
      histos.fill(HIST("hSpFT0AvsNsigmaHe3VsPtvsCent"), spFT0A, nSigmaTPC, pt, ref_cent);
      histos.fill(HIST("hSpFT0CvsNsigmaHe3VsPtvsCent"), spFT0C, nSigmaTPC, pt, ref_cent);
      histos.fill(HIST("hSpFV0AvsNsigmaHe3VsPtvsCent"), spFV0A, nSigmaTPC, pt, ref_cent);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<nucleiFlow>(cfgc, TaskName{"nucleiFlow"})};
}
