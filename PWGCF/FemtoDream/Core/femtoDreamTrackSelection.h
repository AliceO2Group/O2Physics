// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file FemtoDreamTrackCuts.h
/// \brief Definition of the FemtoDreamTrackCuts
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch

#ifndef PWGCF_FEMTODREAM_CORE_FEMTODREAMTRACKSELECTION_H_
#define PWGCF_FEMTODREAM_CORE_FEMTODREAMTRACKSELECTION_H_

#include <string>
#include <vector>
#include <cmath>
#include <iostream>

#include "PWGCF/DataModel/FemtoDerived.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "PWGCF/FemtoDream/Core/femtoDreamObjectSelection.h"
#include "ReconstructionDataFormats/PID.h"
#include "Framework/HistogramRegistry.h"

using namespace o2::framework;

namespace o2::analysis::femtoDream
{
namespace femtoDreamTrackSelection
{
/// The different selections this task is capable of doing
enum TrackSel { kSign,         ///< Sign of the track
                kpTMin,        ///< Min. p_T (GeV/c)
                kpTMax,        ///< Max. p_T (GeV/c)
                kEtaMax,       ///< Max. |eta|
                kTPCnClsMin,   ///< Min. number of TPC clusters
                kTPCfClsMin,   ///< Min. fraction of crossed rows/findable TPC clusters
                kTPCcRowsMin,  ///< Min. number of crossed TPC rows
                kTPCsClsMax,   ///< Max. number of shared TPC clusters
                kITSnClsMin,   ///< Min. number of ITS clusters
                kITSnClsIbMin, ///< Min. number of ITS clusters in the inner barrel
                kDCAxyMax,     ///< Max. DCA_xy (cm)
                kDCAzMax,      ///< Max. DCA_z (cm)
                kDCAMin,       ///< Min. DCA_xyz (cm)
                kPIDnSigmaMax  ///< Max. |n_sigma| for PID
};

enum TrackContainerPosition {
  kCuts,
  kPID
}; /// Position in the full track cut container

} // namespace femtoDreamTrackSelection

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
class FemtoDreamTrackSelection : public FemtoDreamObjectSelection<float, femtoDreamTrackSelection::TrackSel>
{
 public:
  FemtoDreamTrackSelection() : nRejectNotPropagatedTracks(false),
                               nPtMinSel(0),
                               nPtMaxSel(0),
                               nEtaSel(0),
                               nTPCnMinSel(0),
                               nTPCfMinSel(0),
                               nTPCcMinSel(0),
                               nTPCsMaxSel(0),
                               nITScMinSel(0),
                               nITScIbMinSel(0),
                               nDCAxyMaxSel(0),
                               nDCAzMaxSel(0),
                               nDCAMinSel(0),
                               nPIDnSigmaSel(0),
                               pTMin(9999999.),
                               pTMax(-9999999.),
                               etaMax(-9999999.),
                               nClsMin(9999999.),
                               fClsMin(9999999.),
                               cTPCMin(9999999.),
                               sTPCMax(-9999999.),
                               dcaXYMax(-9999999.),
                               dcaZMax(-9999999.),
                               dcaMin(9999999.),
                               nSigmaPIDMax(9999999.),
                               nSigmaPIDOffsetTPC(0.),
                               nSigmaPIDOffsetTOF(0.) {}

  /// Initializes histograms for the task
  /// \tparam part Type of the particle for proper naming of the folders for QA
  /// \tparam tracktype Type of track (track, positive child, negative child) for proper naming of the folders for QA
  /// \tparam cutContainerType Data type of the bit-wise container for the selections
  /// \param registry HistogramRegistry for QA output
  template <o2::aod::femtodreamparticle::ParticleType part, o2::aod::femtodreamparticle::TrackType tracktype, typename cutContainerType>
  void init(HistogramRegistry* QAregistry, HistogramRegistry* Registry);

  /// Passes the species to the task for which PID needs to be stored
  /// \tparam T Data type of the configurable passed to the functions
  /// \param pids Configurable with the species
  template <typename T>
  void setPIDSpecies(T& pids)
  {
    std::vector<int> tmpPids = pids; /// necessary due to some features of the configurable
    for (o2::track::PID pid : tmpPids) {
      mPIDspecies.push_back(pid);
    }
  }

  /// Computes the n_sigma for a track and a particle-type hypothesis in the TPC
  /// \tparam T Data type of the track
  /// \param track Track for which PID is evaluated
  /// \param pid Particle species for which PID is evaluated
  /// \return Value of n_{sigma, TPC}
  template <typename T>
  auto getNsigmaTPC(T const& track, o2::track::PID pid);

  /// Computes the n_sigma for a track and a particle-type hypothesis in the TOF
  /// \tparam T Data type of the track
  /// \param track Track for which PID is evaluated
  /// \param pid Particle species for which PID is evaluated
  /// \return Value of n_{sigma, TOF}
  template <typename T>
  auto getNsigmaTOF(T const& track, o2::track::PID pid);

  /// Computes the n_sigma for a track and a particle-type hypothesis in the ITS
  /// \tparam T Data type of the track
  /// \param track Track for which PID is evaluated
  /// \param pid Particle species for which PID is evaluated
  /// \return Value of n_{sigma, ITS}
  template <typename T>
  auto getNsigmaITS(T const& track, o2::track::PID pid);

  /// Checks whether the most open combination of all selection criteria is fulfilled
  /// \tparam T Data type of the track
  /// \param track Track
  /// \return Whether the most open combination of all selection criteria is fulfilled
  template <typename T>
  bool isSelectedMinimal(T const& track);

  /// Obtain the bit-wise container for the selections
  /// Pt, eta and dca are not necessarily taken from the track table. For example, for V0 daughters they are recaluated and stored in the V0 table
  /// \todo For the moment, PID is separated from the other selections, hence instead of a single value an std::array of size two is returned
  /// \tparam cutContainerType Data type of the bit-wise container for the selections
  /// \tparam T Data type of the track
  /// \param track Track
  /// \param Pt pt of the track
  /// \param Eta eta of the track
  /// \param Dca dca of the track with respect to primary vertex
  /// \return The bit-wise container for the selections, separately with all selection criteria, and the PID
  template <bool useItsPid = false, typename cutContainerType, typename T, typename R>
  std::array<cutContainerType, 2> getCutContainer(T const& track, R Pt, R Eta, R Dcaxy);

  /// Some basic QA histograms
  /// \tparam part Type of the particle for proper naming of the folders for QA
  /// \tparam tracktype Type of track (track, positive child, negative child) for proper naming of the folders for QA
  /// \tparam T Data type of the track
  /// \param track Track
  template <o2::aod::femtodreamparticle::ParticleType part, o2::aod::femtodreamparticle::TrackType tracktype, bool isHF = false, int cutstage = 1, typename T>
  void fillQA(T const& track);

  /// Helper function to obtain the name of a given selection criterion for consistent naming of the configurables
  /// \param iSel Track selection variable to be examined
  /// \param prefix Additional prefix for the name of the configurable
  /// \param suffix Additional suffix for the name of the configurable
  static std::string getSelectionName(femtoDreamTrackSelection::TrackSel iSel, std::string_view prefix = "", std::string_view suffix = "")
  {
    std::string outString = static_cast<std::string>(prefix);
    outString += static_cast<std::string>(mSelectionNames[iSel]);
    outString += suffix;
    return outString;
  }

  /// Helper function to obtain the index of a given selection variable for consistent naming of the configurables
  /// \param obs Track selection variable (together with prefix) got from file
  /// \param prefix Additional prefix for the output of the configurable
  static int findSelectionIndex(const std::string_view& obs, std::string_view prefix = "")
  {
    for (int index = 0; index < kNtrackSelection; index++) {
      std::string comp = static_cast<std::string>(prefix) + static_cast<std::string>(mSelectionNames[index]);
      std::string_view cmp{comp};
      if (obs.compare(cmp) == 0)
        return index;
    }
    return -1;
  }

  /// Helper function to obtain the type of a given selection variable for consistent naming of the configurables
  /// \param iSel Track selection variable whose type is returned
  static femtoDreamSelection::SelectionType getSelectionType(femtoDreamTrackSelection::TrackSel iSel)
  {
    return mSelectionTypes[iSel];
  }

  /// Helper function to obtain the helper string of a given selection criterion for consistent description of the configurables
  /// \param iSel Track selection variable to be examined
  /// \param prefix Additional prefix for the output of the configurable
  static std::string getSelectionHelper(femtoDreamTrackSelection::TrackSel iSel, std::string_view prefix = "")
  {
    std::string outString = static_cast<std::string>(prefix);
    outString += static_cast<std::string>(mSelectionHelper[iSel]);
    return outString;
  }

  int getSigmaPIDMax()
  {
    return nSigmaPIDMax;
  }

  void setRejectNotPropagatedTracks(bool reject)
  {
    nRejectNotPropagatedTracks = reject;
  }
  void setnSigmaPIDOffset(float offsetTPC, float offsetTOF)
  {
    nSigmaPIDOffsetTPC = offsetTPC;
    nSigmaPIDOffsetTOF = offsetTOF;
  }

 private:
  bool nRejectNotPropagatedTracks;
  int nPtMinSel;
  int nPtMaxSel;
  int nEtaSel;
  int nTPCnMinSel;
  int nTPCfMinSel;
  int nTPCcMinSel;
  int nTPCsMaxSel;
  int nITScMinSel;
  int nITScIbMinSel;
  int nDCAxyMaxSel;
  int nDCAzMaxSel;
  int nDCAMinSel;
  int nPIDnSigmaSel;
  float pTMin;
  float pTMax;
  float etaMax;
  float nClsMin;
  float fClsMin;
  float cTPCMin;
  float sTPCMax;
  float nITSclsMin;
  float nITSclsIbMin;
  float dcaXYMax;
  float dcaZMax;
  float dcaMin;
  float nSigmaPIDMax;
  float nSigmaPIDOffsetTPC;
  float nSigmaPIDOffsetTOF;
  std::vector<o2::track::PID> mPIDspecies; ///< All the particle species for which the n_sigma values need to be stored
  static constexpr int kNtrackSelection = 14;
  static constexpr std::string_view mSelectionNames[kNtrackSelection] = {"Sign",
                                                                         "PtMin",
                                                                         "PtMax",
                                                                         "EtaMax",
                                                                         "TPCnClsMin",
                                                                         "TPCfClsMin",
                                                                         "TPCcRowsMin",
                                                                         "TPCsClsMax",
                                                                         "ITSnClsMin",
                                                                         "ITSnClsIbMin",
                                                                         "DCAxyMax",
                                                                         "DCAzMax",
                                                                         "DCAMin",
                                                                         "PIDnSigmaMax"}; ///< Name of the different selections

  static constexpr femtoDreamSelection::SelectionType mSelectionTypes[kNtrackSelection]{femtoDreamSelection::kEqual,
                                                                                        femtoDreamSelection::kLowerLimit,
                                                                                        femtoDreamSelection::kUpperLimit,
                                                                                        femtoDreamSelection::kAbsUpperLimit,
                                                                                        femtoDreamSelection::kLowerLimit,
                                                                                        femtoDreamSelection::kLowerLimit,
                                                                                        femtoDreamSelection::kLowerLimit,
                                                                                        femtoDreamSelection::kUpperLimit,
                                                                                        femtoDreamSelection::kLowerLimit,
                                                                                        femtoDreamSelection::kLowerLimit,
                                                                                        femtoDreamSelection::kAbsUpperLimit,
                                                                                        femtoDreamSelection::kAbsUpperLimit,
                                                                                        femtoDreamSelection::kAbsUpperLimit,  // <-----TODO this should be a lower limit, no??
                                                                                        femtoDreamSelection::kAbsUpperLimit}; ///< Map to match a variable with its type

  static constexpr std::string_view mSelectionHelper[kNtrackSelection] = {"Sign of the track",
                                                                          "Minimal pT (GeV/c)",
                                                                          "Maximal pT (GeV/c)",
                                                                          "Maximal eta",
                                                                          "Minimum number of TPC clusters",
                                                                          "Minimum fraction of crossed rows/findable clusters",
                                                                          "Minimum number of crossed TPC rows",
                                                                          "Maximal number of shared TPC cluster",
                                                                          "Minimum number of ITS clusters",
                                                                          "Minimum number of ITS clusters in the inner barrel",
                                                                          "Maximal DCA_xy (cm)",
                                                                          "Maximal DCA_z (cm)",
                                                                          "Minimal DCA (cm)",
                                                                          "Maximal PID (nSigma)"}; ///< Helper information for the different selections
  static constexpr int kNcutStages = 2;
  static constexpr std::string_view mCutStage[kNcutStages] = {"BeforeSel", "AfterSel"};
}; // namespace femtoDream

template <o2::aod::femtodreamparticle::ParticleType part, o2::aod::femtodreamparticle::TrackType tracktype, typename cutContainerType>
void FemtoDreamTrackSelection::init(HistogramRegistry* QAregistry, HistogramRegistry* Registry)
{
  if (QAregistry && Registry) {
    mHistogramRegistry = Registry;
    mQAHistogramRegistry = QAregistry;
    std::string folderName = static_cast<std::string>(o2::aod::femtodreamparticle::ParticleTypeName[part]) + "/" + static_cast<std::string>(o2::aod::femtodreamparticle::TrackTypeName[tracktype]);

    /// check whether the number of selection exceeds the bitmap size
    unsigned int nSelections = getNSelections() - getNSelections(femtoDreamTrackSelection::kPIDnSigmaMax);
    if (nSelections > 8 * sizeof(cutContainerType)) {
      LOG(fatal) << "FemtoDreamTrackCuts: Number of selections too large for your container - quitting!";
    }

    for (int istage = 0; istage < kNcutStages; istage++) {
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hPt").c_str(), "; #it{p}_{T} (GeV/#it{c}); Entries", kTH1F, {{240, 0, 6}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hEta").c_str(), "; #eta; Entries", kTH1F, {{200, -1.5, 1.5}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hPhi").c_str(), "; #phi; Entries", kTH1F, {{200, 0, 2. * M_PI}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hTPCfindable").c_str(), "; TPC findable clusters; Entries", kTH1F, {{163, -0.5, 162.5}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hTPCfound").c_str(), "; TPC found clusters; Entries", kTH1F, {{163, -0.5, 162.5}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hTPCcrossedOverFindalbe").c_str(), "; TPC ratio findable; Entries", kTH1F, {{100, 0.5, 1.5}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hTPCcrossedRows").c_str(), "; TPC crossed rows; Entries", kTH1F, {{163, 0, 163}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hTPCfindableVsCrossed").c_str(), ";TPC findable clusters ; TPC crossed rows;", kTH2F, {{163, 0, 163}, {163, 0, 163}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hTPCshared").c_str(), "; TPC shared clusters; Entries", kTH1F, {{163, -0.5, 162.5}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hITSclusters").c_str(), "; ITS clusters; Entries", kTH1F, {{10, -0.5, 9.5}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hITSclustersIB").c_str(), "; ITS clusters in IB; Entries", kTH1F, {{10, -0.5, 9.5}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hDCAxy").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {{100, 0, 10}, {500, -5, 5}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hDCAz").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{z} (cm)", kTH2F, {{100, 0, 10}, {500, -5, 5}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hDCA").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA (cm)", kTH2F, {{100, 0, 10}, {301, 0., 1.5}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/hTPCdEdX").c_str(), "; #it{p} (GeV/#it{c}); TPC Signal", kTH2F, {{100, 0, 10}, {1000, 0, 1000}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/nSigmaTPC_el").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{e}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/nSigmaTPC_pi").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{#pi}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/nSigmaTPC_K").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{K}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/nSigmaTPC_p").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{p}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/nSigmaTPC_d").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{d}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/nSigmaTOF_el").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{e}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/nSigmaTOF_pi").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{#pi}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/nSigmaTOF_K").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{K}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/nSigmaTOF_p").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{p}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/nSigmaTOF_d").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{d}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/nSigmaComb_el").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{e}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/nSigmaComb_pi").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{#pi}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/nSigmaComb_K").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{K}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/nSigmaComb_p").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{p}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
      mQAHistogramRegistry->add((folderName + "/" + static_cast<std::string>(mCutStage[istage]) + "/nSigmaComb_d").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{d}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    }
  }

  /// set cuts
  nPtMinSel = getNSelections(femtoDreamTrackSelection::kpTMin);
  nPtMaxSel = getNSelections(femtoDreamTrackSelection::kpTMax);
  nEtaSel = getNSelections(femtoDreamTrackSelection::kEtaMax);
  nTPCnMinSel = getNSelections(femtoDreamTrackSelection::kTPCnClsMin);
  nTPCfMinSel = getNSelections(femtoDreamTrackSelection::kTPCfClsMin);
  nTPCcMinSel = getNSelections(femtoDreamTrackSelection::kTPCcRowsMin);
  nTPCsMaxSel = getNSelections(femtoDreamTrackSelection::kTPCsClsMax);
  nITScMinSel = getNSelections(femtoDreamTrackSelection::kITSnClsMin);
  nITScIbMinSel = getNSelections(femtoDreamTrackSelection::kITSnClsIbMin);
  nDCAxyMaxSel = getNSelections(femtoDreamTrackSelection::kDCAxyMax);
  nDCAzMaxSel = getNSelections(femtoDreamTrackSelection::kDCAzMax);
  nDCAMinSel = getNSelections(femtoDreamTrackSelection::kDCAMin);
  nPIDnSigmaSel = getNSelections(femtoDreamTrackSelection::kPIDnSigmaMax);

  pTMin = getMinimalSelection(femtoDreamTrackSelection::kpTMin, femtoDreamSelection::kLowerLimit);
  pTMax = getMinimalSelection(femtoDreamTrackSelection::kpTMax, femtoDreamSelection::kUpperLimit);
  etaMax = getMinimalSelection(femtoDreamTrackSelection::kEtaMax, femtoDreamSelection::kAbsUpperLimit);
  nClsMin = getMinimalSelection(femtoDreamTrackSelection::kTPCnClsMin, femtoDreamSelection::kLowerLimit);
  fClsMin = getMinimalSelection(femtoDreamTrackSelection::kTPCfClsMin, femtoDreamSelection::kLowerLimit);
  cTPCMin = getMinimalSelection(femtoDreamTrackSelection::kTPCcRowsMin, femtoDreamSelection::kLowerLimit);
  sTPCMax = getMinimalSelection(femtoDreamTrackSelection::kTPCsClsMax, femtoDreamSelection::kUpperLimit);
  nITSclsMin = getMinimalSelection(femtoDreamTrackSelection::kITSnClsMin, femtoDreamSelection::kLowerLimit);
  nITSclsIbMin = getMinimalSelection(femtoDreamTrackSelection::kITSnClsIbMin, femtoDreamSelection::kLowerLimit);
  dcaXYMax = getMinimalSelection(femtoDreamTrackSelection::kDCAxyMax, femtoDreamSelection::kAbsUpperLimit);
  dcaZMax = getMinimalSelection(femtoDreamTrackSelection::kDCAzMax, femtoDreamSelection::kAbsUpperLimit);
  dcaMin = getMinimalSelection(femtoDreamTrackSelection::kDCAMin, femtoDreamSelection::kAbsLowerLimit);
  nSigmaPIDMax = getMinimalSelection(femtoDreamTrackSelection::kPIDnSigmaMax, femtoDreamSelection::kAbsUpperLimit);
}

template <typename T>
auto FemtoDreamTrackSelection::getNsigmaTPC(T const& track, o2::track::PID pid)
{
  return o2::aod::pidutils::tpcNSigma(pid, track);
}

template <typename T>
auto FemtoDreamTrackSelection::getNsigmaTOF(T const& track, o2::track::PID pid)
{
  /// skip tracks without TOF signal
  if (!track.hasTOF()) {
    return 999.f;
  }
  return o2::aod::pidutils::tofNSigma(pid, track);
}

template <typename T>
auto FemtoDreamTrackSelection::getNsigmaITS(T const& track, o2::track::PID pid)
{
  if (pid == o2::track::PID::Electron) {
    return track.itsNSigmaEl();
  } else if (pid == o2::track::PID::Pion) {
    return track.itsNSigmaPi();
  } else if (pid == o2::track::PID::Kaon) {
    return track.itsNSigmaKa();
  } else if (pid == o2::track::PID::Proton) {
    return track.itsNSigmaPr();
  } else if (pid == o2::track::PID::Deuteron) {
    return track.itsNSigmaDe();
  } else if (pid == o2::track::PID::Triton) {
    return track.itsNSigmaTr();
  } else if (pid == o2::track::PID::Helium3) {
    return track.itsNSigmaHe();
  }
  // if nothing matched, return default value
  return -999.f;
}

template <typename T>
bool FemtoDreamTrackSelection::isSelectedMinimal(T const& track)
{
  const auto pT = track.pt();
  const auto eta = track.eta();
  const auto tpcNClsF = track.tpcNClsFound();
  const auto tpcRClsC = track.tpcCrossedRowsOverFindableCls();
  const auto tpcNClsC = track.tpcNClsCrossedRows();
  const auto tpcNClsS = track.tpcNClsShared();
  const auto itsNCls = track.itsNCls();
  const auto itsNClsIB = track.itsNClsInnerBarrel();
  const auto dcaXY = track.dcaXY();
  const auto dcaZ = track.dcaZ();
  const auto dca = track.dcaXY(); // Accordingly to FemtoDream in AliPhysics  as well as LF analysis,
                                  // only dcaXY should be checked; NOT std::sqrt(pow(dcaXY, 2.) + pow(dcaZ, 2.))
  std::vector<float> pidTPC, pidTOF;
  for (auto it : mPIDspecies) {
    pidTPC.push_back(getNsigmaTPC(track, it));
    pidTOF.push_back(getNsigmaTOF(track, it));
  }

  if (nPtMinSel > 0 && pT < pTMin) {
    return false;
  }
  if (nPtMaxSel > 0 && pT > pTMax) {
    return false;
  }
  if (nEtaSel > 0 && std::fabs(eta) > etaMax) {
    return false;
  }
  if (nTPCnMinSel > 0 && tpcNClsF < nClsMin) {
    return false;
  }
  if (nTPCfMinSel > 0 && tpcRClsC < fClsMin) {
    return false;
  }
  if (nTPCcMinSel > 0 && tpcNClsC < cTPCMin) {
    return false;
  }
  if (nTPCsMaxSel > 0 && tpcNClsS > sTPCMax) {
    return false;
  }
  if (nITScMinSel > 0 && itsNCls < nITSclsMin) {
    return false;
  }
  if (nITScIbMinSel > 0 && itsNClsIB < nITSclsIbMin) {
    return false;
  }
  if (nDCAxyMaxSel > 0 && std::fabs(dcaXY) > dcaXYMax) {
    return false;
  }
  if (nDCAzMaxSel > 0 && std::fabs(dcaZ) > dcaZMax) {
    return false;
  }
  if (nDCAMinSel > 0 && std::fabs(dca) < dcaMin) {
    return false;
  }
  if (nRejectNotPropagatedTracks && std::fabs(dca) > 1e3) {
    return false;
  }

  if (nPIDnSigmaSel > 0) {
    bool isFulfilled = false;
    for (size_t i = 0; i < pidTPC.size(); ++i) {
      auto pidTPCVal = pidTPC.at(i);
      if (std::fabs(pidTPCVal - nSigmaPIDOffsetTPC) < nSigmaPIDMax) {
        isFulfilled = true;
      }
    }
    if (!isFulfilled) {
      return isFulfilled;
    }
  }
  return true;
}

template <bool useItsPid, typename cutContainerType, typename T, typename R>
std::array<cutContainerType, 2> FemtoDreamTrackSelection::getCutContainer(T const& track, R Pt, R Eta, R Dca)
{
  cutContainerType output = 0;
  size_t counter = 0;
  cutContainerType outputPID = 0;
  const auto sign = track.sign();
  const auto pt = Pt;
  const auto eta = Eta;
  const auto tpcNClsF = track.tpcNClsFound();
  const auto tpcRClsC = track.tpcCrossedRowsOverFindableCls();
  const auto tpcNClsC = track.tpcNClsCrossedRows();
  const auto tpcNClsS = track.tpcNClsShared();
  const auto itsNCls = track.itsNCls();
  const auto itsNClsIB = track.itsNClsInnerBarrel();
  const auto dcaXY = track.dcaXY();
  const auto dcaZ = track.dcaZ();
  const auto dca = Dca;

  std::vector<float> pidTPC, pidTOF, pidITS;
  for (auto it : mPIDspecies) {
    pidTPC.push_back(getNsigmaTPC(track, it));
    pidTOF.push_back(getNsigmaTOF(track, it));
    if constexpr (useItsPid) {
      pidITS.push_back(getNsigmaITS(track, it));
    }
  }

  float observable = 0.;
  for (auto& sel : mSelections) {
    const auto selVariable = sel.getSelectionVariable();
    if (selVariable == femtoDreamTrackSelection::kPIDnSigmaMax) {
      /// PID needsgetNsigmaITSto be handled a bit differently since we may need more than one species
      for (size_t i = 0; i < mPIDspecies.size(); ++i) {
        auto pidTPCVal = pidTPC.at(i) - nSigmaPIDOffsetTPC;
        auto pidTOFVal = pidTOF.at(i) - nSigmaPIDOffsetTOF;
        auto pidComb = std::sqrt(pidTPCVal * pidTPCVal + pidTOFVal * pidTOFVal);
        sel.checkSelectionSetBitPID(pidTPCVal, outputPID);
        sel.checkSelectionSetBitPID(pidComb, outputPID);
        if constexpr (useItsPid) {
          auto pidITSVal = pidITS.at(i);
          sel.checkSelectionSetBitPID(pidITSVal, outputPID);
        }
      }
    } else {
      /// for the rest it's all the same
      switch (selVariable) {
        case (femtoDreamTrackSelection::kSign):
          observable = sign;
          break;
        case (femtoDreamTrackSelection::kpTMin):
        case (femtoDreamTrackSelection::kpTMax):
          observable = pt;
          break;
        case (femtoDreamTrackSelection::kEtaMax):
          observable = eta;
          break;
        case (femtoDreamTrackSelection::kTPCnClsMin):
          observable = tpcNClsF;
          break;
        case (femtoDreamTrackSelection::kTPCfClsMin):
          observable = tpcRClsC;
          break;
        case (femtoDreamTrackSelection::kTPCcRowsMin):
          observable = tpcNClsC;
          break;
        case (femtoDreamTrackSelection::kTPCsClsMax):
          observable = tpcNClsS;
          break;
        case (femtoDreamTrackSelection::kITSnClsMin):
          observable = itsNCls;
          break;
        case (femtoDreamTrackSelection::kITSnClsIbMin):
          observable = itsNClsIB;
          break;
        case (femtoDreamTrackSelection::kDCAxyMax):
          observable = dcaXY;
          break;
        case (femtoDreamTrackSelection::kDCAzMax):
          observable = dcaZ;
          break;
        case (femtoDreamTrackSelection::kDCAMin):
          observable = dca;
          break;
        case (femtoDreamTrackSelection::kPIDnSigmaMax):
          break;
      }
      sel.checkSelectionSetBit(observable, output, counter, mHistogramRegistry);
    }
  }
  return {output, outputPID};
}

template <o2::aod::femtodreamparticle::ParticleType part, o2::aod::femtodreamparticle::TrackType tracktype, bool isHF, int cutstage, typename T>
void FemtoDreamTrackSelection::fillQA(T const& track)
{
  if (mQAHistogramRegistry) {
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hPt"), track.pt());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hEta"), track.eta());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hPhi"), track.phi());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hTPCfindable"), track.tpcNClsFindable());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hTPCfound"), track.tpcNClsFound());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hTPCcrossedOverFindalbe"), track.tpcCrossedRowsOverFindableCls());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hTPCcrossedRows"), track.tpcNClsCrossedRows());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hTPCfindableVsCrossed"), track.tpcNClsFindable(), track.tpcNClsCrossedRows());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hTPCshared"), track.tpcNClsShared());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hITSclusters"), track.itsNCls());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hITSclustersIB"), track.itsNClsInnerBarrel());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hDCAxy"), track.pt(), track.dcaXY());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hDCAz"), track.pt(), track.dcaZ());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hDCA"), track.pt(), std::sqrt(pow(track.dcaXY(), 2.) + pow(track.dcaZ(), 2.)));
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/hTPCdEdX"), track.p(), track.tpcSignal());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/nSigmaTPC_pi"), track.p(), track.tpcNSigmaPi());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/nSigmaTPC_K"), track.p(), track.tpcNSigmaKa());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/nSigmaTPC_p"), track.p(), track.tpcNSigmaPr());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/nSigmaTOF_pi"), track.p(), track.tofNSigmaPi());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/nSigmaTOF_K"), track.p(), track.tofNSigmaKa());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/nSigmaTOF_p"), track.p(), track.tofNSigmaPr());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/nSigmaComb_pi"), track.p(), std::sqrt(track.tpcNSigmaPi() * track.tpcNSigmaPi() + track.tofNSigmaPi() * track.tofNSigmaPi()));
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/nSigmaComb_K"), track.p(), std::sqrt(track.tpcNSigmaKa() * track.tpcNSigmaKa() + track.tofNSigmaKa() * track.tofNSigmaKa()));
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/nSigmaComb_p"), track.p(), std::sqrt(track.tpcNSigmaPr() * track.tpcNSigmaPr() + track.tofNSigmaPr() * track.tofNSigmaPr()));
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/nSigmaTPC_d"), track.p(), track.tpcNSigmaDe());
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/nSigmaComb_d"), track.p(), std::sqrt(track.tpcNSigmaDe() * track.tpcNSigmaDe() + track.tofNSigmaDe() * track.tofNSigmaDe()));
    mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/nSigmaTOF_d"), track.p(), track.tofNSigmaDe());

    if constexpr (!isHF) {
      mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/nSigmaComb_el"), track.p(), std::sqrt(track.tpcNSigmaEl() * track.tpcNSigmaEl() + track.tofNSigmaEl() * track.tofNSigmaEl()));
      mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/nSigmaTOF_el"), track.p(), track.tofNSigmaEl());
      mQAHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtodreamparticle::TrackTypeName[tracktype]) + HIST("/") + HIST(mCutStage[cutstage]) + HIST("/nSigmaTPC_el"), track.p(), track.tpcNSigmaEl());
    }
  }
}

} // namespace o2::analysis::femtoDream

#endif // PWGCF_FEMTODREAM_CORE_FEMTODREAMTRACKSELECTION_H_
