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

/// \file FemtoWorldTrackCuts.h
/// \brief Definition of the FemtoWorldTrackCuts
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Luca Barioglio, TU München, luca.barioglio@cern.ch
/// \author Zuzanna Chochulska, WUT Warsaw, zchochul@cern.ch

#ifndef FEMTOWORLDTRACKSELECTION_H_
#define FEMTOWORLDTRACKSELECTION_H_

#include "PWGCF/FemtoWorld/Core/FemtoWorldObjectSelection.h"
#include "PWGCF/FemtoWorld/DataModel/FemtoWorldDerived.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/PID.h"

#include <cmath>
#include <iostream>

using namespace o2::framework;

namespace o2::analysis::femtoWorld
{
namespace femtoWorldTrackSelection
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

} // namespace femtoWorldTrackSelection

/// \class FemtoWorldTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
class FemtoWorldTrackSelection : public FemtoWorldObjectSelection<float, femtoWorldTrackSelection::TrackSel>
{
 public:
  FemtoWorldTrackSelection() : nRejectNotPropagatedTracks(false),
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
                               nSigmaPIDMax(9999999.) {};

  /// Initializes histograms for the task
  /// \tparam part Type of the particle for proper naming of the folders for QA
  /// \tparam tracktype Type of track (track, positive child, negative child) for proper naming of the folders for QA
  /// \tparam cutContainerType Data type of the bit-wise container for the selections
  /// \param registry HistogramRegistry for QA output
  template <o2::aod::femtoworldparticle::ParticleType part, o2::aod::femtoworldparticle::TrackType tracktype, typename cutContainerType>
  void init(HistogramRegistry* registry);

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

  /// Checks whether the most open combination of all selection criteria is fulfilled
  /// \tparam T Data type of the track
  /// \param track Track
  /// \return Whether the most open combination of all selection criteria is fulfilled
  template <typename T>
  bool isSelectedMinimal(T const& track);

  /// Obtain the bit-wise container for the selections
  /// \todo For the moment, PID is separated from the other selections, hence instead of a single value an std::array of size two is returned
  /// \tparam cutContainerType Data type of the bit-wise container for the selections
  /// \tparam T Data type of the track
  /// \param track Track
  /// \return The bit-wise container for the selections, separately with all selection criteria, and the PID
  template <typename cutContainerType, typename T>
  std::array<cutContainerType, 2> getCutContainer(T const& track);

  /// Some basic QA histograms
  /// \tparam part Type of the particle for proper naming of the folders for QA
  /// \tparam tracktype Type of track (track, positive child, negative child) for proper naming of the folders for QA
  /// \tparam T Data type of the track
  /// \param track Track
  template <o2::aod::femtoworldparticle::ParticleType part, o2::aod::femtoworldparticle::TrackType tracktype, typename T>
  void fillQA(T const& track);

  /// Helper function to obtain the name of a given selection criterion for consistent naming of the configurables
  /// \param iSel Track selection variable to be examined
  /// \param prefix Additional prefix for the name of the configurable
  /// \param suffix Additional suffix for the name of the configurable
  static std::string getSelectionName(femtoWorldTrackSelection::TrackSel iSel, std::string_view prefix = "", std::string_view suffix = "")
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
    LOGF(info, "Variable %s not found", obs);
    return -1;
  }

  /// Helper function to obtain the type of a given selection variable for consistent naming of the configurables
  /// \param iSel Track selection variable whose type is returned
  static femtoWorldSelection::SelectionType getSelectionType(femtoWorldTrackSelection::TrackSel iSel)
  {
    return mSelectionTypes[iSel];
  }

  /// Helper function to obtain the helper string of a given selection criterion for consistent description of the configurables
  /// \param iSel Track selection variable to be examined
  /// \param prefix Additional prefix for the output of the configurable
  static std::string getSelectionHelper(femtoWorldTrackSelection::TrackSel iSel, std::string_view prefix = "")
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

  static constexpr femtoWorldSelection::SelectionType mSelectionTypes[kNtrackSelection]{femtoWorldSelection::kEqual,
                                                                                        femtoWorldSelection::kLowerLimit,
                                                                                        femtoWorldSelection::kUpperLimit,
                                                                                        femtoWorldSelection::kAbsUpperLimit,
                                                                                        femtoWorldSelection::kLowerLimit,
                                                                                        femtoWorldSelection::kLowerLimit,
                                                                                        femtoWorldSelection::kLowerLimit,
                                                                                        femtoWorldSelection::kUpperLimit,
                                                                                        femtoWorldSelection::kLowerLimit,
                                                                                        femtoWorldSelection::kLowerLimit,
                                                                                        femtoWorldSelection::kAbsUpperLimit,
                                                                                        femtoWorldSelection::kAbsUpperLimit,
                                                                                        femtoWorldSelection::kAbsUpperLimit,
                                                                                        femtoWorldSelection::kAbsUpperLimit}; ///< Map to match a variable with its type

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
}; // namespace femtoWorld

template <o2::aod::femtoworldparticle::ParticleType part, o2::aod::femtoworldparticle::TrackType tracktype, typename cutContainerType>
void FemtoWorldTrackSelection::init(HistogramRegistry* registry)
{
  if (registry) {
    mHistogramRegistry = registry;
    std::string folderName = static_cast<std::string>(o2::aod::femtoworldparticle::ParticleTypeName[part]) + "/" + static_cast<std::string>(o2::aod::femtoworldparticle::TrackTypeName[tracktype]);

    /// check whether the number of selection exceeds the bitmap size
    unsigned int nSelections = getNSelections() - getNSelections(femtoWorldTrackSelection::kPIDnSigmaMax);
    if (nSelections > 8 * sizeof(cutContainerType)) {
      LOG(fatal) << "FemtoWorldTrackCuts: Number of selections too large for your container - quitting!";
    }

    mHistogramRegistry->add((folderName + "/hPt").c_str(), "; #it{p}_{T} (GeV/#it{c}); Entries", kTH1F, {{240, 0, 6}});
    mHistogramRegistry->add((folderName + "/hEta").c_str(), "; #eta; Entries", kTH1F, {{200, -1.5, 1.5}});
    mHistogramRegistry->add((folderName + "/hPhi").c_str(), "; #phi; Entries", kTH1F, {{200, 0, 2. * M_PI}});
    mHistogramRegistry->add((folderName + "/hTPCfindable").c_str(), "; TPC findable clusters; Entries", kTH1F, {{163, -0.5, 162.5}});
    mHistogramRegistry->add((folderName + "/hTPCfound").c_str(), "; TPC found clusters; Entries", kTH1F, {{163, -0.5, 162.5}});
    mHistogramRegistry->add((folderName + "/hTPCcrossedOverFindalbe").c_str(), "; TPC ratio findable; Entries", kTH1F, {{100, 0.5, 1.5}});
    mHistogramRegistry->add((folderName + "/hTPCcrossedRows").c_str(), "; TPC crossed rows; Entries", kTH1F, {{163, 0, 163}});
    mHistogramRegistry->add((folderName + "/hTPCfindableVsCrossed").c_str(), ";TPC findable clusters ; TPC crossed rows;", kTH2F, {{163, 0, 163}, {163, 0, 163}});
    mHistogramRegistry->add((folderName + "/hTPCshared").c_str(), "; TPC shared clusters; Entries", kTH1F, {{163, -0.5, 162.5}});
    mHistogramRegistry->add((folderName + "/hITSclusters").c_str(), "; ITS clusters; Entries", kTH1F, {{10, -0.5, 9.5}});
    mHistogramRegistry->add((folderName + "/hITSclustersIB").c_str(), "; ITS clusters in IB; Entries", kTH1F, {{10, -0.5, 9.5}});
    mHistogramRegistry->add((folderName + "/hDCAxy").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {{100, 0, 10}, {500, -5, 5}});
    mHistogramRegistry->add((folderName + "/hDCAz").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{z} (cm)", kTH2F, {{100, 0, 10}, {500, -5, 5}});
    mHistogramRegistry->add((folderName + "/hDCA").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA (cm)", kTH2F, {{100, 0, 10}, {301, 0., 1.5}});
    mHistogramRegistry->add((folderName + "/hTPCdEdX").c_str(), "; #it{p} (GeV/#it{c}); TPC Signal", kTH2F, {{100, 0, 10}, {1000, 0, 1000}});
    mHistogramRegistry->add((folderName + "/nSigmaTPC_el").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{e}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    mHistogramRegistry->add((folderName + "/nSigmaTPC_pi").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{#pi}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    mHistogramRegistry->add((folderName + "/nSigmaTPC_K").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{K}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    mHistogramRegistry->add((folderName + "/nSigmaTPC_p").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{p}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    mHistogramRegistry->add((folderName + "/nSigmaTPC_d").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TPC}^{d}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    mHistogramRegistry->add((folderName + "/nSigmaTOF_el").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{e}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    mHistogramRegistry->add((folderName + "/nSigmaTOF_pi").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{#pi}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    mHistogramRegistry->add((folderName + "/nSigmaTOF_K").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{K}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    mHistogramRegistry->add((folderName + "/nSigmaTOF_p").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{p}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    mHistogramRegistry->add((folderName + "/nSigmaTOF_d").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{TOF}^{d}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    mHistogramRegistry->add((folderName + "/nSigmaComb_el").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{e}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    mHistogramRegistry->add((folderName + "/nSigmaComb_pi").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{#pi}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    mHistogramRegistry->add((folderName + "/nSigmaComb_K").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{K}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    mHistogramRegistry->add((folderName + "/nSigmaComb_p").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{p}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
    mHistogramRegistry->add((folderName + "/nSigmaComb_d").c_str(), "; #it{p} (GeV/#it{c}); n#sigma_{comb}^{d}", kTH2F, {{100, 0, 10}, {100, -5, 5}});
  }
  /// set cuts
  nPtMinSel = getNSelections(femtoWorldTrackSelection::kpTMin);
  nPtMaxSel = getNSelections(femtoWorldTrackSelection::kpTMax);
  nEtaSel = getNSelections(femtoWorldTrackSelection::kEtaMax);
  nTPCnMinSel = getNSelections(femtoWorldTrackSelection::kTPCnClsMin);
  nTPCfMinSel = getNSelections(femtoWorldTrackSelection::kTPCfClsMin);
  nTPCcMinSel = getNSelections(femtoWorldTrackSelection::kTPCcRowsMin);
  nTPCsMaxSel = getNSelections(femtoWorldTrackSelection::kTPCsClsMax);
  nITScMinSel = getNSelections(femtoWorldTrackSelection::kITSnClsMin);
  nITScIbMinSel = getNSelections(femtoWorldTrackSelection::kITSnClsIbMin);
  nDCAxyMaxSel = getNSelections(femtoWorldTrackSelection::kDCAxyMax);
  nDCAzMaxSel = getNSelections(femtoWorldTrackSelection::kDCAzMax);
  nDCAMinSel = getNSelections(femtoWorldTrackSelection::kDCAMin);
  nPIDnSigmaSel = getNSelections(femtoWorldTrackSelection::kPIDnSigmaMax);

  pTMin = getMinimalSelection(femtoWorldTrackSelection::kpTMin, femtoWorldSelection::kLowerLimit);
  pTMax = getMinimalSelection(femtoWorldTrackSelection::kpTMax, femtoWorldSelection::kUpperLimit);
  etaMax = getMinimalSelection(femtoWorldTrackSelection::kEtaMax, femtoWorldSelection::kAbsUpperLimit);
  nClsMin = getMinimalSelection(femtoWorldTrackSelection::kTPCnClsMin, femtoWorldSelection::kLowerLimit);
  fClsMin = getMinimalSelection(femtoWorldTrackSelection::kTPCfClsMin, femtoWorldSelection::kLowerLimit);
  cTPCMin = getMinimalSelection(femtoWorldTrackSelection::kTPCcRowsMin, femtoWorldSelection::kLowerLimit);
  sTPCMax = getMinimalSelection(femtoWorldTrackSelection::kTPCsClsMax, femtoWorldSelection::kUpperLimit);
  nITSclsMin = getMinimalSelection(femtoWorldTrackSelection::kITSnClsMin, femtoWorldSelection::kLowerLimit);
  nITSclsIbMin = getMinimalSelection(femtoWorldTrackSelection::kITSnClsIbMin, femtoWorldSelection::kLowerLimit);
  dcaXYMax = getMinimalSelection(femtoWorldTrackSelection::kDCAxyMax, femtoWorldSelection::kAbsUpperLimit);
  dcaZMax = getMinimalSelection(femtoWorldTrackSelection::kDCAzMax, femtoWorldSelection::kAbsUpperLimit);
  dcaMin = getMinimalSelection(femtoWorldTrackSelection::kDCAMin, femtoWorldSelection::kAbsLowerLimit);
  nSigmaPIDMax = getMinimalSelection(femtoWorldTrackSelection::kPIDnSigmaMax, femtoWorldSelection::kAbsUpperLimit);
}

template <typename T>
auto FemtoWorldTrackSelection::getNsigmaTPC(T const& track, o2::track::PID pid)
{
  return o2::aod::pidutils::tpcNSigma(pid, track);
}

template <typename T>
auto FemtoWorldTrackSelection::getNsigmaTOF(T const& track, o2::track::PID pid)
{
  /// skip tracks without TOF signal
  if (!track.hasTOF()) {
    return 999.f;
  }

  return o2::aod::pidutils::tofNSigma(pid, track);
}

template <typename T>
bool FemtoWorldTrackSelection::isSelectedMinimal(T const& track)
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
  const auto dca = track.dcaXY(); // Accordingly to FemtoWorld in AliPhysics  as well as LF analysis,
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
  if (nEtaSel > 0 && std::abs(eta) > etaMax) {
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
  if (nDCAxyMaxSel > 0 && std::abs(dcaXY) > dcaXYMax) {
    return false;
  }
  if (nDCAzMaxSel > 0 && std::abs(dcaZ) > dcaZMax) {
    return false;
  }
  if (nDCAMinSel > 0 && std::abs(dca) < dcaMin) {
    return false;
  }
  if (nRejectNotPropagatedTracks && std::abs(dca) > 1e3) {
    return false;
  }

  if (nPIDnSigmaSel > 0) {
    bool isFulfilled = false;
    for (size_t i = 0; i < pidTPC.size(); ++i) {
      auto pidTPCVal = pidTPC.at(i);
      if (std::abs(pidTPCVal) < nSigmaPIDMax) {
        isFulfilled = true;
      }
    }
    if (!isFulfilled) {
      return isFulfilled;
    }
  }
  return true;
}

template <typename cutContainerType, typename T>
std::array<cutContainerType, 2> FemtoWorldTrackSelection::getCutContainer(T const& track)
{
  cutContainerType output = 0;
  size_t counter = 0;
  cutContainerType outputPID = 0;
  size_t counterPID = 0;
  const auto sign = track.sign();
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
  const auto dca = std::sqrt(pow(dcaXY, 2.) + pow(dcaZ, 2.));

  std::vector<float> pidTPC, pidTOF;
  for (auto it : mPIDspecies) {
    pidTPC.push_back(getNsigmaTPC(track, it));
    pidTOF.push_back(getNsigmaTOF(track, it));
  }

  float observable = 0.;
  for (auto& sel : mSelections) {
    const auto selVariable = sel.getSelectionVariable();
    if (selVariable == femtoWorldTrackSelection::kPIDnSigmaMax) {
      /// PID needs to be handled a bit differently since we may need more than one species
      for (size_t i = 0; i < pidTPC.size(); ++i) {
        auto pidTPCVal = pidTPC.at(i);
        auto pidTOFVal = pidTOF.at(i);
        sel.checkSelectionSetBit(pidTPCVal, outputPID, counterPID);
        auto pidComb = std::sqrt(pidTPCVal * pidTPCVal + pidTOFVal * pidTOFVal);
        sel.checkSelectionSetBit(pidComb, outputPID, counterPID);
      }

    } else {
      /// for the rest it's all the same
      switch (selVariable) {
        case (femtoWorldTrackSelection::kSign):
          observable = sign;
          break;
        case (femtoWorldTrackSelection::kpTMin):
        case (femtoWorldTrackSelection::kpTMax):
          observable = pT;
          break;
        case (femtoWorldTrackSelection::kEtaMax):
          observable = eta;
          break;
        case (femtoWorldTrackSelection::kTPCnClsMin):
          observable = tpcNClsF;
          break;
        case (femtoWorldTrackSelection::kTPCfClsMin):
          observable = tpcRClsC;
          break;
        case (femtoWorldTrackSelection::kTPCcRowsMin):
          observable = tpcNClsC;
          break;
        case (femtoWorldTrackSelection::kTPCsClsMax):
          observable = tpcNClsS;
          break;
        case (femtoWorldTrackSelection::kITSnClsMin):
          observable = itsNCls;
          break;
        case (femtoWorldTrackSelection::kITSnClsIbMin):
          observable = itsNClsIB;
          break;
        case (femtoWorldTrackSelection::kDCAxyMax):
          observable = dcaXY;
          break;
        case (femtoWorldTrackSelection::kDCAzMax):
          observable = dcaZ;
          break;
        case (femtoWorldTrackSelection::kDCAMin):
          observable = dca;
          break;
        case (femtoWorldTrackSelection::kPIDnSigmaMax):
          break;
      }
      sel.checkSelectionSetBit(observable, output, counter);
    }
  }
  return {output, outputPID};
}

template <o2::aod::femtoworldparticle::ParticleType part, o2::aod::femtoworldparticle::TrackType tracktype, typename T>
void FemtoWorldTrackSelection::fillQA(T const& track)
{
  if (mHistogramRegistry) {
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/hPt"), track.pt());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/hEta"), track.eta());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/hPhi"), track.phi());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/hTPCfindable"), track.tpcNClsFindable());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/hTPCfound"), track.tpcNClsFound());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/hTPCcrossedOverFindalbe"), track.tpcCrossedRowsOverFindableCls());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/hTPCcrossedRows"), track.tpcNClsCrossedRows());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/hTPCfindableVsCrossed"), track.tpcNClsFindable(), track.tpcNClsCrossedRows());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/hTPCshared"), track.tpcNClsShared());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/hITSclusters"), track.itsNCls());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/hITSclustersIB"), track.itsNClsInnerBarrel());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/hDCAxy"), track.pt(), track.dcaXY());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/hDCAz"), track.pt(), track.dcaZ());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/hDCA"), track.pt(), std::sqrt(pow(track.dcaXY(), 2.) + pow(track.dcaZ(), 2.)));
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/hTPCdEdX"), track.p(), track.tpcSignal());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/nSigmaTPC_el"), track.p(), track.tpcNSigmaEl());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/nSigmaTPC_pi"), track.p(), track.tpcNSigmaPi());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/nSigmaTPC_K"), track.p(), track.tpcNSigmaKa());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/nSigmaTPC_p"), track.p(), track.tpcNSigmaPr());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/nSigmaTPC_d"), track.p(), track.tpcNSigmaDe());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/nSigmaTOF_el"), track.p(), track.tofNSigmaEl());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/nSigmaTOF_pi"), track.p(), track.tofNSigmaPi());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/nSigmaTOF_K"), track.p(), track.tofNSigmaKa());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/nSigmaTOF_p"), track.p(), track.tofNSigmaPr());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/nSigmaTOF_d"), track.p(), track.tofNSigmaDe());
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/nSigmaComb_el"), track.p(), std::sqrt(track.tpcNSigmaEl() * track.tpcNSigmaEl() + track.tofNSigmaEl() * track.tofNSigmaEl()));
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/nSigmaComb_pi"), track.p(), std::sqrt(track.tpcNSigmaPi() * track.tpcNSigmaPi() + track.tofNSigmaPi() * track.tofNSigmaPi()));
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/nSigmaComb_K"), track.p(), std::sqrt(track.tpcNSigmaKa() * track.tpcNSigmaKa() + track.tofNSigmaKa() * track.tofNSigmaKa()));
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/nSigmaComb_p"), track.p(), std::sqrt(track.tpcNSigmaPr() * track.tpcNSigmaPr() + track.tofNSigmaPr() * track.tofNSigmaPr()));
    mHistogramRegistry->fill(HIST(o2::aod::femtoworldparticle::ParticleTypeName[part]) + HIST("/") + HIST(o2::aod::femtoworldparticle::TrackTypeName[tracktype]) + HIST("/nSigmaComb_d"), track.p(), std::sqrt(track.tpcNSigmaDe() * track.tpcNSigmaDe() + track.tofNSigmaDe() * track.tofNSigmaDe()));
  }
}

} // namespace o2::analysis::femtoWorld

#endif /* FEMTOWORLDTRACKSELECTION_H_ */
