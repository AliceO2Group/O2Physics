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

/// \file gammaConversions.h
/// \brief contains type definitions for gammaConversions.cxx
/// \author stephan.friedrich.stiefelmaier@cern.ch

#ifndef PWGEM_PHOTONMESON_TASKS_GAMMACONVERSIONS_H_
#define PWGEM_PHOTONMESON_TASKS_GAMMACONVERSIONS_H_

#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/Logger.h>

#include <TH1.h>
#include <TH2.h>

#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

typedef std::map<std::string, o2::framework::HistPtr> mapStringHistPtr;

// define this in order to have a constructor of the HistogramSpec which copies the name into the title and to have dataOnly flag
struct MyHistogramSpec {
  MyHistogramSpec(char const* const name_, char const* const title_, o2::framework::HistogramConfigSpec config_, bool callSumw2_ = false, bool dataOnly_ = false)
    : m{name_, title_, config_, callSumw2_}, fDataOnly(dataOnly_) {}
  MyHistogramSpec(char const* const name_, o2::framework::HistogramConfigSpec config_, bool callSumw2_ = false, bool dataOnly_ = false)
    : m{name_, name_, config_, callSumw2_}, fDataOnly(dataOnly_) {}
  o2::framework::HistogramSpec m{};
  bool fDataOnly{false};
};

enum class ePhotonCuts {
  kV0In,
  kTrackEta,
  kTrackPt,
  kElectronPID,
  kPionRejLowMom,
  kPionRejHighMom,
  kTPCFoundOverFindableCls,
  kTPCCrossedRowsOverFindableCls,
  kV0Radius,
  kArmenteros,
  kPsiPair,
  kCosinePA,
  kRZLine,
  kV0Out
};

/* These are the distinguished cases why a reconstructed V0, may not become a MC validated photon.
 * This is done before reconstruction cuts. (except kMcValAfterRecCuts, see comment.) */
enum class eV0McValidation {
  kV0in,
  kFakeV0,
  kMcMotherIn,
  kNoPhysicalPrimary,
  kNoPhoton,
  kOutsideMCEtaAcc,
  kMcValidatedPhotonOut, // meaning the V0 comes from a true, primary, mc photon that has a true |eta| < fTruePhotonEtaMax.
  kMcValAfterRecCuts     // the kMcValidatedPhotonOut which also pass reconstruction cuts.
};

enum class eV0Decays {
  ee1,         // Electron - Positron with same mother (true V0)
  ee2,         // Electron - Positron with different mother
  epi,         // Electron/Positron - Pion
  ek,          // Electron/Positron - Kaon
  ep,          // Electron/Positron - Proton/Antiproton
  emu,         // Electron/Positron - Muon
  pipi,        // Pion - Pion
  pik,         // Pion - Kaon
  pip,         // Pion - Proton/Antiproton
  pimu,        // Pion - Muon
  pKmu,        // Proton/Antiproton - Kaon/Muon
  other,       // other
  nomcparticle // no mc particle was found for this track
};

enum eBeforeAfterRecCuts { kBeforeRecCuts,
                           kAfterRecCuts };
enum eV0HistoFlavor { kRec,
                      kMCTrue,
                      kMCVal,
                      kRes };

struct tV0Kind {
  explicit tV0Kind(std::string const& thePath) : mPath{thePath} {}

  // todo: remove
  template <typename T>
  void appendSuffixToTitleI(o2::framework::HistPtr& theHistPtr, std::string const* theSuffix)
  {
    auto lHisto = std::get<std::shared_ptr<T>>(theHistPtr);
    if (lHisto) {
      std::string lTitle(lHisto->GetTitle());
      lHisto->SetTitle(lTitle.append(*theSuffix).data());
    } else {
      LOGF(info, "SFS appendSuffixToTitle(): %s could not be obtained in order to append suffix to title.");
    }
  }

  // todo: bind tV0Kind to gammaConverions such that theOfficialRegistry and theCheckDataOnly are already known
  void addHistosToOfficalRegistry(o2::framework::HistogramRegistry& theOfficialRegistry,
                                  std::vector<MyHistogramSpec> const& theHistoDefinitions,
                                  std::string const* theSuffix = nullptr,
                                  bool theCheckDataOnly = false)
  {
    for (auto& tHisto : theHistoDefinitions) {
      if (theCheckDataOnly && tHisto.fDataOnly) {
        continue;
      }
      std::string lFullName(mPath + tHisto.m.name + (theSuffix ? *theSuffix : std::string("")));
      LOGF(info, "adding %s %d", lFullName, tHisto.fDataOnly);
      o2::framework::HistPtr lHistPtr = theOfficialRegistry.add(lFullName.data(), tHisto.m.title.data(), tHisto.m.config);
      mContainer.insert(std::pair{tHisto.m.name, lHistPtr});

      // todo ugly: remove
      if (theSuffix) {
        if (tHisto.m.config.type == o2::framework::kTH1F) {
          appendSuffixToTitleI<TH1>(lHistPtr, theSuffix);
        } else if (tHisto.m.config.type == o2::framework::kTH2F) {
          appendSuffixToTitleI<TH2>(lHistPtr, theSuffix);
        }
      }
    }
  }

  std::string mPath{};
  mapStringHistPtr mContainer{};
};

struct tMotherDirV0Kinds {
  explicit tMotherDirV0Kinds(std::string const& thePath) : mPath{thePath} {}
  std::string mPath{""};
  tV0Kind mV0Kind[4]{tV0Kind{mPath + "Rec/"}, tV0Kind{mPath + "MCTrue/"}, tV0Kind{mPath + "MCVal/"}, tV0Kind{mPath + "Res/"}};
};

struct tMotherDirRejMc {
  explicit tMotherDirRejMc(std::string const& thePath) : mPath{thePath} {}
  std::string mPath{""};
  tMotherDirV0Kinds mBeforeAfterRecCuts[2]{tMotherDirV0Kinds{mPath + "beforeRecCuts/"}, tMotherDirV0Kinds{mPath + "afterRecCuts/"}};
};

// for collision track v0
struct tHistoFolderCTV {
  explicit tHistoFolderCTV(std::string const& thePath) : mPath{thePath} {}

  std::string mPath{""};
  tV0Kind mSpecialHistos{mPath};
  tMotherDirV0Kinds mBeforeAfterRecCuts[2]{tMotherDirV0Kinds{mPath + "beforeRecCuts/"},
                                           tMotherDirV0Kinds{mPath + "afterRecCuts/"}};

  tMotherDirRejMc mRejectedByMc[2]{tMotherDirRejMc{mPath + "rejectedByMc/kNoPhysicalPrimary/"},
                                   tMotherDirRejMc{mPath + "rejectedByMc/kOutsideMCEtaAcc/"}};
};

enum class eMcRejectedSaved { kNoPhysicalPrimary,
                              kOutsideMCEtaAcc };

struct tHistoRegistry {
  explicit tHistoRegistry(std::string const& thePath) : mPath{thePath} {}

  std::string mPath{""};
  tHistoFolderCTV mCollision{mPath + "Collision/"};
  tHistoFolderCTV mTrack{mPath + "Track/"};
  tHistoFolderCTV mV0{mPath + "V0/"};
};

#endif // PWGEM_PHOTONMESON_TASKS_GAMMACONVERSIONS_H_
