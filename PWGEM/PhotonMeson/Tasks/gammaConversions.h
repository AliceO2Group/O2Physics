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

/// \brief contains type definitions for gammaConversions.cxx
/// \author stephan.friedrich.stiefelmaier@cern.ch

#include "Framework/AnalysisTask.h"

using namespace o2::framework;

typedef std::map<std::string, HistPtr> mapStringHistPtr;

// define this in order to have a constructor of the HistogramSpec which copies the name into the title and to have dataOnly flag
struct MyHistogramSpec {
  MyHistogramSpec(char const* const name_, char const* const title_, HistogramConfigSpec config_, bool callSumw2_ = false, bool dataOnly_ = false)
    : m{name_, title_, config_, callSumw2_}, fDataOnly(dataOnly_) {}
  MyHistogramSpec(char const* const name_, HistogramConfigSpec config_, bool callSumw2_ = false, bool dataOnly_ = false)
    : m{name_, name_, config_, callSumw2_}, fDataOnly(dataOnly_) {}
  HistogramSpec m{};
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

enum eBeforeAfterRecCuts { kBeforeRecCuts,
                           kAfterRecCuts };
enum eV0HistoFlavor { kRec,
                      kMCTrue,
                      kMCVal,
                      kRes };

struct tV0Kind {
  tV0Kind(std::string const& thePath) : mPath{thePath} {}

  // todo: remove
  template <typename T>
  void appendSuffixToTitleI(HistPtr& theHistPtr, std::string const* theSuffix)
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
  void addHistosToOfficalRegistry(HistogramRegistry& theOfficialRegistry,
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
      HistPtr lHistPtr = theOfficialRegistry.add(lFullName.data(), tHisto.m.title.data(), tHisto.m.config);
      mContainer.insert(std::pair{tHisto.m.name, lHistPtr});

      // todo ugly: remove
      if (theSuffix) {
        if (tHisto.m.config.type == kTH1F) {
          appendSuffixToTitleI<TH1>(lHistPtr, theSuffix);
        } else if (tHisto.m.config.type == kTH2F) {
          appendSuffixToTitleI<TH2>(lHistPtr, theSuffix);
        }
      }
    }
  }

  std::string mPath{};
  mapStringHistPtr mContainer{};
};

struct tMotherDirV0Kinds {
  tMotherDirV0Kinds(std::string const& thePath) : mPath{thePath} {}
  std::string mPath{""};
  tV0Kind mV0Kind[4]{mPath + "Rec/", mPath + "MCTrue/", mPath + "MCVal/", mPath + "Res/"};
};

struct tMotherDirRejMc {
  tMotherDirRejMc(std::string const& thePath) : mPath{thePath} {}
  std::string mPath{""};
  tMotherDirV0Kinds mBeforeAfterRecCuts[2]{mPath + "beforeRecCuts/", mPath + "afterRecCuts/"};
};

// for collision track v0
struct tHistoFolderCTV {
  tHistoFolderCTV(std::string const& thePath) : mPath{thePath} {}

  std::string mPath{""};
  tV0Kind mSpecialHistos{mPath};
  tMotherDirV0Kinds mBeforeAfterRecCuts[2]{mPath + "beforeRecCuts/",
                                           mPath + "afterRecCuts/"};

  tMotherDirRejMc mRejectedByMc[2]{mPath + "rejectedByMc/kNoPhysicalPrimary/",
                                   mPath + "rejectedByMc/kOutsideMCEtaAcc/"};
};

enum class eMcRejectedSaved { kNoPhysicalPrimary,
                              kOutsideMCEtaAcc };

struct tHistoRegistry {
  tHistoRegistry(std::string const& thePath) : mPath{thePath} {}

  std::string mPath{""};
  tHistoFolderCTV mCollision{mPath + "Collision/"};
  tHistoFolderCTV mTrack{mPath + "Track/"};
  tHistoFolderCTV mV0{mPath + "V0/"};
};
