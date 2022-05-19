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

/// \brief perform photon conversion analysis on V0 candidates from aod::StoredV0Datas
/// dependencies: o2-analysis-lf-lambdakzerobuilder
/// \author stephan.friedrich.stiefelmaier@cern.ch

#include "gammaTables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

#include <TVector3.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// using collisionEvSelIt = soa::Join<aod::Collisions, aod::EvSels>::iterator;
struct GammaConversions {

  // event cuts
  /*
  Configurable<bool> fDoEventSel{"fDoEventSel", 0, "demand sel7 for events"};
  Configurable<float> fCentMin{"fCentMin", 0.0, "lower bound of centrality selection"};
  Configurable<float> fCentMax{"fCentMax", 100.0, "upper bound of centrality selection"};
  */

  // track cuts
  Configurable<float> fTrackEtaMax{"fTrackEtaMax", 0.8, "accepted track eta range"};
  Configurable<float> fTruePhotonEtaMax{"fTruePhotonEtaMax", 0.8, "eta range for true mc photons that serve for validation"};
  Configurable<float> fTrackPtMin{"fTrackPtMin", 0.04, "minimum daughter track pt"};
  Configurable<float> fPIDnSigmaElectronMin{"fPIDnSigmaElectronMin", -3., "minimum sigma electron PID for V0 daughter tracks. Set 0 to disable."};
  Configurable<float> fPIDnSigmaElectronMax{"fPIDnSigmaElectronMax", 3., "maximum sigma electron PID for V0 daughter tracks"};
  Configurable<float> fPIDPionRejectionPMin{"fPIDPionRejectionPMin", 0.4, "minimum track momentum to apply any pion rejection"};                              // case 7:  // 0.4 GeV
  Configurable<float> fPIDPionRejectionPBoarder{"fPIDPionRejectionPBoarder", 8., "border between low and high momentum pion rejection"};                      // case 7:  // 8. GeV
  Configurable<float> fPIDnSigmaAbovePionLineLowPMin{"fPIDnSigmaAbovePionLineLowPMin", 3., "minimum sigma to be over the pion line for low momentum tracks"}; // case 4: 3.0sigma, 1.0 sigma at high momentum
  Configurable<float> fPIDnSigmaAbovePionLineHighPMin{"fPIDnSigmaAbovePionLineHighPMin", 1., "minimum sigma to be over the pion line for high momentum tracks"};
  Configurable<float> fMinTPCFoundOverFindableCls{"fMinTPCNClsFoundOverFindable", 0.6, "minimum ratio found tpc clusters over findable"}; // case 9:  // 0.6
  Configurable<float> fMinTPCCrossedRowsOverFindableCls{"fMinTPCCrossedRowsOverFindableCls", 0.0, "minimum ratio TPC crossed rows over findable clusters"};

  // V0 cuts
  Configurable<float> fV0CosPAngleMin{"fV0CosPAngleMin", 0.85, "Set negative to disable. Minimum cosinus of the pointing angle"}; // case 4
  Configurable<float> fV0RMin{"fV0RMin", 5., "minimum conversion radius of the V0s"};
  Configurable<float> fV0RMax{"fV0RMax", 180., "maximum conversion radius of the V0s"};
  Configurable<float> fV0PhotonAsymmetryMax{"fV0PhotonAsymmetryMax", 0.95, "maximum photon asymetry. Set negative do disable cut."};
  Configurable<float> fV0PsiPairMax{"fV0PsiPairMax", 0.1, "maximum psi angle of the track pair. Set negative do disable cut. "};
  Configurable<float> fV0QtPtMultiplicator{"fV0QtPtMultiplicator", 0.11, "Multiply pt of V0s by this value to get the 2nd denominator in the armenteros cut. The products maximum value is fV0QtMax."};
  Configurable<float> fV0QtMax{"fV0QtMax", 0.040, "the maximum value of the product, that is the maximum qt"};

  // define this in order to have a constructor of the HistogramSpec which copies the name into the title
  struct MyHistogramSpec {
    MyHistogramSpec(char const* const name_, char const* const title_, HistogramConfigSpec config_, bool callSumw2_ = false, bool dataOnly_ = false)
      : m{name_, title_, config_, callSumw2_}, fDataOnly(dataOnly_) {}
    MyHistogramSpec(char const* const name_, HistogramConfigSpec config_, bool callSumw2_ = false, bool dataOnly_ = false)
      : m{name_, name_, config_, callSumw2_}, fDataOnly(dataOnly_) {}
    HistogramSpec m{};
    bool fDataOnly{false};
  };

  // collision histograms
  std::vector<MyHistogramSpec> fCollisionHistoDefinitions{
    {"hCollisionZ", "hCollisionZ;z (cm);counts", {HistType::kTH1F, {{800, -50.f, 50.f}}}}};

  // track histograms
  std::vector<MyHistogramSpec> fTrackHistoDefinitions{
    {"hTrackPt", "hTrackPt;p_{T} (GeV/c);counts", {HistType::kTH1F, {{800, 0.0f, 25.0f}}}},
    {"hTrackEta", "hTrackEta;#eta;counts", {HistType::kTH1F, {{800, -2.f, 2.f}}}},
    {"hTrackPhi", "hTrackPhi;#phi (rad);counts", {HistType::kTH1F, {{800, 0.f, 2.f * M_PI}}}},
    {"hTPCFoundOverFindableCls", "hTPCFoundOverFindableCls;TPCFoundOverFindableCls;counts", {HistType::kTH1F, {{800, 0.9f, 1.01f}}}},
    {"hTPCCrossedRowsOverFindableCls", "hTPCCrossedRowsOverFindableCls;TPCCrossedRowsOverFindableCls;counts", {HistType::kTH1F, {{800, 0.8f, 1.5f}}}},
    {"hTPCdEdx", "hTPCdEdx;p (GeV/c);TPCdEdx", {HistType::kTH2F, {{800, 0.03f, 20.f}, {800, 0.f, 200.f}}}},
    {"hTPCdEdxSigEl", "hTPCdEdxSigEl;p (GeV/c);TPCdEdxSigEl", {HistType::kTH2F, {{800, 0.03f, 20.f}, {800, -10.f, 10.f}}}},
    {"hTPCdEdxSigPi", "hTPCdEdxSigPi;p (GeV/c);TPCdEdxSigPi", {HistType::kTH2F, {{800, 0.03f, 20.f}, {800, -10.f, 10.f}}}}};

  // v0 histograms
  std::vector<MyHistogramSpec> fV0HistoDefinitions{
    {"hPt", "hPt;p_{T} (GeV/c);counts", {HistType::kTH1F, {{800, 0.0f, 25.0f}}}},
    {"hEta", "hEta;#eta;counts", {HistType::kTH1F, {{800, -2.f, 2.f}}}},
    {"hPhi", "hPhi;#phi (rad);counts", {HistType::kTH1F, {{800, 0.f, 2.f * M_PI}}}},
    {"hConvPointR", "hConvPointR;conversion radius (cm);counts", {HistType::kTH1F, {{800, 0.f, 200.f}}}},
    {"hArmenteros", "hArmenteros;#alpha;q_{T} (GeV/c)", {HistType::kTH2F, {{800, -1.f, 1.f}, {800, 0.f, 0.25f}}}},
    {"hPsiPt", "hPsiPt;#Psi;p_{T} (GeV/c)", {HistType::kTH2F, {{800, -2.f, 2.f}, {800, 0.f, 10.f}}}},
    {"hCosPAngle", "hCosPAngle;CosPAngle;counts", {HistType::kTH1F, {{800, 0.99f, 1.005f}}}},
  };

  // only in mc
  // resolution histos
  std::vector<MyHistogramSpec> fV0ResolutionHistoDefinitions{
    {"hPtRes", "hPtRes_Rec-MC;#Delta p_T (GeV/c);counts", {HistType::kTH1F, {{800, -5.f, 5.f}}}},
    {"hEtaRes", "hEtaRes_Rec-MC;#Delta #eta;counts", {HistType::kTH1F, {{800, -3.145f, 3.145f}}}},
    {"hPhiRes", "hPhiRes_Rec-MC;#Delta #phi (rad);counts", {HistType::kTH1F, {{800, -3.145f, 3.145f}}}},
    {"hConvPointRRes", "hConvPointRRes_Rec-MC;#Delta R (cm);counts", {HistType::kTH1F, {{800, -200.f, 200.f}}}},
    {"hConvPointAbsoluteDistanceRes", "hConvPointAbsoluteDistanceRes;euclidean distance (cm);counts", {HistType::kTH1F, {{800, -0.0f, 200.f}}}},
  };

  // think of better name
  std::vector<MyHistogramSpec> fSpecialHistoDefinitions{
    {"hIsPhotonSelected", "hIsPhotonSelected;cut categories;counts", {HistType::kTH1F, {{13, -0.0f, 12.5f}}}},
    {"hTruePhotonSelection", "hTruePhotonSelection;cut categories;counts", {HistType::kTH1F, {{13, -0.0f, 12.5f}}}, true /*dataOnly_*/}};

  std::map<std::string, size_t> fPhotonCutIndeces{
    {"kV0In", 0},
    {"kTrackEta", 1},
    {"kTrackPt", 2},
    {"kElectronPID", 3},
    {"kPionRejLowMom", 4},
    {"kPionRejHighMom", 5},
    {"kTPCFoundOverFindableCls", 6},
    {"kTPCCrossedRowsOverFindableCls", 7},
    {"kV0Radius", 8},
    {"kArmenteros", 9},
    {"kPsiPair", 10},
    {"kCosinePA", 11},
    {"kV0Out", 12}};

  std::map<std::string, size_t> fMcPhotonCutIndeces{
    {"kV0in", 0},
    {"kFakeV0", 1},
    {"kMcMotherIn", 2},
    {"kNoPhysicalPrimary", 3},
    {"kNoGammaMother", 4},
    {"kOutsideMCEtaAcc", 5},
    {"kGoodMcPhotonOut", 12}};

  std::vector<std::string> fHistoSuffixes{"_MCRec", "_MCTrue", "_MCVal", "_Res"};
  enum eBeforeAfterRecCuts { kBeforeRecCuts,
                             kAfterRecCuts };
  enum eV0HistoFlavor { kRec,
                        kMCTrue,
                        kMCVal,
                        kRes };

  typedef std::map<std::string, HistPtr> mapStringHistPtr;
  struct tV0Kind {
    tV0Kind(std::string thePath) : mPath{thePath} {}

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

    void addHistosToOfficalRegistry(HistogramRegistry& theOfficialRegistry,
                                    std::vector<MyHistogramSpec> const& theHistoDefinitions,
                                    std::string const* theSuffix = nullptr)
    {
      for (auto& tHisto : theHistoDefinitions) {
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
    tMotherDirV0Kinds(std::string thePath) : mPath{thePath} {}
    std::string mPath{""};
    tV0Kind mV0Kind[4]{mPath + "Rec/", mPath + "MCTrue/", mPath + "MCVal/", mPath + "Res/"};
  };

  struct tMotherDirRejMc {
    tMotherDirRejMc(std::string thePath) : mPath{thePath} {}
    std::string mPath{""};
    tMotherDirV0Kinds mBeforeAfterRecCuts[2]{mPath + "beforeRecCuts/", mPath + "afterRecCuts/"};
  };

  // for collision track v0
  struct tHistoFolderCTV {
    tHistoFolderCTV(std::string thePath) : mPath{thePath} {}

    std::string mPath{""};
    tV0Kind mSpecialHistos{mPath};
    tMotherDirV0Kinds mBeforeAfterRecCuts[2]{mPath + "beforeRecCuts/",
                                             mPath + "afterRecCuts/"};

    tMotherDirRejMc mRejectedByMc[2]{mPath + "rejectedByMc/kNoPhysicalPrimary/",
                                     mPath + "rejectedByMc/kOutsideMCEtaAcc/"};
  };
  enum eMcRejectionReason { kNoPhysicalPrimary,
                            kOutsideMCEtaAcc };

  struct tHistoRegistry {
    tHistoRegistry(std::string thePath) : mPath{thePath} {}

    std::string mPath{""};
    tHistoFolderCTV mCollision{mPath + "Collision/"};
    tHistoFolderCTV mTrack{mPath + "Track/"};
    tHistoFolderCTV mV0{mPath + "V0/"};
  };
  tHistoRegistry fMyRegistry{""};
  HistogramRegistry fHistogramRegistry{"fHistogramRegistry"};

  void init(InitContext const&)
  {
    auto addLablesToHisto = [](auto const& theContainer, std::string const& theHistoName, auto const& theLables) {
      auto lHisto = theContainer.find(theHistoName);
      if (lHisto != theContainer.end()) {
        TAxis* lXaxis = std::get<std::shared_ptr<TH1>>(lHisto->second)->GetXaxis();
        for (auto& lPairIt : theLables) {
          lXaxis->SetBinLabel(lPairIt.second + 1, lPairIt.first.data());
        }
      }
    };

    if (doprocessRec) {
      fHistoSuffixes[0] = "Rec";
      fSpecialHistoDefinitions.pop_back();
    }

    fMyRegistry.mV0.mSpecialHistos.addHistosToOfficalRegistry(
      fHistogramRegistry,
      fSpecialHistoDefinitions);

    // do some labeling
    addLablesToHisto(fMyRegistry.mV0.mSpecialHistos.mContainer, "hIsPhotonSelected", fPhotonCutIndeces);
    if (doprocessMc) {
      addLablesToHisto(fMyRegistry.mV0.mSpecialHistos.mContainer, "hTruePhotonSelection", fMcPhotonCutIndeces);
    }

    for (size_t iBARecCuts = 0; iBARecCuts < 2; ++iBARecCuts) {
      for (size_t iMcKind = 0; iMcKind < (doprocessMc ? 4 : 1); ++iMcKind) {
        std::string const* lMcSuffix = &fHistoSuffixes[iMcKind];

        // for track and collision histos we only plot reconstructed quantities at the moment
        if (iMcKind == kRec) {
          // collision histograms
          fMyRegistry.mCollision.mBeforeAfterRecCuts[iBARecCuts].mV0Kind[iMcKind].addHistosToOfficalRegistry(fHistogramRegistry,
                                                                                                             fCollisionHistoDefinitions,
                                                                                                             lMcSuffix);

          // track histograms
          fMyRegistry.mTrack.mBeforeAfterRecCuts[iBARecCuts].mV0Kind[iMcKind].addHistosToOfficalRegistry(fHistogramRegistry,
                                                                                                         fTrackHistoDefinitions,
                                                                                                         lMcSuffix);
        }

        // v0 histograms
        fMyRegistry.mV0.mBeforeAfterRecCuts[iBARecCuts].mV0Kind[iMcKind].addHistosToOfficalRegistry(fHistogramRegistry,
                                                                                                    (iMcKind < 3) ? fV0HistoDefinitions : fV0ResolutionHistoDefinitions,
                                                                                                    lMcSuffix);

        if (doprocessMc) {
          // v0 mc rejection histos
          for (size_t iRejReason = 0; iRejReason < 2; ++iRejReason) {
            fMyRegistry.mV0.mRejectedByMc[iRejReason].mBeforeAfterRecCuts[iBARecCuts].mV0Kind[iMcKind].addHistosToOfficalRegistry(fHistogramRegistry,
                                                                                                                                  (iMcKind < 3) ? fV0HistoDefinitions : fV0ResolutionHistoDefinitions,
                                                                                                                                  lMcSuffix);
          }
        }
      }
    }
  }

  // SFS todo: think about if this is actually too expensive. Going the other way round with the indices as keys wouldnt require lookups at inserting but pbly produce a but of code duplication at the definition of the cut names
  size_t gMax_size = (size_t)-1;
  size_t getPhotonCutIndex(std::string const& theKey, bool theMcCuts = false)
  {
    auto const& lMap = theMcCuts ? fMcPhotonCutIndeces : fPhotonCutIndeces;
    auto lPairIt = lMap.find(theKey);
    if (lPairIt != lMap.end()) {
      return lPairIt->second;
    }
    return gMax_size;
  }

  void fillIsPhotonSelected(std::string theCase)
  {
    fillTH1(fMyRegistry.mV0.mSpecialHistos.mContainer, "hIsPhotonSelected", getPhotonCutIndex(theCase));
  }

  void fillTruePhotonSelection(std::string theCase)
  {
    fillTH1(fMyRegistry.mV0.mSpecialHistos.mContainer, "hTruePhotonSelection", getPhotonCutIndex(theCase, true /*theMcCuts*/));
  }

  template <typename TV0, typename TTRACKS>
  bool processV0(TV0 const& theV0, const float& theV0CosinePA, TTRACKS const& theTwoV0Daughters)
  {
    fillReconstructedInfoHistograms(
      kBeforeRecCuts,
      theV0,
      theTwoV0Daughters,
      theV0CosinePA);

    // apply track cuts and photon cuts
    if (!(v0PassesTrackCuts(theTwoV0Daughters) &&
          v0PassesPhotonCuts(theV0, theV0CosinePA))) {
      return false;
    }

    fillReconstructedInfoHistograms(
      kAfterRecCuts,
      theV0,
      theTwoV0Daughters,
      theV0CosinePA);

    return kTRUE;
  }

  template <typename TV0, typename TMCGAMMA>
  void fillAllV0HistogramsForRejectedByMc(int theRejReason,
                                          int theBefAftRec,
                                          TMCGAMMA const& theMcPhoton,
                                          TV0 const& theV0,
                                          float const& theV0CosinePA)
  {
    fillV0Histograms(
      fMyRegistry.mV0.mRejectedByMc[theRejReason].mBeforeAfterRecCuts[theBefAftRec].mV0Kind[kRec].mContainer,
      theV0,
      theV0CosinePA);

    fillTruePhotonHistogramsForRejectedByMc(theRejReason,
                                            theBefAftRec,
                                            theMcPhoton,
                                            theV0,
                                            theV0CosinePA);
  }

  template <typename TV0, typename TMCGAMMA>
  bool v0IsGoodValidatedMcPhoton(TMCGAMMA const& theMcPhoton, TV0 const& theV0, float const& theV0CosinePA, bool theV0PassesRecCuts)
  {
    auto fillRejectionHistos = [&](int theRejReason) {
      fillAllV0HistogramsForRejectedByMc(theRejReason,
                                         kBeforeRecCuts,
                                         theMcPhoton,
                                         theV0,
                                         theV0CosinePA);
      if (theV0PassesRecCuts) {
        fillAllV0HistogramsForRejectedByMc(theRejReason,
                                           kAfterRecCuts,
                                           theMcPhoton,
                                           theV0,
                                           theV0CosinePA);
      }
    };

    fillTruePhotonSelection("kMcMotherIn");
    if (!theMcPhoton.isPhysicalPrimary()) {
      fillTruePhotonSelection("kNoPhysicalPrimary");
      fillRejectionHistos(kNoPhysicalPrimary);
      return false;
    }

    if (std::abs(theMcPhoton.eta()) > fTruePhotonEtaMax) {
      fillTruePhotonSelection("kOutsideMCEtaAcc");
      fillRejectionHistos(kOutsideMCEtaAcc);
      return false;
    }
    fillTruePhotonSelection("kGoodMcPhotonOut");
    return true;
  }

  template <typename TV0, typename TMCGAMMATABLE, typename TTRACKS>
  void processMcPhoton(TMCGAMMATABLE const& theMcPhotonForThisV0AsTable,
                       TV0 const& theV0,
                       float const& theV0CosinePA,
                       bool theV0PassesRecCuts,
                       TTRACKS const& theTwoV0Daughters)
  {
    fillTruePhotonSelection("kV0in");

    // is a table that might be empty
    if (theMcPhotonForThisV0AsTable.begin() == theMcPhotonForThisV0AsTable.end()) {
      fillTruePhotonSelection("kFakeV0");
      return;
    }
    auto const lMcPhoton = theMcPhotonForThisV0AsTable.begin();

    if (!v0IsGoodValidatedMcPhoton(lMcPhoton,
                                   theV0,
                                   theV0CosinePA,
                                   theV0PassesRecCuts)) {
      return;
    }

    fillTruePhotonHistograms(kBeforeRecCuts,
                             lMcPhoton,
                             theV0,
                             theV0CosinePA);

    if (theV0PassesRecCuts) {
      fillTruePhotonHistograms(kAfterRecCuts,
                               lMcPhoton,
                               theV0,
                               theV0CosinePA);
    }
  }

  template <typename TV0, typename TMCGAMMA>
  void fillV0ResolutionHistograms(mapStringHistPtr& theContainer,
                                  TMCGAMMA const& theMcPhoton,
                                  TV0 const& theV0)
  {
    TVector3 lConvPointRec(theV0.x(), theV0.y(), theV0.z());
    TVector3 lConvPointTrue(theMcPhoton.conversionX(), theMcPhoton.conversionY(), theMcPhoton.conversionZ());

    fillTH1(theContainer, "hPtRes", theV0.pt() - theMcPhoton.pt());
    fillTH1(theContainer, "hEtaRes", theV0.eta() - theMcPhoton.eta());
    fillTH1(theContainer, "hPhiRes", theV0.phi() - theMcPhoton.phi());
    fillTH1(theContainer, "hConvPointRRes", theV0.v0radius() - lConvPointTrue.Perp());
    fillTH1(theContainer, "hConvPointAbsoluteDistanceRes", TVector3(lConvPointRec - lConvPointTrue).Mag());
  }

  template <typename TV0, typename TMCGAMMA> // use enumtypes?
  void fillTruePhotonHistograms(int theBefAftRec, TMCGAMMA const& theMcPhoton, TV0 const& theV0, float const& theV0CosinePA)
  {
    fillV0HistogramsMcGamma(
      fMyRegistry.mV0.mBeforeAfterRecCuts[theBefAftRec].mV0Kind[kMCTrue].mContainer,
      theMcPhoton);

    fillV0Histograms(
      fMyRegistry.mV0.mBeforeAfterRecCuts[theBefAftRec].mV0Kind[kMCVal].mContainer,
      theV0,
      theV0CosinePA);

    fillV0ResolutionHistograms(
      fMyRegistry.mV0.mBeforeAfterRecCuts[theBefAftRec].mV0Kind[kRes].mContainer,
      theMcPhoton,
      theV0);
  }

  template <typename TV0, typename TMCGAMMA> // use enumtypes?
  void fillTruePhotonHistogramsForRejectedByMc(int theRejReason,
                                               int theBefAftRec,
                                               TMCGAMMA const& theMcPhoton,
                                               TV0 const& theV0,
                                               float const& theV0CosinePA)
  {
    fillV0HistogramsMcGamma(
      fMyRegistry.mV0.mRejectedByMc[theRejReason].mBeforeAfterRecCuts[theBefAftRec].mV0Kind[kMCTrue].mContainer,
      theMcPhoton);

    fillV0Histograms(
      fMyRegistry.mV0.mRejectedByMc[theRejReason].mBeforeAfterRecCuts[theBefAftRec].mV0Kind[kMCVal].mContainer,
      theV0,
      theV0CosinePA);

    fillV0ResolutionHistograms(
      fMyRegistry.mV0.mRejectedByMc[theRejReason].mBeforeAfterRecCuts[theBefAftRec].mV0Kind[kRes].mContainer,
      theMcPhoton,
      theV0);
  }

  void processRec(aod::Collisions::iterator const& theCollision,
                  aod::V0Datas const& theV0s,
                  aod::V0DaughterTracks const& theAllTracks)
  {
    fillTH1(fMyRegistry.mCollision.mBeforeAfterRecCuts[kBeforeRecCuts].mV0Kind[kRec].mContainer,
            "hCollisionZ",
            theCollision.posZ());

    for (auto& lV0 : theV0s) {

      auto lTwoV0Daughters = theAllTracks.sliceBy(aod::v0data::v0Id, lV0.v0Id());
      float lV0CosinePA = lV0.v0cosPA(theCollision.posX(), theCollision.posY(), theCollision.posZ());

      if (!processV0(lV0, lV0CosinePA, lTwoV0Daughters)) {
        continue;
      }
    }
  }
  PROCESS_SWITCH(GammaConversions, processRec, "process reconstructed info", true);

  void processMc(aod::Collisions::iterator const& theCollision,
                 aod::V0Datas const& theV0s,
                 aod::V0DaughterTracks const& theAllTracks,
                 aod::McGammasTrue const& theV0sTrue)
  {
    fillTH1(fMyRegistry.mCollision.mBeforeAfterRecCuts[kBeforeRecCuts].mV0Kind[kRec].mContainer,
            "hCollisionZ",
            theCollision.posZ());

    for (auto& lV0 : theV0s) {

      auto lTwoV0Daughters = theAllTracks.sliceBy(aod::v0data::v0Id, lV0.v0Id());
      float lV0CosinePA = lV0.v0cosPA(theCollision.posX(), theCollision.posY(), theCollision.posZ());

      // check if V0 passes rec cuts and fill beforeRecCuts,afterRecCuts [kRec]
      bool lV0PassesRecCuts = processV0(lV0, lV0CosinePA, lTwoV0Daughters);

      // check if it comes from a true photon (lMcPhotonForThisV0AsTable is a table that might be empty)
      auto lMcPhotonForThisV0AsTable = theV0sTrue.sliceBy(aod::v0data::v0Id, lV0.v0Id());
      processMcPhoton(lMcPhotonForThisV0AsTable,
                      lV0,
                      lV0CosinePA,
                      lV0PassesRecCuts,
                      lTwoV0Daughters);
    }
  }
  PROCESS_SWITCH(GammaConversions, processMc, "process reconstructed info and mc", false);

  template <typename T>
  std::shared_ptr<T> getTH(mapStringHistPtr const& theMap, std::string const& theName)
  {
    auto lPairIt = theMap.find(theName);
    if (lPairIt == theMap.end()) {
      LOGF(warning, "SFS getTH(): No element with key %s in map.", theName.data());
      return std::shared_ptr<T>(); // SFS todo verify this is correct
    }
    if (!std::holds_alternative<std::shared_ptr<T>>(lPairIt->second)) {
      LOGF(warning, "SFS getTH(): No shared_ptr<T> in map for key %s.", theName.data());
      return std::shared_ptr<T>();
    }
    std::shared_ptr<T> lHisto = std::get<std::shared_ptr<T>>(lPairIt->second);
    if (lHisto == nullptr) {
      LOGF(warning, "SFS getTH(): Found shared_ptr<T> for key %s but it is a nullptr.", theName.data());
      return std::shared_ptr<T>();
    }
    return lHisto;
  }

  void fillTH1(mapStringHistPtr const& theMap, std::string const& theName, float theValue)
  {
    std::shared_ptr<TH1> lHisto = getTH<TH1>(theMap, theName);
    if (lHisto != nullptr) {
      lHisto->Fill(theValue);
    }
  }

  void fillTH2(mapStringHistPtr const& theMap, std::string const& theName, float theValueX, float theValueY)
  {
    std::shared_ptr<TH2> lHisto = getTH<TH2>(theMap, theName);
    if (lHisto != nullptr) {
      lHisto->Fill(theValueX, theValueY);
    }
  }

  template <typename TTRACKS>
  void fillTrackHistograms(mapStringHistPtr& theContainer, TTRACKS const& theTwoV0Daughters)
  {
    auto fillTrackHistogramsI = [&](auto const& theTrack) {
      fillTH1(theContainer, "hTrackEta", theTrack.eta());
      fillTH1(theContainer, "hTrackPhi", theTrack.phi());
      fillTH1(theContainer, "hTrackPt", theTrack.pt());
      fillTH1(theContainer, "hTPCFoundOverFindableCls", theTrack.tpcFoundOverFindableCls());
      fillTH1(theContainer, "hTPCCrossedRowsOverFindableCls", theTrack.tpcCrossedRowsOverFindableCls());
      fillTH2(theContainer, "hTPCdEdxSigEl", theTrack.p(), theTrack.tpcNSigmaEl());
      fillTH2(theContainer, "hTPCdEdxSigPi", theTrack.p(), theTrack.tpcNSigmaPi());
      fillTH2(theContainer, "hTPCdEdx", theTrack.p(), theTrack.tpcSignal());
    };

    for (auto& lTrack : theTwoV0Daughters) {
      fillTrackHistogramsI(lTrack);
    }
  }

  template <typename TV0>
  void fillV0Histograms(mapStringHistPtr& theContainer, TV0 const& theV0, float const& theV0CosinePA)
  {
    fillTH1(theContainer, "hEta", theV0.eta());
    fillTH1(theContainer, "hPhi", theV0.phi());
    fillTH1(theContainer, "hPt", theV0.pt());
    fillTH1(theContainer, "hConvPointR", theV0.v0radius());
    fillTH1(theContainer, "hCosPAngle", theV0CosinePA);
    fillTH2(theContainer, "hArmenteros", theV0.alpha(), theV0.qtarm());
    fillTH2(theContainer, "hPsiPt", theV0.psipair(), theV0.pt());
  }

  // SFS todo: combine fillV0Histograms and fillV0HistogramsMcGamma
  template <typename TMCGAMMA>
  void fillV0HistogramsMcGamma(mapStringHistPtr& theContainer, TMCGAMMA const& theMcGamma)
  {
    fillTH1(theContainer, "hEta", theMcGamma.eta());
    fillTH1(theContainer, "hPhi", theMcGamma.phi());
    fillTH1(theContainer, "hPt", theMcGamma.pt());
    fillTH1(theContainer, "hConvPointR", theMcGamma.v0Radius());
  }

  template <typename T>
  bool trackPassesCuts(const T& theTrack)
  {
    // single track eta cut
    if (TMath::Abs(theTrack.eta()) > fTrackEtaMax) {
      fillIsPhotonSelected("kTrackEta");
      return kFALSE;
    }

    // single track pt cut
    if (theTrack.pt() < fTrackPtMin) {
      fillIsPhotonSelected("kTrackPt");
      return kFALSE;
    }

    if (!(selectionPIDTPC_track(theTrack))) {
      return kFALSE;
    }

    if (theTrack.tpcFoundOverFindableCls() < fMinTPCFoundOverFindableCls) {
      fillIsPhotonSelected("kTPCFoundOverFindableCls");
      return kFALSE;
    }

    if (theTrack.tpcCrossedRowsOverFindableCls() < fMinTPCCrossedRowsOverFindableCls) {
      fillIsPhotonSelected("kTPCCrossedRowsOverFindableCls");
      return kFALSE;
    }
    return kTRUE;
  }

  template <typename TTRACKS>
  bool v0PassesTrackCuts(TTRACKS const& theTwoV0Daughters)
  {
    for (auto& lTrack : theTwoV0Daughters) {
      if (!trackPassesCuts(lTrack)) {
        return kFALSE;
      }
    }
    return true;
  }

  template <typename TV0>
  bool v0PassesPhotonCuts(const TV0& theV0, float theV0CosinePA)
  {
    if (theV0.v0radius() < fV0RMin || theV0.v0radius() > fV0RMax) {
      fillIsPhotonSelected("kV0Radius");
      return kFALSE;
    }

    if (fV0PhotonAsymmetryMax > 0. && !ArmenterosQtCut(theV0.alpha(), theV0.qtarm(), theV0.pt())) {
      fillIsPhotonSelected("kArmenteros");
      return kFALSE;
    }

    if (fV0PsiPairMax > 0. && TMath::Abs(theV0.psipair()) > fV0PsiPairMax) {
      fillIsPhotonSelected("kPsiPair");
      return kFALSE;
    }

    if (fV0CosPAngleMin > 0. && theV0CosinePA < fV0CosPAngleMin) {
      fillIsPhotonSelected("kCosinePA");
      return kFALSE;
    }
    return kTRUE;
  }

  template <typename TV0, typename TTRACKS>
  void fillReconstructedInfoHistograms(int theBefAftRec, TV0 const& theV0, TTRACKS const& theTwoV0Daughters, float const& theV0CosinePA)
  {
    fillIsPhotonSelected(!theBefAftRec ? "kV0In" : "kV0Out");

    fillTrackHistograms(
      fMyRegistry.mTrack.mBeforeAfterRecCuts[theBefAftRec].mV0Kind[kRec].mContainer,
      theTwoV0Daughters);

    fillV0Histograms(
      fMyRegistry.mV0.mBeforeAfterRecCuts[theBefAftRec].mV0Kind[kRec].mContainer,
      theV0,
      theV0CosinePA);
  }

  Bool_t ArmenterosQtCut(Double_t theAlpha, Double_t theQt, Double_t thePt)
  {
    // in AliPhysics this is the cut for if fDo2DQt && fDoQtGammaSelection == 2
    Float_t lQtMaxPtDep = fV0QtPtMultiplicator * thePt;
    if (lQtMaxPtDep > fV0QtMax) {
      lQtMaxPtDep = fV0QtMax;
    }
    if ((TMath::Power(theAlpha / fV0PhotonAsymmetryMax, 2) + TMath::Power(theQt / lQtMaxPtDep, 2)) >= 1) {
      return kFALSE;
    }
    return kTRUE;
  }

  template <typename T>
  bool selectionPIDTPC_track(const T& theTrack)
  {
    // TPC Electron Line
    if (fPIDnSigmaElectronMin && (theTrack.tpcNSigmaEl() < fPIDnSigmaElectronMin || theTrack.tpcNSigmaEl() > fPIDnSigmaElectronMax)) {
      fillIsPhotonSelected("kElectronPID");
      return kFALSE;
    }

    // TPC Pion Line
    if (theTrack.p() > fPIDPionRejectionPMin) {
      // low pt Pion rej
      if (theTrack.p() < fPIDPionRejectionPBoarder) {
        if (theTrack.tpcNSigmaEl() > fPIDnSigmaElectronMin && theTrack.tpcNSigmaEl() < fPIDnSigmaElectronMax && theTrack.tpcNSigmaPi() < fPIDnSigmaAbovePionLineLowPMin) {
          fillIsPhotonSelected("kPionRejLowMom");
          return kFALSE;
        }
      }
      // High Pt Pion rej
      else {
        if (theTrack.tpcNSigmaEl() > fPIDnSigmaElectronMin && theTrack.tpcNSigmaEl() < fPIDnSigmaElectronMax && theTrack.tpcNSigmaPi() < fPIDnSigmaAbovePionLineHighPMin) {
          fillIsPhotonSelected("kPionRejHighMom");
          return kFALSE;
        }
      }
    }
    return kTRUE;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<GammaConversions>(cfgc)};
}
