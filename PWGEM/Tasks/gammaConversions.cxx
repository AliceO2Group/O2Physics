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
      : m{name_, title_, config_, callSumw2_}, fDataOnly{dataOnly_} {}
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

  enum eRecTrueEnum { kRec,
                      kMCTrue,
                      kMCVal };

  typedef std::map<std::string, HistPtr> mapStringHistPtr;
  mapStringHistPtr fCollisionHistos;
  mapStringHistPtr fTrackHistos;
  mapStringHistPtr fV0ResolutionHistos;
  mapStringHistPtr fSpecialHistos;
  std::vector<mapStringHistPtr> fRecTrueV0Histos{3};

  std::vector<std::string> fRecTrueStrings{"_MCRec", "_MCTrue", "_MCVal"};
  std::string fPathRec{"Rec/"};
  std::string fPathMc{"MC/"};

  std::string fPathRecCollisionHistos{fPathRec + "collision/"};
  std::string fPathRecTrackHistos{fPathRec + "track/"};

  std::vector<std::string> fPathsV0Histos{
    fPathRec + "v0/",
    fPathMc + "v0_afterMcCuts/MCTrue/",
    fPathMc + "v0_afterMcCuts/MCVal/"};
  std::string fPathResolutions{fPathMc + "v0_afterMcCuts/resolutions/"};

  std::string fFullNameIsPhotonSelectedHisto{fPathRec + "hIsPhotonSelected"};
  std::string fFullNameCutsOnMcTruthInfoHisto{fPathMc + "hCutsOnMcTruthInfo"};

  // think of better name
  std::vector<MyHistogramSpec> fSpecialHistoDefinitions{
    {fFullNameIsPhotonSelectedHisto.data(), "hIsPhotonSelected;cut categories;counts", {HistType::kTH1F, {{13, -0.0f, 12.5f}}}},
    {fFullNameCutsOnMcTruthInfoHisto.data(), "hCutsOnMcTruthInfo;cut categories;counts", {HistType::kTH1F, {{13, -0.0f, 12.5f}}}, true /*dataOnly_*/}};

  HistogramRegistry fHistogramRegistry{"fHistogramRegistry"};

  std::map<std::string, size_t> fPhotonCutIndeces{
    {"kPhotonIn", 0},
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
    {"kPhotonOut", 12}};

  std::map<std::string, size_t> fMcPhotonCutIndeces{
    {"kMcPhotonIn", 0},
    {"kMcPhysicalPrimary", 1},
    {"kMcEta", 2},
    {"kMcPhotonOut", 12}};

  // helper function for init::addHistosToRegistry. Can't be a lambda itself because it's templated. (possible in c++20)
  template <typename T>
  void appendSuffixToTitle(HistPtr& theHistPtr, std::string* theSuffix)
  {
    auto lHisto = std::get<std::shared_ptr<T>>(theHistPtr);
    if (lHisto) {
      std::string lTitle(lHisto->GetTitle());
      lHisto->SetTitle(lTitle.append(*theSuffix).data());
    } else {
      LOGF(info, "SFS appendSuffixToTitle(): %s could not be obtained in order to append suffix to title.");
    }
  }

  void init(InitContext const&)
  {
    auto addHistosToRegistry = [&](auto& theContainer, auto const& theHistoDefinitions, std::string const& thePath, std::string* theSuffix = nullptr) {
      for (auto& tHisto : theHistoDefinitions) {
        if (!doprocessMc && tHisto.fDataOnly) {
          continue;
        }
        std::string lFullName(thePath + tHisto.m.name + (theSuffix ? *theSuffix : std::string("")));
        LOGF(info, "adding %s %d", lFullName, tHisto.fDataOnly);
        HistPtr lHistPtr = fHistogramRegistry.add(lFullName.data(), tHisto.m.title.data(), tHisto.m.config);
        theContainer.insert(std::pair{lFullName, lHistPtr});

        if (theSuffix) {
          if (tHisto.m.config.type == kTH1F) {
            appendSuffixToTitle<TH1>(lHistPtr, theSuffix);
          } else if (tHisto.m.config.type == kTH2F) {
            appendSuffixToTitle<TH2>(lHistPtr, theSuffix);
          }
        }
      }
    };

    auto addLablesToHisto = [](auto const& theContainer, std::string const& theHistoName, auto const& theLables) {
      auto lHisto = theContainer.find(theHistoName);
      if (lHisto != theContainer.end()) {
        TAxis* lXaxis = std::get<std::shared_ptr<TH1>>(lHisto->second)->GetXaxis();
        for (auto& lPairIt : theLables) {
          lXaxis->SetBinLabel(lPairIt.second + 1, lPairIt.first.data());
        }
      }
    };

    // todo: clean up
    if (doprocessRec) {
      for (auto& lString : fRecTrueStrings) {
        lString.replace(1, 2, std::string("")); // remove 'MC' in '_MC*'
      }
    }

    // add single Histograms
    addHistosToRegistry(fSpecialHistos, fSpecialHistoDefinitions, std::string(""));

    // do some labeling
    addLablesToHisto(fSpecialHistos, fFullNameIsPhotonSelectedHisto, fPhotonCutIndeces);

    for (auto bac : std::vector<std::string>{"beforeRecCuts/", "afterRecCuts/"}) {

      // collision histograms
      addHistosToRegistry(fCollisionHistos,
                          fCollisionHistoDefinitions,
                          fPathRecCollisionHistos + bac,
                          &fRecTrueStrings[kRec]);

      // track histograms
      addHistosToRegistry(fTrackHistos,
                          fTrackHistoDefinitions,
                          fPathRecTrackHistos + bac,
                          &fRecTrueStrings[kRec]);

      // v0 histograms
      std::vector<eRecTrueEnum> lRecTrue{kRec};
      if (doprocessMc) {
        lRecTrue.push_back(kMCTrue);
        lRecTrue.push_back(kMCVal);

        // v0 Resolution histos
        addHistosToRegistry(fV0ResolutionHistos,
                            fV0ResolutionHistoDefinitions,
                            fPathResolutions + bac);

        addLablesToHisto(fSpecialHistos, fFullNameCutsOnMcTruthInfoHisto, fMcPhotonCutIndeces);
      }

      for (auto iRecTrue : lRecTrue) {
        addHistosToRegistry(fRecTrueV0Histos[iRecTrue],
                            fV0HistoDefinitions,
                            fPathsV0Histos[iRecTrue] + bac,
                            &fRecTrueStrings[iRecTrue]);
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

  template <typename TV0, typename TTRACKS>
  bool processPhoton(TV0 const& theV0, const float& theV0CosinePA, TTRACKS const& theTracks)
  {
    auto lV0Tracks = theTracks.sliceBy(aod::v0data::v0Id, theV0.v0Id());

    fillReconstructedInfoHistograms("beforeRecCuts/",
                                    theV0,
                                    lV0Tracks,
                                    theV0CosinePA);

    for (auto& lTrack : lV0Tracks) {
      if (!trackPassesCuts(lTrack)) {
        return kFALSE;
      }
    }

    // apply photon cuts
    if (!photonPassesCuts(theV0, theV0CosinePA)) {
      return kFALSE;
    }

    fillReconstructedInfoHistograms("afterRecCuts/",
                                    theV0,
                                    lV0Tracks,
                                    theV0CosinePA);
    return kTRUE;
  }

  template <typename TV0, typename TMCGAMMATABLE>
  bool processMcPhoton(TMCGAMMATABLE const& theMcPhotonForThisV0AsTable, TV0 const& theV0, float const& theV0CosinePA)
  {
    // is a table that might be empty
    if (theMcPhotonForThisV0AsTable.begin() == theMcPhotonForThisV0AsTable.end()) {
      return false;
    }

    auto const lMcPhoton = theMcPhotonForThisV0AsTable.begin();
    fillTH1(fSpecialHistos, fFullNameCutsOnMcTruthInfoHisto, getPhotonCutIndex("kMcPhotonIn", true /*theMcCuts*/));

    if (!lMcPhoton.isPhysicalPrimary()) {
      fillTH1(fSpecialHistos, fFullNameCutsOnMcTruthInfoHisto, getPhotonCutIndex("kMcPhysicalPrimary", true /*theMcCuts*/));
      return false;
    }

    if (std::abs(lMcPhoton.eta()) > fTruePhotonEtaMax) {
      fillTH1(fSpecialHistos, fFullNameCutsOnMcTruthInfoHisto, getPhotonCutIndex("kMcEta", true /*theMcCuts*/));
      return false;
    }
    fillTH1(fSpecialHistos, fFullNameCutsOnMcTruthInfoHisto, getPhotonCutIndex("kMcPhotonOut", true /*theMcCuts*/));

    fillTruePhotonHistograms(lMcPhoton,
                             theV0,
                             theV0CosinePA,
                             "beforeRecCuts/");
    return true;
  }

  template <typename TV0, typename TMCGAMMA>
  void fillTruePhotonHistograms(TMCGAMMA const& theMcPhoton, TV0 const& theV0, float const& theV0CosinePA, std::string theBAC)
  {
    fillV0HistogramsMcGamma(kMCTrue, theBAC, theMcPhoton);
    fillV0Histograms(kMCVal, theBAC, theV0, theV0CosinePA);

    // v0 resolution histos
    {
      TVector3 lConvPointRec(theV0.x(), theV0.y(), theV0.z());
      TVector3 lConvPointTrue(theMcPhoton.conversionX(), theMcPhoton.conversionY(), theMcPhoton.conversionZ());

      std::string lPath(fPathResolutions + theBAC);
      fillTH1(fV0ResolutionHistos, lPath + "hPtRes", theV0.pt() - theMcPhoton.pt());
      fillTH1(fV0ResolutionHistos, lPath + "hEtaRes", theV0.eta() - theMcPhoton.eta());
      fillTH1(fV0ResolutionHistos, lPath + "hPhiRes", theV0.phi() - theMcPhoton.phi());
      fillTH1(fV0ResolutionHistos, lPath + "hConvPointRRes", theV0.v0radius() - lConvPointTrue.Perp());
      fillTH1(fV0ResolutionHistos, lPath + "hConvPointAbsoluteDistanceRes", TVector3(lConvPointRec - lConvPointTrue).Mag());
    }
  }

  void processRec(aod::Collisions::iterator const& theCollision,
                  aod::V0Datas const& theV0s,
                  aod::V0DaughterTracks const& theTracks)
  {
    fillTH1(fCollisionHistos, fPathRecCollisionHistos + "beforeRecCuts/hCollisionZ" + fRecTrueStrings[kRec], theCollision.posZ());

    for (auto& lV0 : theV0s) {

      float lV0CosinePA = lV0.v0cosPA(theCollision.posX(), theCollision.posY(), theCollision.posZ());

      if (!processPhoton(lV0, lV0CosinePA, theTracks)) {
        continue;
      }
    }
  }
  PROCESS_SWITCH(GammaConversions, processRec, "process reconstructed info", true);

  void processMc(aod::Collisions::iterator const& theCollision,
                 aod::V0Datas const& theV0s,
                 aod::V0DaughterTracks const& theTracks,
                 aod::McGammasTrue const& theV0sTrue)
  {
    fillTH1(fCollisionHistos, fPathRecCollisionHistos + "beforeRecCuts/hCollisionZ" + fRecTrueStrings[kRec], theCollision.posZ());

    for (auto& lV0 : theV0s) {

      float lV0CosinePA = lV0.v0cosPA(theCollision.posX(), theCollision.posY(), theCollision.posZ());

      // is a table that might be empty
      auto lMcPhotonForThisV0AsTable = theV0sTrue.sliceBy(aod::v0data::v0Id, lV0.v0Id());
      bool lComesFromMcPhotonInAcceptance = processMcPhoton(lMcPhotonForThisV0AsTable, lV0, lV0CosinePA);

      if (!processPhoton(lV0, lV0CosinePA, theTracks)) {
        continue;
      }

      if (lComesFromMcPhotonInAcceptance) {
        fillTruePhotonHistograms(lMcPhotonForThisV0AsTable.begin(),
                                 lV0,
                                 lV0CosinePA,
                                 "afterRecCuts/");
      }
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
  void fillTrackHistograms(std::string const& theBAC, TTRACKS const& theV0Tracks)
  {
    std::string lPath = fPathRecTrackHistos + theBAC;
    std::string& lSuffix = fRecTrueStrings[kRec];
    auto fillTrackHistogramsI = [&](auto const& theTrack) {
      fillTH1(fTrackHistos, lPath + "hTrackEta" + lSuffix, theTrack.eta());
      fillTH1(fTrackHistos, lPath + "hTrackPhi" + lSuffix, theTrack.phi());
      fillTH1(fTrackHistos, lPath + "hTrackPt" + lSuffix, theTrack.pt());
      fillTH1(fTrackHistos, lPath + "hTPCFoundOverFindableCls" + lSuffix, theTrack.tpcFoundOverFindableCls());
      fillTH1(fTrackHistos, lPath + "hTPCCrossedRowsOverFindableCls" + lSuffix, theTrack.tpcCrossedRowsOverFindableCls());
      fillTH2(fTrackHistos, lPath + "hTPCdEdxSigEl" + lSuffix, theTrack.p(), theTrack.tpcNSigmaEl());
      fillTH2(fTrackHistos, lPath + "hTPCdEdxSigPi" + lSuffix, theTrack.p(), theTrack.tpcNSigmaPi());
      fillTH2(fTrackHistos, lPath + "hTPCdEdx" + lSuffix, theTrack.p(), theTrack.tpcSignal());
    };

    for (auto& lTrack : theV0Tracks) {
      fillTrackHistogramsI(lTrack);
    }
  }

  template <typename TV0>
  void fillV0Histograms(eRecTrueEnum theRecTrue, std::string const& theBAC, TV0 const& theV0, float const& theV0CosinePA)
  {
    mapStringHistPtr& lContainer = fRecTrueV0Histos[theRecTrue];
    std::string lPath = fPathsV0Histos[theRecTrue] + theBAC;
    std::string& lSuffix = fRecTrueStrings[theRecTrue];
    auto fullName = [&lPath, &lSuffix](std::string theName) {
      return lPath + theName + lSuffix;
    };
    fillTH1(lContainer, fullName("hEta"), theV0.eta());
    fillTH1(lContainer, fullName("hPhi"), theV0.phi());
    fillTH1(lContainer, fullName("hPt"), theV0.pt());
    fillTH1(lContainer, fullName("hConvPointR"), theV0.v0radius());
    fillTH1(lContainer, fullName("hCosPAngle"), theV0CosinePA);
    fillTH2(lContainer, fullName("hArmenteros"), theV0.alpha(), theV0.qtarm());
    fillTH2(lContainer, fullName("hPsiPt"), theV0.psipair(), theV0.pt());
  }

  // SFS todo: combine fillV0Histograms and fillV0HistogramsMcGamma
  template <typename TMCGAMMA>
  void fillV0HistogramsMcGamma(eRecTrueEnum theRecTrue, std::string const& theBAC, TMCGAMMA const& theMcGamma)
  {
    mapStringHistPtr& lContainer = fRecTrueV0Histos[theRecTrue];
    std::string lPath = fPathsV0Histos[theRecTrue] + theBAC;
    std::string& lSuffix = fRecTrueStrings[theRecTrue];
    auto fullName = [&lPath, &lSuffix](std::string theName) {
      return lPath + theName + lSuffix;
    };
    fillTH1(lContainer, fullName("hEta"), theMcGamma.eta());
    fillTH1(lContainer, fullName("hPhi"), theMcGamma.phi());
    fillTH1(lContainer, fullName("hPt"), theMcGamma.pt());
    fillTH1(lContainer, fullName("hConvPointR"), theMcGamma.v0Radius());
  }

  template <typename T>
  bool trackPassesCuts(const T& theTrack)
  {
    // single track eta cut
    if (TMath::Abs(theTrack.eta()) > fTrackEtaMax) {
      fillTH1(fSpecialHistos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kTrackEta"));
      return kFALSE;
    }

    // single track pt cut
    if (theTrack.pt() < fTrackPtMin) {
      fillTH1(fSpecialHistos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kTrackPt"));
      return kFALSE;
    }

    if (!(selectionPIDTPC_track(theTrack))) {
      return kFALSE;
    }

    if (theTrack.tpcFoundOverFindableCls() < fMinTPCFoundOverFindableCls) {
      fillTH1(fSpecialHistos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kTPCFoundOverFindableCls"));
      return kFALSE;
    }

    if (theTrack.tpcCrossedRowsOverFindableCls() < fMinTPCCrossedRowsOverFindableCls) {
      fillTH1(fSpecialHistos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kTPCCrossedRowsOverFindableCls"));
      return kFALSE;
    }
    return kTRUE;
  }

  template <typename T>
  bool photonPassesCuts(const T& theV0, float theV0CosinePA)
  {
    if (theV0.v0radius() < fV0RMin || theV0.v0radius() > fV0RMax) {
      fillTH1(fSpecialHistos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kV0Radius"));
      return kFALSE;
    }

    if (fV0PhotonAsymmetryMax > 0. && !ArmenterosQtCut(theV0.alpha(), theV0.qtarm(), theV0.pt())) {
      fillTH1(fSpecialHistos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kArmenteros"));
      return kFALSE;
    }

    if (fV0PsiPairMax > 0. && TMath::Abs(theV0.psipair()) > fV0PsiPairMax) {
      fillTH1(fSpecialHistos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kPsiPair"));
      return kFALSE;
    }

    if (fV0CosPAngleMin > 0. && theV0CosinePA < fV0CosPAngleMin) {
      fillTH1(fSpecialHistos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kCosinePA"));
      return kFALSE;
    }

    return kTRUE;
  }

  template <typename TV0, typename TTRACKS>
  void fillReconstructedInfoHistograms(std::string theBAC, TV0 const& theV0, TTRACKS const& theV0Tracks, float const& theV0CosinePA)
  {
    fillTrackHistograms(theBAC, theV0Tracks);
    fillV0Histograms(kRec, theBAC, theV0, theV0CosinePA);

    if (theBAC == "beforeRecCuts/") {
      fillTH1(fSpecialHistos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kPhotonIn"));
    } else if (theBAC == "afterRecCuts/") {
      fillTH1(fSpecialHistos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kPhotonOut"));
    }
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
      fillTH1(fSpecialHistos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kElectronPID"));
      return kFALSE;
    }

    // TPC Pion Line
    if (theTrack.p() > fPIDPionRejectionPMin) {
      // low pt Pion rej
      if (theTrack.p() < fPIDPionRejectionPBoarder) {
        if (theTrack.tpcNSigmaEl() > fPIDnSigmaElectronMin && theTrack.tpcNSigmaEl() < fPIDnSigmaElectronMax && theTrack.tpcNSigmaPi() < fPIDnSigmaAbovePionLineLowPMin) {
          fillTH1(fSpecialHistos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kPionRejLowMom"));
          return kFALSE;
        }
      }
      // High Pt Pion rej
      else {
        if (theTrack.tpcNSigmaEl() > fPIDnSigmaElectronMin && theTrack.tpcNSigmaEl() < fPIDnSigmaElectronMax && theTrack.tpcNSigmaPi() < fPIDnSigmaAbovePionLineHighPMin) {
          fillTH1(fSpecialHistos, fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kPionRejHighMom"));
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
