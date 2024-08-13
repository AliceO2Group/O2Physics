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

/// \file bjetTreeMerger.cxx
/// \brief Task for merging time frames with b-jet information. Can not be used on hyperloop
///
/// \author Hadi Hassan <hadi.hassan@cern.ch>, University of Jyväskylä

#include "TTree.h"
#include "TFile.h"

#include <iostream>
#include <chrono>
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/runDataProcessing.h"
#include "Framework/ASoAHelpers.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace jetInfo
{
// DECLARE_SOA_INDEX_COLUMN(JetIndex, jetindex); //! The jet index
DECLARE_SOA_COLUMN(JetpT, jetpt, float);       //! jet pT
DECLARE_SOA_COLUMN(JetEta, jeteta, float);     //! jet eta
DECLARE_SOA_COLUMN(JetPhi, jetphi, float);     //! jet phi
DECLARE_SOA_COLUMN(NTracks, nTracks, int16_t); //! number of charged tracks inside the jet
DECLARE_SOA_COLUMN(NSV, nSV, int16_t);         //! Number of secondary vertices in the jet
DECLARE_SOA_COLUMN(JetMass, mass, float);      //! The jet mass
DECLARE_SOA_COLUMN(JetFlavor, jetFl, int16_t); //! The jet flavor (b, c, or lf)
DECLARE_SOA_COLUMN(JetR, jetR, int16_t);       //! The jet radius
} // namespace jetInfo

DECLARE_SOA_TABLE(bjetParams, "AOD", "BJETPARAM",
                  o2::soa::Index<>,
                  jetInfo::JetpT,
                  jetInfo::JetEta,
                  jetInfo::JetPhi,
                  jetInfo::NTracks,
                  jetInfo::NSV,
                  jetInfo::JetMass,
                  jetInfo::JetFlavor,
                  jetInfo::JetR);

using bjetParam = bjetParams::iterator;

namespace trackInfo
{
DECLARE_SOA_INDEX_COLUMN(bjetParam, jetindex);                         //! The jet index
DECLARE_SOA_COLUMN(TrackpT, trackpt, float);                           //! The track pT
DECLARE_SOA_COLUMN(TrackEta, tracketa, float);                         //! The track eta
DECLARE_SOA_COLUMN(DotProdTrackJet, trackdotjet, float);               //! The dot product between the track and the jet
DECLARE_SOA_COLUMN(DotProdTrackJetOverJet, trackdotjetoverjet, float); //! The dot product between the track and the jet over the jet momentum
DECLARE_SOA_COLUMN(DeltaRJetTrack, rjettrack, float);                  //! The DR jet-track
DECLARE_SOA_COLUMN(SignedIP2D, ip2d, float);                           //! The track signed 2D IP
DECLARE_SOA_COLUMN(SignedIP2DSign, ip2dsigma, float);                  //! The track signed 2D IP significance
DECLARE_SOA_COLUMN(SignedIP3D, ip3d, float);                           //! The track signed 3D IP
DECLARE_SOA_COLUMN(SignedIP3DSign, ip3dsigma, float);                  //! The track signed 3D IP significance
DECLARE_SOA_COLUMN(MomFraction, momfraction, float);                   //! The track momentum fraction of the jets
DECLARE_SOA_COLUMN(DeltaRTrackVertex, rtrackvertex, float);            //! DR between the track and the closest SV, to be decided whether to add to or not
// DECLARE_SOA_COLUMN(DCATrackJet, dcatrackjet, float);                              //! The distance between track and jet, unfortunately it cannot be calculated in O2
} // namespace trackInfo

DECLARE_SOA_TABLE(bjetTracksParams, "AOD", "BJETTRACKSPARAM",
                  o2::soa::Index<>,
                  trackInfo::bjetParamId,
                  trackInfo::TrackpT,
                  trackInfo::TrackEta,
                  trackInfo::DotProdTrackJet,
                  trackInfo::DotProdTrackJetOverJet,
                  trackInfo::DeltaRJetTrack,
                  trackInfo::SignedIP2D,
                  trackInfo::SignedIP2DSign,
                  trackInfo::SignedIP3D,
                  trackInfo::SignedIP3DSign,
                  trackInfo::MomFraction,
                  trackInfo::DeltaRTrackVertex);

using bjetTracksParam = bjetTracksParams::iterator;

namespace SVInfo
{
DECLARE_SOA_INDEX_COLUMN(bjetParam, jetindex);            //! The jet index
DECLARE_SOA_COLUMN(SVpT, svpt, float);                    //! The SV pT
DECLARE_SOA_COLUMN(DeltaRSVJet, rsvjet, float);           //! The DR jet-SV
DECLARE_SOA_COLUMN(SVMass, mass, float);                  //! The SV mass
DECLARE_SOA_COLUMN(SVfE, svfe, float);                    //! The SV energy fraction
DECLARE_SOA_COLUMN(IPXY, ipxy, float);                    //! The SV 2D IP
DECLARE_SOA_COLUMN(CPA, cpa, float);                      //! Cosine pointing angle between the SV direction and momentum
DECLARE_SOA_COLUMN(Chi2PCA, chi2pca, float);              //! Sum of (non-weighted) distances of the secondary vertex to its prongsm
DECLARE_SOA_COLUMN(DecayLength2D, lxy, float);            //! The decay length of the SV in XY
DECLARE_SOA_COLUMN(DecayLength2DError, lxysigma, float);  //! The decay length of the SV in XY significance
DECLARE_SOA_COLUMN(DecayLength3D, lxyz, float);           //! The decay length of the SV in 3D
DECLARE_SOA_COLUMN(DecayLength3DError, lxyzsigma, float); //! The decay length of the SV in 3d significance
// DECLARE_SOA_COLUMN(SVDispersion, svdispersion, float);                              //! The SV dispersion, unfortunately it cannot be calculated in O2
} // namespace SVInfo

DECLARE_SOA_TABLE(bjetSVParams, "AOD", "BJETSVPARAM",
                  o2::soa::Index<>,
                  SVInfo::bjetParamId,
                  SVInfo::SVpT,
                  SVInfo::DeltaRSVJet,
                  SVInfo::SVMass,
                  SVInfo::SVfE,
                  SVInfo::IPXY,
                  SVInfo::CPA,
                  SVInfo::Chi2PCA,
                  SVInfo::DecayLength2D,
                  SVInfo::DecayLength2DError,
                  SVInfo::DecayLength3D,
                  SVInfo::DecayLength3DError);

using bjetSVParam = bjetSVParams::iterator;

namespace constituents
{
DECLARE_SOA_INDEX_COLUMN(bjetParam, jetindex);
DECLARE_SOA_ARRAY_INDEX_COLUMN(bjetTracksParam, tracks);
DECLARE_SOA_ARRAY_INDEX_COLUMN(bjetSVParam, svs);
} // namespace constituents

DECLARE_SOA_TABLE(bjetConstituents, "AOD", "BJETCONSTIT",
                  constituents::bjetParamId,
                  constituents::bjetTracksParamIds,
                  constituents::bjetSVParamIds);

} // namespace o2::aod

struct bjetTree {
  float mJetpT;
  float mJetEta;
  float mJetPhi;
  int mNTracks;
  int mNSV;
  float mJetMass;
  int mJetFlavor;

  std::array<float, 10> mTrackpT;
  std::array<float, 10> mTrackEta;
  std::array<float, 10> mDotProdTrackJet;
  std::array<float, 10> mDotProdTrackJetOverJet;
  std::array<float, 10> mDeltaRJetTrack;
  std::array<float, 10> mSignedIP2D;
  std::array<float, 10> mSignedIP2DSign;
  std::array<float, 10> mSignedIP3D;
  std::array<float, 10> mSignedIP3DSign;
  std::array<float, 10> mMomFraction;
  std::array<float, 10> mDeltaRTrackVertex;

  std::array<float, 10> mSVpT;
  std::array<float, 10> mDeltaRSVJet;
  std::array<float, 10> mSVMass;
  std::array<float, 10> mSVfE;
  std::array<float, 10> mIPXY;
  std::array<float, 10> mCPA;
  std::array<float, 10> mChi2PCA;
  std::array<float, 10> mDecayLength2D;
  std::array<float, 10> mDecayLength2DError;
  std::array<float, 10> mDecayLength3D;
  std::array<float, 10> mDecayLength3DError;
};

struct BJetMerger {

  bjetTree treeWords;

  OutputObj<TTree> myTree{"myTree"};

  // TFile* myFile;
  // TTree* myTree;

  std::chrono::steady_clock::time_point start;

  void init(InitContext const&)
  {
    myTree.setObject(new TTree("myTree", "Test tree"));
    start = std::chrono::steady_clock::now();

    // myFile = new TFile("myFile.root", "RECREATE");
    // myTree = new TTree("myTree", "Test tree");

    // myTree->SetMaxTreeSize(Long64_t(1e9));

    // Jet info branches
    myTree->Branch("mJetpT", &treeWords.mJetpT, "mJetpT/F");
    myTree->Branch("mJetEta", &treeWords.mJetEta, "mJetEta/F");
    myTree->Branch("mJetPhi", &treeWords.mJetPhi, "mJetPhi/F");
    myTree->Branch("mNTracks", &treeWords.mNTracks, "mNTracks/I");
    myTree->Branch("mNSV", &treeWords.mNSV, "mNSV/I");
    myTree->Branch("mJetMass", &treeWords.mJetMass, "mJetMass/F");
    myTree->Branch("mJetFlavor", &treeWords.mJetFlavor, "mJetFlavor/I");
    // Track info branches
    myTree->Branch("mTrackpT", &treeWords.mTrackpT, "mTrackpT[10]/F");
    myTree->Branch("mTrackEta", &treeWords.mTrackEta, "mTrackEta[10]/F");
    myTree->Branch("mDotProdTrackJet", &treeWords.mDotProdTrackJet, "mDotProdTrackJet[10]/F");
    myTree->Branch("mDotProdTrackJetOverJet", &treeWords.mDotProdTrackJetOverJet, "mDotProdTrackJetOverJet[10]/F");
    myTree->Branch("mDeltaRJetTrack", &treeWords.mDeltaRJetTrack, "mDeltaRJetTrack[10]/F");
    myTree->Branch("mSignedIP2D", &treeWords.mSignedIP2D, "mSignedIP2D[10]/F");
    myTree->Branch("mSignedIP2DSign", &treeWords.mSignedIP2DSign, "mSignedIP2DSign[10]/F");
    myTree->Branch("mSignedIP3D", &treeWords.mSignedIP3D, "mSignedIP3D[10]/F");
    myTree->Branch("mSignedIP3DSign", &treeWords.mSignedIP3DSign, "mSignedIP3DSign[10]/F");
    myTree->Branch("mMomFraction", &treeWords.mMomFraction, "mMomFraction[10]/F");
    myTree->Branch("mDeltaRTrackVertex", &treeWords.mDeltaRTrackVertex, "mDeltaRTrackVertex[10]/F");
    // Secondary vertex info branches
    myTree->Branch("mSVpT", &treeWords.mSVpT, "mSVpT[10]/F");
    myTree->Branch("mDeltaRSVJet", &treeWords.mDeltaRSVJet, "mDeltaRSVJet[10]/F");
    myTree->Branch("mSVMass", &treeWords.mSVMass, "mSVMass[10]/F");
    myTree->Branch("mSVfE", &treeWords.mSVfE, "mSVfE[10]/F");
    myTree->Branch("mIPXY", &treeWords.mIPXY, "mIPXY[10]/F");
    myTree->Branch("mCPA", &treeWords.mCPA, "mCPA[10]/F");
    myTree->Branch("mChi2PCA", &treeWords.mChi2PCA, "mChi2PCA[10]/F");
    myTree->Branch("mDecayLength2D", &treeWords.mDecayLength2D, "mDecayLength2D[10]/F");
    myTree->Branch("mDecayLength2DError", &treeWords.mDecayLength2DError, "mDecayLength2DError[10]/F");
    myTree->Branch("mDecayLength3D", &treeWords.mDecayLength3D, "mDecayLength3D[10]/F");
    myTree->Branch("mDecayLength3DError", &treeWords.mDecayLength3DError, "mDecayLength3DError[10]/F");
  }

  using bjetParamsConstit = soa::Join<aod::bjetParams, aod::bjetConstituents>;

  template <typename T, std::size_t N>
  void printArray(const std::array<T, N>& arr)
  {
    std::cout << "[";
    for (std::size_t i = 0; i < N; ++i) {
      std::cout << arr[i];
      if (i < N - 1) {
        std::cout << ", ";
      }
    }
    std::cout << "]" << std::endl;
  }

  template <typename T>
  void printVector(const std::vector<T>& vec)
  {
    std::cout << "[";
    for (std::size_t i = 0; i < vec.size(); ++i) {
      std::cout << vec[i];
      if (i < vec.size() - 1) {
        std::cout << ", ";
      }
    }
    std::cout << "]" << std::endl;
  }

  void process(bjetParamsConstit const& allJets, aod::bjetTracksParams const& /*allTracks*/, aod::bjetSVParams const& /*allSVs*/)
  {

    for (auto& analysisJet : allJets) {

      treeWords.mJetpT = analysisJet.jetpt();
      treeWords.mJetEta = analysisJet.jeteta();
      treeWords.mJetPhi = analysisJet.jetphi();
      treeWords.mJetFlavor = analysisJet.jetFl();
      treeWords.mJetMass = analysisJet.mass();
      treeWords.mNTracks = analysisJet.nTracks();
      treeWords.mNSV = analysisJet.nSV();

      /*if (analysisJet.globalIndex() % 10000 == 0) {
        auto now = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
        std::cout << "Index: " << analysisJet.globalIndex() << " jet pT: " << analysisJet.jetpt() << " with nSV: " << analysisJet.template svs_as<aod::bjetSVParams>().size() << " and nTracks: " << analysisJet.template tracks_as<aod::bjetTracksParams>().size() << ", Elapsed time: " << elapsed << std::endl;
        start = now;
      }*/

      auto compareDecayLength = [](const aod::bjetSVParams::iterator& sv1, const aod::bjetSVParams::iterator& sv2) {
        return (sv1.lxy() / sv1.lxysigma()) > (sv2.lxy() / sv2.lxysigma());
      };

      auto secondaryVertices = analysisJet.template svs_as<aod::bjetSVParams>();

      std::sort(secondaryVertices.begin(), secondaryVertices.end(), compareDecayLength);

      int countSV = 0;
      for (auto SecondVertex : secondaryVertices) {

        treeWords.mSVpT[countSV] = SecondVertex.svpt();
        treeWords.mDeltaRSVJet[countSV] = SecondVertex.rsvjet();
        treeWords.mSVMass[countSV] = SecondVertex.mass();
        treeWords.mSVfE[countSV] = SecondVertex.svfe();
        treeWords.mIPXY[countSV] = SecondVertex.ipxy();
        treeWords.mCPA[countSV] = SecondVertex.cpa();
        treeWords.mChi2PCA[countSV] = SecondVertex.chi2pca();
        treeWords.mDecayLength2D[countSV] = SecondVertex.lxy();
        treeWords.mDecayLength2DError[countSV] = SecondVertex.lxysigma();
        treeWords.mDecayLength3D[countSV] = SecondVertex.lxyz();
        treeWords.mDecayLength3DError[countSV] = SecondVertex.lxyzsigma();

        countSV++;
        if (countSV >= 10) {
          break;
        }
      }

      auto compareIP = [](const aod::bjetTracksParams::iterator& track1, const aod::bjetTracksParams::iterator& track2) {
        return (track1.ip2d() / track1.ip2dsigma()) > (track2.ip2d() / track2.ip2dsigma());
      };

      auto jetTracks = analysisJet.template tracks_as<aod::bjetTracksParams>();

      std::sort(jetTracks.begin(), jetTracks.end(), compareIP);

      int coutTrack = 0;
      for (auto mytrack : jetTracks) {

        treeWords.mTrackpT[coutTrack] = mytrack.trackpt();
        treeWords.mTrackEta[coutTrack] = mytrack.tracketa();
        treeWords.mDotProdTrackJet[coutTrack] = mytrack.trackdotjet();
        treeWords.mDotProdTrackJetOverJet[coutTrack] = mytrack.trackdotjetoverjet();
        treeWords.mDeltaRJetTrack[coutTrack] = mytrack.rjettrack();
        treeWords.mSignedIP2D[coutTrack] = mytrack.ip2d();
        treeWords.mSignedIP2DSign[coutTrack] = mytrack.ip2dsigma();
        treeWords.mSignedIP3D[coutTrack] = mytrack.ip3d();
        treeWords.mSignedIP3DSign[coutTrack] = mytrack.ip3dsigma();
        treeWords.mMomFraction[coutTrack] = mytrack.momfraction();
        treeWords.mDeltaRTrackVertex[coutTrack] = mytrack.rtrackvertex();

        coutTrack++;
        if (coutTrack >= 10) {
          break;
        }
      }

      myTree->Fill();
      treeWords = bjetTree();
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<BJetMerger>(cfgc, TaskName{"bjet-tree-merger"})};
}
