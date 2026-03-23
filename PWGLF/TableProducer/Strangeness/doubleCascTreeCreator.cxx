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

#include "PWGLF/DataModel/LFDoubleCascTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "TDatabasePDG.h"

#include <random>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::MultZeqs, aod::FT0Mults>::iterator;
using FullCascades = aod::CascDataExt;
using TracksFull = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCKa>;

struct doubleCascCand {
  float ptCasc1 = -999.f; // signed pt of the cascade
  float etaCasc1 = -999.f;
  float phiCasc1 = -999.f;
  float cascDecLength1 = -999.f;
  float omegaMassCasc1 = -999.f;
  float xiMassCasc1 = -999.f;
  float cosPACasc1 = -999.f;
  float dcaBachPVCasc1 = -999.f;
  float dcaV0BachCasc1 = -999.f;
  float nSigmaKBach1 = -999.f;

  float ptCasc2 = -999.f;
  float etaCasc2 = -999.f;
  float phiCasc2 = -999.f;
  float cascDecLength2 = -999.f;
  float omegaMassCasc2 = -999.f;
  float xiMassCasc2 = -999.f;
  float cosPACasc2 = -999.f;
  float dcaBachPVCasc2 = -999.f;
  float dcaV0BachCasc2 = -999.f;
  float nSigmaKBach2 = -999.f;
  float doubleOmegaMass = -999.f;
};

struct doubleCascTreeCreator {
  Produces<o2::aod::DoubleCascTable> doubleCascTable;
  std::vector<doubleCascCand> doubleCascCands;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::vertexing::DCAFitterN<2> fitter;

  int mRunNumber;

  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  Configurable<bool> cfgSkimmedProcessing{"cfgSkimmedProcessing", false, "Skimmed dataset processing"};
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  ConfigurableAxis centAxis{"centAxis", {106, 0, 106}, "binning for the centrality"};
  ConfigurableAxis zVtxAxis{"zVtxBins", {100, -20.f, 20.f}, "Binning for the vertex z in cm"};

  // binning of (anti)lambda mass QA histograms
  ConfigurableAxis massOmegaAxis{"massOmegaAxis", {400, o2::constants::physics::MassOmegaMinus - 0.05f, o2::constants::physics::MassOmegaMinus + 0.05}, "binning for the Omega invariant-mass"};
  ConfigurableAxis massXiAxis{"massXiAxis", {400, o2::constants::physics::MassXiMinus - 0.05f, o2::constants::physics::MassXiMinus + 0.05f}, "binning for the Xi invariant-mass"};

  Configurable<float> zVtxMax{"zVtxMax", 10.0f, "maximum z position of the primary vertex"};
  Configurable<float> etaMax{"etaMax", 0.9f, "maximum eta"};
  ConfigurableAxis momAxis{"momAxisFine", {5.e2, 0.f, 5.f}, "momentum axis binning"};

  Configurable<float> cascPtMin{"cascPtMin", 1.f, "minimum (anti)casc pT (GeV/c)"};
  Configurable<float> cascPtMax{"cascPtMax", 5.f, "maximum (anti)casc pT (GeV/c)"};

  Configurable<float> minNCrossedRows{"minNCrossedRows", 100, "Minimum number of crossed TPC rows"};
  Configurable<float> minNITSClus{"minNITSClus", 0., "Minimum number of ITS clusters"};
  Configurable<float> minNTPCClus{"minNTPCClus", 80, "Minimum number of TPC clusters"};
  Configurable<float> maxNSharedTPCClus{"maxNSharedTPCClus", 5, "Maximum number of shared TPC clusters"};

  Configurable<double> minCascCosPA{"minCascCosPA", 0.99f, "Minimum cosine of the pointing angle of the cascade"};
  Configurable<float> nSigmaTPCCut{"nSigmaTPCCut", 3.f, "Number of sigmas for the TPC PID"};
  Configurable<float> dcaBachToPV{"dcaBachToPV", 0.05f, "DCA of the bachelor to the primary vertex"};
  Configurable<float> dcaV0Bach{"dcaV0Bach", 1.f, "DCA between the V0 daughters"};
  Configurable<float> mXiWindow{"mXiWindow", 0.02f, "mXiWindow"};
  Configurable<float> mOmegaWindow{"mOmegaWindow", 0.01f, "mOmegaWindow"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  template <class T>
  bool selectTrack(T const& track)
  {
    if (std::abs(track.eta()) > etaMax) {
      return false;
    }
    if (track.itsNCls() < minNITSClus ||
        track.tpcNClsFound() < minNTPCClus ||
        track.tpcNClsCrossedRows() < minNCrossedRows ||
        track.tpcNClsCrossedRows() < 0.8 * track.tpcNClsFindable() ||
        track.tpcNClsShared() > maxNSharedTPCClus) {
      return false;
    }
    return true;
  }

  template <class Bc>
  void initCCDB(Bc const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    auto timestamp = bc.timestamp();
    LOG(info) << "Retrieved GRP for timestamp " << timestamp;
    mRunNumber = bc.runNumber();
    if (cfgSkimmedProcessing) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), "fDoubleOmega,fOmegaXi");
      zorro.populateHistRegistry(histos, bc.runNumber());
    }
  }

  void init(o2::framework::InitContext&)
  {
    mRunNumber = 0;
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    zorroSummary.setObject(zorro.getZorroSummary());

    // event QA
    histos.add<TH1>("QA/zVtx", ";#it{z}_{vtx} (cm);Entries", HistType::kTH1F, {zVtxAxis});
    histos.add<TH2>("QA/massXi1", ";#it{p}_{T} (GeV/#it{c});#it{M}(#Lambda + #pi^{-}) (GeV/#it{c}^{2});Entries", HistType::kTH2F, {momAxis, massXiAxis});
    histos.add<TH2>("QA/massOmega1", ";#it{p}_{T} (GeV/#it{c});#it{M}(#Omega + #pi^{-}) (GeV/#it{c}^{2});Entries", HistType::kTH2F, {momAxis, massOmegaAxis});
    histos.add<TH2>("QA/massXi2", ";#it{p}_{T} (GeV/#it{c});#it{M}(#Lambda + #pi^{-}) (GeV/#it{c}^{2});Entries", HistType::kTH2F, {momAxis, massXiAxis});
    histos.add<TH2>("QA/massOmega2", ";#it{p}_{T} (GeV/#it{c});#it{M}(#Omega + #pi^{-}) (GeV/#it{c}^{2});Entries", HistType::kTH2F, {momAxis, massOmegaAxis});
  }

  template <class C, class T>
  bool isSelectedCasc(C const& collision, T const&, FullCascades::iterator const& casc)
  {

    auto bachelor = casc.bachelor_as<T>();
    auto posDau = casc.posTrack_as<T>();
    auto negDau = casc.negTrack_as<T>();

    if (!selectTrack(bachelor) || !selectTrack(posDau) || !selectTrack(negDau)) {
      return false;
    }
    if (casc.sign() > 0) {
      if (TMath::Abs(posDau.tpcNSigmaPi()) > nSigmaTPCCut || TMath::Abs(negDau.tpcNSigmaPr()) > nSigmaTPCCut) {
        return false;
      }
    } else if (casc.sign() < 0) {
      if (TMath::Abs(negDau.tpcNSigmaPi()) > nSigmaTPCCut || TMath::Abs(posDau.tpcNSigmaPr()) > nSigmaTPCCut) {
        return false;
      }
    }
    if (TMath::Abs(casc.dcabachtopv()) < dcaBachToPV) {
      return false;
    }
    if (TMath::Abs(casc.dcacascdaughters()) > dcaV0Bach) {
      return false;
    }
    if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < minCascCosPA) {
      return false;
    }
    if (TMath::Abs(casc.eta()) > etaMax) {
      return false;
    }
    // mass cuts
    bool massInWindow = false;
    if (casc.mOmega() > o2::constants::physics::MassOmegaMinus - mOmegaWindow && casc.mOmega() < o2::constants::physics::MassOmegaMinus + mOmegaWindow) {
      massInWindow = true;
    }
    if (casc.mXi() > o2::constants::physics::MassXiMinus - mXiWindow && casc.mXi() < o2::constants::physics::MassXiMinus + mXiWindow) {
      massInWindow = true;
    }
    if (!massInWindow) {
      return false;
    }
    return true;
  };

  template <class T>
  float doubleOmegaMass(T const&, FullCascades::iterator const& casc1, FullCascades::iterator const& casc2)
  {
    // the fake omega decay is the one with the smaller radius
    auto& fakeOmega = casc1.cascradius() < casc2.cascradius() ? casc1 : casc2;
    auto& realOmega = casc1.cascradius() < casc2.cascradius() ? casc2 : casc1;
    auto kaon = fakeOmega.bachelor_as<T>();
    float momKaon[3] = {kaon.px(), kaon.py(), kaon.pz()};
    float momLambda[3] = {fakeOmega.pxlambda(), fakeOmega.pylambda(), fakeOmega.pzlambda()};
    // now compute real Omega-lambda-kaon mass
    float momTot[3] = {momKaon[0] + momLambda[0] + realOmega.px(), momKaon[1] + momLambda[1] + realOmega.py(), momKaon[2] + momLambda[2] + realOmega.pz()};
    float eK = std::sqrt(o2::constants::physics::MassKaonCharged * o2::constants::physics::MassKaonCharged + momKaon[0] * momKaon[0] + momKaon[1] * momKaon[1] + momKaon[2] * momKaon[2]);
    float eL = std::sqrt(o2::constants::physics::MassLambda0 * o2::constants::physics::MassLambda0 + momLambda[0] * momLambda[0] + momLambda[1] * momLambda[1] + momLambda[2] * momLambda[2]);
    float eO = std::sqrt(o2::constants::physics::MassOmegaMinus * o2::constants::physics::MassOmegaMinus + realOmega.px() * realOmega.px() + realOmega.py() * realOmega.py() + realOmega.pz() * realOmega.pz());
    float eTot = eK + eL + eO;
    float mass = std::sqrt(eTot * eTot - momTot[0] * momTot[0] - momTot[1] * momTot[1] - momTot[2] * momTot[2]);
    return mass;
  }

  template <class C, class T>
  void fillDoubleCasc(C const& collision, T const& tracks, FullCascades const& cascades)
  {
    doubleCascCands.clear();

    for (auto& casc1 : cascades) {
      if (!isSelectedCasc(collision, tracks, casc1)) {
        continue;
      }
      histos.fill(HIST("QA/massXi1"), casc1.pt(), casc1.mXi());
      histos.fill(HIST("QA/massOmega1"), casc1.pt(), casc1.mOmega());
      for (auto& casc2 : cascades) {
        if (!isSelectedCasc(collision, tracks, casc2)) {
          continue;
        }
        histos.fill(HIST("QA/massXi2"), casc2.pt(), casc2.mXi());
        histos.fill(HIST("QA/massOmega2"), casc2.pt(), casc2.mOmega());

        // check that the cascades do not share any track
        std::vector<int> trackIdsCasc1 = {casc1.posTrackId(), casc1.negTrackId(), casc1.bachelorId()};
        std::vector<int> trackIdsCasc2 = {casc2.posTrackId(), casc2.negTrackId(), casc2.bachelorId()};
        bool shareTrack = false;
        for (auto id1 : trackIdsCasc1) {
          for (auto id2 : trackIdsCasc2) {
            if (id1 == id2) {
              shareTrack = true;
              break;
            }
          }
          if (shareTrack) {
            break;
          }
        }
        if (shareTrack) {
          continue;
        }

        auto bach1 = casc1.bachelor_as<T>();
        auto bach2 = casc2.bachelor_as<T>();

        doubleCascCand cand;
        cand.ptCasc1 = casc1.pt();
        cand.etaCasc1 = casc1.eta();
        cand.phiCasc1 = casc1.phi();
        cand.cascDecLength1 = std::hypot(casc1.x() - collision.posX(), casc1.y() - collision.posY(), casc1.z() - collision.posZ());
        cand.omegaMassCasc1 = casc1.mOmega();
        cand.xiMassCasc1 = casc1.mXi();
        cand.cosPACasc1 = casc1.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
        cand.dcaBachPVCasc1 = casc1.dcabachtopv();
        cand.dcaV0BachCasc1 = casc1.dcacascdaughters();
        cand.nSigmaKBach1 = bach1.tpcNSigmaKa();

        cand.ptCasc2 = casc2.pt();
        cand.etaCasc2 = casc2.eta();
        cand.phiCasc2 = casc2.phi();
        cand.cascDecLength2 = std::hypot(casc2.x() - collision.posX(), casc2.y() - collision.posY(), casc2.z() - collision.posZ());
        cand.omegaMassCasc2 = casc2.mOmega();
        cand.xiMassCasc2 = casc2.mXi();
        cand.cosPACasc2 = casc2.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
        cand.dcaBachPVCasc2 = casc2.dcabachtopv();
        cand.dcaV0BachCasc2 = casc2.dcacascdaughters();
        cand.nSigmaKBach2 = bach2.tpcNSigmaKa();

        cand.doubleOmegaMass = doubleOmegaMass(tracks, casc1, casc2);

        doubleCascCands.push_back(cand);
      }
    }
  };

  void processData(Collisions const& collision, TracksFull const& tracks, FullCascades const& cascades, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    if (!collision.sel8())
      return;

    if (std::abs(collision.posZ()) > zVtxMax)
      return;

    if (!collision.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
      return;

    if (cfgSkimmedProcessing) {
      zorro.isSelected(collision.bc_as<aod::BCsWithTimestamps>().globalBC()); /// Just let Zorro do the accounting
    }
    histos.fill(HIST("QA/zVtx"), collision.posZ());
    fillDoubleCasc(collision, tracks, cascades);

    for (auto& cand : doubleCascCands) {
      doubleCascTable(
        cand.ptCasc1,
        cand.etaCasc1,
        cand.phiCasc1,
        cand.cascDecLength1,
        cand.omegaMassCasc1,
        cand.xiMassCasc1,
        cand.cosPACasc1,
        cand.dcaBachPVCasc1,
        cand.dcaV0BachCasc1,
        cand.nSigmaKBach1,
        cand.ptCasc2,
        cand.etaCasc2,
        cand.phiCasc2,
        cand.cascDecLength2,
        cand.omegaMassCasc2,
        cand.xiMassCasc2,
        cand.cosPACasc2,
        cand.dcaBachPVCasc2,
        cand.dcaV0BachCasc2,
        cand.nSigmaKBach2,
        cand.doubleOmegaMass);
    }
  }
  PROCESS_SWITCH(doubleCascTreeCreator, processData, "process (Run 3)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<doubleCascTreeCreator>(cfgc)};
}
