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

/// \file stronglyIntensiveCorr.cxx
/// \brief Forward-backward multiplicity correlations for inclusive charged particles.
/// \author Iwona Sputowska

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TAxis.h>
#include <TH1.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct StronglyIntensiveCorr {
  // ------------------------------------------------------------------
  // Configurables
  // ------------------------------------------------------------------

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgEtaWidth{"cfgEtaWidth", 0.2f, "FB window width"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "absolute eta cut"};
  Configurable<float> cfgCutPtLower{"cfgCutPtLower", 0.2f, "Lower pT cut"};
  Configurable<float> cfgCutPtUpper{"cfgCutPtUpper", 5.0f, "Upper pT cut"};
  Configurable<float> cfgCutTrackDcaZ{"cfgCutTrackDcaZ", 2.0f, "Maximum DCAz"};

  Configurable<int> cfgITScluster{"cfgITScluster", 5, "Minimum number of ITS clusters"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 50, "Minimum number of TPC clusters"};
  Configurable<int> cfgTPCnCrossedRows{"cfgTPCnCrossedRows", 70, "Minimum number of TPC crossed rows"};
  Configurable<float> cfgCutTpcChi2NCl{"cfgCutTpcChi2NCl", 2.5f, "Maximum TPC chi2/NCl"};
  Configurable<float> cfgCutItsChi2NCl{"cfgCutItsChi2NCl", 36.0f, "Maximum ITS chi2/NCl"};

  Configurable<bool> cfgRequireGlobalTrack{"cfgRequireGlobalTrack", true, "Require isGlobalTrack OR isGlobalTrackSDD"};
  Configurable<bool> cfgUseDCAzCut{"cfgUseDCAzCut", true, "Use DCAz cut"};
  Configurable<bool> cfgSel8Trig{"cfgSel8Trig", true, "Require sel8"};

  Configurable<float> cfgCentMin{"cfgCentMin", 0.0f, "Minimum centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 90.0f, "Maximum centrality"};
  Configurable<int> cfgCentralityChoice{"cfgCentralityChoice", 0, "Centrality estimator: 0=FT0C, 1=FT0M"};
  Configurable<float> cfgCentWindowWidth{"cfgCentWindowWidth", 10.0f, "Centrality window width around class centers"};

  Configurable<int> cfgNSubsamples{"cfgNSubsamples", 20, "Number of subsamples; max is NSubsamplesMax"};

  Configurable<bool> cfgEvSelkNoSameBunchPileup{"cfgEvSelkNoSameBunchPileup", true, "Pileup removal"};
  Configurable<bool> cfgUseGoodITSLayerAllCut{"cfgUseGoodITSLayerAllCut", true, "Remove time intervals with bad ITS layers"};
  Configurable<bool> cfgEvSelUseGoodZvtxFT0vsPV{"cfgEvSelUseGoodZvtxFT0vsPV", true, "Good z-vertex FT0 vs PV cut"};

  // Eta-vertex-phi QA and occupancy selection.
  Configurable<bool> cfgUseNoCollInTimeRangeStandard{"cfgUseNoCollInTimeRangeStandard", true, "Reject nearby collisions"};
  Configurable<bool> cfgUseOccupancyCut{"cfgUseOccupancyCut", true, "Apply track occupancy cut"};
  Configurable<int> cfgMinOccupancy{"cfgMinOccupancy", 500, "Minimum track occupancy (inclusive)"};
  Configurable<int> cfgMaxOccupancy{"cfgMaxOccupancy", 3000, "Maximum track occupancy (inclusive)"};
  Configurable<bool> cfgUseNoITSROFrameBorder{"cfgUseNoITSROFrameBorder", true, "Reject ITS readout frame borders"};
  Configurable<bool> cfgUseNoTimeFrameBorder{"cfgUseNoTimeFrameBorder", true, "Reject time frame borders"};
  Configurable<bool> cfgFillEtaVtxQA{"cfgFillEtaVtxQA", true, "Fill centrality-resolved eta versus vertex QA for data"};
  Configurable<bool> cfgFillEtaVtxPhiQA{"cfgFillEtaVtxPhiQA", true, "Fill centrality-resolved eta, vertex and phi QA for data"};

  Configurable<bool> cfgUseOnlyPhysicalPrimary{"cfgUseOnlyPhysicalPrimary", true, "MC truth: require physical primary"};
  Configurable<bool> cfgUseOnlyChargedMC{"cfgUseOnlyChargedMC", true, "MC truth: require charged particle"};

  // ------------------------------------------------------------------
  // FB binning
  // ------------------------------------------------------------------
  static constexpr int NEtaGaps = 8;
  static constexpr int NPtBins = 6;
  static constexpr int NPhiBins = 16;
  static constexpr int NSubsamplesMax = 20;
  static constexpr int NCentClasses = 8;
  static constexpr double TwoPi = 6.28318530717958647692;

  std::array<double, NEtaGaps> etaCenter = {0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0};

  std::array<double, NPtBins + 1> ptEdges = {0.2, 0.5, 0.8, 1.0, 1.5, 2.0, 5.0};
  std::array<double, NCentClasses + 1> centEdges = {0., 10., 20., 30., 40., 50., 60., 70., 80.};

  using EtaPtPhiArray = std::array<std::array<std::array<double, NPhiBins>, NPtBins>, NEtaGaps>;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Service<o2::framework::O2DatabasePDG> pdg{};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;

  // Intentionally loose track filter: manual cuts are done in selTrack(), so QA before/after sees tracks.
  Filter trackFilter = (aod::track::pt > -1.0f);

  using AodCollisions = soa::Filtered<soa::Join<
    aod::Collisions,
    aod::EvSels,
    aod::CentFT0Cs,
    aod::CentFT0Ms,
    aod::Mults>>;
  using AodCollision = AodCollisions::iterator;

  using AodTracks = soa::Filtered<soa::Join<aod::Tracks,
                                            aod::TrackSelection,
                                            aod::TracksExtra,
                                            aod::TracksDCA>>;

  using MyMCRecCollisions = soa::Filtered<soa::Join<
    aod::Collisions,
    aod::EvSels,
    aod::CentFT0Cs,
    aod::CentFT0Ms,
    aod::Mults,
    aod::McCollisionLabels>>;
  using MyMCRecCollision = MyMCRecCollisions::iterator;

  using MyMCTracks = soa::Filtered<soa::Join<aod::Tracks,
                                             aod::TracksExtra,
                                             aod::TracksDCA,
                                             aod::TrackSelection,
                                             aod::McTrackLabels>>;

  Preslice<MyMCTracks> perCollision = aod::track::collisionId;
  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;

  // ------------------------------------------------------------------
  // Init
  // ------------------------------------------------------------------
  void init(InitContext&)
  {
    if (cfgUseOccupancyCut.value && (cfgMinOccupancy.value < 0 || cfgMaxOccupancy.value < cfgMinOccupancy.value)) {
      LOGF(fatal, "Invalid occupancy interval: require 0 <= cfgMinOccupancy <= cfgMaxOccupancy");
    }

    AxisSpec centAxis = {100, 0.0, 100.0, "centrality (%)"};
    std::vector<double> centClassEdges(centEdges.begin(), centEdges.end());
    AxisSpec centClassAxis = {centClassEdges, "centrality class (%)"};
    AxisSpec etaGapAxis = {NEtaGaps, -0.05, 1.55, "#Delta#eta"};

    // Seven physical, non-overlapping FB gaps:
    // DeltaEta = 0.0, 0.2, ..., 1.2.
    std::vector<double> physicalGapEdges = {-0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3};

    AxisSpec physicalGapAxis = {physicalGapEdges, "#Delta#eta"};
    AxisSpec nAxis = {500, 0.0, 500.0, "N"};
    AxisSpec counterAxis = {20, 0.5, 20.5, "counter"};
    AxisSpec ptAxis = {{0.2, 0.5, 0.8, 1.0, 1.5, 2.0, 5.0}, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec phiAxis = {NPhiBins, 0.0, TwoPi, "#varphi"};
    AxisSpec subsampleAxis = {NSubsamplesMax, -0.5, static_cast<double>(NSubsamplesMax) - 0.5, "subsample"};

    AxisSpec nchAxis = {500, 0.0, 500.0, "N_{ch}"};
    AxisSpec nchPvAxis = {1000, 0.0, 1000.0, "N_{ch}^{PV}"};
    AxisSpec itsClsAxis = {21, -0.5, 20.5, "N ITS clusters"};
    AxisSpec tpcClsAxis = {201, -0.5, 200.5, "N TPC clusters / crossed rows"};
    AxisSpec chi2ITSAxis = {100, 0.0, 50.0, "#chi^{2}_{ITS}/N_{cls}"};
    AxisSpec chi2TPCAxis = {100, 0.0, 5.0, "#chi^{2}_{TPC}/N_{cls}"};
    AxisSpec dcaXYAxis = {200, -0.5, 0.5, "DCA_{xy} (cm)"};
    AxisSpec dcaZAxis = {200, -2.0, 2.0, "DCA_{z} (cm)"};
    AxisSpec etaAxis = {160, -1.6, 1.6, "#eta"};
    AxisSpec ptFullAxis = {200, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec phiAxisFull = {64, 0.0, TwoPi, "#varphi"};
    AxisSpec chargeAxis = {5, -2.5, 2.5, "charge"};
    AxisSpec vtxZAxis = {240, -30.0, 30.0, "V_{z} (cm)"};

    // Event QA
    histos.add("QA/eventCounter", "event counter", kTH1D, {counterAxis});
    histos.add("QA/trackCounter", "track counter", kTH1D, {counterAxis});
    histos.add("QA/hTrackSliceSize", "number of tracks in collision-associated slice before manual cuts", kTH1D, {{500, -0.5, 499.5}});
    histos.add("QA/hGlobalTrackFlags", "global-track flags: 1=all, 2=isGlobalTrack, 3=isGlobalTrackSDD, 4=global OR SDD, 5=fail both", kTH1D, {counterAxis});
    histos.add("QA/hCentrality", "centrality", kTH1D, {centAxis});
    histos.add("QA/hNF", "nF", kTH1D, {nAxis});
    histos.add("QA/hNB", "nB", kTH1D, {nAxis});
    histos.add("QA/hNFvsNB", "nF vs nB", kTH2D, {nAxis, nAxis});
    histos.add("QA/fb/nFvsCent", "nF vs centrality;centrality (%);nF", kTH2D, {centAxis, nAxis});
    histos.add("QA/fb/nBvsCent", "nB vs centrality;centrality (%);nB", kTH2D, {centAxis, nAxis});
    histos.add("QA/fb/nFvsNB", "nF vs nB;nF;nB", kTH2D, {nAxis, nAxis});
    histos.add("QA/fb/nFminusNB", "nF-nB;nF-nB;Counts", kTH1D, {{500, -250.0, 250.0}});
    histos.add("QA/fb/nFplusNB", "nF+nB;nF+nB;Counts", kTH1D, {{500, 0.0, 500.0}});

    // Event-level multiplicity QA for every physical gap.
    // Axes: centrality x DeltaEta x multiplicity.
    histos.add("QA/fbAllGaps/nFvsCentGap", "Forward multiplicity;centrality (%);#Delta#eta;n_{F}", kTH3D, {centAxis, physicalGapAxis, nAxis});

    histos.add("QA/fbAllGaps/nBvsCentGap", "Backward multiplicity;centrality (%);#Delta#eta;n_{B}", kTH3D, {centAxis, physicalGapAxis, nAxis});

    histos.add("QA/mcReco/nRecoAll", "MC reco all selected;n;Counts", kTH1D, {nAxis});
    histos.add("QA/mcReco/nRecoPrimary", "MC reco primary selected;n;Counts", kTH1D, {nAxis});
    histos.add("QA/mcReco/nRecoSecondary", "MC reco secondary selected;n;Counts", kTH1D, {nAxis});
    histos.add("QA/mcReco/nRecoFake", "MC reco fake/unmatched selected;n;Counts", kTH1D, {nAxis});

    histos.add("eventQA/before/nPVeta1TracksVsCent", "before;centrality (%);N_{ch}^{PV, |#eta|<1}", kTH2D, {centAxis, nchAxis});
    histos.add("eventQA/after/nPVeta1TracksVsCent", "after;centrality (%);N_{ch}^{PV, |#eta|<1}", kTH2D, {centAxis, nchAxis});
    histos.add("eventQA/before/nPVTracksVsCent", "before;centrality (%);N_{ch}^{PV}", kTH2D, {centAxis, nchPvAxis});
    histos.add("eventQA/after/nPVTracksVsCent", "after;centrality (%);N_{ch}^{PV}", kTH2D, {centAxis, nchPvAxis});
    histos.add("eventQA/before/vtxZ", "before event selection;V_{z} (cm);Counts", kTH1D, {vtxZAxis});
    histos.add("eventQA/after/vtxZ", "after event selection;V_{z} (cm);Counts", kTH1D, {vtxZAxis});
    histos.add("eventQA/before/centFT0C", "before;FT0C centrality (%);Counts", kTH1D, {centAxis});
    histos.add("eventQA/after/centFT0C", "after;FT0C centrality (%);Counts", kTH1D, {centAxis});
    histos.add("eventQA/before/centFT0M", "before;FT0M centrality (%);Counts", kTH1D, {centAxis});
    histos.add("eventQA/after/centFT0M", "after;FT0M centrality (%);Counts", kTH1D, {centAxis});
    histos.add("eventQA/before/centFT0CvsFT0M", "before;FT0C;FT0M", kTH2D, {centAxis, centAxis});
    histos.add("eventQA/after/centFT0CvsFT0M", "after;FT0C;FT0M", kTH2D, {centAxis, centAxis});

    // The centrality coordinate is the percentile, with full 10%-wide classes.
    // The 4D QA has coarser z/eta bins to limit memory; its axes are cent,z,eta,phi.
    if (cfgFillEtaVtxQA.value || cfgFillEtaVtxPhiQA.value) {
      AxisSpec etaVtxZAxis = {400, -10.0, 10.0, "z_{vtx} (cm)"};
      histos.add("eventQA/after/vtxZvsCentClass", "Selected data events;centrality (%);z_{vtx} (cm)", kTH2D, {centClassAxis, etaVtxZAxis});
      histos.add("eventQA/after/centralityForEtaVtx", "Selected data events;centrality (%);events", kTH1D, {centAxis});
      if (cfgFillEtaVtxQA.value) {
        AxisSpec etaVtxEtaBeforeAxis = {200, -1.0, 1.0, "#eta"};
        AxisSpec etaVtxEtaAfterAxis = {160, -0.8, 0.8, "#eta"};
        histos.add("trackQA/before/etaVsVtxZCent", "Tracks before manual selection in selected data events;centrality (%);z_{vtx} (cm);#eta", kTH3D, {centClassAxis, etaVtxZAxis, etaVtxEtaBeforeAxis});
        histos.add("trackQA/after/etaVsVtxZCent", "Selected data tracks;centrality (%);z_{vtx} (cm);#eta", kTH3D, {centClassAxis, etaVtxZAxis, etaVtxEtaAfterAxis});
      }
    }
    if (cfgFillEtaVtxPhiQA.value) {
      AxisSpec etaVtxPhiZAxis = {200, -10.0, 10.0, "z_{vtx} (cm)"};
      AxisSpec etaVtxPhiEtaBeforeAxis = {100, -1.0, 1.0, "#eta"};
      AxisSpec etaVtxPhiEtaAfterAxis = {80, -0.8, 0.8, "#eta"};
      AxisSpec etaVtxPhiAxis = {36, 0.0, TwoPi, "#varphi (rad)"};
      // Dense THnD: about 108 MiB for both histograms, including flow bins.
      histos.add("trackQA/before/etaVsVtxZPhiCent", "Tracks before manual selection in selected data events;centrality (%);z_{vtx} (cm);#eta;#varphi (rad)", kTHnD, {centClassAxis, etaVtxPhiZAxis, etaVtxPhiEtaBeforeAxis, etaVtxPhiAxis});
      histos.add("trackQA/after/etaVsVtxZPhiCent", "Selected data tracks;centrality (%);z_{vtx} (cm);#eta;#varphi (rad)", kTHnD, {centClassAxis, etaVtxPhiZAxis, etaVtxPhiEtaAfterAxis, etaVtxPhiAxis});
    }

    // Fine 1D occupancy retains the exact integer cut boundaries and -1.
    // Coarser 2D occupancy keeps -1 in its own bin; >=20000 is overflow there.
    AxisSpec occupancyAxis = {20002, -1.5, 20000.5, "track occupancy"};
    std::vector<double> occupancyEdges = {-1.5, -0.5};
    for (int upper = 10; upper <= 20000; upper += 10) {
      occupancyEdges.push_back(upper - 0.5);
    }
    AxisSpec occupancyVsCentAxis = {occupancyEdges, "track occupancy"};
    histos.add("eventQA/occupancy/beforeExtraCuts", "After original event-quality cuts, before extra cuts;track occupancy;events", kTH1D, {occupancyAxis});
    histos.add("eventQA/occupancy/beforeExtraCutsVsCent", "After original event-quality cuts, before extra cuts;centrality (%);track occupancy", kTH2D, {centAxis, occupancyVsCentAxis});
    histos.add("eventQA/occupancy/beforeOccupancyCut", "After nearby-collision rejection, before occupancy cut;track occupancy;events", kTH1D, {occupancyAxis});
    histos.add("eventQA/occupancy/beforeOccupancyCutVsCent", "After nearby-collision rejection, before occupancy cut;centrality (%);track occupancy", kTH2D, {centAxis, occupancyVsCentAxis});
    histos.add("eventQA/occupancy/after", "After all event selections including centrality;track occupancy;events", kTH1D, {occupancyAxis});
    histos.add("eventQA/occupancy/afterVsCent", "After all event selections including centrality;centrality (%);track occupancy", kTH2D, {centAxis, occupancyVsCentAxis});

    const std::array<const char*, 11> eventSelectionLabels = {
      "Input after collision filter", "Vertex z", "sel8", "Good ITS layers", "No same-bunch pileup", "Good FT0-PV z", "No ITS ROF border", "No TF border", "No nearby collision", "Occupancy interval", "Centrality / selected"};
    auto eventCounter = histos.get<TH1>(HIST("QA/eventCounter"));
    eventCounter->SetTitle("Sequential event selection;passed or disabled cut;collisions");
    for (unsigned int i = 0; i < eventSelectionLabels.size(); ++i) {
      eventCounter->GetXaxis()->SetBinLabel(i + 1, eventSelectionLabels[i]);
    }

    // Track QA
    histos.add("trackQA/before/nITSClusters", "before;N ITS clusters;Counts", kTH1D, {itsClsAxis});
    histos.add("trackQA/after/nITSClusters", "after;N ITS clusters;Counts", kTH1D, {itsClsAxis});
    histos.add("trackQA/before/chi2prITScls", "before;#chi^{2}_{ITS}/N_{cls};Counts", kTH1D, {chi2ITSAxis});
    histos.add("trackQA/after/chi2prITScls", "after;#chi^{2}_{ITS}/N_{cls};Counts", kTH1D, {chi2ITSAxis});
    histos.add("trackQA/before/nTPCClusters", "before;N TPC clusters;Counts", kTH1D, {tpcClsAxis});
    histos.add("trackQA/after/nTPCClusters", "after;N TPC clusters;Counts", kTH1D, {tpcClsAxis});
    histos.add("trackQA/before/chi2prTPCcls", "before;#chi^{2}_{TPC}/N_{cls};Counts", kTH1D, {chi2TPCAxis});
    histos.add("trackQA/after/chi2prTPCcls", "after;#chi^{2}_{TPC}/N_{cls};Counts", kTH1D, {chi2TPCAxis});
    histos.add("trackQA/before/nTPCCrossedRows", "before;N TPC crossed rows;Counts", kTH1D, {tpcClsAxis});
    histos.add("trackQA/after/nTPCCrossedRows", "after;N TPC crossed rows;Counts", kTH1D, {tpcClsAxis});
    histos.add("trackQA/before/eta", "before;#eta;Counts", kTH1D, {etaAxis});
    histos.add("trackQA/after/eta", "after;#eta;Counts", kTH1D, {etaAxis});
    histos.add("trackQA/before/pt", "before;#it{p}_{T} (GeV/#it{c});Counts", kTH1D, {ptFullAxis});
    histos.add("trackQA/after/pt", "after;#it{p}_{T} (GeV/#it{c});Counts", kTH1D, {ptFullAxis});
    histos.add("trackQA/before/dcaXY", "before;DCA_{xy} (cm);Counts", kTH1D, {dcaXYAxis});
    histos.add("trackQA/after/dcaXY", "after;DCA_{xy} (cm);Counts", kTH1D, {dcaXYAxis});
    histos.add("trackQA/before/dcaZ", "before;DCA_{z} (cm);Counts", kTH1D, {dcaZAxis});
    histos.add("trackQA/after/dcaZ", "after;DCA_{z} (cm);Counts", kTH1D, {dcaZAxis});
    histos.add("trackQA/before/phi", "before;#varphi;Counts", kTH1D, {phiAxisFull});
    histos.add("trackQA/after/phi", "after;#varphi;Counts", kTH1D, {phiAxisFull});
    histos.add("trackQA/before/charge", "before;charge;Counts", kTH1D, {chargeAxis});
    histos.add("trackQA/after/charge", "after;charge;Counts", kTH1D, {chargeAxis});
    histos.add("trackQA/before/ptVsEta", "before;#eta;#it{p}_{T}", kTH2D, {etaAxis, ptFullAxis});
    histos.add("trackQA/after/ptVsEta", "after;#eta;#it{p}_{T}", kTH2D, {etaAxis, ptFullAxis});
    histos.add("trackQA/before/dcaXYvsPt", "before;DCA_{xy};#it{p}_{T}", kTH2D, {dcaXYAxis, ptFullAxis});
    histos.add("trackQA/after/dcaXYvsPt", "after;DCA_{xy};#it{p}_{T}", kTH2D, {dcaXYAxis, ptFullAxis});
    histos.add("trackQA/before/dcaZvsPt", "before;DCA_{z};#it{p}_{T}", kTH2D, {dcaZAxis, ptFullAxis});
    histos.add("trackQA/after/dcaZvsPt", "after;DCA_{z};#it{p}_{T}", kTH2D, {dcaZAxis, ptFullAxis});

    // Chosen centrality binning, default 0-10 ... 70-80.

    histos.add("SIcentClass/pNF_cent_etaGap", "<nF>;centrality class (%);#Delta#eta", kTProfile2D, {centClassAxis, etaGapAxis});
    histos.add("SIcentClass/pNB_cent_etaGap", "<nB>;centrality class (%);#Delta#eta", kTProfile2D, {centClassAxis, etaGapAxis});
    histos.add("SIcentClass/pNF2_cent_etaGap", "<nF^{2}>;centrality class (%);#Delta#eta", kTProfile2D, {centClassAxis, etaGapAxis});
    histos.add("SIcentClass/pNB2_cent_etaGap", "<nB^{2}>;centrality class (%);#Delta#eta", kTProfile2D, {centClassAxis, etaGapAxis});
    histos.add("SIcentClass/pNFNB_cent_etaGap", "<nF nB>;centrality class (%);#Delta#eta", kTProfile2D, {centClassAxis, etaGapAxis});

    histos.add("SIcentWindow/pNF_cent_etaGap", "<nF> in centrality window;centrality class (%);#Delta#eta", kTProfile2D, {centClassAxis, etaGapAxis});
    histos.add("SIcentWindow/pNB_cent_etaGap", "<nB> in centrality window;centrality class (%);#Delta#eta", kTProfile2D, {centClassAxis, etaGapAxis});
    histos.add("SIcentWindow/pNF2_cent_etaGap", "<nF^{2}> in centrality window;centrality class (%);#Delta#eta", kTProfile2D, {centClassAxis, etaGapAxis});
    histos.add("SIcentWindow/pNB2_cent_etaGap", "<nB^{2}> in centrality window;centrality class (%);#Delta#eta", kTProfile2D, {centClassAxis, etaGapAxis});
    histos.add("SIcentWindow/pNFNB_cent_etaGap", "<nF nB> in centrality window;centrality class (%);#Delta#eta", kTProfile2D, {centClassAxis, etaGapAxis});

    // Differential eta-gap x pT x phi moments.
    histos.add("SI3D/pNF_etaGapPtPhi", "<nF>;#Delta#eta;#it{p}_{T};#varphi", kTProfile3D, {etaGapAxis, ptAxis, phiAxis});
    histos.add("SI3D/pNB_etaGapPtPhi", "<nB>;#Delta#eta;#it{p}_{T};#varphi", kTProfile3D, {etaGapAxis, ptAxis, phiAxis});
    histos.add("SI3D/pNF2_etaGapPtPhi", "<nF^{2}>;#Delta#eta;#it{p}_{T};#varphi", kTProfile3D, {etaGapAxis, ptAxis, phiAxis});
    histos.add("SI3D/pNB2_etaGapPtPhi", "<nB^{2}>;#Delta#eta;#it{p}_{T};#varphi", kTProfile3D, {etaGapAxis, ptAxis, phiAxis});
    histos.add("SI3D/pNFNB_etaGapPtPhi", "<nF nB>;#Delta#eta;#it{p}_{T};#varphi", kTProfile3D, {etaGapAxis, ptAxis, phiAxis});
    histos.add("SI3D/pF20_etaGapPtPhi", "<nF(nF-1)>;#Delta#eta;#it{p}_{T};#varphi", kTProfile3D, {etaGapAxis, ptAxis, phiAxis});
    histos.add("SI3D/pB20_etaGapPtPhi", "<nB(nB-1)>;#Delta#eta;#it{p}_{T};#varphi", kTProfile3D, {etaGapAxis, ptAxis, phiAxis});
    histos.add("SI3D/pF11_etaGapPtPhi", "<nF nB>;#Delta#eta;#it{p}_{T};#varphi", kTProfile3D, {etaGapAxis, ptAxis, phiAxis});

    // Subsample centrality-window output.
    histos.add("SubsampleCentWindow/pNF_sub_cent_etaGap", "<nF> subsamples;subsample;centrality class (%);#Delta#eta", kTProfile3D, {subsampleAxis, centClassAxis, etaGapAxis});
    histos.add("SubsampleCentWindow/pNB_sub_cent_etaGap", "<nB> subsamples;subsample;centrality class (%);#Delta#eta", kTProfile3D, {subsampleAxis, centClassAxis, etaGapAxis});
    histos.add("SubsampleCentWindow/pNF2_sub_cent_etaGap", "<nF^{2}> subsamples;subsample;centrality class (%);#Delta#eta", kTProfile3D, {subsampleAxis, centClassAxis, etaGapAxis});
    histos.add("SubsampleCentWindow/pNB2_sub_cent_etaGap", "<nB^{2}> subsamples;subsample;centrality class (%);#Delta#eta", kTProfile3D, {subsampleAxis, centClassAxis, etaGapAxis});
    histos.add("SubsampleCentWindow/pNFNB_sub_cent_etaGap", "<nF nB> subsamples;subsample;centrality class (%);#Delta#eta", kTProfile3D, {subsampleAxis, centClassAxis, etaGapAxis});

    // THnSparse sums per subsample for differential studies.
    histos.add("SubsampleSI3D/EventCount_sub_etaGapPtPhi", "event count;subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});
    histos.add("SubsampleSI3D/SumNF_sub_etaGapPtPhi", "sum nF;subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});
    histos.add("SubsampleSI3D/SumNB_sub_etaGapPtPhi", "sum nB;subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});
    histos.add("SubsampleSI3D/SumNF2_sub_etaGapPtPhi", "sum nF^{2};subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});
    histos.add("SubsampleSI3D/SumNB2_sub_etaGapPtPhi", "sum nB^{2};subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});
    histos.add("SubsampleSI3D/SumNFNB_sub_etaGapPtPhi", "sum nF nB;subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});
    histos.add("SubsampleSI3D/SumF20_sub_etaGapPtPhi", "sum nF(nF-1);subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});
    histos.add("SubsampleSI3D/SumB20_sub_etaGapPtPhi", "sum nB(nB-1);subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});
    histos.add("SubsampleSI3D/SumF11_sub_etaGapPtPhi", "sum nF nB;subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});

    histos.add("Reco/SI3D/pNF_etaGapPtPhi", "<nF> reco;#Delta#eta;#it{p}_{T};#varphi", kTProfile3D, {etaGapAxis, ptAxis, phiAxis});
    histos.add("Reco/SI3D/pNB_etaGapPtPhi", "<nB> reco;#Delta#eta;#it{p}_{T};#varphi", kTProfile3D, {etaGapAxis, ptAxis, phiAxis});
    histos.add("Reco/SI3D/pNF2_etaGapPtPhi", "<nF^{2}> reco;#Delta#eta;#it{p}_{T};#varphi", kTProfile3D, {etaGapAxis, ptAxis, phiAxis});
    histos.add("Reco/SI3D/pNB2_etaGapPtPhi", "<nB^{2}> reco;#Delta#eta;#it{p}_{T};#varphi", kTProfile3D, {etaGapAxis, ptAxis, phiAxis});
    histos.add("Reco/SI3D/pNFNB_etaGapPtPhi", "<nF nB> reco;#Delta#eta;#it{p}_{T};#varphi", kTProfile3D, {etaGapAxis, ptAxis, phiAxis});

    histos.add("Prim/SI3D/pNF_etaGapPtPhi", "<nF> primary reco;#Delta#eta;#it{p}_{T};#varphi", kTProfile3D, {etaGapAxis, ptAxis, phiAxis});
    histos.add("Prim/SI3D/pNB_etaGapPtPhi", "<nB> primary reco;#Delta#eta;#it{p}_{T};#varphi", kTProfile3D, {etaGapAxis, ptAxis, phiAxis});
    histos.add("Prim/SI3D/pNF2_etaGapPtPhi", "<nF^{2}> primary reco;#Delta#eta;#it{p}_{T};#varphi", kTProfile3D, {etaGapAxis, ptAxis, phiAxis});
    histos.add("Prim/SI3D/pNB2_etaGapPtPhi", "<nB^{2}> primary reco;#Delta#eta;#it{p}_{T};#varphi", kTProfile3D, {etaGapAxis, ptAxis, phiAxis});
    histos.add("Prim/SI3D/pNFNB_etaGapPtPhi", "<nF nB> primary reco;#Delta#eta;#it{p}_{T};#varphi", kTProfile3D, {etaGapAxis, ptAxis, phiAxis});

    histos.add("Gen/SI3D/pNF_etaGapPtPhi", "<nF> generated;#Delta#eta;#it{p}_{T};#varphi", kTProfile3D, {etaGapAxis, ptAxis, phiAxis});
    histos.add("Gen/SI3D/pNB_etaGapPtPhi", "<nB> generated;#Delta#eta;#it{p}_{T};#varphi", kTProfile3D, {etaGapAxis, ptAxis, phiAxis});
    histos.add("Gen/SI3D/pNF2_etaGapPtPhi", "<nF^{2}> generated;#Delta#eta;#it{p}_{T};#varphi", kTProfile3D, {etaGapAxis, ptAxis, phiAxis});
    histos.add("Gen/SI3D/pNB2_etaGapPtPhi", "<nB^{2}> generated;#Delta#eta;#it{p}_{T};#varphi", kTProfile3D, {etaGapAxis, ptAxis, phiAxis});
    histos.add("Gen/SI3D/pNFNB_etaGapPtPhi", "<nF nB> generated;#Delta#eta;#it{p}_{T};#varphi", kTProfile3D, {etaGapAxis, ptAxis, phiAxis});

    // Sample-separated MC subsample sums for SI3D.
    // Axes: subsample x eta gap x pT x phi. EventCount is the denominator.
    histos.add("SubsampleRecoSI3D/EventCount_sub_etaGapPtPhi", "event count reco;subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});
    histos.add("SubsampleRecoSI3D/SumNF_sub_etaGapPtPhi", "sum nF reco;subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});
    histos.add("SubsampleRecoSI3D/SumNB_sub_etaGapPtPhi", "sum nB reco;subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});
    histos.add("SubsampleRecoSI3D/SumNF2_sub_etaGapPtPhi", "sum nF^{2} reco;subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});
    histos.add("SubsampleRecoSI3D/SumNB2_sub_etaGapPtPhi", "sum nB^{2} reco;subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});
    histos.add("SubsampleRecoSI3D/SumNFNB_sub_etaGapPtPhi", "sum nF nB reco;subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});

    histos.add("SubsamplePrimSI3D/EventCount_sub_etaGapPtPhi", "event count primary reco;subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});
    histos.add("SubsamplePrimSI3D/SumNF_sub_etaGapPtPhi", "sum nF primary reco;subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});
    histos.add("SubsamplePrimSI3D/SumNB_sub_etaGapPtPhi", "sum nB primary reco;subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});
    histos.add("SubsamplePrimSI3D/SumNF2_sub_etaGapPtPhi", "sum nF^{2} primary reco;subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});
    histos.add("SubsamplePrimSI3D/SumNB2_sub_etaGapPtPhi", "sum nB^{2} primary reco;subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});
    histos.add("SubsamplePrimSI3D/SumNFNB_sub_etaGapPtPhi", "sum nF nB primary reco;subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});

    histos.add("SubsampleGenSI3D/EventCount_sub_etaGapPtPhi", "event count generated;subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});
    histos.add("SubsampleGenSI3D/SumNF_sub_etaGapPtPhi", "sum nF generated;subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});
    histos.add("SubsampleGenSI3D/SumNB_sub_etaGapPtPhi", "sum nB generated;subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});
    histos.add("SubsampleGenSI3D/SumNF2_sub_etaGapPtPhi", "sum nF^{2} generated;subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});
    histos.add("SubsampleGenSI3D/SumNB2_sub_etaGapPtPhi", "sum nB^{2} generated;subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});
    histos.add("SubsampleGenSI3D/SumNFNB_sub_etaGapPtPhi", "sum nF nB generated;subsample;#Delta#eta;#it{p}_{T};#varphi", kTHnSparseD, {subsampleAxis, etaGapAxis, ptAxis, phiAxis});
  }

  // ------------------------------------------------------------------
  // Helpers
  // ------------------------------------------------------------------
  template <typename C>
  float getCentrality(C const& coll) const
  {
    switch (cfgCentralityChoice.value) {
      case 0:
        return coll.centFT0C();
      case 1:
        return coll.centFT0M();
      default:
        return -999.f;
    }
  }

  int getCentClass(float cent) const
  {
    for (int ic = 0; ic < NCentClasses; ++ic) {
      if (cent >= centEdges[ic] && cent < centEdges[ic + 1]) {
        return ic;
      }
    }
    return -1;
  }

  double getCentClassCenter(int ic) const
  {
    if (ic < 0 || ic >= NCentClasses) {
      return -999.0;
    }
    return 0.5 * (centEdges[ic] + centEdges[ic + 1]);
  }

  int getCentWindowClass(float cent) const
  {
    const double halfWidth = 0.5 * static_cast<double>(cfgCentWindowWidth.value);
    if (halfWidth <= 0.0) {
      return -1;
    }

    for (int ic = 0; ic < NCentClasses; ++ic) {
      const double center = getCentClassCenter(ic);
      if (cent >= center - halfWidth && cent < center + halfWidth) {
        return ic;
      }
    }
    return -1;
  }

  template <typename C>
  int getSubsampleIndex(C const& coll) const
  {
    const int nSubsamples = std::min(std::max(cfgNSubsamples.value, 0), NSubsamplesMax);
    if (nSubsamples <= 0) {
      return -1;
    }
    return static_cast<int>(coll.globalIndex() % nSubsamples);
  }

  int getPtBin(double pt) const
  {
    for (int ipt = 0; ipt < NPtBins; ++ipt) {
      if (pt >= ptEdges[ipt] && pt < ptEdges[ipt + 1]) {
        return ipt;
      }
    }
    return -1;
  }

  int getPhiBin(double phi) const
  {
    while (phi < 0.0) {
      phi += TwoPi;
    }
    while (phi >= TwoPi) {
      phi -= TwoPi;
    }

    const int iphi = static_cast<int>(phi / TwoPi * NPhiBins);
    if (iphi < 0 || iphi >= NPhiBins) {
      return -1;
    }
    return iphi;
  }

  bool isChargedMCParticle(int pdgCode) const
  {
    auto* particle = pdg->GetParticle(pdgCode);
    if (!particle) {
      return false;
    }

    constexpr double MinAbsCharge = 0.0;
    return std::abs(particle->Charge()) > MinAbsCharge;
  }

  template <typename C>
  bool selCollision(C const& coll, float& cent)
  {
    histos.fill(HIST("QA/eventCounter"), 1);

    if (std::abs(coll.posZ()) > cfgCutVertex.value) {
      return false;
    }
    histos.fill(HIST("QA/eventCounter"), 2);

    if (cfgSel8Trig.value && !coll.sel8()) {
      return false;
    }
    histos.fill(HIST("QA/eventCounter"), 3);

    if (cfgUseGoodITSLayerAllCut.value && !coll.selection_bit(aod::evsel::kIsGoodITSLayersAll)) {
      return false;
    }
    histos.fill(HIST("QA/eventCounter"), 4);

    if (cfgEvSelkNoSameBunchPileup.value && !coll.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    histos.fill(HIST("QA/eventCounter"), 5);

    if (cfgEvSelUseGoodZvtxFT0vsPV.value && !coll.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    histos.fill(HIST("QA/eventCounter"), 6);

    // Occupancy is supplied by EvSels, not the current event's nF or nB.
    const int occupancy = coll.trackOccupancyInTimeRange();
    cent = getCentrality(coll);
    const bool hasValidCentrality = std::isfinite(cent) && cent >= 0.f && cent < 100.f;
    histos.fill(HIST("eventQA/occupancy/beforeExtraCuts"), occupancy);
    if (hasValidCentrality) {
      histos.fill(HIST("eventQA/occupancy/beforeExtraCutsVsCent"), cent, occupancy);
    }

    if (cfgUseNoITSROFrameBorder.value && !coll.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    histos.fill(HIST("QA/eventCounter"), 7);

    if (cfgUseNoTimeFrameBorder.value && !coll.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    histos.fill(HIST("QA/eventCounter"), 8);

    if (cfgUseNoCollInTimeRangeStandard.value && !coll.selection_bit(aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    histos.fill(HIST("QA/eventCounter"), 9);

    histos.fill(HIST("eventQA/occupancy/beforeOccupancyCut"), occupancy);
    if (hasValidCentrality) {
      histos.fill(HIST("eventQA/occupancy/beforeOccupancyCutVsCent"), cent, occupancy);
    }
    if (cfgUseOccupancyCut.value && (occupancy < cfgMinOccupancy.value || occupancy > cfgMaxOccupancy.value)) {
      return false;
    }
    histos.fill(HIST("QA/eventCounter"), 10);

    constexpr float MinCentrality = 0.f;
    constexpr float MaxCentrality = 100.f;

    if (!std::isfinite(cent) || cent < MinCentrality || cent >= MaxCentrality) {
      return false;
    }
    if (cent < cfgCentMin.value || cent >= cfgCentMax.value) {
      return false;
    }
    histos.fill(HIST("QA/eventCounter"), 11);
    histos.fill(HIST("eventQA/occupancy/after"), occupancy);
    histos.fill(HIST("eventQA/occupancy/afterVsCent"), cent, occupancy);

    return true;
  }

  template <typename T>
  void fillTrackQABefore(T const& track)
  {
    histos.fill(HIST("trackQA/before/nITSClusters"), track.itsNCls());
    histos.fill(HIST("trackQA/before/chi2prITScls"), track.itsChi2NCl());
    histos.fill(HIST("trackQA/before/nTPCClusters"), track.tpcNClsFound());
    histos.fill(HIST("trackQA/before/chi2prTPCcls"), track.tpcChi2NCl());
    histos.fill(HIST("trackQA/before/nTPCCrossedRows"), track.tpcNClsCrossedRows());
    histos.fill(HIST("trackQA/before/eta"), track.eta());
    histos.fill(HIST("trackQA/before/pt"), track.pt());
    histos.fill(HIST("trackQA/before/dcaXY"), track.dcaXY());
    histos.fill(HIST("trackQA/before/dcaZ"), track.dcaZ());
    histos.fill(HIST("trackQA/before/phi"), track.phi());
    histos.fill(HIST("trackQA/before/charge"), track.sign());
    histos.fill(HIST("trackQA/before/ptVsEta"), track.eta(), track.pt());
    histos.fill(HIST("trackQA/before/dcaXYvsPt"), track.dcaXY(), track.pt());
    histos.fill(HIST("trackQA/before/dcaZvsPt"), track.dcaZ(), track.pt());
  }

  template <typename T>
  void fillTrackQAAfter(T const& track)
  {
    histos.fill(HIST("trackQA/after/nITSClusters"), track.itsNCls());
    histos.fill(HIST("trackQA/after/chi2prITScls"), track.itsChi2NCl());
    histos.fill(HIST("trackQA/after/nTPCClusters"), track.tpcNClsFound());
    histos.fill(HIST("trackQA/after/chi2prTPCcls"), track.tpcChi2NCl());
    histos.fill(HIST("trackQA/after/nTPCCrossedRows"), track.tpcNClsCrossedRows());
    histos.fill(HIST("trackQA/after/dcaXYvsPt"), track.dcaXY(), track.pt());
    histos.fill(HIST("trackQA/after/dcaZvsPt"), track.dcaZ(), track.pt());
    histos.fill(HIST("trackQA/after/eta"), track.eta());
    histos.fill(HIST("trackQA/after/pt"), track.pt());
    histos.fill(HIST("trackQA/after/dcaXY"), track.dcaXY());
    histos.fill(HIST("trackQA/after/dcaZ"), track.dcaZ());
    histos.fill(HIST("trackQA/after/phi"), track.phi());
    histos.fill(HIST("trackQA/after/charge"), track.sign());
    histos.fill(HIST("trackQA/after/ptVsEta"), track.eta(), track.pt());
  }

  template <typename T>
  bool selTrack(T const& track)
  {
    histos.fill(HIST("QA/trackCounter"), 1);

    histos.fill(HIST("QA/hGlobalTrackFlags"), 1);
    const bool isGlobal = track.isGlobalTrack();
    const bool isGlobalSDD = track.isGlobalTrackSDD();
    if (isGlobal) {
      histos.fill(HIST("QA/hGlobalTrackFlags"), 2);
    }
    if (isGlobalSDD) {
      histos.fill(HIST("QA/hGlobalTrackFlags"), 3);
    }
    if (isGlobal || isGlobalSDD) {
      histos.fill(HIST("QA/hGlobalTrackFlags"), 4);
    } else {
      histos.fill(HIST("QA/hGlobalTrackFlags"), 5);
    }

    if (cfgRequireGlobalTrack.value && !(isGlobal || isGlobalSDD)) {
      return false;
    }
    histos.fill(HIST("QA/trackCounter"), 2);

    if (track.tpcNClsFound() < cfgTPCcluster.value) {
      return false;
    }
    histos.fill(HIST("QA/trackCounter"), 3);

    if (track.tpcNClsCrossedRows() < cfgTPCnCrossedRows.value) {
      return false;
    }
    histos.fill(HIST("QA/trackCounter"), 4);

    if (track.itsNCls() < cfgITScluster.value) {
      return false;
    }
    histos.fill(HIST("QA/trackCounter"), 5);

    if (track.tpcChi2NCl() >= cfgCutTpcChi2NCl.value) {
      return false;
    }
    histos.fill(HIST("QA/trackCounter"), 6);

    if (track.itsChi2NCl() >= cfgCutItsChi2NCl.value) {
      return false;
    }
    histos.fill(HIST("QA/trackCounter"), 7);

    if (track.pt() < cfgCutPtLower.value || track.pt() >= cfgCutPtUpper.value) {
      return false;
    }
    histos.fill(HIST("QA/trackCounter"), 8);

    if (std::abs(track.eta()) >= cfgCutEta.value) {
      return false;
    }
    histos.fill(HIST("QA/trackCounter"), 9);

    // No DCAxy cut here. It is kept only for QA.
    histos.fill(HIST("QA/trackCounter"), 10);

    if (cfgUseDCAzCut.value && std::abs(track.dcaZ()) > cfgCutTrackDcaZ.value) {
      return false;
    }
    histos.fill(HIST("QA/trackCounter"), 11);

    if (track.sign() == 0) {
      return false;
    }
    histos.fill(HIST("QA/trackCounter"), 12);

    return true;
  }

  template <typename P>
  bool selMCParticle(P const& particle)
  {
    if (cfgUseOnlyPhysicalPrimary.value && !particle.isPhysicalPrimary()) {
      return false;
    }
    if (cfgUseOnlyChargedMC.value && !isChargedMCParticle(particle.pdgCode())) {
      return false;
    }
    if (particle.pt() < cfgCutPtLower.value || particle.pt() >= cfgCutPtUpper.value) {
      return false;
    }
    if (std::abs(particle.eta()) >= cfgCutEta.value) {
      return false;
    }
    return true;
  }

  void addToFBCounters(double eta, double pt, double phi,
                       std::array<double, NEtaGaps>& nF,
                       std::array<double, NEtaGaps>& nB,
                       EtaPtPhiArray& nFPtPhi,
                       EtaPtPhiArray& nBPtPhi)
  {
    const int ipt = getPtBin(pt);
    const int iphi = getPhiBin(phi);
    const bool fillPtPhi = (ipt >= 0 && iphi >= 0);

    for (int i = 0; i < NEtaGaps; ++i) {
      const double etaLow = etaCenter[i] - cfgEtaWidth.value / 2.0;
      const double etaHigh = etaCenter[i] + cfgEtaWidth.value / 2.0;

      if (eta > etaLow && eta < etaHigh) {
        nF[i] += 1.0;
        if (fillPtPhi) {
          nFPtPhi[i][ipt][iphi] += 1.0;
        }
      }

      if (eta > -etaHigh && eta < -etaLow) {
        nB[i] += 1.0;
        if (fillPtPhi) {
          nBPtPhi[i][ipt][iphi] += 1.0;
        }
      }
    }
  }

  void fillFBQA(float cent,
                std::array<double, NEtaGaps> const& nF,
                std::array<double, NEtaGaps> const& nB)
  {
    constexpr int DefaultGapBin = 0;
    const double nf = nF[DefaultGapBin];
    const double nb = nB[DefaultGapBin];

    histos.fill(HIST("QA/fb/nFvsCent"), cent, nf);
    histos.fill(HIST("QA/fb/nBvsCent"), cent, nb);
    histos.fill(HIST("QA/fb/nFvsNB"), nf, nb);
    histos.fill(HIST("QA/fb/nFminusNB"), nf - nb);
    histos.fill(HIST("QA/fb/nFplusNB"), nf + nb);
  }

  void fillAllGapMultiplicityQA(float cent,
                                std::array<double, NEtaGaps> const& nF,
                                std::array<double, NEtaGaps> const& nB)
  {
    for (int i = 0; i < NEtaGaps; ++i) {
      const double gap = 2.0 * (etaCenter[i] - cfgEtaWidth.value / 2.0);

      // etaCenter = 0 corresponds to overlapping windows and
      // gives DeltaEta = -0.2. Keep only the seven physical,
      // non-overlapping values from 0.0 to 1.2.
      if (gap < -1.e-6 || gap > 1.200001) {
        continue;
      }

      histos.fill(HIST("QA/fbAllGaps/nFvsCentGap"), cent, gap, nF[i]);

      histos.fill(HIST("QA/fbAllGaps/nBvsCentGap"), cent, gap, nB[i]);
    }
  }

  void fillInclusiveEtaGapMoments(std::array<double, NEtaGaps> const& nF,
                                  std::array<double, NEtaGaps> const& nB)
  {
    for (int i = 0; i < NEtaGaps; ++i) {

      histos.fill(HIST("QA/hNF"), nF[i]);
      histos.fill(HIST("QA/hNB"), nB[i]);
      histos.fill(HIST("QA/hNFvsNB"), nF[i], nB[i]);
    }
  }

  void fillCentClassEtaGapMoments(int centClass,
                                  std::array<double, NEtaGaps> const& nF,
                                  std::array<double, NEtaGaps> const& nB)
  {
    if (centClass < 0 || centClass >= NCentClasses) {
      return;
    }

    const double centCenter = getCentClassCenter(centClass);
    for (int i = 0; i < NEtaGaps; ++i) {
      const double gap = 2.0 * (etaCenter[i] - cfgEtaWidth.value / 2.0);
      histos.fill(HIST("SIcentClass/pNF_cent_etaGap"), centCenter, gap, nF[i]);
      histos.fill(HIST("SIcentClass/pNB_cent_etaGap"), centCenter, gap, nB[i]);
      histos.fill(HIST("SIcentClass/pNF2_cent_etaGap"), centCenter, gap, nF[i] * nF[i]);
      histos.fill(HIST("SIcentClass/pNB2_cent_etaGap"), centCenter, gap, nB[i] * nB[i]);
      histos.fill(HIST("SIcentClass/pNFNB_cent_etaGap"), centCenter, gap, nF[i] * nB[i]);
    }
  }

  void fillCentWindowEtaGapMoments(int centWindowClass,
                                   std::array<double, NEtaGaps> const& nF,
                                   std::array<double, NEtaGaps> const& nB)
  {
    if (centWindowClass < 0 || centWindowClass >= NCentClasses) {
      return;
    }

    const double centCenter = getCentClassCenter(centWindowClass);
    for (int i = 0; i < NEtaGaps; ++i) {
      const double gap = 2.0 * (etaCenter[i] - cfgEtaWidth.value / 2.0);
      histos.fill(HIST("SIcentWindow/pNF_cent_etaGap"), centCenter, gap, nF[i]);
      histos.fill(HIST("SIcentWindow/pNB_cent_etaGap"), centCenter, gap, nB[i]);
      histos.fill(HIST("SIcentWindow/pNF2_cent_etaGap"), centCenter, gap, nF[i] * nF[i]);
      histos.fill(HIST("SIcentWindow/pNB2_cent_etaGap"), centCenter, gap, nB[i] * nB[i]);
      histos.fill(HIST("SIcentWindow/pNFNB_cent_etaGap"), centCenter, gap, nF[i] * nB[i]);
    }
  }

  void fillSubsampleCentWindowEtaGapMoments(int isub,
                                            int centWindowClass,
                                            std::array<double, NEtaGaps> const& nF,
                                            std::array<double, NEtaGaps> const& nB)
  {
    if (isub < 0 || isub >= NSubsamplesMax) {
      return;
    }
    if (centWindowClass < 0 || centWindowClass >= NCentClasses) {
      return;
    }

    const double centCenter = getCentClassCenter(centWindowClass);
    for (int i = 0; i < NEtaGaps; ++i) {
      const double gap = 2.0 * (etaCenter[i] - cfgEtaWidth.value / 2.0);
      histos.fill(HIST("SubsampleCentWindow/pNF_sub_cent_etaGap"), isub, centCenter, gap, nF[i]);
      histos.fill(HIST("SubsampleCentWindow/pNB_sub_cent_etaGap"), isub, centCenter, gap, nB[i]);
      histos.fill(HIST("SubsampleCentWindow/pNF2_sub_cent_etaGap"), isub, centCenter, gap, nF[i] * nF[i]);
      histos.fill(HIST("SubsampleCentWindow/pNB2_sub_cent_etaGap"), isub, centCenter, gap, nB[i] * nB[i]);
      histos.fill(HIST("SubsampleCentWindow/pNFNB_sub_cent_etaGap"), isub, centCenter, gap, nF[i] * nB[i]);
    }
  }

  void fillDifferentialEtaPtPhiMoments(EtaPtPhiArray const& nF,
                                       EtaPtPhiArray const& nB)
  {
    for (int igap = 0; igap < NEtaGaps; ++igap) {
      const double gap = 2.0 * (etaCenter[igap] - cfgEtaWidth.value / 2.0);
      for (int ipt = 0; ipt < NPtBins; ++ipt) {
        const double ptCenter = 0.5 * (ptEdges[ipt] + ptEdges[ipt + 1]);
        for (int iphi = 0; iphi < NPhiBins; ++iphi) {
          const double phiCenter = (iphi + 0.5) * TwoPi / NPhiBins;
          const double nf = nF[igap][ipt][iphi];
          const double nb = nB[igap][ipt][iphi];

          histos.fill(HIST("SI3D/pNF_etaGapPtPhi"), gap, ptCenter, phiCenter, nf);
          histos.fill(HIST("SI3D/pNB_etaGapPtPhi"), gap, ptCenter, phiCenter, nb);
          histos.fill(HIST("SI3D/pNF2_etaGapPtPhi"), gap, ptCenter, phiCenter, nf * nf);
          histos.fill(HIST("SI3D/pNB2_etaGapPtPhi"), gap, ptCenter, phiCenter, nb * nb);
          histos.fill(HIST("SI3D/pNFNB_etaGapPtPhi"), gap, ptCenter, phiCenter, nf * nb);
          histos.fill(HIST("SI3D/pF20_etaGapPtPhi"), gap, ptCenter, phiCenter, nf * (nf - 1.0));
          histos.fill(HIST("SI3D/pB20_etaGapPtPhi"), gap, ptCenter, phiCenter, nb * (nb - 1.0));
          histos.fill(HIST("SI3D/pF11_etaGapPtPhi"), gap, ptCenter, phiCenter, nf * nb);
        }
      }
    }
  }

  void fillRecoDifferentialEtaPtPhiMoments(EtaPtPhiArray const& nF,
                                           EtaPtPhiArray const& nB)
  {
    for (int igap = 0; igap < NEtaGaps; ++igap) {
      const double gap = 2.0 * (etaCenter[igap] - cfgEtaWidth.value / 2.0);
      for (int ipt = 0; ipt < NPtBins; ++ipt) {
        const double ptCenter = 0.5 * (ptEdges[ipt] + ptEdges[ipt + 1]);
        for (int iphi = 0; iphi < NPhiBins; ++iphi) {
          const double phiCenter = (iphi + 0.5) * TwoPi / NPhiBins;
          const double nf = nF[igap][ipt][iphi];
          const double nb = nB[igap][ipt][iphi];

          histos.fill(HIST("Reco/SI3D/pNF_etaGapPtPhi"), gap, ptCenter, phiCenter, nf);
          histos.fill(HIST("Reco/SI3D/pNB_etaGapPtPhi"), gap, ptCenter, phiCenter, nb);
          histos.fill(HIST("Reco/SI3D/pNF2_etaGapPtPhi"), gap, ptCenter, phiCenter, nf * nf);
          histos.fill(HIST("Reco/SI3D/pNB2_etaGapPtPhi"), gap, ptCenter, phiCenter, nb * nb);
          histos.fill(HIST("Reco/SI3D/pNFNB_etaGapPtPhi"), gap, ptCenter, phiCenter, nf * nb);
        }
      }
    }
  }
  void fillPrimDifferentialEtaPtPhiMoments(EtaPtPhiArray const& nF,
                                           EtaPtPhiArray const& nB)
  {
    for (int igap = 0; igap < NEtaGaps; ++igap) {
      const double gap = 2.0 * (etaCenter[igap] - cfgEtaWidth.value / 2.0);

      for (int ipt = 0; ipt < NPtBins; ++ipt) {
        const double ptCenter = 0.5 * (ptEdges[ipt] + ptEdges[ipt + 1]);

        for (int iphi = 0; iphi < NPhiBins; ++iphi) {
          const double phiCenter = (iphi + 0.5) * TwoPi / NPhiBins;

          const double nf = nF[igap][ipt][iphi];
          const double nb = nB[igap][ipt][iphi];

          histos.fill(HIST("Prim/SI3D/pNF_etaGapPtPhi"), gap, ptCenter, phiCenter, nf);
          histos.fill(HIST("Prim/SI3D/pNB_etaGapPtPhi"), gap, ptCenter, phiCenter, nb);
          histos.fill(HIST("Prim/SI3D/pNF2_etaGapPtPhi"), gap, ptCenter, phiCenter, nf * nf);
          histos.fill(HIST("Prim/SI3D/pNB2_etaGapPtPhi"), gap, ptCenter, phiCenter, nb * nb);
          histos.fill(HIST("Prim/SI3D/pNFNB_etaGapPtPhi"), gap, ptCenter, phiCenter, nf * nb);
        }
      }
    }
  }
  void fillGenDifferentialEtaPtPhiMoments(EtaPtPhiArray const& nF,
                                          EtaPtPhiArray const& nB)
  {
    for (int igap = 0; igap < NEtaGaps; ++igap) {
      const double gap = 2.0 * (etaCenter[igap] - cfgEtaWidth.value / 2.0);

      for (int ipt = 0; ipt < NPtBins; ++ipt) {
        const double ptCenter = 0.5 * (ptEdges[ipt] + ptEdges[ipt + 1]);

        for (int iphi = 0; iphi < NPhiBins; ++iphi) {
          const double phiCenter = (iphi + 0.5) * TwoPi / NPhiBins;

          const double nf = nF[igap][ipt][iphi];
          const double nb = nB[igap][ipt][iphi];

          histos.fill(HIST("Gen/SI3D/pNF_etaGapPtPhi"), gap, ptCenter, phiCenter, nf);
          histos.fill(HIST("Gen/SI3D/pNB_etaGapPtPhi"), gap, ptCenter, phiCenter, nb);
          histos.fill(HIST("Gen/SI3D/pNF2_etaGapPtPhi"), gap, ptCenter, phiCenter, nf * nf);
          histos.fill(HIST("Gen/SI3D/pNB2_etaGapPtPhi"), gap, ptCenter, phiCenter, nb * nb);
          histos.fill(HIST("Gen/SI3D/pNFNB_etaGapPtPhi"), gap, ptCenter, phiCenter, nf * nb);
        }
      }
    }
  }

  void fillRecoSubsampleDifferentialEtaPtPhiMoments(int isub,
                                                    EtaPtPhiArray const& nF,
                                                    EtaPtPhiArray const& nB)
  {
    if (isub < 0 || isub >= NSubsamplesMax) {
      return;
    }

    for (int igap = 0; igap < NEtaGaps; ++igap) {
      const double gap = 2.0 * (etaCenter[igap] - cfgEtaWidth.value / 2.0);
      for (int ipt = 0; ipt < NPtBins; ++ipt) {
        const double ptCenter = 0.5 * (ptEdges[ipt] + ptEdges[ipt + 1]);
        for (int iphi = 0; iphi < NPhiBins; ++iphi) {
          const double phiCenter = (iphi + 0.5) * TwoPi / NPhiBins;
          const double nf = nF[igap][ipt][iphi];
          const double nb = nB[igap][ipt][iphi];

          histos.fill(HIST("SubsampleRecoSI3D/EventCount_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, 1.0);
          histos.fill(HIST("SubsampleRecoSI3D/SumNF_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, nf);
          histos.fill(HIST("SubsampleRecoSI3D/SumNB_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, nb);
          histos.fill(HIST("SubsampleRecoSI3D/SumNF2_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, nf * nf);
          histos.fill(HIST("SubsampleRecoSI3D/SumNB2_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, nb * nb);
          histos.fill(HIST("SubsampleRecoSI3D/SumNFNB_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, nf * nb);
        }
      }
    }
  }

  void fillPrimSubsampleDifferentialEtaPtPhiMoments(int isub,
                                                    EtaPtPhiArray const& nF,
                                                    EtaPtPhiArray const& nB)
  {
    if (isub < 0 || isub >= NSubsamplesMax) {
      return;
    }

    for (int igap = 0; igap < NEtaGaps; ++igap) {
      const double gap = 2.0 * (etaCenter[igap] - cfgEtaWidth.value / 2.0);
      for (int ipt = 0; ipt < NPtBins; ++ipt) {
        const double ptCenter = 0.5 * (ptEdges[ipt] + ptEdges[ipt + 1]);
        for (int iphi = 0; iphi < NPhiBins; ++iphi) {
          const double phiCenter = (iphi + 0.5) * TwoPi / NPhiBins;
          const double nf = nF[igap][ipt][iphi];
          const double nb = nB[igap][ipt][iphi];

          histos.fill(HIST("SubsamplePrimSI3D/EventCount_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, 1.0);
          histos.fill(HIST("SubsamplePrimSI3D/SumNF_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, nf);
          histos.fill(HIST("SubsamplePrimSI3D/SumNB_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, nb);
          histos.fill(HIST("SubsamplePrimSI3D/SumNF2_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, nf * nf);
          histos.fill(HIST("SubsamplePrimSI3D/SumNB2_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, nb * nb);
          histos.fill(HIST("SubsamplePrimSI3D/SumNFNB_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, nf * nb);
        }
      }
    }
  }

  void fillGenSubsampleDifferentialEtaPtPhiMoments(int isub,
                                                   EtaPtPhiArray const& nF,
                                                   EtaPtPhiArray const& nB)
  {
    if (isub < 0 || isub >= NSubsamplesMax) {
      return;
    }

    for (int igap = 0; igap < NEtaGaps; ++igap) {
      const double gap = 2.0 * (etaCenter[igap] - cfgEtaWidth.value / 2.0);
      for (int ipt = 0; ipt < NPtBins; ++ipt) {
        const double ptCenter = 0.5 * (ptEdges[ipt] + ptEdges[ipt + 1]);
        for (int iphi = 0; iphi < NPhiBins; ++iphi) {
          const double phiCenter = (iphi + 0.5) * TwoPi / NPhiBins;
          const double nf = nF[igap][ipt][iphi];
          const double nb = nB[igap][ipt][iphi];

          histos.fill(HIST("SubsampleGenSI3D/EventCount_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, 1.0);
          histos.fill(HIST("SubsampleGenSI3D/SumNF_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, nf);
          histos.fill(HIST("SubsampleGenSI3D/SumNB_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, nb);
          histos.fill(HIST("SubsampleGenSI3D/SumNF2_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, nf * nf);
          histos.fill(HIST("SubsampleGenSI3D/SumNB2_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, nb * nb);
          histos.fill(HIST("SubsampleGenSI3D/SumNFNB_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, nf * nb);
        }
      }
    }
  }

  void fillSubsampleDifferentialEtaPtPhiMoments(int isub,
                                                EtaPtPhiArray const& nF,
                                                EtaPtPhiArray const& nB)
  {
    if (isub < 0 || isub >= NSubsamplesMax) {
      return;
    }

    for (int igap = 0; igap < NEtaGaps; ++igap) {
      const double gap = 2.0 * (etaCenter[igap] - cfgEtaWidth.value / 2.0);
      for (int ipt = 0; ipt < NPtBins; ++ipt) {
        const double ptCenter = 0.5 * (ptEdges[ipt] + ptEdges[ipt + 1]);
        for (int iphi = 0; iphi < NPhiBins; ++iphi) {
          const double phiCenter = (iphi + 0.5) * TwoPi / NPhiBins;
          const double nf = nF[igap][ipt][iphi];
          const double nb = nB[igap][ipt][iphi];

          histos.fill(HIST("SubsampleSI3D/EventCount_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, 1.0);
          histos.fill(HIST("SubsampleSI3D/SumNF_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, nf);
          histos.fill(HIST("SubsampleSI3D/SumNB_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, nb);
          histos.fill(HIST("SubsampleSI3D/SumNF2_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, nf * nf);
          histos.fill(HIST("SubsampleSI3D/SumNB2_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, nb * nb);
          histos.fill(HIST("SubsampleSI3D/SumNFNB_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, nf * nb);
          histos.fill(HIST("SubsampleSI3D/SumF20_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, nf * (nf - 1.0));
          histos.fill(HIST("SubsampleSI3D/SumB20_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, nb * (nb - 1.0));
          histos.fill(HIST("SubsampleSI3D/SumF11_sub_etaGapPtPhi"), isub, gap, ptCenter, phiCenter, nf * nb);
        }
      }
    }
  }

  template <typename C>
  void fillEventQABefore(C const& coll, float centBefore)
  {
    histos.fill(HIST("eventQA/before/vtxZ"), coll.posZ());
    histos.fill(HIST("eventQA/before/nPVeta1TracksVsCent"), centBefore, coll.multNTracksPVeta1());
    histos.fill(HIST("eventQA/before/nPVTracksVsCent"), centBefore, coll.multNTracksPV());
    histos.fill(HIST("eventQA/before/centFT0C"), coll.centFT0C());
    histos.fill(HIST("eventQA/before/centFT0M"), coll.centFT0M());
    histos.fill(HIST("eventQA/before/centFT0CvsFT0M"), coll.centFT0C(), coll.centFT0M());
  }

  template <typename C>
  void fillEventQAAfter(C const& coll, float cent)
  {
    histos.fill(HIST("eventQA/after/vtxZ"), coll.posZ());
    histos.fill(HIST("eventQA/after/nPVeta1TracksVsCent"), cent, coll.multNTracksPVeta1());
    histos.fill(HIST("eventQA/after/nPVTracksVsCent"), cent, coll.multNTracksPV());
    histos.fill(HIST("eventQA/after/centFT0C"), coll.centFT0C());
    histos.fill(HIST("eventQA/after/centFT0M"), coll.centFT0M());
    histos.fill(HIST("eventQA/after/centFT0CvsFT0M"), coll.centFT0C(), coll.centFT0M());
  }

  template <typename C>
  void fillAllFBOutputs(C const& coll,
                        float cent,
                        std::array<double, NEtaGaps> const& nF,
                        std::array<double, NEtaGaps> const& nB,
                        EtaPtPhiArray const& nFPtPhi,
                        EtaPtPhiArray const& nBPtPhi)
  {
    const int centClass = getCentClass(cent);
    const int centWindowClass = getCentWindowClass(cent);
    const int subsampleIndex = getSubsampleIndex(coll);

    histos.fill(HIST("QA/hCentrality"), cent);
    fillFBQA(cent, nF, nB);
    fillAllGapMultiplicityQA(cent, nF, nB);
    fillInclusiveEtaGapMoments(nF, nB);
    fillCentClassEtaGapMoments(centClass, nF, nB);
    fillCentWindowEtaGapMoments(centWindowClass, nF, nB);
    fillSubsampleCentWindowEtaGapMoments(subsampleIndex, centWindowClass, nF, nB);
    fillDifferentialEtaPtPhiMoments(nFPtPhi, nBPtPhi);
    fillSubsampleDifferentialEtaPtPhiMoments(subsampleIndex, nFPtPhi, nBPtPhi);
  }

  template <bool AfterSelection, typename T>
  void fillEtaVertexQA(float cent, float vertexZ, T const& track)
  {
    // Both stages use selected data events; "before" precedes selTrack().
    if (cfgFillEtaVtxQA.value) {
      if constexpr (AfterSelection) {
        histos.fill(HIST("trackQA/after/etaVsVtxZCent"), cent, vertexZ, track.eta());
      } else {
        histos.fill(HIST("trackQA/before/etaVsVtxZCent"), cent, vertexZ, track.eta());
      }
    }
    if (cfgFillEtaVtxPhiQA.value && std::isfinite(track.phi())) {
      double phi = std::fmod(static_cast<double>(track.phi()), TwoPi);
      if (phi < 0.0) {
        phi += TwoPi;
      }
      if constexpr (AfterSelection) {
        histos.fill(HIST("trackQA/after/etaVsVtxZPhiCent"), cent, vertexZ, track.eta(), phi);
      } else {
        histos.fill(HIST("trackQA/before/etaVsVtxZPhiCent"), cent, vertexZ, track.eta(), phi);
      }
    }
  }

  // ------------------------------------------------------------------
  // Data
  // ------------------------------------------------------------------
  void processData(AodCollision const& coll, AodTracks const& tracks)
  {
    float cent = -999.0f;
    const float centBefore = getCentrality(coll);
    fillEventQABefore(coll, centBefore);

    if (!selCollision(coll, cent)) {
      return;
    }
    fillEventQAAfter(coll, cent);
    if (cfgFillEtaVtxQA.value || cfgFillEtaVtxPhiQA.value) {
      histos.fill(HIST("eventQA/after/vtxZvsCentClass"), cent, coll.posZ());
      histos.fill(HIST("eventQA/after/centralityForEtaVtx"), cent);
    }

    std::array<double, NEtaGaps> nF{};
    std::array<double, NEtaGaps> nB{};
    EtaPtPhiArray nFPtPhi{};
    EtaPtPhiArray nBPtPhi{};

    histos.fill(HIST("QA/hTrackSliceSize"), static_cast<double>(tracks.size()));

    for (auto const& track : tracks) {
      fillTrackQABefore(track);
      fillEtaVertexQA<false>(cent, coll.posZ(), track);
      if (!selTrack(track)) {
        continue;
      }
      fillTrackQAAfter(track);
      fillEtaVertexQA<true>(cent, coll.posZ(), track);

      addToFBCounters(track.eta(), track.pt(), track.phi(), nF, nB, nFPtPhi, nBPtPhi);
    }

    fillAllFBOutputs(coll, cent, nF, nB, nFPtPhi, nBPtPhi);
  }
  PROCESS_SWITCH(StronglyIntensiveCorr, processData, "process real data", true);

  // ------------------------------------------------------------------
  // MC reconstructed: selected reco tracks with MC label.
  // This is still reconstructed-level kinematics/cuts.
  // ------------------------------------------------------------------
  void processMCRec(MyMCRecCollision const& coll,
                    MyMCTracks const& tracks,
                    aod::McParticles const&)
  {
    float cent = -999.0f;
    const float centBefore = getCentrality(coll);
    fillEventQABefore(coll, centBefore);

    if (!selCollision(coll, cent)) {
      return;
    }
    fillEventQAAfter(coll, cent);

    std::array<double, NEtaGaps> nF{};
    std::array<double, NEtaGaps> nB{};
    EtaPtPhiArray nFPtPhi{};
    EtaPtPhiArray nBPtPhi{};
    int nRecoAll = 0;
    int nRecoPrimary = 0;
    int nRecoSecondary = 0;
    int nRecoFake = 0;

    histos.fill(HIST("QA/hTrackSliceSize"), static_cast<double>(tracks.size()));

    for (auto const& track : tracks) {
      fillTrackQABefore(track);

      if (!selTrack(track)) {
        continue;
      }

      ++nRecoAll;

      if (!track.has_mcParticle()) {
        ++nRecoFake;
        continue;
      }

      auto mcParticle = track.template mcParticle_as<aod::McParticles>();

      if (mcParticle.isPhysicalPrimary() && isChargedMCParticle(mcParticle.pdgCode())) {
        ++nRecoPrimary;
      } else {
        ++nRecoSecondary;
      }

      fillTrackQAAfter(track);

      addToFBCounters(track.eta(), track.pt(), track.phi(), nF, nB, nFPtPhi, nBPtPhi);
    }
    histos.fill(HIST("QA/mcReco/nRecoAll"), nRecoAll);
    histos.fill(HIST("QA/mcReco/nRecoPrimary"), nRecoPrimary);
    histos.fill(HIST("QA/mcReco/nRecoSecondary"), nRecoSecondary);
    histos.fill(HIST("QA/mcReco/nRecoFake"), nRecoFake);
    const int subsampleIndex = getSubsampleIndex(coll);

    // For MC running with all switches enabled, do not fill common SI histograms here.
    // Keep Reco/Prim/Gen separated in multidimensional outputs.
    fillRecoDifferentialEtaPtPhiMoments(nFPtPhi, nBPtPhi);
    fillRecoSubsampleDifferentialEtaPtPhiMoments(subsampleIndex, nFPtPhi, nBPtPhi);
  }
  PROCESS_SWITCH(StronglyIntensiveCorr, processMCRec, "process MC reconstructed tracks", false);

  // ------------------------------------------------------------------
  // MC primary reconstructed: selected reco tracks matched to physical-primary charged MC particles.
  // ------------------------------------------------------------------
  void processMCPrimaryReco(MyMCRecCollision const& coll,
                            MyMCTracks const& tracks,
                            aod::McParticles const&)
  {
    float cent = -999.0f;
    const float centBefore = getCentrality(coll);
    fillEventQABefore(coll, centBefore);

    if (!selCollision(coll, cent)) {
      return;
    }
    fillEventQAAfter(coll, cent);

    std::array<double, NEtaGaps> nF{};
    std::array<double, NEtaGaps> nB{};
    EtaPtPhiArray nFPtPhi{};
    EtaPtPhiArray nBPtPhi{};

    histos.fill(HIST("QA/hTrackSliceSize"), static_cast<double>(tracks.size()));

    for (auto const& track : tracks) {
      fillTrackQABefore(track);
      if (!selTrack(track)) {
        continue;
      }
      if (!track.has_mcParticle()) {
        continue;
      }

      auto mcParticle = track.template mcParticle_as<aod::McParticles>();
      if (cfgUseOnlyPhysicalPrimary.value && !mcParticle.isPhysicalPrimary()) {
        continue;
      }
      if (cfgUseOnlyChargedMC.value && !isChargedMCParticle(mcParticle.pdgCode())) {
        continue;
      }

      fillTrackQAAfter(track);
      addToFBCounters(track.eta(), track.pt(), track.phi(), nF, nB, nFPtPhi, nBPtPhi);
    }

    const int subsampleIndex = getSubsampleIndex(coll);

    // For MC running with all switches enabled, do not fill common SI histograms here.
    // Keep Reco/Prim/Gen separated in multidimensional outputs.
    fillPrimDifferentialEtaPtPhiMoments(nFPtPhi, nBPtPhi);
    fillPrimSubsampleDifferentialEtaPtPhiMoments(subsampleIndex, nFPtPhi, nBPtPhi);
  }
  PROCESS_SWITCH(StronglyIntensiveCorr, processMCPrimaryReco, "process MC primary reconstructed tracks", false);

  // ------------------------------------------------------------------
  // MC generated: generator-level charged physical primaries from aod::McParticles.
  // Uses reconstructed collision only to get event selection + centrality and then maps to mcCollisionId.
  // ------------------------------------------------------------------
  void processMCGenerated(MyMCRecCollision const& coll,
                          aod::McParticles const& mcParticles)
  {
    if (!coll.has_mcCollision()) {
      return;
    }

    float cent = -999.0f;
    const float centBefore = getCentrality(coll);
    fillEventQABefore(coll, centBefore);

    if (!selCollision(coll, cent)) {
      return;
    }
    fillEventQAAfter(coll, cent);

    std::array<double, NEtaGaps> nF{};
    std::array<double, NEtaGaps> nB{};
    EtaPtPhiArray nFPtPhi{};
    EtaPtPhiArray nBPtPhi{};

    auto particlesThisCollision = mcParticles.sliceBy(perMcCollision, coll.mcCollisionId());

    for (auto const& particle : particlesThisCollision) {
      if (!selMCParticle(particle)) {
        continue;
      }
      addToFBCounters(particle.eta(), particle.pt(), particle.phi(), nF, nB, nFPtPhi, nBPtPhi);
    }

    const int subsampleIndex = getSubsampleIndex(coll);

    // For MC running with all switches enabled, do not fill common SI histograms here.
    // Keep Reco/Prim/Gen separated in multidimensional outputs.
    fillGenDifferentialEtaPtPhiMoments(nFPtPhi, nBPtPhi);
    fillGenSubsampleDifferentialEtaPtPhiMoments(subsampleIndex, nFPtPhi, nBPtPhi);
  }
  PROCESS_SWITCH(StronglyIntensiveCorr, processMCGenerated, "process MC generated particles", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<StronglyIntensiveCorr>(cfgc)};
}
