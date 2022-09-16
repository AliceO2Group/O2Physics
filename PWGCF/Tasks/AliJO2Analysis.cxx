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
/// \author Jasper Parkkila (jparkkil@cern.ch)
/// \author Dong Jo Kim (djkim@jyu.fi)
/// \since Sep 2022

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include <CCDB/BasicCCDBManager.h>
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "ReconstructionDataFormats/V0.h"

#include "AliJO2Catalyst.h"
#include "AliJFFlucAnalysis.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

class TrackIter : public JTrackIterInterface{
 public:
  TrackIter(aod::ParticleTrack::iterator _m1) : m1(_m1){};
  ~TrackIter(){};
  aod::ParticleTrack::iterator m1;
  PtEtaPhiEVector v;
  PtEtaPhiEVector & deref(){
    v.SetPt((*m1).pt());
    v.SetEta((*m1).eta());
    v.SetPhi((*m1).phi());
    return v;
  }
  void increment(){
    ++m1;
  }
  bool equals(const JTrackIterInterface *prhs){
    return m1 == static_cast<const TrackIter *>(prhs)->m1;
  }
};

class Tracks : public TracksBase{
 public:
  Tracks(aod::ParticleTrack const *_ptracks) : ptracks(_ptracks){};
  ~Tracks(){};
  aod::ParticleTrack const *ptracks;
  JTrackIter begin(){
    return JTrackIter(
      std::unique_ptr<JTrackIterInterface>(new TrackIter(ptracks->begin())));
  }
  JTrackIter end(){
    /*return JTrackIter(
      std::unique_ptr<JTrackIterInterface>(new TrackIterEnd(ptracks->end())));*/
    return JTrackIter(
      std::unique_ptr<JTrackIterInterface>(new TrackIter(ptracks->begin()+ptracks->size()))); //return iterator at end()
  }
  size_t size() const{
    return ptracks->size();
  }
};

struct JFlucAnalysis
{
 public:
  ~JFlucAnalysis()
  {
    delete pcf;
  }

  O2_DEFINE_CONFIGURABLE(etamin, double, 0.4, "Minimal eta for tracks");
  O2_DEFINE_CONFIGURABLE(etamax, double, 0.8, "Maximal eta for tracks");

  // OutputObj<TDirectory> output{TDirectory("jflucO2","jflucO2","",0)};
  OutputObj<TDirectory> output{"jflucO2"};

  void init(InitContext const& ic)
  {
    //
    pcf = new AliJFFlucAnalysis("jflucAnalysis");
    pcf->SetNumBins(sizeof(jflucCentBins) / sizeof(jflucCentBins[0]));
    pcf->AddFlags(AliJFFlucAnalysis::FLUC_EBE_WEIGHTING);

    output->cd();
    pcf->UserCreateOutputObjects();
  }

  void process(soa::Join<aod::Collisions, aod::CollisionData>::iterator const& collision, aod::ParticleTrack const& tracks)
  {
    if (tracks.size() == 0)
      return; // rejected event

    const double fVertex[3] = {collision.posX(), collision.posY(), collision.posZ()};

    Tracks tracksInt(&tracks);

    pcf->Init();
    pcf->SetInputList(&tracksInt);
    pcf->SetEventCentralityAndBin(collision.cent(), collision.cbin());
    pcf->SetEventVertex(fVertex);
    pcf->SetEtaRange(etamin, etamax);
    pcf->UserExec("");
  }
  AliJFFlucAnalysis* pcf;
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JFlucAnalysis>(cfgc)};
}
