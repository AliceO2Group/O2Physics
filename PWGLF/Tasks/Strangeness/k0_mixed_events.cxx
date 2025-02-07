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
//
/// \file k0_mixed_events.cxx
/// \brief Femto3D pair mixing task
/// \author Sofia Tomassini, Gleb Romanenko, Nicolò Jacazio
/// \since 31 May 2023

#include <vector>
#include <map>
#include <memory>
#include <string>
#include <utility>

#include <TParameter.h>
#include <TH1F.h>
#include <TLorentzVector.h>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoA.h"
#include "Framework/DataTypes.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/Expressions.h"
#include "Framework/StaticFor.h"

#include "MathUtils/Utils.h"
#include "Common/DataModel/Multiplicity.h"

#include "PWGCF/Femto3D/DataModel/singletrackselector.h"
#include "PWGCF/Femto3D/Core/femto3dPairTask.h"
#include "Common/DataModel/Centrality.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "PWGLF/Utils/inelGt.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "CCDB/BasicCCDBManager.h"
#include "DetectorsBase/Propagator.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

using FilteredCollisions = soa::Filtered<aod::SingleCollSels>;
using FilteredTracks = soa::Filtered<soa::Join<aod::SingleTrackSels, aod::SinglePIDEls, aod::SinglePIDPis, aod::SinglePIDKas, aod::SinglePIDPrs, aod::SinglePIDDes, aod::SinglePIDTrs, aod::SinglePIDHes>>;

using RecoTracks = soa::Join<aod::Tracks, aod::TracksExtra,
                             aod::TracksDCA,
                             aod::pidTPCFullEl,
                             aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                             aod::pidTPCFullDe, aod::pidTPCFullTr, aod::pidTPCFullHe,
                             aod::pidTOFFullEl, aod::pidTOFFullMu,
                             aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr,
                             aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe,
                             aod::TrackSelection>;

typedef std::shared_ptr<FilteredTracks::iterator> trkType;
typedef std::shared_ptr<RecoTracks::iterator> trkTypeData;

typedef std::shared_ptr<FilteredCollisions::iterator> colType;

using MyFemtoPair = o2::aod::singletrackselector::FemtoPair<trkType>;

class ResoPair : public MyFemtoPair
{
 public:
  ResoPair() {}
  ResoPair(trkType const& first, trkType const& second) : MyFemtoPair(first, second)
  {
    setPair(first, second);
  }
  ResoPair(trkType const& first, trkType const& second, const bool& isidentical) : MyFemtoPair(first, second, isidentical) {}
  bool isClosePair() const { return MyFemtoPair::IsClosePair(mDeltaEta, mDeltaPhi, mRadius); }
  void setEtaDiff(const float deta) { mDeltaEta = deta; }
  void setPhiStarDiff(const float dphi) { mDeltaPhi = dphi; }
  void setRadius(const float r) { mRadius = r; }
  void setPair(trkType const& first, trkType const& second)
  {
    MyFemtoPair::SetPair(first, second);
    lDecayDaughter1.SetPtEtaPhiM(first->pt(), first->eta(), first->phi(), particle_mass(GetPDG1()));
    lDecayDaughter2.SetPtEtaPhiM(second->pt(), second->eta(), second->phi(), particle_mass(GetPDG2()));
    lResonance = lDecayDaughter1 + lDecayDaughter2;
  }
  void setPair(trkTypeData const& first, trkTypeData const& second)
  {
    // MyFemtoPair::SetPair(first, second);
    lDecayDaughter1.SetPtEtaPhiM(first->pt(), first->eta(), first->phi(), particle_mass(GetPDG1()));
    lDecayDaughter2.SetPtEtaPhiM(second->pt(), second->eta(), second->phi(), particle_mass(GetPDG2()));
    lResonance = lDecayDaughter1 + lDecayDaughter2;
  }
  float getInvMass() const
  {
    // LOG(info) << "Mass = " << lResonance.M() << " 1 " << lDecayDaughter1.M() << " 2 " << lDecayDaughter2.M();
    return lResonance.M();
  }
  float getPt() const { return lResonance.Pt(); }
  float getRapidity() const { return lResonance.Rapidity(); }

 private:
  TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
  float mDeltaEta = 0.01;
  float mDeltaPhi = 0.01;
  float mRadius = 1.2;
};

struct K0MixedEvents {
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<std::pair<float, float>> multPercentileCut{"multPercentileCut", std::pair<float, float>{-100.f, 1000.f}, "[min., max.] centrality range to keep events within"};
  Configurable<std::pair<float, float>> momentumCut{"momentumCut", std::pair<float, float>{0.f, 100.f}, "[min., max.] momentum range to keep candidates within"};
  Configurable<float> dcaxyCut{"dcaxyCut", -100.f, "dcaXY range to keep candidates within"};
  Configurable<float> dcazCut{"dcazCut", -100.f, "dcaZ range to keep candidates within"};
  Configurable<float> dcaxyExclusionCut{"dcaxyExclusionCut", 100.f, "dcaXY range to discard candidates within"};
  Configurable<float> dcazExclusionCut{"dcazExclusionCut", 100.f, "dcaZ range to discard candidates within"};

  Configurable<float> _eta{"eta", 100.0, "abs eta value limit"};
  Configurable<int16_t> _tpcNClsFound{"minTpcNClsFound", 0, "minimum allowed number of TPC clasters"};
  Configurable<float> _tpcChi2NCl{"tpcChi2NCl", 100.0, "upper limit for chi2 value of a fit over TPC clasters"};
  Configurable<float> _tpcCrossedRowsOverFindableCls{"tpcCrossedRowsOverFindableCls", 0, "lower limit of TPC CrossedRows/FindableCls value"};
  Configurable<int> _tpcNClsShared{"maxTpcNClsShared", 100, "maximum allowed number of TPC shared clasters"};
  Configurable<int> _itsNCls{"minItsNCls", 0, "minimum allowed number of ITS clasters"};
  Configurable<float> _itsChi2NCl{"itsChi2NCl", 100.0, "upper limit for chi2 value of a fit over ITS clasters"};
  Configurable<float> _vertexZ{"VertexZ", 10.0, "abs vertexZ value limit"};
  Configurable<float> _maxy{"_maxy", 100.0, "maximum y of both particles in a pair"};

  Configurable<int> _sign_1{"sign_1", 1, "sign of the first particle in a pair"};
  Configurable<int> _particlePDG_1{"particlePDG_1", 2212, "PDG code of the first particle in a pair to perform PID for (only proton and deurton are supported now)"};
  Configurable<std::vector<float>> _tpcNSigma_1{"tpcNSigma_1", std::vector<float>{-3.0f, 3.0f}, "first particle PID: Nsigma range in TPC before the TOF is used"};
  Configurable<float> _PIDtrshld_1{"PIDtrshld_1", 10.0, "first particle PID: value of momentum from which the PID is done with TOF (before that only TPC is used)"};
  Configurable<std::vector<float>> _tofNSigma_1{"tofNSigma_1", std::vector<float>{-3.0f, 3.0f}, "first particle PID: Nsigma range in TOF"};

  Configurable<int> _sign_2{"sign_2", 1, "sign of the second particle in a pair"};
  Configurable<int> _particlePDG_2{"particlePDG_2", 2212, "PDG code of the second particle in a pair to perform PID for (only proton and deurton are supported now)"};
  Configurable<std::vector<float>> _tpcNSigma_2{"tpcNSigma_2", std::vector<float>{-3.0f, 3.0f}, "second particle PID: Nsigma range in TPC before the TOF is used"};
  Configurable<float> _PIDtrshld_2{"PIDtrshld_2", 10.0, "second particle PID: value of momentum from which the PID is done with TOF (before that only TPC is used)"};
  Configurable<std::vector<float>> _tofNSigma_2{"tofNSigma_2", std::vector<float>{-3.0f, 3.0f}, "second particle PID: Nsigma range in TOF"};

  Configurable<int> _particlePDGtoReject{"particlePDGtoRejectFromSecond", 0, "applied only if the particles are non-identical and only to the second particle in the pair!!!"};
  Configurable<std::vector<float>> _rejectWithinNsigmaTOF{"rejectWithinNsigmaTOF", std::vector<float>{-0.0f, 0.0f}, "TOF rejection Nsigma range for the particle specified with PDG to be rejected"};

  Configurable<float> _deta{"deta", 1, "minimum allowed defference in eta between two tracks in a pair"};
  Configurable<float> _dphi{"dphi", 1, "minimum allowed defference in phi_star between two tracks in a pair"};
  Configurable<float> _radiusTPC{"radiusTPC", 1.2, "TPC radius to calculate phi_star for"};

  Configurable<bool> doMixedEvent{"doMixedEvent", false, "Do the mixed event"};
  Configurable<int> _multbinwidth{"multbinwidth", 50, "width of multiplicity bins within which the mixing is done"};
  Configurable<int> _vertexbinwidth{"vertexbinwidth", 2, "width of vertexZ bins within which the mixing is done"};

  // Mag field
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  // Event selections
  Configurable<bool> sel8{"evSelsel8", 0, "Apply sel8 event selection"};
  Configurable<bool> isTriggerTVX{"evSelisTriggerTVX", 1, "Is Trigger TVX"};
  Configurable<bool> isNoTimeFrameBorder{"evSelisNoTimeFrameBorder", 1, "Is No Time Frame Border"};
  Configurable<bool> isNoITSROFrameBorder{"evSelisNoITSROFrameBorder", 1, "Is No ITS Readout Frame Border"};
  Configurable<bool> isVertexTOFmatched{"evSelisVertexTOFmatched", 0, "Is Vertex TOF matched"};
  Configurable<bool> isGoodZvtxFT0vsPV{"evSelisGoodZvtxFT0vsPV", 0, "isGoodZvtxFT0vsPV"};
  Configurable<bool> isInelGt0{"evSelisInelGt0", 0, "isInelGt0"};

  // Binnings
  ConfigurableAxis invMassBinning{"invMassBinning", {500, 0.4, 0.6}, "k* binning of the CF (Nbins, lowlimit, uplimit)"};
  ConfigurableAxis ptBinning{"ptBinning", {1000, 0.f, 10.f}, "pT binning (Nbins, lowlimit, uplimit)"};
  ConfigurableAxis dcaXyBinning{"dcaXyBinning", {100, -1.f, 1.f}, "dcaXY binning (Nbins, lowlimit, uplimit)"};
  ConfigurableAxis multPercentileBinning{"multPercentileBinning", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0}, "Binning in multiplicity percentile"};

  bool IsIdentical;

  std::pair<int, std::vector<float>> TPCcuts_1;
  std::pair<int, std::vector<float>> TOFcuts_1;

  std::pair<int, std::vector<float>> TPCcuts_2;
  std::pair<int, std::vector<float>> TOFcuts_2;

  std::unique_ptr<ResoPair> Pair = std::make_unique<ResoPair>();

  Filter pFilter = o2::aod::singletrackselector::p > momentumCut.value.first&& o2::aod::singletrackselector::p < momentumCut.value.second;
  Filter etaFilter = nabs(o2::aod::singletrackselector::eta) < _eta;
  Filter tpcTrkFilter = o2::aod::singletrackselector::tpcNClsFound >= _tpcNClsFound && o2::aod::singletrackselector::tpcNClsShared <= (uint8_t)_tpcNClsShared;

  // Filter itsNClsFilter = o2::aod::singletrackselector::itsNCls >= (uint8_t)_itsNCls;

  Filter vertexFilter = nabs(o2::aod::singletrackselector::posZ) < _vertexZ;
  Filter multPercentileFilter = o2::aod::singletrackselector::multPerc > multPercentileCut.value.first&& o2::aod::singletrackselector::multPerc < multPercentileCut.value.second;

  const char* pdgToSymbol(const int pdg)
  {
    switch (std::abs(pdg)) {
      case 211:
        return "#pi";
      case 321:
        return "K";
      case 2212:
        return "p";
      case 1000010020:
        return "d";
    }
    return "X";
  }

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    IsIdentical = (_sign_1 * _particlePDG_1 == _sign_2 * _particlePDG_2);
    LOG(info) << "IsIdentical=" << IsIdentical << "; sign1=" << _sign_1 << "; Pdg1=" << _particlePDG_1 << "; total1=" << _sign_1 * _particlePDG_1 << " -- Pdg2=" << _particlePDG_2 << "; sign2=" << _sign_2 << "; total2=" << _sign_2 * _particlePDG_2;

    Pair->SetIdentical(IsIdentical);
    Pair->SetPDG1(_particlePDG_1);
    Pair->SetPDG2(_particlePDG_2);
    Pair->setEtaDiff(_deta);
    Pair->setPhiStarDiff(_dphi);
    Pair->setRadius(_radiusTPC);

    TPCcuts_1 = std::make_pair(_particlePDG_1, _tpcNSigma_1);
    TOFcuts_1 = std::make_pair(_particlePDG_1, _tofNSigma_1);
    TPCcuts_2 = std::make_pair(_particlePDG_2, _tpcNSigma_2);
    TOFcuts_2 = std::make_pair(_particlePDG_2, _tofNSigma_2);

    const AxisSpec invMassAxis{invMassBinning, "Inv. mass (GeV/c^{2})"};
    const AxisSpec ptAxis{ptBinning, "#it{p}_{T} (GeV/c)"};
    const AxisSpec dcaXyAxis{dcaXyBinning, "DCA_{xy} (cm)"};
    const AxisSpec dcaZAxis{dcaXyBinning, "DCA_{z} (cm)"};
    const AxisSpec multPercentileAxis{multPercentileBinning, "Mult. Perc."};

    registry.add("hNEvents", "hNEvents", {HistType::kTH1I, {{11, 0.f, 11.f}}});
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(1, "all");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(2, "sel8");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(3, "TVX");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(4, "zvertex");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(5, "TFBorder");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(6, "ITSROFBorder");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(7, "isTOFVertexMatched");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(8, "isGoodZvtxFT0vsPV");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(9, "InelGT0");
    registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(10, "Applied selection");

    registry.add("Trks", "Trks", kTH1D, {{2, 0.5, 2.5, "Tracks"}});
    registry.add("VTXc", "VTXc", kTH1D, {{100, -20., 20., "vtx"}});
    registry.add("VTX", "VTX", kTH1D, {{100, -20., 20., "vtx"}});
    registry.add("multPerc", "multPerc", kTH1D, {multPercentileAxis});

    registry.add("SEcand", "SEcand", kTH1D, {{2, 0.5, 2.5}});
    registry.add("SE", "SE", kTH1D, {invMassAxis});
    registry.add("ME", "ME", kTH1D, {invMassAxis});
    registry.add("SEvsPt", "SEvsPt", kTH3F, {invMassAxis, ptAxis, multPercentileAxis});
    if (doMixedEvent) {
      registry.add("MEvsPt", "MEvsPt", kTH3F, {invMassAxis, ptAxis, multPercentileAxis});
    }
    registry.add("eta_first", Form("eta_%i", _particlePDG_1.value), kTH2F, {ptAxis, {100, -10., 10., "#eta"}});
    registry.add("p_first", Form("p_%i", _particlePDG_1.value), kTH1D, {ptAxis});
    registry.add("dcaXY_first", Form("dcaXY_%i", _particlePDG_1.value), kTH2F, {ptAxis, dcaXyAxis});
    registry.add("dcaZ_first", Form("dcaZ_%i", _particlePDG_1.value), kTH2F, {ptAxis, dcaZAxis});
    registry.add("nsigmaTOF_first", Form("nsigmaTOF_%i", _particlePDG_1.value), kTH2F, {ptAxis, {100, -10., 10., Form("N#sigma_{TOF}(%s))", pdgToSymbol(_particlePDG_1))}});
    registry.add("nsigmaTPC_first", Form("nsigmaTPC_%i", _particlePDG_1.value), kTH2F, {ptAxis, {100, -10., 10., Form("N#sigma_{TPC}(%s))", pdgToSymbol(_particlePDG_1))}});
    registry.add("rapidity_first", Form("rapidity_%i", _particlePDG_1.value), kTH2F, {ptAxis, {100, -10., 10., Form("y(%s)", pdgToSymbol(_particlePDG_1))}});

    if (!IsIdentical) {
      registry.add("p_second", Form("p_%i", _particlePDG_2.value), kTH1D, {ptAxis});
      registry.add("dcaXY_second", Form("dcaXY_%i", _particlePDG_2.value), kTH2F, {ptAxis, dcaXyAxis});
      registry.add("dcaZ_second", Form("dcaZ_%i", _particlePDG_1.value), kTH2F, {ptAxis, dcaZAxis});
      registry.add("nsigmaTOF_second", Form("nsigmaTOF_%i", _particlePDG_2.value), kTH2F, {ptAxis, {100, -10., 10., Form("N#sigma_{TOF}(%s))", pdgToSymbol(_particlePDG_2))}});
      registry.add("nsigmaTPC_second", Form("nsigmaTPC_%i", _particlePDG_2.value), kTH2F, {ptAxis, {100, -10., 10., Form("N#sigma_{TPC}(%s))", pdgToSymbol(_particlePDG_2))}});
      registry.add("rapidity_second", Form("rapidity_%i", _particlePDG_2.value), kTH2F, {ptAxis, {100, -10., 10., Form("y(%s)", pdgToSymbol(_particlePDG_2))}});
    }

    if (!doprocessMC) {
      return;
    }
    registry.add("MC/multPerc", "multPerc", kTH1D, {multPercentileAxis});
    registry.add("MC/multPercWMcCol", "multPercWMcCol", kTH1D, {multPercentileAxis});
    registry.add("MC/generatedInRecoEvs", "generatedInRecoEvs", kTH2D, {ptAxis, multPercentileAxis});
    registry.add("MC/SE", "SE", kTH1D, {invMassAxis});
    registry.add("MC/SEvsPt", "SEvsPt", kTH3F, {invMassAxis, ptAxis, multPercentileAxis});
    registry.addClone("MC/", "MCCent/");
    registry.add("MCCent/generatedInGenEvs", "generatedInGenEvs", kTH2D, {ptAxis, multPercentileAxis});
  }

  int mRunNumber = 0;
  float d_bz = 0.f;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  void initCCDB(aod::BCsWithTimestamps::iterator const& bc) // inspired by PWGLF/TableProducer/lambdakzerobuilder.cxx
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    d_bz = 0.f;

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();
    d_bz = 0.1 * d_bz;
  }

  template <typename Type>
  void mixTracks(Type const& tracks, const float centrality)
  { // template for identical particles from the same collision

    LOG(debug) << "Mixing tracks of the same event";
    for (uint32_t trk1 = 0; trk1 < tracks.size(); trk1++) { // nested loop for all the combinations
      for (uint32_t trk2 = trk1 + 1; trk2 < tracks.size(); trk2++) {

        Pair->setPair(tracks[trk1], tracks[trk2]);

        registry.fill(HIST("SEcand"), 1.f);
        // if (!Pair->isClosePair()) {
        //   continue;
        // }
        if (std::abs(Pair->getRapidity()) > 0.5f) {
          continue;
        }
        registry.fill(HIST("SEcand"), 2.f);
        registry.fill(HIST("SE"), Pair->getInvMass());                                // close pair rejection and fillig the SE histo
        registry.fill(HIST("SEvsPt"), Pair->getInvMass(), Pair->getPt(), centrality); // close pair rejection and fillig the SE histo
      }
    }
  }

  template <bool isSameEvent = false, typename Type>
  void mixTracks(Type const& tracks1, Type const& tracks2, const float centrality)
  {
    LOG(debug) << "Mixing tracks of two different events";
    for (auto trk1 : tracks1) {
      for (auto trk2 : tracks2) {

        Pair->setPair(trk1, trk2);

        if constexpr (isSameEvent) {
          registry.fill(HIST("SEcand"), 1.f);
        }
        // if (!Pair->isClosePair()) {
        //   continue;
        // }
        if (std::abs(Pair->getRapidity()) > 0.5f) {
          continue;
        }
        if constexpr (isSameEvent) {
          registry.fill(HIST("SEcand"), 2.f);
          registry.fill(HIST("SE"), Pair->getInvMass());
          registry.fill(HIST("SEvsPt"), Pair->getInvMass(), Pair->getPt(), centrality);
        } else {
          registry.fill(HIST("ME"), Pair->getInvMass());
          registry.fill(HIST("MEvsPt"), Pair->getInvMass(), Pair->getPt(), centrality);
        }
      }
    }
  }

  template <typename TrkType>
  bool isTrackSelected(TrkType const& track)
  {
    if (track.itsChi2NCl() > 36.f)
      return false;
    if (track.itsChi2NCl() < 0.f)
      return false;
    if (track.tpcChi2NCl() < 0.f)
      return false;
    if (track.tpcChi2NCl() > 4.f)
      return false;
    if (track.itsNCls() < _itsNCls) {
      return false;
    }
    if (track.itsChi2NCl() > _itsChi2NCl) {
      return false;
    }
    if (track.tpcChi2NCl() > _tpcChi2NCl) {
      return false;
    }
    if (track.tpcCrossedRowsOverFindableCls() < _tpcCrossedRowsOverFindableCls) {
      return false;
    }
    if (std::abs(track.dcaXY()) > dcaxyCut) {
      return false;
    }
    if (std::abs(track.dcaXY()) < dcaxyExclusionCut) {
      return false;
    }
    if (std::abs(track.dcaZ()) > dcazCut) {
      return false;
    }
    if (std::abs(track.dcaZ()) < dcazExclusionCut) {
      return false;
    }
    if (track.p() < momentumCut.value.first) {
      return false;
    }
    if (track.p() > momentumCut.value.second) {
      return false;
    }
    if (std::abs(track.eta()) >= _eta) {
      return false;
    }
    if (track.tpcNClsFound() < _tpcNClsFound) {
      return false;
    }
    if (track.tpcNClsShared() > _tpcNClsShared) {
      return false;
    }

    return true;
  }

  // Event selection
  template <typename TCollision>
  bool acceptEvent(TCollision const& collision, bool fill = true)
  {
    if (fill) {
      registry.fill(HIST("hNEvents"), 0.5);
    }
    if (sel8 && !collision.sel8()) {
      return false;
    }
    if (fill) {
      registry.fill(HIST("hNEvents"), 1.5);
    }
    if (isTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    if (fill) {
      registry.fill(HIST("hNEvents"), 2.5);
    }
    if (TMath::Abs(collision.posZ()) > _vertexZ) {
      return false;
    }
    if (fill) {
      registry.fill(HIST("hNEvents"), 3.5);
    }
    if (isNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (fill) {
      registry.fill(HIST("hNEvents"), 4.5);
    }
    if (isNoITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (fill) {
      registry.fill(HIST("hNEvents"), 5.5);
    }
    if (isVertexTOFmatched && !collision.selection_bit(aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    if (fill) {
      registry.fill(HIST("hNEvents"), 6.5);
    }
    if (isGoodZvtxFT0vsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (fill) {
      registry.fill(HIST("hNEvents"), 7.5);
    }
    if (isInelGt0 && !collision.isInelGt0()) {
      return false;
    }
    if (fill) {
      registry.fill(HIST("hNEvents"), 8.5);
    }
    if (collision.centFT0M() > multPercentileCut.value.second)
      return false;
    if (collision.centFT0M() < multPercentileCut.value.first)
      return false;
    return true;
  }

  void processDerived(FilteredTracks const& tracks, FilteredCollisions const& collisions)
  {
    LOG(debug) << "Processing " << collisions.size() << " collisions and " << tracks.size() << " tracks";
    std::map<int64_t, std::vector<trkType>> selectedtracks_1;
    std::map<int64_t, std::vector<trkType>> selectedtracks_2;
    std::map<std::pair<int, int>, std::vector<colType>> mixbins;
    if (_particlePDG_1 == 0 || _particlePDG_2 == 0) {
      LOGF(fatal, "One of passed PDG is 0!!!");
    }
    registry.fill(HIST("Trks"), 2.f, tracks.size());
    for (auto collision : collisions) {
      LOG(debug) << "Collision index " << collision.globalIndex();
      registry.fill(HIST("VTXc"), collision.posZ());
      registry.fill(HIST("multPerc"), collision.multPerc());
    }

    for (auto track : tracks) {
      LOG(debug) << "Track index " << track.singleCollSelId();
      if (!isTrackSelected(track)) {
        continue;
      }
      registry.fill(HIST("Trks"), 1);
      const auto& col = track.singleCollSel_as<FilteredCollisions>();
      if (std::abs(col.posZ()) > _vertexZ)
        continue;
      if (col.multPerc() > multPercentileCut.value.second || col.multPerc() < multPercentileCut.value.first)
        continue;
      registry.fill(HIST("VTX"), col.posZ());
      registry.fill(HIST("eta_first"), track.pt(), track.eta());
      registry.fill(HIST("rapidity_first"), track.pt(), track.rapidity(particle_mass(_particlePDG_1)));

      if ((track.sign() == _sign_1) &&
          (track.p() < _PIDtrshld_1 ? o2::aod::singletrackselector::TPCselection(track, TPCcuts_1) : o2::aod::singletrackselector::TOFselection(track, TOFcuts_1))) { // filling the map: eventID <-> selected particles1
        selectedtracks_1[track.singleCollSelId()].push_back(std::make_shared<decltype(track)>(track));

        registry.fill(HIST("p_first"), track.p());
        registry.fill(HIST("dcaXY_first"), track.pt(), track.dcaXY());
        registry.fill(HIST("dcaZ_first"), track.pt(), track.dcaZ());
        switch (_particlePDG_1) {
          case 211:
            registry.fill(HIST("nsigmaTOF_first"), track.p(), track.tofNSigmaPi());
            registry.fill(HIST("nsigmaTPC_first"), track.p(), track.tpcNSigmaPi());
            break;
          case 321:
            registry.fill(HIST("nsigmaTOF_first"), track.p(), track.tofNSigmaKa());
            registry.fill(HIST("nsigmaTPC_first"), track.p(), track.tpcNSigmaKa());
            break;
          case 2212:
            registry.fill(HIST("nsigmaTOF_first"), track.p(), track.tofNSigmaPr());
            registry.fill(HIST("nsigmaTPC_first"), track.p(), track.tpcNSigmaPr());
            break;
          case 1000010020:
            registry.fill(HIST("nsigmaTOF_first"), track.p(), track.tofNSigmaDe());
            registry.fill(HIST("nsigmaTPC_first"), track.p(), track.tpcNSigmaDe());
            break;
          default:
            LOG(fatal) << "PDG code 1: " << _particlePDG_1 << " is not supported!!!";
        }
      }

      if (IsIdentical) {
        continue;
      } else if ((track.sign() == _sign_2) &&
                 (_particlePDGtoReject != 0 || !o2::aod::singletrackselector::TOFselection(track, std::make_pair(_particlePDGtoReject, _rejectWithinNsigmaTOF))) &&
                 (track.p() < _PIDtrshld_2 ? o2::aod::singletrackselector::TPCselection(track, TPCcuts_2) : o2::aod::singletrackselector::TOFselection(track, TOFcuts_2))) { // filling the map: eventID <-> selected particles2 if (see condition above ^)
        selectedtracks_2[track.singleCollSelId()].push_back(std::make_shared<decltype(track)>(track));

        registry.fill(HIST("p_second"), track.p());
        registry.fill(HIST("dcaXY_second"), track.pt(), track.dcaXY());
        registry.fill(HIST("dcaZ_second"), track.pt(), track.dcaZ());
        switch (_particlePDG_2) {
          case 211:
            registry.fill(HIST("nsigmaTOF_second"), track.p(), track.tofNSigmaPi());
            registry.fill(HIST("nsigmaTPC_second"), track.p(), track.tpcNSigmaPi());
            break;
          case 321:
            registry.fill(HIST("nsigmaTOF_second"), track.p(), track.tofNSigmaKa());
            registry.fill(HIST("nsigmaTPC_second"), track.p(), track.tpcNSigmaKa());
            break;
          case 2212:
            registry.fill(HIST("nsigmaTOF_second"), track.p(), track.tofNSigmaPr());
            registry.fill(HIST("nsigmaTPC_second"), track.p(), track.tpcNSigmaPr());
            break;
          case 1000010020:
            registry.fill(HIST("nsigmaTOF_second"), track.p(), track.tofNSigmaDe());
            registry.fill(HIST("nsigmaTPC_second"), track.p(), track.tpcNSigmaDe());
            break;
          default:
            LOG(fatal) << "PDG code 2: " << _particlePDG_2 << " is not supported!!!";
        }
      }
    }

    for (auto collision : collisions) {
      if (selectedtracks_1.find(collision.globalIndex()) == selectedtracks_1.end()) {
        if (IsIdentical)
          continue;
        else if (selectedtracks_2.find(collision.globalIndex()) == selectedtracks_2.end())
          continue;
      }

      mixbins[std::pair<int, int>{round(collision.posZ() / _vertexbinwidth), floor(collision.mult() / _multbinwidth)}].push_back(std::make_shared<decltype(collision)>(collision));
    }

    //====================================== mixing starts here ======================================

    if (IsIdentical) { //====================================== mixing identical ======================================

      for (auto i = mixbins.begin(); i != mixbins.end(); i++) { // iterating over all vertex&mult bins

        for (uint32_t indx1 = 0; indx1 < (i->second).size(); indx1++) { // loop over all the events in each vertex&mult bin

          auto col1 = (i->second)[indx1];

          Pair->SetMagField1(col1->magField());
          Pair->SetMagField2(col1->magField());

          mixTracks(selectedtracks_1[col1->index()], col1->multPerc()); // mixing SE identical
          if (!doMixedEvent) {
            continue;
          }

          for (uint32_t indx2 = indx1 + 1; indx2 < (i->second).size(); indx2++) { // nested loop for all the combinations of collisions in a chosen mult/vertex bin

            auto col2 = (i->second)[indx2];

            Pair->SetMagField2(col2->magField());
            mixTracks(selectedtracks_1[col1->index()], selectedtracks_1[col2->index()], col1->multPerc()); // mixing ME identical
          }
        }
      }

      //====================================== end of mixing identical ======================================
    } else {
      //====================================== mixing non-identical ======================================

      for (auto i = mixbins.begin(); i != mixbins.end(); i++) { // iterating over all vertex&mult bins

        for (uint32_t indx1 = 0; indx1 < (i->second).size(); indx1++) { // loop over all the events in each vertex&mult bin

          auto col1 = (i->second)[indx1];

          Pair->SetMagField1(col1->magField());
          Pair->SetMagField2(col1->magField());

          mixTracks<true>(selectedtracks_1[col1->index()], selectedtracks_2[col1->index()], col1->multPerc()); // mixing SE non-identical
          if (!doMixedEvent) {
            continue;
          }

          for (uint32_t indx2 = indx1 + 1; indx2 < (i->second).size(); indx2++) { // nested loop for all the combinations of collisions in a chosen mult/vertex bin

            auto col2 = (i->second)[indx2];

            Pair->SetMagField2(col2->magField());
            mixTracks(selectedtracks_1[col1->index()], selectedtracks_2[col2->index()], col1->multPerc()); // mixing ME non-identical
          }
        }
      }

    } //====================================== end of mixing non-identical ======================================

    // clearing up
    for (auto i = selectedtracks_1.begin(); i != selectedtracks_1.end(); i++)
      (i->second).clear();
    selectedtracks_1.clear();

    if (!IsIdentical) {
      for (auto i = selectedtracks_2.begin(); i != selectedtracks_2.end(); i++)
        (i->second).clear();
      selectedtracks_2.clear();
    }

    for (auto i = mixbins.begin(); i != mixbins.end(); i++)
      (i->second).clear();
    mixbins.clear();
  }
  PROCESS_SWITCH(K0MixedEvents, processDerived, "process derived", true);

  using RecoCollisions = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Ms, aod::PVMults>;

  void processData(RecoTracks const& tracks, RecoCollisions const& collisions, BCsWithTimestamps const& bcs)
  {
    initCCDB(bcs.iteratorAt(0));
    LOG(debug) << "Processing " << collisions.size() << " collisions and " << tracks.size() << " tracks";
    std::map<int64_t, std::vector<trkTypeData>> selectedtracks_1;
    std::map<int64_t, std::vector<trkTypeData>> selectedtracks_2;
    std::map<std::pair<int, int>, std::vector<std::shared_ptr<RecoCollisions::iterator>>> mixbins;
    if (_particlePDG_1 == 0 || _particlePDG_2 == 0) {
      LOGF(fatal, "One of passed PDG is 0!!!");
    }

    registry.fill(HIST("Trks"), 2.f, tracks.size());
    for (auto collision : collisions) {
      if (!acceptEvent(collision))
        continue;
      LOG(debug) << "Collision index " << collision.globalIndex();
      registry.fill(HIST("VTXc"), collision.posZ());
      registry.fill(HIST("multPerc"), collision.centFT0M());
    }

    for (auto track : tracks) {
      if (!isTrackSelected(track)) {
        continue;
      }
      // if (!track.isGlobalTrackWoDCA()) {
      //   continue;
      // }
      if (track.trackType() != aod::track::Track) {
        continue;
      }
      if (track.tofChi2() >= 10.f) {
        continue;
      }
      if (!track.has_collision()) {
        continue;
      }
      registry.fill(HIST("Trks"), 1);
      const auto& col = track.collision_as<RecoCollisions>();
      if (!acceptEvent(col, false))
        continue;
      if (std::abs(col.posZ()) > _vertexZ)
        continue;
      if (col.centFT0M() > multPercentileCut.value.second || col.centFT0M() < multPercentileCut.value.first)
        continue;
      registry.fill(HIST("VTX"), col.posZ());
      registry.fill(HIST("eta_first"), track.pt(), track.eta());
      registry.fill(HIST("rapidity_first"), track.pt(), track.rapidity(particle_mass(_particlePDG_1)));

      if ((track.sign() == _sign_1) &&
          (track.p() < _PIDtrshld_1 ? o2::aod::singletrackselector::TPCselection(track, TPCcuts_1) : o2::aod::singletrackselector::TOFselection(track, TOFcuts_1))) { // filling the map: eventID <-> selected particles1
        selectedtracks_1[track.collisionId()].push_back(std::make_shared<decltype(track)>(track));

        registry.fill(HIST("p_first"), track.p());
        registry.fill(HIST("dcaXY_first"), track.pt(), track.dcaXY());
        registry.fill(HIST("dcaZ_first"), track.pt(), track.dcaZ());
        switch (_particlePDG_1) {
          case 211:
            registry.fill(HIST("nsigmaTOF_first"), track.p(), track.tofNSigmaPi());
            registry.fill(HIST("nsigmaTPC_first"), track.p(), track.tpcNSigmaPi());
            break;
          case 321:
            registry.fill(HIST("nsigmaTOF_first"), track.p(), track.tofNSigmaKa());
            registry.fill(HIST("nsigmaTPC_first"), track.p(), track.tpcNSigmaKa());
            break;
          case 2212:
            registry.fill(HIST("nsigmaTOF_first"), track.p(), track.tofNSigmaPr());
            registry.fill(HIST("nsigmaTPC_first"), track.p(), track.tpcNSigmaPr());
            break;
          case 1000010020:
            registry.fill(HIST("nsigmaTOF_first"), track.p(), track.tofNSigmaDe());
            registry.fill(HIST("nsigmaTPC_first"), track.p(), track.tpcNSigmaDe());
            break;
          default:
            LOG(fatal) << "PDG code 1: " << _particlePDG_1 << " is not supported!!!";
        }
      }

      if (IsIdentical) {
        continue;
      } else if ((track.sign() == _sign_2) &&
                 (_particlePDGtoReject != 0 || !o2::aod::singletrackselector::TOFselection(track, std::make_pair(_particlePDGtoReject, _rejectWithinNsigmaTOF))) &&
                 (track.p() < _PIDtrshld_2 ? o2::aod::singletrackselector::TPCselection(track, TPCcuts_2) : o2::aod::singletrackselector::TOFselection(track, TOFcuts_2))) { // filling the map: eventID <-> selected particles2 if (see condition above ^)
        selectedtracks_2[track.collisionId()].push_back(std::make_shared<decltype(track)>(track));

        registry.fill(HIST("p_second"), track.p());
        registry.fill(HIST("dcaXY_second"), track.pt(), track.dcaXY());
        registry.fill(HIST("dcaZ_second"), track.pt(), track.dcaZ());
        switch (_particlePDG_2) {
          case 211:
            registry.fill(HIST("nsigmaTOF_second"), track.p(), track.tofNSigmaPi());
            registry.fill(HIST("nsigmaTPC_second"), track.p(), track.tpcNSigmaPi());
            break;
          case 321:
            registry.fill(HIST("nsigmaTOF_second"), track.p(), track.tofNSigmaKa());
            registry.fill(HIST("nsigmaTPC_second"), track.p(), track.tpcNSigmaKa());
            break;
          case 2212:
            registry.fill(HIST("nsigmaTOF_second"), track.p(), track.tofNSigmaPr());
            registry.fill(HIST("nsigmaTPC_second"), track.p(), track.tpcNSigmaPr());
            break;
          case 1000010020:
            registry.fill(HIST("nsigmaTOF_second"), track.p(), track.tofNSigmaDe());
            registry.fill(HIST("nsigmaTPC_second"), track.p(), track.tpcNSigmaDe());
            break;
          default:
            LOG(fatal) << "PDG code 2: " << _particlePDG_2 << " is not supported!!!";
        }
      }
    }

    for (auto collision : collisions) {
      if (selectedtracks_1.find(collision.globalIndex()) == selectedtracks_1.end()) {
        if (IsIdentical)
          continue;
        else if (selectedtracks_2.find(collision.globalIndex()) == selectedtracks_2.end())
          continue;
      }

      mixbins[std::pair<int, int>{round(collision.posZ() / _vertexbinwidth), floor(collision.multNTracksPVeta1() / _multbinwidth)}].push_back(std::make_shared<decltype(collision)>(collision));
    }

    //====================================== mixing starts here ======================================

    if (IsIdentical) { //====================================== mixing identical ======================================

      for (auto i = mixbins.begin(); i != mixbins.end(); i++) { // iterating over all vertex&mult bins

        for (uint32_t indx1 = 0; indx1 < (i->second).size(); indx1++) { // loop over all the events in each vertex&mult bin

          auto col1 = (i->second)[indx1];

          Pair->SetMagField1(d_bz);
          Pair->SetMagField2(d_bz);

          mixTracks(selectedtracks_1[col1->index()], col1->centFT0M()); // mixing SE identical
          if (!doMixedEvent) {
            continue;
          }

          for (uint32_t indx2 = indx1 + 1; indx2 < (i->second).size(); indx2++) { // nested loop for all the combinations of collisions in a chosen mult/vertex bin

            auto col2 = (i->second)[indx2];

            Pair->SetMagField2(d_bz);
            mixTracks(selectedtracks_1[col1->index()], selectedtracks_1[col2->index()], col1->centFT0M()); // mixing ME identical
          }
        }
      }

      //====================================== end of mixing identical ======================================
    } else {
      //====================================== mixing non-identical ======================================

      for (auto i = mixbins.begin(); i != mixbins.end(); i++) { // iterating over all vertex&mult bins

        for (uint32_t indx1 = 0; indx1 < (i->second).size(); indx1++) { // loop over all the events in each vertex&mult bin

          auto col1 = (i->second)[indx1];

          Pair->SetMagField1(d_bz);
          Pair->SetMagField2(d_bz);

          mixTracks<true>(selectedtracks_1[col1->index()], selectedtracks_2[col1->index()], col1->centFT0M()); // mixing SE non-identical
          if (!doMixedEvent) {
            continue;
          }

          for (uint32_t indx2 = indx1 + 1; indx2 < (i->second).size(); indx2++) { // nested loop for all the combinations of collisions in a chosen mult/vertex bin

            auto col2 = (i->second)[indx2];

            Pair->SetMagField2(d_bz);
            mixTracks(selectedtracks_1[col1->index()], selectedtracks_2[col2->index()], col1->centFT0M()); // mixing ME non-identical
          }
        }
      }

    } //====================================== end of mixing non-identical ======================================

    // clearing up
    for (auto i = selectedtracks_1.begin(); i != selectedtracks_1.end(); i++)
      (i->second).clear();
    selectedtracks_1.clear();

    if (!IsIdentical) {
      for (auto i = selectedtracks_2.begin(); i != selectedtracks_2.end(); i++)
        (i->second).clear();
      selectedtracks_2.clear();
    }

    for (auto i = mixbins.begin(); i != mixbins.end(); i++)
      (i->second).clear();
    mixbins.clear();
  }
  PROCESS_SWITCH(K0MixedEvents, processData, "process data", false);

  using RecoMCCollisions = soa::Join<RecoCollisions, aod::McCollisionLabels>;
  using RecoMCTracks = soa::Join<RecoTracks, aod::McTrackLabels>;
  using GenMCCollisions = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;

  Service<o2::framework::O2DatabasePDG> pdgDB;
  Preslice<aod::McParticles> perMCCol = aod::mcparticle::mcCollisionId;
  Preslice<RecoMCTracks> perCollision = aod::track::collisionId;
  SliceCache cache;
  void processMC(RecoMCCollisions const& collisions,
                 RecoMCTracks const& tracks,
                 GenMCCollisions const& mcCollisions,
                 aod::McParticles const& mcParticles)
  {
    // Loop on reconstructed tracks
    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
    std::vector<std::shared_ptr<RecoMCTracks::iterator>> trkPool1;
    std::vector<std::shared_ptr<RecoMCTracks::iterator>> trkPool2;
    // Loop on reconstructed collisions
    for (const auto& col : collisions) {
      if (!col.sel8()) {
        continue;
      }
      if (std::abs(col.posZ()) > _vertexZ) {
        continue;
      }
      // Loop on tracks
      const auto& tracksInCollision = tracks.sliceByCached(aod::track::collisionId, col.globalIndex(), cache);
      for (const auto& trk : tracksInCollision) {
        if (!trk.has_mcParticle()) {
          continue;
        }
        if (!isTrackSelected(trk)) {
          continue;
        }
        // if (!trk.isGlobalTrackWoDCA()) {
        //   continue;
        // }
        if (trk.trackType() != aod::track::Track) {
          continue;
        }
        if (trk.tofChi2() >= 10.f) {
          continue;
        }
        const auto& part = trk.mcParticle();
        switch (part.pdgCode()) {
          case 211:
            trkPool1.push_back(std::make_shared<RecoMCTracks::iterator>(trk));
            break;
          case -211:
            trkPool2.push_back(std::make_shared<RecoMCTracks::iterator>(trk));
            break;
          default:
            continue;
        }
      }

      for (uint32_t trk1 = 0; trk1 < trkPool1.size(); trk1++) { // nested loop for all the combinations
        lDecayDaughter1.SetPtEtaPhiM(trkPool1[trk1]->pt(), trkPool1[trk1]->eta(), trkPool1[trk1]->phi(), particle_mass(_particlePDG_1));
        for (uint32_t trk2 = 0; trk2 < trkPool2.size(); trk2++) {
          lDecayDaughter2.SetPtEtaPhiM(trkPool2[trk2]->pt(), trkPool2[trk2]->eta(), trkPool2[trk2]->phi(), particle_mass(_particlePDG_2));
          // if (!Pair->isClosePair()) {
          //   continue;
          // }
          lResonance = lDecayDaughter1 + lDecayDaughter2;
          if (std::abs(lResonance.Rapidity()) > 0.5f) {
            continue;
          }
          registry.fill(HIST("MC/SE"), lResonance.M());                                      // close pair rejection and fillig the SE histo
          registry.fill(HIST("MC/SEvsPt"), lResonance.M(), lResonance.Pt(), col.centFT0M()); // close pair rejection and fillig the SE histo
          if (col.has_mcCollision()) {
            registry.fill(HIST("MCCent/SE"), lResonance.M());                                                                        // close pair rejection and fillig the SE histo
            registry.fill(HIST("MCCent/SEvsPt"), lResonance.M(), lResonance.Pt(), col.mcCollision_as<GenMCCollisions>().centFT0M()); // close pair rejection and fillig the SE histo
          }
        }
      }
      trkPool1.clear();
      trkPool2.clear();

      registry.fill(HIST("MC/multPerc"), col.centFT0M());
      if (!col.has_mcCollision()) {
        continue;
      }
      const auto& mcCollision = col.mcCollision_as<GenMCCollisions>();
      registry.fill(HIST("MC/multPercWMcCol"), col.centFT0M());
      registry.fill(HIST("MCCent/multPercWMcCol"), mcCollision.centFT0M());

      // Loop on particles
      const auto& particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
      for (const auto& mcParticle : particlesInCollision) {
        switch (mcParticle.pdgCode()) {
          case 310:
            break;
          default:
            continue;
        }
        if (mcParticle.pdgCode() != 310) {
          LOG(fatal) << "Fatal in PDG";
        }
        if (std::abs(mcParticle.y()) > 0.5) {
          continue;
        }
        registry.fill(HIST("MC/generatedInRecoEvs"), mcParticle.pt(), col.centFT0M());
        registry.fill(HIST("MCCent/generatedInRecoEvs"), mcParticle.pt(), mcCollision.centFT0M());
      }
    }

    // Loop on generated collisions
    for (const auto& mcCollision : mcCollisions) {
      if (std::abs(mcCollision.posZ()) > _vertexZ) {
        continue;
      }
      const auto& particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
      if (!o2::pwglf::isINELgt0mc(particlesInCollision, pdgDB)) {
        continue;
      }
      registry.fill(HIST("MCCent/multPerc"), mcCollision.centFT0M());
      for (const auto& mcParticle : particlesInCollision) {
        switch (mcParticle.pdgCode()) {
          case 310:
            break;
          default:
            continue;
        }
        if (mcParticle.pdgCode() != 310) {
          LOG(fatal) << "Fatal in PDG";
        }
        if (std::abs(mcParticle.y()) > 0.5) {
          continue;
        }
        registry.fill(HIST("MCCent/generatedInGenEvs"), mcParticle.pt(), mcCollision.centFT0M());
      }
    }
  }

  PROCESS_SWITCH(K0MixedEvents, processMC, "process mc", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<K0MixedEvents>(cfgc)}; }
