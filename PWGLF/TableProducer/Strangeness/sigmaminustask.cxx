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

/// \file   sigmaminustask.cxx
/// \brief Example of a simple task for the analysis of the Sigma-minus
/// \author Francesco Mazzaschi <francesco.mazzaschi@cern.ch>

#include "PWGLF/DataModel/LFKinkDecayTables.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU,
                             aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullKa,
                             aod::pidTOFFullPi, aod::pidTOFFullPr, aod::pidTOFFullKa>;
using CollisionsFull = soa::Join<aod::Collisions, aod::EvSel>;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSel>;

struct sigmaminustask {

  // Output Tables
  Produces<aod::SlimKinkCands> outputDataTable;
  Produces<aod::SlimKinkCandsMC> outputDataTableMC;

  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rSigmaMinus{"sigmaminus", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rFindable{"findable", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for event selection
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<float> cutNSigmaPi{"cutNSigmaPi", 4, "NSigmaTPCPion"};
  Configurable<float> cutRapMotherMC{"cutRapMotherMC", 1.0f, "Rapidity cut for mother Sigma in MC"};
  Configurable<float> cutMinQtAP{"cutMinQtAP", 0.15f, "Minimum Qt for Armenteros-Podolanski cut"};
  Configurable<float> cutMaxQtAP{"cutMaxQtAP", 0.20f, "Maximum Qt for Armenteros-Podolanski cut"};
  Configurable<float> cutPtGen{"cutPtGen", 0.5f, "Minimum pT for generated sigma particles"};

  Configurable<bool> fillOutputTree{"fillOutputTree", true, "If true, fill the output tree with Kink candidates"};

  // Configurables for findable tracks (kinkBuilder.cxx efficiency)
  Configurable<float> minPtMothKB{"minPtMothKB", 0.5f, "Minimum pT of the mother"};
  Configurable<float> maxPhiDiffKB{"maxPhiDiffKB", 100.0f, "Max phi difference between the kink daughter and the mother"};
  Configurable<float> maxZDiffKB{"maxZDiffKB", 20.0f, "Max z difference between the kink daughter and the mother"};
  Configurable<float> etaMaxKB{"etaMaxKB", 1.0f, "Max eta for both mother and daughter"};
  Configurable<float> nTPCClusMinDaugKB{"nTPCClusMinDaugKB", 80, "Min daug NTPC clusters"};
  Configurable<float> radiusCutKB{"radiusCutKB", 19.6213f, "Min reconstructed decay radius of the mother"};
  Configurable<float> maxDcaMothPvKB{"maxDcaMothPvKB", 0.1f, "Max DCA of the mother to PV"};
  Configurable<float> minDcaDaugPvKB{"minDcaDaugPvKB", 0.1f, "Min DCA of the daughter to PV"};

  Preslice<aod::KinkCands> mPerCol = aod::track::collisionId;

  // Constants
  float radToDeg = o2::constants::math::Rad2Deg;

  // Services CCDB
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> ccdbPath{"ccdbPath", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrLUT), "Type of material correction"};

  // Runtime variables
  int mRunNumber = 0;
  float mBz = 0.0f;
  o2::base::MatLayerCylSet* matLUT = nullptr;

  void init(InitContext const&)
  {
    // Initialize CCDB
    ccdb->setURL(ccdbPath);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(true);

    matLUT = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));

    // Axes
    const AxisSpec ptAxis{100, -10, 10, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec ptUnsignedAxis{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec nSigmaPiAxis{60, -30, 30, "n#sigma_{#pi}"};
    const AxisSpec nSigmaPrAxis{60, -30, 30, "n#sigma_{p}"};
    const AxisSpec sigmaMassAxis{100, 1.1, 1.4, "m (GeV/#it{c}^{2})"};
    const AxisSpec xiMassAxis{100, 1.2, 1.6, "m_{#Xi} (GeV/#it{c}^{2})"};
    const AxisSpec pdgAxis{10001, -5000, 5000, "PDG code"};
    const AxisSpec vertexZAxis{100, -15., 15., "vrtx_{Z} [cm]"};
    const AxisSpec dcaMothAxis{100, 0, 0.03, "DCA [cm]"};
    const AxisSpec dcaDaugAxis{200, 0, 20, "DCA [cm]"};
    const AxisSpec radiusAxis{100, -1, 40, "Decay radius [cm]"};
    const AxisSpec alphaAPAxis{200, -1.0, 1.0, "#alpha_{AP}"};
    const AxisSpec qtAPAxis{200, 0.0, 0.5, "q_{T,AP}"};
    const AxisSpec cosPointingAngleAxis{100, -1.0, 1.0, "Cos#theta_{PA}"};

    const AxisSpec ptResolutionAxis{100, -1.0, 1.0, "(#it{p}_{T}^{rec} - #it{p}_{T}^{gen}) / #it{p}_{T}^{gen}"};
    const AxisSpec massResolutionAxis{100, -0.5, 0.5, "(m_{rec} - m_{gen}) / m_{gen}"};
    const AxisSpec radiusResolutionAxis{100, -0.5, 0.5, "(r_{rec} - r_{gen}) / r_{gen}"};

    const AxisSpec boolAxis{2, -0.5, 1.5, "Boolean value"};
    const AxisSpec filtersAxis{12, -0.5, 11.5, "Filter index"};
    const AxisSpec fakeITSAxis{8, -1.5, 6.5, "Fake ITS cluster layer"};

    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
    // Sigma-minus reconstruction
    rSigmaMinus.add("h2MassSigmaMinusPt", "h2MassSigmaMinusPt", {HistType::kTH2F, {ptAxis, sigmaMassAxis}});
    rSigmaMinus.add("h2SigmaMassVsXiMass", "h2SigmaMassVsXiMass", {HistType::kTH2F, {xiMassAxis, sigmaMassAxis}});
    rSigmaMinus.add("h2NSigmaTPCPiPt", "h2NSigmaTPCPiPt", {HistType::kTH2F, {ptAxis, nSigmaPiAxis}});
    rSigmaMinus.add("h2DCAMothPt", "h2DCAMothPt", {HistType::kTH2F, {ptAxis, dcaMothAxis}});
    rSigmaMinus.add("h2DCADaugPt", "h2DCADaugPt", {HistType::kTH2F, {ptAxis, dcaDaugAxis}});
    rSigmaMinus.add("h2ArmenterosPreCuts", "h2ArmenterosPreCuts", {HistType::kTH2F, {alphaAPAxis, qtAPAxis}});
    rSigmaMinus.add("h2ArmenterosPostCuts", "h2ArmenterosPostCuts", {HistType::kTH2F, {alphaAPAxis, qtAPAxis}});
    rSigmaMinus.add("h2CosPointingAnglePt", "h2CosPointingAnglePt", {HistType::kTH2F, {ptAxis, cosPointingAngleAxis}});

    if (doprocessMC) {
      // Add MC histograms if needed
      rSigmaMinus.add("h2MassPtMCRec", "h2MassPtMCRec", {HistType::kTH2F, {ptAxis, sigmaMassAxis}});
      rSigmaMinus.add("h2MassPtMCGen", "h2MassPtMCGen", {HistType::kTH2F, {ptAxis, sigmaMassAxis}});

      rSigmaMinus.add("h2MassResolution", "h2MassResolution", {HistType::kTH2F, {ptAxis, massResolutionAxis}});
      rSigmaMinus.add("h2PtResolution", "h2PtResolution", {HistType::kTH2F, {ptAxis, ptResolutionAxis}});
      rSigmaMinus.add("h2RadiusResolution", "h2RadiusResolution", {HistType::kTH2F, {ptAxis, radiusResolutionAxis}});

      rSigmaMinus.add("h2NSigmaTOFPiPt", "h2NSigmaTOFPiPt", {HistType::kTH2F, {ptAxis, nSigmaPiAxis}});
      rSigmaMinus.add("h2NSigmaTOFPrPt", "h2NSigmaTOFPrPt", {HistType::kTH2F, {ptAxis, nSigmaPrAxis}});

      // BC ID comparison histograms
      rSigmaMinus.add("hMcCollIdCoherence", "McCollId (coll == daug)", {HistType::kTH1F, {boolAxis}});
      rSigmaMinus.add("h2CollId_BCId", "(McCollId coherence) vs (EvSelBC == McBC)", {HistType::kTH2F, {boolAxis, boolAxis}});
      rSigmaMinus.add("h2BCId_comp", "(McBC == EvSelBC) vs (BC == EvSelBC)", {HistType::kTH2F, {boolAxis, boolAxis}});
    }

    if (doprocessFindable) {
      std::vector<std::string> filterLabels = {"Initial", "ITS/TPC present", "ITS/TPC quality", "Moth p_{T}", "max #eta", "max #Delta#phi", "max #Delta Z", "max DCAmoth", "min DCAdaug", "min Radius", "sel8 coll", "Daug TOF"};

      // Add findable Sigma histograms
      rFindable.add("hfakeITSfindable", "hfakeITSfindable", {HistType::kTH1F, {fakeITSAxis}});
      rFindable.add("hFilterIndex", "hFilterIndex", {HistType::kTH1F, {filtersAxis}});
      rFindable.add("h2MCRadiusFilterIndex", "h2MCRadiusFilterIndex", {HistType::kTH2F, {filtersAxis, radiusAxis}});
      rFindable.add("h2RecRadiusFilterIndex", "h2RecRadiusFilterIndex", {HistType::kTH2F, {filtersAxis, radiusAxis}});

      auto hFilterIndex = rFindable.get<TH1>(HIST("hFilterIndex"));
      auto h2MCRadiusFilterIndex = rFindable.get<TH2>(HIST("h2MCRadiusFilterIndex"));
      auto h2RecRadiusFilterIndex = rFindable.get<TH2>(HIST("h2RecRadiusFilterIndex"));
      for (size_t i = 0; i < filterLabels.size(); ++i) {
        hFilterIndex->GetXaxis()->SetBinLabel(i + 1, filterLabels[i].c_str());
        h2MCRadiusFilterIndex->GetXaxis()->SetBinLabel(i + 1, filterLabels[i].c_str());
        h2RecRadiusFilterIndex->GetXaxis()->SetBinLabel(i + 1, filterLabels[i].c_str());
      }

      rFindable.add("h2MCRadiusFilter_plus_protonkink", "h2MCRadiusFilter_plus_protonkink", {HistType::kTH2F, {filtersAxis, radiusAxis}});
      rFindable.add("h2MCRadiusFilter_plus_pikink", "h2MCRadiusFilter_plus_pikink", {HistType::kTH2F, {filtersAxis, radiusAxis}});
      rFindable.add("h2MCRadiusFilter_minus_pikink", "h2MCRadiusFilter_minus_pikink", {HistType::kTH2F, {filtersAxis, radiusAxis}});

      rFindable.add("h2PtFilter_plus_protonkink", "h2PtFilter_plus_protonkink", {HistType::kTH2F, {filtersAxis, ptUnsignedAxis}});
      rFindable.add("h2PtFilter_plus_pikink", "h2PtFilter_plus_pikink", {HistType::kTH2F, {filtersAxis, ptUnsignedAxis}});
      rFindable.add("h2PtFilter_minus_pikink", "h2PtFilter_minus_pikink", {HistType::kTH2F, {filtersAxis, ptUnsignedAxis}});

      rFindable.add("h2PtDaugFilter_plus_protonkink", "h2PtDaugFilter_plus_protonkink", {HistType::kTH2F, {filtersAxis, ptUnsignedAxis}});
      rFindable.add("h2PtDaugFilter_plus_pikink", "h2PtDaugFilter_plus_pikink", {HistType::kTH2F, {filtersAxis, ptUnsignedAxis}});
      rFindable.add("h2PtDaugFilter_minus_pikink", "h2PtDaugFilter_minus_pikink", {HistType::kTH2F, {filtersAxis, ptUnsignedAxis}});

      rFindable.add("h2DCAMothPt_protonkink", "h2DCAMothPt_protonkink", {HistType::kTH2F, {ptUnsignedAxis, dcaMothAxis}});
      rFindable.add("h2DCADaugPt_protonkink", "h2DCADaugPt_protonkink", {HistType::kTH2F, {ptUnsignedAxis, dcaDaugAxis}});
      rFindable.add("h2DCAMothPt_pikink", "h2DCAMothPt_pikink", {HistType::kTH2F, {ptUnsignedAxis, dcaMothAxis}});
      rFindable.add("h2DCADaugPt_pikink", "h2DCADaugPt_pikink", {HistType::kTH2F, {ptUnsignedAxis, dcaDaugAxis}});
    }
  }

  void initCCDB(aod::BCs::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
    LOG(info) << "Initializing CCDB for run " << mRunNumber;
    o2::parameters::GRPMagField* grpmag = ccdb->getForRun<o2::parameters::GRPMagField>(grpmagPath, mRunNumber);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    mBz = grpmag->getNominalL3Field();
    if (!matLUT) {
      matLUT = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
    }
    o2::base::Propagator::Instance()->setMatLUT(matLUT);
    LOG(info) << "Task initialized for run " << mRunNumber << " with magnetic field " << mBz << " kZG";
  }

  float alphaAP(const std::array<float, 3>& momMother, const std::array<float, 3>& momKink)
  {
    std::array<float, 3> momMissing = {momMother[0] - momKink[0], momMother[1] - momKink[1], momMother[2] - momKink[2]};
    float lQlP = std::inner_product(momMother.begin(), momMother.end(), momKink.begin(), 0.f);
    float lQlN = std::inner_product(momMother.begin(), momMother.end(), momMissing.begin(), 0.f);
    return (lQlP - lQlN) / (lQlP + lQlN);
  }

  float qtAP(const std::array<float, 3>& momMother, const std::array<float, 3>& momKink)
  {
    float dp = std::inner_product(momMother.begin(), momMother.end(), momKink.begin(), 0.f);
    float p2V0 = std::inner_product(momMother.begin(), momMother.end(), momMother.begin(), 0.f);
    float p2A = std::inner_product(momKink.begin(), momKink.end(), momKink.begin(), 0.f);
    return std::sqrt(p2A - dp * dp / p2V0);
  }

  float cosPAngle(const std::array<float, 3>& momMother, const std::array<float, 3>& posMother, const std::array<float, 3>& posKink)
  {
    std::array<float, 3> vMother = {posKink[0] - posMother[0], posKink[1] - posMother[1], posKink[2] - posMother[2]};
    float pMother = std::sqrt(std::inner_product(momMother.begin(), momMother.end(), momMother.begin(), 0.f));
    float vMotherNorm = std::sqrt(std::inner_product(vMother.begin(), vMother.end(), vMother.begin(), 0.f));
    return (std::inner_product(momMother.begin(), momMother.end(), vMother.begin(), 0.f)) / (pMother * vMotherNorm);
  }

  void processData(CollisionsFull::iterator const& collision, aod::KinkCands const& KinkCands, TracksFull const&)
  {
    if (std::abs(collision.posZ()) > cutzvertex || !collision.sel8()) {
      return;
    }
    rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());

    for (const auto& kinkCand : KinkCands) {
      auto dauTrack = kinkCand.trackDaug_as<TracksFull>();

      if (std::abs(dauTrack.tpcNSigmaPi()) > cutNSigmaPi) {
        continue;
      }

      float alphaAPValue = alphaAP(std::array{kinkCand.pxMoth(), kinkCand.pyMoth(), kinkCand.pzMoth()}, std::array{kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug()});
      float qtValue = qtAP(std::array{kinkCand.pxMoth(), kinkCand.pyMoth(), kinkCand.pzMoth()}, std::array{kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug()});
      rSigmaMinus.fill(HIST("h2ArmenterosPreCuts"), alphaAPValue, qtValue);

      if (qtValue < cutMinQtAP || qtValue > cutMaxQtAP) {
        continue;
      }
      float cosPointingAngleRec = cosPAngle(std::array{kinkCand.pxMoth(), kinkCand.pyMoth(), kinkCand.pzMoth()},
                                            std::array{0.0f, 0.0f, 0.0f},
                                            std::array{kinkCand.xDecVtx(), kinkCand.yDecVtx(), kinkCand.zDecVtx()});
      rSigmaMinus.fill(HIST("h2CosPointingAnglePt"), kinkCand.mothSign() * kinkCand.ptMoth(), cosPointingAngleRec);
      rSigmaMinus.fill(HIST("h2ArmenterosPostCuts"), alphaAPValue, qtValue);
      rSigmaMinus.fill(HIST("h2MassSigmaMinusPt"), kinkCand.mothSign() * kinkCand.ptMoth(), kinkCand.mSigmaMinus());
      rSigmaMinus.fill(HIST("h2SigmaMassVsXiMass"), kinkCand.mXiMinus(), kinkCand.mSigmaMinus());
      rSigmaMinus.fill(HIST("h2NSigmaTPCPiPt"), kinkCand.mothSign() * kinkCand.ptMoth(), dauTrack.tpcNSigmaPi());
      rSigmaMinus.fill(HIST("h2DCAMothPt"), kinkCand.mothSign() * kinkCand.ptMoth(), kinkCand.dcaMothPv());
      rSigmaMinus.fill(HIST("h2DCADaugPt"), kinkCand.mothSign() * kinkCand.ptDaug(), kinkCand.dcaDaugPv());

      if (fillOutputTree) {
        outputDataTable(kinkCand.xDecVtx(), kinkCand.yDecVtx(), kinkCand.zDecVtx(),
                        kinkCand.pxMoth(), kinkCand.pyMoth(), kinkCand.pzMoth(),
                        kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug(),
                        kinkCand.dcaMothPv(), kinkCand.dcaDaugPv(), kinkCand.dcaKinkTopo(),
                        kinkCand.mothSign(),
                        dauTrack.tpcNSigmaPi(), dauTrack.tpcNSigmaPr(), dauTrack.tpcNSigmaKa(),
                        dauTrack.tofNSigmaPi(), dauTrack.tofNSigmaPr(), dauTrack.tofNSigmaKa());
      }
    }
  }
  PROCESS_SWITCH(sigmaminustask, processData, "Data processing", true);

  void processMC(CollisionsFullMC const& collisions, aod::KinkCands const& KinkCands, aod::McTrackLabels const& trackLabelsMC, aod::McParticles const& particlesMC, TracksFull const&, aod::McCollisions const&)
  {
    for (const auto& collision : collisions) {
      if (std::abs(collision.posZ()) > cutzvertex || !collision.sel8()) {
        continue;
      }

      rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
      auto kinkCandPerColl = KinkCands.sliceBy(mPerCol, collision.globalIndex());

      for (const auto& kinkCand : kinkCandPerColl) {

        auto dauTrack = kinkCand.trackDaug_as<TracksFull>();
        auto mothTrack = kinkCand.trackMoth_as<TracksFull>();
        if (dauTrack.sign() != mothTrack.sign()) {
          LOG(info) << "Skipping kink candidate with opposite sign daughter and mother: " << kinkCand.globalIndex();
          continue; // Skip if the daughter has the opposite sign as the mother
        }
        if (std::abs(dauTrack.tpcNSigmaPi()) > cutNSigmaPi) {
          continue;
        }

        // histograms filled with all kink candidates
        float alphaAPValue = alphaAP(std::array{kinkCand.pxMoth(), kinkCand.pyMoth(), kinkCand.pzMoth()}, std::array{kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug()});
        float qtValue = qtAP(std::array{kinkCand.pxMoth(), kinkCand.pyMoth(), kinkCand.pzMoth()}, std::array{kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug()});
        rSigmaMinus.fill(HIST("h2ArmenterosPreCuts"), alphaAPValue, qtValue);
        rSigmaMinus.fill(HIST("h2MassSigmaMinusPt"), kinkCand.mothSign() * kinkCand.ptMoth(), kinkCand.mSigmaMinus());
        rSigmaMinus.fill(HIST("h2SigmaMassVsXiMass"), kinkCand.mXiMinus(), kinkCand.mSigmaMinus());
        rSigmaMinus.fill(HIST("h2NSigmaTPCPiPt"), kinkCand.mothSign() * kinkCand.ptMoth(), dauTrack.tpcNSigmaPi());

        // do MC association
        auto mcLabSigma = trackLabelsMC.rawIteratorAt(mothTrack.globalIndex());
        auto mcLabPiDau = trackLabelsMC.rawIteratorAt(dauTrack.globalIndex());
        if (mcLabSigma.has_mcParticle() && mcLabPiDau.has_mcParticle()) {
          auto mcTrackSigma = mcLabSigma.mcParticle_as<aod::McParticles>();
          auto mcTrackPiDau = mcLabPiDau.mcParticle_as<aod::McParticles>();
          if (!mcTrackPiDau.has_mothers()) {
            continue;
          }
          for (auto& piMother : mcTrackPiDau.mothers_as<aod::McParticles>()) {
            if (piMother.globalIndex() != mcTrackSigma.globalIndex()) {
              continue;
            }
            if (std::abs(mcTrackSigma.pdgCode()) != 3112 && std::abs(mcTrackSigma.pdgCode()) != 3222) {
              continue;
            }
            if (std::abs(mcTrackPiDau.pdgCode()) != 211 && std::abs(mcTrackPiDau.pdgCode()) != 2212) {
              continue;
            }

            float MotherMassMC = std::sqrt(piMother.e() * piMother.e() - piMother.p() * piMother.p());
            float MotherpTMC = piMother.pt();
            float deltaXMother = mcTrackPiDau.vx() - piMother.vx();
            float deltaYMother = mcTrackPiDau.vy() - piMother.vy();
            float decayRadiusMC = std::sqrt(deltaXMother * deltaXMother + deltaYMother * deltaYMother);
            float decayRadiusRec = std::sqrt(kinkCand.xDecVtx() * kinkCand.xDecVtx() + kinkCand.yDecVtx() * kinkCand.yDecVtx());
            float cosPointingAngleRec = cosPAngle(std::array{kinkCand.pxMoth(), kinkCand.pyMoth(), kinkCand.pzMoth()},
                                                  std::array{0.0f, 0.0f, 0.0f},
                                                  std::array{kinkCand.xDecVtx(), kinkCand.yDecVtx(), kinkCand.zDecVtx()});

            // Check coherence of MCcollision Id for daughter MCparticle and reconstructed collision
            bool mcCollisionIdCheck = false;
            if (collision.has_mcCollision()) {
              mcCollisionIdCheck = collision.mcCollision().globalIndex() == mcTrackPiDau.mcCollisionId();
            }
            // Check bunch crossing ID coherence
            auto mcCollision = mcTrackPiDau.template mcCollision_as<aod::McCollisions>();
            // bool BCId_vs_MCBCId = collision.bcId() == mcCollision.bcId();
            bool BCId_vs_EvSel = collision.bcId() == collision.foundBCId();
            bool EvSel_vs_MCBCId = collision.foundBCId() == mcCollision.bcId();

            rSigmaMinus.fill(HIST("hMcCollIdCoherence"), static_cast<int>(mcCollisionIdCheck));
            rSigmaMinus.fill(HIST("h2CollId_BCId"), static_cast<int>(mcCollisionIdCheck), static_cast<int>(EvSel_vs_MCBCId));
            rSigmaMinus.fill(HIST("h2BCId_comp"), static_cast<int>(EvSel_vs_MCBCId), static_cast<int>(BCId_vs_EvSel));

            rSigmaMinus.fill(HIST("h2MassPtMCRec"), kinkCand.mothSign() * kinkCand.ptMoth(), kinkCand.mSigmaMinus());
            rSigmaMinus.fill(HIST("h2MassResolution"), kinkCand.mothSign() * kinkCand.ptMoth(), (kinkCand.mSigmaMinus() - MotherMassMC) / MotherMassMC);
            rSigmaMinus.fill(HIST("h2PtResolution"), kinkCand.mothSign() * kinkCand.ptMoth(), (kinkCand.ptMoth() - MotherpTMC) / MotherpTMC);
            rSigmaMinus.fill(HIST("h2RadiusResolution"), kinkCand.mothSign() * kinkCand.ptMoth(), (decayRadiusRec - decayRadiusMC) / decayRadiusMC);
            rSigmaMinus.fill(HIST("h2DCAMothPt"), kinkCand.mothSign() * kinkCand.ptMoth(), kinkCand.dcaMothPv());
            rSigmaMinus.fill(HIST("h2DCADaugPt"), kinkCand.mothSign() * kinkCand.ptMoth(), kinkCand.dcaDaugPv());
            rSigmaMinus.fill(HIST("h2CosPointingAnglePt"), kinkCand.mothSign() * kinkCand.ptMoth(), cosPointingAngleRec);
            rSigmaMinus.fill(HIST("h2ArmenterosPostCuts"), alphaAPValue, qtValue);
            if (std::abs(mcTrackPiDau.pdgCode()) == 211) {
              rSigmaMinus.fill(HIST("h2NSigmaTOFPiPt"), kinkCand.mothSign() * kinkCand.ptMoth(), dauTrack.tofNSigmaPi());
            } else if (std::abs(mcTrackPiDau.pdgCode()) == 2212) {
              rSigmaMinus.fill(HIST("h2NSigmaTOFPrPt"), kinkCand.mothSign() * kinkCand.ptMoth(), dauTrack.tofNSigmaPr());
            }

            // fill the output table with Mc information
            if (fillOutputTree) {
              outputDataTableMC(kinkCand.xDecVtx(), kinkCand.yDecVtx(), kinkCand.zDecVtx(),
                                kinkCand.pxMoth(), kinkCand.pyMoth(), kinkCand.pzMoth(),
                                kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug(),
                                kinkCand.dcaMothPv(), kinkCand.dcaDaugPv(), kinkCand.dcaKinkTopo(),
                                kinkCand.mothSign(),
                                dauTrack.tpcNSigmaPi(), dauTrack.tpcNSigmaPr(), dauTrack.tpcNSigmaKa(),
                                dauTrack.tofNSigmaPi(), dauTrack.tofNSigmaPr(), dauTrack.tofNSigmaKa(),
                                mcTrackSigma.pdgCode(), mcTrackPiDau.pdgCode(),
                                MotherpTMC, MotherMassMC, decayRadiusMC, mcCollisionIdCheck);
            }
          }
        } // MC association and selection
      } // kink cand loop
    } // collision loop

    // Loop over all generated particles to fill MC histograms
    for (const auto& mcPart : particlesMC) {
      if ((std::abs(mcPart.pdgCode()) != 3112 && std::abs(mcPart.pdgCode()) != 3222) || std::abs(mcPart.y()) > cutRapMotherMC) { // only sigma mothers and rapidity cut
        continue;
      }

      if (mcPart.pt() < cutPtGen) {
        continue; // Skip if pT is below threshold
      }

      if (!mcPart.has_daughters()) {
        continue; // Skip if no daughters
      }
      bool hasSigmaDaughter = false;
      int daug_pdg = 0;
      std::array<float, 3> secVtx;
      std::array<float, 3> momDaug;
      for (const auto& daughter : mcPart.daughters_as<aod::McParticles>()) {
        if (std::abs(daughter.pdgCode()) == 211 || std::abs(daughter.pdgCode()) == 2212) { // Pi or proton daughter
          hasSigmaDaughter = true;
          secVtx = {daughter.vx(), daughter.vy(), daughter.vz()};
          momDaug = {daughter.px(), daughter.py(), daughter.pz()};
          daug_pdg = daughter.pdgCode();
          break; // Found a daughter, exit loop
        }
      }
      if (!hasSigmaDaughter) {
        continue; // Skip if no pi/proton daughter found
      }
      float mcMass = std::sqrt(mcPart.e() * mcPart.e() - mcPart.p() * mcPart.p());
      float mcDecayRadius = std::sqrt((secVtx[0] - mcPart.vx()) * (secVtx[0] - mcPart.vx()) + (secVtx[1] - mcPart.vy()) * (secVtx[1] - mcPart.vy()));
      int sigmaSign = mcPart.pdgCode() > 0 ? 1 : -1; // Determine the sign of the Sigma
      rSigmaMinus.fill(HIST("h2MassPtMCGen"), sigmaSign * mcPart.pt(), mcMass);

      // Fill output table with non reconstructed MC candidates
      if (fillOutputTree) {
        outputDataTableMC(-999, -999, -999,
                          -999, -999, -999,
                          -999, -999, -999,
                          -999, -999, -999,
                          sigmaSign,
                          -999, -999, -999,
                          -999, -999, -999,
                          mcPart.pdgCode(), daug_pdg,
                          mcPart.pt(), mcMass, mcDecayRadius, false);
      }
    }
  }
  PROCESS_SWITCH(sigmaminustask, processMC, "MC processing", false);

  void fillFindableHistograms(int filterIndex, float mcRadius, float recRadius, float ptMoth, float ptDaug, bool isSigmaMinus, bool isPiDaughter)
  {
    rFindable.fill(HIST("hFilterIndex"), filterIndex);
    rFindable.fill(HIST("h2MCRadiusFilterIndex"), filterIndex, mcRadius);
    rFindable.fill(HIST("h2RecRadiusFilterIndex"), filterIndex, recRadius);

    if (isPiDaughter) {
      if (isSigmaMinus) {
        rFindable.fill(HIST("h2MCRadiusFilter_minus_pikink"), filterIndex, mcRadius);
        rFindable.fill(HIST("h2PtFilter_minus_pikink"), filterIndex, ptMoth);
        rFindable.fill(HIST("h2PtDaugFilter_minus_pikink"), filterIndex, ptDaug);
      } else {
        rFindable.fill(HIST("h2MCRadiusFilter_plus_pikink"), filterIndex, mcRadius);
        rFindable.fill(HIST("h2PtFilter_plus_pikink"), filterIndex, ptMoth);
        rFindable.fill(HIST("h2PtDaugFilter_plus_pikink"), filterIndex, ptDaug);
      }
    } else {
      if (!isSigmaMinus) {
        rFindable.fill(HIST("h2MCRadiusFilter_plus_protonkink"), filterIndex, mcRadius);
        rFindable.fill(HIST("h2PtFilter_plus_protonkink"), filterIndex, ptMoth);
        rFindable.fill(HIST("h2PtDaugFilter_plus_protonkink"), filterIndex, ptDaug);
      }
    }
  }

  void processFindable(aod::KinkCands const& kinkCands, aod::McTrackLabels const& trackLabelsMC,
                       TracksFull const& tracks, aod::McParticles const&, CollisionsFullMC const&, aod::BCs const&)
  {
    // A - generated findable track pairs map: mcMother.globalIndex() -> (motherTrack.globalIndex(), daughterTrack.globalIndex())
    std::unordered_map<int64_t, std::pair<int64_t, int64_t>> allCandsIndices;

    for (const auto& track : tracks) {
      auto mcLabel = trackLabelsMC.rawIteratorAt(track.globalIndex());
      if (!mcLabel.has_mcParticle()) {
        continue;
      }
      auto mcParticle = mcLabel.mcParticle_as<aod::McParticles>();

      if (mcParticle.has_daughters() && (std::abs(mcParticle.pdgCode()) == 3112 || std::abs(mcParticle.pdgCode()) == 3222)) {
        allCandsIndices[mcParticle.globalIndex()] = {track.globalIndex(), -1};
      }
    }

    for (const auto& track : tracks) {
      auto mcLabel = trackLabelsMC.rawIteratorAt(track.globalIndex());
      if (!mcLabel.has_mcParticle()) {
        continue;
      }
      auto mcParticle = mcLabel.mcParticle_as<aod::McParticles>();

      if (mcParticle.has_mothers() && (std::abs(mcParticle.pdgCode()) == 211 || std::abs(mcParticle.pdgCode()) == 2212)) {
        for (const auto& mother : mcParticle.mothers_as<aod::McParticles>()) {
          auto it = allCandsIndices.find(mother.globalIndex());
          if (it != allCandsIndices.end()) {
            it->second.second = track.globalIndex();
            break;
          }
        }
      }
    }

    // B - reconstructed kinkcands map: mcMother.globalIndex() -> kinkCand.globalIndex()
    std::unordered_map<int64_t, int64_t> findableToKinkCand;
    for (const auto& kinkCand : kinkCands) {
      auto motherTrack = kinkCand.trackMoth_as<TracksFull>();
      auto daughterTrack = kinkCand.trackDaug_as<TracksFull>();

      auto mcLabMoth = trackLabelsMC.rawIteratorAt(motherTrack.globalIndex());
      auto mcLabDaug = trackLabelsMC.rawIteratorAt(daughterTrack.globalIndex());
      if (!mcLabMoth.has_mcParticle() || !mcLabDaug.has_mcParticle()) {
        continue;
      }
      auto mcMother = mcLabMoth.mcParticle_as<aod::McParticles>();
      auto mcDaughter = mcLabDaug.mcParticle_as<aod::McParticles>();

      if (std::abs(mcMother.pdgCode()) != 3112 && std::abs(mcMother.pdgCode()) != 3222) {
        continue;
      }
      if (std::abs(mcDaughter.pdgCode()) != 211 && std::abs(mcDaughter.pdgCode()) != 2212) {
        continue;
      }

      auto findableIt = allCandsIndices.find(mcMother.globalIndex());
      if (findableIt != allCandsIndices.end() &&
          findableIt->second.first == motherTrack.globalIndex() &&
          findableIt->second.second == daughterTrack.globalIndex()) {

        findableToKinkCand[mcMother.globalIndex()] = kinkCand.globalIndex();
      }
    }

    // C - loop on valid pairs for findable analysis
    for (const auto& [mcMotherIndex, trackIndices] : allCandsIndices) {
      if (trackIndices.second == -1 || trackIndices.first == -1) {
        continue;
      }

      // Retrieve mother and daughter tracks and mcParticles
      auto motherTrack = tracks.rawIteratorAt(trackIndices.first);
      auto daughterTrack = tracks.rawIteratorAt(trackIndices.second);
      auto mcLabMoth = trackLabelsMC.rawIteratorAt(motherTrack.globalIndex());
      auto mcLabDaug = trackLabelsMC.rawIteratorAt(daughterTrack.globalIndex());
      auto mcMother = mcLabMoth.mcParticle_as<aod::McParticles>();
      auto mcDaughter = mcLabDaug.mcParticle_as<aod::McParticles>();

      // Compute useful quantities for histograms
      bool isSigmaMinus = (std::abs(mcMother.pdgCode()) == 3112);
      bool isPiDaughter = (std::abs(mcDaughter.pdgCode()) == 211);

      float recPtDaughter = daughterTrack.pt();
      float recPtMother = motherTrack.pt();
      float mcRadius = std::sqrt((mcMother.vx() - mcDaughter.vx()) * (mcMother.vx() - mcDaughter.vx()) + (mcMother.vy() - mcDaughter.vy()) * (mcMother.vy() - mcDaughter.vy()));
      float recRadius = -1.0;
      if (findableToKinkCand.find(mcMother.globalIndex()) != findableToKinkCand.end()) {
        auto kinkCand = kinkCands.rawIteratorAt(findableToKinkCand[mcMother.globalIndex()]);
        recRadius = std::sqrt(kinkCand.xDecVtx() * kinkCand.xDecVtx() + kinkCand.yDecVtx() * kinkCand.yDecVtx());
      }

      // Check for detector mismatches in ITS mother tracks
      auto mask_value = mcLabMoth.mcMask();
      int mismatchITS_index = -1;

      for (int i = 0; i < 7; ++i) { // ITS has layers 0-6, bit ON means mismatch
        if ((mask_value & (1 << i)) != 0) {
          mismatchITS_index = i;
          break;
        }
      }

      // Define filter index and progressively apply kinkbuilder cuts to track pairs
      int filterIndex = 0;
      fillFindableHistograms(filterIndex, mcRadius, recRadius, recPtMother, recPtDaughter, isSigmaMinus, isPiDaughter);

      // 1 - tracks with right ITS, TPC, TOF signals
      if (motherTrack.has_collision() && motherTrack.hasITS() && !motherTrack.hasTPC() && !motherTrack.hasTOF() &&
          daughterTrack.hasITS() && daughterTrack.hasTPC()) {
        filterIndex += 1;
        fillFindableHistograms(filterIndex, mcRadius, recRadius, recPtMother, recPtDaughter, isSigmaMinus, isPiDaughter);
        rFindable.fill(HIST("hfakeITSfindable"), mismatchITS_index);
      } else {
        continue;
      }

      // 2 - moth+daug track quality cuts
      bool motherGoodITS = motherTrack.hasITS() && motherTrack.itsNCls() < 6 && motherTrack.itsNClsInnerBarrel() == 3 && motherTrack.itsChi2NCl() < 36;
      bool daughterGoodITSTPC = daughterTrack.hasITS() && daughterTrack.hasTPC() && daughterTrack.itsNClsInnerBarrel() == 0 &&
                                daughterTrack.itsNCls() < 4 && daughterTrack.tpcNClsCrossedRows() > 0.8 * daughterTrack.tpcNClsFindable() && daughterTrack.tpcNClsFound() > nTPCClusMinDaugKB;
      if (motherGoodITS && daughterGoodITSTPC) {
        filterIndex += 1;
        fillFindableHistograms(filterIndex, mcRadius, recRadius, recPtMother, recPtDaughter, isSigmaMinus, isPiDaughter);
      } else {
        continue;
      }

      // 3 - mother track min pT
      if (motherTrack.pt() > minPtMothKB) {
        filterIndex += 1;
        fillFindableHistograms(filterIndex, mcRadius, recRadius, recPtMother, recPtDaughter, isSigmaMinus, isPiDaughter);
      } else {
        continue;
      }

      // 4 - geometric cuts: eta
      if (std::abs(motherTrack.eta()) < etaMaxKB && std::abs(daughterTrack.eta()) < etaMaxKB) {
        filterIndex += 1;
        fillFindableHistograms(filterIndex, mcRadius, recRadius, recPtMother, recPtDaughter, isSigmaMinus, isPiDaughter);
      } else {
        continue;
      }

      // 5 - geometric cuts: phi difference
      if (std::abs(motherTrack.phi() - daughterTrack.phi()) * radToDeg < maxPhiDiffKB) {
        filterIndex += 1;
        fillFindableHistograms(filterIndex, mcRadius, recRadius, recPtMother, recPtDaughter, isSigmaMinus, isPiDaughter);
      } else {
        continue;
      }

      // DCA calculation: initialization CCDB
      auto collision = motherTrack.template collision_as<CollisionsFullMC>();
      auto bc = collision.template bc_as<aod::BCs>();
      initCCDB(bc);
      const o2::math_utils::Point3D<float> collVtx{collision.posX(), collision.posY(), collision.posZ()};
      o2::track::TrackParCov trackParCovMoth = getTrackParCov(motherTrack);
      o2::track::TrackParCov trackParCovDaug = getTrackParCov(daughterTrack);

      // get DCA to PV for mother and daughter tracks
      std::array<float, 2> dcaInfoMoth;
      o2::base::Propagator::Instance()->propagateToDCA(collVtx, trackParCovMoth, mBz, 2.f, static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value), &dcaInfoMoth);
      std::array<float, 2> dcaInfoDaug;
      o2::base::Propagator::Instance()->propagateToDCA(collVtx, trackParCovDaug, mBz, 2.f, static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value), &dcaInfoDaug);
      float dcaXYMother = std::abs(dcaInfoMoth[0]);
      float dcaXYDaughter = std::abs(dcaInfoDaug[0]);
      if (isPiDaughter) {
        rFindable.fill(HIST("h2DCAMothPt_pikink"), recPtMother, dcaXYMother);
        rFindable.fill(HIST("h2DCADaugPt_pikink"), recPtDaughter, dcaXYDaughter);
      } else {
        rFindable.fill(HIST("h2DCAMothPt_protonkink"), recPtMother, dcaXYMother);
        rFindable.fill(HIST("h2DCADaugPt_protonkink"), recPtDaughter, dcaXYDaughter);
      }

      // 6 - max Z difference
      if (std::abs(trackParCovMoth.getZ() - trackParCovDaug.getZ()) < maxZDiffKB) {
        filterIndex += 1;
        fillFindableHistograms(filterIndex, mcRadius, recRadius, recPtMother, recPtDaughter, isSigmaMinus, isPiDaughter);
      } else {
        continue;
      }

      // 7 - DCA mother
      if (dcaXYMother < maxDcaMothPvKB) {
        filterIndex += 1;
        fillFindableHistograms(filterIndex, mcRadius, recRadius, recPtMother, recPtDaughter, isSigmaMinus, isPiDaughter);
      } else {
        continue;
      }

      // 8 - DCA daughter
      if (dcaXYDaughter > minDcaDaugPvKB) {
        filterIndex += 1;
        fillFindableHistograms(filterIndex, mcRadius, recRadius, recPtMother, recPtDaughter, isSigmaMinus, isPiDaughter);
      } else {
        continue;
      }

      // 9 - radius cut
      if (recRadius > radiusCutKB) {
        filterIndex += 1;
        fillFindableHistograms(filterIndex, mcRadius, recRadius, recPtMother, recPtDaughter, isSigmaMinus, isPiDaughter);
      } else {
        continue;
      }

      // 10 - collision selection
      if (!(std::abs(collision.posZ()) > cutzvertex || !collision.sel8())) {
        filterIndex += 1;
        fillFindableHistograms(filterIndex, mcRadius, recRadius, recPtMother, recPtDaughter, isSigmaMinus, isPiDaughter);
      } else {
        continue;
      }

      // 11 - TOF daughter presence
      if (daughterTrack.hasTOF()) {
        filterIndex += 1;
        fillFindableHistograms(filterIndex, mcRadius, recRadius, recPtMother, recPtDaughter, isSigmaMinus, isPiDaughter);
      }
    }
  }

  PROCESS_SWITCH(sigmaminustask, processFindable, "Findable Sigma processing", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<sigmaminustask>(cfgc)};
}
