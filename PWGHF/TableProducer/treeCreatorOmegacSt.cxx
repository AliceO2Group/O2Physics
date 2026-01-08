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

/// \file treeCreatorOmegacSt.cxx
/// \brief Task to reconstruct Ωc from strangeness-tracked Ω and pion/kaon
///
/// \author Jochen Klein
/// \author Tiantian Cheng
/// \author Ruiqi Yin

#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGHF/Utils/utilsTrkCandHf.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/Track.h>

#include <TH1.h>
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_trkcandsel;
using namespace o2::hf_evsel;
using namespace o2::hf_centrality;

namespace o2::aod
{
namespace hf_st_charmed_baryon_gen
{
DECLARE_SOA_COLUMN(PxCharmedBaryon, pxCharmedBaryon, float);
DECLARE_SOA_COLUMN(PyCharmedBaryon, pyCharmedBaryon, float);
DECLARE_SOA_COLUMN(PzCharmedBaryon, pzCharmedBaryon, float);
DECLARE_SOA_COLUMN(PdgCodeCharmedBaryon, pdgCodeCharmedBaryon, int);
DECLARE_SOA_COLUMN(PxCasc, pxCasc, float);
DECLARE_SOA_COLUMN(PyCasc, pyCasc, float);
DECLARE_SOA_COLUMN(PzCasc, pzCasc, float);
DECLARE_SOA_COLUMN(PdgCodeCasc, pdgCodeCasc, int);
DECLARE_SOA_COLUMN(DecayLengthCharmedBaryon, decayLengthCharmedBaryon, float);
DECLARE_SOA_COLUMN(DecayLengthXYCharmedBaryon, decayLengthXYCharmedBaryon, float);
DECLARE_SOA_COLUMN(DecayLengthCasc, decayLengthCasc, float);
DECLARE_SOA_COLUMN(DecayLengthXYCasc, decayLengthXYCasc, float);
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int);
DECLARE_SOA_COLUMN(DecayChannel, decayChannel, int);
} // namespace hf_st_charmed_baryon_gen

DECLARE_SOA_TABLE(HfStChBarGens, "AOD", "HFSTCHBARGEN",
                  hf_st_charmed_baryon_gen::PxCharmedBaryon,
                  hf_st_charmed_baryon_gen::PyCharmedBaryon,
                  hf_st_charmed_baryon_gen::PzCharmedBaryon,
                  hf_st_charmed_baryon_gen::PdgCodeCharmedBaryon,
                  hf_st_charmed_baryon_gen::PxCasc,
                  hf_st_charmed_baryon_gen::PyCasc,
                  hf_st_charmed_baryon_gen::PzCasc,
                  hf_st_charmed_baryon_gen::PdgCodeCasc,
                  hf_st_charmed_baryon_gen::DecayLengthCharmedBaryon,
                  hf_st_charmed_baryon_gen::DecayLengthXYCharmedBaryon,
                  hf_st_charmed_baryon_gen::DecayLengthCasc,
                  hf_st_charmed_baryon_gen::DecayLengthXYCasc,
                  hf_st_charmed_baryon_gen::OriginMcGen,
                  hf_st_charmed_baryon_gen::DecayChannel);

// CharmedBaryon -> Casc + Pion/Kaon
//                   -> Lambda + BachPi/BachKa
//                        -> Pr + Pi
namespace hf_st_charmed_baryon
{
DECLARE_SOA_COLUMN(Centrality, centrality, float);
DECLARE_SOA_COLUMN(MassOmega, massOmega, float);
DECLARE_SOA_COLUMN(MassXi, massXi, float);
DECLARE_SOA_COLUMN(MassLambda, massLambda, float);
DECLARE_SOA_COLUMN(NSigmaTpcPion, nSigmaTpcPion, float);
DECLARE_SOA_COLUMN(NSigmaTofPion, nSigmaTofPion, float);
DECLARE_SOA_COLUMN(NSigmaTpcKaon, nSigmaTpcKaon, float);
DECLARE_SOA_COLUMN(NSigmaTofKaon, nSigmaTofKaon, float);
DECLARE_SOA_COLUMN(NSigmaTpcV0Pr, nSigmaTpcV0Pr, float);
DECLARE_SOA_COLUMN(NSigmaTofV0Pr, nSigmaTofV0Pr, float);
DECLARE_SOA_COLUMN(NSigmaTpcV0Pi, nSigmaTpcV0Pi, float);
DECLARE_SOA_COLUMN(NSigmaTofV0Pi, nSigmaTofV0Pi, float);
DECLARE_SOA_COLUMN(NSigmaTpcBachPi, nSigmaTpcBachPi, float);
DECLARE_SOA_COLUMN(NSigmaTofBachPi, nSigmaTofBachPi, float);
DECLARE_SOA_COLUMN(NSigmaTpcBachKa, nSigmaTpcBachKa, float);
DECLARE_SOA_COLUMN(NSigmaTofBachKa, nSigmaTofBachKa, float);
DECLARE_SOA_COLUMN(PxCasc, pxCasc, float);
DECLARE_SOA_COLUMN(PyCasc, pyCasc, float);
DECLARE_SOA_COLUMN(PzCasc, pzCasc, float);
DECLARE_SOA_COLUMN(IsPositiveCasc, isPositiveCasc, bool);
DECLARE_SOA_COLUMN(PxPionOrKaon, pxPionOrKaon, float);
DECLARE_SOA_COLUMN(PyPionOrKaon, pyPionOrKaon, float);
DECLARE_SOA_COLUMN(PzPionOrKaon, pzPionOrKaon, float);
DECLARE_SOA_COLUMN(IsPositivePionOrKaon, isPositivePionOrKaon, bool);
DECLARE_SOA_COLUMN(ItsClusterMapPionOrKaon, itsClusterMapPionOrKaon, uint8_t);
DECLARE_SOA_COLUMN(CpaCharmedBaryon, cpaCharmedBaryon, float);
DECLARE_SOA_COLUMN(CpaXYCharmedBaryon, cpaXYCharmedBaryon, float);
DECLARE_SOA_COLUMN(CpaCasc, cpaCasc, float);
DECLARE_SOA_COLUMN(CpaXYCasc, cpaXYCasc, float);
DECLARE_SOA_COLUMN(DcaXYCasc, dcaXYCasc, float);
DECLARE_SOA_COLUMN(DcaXYUncCasc, dcaXYUncCasc, float);
DECLARE_SOA_COLUMN(DcaZCasc, dcaZCasc, float);
DECLARE_SOA_COLUMN(DcaZUncCasc, dcaZUncCasc, float);
DECLARE_SOA_COLUMN(DcaXYPionOrKaon, dcaXYPionOrKaon, float);
DECLARE_SOA_COLUMN(DcaXYUncPionOrKaon, dcaXYUncPionOrKaon, float);
DECLARE_SOA_COLUMN(DcaZPionOrKaon, dcaZPionOrKaon, float);
DECLARE_SOA_COLUMN(DcaZUncPionOrKaon, dcaZUncPionOrKaon, float);
DECLARE_SOA_COLUMN(DcaXYPr, dcaXYPr, float);
DECLARE_SOA_COLUMN(DcaZPr, dcaZPr, float);
DECLARE_SOA_COLUMN(DcaXYKa, dcaXYKa, float);
DECLARE_SOA_COLUMN(DcaZKa, dcaZKa, float);
DECLARE_SOA_COLUMN(DcaXYPi, dcaXYPi, float);
DECLARE_SOA_COLUMN(DcaZPi, dcaZPi, float);
DECLARE_SOA_COLUMN(Chi2TopologicalCharmedBaryon, chi2TopologicalCharmedBaryon, float);
DECLARE_SOA_COLUMN(Chi2TopologicalCasc, chi2TopologicalCasc, float);
DECLARE_SOA_COLUMN(DecayLengthCharmedBaryon, decayLengthCharmedBaryon, float);
DECLARE_SOA_COLUMN(DecayLengthXYCharmedBaryon, decayLengthXYCharmedBaryon, float);
DECLARE_SOA_COLUMN(DecayLengthCharmedBaryonUntracked, decayLengthCharmedBaryonUntracked, float);
DECLARE_SOA_COLUMN(DecayLengthXYCharmedBaryonUntracked, decayLengthXYCharmedBaryonUntracked, float);
DECLARE_SOA_COLUMN(DecayLengthCasc, decayLengthCasc, float);
DECLARE_SOA_COLUMN(DecayLengthXYCasc, decayLengthXYCasc, float);
DECLARE_SOA_INDEX_COLUMN_FULL(MotherCasc, motherCasc, int, HfStChBarGens, "_Casc");
DECLARE_SOA_INDEX_COLUMN_FULL(MotherPionOrKaon, motherPionOrKaon, int, HfStChBarGens, "_PionOrKaon");
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int);
DECLARE_SOA_COLUMN(ToiMask, toiMask, uint32_t);
} // namespace hf_st_charmed_baryon

DECLARE_SOA_TABLE(HfStChBars, "AOD", "HFSTCHBAR",
                  hf_st_charmed_baryon::Centrality,
                  hf_st_charmed_baryon::MassOmega,
                  hf_st_charmed_baryon::MassXi,
                  hf_st_charmed_baryon::MassLambda,
                  hf_st_charmed_baryon::NSigmaTpcPion,
                  hf_st_charmed_baryon::NSigmaTofPion,
                  hf_st_charmed_baryon::NSigmaTpcKaon,
                  hf_st_charmed_baryon::NSigmaTofKaon,
                  hf_st_charmed_baryon::NSigmaTpcV0Pr,
                  hf_st_charmed_baryon::NSigmaTofV0Pr,
                  hf_st_charmed_baryon::NSigmaTpcV0Pi,
                  hf_st_charmed_baryon::NSigmaTofV0Pi,
                  hf_st_charmed_baryon::NSigmaTpcBachPi,
                  hf_st_charmed_baryon::NSigmaTofBachPi,
                  hf_st_charmed_baryon::NSigmaTpcBachKa,
                  hf_st_charmed_baryon::NSigmaTofBachKa,
                  hf_st_charmed_baryon::PxCasc,
                  hf_st_charmed_baryon::PyCasc,
                  hf_st_charmed_baryon::PzCasc,
                  hf_st_charmed_baryon::IsPositiveCasc,
                  hf_st_charmed_baryon::PxPionOrKaon,
                  hf_st_charmed_baryon::PyPionOrKaon,
                  hf_st_charmed_baryon::PzPionOrKaon,
                  hf_st_charmed_baryon::IsPositivePionOrKaon,
                  hf_st_charmed_baryon::ItsClusterMapPionOrKaon,
                  hf_st_charmed_baryon::CpaCharmedBaryon,
                  hf_st_charmed_baryon::CpaXYCharmedBaryon,
                  hf_st_charmed_baryon::CpaCasc,
                  hf_st_charmed_baryon::CpaXYCasc,
                  hf_st_charmed_baryon::DcaXYCasc,
                  hf_st_charmed_baryon::DcaXYUncCasc,
                  hf_st_charmed_baryon::DcaZCasc,
                  hf_st_charmed_baryon::DcaZUncCasc,
                  hf_st_charmed_baryon::DcaXYPionOrKaon,
                  hf_st_charmed_baryon::DcaXYUncPionOrKaon,
                  hf_st_charmed_baryon::DcaZPionOrKaon,
                  hf_st_charmed_baryon::DcaZUncPionOrKaon,
                  hf_st_charmed_baryon::DcaXYPr,
                  hf_st_charmed_baryon::DcaZPr,
                  hf_st_charmed_baryon::DcaXYKa,
                  hf_st_charmed_baryon::DcaZKa,
                  hf_st_charmed_baryon::DcaXYPi,
                  hf_st_charmed_baryon::DcaZPi,
                  hf_st_charmed_baryon::Chi2TopologicalCharmedBaryon,
                  hf_st_charmed_baryon::Chi2TopologicalCasc,
                  hf_st_charmed_baryon::DecayLengthCharmedBaryon,
                  hf_st_charmed_baryon::DecayLengthXYCharmedBaryon,
                  hf_st_charmed_baryon::DecayLengthCharmedBaryonUntracked,
                  hf_st_charmed_baryon::DecayLengthXYCharmedBaryonUntracked,
                  hf_st_charmed_baryon::DecayLengthCasc,
                  hf_st_charmed_baryon::DecayLengthXYCasc,
                  hf_st_charmed_baryon::MotherCascId,
                  hf_st_charmed_baryon::MotherPionOrKaonId,
                  hf_st_charmed_baryon::OriginMcRec,
                  hf_st_charmed_baryon::ToiMask);
} // namespace o2::aod

struct HfTreeCreatorOmegacSt {
  Produces<aod::HfStChBars> outputTable;
  Produces<aod::HfStChBarGens> outputTableGen;

  Configurable<int> materialCorrectionType{"materialCorrectionType", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrLUT), "Type of material correction"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpMagPath{"grpMagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> matLutPath{"matLutPath", "GLO/Param/MatLUT", "Path of the material LUT"};
  Configurable<bool> propToDCA{"propToDCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<bool> skimmedProcessing{"skimmedProcessing", false, "Put true if you are processing apass*_skimmed datasets"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations if chi2/chi2old > this"};
  Configurable<int> minNoClsTrackedCascade{"minNoClsTrackedCascade", 70, "Minimum number of clusters required for daughters of tracked cascades"};
  Configurable<int> minNoClsTrackedPionOrKaon{"minNoClsTrackedPionOrKaon", 70, "Minimum number of clusters required for associated pions/kaons"};
  Configurable<bool> useSel8Trigger{"useSel8Trigger", true, "filter collisions on sel 8 trigger"};
  Configurable<float> massWindowTrackedOmega{"massWindowTrackedOmega", 0.05, "Inv. mass window for tracked Omega"};
  Configurable<float> massWindowXiExclTrackedOmega{"massWindowXiExclTrackedOmega", 0.005, "Inv. mass window for exclusion of Xi for tracked Omega-"};
  Configurable<float> massWindowTrackedXi{"massWindowTrackedXi", 0., "Inv. mass window for tracked Xi"};
  Configurable<float> massWindowLambda{"massWindowLambda", 0.05, "Inv. mass window for Lambda"};
  Configurable<float> massWindowXiC{"massWindowXiC", 0.1, "Inv. mass window for Xic"};
  Configurable<float> massWindowOmegaC{"massWindowOmegaC", 0.1, "Inv. mass window for Omegac"};
  Configurable<int> selectedChannel{"selectedChannel", 0, "Selected decay channel for Xic or Omegac: 0=all, 1=Xic->Xi+Pi, 2=Xic->Omega+K, 3=Omegac->Omega+Pi, 4=Omegac->Omega+K, 5=Omegac->Xi+Pi"};
  Configurable<float> maxMatchingChi2TrackedCascade{"maxMatchingChi2TrackedCascade", 2000., "Max matching chi2 for tracked cascades"};
  Configurable<bool> recalculateMasses{"recalculateMasses", true, "Recalculate Xi/Omega masses"};
  Configurable<float> nSigmaBachelorMax{"nSigmaBachelorMax", 5., "Max Nsigma for bachelor of tracked cascade"};
  Configurable<float> nSigmaV0PrMax{"nSigmaV0PrMax", 5., "Max Nsigma for proton from V0 from tracked cascade"};
  Configurable<float> nSigmaV0PiMax{"nSigmaV0PiMax", 5., "Max Nsigma for pion from V0 from tracked cascade"};
  Configurable<float> nSigmaPionMax{"nSigmaPionMax", 5., "Max Nsigma for pion to be paired with Omega"};
  Configurable<float> nSigmaKaonMax{"nSigmaKaonMax", 5., "Max Nsigma for kaon to be paired with Omega"};
  Configurable<bool> bzOnly{"bzOnly", true, "Use B_z instead of full field map"};
  Configurable<std::string> cfgTriggersOfInterest{"cfgTriggersOfInterest", "fTrackedOmega,fHfCharmBarToXiBach", "Triggers of interest, comma separated for Zorro"};

  // TOF nsigma configurable cuts (defaults: very loose = 9999)
  Configurable<float> nSigmaAbsTofTrackPiMax{"nSigmaAbsTofTrackPiMax", 9999.f, "Max abs TOF nsigma for charm bachelor track as pion"};
  Configurable<float> nSigmaAbsTofTrackKaMax{"nSigmaAbsTofTrackKaMax", 9999.f, "Max abs TOF nsigma for charm bachelor track as kaon"};
  Configurable<float> nSigmaAbsTofV0PrMax{"nSigmaAbsTofV0PrMax", 9999.f, "Max abs TOF nsigma for V0 proton"};
  Configurable<float> nSigmaAbsTofV0PiMax{"nSigmaAbsTofV0PiMax", 9999.f, "Max abs TOF nsigma for V0 pion"};
  Configurable<float> nSigmaAbsTofBachelorPiMax{"nSigmaAbsTofBachelorPiMax", 9999.f, "Max abs TOF nsigma for cascade bachelor as pion"};
  Configurable<float> nSigmaAbsTofBachelorKaMax{"nSigmaAbsTofBachelorKaMax", 9999.f, "Max abs TOF nsigma for cascade bachelor as kaon"};
  // whether to require TOF in addition to TPC
  Configurable<bool> useTofPid{"useTofPid", false, "Require TOF PID together with TPC PID when true; otherwise use TPC only"};

  // pT cut for CharmBaryon
  Configurable<float> ptCharmedBaryonMin{"ptCharmedBaryonMin", 0.f, "Minimum pT for the charmed baryon"};
  Configurable<float> ptCharmedBaryonMax{"ptCharmedBaryonMax", 50.f, "Maximum pT for the charmed baryon"};

  // CPA cuts (defaults: >= 0.9)
  Configurable<float> cpaCharmedBaryonMin{"cpaCharmedBaryonMin", 0.9f, "Minimum CPA for charmed baryon"};
  Configurable<float> cpaXYCharmedBaryonMin{"cpaXYCharmedBaryonMin", 0.9f, "Minimum CPA XY for charmed baryon"}; // remove maybe???
  Configurable<float> cpaCascMin{"cpaCascMin", 0.9f, "Minimum CPA for cascade"};
  Configurable<float> cpaXYCascMin{"cpaXYCascMin", 0.9f, "Minimum CPA XY for cascade"};

  // Chi2 cuts (defaults: <= 10)
  Configurable<float> chi2TopCharmedBaryonMax{"chi2TopCharmedBaryonMax", 10.f, "Maximum chi2 for top charmed baryon fit"};
  Configurable<float> chi2TopCascMax{"chi2TopCascMax", 10.f, "Maximum topologyChi2 of the tracked cascade"};

  // Impact parameter cuts (defaults: absolute value < 10)
  Configurable<float> absImpactParameterCascYMax{"absImpactParameterCascYMax", 10.f, "Max absolute value of impact parameter of cascade (Y)"};
  Configurable<float> absDcaCascZMax{"absDcaCascZMax", 10.f, "Max absolute value of DCA Z of cascade"};
  Configurable<float> absImpactParameterCharmBachelorYMax{"absImpactParameterCharmBachelorYMax", 10.f, "Max absolute value of impact parameter of charm bachelor (Y)"};
  Configurable<float> absDcaCharmBachelorZMax{"absDcaCharmBachelorZMax", 10.f, "Max absolute value of DCA Z of charm bachelor"};
  Configurable<float> ptCascadeMin{"ptCascadeMin", 0.f, "Minimum pT requirement for the cascade track"};
  Configurable<float> ptPionOrKaonMin{"ptPionOrKaonMin", 0.f, "Minimum pT requirement for the charm baryon daughter pion/kaon track"};
  Configurable<float> impProductMin{"impProductMin", -1.f, "Lower bound for DCA impact product"};
  Configurable<float> impProductMax{"impProductMax", 0.f, "Upper bound for DCA impact product"};

  SliceCache cache;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::vertexing::DCAFitterN<2> df2;
  HfEventSelection hfEvSel;

  static constexpr int ItsNClsMin{4};
  static constexpr float TpcNclsFindableFraction{0.8f};
  static constexpr float TpcChi2NclMax{4.f};
  static constexpr float ItsChi2NclMax{36.f};
  static constexpr std::array<float, 2> TofWoSignalRange{998.f, 1000.f};

  float bz{0.f};
  int runNumber{0};
  std::map<int, int> mapMcPartToGenTable;

  using Collisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>;
  using CollisionsFT0C = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>;
  using CollisionsFT0M = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>>;
  using TracksExt = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr>;
  using TracksExtMc = soa::Join<TracksExt, aod::McTrackLabels>;

  Filter collisionFilter = (useSel8Trigger.node() == false) || (o2::aod::evsel::sel8 == true);

  // Preslice<aod::Tracks> perCol = aod::track::collisionId;
  PresliceUnsorted<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;
  PresliceUnsorted<aod::AssignedTrackedCascades> assignedTrackedCascadesPerCollision = aod::track::collisionId;

  std::shared_ptr<TH1> hCandidatesPrPi, hCandidatesV0Pi, hCandidatesCascPiOrK;
  HistogramRegistry registry{
    "registry",
    {
      {"hDca", "DCA;DCA (cm)", {HistType::kTH1D, {{200, 0., .5}}}},
      {"hDcaXY", "DCA;DCA_{xy} (cm)", {HistType::kTH1D, {{200, -.5, .5}}}},
      {"hDcaXYVsPt", "DCA;p_{T} (GeV/#it{c};DCA_{xy} (cm)", {HistType::kTH2D, {{200, 0., 10.}, {200, -.5, .5}}}},
      {"hDcaZ", "DCA;DCA_{z} (cm)", {HistType::kTH1D, {{200, -.5, .5}}}},
      {"hDcaZVsPt", "DCA;p_{T} (GeV/#it{c});DCA_{z} (cm)", {HistType::kTH2D, {{200, 0., 10.}, {200, -.5, .5}}}},
      {"hDcaVsPt", "DCA;DCA (cm);p_{T} (GeV/#it{c})", {HistType::kTH2D, {{200, 0., .5}, {200, 0., 10.}}}},
      {"hDcaVsR", "DCA;DCA (cm);R (cm)", {HistType::kTH2D, {{200, 0., .5}, {200, 0., 10.}}}},
      {"hDecayLength", "Decay length;L (#mum)", {HistType::kTH1D, {{200, 0., 500.}}}},
      {"hDecayLengthId", "Decay length (true #Omega_{c});L (#mum)", {HistType::kTH1D, {{200, 0., 500.}}}},
      {"hDecayLengthGen", "Decay length (gen);L (#mum)", {HistType::kTH1D, {{200, 0., 500.}}}},
      {"hDeltaDecayLength", "#Delta decay length (gen);#Delta L (#mum)", {HistType::kTH1D, {{200, -250., 250.}}}},
      {"hDecayLengthScaled", "Decay length * M/p;L (#mum / #it{c})", {HistType::kTH1D, {{200, 0., 500.}}}},
      {"hDecayLengthScaledId", "Decay length * M/p (true #Omega_{c});L (#mum / #it{c})", {HistType::kTH1D, {{200, 0., 500.}}}},
      {"hDecayLengthScaledGen", "Decay length * M/p (MC id);L (#mum / #it{c})", {HistType::kTH1D, {{200, 0., 500.}}}},
      {"hDecayLengthScaledMc", "Decay length * M/p (MC);L (#mum / #it{c})", {HistType::kTH1D, {{200, 0., 500.}}}},
      {"hMassOmegaPi", "inv. mass #Omega + #pi;inv. mass (GeV/#it{c}^{2})", {HistType::kTH1D, {{400, 1.5, 3.}}}},
      {"hMassOmegaPiVsPt", "inv. mass #Omega + #pi;inv. mass (GeV/#it{c}^{2});p_{T} (GeV/#it{c})", {HistType::kTH2D, {{400, 1.5, 3.}, {10, 0., 10.}}}},
      {"hMassOmegaK", "inv. mass #Omega + K;inv. mass (GeV/#it{c}^{2})", {HistType::kTH1D, {{400, 1.5, 3.}}}},
      {"hMassOmegaKVsPt", "inv. mass #Omega + K;inv. mass (GeV/#it{c}^{2});p_{T} (GeV/#it{c})", {HistType::kTH2D, {{400, 1.5, 3.}, {10, 0., 10.}}}},
      {"hMassXiPiVsPt", "inv. mass #Xi + #pi;inv. mass (GeV/#it{c}^{2});p_{T} (GeV/#it{c})", {HistType::kTH2D, {{400, 1.5, 3.}, {10, 0., 10.}}}},
      {"hMassOmegacId", "inv. mass #Omega + #pi (MC ID);inv. mass (GeV/#it{c}^{2})", {HistType::kTH1D, {{400, 1.5, 3.}}}},
      {"hMassOmegacGen", "inv. mass #Omega + #pi (from MC);inv. mass (GeV/#it{c}^{2})", {HistType::kTH1D, {{400, 1.5, 3.}}}},
      {"hPtVsMassOmega", "#Omega mass;p_{T} (GeV/#it{c});m (GeV/#it{c}^3)", {HistType::kTH2D, {{200, 0., 10.}, {1000, 1., 3.}}}},
      {"hDeltaPtVsPt", "Delta pt;p_{T} (GeV/#it{c});#Delta p_{T} / p_{T}", {HistType::kTH2D, {{200, 0., 10.}, {200, -1., 1.}}}},
    }};

  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  void init(InitContext const&)
  {
    df2.setPropagateToPCA(propToDCA);
    df2.setMaxR(maxR);
    df2.setMaxDZIni(maxDZIni);
    df2.setMinParamChange(minParamChange);
    df2.setMinRelChi2Change(minRelChi2Change);
    df2.setUseAbsDCA(useAbsDCA);

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    if (static_cast<o2::base::Propagator::MatCorrType>(materialCorrectionType.value) == o2::base::Propagator::MatCorrType::USEMatCorrLUT) {
      auto* lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));
      o2::base::Propagator::Instance(true)->setMatLUT(lut);
    }

    /// candidate monitoring
    hCandidatesPrPi = registry.add<TH1>("hCandidatesPrPi", "Pr-Pi candidates counter", {HistType::kTH1D, {axisCands}});
    hCandidatesV0Pi = registry.add<TH1>("hCandidatesV0Pi", "V0-Pi candidates counter", {HistType::kTH1D, {axisCands}});
    hCandidatesCascPiOrK = registry.add<TH1>("hCandidatesCascPiOrK", "Casc-Pi/K candidates counter", {HistType::kTH1D, {axisCands}});
    setLabelHistoCands(hCandidatesPrPi);
    setLabelHistoCands(hCandidatesV0Pi);
    setLabelHistoCands(hCandidatesCascPiOrK);

    // init HF event selection helper (centrality, event cuts, monitoring)
    hfEvSel.init(registry);
  }

  // processMC: loop over MC objects
  // processData: loop over reconstructed objects, no MC information
  // processGen: loop over reconstructed objects, use MC information (mutually exclusive? combine?)
  int indexRec = -1;
  int indexRecCharmBaryon = -1;
  int8_t sign = -9;
  int8_t signCasc = -9;
  int8_t signV0 = -9;
  int8_t origin = 0; // to be used for prompt/non prompt
  int8_t nPiToMuV0{0}, nPiToMuCasc{0}, nPiToMuOmegac0{0};
  int8_t nKaToPiCasc{0}, nKaToPiOmegac0{0};
  std::vector<int> idxBhadMothers;
  int decayChannel = -1; // flag for different decay channels
  bool isMatched = false;
  static constexpr std::size_t NDaughters{2u};
  enum ChannelType {
    AnyChannel = 0,
    XicToXiPi = 1,
    XicToOmegaKa = 2,
    OmegacToOmegaPi = 3,
    OmegacToOmegaKa = 4,
    OmegacToXiPi = 5
  };

  void processMc(aod::McCollisions const&,
                 aod::McParticles const& mcParticles)
  {
    mapMcPartToGenTable.clear();
    for (const auto& mcParticle : mcParticles) {
      const bool isOmegaC = std::abs(mcParticle.pdgCode()) == constants::physics::Pdg::kOmegaC0;
      const bool isXiC = std::abs(mcParticle.pdgCode()) == constants::physics::Pdg::kXiC0;
      if (isOmegaC || isXiC) {
        const auto daughters = mcParticle.daughters_as<aod::McParticles>();
        if (daughters.size() == NDaughters) {
          int idxPionDaughter = -1;
          int idxCascDaughter = -1;
          int idxKaonDaughter = -1;
          const auto daughters = mcParticle.daughters_as<aod::McParticles>();
          for (const auto& daughter : daughters) {
            if (idxCascDaughter < 0 && (std::abs(daughter.pdgCode()) == (isOmegaC ? kOmegaMinus : kXiMinus))) {
              idxCascDaughter = daughter.globalIndex();
            }
            if (idxPionDaughter < 0 && (std::abs(daughter.pdgCode()) == kPiPlus)) {
              idxPionDaughter = daughter.globalIndex();
            }
            if (idxKaonDaughter < 0 && (std::abs(daughter.pdgCode()) == kKPlus)) {
              idxKaonDaughter = daughter.globalIndex();
            }
          }
          if (idxPionDaughter >= 0 && idxCascDaughter >= 0) {
            decayChannel = o2::aod::hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi; // OmegaC -> Omega + Pi
          } else if (idxKaonDaughter >= 0 && idxCascDaughter >= 0) {
            decayChannel = o2::aod::hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK; // OmegaC -> Omega + K
          } else {
            decayChannel = -1; // LOG(warning) << "Decay channel not recognized!";
          }
          if (decayChannel != -1) {
            int const idxDaughter = (decayChannel == o2::aod::hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi) ? idxPionDaughter : idxKaonDaughter;
            auto particle = mcParticles.rawIteratorAt(idxDaughter);
            origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, false, &idxBhadMothers);

            const auto& cascDaughter = mcParticles.iteratorAt(idxCascDaughter);
            const auto& mcColl = mcParticle.mcCollision();
            std::array<double, 3> const primaryVertexPosGen = {mcColl.posX(), mcColl.posY(), mcColl.posZ()};
            std::array<double, 3> const secondaryVertexGen = {cascDaughter.vx(), cascDaughter.vy(), cascDaughter.vz()};
            float decayLengthCascGen = -1.;
            float decayLengthXYCascGen = -1.;
            if (cascDaughter.has_daughters()) {
              const auto& cascDecayDaughter = cascDaughter.daughters_as<aod::McParticles>().iteratorAt(0);
              std::array<double, 3> const tertiaryVertexGen = {cascDecayDaughter.vx(), cascDecayDaughter.vy(), cascDecayDaughter.vz()};
              decayLengthCascGen = RecoDecay::distance(tertiaryVertexGen, primaryVertexPosGen);
              decayLengthXYCascGen = RecoDecay::distanceXY(tertiaryVertexGen, primaryVertexPosGen);
            }
            const auto decayLengthGen = RecoDecay::distance(secondaryVertexGen, primaryVertexPosGen);
            const auto decayLengthXYGen = RecoDecay::distanceXY(secondaryVertexGen, primaryVertexPosGen);
            registry.fill(HIST("hDecayLengthScaledMc"), decayLengthGen * o2::constants::physics::MassOmegaC0 / mcParticle.mothers_first_as<aod::McParticles>().p() * 1e4);
            outputTableGen(
              mcParticle.px(),
              mcParticle.py(),
              mcParticle.pz(),
              mcParticle.pdgCode(),
              cascDaughter.px(),
              cascDaughter.py(),
              cascDaughter.pz(),
              cascDaughter.pdgCode(),
              decayLengthGen,
              decayLengthXYGen,
              decayLengthCascGen,
              decayLengthXYCascGen,
              origin,
              decayChannel);
            mapMcPartToGenTable[mcParticle.globalIndex()] = outputTableGen.lastIndex();
          }
        }
      }
    }
  }
  PROCESS_SWITCH(HfTreeCreatorOmegacSt, processMc, "Process MC", true);

  template <o2::hf_centrality::CentralityEstimator CentEstimator, typename TracksType, typename TColl>
  void fillTable(TColl const& collisions,
                 aod::AssignedTrackedCascades const& trackedCascades,
                 aod::TrackAssoc const& trackIndices,
                 std::optional<std::reference_wrapper<aod::McParticles const>> mcParticles = std::nullopt)
  {
    const auto matCorr = static_cast<o2::base::Propagator::MatCorrType>(materialCorrectionType.value);

    for (const auto& collision : collisions) {
      // Event selection & centrality using HF helper
      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      if (rejectionMask != 0) {
        continue; // reject collisions failing HF event selection
      }
      hfEvSel.fillHistograms(collision, rejectionMask, centrality);
      const auto bc = collision.template bc_as<aod::BCsWithTimestamps>();

      // Update CCDB-based magnetic field only when run changes to avoid redundant lookups
      if (runNumber != bc.runNumber()) {
        if (skimmedProcessing) {
          if (runNumber == 0) {
            zorroSummary.setObject(zorro.getZorroSummary());
          }
          zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), cfgTriggersOfInterest.value);
          zorro.populateHistRegistry(registry, bc.runNumber());
        }
        runNumber = bc.runNumber();
        auto timestamp = bc.timestamp();

        bool fieldLoaded = false;
        if (auto* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, timestamp)) {
          o2::base::Propagator::initFieldFromGRP(grpo);
          bz = grpo->getNominalL3Field();
          fieldLoaded = true;
        } else if (auto* grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpMagPath, timestamp)) {
          o2::base::Propagator::initFieldFromGRP(grpmag);
          bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
          fieldLoaded = true;
        } else {
          LOG(error) << "CCDB miss for GRP/field (" << grpMagPath << ", " << grpPath << ") at timestamp " << timestamp << "; skipping collision to avoid bad B field.";
        }
        if (!fieldLoaded) {
          continue; // skip this collision if CCDB does not provide field
        }
        df2.setBz(bz);
      }

      uint32_t toiMask = 0;
      if (skimmedProcessing) {
        bool const sel = zorro.isSelected(bc.globalBC());
        if (sel) {
          std::vector<bool> toivect = zorro.getTriggerOfInterestResults();
          for (size_t i{0}; i < toivect.size(); i++) {
            toiMask |= static_cast<uint32_t>(toivect[i]) << i;
          }
        }
      }

      const auto primaryVertex = getPrimaryVertex(collision);
      const std::array<double, 3> primaryVertexPos = {primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()};

      const auto collId = collision.globalIndex();
      auto groupedTrackIds = trackIndices.sliceBy(trackIndicesPerCollision, collId);
      auto groupedTrackedCascades = trackedCascades.sliceBy(assignedTrackedCascadesPerCollision, collId);

      o2::dataformats::DCA impactParameterCasc;
      for (const auto& trackedCascade : groupedTrackedCascades) {
        const auto trackCasc = trackedCascade.template track_as<TracksType>();
        int trackCascMotherId = -1;
        if constexpr (std::is_same<TracksType, TracksExtMc>::value) {
          if (trackCasc.has_mcParticle() && trackCasc.mcParticle().has_mothers()) {
            if (auto res = mapMcPartToGenTable.find(trackCasc.mcParticle().mothersIds()[0]); res != mapMcPartToGenTable.end()) {
              trackCascMotherId = res->second;
            }
          }
        }
        auto trackParCovCasc = getTrackParCov(trackCasc);
        if (bzOnly) {
          o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackParCovCasc, bz, 2.f, matCorr, &impactParameterCasc);
        } else {
          o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovCasc, 2.f, matCorr, &impactParameterCasc);
        }

        const auto& casc = trackedCascade.cascade();
        const auto& bachelor = casc.template bachelor_as<TracksType>();
        const auto& v0 = casc.v0();
        const auto& v0TrackPos = v0.template posTrack_as<TracksType>();
        const auto& v0TrackNeg = v0.template negTrack_as<TracksType>();

        if (!v0TrackPos.hasTPC() || !v0TrackNeg.hasTPC() || !bachelor.hasTPC() ||
            v0TrackPos.tpcNClsFindable() < minNoClsTrackedCascade ||
            v0TrackNeg.tpcNClsFindable() < minNoClsTrackedCascade ||
            bachelor.tpcNClsFindable() < minNoClsTrackedCascade) {
          continue;
        }

        const auto& v0TrackPr = trackCasc.sign() < 0 ? v0TrackPos : v0TrackNeg;
        const auto& v0TrackPi = trackCasc.sign() < 0 ? v0TrackNeg : v0TrackPos;

        // track propagation
        hCandidatesPrPi->Fill(SVFitting::BeforeFit);
        try {
          if (!df2.process(getTrackParCov(v0TrackPr), getTrackParCov(v0TrackPi))) {
            continue;
          }
        } catch (const std::runtime_error& error) {
          LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN for Pr-Pi cannot work, skipping the candidate.";
          hCandidatesPrPi->Fill(SVFitting::Fail);
          continue;
        }
        hCandidatesPrPi->Fill(SVFitting::FitOk);

        std::array<double, NDaughters> const massesV0Daughters{o2::constants::physics::MassProton, o2::constants::physics::MassPiMinus};
        std::array<std::array<float, 3>, NDaughters> momentaV0Daughters{};
        o2::track::TrackPar const trackParV0Pr = df2.getTrackParamAtPCA(0);
        trackParV0Pr.getPxPyPzGlo(momentaV0Daughters[0]);
        o2::track::TrackPar const trackParV0Pi = df2.getTrackParamAtPCA(1);
        trackParV0Pi.getPxPyPzGlo(momentaV0Daughters[1]);
        const auto massV0 = RecoDecay::m(momentaV0Daughters, massesV0Daughters);

        o2::track::TrackParCov const trackParCovV0 = df2.createParentTrackParCov(0);
        hCandidatesV0Pi->Fill(SVFitting::BeforeFit);
        try {
          if (!df2.process(trackParCovV0, getTrackParCov(bachelor))) {
            continue;
          }
        } catch (const std::runtime_error& error) {
          LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN for V0-bachelor cannot work, skipping the candidate.";
          hCandidatesV0Pi->Fill(SVFitting::Fail);
          continue;
        }
        hCandidatesV0Pi->Fill(SVFitting::FitOk);

        const auto& secondaryVertex = df2.getPCACandidate();
        const auto decayLengthCasc = RecoDecay::distance(secondaryVertex, primaryVertexPos);
        const auto decayLengthCascXY = RecoDecay::distanceXY(secondaryVertex, primaryVertexPos);
        o2::track::TrackPar const trackParV0 = df2.getTrackParamAtPCA(0);
        o2::track::TrackPar const trackParBachelor = df2.getTrackParamAtPCA(1);
        std::array<std::array<float, 3>, NDaughters> momentaCascDaughters{};
        trackParV0.getPxPyPzGlo(momentaCascDaughters[0]);
        trackParBachelor.getPxPyPzGlo(momentaCascDaughters[1]);
        o2::track::TrackParCov const trackParCovCascUntracked = df2.createParentTrackParCov(0);
        std::array<float, 3> pCasc{};
        trackParCovCascUntracked.getPxPyPzGlo(pCasc);
        const auto cpaCasc = RecoDecay::cpa(primaryVertexPos, df2.getPCACandidate(), pCasc);
        const auto cpaXYCasc = RecoDecay::cpaXY(primaryVertexPos, df2.getPCACandidate(), pCasc);

        std::array<double, NDaughters> const massesXiDaughters = {o2::constants::physics::MassLambda0, o2::constants::physics::MassPiPlus};
        const auto massXi = RecoDecay::m(momentaCascDaughters, massesXiDaughters);
        std::array<double, NDaughters> const massesOmegaDaughters = {o2::constants::physics::MassLambda0, o2::constants::physics::MassKPlus};
        const auto massOmega = RecoDecay::m(momentaCascDaughters, massesOmegaDaughters);

        registry.fill(HIST("hDca"), std::sqrt(impactParameterCasc.getR2()));
        registry.fill(HIST("hDcaXY"), impactParameterCasc.getY());
        registry.fill(HIST("hDcaXYVsPt"), trackParCovCasc.getPt(), impactParameterCasc.getY());
        registry.fill(HIST("hDcaZ"), impactParameterCasc.getZ());
        registry.fill(HIST("hDcaZVsPt"), trackParCovCasc.getPt(), impactParameterCasc.getZ());
        registry.fill(HIST("hDcaVsPt"), impactParameterCasc.getY(), trackCasc.pt());
        registry.fill(HIST("hDcaVsR"), impactParameterCasc.getY(), RecoDecay::sqrtSumOfSquares(trackCasc.x(), trackCasc.y()));
        registry.fill(HIST("hPtVsMassOmega"), trackCasc.pt(), massOmega);

        if ((std::abs(massOmega - o2::constants::physics::MassOmegaMinus) < massWindowTrackedOmega) ||
            (std::abs(massXi - o2::constants::physics::MassXiMinus) < massWindowTrackedXi)) {
          // combine configurable TPC and TOF nsigma requirements for bachelor and V0 daughters
          // PID: decide based on `useTofPid` switch
          const bool tpcBachelor = (std::abs(bachelor.tpcNSigmaKa()) < nSigmaBachelorMax) || (std::abs(bachelor.tpcNSigmaPi()) < nSigmaBachelorMax);
          const float tofBachelorPiAbs = std::abs(bachelor.tofNSigmaPi());
          const float tofBachelorKaAbs = std::abs(bachelor.tofNSigmaKa());
          const bool tofBachelorPiPass = (tofBachelorPiAbs > TofWoSignalRange[0] && tofBachelorPiAbs < TofWoSignalRange[1]) || (tofBachelorPiAbs < nSigmaAbsTofBachelorPiMax);
          const bool tofBachelorKaPass = (tofBachelorKaAbs > TofWoSignalRange[0] && tofBachelorKaAbs < TofWoSignalRange[1]) || (tofBachelorKaAbs < nSigmaAbsTofBachelorKaMax);
          const bool tofBachelorPass = tofBachelorPiPass || tofBachelorKaPass;
          const bool bachelorPass = useTofPid.value ? (tpcBachelor && tofBachelorPass) : tpcBachelor;

          const bool tpcV0Pr = (std::abs(v0TrackPr.tpcNSigmaPr()) < nSigmaV0PrMax);
          const bool tpcV0Pi = (std::abs(v0TrackPi.tpcNSigmaPi()) < nSigmaV0PiMax);
          const float tofV0PrAbs = std::abs(v0TrackPr.tofNSigmaPr());
          const float tofV0PiAbs = std::abs(v0TrackPi.tofNSigmaPi());
          const bool tofV0PrPass = (tofV0PrAbs > TofWoSignalRange[0] && tofV0PrAbs < TofWoSignalRange[1]) || (tofV0PrAbs < nSigmaAbsTofV0PrMax);
          const bool tofV0PiPass = (tofV0PiAbs > TofWoSignalRange[0] && tofV0PiAbs < TofWoSignalRange[1]) || (tofV0PiAbs < nSigmaAbsTofV0PiMax);
          const bool v0PrPass = useTofPid.value ? (tpcV0Pr && tofV0PrPass) : tpcV0Pr;
          const bool v0PiPass = useTofPid.value ? (tpcV0Pi && tofV0PiPass) : tpcV0Pi;

          if (bachelorPass && v0PrPass && v0PiPass) {

            std::array<double, NDaughters> const massesOmegacToOmegaPi{o2::constants::physics::MassOmegaMinus, o2::constants::physics::MassPiPlus};
            std::array<double, NDaughters> const massOmegaK{o2::constants::physics::MassOmegaMinus, o2::constants::physics::MassKPlus};
            std::array<double, NDaughters> const massesXicToXiPi{o2::constants::physics::MassXiMinus, o2::constants::physics::MassPiPlus};

            std::array<std::array<float, 3>, NDaughters> momenta{};

            auto trackParCovPr = getTrackParCov(v0TrackPr);
            auto trackParCovKa = getTrackParCov(v0TrackPi);
            auto trackParCovPi = getTrackParCov(bachelor);
            o2::dataformats::DCA impactParameterPr;
            o2::dataformats::DCA impactParameterKa;
            o2::dataformats::DCA impactParameterPi;
            if (bzOnly) {
              o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackParCovPr, bz, 2.f, matCorr, &impactParameterPr);
              o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackParCovKa, bz, 2.f, matCorr, &impactParameterKa);
              o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackParCovPi, bz, 2.f, matCorr, &impactParameterPi);
            } else {
              o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovPr, 2.f, matCorr, &impactParameterPr);
              o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovKa, 2.f, matCorr, &impactParameterKa);
              o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovPi, 2.f, matCorr, &impactParameterPi);
            }
            for (const auto& trackId : groupedTrackIds) {
              const auto track = trackId.template track_as<TracksType>();
              if (track.globalIndex() == v0TrackPr.globalIndex() ||
                  track.globalIndex() == v0TrackPi.globalIndex() ||
                  track.globalIndex() == bachelor.globalIndex()) {
                continue;
              }
              if ((track.itsNCls() >= ItsNClsMin) &&
                  (track.tpcNClsFound() >= minNoClsTrackedPionOrKaon) &&
                  (track.tpcNClsCrossedRows() >= minNoClsTrackedPionOrKaon) &&
                  (track.tpcNClsCrossedRows() >= TpcNclsFindableFraction * track.tpcNClsFindable()) &&
                  (track.tpcChi2NCl() <= TpcChi2NclMax) &&
                  (track.itsChi2NCl() <= ItsChi2NclMax)) {
                // TPC and TOF PID for pion/kaon from chambaryon
                const bool passTPCpid = (std::abs(track.tpcNSigmaPi()) < nSigmaPionMax) || (std::abs(track.tpcNSigmaKa()) < nSigmaKaonMax);
                const float tofPiAbs = std::abs(track.tofNSigmaPi());
                const float tofKaAbs = std::abs(track.tofNSigmaKa());
                const bool tofPiPass = (tofPiAbs > TofWoSignalRange[0] && tofPiAbs < TofWoSignalRange[1]) || (tofPiAbs < nSigmaAbsTofTrackPiMax);
                const bool tofKaPass = (tofKaAbs > TofWoSignalRange[0] && tofKaAbs < TofWoSignalRange[1]) || (tofKaAbs < nSigmaAbsTofTrackKaMax);
                const bool passTOFpid = tofPiPass || tofKaPass;
                if (useTofPid.value) {
                  if (!(passTPCpid && passTOFpid)) {
                    continue;
                  }
                } else {
                  if (!passTPCpid) {
                    continue;
                  }
                }
                LOGF(debug, "  .. combining with pion/kaon candidate %d", track.globalIndex());
                int trackMotherId = -1;
                if constexpr (std::is_same<TracksType, TracksExtMc>::value) {
                  if (track.has_mcParticle() && track.mcParticle().has_mothers()) {
                    if (auto res = mapMcPartToGenTable.find(track.mcParticle().mothersIds()[0]); res != mapMcPartToGenTable.end()) {
                      trackMotherId = res->second;
                    }
                  }
                }
                auto trackParCovCasc = getTrackParCov(trackCasc);
                auto trackParCovPionOrKaon = getTrackParCov(track);
                o2::dataformats::DCA impactParameterPion;
                if (bzOnly) {
                  o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackParCovPionOrKaon, bz, 2.f, matCorr, &impactParameterPion);
                } else {
                  o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovPionOrKaon, 2.f, matCorr, &impactParameterPion);
                }
                const float impactProduct = impactParameterCasc.getY() * impactParameterPion.getY();
                hCandidatesCascPiOrK->Fill(SVFitting::BeforeFit);
                try {
                  auto decayLengthUntracked = -1.;
                  auto decayLengthXYUntracked = -1.;
                  if (df2.process(trackParCovCascUntracked, trackParCovPionOrKaon)) {
                    const auto& secondaryVertexUntracked = df2.getPCACandidate();
                    decayLengthUntracked = RecoDecay::distance(secondaryVertexUntracked, primaryVertexPos);
                    decayLengthXYUntracked = RecoDecay::distanceXY(secondaryVertexUntracked, primaryVertexPos);
                  }
                  if (df2.process(trackParCovCasc, trackParCovPionOrKaon)) {
                    const auto& secondaryVertex = df2.getPCACandidate();
                    const auto decayLength = RecoDecay::distance(secondaryVertex, primaryVertexPos);
                    const auto decayLengthXY = RecoDecay::distanceXY(secondaryVertex, primaryVertexPos);
                    const auto chi2TopCharmedBaryon = df2.getChi2AtPCACandidate();
                    std::array<float, 3> pCharmedBaryon{};
                    df2.createParentTrackParCov().getPxPyPzGlo(pCharmedBaryon);
                    const auto cpaCharmedBaryon = RecoDecay::cpa(primaryVertexPos, df2.getPCACandidate(), pCharmedBaryon);
                    const auto cpaXYCharmedBaryon = RecoDecay::cpaXY(primaryVertexPos, df2.getPCACandidate(), pCharmedBaryon);

                    df2.getTrackParamAtPCA(0).getPxPyPzGlo(momenta[0]);
                    df2.getTrackParamAtPCA(1).getPxPyPzGlo(momenta[1]);
                    const auto ptCascade = std::hypot(momenta[0][0], momenta[0][1]);
                    const auto ptPionOrKaon = std::hypot(momenta[1][0], momenta[1][1]);
                    const auto ptCharmedBaryon = RecoDecay::pt(momenta[0], momenta[1]);
                    const auto massOmegaPion = RecoDecay::m(momenta, massesOmegacToOmegaPi);
                    const auto massOmegaKaon = RecoDecay::m(momenta, massOmegaK);
                    const auto massXiPion = RecoDecay::m(momenta, massesXicToXiPi);
                    registry.fill(HIST("hMassOmegaPi"), massOmegaPion);
                    registry.fill(HIST("hMassOmegaPiVsPt"), massOmegaPion, ptCharmedBaryon);
                    registry.fill(HIST("hMassOmegaK"), massOmegaKaon);
                    registry.fill(HIST("hMassOmegaKVsPt"), massOmegaKaon, ptCharmedBaryon);
                    registry.fill(HIST("hMassXiPiVsPt"), massXiPion, ptCharmedBaryon);
                    const bool massOmegacToOmegaPiPass = std::abs(massOmegaPion - o2::constants::physics::MassOmegaC0) < massWindowOmegaC;
                    const bool massOmegacToOmegaKPass = std::abs(massOmegaKaon - o2::constants::physics::MassOmegaC0) < massWindowOmegaC;
                    const bool massOmegacToXiPiPass = std::abs(massXiPion - o2::constants::physics::MassOmegaC0) < massWindowOmegaC;
                    const bool massXicToOmegaKaPass = std::abs(massOmegaKaon - o2::constants::physics::MassXiC0) < massWindowXiC;
                    const bool massXicToXiPiPass = std::abs(massXiPion - o2::constants::physics::MassXiC0) < massWindowXiC;
                    const bool anyChannelPass = massOmegacToOmegaKPass || massOmegacToOmegaPiPass || massXicToXiPiPass || massXicToOmegaKaPass || massOmegacToXiPiPass;
                    bool passSelectedChannel = true;
                    switch (selectedChannel.value) {
                      case AnyChannel:
                        passSelectedChannel = anyChannelPass;
                        break;
                      case XicToXiPi:
                        passSelectedChannel = massXicToXiPiPass;
                        break;
                      case XicToOmegaKa:
                        passSelectedChannel = massXicToOmegaKaPass;
                        break;
                      case OmegacToOmegaPi:
                        passSelectedChannel = massOmegacToOmegaPiPass;
                        break;
                      case OmegacToOmegaKa:
                        passSelectedChannel = massOmegacToOmegaKPass;
                        break;
                      case OmegacToXiPi:
                        passSelectedChannel = massOmegacToXiPiPass;
                        break;
                      default:
                        LOG(warning) << "Unexpected selectedChannel value: " << selectedChannel.value << ". Rejecting candidate.";
                        passSelectedChannel = false;
                        break;
                    }
                    if (!passSelectedChannel) {
                      continue;
                    }
                    //--- do the MC Rec match
                    if (mcParticles) {
                      auto arrayDaughters = std::array{
                        trackId.template track_as<TracksExtMc>(), // bachelor <- charm baryon
                        casc.template bachelor_as<TracksExtMc>(), // bachelor <- cascade
                        v0.template posTrack_as<TracksExtMc>(),   // p <- lambda
                        v0.template negTrack_as<TracksExtMc>()};  // pi <- lambda

                      auto arrayDaughtersCasc = std::array{
                        casc.template bachelor_as<TracksExtMc>(), // bachelor <- cascade
                        v0.template posTrack_as<TracksExtMc>(),   // p <- lambda
                        v0.template negTrack_as<TracksExtMc>()};  // pi <- lambda
                      auto arrayDaughtersV0 = std::array{
                        v0.template posTrack_as<TracksExtMc>(),  // p <- lambda
                        v0.template negTrack_as<TracksExtMc>()}; // pi <- lambda

                      if (decayChannel == o2::aod::hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaPi) {
                        // Match Omegac0 → Omega- + Pi+
                        indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, false>(mcParticles->get(), arrayDaughters, o2::constants::physics::kOmegaC0,
                                                                                               std::array{+kPiPlus, +kKMinus, +kProton, +kPiMinus}, true, &sign, 3, &nPiToMuOmegac0, &nKaToPiOmegac0);
                        indexRecCharmBaryon = indexRec;
                        if (indexRec > -1) {
                          // Omega- → K pi p (Cascade match)
                          indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(mcParticles->get(), arrayDaughtersCasc, +kOmegaMinus, std::array{+kKMinus, +kProton, +kPiMinus}, true, &signCasc, 2, &nPiToMuCasc, &nKaToPiCasc);
                          if (indexRec > -1) {
                            // Lambda → p pi (Lambda match)
                            indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(mcParticles->get(), arrayDaughtersV0, +kLambda0, std::array{+kProton, +kPiMinus}, true, &signV0, 1, &nPiToMuV0);
                            if (indexRec > -1) {
                              isMatched = true;
                            }
                          }
                        }
                      } else if (decayChannel == o2::aod::hf_cand_casc_lf::DecayType2Prong::OmegaczeroToOmegaK) {
                        // Match Omegac0 → Omega- + K+
                        indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, false>(mcParticles->get(), arrayDaughters, o2::constants::physics::kOmegaC0,
                                                                                               std::array{+kKPlus, +kKMinus, +kProton, +kPiMinus}, true, &sign, 3, &nPiToMuOmegac0, &nKaToPiOmegac0);
                        indexRecCharmBaryon = indexRec;
                        if (indexRec > -1) {
                          // Omega- → K pi p (Cascade match)
                          indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(mcParticles->get(), arrayDaughtersCasc, +kOmegaMinus, std::array{+kKMinus, +kProton, +kPiMinus}, true, &signCasc, 2, &nPiToMuCasc, &nKaToPiCasc);
                          if (indexRec > -1) {
                            // Lambda → p pi (Lambda match)
                            indexRec = RecoDecay::getMatchedMCRec<false, true, false, true, true>(mcParticles->get(), arrayDaughtersV0, +kLambda0, std::array{+kProton, +kPiMinus}, true, &signV0, 1, &nPiToMuV0);
                            if (indexRec > -1) {
                              isMatched = true;
                            }
                          }
                        }
                      }
                      if (isMatched && indexRecCharmBaryon > -1) {
                        auto particle = mcParticles->get().rawIteratorAt(indexRecCharmBaryon);
                        origin = RecoDecay::getCharmHadronOrigin(mcParticles->get(), particle, false, &idxBhadMothers);
                      }
                    }
                    // apply configurable quality cuts (CPA, chi2, impact parameters)
                    if (cpaCharmedBaryon < cpaCharmedBaryonMin)
                      continue;
                    if (cpaXYCharmedBaryon < cpaXYCharmedBaryonMin)
                      continue;
                    if (cpaCasc < cpaCascMin)
                      continue;
                    if (cpaXYCasc < cpaXYCascMin)
                      continue;
                    if (chi2TopCharmedBaryon > chi2TopCharmedBaryonMax)
                      continue;
                    if (trackedCascade.topologyChi2() > chi2TopCascMax)
                      continue;
                    if (std::abs(impactParameterCasc.getY()) > absImpactParameterCascYMax || std::abs(impactParameterCasc.getZ()) > absDcaCascZMax)
                      continue;
                    if (std::abs(impactParameterPion.getY()) > absImpactParameterCharmBachelorYMax || std::abs(impactParameterPion.getZ()) > absDcaCharmBachelorZMax)
                      continue;
                    if (ptCascade < ptCascadeMin)
                      continue;
                    if (ptPionOrKaon < ptPionOrKaonMin)
                      continue;
                    if (ptCharmedBaryon < ptCharmedBaryonMin)
                      continue;
                    if (ptCharmedBaryon > ptCharmedBaryonMax)
                      continue;
                    if (impactProduct < impProductMin || impactProduct > impProductMax)
                      continue;
                    registry.fill(HIST("hDecayLength"), decayLength * 1e4);
                    registry.fill(HIST("hDecayLengthScaled"), decayLength * o2::constants::physics::MassOmegaC0 / RecoDecay::p(momenta[0], momenta[1]) * 1e4);
                    outputTable(centrality,
                                massOmega,
                                massXi,
                                massV0,
                                track.tpcNSigmaPi(),
                                track.tofNSigmaPi(),
                                track.tpcNSigmaKa(),
                                track.tofNSigmaKa(),
                                v0TrackPr.tpcNSigmaPr(),
                                v0TrackPr.tofNSigmaPr(),
                                v0TrackPi.tpcNSigmaPi(),
                                v0TrackPi.tofNSigmaPi(),
                                bachelor.tpcNSigmaPi(),
                                bachelor.tofNSigmaPi(),
                                bachelor.tpcNSigmaKa(),
                                bachelor.tofNSigmaKa(),
                                momenta[0][0], // cascade momentum
                                momenta[0][1],
                                momenta[0][2],
                                static_cast<bool>(trackCasc.sign() > 0),
                                momenta[1][0], // pion/kaon momentum
                                momenta[1][1],
                                momenta[1][2],
                                static_cast<bool>(track.sign() > 0),
                                static_cast<uint8_t>(track.itsClusterMap()),
                                cpaCharmedBaryon,
                                cpaXYCharmedBaryon,
                                cpaCasc,
                                cpaXYCasc,
                                impactParameterCasc.getY(),
                                std::sqrt(impactParameterCasc.getSigmaY2()),
                                impactParameterCasc.getZ(),
                                std::sqrt(impactParameterCasc.getSigmaZ2()),
                                impactParameterPion.getY(),
                                std::sqrt(impactParameterPion.getSigmaY2()),
                                impactParameterPion.getZ(),
                                std::sqrt(impactParameterPion.getSigmaZ2()),
                                impactParameterPr.getY(),
                                impactParameterPr.getZ(),
                                impactParameterKa.getY(),
                                impactParameterKa.getZ(),
                                impactParameterPi.getY(),
                                impactParameterPi.getZ(),
                                chi2TopCharmedBaryon,
                                trackedCascade.topologyChi2(),
                                decayLength,
                                decayLengthXY,
                                decayLengthUntracked,
                                decayLengthXYUntracked,
                                decayLengthCasc,
                                decayLengthCascXY,
                                trackCascMotherId,
                                trackMotherId,
                                origin,
                                toiMask);
                  } else {
                    continue;
                  }
                } catch (const std::runtime_error& error) {
                  LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN for Casc-Pi/K cannot work, skipping the candidate.";
                  hCandidatesCascPiOrK->Fill(SVFitting::Fail);
                  continue;
                }
                hCandidatesCascPiOrK->Fill(SVFitting::FitOk);
              }
            }
          }
        }
      }
    }
  }

  // No centrality selection
  void processDataNocent(Collisions const& collisions,
                         soa::SmallGroups<aod::AssignedTrackedCascades> const& trackedCascades,
                         aod::TrackAssoc const& trackIndices,
                         aod::Cascades const&,
                         aod::V0s const&,
                         TracksExt const&,
                         aod::BCsWithTimestamps const&)
  {
    fillTable<CentralityEstimator::None, TracksExt>(collisions, trackedCascades, trackIndices);
  }
  PROCESS_SWITCH(HfTreeCreatorOmegacSt, processDataNocent, "Process data (No centrality)", true);

  // FT0C centrality selection
  void processDataFT0C(CollisionsFT0C const& collisions,
                       soa::SmallGroups<aod::AssignedTrackedCascades> const& trackedCascades,
                       aod::TrackAssoc const& trackIndices,
                       aod::Cascades const&,
                       aod::V0s const&,
                       TracksExt const&,
                       aod::BCsWithTimestamps const&)
  {
    fillTable<CentralityEstimator::FT0C, TracksExt>(collisions, trackedCascades, trackIndices);
  }
  PROCESS_SWITCH(HfTreeCreatorOmegacSt, processDataFT0C, "Process data (FT0C centrality)", false);

  // FT0M centrality selection
  void processDataFT0M(CollisionsFT0M const& collisions,
                       soa::SmallGroups<aod::AssignedTrackedCascades> const& trackedCascades,
                       aod::TrackAssoc const& trackIndices,
                       aod::Cascades const&,
                       aod::V0s const&,
                       TracksExt const&,
                       aod::BCsWithTimestamps const&)
  {
    fillTable<CentralityEstimator::FT0M, TracksExt>(collisions, trackedCascades, trackIndices);
  }
  PROCESS_SWITCH(HfTreeCreatorOmegacSt, processDataFT0M, "Process data (FT0M centrality)", false);

  // MC reco: No centrality selection
  void processMcRecNocent(Collisions const& collisions,
                          aod::McCollisions const&,
                          soa::SmallGroups<aod::AssignedTrackedCascades> const& trackedCascades,
                          aod::TrackAssoc const& trackIndices,
                          aod::Cascades const&,
                          aod::V0s const&,
                          TracksExtMc const&,
                          aod::McParticles const& mcParticles,
                          aod::BCsWithTimestamps const&)
  {
    fillTable<CentralityEstimator::None, TracksExtMc>(collisions, trackedCascades, trackIndices, mcParticles);
  }
  PROCESS_SWITCH(HfTreeCreatorOmegacSt, processMcRecNocent, "Process MC reco (No centrality)", false);

  // MC reco: FT0C centrality selection
  void processMcRecFT0C(CollisionsFT0C const& collisions,
                        aod::McCollisions const&,
                        soa::SmallGroups<aod::AssignedTrackedCascades> const& trackedCascades,
                        aod::TrackAssoc const& trackIndices,
                        aod::Cascades const&,
                        aod::V0s const&,
                        TracksExtMc const&,
                        aod::McParticles const& mcParticles,
                        aod::BCsWithTimestamps const&)
  {
    fillTable<CentralityEstimator::FT0C, TracksExtMc>(collisions, trackedCascades, trackIndices, mcParticles);
  }
  PROCESS_SWITCH(HfTreeCreatorOmegacSt, processMcRecFT0C, "Process MC reco (FT0C centrality)", false);

  // MC reco: FT0M centrality selection
  void processMcRecFT0M(CollisionsFT0M const& collisions,
                        aod::McCollisions const&,
                        soa::SmallGroups<aod::AssignedTrackedCascades> const& trackedCascades,
                        aod::TrackAssoc const& trackIndices,
                        aod::Cascades const&,
                        aod::V0s const&,
                        TracksExtMc const&,
                        aod::McParticles const& mcParticles,
                        aod::BCsWithTimestamps const&)
  {
    fillTable<CentralityEstimator::FT0M, TracksExtMc>(collisions, trackedCascades, trackIndices, mcParticles);
  }
  PROCESS_SWITCH(HfTreeCreatorOmegacSt, processMcRecFT0M, "Process MC reco (FT0M centrality)", false);

  void processMcGen(aod::Collision const& collision,
                    aod::McCollisions const&,
                    soa::SmallGroups<aod::AssignedTrackedCascades> const& trackedCascades,
                    aod::Cascades const&,
                    aod::V0s const&,
                    TracksExtMc const& tracks,
                    aod::McParticles const&,
                    aod::BCsWithTimestamps const&)
  {
    const auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    if (runNumber != bc.runNumber()) {
      runNumber = bc.runNumber();
      auto timestamp = bc.timestamp();

      if (auto* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, timestamp)) {
        o2::base::Propagator::initFieldFromGRP(grpo);
      } else if (auto* grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpMagPath, timestamp)) {
        o2::base::Propagator::initFieldFromGRP(grpmag);
      } else {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpMagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << timestamp;
      }
    }

    const auto matCorr = static_cast<o2::base::Propagator::MatCorrType>(materialCorrectionType.value);
    const auto primaryVertex = getPrimaryVertex(collision);
    o2::dataformats::DCA impactParameterCasc;
    for (const auto& trackedCascade : trackedCascades) {
      const auto trackCasc = trackedCascade.template track_as<TracksExtMc>();
      auto trackParCovCasc = getTrackParCov(trackCasc);
      if (bzOnly) {
        o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackParCovCasc, bz, 2.f, matCorr, &impactParameterCasc);
      } else {
        o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovCasc, 2.f, matCorr, &impactParameterCasc);
      }

      const auto& casc = trackedCascade.cascade();
      const auto& bachelor = casc.template bachelor_as<TracksExtMc>();
      const auto& v0 = casc.v0();
      const auto& v0TrackPos = v0.template posTrack_as<TracksExtMc>();
      const auto& v0TrackNeg = v0.template negTrack_as<TracksExtMc>();

      if (!v0TrackPos.has_mcParticle() || !v0TrackNeg.has_mcParticle() || !bachelor.has_mcParticle()) {
        continue;
      }

      LOGF(debug, "v0TrackPos (id: %d, pdg: %d) has mother %d", v0TrackPos.mcParticleId(),
           v0TrackPos.mcParticle().pdgCode(), v0TrackPos.mcParticle().has_mothers() ? v0TrackPos.mcParticle().mothersIds()[0] : -1);
      LOGF(debug, "v0TrackNeg (id: %d, pdg: %d) has mother %d", v0TrackNeg.mcParticleId(),
           v0TrackNeg.mcParticle().pdgCode(), v0TrackNeg.mcParticle().has_mothers() ? v0TrackNeg.mcParticle().mothersIds()[0] : -1);

      LOG(debug) << "bachelor with PDG code: " << bachelor.mcParticle().pdgCode();
      if (v0TrackPos.mcParticle().has_mothers() && v0TrackNeg.mcParticle().has_mothers() &&
          v0TrackPos.mcParticle().mothersIds()[0] == v0TrackNeg.mcParticle().mothersIds()[0]) {
        const auto v0part = v0TrackPos.mcParticle().mothers_first_as<aod::McParticles>();
        LOG(debug) << "v0 with PDG code: " << v0part.pdgCode();
        if (v0part.has_mothers() && bachelor.mcParticle().has_mothers() &&
            v0part.mothersIds()[0] == bachelor.mcParticle().mothersIds()[0]) {
          const auto mother = v0part.mothers_as<aod::McParticles>()[0];
          const auto pdgCode = mother.pdgCode();
          LOG(debug) << "cascade with PDG code: " << pdgCode;
          if (std::abs(pdgCode) == kOmegaMinus) {
            LOG(debug) << "found Omega, looking for pions";
            std::array<double, NDaughters> const masses{o2::constants::physics::MassOmegaMinus, o2::constants::physics::MassPiPlus};
            std::array<std::array<float, 3>, NDaughters> momenta{};
            std::array<double, 3> const primaryVertexPos = {primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()};
            const auto& mcColl = mother.mcCollision();
            std::array<double, 3> const primaryVertexPosGen = {mcColl.posX(), mcColl.posY(), mcColl.posZ()};

            for (const auto& track : tracks) {
              if (!track.has_mcParticle()) {
                continue;
              }
              const auto mcpart = track.mcParticle();
              if (mcpart.pdgCode() == (pdgCode > 0 ? kPiPlus : -kPiPlus)) {
                LOGF(debug, "combining Omega with pion %d", track.globalIndex());
                auto trackParCovPion = getTrackParCov(track);
                o2::dataformats::DCA impactParameterPion;
                if (bzOnly) {
                  o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackParCovPion, bz, 2.f, matCorr, &impactParameterPion);
                } else {
                  o2::base::Propagator::Instance()->propagateToDCABxByBz(primaryVertex, trackParCovPion, 2.f, matCorr, &impactParameterPion);
                }

                trackParCovCasc.getPxPyPzGlo(momenta[0]);
                trackParCovPion.getPxPyPzGlo(momenta[1]);
                registry.fill(HIST("hDeltaPtVsPt"), mcpart.pt(), (trackParCovPion.getPt() - mcpart.pt()) / mcpart.pt());
                registry.fill(HIST("hMassOmegacId"), RecoDecay::m(momenta, masses));

                hCandidatesCascPiOrK->Fill(SVFitting::BeforeFit);
                try {
                  if (df2.process(trackParCovCasc, trackParCovPion) != 0) {
                    const auto& secondaryVertex = df2.getPCACandidate();
                    const auto decayLength = RecoDecay::distance(secondaryVertex, primaryVertexPos);
                    if (mother.has_mothers()) {
                      const auto& cand = mother.template mothers_first_as<aod::McParticles>();
                      if (std::abs(cand.pdgCode()) == constants::physics::Pdg::kOmegaC0 && mcpart.has_mothers()) {
                        if (mcpart.mothersIds()[0] == cand.globalIndex()) {
                          registry.fill(HIST("hDecayLengthId"), decayLength * 1e4);
                          registry.fill(HIST("hDecayLengthScaledId"), decayLength * o2::constants::physics::MassOmegaC0 / RecoDecay::p(momenta[0], momenta[1]) * 1e4);

                          std::array<double, 3> const secondaryVertexGen = {mother.vx(), mother.vy(), mother.vz()};
                          const auto decayLengthGen = RecoDecay::distance(secondaryVertexGen, primaryVertexPosGen);
                          registry.fill(HIST("hDecayLengthGen"), decayLengthGen * 1e4);
                          registry.fill(HIST("hDecayLengthScaledGen"), decayLengthGen * o2::constants::physics::MassOmegaC0 / RecoDecay::p(momenta[0], momenta[1]) * 1e4);

                          registry.fill(HIST("hDeltaDecayLength"), (decayLength - decayLengthGen) * 1e4);
                        }
                      }
                    }
                    hCandidatesCascPiOrK->Fill(SVFitting::FitOk);
                  }
                } catch (const std::runtime_error& error) {
                  LOG(info) << "Run time error found: " << error.what() << ". DCAFitterN for Casc-Pi cannot work, skipping the candidate.";
                  hCandidatesCascPiOrK->Fill(SVFitting::Fail);
                  continue;
                }

                // MC-based mass
                momenta[0] = mother.pVector();
                momenta[1] = mcpart.pVector();
                registry.fill(HIST("hMassOmegacGen"), RecoDecay::m(momenta, masses));
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(HfTreeCreatorOmegacSt, processMcGen, "Process using MC information", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTreeCreatorOmegacSt>(cfgc)};
}
