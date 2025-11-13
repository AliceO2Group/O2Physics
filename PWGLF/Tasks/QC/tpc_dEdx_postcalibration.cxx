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
///
/// \author Alberto Caliva (alberto.caliva@cern.ch)
/// \since June 27, 2023

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

using namespace o2;
using namespace o2::framework;

using PIDTracks = soa::Join<
  aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTOFbeta,
  aod::pidTOFmass, aod::TrackSelection, aod::TrackSelectionExtension,
  aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
  aod::pidTPCFullTr, aod::pidTPCFullHe, aod::pidTOFFullPi, aod::pidTOFFullKa,
  aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe>;

using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels>;

struct tpc_dEdx_postcalibration {

  // dE/dx for all charged particles
  HistogramRegistry registryCh{
    "registryCh",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  // dE/dx and nsigma_{TPC} for different hadron species
  HistogramRegistry registryPi{
    "registryPi",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  HistogramRegistry registryKa{
    "registryKa",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  HistogramRegistry registryPr{
    "registryPr",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  HistogramRegistry registryDe{
    "registryDe",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  HistogramRegistry registryTr{
    "registryTr",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};
  HistogramRegistry registryHe{
    "registryHe",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  // External Parameters
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 70.0f,
                                      "min number of found TPC clusters"};
  Configurable<float> minNCrossedRowsTPC{
    "minNCrossedRowsTPC", 70.0f, "min number of found TPC crossed rows"};
  Configurable<float> minNClsTPCdEdx{
    "minNClsTPCdEdx", 50.0f, "min number of TPC clusters for PID"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f,
                                 "max chi2 per cluster TPC"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f,
                                 "max chi2 per cluster ITS"};
  Configurable<float> etaMin{"etaMin", -0.8f, "etaMin"};
  Configurable<float> etaMax{"etaMax", +0.8f, "etaMax"};
  Configurable<float> v0cospaMin{"v0cospaMin", 0.998, "Minimum V0 CosPA"};
  Configurable<float> minimumV0Radius{"minimumV0Radius", 0.5f,
                                      "Minimum V0 Radius"};
  Configurable<float> maximumV0Radius{"maximumV0Radius", 100.0f,
                                      "Maximum V0 Radius"};
  Configurable<float> dcaV0DaughtersMax{"dcaV0DaughtersMax", 0.5f,
                                        "Maximum DCA Daughters"};
  Configurable<float> nsigmaTOFmax{"nsigmaTOFmax", 3.0f, "Maximum nsigma TOF"};
  Configurable<float> minMassK0s{"minMassK0s", 0.4f, "Minimum Mass K0s"};
  Configurable<float> maxMassK0s{"maxMassK0s", 0.6f, "Maximum Mass K0s"};
  Configurable<float> minMassLambda{"minMassLambda", 1.1f,
                                    "Minimum Mass Lambda"};
  Configurable<float> maxMassLambda{"maxMassLambda", 1.2f,
                                    "Maximum Mass Lambda"};
  Configurable<float> minReqClusterITS{
    "minReqClusterITS", 4.0f, "min number of clusters required in ITS"};
  Configurable<float> maxDCAxy{"maxDCAxy", 0.1f, "maxDCAxy"};
  Configurable<float> maxDCAz{"maxDCAz", 0.1f, "maxDCAz"};
  Configurable<bool> eventSelection{"eventSelection", true, "event selection"};
  Configurable<bool> useTOFpi{"useTOFpi", true, "use TOF for pion ID"};
  Configurable<bool> useTOFpr{"useTOFpr", true, "use TOF for proton ID"};
  Configurable<bool> doContaminations{"doContaminations", false, "Flag to produce the plots for the contaminations"};
  Configurable<bool> addTOF{"addTOF", false, "Flag to produce the TOF plots"};
  Configurable<bool> usePt{"usePt", false, "Flag to use PT instead of the TPCInnerParam"};
  ConfigurableAxis pBins{"pBins", {VARIABLE_WIDTH, -10.0, -9.660508789898131, -9.332543007969905, -9.01571137605957, -8.709635899560805, -8.413951416451948, -8.128305161640986, -7.852356346100718, -7.585775750291836, -7.328245331389037, -7.07945784384138, -6.839116472814292, -6.6069344800759575, -6.3826348619054825, -6.165950018614822, -5.956621435290104, -5.754399373371567, -5.559042572704037, -5.370317963702527, -5.188000389289609, -5.011872336272719, -4.841723675840994, -4.677351412871981, -4.518559443749222, -4.365158322401657, -4.216965034285822, -4.073802778041126, -3.9355007545577725, -3.801893963205613, -3.672823004980846, -3.5481338923357533, -3.427677865464501, -3.311311214825911, -3.198895109691397, -3.0902954325135887, -2.9853826189179604, -2.8840315031266055, -2.7861211686297693, -2.691534803926914, -2.6001595631652723, -2.5118864315095797, -2.4266100950824145, -2.344228815319923, -2.2646443075930596, -2.1877616239495516, -2.1134890398366455, -2.0417379446695296, -1.9724227361148534, -1.9054607179632463, -1.8407720014689564, -1.7782794100389228, -1.7179083871575875, -1.6595869074375598, -1.6032453906900417, -1.5488166189124812, -1.4962356560944328, -1.4454397707459279, -1.3963683610559376, -1.3489628825916533, -1.3031667784522987, -1.2589254117941675, -1.2161860006463678, -1.174897554939529, -1.1350108156723142, -1.096478196143185, -1.0592537251772887, -1.0232929922807537, -0.9885530946569385, -0.9549925860214359, -0.9225714271547628, -0.8912509381337455, -0.8609937521846003, -0.8317637711026709, -0.8035261221856173, -0.7762471166286916, -0.7498942093324559, -0.7244359600749899, -0.6998419960022735, -0.6760829753919816, -0.6531305526474723, -0.630957344480193, -0.609536897240169, -0.588843655355589, -0.5688529308438414, -0.5495408738576245, -0.5308844442309882, -0.5128613839913648, -0.49545019080479, -0.47863009232263826, -0.46238102139926035, -0.44668359215096304, -0.43151907682776525, -0.4168693834703353, -0.4027170343254591, -0.3890451449942805, -0.37583740428844414, -0.3630780547701014, -0.35075187395256796, -0.33884415613920255, -0.3273406948788381, -0.31622776601683794, -0.30549211132155124, -0.29512092266663853, -0.2851018267503908, -0.2754228703338166, -0.26607250597988097, -0.25703957827688634, -0.24831331052955705, -0.239883291901949, -0.23173946499684786, -0.2238721138568339, -0.21627185237270202, -0.20892961308540386, -0.20183663636815607, -0.19498445997580455, -0.18836490894898, -0.18197008586099836, -0.1757923613958692, -0.16982436524617442, -0.16405897731995386, -0.15848931924611134, -0.15310874616820302, -0.1479108388168207, -0.1428893958511103, -0.13803842646028847, -0.1333521432163324, -0.12882495516931336, -0.1244514611771385, -0.12022644346174131, -0.11614486138403426, -0.11220184543019636, -0.10839269140212034, -0.10471285480508996, -0.10115794542598983, -0.09772372209558107, -0.09440608762859236, -0.09120108393559097, -0.08810488730080138, -0.08511380382023763, -0.08222426499470713, -0.07943282347242814, -0.0767361489361819, -0.07413102413009177, -0.0716143410212902, -0.06918309709189363, -0.06683439175686146, -0.06456542290346556, -0.06237348354824192, -0.06025595860743578, -0.05821032177708716, -0.05623413251903491, -0.05432503314924331, -0.05248074602497726, -0.050699070827470445, -0.04897788193684462, -0.04731512589614803, -0.04570881896148749, -0.044157044735331254, -0.04265795188015926, -0.04120975190973302, -0.039810717055349734, -0.03845917820453535, -0.03715352290971724, -0.03589219346450052, -0.034673685045253165, -0.03349654391578276, -0.03235936569296283, -0.03126079367123956, -0.03019951720402016, -0.02917427014001166, -0.028183829312644536, -0.027227013080779128, -0.026302679918953815, -0.02540972705549305, -0.02454708915685031, -0.023713737056616554, -0.022908676527677724, -0.022130947096056376, -0.021379620895022326, -0.02065380155810529, -0.0199526231496888, -0.019275249131909356, -0.018620871366628676, -0.017988709151287873, -0.017378008287493755, -0.016788040181225608, -0.0162181009735893, -0.015667510701081494, -0.01513561248436208, -0.014621771744567183, -0.01412537544622754, -0.013645831365889245, -0.013182567385564075, -0.012735030810166616, -0.012302687708123818, -0.011885022274370183, -0.01148153621496883, -0.01109174815262401, -0.010715193052376065, -0.010351421666793436, -0.01, 0.01, 0.010351421666793436, 0.010715193052376065, 0.01109174815262401, 0.01148153621496883, 0.011885022274370183, 0.012302687708123818, 0.012735030810166616, 0.013182567385564075, 0.013645831365889245, 0.01412537544622754, 0.014621771744567183, 0.01513561248436208, 0.015667510701081494, 0.0162181009735893, 0.016788040181225608, 0.017378008287493755, 0.017988709151287873, 0.018620871366628676, 0.019275249131909356, 0.0199526231496888, 0.02065380155810529, 0.021379620895022326, 0.022130947096056376, 0.022908676527677724, 0.023713737056616554, 0.02454708915685031, 0.02540972705549305, 0.026302679918953815, 0.027227013080779128, 0.028183829312644536, 0.02917427014001166, 0.03019951720402016, 0.03126079367123956, 0.03235936569296283, 0.03349654391578276, 0.034673685045253165, 0.03589219346450052, 0.03715352290971724, 0.03845917820453535, 0.039810717055349734, 0.04120975190973302, 0.04265795188015926, 0.044157044735331254, 0.04570881896148749, 0.04731512589614803, 0.04897788193684462, 0.050699070827470445, 0.05248074602497726, 0.05432503314924331, 0.05623413251903491, 0.05821032177708716, 0.06025595860743578, 0.06237348354824192, 0.06456542290346556, 0.06683439175686146, 0.06918309709189363, 0.0716143410212902, 0.07413102413009177, 0.0767361489361819, 0.07943282347242814, 0.08222426499470713, 0.08511380382023763, 0.08810488730080138, 0.09120108393559097, 0.09440608762859236, 0.09772372209558107, 0.10115794542598983, 0.10471285480508996, 0.10839269140212034, 0.11220184543019636, 0.11614486138403426, 0.12022644346174131, 0.1244514611771385, 0.12882495516931336, 0.1333521432163324, 0.13803842646028847, 0.1428893958511103, 0.1479108388168207, 0.15310874616820302, 0.15848931924611134, 0.16405897731995386, 0.16982436524617442, 0.1757923613958692, 0.18197008586099836, 0.18836490894898, 0.19498445997580455, 0.20183663636815607, 0.20892961308540386, 0.21627185237270202, 0.2238721138568339, 0.23173946499684786, 0.239883291901949, 0.24831331052955705, 0.25703957827688634, 0.26607250597988097, 0.2754228703338166, 0.2851018267503908, 0.29512092266663853, 0.30549211132155124, 0.31622776601683794, 0.3273406948788381, 0.33884415613920255, 0.35075187395256796, 0.3630780547701014, 0.37583740428844414, 0.3890451449942805, 0.4027170343254591, 0.4168693834703353, 0.43151907682776525, 0.44668359215096304, 0.46238102139926035, 0.47863009232263826, 0.49545019080479, 0.5128613839913648, 0.5308844442309882, 0.5495408738576245, 0.5688529308438414, 0.588843655355589, 0.609536897240169, 0.630957344480193, 0.6531305526474723, 0.6760829753919816, 0.6998419960022735, 0.7244359600749899, 0.7498942093324559, 0.7762471166286916, 0.8035261221856173, 0.8317637711026709, 0.8609937521846003, 0.8912509381337455, 0.9225714271547628, 0.9549925860214359, 0.9885530946569385, 1.0232929922807537, 1.0592537251772887, 1.096478196143185, 1.1350108156723142, 1.174897554939529, 1.2161860006463678, 1.2589254117941675, 1.3031667784522987, 1.3489628825916533, 1.3963683610559376, 1.4454397707459279, 1.4962356560944328, 1.5488166189124812, 1.6032453906900417, 1.6595869074375598, 1.7179083871575875, 1.7782794100389228, 1.8407720014689564, 1.9054607179632463, 1.9724227361148534, 2.0417379446695296, 2.1134890398366455, 2.1877616239495516, 2.2646443075930596, 2.344228815319923, 2.4266100950824145, 2.5118864315095797, 2.6001595631652723, 2.691534803926914, 2.7861211686297693, 2.8840315031266055, 2.9853826189179604, 3.0902954325135887, 3.198895109691397, 3.311311214825911, 3.427677865464501, 3.5481338923357533, 3.672823004980846, 3.801893963205613, 3.9355007545577725, 4.073802778041126, 4.216965034285822, 4.365158322401657, 4.518559443749222, 4.677351412871981, 4.841723675840994, 5.011872336272719, 5.188000389289609, 5.370317963702527, 5.559042572704037, 5.754399373371567, 5.956621435290104, 6.165950018614822, 6.3826348619054825, 6.6069344800759575, 6.839116472814292, 7.07945784384138, 7.328245331389037, 7.585775750291836, 7.852356346100718, 8.128305161640986, 8.413951416451948, 8.709635899560805, 9.01571137605957, 9.332543007969905, 9.660508789898131, 10.0}, "Binning in TPC inner param. or pT"};
  ConfigurableAxis dEdxBins{"dEdxBins", {3000, 0.f, 1500.f}, "Binning in dE/dx"};
  ConfigurableAxis nsigmaBins{"nsigmaBins", {200, -5, 5}, "Binning in nsigma"};

  void init(InitContext const&)
  {
    AxisSpec pAxis{pBins, "z#upoint p (GeV/c)"};
    if (usePt) {
      pAxis.title = "#it{p}_{T} (GeV/c)";
    }
    const AxisSpec dedxAxis{dEdxBins, "d#it{E}/d#it{x} Arb. units"};
    const AxisSpec nsigmaAxis{nsigmaBins, "n#sigma_{TPC}"};

    // Raw dE/dx vs. TPC momentum
    registryCh.add(
      "dEdx_vs_Momentum", "dE/dx", HistType::kTH2F,
      {pAxis, dEdxBins});
    registryPi.add(
      "dEdx_vs_Momentum_Pi", "dE/dx", HistType::kTH2F,
      {pAxis, dEdxBins});
    registryKa.add(
      "dEdx_vs_Momentum_Ka", "dE/dx", HistType::kTH2F,
      {pAxis, dEdxBins});
    registryPr.add(
      "dEdx_vs_Momentum_Pr", "dE/dx", HistType::kTH2F,
      {pAxis, dEdxBins});
    registryDe.add(
      "dEdx_vs_Momentum_De", "dE/dx", HistType::kTH2F,
      {pAxis, dEdxBins});
    registryTr.add(
      "dEdx_vs_Momentum_Tr", "dE/dx", HistType::kTH2F,
      {pAxis, dEdxBins});
    registryHe.add(
      "dEdx_vs_Momentum_He", "dE/dx", HistType::kTH2F,
      {pAxis, dEdxBins});

    // nsigma_(TPC) vs. TPC momentum
    registryPi.add(
      "nsigmaTPC_vs_Momentum_Pi", "nsigmaTPC", HistType::kTH2F,
      {pAxis, nsigmaAxis});
    registryKa.add(
      "nsigmaTPC_vs_Momentum_Ka", "nsigmaTPC", HistType::kTH2F,
      {pAxis, nsigmaAxis});
    registryPr.add(
      "nsigmaTPC_vs_Momentum_Pr", "nsigmaTPC", HistType::kTH2F,
      {pAxis, nsigmaAxis});
    registryDe.add(
      "nsigmaTPC_vs_Momentum_De", "nsigmaTPC", HistType::kTH2F,
      {pAxis, nsigmaAxis});
    registryTr.add(
      "nsigmaTPC_vs_Momentum_Tr", "nsigmaTPC", HistType::kTH2F,
      {pAxis, nsigmaAxis});
    registryHe.add(
      "nsigmaTPC_vs_Momentum_He", "nsigmaTPC", HistType::kTH2F,
      {pAxis, nsigmaAxis});
    if (addTOF) {
      registryPi.add("nsigmaTOF_vs_Momentum_Pi", "nsigmaTOF", HistType::kTH2F, {pAxis, nsigmaAxis});
      registryKa.add("nsigmaTOF_vs_Momentum_Ka", "nsigmaTOF", HistType::kTH2F, {pAxis, nsigmaAxis});
      registryPr.add("nsigmaTOF_vs_Momentum_Pr", "nsigmaTOF", HistType::kTH2F, {pAxis, nsigmaAxis});
    }

    if (doContaminations) { // If requested, also produce the contaminations plots
      registryPi.add("nsigmaTPC_vs_Momentum_Ka", "nsigmaTPC", HistType::kTH2F, {pAxis, nsigmaAxis});
      registryPi.add("nsigmaTPC_vs_Momentum_Pr", "nsigmaTPC", HistType::kTH2F, {pAxis, nsigmaAxis});

      registryKa.add("nsigmaTPC_vs_Momentum_Pi", "nsigmaTPC", HistType::kTH2F, {pAxis, nsigmaAxis});
      registryKa.add("nsigmaTPC_vs_Momentum_Pr", "nsigmaTPC", HistType::kTH2F, {pAxis, nsigmaAxis});

      registryPr.add("nsigmaTPC_vs_Momentum_Pi", "nsigmaTPC", HistType::kTH2F, {pAxis, nsigmaAxis});
      registryPr.add("nsigmaTPC_vs_Momentum_Ka", "nsigmaTPC", HistType::kTH2F, {pAxis, nsigmaAxis});

      if (addTOF) {
        registryPi.add("nsigmaTOF_vs_Momentum_Ka", "nsigmaTOF", HistType::kTH2F, {pAxis, nsigmaAxis});
        registryPi.add("nsigmaTOF_vs_Momentum_Pr", "nsigmaTOF", HistType::kTH2F, {pAxis, nsigmaAxis});

        // Not ok until we remove ID with TOF
        // registryKa.add("nsigmaTOF_vs_Momentum_Pi", "nsigmaTOF", HistType::kTH2F, {pAxis, nsigmaAxis});
        // registryKa.add("nsigmaTOF_vs_Momentum_Pr", "nsigmaTOF", HistType::kTH2F, {pAxis, nsigmaAxis});

        registryPr.add("nsigmaTOF_vs_Momentum_Pi", "nsigmaTOF", HistType::kTH2F, {pAxis, nsigmaAxis});
        registryPr.add("nsigmaTOF_vs_Momentum_Ka", "nsigmaTOF", HistType::kTH2F, {pAxis, nsigmaAxis});
      }
    }
    // Event Counter
    registryCh.add("histRecVtxZData", "collision z position", HistType::kTH1F, {{200, -20.0, +20.0, "z_{vtx} (cm)"}});
  }

  // Single-Track Selection
  template <typename T1, typename C>
  bool passedSingleTrackSelection(const T1& track, const C& /*collision*/)
  {
    // Single-Track Selections
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsFound() < minTPCnClsFound)
      return false;
    if (track.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    // if (track.tpcSignalN() < minNClsTPCdEdx)
    // return false;
    if (track.tpcChi2NCl() > maxChi2TPC)
      return false;
    if (track.eta() < etaMin || track.eta() > etaMax)
      return false;

    return true;
  }

  // General V0 Selections
  template <typename T1, typename C>
  bool passedV0Selection(const T1& v0, const C& /*collision*/)
  {
    if (v0.v0cosPA() < v0cospaMin)
      return false;
    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;

    return true;
  }

  // K0s Selections
  template <typename T1, typename T2, typename C>
  bool passedK0Selection(const T1& v0, const T2& ntrack, const T2& ptrack,
                         const C& collision)
  {
    // Single-Track Selections
    if (!passedSingleTrackSelection(ptrack, collision))
      return false;
    if (!passedSingleTrackSelection(ntrack, collision))
      return false;

    if (ptrack.tpcInnerParam() > 0.6) {
      if (!ptrack.hasTOF())
        return false;
      if (TMath::Abs(ptrack.tofNSigmaPi()) > nsigmaTOFmax)
        return false;
    }

    if (ntrack.tpcInnerParam() > 0.6) {
      if (!ntrack.hasTOF())
        return false;
      if (TMath::Abs(ntrack.tofNSigmaPi()) > nsigmaTOFmax)
        return false;
    }

    // Invariant-Mass Selection
    if (v0.mK0Short() < minMassK0s || v0.mK0Short() > maxMassK0s)
      return false;

    return true;
  }

  // Lambda Selections
  template <typename T1, typename T2, typename C>
  bool passedLambdaSelection(const T1& v0, const T2& ntrack, const T2& ptrack,
                             const C& collision)
  {
    // Single-Track Selections
    if (!passedSingleTrackSelection(ptrack, collision))
      return false;
    if (!passedSingleTrackSelection(ntrack, collision))
      return false;

    if (ptrack.tpcInnerParam() > 0.6) {
      if (!ptrack.hasTOF())
        return false;
      if (TMath::Abs(ptrack.tofNSigmaPr()) > nsigmaTOFmax)
        return false;
    }

    if (ntrack.tpcInnerParam() > 0.6) {
      if (!ntrack.hasTOF())
        return false;
      if (TMath::Abs(ntrack.tofNSigmaPi()) > nsigmaTOFmax)
        return false;
    }

    // Invariant-Mass Selection
    if (v0.mLambda() < minMassLambda || v0.mLambda() > maxMassLambda)
      return false;

    return true;
  }

  // AntiLambda Selections
  template <typename T1, typename T2, typename C>
  bool passedAntiLambdaSelection(const T1& v0, const T2& ntrack,
                                 const T2& ptrack, const C& collision)
  {

    // Single-Track Selections
    if (!passedSingleTrackSelection(ptrack, collision))
      return false;
    if (!passedSingleTrackSelection(ntrack, collision))
      return false;

    if (ptrack.tpcInnerParam() > 0.6) {
      if (!ptrack.hasTOF())
        return false;
      if (TMath::Abs(ptrack.tofNSigmaPi()) > nsigmaTOFmax)
        return false;
    }

    if (ntrack.tpcInnerParam() > 0.6) {
      if (!ntrack.hasTOF())
        return false;
      if (TMath::Abs(ntrack.tofNSigmaPr()) > nsigmaTOFmax)
        return false;
    }

    // Invariant-Mass Selection
    if (v0.mAntiLambda() < minMassLambda || v0.mAntiLambda() > maxMassLambda)
      return false;

    return true;
  }

  // Process Data
  void process(SelectedCollisions::iterator const& collision,
               aod::V0Datas const& fullV0s, PIDTracks const& tracks)
  {
    // Event Selection
    if (!collision.sel8())
      return;

    // Event Counter
    registryCh.fill(HIST("histRecVtxZData"), collision.posZ());

    // Kaons and nuclei
    for (auto& trk : tracks) {

      if (!passedSingleTrackSelection(trk, collision))
        continue;
      if (!trk.passedTPCRefit())
        continue;
      float signedP = trk.sign() * trk.tpcInnerParam();
      if (usePt) {
        signedP = trk.sign() * trk.pt();
      }

      // Charged Particles
      registryCh.fill(HIST("dEdx_vs_Momentum"), signedP,
                      trk.tpcSignal());

      // Kaons
      if (trk.tpcInnerParam() > 0.4 && trk.hasTOF() && TMath::Abs(trk.tofNSigmaKa()) < 2.0) {
        registryKa.fill(HIST("dEdx_vs_Momentum_Ka"), signedP, trk.tpcSignal());
        registryKa.fill(HIST("nsigmaTPC_vs_Momentum_Ka"), signedP, trk.tpcNSigmaKa());

        if (doContaminations) {
          registryKa.fill(HIST("nsigmaTPC_vs_Momentum_Pi"), signedP, trk.tpcNSigmaPi());
          registryKa.fill(HIST("nsigmaTPC_vs_Momentum_Pr"), signedP, trk.tpcNSigmaPr());
        }
      }

      if (trk.tpcInnerParam() < 0.4) {
        registryKa.fill(HIST("dEdx_vs_Momentum_Ka"), signedP,
                        trk.tpcSignal());
        registryKa.fill(HIST("nsigmaTPC_vs_Momentum_Ka"), signedP,
                        trk.tpcNSigmaKa());

        if (doContaminations) {
          registryKa.fill(HIST("nsigmaTPC_vs_Momentum_Pi"), signedP,
                          trk.tpcNSigmaPi());
          registryKa.fill(HIST("nsigmaTPC_vs_Momentum_Pr"), signedP,
                          trk.tpcNSigmaPr());
        }
      }

      // Selection of high dE/dx Objects
      if (trk.tpcNSigmaDe() < -4.0)
        continue;

      // Deuterons
      if (trk.tpcInnerParam() > 1.0 && trk.hasTOF() && TMath::Abs(trk.tofNSigmaDe()) < 3.0) {
        registryDe.fill(HIST("dEdx_vs_Momentum_De"), signedP,
                        trk.tpcSignal());
        registryDe.fill(HIST("nsigmaTPC_vs_Momentum_De"), signedP,
                        trk.tpcNSigmaDe());
      }

      if (trk.tpcInnerParam() < 1.0) {
        registryDe.fill(HIST("dEdx_vs_Momentum_De"), signedP,
                        trk.tpcSignal());
        registryDe.fill(HIST("nsigmaTPC_vs_Momentum_De"), signedP,
                        trk.tpcNSigmaDe());
      }

      // Heavier Nuclei
      registryTr.fill(HIST("dEdx_vs_Momentum_Tr"), signedP,
                      trk.tpcSignal());
      registryTr.fill(HIST("nsigmaTPC_vs_Momentum_Tr"), signedP,
                      trk.tpcNSigmaTr());
      registryHe.fill(HIST("dEdx_vs_Momentum_He"), 2.0 * signedP,
                      trk.tpcSignal());
      registryHe.fill(HIST("nsigmaTPC_vs_Momentum_He"),
                      2.0 * signedP, trk.tpcNSigmaHe());
    }

    // Loop over Reconstructed V0s
    for (auto& v0 : fullV0s) {

      // Standard V0 Selections
      if (!passedV0Selection(v0, collision)) {
        continue;
      }

      if (v0.dcaV0daughters() > dcaV0DaughtersMax) {
        continue;
      }

      // Positive and Negative Tracks
      const auto& posTrack = v0.posTrack_as<PIDTracks>();
      const auto& negTrack = v0.negTrack_as<PIDTracks>();

      if (!posTrack.passedTPCRefit())
        continue;
      if (!negTrack.passedTPCRefit())
        continue;

      float signedPpos = posTrack.sign() * posTrack.tpcInnerParam();
      float signedPneg = negTrack.sign() * negTrack.tpcInnerParam();
      if (usePt) {
        signedPpos = posTrack.sign() * posTrack.pt();
        signedPneg = negTrack.sign() * negTrack.pt();
      }

      // K0s Selection
      if (passedK0Selection(v0, negTrack, posTrack, collision)) {
        registryPi.fill(HIST("dEdx_vs_Momentum_Pi"), signedPneg, negTrack.tpcSignal());
        registryPi.fill(HIST("dEdx_vs_Momentum_Pi"), signedPpos, posTrack.tpcSignal());
        registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Pi"),
                        signedPneg, negTrack.tpcNSigmaPi());
        registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Pi"),
                        signedPpos, posTrack.tpcNSigmaPi());
        if (addTOF) {
          registryPi.fill(HIST("nsigmaTOF_vs_Momentum_Pi"), signedPneg, negTrack.tofNSigmaPi());
          registryPi.fill(HIST("nsigmaTOF_vs_Momentum_Pi"), signedPpos, posTrack.tofNSigmaPi());
        }

        if (doContaminations) {
          registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Ka"), signedPneg, negTrack.tpcNSigmaKa());
          registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Pr"), signedPneg, negTrack.tpcNSigmaPr());
          registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Ka"), signedPpos, posTrack.tpcNSigmaKa());
          registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Pr"), signedPpos, posTrack.tpcNSigmaPr());
          if (addTOF) {
            registryPi.fill(HIST("nsigmaTOF_vs_Momentum_Ka"), signedPneg, negTrack.tofNSigmaKa());
            registryPi.fill(HIST("nsigmaTOF_vs_Momentum_Pr"), signedPneg, negTrack.tofNSigmaPr());
            registryPi.fill(HIST("nsigmaTOF_vs_Momentum_Ka"), signedPpos, posTrack.tofNSigmaKa());
            registryPi.fill(HIST("nsigmaTOF_vs_Momentum_Pr"), signedPpos, posTrack.tofNSigmaPr());
          }
        }
      }

      // Lambda Selection
      if (passedLambdaSelection(v0, negTrack, posTrack, collision)) {
        registryPr.fill(HIST("dEdx_vs_Momentum_Pr"), signedPpos, posTrack.tpcSignal());
        registryPr.fill(HIST("nsigmaTPC_vs_Momentum_Pr"),
                        signedPpos, posTrack.tpcNSigmaPr());
        registryPi.fill(HIST("dEdx_vs_Momentum_Pi"), signedPneg, negTrack.tpcSignal());
        registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Pi"),
                        signedPneg, negTrack.tpcNSigmaPi());
        if (addTOF) {
          registryPr.fill(HIST("nsigmaTOF_vs_Momentum_Pr"), signedPpos, posTrack.tofNSigmaPr());
          registryPi.fill(HIST("nsigmaTOF_vs_Momentum_Pi"), signedPneg, negTrack.tofNSigmaPi());
        }
        if (doContaminations) {
          registryPr.fill(HIST("nsigmaTPC_vs_Momentum_Pi"), signedPpos, posTrack.tpcNSigmaPi());
          registryPr.fill(HIST("nsigmaTPC_vs_Momentum_Ka"), signedPpos, posTrack.tpcNSigmaKa());

          registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Ka"), signedPneg, negTrack.tpcNSigmaKa());
          registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Pr"), signedPneg, negTrack.tpcNSigmaPr());
          if (addTOF) {
            registryPr.fill(HIST("nsigmaTOF_vs_Momentum_Pi"), signedPpos, posTrack.tofNSigmaPi());
            registryPr.fill(HIST("nsigmaTOF_vs_Momentum_Ka"), signedPpos, posTrack.tofNSigmaKa());

            registryPi.fill(HIST("nsigmaTOF_vs_Momentum_Ka"), signedPneg, negTrack.tofNSigmaKa());
            registryPi.fill(HIST("nsigmaTOF_vs_Momentum_Pr"), signedPneg, negTrack.tofNSigmaPr());
          }
        }
      }

      // AntiLambda Selection
      if (passedAntiLambdaSelection(v0, negTrack, posTrack, collision)) {
        registryPi.fill(HIST("dEdx_vs_Momentum_Pi"), signedPpos, posTrack.tpcSignal());
        registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Pi"),
                        signedPpos, posTrack.tpcNSigmaPi());
        registryPr.fill(HIST("dEdx_vs_Momentum_Pr"), signedPneg, negTrack.tpcSignal());
        registryPr.fill(HIST("nsigmaTPC_vs_Momentum_Pr"),
                        signedPneg, negTrack.tpcNSigmaPr());
        if (addTOF) {
          registryPi.fill(HIST("nsigmaTOF_vs_Momentum_Pi"), signedPpos, posTrack.tofNSigmaPi());
          registryPr.fill(HIST("nsigmaTOF_vs_Momentum_Pr"), signedPneg, negTrack.tofNSigmaPr());
        }
        if (doContaminations) {
          registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Ka"), signedPpos, posTrack.tpcNSigmaKa());
          registryPi.fill(HIST("nsigmaTPC_vs_Momentum_Pr"), signedPpos, posTrack.tpcNSigmaPr());

          registryPr.fill(HIST("nsigmaTPC_vs_Momentum_Pi"), signedPneg, negTrack.tpcNSigmaPi());
          registryPr.fill(HIST("nsigmaTPC_vs_Momentum_Ka"), signedPneg, negTrack.tpcNSigmaKa());
          if (addTOF) {
            registryPi.fill(HIST("nsigmaTOF_vs_Momentum_Ka"), signedPpos, posTrack.tofNSigmaKa());
            registryPi.fill(HIST("nsigmaTOF_vs_Momentum_Pr"), signedPpos, posTrack.tofNSigmaPr());

            registryPr.fill(HIST("nsigmaTOF_vs_Momentum_Pi"), signedPneg, negTrack.tofNSigmaPi());
            registryPr.fill(HIST("nsigmaTOF_vs_Momentum_Ka"), signedPneg, negTrack.tofNSigmaKa());
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<tpc_dEdx_postcalibration>(cfgc)};
}
