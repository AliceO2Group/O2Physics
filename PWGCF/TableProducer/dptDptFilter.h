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

/// \file dptDptFilter.h
/// \brief Filters collisions and tracks according to selection criteria
/// \author victor.gonzalez.sebastian@gmail.com

#ifndef PWGCF_TABLEPRODUCER_DPTDPTFILTER_H_
#define PWGCF_TABLEPRODUCER_DPTDPTFILTER_H_

#include "PWGCF/Core/AnalysisConfigurableCuts.h"

#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/Core/MetadataHelper.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "ReconstructionDataFormats/PID.h"
#include <CCDB/BasicCCDBManager.h>

#include <TF1.h>
#include <TFormula.h>
#include <TList.h>
#include <TMCProcess.h>
#include <TPDGCode.h>

#include <Rtypes.h>

#include <bitset>
#include <fstream>
#include <functional>
#include <iomanip>
#include <locale>
#include <map>
#include <ranges>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace o2
{
namespace aod
{
using CollisionsEvSelCent = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentNTPVs>;
using CollisionEvSelCent = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentNTPVs>::iterator;
using CollisionsEvSelRun2Cent = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::CentRun2V0Ms, aod::CentRun2CL0s, aod::CentRun2CL1s>;
using CollisionEvSelRun2Cent = soa::Join<aod::Collisions, aod::Mults, aod::EvSels, aod::CentRun2V0Ms, aod::CentRun2CL0s, aod::CentRun2CL1s>::iterator;
using CollisionsEvSel = soa::Join<aod::Collisions, aod::Mults, aod::EvSels>;
using CollisionEvSel = soa::Join<aod::Collisions, aod::Mults, aod::EvSels>::iterator;
using TrackData = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>::iterator;
} // namespace aod
namespace analysis
{
namespace dptdptfilter
{
/// \enum SystemType
/// \brief The type of the system under analysis
enum SystemType {
  SystemNoSystem = 0, ///< no system defined
  SystemPp,           ///< **p-p** system
  SystemPPb,          ///< **p-Pb** system
  SystemPbp,          ///< **Pb-p** system
  SystemPbPb,         ///< **Pb-Pb** system
  SystemXeXe,         ///< **Xe-Xe** system
  SystemPpRun3,       ///< **p-p Run 3** system
  SystemPbPbRun3,     ///< **Pb-Pb Run 3** system
  SystemNeNeRun3,     ///< **Ne-Ne Run 3** system
  SystemOORun3,       ///< **O-O Run 3** system
  SystemPORun3,       ///< **p-O Run 3** system
  SystemNoOfSystems   ///< number of handled systems
};

/// @brief SystemType prefix increment operator
/// @param ipar value
/// @return the incremented value
inline SystemType& operator++(SystemType& ipar)
{
  return ipar = static_cast<SystemType>(static_cast<int>(ipar) + 1);
}

/// @brief SystemType postfix increment operator
/// @param ipar the value
/// @param empty
/// @return the same value
inline SystemType operator++(SystemType& ipar, int)
{
  SystemType iparTmp(ipar);
  ++ipar;
  return iparTmp;
}

/// \std::map systemInternalCodesMap
/// \brief maps system names to internal system codes
static const std::map<std::string_view, int> systemInternalCodesMap{
  {"", SystemNoSystem},
  {"pp", SystemPp},
  {"pPb", SystemPPb},
  {"Pbp", SystemPbp},
  {"PbPb", SystemPbPb},
  {"XeXe", SystemXeXe},
  {"ppRun3", SystemPpRun3},
  {"PbPbRun3", SystemPbPbRun3},
  {"NeNeRun3", SystemNeNeRun3},
  {"OORun3", SystemOORun3},
  {"pORun3", SystemPORun3}};

/// \std::map systemExternalNamesMap
/// \brief maps system internal codes to system external names
static const std::map<int, std::string_view> systemExternalNamesMap{
  {SystemNoSystem, ""},
  {SystemPp, "pp"},
  {SystemPPb, "pPb"},
  {SystemPbp, "Pbp"},
  {SystemPbPb, "PbPb"},
  {SystemXeXe, "XeXe"},
  {SystemPpRun3, "ppRun3"},
  {SystemPbPbRun3, "PbPbRun3"},
  {SystemNeNeRun3, "NeNeRun3"},
  {SystemOORun3, "OORun3"},
  {SystemPORun3, "pORun3"}};

/// \enum DataType
/// \brief Which kind of data is the task addressing
enum DataType {
  kData = 0,     ///< actual data, not generated
  kMC,           ///< Generator level and detector level
  kFastMC,       ///< Gererator level but stored dataset
  kOnTheFly,     ///< On the fly generator level data
  kDataNoEvtSel, ///< actual data but not event selection available yet
  knGenData      ///< number of different generator data types
};

/// \enum CentMultEstimatorType
/// \brief The detector used to estimate centrality/multiplicity
enum CentMultEstimatorType {
  CentMultNOCM = 0,      ///< do not use centrality/multiplicity estimator
  CentMultV0M,           ///< V0M centrality/multiplicity estimator Run 1/2
  CentMultCL0,           ///< CL0 centrality/multiplicity estimator Run 1/2
  CentMultCL1,           ///< CL1 centrality/multiplicity estimator Run 1/2
  CentMultFV0A,          ///< FV0A centrality/multiplicity estimator Run 3
  CentMultFT0M,          ///< FT0M centrality/multiplicity estimator Run 3
  CentMultFT0A,          ///< FT0A centrality/multiplicity estimator Run 3
  CentMultFT0C,          ///< FT0C centrality/multiplicity estimator Run 3
  CentMultNTPV,          ///< NTPV centrality/multiplicity estimator Run 3
  CentMultNOOFESTIMATORS ///< number of centrality/mutiplicity estimator
};

/// \std::map estimatorInternalCodesMap
/// \brief maps centrality/multiplicity estimator names to internal estimator codes
static const std::map<std::string_view, int> estimatorInternalCodesMap{
  {"NOCM", CentMultNOCM},
  {"V0M", CentMultV0M},
  {"CL0", CentMultCL0},
  {"CL1", CentMultCL1},
  {"FV0A", CentMultFV0A},
  {"FT0M", CentMultFT0M},
  {"FT0A", CentMultFT0A},
  {"FT0C", CentMultFT0C},
  {"NTPV", CentMultNTPV}};

/// \std::map estimatorExternalNamesMap
/// \brief maps internal estimator codes to centrality/multiplicity estimator external names
static const std::map<int, std::string_view> estimatorExternalNamesMap{
  {CentMultNOCM, "NOCM"},
  {CentMultV0M, "V0M"},
  {CentMultCL0, "CL0"},
  {CentMultCL1, "CL1"},
  {CentMultFV0A, "FV0A"},
  {CentMultFT0M, "FT0M"},
  {CentMultFT0A, "FT0A"},
  {CentMultFT0C, "FT0C"},
  {CentMultNTPV, "NTPV"}};

/// \enum MultSourceType
/// \brief The multiplicity source
enum MultSourceType {
  MultSourceT0A = 0,        ///< T0A multiplicity
  MultSourceT0C,            ///< T0C multiplicity
  MultSourceT0M,            ///< T0M multiplicity
  MultSourceV0A,            ///< V0A multiplicity
  MultSourceV0C,            ///< V0C multiplicity
  MultSourceV0M,            ///< V0M multiplicity
  MultSourceNtracks,        ///< number of tracks multiplicity
  MultSourcePvContributors, ///< number of primary vertex contributors
  MultSourceNOOFSOURCES     ///< number multiplicity sources
};

/// \enum MultRunType
/// \brief The multiplicity LHC run
enum MultRunType {
  MultRunRUN1RUN2 = 0, ///< LHC Run 1 or Run 2
  MultRunRUN3,         ///< LHC Run 3
  MultRunNOOFRUNS      ///< number of runs for multiplicity
};

/// \std::map multRunForSystemMap
/// \brief maps the system to the lhc Run for multiplicity
static const std::map<int, MultRunType> multRunForSystemMap{
  {SystemNoSystem, MultRunRUN1RUN2},
  {SystemPp, MultRunRUN1RUN2},
  {SystemPPb, MultRunRUN1RUN2},
  {SystemPbp, MultRunRUN1RUN2},
  {SystemPbPb, MultRunRUN1RUN2},
  {SystemXeXe, MultRunRUN1RUN2},
  {SystemPpRun3, MultRunRUN3},
  {SystemPbPbRun3, MultRunRUN3},
  {SystemNeNeRun3, MultRunRUN3},
  {SystemOORun3, MultRunRUN3},
  {SystemPORun3, MultRunRUN3}};

/// \std::map estimatorMultiplicitySourceMap
/// \brief maps internal estimator codes internal multiplicity sources
static const std::map<int, int> estimatorMultiplicitySourceMap{
  {CentMultNOCM, MultSourceT0C},
  {CentMultV0M, MultSourceV0M},
  {CentMultCL0, MultSourceT0C}, /* TODO: for Run1,2 */
  {CentMultCL1, MultSourceT0C}, /* TODO: for Run1,2 */
  {CentMultFV0A, MultSourceV0A},
  {CentMultFT0M, MultSourceT0M},
  {CentMultFT0A, MultSourceT0A},
  {CentMultFT0C, MultSourceT0C},
  {CentMultNTPV, MultSourcePvContributors}};

/// \std::vector<std::map> multiplicitySourceExternalNamesMap
/// \brief maps internal multiplicity source to external names for the LHC runs
static const std::vector<std::map<int, std::string_view>> multiplicitySourceExternalNamesMap{
  /* Run 1 and Run 2 */
  {
    {MultSourceT0A, "T0A multiplicity"},
    {MultSourceT0C, "T0C multiplicity"},
    {MultSourceT0M, "T0M multiplicity"},
    {MultSourceV0A, "V0A multiplicity"},
    {MultSourceV0C, "V0C multiplicity"},
    {MultSourceV0M, "V0M multiplicity"},
    {MultSourceNtracks, "Global tracks"},
    {MultSourcePvContributors, "PV contributors"}},
  /* Run 3 */
  {
    {MultSourceT0A, "FT0A multiplicity"},
    {MultSourceT0C, "FT0C multiplicity"},
    {MultSourceT0M, "FT0M multiplicity"},
    {MultSourceV0A, "FV0A multiplicity"},
    {MultSourceV0C, "WRONG SOURCE"},
    {MultSourceV0M, "FV0M multiplicity"},
    {MultSourceNtracks, "Global tracks"},
    {MultSourcePvContributors, "PV contributors"}}};

/// \std::map multiplicitySourceConfigNamesMap
/// \brief maps internal multiplicity source to external configuration names
/// At configuration time neither the system nor the lhc run is known
static const std::map<int, std::string_view> multiplicitySourceConfigNamesMap{
  {MultSourceT0A, "FT0A (T0A)"},
  {MultSourceT0C, "FT0C (T0C)"},
  {MultSourceT0M, "FT0M (T0M)"},
  {MultSourceV0A, "FV0A (V0A)"},
  {MultSourceV0C, "WRONG (V0C)"},
  {MultSourceV0M, "FV0M (V0M)"},
  {MultSourceNtracks, "Global tracks"},
  {MultSourcePvContributors, "PV contributors"}};

/// \enum CentMultCorrelationsParams
/// \brief internal codes for the supported parameters for centrality/multiplicity correlations exclusion cuts
enum CentMultCorrelationsParams {
  CentMultCorrelationsMT0A = 0,  ///< multiplicity from T0A
  CentMultCorrelationsMT0C,      ///< multiplicity from T0C
  CentMultCorrelationsMT0M,      ///< multiplicity from T0M
  CentMultCorrelationsMV0A,      ///< multiplicity from V0A
  CentMultCorrelationsMV0C,      ///< multiplicity from V0C (only Run 1 & Run 2)
  CentMultCorrelationsMV0M,      ///< multiplicity from V0M
  CentMultCorrelationsMNGLTRK,   ///< multiplicity from number of global tracks
  CentMultCorrelationsMNPVC,     ///< multiplicity from number of PV contributors
  CentMultCorrelationsCT0A,      ///< centrality from T0A
  CentMultCorrelationsCT0C,      ///< centrality from T0C
  CentMultCorrelationsCT0M,      ///< centrality from T0M
  CentMultCorrelationsCV0A,      ///< centrality from V0A
  CentMultCorrelationsCNTPV,     ///< centrality from number of PV contributors
  CentMultCorrelationsNOOFPARAMS ///< the number of parameters supported
};

/// @brief prefix increment operator
/// @param ipar value
/// @return the incremented value
inline CentMultCorrelationsParams& operator++(CentMultCorrelationsParams& ipar)
{
  return ipar = static_cast<CentMultCorrelationsParams>(static_cast<int>(ipar) + 1);
}

/// @brief postfix increment operator
/// @param ipar the value
/// @param empty
/// @return the same value
inline CentMultCorrelationsParams operator++(CentMultCorrelationsParams& ipar, int)
{
  CentMultCorrelationsParams iparTmp(ipar);
  ++ipar;
  return iparTmp;
}

/// \std::map centMultCorrelationsParamsMap
/// \brief maps centrality/multiplicity correlations parameters names to internal codes
static const std::map<std::string_view, CentMultCorrelationsParams> centMultCorrelationsParamsMap{
  {"MT0A", CentMultCorrelationsMT0A},
  {"MT0C", CentMultCorrelationsMT0C},
  {"MT0M", CentMultCorrelationsMT0M},
  {"MV0A", CentMultCorrelationsMV0A},
  {"MV0C", CentMultCorrelationsMV0C},
  {"MV0M", CentMultCorrelationsMV0M},
  {"MNGLTRK", CentMultCorrelationsMNGLTRK},
  {"MNPVC", CentMultCorrelationsMNPVC},
  {"CT0A", CentMultCorrelationsCT0A},
  {"CT0C", CentMultCorrelationsCT0C},
  {"CT0M", CentMultCorrelationsCT0M},
  {"CV0A", CentMultCorrelationsCV0A},
  {"CNTPV", CentMultCorrelationsCNTPV}};

/// \std::map centMultCorrelationsParamsNamesMap
/// \brief maps centrality/multiplicity correlations parameters internal codes to their external names
static const std::map<CentMultCorrelationsParams, std::string_view> centMultCorrelationsParamsNamesMap{
  {CentMultCorrelationsMT0A, "MT0A"},
  {CentMultCorrelationsMT0C, "MT0C"},
  {CentMultCorrelationsMT0M, "MT0M"},
  {CentMultCorrelationsMV0A, "MV0A"},
  {CentMultCorrelationsMV0C, "MV0C"},
  {CentMultCorrelationsMV0M, "MV0M"},
  {CentMultCorrelationsMNGLTRK, "MNGLTRK"},
  {CentMultCorrelationsMNPVC, "MNPVC"},
  {CentMultCorrelationsCT0A, "CT0A"},
  {CentMultCorrelationsCT0C, "CT0C"},
  {CentMultCorrelationsCT0M, "CT0M"},
  {CentMultCorrelationsCV0A, "CV0A"},
  {CentMultCorrelationsCNTPV, "CNTPV"}};

/// \enum TriggerSelectionTags
/// \brief The potential trigger tags to apply for event selection
enum TriggerSelectionTags {
  TriggSelNONE = 0,           ///< do not use trigger selection
  TriggSelMB,                 ///< Minimum bias trigger
  TriggSelNOSAMEBUNCHPUP,     ///< No same bunch pile up
  TriggSelNUMPVCONTRIBUTORS,  ///< Number of primary vertex contributors
  TriggSelVTXTOFMATCHED,      ///< at least one primary vertex contributor is matched to TOF
  TriggSelVTXTRDMATCHED,      ///< at least one primary vertex contributor is matched to TRD
  TriggSelNOCOLLINTRSTD,      ///< no other collision in standard time range gap
  TriggSelNOCOLLINROFSTD,     ///< no other collision in standard readout frame gap
  TriggSelISVERTEXITSTPC,     ///< primary vertex contributors are matched tracks ITS+TPC
  TriggSelISGOODZVTXFT0VSPV,  ///< vertex extracted from FT0 is compatible with the one from primary vertex contributors
  TriggSelGOODITSLAYER3,      ///< good the 3 ITS layer
  TriggSelGOODITSLAYER0123,   ///< check good the 0,1,2,and 3 ITS layers
  TriggSelGOODITSLAYERALL,    ///< check good all ITS layers
  TriggSelNOGOODITSLAYER3,    ///< check no good the  3 ITS layer
  TriggSelNOGOODITSLAYER0123, ///< check no good the 0,1,2,and 3 ITS layers
  TriggSelNOGOODITSLAYERALL,  ///< check no good all ITS layers
  TriggSelNOOFTRIGGERS        ///< number of triggers for event selection
};

/// \std::map triggerSelectionBitsMap
/// \brief maps trigger selection tags to internal trigger selection bits
static const std::map<std::string_view, int> triggerSelectionBitsMap{
  {"none", TriggSelNONE},
  {"mb", TriggSelMB},
  {"nosamebunchpup", TriggSelNOSAMEBUNCHPUP},
  {"numpvcontr", TriggSelNUMPVCONTRIBUTORS},
  {"vtxtofmatched", TriggSelVTXTOFMATCHED},
  {"vtxtrdmatched", TriggSelVTXTRDMATCHED},
  {"nocollintrstd", TriggSelNOCOLLINTRSTD},
  {"nocollinrofstd", TriggSelNOCOLLINROFSTD},
  {"isvtxitstpc", TriggSelISVERTEXITSTPC},
  {"isgoodvtxft0vspv", TriggSelISGOODZVTXFT0VSPV},
  {"gooditslayer3", TriggSelGOODITSLAYER3},
  {"gooditslayer0123", TriggSelGOODITSLAYER0123},
  {"gooditslayerall", TriggSelGOODITSLAYERALL},
  {"nogooditslayer3", TriggSelNOGOODITSLAYER3},
  {"nogooditslayer0123", TriggSelNOGOODITSLAYER0123},
  {"nogooditslayerall", TriggSelNOGOODITSLAYERALL}};

/// \std::map triggerSelectionExternalNamesMap
/// \brief maps trigger selection bits to external names
static const std::map<int, std::string_view> triggerSelectionExternalNamesMap{
  {TriggSelNONE, "none"},
  {TriggSelMB, "Sel8"}, ///< Sel8 includes kIsTriggerTVX, kNoTimeFrameBorder, and kNoITSROFrameBorder
  {TriggSelNOSAMEBUNCHPUP, "No same bunch pileup"},
  {TriggSelNUMPVCONTRIBUTORS, "Number PV contributors"},
  {TriggSelVTXTOFMATCHED, "PV contributor TOF matched"},
  {TriggSelVTXTRDMATCHED, "PV contributor TRD matched"},
  {TriggSelNOCOLLINTRSTD, "No coll in TR standard"},
  {TriggSelNOCOLLINROFSTD, "No coll in ROF standard"},
  {TriggSelISVERTEXITSTPC, "Vertex from ITS and TPC"},
  {TriggSelISGOODZVTXFT0VSPV, "Good vtx FT0 vs PV"},
  {TriggSelGOODITSLAYER3, "Good ITS layer 3"},
  {TriggSelGOODITSLAYER0123, "Good ITS layers 0,1,2,3"},
  {TriggSelGOODITSLAYERALL, "Good ITS layer all"},
  {TriggSelNOGOODITSLAYER3, "No good ITS layer 3"},
  {TriggSelNOGOODITSLAYER0123, "No good ITS layer 0,1,2,3"},
  {TriggSelNOGOODITSLAYERALL, "No good ITS layer all"}};

/// \enum OccupancyEstimationType
/// \brief The type of occupancy estimation
enum OccupancyEstimationType {
  OccupancyNOOCC = 0,     ///< do not use occupancy estimation
  OccupancyTRACKSOCC,     ///< occupancy estimated using tracks
  OccupancyFT0COCC,       ///< occupancy estimated using the FT0C
  OccupancyNOOFESTIMATORS ///< the number of occupancy estimators
};

/// \enum CollisionSelectionFlags
/// \brief The different criteria for selecting/rejecting collisions
enum CollisionSelectionFlags {
  CollSelIN = 0,           ///< new unhandled, yet, event
  CollSelMBBIT,            ///< minimum  bias
  CollSelINT7BIT,          ///< INT7 Run 1/2
  CollSelSEL7BIT,          ///< Sel7 Run 1/2
  CollSelTRIGGSELBIT,      ///< Accepted by trigger selection
  CollSelRCTBIT,           ///< Accetped by the RCT information
  CollSelOCCUPANCYBIT,     ///< occupancy within limits
  CollSelCENTRALITYBIT,    ///< centrality cut passed
  CollSelZVERTEXBIT,       ///< zvtx cut passed
  CollSelMULTCORRELATIONS, ///< multiplicities correlations passed
  CollSelSELECTED,         ///< the event has passed all selections
  CollSelNOOFFLAGS         ///< number of flags
};

constexpr std::bitset<32> CollSelACCEPTEDRUN3 = BIT(CollSelTRIGGSELBIT) | BIT(CollSelRCTBIT) | BIT(CollSelOCCUPANCYBIT) | BIT(CollSelCENTRALITYBIT) | BIT(CollSelZVERTEXBIT) | BIT(CollSelMULTCORRELATIONS);
constexpr std::bitset<32> CollSelPREMULTACCEPTEDRUN3 = BIT(CollSelTRIGGSELBIT) | BIT(CollSelRCTBIT) | BIT(CollSelOCCUPANCYBIT) | BIT(CollSelCENTRALITYBIT) | BIT(CollSelZVERTEXBIT);

/// \std::mag collisionSelectionExternalNamesMap
/// \brief maps collision selection bits to external names
static const std::map<int, std::string_view> collisionSelectionExternalNamesMap{
  {CollSelIN, "In"},
  {CollSelMBBIT, "MB"},
  {CollSelINT7BIT, "INT7"},
  {CollSelSEL7BIT, "Sel7"},
  {CollSelTRIGGSELBIT, "Trigger selection"},
  {CollSelRCTBIT, "RCT accepted"},
  {CollSelOCCUPANCYBIT, "Occupancy"},
  {CollSelCENTRALITYBIT, "Centrality"},
  {CollSelZVERTEXBIT, "z vertex"},
  {CollSelMULTCORRELATIONS, "Multiplicities correlations"},
  {CollSelSELECTED, "Selected"}};

/// \enum StrongDebugging
/// \brief Enable a per track information debugging. Only for local analyses
enum StrongDebugging {
  kNODEBUG = 0, ///< do not debug
  kDEBUG        ///< output debugging information on a per track basis to a text file
};

/// \enum TpcExclusionMethod
/// \brief Methods for excluding tracks witin the TPC
enum TpcExclusionMethod {
  kNOEXCLUSION = 0, ///< do not exclude tracks within the TPC
  kSTATIC,          ///< exclude tracks statically on the bins of the TPC sector borders; only valid if 72 bins and origin shifted by 0.5
  kDYNAMIC          ///< pT dependent exclusion matching the sector borders a la Alex Dobrin
};

//============================================================================================
// The debug output stream
//============================================================================================
std::ofstream debugstream;

//============================================================================================
// The overall minimum momentum
//============================================================================================
float overallminp = 0.0f;

//============================================================================================
// The collision selection and trigger selection flags and configuration objects
//============================================================================================
std::bitset<32> collisionFlags;
std::bitset<32> collisionSelectionFlags;
std::bitset<32> triggerFlags;
std::bitset<32> triggerSelectionFlags;

//============================================================================================
// The collision exclusion using correlations between multiplicity/centrality observables
//============================================================================================
bool useCentralityMultiplicityCorrelationsExclusion = false;
std::vector<float> collisionMultiplicityCentralityObservables = {};
std::vector<int> observableIndexForCentralityMultiplicityParameter = {};
TFormula* multiplicityCentralityCorrelationsExclusion = nullptr;

//============================================================================================
// The input data metadata access helper
//============================================================================================
o2::common::core::MetadataHelper metadataInfo;

//============================================================================================
// The RCT information access
//============================================================================================
bool useRctInformation = false;
o2::aod::rctsel::RCTFlagsChecker rctChecker;

//============================================================================================
// The DptDptFilter configuration objects
//============================================================================================
int ptbins = 18;
float ptlow = 0.2, ptup = 2.0;
int etabins = 16;
float etalow = -0.8, etaup = 0.8;
int zvtxbins = 40;
float zvtxlow = -10.0, zvtxup = 10.0;
int phibins = 72;
float philow = 0.0f;
float phiup = constants::math::TwoPI;
float phibinshift = 0.0f;

struct TpcExcludeTrack;             ///< forward declaration of the excluder object
bool onlyInOneSide = false;         ///< select only tracks that don't cross the TPC central membrane
extern TpcExcludeTrack tpcExcluder; ///< the TPC excluder object instance

/* selection criteria from PWGMM */
static constexpr int TrackTypePWGMM = 4;
// default quality criteria for tracks with ITS contribution
static constexpr o2::aod::track::TrackSelectionFlags::flagtype TrackSelectionITS =
  o2::aod::track::TrackSelectionFlags::kITSNCls | o2::aod::track::TrackSelectionFlags::kITSChi2NDF |
  o2::aod::track::TrackSelectionFlags::kITSHits;

// default quality criteria for tracks with TPC contribution
static constexpr o2::aod::track::TrackSelectionFlags::flagtype TrackSelectionTPC =
  o2::aod::track::TrackSelectionFlags::kTPCNCls |
  o2::aod::track::TrackSelectionFlags::kTPCCrossedRowsOverNCls |
  o2::aod::track::TrackSelectionFlags::kTPCChi2NDF;

// default standard DCA cuts
static constexpr o2::aod::track::TrackSelectionFlags::flagtype TrackSelectionDCA =
  o2::aod::track::TrackSelectionFlags::kDCAz | o2::aod::track::TrackSelectionFlags::kDCAxy;

struct DptDptTrackSelection; // forward struct declaration
int tracktype = 1;
std::vector<DptDptTrackSelection*> trackFilters = {}; // the vector of track selectors

struct DptDptTrackSelection {
  DptDptTrackSelection(TrackSelection* stdTs, TList* outputList, const char* name) : stdTrackSelection(stdTs)
  {
    passedHistogram = new TH1F(name, name, ptbins, ptlow, ptup);
    outputList->Add(passedHistogram);
  }
  DptDptTrackSelection(TrackSelection* stdTs, std::function<float(float)> ptDepCut, TList* outputList, const char* name)
    : stdTrackSelection(stdTs),
      maxDcazPtDep(ptDepCut)
  {
    passedHistogram = new TH1F(name, name, ptbins, ptlow, ptup);
    outputList->Add(passedHistogram);
  }
  void setMaxDcaXY(float max)
  {
    maxDCAxy = max;
    stdTrackSelection->SetMaxDcaXY(max);
  }
  void setMaxDcaZ(float max)
  {
    maxDCAz = max;
    stdTrackSelection->SetMaxDcaZ(max);
  }
  void setMaxDcazPtDep(std::function<float(float)> ptDepCut)
  {
    maxDcazPtDep = ptDepCut;
  }
  void setRequirePvContributor(bool pvc = true)
  {
    requirePvContributor = pvc;
  }

  template <typename TrackObject>
  bool isSelected(TrackObject const& track) const
  {
    if (stdTrackSelection->IsSelected(track)) {
      auto checkDca2Dcut = [&](auto const& track) {
        if (dca2Dcut) {
          if (track.dcaXY() * track.dcaXY() / maxDCAxy / maxDCAxy + track.dcaZ() * track.dcaZ() / maxDCAz / maxDCAz > 1) {
            return false;
          } else {
            return true;
          }
        } else {
          return true;
        }
      };
      auto checkDcaZcut = [&](auto const& track) {
        return ((maxDcazPtDep) ? std::fabs(track.dcaZ()) <= maxDcazPtDep(track.pt()) : true);
      };

      /* tight pT dependent DCAz cut */
      if (!checkDcaZcut(track)) {
        return false;
      }
      /* 2D DCA xy-o-z cut */
      if (!checkDca2Dcut(track)) {
        return false;
      }
      /* primary vertex contributor */
      if (requirePvContributor) {
        if (!track.isPVContributor()) {
          return false;
        }
      }
      passedHistogram->Fill(track.pt());
      return true;
    } else {
      return false;
    }
  }

  static void initializeTrackSelection(TrackSelectionTuneCfg& tune, TList* outputList)
  {
    auto addTrackFilter = [](auto filter) {
      trackFilters.push_back(filter);
    };
    auto highQualityTpcTrack = [](TList* outList, const char* name) {
      DptDptTrackSelection* tpcTrack = new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelection()), outList, name);
      tpcTrack->stdTrackSelection->ResetITSRequirements();
      tpcTrack->stdTrackSelection->SetRequireITSRefit(false);
      tpcTrack->stdTrackSelection->SetMinNClustersTPC(120);
      tpcTrack->stdTrackSelection->SetMaxTPCFractionSharedCls(0.2f);
      return tpcTrack;
    };
    auto highQualityItsOnlyTrack = [](TList* outList, const char* name) {
      DptDptTrackSelection* itsTrack = new DptDptTrackSelection(new TrackSelection(), [](float pt) { return 0.004f + 0.013f / pt; }, outList, name);
      itsTrack->stdTrackSelection->SetTrackType(o2::aod::track::TrackTypeEnum::Track);
      itsTrack->stdTrackSelection->SetRequireITSRefit(true);
      itsTrack->stdTrackSelection->SetRequireHitsInITSLayers(2, {0, 1, 2});
      itsTrack->stdTrackSelection->SetMaxChi2PerClusterITS(36.0f);
      itsTrack->stdTrackSelection->SetMaxDcaXYPtDep([](float pt) { return 0.004f + 0.013f / pt; });
      return itsTrack;
    };
    switch (tracktype) {
      case 1: { /* Run2 global track */
        DptDptTrackSelection* globalRun2 = new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelection()), outputList, "TType1Global");
        globalRun2->stdTrackSelection->SetTrackType(o2::aod::track::Run2Track); // Run 2 track asked by default
        globalRun2->stdTrackSelection->SetMaxChi2PerClusterTPC(2.5f);
        DptDptTrackSelection* globalSDDRun2 = new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelectionSDD()), outputList, "TType1Sdd");
        globalSDDRun2->stdTrackSelection->SetTrackType(o2::aod::track::Run2Track); // Run 2 track asked by default
        globalSDDRun2->stdTrackSelection->SetMaxChi2PerClusterTPC(2.5f);
        addTrackFilter(globalRun2);
        addTrackFilter(globalSDDRun2);
      } break;
      case 3: { /* Run3 track */
        DptDptTrackSelection* globalRun3 = new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelection()), outputList, "TType3Global");
        globalRun3->stdTrackSelection->SetTrackType(o2::aod::track::TrackTypeEnum::Track);
        globalRun3->stdTrackSelection->ResetITSRequirements();
        globalRun3->stdTrackSelection->SetRequireHitsInITSLayers(1, {0, 1, 2});
        DptDptTrackSelection* globalSDDRun3 = new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelection()), outputList, "TType3Sdd");
        globalSDDRun3->stdTrackSelection->SetTrackType(o2::aod::track::TrackTypeEnum::Track);
        globalSDDRun3->stdTrackSelection->ResetITSRequirements();
        globalSDDRun3->stdTrackSelection->SetRequireNoHitsInITSLayers({0, 1, 2});
        globalSDDRun3->stdTrackSelection->SetRequireHitsInITSLayers(1, {3});
        addTrackFilter(globalRun3);
        addTrackFilter(globalSDDRun3);
      } break;
      case 5: { /* Run2 TPC only track */
        DptDptTrackSelection* tpcOnly = new DptDptTrackSelection(new TrackSelection, outputList, "TType5");
        tpcOnly->stdTrackSelection->SetTrackType(o2::aod::track::Run2Track); // Run 2 track asked by default
        tpcOnly->stdTrackSelection->SetMinNClustersTPC(50);
        tpcOnly->stdTrackSelection->SetMaxChi2PerClusterTPC(4);
        tpcOnly->setMaxDcaZ(3.2f);
        tpcOnly->setMaxDcaXY(2.4f);
        tpcOnly->dca2Dcut = true;
        addTrackFilter(tpcOnly);
      } break;
      case 7: { /* Run3 TPC only track */
        DptDptTrackSelection* tpcOnly = new DptDptTrackSelection(new TrackSelection, outputList, "TType7");
        tpcOnly->stdTrackSelection->SetTrackType(o2::aod::track::TrackTypeEnum::Track);
        tpcOnly->stdTrackSelection->SetMinNClustersTPC(50);
        tpcOnly->stdTrackSelection->SetMaxChi2PerClusterTPC(4);
        tpcOnly->setMaxDcaZ(3.2f);
        tpcOnly->setMaxDcaXY(2.4f);
        tpcOnly->dca2Dcut = true;
        addTrackFilter(tpcOnly);
      } break;
      case 10: { /* Run3 track primary vertex contributor */
        DptDptTrackSelection* globalRun3 = new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelection()), outputList, "TType10Global");
        globalRun3->stdTrackSelection->SetTrackType(o2::aod::track::TrackTypeEnum::Track);
        globalRun3->stdTrackSelection->ResetITSRequirements();
        globalRun3->stdTrackSelection->SetRequireHitsInITSLayers(1, {0, 1, 2});
        globalRun3->setRequirePvContributor(true);
        DptDptTrackSelection* globalSDDRun3 = new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelection()), outputList, "TType10Sdd");
        globalSDDRun3->stdTrackSelection->SetTrackType(o2::aod::track::TrackTypeEnum::Track);
        globalSDDRun3->stdTrackSelection->ResetITSRequirements();
        globalSDDRun3->stdTrackSelection->SetRequireNoHitsInITSLayers({0, 1, 2});
        globalSDDRun3->stdTrackSelection->SetRequireHitsInITSLayers(1, {3});
        globalSDDRun3->setRequirePvContributor(true);
        addTrackFilter(globalRun3);
        addTrackFilter(globalSDDRun3);
      } break;
      case 30: { /* Run 3 default global track: kAny on 3 IB layers of ITS */
        addTrackFilter(new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny, TrackSelection::GlobalTrackRun3DCAxyCut::Default)), outputList, "TType30"));
      } break;
      case 31: { /* Run 3 global track: kTwo on 3 IB layers of ITS */
        addTrackFilter(new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibTwo, TrackSelection::GlobalTrackRun3DCAxyCut::Default)), outputList, "TType31"));
      } break;
      case 32: { /* Run 3 global track: kAny on all 7 layers of ITS */
        addTrackFilter(new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSallAny, TrackSelection::GlobalTrackRun3DCAxyCut::Default)), outputList, "TType32"));
      } break;
      case 33: { /* Run 3 global track: kAll on all 7 layers of ITS */
        addTrackFilter(new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSall7Layers, TrackSelection::GlobalTrackRun3DCAxyCut::Default)), outputList, "TType33"));
      } break;
      case 40: { /* Run 3 global track: kAny on 3 IB layers of ITS, tighter DCAxy */
        addTrackFilter(new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)), outputList, "TType40"));
      } break;
      case 41: { /* Run 3 global track: kTwo on 3 IB layers of ITS, tighter DCAxy */
        addTrackFilter(new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibTwo, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)), outputList, "TType41"));
      } break;
      case 42: { /* Run 3 global track: kAny on all 7 layers of ITS, tighter DCAxy */
        addTrackFilter(new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSallAny, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)), outputList, "TType42"));
      } break;
      case 43: { /* Run 3 global track: kAll on all 7 layers of ITS, tighter DCAxy */
        addTrackFilter(new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSall7Layers, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)), outputList, "TType43"));
      } break;
      case 50: { /* Run 3 global track: kAny on 3 IB layers of ITS, tighter DCAxy, tighter pT dep DCAz */
        addTrackFilter(new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)), [](float pt) { return 0.004f + 0.013f / pt; }, outputList, "TType50"));
      } break;
      case 51: { /* Run 3 global track: kTwo on 3 IB layers of ITS, tighter DCAxy, tighter pT dep DCAz */
        addTrackFilter(new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibTwo, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)), [](float pt) { return 0.004f + 0.013f / pt; }, outputList, "TType51"));
      } break;
      case 52: { /* Run 3 global track: kAny on all 7 layers of ITS, tighter DCAxy, tighter pT dep DCAz */
        addTrackFilter(new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSallAny, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)), [](float pt) { return 0.004f + 0.013f / pt; }, outputList, "TType52"));
      } break;
      case 53: { /* Run 3 global track: kAll on all 7 layers of ITS, tighter DCAxy, tighter pT dep DCAz */
        addTrackFilter(new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSall7Layers, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)), [](float pt) { return 0.004f + 0.013f / pt; }, outputList, "TType53"));
      } break;
      case 60: { /* Run 3 global track: kAny on 3 IB layers of ITS, tighter DCAxy, tighter pT dep DCAz, plus TPC+TOF only tracks */
        addTrackFilter(new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)), [](float pt) { return 0.004f + 0.013f / pt; }, outputList, "TType60Global"));
        addTrackFilter(highQualityTpcTrack(outputList, "TType60Tpc"));
      } break;
      case 61: { /* Run 3 global track: kTwo on 3 IB layers of ITS, tighter DCAxy, tighter pT dep DCAz, plus TPC+TOF only tracks */
        addTrackFilter(new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibTwo, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)), [](float pt) { return 0.004f + 0.013f / pt; }, outputList, "TType61Global"));
        addTrackFilter(highQualityTpcTrack(outputList, "TType61Tpc"));
      } break;
      case 62: { /* Run 3 global track: kAny on all 7 layers of ITS, tighter DCAxy, tighter pT dep DCAz, plus TPC+TOF only tracks */
        addTrackFilter(new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSallAny, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)), [](float pt) { return 0.004f + 0.013f / pt; }, outputList, "TType62Global"));
        addTrackFilter(highQualityTpcTrack(outputList, "TType62Tpc"));
      } break;
      case 63: { /* Run 3 global track: kAll on all 7 layers of ITS, tighter DCAxy, tighter pT dep DCAz, plus TPC+TOF only tracks */
        addTrackFilter(new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSall7Layers, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)), [](float pt) { return 0.004f + 0.013f / pt; }, outputList, "TType63Global"));
        addTrackFilter(highQualityTpcTrack(outputList, "TType63Tpc"));
      } break;
      case 70: { /* Run 3 global track: kAny on 3 IB layers of ITS, tighter DCAxy, tighter pT dep DCAz, plus ITS only tracks */
        addTrackFilter(new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibAny, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)), [](float pt) { return 0.004f + 0.013f / pt; }, outputList, "TType70Global"));
        addTrackFilter(highQualityItsOnlyTrack(outputList, "TType70Its"));
      } break;
      case 71: { /* Run 3 global track: kTwo on 3 IB layers of ITS, tighter DCAxy, tighter pT dep DCAz, plus ITS only tracks */
        addTrackFilter(new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSibTwo, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)), [](float pt) { return 0.004f + 0.013f / pt; }, outputList, "TType71Global"));
        addTrackFilter(highQualityItsOnlyTrack(outputList, "TType71Its"));
      } break;
      case 72: { /* Run 3 global track: kAny on all 7 layers of ITS, tighter DCAxy, tighter pT dep DCAz, plus ITS only tracks */
        addTrackFilter(new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSallAny, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)), [](float pt) { return 0.004f + 0.013f / pt; }, outputList, "TType72Global"));
        addTrackFilter(highQualityItsOnlyTrack(outputList, "TType72Its"));
      } break;
      case 73: { /* Run 3 global track: kAll on all 7 layers of ITS, tighter DCAxy, tighter pT dep DCAz, plus ITS only tracks */
        addTrackFilter(new DptDptTrackSelection(new TrackSelection(getGlobalTrackSelectionRun3ITSMatch(TrackSelection::GlobalTrackRun3ITSMatching::Run3ITSall7Layers, TrackSelection::GlobalTrackRun3DCAxyCut::ppPass3)), [](float pt) { return 0.004f + 0.013f / pt; }, outputList, "TType73Global"));
        addTrackFilter(highQualityItsOnlyTrack(outputList, "TType73Its"));
      } break;
      default:
        break;
    }
    if (tune.mUseIt) {
      for (auto const& filter : trackFilters) {
        if (tune.mUseITSclusters) {
          filter->stdTrackSelection->ResetITSRequirements();
          filter->stdTrackSelection->SetRequireHitsInITSLayers(1, {0, 1, 2});
          filter->stdTrackSelection->SetMinNClustersITS(tune.mITSclusters);
        }
        if (tune.mUseTPCclusters) {
          filter->stdTrackSelection->SetMinNClustersTPC(tune.mTPCclusters);
        }
        if (tune.mUseTPCxRows) {
          filter->stdTrackSelection->SetMinNCrossedRowsTPC(tune.mTPCxRows);
        }
        if (tune.mUseTPCXRoFClusters) {
          filter->stdTrackSelection->SetMinNCrossedRowsOverFindableClustersTPC(tune.mTPCXRoFClusters);
        }
        if (tune.mUseDCAxy) {
          /* DCAxy is tricky due to how the pT dependence is implemented */
          filter->stdTrackSelection->SetMaxDcaXYPtDep([&tune](float) { return tune.mDCAxy; });
          filter->setMaxDcaXY(tune.mDCAxy);
        }
        if (tune.mUseDCAz) {
          /* DCAz is tricky due to how the pT dependence is implemented */
          filter->setMaxDcazPtDep([&tune](float) { return tune.mDCAz; });
          filter->setMaxDcaZ(tune.mDCAz);
        }
        if (tune.mUseFractionTpcSharedClusters) {
          filter->stdTrackSelection->SetMaxTPCFractionSharedCls(tune.mFractionTpcSharedClusters);
        }
      }
    }
  }

  float maxDCAxy = 1e6;
  float maxDCAz = 1e6;
  TrackSelection* stdTrackSelection = nullptr;
  std::function<float(float)> maxDcazPtDep = {};
  TH1* passedHistogram = nullptr;
  bool dca2Dcut = false;
  bool requirePvContributor = false;
};

inline TList* getCCDBInput(auto& ccdb, const char* ccdbpath, const char* ccdbdate, bool periodInPath = false, const std::string& suffix = "")
{
  std::tm cfgtm = {};
  std::stringstream ss(ccdbdate);
  ss >> std::get_time(&cfgtm, "%Y%m%d");
  cfgtm.tm_hour = 12;
  int64_t timestamp = std::mktime(&cfgtm) * 1000;

  auto cleanPeriod = [](const auto& str) {
    std::string tmpStr = str;
    size_t pos = tmpStr.find('_');
    if (pos != std::string::npos) {
      tmpStr.erase(pos);
    }
    return tmpStr;
  };

  std::string actualPeriod = cleanPeriod(metadataInfo.get("LPMProductionTag"));
  std::string actualPath = ccdbpath;
  if (periodInPath) {
    actualPath = actualPath + "/" + actualPeriod;
  }
  if (suffix.length() > 0) {
    actualPeriod = actualPeriod + "_" + suffix;
  }

  TList* lst = nullptr;
  std::map<std::string, std::string> metadata{{"Period", actualPeriod}};
  lst = ccdb->template getSpecific<TList>(actualPath, timestamp, metadata);
  if (lst != nullptr) {
    LOGF(info, "Correctly loaded CCDB input object");
  } else {
    LOGF(error, "CCDB input object could not be loaded");
  }
  return lst;
}

SystemType fSystem = SystemNoSystem;
MultRunType fLhcRun = MultRunRUN1RUN2;
DataType fDataType = kData;
CentMultEstimatorType fCentMultEstimator = CentMultV0M;
OccupancyEstimationType fOccupancyEstimation = OccupancyNOOCC; /* the occupancy estimator to use */

float fMinOccupancy = 0.0f; /* the minimum allowed occupancy */
float fMaxOccupancy = 1e6f; /* the maximum allowed occupancy */

/* adaptations for the pp nightly checks */
analysis::CheckRangeCfg traceDCAOutliers;
bool traceOutOfSpeciesParticles = false;
int recoIdMethod = 0;
float particleMaxDCAxy = 999.9f;
float particleMaxDCAZ = 999.9f;
bool traceCollId0 = false;

inline std::bitset<32> getTriggerSelection(std::string_view const& triggstr)
{
  std::bitset<32> flags;

  auto split = [](const auto s) {
    std::vector<std::string_view> tokens;
    std::string_view token;

    size_t posStart = 0;
    size_t posEnd;
    while ((posEnd = s.find("+", posStart)) != std::string_view::npos) {
      token = s.substr(posStart, posEnd - posStart);
      posStart = posEnd + 1;
      tokens.push_back(token);
    }
    tokens.push_back(s.substr(posStart));
    return tokens;
  };

  std::vector<std::string_view> tags = split(triggstr);

  for (const auto& tag : tags) {
    if (triggerSelectionBitsMap.contains(tag)) {
      flags.set(triggerSelectionBitsMap.at(tag), true);
    } else {
      LOGF(fatal, "Wrong trigger selection tag: %s", tag.data());
    }
  }
  return flags;
}

inline SystemType getSystemType(auto const& periodsForSysType)
{
  auto period = metadataInfo.get("LPMProductionTag");
  auto anchoredPeriod = metadataInfo.get("AnchorProduction");
  bool checkAnchor = anchoredPeriod.length() > 0;

  for (SystemType sT = SystemNoSystem; sT < SystemNoOfSystems; ++sT) {
    const std::string& periods = periodsForSysType[static_cast<int>(sT)][0];
    auto contains = [periods](auto const& period) {
      if (periods.find(period) != std::string::npos) {
        return true;
      }
      return false;
    };
    if (periods.length() > 0) {
      if (contains(period) || (checkAnchor && contains(anchoredPeriod))) {
        LOGF(info, "DptDptCorrelations::getSystemType(). Assigned system type %s for period %s", systemExternalNamesMap.at(static_cast<int>(sT)).data(), period.c_str());
        return sT;
      }
    }
  }
  LOGF(fatal, "DptDptCorrelations::getSystemType(). No system type for period: %s", period.c_str());
  return SystemPbPb;
}

/// \brief Type of data according to the configuration string
/// \param datastr The data type configuration string
/// \return Internal code for the passed kind of data string
inline DataType getDataType(std::string const& datastr)
{
  /* we have to figure out how extract the type of data*/
  if (datastr.empty() || (datastr == "data")) {
    return kData;
  } else if (datastr == "datanoevsel") {
    return kDataNoEvtSel;
  } else if (datastr == "MC") {
    return kMC;
  } else if (datastr == "FastMC") {
    return kFastMC;
  } else if (datastr == "OnTheFlyMC") {
    return kOnTheFly;
  } else {
    LOGF(fatal, "DptDptCorrelations::getDataType(). Wrong type of dat: %d", datastr.c_str());
  }
  return kData;
}

inline CentMultEstimatorType getCentMultEstimator(std::string_view const& datastr)
{
  if (estimatorInternalCodesMap.contains(datastr)) {
    return static_cast<CentMultEstimatorType>(estimatorInternalCodesMap.at(datastr));
  } else {
    LOGF(fatal, "Centrality/Multiplicity estimator %s not supported yet", datastr.data());
  }
  return CentMultNOCM;
}

inline std::string_view getCentMultEstimatorName(CentMultEstimatorType est)
{
  if (estimatorExternalNamesMap.contains(est)) {
    return estimatorExternalNamesMap.at(est);
  } else {
    LOGF(fatal, "Centrality/Multiplicity estimator %d not supported yet", static_cast<int>(est));
  }
  return "WRONG";
}

inline OccupancyEstimationType getOccupancyEstimator(const std::string_view& estimator)
{
  if (estimator == "None") {
    return OccupancyNOOCC;
  } else if (estimator == "Tracks") {
    return OccupancyTRACKSOCC;
  } else if (estimator == "FT0C") {
    return OccupancyFT0COCC;
  } else {
    LOGF(fatal, "Occupancy estimator %s not supported yet", estimator.data());
    return OccupancyNOOCC;
  }
}

/// @brief  gets the exclusion formula from the corresponding expression and initializes the exclusion machinery
/// @param  std::string_view with the formula expression
/// @return the expression TFormula
inline TFormula* getExclusionFormula(std::string_view formula)
{
  collisionMultiplicityCentralityObservables.resize(CentMultCorrelationsNOOFPARAMS);
  if (formula.length() != 0) {
    useCentralityMultiplicityCorrelationsExclusion = true;
    TFormula* f = new TFormula("Exclussion expression", formula.data());
    int nParameters = f->GetNpar();
    observableIndexForCentralityMultiplicityParameter.resize(nParameters);
    LOGF(info, "Configuring outliers exclusion with the formula %s which has %d parameters", formula.data(), nParameters);
    for (int iPar = 0; iPar < nParameters; ++iPar) {
      if (centMultCorrelationsParamsMap.contains(std::string(f->GetParName(iPar)))) {
        observableIndexForCentralityMultiplicityParameter[iPar] = centMultCorrelationsParamsMap.at(std::string(f->GetParName(iPar)));
        LOGF(info, "\tAssigned observable %s with index %d to the parameter %s with parameter index %d", centMultCorrelationsParamsNamesMap.at(centMultCorrelationsParamsMap.at(std::string(f->GetParName(iPar)))).data(), static_cast<int>(centMultCorrelationsParamsMap.at(std::string(f->GetParName(iPar)))), f->GetParName(iPar), iPar);
      } else {
        LOGF(fatal, "Exclusion expression contains parameter %s which is still not supported. Please, fix it!", f->GetParName(iPar));
      }
    }
    return f;
  } else {
    useCentralityMultiplicityCorrelationsExclusion = false;
    return nullptr;
  }
}

//////////////////////////////////////////////////////////////////////////////////
/// Trigger selection
//////////////////////////////////////////////////////////////////////////////////

/// \brief Trigger selection for reconstructed and detector level collision tables
template <typename CollisionObject>
inline bool triggerSelectionReco(CollisionObject const& collision)
{
  bool trigsel = false;
  switch (fSystem) {
    case SystemPp:
    case SystemPPb:
    case SystemPbp:
    case SystemPbPb:
    case SystemXeXe:
      if (triggerSelectionFlags.test(TriggSelMB)) {
        switch (fDataType) {
          case kData:
            if (collision.alias_bit(kINT7)) {
              collisionFlags.set(CollSelINT7BIT);
              if (collision.sel7()) {
                trigsel = true;
                collisionFlags.set(CollSelSEL7BIT);
              }
            }
            break;
          case kMC:
            if (collision.sel7()) {
              trigsel = true;
              collisionFlags.set(CollSelSEL7BIT);
            }
            break;
          default:
            collisionFlags.set(CollSelMBBIT);
            trigsel = true;
            break;
        }
      } else if (triggerSelectionFlags.test(TriggSelNONE)) {
        trigsel = true;
      }
      break;
    case SystemPpRun3:
    case SystemPbPbRun3:
    case SystemNeNeRun3:
    case SystemPORun3:
    case SystemOORun3: {
      auto setTriggerFlags = [](auto& flags, auto const& coll) {
        flags.set(TriggSelMB, coll.sel8() != 0);
        flags.set(TriggSelNOSAMEBUNCHPUP, coll.selection_bit(aod::evsel::kNoSameBunchPileup));
        flags.set(TriggSelNUMPVCONTRIBUTORS, coll.numContrib() > 1); /* TODO: make it configurable */
        flags.set(TriggSelVTXTOFMATCHED, coll.selection_bit(aod::evsel::kIsVertexTOFmatched));
        flags.set(TriggSelVTXTRDMATCHED, coll.selection_bit(aod::evsel::kIsVertexTRDmatched));
        flags.set(TriggSelNOCOLLINTRSTD, coll.selection_bit(aod::evsel::kNoCollInTimeRangeStandard));
        flags.set(TriggSelNOCOLLINROFSTD, coll.selection_bit(aod::evsel::kNoCollInRofStandard));
        flags.set(TriggSelISVERTEXITSTPC, coll.selection_bit(aod::evsel::kIsVertexITSTPC));
        flags.set(TriggSelISGOODZVTXFT0VSPV, coll.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV));
        flags.set(TriggSelGOODITSLAYER3, coll.selection_bit(aod::evsel::kIsGoodITSLayer3));
        flags.set(TriggSelGOODITSLAYER0123, coll.selection_bit(aod::evsel::kIsGoodITSLayer0123));
        flags.set(TriggSelGOODITSLAYERALL, coll.selection_bit(aod::evsel::kIsGoodITSLayersAll));
        flags.set(TriggSelNOGOODITSLAYER3, !coll.selection_bit(aod::evsel::kIsGoodITSLayer3));
        flags.set(TriggSelNOGOODITSLAYER0123, !coll.selection_bit(aod::evsel::kIsGoodITSLayer0123));
        flags.set(TriggSelNOGOODITSLAYERALL, !coll.selection_bit(aod::evsel::kIsGoodITSLayersAll));
      };
      setTriggerFlags(triggerFlags, collision);

      /* special treatment for no trigger selection */
      if (triggerSelectionFlags.test(TriggSelNONE)) {
        trigsel = true;
      } else {
        trigsel = ((triggerSelectionFlags & triggerFlags) == triggerSelectionFlags);
      }
    } break;
    default:
      break;
  }
  collisionFlags.set(CollSelTRIGGSELBIT, trigsel);
  return trigsel;
}

/// \brief Trigger selection by default: unknow subscribed collision table
template <typename CollisionObject>
inline bool triggerSelection(CollisionObject const&)
{
  LOGF(fatal, "Trigger selection not implemented for this kind of collisions");
  return false;
}

/// \brief Trigger selection for reconstructed collision tables without centrality/multiplicity
template <>
inline bool triggerSelection<aod::CollisionEvSel>(aod::CollisionEvSel const& collision)
{
  return triggerSelectionReco(collision);
}

/// \brief Trigger selection for reconstructed collision tables with Run 2 centrality/multiplicity
template <>
inline bool triggerSelection<aod::CollisionEvSelRun2Cent>(aod::CollisionEvSelRun2Cent const& collision)
{
  return triggerSelectionReco(collision);
}

/// \brief Trigger selection for reconstructed collision tables with centrality/multiplicity
template <>
inline bool triggerSelection<aod::CollisionEvSelCent>(aod::CollisionEvSelCent const& collision)
{
  return triggerSelectionReco(collision);
}

/// \brief Trigger selection for detector level collision tables without centrality/multiplicity
template <>
inline bool triggerSelection<soa::Join<aod::CollisionsEvSel, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSel, aod::McCollisionLabels>::iterator const& collision)
{
  return triggerSelectionReco(collision);
}

/// \brief Trigger selection for detector level collision tables with Run 2 centrality/multiplicity
template <>
inline bool triggerSelection<soa::Join<aod::CollisionsEvSelRun2Cent, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSelRun2Cent, aod::McCollisionLabels>::iterator const& collision)
{
  return triggerSelectionReco(collision);
}

/// \brief Trigger selection for detector level collision tables with centrality/multiplicity
template <>
inline bool triggerSelection<soa::Join<aod::CollisionsEvSelCent, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSelCent, aod::McCollisionLabels>::iterator const& collision)
{
  return triggerSelectionReco(collision);
}

/// \brief Trigger selection for generator level collison table
template <>
inline bool triggerSelection<aod::McCollision>(aod::McCollision const&)
{
  return true;
}

//////////////////////////////////////////////////////////////////////////////////
/// Multiplicity extraction
//////////////////////////////////////////////////////////////////////////////////
static constexpr float ValidPercentileLowLimit = 0.0f;
static constexpr float ValidPercentileUpLimit = 100.0f;

/// \brief Extract the collision multiplicity from the event selection information
template <typename CollisionObject>
inline float extractMultiplicity(CollisionObject const& collision, CentMultEstimatorType est)
{
  switch (est) {
    case CentMultV0M:
      return collision.multFV0M();
      break;
    case CentMultCL0:
      return collision.multTracklets();
      break;
    case CentMultCL1:
      return collision.multTracklets();
      break;
    case CentMultFV0A:
      return collision.multFV0A();
      break;
    case CentMultFT0M:
      return collision.multFT0M();
      break;
    case CentMultFT0A:
      return collision.multFT0A();
      break;
    case CentMultFT0C:
      return collision.multFT0C();
      break;
    case CentMultNTPV:
      return collision.multNTracksPV();
      break;
    case CentMultNOCM:
      return collision.multFT0M();
      break;
    default:
      LOGF(fatal, "Centrality/Multiplicity estimator %d not supported yet", (int)est);
      return collision.multFT0M();
      break;
  }
}

//////////////////////////////////////////////////////////////////////////////////
/// Centrality selection
//////////////////////////////////////////////////////////////////////////////////

/// \brief Centrality/multiplicity percentile
template <typename CollisionObject>
  requires(o2::aod::HasRun2Centrality<CollisionObject>)
float getCentMultPercentile(CollisionObject collision)
{
  switch (fCentMultEstimator) {
    case CentMultV0M:
      return collision.centRun2V0M();
    case CentMultCL0:
      return collision.centRun2CL0();
    case CentMultCL1:
      return collision.centRun2CL1();
    default:
      return 105.0;
  }
}

template <typename CollisionObject>
  requires(o2::aod::HasCentrality<CollisionObject>)
float getCentMultPercentile(CollisionObject collision)
{
  switch (fCentMultEstimator) {
    case CentMultFV0A:
      return collision.centFV0A();
    case CentMultFT0M:
      return collision.centFT0M();
    case CentMultFT0A:
      return collision.centFT0A();
    case CentMultFT0C:
      return collision.centFT0C();
    case CentMultNTPV:
      return collision.centNTPV();
    default:
      return 105.0;
  }
}

/// \brief Centrality selection when there is centrality/multiplicity information
template <typename CollisionObject>
inline bool centralitySelectionMult(CollisionObject collision, float& centmult)
{
  float mult = getCentMultPercentile(collision);
  if (mult < ValidPercentileUpLimit && ValidPercentileLowLimit < mult) {
    centmult = mult;
    collisionFlags.set(CollSelCENTRALITYBIT);
    return true;
  }
  return false;
}

/// \brief Centrality selection when there is not centrality/multiplicity information
template <typename CollisionObject>
inline bool centralitySelectionNoMult(CollisionObject const&, float& centmult)
{
  bool centmultsel = false;
  switch (fCentMultEstimator) {
    case CentMultNOCM:
      centmult = 50.0;
      centmultsel = true;
      collisionFlags.set(CollSelCENTRALITYBIT);
      break;
    default:
      break;
  }
  return centmultsel;
}

/// \brief Centrality selection by default: unknown subscribed collision table
template <typename CollisionObject>
inline bool centralitySelection(CollisionObject const&, float&)
{
  LOGF(fatal, "Centrality selection not implemented for this kind of collisions");
  return false;
}

/// \brief Centrality selection for reconstructed and detector level collision tables with centrality/multiplicity information
template <>
inline bool centralitySelection<aod::CollisionEvSelCent>(aod::CollisionEvSelCent const& collision, float& centmult)
{
  return centralitySelectionMult(collision, centmult);
}

/// \brief Centrality selection for reconstructed and detector level collision tables with Run 2 centrality/multiplicity information
template <>
inline bool centralitySelection<aod::CollisionEvSelRun2Cent>(aod::CollisionEvSelRun2Cent const& collision, float& centmult)
{
  return centralitySelectionMult(collision, centmult);
}

/// \brief Centrality selection for reconstructed and detector level collision tables without centrality/multiplicity information
template <>
inline bool centralitySelection<aod::CollisionEvSel>(aod::CollisionEvSel const& collision, float& centmult)
{
  return centralitySelectionNoMult(collision, centmult);
}

/// \brief Centrality selection for detector level collision tables without centrality/multiplicity
template <>
inline bool centralitySelection<soa::Join<aod::CollisionsEvSel, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSel, aod::McCollisionLabels>::iterator const& collision, float& centmult)
{
  return centralitySelectionNoMult(collision, centmult);
}

/// \brief Centrality selection for detector level collision tables with centrality/multiplicity
template <>
inline bool centralitySelection<soa::Join<aod::CollisionsEvSelCent, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSelCent, aod::McCollisionLabels>::iterator const& collision, float& centmult)
{
  return centralitySelectionMult(collision, centmult);
}

/// \brief Centrality selection for detector level collision tables with Run 2 centrality/multiplicity
template <>
inline bool centralitySelection<soa::Join<aod::CollisionsEvSelRun2Cent, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSelRun2Cent, aod::McCollisionLabels>::iterator const& collision, float& centmult)
{
  return centralitySelectionMult(collision, centmult);
}

/// \brief Centrality selection for generator level collision table
template <>
inline bool centralitySelection<aod::McCollision>(aod::McCollision const&, float& centmult)
{
  if (centmult < ValidPercentileUpLimit && ValidPercentileLowLimit < centmult) {
    return true;
  } else {
    return false;
  }
}

/// @brief evalues the exclusion formula for the current multiplicity / centrality parameters
/// @return true if the collision is not excluded according to the formula false otherwise
/// WARNING: it has always to be called after filling the multiplicity / centrality observables
/// NOTE: for MC generator level does nothing, as expected
template <typename CollisionObject>
inline bool isCollisionNotExcluded()
{
  bool notExcluded = true;
  if constexpr (!framework::has_type_v<aod::mccollision::GeneratorsID, typename CollisionObject::all_columns>) {
    if (useCentralityMultiplicityCorrelationsExclusion) {
      [&]() {
        /* set the formula parameter values */
        for (size_t iPar = 0; iPar < observableIndexForCentralityMultiplicityParameter.size(); ++iPar) {
          multiplicityCentralityCorrelationsExclusion->SetParameter(iPar, collisionMultiplicityCentralityObservables[observableIndexForCentralityMultiplicityParameter[iPar]]);
        }
      }();
      notExcluded = multiplicityCentralityCorrelationsExclusion->Eval() != 0;
    }
  }
  collisionFlags.set(CollSelMULTCORRELATIONS, notExcluded);
  return notExcluded;
}

//////////////////////////////////////////////////////////////////////////////////
/// Occupancy selection
//////////////////////////////////////////////////////////////////////////////////

/// \brief select on the collision occupancy
/// \return true if collison passes the occupancy cut false otherwise
template <typename CollisionObject>
inline bool selectOnOccupancy(CollisionObject collision)
{
  switch (fOccupancyEstimation) {
    case OccupancyNOOCC:
      collisionFlags.set(CollSelOCCUPANCYBIT);
      return true;
    case OccupancyTRACKSOCC:
      if ((fMinOccupancy <= collision.trackOccupancyInTimeRange()) && (collision.trackOccupancyInTimeRange() < fMaxOccupancy)) {
        collisionFlags.set(CollSelOCCUPANCYBIT);
        return true;
      } else {
        return false;
      }
    case OccupancyFT0COCC:
      if ((fMinOccupancy <= collision.ft0cOccupancyInTimeRange()) && (collision.ft0cOccupancyInTimeRange() < fMaxOccupancy)) {
        collisionFlags.set(CollSelOCCUPANCYBIT);
        return true;
      } else {
        return false;
      }
    default:
      return false;
  }
}

/// \brief Occupancy selection by default: unknown subscribed collision table
template <typename CollisionObject>
inline bool occupancySelection(CollisionObject const&)
{
  LOGF(fatal, "Occupancy selection not implemented for this kind of collisions");
  return false;
}

/// \brief Occupancy selection for reconstructed and detector level collision tables with centrality/multiplicity information
template <>
inline bool occupancySelection<aod::CollisionEvSelCent>(aod::CollisionEvSelCent const& collision)
{
  return selectOnOccupancy(collision);
}

/// \brief Occupancy selection for reconstructed and detector level collision tables with Run 2 centrality/multiplicity information
template <>
inline bool occupancySelection<aod::CollisionEvSelRun2Cent>(aod::CollisionEvSelRun2Cent const&)
{
  return true;
}

/// \brief Occupancy selection for reconstructed and detector level collision tables without centrality/multiplicity information
template <>
inline bool occupancySelection<aod::CollisionEvSel>(aod::CollisionEvSel const& collision)
{
  return selectOnOccupancy(collision);
}

/// \brief Occupancy selection for detector level collision tables without centrality/multiplicity
template <>
inline bool occupancySelection<soa::Join<aod::CollisionsEvSel, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSel, aod::McCollisionLabels>::iterator const& collision)
{
  return selectOnOccupancy(collision);
}

/// \brief Occupancy selection for detector level collision tables with centrality/multiplicity
template <>
inline bool occupancySelection<soa::Join<aod::CollisionsEvSelCent, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSelCent, aod::McCollisionLabels>::iterator const& collision)
{
  return selectOnOccupancy(collision);
}

/// \brief Occupancy selection for detector level collision tables with Run 2 centrality/multiplicity
template <>
inline bool occupancySelection<soa::Join<aod::CollisionsEvSelRun2Cent, aod::McCollisionLabels>::iterator>(soa::Join<aod::CollisionsEvSelRun2Cent, aod::McCollisionLabels>::iterator const&)
{
  return true;
}

/// \brief Occupancy selection for generator level collision table
template <>
inline bool occupancySelection<aod::McCollision>(aod::McCollision const&)
{
  return true;
}

//////////////////////////////////////////////////////////////////////////////////
/// Event selection
//////////////////////////////////////////////////////////////////////////////////

template <typename CollisionObject>
inline bool isEventSelected(CollisionObject const& collision, float& centormult)
{
  triggerFlags.reset();
  collisionFlags.reset();
  collisionFlags.set(CollSelIN);

  bool trigsel = triggerSelection(collision);

  /* run condition table information */
  bool rctsel = true;
  if (useRctInformation) {
    if constexpr (framework::has_type_v<aod::mccollision::GeneratorsID, typename CollisionObject::all_columns>) {
      /* RCT condition not considered for generator level */
    } else {
      if constexpr (framework::has_type_v<aod::evsel::Rct, typename CollisionObject::all_columns>) {
        if (!rctChecker.checkTable(collision)) {
          rctsel = false;
        }
      } else {
        LOGF(fatal, "RCT check required but the dataset does not have RCT information associated. Please, fix it");
      }
    }
  }
  collisionFlags.set(CollSelRCTBIT, rctsel);

  bool occupancysel = occupancySelection(collision);

  bool zvtxsel = false;
  /* TODO: vertex quality checks */
  if (zvtxlow < collision.posZ() && collision.posZ() < zvtxup) {
    if (onlyInOneSide) {
      if (collision.posZ() != 0.0) {
        /* if only one side, we accept collisions which have zvtx different than zero */
        zvtxsel = true;
        collisionFlags.set(CollSelZVERTEXBIT);
      }
    } else {
      zvtxsel = true;
      collisionFlags.set(CollSelZVERTEXBIT);
    }
  }

  bool centmultsel = centralitySelection(collision, centormult);

  bool centmultexclusion = isCollisionNotExcluded<CollisionObject>();

  bool accepted = trigsel && rctsel && occupancysel && zvtxsel && centmultsel && centmultexclusion;

  if (accepted) {
    collisionFlags.set(CollSelSELECTED);
  }

  return accepted;
}

//////////////////////////////////////////////////////////////////////////////////
/// Track selection
//////////////////////////////////////////////////////////////////////////////////

struct TpcExcludeTrack {
  TpcExcludeTrack()
  {
    method = kNOEXCLUSION;
  }
  explicit TpcExcludeTrack(TpcExclusionMethod m)
  {
    static constexpr float DefaultPhiBinShift = 0.5f;
    static constexpr int DefaultNoOfPhiBins = 72;
    switch (m) {
      case kNOEXCLUSION:
        method = m;
        break;
      case kSTATIC:
        if (phibinshift == DefaultPhiBinShift && phibins == DefaultNoOfPhiBins) {
          method = m;
        } else {
          LOGF(fatal, "Static TPC exclusion method with bin shift: %.2f and number of bins %d. Please fix it", phibinshift, phibins);
        }
        break;
      case kDYNAMIC:
        method = m;
        break;
      default:
        LOGF(fatal, "Wrong TPC tracks exclusion method %d. Please, fix it", static_cast<int>(m));
    }
    philow = 0.0f;
    phiup = constants::math::TwoPI;
    phibinwidth = (phiup - philow) / static_cast<float>(phibins);
    phiup = phiup - phibinwidth * phibinshift;
    philow = philow - phibinwidth * phibinshift;
  }

  template <typename TrackObject>
  int getPhiBinIx(TrackObject const& track)
  {
    float phi = RecoDecay::constrainAngle(track.phi(), philow);
    return static_cast<int>((phi - philow) / phibinwidth);
  }

  template <typename TrackObject>
  bool exclude(TrackObject const& track)
  {
    constexpr int NoOfTpcSectors = 18;
    constexpr float TpcPhiSectorWidth = (constants::math::TwoPI) / NoOfTpcSectors;

    switch (method) {
      case kNOEXCLUSION: {
        return false;
      } break;
      case kSTATIC: {
        int phiBinIx = getPhiBinIx(track);
        /* bins multiple of four have got sector border */
        if ((phiBinIx % 4) != 0) {
          return false;
        } else {
          return true;
        }
      } break;
      case kDYNAMIC: {
        float phiInTpcSector = std::fmod(track.phi(), TpcPhiSectorWidth);
        if (track.sign() > 0) {
          return (phiInTpcSector < positiveUpCut->Eval(track.pt())) && (positiveLowCut->Eval(track.pt()) < phiInTpcSector);
        } else {
          return (phiInTpcSector < negativeUpCut->Eval(track.pt())) && (negativeLowCut->Eval(track.pt()) < phiInTpcSector);
        }
      } break;
      default:
        return false;
    }
  }

  void setCuts(std::string pLowCut, std::string pUpCut, std::string nLowCut, std::string nUpCut)
  {
    LOGF(info, "Setting the TPC exclusion cuts: pLow=%s, pUp=%s, nLow=%s, nUp=%s", pLowCut, pUpCut, nLowCut, nUpCut);
    positiveLowCut = new TF1("posLowCut", pLowCut.c_str(), ptlow, ptup);
    positiveUpCut = new TF1("posUpCut", pUpCut.c_str(), ptlow, ptup);
    negativeLowCut = new TF1("negLowCut", nLowCut.c_str(), ptlow, ptup);
    negativeUpCut = new TF1("negUpCut", nUpCut.c_str(), ptlow, ptup);
  }

  TpcExclusionMethod method = kNOEXCLUSION;
  float phibinwidth = 0.0;
  TF1* positiveLowCut = nullptr;
  TF1* positiveUpCut = nullptr;
  TF1* negativeLowCut = nullptr;
  TF1* negativeUpCut = nullptr;
};

template <typename TrackObject>
inline bool matchTrackType(TrackObject const& track)
{
  using namespace o2::aod::track;

  if (tracktype == TrackTypePWGMM) {
    // under tests MM track selection
    // see: https://indico.cern.ch/event/1383788/contributions/5816953/attachments/2805905/4896281/TrackSel_GlobalTracks_vs_MMTrackSel.pdf
    // it should be equivalent to this
    //       (track.passedDCAxy && track.passedDCAz && track.passedGoldenChi2) &&
    //       (track.passedITSNCls && track.passedITSChi2NDF && track.passedITSHits) &&
    //       (!track.hasTPC || (track.passedTPCNCls && track.passedTPCChi2NDF && track.passedTPCCrossedRowsOverNCls));
    return track.hasITS() && ((track.trackCutFlag() & TrackSelectionITS) == TrackSelectionITS) &&
           (!track.hasTPC() || ((track.trackCutFlag() & TrackSelectionTPC) == TrackSelectionTPC)) &&
           ((track.trackCutFlag() & TrackSelectionDCA) == TrackSelectionDCA);
  } else {
    for (auto const& filter : trackFilters) {
      if (filter->isSelected(track)) {
        return true;
      }
    }
    return false;
  }
}

/// \brief Checks if the passed track is within the acceptance conditions of the analysis
/// \param track the track of interest
/// \return true if the track is in the acceptance, otherwise false
template <typename CollisionsObject, typename TrackObject>
inline bool inTheAcceptance(TrackObject const& track)
{
  /* the side on which the collision happened */
  float side = track.template collision_as<CollisionsObject>().posZ();

  /* overall minimum momentum cut for the analysis */
  if (!(overallminp < track.p())) {
    return false;
  }

  if constexpr (framework::has_type_v<aod::mctracklabel::McParticleId, typename TrackObject::all_columns>) {
    if (track.mcParticleId() < 0) {
      return false;
    }
  }

  /* check the side of the collision and decide if one side */
  if (onlyInOneSide) {
    if (side < 0) {
      if (track.eta() >= 0.0) {
        return false;
      }
    } else if (side > 0) {
      if (track.eta() <= 0.0) {
        return false;
      }
    } else {
      /* if only one side we should not have collisions with zvtx = 0 */
      LOGF(fatal, "Selecting one side tracks with a zvtx zero collision");
    }
  }

  if (ptlow < track.pt() && track.pt() < ptup && etalow < track.eta() && track.eta() < etaup) {
    return !tpcExcluder.exclude(track);
  }
  return false;
}

/// \brief Accepts or not the passed track
/// \param track the track of interest
/// \return true if the track is accepted, otherwise false
template <typename CollisionsObject, typename TrackObject>
inline bool acceptTrack(TrackObject const& track)
{
  if (inTheAcceptance<CollisionsObject>(track)) {
    if (matchTrackType(track)) {
      return true;
    }
  }
  return false;
}

template <typename ParticleObject, typename MCCollisionObject>
void exploreMothers(ParticleObject& particle, MCCollisionObject& collision)
{
  for (const auto& m : particle.template mothers_as<aod::McParticles>()) {
    LOGF(info, "   mother index: %d", m.globalIndex());
    LOGF(info, "   Tracking back mother");
    LOGF(info, "   assigned collision Id: %d, looping on collision Id: %d", m.mcCollisionId(), collision.globalIndex());
    LOGF(info, "   index: %d, pdg code: %d", m.globalIndex(), m.pdgCode());
    LOGF(info, "   Passed  isPhysicalPrimary(): %s", m.isPhysicalPrimary() ? "YES" : "NO");

    exploreMothers(m, collision);
  }
}

inline float getCharge(float pdgCharge)
{
  static constexpr int NoOfBasicChargesPerUnitCharge = 3;
  float charge = (pdgCharge / NoOfBasicChargesPerUnitCharge >= 1) ? 1.0 : ((pdgCharge / NoOfBasicChargesPerUnitCharge <= -1) ? -1.0 : 0);
  return charge;
}

/// \brief Accepts or not the passed generated particle
/// \param track the particle of interest
/// \return `true` if the particle is accepted, `false` otherwise
template <typename ParticleObject, typename MCCollisionObject>
inline bool acceptParticle(ParticleObject& particle, MCCollisionObject const&)
{
  /* overall momentum cut */
  if (!(overallminp < particle.p())) {
    return false;
  }

  if (particle.isPhysicalPrimary()) {
    if ((particle.mcCollisionId() == 0) && traceCollId0) {
      LOGF(info, "Particle %d passed isPhysicalPrimary", particle.globalIndex());
    }

    if (ptlow < particle.pt() && particle.pt() < ptup && etalow < particle.eta() && particle.eta() < etaup) {
      return true;
    }
  } else {
    if ((particle.mcCollisionId() == 0) && traceCollId0) {
      LOGF(info, "Particle %d NOT passed isPhysicalPrimary", particle.globalIndex());
    }
  }
  return false;
}

//////////////////////////////////////////////////////////////////////////////////
/// PID
//////////////////////////////////////////////////////////////////////////////////

struct PIDSpeciesSelection {
  const std::vector<int> pdgcodes = {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton};
  const std::vector<std::string_view> spnames = {"e", "mu", "pi", "ka", "p"};
  const std::vector<std::string_view> sptitles = {"e", "#mu", "#pi", "K", "p"};
  const std::vector<std::string_view> spfnames = {"E", "Mu", "Pi", "Ka", "Pr"};
  const std::vector<std::string_view> spadjnames = {"Electron", "Muon", "Pion", "Kaon", "Proton"};
  const std::vector<double> spmasses = {o2::constants::physics::MassElectron, o2::constants::physics::MassMuon, o2::constants::physics::MassPionCharged, o2::constants::physics::MassKaonCharged, o2::constants::physics::MassProton};
  const std::vector<std::string_view> chadjnames = {"P", "M"};
  const char* hadname = "h";
  const char* hadtitle = "h";
  const char* hadfname = "Ha";
  uint getNSpecies() { return config.size(); }
  const char* getSpeciesName(uint8_t ix) { return spnames[species[ix]].data(); }
  const char* getSpeciesTitle(uint8_t ix) { return sptitles[species[ix]].data(); }
  const char* getSpeciesFName(uint8_t ix) { return spfnames[species[ix]].data(); }
  double getSpeciesMass(uint8_t ix) { return spmasses[species[ix]]; }
  const char* getHadName() { return hadname; }
  const char* getHadTitle() { return hadtitle; }
  const char* getHadFName() { return hadfname; }
  bool isSpeciesBeingSelected(uint8_t sp) { return std::find(species.begin(), species.end(), sp) != species.end(); }
  bool isGlobalSpecies(uint8_t isp, o2::track::PID::ID glsp) { return species[isp] == glsp; }
  void storePIDAdjustments(TList* lst)
  {
    auto storedetectorwithcharge = [&](auto& detectorstore, auto detectorname, auto charge) {
      for (uint isp = 0; isp < spadjnames.size(); ++isp) {
        TString fullhname = TString::Format("%s%s%s_Difference", detectorname, spadjnames[isp].data(), charge);
        detectorstore[isp] = static_cast<TH1*>(lst->FindObject(fullhname.Data()));
      }
    };
    auto reportadjdetectorwithcharge = [&](auto& detectorstore, auto detectorname, auto charge) {
      for (uint isp = 0; isp < spadjnames.size(); ++isp) {
        if (detectorstore[isp] != nullptr) {
          LOGF(info, "Stored nsigmas adjust for detector %s and species %s%s in histogram %s", detectorname, spadjnames[isp].data(), charge, detectorstore[isp]->GetName());
        }
      }
    };
    storedetectorwithcharge(tpcnsigmasshiftpos, "TPC", "P");
    storedetectorwithcharge(tofnsigmasshiftpos, "TOF", "P");
    storedetectorwithcharge(tpcnsigmasshiftneg, "TPC", "M");
    storedetectorwithcharge(tofnsigmasshiftneg, "TOF", "M");

    reportadjdetectorwithcharge(tpcnsigmasshiftpos, "TPC", "P");
    reportadjdetectorwithcharge(tofnsigmasshiftpos, "TOF", "P");
    reportadjdetectorwithcharge(tpcnsigmasshiftneg, "TPC", "M");
    reportadjdetectorwithcharge(tofnsigmasshiftneg, "TOF", "M");
  }
  void addSpecies(uint8_t sp, o2::analysis::TrackSelectionPIDCfg* incfg)
  {
    o2::analysis::TrackSelectionPIDCfg* cfg = new o2::analysis::TrackSelectionPIDCfg(*incfg);
    config.push_back(cfg);
    species.push_back(sp);

    auto last = config[config.size() - 1];
    uint8_t lastsp = species[config.size() - 1];
    LOGF(info, "Inserted species %d with", lastsp);
    LOGF(info, "  minTPC nsigmas: el: %.2f, mu: %.2f, pi: %.2f, ka: %.2f, pr: %.2f", last->mMinNSigmasTPC[0], last->mMinNSigmasTPC[1], last->mMinNSigmasTPC[2], last->mMinNSigmasTPC[3], last->mMinNSigmasTPC[4]);
    LOGF(info, "  maxTPC nsigmas: el: %.2f, mu: %.2f, pi: %.2f, ka: %.2f, pr: %.2f", last->mMaxNSigmasTPC[0], last->mMaxNSigmasTPC[1], last->mMaxNSigmasTPC[2], last->mMaxNSigmasTPC[3], last->mMaxNSigmasTPC[4]);
    LOGF(info, "  minTOF nsigmas: el: %.2f, mu: %.2f, pi: %.2f, ka: %.2f, pr: %.2f", last->mMinNSigmasTOF[0], last->mMinNSigmasTOF[1], last->mMinNSigmasTOF[2], last->mMinNSigmasTOF[3], last->mMinNSigmasTOF[4]);
    LOGF(info, "  maxTOF nsigmas: el: %.2f, mu: %.2f, pi: %.2f, ka: %.2f, pr: %.2f", last->mMaxNSigmasTOF[0], last->mMaxNSigmasTOF[1], last->mMaxNSigmasTOF[2], last->mMaxNSigmasTOF[3], last->mMaxNSigmasTOF[4]);
    LOGF(info, "  %.1f < pT < %.1f", last->mPtMin, last->mPtMax);
  }
  void addExcludedSpecies(uint8_t sp, const o2::analysis::TrackSelectionPIDCfg* incfg)
  {
    o2::analysis::TrackSelectionPIDCfg* cfg = new o2::analysis::TrackSelectionPIDCfg(*incfg);
    configexclude.push_back(cfg);
    speciesexclude.push_back(sp);
    auto last = configexclude[configexclude.size() - 1];
    uint8_t lastsp = speciesexclude[configexclude.size() - 1];

    LOGF(info, "Inserted species %d for exclusion with", lastsp);
    LOGF(info, "  minTPC nsigmas: el: %.2f, mu: %.2f, pi: %.2f, ka: %.2f, pr: %.2f", last->mMinNSigmasTPC[0], last->mMinNSigmasTPC[1], last->mMinNSigmasTPC[2], last->mMinNSigmasTPC[3], last->mMinNSigmasTPC[4]);
    LOGF(info, "  maxTPC nsigmas: el: %.2f, mu: %.2f, pi: %.2f, ka: %.2f, pr: %.2f", last->mMaxNSigmasTPC[0], last->mMaxNSigmasTPC[1], last->mMaxNSigmasTPC[2], last->mMaxNSigmasTPC[3], last->mMaxNSigmasTPC[4]);
    LOGF(info, "  minTOF nsigmas: el: %.2f, mu: %.2f, pi: %.2f, ka: %.2f, pr: %.2f", last->mMinNSigmasTOF[0], last->mMinNSigmasTOF[1], last->mMinNSigmasTOF[2], last->mMinNSigmasTOF[3], last->mMinNSigmasTOF[4]);
    LOGF(info, "  maxTOF nsigmas: el: %.2f, mu: %.2f, pi: %.2f, ka: %.2f, pr: %.2f", last->mMaxNSigmasTOF[0], last->mMaxNSigmasTOF[1], last->mMaxNSigmasTOF[2], last->mMaxNSigmasTOF[3], last->mMaxNSigmasTOF[4]);
    LOGF(info, "  %.1f < pT < %.1f", last->mPtMin, last->mPtMax);
  }
  template <StrongDebugging outdebug, typename TrackObject>
  int8_t whichSpecies(TrackObject const& track)
  {
    TString debuginfo;
    std::vector<float> tpcnsigmas = {track.tpcNSigmaEl(), track.tpcNSigmaMu(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
    std::vector<float> tofnsigmas = {track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr()};

    auto outmomentumdebug = [&]() {
      if constexpr (outdebug != 0) {
        debuginfo += TString::Format("%.5f,%.5f,%.5f,%d,%.2f,%.4f,", track.p(), track.tpcInnerParam(), track.pt(), track.hasTOF() ? 1 : 0, track.tpcSignal(), track.beta());
      }
    };
    auto outnsigmasdebug = [&]() {
      if constexpr (outdebug != 0) {
        for (auto const& tpcn : tpcnsigmas) {
          debuginfo += TString::Format("%.4f,", tpcn);
        }
        for (auto const& tofn : tofnsigmas) {
          debuginfo += TString::Format("%.4f,", tofn);
        }
      }
    };

    /* out debug if needed */
    outmomentumdebug();
    /* out debug if needed */
    outnsigmasdebug();

    auto closeTo = [](auto& values, auto& mindet, auto& maxdet, uint8_t sp) {
      if (mindet[sp] <= values[sp] && values[sp] < maxdet[sp]) {
        return true;
      } else {
        return false;
      }
    };
    auto awayFrom = [&](auto& values, auto& mindet, auto& maxdet, uint8_t sp) {
      for (size_t ix = 0; ix < pdgcodes.size(); ix++) {
        if (ix != sp) {
          if (mindet[ix] <= values[ix] && values[ix] < maxdet[ix]) {
            return false;
          }
        } else {
          continue;
        }
      }
      return true;
    };
    auto closeToTPC = [&](auto& config, uint8_t sp) {
      return closeTo(tpcnsigmas, config->mMinNSigmasTPC, config->mMaxNSigmasTPC, sp);
    };
    auto awayFromTPC = [&](auto& config, uint8_t sp) {
      return awayFrom(tpcnsigmas, config->mMinNSigmasTPC, config->mMaxNSigmasTPC, sp);
    };
    auto closeToTOF = [&](auto& config, uint8_t sp) {
      return closeTo(tofnsigmas, config->mMinNSigmasTOF, config->mMaxNSigmasTOF, sp);
    };
    auto closeToTPCTOF = [&](auto& config, uint8_t sp) {
      if (config->m2Dcut) {
        float a = (config->mMaxNSigmasTPC[sp] - config->mMinNSigmasTPC[sp]) / 2.0;
        float b = (config->mMaxNSigmasTOF[sp] - config->mMinNSigmasTOF[sp]) / 2.0;
        float oa = (config->mMaxNSigmasTPC[sp] + config->mMinNSigmasTPC[sp]) / 2.0;
        float ob = (config->mMaxNSigmasTOF[sp] + config->mMinNSigmasTOF[sp]) / 2.0;
        float vtpc = tpcnsigmas[sp] - oa;
        float vtof = tofnsigmas[sp] - ob;
        return ((vtpc * vtpc / a / a + vtof * vtof / b / b) < 1);
      } else {
        return closeToTPC(config, sp) && closeToTOF(config, sp);
      }
    };
    auto awayFromTPCTOF = [&](auto& config, uint8_t sp) {
      for (uint8_t ix = 0; ix < pdgcodes.size(); ++ix) {
        if (ix != sp) {
          if (closeToTPCTOF(config, ix)) {
            return false;
          }
        } else {
          continue;
        }
      }
      return true;
    };
    auto aboveThreshold = [&](auto& config) {
      return ((config->mPThreshold > 0.0) && (config->mPThreshold < track.p()));
    };
    auto isA = [&](auto& config, uint8_t sp) {
      if (track.hasTOF()) {
        return closeToTPCTOF(config, sp) && awayFromTPCTOF(config, sp);
      } else {
        if (aboveThreshold(config)) {
          if (config->mRequireTOF) {
            return false;
          }
        }
      }
      /* we are here                                                */
      /* - below the threshold                                      */
      /* - above the threshold without TOF and without requiring it */
      /* so we check only the TPC information                       */
      return closeToTPC(config, sp) && awayFromTPC(config, sp);
    };

    auto adjustnsigmas = [&]() {
      if (track.sign() > 0) {
        for (uint isp = 0; isp < spnames.size(); ++isp) {
          if (tpcnsigmasshiftpos[isp] != nullptr) {
            TH1* h = tpcnsigmasshiftpos[isp];
            tpcnsigmas[isp] -= h->GetBinContent(h->GetXaxis()->FindFixBin(track.p()));
          }
          if (tofnsigmasshiftpos[isp] != nullptr) {
            TH1* h = tofnsigmasshiftpos[isp];
            tofnsigmas[isp] -= h->GetBinContent(h->GetXaxis()->FindFixBin(track.p()));
          }
        }
      } else {
        for (uint isp = 0; isp < spnames.size(); ++isp) {
          if (tpcnsigmasshiftneg[isp] != nullptr) {
            TH1* h = tpcnsigmasshiftneg[isp];
            tpcnsigmas[isp] -= h->GetBinContent(h->GetXaxis()->FindFixBin(track.p()));
          }
          if (tofnsigmasshiftneg[isp] != nullptr) {
            TH1* h = tofnsigmasshiftneg[isp];
            tofnsigmas[isp] -= h->GetBinContent(h->GetXaxis()->FindFixBin(track.p()));
          }
        }
      }
    };

    auto outpiddebug = [&](int code, int pid, int pid2) {
      if constexpr (outdebug != 0) {
        int truepid = -1;
        int isphysicalprimary = -1;
        int process = -1;
        if constexpr (framework::has_type_v<aod::mctracklabel::McParticleId, typename TrackObject::all_columns>) {
          if (!(track.mcParticleId() < 0)) {
            auto particle = track.template mcParticle_as<aod::McParticles>();
            truepid = particle.pdgCode();
            isphysicalprimary = particle.isPhysicalPrimary() ? 1 : 0;
            process = particle.getProcess();
          }
        }
        debuginfo += TString::Format("%d,%d,%d,%d,%d,%d\n", code, pid, pid2, truepid, isphysicalprimary, process);
        debugstream << debuginfo;
      }
    };

    /* adjust the nsigmas values if appropriate */
    adjustnsigmas();
    /* out debug info if needed */
    outnsigmasdebug();

    /* let's first check the exclusion from the analysis */
    for (uint8_t ix = 0; ix < configexclude.size(); ++ix) {
      if (isA(configexclude[ix], speciesexclude[ix])) {
        /* out debug info if needed */
        outpiddebug(1, speciesexclude[ix], -1);
        return -(ix + 1);
      }
    }
    /* we don't exclude it so check which species if any required */
    if (config.size() > 0) {
      int8_t id = -127;
      for (uint8_t ix = 0; ix < config.size(); ++ix) {
        if (isA(config[ix], species[ix])) {
          if (id < 0) {
            id = ix;
          } else {
            /* out debug info if needed */
            outpiddebug(2, species[id], species[ix]);
            /* already identified once */
            return -127;
          }
        }
      }
      /* check a species transverse momentum cut */
      if (!(id < 0)) {
        /* identified */
        if ((track.pt() < config[id]->mPtMin) || config[id]->mPtMax < track.pt()) {
          /* rejected */
          outpiddebug(4, species[id], -1);
          return -127;
        }
      }
      /* out debug info if needed */
      if (id < 0) {
        /* not identified */
        outpiddebug(3, -1, -1);
      } else {
        /* identified */
        outpiddebug(0, species[id], -1);
      }
      return id;
    } else {
      outpiddebug(0, 0, -1);
      /* charged hadron */
      return 0;
    }
  }
  template <typename ParticleObject>
  int8_t whichTruthSpecies(ParticleObject part)
  {
    int pdgcode = std::abs(part.pdgCode());
    /* let's first check the exclusion from the analysis */
    for (uint8_t ix = 0; ix < configexclude.size(); ++ix) {
      if (pdgcode == pdgcodes[speciesexclude[ix]]) {
        return -(ix + 1);
      }
    }
    /* we don't exclude it so check which species if any required */
    if (config.size() > 0) {
      for (uint8_t ix = 0; ix < config.size(); ++ix) {
        if (pdgcode == pdgcodes[species[ix]]) {
          /* check a species transverse momentum cut */
          if ((part.pt() < config[ix]->mPtMin) || config[ix]->mPtMax < part.pt()) {
            /* rejected */
            return -127;
          }
          return ix;
        }
      }
      return -127;
    } else {
      return 0;
    }
  }

  template <typename ParticleObject>
  int8_t whichTruthPrimarySpecies(ParticleObject part)
  {
    if (part.isPhysicalPrimary()) {
      return whichTruthSpecies(part);
    } else {
      return -127;
    }
  }

  template <typename ParticleObject>
  int8_t whichTruthSecondarySpecies(ParticleObject part)
  {
    if (!(part.isPhysicalPrimary())) {
      return whichTruthSpecies(part);
    } else {
      return -127;
    }
  }

  template <typename ParticleObject>
  int8_t whichTruthSecFromDecaySpecies(ParticleObject part)
  {
    if (!(part.isPhysicalPrimary()) && (part.getProcess() == TMCProcess::kPDecay)) {
      return whichTruthSpecies(part);
    } else {
      return -127;
    }
  }

  template <typename ParticleObject>
  int8_t whichTruthSecFromMaterialSpecies(ParticleObject part)
  {
    if (!(part.isPhysicalPrimary()) && (part.getProcess() != TMCProcess::kPDecay)) {
      return whichTruthSpecies(part);
    } else {
      return -127;
    }
  }

  std::vector<const o2::analysis::TrackSelectionPIDCfg*> config;        ///< the PID selection configuration of the species to include in the analysis
  std::vector<uint8_t> species;                                         ///< the species index of the species to include in the analysis
  std::vector<const o2::analysis::TrackSelectionPIDCfg*> configexclude; ///< the PID selection configuration of the species to exclude from the analysis
  std::vector<uint8_t> speciesexclude;                                  ///< the species index of teh species to exclude from the analysis
  std::vector<TH1*> tpcnsigmasshiftpos{spnames.size(), nullptr};
  std::vector<TH1*> tpcnsigmasshiftneg{spnames.size(), nullptr};
  std::vector<TH1*> tofnsigmasshiftpos{spnames.size(), nullptr};
  std::vector<TH1*> tofnsigmasshiftneg{spnames.size(), nullptr};
};

} // namespace dptdptfilter
} // namespace analysis
} // namespace o2

#endif // PWGCF_TABLEPRODUCER_DPTDPTFILTER_H_
