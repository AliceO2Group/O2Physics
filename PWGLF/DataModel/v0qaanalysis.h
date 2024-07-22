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
/// \brief QA task for V0 analysis using derived data
///
/// \author Francesca Ercolessi (francesca.ercolessi@cern.ch)

#ifndef PWGLF_DATAMODEL_V0QAANALYSIS_H_
#define PWGLF_DATAMODEL_V0QAANALYSIS_H_

namespace o2::aod
{

namespace myv0candidates
{

DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(V0Pt, v0pt, float);
DECLARE_SOA_COLUMN(RapLambda, raplambda, float);
DECLARE_SOA_COLUMN(RapK0Short, rapk0short, float);
DECLARE_SOA_COLUMN(MassLambda, masslambda, float);
DECLARE_SOA_COLUMN(MassAntiLambda, massantilambda, float);
DECLARE_SOA_COLUMN(MassK0Short, massk0short, float);
DECLARE_SOA_COLUMN(V0Radius, v0radius, float);
DECLARE_SOA_COLUMN(V0CosPA, v0cospa, float);
DECLARE_SOA_COLUMN(V0DCAPosToPV, v0dcapostopv, float);
DECLARE_SOA_COLUMN(V0DCANegToPV, v0dcanegtopv, float);
DECLARE_SOA_COLUMN(V0DCAV0Daughters, v0dcav0daughters, float);
DECLARE_SOA_COLUMN(V0PosEta, v0poseta, float);
DECLARE_SOA_COLUMN(V0NegEta, v0negeta, float);
DECLARE_SOA_COLUMN(V0PosPhi, v0posphi, float);
DECLARE_SOA_COLUMN(V0NegPhi, v0negphi, float);
DECLARE_SOA_COLUMN(V0PosITSHits, v0positshits, float);
DECLARE_SOA_COLUMN(V0NegITSHits, v0negitshits, float);
DECLARE_SOA_COLUMN(CtauLambda, ctaulambda, float);
DECLARE_SOA_COLUMN(CtauAntiLambda, ctauantilambda, float);
DECLARE_SOA_COLUMN(CtauK0Short, ctauk0short, float);
DECLARE_SOA_COLUMN(NTPCSigmaNegPr, ntpcsigmanegpr, float);
DECLARE_SOA_COLUMN(NTPCSigmaPosPr, ntpcsigmapospr, float);
DECLARE_SOA_COLUMN(NTPCSigmaNegPi, ntpcsigmanegpi, float);
DECLARE_SOA_COLUMN(NTPCSigmaPosPi, ntpcsigmapospi, float);
DECLARE_SOA_COLUMN(NTOFSigmaNegPr, ntofsigmanegpr, float);
DECLARE_SOA_COLUMN(NTOFSigmaPosPr, ntofsigmapospr, float);
DECLARE_SOA_COLUMN(NTOFSigmaNegPi, ntofsigmanegpi, float);
DECLARE_SOA_COLUMN(NTOFSigmaPosPi, ntofsigmapospi, float);
DECLARE_SOA_COLUMN(PosHasTOF, poshastof, float);
DECLARE_SOA_COLUMN(NegHasTOF, neghastof, float);
DECLARE_SOA_COLUMN(PDGCode, pdgcode, int);
DECLARE_SOA_COLUMN(IsPhysicalPrimary, isphysprimary, bool);
DECLARE_SOA_COLUMN(MultFT0M, multft0m, float);
DECLARE_SOA_COLUMN(MultFV0A, multfv0a, float);
DECLARE_SOA_COLUMN(EvFlag, evflag, int);

} // namespace myv0candidates

DECLARE_SOA_TABLE(MyV0Candidates, "AOD", "MYV0CANDIDATES", o2::soa::Index<>,
                  myv0candidates::CollisionId, myv0candidates::V0Pt, myv0candidates::RapLambda, myv0candidates::RapK0Short,
                  myv0candidates::MassLambda, myv0candidates::MassAntiLambda, myv0candidates::MassK0Short,
                  myv0candidates::V0Radius, myv0candidates::V0CosPA, myv0candidates::V0DCAPosToPV,
                  myv0candidates::V0DCANegToPV, myv0candidates::V0DCAV0Daughters,
                  myv0candidates::V0PosEta, myv0candidates::V0NegEta, myv0candidates::V0PosPhi, myv0candidates::V0NegPhi,
                  myv0candidates::V0PosITSHits, myv0candidates::V0NegITSHits, myv0candidates::CtauLambda, myv0candidates::CtauAntiLambda, myv0candidates::CtauK0Short,
                  myv0candidates::NTPCSigmaNegPr, myv0candidates::NTPCSigmaPosPr, myv0candidates::NTPCSigmaNegPi, myv0candidates::NTPCSigmaPosPi,
                  myv0candidates::NTOFSigmaNegPr, myv0candidates::NTOFSigmaPosPr, myv0candidates::NTOFSigmaNegPi, myv0candidates::NTOFSigmaPosPi,
                  myv0candidates::PosHasTOF, myv0candidates::NegHasTOF,
                  myv0candidates::PDGCode, myv0candidates::IsPhysicalPrimary,
                  myv0candidates::MultFT0M, myv0candidates::MultFV0A,
                  myv0candidates::EvFlag);

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_V0QAANALYSIS_H_
