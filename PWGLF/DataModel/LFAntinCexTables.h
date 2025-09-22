#ifndef PWGLF_DATAMODEL_LFANTINCEXTABLES_H_
#define PWGLF_DATAMODEL_LFANTINCEXTABLES_H_

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"

namespace o2::aod
{
namespace antinCexNS
{
// Etiqueta/metadata básica
DECLARE_SOA_COLUMN(IsCex,           isCex,           bool);     // 1=CEX (from antin), 0=BG
DECLARE_SOA_COLUMN(MotherPdg,       motherPdg,       int32_t);  // mother PDG
DECLARE_SOA_COLUMN(ColId,           colId,           int32_t);  // mcCollisionId
DECLARE_SOA_COLUMN(PId,             pId,             int32_t);  // proton MC id
DECLARE_SOA_COLUMN(AntipId,         antipId,         int32_t);  // antiproton MC id

// MC (par)
DECLARE_SOA_COLUMN(McPairP,         mcPairP,         float);
DECLARE_SOA_COLUMN(McPairPt,        mcPairPt,        float);
DECLARE_SOA_COLUMN(McPairPz,        mcPairPz,        float);
DECLARE_SOA_COLUMN(McDplane,        mcDplane,        float);
DECLARE_SOA_COLUMN(McAngleDeg,      mcAngleDeg,      float);
DECLARE_SOA_COLUMN(McVtxX,          mcVtxX,          float);
DECLARE_SOA_COLUMN(McVtxY,          mcVtxY,          float);
DECLARE_SOA_COLUMN(McVtxZ,          mcVtxZ,          float);

// Tracks (par, del fitter)
DECLARE_SOA_COLUMN(TrkPairP,        trkPairP,        float);
DECLARE_SOA_COLUMN(TrkPairPt,       trkPairPt,       float);
DECLARE_SOA_COLUMN(TrkPairPz,       trkPairPz,       float);
DECLARE_SOA_COLUMN(TrkAngleDeg,     trkAngleDeg,     float);
DECLARE_SOA_COLUMN(TrkVtxfitDcaPair,trkVtxfitDcaPair,float);
DECLARE_SOA_COLUMN(TrkVtxfitR,      trkVtxfitR,      float);
DECLARE_SOA_COLUMN(TrkVtxfitDistToPv,trkVtxfitDistToPv,float);
DECLARE_SOA_COLUMN(TrkVtxfitSecVtxX,trkVtxfitSecVtxX,float);
DECLARE_SOA_COLUMN(TrkVtxfitSecVtxY,trkVtxfitSecVtxY,float);
DECLARE_SOA_COLUMN(TrkVtxfitSecVtxZ,trkVtxfitSecVtxZ,float);

// Calidad del fit y residuales (fit − MC)
DECLARE_SOA_COLUMN(VtxfitChi2,      vtxfitChi2,      float);
DECLARE_SOA_COLUMN(VtxfitStatus,    vtxfitStatus,    int32_t);
DECLARE_SOA_COLUMN(NCand,           nCand,           int32_t);
DECLARE_SOA_COLUMN(VtxfitDX,        vtxfitDX,        float);
DECLARE_SOA_COLUMN(VtxfitDY,        vtxfitDY,        float);
DECLARE_SOA_COLUMN(VtxfitDZ,        vtxfitDZ,        float);
DECLARE_SOA_COLUMN(VtxfitD3D,       vtxfitD3D,       float);

// Track del protón
DECLARE_SOA_COLUMN(PTrkP,           pTrkP,           float);
DECLARE_SOA_COLUMN(PTrkPx,          pTrkPx,          float);
DECLARE_SOA_COLUMN(PTrkPy,          pTrkPy,          float);
DECLARE_SOA_COLUMN(PTrkPz,          pTrkPz,          float);
DECLARE_SOA_COLUMN(PTrkEta,         pTrkEta,         float);
DECLARE_SOA_COLUMN(PTrkTpcSignal,   pTrkTpcSignal,   float);
DECLARE_SOA_COLUMN(PTrkNClsIts,     pTrkNClsIts,     int16_t);

// Track del antiprotón
DECLARE_SOA_COLUMN(AntipTrkP,         antipTrkP,         float);
DECLARE_SOA_COLUMN(AntipTrkPx,        antipTrkPx,        float);
DECLARE_SOA_COLUMN(AntipTrkPy,        antipTrkPy,        float);
DECLARE_SOA_COLUMN(AntipTrkPz,        antipTrkPz,        float);
DECLARE_SOA_COLUMN(AntipTrkEta,       antipTrkEta,       float);
DECLARE_SOA_COLUMN(AntipTrkTpcSignal, antipTrkTpcSignal, float);
DECLARE_SOA_COLUMN(AntipTrkNClsIts,   antipTrkNClsIts,   int16_t);

// (Opcional pero útil) Máscara de selección por cortes
DECLARE_SOA_COLUMN(SelMask,         selMask,         uint32_t); // bits de cortes superados
} // namespace antinCexNS

// Tabla única para CEX y BG (flag isCex)
DECLARE_SOA_TABLE(antinCexPairs, "AOD", "ANTINCEX",
  antinCexNS::IsCex,
  antinCexNS::MotherPdg, antinCexNS::ColId, antinCexNS::PId, antinCexNS::AntipId,
  antinCexNS::McPairP, antinCexNS::McPairPt, antinCexNS::McPairPz,
  antinCexNS::McDplane, antinCexNS::McAngleDeg, antinCexNS::McVtxX, antinCexNS::McVtxY, antinCexNS::McVtxZ,
  antinCexNS::TrkPairP, antinCexNS::TrkPairPt, antinCexNS::TrkPairPz, antinCexNS::TrkAngleDeg,
  antinCexNS::TrkVtxfitDcaPair, antinCexNS::TrkVtxfitR, antinCexNS::TrkVtxfitDistToPv,
  antinCexNS::TrkVtxfitSecVtxX, antinCexNS::TrkVtxfitSecVtxY, antinCexNS::TrkVtxfitSecVtxZ,
  antinCexNS::VtxfitChi2, antinCexNS::VtxfitStatus, antinCexNS::NCand,
  antinCexNS::VtxfitDX, antinCexNS::VtxfitDY, antinCexNS::VtxfitDZ, antinCexNS::VtxfitD3D,
  antinCexNS::PTrkP, antinCexNS::PTrkPx, antinCexNS::PTrkPy, antinCexNS::PTrkPz, antinCexNS::PTrkEta, antinCexNS::PTrkTpcSignal, antinCexNS::PTrkNClsIts,
  antinCexNS::AntipTrkP, antinCexNS::AntipTrkPx, antinCexNS::AntipTrkPy, antinCexNS::AntipTrkPz, antinCexNS::AntipTrkEta, antinCexNS::AntipTrkTpcSignal, antinCexNS::AntipTrkNClsIts,
  antinCexNS::SelMask
)

} // namespace o2::aod

#endif // PWGLF_DATAMODEL_LFANTINCEXTABLES_H_

