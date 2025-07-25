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

#include "ProfileSubset.h"
#include "TProfile2D.h"

TProfile* ProfileSubset::GetSubset(bool onX, const char* name, int firstbin, int lastbin, int l_nbins, double* l_binarray)
{
  TString expectedName = (onX ? "_pfx" : "_pfy");
  TString pname(name);
  if (pname.IsNull() || name == expectedName)
    pname = TString(GetName()) + expectedName;
  const TAxis& outAxis = (onX ? fXaxis : fYaxis);
  const TArrayD* bins = outAxis.GetXbins();
  int firstOutBin = outAxis.GetFirst();
  int lastOutBin = outAxis.GetLast();
  TProfile* p1 = 0;
  if (l_nbins)
    p1 = new TProfile(pname, GetTitle(), l_nbins, l_binarray);
  else
    p1 = new TProfile(pname, GetTitle(), outAxis.GetNbins(), bins->fArray);
  if (fBinSumw2.fN)
    p1->Sumw2();
  TH2D* h2dW = ProjectionXY("h2temp-W", "W");
  TH2D* h2dN = ProjectionXY("h2temp-N", "B");
  h2dW->SetDirectory(0);
  h2dN->SetDirectory(0);
  if (onX) {
    h2dW->GetXaxis()->SetRange(firstOutBin, lastOutBin);
    h2dN->GetXaxis()->SetRange(firstOutBin, lastOutBin);
  } else {
    h2dW->GetYaxis()->SetRange(firstOutBin, lastOutBin);
    h2dN->GetYaxis()->SetRange(firstOutBin, lastOutBin);
  }
  TH1D* h1W = (onX) ? h2dW->ProjectionX("h1temp-W", firstbin, lastbin) : h2dW->ProjectionY("h1temp-W", firstbin, lastbin);
  TH1D* h1N = (onX) ? h2dN->ProjectionX("h1temp-N", firstbin, lastbin) : h2dN->ProjectionY("h1temp-N", firstbin, lastbin);
  h1W->SetDirectory(0);
  h1N->SetDirectory(0);
  R__ASSERT(h1W->fN == p1->fN);
  R__ASSERT(h1N->fN == p1->fN);
  R__ASSERT(h1W->GetSumw2()->fN != 0);
  for (int i = 0; i < p1->fN; ++i) {
    p1->fArray[i] = h1W->GetBinContent(i);                  // array of profile is sum of all values
    p1->GetSumw2()->fArray[i] = h1W->GetSumw2()->fArray[i]; // array of content square of profile is weight square of the W projected histogram
    p1->SetBinEntries(i, h1N->GetBinContent(i));
    if (fBinSumw2.fN)
      p1->GetBinSumw2()->fArray[i] = h1N->GetSumw2()->fArray[i]; // sum of weight squares are stored to compute errors in h1N histogram
  }
  delete h2dW;
  delete h2dN;
  delete h1W;
  delete h1N;
  p1->SetEntries(p1->GetEffectiveEntries());
  return p1;
};
void ProfileSubset::OverrideBinContent(double x, double y, double x2, double y2, double /*val*/)
{
  if (!fBinSumw2.fN)
    Sumw2();
  TH2D* h2dW = ProjectionXY("h2temp-W", "W");
  TH2D* h2dN = ProjectionXY("h2temp-N", "B");
  int binIndex = FindBin(x, y);
  int binIndex2 = FindBin(x2, y2);
  fArray[binIndex] = h2dW->GetBinContent(binIndex2);
  GetSumw2()->fArray[binIndex] = h2dW->GetSumw2()->fArray[binIndex2];
  SetBinEntries(binIndex, h2dN->GetBinContent(binIndex2));
  if (fBinSumw2.fN)
    GetBinSumw2()->fArray[binIndex] = h2dN->GetSumw2()->fArray[binIndex2];
}
void ProfileSubset::OverrideBinContent(double x, double y, double x2, double y2, TProfile2D* sourceProf)
{
  if (!fBinSumw2.fN)
    Sumw2();
  if (!sourceProf->fN)
    sourceProf->Sumw2();
  TH2D* h2dW = sourceProf->ProjectionXY("h2temp-W", "W");
  TH2D* h2dN = sourceProf->ProjectionXY("h2temp-N", "B");
  int binIndex = FindBin(x, y);
  int binIndex2 = sourceProf->FindBin(x2, y2);
  fArray[binIndex] = h2dW->GetBinContent(binIndex2);
  GetSumw2()->fArray[binIndex] = h2dW->GetSumw2()->fArray[binIndex2];
  SetBinEntries(binIndex, h2dN->GetBinContent(binIndex2));
  if (fBinSumw2.fN)
    GetBinSumw2()->fArray[binIndex] = h2dN->GetSumw2()->fArray[binIndex2];
}
bool ProfileSubset::OverrideBinsWithZero(int xb1, int yb1, int xb2, int yb2)
{
  bool lHaveToQuit = kFALSE;
  if (GetNbinsX() < xb1 || GetNbinsX() < xb2) {
    lHaveToQuit = kTRUE;
    printf("xBins out of range! (%i-%i vs %i)\n", xb1, xb2, GetNbinsX());
  }
  if (GetNbinsY() < yb1 || GetNbinsY() < yb2) {
    lHaveToQuit = kTRUE;
    printf("yBins out of range! (%i-%i vs %i)\n", yb1, yb2, GetNbinsY());
  }
  if (lHaveToQuit)
    return kFALSE;
  for (int ix = xb1; ix <= xb2; ix++) {
    for (int iy = yb1; iy <= yb2; iy++) {
      int bind = FindBin(GetXaxis()->GetBinCenter(ix), GetYaxis()->GetBinCenter(iy));
      fArray[bind] = 0;
      GetSumw2()->fArray[bind] = 0;
      SetBinEntries(bind, 0);
      if (fBinSumw2.fN)
        GetBinSumw2()->fArray[bind] = 0;
    }
  }
  return kTRUE;
}
