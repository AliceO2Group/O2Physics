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

/// \author Luca Micheletti <luca.micheletti@cern.ch>, CERN

#include "CB2Pdf.h"

#include "Riostream.h"

#include "RooAbsCategory.h"
#include "RooAbsReal.h"
#include "TMath.h"

#include <math.h>

ClassImp(CB2Pdf);

CB2Pdf::CB2Pdf(const char* name, const char* title,
               RooAbsReal& _x,
               RooAbsReal& _A,
               RooAbsReal& _B,
               RooAbsReal& _C,
               RooAbsReal& _D,
               RooAbsReal& _E,
               RooAbsReal& _F) : RooAbsPdf(name, title),
                                 x("x", "x", this, _x),
                                 A("A", "A", this, _A),
                                 B("B", "B", this, _B),
                                 C("C", "C", this, _C),
                                 D("D", "D", this, _D),
                                 E("E", "E", this, _E),
                                 F("F", "F", this, _F)
{
}

CB2Pdf::CB2Pdf(const CB2Pdf& other, const char* name) : RooAbsPdf(other, name),
                                                        x("x", this, other.x),
                                                        A("A", this, other.A),
                                                        B("B", this, other.B),
                                                        C("C", this, other.C),
                                                        D("D", this, other.D),
                                                        E("E", this, other.E),
                                                        F("F", this, other.F)
{
}

Double_t CB2Pdf::evaluate() const
{
  Double_t t = (x - A) / B;
  if (C < 0) {
    t = -t;
  }

  Double_t absAlpha = TMath::Abs(C);
  Double_t absAlpha2 = TMath::Abs(E);

  if (t >= -absAlpha and t < absAlpha2) {
    return TMath::Exp(-0.5 * t * t);
  }

  if (t < -absAlpha) {
    Double_t a = TMath::Power(D / absAlpha, D) * TMath::Exp(-0.5 * absAlpha * absAlpha);
    Double_t b = D / absAlpha - absAlpha;
    return a / TMath::Power(b - t, D);
  }

  if (t >= absAlpha2) {
    Double_t c = TMath::Power(F / absAlpha2, F) * TMath::Exp(-0.5 * absAlpha2 * absAlpha2);
    Double_t d = F / absAlpha2 - absAlpha2;
    return c / TMath::Power(d + t, F);
  }

  return 0.;
}
