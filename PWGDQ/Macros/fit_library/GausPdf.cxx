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

#include "Riostream.h"

#include "GausPdf.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include <math.h>
#include "TMath.h"

ClassImp(GausPdf);

GausPdf::GausPdf(const char* name, const char* title,
                 RooAbsReal& _x,
                 RooAbsReal& _A,
                 RooAbsReal& _B) : RooAbsPdf(name, title),
                                   x("x", "x", this, _x),
                                   A("A", "A", this, _A),
                                   B("B", "B", this, _B)
{
}

GausPdf::GausPdf(const GausPdf& other, const char* name) : RooAbsPdf(other, name),
                                                           x("x", this, other.x),
                                                           A("A", this, other.A),
                                                           B("B", this, other.B)
{
}

Double_t GausPdf::evaluate() const
{
  return TMath::Exp(-0.5 * ((x - A) / B) * ((x - A) / B));
}
