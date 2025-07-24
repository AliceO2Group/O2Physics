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

#ifndef GAUSPDF
#define GAUSPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class GausPdf : public RooAbsPdf
{
 public:
  GausPdf(){};
  GausPdf(const char* name, const char* title,
          RooAbsReal& _x,
          RooAbsReal& _A,
          RooAbsReal& _B);
  GausPdf(const GausPdf& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new GausPdf(*this, newname); }
  inline virtual ~GausPdf() {}

 protected:
  RooRealProxy x;
  RooRealProxy A;
  RooRealProxy B;

  Double_t evaluate() const;

 private:
  ClassDef(GausPdf, 1) // Your description goes here...
};

#endif
