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

#ifndef EXPPDF
#define EXPPDF

#include "RooAbsCategory.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooCategoryProxy.h"
#include "RooRealProxy.h"

class ExpPdf : public RooAbsPdf
{
 public:
  ExpPdf(){};
  ExpPdf(const char* name, const char* title,
         RooAbsReal& _x,
         RooAbsReal& _A,
         RooAbsReal& _B);
  ExpPdf(const ExpPdf& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new ExpPdf(*this, newname); }
  inline virtual ~ExpPdf() {}

 protected:
  RooRealProxy x;
  RooRealProxy A;
  RooRealProxy B;

  Double_t evaluate() const;

 private:
  ClassDef(ExpPdf, 1) // Your description goes here...
};

#endif
