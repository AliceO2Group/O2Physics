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
//
/// \author Contact: iarsene@cern.ch, i.c.arsene@fys.uio.no
/// \brief Implementation of AnalysisCut
/// \file AnalysisCut.cxx

#include "PWGDQ/Core/AnalysisCut.h"

#include <Framework/Logger.h>

#include <TNamed.h>

#include <Rtypes.h>

#include <vector>

std::vector<int> AnalysisCut::fgUsedVars = {};

//____________________________________________________________________________
AnalysisCut::AnalysisCut(const char* name, const char* title) : TNamed(name, title),
                                                                fCuts()
{
  //
  // named constructor
  //
}

//____________________________________________________________________________
AnalysisCut& AnalysisCut::operator=(const AnalysisCut& c)
{
  //
  // assignment
  //
  if (this != &c) {
    TNamed::operator=(c);
    fCuts = c.fCuts;
  }
  return (*this);
}

//____________________________________________________________________________
AnalysisCut::AnalysisCut(const AnalysisCut& c) : TNamed(c)
{
  //
  // copy constructor
  //
  if (this != &c) {
    fCuts = c.fCuts;
  }
}

//____________________________________________________________________________
AnalysisCut::~AnalysisCut() = default;

//____________________________________________________________________________
void AnalysisCut::printCuts()
{
  LOGF(info, "**************** AnalysisCut::printCuts");
}
