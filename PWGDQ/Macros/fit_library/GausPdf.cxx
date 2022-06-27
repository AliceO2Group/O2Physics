/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

// Your description goes here...

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
