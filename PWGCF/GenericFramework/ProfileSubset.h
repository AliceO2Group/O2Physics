/*
Author: Vytautas Vislavicius
Ported to O2: Emil Gorm Nielsen
Extention of Generic Flow (https://arxiv.org/abs/1312.3572)
*/
#ifndef ProfileSubset__H
#define ProfileSubset__H
// Helper function to select a subrange of a TProfile
#include "TProfile.h"
#include "TProfile2D.h"
#include "TError.h"

class ProfileSubset : public TProfile2D
{
 public:
  ProfileSubset(TProfile2D& inpf) : TProfile2D(inpf){};
  ~ProfileSubset(){};
  TProfile* GetSubset(bool onx, const char* name, int fb, int lb, int l_nbins = 0, double* l_binarray = 0);
  void OverrideBinContent(double x, double y, double x2, double y2, double val);
  void OverrideBinContent(double x, double y, double x2, double y2, TProfile2D* sourceProf);
  bool OverrideBinsWithZero(int xb1, int yb1, int xb2, int yb2);

  ClassDef(ProfileSubset, 2);
};
#endif