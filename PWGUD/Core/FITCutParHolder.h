#ifndef PWGUD_CORE_FITCUTPARHOLDER_H_
#define PWGUD_CORE_FITCUTPARHOLDER_H_

#include <Rtypes.h>

// object to hold customizable FIT bit thresholds
class FITCutParHolder
{
 public:
  // constructor
  FITCutParHolder(bool saveFITbitsets = true,
                  float thr1_FV0A = 50.,
                  float thr1_FT0A = 50.,
                  float thr1_FT0C = 50.,
                  float thr2_FV0A = 100.,
                  float thr2_FT0A = 100.,
                  float thr2_FT0C = 100.)
    : mSaveFITbitsets{saveFITbitsets},
      mThr1FV0A{thr1_FV0A},
      mThr1FT0A{thr1_FT0A},
      mThr1FT0C{thr1_FT0C},
      mThr2FV0A{thr2_FV0A},
      mThr2FT0A{thr2_FT0A},
      mThr2FT0C{thr2_FT0C}
  {
  }

  // setters
  void SetSaveFITbitsets(bool);
  void SetThr1FV0A(float);
  void SetThr1FT0A(float);
  void SetThr1FT0C(float);
  void SetThr2FV0A(float);
  void SetThr2FT0A(float);
  void SetThr2FT0C(float);

  // getters
  bool saveFITbitsets() const;
  float thr1_FV0A() const;
  float thr1_FT0A() const;
  float thr1_FT0C() const;
  float thr2_FV0A() const;
  float thr2_FT0A() const;
  float thr2_FT0C() const;

 private:
  bool mSaveFITbitsets;

  float mThr1FV0A;
  float mThr1FT0A;
  float mThr1FT0C;

  float mThr2FV0A;
  float mThr2FT0A;
  float mThr2FT0C;

  ClassDefNV(FITCutParHolder, 1);
};

#endif // PWGUD_CORE_FITCUTPARHOLDER_H_