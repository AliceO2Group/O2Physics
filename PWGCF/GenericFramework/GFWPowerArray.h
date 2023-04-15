#ifndef GFWPOWERARRAY__C
#define GFWPOWERARRAY__C
#include <vector>
#include <cmath>
#include <string>
using namespace std;
using std::string;
using std::vector;
typedef vector<int> HarSet;
class GFWPowerArray
{
 public:
  static HarSet GetPowerArray(vector<HarSet> inHarmonics);
  static void PowerArrayTest();

 private:
  static int getHighestHarmonic(const HarSet& inhar);
  static HarSet TrimVec(HarSet hars, int ind);
  static HarSet AddConstant(HarSet hars, int offset);
  static void FlushVectorToMaster(HarSet& masterVector, HarSet& comVec, const int& MaxPower);
  static void RecursiveFunction(HarSet& masterVector, HarSet hars, int offset, const int& MaxPower);
  static void PrintVector(const HarSet& singleSet);
};
#endif
