// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <CCDB/BasicCCDBManager.h>
#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "tables.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

float fReal_fTPCSignalN(float mbb0R, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((0.019869f * mbb0R) + (0.0012031f * a1pt) + (-0.0031766f * atgl) + (-0.0058023f * atgl * mbb0R) + (0.00087494f * a1pt * mbb0R) + (0.0020074f * side) + (-0.0010434f * a1pt * a1pt) + (0.011812f)) * occ / 1.e3f +       //
         ((0.009032f * mbb0R) + (0.0011737f * a1pt) + (-0.0010241f * atgl) + (-0.0075789f * atgl * mbb0R) + (0.00029324f * a1pt * mbb0R) + (0.00052475f * side) + (-0.00045413f * a1pt * a1pt) + (0.0024879f)) * fOccTPCN +       //
         ((0.004255f * mbb0R) + (0.0011954f * a1pt) + (0.0054092f * atgl) + (-0.0033655f * atgl * mbb0R) + (0.00052243f * a1pt * mbb0R) + (-0.0002969f * side) + (-0.00074909f * a1pt * a1pt) + (-0.0075754f)) * fTrackOccMeanN + //
         ((-0.07925f * mbb0R) + (-0.03737f * a1pt) + (0.0017054f * atgl) + (0.093686f * atgl * mbb0R) + (0.023925f * a1pt * mbb0R) + (-0.0083407f * side) + (0.00336f * a1pt * a1pt) + (1.0461f));
};

float fReal_qMax0_R0(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((0.036419f * mbb0R1) + (0.0026895f * a1pt) + (-0.0060575f * atgl) + (0.0014132f * atgl * mbb0R1) + (-0.0019824f * a1pt * mbb0R1) + (0.0027807f * side) + (-0.00035062f * a1pt * a1pt) + (0.00035707f)) * occ / 1.e3f         //
         + ((0.018872f * mbb0R1) + (0.0025347f * a1pt) + (0.00042624f * atgl) + (-0.009566f * atgl * mbb0R1) + (-0.0013989f * a1pt * mbb0R1) + (0.00046257f * side) + (-0.00035991f * a1pt * a1pt) + (-0.0035327f)) * fOccTPCN         //
         + ((-0.0013489f * mbb0R1) + (-0.0012116f * a1pt) + (0.00095043f * atgl) + (0.0012604f * atgl * mbb0R1) + (0.0010987f * a1pt * mbb0R1) + (-0.0015823f * side) + (-2.0682e-05f * a1pt * a1pt) + (-0.0034631f)) * fTrackOccMeanN //
         + ((0.14845f * mbb0R1) + (-0.059922f * a1pt) + (0.088597f * atgl) + (-0.046773f * atgl * mbb0R1) + (0.047535f * a1pt * mbb0R1) + (-0.0023642f * side) + (0.0041058f * a1pt * a1pt) + (0.81337f));
};

float fReal_qTot0_R0(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((0.05379f * mbb0R1) + (0.0078357f * a1pt) + (-0.0045679f * atgl) + (-0.0055742f * atgl * mbb0R1) + (-0.0045914f * a1pt * mbb0R1) + (0.0034146f * side) + (-0.001886f * a1pt * a1pt) + (0.0015936f)) * occ / 1.e3f            //
         + ((0.019984f * mbb0R1) + (0.00074509f * a1pt) + (-0.0014629f * atgl) + (-0.014526f * atgl * mbb0R1) + (-0.00027911f * a1pt * mbb0R1) + (0.00076815f * side) + (-0.00027441f * a1pt * a1pt) + (0.003094f)) * fOccTPCN         //
         + ((0.0074452f * mbb0R1) + (-0.0012576f * a1pt) + (0.0076864f * atgl) + (-0.0057527f * atgl * mbb0R1) + (0.0010178f * a1pt * mbb0R1) + (-0.00060637f * side) + (-0.00021891f * a1pt * a1pt) + (-0.0079396f)) * fTrackOccMeanN //
         + ((-0.14994f * mbb0R1) + (-0.025556f * a1pt) + (-0.0027218f * atgl) + (0.094753f * atgl * mbb0R1) + (0.024869f * a1pt * mbb0R1) + (-0.0059085f * side) + (-0.00031196f * a1pt * a1pt) + (1.0614f));
};

float fReal_qMaxTot0(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((-0.011295f * mbb0R1) + (-0.001087f * a1pt) + (-0.0031759f * atgl) + (0.007727f * atgl * mbb0R1) + (-6.3106e-06f * a1pt * mbb0R1) + (-0.00051982f * side) + (0.00053087f * a1pt * a1pt) + (-0.00089624f)) * occ / 1.e3f     //
         + ((-0.0011264f * mbb0R1) + (0.0015426f * a1pt) + (0.0018148f * atgl) + (0.0034001f * atgl * mbb0R1) + (-0.0014034f * a1pt * mbb0R1) + (-0.00026979f * side) + (-0.00017022f * a1pt * a1pt) + (-0.0042007f)) * fOccTPCN      //
         + ((-0.0052403f * mbb0R1) + (-0.0014282f * a1pt) + (-0.0029884f * atgl) + (0.0036363f * atgl * mbb0R1) + (0.00052755f * a1pt * mbb0R1) + (-0.0011482f * side) + (0.00036221f * a1pt * a1pt) + (0.0021188f)) * fTrackOccMeanN //
         + ((0.29443f * mbb0R1) + (-0.016326f * a1pt) + (0.084094f * atgl) + (-0.13113f * atgl * mbb0R1) + (0.021916f * a1pt * mbb0R1) + (0.0053678f * side) + (0.00047673f * a1pt * a1pt) + (0.73264f));
};

float fReal_qMax1_R0(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((0.021162f * mbb0R1) + (0.0047311f * a1pt) + (-0.0050134f * atgl) + (-0.0020799f * atgl * mbb0R1) + (-0.0021859f * a1pt * mbb0R1) + (0.0013058f * side) + (-0.00097268f * a1pt * a1pt) + (-0.00027099f)) * occ / 1.e3f        //
         + ((0.010442f * mbb0R1) + (0.0032625f * a1pt) + (0.00035052f * atgl) + (-0.0061975f * atgl * mbb0R1) + (-0.00095599f * a1pt * mbb0R1) + (0.00028564f * side) + (-0.00060877f * a1pt * a1pt) + (-0.0034378f)) * fOccTPCN        //
         + ((0.00072113f * mbb0R1) + (0.0012692f * a1pt) + (0.0019225f * atgl) + (-4.9502e-05f * atgl * mbb0R1) + (0.0002153f * a1pt * mbb0R1) + (-0.00078837f * side) + (-0.00062101f * a1pt * a1pt) + (-0.0054873f)) * fTrackOccMeanN //
         + ((0.057368f * mbb0R1) + (-0.078491f * a1pt) + (0.030685f * atgl) + (-0.026226f * atgl * mbb0R1) + (0.053421f * a1pt * mbb0R1) + (-0.0068194f * side) + (0.012138f * a1pt * a1pt) + (0.93283f));
};

float fReal_qTot1_R0(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((0.022051f * mbb0R1) + (0.0054705f * a1pt) + (-0.0027645f * atgl) + (-0.0080202f * atgl * mbb0R1) + (-0.0010462f * a1pt * mbb0R1) + (0.0019932f * side) + (-0.0021057f * a1pt * a1pt) + (0.0095212f)) * occ / 1.e3f         //
         + ((0.0091488f * mbb0R1) + (0.0033904f * a1pt) + (-0.0012161f * atgl) + (-0.0071031f * atgl * mbb0R1) + (-0.00043533f * a1pt * mbb0R1) + (0.00059263f * side) + (-0.00097877f * a1pt * a1pt) + (0.00075281f)) * fOccTPCN     //
         + ((0.0065562f * mbb0R1) + (0.0035568f * a1pt) + (0.0055333f * atgl) + (-0.0034575f * atgl * mbb0R1) + (-0.00084672f * a1pt * mbb0R1) + (-0.00021743f * side) + (-0.0011436f * a1pt * a1pt) + (-0.010477f)) * fTrackOccMeanN //
         + ((-0.076084f * mbb0R1) + (-0.047299f * a1pt) + (-0.0033743f * atgl) + (0.078502f * atgl * mbb0R1) + (0.029237f * a1pt * mbb0R1) + (-0.011353f * side) + (0.0062695f * a1pt * a1pt) + (1.0308f));
};

float fReal_qMaxTot1(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((0.0011394f * mbb0R1) + (0.0028785f * a1pt) + (-0.0026975f * atgl) + (0.0068239f * atgl * mbb0R1) + (-0.0032394f * a1pt * mbb0R1) + (-0.00024361f * side) + (0.0001246f * a1pt * a1pt) + (-0.010253f)) * occ / 1.e3f      //
         + ((0.00055046f * mbb0R1) + (-0.00051361f * a1pt) + (0.0015337f * atgl) + (0.001474f * atgl * mbb0R1) + (-0.00059542f * a1pt * mbb0R1) + (-0.00033977f * side) + (0.00042607f * a1pt * a1pt) + (-0.0031188f)) * fOccTPCN   //
         + ((-0.004931f * mbb0R1) + (-0.0036774f * a1pt) + (-0.0017014f * atgl) + (0.0023117f * atgl * mbb0R1) + (0.0018169f * a1pt * mbb0R1) + (-0.00049096f * side) + (0.0008295f * a1pt * a1pt) + (0.0040567f)) * fTrackOccMeanN //
         + ((0.15081f * mbb0R1) + (-0.016878f * a1pt) + (0.029106f * atgl) + (-0.10469f * atgl * mbb0R1) + (0.017923f * a1pt * mbb0R1) + (0.003732f * side) + (0.0033021f * a1pt * a1pt) + (0.87985f));
};

float fReal_qMax2_R0(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((0.011954f * mbb0R1) + (0.0046224f * a1pt) + (-0.0051857f * atgl) + (-0.00081821f * atgl * mbb0R1) + (-0.0023337f * a1pt * mbb0R1) + (0.0013618f * side) + (-0.00093544f * a1pt * a1pt) + (-6.2866e-05f)) * occ / 1.e3f      //
         + ((0.0072282f * mbb0R1) + (0.0028697f * a1pt) + (0.00079917f * atgl) + (-0.0042087f * atgl * mbb0R1) + (-0.0015179f * a1pt * mbb0R1) + (-2.2e-05f * side) + (-0.0004591f * a1pt * a1pt) + (-0.0032937f)) * fOccTPCN          //
         + ((0.0011864f * mbb0R1) + (0.0017477f * a1pt) + (0.0019628f * atgl) + (-0.00099482f * atgl * mbb0R1) + (3.9852e-05f * a1pt * mbb0R1) + (-0.0006403f * side) + (-0.00070453f * a1pt * a1pt) + (-0.0045286f)) * fTrackOccMeanN //
         + ((0.019613f * mbb0R1) + (-0.12656f * a1pt) + (0.019093f * atgl) + (-0.035163f * atgl * mbb0R1) + (0.081274f * a1pt * mbb0R1) + (-0.0038576f * side) + (0.026742f * a1pt * a1pt) + (0.97812f));
};

float fReal_qTot2_R0(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((0.013053f * mbb0R1) + (0.0052803f * a1pt) + (-0.002833f * atgl) + (-0.0067229f * atgl * mbb0R1) + (-0.0015083f * a1pt * mbb0R1) + (0.0014132f * side) + (-0.0017529f * a1pt * a1pt) + (0.0045371f)) * occ / 1.e3f           //
         + ((0.0046757f * mbb0R1) + (0.0019132f * a1pt) + (-0.00092483f * atgl) + (-0.0040224f * atgl * mbb0R1) + (-0.00017565f * a1pt * mbb0R1) + (0.00014709f * side) + (-0.00054281f * a1pt * a1pt) + (0.00059496f)) * fOccTPCN     //
         + ((0.0024731f * mbb0R1) + (0.0012129f * a1pt) + (0.0041885f * atgl) + (-0.0030637f * atgl * mbb0R1) + (0.00053027f * a1pt * mbb0R1) + (-0.00064605f * side) + (-0.00073676f * a1pt * a1pt) + (-0.0048447f)) * fTrackOccMeanN //
         + ((-0.03766f * mbb0R1) + (-0.017712f * a1pt) + (-0.0048423f * atgl) + (0.094793f * atgl * mbb0R1) + (0.016855f * a1pt * mbb0R1) + (-0.0063542f * side) + (0.00066842f * a1pt * a1pt) + (0.9918f));
};

float fReal_qMaxTot2(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((0.0038681f * mbb0R1) + (0.0054964f * a1pt) + (-0.0011662f * atgl) + (0.004913f * atgl * mbb0R1) + (-0.0037561f * a1pt * mbb0R1) + (-1.4101e-05f * side) + (-0.00061975f * a1pt * a1pt) + (-0.010448f)) * occ / 1.e3f       //
         + ((0.0019623f * mbb0R1) + (0.00028139f * a1pt) + (0.0021717f * atgl) + (-0.00046209f * atgl * mbb0R1) + (-0.0010646f * a1pt * mbb0R1) + (-0.00023386f * side) + (0.00022279f * a1pt * a1pt) + (-0.0031211f)) * fOccTPCN     //
         + ((-0.001751f * mbb0R1) + (-0.0019385f * a1pt) + (-0.0013717f * atgl) + (0.0012036f * atgl * mbb0R1) + (0.00056936f * a1pt * mbb0R1) + (-3.1006e-06f * side) + (0.00060137f * a1pt * a1pt) + (0.0014547f)) * fTrackOccMeanN //
         + ((0.088494f * mbb0R1) + (-0.078423f * a1pt) + (0.011554f * atgl) + (-0.1185f * atgl * mbb0R1) + (0.0446f * a1pt * mbb0R1) + (9.0607e-05f * side) + (0.020254f * a1pt * a1pt) + (0.95417f));
};

float fReal_qMax3_R0(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((0.006657f * mbb0R1) + (0.0027675f * a1pt) + (-0.0075368f * atgl) + (-2.647e-05f * atgl * mbb0R1) + (-0.00045219f * a1pt * mbb0R1) + (0.00022622f * side) + (-0.00059026f * a1pt * a1pt) + (0.002745f)) * occ / 1.e3f  //
         + ((0.0072501f * mbb0R1) + (0.0028788f * a1pt) + (0.0018944f * atgl) + (-0.0050865f * atgl * mbb0R1) + (-0.0012857f * a1pt * mbb0R1) + (-1.9644e-05f * side) + (-0.00039888f * a1pt * a1pt) + (-0.004457f)) * fOccTPCN  //
         + ((0.0093193f * mbb0R1) + (0.01292f * a1pt) + (-0.001179f * atgl) + (0.0034663f * atgl * mbb0R1) + (-0.005494f * a1pt * mbb0R1) + (-0.00023943f * side) + (-0.0031454f * a1pt * a1pt) + (-0.015039f)) * fTrackOccMeanN //
         + ((0.031517f * mbb0R1) + (-0.12968f * a1pt) + (0.013073f * atgl) + (-0.013104f * atgl * mbb0R1) + (0.064576f * a1pt * mbb0R1) + (0.0037192f * side) + (0.035161f * a1pt * a1pt) + (0.94643f));
};

float fReal_qTot3_R0(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((0.01012f * mbb0R1) + (0.0045905f * a1pt) + (-0.0060295f * atgl) + (-0.0045184f * atgl * mbb0R1) + (-0.0011767f * a1pt * mbb0R1) + (0.00018135f * side) + (-0.0013776f * a1pt * a1pt) + (0.0030707f)) * occ / 1.e3f    //
         + ((0.0034954f * mbb0R1) + (0.0005948f * a1pt) + (-0.0006857f * atgl) + (-0.0032487f * atgl * mbb0R1) + (0.00063149f * a1pt * mbb0R1) + (0.00010616f * side) + (-0.00014182f * a1pt * a1pt) + (0.00016516f)) * fOccTPCN //
         + ((0.012009f * mbb0R1) + (0.01333f * a1pt) + (0.0026538f * atgl) + (-0.0024246f * atgl * mbb0R1) + (-0.0051234f * a1pt * mbb0R1) + (9.6145e-05f * side) + (-0.0033937f * a1pt * a1pt) + (-0.016958f)) * fTrackOccMeanN //
         + ((-0.055659f * mbb0R1) + (0.0060279f * a1pt) + (-0.024762f * atgl) + (0.1707f * atgl * mbb0R1) + (0.010545f * a1pt * mbb0R1) + (0.0053018f * side) + (-0.00751f * a1pt * a1pt) + (0.96645f));
};

float fReal_qMaxTot3(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((9.8263e-05f * mbb0R1) + (0.0019565f * a1pt) + (-0.0020669f * atgl) + (0.0051386f * atgl * mbb0R1) + (-0.0015196f * a1pt * mbb0R1) + (0.00022648f * side) + (0.0001314f * a1pt * a1pt) + (-0.0044507f)) * occ / 1.e3f       //
         + ((0.0039512f * mbb0R1) + (0.0027043f * a1pt) + (0.0027183f * atgl) + (-0.0021564f * atgl * mbb0R1) + (-0.0019754f * a1pt * mbb0R1) + (-9.8236e-05f * side) + (-0.0003849f * a1pt * a1pt) + (-0.0048656f)) * fOccTPCN       //
         + ((0.0015951f * mbb0R1) + (0.0030524f * a1pt) + (-0.0021155f * atgl) + (0.004186f * atgl * mbb0R1) + (-0.0020921f * a1pt * mbb0R1) + (-0.00049161f * side) + (-0.00043331f * a1pt * a1pt) + (-0.0032405f)) * fTrackOccMeanN //
         + ((0.078899f * mbb0R1) + (-0.14603f * a1pt) + (0.034436f * atgl) + (-0.17042f * atgl * mbb0R1) + (0.05807f * a1pt * mbb0R1) + (-0.00053445f * side) + (0.045416f * a1pt * a1pt) + (0.98711f));
};

float fReal_fTPCSignalN_mad(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((0.0022684f * mbb0R1) + (0.0010606f * a1pt) + (-6.7709e-06f * atgl) + (-0.0012901f * atgl * mbb0R1) + (1.809e-05f * a1pt * mbb0R1) + (7.4229e-05f * side) + (-0.00031708f * a1pt * a1pt) + (8.2362e-05f)) * occ / 1.e3f         //
         + ((-0.00015603f * mbb0R1) + (-0.0007137f * a1pt) + (-0.00029446f * atgl) + (-3.6583e-05f * atgl * mbb0R1) + (0.00065671f * a1pt * mbb0R1) + (0.00013936f * side) + (0.00012159f * a1pt * a1pt) + (0.0011442f)) * fOccTPCN       //
         + ((-0.00024348f * mbb0R1) + (-0.00082385f * a1pt) + (0.00053475f * atgl) + (-0.00024263f * atgl * mbb0R1) + (0.00033506f * a1pt * mbb0R1) + (1.009e-05f * side) + (0.00016525f * a1pt * a1pt) + (0.00071934f)) * fTrackOccMeanN //
         + ((0.026912f * mbb0R1) + (-0.0022536f * a1pt) + (-0.0018352f * atgl) + (-0.0050977f * atgl * mbb0R1) + (-0.0041775f * a1pt * mbb0R1) + (0.0010817f * side) + (0.00032f * a1pt * a1pt) + (0.042308f));
};

float fReal_qMax0_R0_mad(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((0.0032739f * mbb0R1) + (-0.00083832f * a1pt) + (-0.0013899f * atgl) + (0.00031148f * atgl * mbb0R1) + (3.3722e-05f * a1pt * mbb0R1) + (0.00035792f * side) + (0.00022816f * a1pt * a1pt) + (0.0020259f)) * occ / 1.e3f        //
         + ((0.0011825f * mbb0R1) + (-0.00087854f * a1pt) + (-0.00060302f * atgl) + (-8.301e-06f * atgl * mbb0R1) + (0.000603f * a1pt * mbb0R1) + (0.00014523f * side) + (0.00021131f * a1pt * a1pt) + (0.00089158f)) * fOccTPCN         //
         + ((-0.00039752f * mbb0R1) + (-0.00067871f * a1pt) + (-0.00042804f * atgl) + (0.0003994f * atgl * mbb0R1) + (0.0003911f * a1pt * mbb0R1) + (7.6255e-05f * side) + (0.00019167f * a1pt * a1pt) + (0.00064082f)) * fTrackOccMeanN //
         + ((0.05757f * mbb0R1) + (-0.00010976f * a1pt) + (-0.0053188f * atgl) + (-0.010717f * atgl * mbb0R1) + (-0.0075574f * a1pt * mbb0R1) + (0.003877f * side) + (0.00013875f * a1pt * a1pt) + (0.049855f));
};

float fReal_qTot0_R0_mad(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((0.010376f * mbb0R1) + (0.0020723f * a1pt) + (0.00023197f * atgl) + (-0.003087f * atgl * mbb0R1) + (-0.0014009f * a1pt * mbb0R1) + (0.00062246f * side) + (-0.00057047f * a1pt * a1pt) + (-0.00074635f)) * occ / 1.e3f          //
         + ((0.0039974f * mbb0R1) + (-5.2525e-05f * a1pt) + (-0.00054452f * atgl) + (-0.0017733f * atgl * mbb0R1) + (0.00038107f * a1pt * mbb0R1) + (0.00018845f * side) + (4.6742e-05f * a1pt * a1pt) + (4.6125e-05f)) * fOccTPCN        //
         + ((0.0011989f * mbb0R1) + (1.3873e-05f * a1pt) + (0.0001106f * atgl) + (-0.00039192f * atgl * mbb0R1) + (-0.00012419f * a1pt * mbb0R1) + (-4.8745e-05f * side) + (2.8361e-05f * a1pt * a1pt) + (-0.00013149f)) * fTrackOccMeanN //
         + ((0.061978f * mbb0R1) + (-2.0672e-05f * a1pt) + (-0.0029687f * atgl) + (-0.020046f * atgl * mbb0R1) + (-0.0068452f * a1pt * mbb0R1) + (0.0060925f * side) + (-0.001742f * a1pt * a1pt) + (0.068078f));
};

float fReal_qMaxTot0_mad(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((0.00081584f * mbb0R1) + (0.00060523f * a1pt) + (-0.00077368f * atgl) + (0.00058048f * atgl * mbb0R1) + (-0.00055331f * a1pt * mbb0R1) + (9.2515e-05f * side) + (-7.2075e-05f * a1pt * a1pt) + (0.00040485f)) * occ / 1.e3f       //
         + ((0.00078366f * mbb0R1) + (0.00032795f * a1pt) + (-0.00035044f * atgl) + (-0.00010332f * atgl * mbb0R1) + (-0.00015074f * a1pt * mbb0R1) + (3.965e-05f * side) + (-6.6383e-05f * a1pt * a1pt) + (-0.00023799f)) * fOccTPCN       //
         + ((0.00046204f * mbb0R1) + (0.00046443f * a1pt) + (-2.0903e-05f * atgl) + (5.8387e-05f * atgl * mbb0R1) + (-0.00025348f * a1pt * mbb0R1) + (4.4914e-05f * side) + (-9.0766e-05f * a1pt * a1pt) + (-0.00033152f)) * fTrackOccMeanN //
         + ((0.027694f * mbb0R1) + (-0.0011263f * a1pt) + (0.0032024f * atgl) + (-0.01575f * atgl * mbb0R1) + (-0.0016355f * a1pt * mbb0R1) + (0.0042554f * side) + (0.00065342f * a1pt * a1pt) + (0.024298f));
};

float fReal_qMax1_R0_mad(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((0.0045549f * mbb0R1) + (0.00087673f * a1pt) + (-0.0008244f * atgl) + (-0.00084722f * atgl * mbb0R1) + (-0.00076807f * a1pt * mbb0R1) + (0.00026f * side) + (-8.4332e-05f * a1pt * a1pt) + (-0.00046328f)) * occ / 1.e3f      //
         + ((0.00077207f * mbb0R1) + (-0.00066463f * a1pt) + (-0.00028896f * atgl) + (-0.00051748f * atgl * mbb0R1) + (0.00052313f * a1pt * mbb0R1) + (0.00020918f * side) + (0.00012531f * a1pt * a1pt) + (0.00077041f)) * fOccTPCN    //
         + ((-0.001964f * mbb0R1) + (-0.0018987f * a1pt) + (-0.00051532f * atgl) + (0.00090779f * atgl * mbb0R1) + (0.00092927f * a1pt * mbb0R1) + (-3.5317e-06f * side) + (0.00038901f * a1pt * a1pt) + (0.0023439f)) * fTrackOccMeanN //
         + ((0.052395f * mbb0R1) + (0.0051874f * a1pt) + (-0.0046362f * atgl) + (-0.010693f * atgl * mbb0R1) + (-0.0031306f * a1pt * mbb0R1) + (-0.00029953f * side) + (-0.0014752f * a1pt * a1pt) + (0.05555f));
};

float fReal_qTot1_R0_mad(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((0.0051282f * mbb0R1) + (-0.00036785f * a1pt) + (-0.00018229f * atgl) + (-0.0021317f * atgl * mbb0R1) + (-3.9574e-05f * a1pt * mbb0R1) + (0.00046985f * side) + (-6.4708e-05f * a1pt * a1pt) + (0.0022778f)) * occ / 1.e3f   //
         + ((0.0011332f * mbb0R1) + (-0.00057264f * a1pt) + (-0.00048459f * atgl) + (-0.00053628f * atgl * mbb0R1) + (0.00054566f * a1pt * mbb0R1) + (0.00020773f * side) + (6.1131e-05f * a1pt * a1pt) + (0.0013186f)) * fOccTPCN     //
         + ((-0.0014294f * mbb0R1) + (-0.0016395f * a1pt) + (0.00050097f * atgl) + (0.00019716f * atgl * mbb0R1) + (0.0011792f * a1pt * mbb0R1) + (-6.9506e-05f * side) + (0.00028175f * a1pt * a1pt) + (0.0018185f)) * fTrackOccMeanN //
         + ((0.049612f * mbb0R1) + (0.00022664f * a1pt) + (0.00069225f * atgl) + (-0.018209f * atgl * mbb0R1) + (0.001862f * a1pt * mbb0R1) + (-0.00013476f * side) + (-0.0017829f * a1pt * a1pt) + (0.067225f));
};

float fReal_qMaxTot1_mad(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((0.00030083f * mbb0R1) + (0.00032633f * a1pt) + (-0.00089213f * atgl) + (0.00050071f * atgl * mbb0R1) + (-0.00038379f * a1pt * mbb0R1) + (0.00013246f * side) + (-3.0354e-05f * a1pt * a1pt) + (0.0011596f)) * occ / 1.e3f    //
         + ((0.0011627f * mbb0R1) + (0.001223f * a1pt) + (-0.00011221f * atgl) + (-0.00038303f * atgl * mbb0R1) + (-0.00041462f * a1pt * mbb0R1) + (7.419e-05f * side) + (-0.00031261f * a1pt * a1pt) + (-0.00089038f)) * fOccTPCN      //
         + ((-0.00081109f * mbb0R1) + (-0.00065714f * a1pt) + (4.2318e-05f * atgl) + (0.00013371f * atgl * mbb0R1) + (0.00027371f * a1pt * mbb0R1) + (5.869e-05f * side) + (0.00013454f * a1pt * a1pt) + (0.0010851f)) * fTrackOccMeanN //
         + ((0.022513f * mbb0R1) + (0.0001983f * a1pt) + (0.0016327f * atgl) + (-0.013841f * atgl * mbb0R1) + (-0.0020134f * a1pt * mbb0R1) + (0.00013836f * side) + (0.0010369f * a1pt * a1pt) + (0.027599f));
};

float fReal_qMax2_R0_mad(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((0.0028109f * mbb0R1) + (0.00154f * a1pt) + (-0.0014824f * atgl) + (0.00076005f * atgl * mbb0R1) + (-0.001141f * a1pt * mbb0R1) + (0.00019105f * side) + (-0.00022128f * a1pt * a1pt) + (-0.00068869f)) * occ / 1.e3f          //
         + ((0.00050553f * mbb0R1) + (-0.00021518f * a1pt) + (-0.00018194f * atgl) + (-0.00051275f * atgl * mbb0R1) + (0.00028182f * a1pt * mbb0R1) + (-8.0412e-06f * side) + (1.7056e-05f * a1pt * a1pt) + (0.00046885f)) * fOccTPCN    //
         + ((0.0014012f * mbb0R1) + (0.00080437f * a1pt) + (4.4211e-05f * atgl) + (-0.00021537f * atgl * mbb0R1) + (-0.00062642f * a1pt * mbb0R1) + (4.8299e-05f * side) + (-0.00012079f * a1pt * a1pt) + (-0.001307f)) * fTrackOccMeanN //
         + ((0.025323f * mbb0R1) + (-0.017098f * a1pt) + (-0.0016086f * atgl) + (-0.015473f * atgl * mbb0R1) + (0.010693f * a1pt * mbb0R1) + (0.0010934f * side) + (0.0034083f * a1pt * a1pt) + (0.079031f));
};

float fReal_qTot2_R0_mad(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((0.00030261f * mbb0R1) + (-0.0023321f * a1pt) + (-0.00052643f * atgl) + (-0.0013095f * atgl * mbb0R1) + (0.0012352f * a1pt * mbb0R1) + (0.00021252f * side) + (0.00039201f * a1pt * a1pt) + (0.0041364f)) * occ / 1.e3f        //
         + ((7.5596e-05f * mbb0R1) + (-0.0009219f * a1pt) + (-0.00017619f * atgl) + (-0.0006923f * atgl * mbb0R1) + (0.00048164f * a1pt * mbb0R1) + (6.911e-05f * side) + (0.00014689f * a1pt * a1pt) + (0.0015109f)) * fOccTPCN         //
         + ((0.00026771f * mbb0R1) + (-7.3137e-05f * a1pt) + (0.0010327f * atgl) + (-0.0010829f * atgl * mbb0R1) + (0.00032905f * a1pt * mbb0R1) + (9.1216e-05f * side) + (-6.3908e-05f * a1pt * a1pt) + (-0.0001052f)) * fTrackOccMeanN //
         + ((0.045023f * mbb0R1) + (-0.0012833f * a1pt) + (-0.0011393f * atgl) + (-0.0128f * atgl * mbb0R1) + (-0.00054145f * a1pt * mbb0R1) + (-0.00012961f * side) + (-0.00057735f * a1pt * a1pt) + (0.063766f));
};

float fReal_qMaxTot2_mad(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((-0.00055972f * mbb0R1) + (-5.4245e-05f * a1pt) + (-0.00086381f * atgl) + (0.00045553f * atgl * mbb0R1) + (-0.00011785f * a1pt * mbb0R1) + (0.00012989f * side) + (2.403e-05f * a1pt * a1pt) + (0.0015619f)) * occ / 1.e3f        //
         + ((0.00036557f * mbb0R1) + (0.00044828f * a1pt) + (-0.00033607f * atgl) + (7.8166e-05f * atgl * mbb0R1) + (-6.3217e-05f * a1pt * mbb0R1) + (-3.2699e-07f * side) + (-0.00011895f * a1pt * a1pt) + (-0.0001779f)) * fOccTPCN       //
         + ((0.00026633f * mbb0R1) + (0.00052126f * a1pt) + (-6.9848e-05f * atgl) + (0.00025286f * atgl * mbb0R1) + (-0.00033007f * a1pt * mbb0R1) + (8.4534e-05f * side) + (-0.00011899f * a1pt * a1pt) + (-0.00024994f)) * fTrackOccMeanN //
         + ((0.028394f * mbb0R1) + (0.0062708f * a1pt) + (0.0097613f * atgl) + (-0.025715f * atgl * mbb0R1) + (-0.0040359f * a1pt * mbb0R1) + (-7.9959e-06f * side) + (7.2706e-05f * a1pt * a1pt) + (0.018673f));
};

float fReal_qMax3_R0_mad(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((-0.0016565f * mbb0R1) + (-0.0019009f * a1pt) + (-0.00018012f * atgl) + (-0.00050733f * atgl * mbb0R1) + (0.00084566f * a1pt * mbb0R1) + (0.00033188f * side) + (0.00034497f * a1pt * a1pt) + (0.003469f)) * occ / 1.e3f    //
         + ((0.0013034f * mbb0R1) + (0.001231f * a1pt) + (-0.00029108f * atgl) + (0.00046219f * atgl * mbb0R1) + (-0.00068638f * a1pt * mbb0R1) + (7.6018e-05f * side) + (-0.00027657f * a1pt * a1pt) + (-0.0010937f)) * fOccTPCN     //
         + ((-0.0017241f * mbb0R1) + (-0.0024624f * a1pt) + (-0.00029956f * atgl) + (0.00051832f * atgl * mbb0R1) + (0.0011205f * a1pt * mbb0R1) + (1.466e-06f * side) + (0.00064987f * a1pt * a1pt) + (0.0022073f)) * fTrackOccMeanN //
         + ((0.027891f * mbb0R1) + (-0.042432f * a1pt) + (0.018144f * atgl) + (-0.036612f * atgl * mbb0R1) + (0.01245f * a1pt * mbb0R1) + (-0.00057199f * side) + (0.013407f * a1pt * a1pt) + (0.086444f));
};

float fReal_qTot3_R0_mad(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((2.4179e-05f * mbb0R1) + (-0.0022156f * a1pt) + (0.00063391f * atgl) + (-0.0025629f * atgl * mbb0R1) + (0.0012503f * a1pt * mbb0R1) + (0.00035431f * side) + (0.0003733f * a1pt * a1pt) + (0.0031568f)) * occ / 1.e3f       //
         + ((0.00064537f * mbb0R1) + (-0.00022771f * a1pt) + (-0.00052572f * atgl) + (-0.00013851f * atgl * mbb0R1) + (0.00017452f * a1pt * mbb0R1) + (5.8377e-05f * side) + (8.395e-05f * a1pt * a1pt) + (0.00032015f)) * fOccTPCN   //
         + ((-0.0019497f * mbb0R1) + (-0.0031993f * a1pt) + (0.0010855f * atgl) + (-0.00057412f * atgl * mbb0R1) + (0.0012335f * a1pt * mbb0R1) + (0.00015677f * side) + (0.00076355f * a1pt * a1pt) + (0.0028201f)) * fTrackOccMeanN //
         + ((0.020865f * mbb0R1) + (-0.058609f * a1pt) + (0.022839f * atgl) + (-0.031806f * atgl * mbb0R1) + (0.016076f * a1pt * mbb0R1) + (-0.0010281f * side) + (0.016455f * a1pt * a1pt) + (0.10016f));
};

float fReal_qMaxTot3_mad(float mbb0R1, float a1pt, float atgl, float side, float occ, float fOccTPCN, float fTrackOccMeanN)
{
  return ((0.00031711f * mbb0R1) + (0.0010899f * a1pt) + (-0.00042976f * atgl) + (-0.00021686f * atgl * mbb0R1) + (-0.00044157f * a1pt * mbb0R1) + (6.7083e-05f * side) + (-0.00026787f * a1pt * a1pt) + (0.00018337f)) * occ / 1.e3f      //
         + ((-0.00049422f * mbb0R1) + (-0.00047902f * a1pt) + (-0.00030815f * atgl) + (2.2409e-05f * atgl * mbb0R1) + (0.00035703f * a1pt * mbb0R1) + (-1.1649e-05f * side) + (8.7419e-05f * a1pt * a1pt) + (0.00080965f)) * fOccTPCN      //
         + ((-9.4331e-05f * mbb0R1) + (-9.6348e-05f * a1pt) + (-0.00057261f * atgl) + (0.00085743f * atgl * mbb0R1) + (-0.00016598f * a1pt * mbb0R1) + (1.233e-06f * side) + (0.00010932f * a1pt * a1pt) + (0.00019994f)) * fTrackOccMeanN //
         + ((0.020742f * mbb0R1) + (-0.023041f * a1pt) + (0.030993f * atgl) + (-0.046481f * atgl * mbb0R1) + (0.0041709f * a1pt * mbb0R1) + (0.00054496f * side) + (0.010326f * a1pt * a1pt) + (0.035667f));
};

float clamp(float value, float lo, float hi)
{
  return value < lo ? lo : (value > hi ? hi : value);
}

struct ProduceDynamicExtensionCalib {
  Produces<aod::TracksQACorrectedFull> qacorf;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  ctpRateFetcher fetcher;

  Preslice<aod::Tracks> perColl = aod::track::collisionId;

  using BCs = soa::Join<aod::BCs, aod::Timestamps>;
  using Collisions = soa::Join<aod::Collisions, aod::Mults, aod::EvSels>;
  using Tracks = soa::Join<aod::Tracks, aod::TracksExtra>;

  int runNumber{0};
  int colId{-100};
  int bcId{-100};
  int trkId{-100};
  Collisions::iterator col;
  BCs::iterator bc;
  Tracks::iterator track;

  void process(Collisions const& collisions, BCs const& bcs, /*aod::FT0s const& ft0s,*/ Tracks const& tracks, aod::TracksQAVersion const& tracksQA)
  {
    col = collisions.begin();
    bc = bcs.begin();
    runNumber = bc.runNumber();
    track = tracks.begin();
    qacorf.reserve(tracksQA.size());
    for (auto& trackqa : tracksQA) {
      if (!trackqa.has_track()) {
        qacorf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        continue;
      }
      if (trackqa.trackId() != trkId) {
        track.setCursor(trackqa.trackId());
      }
      if (!track.has_collision()) {
        qacorf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        continue;
      }
      if (track.collisionId() != colId) {
        colId = track.collisionId();
        col.setCursor(colId);
      }
      if (!col.has_foundBC()) {
        qacorf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        continue;
      }
      if (col.foundBCId() != bcId) {
        bc.setCursor(col.foundBCId());
        if (bc.runNumber() != runNumber) {
          runNumber = bc.runNumber();
        }
      }

      float rate = fetcher.fetch(ccdb.service, bc.timestamp(), runNumber, "ZNC hadronic") * 1.e-3;
      float occ = col.trackOccupancyInTimeRange();
      float fOccTPCN = clamp(col.multTPC() / 1100.f, 0.f, 12.f);
      float abs1pt = std::abs(track.signed1Pt());
      float abstgl = std::abs(track.tgl());
      float side = track.tgl() > 0 ? 1.f : 0.f;

      float correction0 = fReal_fTPCSignalN(clamp(50.f / track.tpcSignal(), 0.05f, 1.05f), abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.f);
      auto mbb0R1 = clamp(correction0 * 50.f / track.tpcSignal(), 0.05f, 1.05f);
      float correction1 = fReal_fTPCSignalN(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.);

      qacorf(
        track.tpcSignal() / correction1,
        fReal_qMax0_R0(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.), fReal_qMax1_R0(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.), fReal_qMax2_R0(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.), fReal_qMax3_R0(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.),
        fReal_qMax0_R0_mad(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.), fReal_qMax1_R0_mad(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.), fReal_qMax2_R0_mad(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.), fReal_qMax3_R0_mad(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.),
        fReal_qTot0_R0(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.), fReal_qTot1_R0(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.), fReal_qTot2_R0(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.), fReal_qTot3_R0(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.),
        fReal_qTot0_R0_mad(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.), fReal_qTot1_R0_mad(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.), fReal_qTot2_R0_mad(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.), fReal_qTot3_R0_mad(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.),
        fReal_qMaxTot0(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.), fReal_qMaxTot1(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.), fReal_qMaxTot2(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.), fReal_qMaxTot3(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.),
        fReal_qMaxTot0_mad(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.), fReal_qMaxTot1_mad(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.), fReal_qMaxTot2_mad(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.), fReal_qMaxTot3_mad(mbb0R1, abs1pt, abstgl, side, occ, fOccTPCN, rate / 5.));
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return {adaptAnalysisTask<ProduceDynamicExtensionCalib>(cfgc)};
}
