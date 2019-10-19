#include "TString.h"
#include <iostream>
#include "../include/ExtrapolationSRCRRunner.h"

using namespace std;

ExtrapolationSRCRRunner::ExtrapolationSRCRRunner(){

  inpathCR = "/nfs/dust/cms/user/ehlersni/LQTopMuRun2/fullselection/controlregionbtag";
  inpathCRA = "/nfs/dust/cms/user/ehlersni/LQTopMuRun2/fullselection/controlregionbtagalpha/";
  inpathSR = "/nfs/dust/cms/user/ehlersni/LQTopMuRun2/fullselection/signalregion/";
  outpath = "/nfs/dust/cms/user/ehlersni/LQTopMuRun2/fullselection/extrapolationsrcr/";
  outpathSUB = "/nfs/dust/cms/user/ehlersni/LQTopMuRun2/fullselection/substractmcfromdata/";


}
