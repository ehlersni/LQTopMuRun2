#include "TString.h"
#include <iostream>
#include "../include/SubtractMCfromDATARunner.h"

using namespace std;

SubtractMCfromDATARunner::SubtractMCfromDATARunner(){

  inpathCR = "/nfs/dust/cms/user/ehlersni/LQTopMuRun2/fullselection/controlregionbtag/";
  inpathSR = "/nfs/dust/cms/user/ehlersni/LQTopMuRun2/fullselection/signalregion/";
  outpath = "/nfs/dust/cms/user/ehlersni/LQTopMuRun2/fullselection/extrapolationsrcr/";
  outpathSUB = "/nfs/dust/cms/user/ehlersni/LQTopMuRun2/fullselection/substractmcfromdata/";


}
