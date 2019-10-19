#include "TString.h"
#include <iostream>
#include "../include/EleTriggerEffRunner.h"
#include "../include/ExtrapolationSRCRRunner.h"
#include "../include/DibosonSFRunner.h"
#include "../include/FakeRateRunner.h"
// #include "../include/SubtractMCfromDATARunner.h"

using namespace std;


int main(){

  // FakeRateRunner FRRunner;
  // FRRunner.EleFakeRate();
  // FRRunner.MuFakeRate();
  //
  // DibosonSFRunner DiSFRunner;
  // DiSFRunner.DibosonSF();

   // EleTriggerEffRunner TriggerRunner;
   // TriggerRunner.EleTriggerEff();


  ExtrapolationSRCRRunner ExtrapolationRunner;
  ExtrapolationRunner.ExtrapolationSRCRFunction();
  // ExtrapolationRunner.SubtractMCfromDATA();

  // SubtractMCfromDATARunner Substractor;
  // Substractor.SubtractMCfromDATA();


}
