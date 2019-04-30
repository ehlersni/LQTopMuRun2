#include "TString.h"
#include <iostream>
#include "../include/EleTriggerEffRunner.h"
#include "../include/ExtrapolationSRCRRunner.h"
// #include "../include/SubtractMCfromDATARunner.h"

using namespace std;


int main(){

   EleTriggerEffRunner TriggerRunner;
   TriggerRunner.EleTriggerEff();


  // ExtrapolationSRCRRunner ExtrapolationRunner;
  //// ExtrapolationRunner.ExtrapolationSRCRFunction();
  // ExtrapolationRunner.SubtractMCfromDATA();

  // SubtractMCfromDATARunner Substractor;
  // Substractor.SubtractMCfromDATA();


}
