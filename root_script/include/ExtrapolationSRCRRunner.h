#pragma once

#include "TString.h"
#include <iostream>

class ExtrapolationSRCRRunner{

public:
  ExtrapolationSRCRRunner();
  ~ExtrapolationSRCRRunner() = default;

  void ExtrapolationSRCRFunction();
  void SubtractMCfromDATA();

private:
  TString inpathCR;
  TString inpathSR;
  TString outpath;
  TString outpathSUB;
  TString inpathCRA;

};
