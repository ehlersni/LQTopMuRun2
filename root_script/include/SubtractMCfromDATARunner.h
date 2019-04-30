#pragma once

#include "TString.h"
#include <iostream>

class SubtractMCfromDATARunner{

public:
  SubtractMCfromDATARunner();
  ~SubtractMCfromDATARunner() = default;

  void SubtractMCfromDATA();

private:
  TString inpathCR;
  TString inpathSR;
  TString outpath;
  TString outpathSUB;

};
