#pragma once

#include "TString.h"
#include <iostream>

class EleTriggerEffRunner{

public:
  EleTriggerEffRunner();
  ~EleTriggerEffRunner() = default;

  void EleTriggerEff();

private:
  TString inpath;
  TString outpath;

};
