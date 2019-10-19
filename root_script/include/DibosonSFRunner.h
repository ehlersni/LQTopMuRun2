#pragma once

#include "TString.h"
#include <iostream>

class DibosonSFRunner{

public:
  DibosonSFRunner();
  ~DibosonSFRunner() = default;

  void DibosonSF();

private:
  TString inpathDi;
  TString inpathDibtag;
  TString outpathDi;
  TString outpathDibtag;



};
