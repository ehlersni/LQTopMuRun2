#pragma once

#include "TString.h"
#include <iostream>

class FakeRateRunner{

public:
  FakeRateRunner();
  ~FakeRateRunner() = default;

  void EleFakeRate();
  void MuFakeRate();

private:

  TString pathele;
  TString pathele_out;
  TString pathmu;
  TString pathmu_out;


};
