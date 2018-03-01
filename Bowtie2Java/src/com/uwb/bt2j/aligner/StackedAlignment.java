package com.uwb.bt2j.aligner;

public class StackedAlignment {
  public StackedAlignment() {
  stackRef_ = 
    reset();
  }
  
  public void reset() {
  inited_ = false;
		trimLS_ = trimLH_ = trimRS_ = trimRH_ = 0;
		stackRef_.clear();
		stackRel_.clear();
		stackRead_.clear();
		cigDistMm_ = cigCalc_ = false;
		cigOp_.clear();
		cigRun_.clear();
		mdzCalc_ = false;
		mdzOp_.clear();
		mdzChr_.clear();
		mdzRun_.clear();
  }
}
