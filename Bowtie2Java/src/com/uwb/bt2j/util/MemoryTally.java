package com.uwb.bt2j.util;

public class MemoryTally {
	protected long tots_[];
	protected long tot_;
	protected long peaks_[];
	protected long peak_;
	
	 public void add(int cat, long amt) {
		 tots_[cat] += amt;
			tot_ += amt;
			if(tots_[cat] > peaks_[cat]) {
				peaks_[cat] = tots_[cat];
			}
			if(tot_ > peak_) {
				peak_ = tot_;
			}
	  }
	  
	  public void del(int cat, long amt) {
		  tots_[cat] -= amt;
			tot_ -= amt;
	  }
	  
	  public long total() {
	    return tot_;
	  }
	  
	  public long total(int cat) {
	    return tots_[cat];
	  }
	  
	  public long peak() {
	    return peak_;
	  }
	  
	  public long peak(int cat) {
	    return peaks_[cat];
	  }
}
