package com.uwb.bt2j.aligner;

public class RedundantAlignment {
  protected EList<ESet<RedundantCell>> cells_;
  public class RedundantCell {
  public long rfid;
  public Boolean fw;
  public long rfoff;
  public double rdoff;
  
    public RedundantCell() {
      rfid = 0;
      fw = true;
      rfoff = 0;
      rdoff = 0;
    }
    
    public RedundantCell(long rfid_, Boolean fw_, long rfoff_, double rdoff_) {
    init(rfid_, fw_, rfoff_, rdoff_);
    }
    
    public void init(long rfid_, Boolean fw_, long rfoff_, double rdoff_) {
      rfid  = rfid_;
		fw    = fw_;
		rfoff = rfoff_;
		rdoff = rdoff_;
    }
  }
  
  public RedundantAlignment(int cat) {
    cells_ = cat;
  }
  
  public void reset() {
    cells_.clear();
  }
  
  public void init(double npos) {
    cells_.resize(npos);
		for(double i = 0; i < npos; i++) {
			cells_[i].clear();
		}
  }
  
  public void add(AlignmentResult res) {
  
  }
  
  public Boolean overlap(AlignmentResult res) {
  
  }
}
