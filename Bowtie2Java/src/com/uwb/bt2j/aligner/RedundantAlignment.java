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
  long left, right = res.refoff();
	double len = res.readExtentRows();
	double alignmentStart = res.trimmedLeft(true);
	if(!res.fw()) {
		const_cast<AlnRes&>(res).invertEdits();
	}
	EList<Edit> ned = res.ned();
	double nedidx = 0;
	// For each row...
	for(double i = alignmentStart; i < alignmentStart + len; i++) {
		double diff = 1;  // amount to shift to right for next round
		right = left + 1;
		while(nedidx < ned.size() && ned[nedidx].pos == i) {
			if(ned[nedidx].isRefGap()) {
				// Next my_left will be in same column as this round
				diff = 0;
			}
			nedidx++;
		}
		if(i < alignmentStart + len - 1) {
			// See how many inserts there are before the next read
			// character
			double nedidx_next = nedidx;
			while(nedidx_next < ned.size() && ned[nedidx_next].pos == i+1)
			{
				if(ned[nedidx_next].isReadGap()) {
					right++;
				}
				nedidx_next++;
			}
		}
		for(long j = left; j < right; j++) {
			// Add to db
			RedundantCell c(res.refid(), res.fw(), j, i);
		}
		left = right + diff - 1;
	}
	if(!res.fw()) {
		const_cast<AlnRes&>(res).invertEdits();
	}
  }
  
  public Boolean overlap(AlignmentResult res) {
  long left = res.refoff(), right;
	double len = res.readExtentRows();
	double alignmentStart = res.trimmedLeft(true);
	if(!res.fw()) {
		const_cast<AlnRes&>(res).invertEdits();
	}
	EList<Edit> ned = res.ned();
	double nedidx = 0;
	// For each row...
	Boolean olap = false;
	for(double i = alignmentStart; i < alignmentStart + len; i++) {
		double diff = 1;  // amount to shift to right for next round
		right = left + 1;
		while(nedidx < ned.size() && ned[nedidx].pos == i) {
			if(ned[nedidx].isRefGap()) {
				// Next my_left will be in same column as this round
				diff = 0;
			}
			nedidx++;
		}
		if(i < alignmentStart + len - 1) {
			// See how many inserts there are before the next read
			// character
			double nedidx_next = nedidx;
			while(nedidx_next < ned.size() && ned[nedidx_next].pos == i+1)
			{
				if(ned[nedidx_next].isReadGap()) {
					right++;
				}
				nedidx_next++;
			}
		}
		for(long j = left; j < right; j++) {
			// Add to db
			RedundantCell c(res.refid(), res.fw(), j, i);
			if(cells_[i].contains(c)) {
				olap = true;
				break;
			}
		}
		if(olap) {
			break;
		}
		left = right + diff - 1;
	}
	if(!res.fw()) {
		const_cast<AlnRes&>(res).invertEdits();
	}
	return olap;
  }
}
