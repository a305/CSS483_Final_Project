package com.uwb.bt2j.aligner;

public class AlignmentScore {
  public long score_;
  public int basesAligned;
  public int edits_;
  public long ns_;
  public long gaps_;
  
  public AlignmentScore() {
    reset();
    invalidate();
  }
  
  public AlignmentScore(long score, int basesAligned, int edits, long ns, long gaps) {
  score_ = score;
		basesAligned_ = basesAligned;
		edits_ = edits;
		ns_ = ns;
		gaps_ = gaps;
  }
  
  public void reset() {
  score_ = basesAligned_ = edits_ = ns_ = gaps_ = 0;
  }
  
  public static AlignmentScore INVALID() {
    AlignmentScore s;
		s.invalidate();
		return s;
  }
  
  public Boolean valid() {
    return score_ != MIN_I64;
  }
  
  public void invalidate() {
  score_ = MIN_I64;
		edits_ = basesAligned_ = Integer.MIN_VALUE;
		ns_ = gaps_ = 0;
  }
  
  public void incNs(int nceil) {
    if(++ns_ > nceil) {
			invalidate();
		}
  }
  
  public long score() {
    return score_;
  }
  
  public long gaps() {
    return gaps_;
  }
  
  public long penalty() {
    return -score_;
  }
  
  public long ns() {
    return ns_;
  }
  
  public int basesAligned() {
    return basesAligned_;
  }
  
  public int nedit() {
    return edits_;
  }
}
