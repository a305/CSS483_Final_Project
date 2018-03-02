package com.uwb.bt2j.aligner.dp;

public class DynamicProgrammingBtCandidate {
	public double row;
	public double col;
	public long score;
	public int fate;
	public enum BtCandidate {
		BT_CAND_FATE_SUCCEEDED(1),
				BT_CAND_FATE_FAILED(2),
				BT_CAND_FATE_FILT_START(3),     // skipped b/c starting cell already explored
				BT_CAND_FATE_FILT_DOMINATED(4), // skipped b/c it was dominated
				BT_CAND_FATE_FILT_SCORE(5);      // skipped b/c score not interesting anymore
	private int y;
	BtCandidate(int x){y = x;}
	}
	
	public DynamicProgrammingBtCandidate() {
		reset();
	}
	
	public void reset() {
		init(0, 0, 0); 
	}
	
	public DynamicProgrammingBtCandidate(double row_, double col_, long score_) {
		init(row_, col_, score_);
	}
	
	public void init(double row_, double col_, long score_) {
		row = row_;
		col = col_;
		score = score_;
		// 0 = invalid; this should be set later according to what happens
		// before / during the backtrace
		fate = 0; 
	}
	
	public boolean dominatedBy(DynamicProgrammingBtCandidate o) {
		double SQ = 40;
		double rowhi = row;
		double rowlo = o.row;
		if(rowhi < rowlo){
			double tmp =rowhi;
			rowhi = rowlo;
			rowlo = tmp;
		}
		double colhi = col;
		double collo = o.col;
		if(colhi < collo){
			double tmp =colhi;
			colhi = collo;
			collo = tmp;
		}
		return (colhi - collo) <= SQ &&
		       (rowhi - rowlo) <= SQ;
	}
}
