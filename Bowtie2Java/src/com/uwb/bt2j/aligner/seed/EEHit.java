package com.uwb.bt2j.aligner.seed;

import com.uwb.bt2j.aligner.Edit;

public class EEHit {
	public long top;
	public long bot;
	public Edit     e1;
	public Edit     e2;
	boolean     fw;
	long  score;
	
	public EEHit() {
		reset();
	}
	
	public void reset() {
		top = bot = 0;
		fw = false;
		e1.reset();
		e2.reset();
		score = Long.MIN_VALUE;
	}
	
	public void init(
			long top_,
			long bot_,
			Edit e1_,
			Edit e2_,
			boolean fw_,
			long score_) {
		top = top_; bot = bot_;
		if(e1_ != null) {
			e1 = e1_;
		} else {
			e1.reset();
		}
		if(e2_ != null) {
			e2 = e2_;
		} else {
			e2.reset();
		}
		fw = fw_;
		score = score_;
	}
	
	public int mms() {
		if     (e2.inited()) return 2;
		else if(e1.inited()) return 1;
		else                 return 0;
	}
	
	public int ns() {
		int ns = 0;
		if(e1.inited() && e1.hasN()) {
			ns++;
			if(e2.inited() && e2.hasN()) {
				ns++;
			}
		}
		return ns;
	}
	
	public int refns() {
		int ns = 0;
		if(e1.inited() && e1.chr == 'N') {
			ns++;
			if(e2.inited() && e2.chr == 'N') {
				ns++;
			}
		}
		return ns;
	}
	
	public boolean empty() {
		return bot <= top;
	}
	
	public long size() {
		return bot - top;
	}
}
