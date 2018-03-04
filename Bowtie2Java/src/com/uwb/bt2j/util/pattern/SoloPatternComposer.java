package com.uwb.bt2j.util.pattern;

import com.uwb.bt2j.aligner.Read;
import com.uwb.bt2j.indexer.EList;

import javafx.util.Pair;

public class SoloPatternComposer extends PatternComposer {
	protected double cur_;
	protected EList<PatternSource> src_;
	
	public SoloPatternComposer(PatternParams p, EList<PatternSource> src) {
		super(p);
		cur_ = 0;
		src_ = src;
	}

	@Override
	public void reset() {
		for(int i = 0; i < src_.size(); i++) {
			src[i].reset();
		}
		cur_ = 0;
	}

	public Pair<Boolean, Integer> nextBatch(PerThreadReadBuf pt) {
		double cur = cur_;
		while(cur < src_.size()) {
			// Patterns from srca_[cur_] are unpaired
			Pair<Boolean, Integer> res;
			do {
				res = src_[cur].nextBatch(
					pt,
					true,  // batch A (or pairs)
					true); // grab lock below
			} while(!res.first && res.second == 0);
			if(res.second == 0) {
				ThreadSafe ts(mutex_m);
				if(cur + 1 > cur_) {
					cur_++;
				}
				cur = cur_;
				continue; // on to next pair of PatternSources
			}
			return res;
		}
		return new Pair(true, 0);
	}

	@Override
	public boolean parse(Read ra, Read rb, long rdid) {
		return src[0].parse(ra, rb, rdid);
	}
	
	
}
