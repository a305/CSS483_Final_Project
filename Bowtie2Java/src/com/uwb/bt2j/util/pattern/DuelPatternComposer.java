package com.uwb.bt2j.util.pattern;

import com.uwb.bt2j.aligner.Read;
import com.uwb.bt2j.indexer.EList;

import javafx.util.Pair;

public class DuelPatternComposer extends PatternComposer{
	protected double cur_;
	protected EList<PatternSource> srca_;
	protected EList<PatternSource> srcb_;
	public DuelPatternComposer(PatternParams p, EList<PatternSource> srca, EList<PatternSource> srcb) {
		super(p);
		srca_ = srca;
		srcb_ = srcb;
	}


	@Override
	public void reset() {
		for(int i = 0; i < srca_.size(); i++) {
			srca_[i].reset();
			if(srcb_[i] != null) {
				srcb_[i].reset();
			}
		}
		cur_ = 0;
	}

	@Override
	public Pair<Boolean, Integer> nextBatch(PerThreadReadBuf pt) {
		// 'cur' indexes the current pair of PatternSources
		double cur = cur_;
		while(cur < srca_.size()) {
			if(srcb_[cur] == null) {
				// Patterns from srca_ are unpaired
				Pair<Boolean, Integer> res = srca_[cur].nextBatch(
					pt,
					true,  // batch A (or pairs)
					true); // grab lock below
				boolean done = res.first;
				if(!done && res.second == 0) {
					ThreadSafe ts(mutex_m);
					if(cur + 1 > cur_) cur_++;
					cur = cur_; // Move on to next PatternSource
					continue; // on to next pair of PatternSources
				}
				return new Pair(done, res.second);
			} else {
				Pair<Boolean, Integer> resa, resb;
				// Lock to ensure that this thread gets parallel reads
				// in the two mate files
				{
					ThreadSafe ts(mutex_m);
					resa = srca_[cur].nextBatch(
						pt,
						true,   // batch A
						false); // don't grab lock below
					resb = srcb_[cur].nextBatch(
						pt,
						false,  // batch B
						false); // don't grab lock below
				}
				if(resa.second < resb.second) {
					System.err.println( "Error, fewer reads in file specified with -1 "
					     + "than in file specified with -2");
				} else if(resa.second == 0 && resb.second == 0) {
					ThreadSafe ts(mutex_m);
					if(cur + 1 > cur_) {
						cur_++;
					}
					cur = cur_; // Move on to next PatternSource
					continue; // on to next pair of PatternSources
				} else if(resb.second < resa.second) {
					System.err.println( "Error, fewer reads in file specified with -2 "
					     + "than in file specified with -1");
				}
				return new Pair(resa.first && cur == srca_.size() - 1, resa.second);
			}
		}
		return new Pair(true, 0);
	}

	@Override
	public boolean parse(Read ra, Read rb, long rdid) {
		return srca_[0].parse(ra, rb, rdid);
	}
}
