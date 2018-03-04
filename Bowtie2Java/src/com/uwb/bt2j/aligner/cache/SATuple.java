package com.uwb.bt2j.aligner.cache;

import com.uwb.bt2j.indexer.IndexTypes;
import com.uwb.bt2j.indexer.PListSlice;

public class SATuple {
	public QKey    key;  // sequence key
	public long topf;  // top in BWT index
	public long topb;  // top in BWT' index
	public PListSlice<Long>   offs; // offsets
	
	public SATuple() {
		reset();
	}
	
	public SATuple(QKey k, long tf, long tb, PListSlice<Long> o) {
		init(k, tf, tb, o);
	}
	
	public void init(QKey k, long tf, long tb, PListSlice<Long> o) {
		key = k; topf = tf; topb = tb; offs = o;
	}
	
	public void init(SATuple src, int first, int last) {
		key = src.key;
		topf = (long)(src.topf + first);
		topb = IndexTypes.OFF_MASK; // unknown!
		offs.init(src.offs, first, last);
	}
	
	public void reset() {
		topf = topb = IndexTypes.OFF_MASK; offs.reset(); 
	}
	
	public void setLength(int nlen) {
		offs.setLength(nlen);
	}
	
	public int size() {
		return offs.size();
	}
}