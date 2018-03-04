package com.uwb.bt2j.aligner.cache;

import com.uwb.bt2j.indexer.IndexTypes;

public class SAVal {
	public int topf;
	public int topb;
	public int i;
	public int len;
	
	public SAVal() {
		len = (int) IndexTypes.OFF_MASK;
	}
	
	public boolean valid() {
		return len != IndexTypes.OFF_MASK;
	}
	
	public void init(int tf, int tb, int ii, int ln) {
		topf = tf;
		topb = tb;
		i = ii;
		len = ln;
	}
}