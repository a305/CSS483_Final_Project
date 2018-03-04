package com.uwb.bt2j.aligner.cache;

import com.uwb.bt2j.indexer.BTDnaString;

public class QKey {
	
	public long seq;
	public int len;
	
	public QKey(){
		reset();
	}
	
	public QKey(BTDnaString s) {
		init(s);
	}
	
	public Boolean init(BTDnaString s) {
		seq = 0;
		len = (double)s.length();
		if(len > 32) {
			len = 0xffffffff;
			return false; // wasn't cacheable
		} else {
			// Rightmost char of 's' goes in the least significant bitpair
			for(int i = 0; i < 32 && i < s.length(); i++) {
				int c = (int)s.get(i);
				if(c == 4) {
					len = 0xffffffff;
					return false;
				}
				seq = (seq << 2) | s.get(i);
			}
			return true; // was cacheable
		}
	}
	
	public void toString(BTDnaString s) {
		s.resize(len);
		long sq = seq;
		for(int i = (int) ((len)-1); i >= 0; i--) {
			s.set((double)(sq & 3), i);
			sq >>= 2;
		}
	}
	
	public Boolean cacheable() {
		return len != 0xffffffff;
	}
	
	public void reset() {
		seq = 0; len = 0xffffffff;
	}
}
