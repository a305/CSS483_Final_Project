package com.uwb.bt2j.aligner.groupwalk;

public class SARangeWithOffs<T> {
	public long topf;
	public int len;
	public T offs;
	
	public SARangeWithOffs() {
		reset();
	}
	
	public SARangeWithOffs(long tf, int len, T o) {
		init(tf, len, o);
	}
	
	public void init(long tf, int len_, T o) {
		topf = tf;
		len = len_;
		offs = o;
	}
	
	public void reset() {
		topf = Long.MAX_VALUE;
	}
	
	public boolean inited() {
		return topf != Long.MAX_VALUE;
	}
	
	public int size() {
		return offs.size();
	}
}
