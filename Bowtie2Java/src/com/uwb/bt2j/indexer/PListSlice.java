package com.uwb.bt2j.indexer;

public class PListSlice <T>{
	protected int i_;
	protected int len_;
	protected PList<T> list_;
	
	public PListSlice() {
		i_ = 0;
		len_ = 0;
	}
	
	public PListSlice(PList<T> list, int i, int len) {
		i_ = i;
		len_ = len;
		list_ = list;
	}
	
	public void init(PListSlice<T> sl, int first, int last) {
		i_ = (sl.i_ + first);
		len_ = (last - first);
		list_ = sl.list_;
	}
	
	public void reset() {
		i_ = len_ = 0;
		list_ = null;
	}
	
	public T get(int i) {
		return list_.get(i + i_);
	}
	
	public boolean valid() {
		return len_ != 0;
	}
	
	public int size() {
		return len_;
	}
	
	public void setLength(int nlen) {
		len_ = nlen;
	}
}
