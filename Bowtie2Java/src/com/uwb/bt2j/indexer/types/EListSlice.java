package com.uwb.bt2j.indexer.types;

public class EListSlice <T>{
	protected int i_;
	protected int len_;
	protected EList<T> list_;
	
	public EListSlice() {
		i_ = 0;
		len_ = 0;
	}
	
	public EListSlice(EList<T> list, int i, int len) {
		i_ = i;
		len_ = len;
		list_ = list;
	}
	
	public void init(EListSlice<T> sl, int first, int last) {
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
