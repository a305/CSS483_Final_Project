package com.uwb.bt2j.indexer;

public class EBitList {
	protected EList<Byte> l_;
	protected int max_;
	
	public EBitList(int isz, int cat) {
		l_ = new EList<Byte>(isz,cat);
		reset();
	}
	
	public EBitList(int cat) {
		l_ = new EList<Byte>(cat);
		reset();
	}
	
	public void clear() {
		reset();
	}
	
	public void reset() {
		l_.clear();
		max_ = Integer.MAX_VALUE;
	}
	
	public void set(int off) {
		resize(off);
		l_[off >> 3] |= (1 << (off & 7));
		if(off > max_ || max_ == Integer.MAX_VALUE) {
			max_ = off;
		}
	}
	
	public boolean test(int off) {
		if((int)(off >> 3) >= l_.size()) {
			return false;
		}
		return (l_[off >> 3] & (1 << (off & 7))) != 0;
	}
	
	public int size() {
		return l_.size();
	}
	
	public void resize(int off) {
		if((int)(off >> 3) >= l_.size()) {
			int oldsz = l_.size();
			l_.resize((off >> 3) + 1);
			for(int i = oldsz; i < l_.size(); i++) {
				l_[i] = 0;
			}
		}
	}
	
	public int max() {
		return max_;
	}
}
