package com.uwb.bt2j.aligner;

public class EList_m128i {
	public int cat_;
	public double sz_;
	public double cur_;
	
	public EList_m128i(int cat) {
		cat_ = cat;
		last_alloc_ = null;
		list_ = null;
		sz_ = 0;
		cur_ = 0;
	}
	
	public double size() {
		return cur_;
	}
	
	public double capacity() {
		return sz_;
	}
	
	public void ensure(double thresh) {
		if(list_ == null) lazyInit();
		expandCopy(cur_ + thresh);
	}
	
	public void reserveExact(double newsz) {
		if(list_ == null) lazyInitExact(newsz);
		expandCopyExact(newsz);
	}
	
	public boolean empty() {
		return cur_ == 0;
	}
	
	public boolean nill() {
		return list_ == null;
	}
	
	public void resize(double sz) {
		if(sz > 0 && list_ == null) lazyInit();
		if(sz <= cur_) {
			cur_ = sz;
			return;
		}
		if(sz_ < sz) {
			expandCopy(sz);
		}
		cur_ = sz;
	}
	
	public void zero() {
		if(cur_ > 0) {
			list_ = 0;
		}
	}
	
	public void resizeNoCopy(double sz) {
		if(sz > 0 && list_ == null) lazyInit();
		if(sz <= cur_) {
			cur_ = sz;
			return;
		}
		if(sz_ < sz) {
			expandNoCopy(sz);
		}
		cur_ = sz;
	}
	
	public void resizeExact(double sz) {
		if(sz > 0 && list_ == null) lazyInitExact(sz);
		if(sz <= cur_) {
			cur_ = sz;
			return;
		}
		if(sz_ < sz) expandCopyExact(sz);
		cur_ = sz;
	}
	
	public void clear() {
		cur_ = 0;
	}
	
	public int cat() {
		return cat_;
	}
	
	public void expandCopy(double thresh) {
		if(thresh <= sz_) return;
		double newsz = (sz_ * 2)+1;
		while(newsz < thresh) newsz *= 2;
		expandCopyExact(newsz);
	}
	
	public void expandNoCopy(double thresh) {
		if(thresh <= sz_) return;
		double newsz = (sz_ * 2)+1;
		while(newsz < thresh) newsz *= 2;
		expandNoCopyExact(newsz);
	}
}
