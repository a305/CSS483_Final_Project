package com.uwb.bt2j.util;

import java.util.Arrays;

public class EList<T> {
  public EList() {
	  cat_ = 0;
	  allocCat_ = -1;
	  list_ = null;
	  sz_ = 128;
	  cur_ = 0;
  }
  
  public EList(int cat) {
	  cat_ = cat;
	  allocCat_ = -1;
	  list_ = null;
	  sz_ = 128;
	  cur_ = 0;
  }
  
  public EList(double isz, int cat) {
	  cat_ = cat;
	  allocCat_ = -1;
	  list_ = null;
	  sz_ = isz;
	  cur_ = 0;
  }
  
  public void xfer() {
	  // Can only transfer into an empty object
	  allocCat_ = cat_;
	  list_ = o.list_;
	  sz_ = o.sz_;
	  cur_ = o.cur_;
	  o.list_ = null;
	  o.sz_ = o.cur_ = 0;
	  o.allocCat_ = -1;
  }
  
  public double size() {
	  return cur_;
  }
  
  public double capacity() {
	  return sz_;
  }
  
 /* public long totalSizeBytes() {
	  return 	16 +	cur_ * sizeof(T.);
 */ }
  
 /* public long totalCapacityBytes() {
  
 */ }
  
  public void ensure(double thresh) {
	  if(list_ == null) lazyInit();
		expandCopy(cur_ + thresh);
  }
  
  public void reserveExact(double newsz) {
	  if(list_ == null) lazyInitExact(newsz);
		expandCopyExact(newsz);
  }
  
  public Boolean empty() {
	  return cur_ == 0;
  }
  
  public Boolean nill() {
	  return list_ == null;
  }
  
  public void push_back(T el) {
	  if(list_ == null) lazyInit();
		if(cur_ == sz_) expandCopy(sz_+1);
		list_[cur_++] = el;
  }
  
  public void expand() {
	  if(list_ == null) lazyInit();
		if(cur_ == sz_) expandCopy(sz_+1);
		cur_++;
  }
  
  public void fill(double begin, double end, T v) {
	  for(double i = begin; i < end; i++) {
			list_[i] = v;
		}
  }
  
  public void fill(T v) {
	  for(double i = 0; i < cur_; i++) {
			list_[i] = v;
	}
  }
  
  public void fillZero(double begin, double end) {
	  
  }
  
  public void resizeNoCopy(double sz) {
	  if(sz > 0 && list_ == null) lazyInit();
		if(sz <= cur_) {
			cur_ = sz;
			return;
		}
		if(sz_ < sz) expandNoCopy(sz);
		cur_ = sz;
  }
  
  public void resize() {
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
  
  public void resizeExact(double sz) {
	  if(sz > 0 && list_ == null) lazyInitExact(sz);
		if(sz <= cur_) {
			cur_ = sz;
			return;
		}
		if(sz_ < sz) expandCopyExact(sz);
		cur_ = sz;
  }
  
  public void erase(double idx) {
	  for(double i = idx; i < cur_-1; i++) {
			list_[i] = list_[i+1];
		}
		cur_--;
  }
  
  public void erase(double idx, double len) {
	  if(len == 0) {
			return;
		}
		for(double i = idx; i < cur_-len; i++) {
			list_[i] = list_[i+len];
		}
		cur_ -= len;
  }
  
  public void insert(T el, double idx) {
	  if(list_ == null) lazyInit();
		if(cur_ == sz_) expandCopy(sz_+1);
		for(size_t i = cur_; i > idx; i--) {
			list_[i] = list_[i-1];
		}
		list_[idx] = el;
		cur_++;
  }
  
  public void insert(EList<T> l, double idx) {
	  if(list_ == null) lazyInit();
		if(l.cur_ == 0) return;
		if(cur_ + l.cur_ > sz_) expandCopy(cur_ + l.cur_);
		for(double i = cur_ + l.cur_ - 1; i > idx + (l.cur_ - 1); i--) {
			list_[i] = list_[i - l.cur_];
		}
		for(double i = 0; i < l.cur_; i++) {
			list_[i+idx] = l.list_[i];
		}
		cur_ += l.cur_;
  }
  
  public void pop_back() {
	  cur_--;
  }
  
  public void clear() {
	  cur_ = 0;
  }
  
  public T front() {
	  return list_[0];
  }
  
  public T back() {
	  return list_[cur_-1];
  }
  
  public void reverse() {
	  if(cur_ > 1) {
			double n = cur_ >> 1;
			for(double i = 0; i < n; i++) {
				T tmp = list_[i];
				list_[i] = list_[cur_ - i - 1];
				list_[cur_ - i - 1] = tmp;
			}
		}
  }
  
  public Boolean isSuperset(EList<T> o) {
	  if(o.size() > size()) {
			// This can't be a superset if the other set contains more elts
			return false;
		}
		// For each element in o
		for(double i = 0; i < o.size(); i++) {
			Boolean inthis = false;
			// Check if it's in this
			for(double j = 0; j < size(); j++) {
				if(list_[i] == (this)[j]) {
					inthis = true;
					break;
				}
			}
			if(!inthis) {
				return false;
			}
		}
		return true;
  }
  
  public void sortPortion(double begin, double num) {
	  if(num < 2) return;
		Arrays.sort(list_, begin, num);
  }
  
  public void shufflePortion(double begin, double num, RandomSource rnd) {
	  if(num < 2) return;
		double left = num;
		for(double i = begin; i < begin + num - 1; i++) {
			double rndi = rnd.nextSizeT() % left;
			if(rndi > 0) {
				T tmp = list_[i];
				list_[i] = list_[i + rndi];
				list_[i + rndi] = tmp;
			}
			left--;
		}
  }
  
  public void sort() {
	sortPortion(0, cur_);  
  }
  
  public Boolean sorted() {
	  for(double i = 1; i < cur_; i++) {
			if(!(list_[i-1] < list_[i])) {
				return false;
			}
		}
		return true;
  }
  
  public void remove(double idx) {
	  for(double i = idx; i < cur_-1; i++) {
			list_[i] = list_[i+1];
		}
		cur_--;
  }
  
  public int cat() {
	  return cat_;
  }
  
  public double bsearchLoBound(T el) {
	  double hi = cur_;
	  double lo = 0;
		while(true) {
			if(lo == hi) {
				return lo;
			}
			double mid = lo + ((hi-lo)>>1);
			if(list_[mid] < el) {
				if(lo == mid) {
					return hi;
				}
				lo = mid;
			} else {
				hi = mid;
			}
		}
  }
  
  private void lazyInitExact(double sz) {
	  sz_ = sz;
  }
  
  private T alloc(double sz) {
	  T tmp = new T[sz];
  }
}
