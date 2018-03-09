package com.uwb.bt2j.indexer.types;

import java.util.Arrays;

public class EList<T> {
	private int cat_;
	private int allocCat_;
	private T[] list_;
	private int sz_;
	private int cur_;
	
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
  
  public EList(int isz, int cat) {
	  cat_ = cat;
	  allocCat_ = -1;
	  list_ = null;
	  sz_ = isz;
	  cur_ = 0;
  }
  
  public EList(EList<T> o) {
	  cat_ = 0;
	  allocCat_ = -1;
	  list_ = null;
	  sz_ = 0;
	  cur_ = 0;
  }
  
  public EList(EList<T> o, int cat) {
	  cat_ = cat;
	  allocCat_ = -1;
	  list_ = null;
	  sz_ = 0;
	  cur_ = 0;
  }
  
  public void xfer(EList<T> o) {
	  // Can only transfer into an empty object
	  allocCat_ = cat_;
	  list_ = o.list_;
	  sz_ = o.sz_;
	  cur_ = o.cur_;
	  o.list_ = null;
	  o.sz_ = o.cur_ = 0;
	  o.allocCat_ = -1;
  }
  
  public int size() {
	  return cur_;
  }
  
  public double capacity() {
	  return sz_;
  }
  
  public void ensure(int thresh) {
		expandCopy(cur_ + thresh);
  }
  
  public void reserveExact(int newsz) {
		expandCopyExact(newsz);
  }
  
  public boolean empty() {
	  return cur_ == 0;
  }
  
  public boolean nill() {
	  return list_ == null;
  }
  
  public void push_back(T el) {
		if(cur_ == sz_) expandCopy(sz_+1);
		list_[cur_++] = el;
  }
  
  public void expand() {
		if(cur_ == sz_) expandCopy(sz_+1);
		cur_++;
  }
  
  public void fill(int begin, double end, T v) {
	  for(int i = begin; i < end; i++) {
			list_[i] = v;
		}
  }
  
  public void fill(T v) {
	  for(int i = 0; i < cur_; i++) {
			list_[i] = v;
	}
  }
  
  public void resizeNoCopy(int sz) {
		if(sz <= cur_) {
			cur_ = sz;
			return;
		}
		if(sz_ < sz) expandNoCopy(sz);
		cur_ = sz;
  }
  
  public void resize(int sz) {
		if(sz <= cur_) {
			cur_ = sz;
			return;
		}
		if(sz_ < sz) {
			expandCopy(sz);
		}
		cur_ = sz;
  }
  
  public void resizeExact(int sz) {
		if(sz <= cur_) {
			cur_ = sz;
			return;
		}
		if(sz_ < sz) expandCopyExact(sz);
		cur_ = sz;
  }
  
  public void erase(int idx) {
	  for(int i = idx; i < cur_-1; i++) {
			list_[i] = list_[i+1];
		}
		cur_--;
  }
  
  public void erase(int idx, int len) {
	  if(len == 0) {
			return;
		}
		for(int i = idx; i < cur_-len; i++) {
			list_[i] = list_[i+len];
		}
		cur_ -= len;
  }
  
  public void insert(T el, int idx) {
		if(cur_ == sz_) expandCopy(sz_+1);
		for(int i = cur_; i > idx; i--) {
			list_[i] = list_[i-1];
		}
		list_[idx] = el;
		cur_++;
  }
  
  public T get(int idx) {
	  return list_[idx];
  }
  
  public void insert(EList<T> l, int idx) {
		if(l.cur_ == 0) return;
		if(cur_ + l.cur_ > sz_) expandCopy(cur_ + l.cur_);
		for(int i = cur_ + l.cur_ - 1; i > idx + (l.cur_ - 1); i--) {
			list_[i] = list_[i - l.cur_];
		}
		for(int i = 0; i < l.cur_; i++) {
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
			for(int i = 0; i < n; i++) {
				T tmp = list_[i];
				list_[i] = list_[cur_ - i - 1];
				list_[cur_ - i - 1] = tmp;
			}
		}
  }
  
  public boolean isSuperset(EList<T> o) {
	  if(o.size() > size()) {
			// This can't be a superset if the other set contains more elts
			return false;
		}
		// For each element in o
		for(int i = 0; i < o.size(); i++) {
			boolean inthis = false;
			// Check if it's in this
			for(int j = 0; j < size(); j++) {
				if(list_[i] == list_[j]) {
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
  
  public void sortPortion(int begin, int num) {
	  if(num < 2) return;
		Arrays.sort(list_, begin, num);
  }
  
  public void shufflePortion(int begin, int num, RandomSource rnd) {
	  if(num < 2) return;
		int left = num;
		for(int i = begin; i < begin + num - 1; i++) {
			int rndi = rnd.nextSizeT() % left;
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
  
  public boolean sorted() {
	  for(int i = 1; i < cur_; i++) {
			if(!(list_[i-1] < list_[i])) {
				return false;
			}
		}
		return true;
  }
  
  public void remove(int idx) {
	  for(int i = idx; i < cur_-1; i++) {
			list_[i] = list_[i+1];
		}
		cur_--;
  }
  
  public int cat() {
	  return cat_;
  }
  
  public void set(T o, int idx) {
	  list_[idx] = o;
  }
  
  public void setCat(int cat) {
	  cat_ = cat;
  }
  
  public double bsearchLoBound(T el) {
	  int hi = cur_;
	  int lo = 0;
		while(true) {
			if(lo == hi) {
				return lo;
			}
			int mid = lo + ((hi-lo)>>1);
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
  
  private void expandCopy(int thresh) {
	  if(thresh <= sz_) return;
		int newsz = (sz_ * 2)+1;
		while(newsz < thresh) newsz *= 2;
		expandCopyExact(newsz);
  }
  
  private void expandCopyExact(int newsz) {
	  if(newsz <= sz_) return;
		int cur = cur_;
		T[] tmp = new T[];
		if(list_ != null) {
			for(int i = 0; i < cur_; i++) {
				// Note: operator= is used
				tmp[i] = list_[i];
			}
		}
		list_ = tmp;
		sz_ = newsz;
		cur_ = cur;
  }
  
  private void expandNoCopy(int thresh) {
	  if(thresh <= sz_) return;
		int newsz = (sz_ * 2)+1;
		while(newsz < thresh) newsz *= 2;
		expandNoCopyExact(newsz);
  }
  
  private void expandNoCopyExact(int newsz) {
	T[] tmp = new T[];
	list_ = tmp;
	sz_ = newsz;
  }
}
