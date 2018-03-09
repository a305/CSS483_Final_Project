package com.uwb.bt2j.indexer.types;

import java.util.Arrays;

public class ESet<T> {
	private int cat_;
	private T[] list_;
	private int sz_;
	private int cur_;
	
  public ESet() {
	  cat_ = 0;
	  list_ = null;
	  sz_ = 0;
	  cur_ = 0;
  }
  
  public ESet(int cat) {
	  cat_ = cat;
	  list_ = null;
	  sz_ = 0;
	  cur_ = 0;
  }
  
  public ESet(int isz, int cat) {
	  cat_ = cat;
	  list_ = null;
	  sz_ = isz;
	  cur_ = 0;
  }
  
  public ESet(ESet<T> o, int cat) {
	  cat_ = cat;
	  list_ = null;
  }
  
  public void xfer(ESet<T> o) {
	  // Can only transfer into an empty object
	  list_ = o.list_;
	  sz_ = o.sz_;
	  cur_ = o.cur_;
	  o.list_ = null;
	  o.sz_ = o.cur_ = 0;
  }
  
  public int size() {
	  return cur_;
  }
  
  public boolean empty() {
	  return cur_ == 0;
  }
  
  public boolean nill() {
	  return list_ == null;
  }
  
  public boolean insert(T el) {
	  int i = 0;
		if(cur_ == 0) {
			insert(el, 0);
			return true;
		}
		if(cur_ < 16) {
			// Linear scan
			i = scanLoBound(el);
		} else {
			// Binary search
			i = bsearchLoBound(el);
		}
		if(i < cur_ && list_[i] == el) return false;
		insert(el, i);
		return true;
  }
  
  public boolean contains(T el) {
	  if(cur_ == 0) {
			return false;
		}
		else if(cur_ == 1) {
			return el == list_[0];
		}
		int i;
		if(cur_ < 16) {
			// Linear scan
			i = scanLoBound(el);
		} else {
			// Binary search
			i = bsearchLoBound(el);
		}
		return i != cur_ && list_[i] == el;
  }
  
  public void remove(T el) {
	  int i;
		if(cur_ < 16) {
			// Linear scan
			i = scanLoBound(el);
		} else {
			// Binary search
			i = bsearchLoBound(el);
		}
		erase(i);
  }
  
  public void resize(int sz) {
	  if(sz <= cur_) return;
		if(sz_ < sz) expandCopy(sz);
  }
  
  private void erase(int idx) {
	  for(int i = idx; i < cur_-1; i++) {
			list_[i] = list_[i+1];
		}
		cur_--;
  }
  
  private void insert(T el, int idx) {
	  if(cur_ == sz_) {
			expandCopy(sz_+1);
		}
		for(int i = cur_; i > idx; i--) {
			list_[i] = list_[i-1];
		}
		list_[idx] = el;
		cur_++;
  }
  
  public T get(int idx) {
	  return list_[idx];
  }
  
  public void clear() {
	  cur_ = 0;
  }
  
  private boolean sorted() {
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
  
  public void setCat(int cat) {
	  cat_ = cat;
  }
  
  private int scanLoBound(T el) {
	  for(int i = 0; i < cur_; i++) {
			if(!(list_.get(i) < el)) {
				// Shouldn't be equal
				return i;
			}
		}
		return cur_;
  }
  
  private int bsearchLoBound(T el) {
	  int hi = cur_;
		int lo = 0;
		while(true) {
			if(lo == hi) {
				return lo;
			}
			int mid = lo + ((hi-lo)>>1);
			if(list_.get(mid) < el) {
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
		T[] tmp = new T[];
		for(int i = 0; i < cur_; i++) {
			tmp[i] = list_[i];
		}
		list_ = tmp;
		sz_ = newsz;
  }
}
