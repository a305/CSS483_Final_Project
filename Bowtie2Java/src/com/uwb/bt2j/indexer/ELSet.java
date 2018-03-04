package com.uwb.bt2j.indexer;

import java.util.Arrays;

public class ELSet<T> {
	private int cat_;
	private ESet<T> list_;
	private int sz_;
	private int cur_;
	
  public ELSet() {
	  cat_ = 0;
	  list_ = null;
	  sz_ = 0;
	  cur_ = 0;
  }
  
  public ELSet(int cat) {
	  cat_ = cat;
	  list_ = null;
	  sz_ = 128;
	  cur_ = 0;
  }
  
  public ELSet(int isz, int cat) {
	  cat_ = cat;
	  list_ = null;
	  sz_ = isz;
	  cur_ = 0;
  }
  
  public ELSet(ELSet<T> o) {
	  cat_ = 0;
	  list_ = null;
	  sz_ = 0;
	  cur_ = 0;
  }
  
  public ELSet(ELSet<T> o, int cat) {
	  cat_ = cat;
	  list_ = null;
	  sz_ = 0;
	  cur_ = 0;
  }
  
  public void xfer(ELSet<T> o) {
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
  
  public void expand() {
	  if(cur_ == sz_) expandCopy(sz_+1);
		cur_++;
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
  
  public void clear() {
	  cur_ = 0;
  }
  
  public ESet<T> back() {
	  return list_[cur_-1];
  }
  
  public ESet<T> front() {
	  return list_[0];
  }
  
  public ESet<T> get(int idx) {
	  return list_[idx];
  }
  
  public void setCat(int cat) {
		cat_ = cat;
		if(cat_ != 0) {
			for(int i = 0; i < sz_; i++) {
				list_[i].setCat(cat_);
			}
		}
  }
  
  public int cat() {
	  return cat_;
  }
  
  private void expandCopy(int thresh) {
	  if(thresh <= sz_) return;
		int newsz = (sz_ * 2)+1;
		while(newsz < thresh) newsz *= 2;
		ESet<T> tmp;
		if(list_ != null) {
			for(int i = 0; i < cur_; i++) {
				tmp[i].xfer(list_[i]);
			}
		}
		list_ = tmp;
		sz_ = newsz;
  }
  
  private void expandNoCopy(int thresh) {
	  if(thresh <= sz_) return;
		int newsz = (sz_ * 2)+1;
		while(newsz < thresh) newsz *= 2;
		ESet<T> tmp;
		list_ = tmp;
		sz_ = newsz;
  }
}
