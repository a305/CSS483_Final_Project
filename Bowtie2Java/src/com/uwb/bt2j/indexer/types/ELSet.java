package com.uwb.bt2j.indexer.types;

import java.util.Arrays;

public class ELSet<T> {
	private int cat_;
	private ESet<T> list_;
	private int sz_;
	private int cur_;
	
	  public ELSet(int cat) {
		  cat_ = cat;
		  sz_ = 128;
		  cur_ = 0;
	  }
	  
	  public ELSet(int isz, int cat) {
		  cat_ = cat;
		  sz_ = isz;
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
	  return (ESet<T>) list_.get(cur_-1);
  }
  
  public ESet<T> front() {
	  return (ESet<T>) list_.get(0);
  }
  
  public ESet<T> get(int idx) {
	  return (ESet<T>) list_.get(idx);
  }
  
  public void setCat(int cat) {
		cat_ = cat;
		if(cat_ != 0) {
			for(int i = 0; i < sz_; i++) {
				list_.get(i).setCat(cat_);
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
				tmp.get(i).xfer(list_.get(i));
			}
		}
		list_ = tmp;
		sz_ = newsz;
  }
  
  private void expandNoCopy(int thresh) {
	  if(thresh <= sz_) return;
		int newsz = (sz_ * 2)+1;
		while(newsz < thresh) newsz *= 2;
		ESet<T> tmp = new ESet(newsz);
		if(cat_ != 0) {
			for(int i = 0; i < newsz; i++) {
				tmp.get(i).setCat(cat_);
			}
		}
		list_ = tmp;
		sz_ = newsz;
  }
}
