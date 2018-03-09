package com.uwb.bt2j.indexer.types;

public class ELList <T>{
	
	private int cat_;
	private EList<EList<T>> list_;
	private int sz_;
	private int cur_;
	
	public ELList() {
		cat_ = 0;
	    list_ = null;
	    sz_ = 128;
	    cur_ = 0;
	}
	
  public ELList(int cat) {
    cat_ = cat;
    list_ = null;
    sz_ = 128;
    cur_ = 0;
  }
  
  public ELList(int isz, int cat) {
	  cat_ = cat;
	    list_ = null;
	    sz_ = isz;
	    cur_ = 0;
  }
  
  public ELList(ELList<T> o, int cat) {
	  cat_ = cat;
	    list_ = null;
	    sz_ = 0;
	    cur_ = 0;
  }
  
  public void xfer(ELList<T> o) {
	  list_ = o.list_; // list_ is an array of EList<T>s
		sz_   = o.sz_;
		cur_  = o.cur_;
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
  
  public EList<T> get(int i){
	  return list_.get(i);
  }
  
  public void clear() {
    cur_ = 0;
  }
  
  public EList<T> back() {
	  return list_.get(cur_-1);
  }
  
  public EList<T> front() {
	  return list_.get(0);
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
		EList<EList<T>> tmp = new EList<EList<T>>();
		if(list_ != null) {
			for(int i = 0; i < cur_; i++) {
				tmp.get(i).xfer(list_.get(i));
			}
		}
		list_ = tmp;
		sz_ = newsz;
  }
  
  public void expandNoCopy(int thresh) {
	  if(thresh <= sz_) return;
		int newsz = (sz_ * 2)+1;
		while(newsz < thresh) newsz *= 2;
		EList<EList<T>> tmp = new EList<EList<T>>();
		list_ = tmp;
		sz_ = newsz;
  }
}
