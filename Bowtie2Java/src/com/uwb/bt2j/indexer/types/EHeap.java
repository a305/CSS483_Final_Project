package com.uwb.bt2j.indexer.types;

public class EHeap<T> {
	protected EList<T> l_;
	
	public void insert(T o) {
		int pos = l_.size();
		l_.push_back(o);
		while(pos > 0) {
			int parent = (pos-1) >> 1;
			if(l_.get(pos) < l_.get(parent)) {
				T tmp = new T(l_.get(pos));
				l_.set(l_.get(parent),pos);
				l_.set(tmp,parent);
				pos = parent;
			} else break;
		}
	}
	
	public T top() {
		return l_.get(0);
	}
	
	public T pop() {
		T ret = l_.get(0);
		l_.set(l_.get(l_.size()-1),0);
		l_.resize(l_.size()-1);
		int cur = 0;
		while(true) {
			int c1 = ((cur+1) << 1) - 1;
			int c2 = c1 + 1;
			if(c2 < l_.size()) {
				if(l_.get(c1) < l_.get(cur) && l_.get(c1) <= l_.get(c2)) {
					T tmp = new T(l_.get(c1));
					l_.set(l_.get(cur),c1);
					l_.set(tmp,cur);
					cur = c1;
				} else if(l_.get(c2) < l_.get(cur)) {
					T tmp = new T(l_.get(c2));
					l_.set(l_.get(cur),c2);
					l_.set(tmp ,cur);
					cur = c2;
				} else {
					break;
				}
			} else if(c1 < l_.size()) {
				if(l_.get(c1) < l_.get(cur)) {
					T tmp = new T(l_.get(c1));
					l_.set(l_.get(cur) ,c1);
					l_.set(tmp, cur);
					cur = c1;
				} else {
					break;
				}
			} else {
				break;
			}
		}
		return ret;
	}
	
	public int size() {
		return l_.size();
	}
	
	public boolean empty() {
		return l_.empty();
	}
	
	public void clear() {
		l_.clear();
	}
}
