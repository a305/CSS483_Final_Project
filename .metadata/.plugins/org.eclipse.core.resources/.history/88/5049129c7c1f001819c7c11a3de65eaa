package com.uwb.bt2j.util.types;

public class EHeap<T> {
	protected EList<T> l_;
	
	public void insert(T o) {
		int pos = l_.size();
		l_.push_back(o);
		while(pos > 0) {
			int parent = (pos-1) >> 1;
			if(l_[pos] < l_[parent]) {
				T tmp = new T(l_[pos]);
				l_[pos] = l_[parent];
				l_[parent] = tmp;
				pos = parent;
			} else break;
		}
	}
	
	public T top() {
		return l_[0];
	}
	
	public T pop() {
		T ret = l_[0];
		l_[0] = l_[l_.size()-1];
		l_.resize(l_.size()-1);
		int cur = 0;
		while(true) {
			int c1 = ((cur+1) << 1) - 1;
			int c2 = c1 + 1;
			if(c2 < l_.size()) {
				if(l_[c1] < l_[cur] && l_[c1] <= l_[c2]) {
					T tmp(l_[c1]);
					l_[c1] = l_[cur];
					l_[cur] = tmp;
					cur = c1;
				} else if(l_[c2] < l_[cur]) {
					T tmp(l_[c2]);
					l_[c2] = l_[cur];
					l_[cur] = tmp;
					cur = c2;
				} else {
					break;
				}
			} else if(c1 < l_.size()) {
				if(l_[c1] < l_[cur]) {
					T tmp(l_[c1]);
					l_[c1] = l_[cur];
					l_[cur] = tmp;
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
