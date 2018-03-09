package com.uwb.bt2j.indexer.types;

import javafx.util.Pair;

public class EMap <K,V>{
	private int cat_;
	private Pair<K,V> list_;
	private int sz_;
	private int cur_;
	
	public EMap(int cat) {
		cat_ = cat;
		sz_ = 128;
		cur_ = 0;
	}
	
	public EMap(int isz, int cat) {
		cat_ = cat;
		sz_ = isz;
		cur_ = 0;
	}
	
	public EMap(EMap<K,V> o) {
	}
	
	public int size() {
		return cur_;
	}
	
	public boolean empty() {
		return cur_ == 0;
	}
	
	public boolean insert(Pair<K,V> el) {
		int i = 0;
		if(cur_ == 0) {
			insert(el, 0);
			return true;
		}
		if(cur_ < 16) {
			// Linear scan
			i = scanLoBound(el.getKey());
		} else {
			// Binary search
			i = bsearchLoBound(el.getKey());
		}
		if(list_[i] == el) return false; // already there
		insert(el, i);
		return true; // not already there
	}
	
	public boolean contains(K el) {
		if(cur_ == 0) return false;
		else if(cur_ == 1) return el == list_[0].getKey();
		int i;
		if(cur_ < 16) {
			// Linear scan
			i = scanLoBound(el);
		} else {
			// Binary search
			i = bsearchLoBound(el);
		}
		return i != cur_ && list_[i].getKey() == el;
	}
	
	public boolean containsEx(K el, int i) {
		if(cur_ == 0) return false;
		else if(cur_ == 1) {
			i = 0;
			return el == list_[0].getKey();
		}
		if(cur_ < 16) {
			// Linear scan
			i = scanLoBound(el);
		} else {
			// Binary search
			i = bsearchLoBound(el);
		}
		return i != cur_ && list_[i].getKey() == el;
	}
	
	public void remove(K el) {
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
	
	public Pair<K,V> get(int i) {
		return list_[i];
	}
	
	public void clear() {
		cur_ = 0;
	}
	
	private int scanLoBound(K el) {
		for(int i = 0; i < cur_; i++) {
			if(!(list_[i].getKey() < el)) {
				// Shouldn't be equal
				return i;
			}
		}
		return cur_;
	}
	
	private int bsearchLoBound(K el) {
		int hi = cur_;
		int lo = 0;
		while(true) {
			if(lo == hi) {
				return lo;
			}
			int mid = lo + ((hi-lo)>>1);
			if(list_[mid].getKey() < el) {
				if(lo == mid) {
					return hi;
				}
				lo = mid;
			} else {
				hi = mid;
			}
		}
	}
	
	public boolean sorted() {
		return true;
	}
	
	public void insert(Pair<K,V> el, int idx) {
		if(cur_ == sz_) {
			expandCopy(sz_+1);
		}
		for(int i = cur_; i > idx; i--) {
			list_[i] = list_[i-1];
		}
		list_[idx] = el;
		cur_++;
	}
	
	public void erase(int idx) {
		for(int i = idx; i < cur_-1; i++) {
			list_[i] = list_[i+1];
		}
		cur_--;
	}
	
	public void expandCopy(int thresh) {
		if(thresh <= sz_) return;
		int newsz = sz_ * 2;
		while(newsz < thresh) newsz *= 2;
		Pair<K, V> tmp;
		for(int i = 0; i < cur_; i++) {
			tmp[i] = list_[i];
		}
		list_ = tmp;
		sz_ = newsz;
	}
	
	
}
