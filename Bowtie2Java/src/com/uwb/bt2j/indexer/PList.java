package com.uwb.bt2j.indexer;

public class PList<T> {
	protected int cur_;
	protected int curPage_;
	protected EList<T> pages_;
	
	public PList(int cat) {
		cur_ = 0;
		curPage_ = 0;
		pages_ = cat;
	}
	
	public boolean add(Pool p, T o) {
		if(!ensure(p, 1)) return false;
		if(cur_ == PLIST_PER_PAGE) {
			cur_ = 0;
			curPage_++;
		}
		pages_[curPage_][cur_++] = o;
		return true;
	}
	
	public boolean add(Pool p, EList<T> os) {
		if(!ensure(p, os.size())) return false;
		for(int i = 0; i < os.size(); i++) {
			if(cur_ == PLIST_PER_PAGE) {
				cur_ = 0;
				curPage_++;
			}
			pages_[curPage_][cur_++] = os[i];
		}
		return true;
	}
	
	public boolean copy(Pool p, PList<T> src, int i, int len) {
		if(!ensure(p, src.size())) return false;
		for(int i = 0; i < src.size(); i++) {
			if(cur_ == PLIST_PER_PAGE) {
				cur_ = 0;
				curPage_++;
			}
			pages_[curPage_][cur_++] = src[i];
		}
		return true;
	}
	
	public boolean addFill(Pool p, int num, T o) {
		if(!ensure(p, num)) return false;
		for(int i = 0; i < num; i++) {
			if(cur_ == PLIST_PER_PAGE) {
				cur_ = 0;
				curPage_++;
			}
			pages_[curPage_][cur_++] = o;
		}
		return true;
	}
	
	public void clear() {
		pages_.clear();
		cur_ = curPage_ = 0;
	}
	
	public int size() {
		return curPage_ * PLIST_PER_PAGE + cur_;
	}
	
	public boolean empty() {
		return size() == 0;
	}
	
	public T get(int i) {
		int page = i / PLIST_PER_PAGE;
		int elt = i % PLIST_PER_PAGE;
		return pages_[page][elt];
	}
	
	public T back() {
		int page = (size()-1) / PLIST_PER_PAGE;
		int elt = (size()-1) % PLIST_PER_PAGE;
		return pages_[page][elt];
	}
	
	public T last() {
		if(cur_ == 0) {
			return pages_[pages_.size()-2][PLIST_PER_PAGE-1];
		} else {
			return pages_.back()[cur_-1];
		}
	}
	
	public boolean ensure(Pool p, int num) {
		if(num == 0) return true;
		// Allocation of the first page
		if(pages_.size() == 0) {
			if(expand(p) == null) {
				return false;
			}
		}
		int cur = cur_;
		int curPage = curPage_;
		while(cur + num > PLIST_PER_PAGE) {
			if(curPage == pages_.size()-1 && expand(p) == null) {
				return false;
			}
			num -= (PLIST_PER_PAGE - cur);
			cur = 0;
			curPage++;
		}
		return true;
	}
	
	protected T expand(Pool p) {
		T newpage = (T)p.alloc();
		if(newpage == null) {
			return null;
		}
		pages_.push_back(newpage);
		return pages_.back();
	}
}