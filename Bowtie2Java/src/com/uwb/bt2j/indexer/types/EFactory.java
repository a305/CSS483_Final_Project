package com.uwb.bt2j.indexer.types;

public class EFactory <T>{
	protected EList<T> l_;
	
	public EFactory(int isz, int cat) {
		l_ = new EList<T>(isz,cat);
	}
	
	public EFactory(int cat) {
		l_ = new EList<T>(cat);
	}
	
	public void clear() {
		l_.clear();
	}
	
	public int alloc() {
		l_.expand();
		return l_.size() - 1;
	}
	
	public int size() {
		return l_.size();
	}
	
	public void resize(int sz) {
		l_.resize(sz);
	}
	
	public boolean empty() {
		return size() == 0;
	}
	
	public void pop() {
		l_.resize(l_.size() - 1);
	}
}
