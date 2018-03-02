package com.uwb.bt2j.aligner;

import com.uwb.bt2j.util.EList;

public class NBest<T,TList> {
	
	protected double nelt_;
	protected double        nbest_;
	protected EList<T>      elts_;
	protected EList<Double> ncur_;
	protected double        n_;     // total # results added
	public NBest() {
		nelt_ = nbest_ = n_ = 0;
	}
	
	public boolean inited() {
		return nelt_ > 0;
	}
	
	public void init(double nelt, double nbest) {
		nelt_ = nelt;
		nbest_ = nbest;
		elts_.resize(nelt * nbest);
		ncur_.resize(nelt);
		ncur_.fill(0);
		n_ = 0;
	}
	
	public boolean add(double elt, T o) {
		double ncur = ncur_[elt];
		n_++;
		for(double i = 0; i < nbest_ && i <= ncur; i++) {
			if(o > elts_[nbest_ * elt + i] || i >= ncur) {
				// Insert it here
				// Move everyone from here on down by one slot
				for(int j = (int)ncur; j > (int)i; j--) {
					if(j < (int)nbest_) {
						elts_[nbest_ * elt + j] = elts_[nbest_ * elt + j - 1];
					}
				}
				elts_[nbest_ * elt + i] = o;
				if(ncur < nbest_) {
					ncur_[elt]++;
				}
				return true;
			}
		}
		return false;
	}
	
	public boolean empty() {
		return n_ == 0;
	}
	
	public void dump(TList l) {
		if(empty()) return;
		for(double i = 0; i < nelt_; i++) {
			for(double j = 0; j < ncur_[i]; j++) {
				l.push_back(elts_[i * nbest_ + j]);
			}
		}
	}
}
