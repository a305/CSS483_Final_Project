package com.uwb.bt2j.aligner.cache;

import com.uwb.bt2j.indexer.BTDnaString;
import com.uwb.bt2j.indexer.BTString;
import com.uwb.bt2j.indexer.EList;

public class AlignmentCacheIFace {
	protected QKey qk_;  // key representation for current read substring
	protected QVal qv_; // pointer to value representation for current read substring
	protected QVal qvbuf_; // buffer for when key is uncacheable but we need a qv
	protected boolean cacheable_; // true iff the read substring currently being aligned is cacheable
	
	protected int rangen_; // number of ranges since last alignment job began
	protected int eltsn_;  // number of elements since last alignment job began

	protected AlignmentCache current_; // cache dedicated to the current read
	protected AlignmentCache local_;   // local, unsynchronized cache
	protected AlignmentCache shared_;  // shared, synchronized cache
	
	public AlignmentCacheIFace(AlignmentCache current, AlignmentCache local, AlignmentCache shared) {
		qv_ = null;
		cacheable_ = false;
		rangen_ = 0;
		eltsn_ = 0;
		current_ = current;
		local_ = local;
		shared_ = shared;
	}
	
	public int beginAlign(BTDnaString seq, BTString qual, QVal qv, boolean getLock) {
		qk_.init(seq);
		if(qk_.cacheable()) {
			qv_ = current_.add(qk_, cacheable_, getLock);
		} else {
			qv_ = qvbuf_;
		}
		if(qv_ == null) {
			resetRead();
 			return -1; // Not in memory
		}
		qv_.reset();
		return 0; // Need to search for it
	}
	
	public QVal finishAlign(boolean getLock) {
		if(!qv_.valid()) {
			qv_.init(0, 0, 0);
		}

		QVal qv = qv_;
		resetRead();
		return qv;
	}
	
	public void nextRead() {
		current_.clear();
		resetRead();
	}
	
	public boolean aligning() {
		return qv_ != null;
	}
	
	public void clear() {
		if(current_ != null) current_.clear();
		if(local_   != null) local_.clear();
		if(shared_  != null) shared_.clear();
	}
	
	public boolean addOnTheFly(
			BTDnaString rfseq, // reference sequence close to read seq
			long topf,            // top in BWT index
			long botf,            // bot in BWT index
			long topb,            // top in BWT' index
			long botb,            // bot in BWT' index
			boolean getLock)      // true -> lock is not held by caller) 
	{
		QKey sak = new QKey(rfseq);
		//assert(sak.cacheable());
		if(current_.addOnTheFly(qv_, sak, topf, botf, topb, botb, getLock)) {
			rangen_++;
			eltsn_ += (botf-topf);
			return true;
		}
		return false;
	}
	
	public void queryQval(
			QVal qv,
			EList<SATuple> satups,
			int nrange,
			int nelt,
			boolean getLock = true)
	{
		current_.queryQval(qv, satups, nrange, nelt, getLock);
	}
	
	public AlignmentCache currentCache() {
		return current_;
	}
	
	public int curNumRanges() {
		return rangen_;
	}
	
	public int curNumElts() {
		return eltsn_;
	}
	
	protected void resetRead() {
		cacheable_ = false;
		rangen_ = eltsn_ = 0;
		qv_ = null;
	}
}
