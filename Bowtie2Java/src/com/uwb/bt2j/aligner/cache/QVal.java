package com.uwb.bt2j.aligner.cache;

import com.uwb.bt2j.util.IndexTypes;

public class QVal {
	protected long i_;
	protected long rangen_;
	protected long eltn_;
	
	public QVal() {
		reset();
	}
	
	public int offset() {
		return i_;
	}
	
	public long numRanges() {
		return rangen_;
	}
	
	public long numElts() {
		return eltn_;
	}
	
	public boolean empty() {
		return numRanges() == 0;
	}
	
	public boolean valid() {
		return rangen_ != IndexTypes.OFF_MASK;
	}
	
	public void reset() {
		i_ = 0;
		rangen_ = eltn_ = IndexTypes.OFF_MASK;
	}
	
	public void init(long i, long ranges, long elts) {
		i_ = i; rangen_ = ranges; eltn_ = elts;
	}
	
	public void addRange(long numElts) {
		rangen_++;
		eltn_ += numElts;
	}
}
