package com.uwb.bt2j.aligner.groupwalk;

import com.uwb.bt2j.indexer.IndexTypes;

public class WalkResult {
	public GroupWalkElt elt;
	public long bwrow;
	public long toff;
	
	public WalkResult() {
		reset();
	}
	
	public void reset() {
		elt.reset();
		bwrow = toff = IndexTypes.OFF_MASK;
	}
	
	public void init(
			long oi,  // seed offset index
			boolean f,       // strand
			long r,   // range
			long e,   // element
			long bwr, // BW row
			long len, // length
			long to)  // text offset
			{
		elt.init(oi, f, r, e, len);
		bwrow = bwr;
		toff = to;
	}
}
