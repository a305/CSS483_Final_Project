package com.uwb.bt2j.aligner.groupwalk;

import com.uwb.bt2j.util.IndexTypes;

public class GroupWalkElt {
	protected long offidx;
	protected boolean fw;
	protected long range;
	protected long elt;
	protected long len;
	
	public GroupWalkElt() {
		reset();
	}
	
	public void reset() {
		offidx = range = elt = len = IndexTypes.OFF_MASK;
		fw = false;
	}
	
	public void init(
			long oi,
			boolean f,
			long r,
			long e,
			long l){
		offidx = oi;
		fw = f;
		range = r;
		elt = e;
		len = l;
	}
}
