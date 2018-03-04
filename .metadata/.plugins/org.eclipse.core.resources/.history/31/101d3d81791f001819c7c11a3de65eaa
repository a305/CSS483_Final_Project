package com.uwb.bt2j.aligner.groupwalk;

import com.uwb.bt2j.indexer.Ebwt;
import com.uwb.bt2j.util.BitPairReference;
import com.uwb.bt2j.util.IndexTypes;
import com.uwb.bt2j.util.types.EList;

public class GroupWalk2S<T> {
	protected int elt_;
	protected int rep_;
	protected TStateV st_;
	protected GWHit<T> hit_;
	
	public GroupWalk2S() {
		st_ = new EList(8, 4);
		reset();
	}
	
	public void reset() {
		elt_ = rep_ = 0;
	}
	
	public void init(
			Ebwt ebwtFw,         // forward Bowtie index for walking left
			BitPairReference ref,// bitpair-encoded reference
			SARangeWithOffs<T> sa,     // SA range with offsets
			RandomSource rnd,          // pseudo-random generator for sampling rows
			WalkMetrics met)           // update metrics here
	{
		reset();
				// Init GWHit
				hit_.init(sa, 0, false, 0);
				// Init corresponding GWState
				st_.resize(1);
				st_.back().reset();
				long top = sa.topf;
				long bot = (long)(top + sa.size());
				st_.back().initMap(bot-top);
				st_.ensure(4);
				st_.back().init(
					ebwtFw,             // Bowtie index
					ref,                // bitpair-encoded reference
					sa,                 // SA range with offsets
					st_,                // EList<GWState>
					hit_,               // GWHit
					0,                  // range 0
					false,              // put resolved elements into res_?
					null,               // put resolved elements here
					top,                // BW row at top
					bot,                // BW row at bot
					0,                  // # steps taken
					met);               // update metrics here
				elt_ += sa.size();
	}
	
	public void resolveAll(WalkMetrics met, PerReadMetrics prm) {
		WalkResult res; // ignore results for now
		for(int i = 0; i < elt_; i++) {
			advanceElement((long)i, res, met, prm);
		}
	}
	
	public boolean advanceElement(
			long elt,                // element within the range
			Ebwt ebwtFw,          // forward Bowtie index for walking left
			BitPairReference ref, // bitpair-encoded reference
			SARangeWithOffs<T> sa,      // SA range with offsets
			GroupWalkState gws,         // GroupWalk state; scratch space
			WalkResult res,             // put the result here
			WalkMetrics met,            // metrics
			PerReadMetrics prm)         // per-read metrics
	{
		while(sa.offs[elt] == IndexTypes.OFF_MASK) {
			// Get the GWState that contains our element of interest
			int range = hit_.fmap[elt].first;
			st_.ensure(4);
			GWState<T> st = st_[range];
			// Returns a pair of numbers, the first being the number of
			// resolved but unreported offsets found during this advance, the
			// second being the number of as-yet-unresolved offsets.
			st.advance(
				ebwtFw,
				ref,
				sa,
				hit_,
				(long)range,
				false,
				null,
				st_,
				gws,
				met,
				prm);
		}
		// Report it!
		if(!hit_.reported(elt)) {
			hit_.setReported(elt);
		}
		met.reports++;
		res.init(
			0,              // seed offset
			false,          // orientation
			0,              // range
			elt,            // element
			sa.topf + elt,  // bw row
			(long)sa.len, // length of hit
			sa.offs[elt]);  // resolved text offset
		rep_++;
		return true;
	}
	
	public boolean done() {
		return rep_ == elt_;
	}
	
	public int numElts() {
		return elt_;
	}
}
