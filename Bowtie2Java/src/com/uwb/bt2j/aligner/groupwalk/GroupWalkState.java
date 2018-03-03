package com.uwb.bt2j.aligner.groupwalk;

import com.uwb.bt2j.indexer.Ebwt;
import com.uwb.bt2j.indexer.SideLocus;
import com.uwb.bt2j.util.BitPairReference;
import com.uwb.bt2j.util.IndexTypes;
import com.uwb.bt2j.util.types.EList;

import javafx.util.Pair;

public class GroupWalkState <T> {
	public EList<Boolean> masks[];
	public EList<Long> map;
	public SideLocus tloc;
	public SideLocus bloc;
	public long top;
	public long bot;
	public long step;
	protected EList<Long> map_;
	protected long mapi_;
	
	public GroupWalkState(int cat) {
		masks[0].setCat(cat);
		masks[1].setCat(cat);
		masks[2].setCat(cat);
		masks[3].setCat(cat);
	}
	
	public GroupWalkState() {
		map_ = new EList(0, 4);
		reset();
	}
	
	public Pair<Long,Long> init(
			Ebwt ebwt,             // index to walk left in
			BitPairReference ref,  // bitpair-encoded reference
			SARangeWithOffs<T> sa,       // SA range with offsets
			EList<GroupWalkState> sts,       // EList of GWStates for range being advanced
			GroupWalkHit<T> hit,                // Corresponding hit structure
			long range,               // which range is this?
			boolean reportList,              // if true, "report" resolved offsets immediately by adding them to 'res' list
			EList<WalkResult> res,   // EList where resolved offsets should be appended
			long tp,                  // top of range at this step
			long bt,                  // bot of range at this step
			long st,                  // # steps taken to get to this step
			WalkMetrics met) {
		top = tp;
		bot = bt;
		step = st;
		return init(ebwt, ref, sa, sts, hit, range, reportList, res, met);
	}
	
	public Pair<Long,Long> init(
			Ebwt ebwt,             // index to walk left in
			BitPairReference ref,  // bitpair-encoded reference
			SARangeWithOffs<T> sa,       // SA range with offsets
			EList<GroupWalkState> sts,       // EList of GWStates for range being advanced
			GroupWalkHit<T> hit,                // Corresponding hit structure
			long range,               // which range is this?
			boolean reportList,              // if true, "report" resolved offsets immediately by adding them to 'res' list
			EList<WalkResult> res,   // EList where resolved offsets should be appended
			WalkMetrics met) {
		Pair<Long, Long> ret = new Pair(0, 0);
		long trimBegin = 0, trimEnd = 0;
		boolean empty = true; // assume all resolved until proven otherwise
		// Commit new information, if any, to the PListSlide.  Also,
		// trim and check if we're done.
		for(int i = mapi_; i < map_.size(); i++) {
			boolean resolved = (off(i, sa) != IndexTypes.OFF_MASK);
			if(!resolved) {
				// Elt not resolved yet; try to resolve it now
				long bwrow = (long)(top - mapi_ + i);
				long toff = ebwt.tryOffset(bwrow);
				if(toff != IndexTypes.OFF_MASK) {
					// Yes, toff was resolvable
					met.resolves++;
					toff += step;
					setOff(i, toff, sa, met);
					if(!reportList) ret.first++;
				}
			}
			// Is the element resolved?  We ask this regardless of how it was
			// resolved (whether this function did it just now, whether it did
			// it a while ago, or whether some other function outside GroupWalk
			// did it).
			if(off(i, sa) != IndexTypes.OFF_MASK) {
				if(reportList && !hit.reported(map(i))) {
					// Report it
					long toff = off(i, sa);
					res.expand();
					long origBwRow = sa.topf + map(i);
					res.back().init(
						hit.offidx, // offset idx
						hit.fw,     // orientation
						hit.range,  // original range index
						map(i),     // original element offset
						origBwRow,  // BW row resolved
						hit.len,    // hit length
						toff);      // text offset
					hit.setReported(map(i));
					met.reports++;
				}
				// Offset resolved
				if(empty) {
					// Haven't seen a non-empty entry yet, so we
					// can trim this from the beginning.
					trimBegin++;
				} else {
					trimEnd++;
				}
			} else {
				// Offset not yet resolved
				ret.second++;
				trimEnd = 0;
				empty = false;
				// Set the forward map in the corresponding GWHit
				// object to point to the appropriate element of our
				// range
				long bmap = map(i);
				hit.fmap[bmap].first = range;
				hit.fmap[bmap].second = (long)i;
			}
		}
		// Trim from beginning
		mapi_ += trimBegin;
		top += trimBegin;
		if(trimEnd > 0) {
			// Trim from end
			map_.resize(map_.size() - trimEnd);
			bot -= trimEnd;
		}
		if(empty) {
			return ret;
		}
		// Is there a dollar sign in the middle of the range?
		if(ebwt._zOff > top && ebwt._zOff < bot-1) {
			// Yes, the dollar sign is in the middle of this range.  We
			// must split it into the two ranges on either side of the
			// dollar.  Let 'bot' and 'top' delimit the portion of the
			// range prior to the dollar.
			long oldbot = bot;
			bot = ebwt._zOff;
			// Note: might be able to do additional trimming off the
			// end.
			// Create a new range for the portion after the dollar.
			st.expand();
			st.back().reset();
			long ztop = ebwt._zOff+1;
			st.back().initMap(oldbot - ztop);
			for(int i = ztop; i < oldbot; i++) {
				st.back().map_[i - ztop] = map(i-top+mapi_);
			}
			map_.resize(bot - top + mapi_);
			st.back().init(
				ebwt,
				ref,
				sa,
				st,
				hit,
				(long)st.size()-1,
				reportList,
				res,
				ztop,
				oldbot,
				step,
				met);
		}
		// Prepare SideLocus's for next step
		if(bot-top > 1) {
			SideLocus.initFromTopBot(top, bot, ebwt.eh(), ebwt.ebwt(), tloc, bloc);
		} else {
			tloc.initFromRow(top, ebwt.eh(), ebwt.ebwt());
			bloc.invalidate();
		}
		return ret;
	}
	
	public long off(int i, SARangeWithOffs<T> sa) {
		return sa.offs.get(map_[i]);
	}
	
	public long map(int i) {
		return map_[i];
	}
	
	public long mapi() {
		return mapi_;
	}
	
	public int size() {
		return map_.size() - mapi_;
	}
	
	public boolean done() {
		return size() == 0;
	}
	
	public void setOff(
			double i,
			long off,
			SARangeWithOffs<T> sa,
			WalkMetrics met){
		sa.offs[saoff] = off;
	}
	
	public Pair<Long,Long> advance(
			Ebwt ebwt,            // the forward Bowtie index, for stepping left
			BitPairReference ref, // bitpair-encoded reference
			SARangeWithOffs<T> sa,      // SA range with offsets
			GroupWalkHit<T> hit,               // the associated GWHit object
			long range,              // which range is this?
			boolean reportList,             // if true, "report" resolved offsets immediately by adding them to 'res' list
			EList<WalkResult> res,  // EList where resolved offsets should be appended
			EList<GroupWalkState> st,       // EList of GWStates for range being advanced
			GroupWalkState gws,         // temporary storage for masks
			WalkMetrics met,
			PerReadMetrics prm){
		Pair<Long, Long> ret = new Pair(0, 0);
		if(bloc.valid()) {
			// Still multiple elements being tracked
			long upto[4], in[4];
			upto[0] = in[0] = upto[1] = in[1] =
			upto[2] = in[2] = upto[3] = in[3] = 0;
			met.bwops++;
			prm.nExFmops++;
			// Assert that there's not a dollar sign in the middle of
			// this range
			ebwt.mapLFRange(tloc, bloc, bot-top, upto, in, gws.masks);

			boolean first = true;
			long newtop = 0, newbot = 0;
			gws.map.clear();
			for(int i = 0; i < 4; i++) {
				if(in[i] > 0) {
					// Non-empty range resulted
					if(first) {
						// For the first one, 
						first = false;
						newtop = upto[i];
						newbot = newtop + in[i];
						// Range narrowed so we have to look at the masks
						for(int j = 0; j < gws.masks[i].size(); j++) {
							if(gws.masks[i][j]) {
								gws.map.push_back(map_[j+mapi_]);
							}
						}
					} else {
						// For each beyond the first, create a new
						// GWState and add it to the GWState list. 
						// NOTE: this can cause the underlying list to
						// be expanded which in turn might leave 'st'
						// pointing to bad memory.
						st.expand();
						st.back().reset();
						long ntop = upto[i];
						long nbot = ntop + in[i];
						st.back().mapi_ = 0;
						st.back().map_.clear();
						met.branches++;
						// Range narrowed so we have to look at the masks
						for(int j = 0; j < gws.masks[i].size(); j++) {
							if(gws.masks[i][j]) st.back().map_.push_back(map_[j+mapi_]);
						}
						Pair<Long, Long> rret =
						st.back().init(
							ebwt,        // forward Bowtie index
							ref,         // bitpair-encodede reference
							sa,          // SA range with offsets
							st,          // EList of all GWStates associated with original range
							hit,         // associated GWHit object
							(long)st.size()-1, // range offset
							reportList,  // if true, report hits to 'res' list
							res,         // report hits here if reportList is true
							ntop,        // BW top of new range
							nbot,        // BW bot of new range
							step+1,      // # steps taken to get to this new range
							met);        // update these metrics
						ret.first += rret.first;
						ret.second += rret.second;
					}
				}
			}
			mapi_ = 0;
			top = newtop;
			bot = newbot;
			if(!gws.map.empty()) {
				map_ = gws.map;
			}
		} else {
			// Down to one element
			// Sets top, returns char walked through (which we ignore)
			met.bwops++;
			prm.nExFmops++;
			ebwt.mapLF1(top, tloc);
			bot = top+1;
			if(mapi_ > 0) {
				map_[0] = map_[mapi_];
				mapi_ = 0;
			}
			map_.resize(1);
		}
		step++;
		Pair<Long, Long> rret =
		init<S>(
			ebwt,       // forward Bowtie index
			ref,        // bitpair-encodede reference
			sa,         // SA range with offsets
			st,         // EList of all GWStates associated with original range
			hit,        // associated GWHit object
			range,      // range offset
			reportList, // if true, report hits to 'res' list
			res,        // report hits here if reportList is true
			met);       // update these metrics
		ret.first += rret.first;
		ret.second += rret.second;
		return ret;
	}
	
	public void reset() {
		top = bot = step = mapi_ = 0;
		tloc.invalidate();
		bloc.invalidate();
		map_.clear();
	}
	
	public void initMap(int newsz) {
		mapi_ = 0;
		map_.resize(newsz);
		for(int i = 0; i < newsz; i++) {
			map_[i] = (long)i;
		}
	}
	
	public boolean doneResolving(SARangeWithOffs<T> sa) {
		for(int i = mapi_; i < map_.size(); i++) {
			if(sa.offs[map(i)] == IndexTypes.OFF_MASK) return false;
		}
		return true;
	}
}
