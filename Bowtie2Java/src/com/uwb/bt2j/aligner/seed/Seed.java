package com.uwb.bt2j.aligner.seed;

public class Seed {
	public int len;
	public int type;
	public Constraint overall;
	public Constraint zones[];
	
	public Seed() {
		init(0,0,null);
	}
	
	public Seed(int ln, int ty, Constraint oc) {
		init(ln, ty, oc);
	}
	
	public void init(int ln, int ty, Constraint oc) {
		len = ln;
		type = ty;
		overall = oc;
	}
	
	public Boolean instantiate(
			Read read,
			BTDnaString seq,
			BTString qual,
			Scoring pens,
			int depth,
			int seedoffidx,
			int seedtypeidx,
			Boolean fw,
			InstantiatedSeed is) {
		int seedlen = len;
		if((int)read.length() < seedlen) {
			// Shrink seed length to fit read if necessary
			seedlen = (int)read.length();
		}
		is.steps.resize(seedlen);
		is.zones.resize(seedlen);
		// Fill in 'steps' and 'zones'
		//
		// The 'steps' list indicates which read character should be
		// incorporated at each step of the search process.  Often we will
		// simply proceed from one end to the other, in which case the
		// 'steps' list is ascending or descending.  In some cases (e.g.
		// the 2mm case), we might want to switch directions at least once
		// during the search, in which case 'steps' will jump in the
		// middle.  When an element of the 'steps' list is negative, this
		// indicates that the next
		//
		// The 'zones' list indicates which zone constraint is active at
		// each step.  Each element of the 'zones' list is a pair; the
		// first pair element indicates the applicable zone when
		// considering either mismatch or delete (ref gap) events, while
		// the second pair element indicates the applicable zone when
		// considering insertion (read gap) events.  When either pair
		// element is a negative number, that indicates that we are about
		// to leave the zone for good, at which point we may need to
		// evaluate whether we have reached the zone's budget.
		//
		switch(type) {
			case SEED_TYPE_EXACT: {
				for(int k = 0; k < seedlen; k++) {
					is.steps[k] = -(seedlen - k);
					// Zone 0 all the way
					is.zones[k].first = is.zones[k].second = 0;
				}
				break;
			}
			case SEED_TYPE_LEFT_TO_RIGHT: {
				for(int k = 0; k < seedlen; k++) {
					is.steps[k] = k+1;
					// Zone 0 from 0 up to ceil(len/2), then 1
					is.zones[k].first = is.zones[k].second = ((k < (seedlen+1)/2) ? 0 : 1);
				}
				// Zone 1 ends at the RHS
				is.zones[seedlen-1].first = is.zones[seedlen-1].second = -1;
				break;
			}
			case SEED_TYPE_RIGHT_TO_LEFT: {
				for(int k = 0; k < seedlen; k++) {
					is.steps[k] = -(seedlen - k);
					// Zone 0 from 0 up to floor(len/2), then 1
					is.zones[k].first  = ((k < seedlen/2) ? 0 : 1);
					// Inserts: Zone 0 from 0 up to ceil(len/2)-1, then 1
					is.zones[k].second = ((k < (seedlen+1)/2+1) ? 0 : 1);
				}
				is.zones[seedlen-1].first = is.zones[seedlen-1].second = -1;
				break;
			}
			case SEED_TYPE_INSIDE_OUT: {
				// Zone 0 from ceil(N/4) up to N-floor(N/4)
				int step = 0;
				for(int k = (seedlen+3)/4; k < seedlen - (seedlen/4); k++) {
					is.zones[step].first = is.zones[step].second = 0;
					is.steps[step++] = k+1;
				}
				// Zone 1 from N-floor(N/4) up
				for(int k = seedlen - (seedlen/4); k < seedlen; k++) {
					is.zones[step].first = is.zones[step].second = 1;
					is.steps[step++] = k+1;
				}
				// No Zone 1 if seedlen is short (like 2)
				//assert_eq(1, is.zones[step-1].first);
				is.zones[step-1].first = is.zones[step-1].second = -1;
				// Zone 2 from ((seedlen+3)/4)-1 down to 0
				for(int k = ((seedlen+3)/4)-1; k >= 0; k--) {
					is.zones[step].first = is.zones[step].second = 2;
					is.steps[step++] = -(k+1);
				}
				is.zones[step-1].first = is.zones[step-1].second = -2;
				break;
			}
		}
		// Instantiate constraints
		for(int i = 0; i < 3; i++) {
			is.cons[i] = zones[i];
			is.cons[i].instantiate(read.length());
		}
		is.overall = overall;
		is.overall.instantiate(read.length());
		// Take a sweep through the seed sequence.  Consider where the Ns
		// occur and how zones are laid out.  Calculate the maximum number
		// of positions we can jump over initially (e.g. with the ftab) and
		// perhaps set this function's return value to false, indicating
		// that the arrangements of Ns prevents the seed from aligning.
		Boolean streak = true;
		is.maxjump = 0;
		Boolean ret = true;
		Boolean ltr = (is.steps[0] > 0); // true -> left-to-right
		for(double i = 0; i < is.steps.size(); i++) {
			int off = is.steps[i];
			off = abs(off)-1;
			Constraint cons = is.cons[abs(is.zones[i].first)];
			int c = seq[off];
			int q = qual[off];
			if(ltr != (is.steps[i] > 0) || // changed direction
			   is.zones[i].first != 0 ||   // changed zone
			   is.zones[i].second != 0)    // changed zone
			{
				streak = false;
			}
			if(c == 4) {
				// Induced mismatch
				if(cons.canN(q, pens)) {
					cons.chargeN(q, pens);
				} else {
					// Seed disqualified due to arrangement of Ns
					return false;
				}
			}
			if(streak) is.maxjump++;
		}
		is.seedoff = depth;
		is.seedoffidx = seedoffidx;
		is.fw = fw;
		is.s = this;
		return ret;
	}
	
	public void oneMmSeeds(int ln, EList<Seed> pols, Constraint oall) {
		oall.init();
		// Seed policy 1: left-to-right search
		pols.expand();
		pols.back().len = ln;
		pols.back().type = SEED_TYPE_LEFT_TO_RIGHT;
		pols.back().zones[0] = Constraint.exact();
		pols.back().zones[1] = Constraint.mmBased(1);
		pols.back().zones[2] = Constraint.exact(); // not used
		pols.back().overall = oall;
		// Seed policy 2: right-to-left search
		pols.expand();
		pols.back().len = ln;
		pols.back().type = SEED_TYPE_RIGHT_TO_LEFT;
		pols.back().zones[0] = Constraint.exact();
		pols.back().zones[1] = Constraint.mmBased(1);
		pols.back().zones[1].mmsCeil = 0;
		pols.back().zones[2] = Constraint.exact(); // not used
		pols.back().overall = oall;
	}
	
	public void twoMmSeeds(int ln, EList<Seed> pols, Constraint oall) {
		oall.init();
		// Seed policy 1: left-to-right search
		pols.expand();
		pols.back().len = ln;
		pols.back().type = SEED_TYPE_LEFT_TO_RIGHT;
		pols.back().zones[0] = Constraint.exact();
		pols.back().zones[1] = Constraint.mmBased(2);
		pols.back().zones[2] = Constraint.exact(); // not used
		pols.back().overall = &oall;
		// Seed policy 2: right-to-left search
		pols.expand();
		pols.back().len = ln;
		pols.back().type = SEED_TYPE_RIGHT_TO_LEFT;
		pols.back().zones[0] = Constraint.exact();
		pols.back().zones[1] = Constraint.mmBased(2);
		pols.back().zones[1].mmsCeil = 1; // Must have used at least 1 mismatch
		pols.back().zones[2] = Constraint.exact(); // not used
		pols.back().overall = oall;
		// Seed policy 3: inside-out search
		pols.expand();
		pols.back().len = ln;
		pols.back().type = SEED_TYPE_INSIDE_OUT;
		pols.back().zones[0] = Constraint.exact();
		pols.back().zones[1] = Constraint.mmBased(1);
		pols.back().zones[1].mmsCeil = 0; // Must have used at least 1 mismatch
		pols.back().zones[2] = Constraint.mmBased(1);
		pols.back().zones[2].mmsCeil = 0; // Must have used at least 1 mismatch
		pols.back().overall = oall;
	}
	
	public boolean acceptable() {
		return zones[0].acceptable() &&
			       zones[1].acceptable() &&
			       zones[2].acceptable() &&
			       overall.acceptable();
	}
	
	public static void mmSeeds(int mms,
			int ln,
			EList<Seed> pols,
			Constraint oall) {
		if(mms == 0) {
			zeroMmSeeds(ln, pols, oall);
		} else if(mms == 1) {
			oneMmSeeds(ln, pols, oall);
		} else if(mms == 2) {
			twoMmSeeds(ln, pols, oall);
		}
	}
}
