/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

/* //ifndef ALIGNER_SEED_H_
//define ALIGNER_SEED_H_

//include <iostream>
//include <utility>
//include <limits>
//include "qual.h"
//include "ds.h"
//include "sstring.h"
//include "alphabet.h"
//include "edit.h"
//include "read.h"
// Threading is necessary to synchronize the classes that dump
// intermediate alignment results to files.  Otherwise, all data herein
// is ant and shared, or per-thread.
//include "threading.h"
//include "aligner_result.h"
//include "aligner_cache.h"
//include "scoring.h"
//include "mem_ids.h"
//include "simple_func.h"
//include "btypes.h" */

package com.uwb.bt2j.aligner;

/**
 * A raint to apply to an alignment zone, or to an overall
 * alignment.
 *
 * The raint can put both caps and ceilings on the number and
 * types of edits allowed.
 */
class Constraint {
	
	public int edits;      // // edits permitted
	public int mms;        // // mismatches permitted
	public int ins;        // // insertions permitted
	public int dels;       // // deletions permitted
	public int penalty;    // penalty total permitted
	public int editsCeil;  // <= this many edits can be left at the end
	public int mmsCeil;    // <= this many mismatches can be left at the end
	public int insCeil;    // <= this many inserts can be left at the end
	public int delsCeil;   // <= this many deletions can be left at the end
	public int penaltyCeil;// <= this much leftover penalty can be left at the end
	public SimpleFunc penFunc;// penalty function; function of read len
	public Boolean instantiated; // whether raint is instantiated w/r/t read len
	
	public Constraint() { init(); }
	
	/**
	 * Initialize Constraint to be fully permissive.
	 */
	public void init() {
		edits = mms = ins = dels = penalty = editsCeil = mmsCeil =
		insCeil = delsCeil = penaltyCeil = MAX_I;
		penFunc.reset();
		instantiated = false;
	}
	
	/**
	 * Return true iff penalities and raints prevent us from
	 * adding any edits.
	 */
	public Boolean mustMatch() {
		assert(instantiated);
		return (mms == 0 && edits == 0) ||
		        penalty == 0 ||
		       (mms == 0 && dels == 0 && ins == 0);
	}
	
	/**
	 * Return true iff a mismatch of the given quality is permitted.
	 */
	public Boolean canMismatch(int q,  Scoring cm) {
		assert(instantiated);
		return (mms > 0 || edits > 0) &&
		       penalty >= cm.mm(q);
	}

	/**
	 * Return true iff a mismatch of the given quality is permitted.
	 */
	public Boolean canN(int q,  Scoring cm) {
		assert(instantiated);
		return (mms > 0 || edits > 0) &&
		       penalty >= cm.n(q);
	}
	
	/**
	 * Return true iff a mismatch of *any* quality (even qual=1) is
	 * permitted.
	 */
	public Boolean canMismatch() {
		assert(instantiated);
		return (mms > 0 || edits > 0) && penalty > 0;
	}

	/**
	 * Return true iff a mismatch of *any* quality (even qual=1) is
	 * permitted.
	 */
	public Boolean canN() {
		assert(instantiated);
		return (mms > 0 || edits > 0);
	}
	
	/**
	 * Return true iff a deletion of the given extension (0=open, 1=1st
	 * extension, etc) is permitted.
	 */
	public Boolean canDelete(int ex,  Scoring cm) {
		assert(instantiated);
		return (dels > 0 && edits > 0) &&
		       penalty >= cm.del(ex);
	}

	/**
	 * Return true iff a deletion of any extension is permitted.
	 */
	public Boolean canDelete() {
		assert(instantiated);
		return (dels > 0 || edits > 0) &&
		       penalty > 0;
	}
	
	/**
	 * Return true iff an insertion of the given extension (0=open,
	 * 1=1st extension, etc) is permitted.
	 */
	public Boolean canInsert(int ex,  Scoring cm) {
		assert(instantiated);
		return (ins > 0 || edits > 0) &&
		       penalty >= cm.ins(ex);
	}

	/**
	 * Return true iff an insertion of any extension is permitted.
	 */
	public Boolean canInsert() {
		assert(instantiated);
		return (ins > 0 || edits > 0) &&
		       penalty > 0;
	}
	
	/**
	 * Return true iff a gap of any extension is permitted
	 */
	public Boolean canGap() {
		assert(instantiated);
		return ((ins > 0 || dels > 0) || edits > 0) && penalty > 0;
	}
	
	/**
	 * Charge a mismatch of the given quality.
	 */
	public void chargeMismatch(int q,  Scoring cm) {
		assert(instantiated);
		if(mms == 0) { ////assert_gt(edits, 0);
		edits--; }
		else mms--;
		penalty -= cm.mm(q);
		//assert_geq(mms, 0);
		//assert_geq(edits, 0);
		//assert_geq(penalty, 0);
	}
	
	/**
	 * Charge an N mismatch of the given quality.
	 */
	public void chargeN(int q,  Scoring cm) {
		assert(instantiated);
		if(mms == 0) { //assert_gt(edits, 0); edits--; }
		else mms--;
		penalty -= cm.n(q);
		//assert_geq(mms, 0);
		//assert_geq(edits, 0);
		//assert_geq(penalty, 0);
	}
	
	/**
	 * Charge a deletion of the given extension.
	 */
	public void chargeDelete(int ex,  Scoring cm) {
		assert(instantiated);
		dels--;
		edits--;
		penalty -= cm.del(ex);
		//assert_geq(dels, 0);
		//assert_geq(edits, 0);
		//assert_geq(penalty, 0);
	}

	/**
	 * Charge an insertion of the given extension.
	 */
	public void chargeInsert(int ex,  Scoring cm) {
		assert(instantiated);
		ins--;
		edits--;
		penalty -= cm.ins(ex);
		//assert_geq(ins, 0);
		//assert_geq(edits, 0);
		//assert_geq(penalty, 0);
	}
	
	/**
	 * Once the rained area is completely explored, call this
	 * function to check whether there were *at least* as many
	 * dissimilarities as required by the raint.  Bounds like this
	 * are helpful to resolve instances where two search roots would
	 * otherwise overlap in what alignments they can find.
	 */
	public Boolean acceptable() {
		assert(instantiated);
		return edits   <= editsCeil &&
		       mms     <= mmsCeil   &&
		       ins     <= insCeil   &&
		       dels    <= delsCeil  &&
		       penalty <= penaltyCeil;
	}
	
	/**
	 * Instantiate a raint w/r/t the read length and the ant
	 * and linear coefficients for the penalty function.
	 */
	public static int instantiate(int rdlen,  SimpleFunc func) {
		return func.f<int>((double)rdlen);
	}
	
	/**
	 * Instantiate this raint w/r/t the read length.
	 */
	public void instantiate(int rdlen) {
		assert(!instantiated);
		if(penFunc.initialized()) {
			penalty = Constraint::instantiate(rdlen, penFunc);
		}
		instantiated = true;
	}
	
	//
	// Some static methods for ructing some standard Constraints
	//

	/**
	 * Construct a raint with no edits of any kind allowed.
	 */
	public static Constraint exact() {
		Constraint c;
		c.edits = c.mms = c.ins = c.dels = c.penalty = 0;
		return c;
	}
	
	/**
	 * Construct a raint where the only raint is a total
	 * penalty raint.
	 */
	public static Constraint penaltyBased(int pen) {
		Constraint c;
		c.penalty = pen;
		return c;
	}

	/**
	 * Construct a raint where the only raint is a total
	 * penalty raint related to the length of the read.
	 */
	public static Constraint penaltyFuncBased( SimpleFunc func) {
		Constraint c;
		c.penFunc = f;
		return c;
	}

	/**
	 * Construct a raint where the only raint is a total
	 * penalty raint.
	 */
	public static Constraint mmBased(int mms) {
		Constraint c;
		c.mms = mms;
		c.edits = c.dels = c.ins = 0;
		return c;
	}

	/**
	 * Construct a raint where the only raint is a total
	 * penalty raint.
	 */
	public static Constraint editBased(int edits) {
		Constraint c;
		c.edits = edits;
		c.dels = c.ins = c.mms = 0;
		return c;
	}
}

/**
 * We divide seed search strategies into three categories:
 *
 * 1. A left-to-right search where the left half of the read is
 *    rained to match exactly and the right half is subject to
 *    some looser raint (e.g. 1mm or 2mm)
 * 2. Same as 1, but going right to left with the exact matching half
 *    on the right.
 * 3. Inside-out search where the center half of the read is
 *    rained to match exactly, and the extreme quarters of the
 *    read are subject to a looser raint.
 */
public enum type{
	SEED_TYPE_EXACT(1),
	SEED_TYPE_LEFT_TO_RIGHT(2),
	SEED_TYPE_RIGHT_TO_LEFT(3),
	SEED_TYPE_INSIDE_OUT(4);
	
	private int final value;
	type(int v) {
		value = v;
	}
};

class InstantiatedSeed;

/**
 * Policy dictating how to size and arrange seeds along the length of
 * the read, and what raints to force on the zones of the seed.
 * We assume that seeds are plopped down at regular intervals from the
 * 5' to 3' ends, with the first seed flush to the 5' end.
 *
 * If the read is shorter than a single seed, one seed is used and it
 * is shrunk to accommodate the read.
 */
class Seed {

	public int len;             // length of a seed
	public int type;            // dictates anchor portion, direction of search
	public Constraint *overall; // for the overall alignment

	public Seed() { init(0, 0, null); }

	/**
	 * Construct and initialize this seed with given length and type.
	 */
	public Seed(int ln, int ty, Constraint oc) {
		init(ln, ty, oc);
	}

	/**
	 * Initialize this seed with given length and type.
	 */
	public void init(int ln, int ty, Constraint oc) {
		len = ln;
		type = ty;
		overall = oc;
	}
	
	// If the seed is split into halves, we just use zones[0] and
	// zones[1]; 0 is the near half and 1 is the far half.  If the seed
	// is split into thirds (i.e. inside-out) then 0 is the center, 1
	// is the far portion on the left, and 2 is the far portion on the
	// right.
	public Constraint zones[3];

	/**
	 * Once the rained seed is completely explored, call this
	 * function to check whether there were *at least* as many
	 * dissimilarities as required by all raints.  Bounds like this
	 * are helpful to resolve instances where two search roots would
	 * otherwise overlap in what alignments they can find.
	 */
	public Boolean acceptable() {
		assert(overall != null);
		return zones[0].acceptable() &&
		       zones[1].acceptable() &&
		       zones[2].acceptable() &&
		       overall.acceptable();
	}

	/**
	 * Given a read, depth and orientation, extract a seed data structure
	 * from the read and fill in the steps & zones arrays.  The Seed
	 * contains the sequence and quality values.
	 */
	public Boolean instantiate(
		 Read read,
		 BTDnaString seq, // already-extracted seed sequence
		 BTString qual,   // already-extracted seed quality sequence
		 Scoring pens,
		int depth,
		int seedoffidx,
		int seedtypeidx,
		Boolean fw,
		InstantiatedSeed si) {
			assert(overall != null);
			int seedlen = len;
			if((int)read.length() < seedlen) {
				// Shrink seed length to fit read if necessary
				seedlen = (int)read.length();
			}
			////assert_gt(seedlen, 0);
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
			// The 'zones' list indicates which zone raint is active at
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
					////assert_eq(1, is.zones[step-1].first);
					is.zones[step-1].first = is.zones[step-1].second = -1;
					// Zone 2 from ((seedlen+3)/4)-1 down to 0
					for(int k = ((seedlen+3)/4)-1; k >= 0; k--) {
						is.zones[step].first = is.zones[step].second = 2;
						is.steps[step++] = -(k+1);
					}
					////assert_eq(2, is.zones[step-1].first);
					is.zones[step-1].first = is.zones[step-1].second = -2;
					////assert_eq(seedlen, step);
					break;
				}
				default:
					throw 1;
			}
					// Instantiate raints
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
			Boolean ltr = (is.steps[0] > 0); // true . left-to-right
			for(int i = 0; i < is.steps.size(); i++) {
				////assert_neq(0, is.steps[i]);
				int off = is.steps[i];
				off = abs(off)-1;
				Constraint cons = is.cons[abs(is.zones[i].first)];
				int c = seq[off];  //assert_range(0, 4, c);
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

	/**
	 * Return a list of Seed objects encapsulating
	 */
	public static void mmSeeds(
		int mms,
		int ln,
		EList<Seed>& pols,
		Constraint& oall)
	{
		if(mms == 0) {
			zeroMmSeeds(ln, pols, oall);
		} else if(mms == 1) {
			oneMmSeeds(ln, pols, oall);
		} else if(mms == 2) {
			twoMmSeeds(ln, pols, oall);
		} else throw 1;
	}
	
	public static void zeroMmSeeds(int ln, EList<Seed>, Constraint) {
		oall.init();
		// Seed policy 1: left-to-right search
		pols.expand();
		pols.back().len = ln;
		pols.back().type = SEED_TYPE_EXACT;
		pols.back().zones[0] = Constraint::exact();
		pols.back().zones[1] = Constraint::exact();
		pols.back().zones[2] = Constraint::exact();
		pols.back().overall = oall;
	}
	
	public static void oneMmSeeds (int ln, EList<Seed>, Constraint) {
		oall.init();
		// Seed policy 1: left-to-right search
		pols.expand();
		pols.back().len = ln;
		pols.back().type = SEED_TYPE_LEFT_TO_RIGHT;
		pols.back().zones[0] = Constraint::exact();
		pols.back().zones[1] = Constraint::mmBased(1);
		pols.back().zones[2] = Constraint::exact(); // not used
		pols.back().overall = oall;
		// Seed policy 2: right-to-left search
		pols.expand();
		pols.back().len = ln;
		pols.back().type = SEED_TYPE_RIGHT_TO_LEFT;
		pols.back().zones[0] = Constraint::exact();
		pols.back().zones[1] = Constraint::mmBased(1);
		pols.back().zones[1].mmsCeil = 0;
		pols.back().zones[2] = Constraint::exact(); // not used
		pols.back().overall = oall;
	}
	
	public static void twoMmSeeds (int ln, EList<Seed>, Constraint) {
		oall.init();
		// Seed policy 1: left-to-right search
		pols.expand();
		pols.back().len = ln;
		pols.back().type = SEED_TYPE_LEFT_TO_RIGHT;
		pols.back().zones[0] = Constraint::exact();
		pols.back().zones[1] = Constraint::mmBased(2);
		pols.back().zones[2] = Constraint::exact(); // not used
		pols.back().overall = oall;
		// Seed policy 2: right-to-left search
		pols.expand();
		pols.back().len = ln;
		pols.back().type = SEED_TYPE_RIGHT_TO_LEFT;
		pols.back().zones[0] = Constraint::exact();
		pols.back().zones[1] = Constraint::mmBased(2);
		pols.back().zones[1].mmsCeil = 1; // Must have used at least 1 mismatch
		pols.back().zones[2] = Constraint::exact(); // not used
		pols.back().overall = oall;
		// Seed policy 3: inside-out search
		pols.expand();
		pols.back().len = ln;
		pols.back().type = SEED_TYPE_INSIDE_OUT;
		pols.back().zones[0] = Constraint::exact();
		pols.back().zones[1] = Constraint::mmBased(1);
		pols.back().zones[1].mmsCeil = 0; // Must have used at least 1 mismatch
		pols.back().zones[2] = Constraint::mmBased(1);
		pols.back().zones[2].mmsCeil = 0; // Must have used at least 1 mismatch
		pols.back().overall = oall;
	}
}

/**
 * Types of actions that can be taken by the SeedAligner.
 */
public enum actions{
	SA_ACTION_TYPE_RESET(1),
	SA_ACTION_TYPE_SEARCH_SEED(2), // 2
	SA_ACTION_TYPE_FTAB(3),        // 3
	SA_ACTION_TYPE_FCHR(4),        // 4
	SA_ACTION_TYPE_MATCH(5),       // 5
	SA_ACTION_TYPE_EDIT(6);         // 6
	
	private int final value;
	actions(int v) {
		value = v;
	}
};

/**
 * An instantiated seed is a seed (perhaps modified to fit the read)
 * plus all data needed to conduct a search of the seed.
 */
class InstantiatedSeed {

	public InstantiatedSeed() : steps(AL_CAT), zones(AL_CAT) { }

	// Steps map.  There are as many steps as there are positions in
	// the seed.  The map is a helpful abstraction because we sometimes
	// visit seed positions in an irregular order (e.g. inside-out
	// search).
	public EList<int> steps;

	// Zones map.  For each step, records what raint to charge an
	// edit to.  The first entry in each pair gives the raint for
	// non-insert edits and the second entry in each pair gives the
	// raint for insert edits.  If the value stored is negative,
	// this indicates that the zone is "closed out" after this
	// position, so zone acceptility should be checked.
	public EList<pair<int, int> > zones;

	// Nucleotide sequence covering the seed, extracted from read
	public BTDnaString seq;
	
	// Quality sequence covering the seed, extracted from read
	public BTString qual;
	
	// Initial raints governing zones 0, 1, 2.  We precalculate
	// the effect of Ns on these.
	public Constraint cons[3];
	
	// Overall raint, tailored to the read length.
	public Constraint overall;
	
	// Maximum number of positions that the aligner may advance before
	// its first step.  This lets the aligner know whether it can use
	// the ftab or not.
	public int maxjump;
	
	// Offset of seed from 5' end of read
	public int seedoff;

	// Id for seed offset; ids are such that the smallest index is the
	// closest to the 5' end and consecutive ids are adjacent (i.e.
	// there are no intervening offsets with seeds)
	public int seedoffidx;
	
	// Type of seed (left-to-right, etc)
	public int seedtypeidx;
	
	// Seed comes from forward-oriented read?
	public Boolean fw;
	
	// Filtered out due to the pattern of Ns present.  If true, this
	// seed should be ignored by searchAllSeeds().
	public Boolean nfiltered;
	
	// Seed this was instantiated from
	public Seed s;
	
//ifndef NDEBUG
	/**
	 * Check that InstantiatedSeed is internally consistent.
	 */
	public Boolean repOk()  {
		return true;
	}
//endif
}

/**
 * Simple struct for holding a end-to-end alignments for the read with at most
 * 2 edits.
 */
class EEHit {
	
	public TIndexOffU top;
	public TIndexOffU bot;
	public Edit     e1;
	public Edit     e2;
	public Boolean     fw;
	public long  score;
	public EEHit() { reset(); }
	
	public void reset() {
		top = bot = 0;
		fw = false;
		e1.reset();
		e2.reset();
		score = MIN_I64;
	}
	
	public void init(
		TIndexOffU top_,
		TIndexOffU bot_,
		 Edit e1_,
		 Edit e2_,
		Boolean fw_,
		long score_)
	{
		top = top_; bot = bot_;
		if(e1_ != null) {
			e1 = e1_;
		} else {
			e1.reset();
		}
		if(e2_ != null) {
			e2 = e2_;
		} else {
			e2.reset();
		}
		fw = fw_;
		score = score_;
	}
	
	/**
	 * Return number of mismatches in the alignment.
	 */
	public int mms()  {
		if     (e2.inited()) return 2;
		else if(e1.inited()) return 1;
		else                 return 0;
	}
	
	/**
	 * Return the number of Ns involved in the alignment.
	 */
	public int ns()  {
		int ns = 0;
		if(e1.inited() && e1.hasN()) {
			ns++;
			if(e2.inited() && e2.hasN()) {
				ns++;
			}
		}
		return ns;
	}

	/**
	 * Return the number of Ns involved in the alignment.
	 */
	public int refns()  {
		int ns = 0;
		if(e1.inited() && e1.chr == 'N') {
			ns++;
			if(e2.inited() && e2.chr == 'N') {
				ns++;
			}
		}
		return ns;
	}
	
	/**
	 * Return true iff there is no hit.
	 */
	public Boolean empty()  {
		return bot <= top;
	}
	
	/**
	 * Higher score = higher priority.
	 */
	public Boolean LessThan( EEHit o)  {
		return score > o.score;
	}
	
	/**
	 * Return the size of the alignments SA range.s
	 */
	public TIndexOffU size()  { return bot - top; }
	
//ifndef NDEBUG
	/**
	 * Check that hit is sane w/r/t read.
	 */
	public Boolean repOk( Read rd)  {
		//assert_gt(bot, top);
		if(e1.inited()) {
			//assert_lt(e1.pos, rd.length());
			if(e2.inited()) {
				//assert_lt(e2.pos, rd.length());
			}
		}
		return true;
	}
//endif
}

/**
 * Data structure for holding all of the seed hits associated with a read.  All
 * the seed hits for a given read are encapsulated in a single QVal object.  A
 * QVal refers to a range of values in the qlist, where each qlist value is a 
 * BW range and a slot to hold the hit's suffix array offset.  QVals are kept
 * in two lists (hitsFw_ and hitsRc_), one for seeds on the forward read strand,
 * one for seeds on the reverse read strand.  The list is indexed by read
 * offset index (e.g. 0=closest-to-5', 1=second-closest, etc).
 *
 * An assumption behind this data structure is that all the seeds are found
 * first, then downstream analyses try to extend them.  In between finding the
 * seed hits and extending them, the sort() member function is called, which
 * ranks QVals according to the order they should be extended.  Right now the
 * policy is that QVals with fewer elements (hits) should be tried first.
 */
class SeedResults {
	
	// As seed hits and edits are added they're sorted into these
	// containers
	protected EList<BTDnaString>  seqFw_;       // seqs for seeds from forward read
	protected EList<BTDnaString>  seqRc_;       // seqs for seeds from revcomp read
	protected EList<BTString>     qualFw_;      // quals for seeds from forward read
	protected EList<BTString>     qualRc_;      // quals for seeds from revcomp read
	protected EList<QVal>         hitsFw_;      // hits for forward read
	protected EList<QVal>         hitsRc_;      // hits for revcomp read
	protected EList<EList<InstantiatedSeed> > isFw_; // hits for forward read
	protected EList<EList<InstantiatedSeed> > isRc_; // hits for revcomp read
	protected EList<Boolean>         sortedFw_;    // true iff fw QVal was sorted/ranked
	protected EList<Boolean>         sortedRc_;    // true iff rc QVal was sorted/ranked
	protected int              nonzTot_;     // // offsets with non-zero size
	protected int              uniTot_;      // // offsets unique hit
	protected int              uniTotS_[2];  // // offsets unique hit on each strand
	protected int              repTot_;      // // offsets repetitive hit
	protected int              repTotS_[2];  // // offsets repetitive hit on each strand
	protected int              nonzFw_;      // // offsets into fw read with non-0 size
	protected int              nonzRc_;      // // offsets into rc read with non-0 size
	protected int              numRanges_;   // // ranges added
	protected int              numElts_;     // // elements added
	protected int              numRangesFw_; // // ranges added for fw seeds
	protected int              numEltsFw_;   // // elements added for fw seeds
	protected int              numRangesRc_; // // ranges added for rc seeds
	protected int              numEltsRc_;   // // elements added for rc seeds

	protected EList<int>     offIdx2off_;// map from offset indexes to offsets from 5' end

	// When the sort routine is called, the seed hits collected so far
	// are sorted into another set of containers that allow easy access
	// to hits from the lowest-ranked offset (the one with the fewest
	// BW elements) to the greatest-ranked offset.  Offsets with 0 hits
	// are ignored.
	protected EList<int>     rankOffs_;  // sorted offests of seeds to try
	protected EList<Boolean>         rankFws_;   // sorted orientations assoc. with rankOffs_
	protected Boolean                sorted_;    // true if sort() called since last reset
	
	// These fields set once per read
	protected int              numOffs_;   // // different seed offsets possible
	protected  Read         read_;      // read from which seeds were extracted
	
	protected EEHit               exactFwHit_; // end-to-end exact hit for fw read
	protected EEHit               exactRcHit_; // end-to-end exact hit for rc read
	protected EList<EEHit>        mm1Hit_;     // 1-mismatch end-to-end hits
	protected int              mm1Elt_;     // number of 1-mismatch hit rows
	protected Boolean                mm1Sorted_;  // true iff we've sorted the mm1Hit_ list

	protected EList<int> tmpMedian_; // temporary storage for calculating median

	public SeedResults() :
		seqFw_(AL_CAT),
		seqRc_(AL_CAT),
		qualFw_(AL_CAT),
		qualRc_(AL_CAT),
		hitsFw_(AL_CAT),
		hitsRc_(AL_CAT),
		isFw_(AL_CAT),
		isRc_(AL_CAT),
		sortedFw_(AL_CAT),
		sortedRc_(AL_CAT),
		offIdx2off_(AL_CAT),
		rankOffs_(AL_CAT),
		rankFws_(AL_CAT),
		mm1Hit_(AL_CAT)
	{
		clear();
	}
	
	/**
	 * Set the current read.
	 */
	public void nextRead( Read read) {
		read_ = read;
	}

	/**
	 * Set the appropriate element of either hitsFw_ or hitsRc_ to the given
	 * QVal.  A QVal encapsulates all the BW ranges for reference substrings 
	 * that are within some distance of the seed string.
	 */
	public void add(
		 QVal qv,           // range of ranges in cache
		 AlignmentCache ac, // cache
		int seedIdx,         // seed index (from 5' end)
		Boolean     seedFw)          // whether seed is from forward read
	{
		assert(qv.repOk(ac));
		assert(repOk(ac));
		//assert_lt(seedIdx, hitsFw_.size());
		//assert_gt(numOffs_, 0); // if this fails, probably failed to call reset
		if(qv.empty()) return;
		if(seedFw) {
			assert(!hitsFw_[seedIdx].valid());
			hitsFw_[seedIdx] = qv;
			numEltsFw_ += qv.numElts();
			numRangesFw_ += qv.numRanges();
			if(qv.numRanges() > 0) nonzFw_++;
		} else {
			assert(!hitsRc_[seedIdx].valid());
			hitsRc_[seedIdx] = qv;
			numEltsRc_ += qv.numElts();
			numRangesRc_ += qv.numRanges();
			if(qv.numRanges() > 0) nonzRc_++;
		}
		numElts_ += qv.numElts();
		numRanges_ += qv.numRanges();
		if(qv.numRanges() > 0) {
			nonzTot_++;
			if(qv.numRanges() == 1 && qv.numElts() == 1) {
				uniTot_++;
				uniTotS_[seedFw ? 0 : 1]++;
			} else {
				repTot_++;
				repTotS_[seedFw ? 0 : 1]++;
			}
		}
		assert(repOk(ac));
	}

	/**
	 * Clear buffered seed hits and state.  Set the number of seed
	 * offsets and the read.
	 */
	public void reset(
		 Read read,
		 EList<int> offIdx2off,
		int numOffs)
	{
		//assert_gt(numOffs, 0);
		clearSeeds();
		numOffs_ = numOffs;
		seqFw_.resize(numOffs_);
		seqRc_.resize(numOffs_);
		qualFw_.resize(numOffs_);
		qualRc_.resize(numOffs_);
		hitsFw_.resize(numOffs_);
		hitsRc_.resize(numOffs_);
		isFw_.resize(numOffs_);
		isRc_.resize(numOffs_);
		sortedFw_.resize(numOffs_);
		sortedRc_.resize(numOffs_);
		offIdx2off_ = offIdx2off;
		for(int i = 0; i < numOffs_; i++) {
			sortedFw_[i] = sortedRc_[i] = false;
			hitsFw_[i].reset();
			hitsRc_[i].reset();
			isFw_[i].clear();
			isRc_[i].clear();
		}
		read_ = &read;
		sorted_ = false;
	}
	
	/**
	 * Clear seed-hit state.
	 */
	public void clearSeeds() {
		sortedFw_.clear();
		sortedRc_.clear();
		rankOffs_.clear();
		rankFws_.clear();
		offIdx2off_.clear();
		hitsFw_.clear();
		hitsRc_.clear();
		isFw_.clear();
		isRc_.clear();
		seqFw_.clear();
		seqRc_.clear();
		nonzTot_ = 0;
		uniTot_ = uniTotS_[0] = uniTotS_[1] = 0;
		repTot_ = repTotS_[0] = repTotS_[1] = 0;
		nonzFw_ = 0;
		nonzRc_ = 0;
		numOffs_ = 0;
		numRanges_ = 0;
		numElts_ = 0;
		numRangesFw_ = 0;
		numEltsFw_ = 0;
		numRangesRc_ = 0;
		numEltsRc_ = 0;
	}
	
	/**
	 * Clear seed-hit state and end-to-end alignment state.
	 */
	public void clear() {
		clearSeeds();
		read_ = null;
		exactFwHit_.reset();
		exactRcHit_.reset();
		mm1Hit_.clear();
		mm1Sorted_ = false;
		mm1Elt_ = 0;
		assert(empty());
	}
	
	/**
	 * Extract key summaries from this SeedResults and put into 'ssum'.
	 */
	public void toSeedAlSumm(SeedAlSumm ssum)  {
		// Number of positions with at least 1 range
		ssum.nonzTot   = nonzTot_;
		ssum.nonzFw    = nonzFw_;
		ssum.nonzRc    = nonzRc_;

		// Number of ranges
		ssum.nrangeTot = numRanges_;
		ssum.nrangeFw  = numRangesFw_;
		ssum.nrangeRc  = numRangesRc_;

		// Number of elements
		ssum.neltTot   = numElts_;
		ssum.neltFw    = numEltsFw_;
		ssum.neltRc    = numEltsRc_;
		
		// Other summaries
		ssum.maxNonzRangeFw = ssum.minNonzRangeFw = 0;
		ssum.maxNonzRangeRc = ssum.minNonzRangeRc = 0;
		ssum.maxNonzEltFw = ssum.minNonzEltFw = 0;
		ssum.maxNonzEltRc = ssum.minNonzEltRc = 0;
		for(int i = 0; i < numOffs_; i++) {
			if(hitsFw_[i].valid()) {
				if(ssum.minNonzEltFw == 0 || hitsFw_[i].numElts() < ssum.minNonzEltFw) {
					ssum.minNonzEltFw = hitsFw_[i].numElts();
				}
				if(ssum.maxNonzEltFw == 0 || hitsFw_[i].numElts() > ssum.maxNonzEltFw) {
					ssum.maxNonzEltFw = hitsFw_[i].numElts();
				}
				if(ssum.minNonzRangeFw == 0 || hitsFw_[i].numRanges() < ssum.minNonzRangeFw) {
					ssum.minNonzRangeFw = hitsFw_[i].numRanges();
				}
				if(ssum.maxNonzRangeFw == 0 || hitsFw_[i].numRanges() > ssum.maxNonzRangeFw) {
					ssum.maxNonzRangeFw = hitsFw_[i].numRanges();
				}
			}
			if(hitsRc_[i].valid()) {
				if(ssum.minNonzEltRc == 0 || hitsRc_[i].numElts() < ssum.minNonzEltRc) {
					ssum.minNonzEltRc = hitsRc_[i].numElts();
				}
				if(ssum.maxNonzEltRc == 0 || hitsRc_[i].numElts() > ssum.maxNonzEltRc) {
					ssum.maxNonzEltRc = hitsRc_[i].numElts();
				}
				if(ssum.minNonzRangeRc == 0 || hitsRc_[i].numRanges() < ssum.minNonzRangeRc) {
					ssum.minNonzRangeRc = hitsRc_[i].numRanges();
				}
				if(ssum.maxNonzRangeRc == 0 || hitsRc_[i].numRanges() > ssum.maxNonzRangeRc) {
					ssum.maxNonzRangeRc = hitsRc_[i].numRanges();
				}
			}
		}
	}
	
	/**
	 * Return average number of hits per seed.
	 */
	public float averageHitsPerSeed()  {
		return nonzTot_ == 0 ? 0 : (float)numElts_ / (float)nonzTot_;
	}

	/**
	 * Return fraction of seeds that align uniquely.
	 */
	public int numUniqueSeeds()  {
		return uniTot_;
	}

	/**
	 * Return fraction of seeds that aligned uniquely on the given strand.
	 */
	public int numUniqueSeedsStrand(Boolean fw)  {
		return uniTotS_[fw ? 0 : 1];
	}

	/**
	 * Return fraction of seeds that align repetitively on the given strand.
	 */
	public int numRepeatSeeds()  {
		return repTot_;
	}
	
	/**
	 * Return fraction of seeds that align repetitively.
	 */
	public int numRepeatSeedsStrand(Boolean fw)  {
		return repTotS_[fw ? 0 : 1];
	}
	
	/**
	 * Return median of all the non-zero per-seed // hits
	 */
	public float medianHitsPerSeed()  {
		EList<int>& median = _cast<EList<int>>(tmpMedian_);
		median.clear();
		for(int i = 0; i < numOffs_; i++) {
			if(hitsFw_[i].valid() && hitsFw_[i].numElts() > 0) {
				median.push_back(hitsFw_[i].numElts());
			}
			if(hitsRc_[i].valid() && hitsRc_[i].numElts() > 0) {
				median.push_back(hitsRc_[i].numElts());
			}
		}
		if(tmpMedian_.empty()) {
			return 0.0f;
		}
		median.sort();
		float med1 = (float)median[tmpMedian_.size() >> 1];
		float med2 = med1;
		if((median.size() & 1) == 0) {
			med2 = (float)median[(tmpMedian_.size() >> 1) - 1];
		}
		return med1 + med2 * 0.5f;
	}
	
	/**
	 * Return a number that's meant to quantify how hopeful we are that this
	 * set of seed hits will lead to good alignments.
	 */
	public double uniquenessFactor()  {
		double result = 0.0;
		for(int i = 0; i < numOffs_; i++) {
			if(hitsFw_[i].valid()) {
				int nelt = hitsFw_[i].numElts();
				result += (1.0 / (double)(nelt * nelt));
			}
			if(hitsRc_[i].valid()) {
				int nelt = hitsRc_[i].numElts();
				result += (1.0 / (double)(nelt * nelt));
			}
		}
		return result;
	}

	/**
	 * Return the number of ranges being held.
	 */
	public int numRanges()  { return numRanges_; }

	/**
	 * Return the number of elements being held.
	 */
	public int numElts()  { return numElts_; }

	/**
	 * Return the number of ranges being held for seeds on the forward
	 * read strand.
	 */
	public int numRangesFw()  { return numRangesFw_; }

	/**
	 * Return the number of elements being held for seeds on the
	 * forward read strand.
	 */
	public int numEltsFw()  { return numEltsFw_; }

	/**
	 * Return the number of ranges being held for seeds on the
	 * reverse-complement read strand.
	 */
	public int numRangesRc()  { return numRangesRc_; }

	/**
	 * Return the number of elements being held for seeds on the
	 * reverse-complement read strand.
	 */
	public int numEltsRc()  { return numEltsRc_; }
	
	/**
	 * Given an offset index, return the offset that has that index.
	 */
	public int idx2off(int off)  {
		return offIdx2off_[off];
	}
	
	/**
	 * Return true iff there are 0 hits being held.
	 */
	public Boolean empty()  { return numRanges() == 0; }
	
	/**
	 * Get the QVal representing all the reference hits for the given
	 * orientation and seed offset index.
	 */
	 public QVal hitsAtOffIdx(Boolean fw, int seedoffidx)  {
		////assert_lt(seedoffidx, numOffs_);
		assert(repOk(null));
		return fw ? hitsFw_[seedoffidx] : hitsRc_[seedoffidx];
	}

	/**
	 * Get the Instantiated seeds for the given orientation and offset.
	 */
	public EList<InstantiatedSeed> instantiatedSeeds(Boolean fw, int seedoffidx) {
		////assert_lt(seedoffidx, numOffs_);
		assert(repOk(null));
		return fw ? isFw_[seedoffidx] : isRc_[seedoffidx];
	}
	
	/**
	 * Return the number of different seed offsets possible.
	 */
	public int numOffs()  { return numOffs_; }
	
	/**
	 * Return the read from which seeds were extracted, aligned.
	 */
	public Read read()  { return read_; }
	
//ifndef NDEBUG
	/**
	 * Check that this SeedResults is internally consistent.
	 */
	public Boolean repOk(
		 AlignmentCache ac,
		Boolean requireInited = false) 
	{
		if(requireInited) {
			assert(read_ != null);
		}
		if(numOffs_ > 0) {
			////assert_eq(numOffs_, hitsFw_.size());
			////assert_eq(numOffs_, hitsRc_.size());
			////assert_leq(numRanges_, numElts_);
			////assert_leq(nonzTot_, numRanges_);
			int nonzs = 0;
			for(int fw = 0; fw <= 1; fw++) {
				 EList<QVal> rrs = (fw ? hitsFw_ : hitsRc_);
				for(int i = 0; i < numOffs_; i++) {
					if(rrs[i].valid()) {
						if(rrs[i].numRanges() > 0) nonzs++;
						if(ac != null) {
							assert(rrs[i].repOk(ac));
						}
					}
				}
			}
			////assert_eq(nonzs, nonzTot_);
			assert(!sorted_ || nonzTot_ == rankFws_.size());
			assert(!sorted_ || nonzTot_ == rankOffs_.size());
		}
		return true;
	}
//endif
	
	/**
	 * Populate rankOffs_ and rankFws_ with the list of QVals that need to be
	 * examined for this SeedResults, in order.  The order is ascending by
	 * number of elements, so QVals with fewer elements (i.e. seed sequences
	 * that are more unique) will be tried first and QVals with more elements
	 * (i.e. seed sequences
	 */
	public void rankSeedHits(RandomSource rnd, Boolean all) {
		if(all) {
			for(int i = 1; i < numOffs_; i++) {
				for(int fwi = 0; fwi <= 1; fwi++) {
					Boolean fw = fwi == 0;
					EList<QVal> rrs = (fw ? hitsFw_ : hitsRc_);
					if(rrs[i].valid() && rrs[i].numElts() > 0) {
						rankOffs_.push_back(i);
						rankFws_.push_back(fw);
					}
				}
			}
			if(hitsFw_[0].valid() && hitsFw_[0].numElts() > 0) {
				rankOffs_.push_back(0);
				rankFws_.push_back(true);
			}
			if(hitsRc_[0].valid() && hitsRc_[0].numElts() > 0) {
				rankOffs_.push_back(0);
				rankFws_.push_back(false);
			}
		} else {
			while(rankOffs_.size() < nonzTot_) {
				TIndexOffU minsz = MAX_U32;
				int minidx = 0;
				Boolean minfw = true;
				// Rank seed-hit positions in ascending order by number of elements
				// in all BW ranges
				Boolean rb = rnd.nextBool();
				assert(rb == 0 || rb == 1);
				for(int fwi = 0; fwi <= 1; fwi++) {
					Boolean fw = (fwi == (rb ? 1 : 0));
					EList<QVal> rrs = (fw ? hitsFw_ : hitsRc_);
					EList<Boolean> sorted = (fw ? sortedFw_ : sortedRc_);
					int i = (rnd.nextU32() % (int)numOffs_);
					for(int ii = 0; ii < numOffs_; ii++) {
						if(rrs[i].valid() &&         // valid QVal
						   rrs[i].numElts() > 0 &&   // non-empty
						   !sorted[i] &&             // not already sorted
						   rrs[i].numElts() < minsz) // least elts so far?
						{
							minsz = rrs[i].numElts();
							minidx = i;
							minfw = (fw == 1);
						}
						if((++i) == numOffs_) {
							i = 0;
						}
					}
				}
				////assert_neq(MAX_U32, minsz);
				if(minfw) {
					sortedFw_[minidx] = true;
				} else {
					sortedRc_[minidx] = true;
				}
				rankOffs_.push_back(minidx);
				rankFws_.push_back(minfw);
			}
		}
		////assert_eq(rankOffs_.size(), rankFws_.size());
		sorted_ = true;
	}

	/**
	 * Return the number of orientation/offsets into the read that have
	 * at least one seed hit.
	 */
	public int nonzeroOffsets()  {
		assert(!sorted_ || nonzTot_ == rankFws_.size());
		assert(!sorted_ || nonzTot_ == rankOffs_.size());
		return nonzTot_;
	}
	
	/**
	 * Return true iff all seeds hit for forward read.
	 */
	public Boolean allFwSeedsHit()  {
		return nonzFw_ == numOffs();
	}

	/**
	 * Return true iff all seeds hit for revcomp read.
	 */
	public Boolean allRcSeedsHit()  {
		return nonzRc_ == numOffs();
	}
	
	/**
	 * Return the minimum number of edits that an end-to-end alignment of the
	 * fw read could have.  Uses knowledge of how many seeds have exact hits
	 * and how the seeds overlap.
	 */
	public int fewestEditsEE(Boolean fw, int seedlen, int per)  {
		////assert_gt(seedlen, 0);
		////assert_gt(per, 0);
		int nonz = fw ? nonzFw_ : nonzRc_;
		if(nonz < numOffs()) {
			int maxdepth = (seedlen + per - 1) / per;
			int missing = (int)(numOffs() - nonz);
			return (missing + maxdepth - 1) / maxdepth;
		} else {
			// Exact hit is possible (not guaranteed)
			return 0;
		}
	}

	/**
	 * Return the number of offsets into the forward read that have at
	 * least one seed hit.
	 */
	public int nonzeroOffsetsFw()  {
		return nonzFw_;
	}
	
	/**
	 * Return the number of offsets into the reverse-complement read
	 * that have at least one seed hit.
	 */
	public int nonzeroOffsetsRc()  {
		return nonzRc_;
	}

	/**
	 * Return a QVal of seed hits of the given rank 'r'.  'offidx' gets the id
	 * of the offset from 5' from which it was extracted (0 for the 5-most
	 * offset, 1 for the next closes to 5', etc).  'off' gets the offset from
	 * the 5' end.  'fw' gets true iff the seed was extracted from the forward
	 * read.
	 */
	public QVal hitsByRank(
		int    r,       // in
		int offidx,  // out
		int off,     // out
		Boolean     fw,      // out
		int seedlen) // out
	{
		assert(sorted_);
		////assert_lt(r, nonzTot_);
		if(rankFws_[r]) {
			fw = true;
			offidx = rankOffs_[r];
			//assert_lt(offidx, offIdx2off_.size());
			off = offIdx2off_[offidx];
			seedlen = (int)seqFw_[rankOffs_[r]].length();
			return hitsFw_[rankOffs_[r]];
		} else {
			fw = false;
			offidx = rankOffs_[r];
			////assert_lt(offidx, offIdx2off_.size());
			off = offIdx2off_[offidx];
			seedlen = (int)seqRc_[rankOffs_[r]].length();
			return hitsRc_[rankOffs_[r]];
		}
	}

	/**
	 * Return an EList of seed hits of the given rank.
	 */
	public BTDnaString seqByRank(int r) {
		assert(sorted_);
		////assert_lt(r, nonzTot_);
		return rankFws_[r] ? seqFw_[rankOffs_[r]] : seqRc_[rankOffs_[r]];
	}

	/**
	 * Return an EList of seed hits of the given rank.
	 */
	 public BTString qualByRank(int r) {
		assert(sorted_);
		//assert_lt(r, nonzTot_);
		return rankFws_[r] ? qualFw_[rankOffs_[r]] : qualRc_[rankOffs_[r]];
	}
	
	/**
	 * Return the list of extracted seed sequences for seeds on either
	 * the forward or reverse strand.
	 */
	public EList<BTDnaString> seqs(Boolean fw) { return fw ? seqFw_ : seqRc_; }

	/**
	 * Return the list of extracted quality sequences for seeds on
	 * either the forward or reverse strand.
	 */
	public EList<BTString> quals(Boolean fw) { return fw ? qualFw_ : qualRc_; }

	/**
	 * Return exact end-to-end alignment of fw read.
	 */
	public EEHit exactFwEEHit()  { return exactFwHit_; }

	/**
	 * Return exact end-to-end alignment of rc read.
	 */
	public EEHit exactRcEEHit()  { return exactRcHit_; }
	
	/**
	 * Return  ref to list of 1-mismatch end-to-end alignments.
	 */
	public EList<EEHit> mm1EEHits()  { return mm1Hit_; }
	
	/**
	 * Sort the end-to-end 1-mismatch alignments, prioritizing by score (higher
	 * score = higher priority).
	 */
	public void sort1mmEe(RandomSource rnd) {
		assert(!mm1Sorted_);
		mm1Hit_.sort();
		int streak = 0;
		for(int i = 1; i < mm1Hit_.size(); i++) {
			if(mm1Hit_[i].score == mm1Hit_[i-1].score) {
				if(streak == 0) { streak = 1; }
				streak++;
			} else {
				if(streak > 1) {
					//assert_geq(i, streak);
					mm1Hit_.shufflePortion(i-streak, streak, rnd);
				}
				streak = 0;
			}
		}
		if(streak > 1) {
			mm1Hit_.shufflePortion(mm1Hit_.size() - streak, streak, rnd);
		}
		mm1Sorted_ = true;
	}
	
	/**
	 * Add an end-to-end 1-mismatch alignment.
	 */
	public void add1mmEe(
		TIndexOffU top,
		TIndexOffU bot,
		 Edit e1,
		 Edit e2,
		Boolean fw,
		long score)
	{
		mm1Hit_.expand();
		mm1Hit_.back().init(top, bot, e1, e2, fw, score);
		mm1Elt_ += (bot - top);
	}

	/**
	 * Add an end-to-end exact alignment.
	 */
	public void addExactEeFw(
		TIndexOffU top,
		TIndexOffU bot,
		 Edit e1,
		 Edit e2,
		Boolean fw,
		long score)
	{
		exactFwHit_.init(top, bot, e1, e2, fw, score);
	}

	/**
	 * Add an end-to-end exact alignment.
	 */
	public void addExactEeRc(
		TIndexOffU top,
		TIndexOffU bot,
		 Edit e1,
		 Edit e2,
		Boolean fw,
		long score)
	{
		exactRcHit_.init(top, bot, e1, e2, fw, score);
	}
	
	/**
	 * Clear out the end-to-end exact alignments.
	 */
	public void clearExactE2eHits() {
		exactFwHit_.reset();
		exactRcHit_.reset();
	}
	
	/**
	 * Clear out the end-to-end 1-mismatch alignments.
	 */
	public void clear1mmE2eHits() {
		mm1Hit_.clear();     // 1-mismatch end-to-end hits
		mm1Elt_ = 0;         // number of 1-mismatch hit rows
		mm1Sorted_ = false;  // true iff we've sorted the mm1Hit_ list
	}

	/**
	 * Return the number of distinct exact and 1-mismatch end-to-end hits
	 * found.
	 */
	public int numE2eHits()  {
		return exactFwHit_.size() + exactRcHit_.size() + mm1Elt_;
	}

	/**
	 * Return the number of distinct exact end-to-end hits found.
	 */
	public int numExactE2eHits()  {
		return exactFwHit_.size() + exactRcHit_.size();
	}

	/**
	 * Return the number of distinct 1-mismatch end-to-end hits found.
	 */
	public int num1mmE2eHits()  {
		return mm1Elt_;
	}
	
	/**
	 * Return the length of the read that yielded the seed hits.
	 */
	public int readLength()  {
		assert(read_ != null);
		return read_.length();
	}
}


// Forward decl
class Ebwt;
class SideLocus;

/**
 * Encapsulates a sumamry of what the searchAllSeeds aligner did.
 */
class SeedSearchMetrics {

	public SeedSearchMetrics() : mutex_m() {
	    reset();
	}

	/**
	 * Merge this metrics object with the given object, i.e., sum each
	 * category.  This is the only safe way to update a
	 * SeedSearchMetrics object shread by multiple threads.
	 */
	public void merge( SeedSearchMetrics m, Boolean getLock = false) {
		seedsearch   += m.seedsearch;
		nrange       += m.nrange;
		nelt         += m.nelt;
		possearch    += m.possearch;
		intrahit     += m.intrahit;
		interhit     += m.interhit;
		filteredseed += m.filteredseed;
		ooms         += m.ooms;
		bwops        += m.bwops;
		bweds        += m.bweds;
		bestmin0     += m.bestmin0;
		bestmin1     += m.bestmin1;
		bestmin2     += m.bestmin2;
	}
	
	/**
	 * Set all counters to 0.
	 */
	public void reset() {
		seedsearch =
		nrange =
		nelt =
		possearch =
		intrahit =
		interhit =
		filteredseed =
		ooms =
		bwops =
		bweds =
		bestmin0 =
		bestmin1 =
		bestmin2 = 0;
	}

	long seedsearch;   // // times we executed strategy in InstantiatedSeed
	long nrange;       // // ranges found
	long nelt;         // // range elements found
	long possearch;    // // offsets where aligner executed >= 1 strategy
	long intrahit;     // // offsets where current-read cache gave answer
	long interhit;     // // offsets where across-read cache gave answer
	long filteredseed; // // seed instantiations skipped due to Ns
	long ooms;         // out-of-memory errors
	long bwops;        // Burrows-Wheeler operations
	long bweds;        // Burrows-Wheeler edits
	long bestmin0;     // // times the best min // edits was 0
	long bestmin1;     // // times the best min // edits was 1
	long bestmin2;     // // times the best min // edits was 2
	MUTEX_T  mutex_m;
}

class Pair {
	public int first;
	public int second;
	
	public Pair(int f, int s) {
		first = f;
		second = s;
	}
}

/**
 * Given an index and a seeding scheme, searches for seed hits.
 */
class SeedAligner {
	
	// Following are set in searchAllSeeds then used by searchSeed()
	// and other protected members.
	 protected Ebwt ebwtFw_;       // forward index (BWT)
	 protected Ebwt ebwtBw_;       // backward/mirror index (BWT')
	 protected Scoring sc_;        // scoring scheme
	 protected InstantiatedSeed s_;// current instantiated seed
	
	 protected Read read_;         // read whose seeds are currently being aligned
	
	// The following are set just before a call to searchSeedBi()
	 protected BTDnaString seq_;   // sequence of current seed
	 protected BTString qual_;     // quality string for current seed
	protected int off_;               // offset of seed currently being searched
	protected Boolean fw_;                  // orientation of seed currently being searched
	
	protected EList<Edit> edits_;        // temporary place to sort edits
	protected AlignmentCacheIface ca_;  // local alignment cache for seed alignments
	protected EList<int> offIdx2off_;// offset idx to read offset map, set up instantiateSeeds()
	protected long bwops_;           // Burrows-Wheeler operations
	protected long bwedits_;         // Burrows-Wheeler edits
	protected BTDnaString tmprfdnastr_;  // used in reportHit
	
	protected ESet<BTDnaString> hits_; // Ref hits so far for seed being aligned
	protected BTDnaString tmpdnastr_;
	
	/**
	 * Initialize with index.
	 */
	public SeedAligner() : edits_(AL_CAT), offIdx2off_(AL_CAT) { }

	/**
	 * Given a read and a few coordinates that describe a substring of the
	 * read (or its reverse complement), fill in 'seq' and 'qual' objects
	 * with the seed sequence and qualities.
	 */
	public void instantiateSeq(
		 Read read, // input read
		BTDnaString seq, // output sequence
		BTString qual,   // output qualities
		int len,          // seed length
		int depth,        // seed's 0-based offset from 5' end
		Boolean fw)   // seed's orientation
	{
		// Fill in 'seq' and 'qual'
		int seedlen = len;
		if((int)read.length() < seedlen) seedlen = (int)read.length();
		seq.resize(len);
		qual.resize(len);
		// If fw is false, we take characters starting at the 3' end of the
		// reverse complement of the read.
		for(int i = 0; i < len; i++) {
			seq.set(read.patFw.windowGetDna(i, fw, depth, len), i);
			qual.set(read.qual.windowGet(i, fw, depth, len), i);
		}
	}

	/**
	 * Iterate through the seeds that cover the read and initiate a
	 * search for each seed.
	 */
	public Pair instantiateSeeds(
		 EList<Seed> seeds,   // search seeds
		int off,                 // offset into read to start extracting
		int per,                    // interval between seeds
		 Read read,           // read to align
		 Scoring pens,        // scoring scheme
		Boolean nofw,                  // don't align forward read
		Boolean norc,                  // don't align revcomp read
		AlignmentCacheIface cache, // holds some seed hits from previous reads
		SeedResults sr,            // holds all the seed hits
		SeedSearchMetrics met,     // metrics
		Pair instFw,
		Pair instRc)
	{
		assert(!seeds.empty());
		////assert_gt(read.length(), 0);
		// Check whether read has too many Ns
		offIdx2off_.clear();
		int len = seeds[0].len; // assume they're all the same length
	////#ifndef NDEBUG
		//for(int i = 1; i < seeds.size(); i++) {
			////assert_eq(len, seeds[i].len);
		//}
	////#endif
		// Calc # seeds within read interval
		int nseeds = 1;
		if((int)read.length() - (int)off > len) {
			nseeds += ((int)read.length() - (int)off - len) / per;
		}
		for(int i = 0; i < nseeds; i++) {
			offIdx2off_.push_back(per * i + (int)off);
		}
		Pair ret;
		ret.first = 0;  // # seeds that require alignment
		ret.second = 0; // # seeds that hit in cache with non-empty results
		sr.reset(read, offIdx2off_, nseeds);
		assert(sr.repOk(cache.current(), true)); // require that SeedResult be initialized
		// For each seed position
		for(int fwi = 0; fwi < 2; fwi++) {
			Boolean fw = (fwi == 0);
			if((fw && nofw) || (!fw && norc)) {
				// Skip this orientation b/c user specified --nofw or --norc
				continue;
			}
			// For each seed position
			for(int i = 0; i < nseeds; i++) {
				int depth = i * per + (int)off;
				int seedlen = seeds[0].len;
				// Extract the seed sequence at this offset
				// If fw == true, we extract the characters from i*per to
				// i*(per-1) (exclusive).  If fw == false, 
				instantiateSeq(
					read,
					sr.seqs(fw)[i],
					sr.quals(fw)[i],
					Math.min((int)seedlen, (int)read.length()),
					depth,
					fw);
				QKey qk(sr.seqs(fw)[i], tmpdnastr_);
				// For each search strategy
				EList<InstantiatedSeed> iss = sr.instantiatedSeeds(fw, i);
				for(int j = 0; j < (int)seeds.size(); j++) {
					iss.expand();
					////assert_eq(seedlen, seeds[j].len);
					InstantiatedSeed is = iss.back();
					if(seeds[j].instantiate(
						read,
						sr.seqs(fw)[i],
						sr.quals(fw)[i],
						pens,
						depth,
						i,
						j,
						fw,
						is))
					{
						// Can we fill this seed hit in from the cache?
						ret.first++;
						if(fwi == 0) { instFw.first++; } else { instRc.first++; }
					} else {
						// Seed may fail to instantiate if there are Ns
						// that prevent it from matching
						met.filteredseed++;
						iss.pop_back();
					}
				}
			}
		}
		return ret;
	}

	/**
	 * Iterate through the seeds that cover the read and initiate a
	 * search for each seed.
	 */
	public void searchAllSeeds(
		 EList<Seed> seeds,   // search seeds
		 Ebwt ebwtFw,         // BWT index
		 Ebwt ebwtBw,         // BWT' index
		 Read read,           // read to align
		 Scoring pens,        // scoring scheme
		AlignmentCacheIface cache, // local seed alignment cache
		SeedResults hits,          // holds all the seed hits
		SeedSearchMetrics met,     // metrics
		PerReadMetrics prm)       // per-read metrics
	{
		assert(!seeds.empty());
		assert(ebwtFw != null);
		assert(ebwtFw.isInMemory());
		assert(sr.repOk(&cache.current()));
		ebwtFw_ = ebwtFw;
		ebwtBw_ = ebwtBw;
		sc_ = pens;
		read_ = read;
		ca_ = cache;
		bwops_ = bwedits_ = 0;
		long possearches = 0, seedsearches = 0, intrahits = 0, interhits = 0, ooms = 0;
		// For each instantiated seed
		for(int i = 0; i < (int)sr.numOffs(); i++) {
			int off = sr.idx2off(i);
			for(int fwi = 0; fwi < 2; fwi++) {
				Boolean fw = (fwi == 0);
				assert(sr.repOk(cache.current()));
				EList<InstantiatedSeed> iss = sr.instantiatedSeeds(fw, i);
				if(iss.empty()) {
					// Cache hit in an across-read cache
					continue;
				}
				QVal qv;
				seq_  = sr.seqs(fw)[i];  // seed sequence
				qual_ = sr.quals(fw)[i]; // seed qualities
				off_  = off;              // seed offset (from 5')
				fw_   = fw;               // seed orientation
				// Tell the cache that we've started aligning, so the cache can
				// expect a series of on-the-fly updates
				int ret = cache.beginAlign(seq_, qual_, qv);
				hits_.clear();
				if(ret == -1) {
					// Out of memory when we tried to add key to map
					ooms++;
					continue;
				}
				Boolean abort = false;
				if(ret == 0) {
					// Not already in cache
					assert(cache.aligning());
					possearches++;
					for(int j = 0; j < iss.size(); j++) {
						// Set seq_ and qual_ appropriately, using the seed sequences
						// and qualities already installed in SeedResults
						//assert_eq(fw, iss[j].fw);
						//assert_eq(i, (int)iss[j].seedoffidx);
						s_ = &iss[j];
						// Do the search with respect to seq_, qual_ and s_.
						if(!searchSeedBi()) {
							// Memory exhausted during search
							ooms++;
							abort = true;
							break;
						}
						seedsearches++;
						assert(cache.aligning());
					}
					if(!abort) {
						qv = cache.finishAlign();
					}
				} else {
					// Already in cache
					////assert_eq(1, ret);
					assert(qv.valid());
					intrahits++;
				}
				assert(abort || !cache.aligning());
				if(qv.valid()) {
					sr.add(
						qv,    // range of ranges in cache
						cache.current(), // cache
						i,     // seed index (from 5' end)
						fw);   // whether seed is from forward read
				}
			}
		}
		prm.nSeedRanges = sr.numRanges();
		prm.nSeedElts = sr.numElts();
		prm.nSeedRangesFw = sr.numRangesFw();
		prm.nSeedRangesRc = sr.numRangesRc();
		prm.nSeedEltsFw = sr.numEltsFw();
		prm.nSeedEltsRc = sr.numEltsRc();
		prm.seedMedian = (long)(sr.medianHitsPerSeed() + 0.5);
		prm.seedMean = (long)sr.averageHitsPerSeed();

		prm.nSdFmops += bwops_;
		met.seedsearch += seedsearches;
		met.nrange += sr.numRanges();
		met.nelt += sr.numElts();
		met.possearch += possearches;
		met.intrahit += intrahits;
		met.interhit += interhits;
		met.ooms += ooms;
		met.bwops += bwops_;
		met.bweds += bwedits_;
	}
	/**
	 * Sanity-check a partial alignment produced during oneMmSearch.
	 */
	public Boolean sanityPartial(
		 Ebwt        ebwtFw, // BWT index
		 Ebwt        ebwtBw, // BWT' index
		 BTDnaString seq,
		int             dep,
		int             len,
		Boolean               do1mm,
		TIndexOffU           topfw,
		TIndexOffU           botfw,
		TIndexOffU           topbw,
		TIndexOffU           botbw)
	{
		tmpdnastr_.clear();
		for(int i = dep; i < len; i++) {
			tmpdnastr_.append(seq[i]);
		}
		TIndexOffU top_fw = 0, bot_fw = 0;
		ebwtFw.contains(tmpdnastr_, top_fw, bot_fw);
		////assert_eq(top_fw, topfw);
		////assert_eq(bot_fw, botfw);
		if(do1mm && ebwtBw != null) {
			tmpdnastr_.reverse();
			TIndexOffU top_bw = 0, bot_bw = 0;
			ebwtBw.contains(tmpdnastr_, top_bw, bot_bw);
			////assert_eq(top_bw, topbw);
			////assert_eq(bot_bw, botbw);
		}
		return true;
	}

	/**
	 * Do an exact-matching sweet to establish a lower bound on number of edits
	 * and to find exact alignments.
	 */
	public int exactSweep(
		 Ebwt        ebwt,    // BWT index
		 Read        read,    // read to align
		 Scoring     sc,      // scoring scheme
		Boolean               nofw,    // don't align forward read
		Boolean               norc,    // don't align revcomp read
		int             mineMax, // don't care about edit bounds > this
		int            mineFw,  // minimum // edits for forward read
		int            mineRc,  // minimum // edits for revcomp read
		Boolean               repex,   // report 0mm hits?
		SeedResults       hits,    // holds all the seed hits (and exact hit)
		SeedSearchMetrics met)    // metrics
	{
		TIndexOffU top = 0, bot = 0;
		SideLocus tloc, bloc;
		 int len = read.length();
		int nelt = 0;
		for(int fwi = 0; fwi < 2; fwi++) {
			Boolean fw = (fwi == 0);
			if( fw && nofw) continue;
			if(!fw && norc) continue;
			 BTDnaString seq = fw ? read.patFw : read.patRc;
			assert(!seq.empty());
			int ftabLen = ebwt.eh().ftabChars();
			int dep = 0;
			int nedit = 0;
			Boolean done = false;
			while(dep < len && !done) {
				top = bot = 0;
				int left = len - dep;
				////assert_gt(left, 0);
				Boolean doFtab = ftabLen > 1 && left >= (int)ftabLen;
				if(doFtab) {
					// Does N interfere with use of Ftab?
					for(int i = 0; i < (int)ftabLen; i++) {
						int c = seq[len-dep-1-i];
						if(c > 3) {
							doFtab = false;
							break;
						}
					}
				}
				if(doFtab) {
					// Use ftab
					ebwt.ftabLoHi(seq, len - dep - ftabLen, false, top, bot);
					dep += (int)ftabLen;
				} else {
					// Use fchr
					int c = seq[len-dep-1];
					if(c < 4) {
						top = ebwt.fchr()[c];
						bot = ebwt.fchr()[c+1];
					}
					dep++;
				}
				if(bot <= top) {
					nedit++;
					if(nedit >= mineMax) {
						if(fw) { mineFw = nedit; } else { mineRc = nedit; }
						break;
					}
					continue;
				}
				INIT_LOCS(top, bot, tloc, bloc, ebwt);
				// Keep going
				while(dep < len) {
					int c = seq[len-dep-1];
					if(c > 3) {
						top = bot = 0;
					} else {
						if(bloc.valid()) {
							bwops_ += 2;
							top = ebwt.mapLF(tloc, c);
							bot = ebwt.mapLF(bloc, c);
						} else {
							bwops_++;
							top = ebwt.mapLF1(top, tloc, c);
							if(top == OFF_MASK) {
								top = bot = 0;
							} else {
								bot = top+1;
							}
						}
					}
					if(bot <= top) {
						nedit++;
						if(nedit >= mineMax) {
							if(fw) { mineFw = nedit; } else { mineRc = nedit; }
							done = true;
						}
						break;
					}
					INIT_LOCS(top, bot, tloc, bloc, ebwt);
					dep++;
				}
				if(done) {
					break;
				}
				if(dep == len) {
					// Set the minimum # edits
					if(fw) { mineFw = nedit; } else { mineRc = nedit; }
					// Done
					if(nedit == 0 && bot > top) {
						if(repex) {
							// This is an exact hit
							long score = len * sc.match();
							if(fw) {
								hits.addExactEeFw(top, bot, null, null, fw, score);
								assert(ebwt.contains(seq, null, null));
							} else {
								hits.addExactEeRc(top, bot, null, null, fw, score);
								assert(ebwt.contains(seq, null, null));
							}
						}
						nelt += (bot - top);
					}
					break;
				}
				dep++;
			}
		}
		return nelt;
	}
	/**
	 * Search for end-to-end alignments with up to 1 mismatch.
	 */
	public Boolean oneMmSearch(
		 Ebwt        ebwtFw, // BWT index
		 Ebwt        ebwtBw, // BWT' index
		 Read        read,   // read to align
		 Scoring     sc,     // scoring
		long            minsc,  // minimum score
		Boolean               nofw,   // don't align forward read
		Boolean               norc,   // don't align revcomp read
		Boolean               local,  // 1mm hits must be legal local alignments
		Boolean               repex,  // report 0mm hits?
		Boolean               rep1mm, // report 1mm hits?
		SeedResults       hits,   // holds all the seed hits (and exact hit)
		SeedSearchMetrics met)   // metrics
	{
		assert(!rep1mm || ebwtBw != null);
		 int len = read.length();
		int nceil = sc.nCeil.f<int>((double)len);
		int ns = read.ns();
		if(ns > 1) {
			// Can't align this with <= 1 mismatches
			return false;
		} else if(ns == 1 && !rep1mm) {
			// Can't align this with 0 mismatches
			return false;
		}
		////assert_geq(len, 2);
		assert(!rep1mm || ebwtBw.eh().ftabChars() == ebwtFw.eh().ftabChars());
	////#ifndef NDEBUG
		if(ebwtBw != null) {
			for(int i = 0; i < 4; i++) {
				//assert_eq(ebwtBw.fchr()[i], ebwtFw.fchr()[i]);
			}
		}
	////#endif
		int halfFw = len >> 1;
		int halfBw = len >> 1;
		if((len & 1) != 0) {
			halfBw++;
		}
		////assert_geq(halfFw, 1);
		////assert_geq(halfBw, 1);
		SideLocus tloc, bloc;
		TIndexOffU t[4], b[4];   // dest BW ranges for BWT
		t[0] = t[1] = t[2] = t[3] = 0;
		b[0] = b[1] = b[2] = b[3] = 0;
		TIndexOffU tp[4], bp[4]; // dest BW ranges for BWT'
		tp[0] = tp[1] = tp[2] = tp[3] = 0;
		bp[0] = bp[1] = bp[2] = bp[3] = 0;
		TIndexOffU top = 0, bot = 0, topp = 0, botp = 0;
		// Align fw read / rc read
		Boolean results = false;
		for(int fwi = 0; fwi < 2; fwi++) {
			Boolean fw = (fwi == 0);
			if( fw && nofw) continue;
			if(!fw && norc) continue;
			// Align going right-to-left, left-to-right
			int lim = rep1mm ? 2 : 1;
			for(int ebwtfwi = 0; ebwtfwi < lim; ebwtfwi++) {
				Boolean ebwtfw = (ebwtfwi == 0);
				 Ebwt ebwt  = (ebwtfw ? ebwtFw : ebwtBw);
				 Ebwt ebwtp = (ebwtfw ? ebwtBw : ebwtFw);
				assert(rep1mm || ebwt.fw());
				 BTDnaString& seq =
					(fw ? (ebwtfw ? read.patFw : read.patFwRev) :
						  (ebwtfw ? read.patRc : read.patRcRev));
				assert(!seq.empty());
				 BTString qual =
					(fw ? (ebwtfw ? read.qual    : read.qualRev) :
						  (ebwtfw ? read.qualRev : read.qual));
				int ftabLen = ebwt.eh().ftabChars();
				int nea = ebwtfw ? halfFw : halfBw;
				// Check if there's an N in the near portion
				Boolean skip = false;
				for(int dep = 0; dep < nea; dep++) {
					if(seq[len-dep-1] > 3) {
						skip = true;
						break;
					}
				}
				if(skip) {
					continue;
				}
				int dep = 0;
				// Align near half
				if(ftabLen > 1 && (int)ftabLen <= nea) {
					// Use ftab to jump partway into near half
					Boolean rev = !ebwtfw;
					ebwt.ftabLoHi(seq, len - ftabLen, rev, top, bot);
					if(rep1mm) {
						ebwtp.ftabLoHi(seq, len - ftabLen, rev, topp, botp);
						////assert_eq(bot - top, botp - topp);
					}
					if(bot - top == 0) {
						continue;
					}
					int c = seq[len - ftabLen];
					t[c] = top; b[c] = bot;
					tp[c] = topp; bp[c] = botp;
					dep = ftabLen;
					// initialize tloc, bloc??
				} else {
					// Use fchr to jump in by 1 pos
					int c = seq[len-1];
					////assert_range(0, 3, c);
					top = topp = tp[c] = ebwt.fchr()[c];
					bot = botp = bp[c] = ebwt.fchr()[c+1];
					if(bot - top == 0) {
						continue;
					}
					dep = 1;
					// initialize tloc, bloc??
				}
				INIT_LOCS(top, bot, tloc, bloc, ebwt);
				assert(sanityPartial(ebwt, ebwtp, seq, len-dep, len, rep1mm, top, bot, topp, botp));
				Boolean do_continue = false;
				for(; dep < nea; dep++) {
					////assert_lt(dep, len);
					int rdc = seq[len - dep - 1];
					tp[0] = tp[1] = tp[2] = tp[3] = topp;
					bp[0] = bp[1] = bp[2] = bp[3] = botp;
					if(bloc.valid()) {
						bwops_++;
						t[0] = t[1] = t[2] = t[3] = b[0] = b[1] = b[2] = b[3] = 0;
						ebwt.mapBiLFEx(tloc, bloc, t, b, tp, bp);
						SANITY_CHECK_4TUP(t, b, tp, bp);
						top = t[rdc]; bot = b[rdc];
						if(bot <= top) {
							do_continue = true;
							break;
						}
						topp = tp[rdc]; botp = bp[rdc];
						assert(!rep1mm || bot - top == botp - topp);
					} else {
						////assert_eq(bot, top+1);
						assert(!rep1mm || botp == topp+1);
						bwops_++;
						top = ebwt.mapLF1(top, tloc, rdc);
						if(top == OFF_MASK) {
							do_continue = true;
							break;
						}
						bot = top + 1;
						t[rdc] = top; b[rdc] = bot;
						tp[rdc] = topp; bp[rdc] = botp;
						assert(!rep1mm || b[rdc] - t[rdc] == bp[rdc] - tp[rdc]);
						// topp/botp stay the same
					}
					INIT_LOCS(top, bot, tloc, bloc, ebwt);
					assert(sanityPartial(ebwt, ebwtp, seq, len - dep - 1, len, rep1mm, top, bot, topp, botp));
				}
				if(do_continue) {
					continue;
				}
				// Align far half
				for(; dep < len; dep++) {
					int rdc = seq[len-dep-1];
					int quc = qual[len-dep-1];
					if(rdc > 3 && nceil == 0) {
						break;
					}
					tp[0] = tp[1] = tp[2] = tp[3] = topp;
					bp[0] = bp[1] = bp[2] = bp[3] = botp;
					int clo = 0, chi = 3;
					Boolean match = true;
					if(bloc.valid()) {
						bwops_++;
						t[0] = t[1] = t[2] = t[3] = b[0] = b[1] = b[2] = b[3] = 0;
						ebwt.mapBiLFEx(tloc, bloc, t, b, tp, bp);
						SANITY_CHECK_4TUP(t, b, tp, bp);
						match = rdc < 4;
						top = t[rdc]; bot = b[rdc];
						topp = tp[rdc]; botp = bp[rdc];
					} else {
						////assert_eq(bot, top+1);
						assert(!rep1mm || botp == topp+1);
						bwops_++;
						clo = ebwt.mapLF1(top, tloc);
						match = (clo == rdc);
						////assert_range(-1, 3, clo);
						if(clo < 0) {
							break; // Hit the $
						} else {
							t[clo] = top;
							b[clo] = bot = top + 1;
						}
						bp[clo] = botp;
						tp[clo] = topp;
						assert(!rep1mm || bot - top == botp - topp);
						assert(!rep1mm || b[clo] - t[clo] == bp[clo] - tp[clo]);
						chi = clo;
					}
					//assert(sanityPartial(ebwt, ebwtp, seq, len - dep - 1, len, rep1mm, top, bot, topp, botp));
					if(rep1mm && (ns == 0 || rdc > 3)) {
						for(int j = clo; j <= chi; j++) {
							if(j == rdc || b[j] == t[j]) {
								// Either matches read or isn't a possibility
								continue;
							}
							// Potential mismatch - next, try
							int depm = dep + 1;
							TIndexOffU topm = t[j], botm = b[j];
							TIndexOffU topmp = tp[j], botmp = bp[j];
							////assert_eq(botm - topm, botmp - topmp);
							TIndexOffU tm[4], bm[4];   // dest BW ranges for BWT
							tm[0] = t[0]; tm[1] = t[1];
							tm[2] = t[2]; tm[3] = t[3];
							bm[0] = b[0]; bm[1] = t[1];
							bm[2] = b[2]; bm[3] = t[3];
							TIndexOffU tmp[4], bmp[4]; // dest BW ranges for BWT'
							tmp[0] = tp[0]; tmp[1] = tp[1];
							tmp[2] = tp[2]; tmp[3] = tp[3];
							bmp[0] = bp[0]; bmp[1] = tp[1];
							bmp[2] = bp[2]; bmp[3] = tp[3];
							SideLocus tlocm, blocm;
							INIT_LOCS(topm, botm, tlocm, blocm, *ebwt);
							for(; depm < len; depm++) {
								int rdcm = seq[len - depm - 1];
								tmp[0] = tmp[1] = tmp[2] = tmp[3] = topmp;
								bmp[0] = bmp[1] = bmp[2] = bmp[3] = botmp;
								if(blocm.valid()) {
									bwops_++;
									tm[0] = tm[1] = tm[2] = tm[3] =
									bm[0] = bm[1] = bm[2] = bm[3] = 0;
									ebwt.mapBiLFEx(tlocm, blocm, tm, bm, tmp, bmp);
									SANITY_CHECK_4TUP(tm, bm, tmp, bmp);
									topm = tm[rdcm]; botm = bm[rdcm];
									topmp = tmp[rdcm]; botmp = bmp[rdcm];
									if(botm <= topm) {
										break;
									}
								} else {
									////assert_eq(botm, topm+1);
									////assert_eq(botmp, topmp+1);
									bwops_++;
									topm = ebwt.mapLF1(topm, tlocm, rdcm);
									if(topm == OFF_MASK) {
										break;
									}
									botm = topm + 1;
									// topp/botp stay the same
								}
								INIT_LOCS(topm, botm, tlocm, blocm, *ebwt);
							}
							if(depm == len) {
								// Success; this is a 1MM hit
								int off5p = dep;  // offset from 5' end of read
								int offstr = dep; // offset into patFw/patRc
								if(fw == ebwtfw) {
									off5p = len - off5p - 1;
								}
								if(!ebwtfw) {
									offstr = len - offstr - 1;
								}
								Edit e((int)off5p, j, rdc, EDIT_TYPE_MM, false);
								results = true;
								long score = (len - 1) * sc.match();
								// In --local mode, need to double-check that
								// end-to-end alignment doesn't violate  local
								// alignment principles.  Specifically, it
								// shouldn't to or below 0 anywhere in the middle.
								int pen = sc.score(rdc, (int)(1 << j), quc - 33);
								score += pen;
								Boolean valid = true;
								if(local) {
									long locscore_fw = 0, locscore_bw = 0;
									for(int i = 0; i < len; i++) {
										if(i == dep) {
											if(locscore_fw + pen <= 0) {
												valid = false;
												break;
											}
											locscore_fw += pen;
										} else {
											locscore_fw += sc.match();
										}
										if(len-i-1 == dep) {
											if(locscore_bw + pen <= 0) {
												valid = false;
												break;
											}
											locscore_bw += pen;
										} else {
											locscore_bw += sc.match();
										}
									}
								}
								if(valid) {
									valid = score >= minsc;
								}
								if(valid) {
	//#ifndef NDEBUG
									BTDnaString& rf = tmprfdnastr_;
									rf.clear();
									edits_.clear();
									edits_.push_back(e);
									if(!fw) Edit::invertPoss(edits_, len, false);
									Edit::toRef(fw ? read.patFw : read.patRc, edits_, rf);
									if(!fw) Edit::invertPoss(edits_, len, false);
									//assert_eq(len, rf.length());
									for(int i = 0; i < len; i++) {
										//assert_lt((int)rf[i], 4);
									}
									TIndexOffU toptmp = 0;
									TIndexOffU bottmp = 0;
									assert(ebwtFw.contains(rf, toptmp, bottmp));
	//#endif
									TIndexOffU toprep = ebwtfw ? topm : topmp;
									TIndexOffU botrep = ebwtfw ? botm : botmp;
									//assert_eq(toprep, toptmp);
									//assert_eq(botrep, bottmp);
									hits.add1mmEe(toprep, botrep, e, null, fw, score);
								}
							}
						}
					}
					if(bot > top && match) {
						////assert_lt(rdc, 4);
						if(dep == len-1) {
							// Success; this is an exact hit
							if(ebwtfw && repex) {
								if(fw) {
									results = true;
									long score = len * sc.match();
									hits.addExactEeFw(
										ebwtfw ? top : topp,
										ebwtfw ? bot : botp,
										null, null, fw, score);
									assert(ebwtFw.contains(seq, null, null));
								} else {
									results = true;
									long score = len * sc.match();
									hits.addExactEeRc(
										ebwtfw ? top : topp,
										ebwtfw ? bot : botp,
										null, null, fw, score);
									assert(ebwtFw.contains(seq, null, null));
								}
							}
							break; // End of far loop
						} else {
							INIT_LOCS(top, bot, tloc, bloc, *ebwt);
							assert(sanityPartial(ebwt, ebwtp, seq, len - dep - 1, len, rep1mm, top, bot, topp, botp));
						}
					} else {
						break; // End of far loop
					}
				} // for(; dep < len; dep++)
			} // for(int ebwtfw = 0; ebwtfw < 2; ebwtfw++)
		} // for(int fw = 0; fw < 2; fw++)
		return results;
	}
	/**
	 * Report a seed hit found by searchSeedBi(), but first try to extend it out in
	 * either direction as far as possible without hitting any edits.  This will
	 * allow us to prioritize the seed hits better later on.  Call reportHit() when
	 * we're done, which actually adds the hit to the cache.  Returns result from
	 * calling reportHit().
	 */
	protected Boolean extendAndReportHit(
		TIndexOffU topf,                     // top in BWT
		TIndexOffU botf,                     // bot in BWT
		TIndexOffU topb,                     // top in BWT'
		TIndexOffU botb,                     // bot in BWT'
		short len,                      // length of hit
		DoublyLinkedList<Edit> prevEdit) // previous edit
	{
		int nlex = 0, nrex = 0;
		TIndexOffU t[4], b[4];
		TIndexOffU tp[4], bp[4];
		SideLocus tloc, bloc;
		if(off_ > 0) {
			 Ebwt ebwt = ebwtFw_;
			assert(ebwt != null);
			// Extend left using forward index
			 BTDnaString& seq = fw_ ? read_.patFw : read_.patRc;
			// See what we get by extending 
			TIndexOffU top = topf, bot = botf;
			t[0] = t[1] = t[2] = t[3] = 0;
			b[0] = b[1] = b[2] = b[3] = 0;
			tp[0] = tp[1] = tp[2] = tp[3] = topb;
			bp[0] = bp[1] = bp[2] = bp[3] = botb;
			SideLocus tloc, bloc;
			INIT_LOCS(top, bot, tloc, bloc, *ebwt);
			for(int ii = off_; ii > 0; ii--) {
				int i = ii-1;
				// Get char from read
				int rdc = seq.get(i);
				// See what we get by extending 
				if(bloc.valid()) {
					bwops_++;
					t[0] = t[1] = t[2] = t[3] =
					b[0] = b[1] = b[2] = b[3] = 0;
					ebwt.mapBiLFEx(tloc, bloc, t, b, tp, bp);
					SANITY_CHECK_4TUP(t, b, tp, bp);
					int nonz = -1;
					Boolean abort = false;
					for(int j = 0; j < 4; j++) {
						if(b[i] > t[i]) {
							if(nonz >= 0) {
								abort = true;
								break;
							}
							nonz = j;
							top = t[i]; bot = b[i];
						}
					}
					if(abort || nonz != rdc) {
						break;
					}
				} else {
					////assert_eq(bot, top+1);
					bwops_++;
					int c = ebwt.mapLF1(top, tloc);
					if(c != rdc) {
						break;
					}
					bot = top + 1;
				}
				if(++nlex == 255) {
					break;
				}
				INIT_LOCS(top, bot, tloc, bloc, ebwt);
			}
		}
		int rdlen = read_.length();
		int nright = rdlen - off_ - len;
		if(nright > 0 && ebwtBw_ != null) {
			 Ebwt ebwt = ebwtBw_;
			assert(ebwt != null);
			// Extend right using backward index
			 BTDnaString seq = fw_ ? read_.patFw : read_.patRc;
			// See what we get by extending 
			TIndexOffU top = topb, bot = botb;
			t[0] = t[1] = t[2] = t[3] = 0;
			b[0] = b[1] = b[2] = b[3] = 0;
			tp[0] = tp[1] = tp[2] = tp[3] = topb;
			bp[0] = bp[1] = bp[2] = bp[3] = botb;
			INIT_LOCS(top, bot, tloc, bloc, ebwt);
			for(int i = off_ + len; i < rdlen; i++) {
				// Get char from read
				int rdc = seq.get(i);
				// See what we get by extending 
				if(bloc.valid()) {
					bwops_++;
					t[0] = t[1] = t[2] = t[3] =
					b[0] = b[1] = b[2] = b[3] = 0;
					ebwt.mapBiLFEx(tloc, bloc, t, b, tp, bp);
					SANITY_CHECK_4TUP(t, b, tp, bp);
					int nonz = -1;
					Boolean abort = false;
					for(int j = 0; j < 4; j++) {
						if(b[i] > t[i]) {
							if(nonz >= 0) {
								abort = true;
								break;
							}
							nonz = j;
							top = t[i]; bot = b[i];
						}
					}
					if(abort || nonz != rdc) {
						break;
					}
				} else {
					////assert_eq(bot, top+1);
					bwops_++;
					int c = ebwt.mapLF1(top, tloc);
					if(c != rdc) {
						break;
					}
					bot = top + 1;
				}
				if(++nrex == 255) {
					break;
				}
				INIT_LOCS(top, bot, tloc, bloc, ebwt);
			}
		}
		////assert_lt(nlex, rdlen);
		////assert_leq(nlex, off_);
		////assert_lt(nrex, rdlen);
		return reportHit(topf, botf, topb, botb, len, prevEdit);
	}
	/**
	 * Report a seed hit found by searchSeedBi() by adding it to the cache.  Return
	 * false if the hit could not be reported because of, e.g., cache exhaustion.
	 */
	protected Boolean reportHit(
		TIndexOffU topf,         // top in BWT
		TIndexOffU botf,         // bot in BWT
		TIndexOffU topb,         // top in BWT'
		TIndexOffU botb,         // bot in BWT'
		short len,          // length of hit
		DoublyLinkedList<Edit> prevEdit)  // previous edit
	{
		// Add information about the seed hit to AlignmentCache.  This
		// information eventually makes its way back to the SeedResults
		// object when we call finishAlign(...).
		BTDnaString rf = tmprfdnastr_;
		rf.clear();
		edits_.clear();
		if(prevEdit != null) {
			prevEdit.toList(edits_);
			Edit::sort(edits_);
			assert(Edit::repOk(edits_, seq_));
			Edit::toRef(seq_, edits_, rf);
		} else {
			rf = seq_;
		}
		// Sanity check: shouldn't add the same hit twice.  If this
		// happens, it may be because our zone Constraints are not set up
		// properly and erroneously return true from acceptable() when they
		// should return false in some cases.
		////assert_eq(hits_.size(), ca_.curNumRanges());
		assert(hits_.insert(rf));
		if(!ca_.addOnTheFly(rf, topf, botf, topb, botb)) {
			return false;
		}
		////assert_eq(hits_.size(), ca_.curNumRanges());
	//#ifndef NDEBUG
		// Sanity check that the topf/botf and topb/botb ranges really
		// correspond to the reference sequence aligned to
		{
			BTDnaString rfr;
			TIndexOffU tpf, btf, tpb, btb;
			tpf = btf = tpb = btb = 0;
			assert(ebwtFw_.contains(rf, tpf, btf));
			if(ebwtBw_ != null) {
				rfr = rf;
				rfr.reverse();
				assert(ebwtBw_.contains(rfr, tpb, btb));
				////assert_eq(tpf, topf);
				////assert_eq(btf, botf);
				////assert_eq(tpb, topb);
				////assert_eq(btb, botb);
			}
		}
	//#endif
		return true;
	}
	
	/**
	 * Given an instantiated seed (in s_ and other fields), search
	 */
	protected Boolean searchSeedBi()
	{
		return searchSeedBi(
		0, 0,
		0, 0, 0, 0,
		SideLocus(), SideLocus(),
		s_.cons[0], s_.cons[1], s_.cons[2], s_.overall,
		null);
	}
	
	/**
	 * Main, recursive implementation of the seed search.
	 */
	protected Boolean searchSeedBi(
		int step,              // depth into steps_[] array
		int depth,             // recursion depth
		TIndexOffU topf,         // top in BWT
		TIndexOffU botf,         // bot in BWT
		TIndexOffU topb,         // top in BWT'
		TIndexOffU botb,         // bot in BWT'
		SideLocus tloc,        // locus for top (perhaps unititialized)
		SideLocus bloc,        // locus for bot (perhaps unititialized)
		Constraint c0,         // raints to enforce in seed zone 0
		Constraint c1,         // raints to enforce in seed zone 1
		Constraint c2,         // raints to enforce in seed zone 2
		Constraint overall,    // overall raints
		DoublyLinkedList<Edit> prevEdit)  // previous edit
	{
		assert(s_ != null);
		 InstantiatedSeed s = s_;
		//assert_gt(s.steps.size(), 0);
		assert(ebwtBw_ == null || ebwtBw_.eh().ftabChars() == ebwtFw_.eh().ftabChars());
	//#ifndef NDEBUG
		for(int i = 0; i < 4; i++) {
			assert(ebwtBw_ == null || ebwtBw_.fchr()[i] == ebwtFw_.fchr()[i]);
		}
	//#endif
		if(step == (int)s.steps.size()) {
			// Finished aligning seed
			assert(c0.acceptable());
			assert(c1.acceptable());
			assert(c2.acceptable());
			if(!reportHit(topf, botf, topb, botb, seq_.length(), prevEdit)) {
				return false; // Memory exhausted
			}
			return true;
		}
	//#ifndef NDEBUG
		if(depth > 0) {
			assert(botf - topf == 1 ||  bloc.valid());
			assert(botf - topf > 1  || !bloc.valid());
		}
	//#endif
		int off;
		TIndexOffU tp[4], bp[4]; // dest BW ranges for "prime" index
		if(step == 0) {
			// Just starting
			assert(prevEdit == null);
			assert(!tloc.valid());
			assert(!bloc.valid());
			off = s.steps[0];
			Boolean ltr = off > 0;
			off = abs(off)-1;
			// Check whether/how far we can jump using ftab or fchr
			int ftabLen = ebwtFw_.eh().ftabChars();
			if(ftabLen > 1 && ftabLen <= s.maxjump) {
				if(!ltr) {
					////assert_geq(off+1, ftabLen-1);
					off = off - ftabLen + 1;
				}
				ebwtFw_.ftabLoHi(seq_, off, false, topf, botf);
				//#ifdef NDEBUG
				if(botf - topf == 0) return true;
				//#endif
				//#ifdef NDEBUG
				if(ebwtBw_ != null) {
					topb = ebwtBw_.ftabHi(seq_, off);
					botb = topb + (botf-topf);
				}
				//#else
				if(ebwtBw_ != null) {
					ebwtBw_.ftabLoHi(seq_, off, false, topb, botb);
					////assert_eq(botf-topf, botb-topb);
				}
				if(botf - topf == 0) return true;
				//#endif
				step += ftabLen;
			} else if(s.maxjump > 0) {
				// Use fchr
				int c = (seq_)[off];
				////assert_range(0, 3, c);
				topf = topb = ebwtFw_.fchr()[c];
				botf = botb = ebwtFw_.fchr()[c+1];
				if(botf - topf == 0) return true;
				step++;
			} else {
				////assert_eq(0, s.maxjump);
				topf = topb = 0;
				botf = botb = ebwtFw_.fchr()[4];
			}
			if(step == (int)s.steps.size()) {
				// Finished aligning seed
				assert(c0.acceptable());
				assert(c1.acceptable());
				assert(c2.acceptable());
				if(!reportHit(topf, botf, topb, botb, seq_.length(), prevEdit)) {
					return false; // Memory exhausted
				}
				return true;
			}
			nextLocsBi(tloc, bloc, topf, botf, topb, botb, step);
			assert(tloc.valid());
		} else assert(prevEdit != null);
		assert(tloc.valid());
		assert(botf - topf == 1 ||  bloc.valid());
		assert(botf - topf > 1  || !bloc.valid());
		////assert_geq(step, 0);
		TIndexOffU t[4], b[4]; // dest BW ranges
		Constraint zones[3] = { c0, c1, c2 };
		TIndexOffU lasttot = botf - topf;
		for(int i = step; i < (int)s.steps.size(); i++) {
			////assert_gt(botf, topf);
			assert(botf - topf == 1 ||  bloc.valid());
			assert(botf - topf > 1  || !bloc.valid());
			assert(ebwtBw_ == null || botf-topf == botb-topb);
			assert(tloc.valid());
			off = s.steps[i];
			Boolean ltr = off > 0;
			 Ebwt ebwt = ltr ? ebwtBw_ : ebwtFw_;
			assert(ebwt != null);
			if(ltr) {
				tp[0] = tp[1] = tp[2] = tp[3] = topf;
				bp[0] = bp[1] = bp[2] = bp[3] = botf;
			} else {
				tp[0] = tp[1] = tp[2] = tp[3] = topb;
				bp[0] = bp[1] = bp[2] = bp[3] = botb;
			}
			t[0] = t[1] = t[2] = t[3] = b[0] = b[1] = b[2] = b[3] = 0;
			if(bloc.valid()) {
				// Range delimited by tloc/bloc has size >1.  If size == 1,
				// we use a simpler query (see if(!bloc.valid()) blocks below)
				bwops_++;
				ebwt.mapBiLFEx(tloc, bloc, t, b, tp, bp);
				TIndexOffU tot = (b[0]-t[0])+(b[1]-t[1])+(b[2]-t[2])+(b[3]-t[3]);
				TIndexOffU totp = (bp[0]-tp[0])+(bp[1]-tp[1])+(bp[2]-tp[2])+(bp[3]-tp[3]);
				////assert_eq(tot, totp);
				////assert_leq(tot, lasttot);
				lasttot = tot;
			}
			TIndexOffU tf = ltr ? tp : t, tb = ltr ? t : tp;
			TIndexOffU bf = ltr ? bp : b, bb = ltr ? b : bp;
			off = abs(off)-1;
			//
			Boolean leaveZone = s.zones[i].first < 0;
			//Boolean leaveZoneIns = zones_[i].second < 0;
			Constraint cons    = zones[abs(s.zones[i].first)];
			//Constraint& insCons = *zones[abs(s.zones[i].second)];
			int c = (seq_)[off];  //assert_range(0, 4, c);
			int q = (qual_)[off];
			// Is it legal for us to advance on characters other than 'c'?
			if(!(cons.mustMatch() && !overall.mustMatch()) || c == 4) {
				// There may be legal edits
				Boolean bail = false;
				if(!bloc.valid()) {
					// Range delimited by tloc/bloc has size 1
					TIndexOffU ntop = ltr ? topb : topf;
					bwops_++;
					int cc = ebwt.mapLF1(ntop, tloc);
					////assert_range(-1, 3, cc);
					if(cc < 0) bail = true;
					else { t[cc] = ntop; b[cc] = ntop+1; }
				}
				if(!bail) {
					if((cons.canMismatch(q, sc_) && overall.canMismatch(q, sc_)) || c == 4) {
						Constraint oldCons = cons, oldOvCons = overall;
						SideLocus oldTloc = tloc, oldBloc = bloc;
						if(c != 4) {
							cons.chargeMismatch(q, sc_);
							overall.chargeMismatch(q, sc_);
						}
						// Can leave the zone as-is
						if(!leaveZone || (cons.acceptable() && overall.acceptable())) {
							for(int j = 0; j < 4; j++) {
								if(j == c || b[j] == t[j]) continue;
								// Potential mismatch
								nextLocsBi(tloc, bloc, tf[j], bf[j], tb[j], bb[j], i+1);
								int loff = off;
								if(!ltr) loff = (int)(s.steps.size() - loff - 1);
								assert(prevEdit == null || prevEdit.next == null);
								Edit edit(off, j, c, EDIT_TYPE_MM, false);
								DoublyLinkedList<Edit> editl;
								editl.payload = edit;
								if(prevEdit != null) {
									prevEdit.next = editl;
									editl.prev = prevEdit;
								}
								assert(editl.next == null);
								bwedits_++;
								if(!searchSeedBi(
									i+1,     // depth into steps_[] array
									depth+1, // recursion depth
									tf[j],   // top in BWT
									bf[j],   // bot in BWT
									tb[j],   // top in BWT'
									bb[j],   // bot in BWT'
									tloc,    // locus for top (perhaps unititialized)
									bloc,    // locus for bot (perhaps unititialized)
									c0,      // raints to enforce in seed zone 0
									c1,      // raints to enforce in seed zone 1
									c2,      // raints to enforce in seed zone 2
									overall, // overall raints to enforce
									editl))  // latest edit
								{
									return false;
								}
								if(prevEdit != null) prevEdit.next = null;
							}
						} else {
							// Not enough edits to make this path
							// non-redundant with other seeds
						}
						cons = oldCons;
						overall = oldOvCons;
						tloc = oldTloc;
						bloc = oldBloc;
					}
					if(cons.canGap() && overall.canGap()) {
						throw 1; // TODO
	//					int delEx = 0;
	//					if(cons.canDelete(delEx, *sc_) && overall.canDelete(delEx, *sc_)) {
	//						// Try delete
	//					}
	//					int insEx = 0;
	//					if(insCons.canInsert(insEx, *sc_) && overall.canInsert(insEx, *sc_)) {
	//						// Try insert
	//					}
					}
				} // if(!bail)
			}
			if(c == 4) {
				return true; // couldn't handle the N
			}
			if(leaveZone && (!cons.acceptable() || !overall.acceptable())) {
				// Not enough edits to make this path non-redundant with
				// other seeds
				return true;
			}
			if(!bloc.valid()) {
				assert(ebwtBw_ == null || bp[c] == tp[c]+1);
				// Range delimited by tloc/bloc has size 1
				TIndexOffU top = ltr ? topb : topf;
				bwops_++;
				t[c] = ebwt.mapLF1(top, tloc, c);
				if(t[c] == OFF_MASK) {
					return true;
				}
				//assert_geq(t[c], ebwt.fchr()[c]);
				//assert_lt(t[c],  ebwt.fchr()[c+1]);
				b[c] = t[c]+1;
				//assert_gt(b[c], 0);
			}
			assert(ebwtBw_ == null || bf[c]-tf[c] == bb[c]-tb[c]);
			//assert_leq(bf[c]-tf[c], lasttot);
			lasttot = bf[c]-tf[c];
			if(b[c] == t[c]) {
				return true;
			}
			topf = tf[c]; botf = bf[c];
			topb = tb[c]; botb = bb[c];
			if(i+1 == (int)s.steps.size()) {
				// Finished aligning seed
				assert(c0.acceptable());
				assert(c1.acceptable());
				assert(c2.acceptable());
				if(!reportHit(topf, botf, topb, botb, seq_.length(), prevEdit)) {
					return false; // Memory exhausted
				}
				return true;
			}
			nextLocsBi(tloc, bloc, tf[c], bf[c], tb[c], bb[c], i+1);
		}
		return true;
	}
	/**
	 * Get tloc and bloc ready for the next step.
	 */
	 // Was inline
	protected void nextLocsBi(
		SideLocus tloc,            // top locus
		SideLocus bloc,            // bot locus
		TIndexOffU topf,              // top in BWT
		TIndexOffU botf,              // bot in BWT
		TIndexOffU topb,              // top in BWT'
		TIndexOffU botb,              // bot in BWT'
		int step)                 // step to get ready for
	{
		assert(ebwtBw_ == null || botb > 0);
		//assert_geq(step, 0); // next step can't be first one
		assert(ebwtBw_ == null || botf-topf == botb-topb);
		if(step == (int)s_.steps.size()) return; // no more steps!
		// Which direction are we going in next?
		if(s_.steps[step] > 0) {
			// Left to right; use BWT'
			if(botb - topb == 1) {
				// Already down to 1 row; just init top locus
				tloc.initFromRow(topb, ebwtBw_.eh(), ebwtBw_.ebwt());
				bloc.invalidate();
			} else {
				SideLocus::initFromTopBot(
					topb, botb, ebwtBw_.eh(), ebwtBw_.ebwt(), tloc, bloc);
				assert(bloc.valid());
			}
		} else {
			// Right to left; use BWT
			if(botf - topf == 1) {
				// Already down to 1 row; just init top locus
				tloc.initFromRow(topf, ebwtFw_.eh(), ebwtFw_.ebwt());
				bloc.invalidate();
			} else {
				SideLocus::initFromTopBot(
					topf, botf, ebwtFw_.eh(), ebwtFw_.ebwt(), tloc, bloc);
				assert(bloc.valid());
			}
		}
		// Check if we should update the tracker with this refinement
	//#if 0
		if(botf-topf <= BW_OFF_TRACK_CEIL) {
			if(ot.size() == 0 && prevOt != null && prevOt.size() > 0) {
				// Inherit state from the predecessor
				ot = prevOt;
			}
			Boolean ltr = s_.steps[step-1] > 0;
			int adj = abs(s_.steps[step-1])-1;
			 Ebwt ebwt = ltr ? ebwtBw_ : ebwtFw_;
			ot.update(
				ltr ? topb : topf,    // top
				ltr ? botb : botf,    // bot
				adj,                  // adj (to be subtracted from offset)
				ebwt.offs(),         // offs array
				ebwt.eh().offRate(), // offrate (sample = every 1 << offrate elts)
				null                  // dead
			);
			//assert_gt(ot.size(), 0);
		}
	//#endif
		assert(botf - topf == 1 ||  bloc.valid());
		assert(botf - topf > 1  || !bloc.valid());
	}
}

//define INIT_LOCS(top, bot, tloc, bloc, e) { \
	if(bot - top == 1) { \
		tloc.initFromRow(top, (e).eh(), (e).ebwt()); \
		bloc.invalidate(); \
	} else { \
		SideLocus::initFromTopBot(top, bot, (e).eh(), (e).ebwt(), tloc, bloc); \
		assert(bloc.valid()); \
	} \
}

//define SANITY_CHECK_4TUP(t, b, tp, bp) { \
	ASSERT_ONLY(TIndexOffU tot = (b[0]-t[0])+(b[1]-t[1])+(b[2]-t[2])+(b[3]-t[3])); \
	ASSERT_ONLY(TIndexOffU totp = (bp[0]-tp[0])+(bp[1]-tp[1])+(bp[2]-tp[2])+(bp[3]-tp[3])); \
	//assert_eq(tot, totp); \
}

//endif /*ALIGNER_SEED_H_*/
