package com.uwb.bt2j.aligner;

import com.uwb.bt2j.indexer.Ebwt;

class AlignerDriver {

  	protected DescentRootSelector sel_;
  	protected DescentAlignmentSelector alsel_;
 	protected DescentDriver dr1_;
	protected DescentDriver dr2_;
  	protected DescentStoppingConditions stop_;
  	protected Boolean paired_;
  	protected SimpleFunc totsz_;
  	protected SimpleFunc totfmops_;
  	protected RedundantAlns  red1_;
  	protected RedundantAlns  red2_;
  
  public enum ALDRIVER {
    ALDRIVER_EXHAUSTED_CANDIDATES(1),
	  ALDRIVER_POLICY_FULFILLED(2),
	  ALDRIVER_EXCEEDED_LIMIT(3);
	private int x;
	  ALDRVER(int y) {x = y;}
  }
  
  public AlignerDriver(
  double consExp,
  Boolean prioritizeRoots,
  SimpleFunc rootIval,
  double landing,
  Boolean veryVerbose,
  SimpleFunc totsz,
  SimpleFunc totfmops,
  ) {
    	  alsel_();
	  dr1_ = veryVerbose;
	  dr2_ = veryVerbose;
	  totsz_ = totsz;
		totfmops_ = totfmops;
		if(prioritizeRoots) {
			// Prioritize roots according the quality info & Ns
			sel_ = new PrioritizedRootSelector(consExp, rootIval, landing);
		} else {
			// Take a root every so many positions
			sel_ = new IntervalRootSelector(consExp, rootIval, landing);
		}
  }
  
  public void initRead(Read q1, Boolean nofw, Boolean norc, long minsc, long maxpen, Read q2) {
    	// Initialize search for mate 1.  This includes instantiating and
		// prioritizing all the search roots.
		dr1_.initRead(q1, nofw, norc, minsc, maxpen, q2, sel_);
		red1_.init(q1.length());
		paired_ = false;
		if(q2 != null) {
			// Initialize search for mate 1.  This includes instantiating and
			// prioritizing all the search roots.
			dr2_.initRead(*q2, nofw, norc, minsc, maxpen, &q1, sel_);
			red2_.init(q2->length());
			paired_ = true;
		} else {
			dr2_.reset();
		}
		// Initialize stopping conditions.  We use two conditions:
		// totsz: when memory footprint exceeds this many bytes
		// totfmops: when we've exceeded this many FM Index ops
		double totsz = totsz_.f<double>(q1.length());
		double totfmops = totfmops_.f<double>(q1.length());
		stop_.init(
			totsz,
			0,
			true,
			totfmops);
  }
  
  public int go(Scoring sc, Ebwt ebwtFw, Ebwt ebwtBw, BitPairReference ref, DescentMetrics met, WalkMetrics wlm, PerReadMetrics prm, RandomSource rnd, AlnSinkWrap sink) {
    if(paired_) {
		// Paired-end - alternate between advancing dr1_ / dr2_ whenever a
		// new full alignment is discovered in the one currently being
		// advanced.  Whenever a new full alignment is found, check to see
		// if it pairs with a previously discovered alignment.
		Boolean first1 = rnd.nextBool();
		Boolean first = true;
		DescentStoppingConditions stopc1 = stop_;
		DescentStoppingConditions stopc2 = stop_;
		double totszIncr = (stop_.totsz + 7) / 8;
		stopc1.totsz = totszIncr;
		stopc2.totsz = totszIncr;
		while(stopc1.totsz <= stop_.totsz && stopc2.totsz <= stop_.totsz) {
			if(first && first1 && stopc1.totsz <= stop_.totsz) {
				dr1_.advance(stop_, sc, ebwtFw, ebwtBw, met, prm);
				stopc1.totsz += totszIncr;
			}
			if(stopc2.totsz <= stop_.totsz) {
				dr2_.advance(stop_, sc, ebwtFw, ebwtBw, met, prm);
				stopc2.totsz += totszIncr;
			}
			first = false;
		}
	} else {
		// Unpaired
		double iter = 1;
		while(true) {
			int ret = dr1_.advance(stop_, sc, ebwtFw, ebwtBw, met, prm);
			if(ret == DESCENT_DRIVER_ALN) {
				System.err.println(iter + ". DESCENT_DRIVER_ALN");
			} else if(ret == DESCENT_DRIVER_MEM) {
				System.err.println(iter + ". DESCENT_DRIVER_MEM");
				break;
			} else if(ret == DESCENT_DRIVER_STRATA) {
				// DESCENT_DRIVER_STRATA is returned by DescentDriver.advance()
				// when it has finished with a "non-empty" stratum: a stratum
				// in which at least one alignment was found.  Here we report
				// the alignments in an arbitrary order.
				AlnRes res;
				// Initialize alignment selector with the DescentDriver's
				// alignment sink
				alsel_.init(
					dr1_.query(),
					dr1_.sink(),
					ebwtFw,
					ref,
					rnd,
					wlm);
				while(!alsel_.done() && !sink.state().doneWithMate(true)) {
					res.reset();
					bool ret2 = alsel_.next(
						dr1_,
						ebwtFw,
						ref,
						rnd,
						res,
						wlm,
						prm);
					if(ret2) {
						// Get reference interval involved in alignment
						Interval refival(res.refid(), 0, res.fw(), res.reflen());
						// Does alignment falls off end of reference?
						if(gReportOverhangs &&
						   !refival.containsIgnoreOrient(res.refival()))
						{
							res.clipOutside(true, 0, res.reflen());
							if(res.refExtent() == 0) {
								continue;
							}
						}
						// Alignment fell entirely outside the reference?
						if(!refival.overlapsIgnoreOrient(res.refival())) {
							continue; // yes, fell outside
						}
						// Alignment redundant with one we've seen previously?
						if(red1_.overlap(res)) {
							continue; // yes, redundant
						}
						red1_.add(res); // so we find subsequent redundancies
						if(sink.report(0, &res, null)) {
							// Short-circuited because a limit, e.g. -k, -m or
							// -M, was exceeded
							return ALDRIVER_POLICY_FULFILLED;
						}
					}
				}
				dr1_.sink().advanceStratum();
			} else if(ret == DESCENT_DRIVER_BWOPS) {
				System.err.println(iter + ". DESCENT_DRIVER_BWOPS");
				break;
			} else if(ret == DESCENT_DRIVER_DONE) {
				System.err.println(iter + ". DESCENT_DRIVER_DONE");
				break;
			}
			iter++;
		}
	}
	return ALDRIVER_EXHAUSTED_CANDIDATES;
  }
  
  public void reset() {
    		dr1_.reset();
		dr2_.reset();
		red1_.reset();
		red2_.reset();
  }
  
  public final DescentDriver dr1() {
    return dr1_;
  }
  
  public final DescentDriver dr2() {
    return dr2_;
  }
  
  public class IntervalRootSelector() {
    protected double consExp_;
    protected SimpleFunc rootIval_;
    protected double landing_;
    
    public IntervalRootSelector(double consExp, SimpleFunc rootIval, double landing) {
      	consExp_ = consExp;
	rootIval_ = rootIval;
	landing_ = landing;
    }
    
    public void select(Read q, Read qo, Boolean nofw, Boolean norc, EList<DescentConfig> confs, EList<DescentRoot> roots) {
      int interval = rootIval_.f<Integer>((double)q.length());
	if(qo != null) {
		// Boost interval length by 20% for paired-end reads
		interval = (int)(interval * 1.2 + 0.5);
	}
	float pri = 0.0f;
	for(int fwi = 0; fwi < 2; fwi++) {
		Boolean fw = (fwi == 0);
		if((fw && nofw) || (!fw && norc)) {
			continue;
		}
		// Put down left-to-right roots w/r/t forward and reverse-complement reads
		{
			Boolean first = true;
			double i = 0;
			while(first || (i + landing_ <= q.length())) {
				confs.expand();
				confs.back().cons.init(landing_, consExp_);
				roots.expand();
				roots.back().init(
					i,          // offset from 5' end
					true,       // left-to-right?
					fw,         // fw?
					1,          // landing
					q.length(), // query length
					pri);       // root priority
				i += interval;
				first = false;
			}
		}
		// Put down right-to-left roots w/r/t forward and reverse-complement reads
		{
			Boolean first = true;
			double i = 0;
			while(first || (i + landing_ <= q.length())) {
				confs.expand();
				confs.back().cons.init(landing_, consExp_);
				roots.expand();
				roots.back().init(
					q.length() - i - 1, // offset from 5' end
					false,              // left-to-right?
					fw,                 // fw?
					1,          // landing
					q.length(),         // query length
					pri);               // root priority
				i += interval;
				first = false;
			}
		}
	}
    }
  }
  
  public class PrioritizedRootSelector() {
    protected double consExp_;
    protected SimpleFunc rootIval_;
    protected double landing_;
    EHeap<DescentRoot> rootHeap_;
    EList<int> scoresOrig_[2];
    EList<int> scores_[2];
    
    public PrioritizedRootSelector(double consExp, SimpleFunc rootIval, double landing) {
      	consExp_ = consExp;
	rootIval_ = rootIval;
	landing_ = landing;
    }
    
    public void select(Read q, Read qo, Boolean nofw, Boolean norc, EList<DescentConfig> confs, EList<DescentRoot> roots) {
      final int nPenalty = 150;
	final int endBonus = 150;
	final double qlen = q.length();
	// Calculate interval length
	int interval = rootIval_.f<Integer>((double)qlen);
	double sizeTarget = qlen - landing_ + 1;
	sizeTarget = (double)(Math.ceil((sizeTarget / (float)interval)));
	sizeTarget *= 4;
	// Set up initial score arrays
	for(int i = 0; i < 2; i++) {
		Boolean fw = (i == 0);
		scoresOrig_[i].resize(qlen);
		scores_[i].resize(qlen);
		for(double j = 0; j < qlen; j++) {
			double off5p = fw ? j : (qlen - j - 1);
			int c = q.getc(off5p, fw);
			int sc = q.getq(off5p) - ((c > 3) ? nPenalty : 0);
			scoresOrig_[i][j] = scores_[i][j] = sc;
		}
	}
	rootHeap_.clear();
	for(int fwi = 0; fwi < 2; fwi++) {
		Boolean fw = (fwi == 0);
		if((fw && nofw) || (!fw && norc)) {
			continue;
		}
		int pri = 0;
		double revi = qlen;
		for(double i = 0; i < qlen; i++) {
			revi--;
			pri += scoresOrig_[fwi][i];
			if(i >= landing_) {
				pri -= scoresOrig_[fwi][i - landing_];
			}
			if(i >= landing_-1 && scoresOrig_[fwi][i] > 0) {
				rootHeap_.insert(DescentRoot(
					fw ? i : revi, // offset from 5' end
					false,         // left-to-right?
					fw,            // fw?
					landing_,      // landing length
					qlen,          // query length
					pri + ((revi == 0) ? endBonus : 0))); // root priority
				// Give priority boost for being flush with one end or the
				// other
			}
		}
		pri = 0;
		double i = qlen - revi;
		for(double revi = 0; revi < qlen; revi++) {
			i--;
			pri += scoresOrig_[fwi][i];
			if(revi >= landing_) {
				pri -= scoresOrig_[fwi][i + landing_];
			}
			if(revi >= landing_-1 && scoresOrig_[fwi][i] > 0) {
				rootHeap_.insert(DescentRoot(
					fw ? i : revi, // offset from 5' end
					true,          // left-to-right?
					fw,            // fw?
					landing_,      // landing length
					qlen,          // query length
					pri + ((i == 0) ? endBonus : 0))); // root priority
				// Give priority boost for being flush with one end or the
				// other
			}
		}
	}
	// Now that all the roots are in a heap, we select them one-by-one.
	// Each time we select a root beyond the first, we check to see if an
	// already-selected root's landing area overlaps.  If so, we take away
	// any benefit associated with the bases/qualities in the landing area
	// and then push it back onto the heap if that changes its priority.
	while(roots.size() < sizeTarget) {
		if(rootHeap_.empty()) {
			break;
		}
		DescentRoot r = rootHeap_.pop();
		final double off = r.fw ? r.off5p : (qlen - r.off5p - 1);
		int fwi = r.fw ? 0 : 1;
		// Re-calculate priority
		int pri = 0;
		if(r.l2r) {
			for(double i = 0; i < landing_; i++) {
				pri += scores_[fwi][off + i];
			}
		} else {
			for(double i = 0; i < landing_; i++) {
				pri += scores_[fwi][off - i];
			}
		}
		// Must take end bonus into account when re-calculating
		if((r.l2r && (off == 0)) || (!r.l2r && (off == qlen - 1))) {
			pri += endBonus;
		}
		if(pri == r.pri) {
			// Update the positions in this root's landing area
			if(r.l2r) {
				for(double i = 0; i < landing_; i++) {
					float frac = ((float)i / (float)landing_);
					scores_[fwi][off + i] = (int)(scores_[fwi][off + i] * frac);
				}
			} else {
				for(double i = 0; i < landing_; i++) {
					float frac = ((float)i / (float)landing_);
					scores_[fwi][off - i] = (int)(scores_[fwi][off - i] * frac);
				}
			}
			confs.expand();
			confs.back().cons.init(landing_, consExp_);
			roots.push_back(r);
		} else {
			// Re-insert the root, its priority now changed
			r.pri = pri;
			rootHeap_.insert(r);
		}
	}
    }
  }
}
