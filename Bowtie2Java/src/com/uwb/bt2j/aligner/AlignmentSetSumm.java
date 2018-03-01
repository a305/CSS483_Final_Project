package com.uwb.bt2j.aligner;

public class AlignmentSetSumm {
	protected long other1_;        // # more alignments within N points of second-best
	protected long other2_;        // # more alignments within N points of second-best
	protected Boolean     paired_;        // results are paired
	protected Boolean     exhausted1_;    // searched exhaustively for mate 1 alignments?
	protected Boolean     exhausted2_;    // searched exhaustively for mate 2 alignments?
	protected long orefid_;
	protected long orefoff_;

	protected AlignmentScore bestUScore_;
	protected AlignmentScore bestUDist_;
	protected AlignmentScore bestP1Score_;
	protected AlignmentScore bestP1Dist_;
	protected AlignmentScore bestP2Score_;
	protected AlignmentScore bestP2Dist_;
	protected AlignmentScore bestCScore_;
	protected AlignmentScore bestCDist_;

	protected AlignmentScore bestUnchosenUScore_;
	protected AlignmentScore bestUnchosenUDist_;
	protected AlignmentScore bestUnchosenP1Score_;
	protected AlignmentScore bestUnchosenP1Dist_;
	protected AlignmentScore bestUnchosenP2Score_;
	protected AlignmentScore bestUnchosenP2Dist_;
	protected AlignmentScore bestUnchosenCScore_;
	protected AlignmentScore bestUnchosenCDist_;
	
	public AlignmentSetSumm() {
		reset();
	}
	
	public AlignmentSetSumm(Read rd1, Read rd2,
			EList<AlignmentResult> rs1,
			EList<AlignmentResult> rs2,
			EList<AlignmentResult> rs1u,
			EList<AlignmentResult> rs2u,
			Boolean exhausted1, Boolean exhausted2,
			long orefid,
			long orefoff) {
		init(rd1, rd2, rs1, rs2, rs1u, rs2u, exhausted1, exhausted2, 
			     orefid, orefoff);
	}
	
	public AlignmentSetSumm(long other1, long other2,
			Boolean paired, Boolean exhausted1,
			Boolean exhausted2, long orefid,
			long orefoff) {
		init(
				other1,
				other2,
				paired,
				exhausted1,
				exhausted2,
				orefid,
				orefoff);
	}
	
	public void reset() {
		bestUScore_.invalidate();
		bestP1Score_.invalidate();
		bestP2Score_.invalidate();
		bestCScore_.invalidate();
		bestUDist_.invalidate();
		bestP1Dist_.invalidate();
		bestP2Dist_.invalidate();
		bestCDist_.invalidate();
		bestUnchosenUScore_.invalidate();
		bestUnchosenP1Score_.invalidate();
		bestUnchosenP2Score_.invalidate();
		bestUnchosenCScore_.invalidate();
		bestUnchosenUDist_.invalidate();
		bestUnchosenP1Dist_.invalidate();
		bestUnchosenP2Dist_.invalidate();
		bestUnchosenCDist_.invalidate();
		other1_ = other2_ = 0;
		paired_ = false;
		exhausted1_ = exhausted2_ = false;
		orefid_ = -1;
		orefoff_ = -1;
	}
	
	public void init(Read rd1, Read rd2,
			EList<AlignmentResult> rs1,
			EList<AlignmentResult> rs2,
			EList<AlignmentResult> rs1u,
			EList<AlignmentResult> rs2u,
			Boolean exhausted1, Boolean exhausted2,
			long orefid,
			long orefoff) {
		
	}
	
	public void init(long other1, long other2,
			Boolean paired, Boolean exhausted1,
			Boolean exhausted2, long orefid,
			long orefoff) {
		other1_        = other1;
		other2_        = other2;
		paired_        = paired;
		exhausted1_    = exhausted1;
		exhausted2_    = exhausted2;
		orefid_        = orefid;
		orefoff_       = orefoff;
	}
	
	public Boolean empty() {
		return !AlignerResult.VALID_AL_SCORE(bestScore(true));
	}
	
	public long other1()         { return other1_;        }
	public long other2()         { return other2_;        }
	public Boolean     paired()         { return paired_;        }
	public Boolean     exhausted1()     { return exhausted1_;    }
	public Boolean     exhausted2()     { return exhausted2_;    }
	public long   orefid()         { return orefid_;        }
	public long  orefoff()        { return orefoff_;       }
	
	public AlignmentScore bestUScore()   { return bestUScore_;  }
	public AlignmentScore bestP1Score()  { return bestP1Score_; }
	public AlignmentScore bestP2Score()  { return bestP2Score_; }
	public AlignmentScore bestCScore()   { return bestCScore_;  }
	public AlignmentScore bestUDist()    { return bestUDist_;  }
	public AlignmentScore bestP1Dist()   { return bestP1Dist_; }
	public AlignmentScore bestP2Dist()   { return bestP2Dist_; }
	public AlignmentScore bestCDist()    { return bestCDist_;  }

	public AlignmentScore bestUnchosenUScore()   { return bestUnchosenUScore_;  }
	public AlignmentScore bestUnchosenP1Score()  { return bestUnchosenP1Score_; }
	public AlignmentScore bestUnchosenP2Score()  { return bestUnchosenP2Score_; }
	public AlignmentScore bestUnchosenCScore()   { return bestUnchosenCScore_;  }
	public AlignmentScore bestUnchosenUDist()    { return bestUnchosenUDist_;  }
	public AlignmentScore bestUnchosenP1Dist()   { return bestUnchosenP1Dist_; }
	public AlignmentScore bestUnchosenP2Dist()   { return bestUnchosenP2Dist_; }
	public AlignmentScore bestUnchosenCDist()    { return bestUnchosenCDist_;  }
	
	/**
	 * Return best unchosen alignment score for end 1 or 2 of a pair.
	 */
	public AlignmentScore bestUnchosenPScore(Boolean mate1)  {
		return mate1 ? bestUnchosenP1Score_ : bestUnchosenP2Score_;
	}

	/**
	 * Return best unchosen edit distance for end 1 or 2 of a pair.
	 */
	public AlignmentScore bestUnchosenPDist(Boolean mate1)  {
		return mate1 ? bestUnchosenP1Dist_ : bestUnchosenP2Dist_;
	}
	
	/**
	 * Return best unchosen alignment score for end 1 or 2 whether
	 * the read is a pair or not.
	 */
	public AlignmentScore bestUnchosenScore(Boolean mate1)  {
		return paired_ ? (mate1 ? bestUnchosenP1Score_ : bestUnchosenP2Score_) : bestUnchosenUScore();
	}

	/**
	 * Return best unchosen edit distance for end 1 or 2 whether
	 * the read is a pair or not.
	 */
	public AlignmentScore bestUnchosenDist(Boolean mate1)  {
		return paired_ ? (mate1 ? bestUnchosenP1Dist_ : bestUnchosenP2Dist_) : bestUnchosenUDist();
	}

	public Boolean exhausted(Boolean mate1)  {
		return mate1 ? exhausted1_ : exhausted2_;
	}
	
	/**
	 * Return best alignment score for end 1 or 2 whether the read is
	 * a pair or not.
	 */
	public AlignmentScore bestScore(Boolean mate1)  {
		return paired_ ? (mate1 ? bestP1Score_ : bestP2Score_) : bestUScore_;
	}

	/**
	 * Return best edit distance for end 1 or 2 whether the read is
	 * a pair or not.
	 */
	public AlignmentScore bestDist(Boolean mate1)  {
		return paired_ ? (mate1 ? bestP1Dist_ : bestP2Dist_) : bestUDist_;
	}
	
	public void setBest(
			AlignmentScore bestUScore,
			AlignmentScore bestUDist,
			AlignmentScore bestP1Score,
			AlignmentScore bestP1Dist,
			AlignmentScore bestP2Score,
			AlignmentScore bestP2Dist,
			AlignmentScore bestCScore,
			AlignmentScore bestCDist,
			AlignmentScore bestUnchosenUScore,
			AlignmentScore bestUnchosenUDist,
			AlignmentScore bestUnchosenP1Score,
			AlignmentScore bestUnchosenP1Dist,
			AlignmentScore bestUnchosenP2Score,
			AlignmentScore bestUnchosenP2Dist,
			AlignmentScore bestUnchosenCScore,
			AlignmentScore bestUnchosenCDist)
		{
			bestUScore_ = bestUScore;
			bestUDist_ = bestUDist;
			bestP1Score_ = bestP1Score;
			bestP1Dist_ = bestP1Dist;
			bestP2Score_ = bestP2Score;
			bestP2Dist_ = bestP2Dist;
			bestCScore_ = bestCScore;
			bestCDist_ = bestCDist;
			bestUnchosenUScore_ = bestUnchosenUScore;
			bestUnchosenUDist_ = bestUnchosenUDist;
			bestUnchosenP1Score_ = bestUnchosenP1Score;
			bestUnchosenP1Dist_ = bestUnchosenP1Dist;
			bestUnchosenP2Score_ = bestUnchosenP2Score;
			bestUnchosenP2Dist_ = bestUnchosenP2Dist;
			bestUnchosenCScore_ = bestUnchosenCScore;
			bestUnchosenCDist_ = bestUnchosenCDist;
		}
}
