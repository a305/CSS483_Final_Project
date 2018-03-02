package com.uwb.bt2j.util.mapq;

import com.uwb.bt2j.aligner.AlignerFlags;
import com.uwb.bt2j.aligner.AlignerResult;
import com.uwb.bt2j.aligner.AlignmentSetSumm;
import com.uwb.bt2j.aligner.Scoring;
import com.uwb.bt2j.aligner.SimpleFunc;

public abstract class Mapq {
	public static final long unp_nosec_perf = 44;
	public static final long unp_nosec[] = {
			43, 42, 41, 36, 32, 27, 20, 11, 4, 1, 0
	};
	public static final long unp_sec_perf[] = {
			2, 16, 23, 30, 31, 32, 34, 36, 38, 40, 42
	};
	public static final long unp_sec[][] = {
			{  2,  2,  2,  1,  1, 0, 0, 0, 0, 0, 0},
			{ 20, 14,  7,  3,  2, 1, 0, 0, 0, 0, 0},
			{ 20, 16, 10,  6,  3, 1, 0, 0, 0, 0, 0},
			{ 20, 17, 13,  9,  3, 1, 1, 0, 0, 0, 0},
			{ 21, 19, 15,  9,  5, 2, 2, 0, 0, 0, 0},
			{ 22, 21, 16, 11, 10, 5, 0, 0, 0, 0, 0},
			{ 23, 22, 19, 16, 11, 0, 0, 0, 0, 0, 0},
			{ 24, 25, 21, 30,  0, 0, 0, 0, 0, 0, 0},
			{ 30, 26, 29,  0,  0, 0, 0, 0, 0, 0, 0},
			{ 30, 27,  0,  0,  0, 0, 0, 0, 0, 0, 0},
			{ 30,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0}
	};
	public static final long pair_nosec_perf = 44;
	protected SimpleFunc scoreMin_;
	protected Scoring sc_;
	
	public Mapq(SimpleFunc scoreMin, Scoring sc) {
		scoreMin_ = scoreMin;
		sc_ = sc;
	}
	
	public abstract long mapq(
			AlignmentSetSumm s,
			AlignerFlags flags,
			boolean mate1,
			long rdlen,
			long ordlen,
			String inps
			);
	
	public static boolean bestIsUnique(
			AlignmentSetSumm s,
			AlignerFlags flags,
			boolean mate1,
			long rdlen,
			long ordlen,
			String inps
			) {
		return !AlignerResult.VALID_AL_SCORE(s.bestUnchosenScore(mate1));
	}
	
	public static Mapq new_mapq(
			int version,
			SimpleFunc scoreMin,
			Scoring sc
			) {
		if(version == 3) {
			return new BowtieMapq3(scoreMin, sc);
		} else if(version == 2) {
			return new BowtieMapq2(scoreMin, sc);
		} else {
			return new BowtieMapq(scoreMin, sc);
		}
	}
}
