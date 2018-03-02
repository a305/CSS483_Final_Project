package com.uwb.bt2j.util.mapq;

import com.uwb.bt2j.aligner.AlignerFlags;
import com.uwb.bt2j.aligner.AlignerResult;
import com.uwb.bt2j.aligner.AlignmentSetSumm;
import com.uwb.bt2j.aligner.Scoring;
import com.uwb.bt2j.aligner.SimpleFunc;

public class BowtieMapq extends Mapq{

	public BowtieMapq(SimpleFunc scoreMin, Scoring sc) {
		super(scoreMin, sc);
	}

	@Override
	public long mapq(AlignmentSetSumm s, AlignerFlags flags, boolean mate1, long rdlen, long ordlen, String inps) {
		boolean hasSecbest = AlignerResult.VALID_AL_SCORE(s.bestUnchosenScore(mate1));
		if(!flags.canMax() && !s.exhausted(mate1) && !hasSecbest) {
			return 255;
		}
		long scPer = (long)sc_.perfectScore(rdlen);
		long scMin = scoreMin_.f<long>((float)rdlen);
		long secbest = scMin-1;
		long diff = (scPer - scMin);
		float sixth_2 = (float)(scPer - diff * (double)0.1666f * 2); 
		float sixth_3 = (float)(scPer - diff * (double)0.1666f * 3);
		long ret = 0;
		long best = s.bestScore(mate1).score();
		if(!hasSecbest) {
			// Top third?
			if(best >= sixth_2) {
				ret = 37;
			}
			// Otherwise in top half?
			else if(best >= sixth_3) {
				ret = 25;
			}
			// Otherwise has no second-best?
			else {
				ret = 10;
			}
		} else {
			secbest = s.bestUnchosenScore(mate1).score();
			long bestdiff = abs(abs(static_cast<long>(best))-abs(static_cast<long>(secbest)));
			if(bestdiff >= diff * 0.1666 * 5) {
				ret = 6;
			} else if(bestdiff >= diff * 0.1666 * 4) {
				ret = 5;
			} else if(bestdiff >= diff * 0.1666 * 3) {
				ret = 4;
			} else if(bestdiff >= diff * 0.1666 * 2) {
				ret = 3;
			} else if(bestdiff >= diff * 0.1666 * 1) {
				ret = 2;
			} else {
				ret = 1;
			}
		}
		// Note: modifications to inps must be synchronized
		//if(inps != NULL) {
		//	inps = itoa10<long>(best, inps);
		//	*inps++ = ',';
		//	inps = itoa10<long>(secbest, inps);
		//	*inps++ = ',';
		//	inps = itoa10<TMapq>(ret, inps);
		//}
		return ret;
	}

}
