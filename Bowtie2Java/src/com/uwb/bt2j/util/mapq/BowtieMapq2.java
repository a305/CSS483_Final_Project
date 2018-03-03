package com.uwb.bt2j.util.mapq;

import com.uwb.bt2j.aligner.AlignerFlags;
import com.uwb.bt2j.aligner.AlignerResult;
import com.uwb.bt2j.aligner.AlignmentSetSumm;
import com.uwb.bt2j.aligner.Scoring;
import com.uwb.bt2j.aligner.SimpleFunc;

public class BowtieMapq2 extends Mapq{

	public BowtieMapq2(SimpleFunc scoreMin, Scoring sc) {
		super(scoreMin, sc);
	}

	@Override
	public long mapq(AlignmentSetSumm s, AlignerFlags flags, boolean mate1, long rdlen, long ordlen, String inps) {
		// Did the read have a second-best alignment?
				boolean hasSecbest = s.paired() ?
					AlignerResult.VALID_AL_SCORE(s.bestUnchosenCScore()) :
					AlignerResult.VALID_AL_SCORE(s.bestUnchosenScore(mate1));
				// This corresponds to a scenario where we found one and only one
				// alignment but didn't really look for a second one
				if(!flags.canMax() && !s.exhausted(mate1) && !hasSecbest) {
					return 255;
				}
				// scPer = score of a perfect match
				long scPer = (long)sc_.perfectScore(rdlen);
				if(s.paired()) {
					scPer += (long)sc_.perfectScore(ordlen);
				}
				// scMin = score of a just barely valid match
				long scMin = scoreMin_.f<long>((float)rdlen);
				if(s.paired()) {
					scMin += scoreMin_.f<long>((float)ordlen);
				}
				long secbest = scMin-1;
				long diff = (scPer - scMin);  // scores can vary by up to this much
				long ret = 0;
				long best = s.paired() ?
					s.bestCScore().score() : s.bestScore(mate1).score();
				// best score but normalized so that 0 = worst valid score
				long bestOver = best - scMin;
				if(sc_.monotone) {
					// End-to-end alignment
					if(!hasSecbest) {
						if     (bestOver >= diff * (double)0.8f) ret = 42;
						else if(bestOver >= diff * (double)0.7f) ret = 40;
						else if(bestOver >= diff * (double)0.6f) ret = 24;
						else if(bestOver >= diff * (double)0.5f) ret = 23;
						else if(bestOver >= diff * (double)0.4f) ret = 8;
						else if(bestOver >= diff * (double)0.3f) ret = 3;
						else                                     ret = 0;
					} else {
						secbest = s.paired() ?
							s.bestUnchosenCScore().score() : s.bestUnchosenScore(mate1).score();
						long bestdiff = abs(abs(static_cast<long>(best))-abs(static_cast<long>(secbest)));
						if(bestdiff >= diff * (double)0.9f) {
							if(bestOver == diff) {
								ret = 39;
							} else {
								ret = 33;
							}
						} else if(bestdiff >= diff * (double)0.8f) {
							if(bestOver == diff) {
								ret = 38;
							} else {
								ret = 27;
							}
						} else if(bestdiff >= diff * (double)0.7f) {
							if(bestOver == diff) {
								ret = 37;
							} else {
								ret = 26;
							}
						} else if(bestdiff >= diff * (double)0.6f) {
							if(bestOver == diff) {
								ret = 36;
							} else {
								ret = 22;
							}
						} else if(bestdiff >= diff * (double)0.5f) {
							// Top third is still pretty good
							if       (bestOver == diff) {
								ret = 35;
							} else if(bestOver >= diff * (double)0.84f) {
								ret = 25;
							} else if(bestOver >= diff * (double)0.68f) {
								ret = 16;
							} else {
								ret = 5;
							}
						} else if(bestdiff >= diff * (double)0.4f) {
							// Top third is still pretty good
							if       (bestOver == diff) {
								ret = 34;
							} else if(bestOver >= diff * (double)0.84f) {
								ret = 21;
							} else if(bestOver >= diff * (double)0.68f) {
								ret = 14;
							} else {
								ret = 4;
							}
						} else if(bestdiff >= diff * (double)0.3f) {
							// Top third is still pretty good
							if       (bestOver == diff) {
								ret = 32;
							} else if(bestOver >= diff * (double)0.88f) {
								ret = 18;
							} else if(bestOver >= diff * (double)0.67f) {
								ret = 15;
							} else {
								ret = 3;
							}
						} else if(bestdiff >= diff * (double)0.2f) {
							// Top third is still pretty good
							if       (bestOver == diff) {
								ret = 31;
							} else if(bestOver >= diff * (double)0.88f) {
								ret = 17;
							} else if(bestOver >= diff * (double)0.67f) {
								ret = 11;
							} else {
								ret = 0;
							}
						} else if(bestdiff >= diff * (double)0.1f) {
							// Top third is still pretty good
							if       (bestOver == diff) {
								ret = 30;
							} else if(bestOver >= diff * (double)0.88f) {
								ret = 12;
							} else if(bestOver >= diff * (double)0.67f) {
								ret = 7;
							} else {
								ret = 0;
							}
						} else if(bestdiff > 0) {
							// Top third is still pretty good
							if(bestOver >= diff * (double)0.67f) {
								ret = 6;
							} else {
								ret = 2;
							}
						} else {
							// Top third is still pretty good
							if(bestOver >= diff * (double)0.67f) {
								ret = 1;
							} else {
								ret = 0;
							}
						}
					}
				} else {
					// Local alignment
					if(!hasSecbest) {
						if     (bestOver >= diff * (double)0.8f) ret = 44;
						else if(bestOver >= diff * (double)0.7f) ret = 42;
						else if(bestOver >= diff * (double)0.6f) ret = 41;
						else if(bestOver >= diff * (double)0.5f) ret = 36;
						else if(bestOver >= diff * (double)0.4f) ret = 28;
						else if(bestOver >= diff * (double)0.3f) ret = 24;
						else                                     ret = 22;
					} else {
						secbest = s.paired() ?
							s.bestUnchosenCScore().score() : s.bestUnchosenScore(mate1).score();
						long bestdiff = abs(abs(static_cast<long>(best))-abs(static_cast<long>(secbest)));
						if     (bestdiff >= diff * (double)0.9f) ret = 40;
						else if(bestdiff >= diff * (double)0.8f) ret = 39;
						else if(bestdiff >= diff * (double)0.7f) ret = 38;
						else if(bestdiff >= diff * (double)0.6f) ret = 37;
						else if(bestdiff >= diff * (double)0.5f) {
							if     (bestOver == diff)                 ret = 35;
							else if(bestOver >= diff * (double)0.50f) ret = 25;
							else                                      ret = 20;
						} else if(bestdiff >= diff * (double)0.4f) {
							if     (bestOver == diff)                 ret = 34;
							else if(bestOver >= diff * (double)0.50f) ret = 21;
							else                                      ret = 19;
						} else if(bestdiff >= diff * (double)0.3f) {
							if     (bestOver == diff)                 ret = 33;
							else if(bestOver >= diff * (double)0.5f)  ret = 18;
							else                                      ret = 16;
						} else if(bestdiff >= diff * (double)0.2f) {
							if     (bestOver == diff)                 ret = 32;
							else if(bestOver >= diff * (double)0.5f)  ret = 17;
							else                                      ret = 12;
						} else if(bestdiff >= diff * (double)0.1f) {
							if     (bestOver == diff)                 ret = 31;
							else if(bestOver >= diff * (double)0.5f)  ret = 14;
							else                                      ret = 9;
						} else if(bestdiff > 0) {
							if(bestOver >= diff * (double)0.5f)       ret = 11;
							else                                      ret = 2;
						} else {
							if(bestOver >= diff * (double)0.5f)       ret = 1;
							else                                      ret = 0;
						}
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