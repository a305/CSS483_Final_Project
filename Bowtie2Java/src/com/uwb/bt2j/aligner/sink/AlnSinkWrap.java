package com.uwb.bt2j.aligner.sink;

import com.uwb.bt2j.aligner.PerReadMetrics;
import com.uwb.bt2j.aligner.RandomSource;
import com.uwb.bt2j.aligner.Read;
import com.uwb.bt2j.aligner.Scoring;
import com.uwb.bt2j.aligner.seed.SeedResults;

public class AlnSinkWrap {
	protected AlnSink g_;
	protected ReportingParams rp_;
	protected double threadid_;
	protected Mapq mapq_;
	protected boolean init_;
	protected boolean maxed1_;
	protected boolean maxed2_;
	protected boolean maxedOverall_;
	protected long bestPair_;
	protected long best2Pair_;
	protected long bestUnp1;
	protected long best2Unp1;
	protected long bestUnp2;
	protected long best2Unp2;
	protected final Read rd1_;
	protected final Read rd2_;
	protected double rdid_;
	protected EList<AlignmentResult> rs1_;
	protected EList<AlignmentResult> rs2_;
	protected EList<AlignmentResult> rs1u_;
	protected EList<AlignmentResult> rs2u_;
	protected EList<double> select1_;
	protected EList<double> select2_;
	protected ReportingState st_;
	protected EList<Pair<AlnScore, double>> selectBuf_;
	protected BString obuf_;
	protected StackedAln staln_;

	public int nextRead(
			Read rd1,      // new mate #1
			Read rd2,      // new mate #2
			long rdid,         // read ID for new pair
			boolean qualitiesMatter) // aln policy distinguishes b/t quals?
	{
		init_ = true;
		// Keep copy of new read, so that we can compare it with the
		// next one
		if(rd1 != null) {
			rd1_ = rd1;
		} else rd1_ = null;
		if(rd2 != null) {
			rd2_ = rd2;
		} else rd2_ = null;
		rdid_ = rdid;
		// Caller must now align the read
		maxed1_ = false;
		maxed2_ = false;
		maxedOverall_ = false;
		bestPair_ = best2Pair_ =
		bestUnp1_ = best2Unp1_ =
		bestUnp2_ = best2Unp2_ = Long.MIN_VALUE;
		rs1_.clear();     // clear out paired-end alignments
		rs2_.clear();     // clear out paired-end alignments
		rs1u_.clear();    // clear out unpaired alignments for mate #1
		rs2u_.clear();    // clear out unpaired alignments for mate #2
		st_.nextRead(readIsPair()); // reset state
		// Start from the first stage
		return 0;
	}
	
	public void finishRead(
			SeedResults sr1,         // seed alignment results for mate 1
			SeedResults sr2,         // seed alignment results for mate 2
			boolean               exhaust1,    // mate 1 exhausted?
			boolean               exhaust2,    // mate 2 exhausted?
			boolean               nfilt1,      // mate 1 N-filtered?
			boolean               nfilt2,      // mate 2 N-filtered?
			boolean               scfilt1,     // mate 1 score-filtered?
			boolean               scfilt2,     // mate 2 score-filtered?
			boolean               lenfilt1,    // mate 1 length-filtered?
			boolean               lenfilt2,    // mate 2 length-filtered?
			boolean               qcfilt1,     // mate 1 qc-filtered?
			boolean               qcfilt2,     // mate 2 qc-filtered?
			RandomSource      rnd,         // pseudo-random generator
			ReportingMetrics  met,         // reporting metrics
			PerReadMetrics prm,      // per-read metrics
			Scoring sc,              // scoring scheme
			boolean suppressSeedSummary,       // = true
			boolean suppressAlignments,        // = false
			boolean scUnMapped,                // = false
			boolean xeq)                       // = false
	{
		obuf_.clear();
		OutputQueueMark qqm(g_.outq(), obuf_, rdid_, threadid_);
		assert(init_);
		if(!suppressSeedSummary) {
			if(sr1 != null) {
				assert(rd1_ != null);
				// Mate exists and has non-empty SeedResults
				g_.reportSeedSummary(obuf_, *rd1_, rdid_, threadid_, *sr1, true);
			} else if(rd1_ != null) {
				// Mate exists but has null SeedResults
				g_.reportEmptySeedSummary(obuf_, *rd1_, rdid_, true);
			}
			if(sr2 != null) {
				assert(rd2_ != null);
				// Mate exists and has non-empty SeedResults
				g_.reportSeedSummary(obuf_, *rd2_, rdid_, threadid_, *sr2, true);
			} else if(rd2_ != null) {
				// Mate exists but has null SeedResults
				g_.reportEmptySeedSummary(obuf_, *rd2_, rdid_, true);
			}
		}
		if(!suppressAlignments) {
			// Ask the ReportingState what to report
			st_.finish();
			uint64_t nconcord = 0, ndiscord = 0, nunpair1 = 0, nunpair2 = 0;
			bool pairMax = false, unpair1Max = false, unpair2Max = false;
			st_.getReport(
				nconcord,
				ndiscord,
				nunpair1,
				nunpair2,
				pairMax,
				unpair1Max,
				unpair2Max);
			assert_leq(nconcord, rs1_.size());
			assert_leq(nunpair1, rs1u_.size());
			assert_leq(nunpair2, rs2u_.size());
			assert_leq(ndiscord, 1);
			assert_gt(rp_.khits, 0);
			assert_gt(rp_.mhits, 0);
			assert(!pairMax    || rs1_.size()  >= (uint64_t)rp_.mhits);
			assert(!unpair1Max || rs1u_.size() >= (uint64_t)rp_.mhits);
			assert(!unpair2Max || rs2u_.size() >= (uint64_t)rp_.mhits);
			met.nread++;
			if(readIsPair()) {
				met.npaired++;
			} else {
				met.nunpaired++;
			}
			// Report concordant paired-end alignments if possible
			if(nconcord > 0) {
				AlnSetSumm concordSumm(
					rd1_, rd2_, &rs1_, &rs2_, &rs1u_, &rs2u_,
					exhaust1, exhaust2, -1, -1);
				// Sort by score then pick from low to high
				AlnScore bestUScore, bestP1Score, bestP2Score, bestCScore;
				AlnScore bestUDist, bestP1Dist, bestP2Dist, bestCDist;
				AlnScore bestUnchosenUScore, bestUnchosenP1Score, bestUnchosenP2Score, bestUnchosenCScore;
				AlnScore bestUnchosenUDist, bestUnchosenP1Dist, bestUnchosenP2Dist, bestUnchosenCDist;
				// TODO: should probably package these variables up so it's not
				// such a pain to pass them around
				size_t off = selectByScore(
					&rs1_, &rs2_,
					nconcord, select1_,
					&rs1u_, &rs2u_,
					bestUScore,
					bestUDist,
					bestP1Score,
					bestP1Dist,
					bestP2Score,
					bestP2Dist,
					bestCScore,
					bestCDist,
					bestUnchosenUScore,
					bestUnchosenUDist,
					bestUnchosenP1Score,
					bestUnchosenP1Dist,
					bestUnchosenP2Score,
					bestUnchosenP2Dist,
					bestUnchosenCScore,
					bestUnchosenCDist,
					rnd);
				concordSumm.setBest(
					bestUScore,
					bestUDist,
					bestP1Score,
					bestP1Dist,
					bestP2Score,
					bestP2Dist,
					bestCScore,
					bestCDist,
					bestUnchosenUScore,
					bestUnchosenUDist,
					bestUnchosenP1Score,
					bestUnchosenP1Dist,
					bestUnchosenP2Score,
					bestUnchosenP2Dist,
					bestUnchosenCScore,
					bestUnchosenCDist);
				assert(concordSumm.bestScore(true).valid());
				assert(concordSumm.bestScore(false).valid());
				assert_lt(off, rs1_.size());
				const AlignmentResult *rs1 = &rs1_[off];
				const AlignmentResult *rs2 = &rs2_[off];
				AlnFlags flags1(
					ALN_FLAG_PAIR_CONCORD_MATE1,
					st_.params().mhitsSet(),
					unpair1Max,
					pairMax,
					nfilt1,
					scfilt1,
					lenfilt1,
					qcfilt1,
					st_.params().mixed,
					true,       // primary
					true,       // opp aligned
					rs2->fw(),  // opp fw
					scUnMapped,
					xeq);
				AlnFlags flags2(
					ALN_FLAG_PAIR_CONCORD_MATE2,
					st_.params().mhitsSet(),
					unpair2Max,
					pairMax,
					nfilt2,
					scfilt2,
					lenfilt2,
					qcfilt2,
					st_.params().mixed,
					false,      // primary
					true,       // opp aligned
					rs1->fw(),  // opp fw
					scUnMapped,
					xeq);
				// Issue: we only set the flags once, but some of the flags might
				// vary from pair to pair among the pairs we're reporting.  For
				// instance, whether the a given mate aligns to the forward strand.
				SeedAlSumm ssm1, ssm2;
				sr1->toSeedAlSumm(ssm1);
				sr2->toSeedAlSumm(ssm2);
				for(size_t i = 0; i < rs1_.size(); i++) {
					rs1_[i].setMateParams(ALN_RES_TYPE_MATE1, &rs2_[i], flags1);
					rs2_[i].setMateParams(ALN_RES_TYPE_MATE2, &rs1_[i], flags2);
					assert_eq(abs(rs1_[i].fragmentLength()), abs(rs2_[i].fragmentLength()));
				}
				assert(!select1_.empty());
				g_.reportHits(
					obuf_,
					staln_,
					threadid_,
					rd1_,
					rd2_,
					rdid_,
					select1_,
					null,
					&rs1_,
					&rs2_,
					pairMax,
					concordSumm,
					ssm1,
					ssm2,
					&flags1,
					&flags2,
					prm,
					mapq_,
					sc,
					false);
				if(pairMax) {
					met.nconcord_rep++;
				} else {
					met.nconcord_uni++;
					assert(!rs1_.empty());
					if(rs1_.size() == 1) {
						met.nconcord_uni1++;
					} else {
						met.nconcord_uni2++;
					}
				}
				init_ = false;
				//g_.outq().finishRead(obuf_, rdid_, threadid_);
				return;
			}
			// Report disconcordant paired-end alignments if possible
			else if(ndiscord > 0) {
				ASSERT_ONLY(bool ret =) prepareDiscordants();
				assert(ret);
				assert_eq(1, rs1_.size());
				assert_eq(1, rs2_.size());
				AlnSetSumm discordSumm(
					rd1_, rd2_, &rs1_, &rs2_, &rs1u_, &rs2u_,
					exhaust1, exhaust2, -1, -1);
				const AlignmentResult *rs1 = &rs1_[0];
				const AlignmentResult *rs2 = &rs2_[0];
				AlnFlags flags1(
					ALN_FLAG_PAIR_DISCORD_MATE1,
					st_.params().mhitsSet(),
					false,
					pairMax,
					nfilt1,
					scfilt1,
					lenfilt1,
					qcfilt1,
					st_.params().mixed,
					true,       // primary
					true,       // opp aligned
					rs2->fw(),  // opp fw
					scUnMapped,
					xeq);
				AlnFlags flags2(
					ALN_FLAG_PAIR_DISCORD_MATE2,
					st_.params().mhitsSet(),
					false,
					pairMax,
					nfilt2,
					scfilt2,
					lenfilt2,
					qcfilt2,
					st_.params().mixed,
					false,      // primary
					true,       // opp aligned
					rs1->fw(),  // opp fw
					scUnMapped,
					xeq);
				SeedAlSumm ssm1, ssm2;
				sr1->toSeedAlSumm(ssm1);
				sr2->toSeedAlSumm(ssm2);
				for(size_t i = 0; i < rs1_.size(); i++) {
					rs1_[i].setMateParams(ALN_RES_TYPE_MATE1, &rs2_[i], flags1);
					rs2_[i].setMateParams(ALN_RES_TYPE_MATE2, &rs1_[i], flags2);
					assert(rs1_[i].isFraglenSet() == rs2_[i].isFraglenSet());
					assert(!rs1_[i].isFraglenSet() || abs(rs1_[i].fragmentLength()) == abs(rs2_[i].fragmentLength()));
				}
				AlnScore bestUScore, bestP1Score, bestP2Score, bestCScore;
				AlnScore bestUDist, bestP1Dist, bestP2Dist, bestCDist;
				AlnScore bestUnchosenUScore, bestUnchosenP1Score, bestUnchosenP2Score, bestUnchosenCScore;
				AlnScore bestUnchosenUDist, bestUnchosenP1Dist, bestUnchosenP2Dist, bestUnchosenCDist;
				ASSERT_ONLY(size_t off =) selectByScore(
					&rs1_, &rs2_,
					ndiscord, select1_,
					&rs1u_, &rs2u_,
					bestUScore,
					bestUDist,
					bestP1Score,
					bestP1Dist,
					bestP2Score,
					bestP2Dist,
					bestCScore,
					bestCDist,
					bestUnchosenUScore,
					bestUnchosenUDist,
					bestUnchosenP1Score,
					bestUnchosenP1Dist,
					bestUnchosenP2Score,
					bestUnchosenP2Dist,
					bestUnchosenCScore,
					bestUnchosenCDist,
					rnd);
				discordSumm.setBest(
					bestUScore,
					bestUDist,
					bestP1Score,
					bestP1Dist,
					bestP2Score,
					bestP2Dist,
					bestCScore,
					bestCDist,
					bestUnchosenUScore,
					bestUnchosenUDist,
					bestUnchosenP1Score,
					bestUnchosenP1Dist,
					bestUnchosenP2Score,
					bestUnchosenP2Dist,
					bestUnchosenCScore,
					bestUnchosenCDist);
				assert_eq(0, off);
				assert(!select1_.empty());
				g_.reportHits(
					obuf_,
					staln_,
					threadid_,
					rd1_,
					rd2_,
					rdid_,
					select1_,
					null,
					&rs1_,
					&rs2_,
					pairMax,
					discordSumm,
					ssm1,
					ssm2,
					&flags1,
					&flags2,
					prm,
					mapq_,
					sc,
					false);
				met.nconcord_0++;
				met.ndiscord++;
				init_ = false;
				//g_.outq().finishRead(obuf_, rdid_, threadid_);
				return;
			}
			// If we're at this point, at least one mate failed to align.
			// BTL: That's not true.  It could be that there are no concordant
			// alignments but both mates have unpaired alignments, with one of
			// the mates having more than one.
			//assert(nunpair1 == 0 || nunpair2 == 0);
			assert(!pairMax);

			// Update counters given that one mate didn't align
			if(readIsPair()) {
				met.nconcord_0++;
			}
			if(rd1_ != null) {
				if(nunpair1 > 0) {
					// Update counters
					if(readIsPair()) {
						if(unpair1Max) met.nunp_0_rep++;
						else {
							met.nunp_0_uni++;
							assert(!rs1u_.empty());
							if(rs1u_.size() == 1) {
								met.nunp_0_uni1++;
							} else {
								met.nunp_0_uni2++;
							}
						}
					} else {
						if(unpair1Max) met.nunp_rep++;
						else {
							met.nunp_uni++;
							assert(!rs1u_.empty());
							if(rs1u_.size() == 1) {
								met.nunp_uni1++;
							} else {
								met.nunp_uni2++;
							}
						}
					}
				} else if(unpair1Max) {
					// Update counters
					if(readIsPair())   met.nunp_0_rep++;
					else               met.nunp_rep++;
				} else {
					// Update counters
					if(readIsPair())   met.nunp_0_0++;
					else               met.nunp_0++;
				}
			}
			if(rd2_ != null) {
				if(nunpair2 > 0) {
					// Update counters
					if(readIsPair()) {
						if(unpair2Max) met.nunp_0_rep++;
						else {
							assert(!rs2u_.empty());
							met.nunp_0_uni++;
							if(rs2u_.size() == 1) {
								met.nunp_0_uni1++;
							} else {
								met.nunp_0_uni2++;
							}
						}
					} else {
						if(unpair2Max) met.nunp_rep++;
						else {
							assert(!rs2u_.empty());
							met.nunp_uni++;
							if(rs2u_.size() == 1) {
								met.nunp_uni1++;
							} else {
								met.nunp_uni2++;
							}
						}
					}
				} else if(unpair2Max) {
					// Update counters
					if(readIsPair())   met.nunp_0_rep++;
					else               met.nunp_rep++;
				} else {
					// Update counters
					if(readIsPair())   met.nunp_0_0++;
					else               met.nunp_0++;
				}
			}
			
			const AlignmentResult *repRs1 = null, *repRs2 = null;
			AlnSetSumm summ1, summ2;
			AlnFlags flags1, flags2;
			TRefId refid = -1; TRefOff refoff = -1;
			bool rep1 = rd1_ != null && nunpair1 > 0;
			bool rep2 = rd2_ != null && nunpair2 > 0;

			// This is the preliminary if statement for mate 1 - here we're
			// gathering some preliminary information, making it possible to call
			// g_.reportHits(...) with information about both mates potentially
			if(rep1) {
				// Mate 1 aligned at least once
				summ1.init(
					rd1_, null, null, null, &rs1u_, null,
					exhaust1, exhaust2, -1, -1);
				// Sort by score then pick from low to high
				AlnScore bestUScore, bestP1Score, bestP2Score, bestCScore;
				AlnScore bestUDist, bestP1Dist, bestP2Dist, bestCDist;
				AlnScore bestUnchosenUScore, bestUnchosenP1Score, bestUnchosenP2Score, bestUnchosenCScore;
				AlnScore bestUnchosenUDist, bestUnchosenP1Dist, bestUnchosenP2Dist, bestUnchosenCDist;
				size_t off = selectByScore(
					&rs1u_, null, nunpair1, select1_, null, null,
					bestUScore,
					bestUDist,
					bestP1Score,
					bestP1Dist,
					bestP2Score,
					bestP2Dist,
					bestCScore,
					bestCDist,
					bestUnchosenUScore,
					bestUnchosenUDist,
					bestUnchosenP1Score,
					bestUnchosenP1Dist,
					bestUnchosenP2Score,
					bestUnchosenP2Dist,
					bestUnchosenCScore,
					bestUnchosenCDist,
					rnd);
				summ1.setBest(
					bestUScore,
					bestUDist,
					bestP1Score,
					bestP1Dist,
					bestP2Score,
					bestP2Dist,
					bestCScore,
					bestCDist,
					bestUnchosenUScore,
					bestUnchosenUDist,
					bestUnchosenP1Score,
					bestUnchosenP1Dist,
					bestUnchosenP2Score,
					bestUnchosenP2Dist,
					bestUnchosenCScore,
					bestUnchosenCDist);
				repRs1 = &rs1u_[off];
			} else if(rd1_ != null) {
				// Mate 1 failed to align - don't do anything yet.  First we want
				// to collect information on mate 2 in case that factors into the
				// summary
				assert(!unpair1Max);
			}
			
			if(rep2) {
				summ2.init(
					null, rd2_, null, null, null, &rs2u_,
					exhaust1, exhaust2, -1, -1);
				// Sort by score then pick from low to high
				AlnScore bestUScore, bestP1Score, bestP2Score, bestCScore;
				AlnScore bestUDist, bestP1Dist, bestP2Dist, bestCDist;
				AlnScore bestUnchosenUScore, bestUnchosenP1Score, bestUnchosenP2Score, bestUnchosenCScore;
				AlnScore bestUnchosenUDist, bestUnchosenP1Dist, bestUnchosenP2Dist, bestUnchosenCDist;
				size_t off = selectByScore(
					&rs2u_, null, nunpair2, select2_, null, null,
					bestUScore,
					bestUDist,
					bestP1Score,
					bestP1Dist,
					bestP2Score,
					bestP2Dist,
					bestCScore,
					bestCDist,
					bestUnchosenUScore,
					bestUnchosenUDist,
					bestUnchosenP1Score,
					bestUnchosenP1Dist,
					bestUnchosenP2Score,
					bestUnchosenP2Dist,
					bestUnchosenCScore,
					bestUnchosenCDist,
					rnd);
				summ2.setBest(
					bestUScore,
					bestUDist,
					bestP1Score,
					bestP1Dist,
					bestP2Score,
					bestP2Dist,
					bestCScore,
					bestCDist,
					bestUnchosenUScore,
					bestUnchosenUDist,
					bestUnchosenP1Score,
					bestUnchosenP1Dist,
					bestUnchosenP2Score,
					bestUnchosenP2Dist,
					bestUnchosenCScore,
					bestUnchosenCDist);
				repRs2 = &rs2u_[off];
			} else if(rd2_ != null) {
				// Mate 2 failed to align - don't do anything yet.  First we want
				// to collect information on mate 1 in case that factors into the
				// summary
				assert(!unpair2Max);
			}

			// Now set up flags
			if(rep1) {
				// Initialize flags.  Note: We want to have information about how
				// the other mate aligned (if it did) at this point
				flags1.init(
					readIsPair() ?
						ALN_FLAG_PAIR_UNPAIRED_MATE1 :
						ALN_FLAG_PAIR_UNPAIRED,
					st_.params().mhitsSet(),
					unpair1Max,
					pairMax,
					nfilt1,
					scfilt1,
					lenfilt1,
					qcfilt1,
					st_.params().mixed,
					true,   // primary
					repRs2 != null,                    // opp aligned
					repRs2 == null || repRs2->fw(),    // opp fw
					scUnMapped,
					xeq);
				for(size_t i = 0; i < rs1u_.size(); i++) {
					rs1u_[i].setMateParams(ALN_RES_TYPE_UNPAIRED_MATE1, null, flags1);
				}
			}
			if(rep2) {
				// Initialize flags.  Note: We want to have information about how
				// the other mate aligned (if it did) at this point
				flags2.init(
					readIsPair() ?
						ALN_FLAG_PAIR_UNPAIRED_MATE2 :
						ALN_FLAG_PAIR_UNPAIRED,
					st_.params().mhitsSet(),
					unpair2Max,
					pairMax,
					nfilt2,
					scfilt2,
					lenfilt2,
					qcfilt2,
					st_.params().mixed,
					true,   // primary
					repRs1 != null,                  // opp aligned
					repRs1 == null || repRs1->fw(),  // opp fw
					scUnMapped,
					xeq);
				for(size_t i = 0; i < rs2u_.size(); i++) {
					rs2u_[i].setMateParams(ALN_RES_TYPE_UNPAIRED_MATE2, null, flags2);
				}
			}
			
			// Now report mate 1
			if(rep1) {
				SeedAlSumm ssm1, ssm2;
				if(sr1 != null) sr1->toSeedAlSumm(ssm1);
				if(sr2 != null) sr2->toSeedAlSumm(ssm2);
				assert(!select1_.empty());
				g_.reportHits(
					obuf_,
					staln_,
					threadid_,
					rd1_,
					repRs2 != null ? rd2_ : null,
					rdid_,
					select1_,
					repRs2 != null ? &select2_ : null,
					&rs1u_,
					repRs2 != null ? &rs2u_ : null,
					unpair1Max,
					summ1,
					ssm1,
					ssm2,
					&flags1,
					repRs2 != null ? &flags2 : null,
					prm,
					mapq_,
					sc,
					false);
				assert_lt(select1_[0], rs1u_.size());
				refid = rs1u_[select1_[0]].refid();
				refoff = rs1u_[select1_[0]].refoff();
			}
			
			// Now report mate 2
			//if(rep2 && !rep1) {
			if(rep2) {
				SeedAlSumm ssm1, ssm2;
				if(sr1 != null) sr1->toSeedAlSumm(ssm1);
				if(sr2 != null) sr2->toSeedAlSumm(ssm2);
				assert(!select2_.empty());
				g_.reportHits(
					obuf_,
					staln_,
					threadid_,
					rd2_,
					repRs1 != null ? rd1_ : null,
					rdid_,
					select2_,
					repRs1 != null ? &select1_ : null,
					&rs2u_,
					repRs1 != null ? &rs1u_ : null,
					unpair2Max,
					summ2,
					ssm1,
					ssm2,
					&flags2,
					repRs1 != null ? &flags1 : null,
					prm,
					mapq_,
					sc,
					false);
				assert_lt(select2_[0], rs2u_.size());
				refid = rs2u_[select2_[0]].refid();
				refoff = rs2u_[select2_[0]].refoff();
			}
			
			if(rd1_ != null && nunpair1 == 0) {
				if(nunpair2 > 0) {
					assert_neq(-1, refid);
					summ1.init(
						rd1_, null, null, null, null, null,
						exhaust1, exhaust2, refid, refoff);
				} else {
					summ1.init(
						rd1_, null, null, null, null, null,
						exhaust1, exhaust2, -1, -1);
				}
				SeedAlSumm ssm1, ssm2;
				if(sr1 != null) sr1->toSeedAlSumm(ssm1);
				if(sr2 != null) sr2->toSeedAlSumm(ssm2);
				flags1.init(
					readIsPair() ?
						ALN_FLAG_PAIR_UNPAIRED_MATE1 :
						ALN_FLAG_PAIR_UNPAIRED,
					st_.params().mhitsSet(),
					false,
					false,
					nfilt1,
					scfilt1,
					lenfilt1,
					qcfilt1,
					st_.params().mixed,
					true,           // primary
					repRs2 != null, // opp aligned
					(repRs2 != null) ? repRs2->fw() : false, // opp fw
					scUnMapped,
					xeq);
				g_.reportUnaligned(
					obuf_,      // string to write output to
					staln_,
					threadid_,
					rd1_,    // read 1
					null,    // read 2
					rdid_,   // read id
					summ1,   // summ
					ssm1,    // 
					ssm2,
					&flags1, // flags 1
					null,    // flags 2
					prm,     // per-read metrics
					mapq_,   // MAPQ calculator
					sc,      // scoring scheme
					true);   // get lock?
			}
			if(rd2_ != null && nunpair2 == 0) {
				if(nunpair1 > 0) {
					assert_neq(-1, refid);
					summ2.init(
						null, rd2_, null, null, null, null,
						exhaust1, exhaust2, refid, refoff);
				} else {
					summ2.init(
						null, rd2_, null, null, null, null,
						exhaust1, exhaust2, -1, -1);
				}
				SeedAlSumm ssm1, ssm2;
				if(sr1 != null) sr1->toSeedAlSumm(ssm1);
				if(sr2 != null) sr2->toSeedAlSumm(ssm2);
				flags2.init(
					readIsPair() ?
						ALN_FLAG_PAIR_UNPAIRED_MATE2 :
						ALN_FLAG_PAIR_UNPAIRED,
					st_.params().mhitsSet(),
					false,
					false,
					nfilt2,
					scfilt2,
					lenfilt2,
					qcfilt2,
					st_.params().mixed,
					true,           // primary
					repRs1 != null, // opp aligned
					(repRs1 != null) ? repRs1->fw() : false, // opp fw
					scUnMapped,
					xeq);
				g_.reportUnaligned(
					obuf_,      // string to write output to
					staln_,
					threadid_,
					rd2_,    // read 1
					null,    // read 2
					rdid_,   // read id
					summ2,   // summ
					ssm1,
					ssm2,
					&flags2, // flags 1
					null,    // flags 2
					prm,     // per-read metrics
					mapq_,   // MAPQ calculator
					sc,      // scoring scheme
					true);   // get lock?
			}
		} // if(suppress alignments)
		init_ = false;
		return;
	}
}
	
	public boolean report(int stage, AlignmentResult rs1, AlignmentResult rs2) {
		bool paired = (rs1 != null && rs2 != null);
		bool one = (rs1 != null);
		const AlignmentResult* rsa = one ? rs1 : rs2;
		const AlignmentResult* rsb = one ? rs2 : rs1;
		if(paired) {
			assert(readIsPair());
			st_.foundConcordant();
			rs1_.push_back(*rs1);
			rs2_.push_back(*rs2);
		} else {
			st_.foundUnpaired(one);
			if(one) {
				rs1u_.push_back(*rs1);
			} else {
				rs2u_.push_back(*rs2);
			}
		}
		// Tally overall alignment score
		TAlScore score = rsa->score().score();
		if(rsb != null) score += rsb->score().score();
		// Update best score so far
		if(paired) {
			if(score > bestPair_) {
				best2Pair_ = bestPair_;
				bestPair_ = score;
			} else if(score > best2Pair_) {
				best2Pair_ = score;
			}
		} else {
			if(one) {
				if(score > bestUnp1_) {
					best2Unp1_ = bestUnp1_;
					bestUnp1_ = score;
				} else if(score > best2Unp1_) {
					best2Unp1_ = score;
				}
			} else {
				if(score > bestUnp2_) {
					best2Unp2_ = bestUnp2_;
					bestUnp2_ = score;
				} else if(score > best2Unp2_) {
					best2Unp2_ = score;
				}
			}
		}
		return st_.done();
	}
	
	public void appendSeedSummary(
			BTString     o,
			Read   rd,
			TReadId rdid,
			size_t        seedsTried,
			size_t        nonzero,
			size_t        ranges,
			size_t        elts,
			size_t        seedsTriedFw,
			size_t        nonzeroFw,
			size_t        rangesFw,
			size_t        eltsFw,
			size_t        seedsTriedRc,
			size_t        nonzeroRc,
			size_t        rangesRc,
			size_t        eltsRc
			) {
		char buf[1024];
		bool firstfield = true;
		//
		// Read name
		//
		BEGIN_FIELD;
		printUptoWs(o, rd.name, true);
		
		//
		// Total number of seeds tried
		//
		BEGIN_FIELD;
		WRITE_NUM(o, seedsTried);

		//
		// Total number of seeds tried where at least one range was found.
		//
		BEGIN_FIELD;
		WRITE_NUM(o, nonzero);

		//
		// Total number of ranges found
		//
		BEGIN_FIELD;
		WRITE_NUM(o, ranges);

		//
		// Total number of elements found
		//
		BEGIN_FIELD;
		WRITE_NUM(o, elts);
		
		//
		// The same four numbers, but only for seeds extracted from the
		// forward read representation.
		//
		BEGIN_FIELD;
		WRITE_NUM(o, seedsTriedFw);

		BEGIN_FIELD;
		WRITE_NUM(o, nonzeroFw);

		BEGIN_FIELD;
		WRITE_NUM(o, rangesFw);

		BEGIN_FIELD;
		WRITE_NUM(o, eltsFw);

		//
		// The same four numbers, but only for seeds extracted from the
		// reverse complement read representation.
		//
		BEGIN_FIELD;
		WRITE_NUM(o, seedsTriedRc);

		BEGIN_FIELD;
		WRITE_NUM(o, nonzeroRc);

		BEGIN_FIELD;
		WRITE_NUM(o, rangesRc);

		BEGIN_FIELD;
		WRITE_NUM(o, eltsRc);

		o.append('\n');
	}
	
	public void reportEmptySeedSummary(
			BTString          o,
			Read        rd,
			long            rdid,
			int             threadId,
			boolean               getLock)
	{
		appendSeedSummary(
				o,                     // string to append to
				rd,                    // read
				rdid,                  // read id
				0,                     // # seeds tried
				0,                     // # seeds with non-empty results
				0,                     // # ranges for all seed hits
				0,                     // # elements for all seed hits
				0,                     // # seeds tried from fw read
				0,                     // # seeds with non-empty results from fw read
				0,                     // # ranges for seed hits from fw read
				0,                     // # elements for seed hits from fw read
				0,                     // # seeds tried from rc read
				0,                     // # seeds with non-empty results from fw read
				0,                     // # ranges for seed hits from fw read
				0);                    // # elements for seed hits from fw read
	}
	
	public void reportSeedSummary(
			BTString          o,
			Read        rd,
			long            rdid,
			int             threadId,
			SeedResults rs,
			boolean               getLock)
	{
		appendSeedSummary(
				o,                     // string to write to
				rd,                    // read
				rdid,                  // read id
				rs.numOffs()*2,        // # seeds tried
				rs.nonzeroOffsets(),   // # seeds with non-empty results
				rs.numRanges(),        // # ranges for all seed hits
				rs.numElts(),          // # elements for all seed hits
				rs.numOffs(),          // # seeds tried from fw read
				rs.nonzeroOffsetsFw(), // # seeds with non-empty results from fw read
				rs.numRangesFw(),      // # ranges for seed hits from fw read
				rs.numEltsFw(),        // # elements for seed hits from fw read
				rs.numOffs(),          // # seeds tried from rc read
				rs.nonzeroOffsetsRc(), // # seeds with non-empty results from fw read
				rs.numRangesRc(),      // # ranges for seed hits from fw read
				rs.numEltsRc());       // # elements for seed hits from fw read
	}
	
	public static void printUptoWs(BTString s, T str, boolean chopws) {
		size_t len = str.length();
		for(size_t i = 0; i < len; i++) {
			if(!chopws || (str[i] != ' ' && str[i] != '\t')) {
				s.append(str[i]);
			} else {
				break;
			}
		}
	}
	
	public boolean empty() {
		return rs1_.empty() && rs1u_.empty() && rs2u_.empty();
	}
	
	public final boolean maxed() {
		return maxedOverall_;
	}
	
	public final boolean readIsPair() {
		return rd1_ != null && rd2_ != null;
	}
	
	public final boolean inited() {
		return init_;
	}
	
	public final ReportingState state() {
		return st_;
	}
	
	public final boolean Mmode() {
		return rp_.mhitsSet();
	}
	
	public final boolean allHits() {
		return rp_.allHits();
	}
	
	public final boolean hasSecondBestUnp1() {
		return best2Unp1_ != Long.MIN_VALUE;
	}
	
	public final boolean hasSecondBestUnp2() {
		return best2Unp2_ != Long.MIN_VALUE;
	}
	
	public final boolean hasSecondBestPair() {
		return best2Pair_ != Long.MIN_VALUE;
	}
	
	public final long bestUnp1() {
		return bestUnp1_;
	}
	
	public final long secondBestUnp1() {
		return best2Unp1;
	}
	
	public final long bestUnp2() {
		return bestUnp2_;
	}
	
	public final long secondBestUnp2() {
		return best2Unp2;
	}
	
	public final long bestPair() {
		return bestPair_;
	}
	
	public final long secondBestPar() {
		return best2Pair;
	}
	
	public boolean sameRead(
			Read rd1,      // new mate #1
			Read rd2,      // new mate #2
			boolean qualitiesMatter) {// aln policy distinguishes b/t quals?{
	boolean same = false;
	if(rd1_ != null || rd2_ != null) {
		// This is not the first time the sink was initialized with
		// a read.  Check if new read/pair is identical to previous
		// read/pair
		if((rd1_ == null) == (rd1 == null) &&
		   (rd2_ == null) == (rd2 == null))
		{
			boolean m1same = (rd1 == null && rd1_ == null);
			if(!m1same) {
				assert(rd1 != null);
				assert(rd1_ != null);
				m1same = Read.same(
					rd1.patFw,  // new seq
					rd1.qual,   // new quals
					rd1_.patFw, // old seq
					rd1_.qual,  // old quals
					qualitiesMatter);
			}
			if(m1same) {
				boolean m2same = (rd2 == null && rd2_ == null);
				if(!m2same) {
					m2same = Read.same(
						rd2.patFw,  // new seq
						rd2.qual,   // new quals
						rd2_.patFw, // old seq
						rd2_.qual,  // old quals
						qualitiesMatter);
				}
				same = m2same;
			}
		}
	}
	return same;
	}
	
	public boolean prepareDiscordants() {
		if(rs1u_.size() == 1 && rs2u_.size() == 1) {
			rs1_.push_back(rs1u_[0]);
			rs2_.push_back(rs2u_[0]);
			return true;
		}
		return false;
	}
	
	public double selectAlnsToReport(EList<AlignmentResult> rs, long num, EList<Double> select, RandomSource rnd) {
		size_t sz = rs.size();
		if(sz < num) {
			num = sz;
		}
		if(sz < 1) {
			return 0;
		}
		select.resize((size_t)num);
		if(sz == 1) {
			assert_eq(1, num);
			select[0] = 0;
			return 0;
		}
		// Select a random offset into the list of alignments
		uint32_t off = rnd.nextU32() % (uint32_t)sz;
		uint32_t offOrig = off;
		// Now take elements starting at that offset, wrapping around to 0 if
		// necessary.  Leave the rest.
		for(size_t i = 0; i < num; i++) {
			select[i] = off;
			off++;
			if(off == sz) {
				off = 0;
			}
		}
		return offOrig;
	}
	
	public final double selectByScore(
			EList<AlignmentResult> rs1,    // alignments to select from (mate 1)
			const EList<AlignmentResult> rs2,    // alignments to select from (mate 2, or NULL)
			long             num,    // number of alignments to select
			EList<Integer>       select, // prioritized list to put results in
			const EList<AlignmentResult> rs1u,   // alignments to select from (mate 1)
			const EList<AlignmentResult> rs2u,   // alignments to select from (mate 2, or NULL)
			AlignmentScore            bestUScore,
			AlignmentScore            bestUDist,
			AlignmentScore            bestP1Score,
			AlignmentScore            bestP1Dist,
			AlignmentScore            bestP2Score,
			AlignmentScore            bestP2Dist,
			AlignmentScore            bestCScore,
			AlignmentScore            bestCDist,
			AlignmentScore            bestUnchosenUScore,
			AlignmentScore            bestUnchosenUDist,
			AlignmentScore            bestUnchosenP1Score,
			AlignmentScore            bestUnchosenP1Dist,
			AlignmentScore            bestUnchosenP2Score,
			AlignmentScore            bestUnchosenP2Dist,
			AlignmentScore            bestUnchosenCScore,
			AlignmentScore            bestUnchosenCDist,
			RandomSource        rnd){
		bestUScore.invalidate();
		bestUDist.invalidate();
		bestUnchosenUScore.invalidate();
		bestUnchosenUDist.invalidate();

		bestCScore.invalidate();
		bestP1Score.invalidate();
		bestP2Score.invalidate();
		bestCDist.invalidate();
		bestP1Dist.invalidate();
		bestP2Dist.invalidate();
		bestUnchosenCScore.invalidate();
		bestUnchosenP1Score.invalidate();
		bestUnchosenP2Score.invalidate();
		bestUnchosenCDist.invalidate();
		bestUnchosenP1Dist.invalidate();
		bestUnchosenP2Dist.invalidate();
		
		size_t sz = rs1->size(); // sz = # alignments found
		assert_leq(num, sz);
		if(sz < num) {
			num = sz;
		}
		// num = # to select
		if(sz == 0) {
			return 0;
		}
		select.resize((size_t)num);
		// Use 'selectBuf_' as a temporary list for sorting purposes
		EList<std::pair<AlnScore, size_t> >& buf =
			const_cast<EList<std::pair<AlnScore, size_t> >& >(selectBuf_);
		buf.resize(sz);
		// Sort by score.  If reads are pairs, sort by sum of mate scores.
		for(size_t i = 0; i < sz; i++) {
			buf[i].first = (*rs1)[i].score();
			if(rs2 != NULL) {
				buf[i].first += (*rs2)[i].score();
			}
			buf[i].second = i; // original offset
		}
		buf.sort(); buf.reverse(); // sort in descending order by score
		
		// Randomize streaks of alignments that are equal by score
		size_t streak = 0;
		for(size_t i = 1; i < buf.size(); i++) {
			if(buf[i].first == buf[i-1].first) {
				if(streak == 0) { streak = 1; }
				streak++;
			} else {
				if(streak > 1) {
					assert_geq(i, streak);
					buf.shufflePortion(i-streak, streak, rnd);
				}
				streak = 0;
			}
		}
		if(streak > 1) {
			buf.shufflePortion(buf.size() - streak, streak, rnd);
		}
		
		// Copy the permutation into the 'select' list
		for(size_t i = 0; i < num; i++) { select[i] = buf[i].second; }
		
		if(rs2 == NULL) {
			bestUScore = bestUDist = (*rs1)[select[0]].score();
		}
		
		// For paired-end read, find best alignment score among end
		// alignments not chosen, for both ends
		if(rs2 != NULL) {
			bestCScore = bestCDist = (*rs1)[select[0]].score() + (*rs2)[select[0]].score();
			bestP1Score = bestP1Dist = (*rs1)[select[0]].score();
			bestP2Score = bestP2Dist = (*rs2)[select[0]].score();
			for(size_t i = 0; i < rs1u->size(); i++) {
				if((*rs1u)[i].refcoord() == (*rs1)[select[0]].refcoord()) {
					continue;
				}
				if((*rs1u)[i].score() > bestUnchosenP1Score) {
					bestUnchosenP1Score = (*rs1u)[i].score();
				}
				if((*rs1u)[i].score().basesAligned() > bestUnchosenP1Dist.basesAligned()) {
					bestUnchosenP1Dist = (*rs1u)[i].score();
				}
			}
			for(size_t i = 0; i < rs2u->size(); i++) {
				if((*rs2u)[i].refcoord() == (*rs2)[select[0]].refcoord()) {
					continue;
				}
				if((*rs2u)[i].score() > bestUnchosenP2Score) {
					bestUnchosenP2Score = (*rs2u)[i].score();
				}
				if((*rs2u)[i].score().basesAligned() > bestUnchosenP2Dist.basesAligned()) {
					bestUnchosenP2Dist = (*rs2u)[i].score();
				}
			}
			if(buf.size() > 1) {
				bestUnchosenCScore = buf[1].first;
				for(size_t i = 1; i < buf.size(); i++) {
					AlnScore dist = (*rs1)[buf[i].second].score() +
					                (*rs2)[buf[i].second].score();
					if(dist.basesAligned() > bestUnchosenCDist.basesAligned()) {
						bestUnchosenCDist = dist;
					}
				}
			}
		} else if(buf.size() > 1) {
			bestUnchosenUScore = (*rs1)[buf[1].second].score();
			for(size_t i = 1; i < buf.size(); i++) {
				if((*rs1)[buf[1].second].score().basesAligned() > bestUnchosenUDist.basesAligned()) {
					bestUnchosenUDist = (*rs1)[buf[1].second].score();
				}
			}
		}
		
		// Returns index of the representative alignment, but in 'select' also
		// returns the indexes of the next best selected alignments in order by
		// score.
		return selectBuf_[0].second;
	}
	
	
}
