package com.uwb.bt2j.aligner.seed;

import com.uwb.bt2j.aligner.Edit;
import com.uwb.bt2j.aligner.Read;
import com.uwb.bt2j.util.strings.BTDnaString;
import com.uwb.bt2j.util.strings.BTString;
import com.uwb.bt2j.util.types.EList;

public class SeedResults {
	protected double              nonzTot_;     // # offsets with non-zero size
	protected double              uniTot_;      // # offsets unique hit
	protected double              uniTotS_[];  // # offsets unique hit on each strand
	protected double              repTot_;      // # offsets repetitive hit
	protected double              repTotS_[];  // # offsets repetitive hit on each strand
	protected double              nonzFw_;      // # offsets into fw read with non-0 size
	protected double              nonzRc_;      // # offsets into rc read with non-0 size
	protected double              numRanges_;   // # ranges added
	protected double              numElts_;     // # elements added
	protected double              numRangesFw_; // # ranges added for fw seeds
	protected double              numEltsFw_;   // # elements added for fw seeds
	protected double              numRangesRc_; // # ranges added for rc seeds
	protected double              numEltsRc_;   // # elements added for rc seeds
	protected EList<BTDnaString>  seqFw_;       // seqs for seeds from forward read
	protected EList<BTDnaString>  seqRc_;       // seqs for seeds from revcomp read
	protected EList<BTString>     qualFw_;      // quals for seeds from forward read
	protected EList<BTString>     qualRc_;      // quals for seeds from revcomp read
	protected EList<QVal>         hitsFw_;      // hits for forward read
	protected EList<QVal>         hitsRc_;      // hits for revcomp read
	protected EList<EList<InstantiatedSeed> > isFw_; // hits for forward read
	protected EList<EList<InstantiatedSeed> > isRc_; // hits for revcomp read
	protected EList<boolean>         sortedFw_;    // true iff fw QVal was sorted/ranked
	protected EList<boolean>         sortedRc_;    // true iff rc QVal was sorted/ranked
	protected EList<double> offIdx2off_;
	protected EList<double> rankOffs_;
	protected EList<boolean> rankFws_;
	protected boolean sorted_;
	protected double numOffs_;
	protected Read read_;
	protected EEHit exactFwHit_;
	protected EEHit exactRcHit_;
	protected EList<EEHit> mm1Hit_;
	protected double mm1Elt_;
	protected boolean mm1Sorted_;
	protected EList<double> tmpMedian_;
	
	public SeedResults() {
		seqFw_ = seqRc_ = qualFw_ = qualRc_ = hitsFw_ = hitsRc_ = isFw_ = isRc_ = sortedFw_ = sortedRc_ = offIdx2off_ = rankOffs_ = rankFws_ = mm1Hit_ = 5;
		clear();
	}
	
	public void nextRead(Read read) {
		read_ = read;
	}
	
	public void add(
			QVal qv,           // range of ranges in cache
			AlignmentCache ac, // cache
			double seedIdx,         // seed index (from 5' end)
			boolean     seedFw) {
		if(qv.empty()) return;
		if(seedFw) {
			hitsFw_[seedIdx] = qv;
			numEltsFw_ += qv.numElts();
			numRangesFw_ += qv.numRanges();
			if(qv.numRanges() > 0) nonzFw_++;
		} else {
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
	}
	
	public void reset(Read read, EList<Double> offIdx2off, double numOffs) {
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
		for(double i = 0; i < numOffs_; i++) {
			sortedFw_[i] = sortedRc_[i] = false;
			hitsFw_[i].reset();
			hitsRc_[i].reset();
			isFw_[i].clear();
			isRc_[i].clear();
		}
		read_ = read;
		sorted_ = false;
	}
	
	public void reset(Read read, EList<Integer> offIdx2off, int numOffs) {
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
		read_ = read;
		sorted_ = false;
	}
	
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
	
	public void clear() {
		clearSeeds();
		read_ = null;
		exactFwHit_.reset();
		exactRcHit_.reset();
		mm1Hit_.clear();
		mm1Sorted_ = false;
		mm1Elt_ = 0;
	}
	
	public void toSeedAlSumm(SeedAlSumm ssum) {
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
				for(double i = 0; i < numOffs_; i++) {
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
	
	public float averageHitsPerSeed() {
		return nonzTot_ == 0 ? 0 : (float)numElts_ / (float)nonzTot_;
	}
	
	public double numUniqueSeeds() {
		return uniTot_;
	}
	
	public double numUniqueSeedsStrand(boolean fw) {
		return uniTotS_[fw ? 0 : 1];
	}
	
	public double numRepeatSeeds() {
		return repTot_;
	}
	
	public double numRepeatSeedsStrand(boolean fw) {
		return repTotS_[fw ? 0 : 1];
	}
	
	public float medianHitsPerSeed() {
		EList<double> median =_cast<EList<double>>(tmpMedian_);
		median.clear();
		for(double i = 0; i < numOffs_; i++) {
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
	
	public double uniquenessFactor() {
		double result = 0.0;
		for(double i = 0; i < numOffs_; i++) {
			if(hitsFw_[i].valid()) {
				double nelt = hitsFw_[i].numElts();
				result += (1.0 / (double)(nelt * nelt));
			}
			if(hitsRc_[i].valid()) {
				double nelt = hitsRc_[i].numElts();
				result += (1.0 / (double)(nelt * nelt));
			}
		}
		return result;
	}
	
	/**
	 * Return the number of ranges being held.
	 */
	public double numRanges() { return numRanges_; }

	/**
	 * Return the number of elements being held.
	 */
	public double numElts() { return numElts_; }

	/**
	 * Return the number of ranges being held for seeds on the forward
	 * read strand.
	 */
	public double numRangesFw() { return numRangesFw_; }

	/**
	 * Return the number of elements being held for seeds on the
	 * forward read strand.
	 */
	public double numEltsFw() { return numEltsFw_; }

	/**
	 * Return the number of ranges being held for seeds on the
	 * reverse-complement read strand.
	 */
	public double numRangesRc() { return numRangesRc_; }

	/**
	 * Return the number of elements being held for seeds on the
	 * reverse-complement read strand.
	 */
	public double numEltsRc() { return numEltsRc_; }
	
	/**
	 * Given an offset index, return the offset that has that index.
	 */
	public double idx2off(double off) {
		return offIdx2off_[off];
	}
	
	/**
	 * Return true iff there are 0 hits being held.
	 */
	public boolean empty() { return numRanges() == 0; }
	
	public QVal hitsAtOffIdx(boolean fw, double seedoffidx) {
		return fw ? hitsFw_[seedoffidx] : hitsRc_[seedoffidx];
	}
	
	public EList<InstantiatedSeed> instantiatedSeeds(boolean fw, double seedoffidx) {
		return fw ? isFw_[seedoffidx] : isRc_[seedoffidx];
	}
	
	public double numOffs() {
		return numOffs_;
	}
	
	public Read read() {
		return read_;
	}
	
	public void rankSeedHits(RandomSource rnd, boolean all) {
		if(all) {
			for(double i = 1; i < numOffs_; i++) {
				for(int fwi = 0; fwi <= 1; fwi++) {
					boolean fw = fwi == 0;
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
				long minsz = Long.MAX_VALUE;
				double minidx = 0;
				boolean minfw = true;
				// Rank seed-hit positions in ascending order by number of elements
				// in all BW ranges
				boolean rb = rnd.nextBool();
				for(int fwi = 0; fwi <= 1; fwi++) {
					boolean fw = (fwi == (rb ? 1 : 0));
					EList<QVal> rrs = (fw ? hitsFw_ : hitsRc_);
					EList<boolean> sorted = (fw ? sortedFw_ : sortedRc_);
					double i = (rnd.nextU32() % (double)numOffs_);
					for(double ii = 0; ii < numOffs_; ii++) {
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
				if(minfw) {
					sortedFw_[minidx] = true;
				} else {
					sortedRc_[minidx] = true;
				}
				rankOffs_.push_back(minidx);
				rankFws_.push_back(minfw);
			}
		}
		sorted_ = true;
	}
	
	public double nonzeroOffsets() {
		return nonzTot_;
	}
	
	public boolean allFwSeedsHit() {
		return nonzFw_ == numOffs();
	}
	
	public boolean allRcSeedsHit() {
		return nonzRc_ == numOffs();
	}
	
	public double fewestEditsEE(boolean fw, int seedlen, int per) {
		double nonz = fw ? nonzFw_ : nonzRc_;
		if(nonz < numOffs()) {
			int maxdepth = (seedlen + per - 1) / per;
			int missing = (int)(numOffs() - nonz);
			return (missing + maxdepth - 1) / maxdepth;
		} else {
			// Exact hit is possible (not guaranteed)
			return 0;
		}
	}
	
	public double nonzeroOffsetsFw() {
		return nonzFw_;
	}
	
	public double nonzeroOffsetsRc() {
		return nonzRc_;
	}
	
	public QVal hitsByRank(
			double    r,       // in
			double offidx,  // out
			double off,     // out
			boolean    fw,      // out
			double seedlen) {
		if(rankFws_[r]) {
			fw = true;
			offidx = rankOffs_[r];
			off = offIdx2off_[offidx];
			seedlen = (double)seqFw_[rankOffs_[r]].length();
			return hitsFw_[rankOffs_[r]];
		} else {
			fw = false;
			offidx = rankOffs_[r];
			off = offIdx2off_[offidx];
			seedlen = (double)seqRc_[rankOffs_[r]].length();
			return hitsRc_[rankOffs_[r]];
		}
	}
	
	public BTDnaString seqByRank(double r) {
		return rankFws_[r] ? seqFw_[rankOffs_[r]] : seqRc_[rankOffs_[r]];
	}
	
	public BTString qualByRank(double r) {
		return rankFws_[r] ? qualFw_[rankOffs_[r]] : qualRc_[rankOffs_[r]];
	}
	
	public EList<BTDnaString> seqs(boolean fw) {
		return fw ? seqFw_ : seqRc_; 
	}
	
	public EList<BTString> quals(boolean fw) {
		return fw ? qualFw_ : qualRc_;
	}
	
	public EEHit exactFwEEHit() {
		return exactFwHit_;
	}
	
	public EEHit exactRcEEHit() {
		return exactRcHit_;
	}
	
	public EList<EEHit> mm1EEHits() {
		return mm1Hit_;
	}
	
	public void sort1mmEe(RandomSource rnd) {
		mm1Hit_.sort();
		double streak = 0;
		for(double i = 1; i < mm1Hit_.size(); i++) {
			if(mm1Hit_[i].score == mm1Hit_[i-1].score) {
				if(streak == 0) { streak = 1; }
				streak++;
			} else {
				if(streak > 1) {
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
	
	public void add1mmEe(
			long top,
			long bot,
			Edit e1,
			Edit e2,
			boolean fw,
			long score) {
		mm1Hit_.expand();
		mm1Hit_.back().init(top, bot, e1, e2, fw, score);
		mm1Elt_ += (bot - top);
	}
	
	public void addExactEeRc(
			long top,
			long bot,
			Edit e1,
			Edit e2,
			boolean fw,
			long score) {
		exactRcHit_.init(top, bot, e1, e2, fw, score);
	}
	
	public void clear1mmE2eHits() {
		mm1Hit_.clear();     // 1-mismatch end-to-end hits
		mm1Elt_ = 0;         // number of 1-mismatch hit rows
		mm1Sorted_ = false;  // true iff we've sorted the mm1Hit_ list
	}
	
	public double numE2eHits() {
		return exactFwHit_.size() + exactRcHit_.size() + mm1Elt_;
	}
	
	public double numExactE2eHits() {
		return exactFwHit_.size() + exactRcHit_.size();
	}
	
	public double num1mmE2eHits() {
		return mm1Elt_;
	}
	
	public double readLength() {
		return read_.length();
	}

	public void addExactEeFw(long top, long bot, Edit e1, Edit e2, boolean fw, long score) {
		exactFwHit_.init(top, bot, e1, e2, fw, score);
	}
}
