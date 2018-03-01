package com.uwb.bt2j.aligner;

public class AlignerMetrics {
	public int curBackTrracks_;
	public int curBwtOps_;
	
	protected Boolean first_;
	protected Boolean curIsLowEntropy_;
	protected Boolean curIsHomoPoly_;
	protected Boolean curHadRanges_;
	protected int curNumNs_;
	
	protected int reads_;
	protected int homoReads_;
	protected int lowEntReads_;
	protected int hiEntReads_;
	protected int alignedReads_;
	protected int unalignedReads_;
	protected int threeOrMoreNReads_;
	protected int lessThanThreeNRreads_;
	
	protected RunningStat bwtOpsPerRead_;
	protected RunningStat backtracksPerRead_;
	protected RunningStat bwtOpsPerHomoRead_;
	protected RunningStat backtracksPerHomoRead_;
	protected RunningStat bwtOpsPerLoEntRead_;
	protected RunningStat backtracksPerLoEntRead_;
	protected RunningStat bwtOpsPerHiEntRead_;
	protected RunningStat backtracksPerHiEntRead_;
	protected RunningStat bwtOpsPerAlignedRead_;
	protected RunningStat backtracksPerAlignedRead_;
	protected RunningStat bwtOpsPerUnalignedRead_;
	protected RunningStat backtracksPerUnalignedRead_;
	protected RunningStat bwtOpsPer0nRead_;
	protected RunningStat backtracksPer0nRead_;
	protected RunningStat bwtOpsPer1nRead_;
	protected RunningStat backtracksPer1nRead_;
	protected RunningStat bwtOpsPer2nRead_;
	protected RunningStat backtracksPer2nRead_;
	protected RunningStat bwtOpsPer3orMoreNRead_;
	protected RunningStat backtracksPer3orMoreNRead_;
	protected Timer timer_;
	
	public AlignerMetrics() {
		curBacktracks_ = 0;
		curBwtOps_ = 0;
		first_ = true;
		curIsLowEntropy_ =false;
		curIsHomoPoly_=false;
		curHadRanges_=false;
		curNumNs_=0;
		reads_=0;
		homoReads_=0;
		lowEntReads_=0;
		hiEntReads_=0;
		alignedReads_=0;
		unalignedReads_=0;
		threeOrMoreNReads_=0;
		lessThanThreeNRreads_=0;
		timer_(cout, "", false)
	}
	
	public void nextRead(BTDnaString read) {
		if(!first_) {
			finishRead();
		}
		first_ = false;
		//float ent = entropyDna5(read);
		float ent = 0.0f;
		curIsLowEntropy_ = (ent < 0.75f);
		curIsHomoPoly_ = (ent < 0.001f);
		curHadRanges_ = false;
		curBwtOps_ = 0;
		curBacktracks_ = 0;
		// Count Ns
		curNumNs_ = 0;
		final double len = read.length();
		for(double i = 0; i < len; i++) {
			if((int)read[i] == 4) curNumNs_++;
		}
	}
	
	public void setReadHasRange() {
		curHadRanges_ = true;
	}
	
	public void finishRead() {
		reads_++;
		if(curIsHomoPoly_) homoReads_++;
		else if(curIsLowEntropy_) lowEntReads_++;
		else hiEntReads_++;
		if(curHadRanges_) alignedReads_++;
		else unalignedReads_++;
		bwtOpsPerRead_.push((float)curBwtOps_);
		backtracksPerRead_.push((float)curBacktracks_);
		// Drill down by entropy
		if(curIsHomoPoly_) {
			bwtOpsPerHomoRead_.push((float)curBwtOps_);
			backtracksPerHomoRead_.push((float)curBacktracks_);
		} else if(curIsLowEntropy_) {
			bwtOpsPerLoEntRead_.push((float)curBwtOps_);
			backtracksPerLoEntRead_.push((float)curBacktracks_);
		} else {
			bwtOpsPerHiEntRead_.push((float)curBwtOps_);
			backtracksPerHiEntRead_.push((float)curBacktracks_);
		}
		// Drill down by whether it aligned
		if(curHadRanges_) {
			bwtOpsPerAlignedRead_.push((float)curBwtOps_);
			backtracksPerAlignedRead_.push((float)curBacktracks_);
		} else {
			bwtOpsPerUnalignedRead_.push((float)curBwtOps_);
			backtracksPerUnalignedRead_.push((float)curBacktracks_);
		}
		if(curNumNs_ == 0) {
			lessThanThreeNRreads_++;
			bwtOpsPer0nRead_.push((float)curBwtOps_);
			backtracksPer0nRead_.push((float)curBacktracks_);
		} else if(curNumNs_ == 1) {
			lessThanThreeNRreads_++;
			bwtOpsPer1nRead_.push((float)curBwtOps_);
			backtracksPer1nRead_.push((float)curBacktracks_);
		} else if(curNumNs_ == 2) {
			lessThanThreeNRreads_++;
			bwtOpsPer2nRead_.push((float)curBwtOps_);
			backtracksPer2nRead_.push((float)curBacktracks_);
		} else {
			threeOrMoreNReads_++;
			bwtOpsPer3orMoreNRead_.push((float)curBwtOps_);
			backtracksPer3orMoreNRead_.push((float)curBacktracks_);
		}
	}
	
	public void printSummary() {
		if(!first_) {
			finishRead();
		}
		System.out.println("AlignerMetrics:");
		System.out.println("  # Reads:             " + reads_);
		float hopct = (reads_ > 0) ? (((float)homoReads_)/((float)reads_)) : (0.0f);
		hopct *= 100.0f;
		System.out.println("  % homo-polymeric:    " + (hopct));
		float lopct = (reads_ > 0) ? ((float)lowEntReads_/(float)(reads_)) : (0.0f);
		lopct *= 100.0f;
		System.out.println("  % low-entropy:       " + (lopct));
		float unpct = (reads_ > 0) ? ((float)unalignedReads_/(float)(reads_)) : (0.0f);
		unpct *= 100.0f;
		System.out.println("  % unaligned:         " + (unpct));
		float npct = (reads_ > 0) ? ((float)threeOrMoreNReads_/(float)(reads_)) : (0.0f);
		npct *= 100.0f;
		System.out.println( "  % with 3 or more Ns: " + (npct));
		System.out.println();
		System.out.println( "  Total BWT ops:    avg: " + bwtOpsPerRead_.mean() + ", stddev: " + bwtOpsPerRead_.stddev());
		System.out.println("  Total Backtracks: avg: " + backtracksPerRead_.mean() + ", stddev: " + backtracksPerRead_.stddev());
		time_t elapsed = timer_.elapsed();
		System.out.println("  BWT ops per second:    " + (bwtOpsPerRead_.tot()/elapsed));
		System.out.println("  Backtracks per second: " + (backtracksPerRead_.tot()/elapsed));
		System.out.println();
		System.out.println("  Homo-poly:" );
		System.out.println("    BWT ops:    avg: " + bwtOpsPerHomoRead_.mean() + ", stddev: " + bwtOpsPerHomoRead_.stddev() );
		System.out.println("    Backtracks: avg: " + backtracksPerHomoRead_.mean() + ", stddev: " + backtracksPerHomoRead_.stddev() );
		System.out.println("  Low-entropy:" );;
		System.out.println("    BWT ops:    avg: " + bwtOpsPerLoEntRead_.mean() + ", stddev: " + bwtOpsPerLoEntRead_.stddev() );
		System.out.println("    Backtracks: avg: " + backtracksPerLoEntRead_.mean() + ", stddev: " + backtracksPerLoEntRead_.stddev() );
		System.out.println("  High-entropy:" );
		System.out.println("    BWT ops:    avg: " + bwtOpsPerHiEntRead_.mean() + ", stddev: " + bwtOpsPerHiEntRead_.stddev() );
		System.out.println("    Backtracks: avg: " + backtracksPerHiEntRead_.mean() + ", stddev: " + backtracksPerHiEntRead_.stddev() );
		System.out.println();
		System.out.println("  Unaligned:" );
		System.out.println("    BWT ops:    avg: " + bwtOpsPerUnalignedRead_.mean() + ", stddev: " + bwtOpsPerUnalignedRead_.stddev() );
		System.out.println("    Backtracks: avg: " + backtracksPerUnalignedRead_.mean() + ", stddev: " + backtracksPerUnalignedRead_.stddev() );;
		System.out.println("  Aligned:" );
		System.out.println("    BWT ops:    avg: " + bwtOpsPerAlignedRead_.mean() + ", stddev: " + bwtOpsPerAlignedRead_.stddev() );
		System.out.println("    Backtracks: avg: " + backtracksPerAlignedRead_.mean() + ", stddev: " + backtracksPerAlignedRead_.stddev() );
		System.out.println();
		System.out.println("  0 Ns:" );
		System.out.println("    BWT ops:    avg: " + bwtOpsPer0nRead_.mean() + ", stddev: " + bwtOpsPer0nRead_.stddev() );
		System.out.println("    Backtracks: avg: " + backtracksPer0nRead_.mean() + ", stddev: " + backtracksPer0nRead_.stddev() );
		System.out.println("  1 N:" );
		System.out.println("    BWT ops:    avg: " + bwtOpsPer1nRead_.mean() + ", stddev: " + bwtOpsPer1nRead_.stddev() );
		System.out.println("    Backtracks: avg: " + backtracksPer1nRead_.mean() + ", stddev: " + backtracksPer1nRead_.stddev() );
		System.out.println("  2 Ns:" );
		System.out.println("    BWT ops:    avg: " + bwtOpsPer2nRead_.mean() + ", stddev: " + bwtOpsPer2nRead_.stddev() );
		System.out.println("    Backtracks: avg: " + backtracksPer2nRead_.mean() + ", stddev: " + backtracksPer2nRead_.stddev() );
		System.out.println("  >2 Ns:" );
		System.out.println("    BWT ops:    avg: " + bwtOpsPer3orMoreNRead_.mean() + ", stddev: " + bwtOpsPer3orMoreNRead_.stddev() );
		System.out.println("    Backtracks: avg: " + backtracksPer3orMoreNRead_.mean() + ", stddev: " + backtracksPer3orMoreNRead_.stddev() );
		System.out.println();
	}
}
