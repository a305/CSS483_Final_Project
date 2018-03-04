package com.uwb.bt2j.aligner.sink;

public class ReportingMetrics {
	long nread;         // # reads
	long npaired;       // # pairs
	long nunpaired;     // # unpaired reads
	
	// Paired
	
	//  Concordant
	long nconcord_uni;  // # pairs with unique concordant alns
	long nconcord_uni1; // # pairs with exactly 1 concordant alns
	long nconcord_uni2; // # pairs with >1 concordant aln, still unique
	long nconcord_rep;  // # pairs with repetitive concordant alns
	long nconcord_0;    // # pairs with 0 concordant alns
	//  Discordant
	long ndiscord;      // # pairs with 1 discordant aln
	
	//  Unpaired from failed pairs
	long nunp_0_uni;    // # unique from nconcord_0_ - ndiscord_
	long nunp_0_uni1;   // # pairs with exactly 1 concordant alns
	long nunp_0_rep;    // # repetitive from 
	long nunp_0_0;      // # with 0 alignments

	//  Unpaired from repetitive pairs
	private long nunp_rep_uni;  // # pairs with unique concordant alns
	private long nunp_rep_uni1; // # pairs with exactly 1 concordant alns
	private long nunp_rep_uni2; // # pairs with >1 concordant aln, still unique
	private long nunp_rep_rep;  // # pairs with repetitive concordant alns
	private long nunp_rep_0;    // # pairs with 0 concordant alns
	
	// Unpaired
	
	long nunp_uni;      // # unique from nconcord_0_ - ndiscord_
	long nunp_uni1;     // # pairs with exactly 1 concordant alns
	long nunp_uni2;     // # pairs with >1 concordant aln, still unique
	long nunp_rep;      // # repetitive from 
	long nunp_0;        // # with 0 alignments

	
	private long sum_best1;     // Sum of all the best alignment scores
	private long sum_best2;     // Sum of all the second-best alignment scores
	private long sum_best;      // Sum of all the best and second-best
	
	public void reset() {
		init(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
	}
	
	public void init(
			long nread_,
			long npaired_,
			long nunpaired_,
			long nconcord_uni_,
			long nconcord_uni1_,
			long nconcord_uni2_,
			long nconcord_rep_,
			long nconcord_0_,
			long ndiscord_,
			long nunp_0_uni_,
			long nunp_0_uni1_,
			long nunp_0_uni2_,
			long nunp_0_rep_,
			long nunp_0_0_,
			long nunp_rep_uni_,
			long nunp_rep_uni1_,
			long nunp_rep_uni2_,
			long nunp_rep_rep_,
			long nunp_rep_0_,
			long nunp_uni_,
			long nunp_uni1_,
			long nunp_uni2_,
			long nunp_rep_,
			long nunp_0_,
			long sum_best1_,
			long sum_best2_,
			long sum_best_) {
		nread         = nread_;
		
		npaired       = npaired_;
		nunpaired     = nunpaired_;
		
		nconcord_uni  = nconcord_uni_;
		nconcord_uni1 = nconcord_uni1_;
		nconcord_uni2 = nconcord_uni2_;
		nconcord_rep  = nconcord_rep_;
		nconcord_0    = nconcord_0_;
		
		ndiscord      = ndiscord_;
		
		nunp_0_uni    = nunp_0_uni_;
		nunp_0_uni1   = nunp_0_uni1_;
		nunp_0_uni2   = nunp_0_uni2_;
		nunp_0_rep    = nunp_0_rep_;
		nunp_0_0      = nunp_0_0_;

		nunp_rep_uni  = nunp_rep_uni_;
		nunp_rep_uni1 = nunp_rep_uni1_;
		nunp_rep_uni2 = nunp_rep_uni2_;
		nunp_rep_rep  = nunp_rep_rep_;
		nunp_rep_0    = nunp_rep_0_;

		nunp_uni      = nunp_uni_;
		nunp_uni1     = nunp_uni1_;
		nunp_uni2     = nunp_uni2_;
		nunp_rep      = nunp_rep_;
		nunp_0        = nunp_0_;

		sum_best1     = sum_best1_;
		sum_best2     = sum_best2_;
		sum_best      = sum_best_;
	}
	
	public void merge(final ReportingMetrics met) {
		nread         += met.nread;

		npaired       += met.npaired;
		nunpaired     += met.nunpaired;

		nconcord_uni  += met.nconcord_uni;
		nconcord_uni1 += met.nconcord_uni1;
		nconcord_uni2 += met.nconcord_uni2;
		nconcord_rep  += met.nconcord_rep;
		nconcord_0    += met.nconcord_0;

		ndiscord      += met.ndiscord;

		nunp_0_uni    += met.nunp_0_uni;
		nunp_0_uni1   += met.nunp_0_uni1;
		nunp_0_uni2   += met.nunp_0_uni2;
		nunp_0_rep    += met.nunp_0_rep;
		nunp_0_0      += met.nunp_0_0;

		nunp_rep_uni  += met.nunp_rep_uni;
		nunp_rep_uni1 += met.nunp_rep_uni1;
		nunp_rep_uni2 += met.nunp_rep_uni2;
		nunp_rep_rep  += met.nunp_rep_rep;
		nunp_rep_0    += met.nunp_rep_0;

		nunp_uni      += met.nunp_uni;
		nunp_uni1     += met.nunp_uni1;
		nunp_uni2     += met.nunp_uni2;
		nunp_rep      += met.nunp_rep;
		nunp_0        += met.nunp_0;

		sum_best1     += met.sum_best1;
		sum_best2     += met.sum_best2;
		sum_best      += met.sum_best;
	}
}
