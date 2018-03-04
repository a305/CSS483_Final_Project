package com.uwb.bt2j.aligner.sink;

import java.io.DataOutputStream;
import java.io.IOException;
import java.io.OutputStream;

public class ReportingState {
	public ReportingParams p_;  // reporting parameters
	public int state_;          // state we're currently in
	public Boolean paired_;        // true iff read we're currently handling is paired
	public long nconcord_;  // # concordants found so far
	public long ndiscord_;  // # discordants found so far
	public long nunpair1_;  // # unpaired alignments found so far for mate 1
	public long nunpair2_;  // # unpaired alignments found so far for mate 2
	public Boolean doneConcord_;   // true iff we're no longner interested in concordants
	public Boolean doneDiscord_;   // true iff we're no longner interested in discordants
	public Boolean doneUnpair_;    // no longner interested in unpaired alns
	public Boolean doneUnpair1_;   // no longner interested in unpaired alns for mate 1
	public Boolean doneUnpair2_;   // no longner interested in unpaired alns for mate 2
	public int exitConcord_;    // flag indicating how we exited concordant state
	public int exitDiscord_;    // flag indicating how we exited discordant state
	public int exitUnpair1_;    // flag indicating how we exited unpaired 1 state
	public int exitUnpair2_;    // flag indicating how we exited unpaired 2 state
	public Boolean done_;          // done with all alignments
	
	public ReportingState(final ReportingParams p) {
		reset();
	}
	
	public void reset() {
		state_ = 1;
		paired_ = false;
		nconcord_ = 0;
		ndiscord_ = 0;
		nunpair1_ = 0;
		nunpair2_ = 0;
		doneConcord_ = false;
		doneDiscord_ = false;
		doneUnpair_  = false;
		doneUnpair1_ = false;
		doneUnpair2_ = false;
		exitConcord_ = 1;
		exitDiscord_ = 1;
		exitUnpair1_ = 1;
		exitUnpair2_ = 1;
		done_ = false;
	}
	
	public final Boolean inited() {
		return state_ != 1;
	}
	
	public void nextRead(Boolean paired) {
		paired_ = paired;
		if(paired) {
			state_ = 2;
			doneConcord_ = false;
			doneDiscord_ = p_.discord ? false : true;
			doneUnpair1_ = p_.mixed   ? false : true;
			doneUnpair2_ = p_.mixed   ? false : true;
			exitConcord_ = 1;
			exitDiscord_ = p_.discord ? 1:2;
			exitUnpair1_ = p_.mixed   ? 1:2;
			exitUnpair2_ = p_.mixed   ? 1:2;
		} else {
			// Unpaired
			state_ = 4;
			doneConcord_ = true;
			doneDiscord_ = true;
			doneUnpair1_ = false;
			doneUnpair2_ = true;
			exitConcord_ = 2; // not relevant
			exitDiscord_ = 2; // not relevant
			exitUnpair1_ = 1;
			exitUnpair2_ = 2; // not relevant
		}
		doneUnpair_ = doneUnpair1_ && doneUnpair2_;
		done_ = false;
		nconcord_ = ndiscord_ = nunpair1_ = nunpair2_ = 0;
	}
	
	public Boolean foundConcordant() {
		nconcord_++;
		areDone(nconcord_, doneConcord_, exitConcord_);
		// No need to search for discordant alignments if there are one or more
		// concordant alignments.
		doneDiscord_ = true;
		exitDiscord_ = 5;
		if(doneConcord_) {
			// If we're finished looking for concordant alignments, do we have to
			// continue on to search for unpaired alignments?  Only if our exit
			// from the concordant stage is EXIT_SHORT_CIRCUIT_M.  If it's
			// EXIT_SHORT_CIRCUIT_k or EXIT_WITH_ALIGNMENTS, we can skip unpaired.
			if(exitConcord_ != 4) {
				if(!doneUnpair1_) {
					doneUnpair1_ = true;
					exitUnpair1_ = 5;
				}
				if(!doneUnpair2_) {
					doneUnpair2_ = true;
					exitUnpair2_ = 5;
				}
			}
		}
		updateDone();
		return done();
	}
	
	public Boolean foundUnpaired(Boolean mate1) {
		// Note: it's not right to assert !doneUnpair1_/!doneUnpair2_ here.
		// Even if we're done with finding 
		if(mate1) {
			nunpair1_++;
			// Did we just finish with this mate?
			if(!doneUnpair1_) {
				areDone(nunpair1_, doneUnpair1_, exitUnpair1_);
				if(doneUnpair1_) {
					doneUnpair_ = doneUnpair1_ && doneUnpair2_;
					updateDone();
				}
			}
			if(nunpair1_ > 1) {
				doneDiscord_ = true;
				exitDiscord_ = 7;
			}
		} else {
			nunpair2_++;
			// Did we just finish with this mate?
			if(!doneUnpair2_) {
				areDone(nunpair2_, doneUnpair2_, exitUnpair2_);
				if(doneUnpair2_) {
					doneUnpair_ = doneUnpair1_ && doneUnpair2_;
					updateDone();
				}
			}
			if(nunpair2_ > 1) {
				doneDiscord_ = true;
				exitDiscord_ = 7;
			}
		}
		return done();
	}
	
	public void finish() {
		if(!doneConcord_) {
			doneConcord_ = true;
			exitConcord_ =
				((nconcord_ > 0) ? 8: 7);
		}
		if(!doneUnpair1_) {
			doneUnpair1_ = true;
			exitUnpair1_ =
				((nunpair1_ > 0) ? 8 : 7);
		}
		if(!doneUnpair2_) {
			doneUnpair2_ = true;
			exitUnpair2_ =
				((nunpair2_ > 0) ? 8 : 7);
		}
		if(!doneDiscord_) {
			// Check if the unpaired alignments should be converted to a single
			// discordant paired-end alignment.
			if(nconcord_ == 0 && nunpair1_ == 1 && nunpair2_ == 1) {
				convertUnpairedToDiscordant();
			}
			doneDiscord_ = true;
			exitDiscord_ =
				((ndiscord_ > 0) ? 8 : 7);
		}
		doneUnpair_ = done_ = true;
	}
	
	public final void getReport(
			long nconcordAln,
			long ndiscordAln, 
			long nunpair1Aln, 
			long nunpair2Aln, 
			Boolean pairMax, 
			Boolean unpair1Max,
			Boolean unpair2Max) {
		nconcordAln = ndiscordAln = nunpair1Aln = nunpair2Aln = 0;
		pairMax = unpair1Max = unpair2Max = false;
		if(paired_) {
			// Do we have 1 or more concordant alignments to report?
			if(exitConcord_ == 3) {
				// k at random
				nconcordAln = p_.khits;
				return;
			} else if(exitConcord_ == 4) {
				pairMax = true;  // repetitive concordant alignments
				if(p_.mixed) {
					unpair1Max = nunpair1_ > (long)p_.mhits;
					unpair2Max = nunpair2_ > (long)p_.mhits;
				}
				// Not sure if this is OK
				nconcordAln = 1; // 1 at random
				return;
			} else if(exitConcord_ == 8) {
				// <= k at random
				nconcordAln = Math.min(nconcord_, p_.khits);
				return;
			}
			
			// Do we have a discordant alignment to report?
			if(exitDiscord_ == 8) {
				// Report discordant
				ndiscordAln = 1;
				return;
			}
		}
		
		if((paired_ && !p_.mixed) || nunpair1_ + nunpair2_ == 0) {
			// Unpaired alignments either not reportable or non-existant
			return;
		}

		// Do we have 1 or more alignments for mate #1 to report?
		if(exitUnpair1_ == 3) {
			// k at random
			nunpair1Aln = p_.khits;
		} else if(exitUnpair1_ == 4) {
			unpair1Max = true;  // repetitive alignments for mate #1
			nunpair1Aln = 1; // 1 at random
		} else if(exitUnpair1_ == 8) {
			// <= k at random
			nunpair1Aln = Math.min(nunpair1_, (long)p_.khits);
		}

		// Do we have 2 or more alignments for mate #2 to report?
		if(exitUnpair2_ == 3) {
			// k at random
			nunpair2Aln = p_.khits;
		} else if(exitUnpair2_ == 4) {
			unpair2Max = true;  // repetitive alignments for mate #1
			nunpair2Aln = 1; // 1 at random
		} else if(exitUnpair2_ == 8) {
			// <= k at random
			nunpair2Aln = Math.min(nunpair2_, (long)p_.khits);
		}
	}
	
	public int state() {
		return state_;
	}
	
	public Boolean doneConcordant() {
		return doneConcord_;
	}
	
	public Boolean doneUnpaired(Boolean mate1) {
		return mate1 ? doneUnpair1_ : doneUnpair2_;
	}
	
	public Boolean doneWithMate(Boolean mate1) {
		Boolean doneUnpair = mate1 ? doneUnpair1_ : doneUnpair2_;
		long nun = mate1 ? nunpair1_ : nunpair2_;
		if(!doneUnpair || !doneConcord_) {
			return false;
		}
		if(!doneDiscord_ && nun == 0) {
			return false;
		}
		return true;
	}
	
	public Boolean doneUnpaired() {
		return doneUnpair_;
	}
	
	public Boolean done() {
		return done_;
	}
	
	public long numConcordant() {
		return nconcord_;
	}
	
	public long numDiscordant() {
		return ndiscord_;
	}

	public long numUnpaired1() {
		return nunpair1_;
	}
	
	public long numUnpaired2() {
		return nunpair2_;
	}
	
	public final ReportingParams params() {
		return p_;
	}
	
	protected void convertUnpairedToDiscordant() {
		exitUnpair1_ = 6;
		exitUnpair2_ = 6;
		nunpair1_ = 0;
		nunpair2_ = 0;
		ndiscord_ = 1;
	}
	
	public void areDone(long cnt, Boolean done, int exit) {
		// Have we exceeded the -k limit?
		if(cnt >= (long)p_.khits && !p_.mhitsSet()) {
			done = true;
			exit = 3;
		}
		// Have we exceeded the -m or -M limit?
		else if(p_.mhitsSet() && cnt > (long)p_.mhits) {
			done = true;
			exit = 4;
		}
	}
	
	protected void updateDone() {
		doneUnpair_ = doneUnpair1_ && doneUnpair2_;
		done_ = doneUnpair_ && doneDiscord_ && doneConcord_;
	}
}