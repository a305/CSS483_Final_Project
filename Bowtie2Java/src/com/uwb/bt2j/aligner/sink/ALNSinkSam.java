package com.uwb.bt2j.aligner.sink;

public class ALNSinkSam {
	protected final SamConfig samc_;
	protected BTDnaString dseq_;
	protected BTString dqual_;
	
	public ALNSinkSam() {
		
	}
	
	public void append(
			BString o,
			StackedAln staln,
			double threadId,
			final Read rd1,
			final Read rd2,
			final long rdid,
			AlnRes rs1,
			AlnRes rs2,
			final AlnSetSumm summ,
			final SeedAlSumm ssm1,
			final SeedAlSumm ssm2,
			final AlnFlags flags1,
			final AlnFlags flags2,
			final PerReadMetrics prm,
			final Mapq mapq,
			final Scoring sc,
			Boolean report2) {
		if(rd1 != null) {
			appendMate(o, staln, rd1, rd2, rdid, rs1, rs2, summ, ssm1, ssm2,
			           flags1, prm, mapq, sc);
		}
		if(rd2 != null && report2) {
			appendMate(o, staln, rd2, rd1, rdid, rs2, rs1, summ, ssm2, ssm1,
			           flags2, prm, mapq, sc);
		}
	}
	
	protected void appendMate() {
		if(rs == NULL && samc_.omitUnalignedReads()) {
			return;
		}
		char buf[1024];
		char mapqInps[1024];
		if(rs != NULL) {
			staln.reset();
			rs->initStacked(rd, staln);
			staln.leftAlign(false /* not past MMs */);
		}
		int offAdj = 0;
		// QNAME
		samc_.printReadName(o, rd.name, flags.partOfPair());
		o.append('\t');
		// FLAG
		int fl = 0;
		if(flags.partOfPair()) {
			fl |= SAM_FLAG_PAIRED;
			if(flags.alignedConcordant()) {
				fl |= SAM_FLAG_MAPPED_PAIRED;
	 		}
			if(!flags.mateAligned()) {
				// Other fragment is unmapped
				fl |= SAM_FLAG_MATE_UNMAPPED;
			}
			fl |= (flags.readMate1() ?
				SAM_FLAG_FIRST_IN_PAIR : SAM_FLAG_SECOND_IN_PAIR);
			if(flags.mateAligned()) {
				if(!flags.isOppFw()) {
					fl |= SAM_FLAG_MATE_STRAND;
				}
			}
		}
		if(!flags.isPrimary()) {
			fl |= SAM_FLAG_NOT_PRIMARY;
		}
		if(rs != NULL && !rs->fw()) {
			fl |= SAM_FLAG_QUERY_STRAND;
		}
		if(rs == NULL) {
			// Failed to align
			fl |= SAM_FLAG_UNMAPPED;
		}
		itoa10<int>(fl, buf);
		o.append(buf);
		o.append('\t');
		// RNAME
		if(rs != NULL) {
			samc_.printRefNameFromIndex(o, (size_t)rs->refid());
			o.append('\t');
		} else {
			if(summ.orefid() != -1) {
				// Opposite mate aligned but this one didn't - print the opposite
				// mate's RNAME and POS as is customary
				assert(flags.partOfPair());
				samc_.printRefNameFromIndex(o, (size_t)summ.orefid());
			} else {		
				// No alignment
				o.append('*');
			}
			o.append('\t');
		}
		// POS
		// Note: POS is *after* soft clipping.  I.e. POS points to the
		// upstream-most character *involved in the clipped alignment*.
		if(rs != NULL) {
			itoa10<int64_t>(rs->refoff()+1+offAdj, buf);
			o.append(buf);
			o.append('\t');
		} else {
			if(summ.orefid() != -1) {
				// Opposite mate aligned but this one didn't - print the opposite
				// mate's RNAME and POS as is customary
				assert(flags.partOfPair());
				itoa10<int64_t>(summ.orefoff()+1+offAdj, buf);
				o.append(buf);
			} else {
				// No alignment
				o.append('0');
			}
			o.append('\t');
		}
		// MAPQ
		mapqInps[0] = '\0';
		if(rs != NULL) {
			itoa10<TMapq>(mapqCalc.mapq(
				summ, flags, rd.mate < 2, rd.length(),
				rdo == NULL ? 0 : rdo->length(), mapqInps), buf);
			o.append(buf);
			o.append('\t');
		} else {
			// No alignment
			o.append("0\t");
		}
		// CIGAR
		if(rs != NULL) {
			staln.buildCigar(flags.xeq());
			staln.writeCigar(&o, NULL);
			o.append('\t');
		} else {
			// No alignment
			o.append("*\t");
		}
		// RNEXT
		if(rs != NULL && flags.partOfPair()) {
			if(rso != NULL && rs->refid() != rso->refid()) {
				samc_.printRefNameFromIndex(o, (size_t)rso->refid());
				o.append('\t');
			} else {
				o.append("=\t");
			}
		} else if(summ.orefid() != -1) {
			// The convention if this mate fails to align but the other doesn't is
			// to copy the mate's details into here
			o.append("=\t");
		} else {
			o.append("*\t");
		}
		// PNEXT
		if(rs != NULL && flags.partOfPair()) {
			if(rso != NULL) {
				itoa10<int64_t>(rso->refoff()+1, buf);
				o.append(buf);
				o.append('\t');
			} else {
				// The convenstion is that if this mate aligns but the opposite
				// doesn't, we print this mate's offset here
				itoa10<int64_t>(rs->refoff()+1, buf);
				o.append(buf);
				o.append('\t');
			}
		} else if(summ.orefid() != -1) {
			// The convention if this mate fails to align but the other doesn't is
			// to copy the mate's details into here
			itoa10<int64_t>(summ.orefoff()+1, buf);
			o.append(buf);
			o.append('\t');
		} else {
			o.append("0\t");
		}
		// ISIZE
		if(rs != NULL && rs->isFraglenSet()) {
			itoa10<int64_t>(rs->fragmentLength(), buf);
			o.append(buf);
			o.append('\t');
		} else {
			// No fragment
			o.append("0\t");
		}
		// SEQ
		if(!flags.isPrimary() && samc_.omitSecondarySeqQual()) {
			o.append('*');
		} else {
			// Print the read
			if(rd.patFw.length() == 0) {
				o.append('*');
			} else {
				if(rs == NULL || rs->fw()) {
					o.append(rd.patFw.toZBuf());
				} else {
					o.append(rd.patRc.toZBuf());
				}
			}
		}
		o.append('\t');
		// QUAL
		if(!flags.isPrimary() && samc_.omitSecondarySeqQual()) {
			o.append('*');
		} else {
			// Print the quals
			if(rd.qual.length() == 0) {
				o.append('*');
			} else {
				if(rs == NULL || rs->fw()) {
					o.append(rd.qual.toZBuf());
				} else {
					o.append(rd.qualRev.toZBuf());
				}
			}
		}
		o.append('\t');
		//
		// Optional fields
		//
		if(rs != NULL) {
			samc_.printAlignedOptFlags(
				o,           // output buffer
				true,        // first opt flag printed is first overall?
				rd,          // read
				rdo,         // opposite read
				*rs,         // individual alignment result
				staln,       // stacked alignment
				flags,       // alignment flags
				summ,        // summary of alignments for this read
				ssm,         // seed alignment summary
				prm,         // per-read metrics
				sc,          // scoring scheme
				mapqInps);   // inputs to MAPQ calculation
		} else {
			samc_.printEmptyOptFlags(
				o,           // output buffer
				true,        // first opt flag printed is first overall?
				rd,          // read
				flags,       // alignment flags
				summ,        // summary of alignments for this read
				ssm,         // seed alignment summary
				prm,         // per-read metrics
				sc);         // scoring scheme
		}
		o.append('\n');
	}
}
