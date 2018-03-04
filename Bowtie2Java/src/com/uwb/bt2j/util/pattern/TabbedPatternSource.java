package com.uwb.bt2j.util.pattern;

import com.uwb.bt2j.aligner.Read;
import com.uwb.bt2j.indexer.EList;

import javafx.util.Pair;

public class TabbedPatternSource extends CFilePatternSource{
	protected boolean secondName_;

	public TabbedPatternSource(PatternParams p, EList<String> infiles, boolean secondName) {
		super(p, infiles);
		secondName_ = secondName;
	}

	
	@Override
	protected Pair<Boolean, Integer> nextBatchFromFile(PerThreadReadBuf pt, boolean batch_a, int readi) {
		int c = getc_wrapper();
		while(c >= 0 && (c == '\n' || c == '\r')) {
			c = getc_wrapper();
		}
		EList<Read> readbuf = batch_a ? pt.bufa_ : pt.bufb_;
		// Read until we run out of input or until we've filled the buffer
		for(; readi < pt.max_buf_ && c >= 0; readi++) {
			readbuf[readi].readOrigBuf.clear();
			while(c >= 0 && c != '\n' && c != '\r') {
				readbuf[readi].readOrigBuf.append(c);
				c = getc_wrapper();
			}
			while(c >= 0 && (c == '\n' || c == '\r') && readi < pt.max_buf_ - 1) {
				c = getc_wrapper();
			}
		}
		return new Pair(c < 0, readi);
	}

	@Override
	protected void resetForNextFile() {
		// TODO Auto-generated method stub
		
	}
	
	public boolean parse(Read ra, Read rb, long rdid) {
		// Light parser (nextBatchFromFile) puts unparsed data
		// into Read& r, even when the read is paired.

		int c = '\t';
		double cur = 0;
		double buflen = ra.readOrigBuf.length();
		
		// Loop over the two ends
		for(int endi = 0; endi < 2 && c == '\t'; endi++) {
			Read r = ((endi == 0) ? ra : rb);
			// Parse name if (a) this is the first end, or
			// (b) this is tab6
			if(endi < 1 || secondName_) {
				// Parse read name
				c = ra.readOrigBuf[cur++];
				while(c != '\t' && cur < buflen) {
					r.name.append(c);
					c = ra.readOrigBuf[cur++];
				}
				if(cur >= buflen) {
					return false; // record ended prematurely
				}
			} else if(endi > 0) {
				// if this is the second end and we're parsing
				// tab5, copy name from first end
				rb.name = ra.name;
			}

			// Parse sequence
			c = ra.readOrigBuf[cur++];
			int nchar = 0;
			while(c != '\t' && cur < buflen) {
				if(isalpha(c)) {
					if(nchar++ >= pp_.trim5) {
						r.patFw.append(asc2dna[c]); // ascii to int
					}
				}
				c = ra.readOrigBuf[cur++];
			}
			if(cur >= buflen) {
				return false; // record ended prematurely
			}
			// record amt trimmed from 5' end due to --trim5
			r.trimmed5 = (int)(nchar - r.patFw.length());
			// record amt trimmed from 3' end due to --trim3
			r.trimmed3 = (int)(r.patFw.trimEnd(pp_.trim3));
			
			// Parse qualities
			c = ra.readOrigBuf[cur++];
			int nqual = 0;
			if (pp_.intQuals) {
				int cur_int = 0;
				while(c != '\t' && c != '\n' && c != '\r' && cur < buflen) {
					cur_int *= 10;
					cur_int += (int)(c - '0');
					c = ra.readOrigBuf[cur++];
					if(c == ' ' || c == '\t' || c == '\n' || c == '\r') {
						char cadd = intToPhred33(cur_int, pp_.solexa64);
						cur_int = 0;
						if(++nqual > pp_.trim5) {
							r.qual.append(cadd);
						}
					}
				}
			} else {
				while(c != '\t' && c != '\n' && c != '\r') {
					if(c == ' ') {
						wrongQualityFormat(r.name);
						return false;
					}
					char cadd = charToPhred33(c, pp_.solexa64, pp_.phred64);
					if(++nqual > pp_.trim5) {
						r.qual.append(cadd);
					}
					if(cur >= buflen) break;
					c = ra.readOrigBuf[cur++];
				}
			}
			if(nchar > nqual) {
				tooFewQualities(r.name);
				return false;
			} else if(nqual > nchar) {
				tooManyQualities(r.name);
				return false;
			}
			r.qual.trimEnd(pp_.trim3);
		}
		return true;
	}	
	

	@Override
	public void nextBatch(PerThreadReadBuf pt, boolean batch_a, boolean lock) {
		// TODO Auto-generated method stub
		
	}

}
