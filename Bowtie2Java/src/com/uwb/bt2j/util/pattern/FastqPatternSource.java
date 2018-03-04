package com.uwb.bt2j.util.pattern;

import com.uwb.bt2j.aligner.Read;
import com.uwb.bt2j.indexer.EList;

import javafx.util.Pair;

public class FastqPatternSource extends CFilePatternSource{

	private boolean first_;
	private boolean interleaved_;
	
	public FastqPatternSource(PatternParams p, EList<String> infiles, boolean interleaved) {
		super(p, infiles);
		first_ = true;
		interleaved_ = interleaved;
	}
	
	public void reset() {
		super.reset();
		first_ = true;
	}

	protected Pair<Boolean, Integer> nextBatchFromFile(PerThreadReadBuf pt, boolean batch_a, int readi) {
		int c = -1;
		EList<Read> readbuf = batch_a ? pt.bufa_ : pt.bufb_;
		if(first_) {
			c = getc_wrapper();
			if (c == -1) {
				return new Pair(true, 0);
			}
			while(c == '\r' || c == '\n') {
				c = getc_wrapper();
			}
			if(c != '@') {
				System.err.println("Error: reads file does not look like a FASTQ file");
			}
			first_ = false;
			readbuf[readi].readOrigBuf.append('@');	
		}
		
		boolean done = false, aborted = false;
		// Read until we run out of input or until we've filled the buffer
		while (readi < pt.max_buf_ && !done) {
			long buf = readbuf[readi].readOrigBuf;
			int newlines = 4;
			while(newlines) {
				c = getc_wrapper();
				done = c < 0;
				if(c == '\n' || (done && newlines == 1)) {
					// Saw newline, or EOF that we're
					// interpreting as final newline
					newlines--;
					c = '\n';
				} else if(done) {
					// account for newline at the end of the file
					if (newlines == 4) {
						newlines = 0;
					}
					else {
						aborted = true; // Unexpected EOF
					}
					break;
				}
				buf.append(c);
			}
			if (c > 0) {
				if (interleaved_) {
					// alternate between read buffers
					batch_a = !batch_a;
					readbuf = batch_a ? pt.bufa_ : pt.bufb_;
					// increment read counter after each pair gets read
					readi = batch_a ? readi+1 : readi;
				}
				else {
					readi++;
				}
			}
		}
		if(aborted) {
			readi--;
		}
		return new Pair(done, readi);
	}

	@Override
	protected void resetForNextFile() {
		first_ = true;
	}

	@Override
	public void nextBatch(PerThreadReadBuf pt, boolean batch_a, boolean lock) {
		// TODO Auto-generated method stub
		
	}

	public boolean parse(Read r, Read rb, long rdid) {
		// We assume the light parser has put the raw data for the separate ends
		// into separate Read objects. That doesn't have to be the case, but
		// that's how we've chosen to do it for FastqPatternSource

		int c;
		double cur = 1;
		double buflen = r.readOrigBuf.length();

		// Parse read name
		while(true) {
			c = r.readOrigBuf[cur++];
			if(c == '\n' || c == '\r') {
				do {
					c = r.readOrigBuf[cur++];
				} while(c == '\n' || c == '\r');
				break;
			}
			r.name.append(c);
		}
		
		// Parse sequence
		int nchar = 0;
		while(c != '+') {
			if(c == '.') {
				c = 'N';
			}
			if(isalpha(c)) {
				// If it's past the 5'-end trim point
				if(nchar++ >= pp_.trim5) {
					r.patFw.append(asc2dna[c]);
				}
			}
			c = r.readOrigBuf[cur++];
		}
		r.trimmed5 = (int)(nchar - r.patFw.length());
		r.trimmed3 = (int)(r.patFw.trimEnd(pp_.trim3));
		
		do {
			c = r.readOrigBuf[cur++];
		} while(c != '\n' && c != '\r');
		while(cur < buflen && (c == '\n' || c == '\r')) {
			c = r.readOrigBuf[cur++];
		}
		
		if(nchar > 0) {
			int nqual = 0;
			if (pp_.intQuals) {
				int cur_int = 0;
				while(c != '\t' && c != '\n' && c != '\r') {
					cur_int *= 10;
					cur_int += (int)(c - '0');
					c = r.readOrigBuf[cur++];
					if(c == ' ' || c == '\t' || c == '\n' || c == '\r') {
						char cadd = intToPhred33(cur_int, pp_.solexa64);
						cur_int = 0;
						if(++nqual > pp_.trim5) {
							r.qual.append(cadd);
						}
					}
				}
			} else {
				c = charToPhred33(c, pp_.solexa64, pp_.phred64);
				if(nqual++ >= r.trimmed5) {
					r.qual.append(c);
				}
				while(cur < r.readOrigBuf.length()) {
					c = r.readOrigBuf[cur++];
					if (c == ' ') {
						wrongQualityFormat(r.name);
						return false;
					}
					if(c == '\r' || c == '\n') {
						break;
					}
					c = charToPhred33(c, pp_.solexa64, pp_.phred64);
					if(nqual++ >= r.trimmed5) {
						r.qual.append(c);
					}
				}
				r.qual.trimEnd(r.trimmed3);
				if(r.qual.length() < r.patFw.length()) {
					tooFewQualities(r.name);
					return false;
				} else if(r.qual.length() > r.patFw.length()) {
					tooManyQualities(r.name);
					return false;
				}
			}
		}
		// Set up a default name if one hasn't been set
		if(r.name.empty()) {
			String cbuf;
			cbuf = String.valueOf((long)rdid);
			r.name.install(cbuf);
		}
		r.parsed = true;
		if(!rb.parsed && !rb.readOrigBuf.empty()) {
			return parse(rb, r, rdid);
		}
		return true;
	}
	
}
