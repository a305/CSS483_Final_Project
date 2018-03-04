package com.uwb.bt2j.util.pattern;

import com.uwb.bt2j.aligner.Read;
import com.uwb.bt2j.indexer.EList;
import com.uwb.bt2j.indexer.FileBuf;

import javafx.util.Pair;

public class FastaPatternSource extends CFilePatternSource {
	protected boolean first_;
	
	public FastaPatternSource(PatternParams p, EList<String> infiles) {
		super(p, infiles);
		first_ = true;
	}
	
	public boolean parse(Read r, Read rb, long rdid) {
		// We assume the light parser has put the raw data for the separate ends
		// into separate Read objects.	That doesn't have to be the case, but
		// that's how we've chosen to do it for FastqPatternSource

		int c = -1;
		double cur = 1;
		double buflen = r.readOrigBuf.length();
		
		// Parse read name
		while(cur < buflen) {
			c = r.readOrigBuf[cur++];
			if(c == '\n' || c == '\r') {
				do {
					c = r.readOrigBuf[cur++];
				} while((c == '\n' || c == '\r') && cur < buflen);
				break;
			}
			r.name.append(c);
		}
		if(cur >= buflen) {
			return false; // FASTA ended prematurely
		}
		
		// Parse sequence
		int nchar = 0;
		while(cur < buflen) {
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
			if ((c == '\n' || c == '\r')
					&& cur < buflen
					&& r.readOrigBuf[cur] != '>') {
				c = r.readOrigBuf[cur++];
			}
		}
		r.trimmed5 = (int)(nchar - r.patFw.length());
		r.trimmed3 = (int)(r.patFw.trimEnd(pp_.trim3));
		
		for(int i = 0; i < r.patFw.length(); i++) {
			r.qual.append('I');
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

	@Override
	protected Pair<Boolean, Integer> nextBatchFromFile(PerThreadReadBuf pt, boolean batch_a, int readi) {
		int c;
		EList<Read> readbuf = batch_a ? pt.bufa_ : pt.bufb_;
		if(first_) {
			c = getc_wrapper();
			if (c == -1) {
				return new Pair(true, 0);
			}
			while(c == '\r' || c == '\n') {
				c = getc_wrapper();
			}
			if(c != '>') {
				System.err.println("Error: reads file does not look like a FASTA file");
			}
			first_ = false;
		}
		boolean done = false;
		// Read until we run out of input or until we've filled the buffer
		for(; readi < pt.max_buf_ && !done; readi++) {
			long buf = readbuf[readi].readOrigBuf;
			buf.clear();
			buf.append('>');
			while(true) {
				c = getc_wrapper();
				if(c < 0 || c == '>') {
					done = c < 0;
					break;
				}
				buf.append(c);
			}
		}
		// Immediate EOF case
		if(done && readbuf[readi-1].readOrigBuf.length() == 1) {
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
	
	protected static int skipToNextFastaRecord(FileBuf in) {
		int c;
		while((c = in.get()) != '>') {
			if(in.eof()) return -1;
		}
		return c;
	}
	
	
	public void reset() {
		first_ = true;
		super.reset();
	}
}
