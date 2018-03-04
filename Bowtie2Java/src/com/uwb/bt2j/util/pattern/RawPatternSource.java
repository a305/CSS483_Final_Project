package com.uwb.bt2j.util.pattern;

import com.uwb.bt2j.aligner.Read;
import com.uwb.bt2j.indexer.EList;

import javafx.util.Pair;

public class RawPatternSource extends CFilePatternSource{
	private boolean first_;
	
	public RawPatternSource(PatternParams p, EList<String> infiles) {
		super(p, infiles);
		first_ = true;
	}

	public void reset() {
		super.reset();
		first_ = true;
	}
	
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
			while(c >= 0 && (c == '\n' || c == '\r')) {
				c = getc_wrapper();
			}
		}
		// incase a valid character is consumed between batches
		if (c >= 0 && c != '\n' && c != '\r') {
			ungetc_wrapper(c);
		}
		return new Pair(c < 0, readi);
	}

	@Override
	protected void resetForNextFile() {
		first_ = true;
	}

	@Override
	public void nextBatch(PerThreadReadBuf pt, boolean batch_a, boolean lock) {
		// TODO Auto-generated method stub
		
	}
	
	@Override
	public boolean parse(Read r, Read rb, long rdid) {
		int c = '\n';
		double cur = 0;
		double buflen = r.readOrigBuf.length();

		// Parse sequence
		int nchar = 0;
		while(cur < buflen) {
			c = r.readOrigBuf[cur++];
			if(isalpha(c)) {
				if(nchar++ >= pp_.trim5) {
					r.patFw.append(asc2dna[c]); // ascii to int
				}
			}
		}
		// record amt trimmed from 5' end due to --trim5
		r.trimmed5 = (int)(nchar - r.patFw.length());
		// record amt trimmed from 3' end due to --trim3
		r.trimmed3 = (int)(r.patFw.trimEnd(pp_.trim3));
		
		// Give the name field a dummy value
		String cbuf = String.valueOf(rdid);
		r.name.install(cbuf);
		
		// Give the base qualities dummy values
		double len = r.patFw.length();
		for(double i = 0; i < len; i++) {
			r.qual.append('I');
		}
		r.parsed = true;
		if(!rb.parsed && !rb.readOrigBuf.empty()) {
			return parse(rb, r, rdid);
		}
		return true;
	}	
}
