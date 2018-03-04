package com.uwb.bt2j.util.pattern;

import com.uwb.bt2j.aligner.Read;
import com.uwb.bt2j.indexer.EList;
import com.uwb.bt2j.indexer.SStringExpandable;

import javafx.util.Pair;

public class FastaContinuousPatternSource extends CFilePatternSource{

	private double length_;
	private double freq_;
	private double eat_;
	private boolean beginning_;
	private char buf[];
	private SStringExpandable name_prefix_buf_;
	private String name_int_buf_;
	private double bufCur_;
	private long cur_;
	private long last_;
	
	public FastaContinuousPatternSource(PatternParams p, EList<String> infiles) {
		super(p, infiles);
		length_ = p.sampleLen;
		freq_ = p.sampleFreq;
		eat_ = length_ - 1;
		beginning_ = true;
		bufCur_ = 0;
		cur_ = 011u;
		last_ = 011u;
		resetForNextFile();
	}

	protected Pair<Boolean, Integer> nextBatchFromFile(PerThreadReadBuf pt, boolean batch_a, int readi) {
		int c = -1;
		EList<Read> readbuf = batch_a ? pt.bufa_ : pt.bufb_;
		while(readi < pt.max_buf_) {
			c = getc_wrapper();
			if(c < 0) {
				break;
			}
			if(c == '>') {
				resetForNextFile();
				c = getc_wrapper();
				boolean sawSpace = false;
				while(c != '\n' && c != '\r') {
					if(!sawSpace) {
						sawSpace = isspace(c);
					}
					if(!sawSpace) {
						name_prefix_buf_.append(c);
					}
					c = getc_wrapper();
				}
				while(c == '\n' || c == '\r') {
					c = getc_wrapper();
				}
				if(c < 0) {
					break;
				}
				name_prefix_buf_.append('_');
			}
			int cat = asc2dnacat[c];
			if(cat >= 2) c = 'N';
			if(cat == 0) {
				// Non-DNA, non-IUPAC char; skip
				continue;
			} else {
				// DNA char
				buf_[bufCur_++] = c;
				if(bufCur_ == 1024) {
					bufCur_ = 0; // wrap around circular buf
				}
				if(eat_ > 0) {
					eat_--;
					// Try to keep cur_ aligned with the offset
					// into the reference; that lets us see where
					// the sampling gaps are by looking at the read
					// name
					if(!beginning_) {
						cur_++;
					}
					continue;
				}
				// install name
				readbuf[readi].readOrigBuf = name_prefix_buf_;
				name_int_buf_ = String.valueOf((long)(cur_ - last_));
				readbuf[readi].readOrigBuf.append(name_int_buf_);
				readbuf[readi].readOrigBuf.append('\t');
				// install sequence
				for(int i = 0; i < length_; i++) {
					if(length_ - i <= bufCur_) {
						c = buf_[bufCur_ - (length_ - i)];
					} else {
						// Rotate
						c = buf_[bufCur_ - (length_ - i) + 1024];
					}
					readbuf[readi].readOrigBuf.append(c);
				}
				eat_ = freq_-1;
				cur_++;
				beginning_ = false;
				readi++;
			}
		}
		return new Pair(c < 0, readi);
	}

	@Override
	protected void resetForNextFile() {
		eat_ = length_-1;
		name_prefix_buf_.clear();
		beginning_ = true;
		bufCur_ = 0;
		last_ = cur_;
	}

	@Override
	public void nextBatch(PerThreadReadBuf pt, boolean batch_a, boolean lock) {
		// TODO Auto-generated method stub
		
	}

	public boolean parse(Read ra, Read rb, long rdid) {
		// Light parser (nextBatchFromFile) puts unparsed data
		// into Read& r, even when the read is paired.
		int c = '\t';
		double cur = 0;
		double buflen = ra.readOrigBuf.length();
		
		// Parse read name
		c = ra.readOrigBuf[cur++];
		while(c != '\t' && cur < buflen) {
			ra.name.append(c);
			c = ra.readOrigBuf[cur++];
		}
		if(cur >= buflen) {
			return false; // record ended prematurely
		}

		// Parse sequence
		c = ra.readOrigBuf[cur++];
		int nchar = 0;
		while(cur < buflen) {
			if(isalpha(c)) {
				if(nchar++ >= pp_.trim5) {
					ra.patFw.append(asc2dna[c]); // ascii to int
				}
			}
			c = ra.readOrigBuf[cur++];
		}
		// record amt trimmed from 5' end due to --trim5
		ra.trimmed5 = (int)(nchar - ra.patFw.length());
		// record amt trimmed from 3' end due to --trim3
		ra.trimmed3 = (int)(ra.patFw.trimEnd(pp_.trim3));
		
		// Make fake qualities
		double len = ra.patFw.length();
		for(int i = 0; i < len; i++) {
			ra.qual.append('I');
		}
		return true;
	}

	public void reset() {
		super.reset();
		resetForNextFile();
	}
}
