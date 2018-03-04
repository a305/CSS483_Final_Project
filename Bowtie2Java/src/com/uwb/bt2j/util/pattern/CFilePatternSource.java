package com.uwb.bt2j.util.pattern;

import java.io.File;

import com.uwb.bt2j.aligner.Read;
import com.uwb.bt2j.indexer.EList;

import javafx.util.Pair;

public abstract class CFilePatternSource extends PatternSource{
	public EList<String> infiles_;
	public ELIst<boolean> errs_;
	public double filecur_;
	public File fp_;
	public gzFile zfp_;
	public boolean is_open_;
	public double skip_;
	public boolean first_;
	public char buf_[];
	public boolean compressed_;
	
	public CFilePatternSource(PatternParams p,EList<String> infiles) {
		super(p);
		infiles_ = infiles;
		filecur_ = 0;
		fp_ = null;
		zfp_ = null;
		is_open_ = false;
		skip_ = p.skip;
		first_ = true;
		compressed_ = false;
		errs_.resize(infiles_.size());
		errs_.fill(0, infiles_.size(), false);
		open(); // open first file in the list
		filecur_++;
	}

	protected abstract Pair<Boolean, Integer> nextBatchFromFile(PerThreadReadBuf pt, boolean batch_a, int read_idx);
	
	private Pair<Boolean, Integer> nextBatchImpl(PerThreadReadBuf pt, boolean batch_a) {
		boolean done = false;
		int nread = 0;
		pt.setReadId(readCnt_);
		while(true) { // loop that moves on to next file when needed
			do {
				Pair<Boolean, Integer> ret = nextBatchFromFile(pt, batch_a, nread);
				done = ret.first;
				nread = ret.second;
			} while(!done && nread == 0); // not sure why this would happen
			if(done && filecur_ < infiles_.size()) { // finished with this file
				open();
				resetForNextFile(); // reset state to handle a fresh file
				filecur_++;
				if(nread == 0 || (nread < pt.max_buf_)) {
					continue;
				}
				done = false;
			}
			break;
		}
		readCnt_ += nread;
		return new Pair(done, nread);
	}
	
	protected abstract void resetForNextFile();
	
	public void open() {
		if(is_open_) {
			is_open_ = false;
			if (compressed_) {
				gzclose(zfp_);
				zfp_ = null;
			}
			else {
	      			fp_ = null;
	      		}
		}
		while(filecur_ < infiles_.size()) {
			if(infiles_[filecur_] == "-") {
				// always assume that data from stdin is compressed
				compressed_ = true;
				int fn = dup(fileno(stdin));
				zfp_ = gzdopen(fn, "rb");
			}
			else {
				compressed_ = false;
				if (is_gzipped_file(infiles_[filecur_])) {
					compressed_ = true;
					zfp_ = gzopen(infiles_[filecur_], "rb");
				}
				else {
					fp_ = fopen(infiles_[filecur_], "rb");
				}
				if((compressed_ && zfp_ == null) || (!compressed_ && fp_ == null)) {
					if(!errs_[filecur_]) {
						System.err.println("Warning: Could not open read file \""
						     + infiles_[filecur_]
						     + "\" for reading; skipping...");
						errs_[filecur_] = true;
	      				}
	      				filecur_++;
	      				continue;
				}
			}
			is_open_ = true;
			return;
		}
		System.err.println("Error: No input read files were valid");
		return;
	}
	
	public int getc_wrapper() {
		return gzgetc(zfp_);
	}
	
	public int ungetc_wrapper(int c) {
		return gzungetc(c, zfp_);
	}
	
	public boolean is_gzipped_file(String filename) {
		int pos = filename.lastIndexOf(".");
		String ext = (pos == -1) ? "" : filename.substring(pos + 1);
		if (ext == "" || ext == "gz" || ext == "Z") {
			return true;
		}
		return false;
	}

	public void nextBatch(PerThreadReadBuf pt, boolean batch_a) {
			// synchronization at this level because both reading and manipulation of
			// current file pointer have to be protected
			ThreadSafe ts(mutex);
			return nextBatchImpl(pt, batch_a);
	}

	@Override
	public abstract boolean parse(Read ra, Read rb, long rdid);

	@Override
	public void reset() {
		readCnt_ = 0;
		filecur_ = 0;
		open();
		filecur_++;
	}
}
