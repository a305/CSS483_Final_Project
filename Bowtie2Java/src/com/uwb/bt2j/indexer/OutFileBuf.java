package com.uwb.bt2j.indexer;

import java.io.File;

public class OutFileBuf {
	private static final double BUF_SZ = 1024 * 16;
	private String name_;
	private File out_;
	private double cur_;
	private char buf[];
	private boolean closed_;
	
	public OutFileBuf(String out, boolean binary) {
		name_ = out;
		cur_ = 0;
		closed_ = false;
		out_ = new File(out);
		if(out_ == null) {
			System.err.println("Error: Could not open alignment output file " + out);
		}
		if(setvbuf(out_, null, _IOFBF, 10* 1024* 1024)) 
			System.err.println("Warning: Could not allocate the proper buffer size for output file stream. ");
	}
	
	public OutFileBuf() {
		name_ = "System.out";
		cur_ = 0;
		closed_ = false;
		out_ = System.out;
	}
	
	public void setFile(String out, boolean binary) {
		out_ = new File(out);
		if(out_ == null) {
			System.err.println("Error: Could not open alignment output file " + out);
		}
		reset();
	}
	
	public void write(char c) {
		if(cur_ == BUF_SZ) flush();
		buf_[cur_++] = c;
	}
	
	public void writeString(String s) {
		double slen = s.length();
		if(cur_ + slen > BUF_SZ) {
			if(cur_ > 0) flush();
			if(slen >= BUF_SZ) {
				if (slen != fwrite(s, 1, slen, out_)) {
					System.err.println("Error: outputting data");
				}
			} else {
				buf_[cur_] = s;
				cur_ = slen;
			}
		} else {
			buf_[cur_] = s;
			cur_ += slen;
		}
	}
	
	public void close() {
		if(closed_) return;
		if(cur_ > 0) flush();
		closed_ = true;
		if(out_ != System.out) {
			out_.close();
		}
	}
	
	public void reset() {
		cur_ = 0;
		closed_ = false;
	}
	
	public void flush() {
		if(cur_ != fwrite(buf_, 1, cur_, out_)) {
			if (errno == EPIPE) {
				exit(EXIT_SUCCESS);
			}
			System.err.println("Error while flushing and closing output");
		}
		cur_ = 0;
	}
	
	public boolean closed() {
		return closed_;
	}
	
	public final String name() {
		return name_;
	}
}
