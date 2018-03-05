package com.uwb.bt2j.indexer.types;

import java.io.File;
import java.io.FileInputStream;
import java.io.OutputStream;

public class RefRecord {
	public long off;
	public long len;
	public boolean first;
	
	public enum ReadDir {
		REF_READ_FORWARD, // don't reverse reference sequence
				REF_READ_REVERSE,     // reverse entire reference sequence
				REF_READ_REVERSE_EACH // reverse each unambiguous stretch of reference
	}
	
	public RefRecord() {}
	public RefRecord(long _off, long _len, boolean _first) {
		off = _off;
		len = _len;
		first = _first;
	}
	
	public RefRecord(File in, boolean swap) {
		if(!in.canRead()) {
			System.err.println("Error reading RefRecord offset from FILE");
		}
		if(!in.canRead()) {
			System.err.println("Error reading RefRecord offset from FILE");
		}
		first = in.canRead() ? true : false;
	}
	
	public void write(OutputStream out, boolean be) {
		WordIO.writeU<Long>(out, off, be);
		WordIO.writeU<Long>(out, len, be);
		out.write(first ? 1 : 0);
	}
}
