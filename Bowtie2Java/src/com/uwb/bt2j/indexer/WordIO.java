package com.uwb.bt2j.indexer;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

public class WordIO <T> {
	public void writeU(OutputStream out, T x) throws IOException {
		out.write((int) x);
	}
	
	public void writeI(OutputStream out, T x) throws IOException {
		out.write((int) x);
	}
	
	public T readU(InputStream in) {
		T x;
		x = in.read();
			return x;
	}
}
