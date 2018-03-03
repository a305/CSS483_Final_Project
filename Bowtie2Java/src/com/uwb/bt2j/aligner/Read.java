package com.uwb.bt2j.aligner;

import java.io.OutputStream;

import com.uwb.bt2j.util.strings.BTDnaString;
import com.uwb.bt2j.util.strings.BTString;
import com.uwb.bt2j.util.strings.SStringExpandable;

import javafx.util.Pair;

public class Read {
	
	public BTDnaString patFw;            // forward-strand sequence
	public BTDnaString patRc;            // reverse-complement sequence
	public BTString    qual;             // quality values

	public BTDnaString patFwRev;
	public BTDnaString patRcRev;
	public BTString    qualRev;

	// For remembering the exact input text used to define a read
	public SStringExpandable<String> readOrigBuf;

	public BTString name;      // read name
	public long  rdid;      // 0-based id based on pair's offset in read file(s)
	public int      mate;      // 0 = single-end, 1 = mate1, 2 = mate2
	public double seed;      // random seed
	public boolean     parsed;    // true iff read has been fully parsed
	int   ns_;       // # Ns
	char     filter;    // if read format permits filter char, set it here
	public int      trimmed5;  // amount actually trimmed off 5' end
	public int      trimmed3;  // amount actually trimmed off 3' end
	HitSet  hitset;    // holds previously-found hits; for chaining
	public Read() {
		reset();
	}
	
	public Read(String nm, String seq, String ql) {
		init(nm, seq, ql);
	}
	
	public void reset() {
		rdid = 0;
		trimmed5 = trimmed3 = 0;
		readOrigBuf.clear();
		patFw.clear();
		patRc.clear();
		qual.clear();
		patFwRev.clear();
		patRcRev.clear();
		qualRev.clear();
		name.clear();
		filter = '?';
		seed = 0;
		parsed = false;
		ns_ = 0;
	}
	
	public void finalize() {
		for(int i = 0; i < patFw.length(); i++) {
			if((int)patFw[i] > 3) {
				ns_++;
			}
		}
		constructRevComps();
		constructReverses();
	}
	
	public void init(String nm, String seq, String ql) {
		reset();
		patFw.installChars(seq);
		qual.install(ql);
		for(int i = 0; i < patFw.length(); i++) {
			if((int)patFw[i] > 3) {
				ns_++;
			}
		}
		constructRevComps();
		constructReverses();
		if(nm != null) name.install(nm);
	}
	
	public boolean empty() {
		return patFw.empty();
	}
	
	public int length() {
		return patFw.length();
	}
	
	public int ns() {
		return ns_;
	}
	
	public void constructRevComps() {
		patRc.installReverseComp(patFw);
	}
	
	public void constructReverses() {
		patFwRev.installReverse(patFw);
		patRcRev.installReverse(patRc);
		qualRev.installReverse(qual);
	}
	
	public void fixMateName(int i) {
		int namelen = name.length();
		boolean append = false;
		if(namelen < 2) {
			// Name is too short to possibly have /1 or /2 on the end
			append = true;
		} else {
			if(i == 1) {
				// append = true iff mate name does not already end in /1
				append =
					name[namelen-2] != '/' ||
					name[namelen-1] != '1';
			} else {
				// append = true iff mate name does not already end in /2
				append =
					name[namelen-2] != '/' ||
					name[namelen-1] != '2';
			}
		}
		if(append) {
			name.append('/');
			name.append("012".charAt(i));
		}
	}
	
	public void dump(OutputStream os) {
		os.write(name + ' ' +
		patFw
		+ ' '
		+ qual.toZBuf() +" ");
	}
	
	public static boolean same(BTDnaString seq1,
			BTString    qual1,
			BTDnaString seq2,
			BTString    qual2,
			boolean qualitiesMatter) {
		if(seq1.length() != seq2.length()) {
			return false;
		}
		for(int i = 0; i < seq1.length(); i++) {
			if(seq1[i] != seq2[i]) return false;
		}
		if(qualitiesMatter) {
			if(qual1.length() != qual2.length()) {
				return false;
			}
			for(int i = 0; i < qual1.length(); i++) {
				if(qual1[i] != qual2[i]) return false;
			}
		}
		return true;
	}
	
	public Pair<Integer, Integer> get(int off5p, boolean fw) {
		int c = (int)patFw[off5p];
        int q = qual[off5p];
		return new Pair((!fw && c < 4) ? (c ^ 3) : c, q - 33);
	}
	
	public int getc(int off5p, boolean fw) {
		int c = (int)patFw[off5p];
		return (!fw && c < 4) ? (c ^ 3) : c;
	}
	
	public int getq(int off5p) {
		int q = qual[off5p];
		return q-33;
	}
}