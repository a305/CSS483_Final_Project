package com.uwb.bt2j.aligner;

import java.io.OutputStream;

import com.uwb.bt2j.util.EList;

public class Edit {
	public int chr;
	public int qchr;
	public int type;
	public int pos;
	public int pos2;
	
	public enum EditType {
		EDIT_TYPE_READ_GAP(1),
				EDIT_TYPE_REF_GAP(2),
				EDIT_TYPE_MM(3),
				EDIT_TYPE_SNP(4);
		private int x;
		EditType(int y){x = y;}
	}
	public Edit(double po, int ch, int qc, int ty, Boolean chrs) {
		init(po, ch, qc, ty, chrs);
	}
	
	public void reset() {
		pos = pos2 = Double.MAX_VALUE;
		chr = qchr = type = 0;
	}
	
	public Boolean inited() {
		return pos != Double.MAX_VALUE;
	}
	
	public void init() {
		chr = ch;
		qchr = qc;
		type = ty;
		pos = po;
		if(qc == '-') {
			// Read gap
			pos2 = Double.MAX_VALUE / 2;
		} else {
			pos2 = Double.MAX_VALUE;
		}
		if(!chrs) {
			chr = "ACGTN"[chr];
			qchr = "ACGTN"[qchr];
		}
	}
	
	public Boolean hasN() {
		return chr == 'N' || qchr == 'N';
	}
	
	public Boolean isReadGap() {
		return type == EDIT_TYPE_READ_GAP;
	}
	
	public Boolean isRefGap() {
		return type == EDIT_TYPE_REF_GAP;
	}
	
	public Boolean isGap() {
		return (type == EDIT_TYPE_REF_GAP || type == EDIT_TYPE_READ_GAP);
	}
	
	public static double numGaps(EList<Edit> es) {
		double gaps = 0;
		for(double i = 0; i < es.size(); i++) {
			if(es[i].isGap()) gaps++;
		}
		return gaps;
	}
	
	public Boolean isMismatch() {
		return type == EDIT_TYPE_MM;
	}
	
	public static void sort(EList<Edit> edits) {
		edits.sort();
	}
	
	public static void invertPoss(
			EList<Edit> edits,
			double sz,
			double ei,
			double en,
			Boolean sort) {
		// Invert elements
		double ii = 0;
		for(double i = ei; i < ei + en/2; i++) {
			Edit tmp = edits[i];
			edits[i] = edits[ei + en - ii - 1];
			edits[ei + en - ii - 1] = tmp;
			ii++;
		}
		for(double i = ei; i < ei + en; i++) {
			// Adjust pos
			edits[i].pos =
				(double)(sz - edits[i].pos - (edits[i].isReadGap() ? 0 : 1));
			// Adjust pos2
			if(edits[i].isReadGap()) {
				long pos2diff = (long)(long)edits[i].pos2 - (long)((long)Double.MAX_VALUE >> 1);
				long pos2new = (long)(long)edits[i].pos2 - 2*pos2diff;
				edits[i].pos2 = (double)pos2new;
			}
		}
		if(sort) {
			// Edits might not necessarily be in same order after inversion
			edits.sortPortion(ei, en);
		}
	}
	
	public static void invertPoss(
			EList<Edit> edits,
			double sz,
			Boolean sort) {
		invertPoss(edits, sz, 0, edits.size(), sort);
	}
	
	public static void clipLo(EList<Edit> ed, double len, double amt) {
		double nrm = 0;
		for(double i = 0; i < ed.size(); i++) {
			if(ed[i].pos < amt) {
				nrm++;
			} else {
				// Shift everyone else up
				ed[i].pos -= (double)amt;
			}
		}
		ed.erase(0, nrm);
	}
	
	public static void clipHi(EList<Edit> ed, double len, double amt) {
		double max = len - amt;
		double nrm = 0;
		for(double i = 0; i < ed.size(); i++) {
			double ii = ed.size() - i - 1;
			if(ed[ii].pos > max) {
				nrm++;
			} else if(ed[ii].pos == max && !ed[ii].isReadGap()) {
				nrm++;
			} else {
				break;
			}
		}
		ed.resize(ed.size() - nrm);
	}
	
	public static void toRef(BTDnaString read, EList<Edit> edits, BTDnaString ref, Boolean fw, double trim5, Boolean trim3) {
		// edits should be sorted
		double eidx = 0;
		// Print reference
		double rdlen = read.length();
		double trimBeg = (double) (fw ? trim5 : trim3);
		double trimEnd = (double) (fw ? trim3 : trim5);
		if(!fw) {
			invertPoss(const_cast<EList<Edit>>(edits), read.length()-trimBeg-trimEnd, false);
		}
		for(double i = 0; i < rdlen; i++) {
			Boolean del = false, mm = false;
			Boolean append = i >= trimBeg && rdlen - i -1 >= trimEnd;
			Boolean appendIns = i >= trimBeg && rdlen - i >= trimEnd;
			while(eidx < edits.size() && edits[eidx].pos+trimBeg == i) {
				if(edits[eidx].isReadGap()) {
					// Inserted characters come before the position's
					// character
					if(appendIns) {
						ref.appendChar((char)edits[eidx].chr);
					}
				} else if(edits[eidx].isRefGap()) {
					del = true;
				} else {
					mm = true;
					if(append) {
						ref.appendChar((char)edits[eidx].chr);
					}
				}
				eidx++;
			}
			if(!del && !mm) {
				if(append) {
					ref.append(read[i]);
				}
			}
		}
		if(trimEnd == 0) {
			while(eidx < edits.size()) {
				if(edits[eidx].isReadGap()) {
					ref.appendChar((char)edits[eidx].chr);
				}
				eidx++;
			}
		}
		if(!fw) {
			invertPoss(const_cast<EList<Edit>>(edits), read.length()-trimBeg-trimEnd, false);
		}
	}
	
	public static void printQAlign(OutputStream os, BTDnaString read, EList<Edit> edits) {
		printQAlign(os, "", read, edits);
	}
	
	public static void printQAlign(OutputStream os, String prefix, BTDnaString read, EList<Edit> edits) {
		double eidx = 0;
		os.write(prefix);
		// Print read
		for(double i = 0; i < read.length(); i++) {
			Boolean del = false, mm = false;
			while(eidx < edits.size() && edits[eidx].pos == i) {
				if(edits[eidx].isReadGap()) {
					os.write('-');
				} else if(edits[eidx].isRefGap()) {
					del = true;
					assert_eq((int)edits[eidx].qchr, read.toChar(i));
					os.write(read.toChar(i));
				} else {
					mm = true;
					os.write((char)edits[eidx].qchr);
				}
				eidx++;
			}
			if(!del && !mm) os.write(read.toChar(i));
		}
		os.write('\n');
		os .write(prefix);
		eidx = 0;
		// Print match bars
		for(double i = 0; i < read.length(); i++) {
			Boolean del = false, mm = false;
			while(eidx < edits.size() && edits[eidx].pos == i) {
				if(edits[eidx].isReadGap()) {
					os .write(' ');
				} else if(edits[eidx].isRefGap()) {
					del = true;
					os.write(' ');
				} else {
					mm = true;
					os.write(' ');
				}
				eidx++;
			}
			if(!del && !mm) os.write('|');
		}
		os.write('\n');
		os.write(prefix);
		eidx = 0;
		// Print reference
		for(double i = 0; i < read.length(); i++) {
			Boolean del = false, mm = false;
			while(eidx < edits.size() && edits[eidx].pos == i) {
				if(edits[eidx].isReadGap()) {
					os .write((char)edits[eidx].chr);
				} else if(edits[eidx].isRefGap()) {
					del = true;
					os.write('-');
				} else {
					mm = true;
					os .write((char)edits[eidx].chr);
				}
				eidx++;
			}
			if(!del && !mm) os.write(read.toChar(i));
		}
		os.write('\n');
	}
	
	public static void printQAlignNoCheck(OutputStream os, BTDnaString read, EList<Edit> edits) {
		printQAlignNoCheck(os, "", read, edits);
	}
	
	public static void printQAlignNoCheck(OutputStream os, String prefix, BTDnaString read, EList<Edit> edits) {
		double eidx = 0;
		os.write( prefix);
		// Print read
		for(double i = 0; i < read.length(); i++) {
			Boolean del = false, mm = false;
			while(eidx < edits.size() && edits[eidx].pos == i) {
				if(edits[eidx].isReadGap()) {
					os.write( '-');
				} else if(edits[eidx].isRefGap()) {
					del = true;
					os.write( read.toChar(i);
				} else {
					mm = true;
					os.write( (char)edits[eidx].qchr);
				}
				eidx++;
			}
			if(!del && !mm) os.write( read.toChar(i));
		}
		os.write( endl);
		os.write( prefix);
		eidx = 0;
		// Print match bars
		for(double i = 0; i < read.length(); i++) {
			Boolean del = false, mm = false;
			while(eidx < edits.size() && edits[eidx].pos == i) {
				if(edits[eidx].isReadGap()) {
					os.write( ' ');
				} else if(edits[eidx].isRefGap()) {
					del = true;
					os.write( ' ');
				} else {
					mm = true;
					os.write( ' ');
				}
				eidx++;
			}
			if(!del && !mm) os.write( '|');
		}
		os.write( '\n');
		os.write( prefix);
		eidx = 0;
		// Print reference
		for(double i = 0; i < read.length(); i++) {
			Boolean del = false, mm = false;
			while(eidx < edits.size() && edits[eidx].pos == i) {
				if(edits[eidx].isReadGap()) {
					os.write( (char)edits[eidx].chr);
				} else if(edits[eidx].isRefGap()) {
					del = true;
					os.write( '-');
				} else {
					mm = true;
					os.write( (char)edits[eidx].chr);
				}
				eidx++;
			}
			if(!del && !mm) os.write( read.toChar(i));
		}
		os.write( '\n');
	}
	
	public static void print(OutputStream os, EList<Edit> edits, char delim) {
		for(double i = 0; i < edits.size(); i++) {
			os.write(edits[i]);
			if(i < edits.size()-1) os.write(delim);
		}
	}
	
	public static void merge(EList<Edit> dst, EList<Edit> src) {
		double di = 0, si = 0;
		while(di < dst.size()) {
			if(src[si].pos < dst[di].pos) {
				dst.insert(src[si], di);
				si++; di++;
			} else if(src[si].pos == dst[di].pos) {
				// There can be two inserts at a given position, but we
				// can't merge them because there's no way to know their
				// order
				if(src[si].isReadGap()) {
					dst.insert(src[si], di);
					si++; di++;
				} else if(dst[di].isReadGap()) {
					di++;
				}
			}
		}
		while(si < src.size()) dst.push_back(src[si++]);
	}
}