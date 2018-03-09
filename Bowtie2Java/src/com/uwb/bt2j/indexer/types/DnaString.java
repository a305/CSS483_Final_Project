package com.uwb.bt2j.indexer.types;
import com.uwb.bt2j.indexer.*;
import com.uwb.bt2j.indexer.util.Alphabet;

public class DnaString {
	protected char cs_[];
	
	public DnaString(String b, boolean chars, boolean colors) {
		if(chars) {
			if(colors) {
				installColors(b, b.length());
			} else {
				installChars(b, b.length());
			}
		}
		cs_ = b.toCharArray();
	}
	
	public void installReverseComp(String b, int sz) {
		for(int i = 0; i < sz; i++)
			cs_[i] = (b.charAt(sz-i-1) == 4 ? (char)4 : (char)(b.charAt(sz-i-1) ^ 3));
	}
	
	public void reverseComp() {
		for(int i = 0; i < (cs_.length >> 1); i++) {
			char tmp1 = (char) (cs_[i] == 4 ? 4 : cs_[i] ^ 3);
			char tmp2 = (char) (cs_[cs_.length-i-1] == 4 ? 4 : cs_[cs_.length-i-1] ^ 3);
			this.cs_[i] = tmp2;
			this.cs_[cs_.length-i-1] = tmp1;
		}
		// Do middle element iff there are an odd number
		if((this.cs_.length & 1) != 0) {
			char tmp = this.cs_[this.cs_.length >> 1];
			tmp = (char) (tmp == 4 ? 4 : tmp ^ 3);
			this.cs_[this.cs_.length >> 1] = tmp;
		}
	}
	
	public void installChars(String b, int sz) {
		for(int i = 0; i < sz; i++) {
			cs_[i] = Alphabet.asc2dna[(int)b.charAt(i)];
		}
	}
	
	public void installColors(String b, int sz) {
		for(int i = 0; i < sz; i++) {
			cs_[i] = Alphabet.asc2col[(int)b.charAt(i)];
		}
	}
	
	public void setChar(char c, int idx) {
		cs_[idx] = Alphabet.asc2dna[(int)c];
	}
	
	public void appendChar(char c) {
		String tmp = new String(cs_) + Alphabet.asc2dna[(int)c];
		cs_ = tmp.toCharArray();
	}
	
	public char toChar(int idx) {
		return "ACGTN".charAt((int)cs_[idx]);
	}
	
	public char toColor(int idx) {
		int c = (int) cs_[idx];
		return "0123".charAt(c);
	}
	
	public int nwords() {
		return (cs_.length + 15) >> 4;
	}
	
	public char windowGetDna(int i, boolean fw, int depth, int len) {
		if(fw)	return cs_[depth+i];
		return Alphabet.compDna(cs_[depth+len-i-1]);
	}
	
	public String randSubstr(String o, RandomSource rnd, int len, boolean watson, boolean crick) {
		char dst[] = o.toCharArray();
		int poss = cs_.length - len + 1;
		int rndoff = (int)(rnd.nextU32() % poss);
		boolean fw;
		if     (watson && !crick) fw = true;
		else if(!watson && crick) fw = false;
		else {
			fw = rnd.nextBool();
		}
		if(fw) {
			// Install Watson substring
			for(int i = 0; i < len; i++) {
				dst[i] = cs_[i + rndoff];
			}
		} else {
			// Install Crick substring
			for(int i = 0; i < len; i++) {
				dst[i] = Alphabet.maskcomp[(int)cs_[i + rndoff + (len - i - 1)]];
			}
		}
		return new String(dst);
	}
}
