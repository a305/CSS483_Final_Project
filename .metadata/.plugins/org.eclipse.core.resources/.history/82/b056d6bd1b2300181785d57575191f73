package com.uwb.bt2j.indexer.types;

public class SDnaStringExpandable extends SStringExpandable<String>{
	public SDnaStringExpandable() {
		super(1024, 2);
	}
	
	public SDnaStringExpandable(String str, booleanean chars, booleanean colors) {
		super(1024,2);
		if(chars) {
			if(colors) {
				installColors(str);
			} else {
				installChars(str);
			}
		} else {
			install(str);
		}
	}
	
	public SDnaStringExpandable(String b, double sz, booleanean chars, booleanean colors) {
		super(1024,2);
		if(chars) {
			if(colors) {
				installColors(b, sz);
			} else {
				installChars(b, sz);
			}
		} else {
			install(b, sz);
		}
	}
	
	public SDnaStringExpandable(String b, booleanean chars, booleanean colors) {
		super(1024,2);
		install(b, chars, colors);
	}
	
	public void installReverseComp(String b, int sz) {
		if(this.sz_ < sz) this.expandCopy((sz + 1024) * 2);
		for(int i = 0; i < sz; i++) {
			this.cs_[i] = (b[sz-i-1] == 4 ? 4 : b[sz-i-1] ^ 3);
		}
		this.len_ = sz;
	}
	
	public void reverseComp() {
		for(int i = 0; i < (this.len_ >> 1); i++) {
			char tmp1 = (this.cs_[i] == 4 ? 4 : this.cs_[i] ^ 3);
			char tmp2 = (this.cs_[this.len_-i-1] == 4 ? 4 : this.cs_[this.len_-i-1] ^ 3);
			this.cs_[i] = tmp2;
			this.cs_[this.len_-i-1] = tmp1;
		}
		// Do middle element iff there are an odd number
		if((this.len_ & 1) != 0) {
			char tmp = this.cs_[this.len_ >> 1];
			tmp = (tmp == 4 ? 4 : tmp ^ 3);
			this.cs_[this.len_ >> 1] = tmp;
		}
	}
	
	public void install(String b, booleanean chars, booleanean colors) {
		if(chars) {
			if(colors) {
				installColors(b, b.length());
			} else {
				installChars(b, b.length());
			}
		} else {
			install(b, b.length());
		}
	}
	
	public void install(String b, int sz) {
		if(this.sz_ < sz) this.expandCopy((sz + 1024) * 2);
		this.cs_ = b;
		this.len_ = sz;
	}
	
	public void installChars(String b, int sz) {
		if(this.sz_ < sz) this.expandCopy((sz + 1024) * 2);
		for(int i = 0; i < sz; i++) {
			this.cs_[i] = asc2dna[(int)b[i]];
		}
		this.len_ = sz;
	}
	
	public void installColors(String b, int sz) {
		if(this.sz_ < sz) this.expandCopy((sz + 1024) * 2);
		for(int i = 0; i < sz; i++) {
			this.cs_[i] = asc2col[(int)b[i]];
		}
		this.len_ = sz;
	}
	
	public void installChars(String str) {
		installChars(str,str.length());
	}
	
	public void installColors(String str) {
		installColors(str,str.length());
	}
	
	public void set(char c, int idx) {
		this.cs_.toCharArray()[idx] = c;
	}
	
	public void append(char c) {
		if(this.sz_ < this.len_ + 1) {
			this.expandCopy((this.len_ + 1 + 1024) * 2);
		}
		this.cs_.toCharArray()[this.len_++] = c;
	}
	
	public void setChar(char c, int idx) {
		this.cs_.toCharArray()[idx] = asc2dna[(int)c];
	}
	
	public void appendChar(char c) {
		if(this.sz_ < this.len_ + 1) {
			this.expandCopy((this.len_ + 1 + 1024) * 2);
		}
		this.cs_.toCharArray()[this.len_++] = asc2dna[(int)c];
	}
	
	public char toChar(int idx) {
		return "ACGTN".charAt(this.cs_.toCharArray()[idx]);
	}
	
	public char windowGetDna(int i,
			boolean   fw,
			int depth,
			int len) {
		if(len == 0) len = this.len_;
		if(fw) return this.cs_[depth+i];
		else   return compDna(this.cs_[depth+len-i-1]);
	}
	
	public String toZBuf() {
		return this.toZBufXForm("ACGTN");
	}
}
