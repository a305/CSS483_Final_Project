package com.uwb.bt2j.aligner.seed;

public class Constraint {
	public Constraint exact() {
		Constraint c;
		c.edits = c.mms = c.ins = c.dels = c.penalty = 0;
		return c;
	}
	
	public Constraint penaltyBased(int pen) {
		Constraint c;
		c.penalty = pen;
		return c;
	}
	
	public Constraint penaltyFuncBased(SimpleFunc f) {
		Constraint c;
		c.penFunc = f;
		return c;
	}
	
	public Constraint mmBased(int mms) {
		Constraint c;
		c.mms = mms;
		c.edits = c.dels = c.ins = 0;
		return c;
	}
	
	public Constraint editBased(int edits) {
		Constraint c;
		c.edits = edits;
		c.dels = c.ins = c.mms = 0;
		return c;
	}
}
