package com.uwb.bt2j.aligner.seed;

public class InstantiatedSeed {
	public EList<int> steps;
	public EList<Pair<int,int>> zones;
	public BTDnaString seq;
	public BTString qual;
	public Constraint cons[];
	public Constraint overall;
	public int maxjump;
	public int seedoff;
	public int seedoffidx;
	public boolean fw;
	public boolean nfiltered;
	public Seed s;
	
	public InstantiatedSeed() {
		steps, zones = 5;
	}
}
