package com.uwb.bt2j.aligner;

public class RandomSource {
	public static final int DEFUALT_A = 1664525;
	public static final int DEFUALT_C = 1013904223;
	private int a;
	private int c;
	private int last;
	private int lastOff;
	private boolean inited_;
	
	public RandomSource() {
		a = DEFUALT_A;
		c = DEFUALT_C;
		inited_ = false;
	}
	
	public RandomSource(int _last) {
		a = DEFUALT_A;
		c = DEFUALT_C;
		last = _last;
		inited_ = true;
	}
	
	public RandomSource(int _a, int _c) {
		a = _a;
		c = _c;
		inited_ = false;
	}
	
	public void init(int seed) {
		last = seed;
		inited_ = true;
		lastOff = 30;
	}
	
	public int nextU32() {
		int ret;
		last = a * last + c;
		ret = last >> 16;
		last = a * last + c;
		ret ^= last;
		lastOff = 0;
		return ret;
	}
	
	public long nextU64() {
		long first = nextU32();
		first = first << 32;
		long second = nextU32();
		return first | second;
	}
	
	public double nextSizeT() {
		return nextU32();
	}
	
	public int nextU32Range(int lo, int hi) {
		int ret = lo;
		if(hi > lo) {
			ret += (nextU32() % (hi-lo+1));
		}
		return ret;
	}
	
	public int nextU2() {
		if(lastOff > 30) {
			nextU32();
		}
		int ret = (last >> lastOff) & 3;
		lastOff += 2;
		return ret;
	}
	
	public int nextBool() {
		if(lastOff > 31) {
			nextU32();
		}
		int ret = (last >> lastOff) & 1;
		lastOff++;
		return ret;
	}
	
	public int nextFromProbs(float[] weights, double numWeights) {
		float f = nextFloat();
		float tot = 0.0f; // total weight seen so far
		for(int i = 0; i < numWeights; i++) {
			tot += weights[i];
			if(f < tot) return i;
		}
		return (int)(numWeights-1);
	}
	
	public float nextFloat() {
		return (float)nextU32() / (float)0xffffffff;
	}
	
	public static int nextU32(int last, int a, int c) {
		return (a * last) + c;
	}
	
	public int currentA() {
		return a;
	}
	
	public int currentC() {
		return c;
	}
	
	public int currentLast() {
		return last;
	}
}
