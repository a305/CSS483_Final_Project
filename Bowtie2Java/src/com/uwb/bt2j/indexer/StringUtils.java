package com.uwb.bt2j.indexer;

public class StringUtils<S, T> {
	public double len(S s){
		return 
	}
	
	public Boolean equals(S s, T t) {
		
	}
	
	public Boolean notEquals(S s, T t) {
		
	}
	
	public Boolean suffixUpToEquals() {
		
	}
	
	public Boolean suffixUpToNotEquals() {
		
	}
	
	public Boolean lessThan() {
		
	}
	
	public Boolean greaterThan() {
		
	}
	
	public Boolean lessThanOrEquals() {
		
	}
	
	public Boolean greaterThanOrEquals() {
		
	}
	
	public Boolean suffixLessThan() {
		
	}
	
	public Boolean suffixGreaterThan() {
		
	}
	
	public Boolean suffixLessThanOrEquals() {
		
	}
	
	public Boolean suffixGreaterThanOrEquals() {
		
	}
	
	public Boolean suffixUpToLessThan() {
		
	}
	
	public Boolean prefixLessThan() {
		
	}
	
	public Boolean prefixGreaterThan() {
		
	}
	
	public Boolean prefixLessThanOrEquals() {
		
	}
	
	public Boolean prefixGreaterThanOrEquals() {
		
	}
	
	public static int hash_string(String s) {
		int ret = 0;
		int a = 63689;
		int b = 378551;
		for(int i = 0; i < s.length(); i++) {
			ret = (ret * a) + (int)s.charAt(i);
			if(a == 0) {
				a += b;
			} else {
				a *= b;
			}
			if(a == 0) {
				a += b;
			}
		}
		return ret;
	}
}