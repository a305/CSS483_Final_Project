package com.uwb.bt2j.indexer.util;

import com.uwb.bt2j.indexer.types.EList;

public class StringUtils<S, T> {
	
	public static EList<String> tokenize(String s, String delims, int max) {
		//string::size_type lastPos = s.find_first_not_of(delims, 0);
		int lastPos = 0;
		EList<String> ss = new EList();
		String[] tokens = s.substring(lastPos).split(delims);
		int pos = tokens.length > 1 ? tokens[0].length() : -1;
		
		while (pos != -1 || lastPos != -1) {
			ss.push_back(s.substring(lastPos, pos - lastPos));
			String[] t = s.substring(pos).split(delims);
			lastPos = t.length > 1 ? s.indexOf(t[0]) : -1;
			pos = s.indexOf(delims, lastPos);
			if(ss.size() == (max - 1)) {
				pos = -1;
			}
		}
		
		return ss;
	}
	
	public static EList<String> tokenize(String s, String delim) {
		String[] tokens = s.split(delim);
		EList<String> ss = new EList();
		for (String u : tokens) {
			ss.push_back(u);
		}
		return ss;
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
	
	public static String reverse(String in) {
		char ret[] = in.toCharArray();
		int back = ret.length - 1;
		for(int i = 0; i < ret.length / 2; i++) {
			char tmp = ret[i];
			back -= i;
			ret[i] = ret[back];
			ret[back] = tmp;
		}
		return new String(ret);
	}
	
	public static String fill(char el, int n) {
		String ret = "";
		for(int i = 0; i < n; i++)
			ret += el;
		return ret;
	}
	
	public static String set(String o, char c, int i) {
		char ret[] = o.toCharArray();
		ret[i] = c;
		return new String(ret);
	}
	
	public static String remove(String o, int i) {
		return o.substring(0,i) + o.substring(i + 1,o.length());
	}
	
	public static String insert(String o, char c, int i) {
		return o.substring(0,i) + c + o.substring(i,o.length());
	}
	
	public static String trimBegin(String o, int n){
		return o.substring(n,o.length());
	}
	
	public static String trimEnd(String o, int n){
		return o.substring(0,o.length() - n);
	}
}
