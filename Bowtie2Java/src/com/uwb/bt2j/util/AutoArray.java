package com.uwb.bt2j.util;

import java.lang.reflect.Array;

public class AutoArray<T>{
	private int cat_;
	private T[] t_;
	private int sz_;
	
  public AutoArray(int sz, int cat) {
    cat_ = cat;
    t_ = (T[])Array.newInstance(T[],T);
    sz_ = sz;
  }
  
  public int size() {
	  return sz_;
  }
}
