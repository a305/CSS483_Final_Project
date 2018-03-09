package com.uwb.bt2j.indexer.util;

public class BitPack {
	public static void pack_2b_in_8b(int two, int eight, int off) {
		eight |= (two << (off*2));
	}
	
	public static int unpack_2b_from_8b(int eight, int off) {
		return ((eight >> (off*2)) & 0x3);
	}
	
	public static void pack_2b_in_32b(int two, int thirty2, int off) {
		thirty2 |= (two << (off*2));
	}
	
	public static int unpack_2b_from_32b(int thirty2, int off) {
		return ((thirty2 >> (off*2)) & 0x3);
	}
}
