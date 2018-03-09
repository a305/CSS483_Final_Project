package com.uwb.bt2j.indexer.util;

public class Alphabet {
  public static char asc2dnacat[] = {
	/*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0,
	       /*                                        - */
	/*  48 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  64 */ 0, 1, 2, 1, 2, 0, 0, 1, 2, 0, 0, 2, 0, 2, 2, 0,
	       /*    A  B  C  D        G  H        K     M  N */
	/*  80 */ 0, 0, 2, 2, 1, 0, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0,
	       /*       R  S  T     V  W  X  Y */
	/*  96 */ 0, 1, 2, 1, 2, 0, 0, 1, 2, 0, 0, 2, 0, 2, 2, 0,
	       /*    a  b  c  d        g  h        k     m  n */
	/* 112 */ 0, 0, 2, 2, 1, 0, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0,
	       /*       r  s  t     v  w  x  y */
	/* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  };
  
  public static char mask2dna[] = {
	'?', // 0
	'A', // 1
	'C', // 2
	'M', // 3
	'G', // 4
	'R', // 5
	'S', // 6
	'V', // 7
	'T', // 8
	'W', // 9
	'Y', // 10
	'H', // 11
	'K', // 12
	'D', // 13
	'B', // 14
	'N', // 15 (inclusive N)
	'N'  // 16 (exclusive N)
  };
  
  public static char mask2iupac[] = {
	'A', // 0001
	'C', // 0010
	'M', // 0011
	'G', // 0100
	'R', // 0101
	'S', // 0110
	'V', // 0111
	'T', // 1000
	'W', // 1001
	'Y', // 1010
	'H', // 1011
	'K', // 1100
	'D', // 1101
	'B', // 1110
	'N', // 1111
  };
  
  public static char asc2dnamask[] = {
	/*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  48 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  64 */ 0, 1,14, 2,13, 0, 0, 4,11, 0, 0,12, 0, 3,15, 0,
	       /*    A  B  C  D        G  H        K     M  N */
	/*  80 */ 0, 0, 5, 6, 8, 0, 7, 9, 0,10, 0, 0, 0, 0, 0, 0,
	       /*       R  S  T     V  W     Y */
	/*  96 */ 0, 1,14, 2,13, 0, 0, 4,11, 0, 0,12, 0, 3,15, 0,
	       /*    a  b  c  d        g  h        k     m  n */
	/* 112 */ 0, 0, 5, 6, 8, 0, 7, 9, 0,10, 0, 0, 0, 0, 0, 0,
	       /*       r  s  t     v  w     y */
	/* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  };
  
  public static int mask2popcnt[] = {
	0, 1, 1, 2, 1, 2, 2, 3,
	1, 2, 2, 3, 2, 3, 3, 4,
	1, 2, 2, 3, 2, 3, 3, 4,
	2, 3, 3, 4, 3, 4, 4, 5
  };
  
  public static char asc2dna[] = {
	/*   0 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	/*  16 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	/*  32 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	/*                                               - */
	/*  48 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	/*  64 */ 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
	       /*    A  B  C  D        G  H        K     M  N */
	/*  80 */ 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	       /*       R  S  T  U  V  W     Y */
	/*  96 */ 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
	       /*    a  b  c  d        g  h        k     m  n */
	/* 112 */ 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	       /*       r  s  t  u  v  w     y */
	/* 128 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	/* 144 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	/* 160 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	/* 176 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	/* 192 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	/* 208 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	/* 224 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	/* 240 */ 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  };
  
  public static char asc2col[] = {
	/*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 4, 0,
	       /*                                        -  . */
	/*  48 */ 0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	       /* 0  1  2  3 */
	/*  64 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  80 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  96 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 112 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  };
  
  
  public static char maskcomp[] = {
			0,  // 0000 (!) -> 0000 (!)
			8,  // 0001 (A) -> 1000 (T)
			4,  // 0010 (C) -> 0100 (G)
			12, // 0011 (M) -> 1100 (K)
			2,  // 0100 (G) -> 0010 (C)
			10, // 0101 (R) -> 1010 (Y)
			6,  // 0110 (S) -> 0110 (S)
			14, // 0111 (V) -> 1110 (B)
			1,  // 1000 (T) -> 0001 (A)
			9,  // 1001 (W) -> 1001 (W)
			5,  // 1010 (Y) -> 0101 (R)
			13, // 1011 (H) -> 1101 (D)
			3,  // 1100 (K) -> 0011 (M)
			11, // 1101 (D) -> 1011 (H)
			7,  // 1110 (B) -> 0111 (V)
			15 // 1111 (N) -> 1111 (N)
		  };
  
  public static char dnacomp[] = {3, 2, 1, 0, 4};
  
  public static byte dinuc2color[][] = {
	/* A */ {0, 1, 2, 3, 4},
	/* C */ {1, 0, 3, 2, 4},
	/* G */ {2, 3, 0, 1, 4},
	/* T */ {3, 2, 1, 0, 4},
	/* N */ {4, 4, 4, 4, 4}
  };
  
  public static final String iupacs = "!ACMGRSVTWYHKDBN!acmgrsvtwyhkdbn";
  
  public static Boolean isUnambigNuc(char c) {
	  return asc2dnacat[(int) c] == 1;
  }
  
  public static char comp(char c) {
	  switch(c) {
		case 'a': return 't';
		case 'A': return 'T';
		case 'c': return 'g';
		case 'C': return 'G';
		case 'g': return 'c';
		case 'G': return 'C';
		case 't': return 'a';
		case 'T': return 'A';
		default: return c;
	  }
  }
  
  public static char compDna(int c) {
		return dnacomp[c];
  }
  
  public static void decodeNuc(char c , int num, int[] alts) {
	  switch(c) {
		case 'A': alts[0] = 0; num = 1; break;
		case 'C': alts[0] = 1; num = 1; break;
		case 'G': alts[0] = 2; num = 1; break;
		case 'T': alts[0] = 3; num = 1; break;
		case 'M': alts[0] = 0; alts[1] = 1; num = 2; break;
		case 'R': alts[0] = 0; alts[1] = 2; num = 2; break;
		case 'W': alts[0] = 0; alts[1] = 3; num = 2; break;
		case 'S': alts[0] = 1; alts[1] = 2; num = 2; break;
		case 'Y': alts[0] = 1; alts[1] = 3; num = 2; break;
		case 'K': alts[0] = 2; alts[1] = 3; num = 2; break;
		case 'V': alts[0] = 0; alts[1] = 1; alts[2] = 2; num = 3; break;
		case 'H': alts[0] = 0; alts[1] = 1; alts[2] = 3; num = 3; break;
		case 'D': alts[0] = 0; alts[1] = 2; alts[2] = 3; num = 3; break;
		case 'B': alts[0] = 1; alts[1] = 2; alts[2] = 3; num = 3; break;
		case 'N': alts[0] = 0; alts[1] = 1; alts[2] = 2; alts[3] = 3; num = 4; break;
		default: {
			System.err.println("Bad IUPAC code: " + c + ", (int: " + (int)c + ")");
		}
		}
  }
  
  public static void setIupacsCat(byte cat) {
		asc2dnacat[(int)'B'] = asc2dnacat[(int)'b'] =
		asc2dnacat[(int)'D'] = asc2dnacat[(int)'d'] =
		asc2dnacat[(int)'H'] = asc2dnacat[(int)'h'] =
		asc2dnacat[(int)'K'] = asc2dnacat[(int)'k'] =
		asc2dnacat[(int)'M'] = asc2dnacat[(int)'m'] =
		asc2dnacat[(int)'N'] = asc2dnacat[(int)'n'] =
		asc2dnacat[(int)'R'] = asc2dnacat[(int)'r'] =
		asc2dnacat[(int)'S'] = asc2dnacat[(int)'s'] =
		asc2dnacat[(int)'V'] = asc2dnacat[(int)'v'] =
		asc2dnacat[(int)'W'] = asc2dnacat[(int)'w'] =
		asc2dnacat[(int)'X'] = asc2dnacat[(int)'x'] =
		asc2dnacat[(int)'Y'] = asc2dnacat[(int)'y'] = cat;
  }
}
