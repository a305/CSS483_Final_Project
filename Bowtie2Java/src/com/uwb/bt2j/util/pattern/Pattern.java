package com.uwb.bt2j.util.pattern;

public class Pattern {
	public static double genRandSeed(BTDnaString qry, BTString qual, BTString name, double seed) {
		double rseed = (seed + 101) * 59 * 61 * 67 * 71 * 73 * 79 * 83;
		double qlen = qry.length();
		
		for(int i = 0; i < qlen; i++) {
			int p = (int)qual[i];
			int off = ((i & 3) << 3);
			rseed ^= (double)(p << off);
		}
		
		for(int i = 0; i < qlen; i++) {
			int p = (int)qual[i];
			int off = ((i & 3) << 3);
			rseed ^= (double)(p << off);
		}
		
		double namelen = name.length();
		for(int i = 0; i < namelen; i++) {
			int p = (int)name[i];
			if(p == '/')	break;
			int off = ((i & 3) << 3);
			rseed ^= (double)(p << off);
		}
		
		return rseed;
	}	
	
	public void wrongQualityFormat(BTString read_name) {
		System.err.println("Error: Encountered one or more spaces while parsing the quality "
				 + "string for read " + read_name + ".  If this is a FASTQ file "
				 + "with integer (non-ASCII-encoded) qualities, try re-running with "
				 + "the --integer-quals option.");
	}
	
	public void tooFewQualities(BTString read_name) {
		System.err.println("Error: Read " + read_name + " has more read characters than "
				 + "quality values.");
	}
	
	public void tooManyQualities(BTString read_name) {
		System.err.println("Error: Read " + read_name + " has more quality values than read "
		 + "characters.");
	}
}