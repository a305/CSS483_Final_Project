package com.uwb.bt2j.indexer;

import java.io.File;

import com.uwb.bt2j.indexer.types.EList;

class FASTAParser<T, U, V> {
  public static void parseFastaLens(
		  File infile,   // filename
			EList<Integer> namelens, // destination for fasta name lengths
			EList<Integer> seqlens
		  ) {
		if(infile == null) {
			System.err.println("Could not open sequence file");
		}
		FileBuf fb = new FileBuf(infile);
		while(!fb.eof()) {
			namelens.expand();
			namelens.insert(0, namelens.back());
			seqlens.insert(0, seqlens.back());
			fb.parseFastaRecordLength(namelens.back(), seqlens.back());
			if(seqlens.back() == 0) {
				// Couldn't read a record.  We're probably done with this file.
				namelens.pop_back();
				seqlens.pop_back();
				continue;
			}
		}
		fb.close();
  }
  
  public void parseFasta(
		  	File infile,   // filename
			EList<U> names,    // destination for fasta names
			EList<Integer>   namelens, // destination for fasta name lengths
			EList<V>  seqs,     // destination for fasta sequences
			EList<Integer>   seqlens
		  ) {
	    int cur = namelens.size();
		parseFastaLens(infile, namelens, seqlens);
		if(infile == null) {
			System.err.println("Could not open sequence file");
		}
		FileBuf fb = new FileBuf(infile);
		while(!fb.eof()) {
			// Add a new empty record to the end
			names.expand();
			seqs.expand();
			names.insert((T)new char[namelens.get(cur)+1], (int)names.back());
			seqs.insert((T)new char[seqlens[cur]+1], (int)seqs.back());
			fb.parseFastaRecord(names.back(), seqs.back());
			if(seqs.empty()) {
				// Couldn't read a record.  We're probably done with this file.
				names.pop_back();
				seqs.pop_back();
				continue;
			}
		}
		fb.close();
  }
  
  public void parseFastas(
		  EList<T>    infiles,   // filename
			EList<U> names,    // destination for fasta names
			EList<Integer>   namelens, // destination for fasta name lengths
			EList<V>  seqs,     // destination for fasta sequences
			EList<Integer>   seqlens
		  ) {
		for(int i = 0; i < infiles.size(); i++) {
			parseFasta(
				(File)infiles.get(i),
				names,
				namelens,
				seqs,
				seqlens);
		}
  }
}
