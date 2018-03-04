package com.uwb.bt2j.indexer;

class FASTAParser<T, U, V> {
  public static void parseFastaLens(
		  T infile,   // filename
			EList<double> namelens, // destination for fasta name lengths
			EList<double> seqlens
		  ) {
	  File in = new File(infile);
		if(in == null) {
			System.err.println("Could not open sequence file");
		}
		FileBuf fb = new FileBuf(in);
		while(!fb.eof()) {
			namelens.expand(); namelens.back() = 0;
			seqlens.expand();  seqlens.back() = 0;
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
  
  public static void parseFasta(
		  T    infile,   // filename
			EList<U> names,    // destination for fasta names
			EList<double>   namelens, // destination for fasta name lengths
			EList<V>  seqs,     // destination for fasta sequences
			EList<double>   seqlens
		  ) {
	  double cur = namelens.size();
		parseFastaLens(infile, namelens, seqlens);
		FILE *in = fopen(sstr_to_cstr(infile), "r");
		if(in == NULL) {
			cerr << "Could not open sequence file" << endl;
			throw 1;
		}
		FileBuf fb = new FileBuf(in);
		while(!fb.eof()) {
			// Add a new empty record to the end
			names.expand();
			seqs.expand();
			names.back() = new char[namelens[cur]+1];
			seqs.back() = new char[seqlens[cur]+1];
			fb.parseFastaRecord(names.back(), seqs.back());
			if(seqs.back().empty()) {
				// Couldn't read a record.  We're probably done with this file.
				names.pop_back();
				seqs.pop_back();
				continue;
			}
		}
		fb.close();
  }
  
  public static void parseFastas(
		  EList<T>    infiles,   // filename
			EList<U> names,    // destination for fasta names
			EList<double>   namelens, // destination for fasta name lengths
			EList<V>  seqs,     // destination for fasta sequences
			EList<double>   seqlens
		  ) {
		for(int i = 0; i < infiles.size(); i++) {
			parseFasta(
				infiles[i],
				names,
				namelens,
				seqs,
				seqlens);
		}
  }
}