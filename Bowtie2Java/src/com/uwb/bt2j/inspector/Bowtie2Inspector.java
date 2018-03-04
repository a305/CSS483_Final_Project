package com.uwb.bt2j.inspector;

import java.io.IOException;
import java.io.OutputStream;

import com.uwb.bt2j.indexer.BitPairReference;
import com.uwb.bt2j.indexer.EList;

class Bowtie2Inspector {
	public static boolean showVersion = false;
	public static int verbose = 0;
	public static int names_only = 0;
	public static int summarize_only = 0;
	public static int across = 60;
	public static boolean refFromEbwt = false;
	public static String wrapper;
	public static final String short_options = "vhnsea:";
	public static String argv0;
	
	public enum ARGS {
		ARG_VERSION(256),
		ARG_WRAPPER(257),
		ARG_USAGE(258);
		
		private int x;
		ARGS(int y){this.x = y;}
	};
	
	public static void printUsage(OutputStream out) {
		out.write(("Bowtie 2 version " + BOWTIE2_VERSION + " by Ben Langmead (langmea@cs.jhu.edu, www.cs.jhu.edu/~langmea)" + "\n"
		+ "Usage: bowtie2-inspect [options]* <bt2_base>" + "\n"
		+ "  <bt2_base>         bt2 filename minus trailing .1." + gEbwt_ext + "/.2." + gEbwt_ext + "\n"
		+ "\n"
		+ "  By default, prints FASTA records of the indexed nucleotide sequences to" + "\n"
		+ "  standard out.  With -n, just prints names.  With -s, just prints a summary of" + "\n"
		+ "  the index parameters and sequences." + "\n"
		+ "\n"
		+ "Options:" + "\n").getBytes());
		if(wrapper == "basic-0") {
			out.write(("  --large-index      force inspection of the 'large' index, even if a" + "\n"
				+ "                     'small' one is present." + "\n").getBytes());
		}
		out.write(("  -a/--across <int>  Number of characters across in FASTA output (default: 60)" + "\n"
		+ "  -n/--names         Print reference sequence names only" + "\n"
		+ "  -s/--summary       Print summary incl. ref names, lengths, index properties" + "\n"
		+ "  -v/--verbose       Verbose output (for debugging)" + "\n"
		+ "  -h/--help          print detailed description of tool and its options" + "\n"
		+ "  --help             print this usage message" + "\n").getBytes());
		;
		if(wrapper.empty()) {
			System.err.print("\n"
			     + "*** Warning ***" + "\n"
				 + "'boowtie2-inspect' was run directly.  It is recommended "
				 + "to use the wrapper script instead."
				 + "\n" + "\n");
		}
	}
	
	public static int parseInt(int lower, String errmsg) {
		long l = Long.parseLong(optarg);

			if (l < lower) {
				System.err.println(errmsg);
				printUsage(System.err);
			}
			return (int)l;

		System.err.println( errmsg);
		printUsage(System.err);
		return -1;
	}
	
	public static void parseOptions(String[] args) {
		int option_index = 0;
		int next_option;
		do {
			next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
			switch (next_option) {
				case ARG_WRAPPER:
					wrapper = optarg;
					break;
				case ARG_USAGE:
				case 'h':
					printUsage(System.out);
					break;
				case 'v': verbose = true; break;
				case ARG_VERSION: showVersion = true; break;
				case 'e': refFromEbwt = true; break;
				case 'n': names_only = true; break;
				case 's': summarize_only = true; break;
				case 'a': across = parseInt(-1, "-a/--across arg must be at least 1"); break;
				case -1: break; /* Done with options. */
				case 0:
					if (long_options[option_index].flag != 0)
						break;
				default:
					printUsage(System.err);
			}
		} while(next_option != -1);
	}
	
	public static void printFastaRecord(OutputStream fout, String defline, String seq) throws IOException {
		fout.write('>');
		fout.write(defline.getBytes());

		if(across > 0) {
			int i = 0;
			while (i + across < seq.length())
			{
				fout.write( seq.subString(i, across).getBytes());
				i += across;
			}
			if (i < seq.length())
				fout.write(seq.subString(i).getBytes());
		} else {
			fout.write(seq.getBytes());
		}
	}
	
	public static void printRefSequence(
			OutputStream fout,
			BitPairReference ref,
			String name,
			int refi,
			int len) {
		boolean newlines = across > 0;
		int myacross = across > 0 ? across : 60;
		int incr = myacross * 1000;
		double buf = (incr + 128)/4;
		fout.write((">" + name + "\n").getBytes());
		for(int i = 0; i < len; i += incr) {
			int amt = Math.min(incr, len-i);
			int off = ref.getStretch(buf, refi, i, amt);
			short cb = (short) (buf + off);
			for(int j = 0; j < amt; j++) {
				if(newlines && j > 0 && (j % myacross) == 0) fout.write('\n');
				fout.write("ACGTN".charAt(cb[j]));
			}
			fout.write('\n');
		}
	}
	
	public static void printRefSequences(
			OutputStream fout,
			boolean color,
			EList<String> refnames,
			long plen,
			String adjustedEbwtFileBase) {
		BitPairReference ref = new BitPairReference(
				adjustedEbwtFileBase, // input basename
				color,                // true -> expect colorspace reference
				false,                // sanity-check reference
				null,                 // infiles
				null,                 // originals
				false,                // infiles are sequences
				false,                // memory-map
				false,                // use shared memory
				false,                // sweep mm-mapped ref
				verbose,              // be talkative
				verbose);             // be talkative at startup
			for(int i = 0; i < ref.numRefs(); i++) {
				print_ref_sequence(
					fout,
					ref,
					refnames[i],
					i,
					plen[i] + (color ? 1 : 0));
			}
	}
	
	public static void printIndexSequences(OutputStream fout, Ebwt ebwt) {
		EList<String> refnames = (ebwt.refnames());

		TStr cat_ref;
		ebwt.restore(cat_ref);

		long curr_ref = OFF_MASK;
		String curr_ref_seq = "";
		long curr_ref_len = OFF_MASK;
		long last_text_off = 0;
		int orig_len = cat_ref.length();
		long tlen = OFF_MASK;
		boolean first = true;
		for(int i = 0; i < orig_len; i++) {
			long tidx = OFF_MASK;
			long textoff = OFF_MASK;
			tlen = OFF_MASK;
			boolean straddled = false;
			ebwt.joinedToTextOff(1 /* qlen */, (long)i, tidx, textoff, tlen, true, straddled);

			if (tidx != OFF_MASK && textoff < tlen)
			{
				if (curr_ref != tidx)
				{
					if (curr_ref != OFF_MASK)
					{
						// Add trailing gaps, if any exist

						print_fasta_record(fout, refnames.get((int) curr_ref), curr_ref_seq);
					}
					curr_ref = tidx;
					curr_ref_seq = "";
					curr_ref_len = tlen;
					last_text_off = 0;
					first = true;
				}

				long textoff_adj = textoff;
				if(first && textoff > 0) textoff_adj++;

				curr_ref_seq.push_back("ACGT".charAt(cat_ref[i]));
				last_text_off = textoff;
				first = false;
			}
		}
		if (curr_ref < refnames.size())
		{
			print_fasta_record(fout, refnames.get(curr_ref), curr_ref_seq);
		}
	}
	
	public static void printIndexSequencesNames(String fname, OutputStream fout) {
		EList<String> p_refnames;
		readEbwtRefnames(fname, p_refnames);
		for(int i = 0; i < p_refnames.size(); i++) {
			System.out.println(p_refnames[i]);
		}
	}
	
	public static void printIndexSummary() {
		
	}
	
	public static void driver(String ebwtFileBase, String query) {
		// Adjust
		String adjustedEbwtFileBase = adjustEbwtBase(argv0, ebwtFileBase, verbose);

		if (names_only) {
			print_index_sequence_names(adjustedEbwtFileBase, System.out);
		} else if(summarize_only) {
			print_index_summary(adjustedEbwtFileBase, System.out);
		} else {
			// Initialize Ebwt object
			boolean color = readEbwtColor(adjustedEbwtFileBase);
			Ebwt ebwt = new Ebwt(
				adjustedEbwtFileBase, 
				color,                // index is colorspace
				-1,                   // don't care about entire-reverse
				true,                 // index is for the forward direction
				-1,                   // offrate (-1 = index default)
				0,                    // offrate-plus (0 = index default)
				false,                // use memory-mapped IO
				false,                // use shared memory
				false,                // sweep memory-mapped memory
				true,                 // load names?
				true,                 // load SA sample?
				true,                 // load ftab?
				true,                 // load rstarts?
				verbose,              // be talkative?
				verbose,              // be talkative at startup?
				false,                // pass up memory exceptions?
				false);               // sanity check?
			// Load whole index into memory
			if(refFromEbwt) {
				ebwt.loadIntoMemory(
					-1,     // color
					-1,     // need entire reverse
					true,   // load SA sample
					true,   // load ftab
					true,   // load rstarts
					true,   // load names
					false); // verbose
				print_index_sequences<SString<Character> >(System.out, ebwt);
			} else {
				EList<String> refnames;
				readEbwtRefnames(adjustedEbwtFileBase, refnames);
				print_ref_sequences(
					System.out,
					readEbwtColor(ebwtFileBase),
					refnames,
					ebwt.plen(),
					adjustedEbwtFileBase);
			}
			// Evict any loaded indexes from memory
			if(ebwt.isInMemory()) {
				ebwt.evictFromMemory();
			}
		}
	}
	
	public static void main(String[] args){
		try {
			String ebwtFile;  // read serialized Ebwt from this file
			String query;   // read query String(s) from this file
			EList<String> queries;
			String outfile; // write query results to this file
			argv0 = args[0];
			parseOptions(args);
			if(showVersion) {
				System.out.println( argv0 + " version " + BOWTIE2_VERSION );
				System.out.println( "Built on " + BUILD_HOST);
				System.out.println( BUILD_TIME);
				System.out.println( "Compiler: " + COMPILER_VERSION);
				System.out.println( "Options: " + COMPILER_OPTIONS);
			}

			// Get input filename
			if(optind >= args.length) {
				System.err.println("No index name given!");
				printUsage(System.err);
			}
			ebwtFile = args[optind++];

			// Optionally summarize
			if(verbose) {
				System.out.println( "Input ebwt file: \"" + ebwtFile + "\"");
				System.out.println( "Output file: \"" + outfile + "\"");
				System.out.println( "Local endianness: " + (currentlyBigEndian()? "big":"little"));
			}
			driver(ebwtFile, query);
		} catch(Exception e) {
			System.err.println("Error: Encountered exception: '" + e + "'");
			System.err.print("Command: ");
			for(int i = 0; i < args.length; i++) System.err.print(args[i] + " ");
			System.err.println();
		}
	}
}