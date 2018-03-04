package com.uwb.bt2j.indexer;

import java.io.OutputStream;

import com.uwb.bt2j.util.Formats.FileFormat;
import com.uwb.bt2j.util.IndexTypes;
import com.uwb.bt2j.util.Tokenize;
import com.uwb.bt2j.util.types.EList;

public class Indexer <T> {
	
	// Build parameters
	public int verbose;
	public static int sanityCheck;
	public static FileFormat format;
	public static long bmax;
	public static long bmaxMultSqrt;
	public static int bmaxDivN;
	public static int dcv;
	public static int noDc;
	public static int entireSA;
	public static int seed;
	public static int showVersion;
	//   Ebwt parameters
	public static int lineRate;
	public static int linesPerSide;
	public static int offRate;
	public static int ftabChars;
	public static int  bigEndian;
	public static boolean nsToAs;    // convert Ns to As
	public static boolean doSaFile;  // make a file with just the suffix array in it
	public static boolean doBwtFile; // make a file with just the BWT string in it
	public static boolean autoMem;
	public static boolean packed;
	public static boolean writeRef;
	public static boolean justRef;
	public static boolean reverseEach;
	public static int nthreads;
	public static String wrapper;
	public static final String short_options = "qraph?nscfl:i:o:t:h:3C";
	public EList<String> filesWritten;
	
	public enum GetOpts {
		ARG_BMAX(256),
		ARG_BMAX_MULT(257),
		ARG_BMAX_DIV(258),
		ARG_DCV(259),
		ARG_SEED(260),
		ARG_CUTOFF(261),
		ARG_PMAP(262),
		ARG_NTOA(263),
		ARG_USAGE(264),
		ARG_REVERSE_EACH(265),
		ARG_SA(266),
		ARG_THREADS(267),
		ARG_WRAPPER(268);
		private int x;
		GetOpts(int y){x = y;}
	}
	
	public static void main(String[] args) {
		if(argc > 2 && strcmp(argv[1], "-A") == 0) {
			String file = args[2];
			ifstream in;
			in.open(file);
			String buf;
			int lastret = -1;
			while(in.getline(buf, 4095)) {
				EList<String> argz = new EList(9);
				argz.push_back(args[0]);
				buf = Tokenize.tokenize(" \t", args);
				for(int i = 0; i < args.length; i++) {
					argz.insert(args[i],i);
				}
				if(args.length == 1) continue;
				lastret = bowtie_build((int)args.length,argz);
			}
			if(lastret == -1) {
				System.err.println( "Warning: No arg strings parsed from " + file);
			}
		}
		string outfile;
		try {
			// Reset all global state, including getopt state
			opterr = optind = 1;
			resetOptions();

			string infile;
			EList<string> infiles(MISC_CAT);

			if(parseOptions(argc, argv)) {
				return 0;
			}
			argv0 = argv[0];
			if(showVersion) {
				cout << argv0 << " version " << string(BOWTIE2_VERSION).c_str() << endl;
				if(sizeof(void*) == 4) {
					cout << "32-bit" << endl;
				} else if(sizeof(void*) == 8) {
					cout << "64-bit" << endl;
				} else {
					cout << "Neither 32- nor 64-bit: sizeof(void*) = " << sizeof(void*) << endl;
				}
				cout << "Built on " << BUILD_HOST << endl;
				cout << BUILD_TIME << endl;
				cout << "Compiler: " << COMPILER_VERSION << endl;
				cout << "Options: " << COMPILER_OPTIONS << endl;
				cout << "Sizeof {int, long, long long, void*, size_t, off_t}: {"
					 << sizeof(int)
					 << ", " << sizeof(long) << ", " << sizeof(long long)
					 << ", " << sizeof(void *) << ", " << sizeof(size_t)
					 << ", " << sizeof(off_t) << "}" << endl;
				return 0;
			}

			// Get input filename
			if(optind >= argc) {
				cerr << "No input sequence or sequence file specified!" << endl;
				printUsage(cerr);
				return 1;
			}
			infile = argv[optind++];

			// Get output filename
			if(optind >= argc) {
				cerr << "No output file specified!" << endl;
				printUsage(cerr);
				return 1;
			}
			outfile = argv[optind++];

			tokenize(infile, ",", infiles);
			if(infiles.size() < 1) {
				cerr << "Tokenized input file list was empty!" << endl;
				printUsage(cerr);
				return 1;
			}

			// Optionally summarize
			if(verbose) {
				cout << "Settings:" << endl
					 << "  Output files: \"" << outfile.c_str() << ".*." + gEbwt_ext + "\"" << endl
					 << "  Line rate: " << lineRate << " (line is " << (1<<lineRate) << " bytes)" << endl
					 << "  Lines per side: " << linesPerSide << " (side is " << ((1<<lineRate)*linesPerSide) << " bytes)" << endl
					 << "  Offset rate: " << offRate << " (one in " << (1<<offRate) << ")" << endl
					 << "  FTable chars: " << ftabChars << endl
					 << "  Strings: " << (packed? "packed" : "unpacked") << endl
					 ;
				if(bmax == OFF_MASK) {
					cout << "  Max bucket size: default" << endl;
				} else {
					cout << "  Max bucket size: " << bmax << endl;
				}
				if(bmaxMultSqrt == OFF_MASK) {
					cout << "  Max bucket size, sqrt multiplier: default" << endl;
				} else {
					cout << "  Max bucket size, sqrt multiplier: " << bmaxMultSqrt << endl;
				}
				if(bmaxDivN == 0xffffffff) {
					cout << "  Max bucket size, len divisor: default" << endl;
				} else {
					cout << "  Max bucket size, len divisor: " << bmaxDivN << endl;
				}
				cout << "  Difference-cover sample period: " << dcv << endl;
				cout << "  Endianness: " << (bigEndian? "big":"little") << endl
					 << "  Actual local endianness: " << (currentlyBigEndian()? "big":"little") << endl
					 << "  Sanity checking: " << (sanityCheck? "enabled":"disabled") << endl;
		#ifdef NDEBUG
				cout << "  Assertions: disabled" << endl;
		#else
				cout << "  Assertions: enabled" << endl;
		#endif
				cout << "  Random seed: " << seed << endl;
				cout << "  Sizeofs: void*:" << sizeof(void*) << ", int:" << sizeof(int) << ", long:" << sizeof(long) << ", size_t:" << sizeof(size_t) << endl;
				cout << "Input files DNA, " << file_format_names[format].c_str() << ":" << endl;
				for(size_t i = 0; i < infiles.size(); i++) {
					cout << "  " << infiles[i].c_str() << endl;
				}
			}
			// Seed random number generator
			srand(seed);
			{
				Timer timer(cout, "Total time for call to driver() for forward index: ", verbose);
				if(!packed) {
					try {
						driver<SString<char> >(infile, infiles, outfile, false, REF_READ_FORWARD);
					} catch(bad_alloc& e) {
						if(autoMem) {
							cerr << "Switching to a packed string representation." << endl;
							packed = true;
						} else {
							throw e;
						}
					}
				}
				if(packed) {
					driver<S2bDnaString>(infile, infiles, outfile, true, REF_READ_FORWARD);
				}
			}
			int reverseType = reverseEach ? REF_READ_REVERSE_EACH : REF_READ_REVERSE;
			srand(seed);
			Timer timer(cout, "Total time for backward call to driver() for mirror index: ", verbose);
			if(!packed) {
				try {
					driver<SString<char> >(infile, infiles, outfile + ".rev", false, reverseType);
				} catch(bad_alloc& e) {
					if(autoMem) {
						cerr << "Switching to a packed string representation." << endl;
						packed = true;
					} else {
						throw e;
					}
				}
			}
			if(packed) {
				driver<S2bDnaString>(infile, infiles, outfile + ".rev", true, reverseType);
			}
			return 0;
		} catch(std::exception& e) {
			cerr << "Error: Encountered exception: '" << e.what() << "'" << endl;
			cerr << "Command: ";
			for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
			cerr << endl;
			deleteIdxFiles(outfile, writeRef || justRef, justRef);
			return 1;
		} catch(int e) {
			if(e != 0) {
				cerr << "Error: Encountered internal Bowtie 2 exception (#" << e << ")" << endl;
				cerr << "Command: ";
				for(int i = 0; i < argc; i++) cerr << argv[i] << " ";
				cerr << endl;
			}
			deleteIdxFiles(outfile, writeRef || justRef, justRef);
			return e;
		}
	}
	
	public static void resetOptions() {
		verbose      = true;  // be talkative (default)
		sanityCheck  = 0;     // do slow sanity checks
		format       = FileFormat.FASTA; // input sequence format
		bmax         = IndexTypes.OFF_MASK; // max blockwise SA bucket size
		bmaxMultSqrt = IndexTypes.OFF_MASK; // same, as multplier of sqrt(n)
		bmaxDivN     = 4;          // same, as divisor of n
		dcv          = 1024;  // bwise SA difference-cover sample sz
		noDc         = 0;     // disable difference-cover sample
		entireSA     = 0;     // 1 = disable blockwise SA
		seed         = 0;     // srandom seed
		showVersion  = 0;     // just print version and quit?
		//   Ebwt parameters
		lineRate     = Ebwt.default_lineRate; // a "line" is 64 or 128 bytes
		linesPerSide = 1;  // 1 64-byte line on a side
		offRate      = 4;  // sample 1 out of 16 SA elts
		ftabChars    = 10; // 10 chars in initial lookup table
		bigEndian    = 0;  // little endian
		nsToAs       = false; // convert reference Ns to As prior to indexing
		doSaFile     = false; // make a file with just the suffix array in it
		doBwtFile    = false; // make a file with just the BWT string in it
		autoMem      = true;  // automatically adjust memory usage parameters
		packed       = false; //
		writeRef     = true;  // write compact reference to .3.gEbwt_ext/.4.gEbwt_ext
		justRef      = false; // *just* write compact reference, don't index
		reverseEach  = false;
	    nthreads     = 1;
		wrapper = "";
	}
	
	public static void printUsage(OutputStream out) {
		out.write(( "Bowtie 2 version " + string(BOWTIE2_VERSION) + " by Ben Langmead (langmea@cs.jhu.edu, www.cs.jhu.edu/~langmea)" + "\n").getBytes());
			String tool_name = "bowtie2-build-s";
			if(wrapper == "basic-0") {
				tool_name = "bowtie2-build";
			}
			
			//               1         2         3         4         5         6         7         8
			//      12345678901234567890123456789012345678901234567890123456789012345678901234567890
			out.write(( "Usage: " + tool_name + " [options]* <reference_in> <bt2_index_base>" + "\n"
			    + "    reference_in            comma-separated list of files with ref sequences" + "\n"
			    + "    bt2_index_base          write " + IndexTypes.gEbwt_ext + " data to files with this dir/basename" + "\n"
				+ "*** Bowtie 2 indexes work only with v2 (not v1).  Likewise for v1 indexes. ***" + "\n"
			    + "Options:" + "\n"
			    + "    -f                      reference files are Fasta (default)" + "\n"
			    + "    -c                      reference sequences given on cmd line (as" + "\n"
				+ "                            <reference_in>)" + "\n").getBytes());
			if(wrapper == "basic-0") {
			out.write(( "    --large-index           force generated index to be 'large', even if ref" + "\n"
				+ "                            has fewer than 4 billion nucleotides" + "\n").getBytes());
			}
			out.write(( "    -a/--noauto             disable automatic -p/--bmax/--dcv memory-fitting" + "\n"
			    + "    -p/--packed             use packed strings internally; slower, less memory" + "\n"
			    + "    --bmax <int>            max bucket sz for blockwise suffix-array builder" + "\n"
			    + "    --bmaxdivn <int>        max bucket sz as divisor of ref len (default: 4)" + "\n"
			    + "    --dcv <int>             diff-cover period for blockwise (default: 1024)" + "\n"
			    + "    --nodc                  disable diff-cover (algorithm becomes quadratic)" + "\n"
			    + "    -r/--noref              don't build .3/.4 index files" + "\n"
			    + "    -3/--justref            just build .3/.4 index files" + "\n"
			    + "    -o/--offrate <int>      SA is sampled every 2^<int> BWT chars (default: 5)" + "\n"
			    + "    -t/--ftabchars <int>    # of chars consumed in initial lookup (default: 10)" + "\n"
		        + "    --threads <int>         # of threads" + "\n"
			    + "    --seed <int>            seed for random number generator" + "\n"
			    + "    -q/--quiet              verbose output (for debugging)" + "\n"
			    + "    -h/--help               print detailed description of tool and its options" + "\n"
			    + "    --usage                 print this usage message" + "\n"
			    + "    --version               print version information and quit" + "\n").getBytes());
			    ;
			if(wrapper.equals("")) {
				System.err.println("\n"
				     + "*** Warning ***" + "\n"
					 + "'" + tool_name + "' was run directly.  It is recommended "
					 + "that you run the wrapper script 'bowtie2-build' instead."
					 + "\n");
			}
	}
	
	public static boolean parseOptions(String args) {
		int option_index = 0;
		int next_option;
		bool bmaxDivNSet = false;
		bool abort = false;
		do {
			next_option = getopt_long(
				argc, const_cast<char**>(argv),
				short_options, long_options, &option_index);
			switch (next_option) {
				case ARG_WRAPPER:
					wrapper = optarg;
					break;
				case 'f': format = FASTA; break;
				case 'c': format = CMDLINE; break;
				case 'p': packed = true; break;
				case 'l':
					lineRate = parseNumber<int>(3, "-l/--lineRate arg must be at least 3");
					break;
				case 'i':
					linesPerSide = parseNumber<int>(1, "-i/--linesPerSide arg must be at least 1");
					break;
				case 'o':
					offRate = parseNumber<int>(0, "-o/--offRate arg must be at least 0");
					break;
				case '3':
					justRef = true;
					break;
				case 't':
					ftabChars = parseNumber<int>(1, "-t/--ftabChars arg must be at least 1");
					break;
				case 'n':
					// all f-s is used to mean "not set", so put 'e' on end
					bmax = 0xfffffffe;
					break;
				case 'h':
				case ARG_USAGE:
					printUsage(cout);
					abort = true;
					break;
				case ARG_BMAX:
					bmax = parseNumber<TIndexOffU>(1, "--bmax arg must be at least 1");
					bmaxMultSqrt = OFF_MASK; // don't use multSqrt
					bmaxDivN = 0xffffffff;     // don't use multSqrt
					break;
				case ARG_BMAX_MULT:
					bmaxMultSqrt = parseNumber<TIndexOffU>(1, "--bmaxmultsqrt arg must be at least 1");
					bmax = OFF_MASK;     // don't use bmax
					bmaxDivN = 0xffffffff; // don't use multSqrt
					break;
				case ARG_BMAX_DIV:
					bmaxDivNSet = true;
					bmaxDivN = parseNumber<uint32_t>(1, "--bmaxdivn arg must be at least 1");
					bmax = OFF_MASK;         // don't use bmax
					bmaxMultSqrt = OFF_MASK; // don't use multSqrt
					break;
				case ARG_DCV:
					dcv = parseNumber<int>(3, "--dcv arg must be at least 3");
					break;
				case ARG_SEED:
					seed = parseNumber<int>(0, "--seed arg must be at least 0");
					break;
				case ARG_REVERSE_EACH:
					reverseEach = true;
					break;
				case ARG_SA:
					doSaFile = true;
					break;
				case ARG_NTOA: nsToAs = true; break;
	            case ARG_THREADS:
	                nthreads = parseNumber<int>(0, "--threads arg must be at least 1");
	                break;
				case 'a': autoMem = false; break;
				case 'q': verbose = false; break;
				case 's': sanityCheck = true; break;
				case 'r': writeRef = false; break;

				case -1: /* Done with options. */
					break;
				case 0:
					if (long_options[option_index].flag != 0)
						break;
				default:
					printUsage(cerr);
					throw 1;
			}
		} while(next_option != -1);
		if(bmax < 40) {
			cerr << "Warning: specified bmax is very small (" << bmax << ").  This can lead to" << endl
			     << "extremely slow performance and memory exhaustion.  Perhaps you meant to specify" << endl
			     << "a small --bmaxdivn?" << endl;
		}
		if (!bmaxDivNSet) {
			bmaxDivN *= nthreads;
		}
		return abort;
	}
	
	public static void deleteIdxFiles(String outfile, boolean doRef, boolean justRef) {
		for(int i = 0; i < filesWritten.size(); i++) {
			System.err.println("Deleting \"" + filesWritten[i]
			     + "\" file written during aborted indexing attempt.");
			remove(filesWritten[i]);
		}
	}
	
	public static void driver(String infile, EList<String> infiles, String outfile, boolean packed, int reverse) {
		EList<FileBuf*> is(MISC_CAT);
		bool bisulfite = false;
		RefReadInParams refparams(false, reverse, nsToAs, bisulfite);
		assert_gt(infiles.size(), 0);
		if(format == CMDLINE) {
			// Adapt sequence strings to stringstreams open for input
			stringstream *ss = new stringstream();
			for(size_t i = 0; i < infiles.size(); i++) {
				(*ss) << ">" << i << endl << infiles[i].c_str() << endl;
			}
			FileBuf *fb = new FileBuf(ss);
			assert(fb != NULL);
			assert(!fb->eof());
			assert(fb->get() == '>');
			ASSERT_ONLY(fb->reset());
			assert(!fb->eof());
			is.push_back(fb);
		} else {
			// Adapt sequence files to ifstreams
			for(size_t i = 0; i < infiles.size(); i++) {
				FileBuf *fb;

				size_t idx = infiles[i].find_last_of(".");
				std::string ext = (idx == std::string::npos) ? "" : infiles[i].substr(idx + 1);
				if (ext == "" || ext == "gz" || ext == "Z") {
					gzFile zFp = gzopen(infiles[i].c_str(), "rb");
					if (zFp == NULL) {
						cerr << "Error: could not open "<< infiles[i].c_str() << endl;
						throw 1;
					}
					fb = new FileBuf(zFp);
				}
				else {
					FILE *f = fopen(infiles[i].c_str(), "rb");
					if (f == NULL) {
						cerr << "Error: could not open "<< infiles[i].c_str() << endl;
						throw 1;
					}
					fb = new FileBuf(f);
				}
				assert(fb != NULL);
				if(fb->peek() == -1 || fb->eof()) {
					cerr << "Warning: Empty fasta file: '" << infile.c_str() << "'" << endl;
					continue;
				}
				assert(!fb->eof());
				assert(fb->get() == '>');
				ASSERT_ONLY(fb->reset());
				assert(!fb->eof());
				is.push_back(fb);
			}
		}
		if(is.empty()) {
			cerr << "Warning: All fasta inputs were empty" << endl;
			throw 1;
		}
		if(!reverse) {
	#ifdef BOWTIE_64BIT_INDEX
	          if (verbose) cerr << "Building a LARGE index" << endl;
	#else
	          if (verbose) cerr << "Building a SMALL index" << endl;
	#endif
		}
		// Vector for the ordered list of "records" comprising the input
		// sequences.  A record represents a stretch of unambiguous
		// characters in one of the input sequences.
		EList<RefRecord> szs(MISC_CAT);
		std::pair<size_t, size_t> sztot;
		{
			if(verbose) cout << "Reading reference sizes" << endl;
			Timer _t(cout, "  Time reading reference sizes: ", verbose);
			if(!reverse && (writeRef || justRef)) {
				filesWritten.push_back(outfile + ".3." + gEbwt_ext);
				filesWritten.push_back(outfile + ".4." + gEbwt_ext);
				sztot = BitPairReference::szsFromFasta(is, outfile, bigEndian, refparams, szs, sanityCheck);
			} else {
				sztot = BitPairReference::szsFromFasta(is, string(), bigEndian, refparams, szs, sanityCheck);
			}
		}
		if(justRef) return;
		assert_gt(sztot.first, 0);
		assert_gt(sztot.second, 0);
		assert_gt(szs.size(), 0);
		// Construct index from input strings and parameters
		filesWritten.push_back(outfile + ".1." + gEbwt_ext);
		filesWritten.push_back(outfile + ".2." + gEbwt_ext);
		Ebwt ebwt(
			TStr(),
			packed,
			0,
			1,  // TODO: maybe not?
			lineRate,
			offRate,      // suffix-array sampling rate
			ftabChars,    // number of chars in initial arrow-pair calc
	              nthreads,
			outfile,      // basename for .?.ebwt files
			reverse == 0, // fw
			!entireSA,    // useBlockwise
			bmax,         // block size for blockwise SA builder
			bmaxMultSqrt, // block size as multiplier of sqrt(len)
			bmaxDivN,     // block size as divisor of len
			noDc? 0 : dcv,// difference-cover period
			is,           // list of input streams
			szs,          // list of reference sizes
			(TIndexOffU)sztot.first,  // total size of all unambiguous ref chars
			refparams,    // reference read-in parameters
			seed,         // pseudo-random number generator seed
			-1,           // override offRate
			doSaFile,     // make a file with just the suffix array in it
			doBwtFile,    // make a file with just the BWT string in it
			verbose,      // be talkative
			autoMem,      // pass exceptions up to the toplevel so that we can adjust memory settings automatically
			sanityCheck); // verify results and internal consistency
		// Note that the Ebwt is *not* resident in memory at this time.  To
		// load it into memory, call ebwt.loadIntoMemory()
		if(verbose) {
			// Print Ebwt's vital stats
			ebwt.eh().print(cout);
		}
		if(sanityCheck) {
			// Try restoring the original string (if there were
			// multiple texts, what we'll get back is the joined,
			// padded string, not a list)
			ebwt.loadIntoMemory(
				0,
				reverse ? (refparams.reverse == REF_READ_REVERSE) : 0,
				true,  // load SA sample?
				true,  // load ftab?
				true,  // load rstarts?
				false,
				false);
			SString<char> s2;
			ebwt.restore(s2);
			ebwt.evictFromMemory();
			{
				SString<char> joinedss = Ebwt::join<SString<char> >(
					is,          // list of input streams
					szs,         // list of reference sizes
					(TIndexOffU)sztot.first, // total size of all unambiguous ref chars
					refparams,   // reference read-in parameters
					seed);       // pseudo-random number generator seed
				if(refparams.reverse == REF_READ_REVERSE) {
					joinedss.reverse();
				}
				assert_eq(joinedss.length(), s2.length());
				assert(sstr_eq(joinedss, s2));
			}
			if(verbose) {
				if(s2.length() < 1000) {
					cout << "Passed restore check: " << s2.toZBuf() << endl;
				} else {
					cout << "Passed restore check: (" << s2.length() << " chars)" << endl;
				}
			}
		}
	}
}
