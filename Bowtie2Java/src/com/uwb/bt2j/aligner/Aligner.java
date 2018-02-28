package com.uwb.bt2j.aligner;

import org.omg.CORBA_2_3.portable.OutputStream;

import com.uwb.bt2j.util.EList;

import javafx.util.Pair;

class Aligner<T> {
	public static int FNAME_SIZE;
	public static int thread_counter;
	public static EList<String> mates1;
	public static EList<String> mates2;
	public static EList<String> mates12;
	public static EList<String> qualities;
	public static EList<String> qualities1;
	public static EList<String> qualities2;
	public static EList<String> queries;
	public static EList<String> presetList;
	public static EList<Pair<int, String>> extraOpts;
	
	public static String short_options = "fF:qbzhcu:rv:s:aP:t3:5:w:p:k:M:1:2:I:X:CQ:N:i:L:U:x:S:g:O:D:R:";
	public static String argstr;
	public static String arg0;
	public static String adjIdxBase;
	public static String origString;
	public static String metricsFile;
	public static String threadStealingDir;
	public static String rgid;
	public static String rgs;
	public static String rgs_optflag;
	public static String polstr;
	public static String defaultPreset;
	public static String wrapper;
	public static String outfile;
	public static String logDps;
	public static String logDpsOpp;
	public static String bt2index;
	
	public static double qUpTo;
	public static double khits;
	public static double mhits;
	public static double cacheLimit;
	public static double cacheSize;
	public static double skipReads;
	public static double fastaContLen;
	public static double fastaContFreq;
	public static double descLanding;
	public static double multiseedOff;
	public static double seedCacheLocalMB;
	public static double seedCacheCurrentMB;
	public static double exactCacheCurrentMB;
	public static double maxHalf;
	public static double maxIters;
	public static double maxUg;
	public static double maxDp;
	public static double maxItersIncr;
	public static double maxEeStreak;
	public static double maxUgStreak;
	public static double maxDpStreak;
	public static double maxStreakIncr;
	public static double maxMateStreak;
	public static double cminlen;
	public static double cpow2;
	public static double do1mmMinLen;
	public static double nSeedRounds;
	public static double extraOptsCur;
	
	public static float bwaSwLikeC;
	public static float bwaSwLikeT;
	
	public static Boolean startVerbose;
	public static Boolean showVersion;
	public static Boolean metricsStderr;
	public static Boolean metricsPerRead;
	public static Boolean allHits;
	public static Boolean solexaQuals;
	public static Boolean phred64Quals;
	public static Boolean integerQuals;
	public static Boolean threadStealing;
	public static Boolean noRefNames;
	public static Boolean fileParallel;
	public static Boolean useShmem;
	public static Boolean useMm;
	public static Boolean mmSweep;
	public static Boolean hadoopOut;
	public static Boolean fullRef;
	public static Boolean samTruncQname;
	public static Boolean samOmitSecSeqQual;
	public static Boolean samNoUnal;
	public static Boolean samNoHead;
	public static Boolean samNoSQ;
	public static Boolean sam_print_as;
	public static Boolean sam_print_xs;
	public static Boolean sam_print_xss;
	public static Boolean sam_print_yn;
	public static Boolean sam_print_xn;
	public static Boolean sam_print_x0;
	public static Boolean sam_print_x1;
	public static Boolean sam_print_xm;
	public static Boolean sam_print_xo;
	public static Boolean sam_print_xg;
	public static Boolean sam_print_nm;
	public static Boolean sam_print_md;
	public static Boolean sam_print_yf;
	public static Boolean sam_print_yi;
	public static Boolean sam_print_ym;
	public static Boolean sam_print_yp;
	public static Boolean sam_print_yt;
	public static Boolean sam_print_ys;
	public static Boolean sam_print_zs;
	public static Boolean sam_print_xr;
	public static Boolean sam_print_xt;
	public static Boolean sam_print_xd;
	public static Boolean sam_print_xu;
	public static Boolean sam_print_yl;
	public static Boolean sam_print_ye;
	public static Boolean sam_print_yu;
	public static Boolean sam_print_xp;
	public static Boolean sam_print_yr;
	public static Boolean sam_print_zb;
	public static Boolean sam_print_zr;
	public static Boolean sam_print_zf;
	public static Boolean sam_print_zm;
	public static Boolean sam_print_zi;
	public static Boolean sam_print_zp;
	public static Boolean sam_print_zu;
	public static Boolean sam_print_zt;
	public static Boolean bwaSwLike;
	public static Boolean gSeedLenIsSet;
	public static Boolean qcFilter;
	public static Boolean msample;
	public static Boolean msNoCache;
	public static Boolean penNCatPair;
	public static Boolean localAlign;
	public static Boolean noisyHpolymer;
	public static Boolean descPrioritizeRoots;
	public static Boolean seedSumm;
	public static Boolean scUnMapped;
	public static Boolean doUngapped;
	public static Boolean xeq;
	public static Boolean doExtend;
	public static Boolean enable8;
	public static Boolean doTri;
	public static Boolean ignoreQuals;
	public static Boolean doExactUpFront;
	public static Boolean do1mmUpFront;
	public static Boolean reorder;
	public static Boolean arbitraryRandom;
	public static Boolean bowtie2p5;
	public static Boolean saw_M;
	public static Boolean saw_a;
	public static Boolean saw_k;
	
	public Boolean gMate1fw;
	public Boolean gMate2fw;
	public Boolean gFlippedMatesOK;
	public Boolean gDovetailMatesOK;
	public Boolean gContainsMatesOK;
	public Boolean gOlapMatesOK;
	public Boolean gExpandToFrag;
	public Boolean gReportDiscordant;
	public Boolean gReportMixed;
	public Boolean gNoFw;
	public Boolean gNorc;
	public Boolean gReportOverhangs;
	
	public static void main(String[] args) {
		try {
			
			// Reset all global state, including getopt state
			int opterr, optind = 1;
			resetOptions();
			
			for(int i = 0; i < args.length; i++) {
				argstr += args[i];
				if(i < args.length-1) argstr += " ";
			}
			
			if(startVerbose) {
				System.err.println("Entered main(): ");
				//logTime(cerr, true);
			}
			
			parseOptions(args);
			arg0 = args[0];
			if(showVersion) {
				cout << argv0 << " version " << BOWTIE2_VERSION << endl;
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
			{
				Timer _t(cerr, "Overall time: ", timing);
				if(startVerbose) {
					cerr << "Parsing index and read arguments: "; logTime(cerr, true);
				}

				// Get index basename (but only if it wasn't specified via --index)
				if(bt2index.empty()) {
					cerr << "No index, query, or output file specified!" << endl;
					printUsage(cerr);
					return 1;
				}
		
				if(thread_stealing && thread_stealing_dir.empty()) {
					cerr << "When --thread-ceiling is specified, must also specify --thread-piddir" << endl;
					printUsage(cerr);
					return 1;
				}

				// Get query filename
				bool got_reads = !queries.empty() || !mates1.empty() || !mates12.empty();
				
				if(optind >= argc) {
					if(!got_reads) {
						printUsage(cerr);
						cerr << "***" << endl
						     << "Error: Must specify at least one read input with -U/-1/-2" << endl;
						return 1;
					}
				} else if(!got_reads) {
					// Tokenize the list of query files
					tokenize(argv[optind++], ",", queries);
					if(queries.empty()) {
						cerr << "Tokenized query file list was empty!" << endl;
						printUsage(cerr);
						return 1;
					}
				}

				// Get output filename
				if(optind < argc && outfile.empty()) {
					outfile = argv[optind++];
					cerr << "Warning: Output file '" << outfile.c_str()
					     << "' was specified without -S.  This will not work in "
						 << "future Bowtie 2 versions.  Please use -S instead."
						 << endl;
				}

				// Extra parametesr?
				if(optind < argc) {
					cerr << "Extra parameter(s) specified: ";
					for(int i = optind; i < argc; i++) {
						cerr << "\"" << argv[i] << "\"";
						if(i < argc-1) cerr << ", ";
					}
					cerr << endl;
					if(mates1.size() > 0) {
						cerr << "Note that if <mates> files are specified using -1/-2, a <singles> file cannot" << endl
							 << "also be specified.  Please run bowtie separately for mates and singles." << endl;
					}
					throw 1;
				}

				// Optionally summarize
				if(gVerbose) {
					cout << "Input " + gEbwt_ext +" file: \"" << bt2index.c_str() << "\"" << endl;
					cout << "Query inputs (DNA, " << file_format_names[format].c_str() << "):" << endl;
					for(size_t i = 0; i < queries.size(); i++) {
						cout << "  " << queries[i].c_str() << endl;
					}
					cout << "Quality inputs:" << endl;
					for(size_t i = 0; i < qualities.size(); i++) {
						cout << "  " << qualities[i].c_str() << endl;
					}
					cout << "Output file: \"" << outfile.c_str() << "\"" << endl;
					cout << "Local endianness: " << (currentlyBigEndian()? "big":"little") << endl;
					cout << "Sanity checking: " << (sanityCheck? "enabled":"disabled") << endl;
				#ifdef NDEBUG
					cout << "Assertions: disabled" << endl;
				#else
					cout << "Assertions: enabled" << endl;
				#endif
				}
				if(ipause) {
					cout << "Press key to continue..." << endl;
					getchar();
				}
				driver<SString<char> >("DNA", bt2index, outfile);
			}
			return 0;
		} catch(Exception e) {
			System.err.println("Error: Encountered exception: '" + e + "'");
			System.err.println("Command: ");
			for(int i = 0; i < argc; i++)
				System.err.println(argv[i] + " ");
			return 1;
		} catch(int e) {
			if(e != 0) {
				System.err.println("Error: Encountered internal Bowtie 2 exception (#" + e + ")");
				System.err.println("Command: ");
				for(int i = 0; i < argc; i++)
					System.err.println(argv[i] + " ");
			}
			return e;
	}
	}
	public static void resetOptions() {
		mates1.clear();
		mates2.clear();
		mates12.clear();
		adjIdxBase	            = "";
		gVerbose                = 0;
		startVerbose			= 0;
		gQuiet					= false;
		sanityCheck				= 0;  // enable expensive sanity checks
		format					= FASTQ; // default read format is FASTQ
		origString				= ""; // reference text, or filename(s)
		seed					= 0; // srandom() seed
		timing					= 0; // whether to report basic timing data
		metricsIval				= 1; // interval between alignment metrics messages (0 = no messages)
		metricsFile             = ""; // output file to put alignment metrics in
		metricsStderr           = false; // print metrics to stderr (in addition to --metrics-file if it's specified
		metricsPerRead          = false; // report a metrics tuple for every read?
		allHits					= false; // for multihits, report just one
		showVersion				= false; // just print version and quit?
		ipause					= 0; // pause before maching?
		qUpto					= 0xffffffff; // max # of queries to read
		gTrim5					= 0; // amount to trim from 5' end
		gTrim3					= 0; // amount to trim from 3' end
		offRate					= -1; // keep default offRate
		solexaQuals				= false; // quality strings are solexa quals, not phred, and subtract 64 (not 33)
		phred64Quals			= false; // quality chars are phred, but must subtract 64 (not 33)
		integerQuals			= false; // quality strings are space-separated strings of integers, not ASCII
		nthreads				= 1;     // number of pthreads operating concurrently
		thread_ceiling			= 0;     // max # threads user asked for
		thread_stealing_dir		= ""; // keep track of pids in this directory
		thread_stealing			= false; // true iff thread stealing is in use
		FNAME_SIZE				= 4096;
		outType					= OUTPUT_SAM;  // style of output
		noRefNames				= false; // true -> print reference indexes; not names
		khits					= 1;     // number of hits per read; >1 is much slower
		mhits					= 50;    // stop after finding this many alignments+1
		partitionSz				= 0;     // output a partitioning key in first field
		readsPerBatch			= 16;    // # reads to read from input file at once
		fileParallel			= false; // separate threads read separate input files in parallel
		useShmem				= false; // use shared memory to hold the index
		useMm					= false; // use memory-mapped files to hold the index
		mmSweep					= false; // sweep through memory-mapped files immediately after mapping
		gMinInsert				= 0;     // minimum insert size
		gMaxInsert				= 500;   // maximum insert size
		gMate1fw				= true;  // -1 mate aligns in fw orientation on fw strand
		gMate2fw				= false; // -2 mate aligns in rc orientation on fw strand
		gFlippedMatesOK         = false; // allow mates to be in wrong order
		gDovetailMatesOK        = false; // allow one mate to extend off the end of the other
		gContainMatesOK         = true;  // allow one mate to contain the other in PE alignment
		gOlapMatesOK            = true;  // allow mates to overlap in PE alignment
		gExpandToFrag           = true;  // incr max frag length to =larger mate len if necessary
		gReportDiscordant       = true;  // find and report discordant paired-end alignments
		gReportMixed            = true;  // find and report unpaired alignments for paired reads

		cacheLimit				= 5;     // ranges w/ size > limit will be cached
		cacheSize				= 0;     // # words per range cache
		skipReads				= 0;     // # reads/read pairs to skip
		gNofw					= false; // don't align fw orientation of read
		gNorc					= false; // don't align rc orientation of read
		fastaContLen			= 0;
		fastaContFreq			= 0;
		hadoopOut				= false; // print Hadoop status and summary messages
		fullRef					= false; // print entire reference name instead of just up to 1st space
		samTruncQname           = true;  // whether to truncate QNAME to 255 chars
		samOmitSecSeqQual       = false; // omit SEQ/QUAL for 2ndary alignments?
		samNoUnal               = false; // omit SAM records for unaligned reads
		samNoHead				= false; // don't print any header lines in SAM output
		samNoSQ					= false; // don't print @SQ header lines
		sam_print_as            = true;
		sam_print_xs            = true;
		sam_print_xss           = false; // Xs:i and Ys:i
		sam_print_yn            = false; // YN:i and Yn:i
		sam_print_xn            = true;
		sam_print_x0            = true;
		sam_print_x1            = true;
		sam_print_xm            = true;
		sam_print_xo            = true;
		sam_print_xg            = true;
		sam_print_nm            = true;
		sam_print_md            = true;
		sam_print_yf            = true;
		sam_print_yi            = false;
		sam_print_ym            = false;
		sam_print_yp            = false;
		sam_print_yt            = true;
		sam_print_ys            = true;
		sam_print_zs            = false;
		sam_print_xr            = false;
		sam_print_xt            = false;
		sam_print_xd            = false;
		sam_print_xu            = false;
		sam_print_yl            = false;
		sam_print_ye            = false;
		sam_print_yu            = false;
		sam_print_xp            = false;
		sam_print_yr            = false;
		sam_print_zb            = false;
		sam_print_zr            = false;
		sam_print_zf            = false;
		sam_print_zm            = false;
		sam_print_zi            = false;
		sam_print_zp            = false;
		sam_print_zu            = false;
		sam_print_zt            = false;
		bwaSwLike               = false;
		gSeedLenIsSet			= false;
		bwaSwLikeC              = 5.5f;
		bwaSwLikeT              = 20.0f;
		gDefaultSeedLen			= DEFAULT_SEEDLEN;
		qcFilter                = false; // don't believe upstream qc by default
		rgid					= "";    // SAM outputs for @RG header line
		rgs						= "";    // SAM outputs for @RG header line
		rgs_optflag				= "";    // SAM optional flag to add corresponding to @RG ID
		msample				    = true;
		gGapBarrier				= 4;     // disallow gaps within this many chars of either end of alignment
		qualities.clear();
		qualities1.clear();
		qualities2.clear();
		polstr.clear();
		msNoCache       = true; // true -> disable local cache
		bonusMatchType  = DEFAULT_MATCH_BONUS_TYPE;
		bonusMatch      = DEFAULT_MATCH_BONUS;
		penMmcType      = DEFAULT_MM_PENALTY_TYPE;
		penMmcMax       = DEFAULT_MM_PENALTY_MAX;
		penMmcMin       = DEFAULT_MM_PENALTY_MIN;
		penNType        = DEFAULT_N_PENALTY_TYPE;
		penN            = DEFAULT_N_PENALTY;
		penNCatPair     = DEFAULT_N_CAT_PAIR; // concatenate mates before N filtering?
		localAlign      = false;     // do local alignment in DP steps
		noisyHpolymer   = false;
		penRdGapConst   = DEFAULT_READ_GAP_CONST;
		penRfGapConst   = DEFAULT_REF_GAP_CONST;
		penRdGapLinear  = DEFAULT_READ_GAP_LINEAR;
		penRfGapLinear  = DEFAULT_REF_GAP_LINEAR;
		scoreMin.init  (SIMPLE_FUNC_LINEAR, DEFAULT_MIN_CONST,   DEFAULT_MIN_LINEAR);
		nCeil.init     (SIMPLE_FUNC_LINEAR, 0.0f, DMAX, 2.0f, 0.1f);
		msIval.init    (SIMPLE_FUNC_LINEAR, 1.0f, DMAX, DEFAULT_IVAL_B, DEFAULT_IVAL_A);
		descConsExp     = 2.0;
		descPrioritizeRoots = false;
		descLanding = 20;
		descentTotSz.init(SIMPLE_FUNC_LINEAR, 1024.0, DMAX, 0.0, 1024.0);
		descentTotFmops.init(SIMPLE_FUNC_LINEAR, 100.0, DMAX, 0.0, 10.0);
		multiseedMms    = DEFAULT_SEEDMMS;
		multiseedLen    = gDefaultSeedLen;
		multiseedOff    = 0;
		seedCacheLocalMB   = 32; // # MB to use for non-shared seed alignment cacheing
		seedCacheCurrentMB = 20; // # MB to use for current-read seed hit cacheing
		exactCacheCurrentMB = 20; // # MB to use for current-read seed hit cacheing
		maxhalf            = 15; // max width on one side of DP table
		seedSumm           = false; // print summary information about seed hits, not alignments
		scUnMapped         = false; // consider soft clipped bases unmapped when calculating TLEN
		xeq                = false; // use =/X instead of M in CIGAR string
		doUngapped         = true;  // do ungapped alignment
		maxIters           = 400;   // max iterations of extend loop
		maxUg              = 300;   // stop after this many ungap extends
		maxDp              = 300;   // stop after this many dp extends
		maxItersIncr       = 20;    // amt to add to maxIters for each -k > 1
		maxEeStreak        = 15;    // stop after this many end-to-end fails in a row
		maxUgStreak        = 15;    // stop after this many ungap fails in a row
		maxDpStreak        = 15;    // stop after this many dp fails in a row
		maxStreakIncr      = 10;    // amt to add to streak for each -k > 1
		maxMateStreak      = 10;    // in PE: abort seed range after N mate-find fails
		doExtend           = true;  // do seed extensions
		enable8            = true;  // use 8-bit SSE where possible?
		cminlen            = 2000;  // longer reads use checkpointing
		cpow2              = 4;     // checkpoint interval log2
		doTri              = false; // do triangular mini-fills?
		defaultPreset      = "sensitive%LOCAL%"; // default preset; applied immediately
		extra_opts.clear();
		extra_opts_cur = 0;
		bt2index.clear();        // read Bowtie 2 index from files with this prefix
		ignoreQuals = false;     // all mms incur same penalty, regardless of qual
		wrapper.clear();         // type of wrapper script, so we can print correct usage
		queries.clear();         // list of query files
		outfile.clear();         // write SAM output to this file
		mapqv = 2;               // MAPQ calculation version
		tighten = 3;             // -M tightening mode
		doExactUpFront = true;   // do exact search up front if seeds seem good enough
		do1mmUpFront = true;    // do 1mm search up front if seeds seem good enough
		seedBoostThresh = 300;   // if average non-zero position has more than this many elements
		nSeedRounds = 2;         // # rounds of seed searches to do for repetitive reads
		do1mmMinLen = 60;        // length below which we disable 1mm search
		reorder = false;         // reorder SAM records with -p > 1
		sampleFrac = 1.1f;       // align all reads
		arbitraryRandom = false; // let pseudo-random seeds be a function of read properties
		bowtie2p5 = false;
		logDps.clear();          // log seed-extend dynamic programming problems
		logDpsOpp.clear();       // log mate-search dynamic programming problems
	}
	
	public static void parseOptions(String[] args) {
		int option_index = 0;
		int next_option;
		saw_M = false;
		saw_a = false;
		saw_k = false;
		presetList.clear();
		if(startVerbose) { cerr << "Parsing options: "; logTime(cerr, true); }
		while(true) {
			next_option = getopt_long(
				argc, const_cast<char**>(argv),
				short_options, long_options, &option_index);
			const char * arg = optarg;
			if(next_option == EOF) {
				if(extra_opts_cur < extra_opts.size()) {
					next_option = extra_opts[extra_opts_cur].first;
					arg = extra_opts[extra_opts_cur].second.c_str();
					extra_opts_cur++;
				} else {
					break;
				}
			}
			parseOption(next_option, arg);
		}
		// Now parse all the presets.  Might want to pick which presets version to
		// use according to other parameters.
		auto_ptr<Presets> presets(new PresetsV0());
		// Apply default preset
		if(!defaultPreset.empty()) {
			polstr = applyPreset(defaultPreset, *presets.get()) + polstr;
		}
		// Apply specified presets
		for(size_t i = 0; i < presetList.size(); i++) {
			polstr += applyPreset(presetList[i], *presets.get());
		}
		for(size_t i = 0; i < extra_opts.size(); i++) {
			next_option = extra_opts[extra_opts_cur].first;
			const char *arg = extra_opts[extra_opts_cur].second.c_str();
			parseOption(next_option, arg);
		}
		// Remove initial semicolons
		while(!polstr.empty() && polstr[0] == ';') {
			polstr = polstr.substr(1);
		}
		if(gVerbose) {
			cerr << "Final policy string: '" << polstr.c_str() << "'" << endl;
		}
		size_t failStreakTmp = 0;
		SeedAlignmentPolicy::parseString(
			polstr,
			localAlign,
			noisyHpolymer,
			ignoreQuals,
			bonusMatchType,
			bonusMatch,
			penMmcType,
			penMmcMax,
			penMmcMin,
			penNType,
			penN,
			penRdGapConst,
			penRfGapConst,
			penRdGapLinear,
			penRfGapLinear,
			scoreMin,
			nCeil,
			penNCatPair,
			multiseedMms,
			multiseedLen,
			msIval,
			failStreakTmp,
			nSeedRounds);
		if(failStreakTmp > 0) {
			maxEeStreak = failStreakTmp;
			maxUgStreak = failStreakTmp;
			maxDpStreak = failStreakTmp;
		}
		if(saw_a || saw_k) {
			msample = false;
			mhits = 0;
		} else {
			assert_gt(mhits, 0);
			msample = true;
		}
		if(mates1.size() != mates2.size()) {
			cerr << "Error: " << mates1.size() << " mate files/sequences were specified with -1, but " << mates2.size() << endl
			     << "mate files/sequences were specified with -2.  The same number of mate files/" << endl
			     << "sequences must be specified with -1 and -2." << endl;
			throw 1;
		}
		if(qualities.size() && format != FASTA) {
			cerr << "Error: one or more quality files were specified with -Q but -f was not" << endl
			     << "enabled.  -Q works only in combination with -f and -C." << endl;
			throw 1;
		}
		if(qualities1.size() && format != FASTA) {
			cerr << "Error: one or more quality files were specified with --Q1 but -f was not" << endl
			     << "enabled.  --Q1 works only in combination with -f and -C." << endl;
			throw 1;
		}
		if(qualities2.size() && format != FASTA) {
			cerr << "Error: one or more quality files were specified with --Q2 but -f was not" << endl
			     << "enabled.  --Q2 works only in combination with -f and -C." << endl;
			throw 1;
		}
		if(qualities1.size() > 0 && mates1.size() != qualities1.size()) {
			cerr << "Error: " << mates1.size() << " mate files/sequences were specified with -1, but " << qualities1.size() << endl
			     << "quality files were specified with --Q1.  The same number of mate and quality" << endl
			     << "files must sequences must be specified with -1 and --Q1." << endl;
			throw 1;
		}
		if(qualities2.size() > 0 && mates2.size() != qualities2.size()) {
			cerr << "Error: " << mates2.size() << " mate files/sequences were specified with -2, but " << qualities2.size() << endl
			     << "quality files were specified with --Q2.  The same number of mate and quality" << endl
			     << "files must sequences must be specified with -2 and --Q2." << endl;
			throw 1;
		}
		if(!rgs.empty() && rgid.empty()) {
			cerr << "Warning: --rg was specified without --rg-id also "
			     << "being specified.  @RG line is not printed unless --rg-id "
				 << "is specified." << endl;
		}
		// Check for duplicate mate input files
		if(format != CMDLINE) {
			for(size_t i = 0; i < mates1.size(); i++) {
				for(size_t j = 0; j < mates2.size(); j++) {
					if(mates1[i] == mates2[j] && !gQuiet) {
						cerr << "Warning: Same mate file \"" << mates1[i].c_str() << "\" appears as argument to both -1 and -2" << endl;
					}
				}
			}
		}
		// If both -s and -u are used, we need to adjust qUpto accordingly
		// since it uses rdid to know if we've reached the -u limit (and
		// rdids are all shifted up by skipReads characters)
		if(qUpto + skipReads > qUpto) {
			qUpto += skipReads;
		}
		if(useShmem && useMm && !gQuiet) {
			cerr << "Warning: --shmem overrides --mm..." << endl;
			useMm = false;
		}
		if(gGapBarrier < 1) {
			cerr << "Warning: --gbar was set less than 1 (=" << gGapBarrier
			     << "); setting to 1 instead" << endl;
			gGapBarrier = 1;
		}
		if(bonusMatch > 0 && !scoreMin.alwaysPositive()) {
			cerr << "Error: the match penalty is greater than 0 (" << bonusMatch
			     << ") but the --score-min function can be less than or equal to "
				 << "zero.  Either let the match penalty be 0 or make --score-min "
				 << "always positive." << endl;
			throw 1;
		}
		if(multiseedMms >= multiseedLen) {
			assert_gt(multiseedLen, 0);
			cerr << "Warning: seed mismatches (" << multiseedMms
			     << ") is less than seed length (" << multiseedLen
				 << "); setting mismatches to " << (multiseedMms-1)
				 << " instead" << endl;
			multiseedMms = multiseedLen-1;
		}
		sam_print_zm = sam_print_zm && bowtie2p5;
	#ifndef NDEBUG
		if(!gQuiet) {
			cerr << "Warning: Running in debug mode.  Please use debug mode only "
				 << "for diagnosing errors, and not for typical use of Bowtie 2."
				 << endl;
		}
	}
	
	public void printArgDesc(OutputStream out) {
		double i = 0;
		while(long_options[i].name != 0) {
			out.write(long_options[i].name + "\t"
			    + (long_options[i].has_arg == no_argument ? 0 : 1));
			i++;
		}
		double solen = strlen(short_options);
		for(i = 0; i < solen; i++) {
			// Has an option?  Does if next char is :
			if(i == solen-1) {
				System.out.println((char)short_options[i] + "\t" + 0);
			} else {
				if(short_options[i+1] == ':') {
					// Option with argument
					System.out.println((char)short_options[i] + "\t" + 1);
					i++; // skip the ':'
				} else {
					// Option with no argument
					System.out.println((char)short_options[i] + "\t" + 0 );
				}
			}
		}
	}
	
	public static void printUsage(OutputStream out) {
		out.write("Bowtie 2 version " + BOWTIE2_VERSION + " by Ben Langmead (langmea@cs.jhu.edu, www.cs.jhu.edu/~langmea)");
		String tool_name = "bowtie2-align";
		if(wrapper == "basic-0") {
			tool_name = "bowtie2";
		}
		out.write("Usage: " + "\n"
		    + "  " + tool_name + " [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i>} [-S <sam>]" + "\n"
		    + "\n"
			+     "  <bt2-idx>  Index filename prefix (minus trailing .X." + gEbwt_ext + ")." + "\n"
			+     "             NOTE: Bowtie 1 and Bowtie 2 indexes are not compatible." + "\n"
		    +     "  <m1>       Files with #1 mates, paired with files in <m2>." + "\n");
		if(wrapper == "basic-0") {
			out.write("             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)." + "\n");
		}
		out.write("  <m2>       Files with #2 mates, paired with files in <m1>." + "\n");
		if(wrapper == "basic-0") {
			out.write("             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)." + "\n");
		}
		out.write("  <r>        Files with unpaired reads." + "\n");
		if(wrapper == "basic-0") {
			out.write("             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)." + "\n");
		}
		out.write("  <i>        Files with interleaved paired-end FASTQ reads" + "\n");
		if(wrapper == "basic-0") {
			out.write("             Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2)." + "\n");
		}
		out.write("  <sam>      File for SAM output (default: stdout)" + "\n"
		    + "\n"
		    + "  <m1>, <m2>, <r> can be comma-separated lists (no whitespace) and can be" + "\n"
			+ "  specified many times.  E.g. '-U file1.fq,file2.fq -U file3.fq'." + "\n"
			// Wrapper script should write <bam> line next
			+ "\n"
		    + "Options (defaults in parentheses):" + "\n"
			+ "\n"
		    + " Input:" + "\n"
		    + "  -q                 query input files are FASTQ .fq/.fastq (default)" + "\n"
			+ "  --tab5             query input files are TAB5 .tab5" + "\n"
			+ "  --tab6             query input files are TAB6 .tab6" + "\n"
		    + "  --qseq             query input files are in Illumina's qseq format" + "\n"
		    + "  -f                 query input files are (multi-)FASTA .fa/.mfa" + "\n"
		    + "  -r                 query input files are raw one-sequence-per-line"+ "\n"
		    + "  -F k:<int>,i:<int> query input files are continuous FASTA where reads" + "\n"
		    + "                     are substrings (k-mers) extracted from a FASTA file <s>" + "\n"
		    + "                     and aligned at offsets 1, 1+i, 1+2i ... end of reference" + "\n"
		    + "  -c                 <m1>, <m2>, <r> are sequences themselves, not files" + "\n"
		    + "  -s/--skip <int>    skip the first <int> reads/pairs in the input (none)" + "\n"
		    + "  -u/--upto <int>    stop after first <int> reads/pairs (no limit)" + "\n"
		    + "  -5/--trim5 <int>   trim <int> bases from 5'/left end of reads (0)" + "\n"
		    + "  -3/--trim3 <int>   trim <int> bases from 3'/right end of reads (0)" + "\n"
		    + "  --phred33          qualities are Phred+33 (default)" + "\n"
		    + "  --phred64          qualities are Phred+64" + "\n"
		    + "  --int-quals        qualities encoded as space-delimited integers"+ "\n"
		    + "\n"
		    + " Presets:                 Same as:" + "\n"
			+ "  For --end-to-end:" + "\n"
			+ "   --very-fast            -D 5 -R 1 -N 0 -L 22 -i S,0,2.50" + "\n"
			+ "   --fast                 -D 10 -R 2 -N 0 -L 22 -i S,0,2.50" + "\n"
			+ "   --sensitive            -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 (default)" + "\n"
			+ "   --very-sensitive       -D 20 -R 3 -N 0 -L 20 -i S,1,0.50" + "\n"
			+ "\n"
			+ "  For --local:" + "\n"
			+ "   --very-fast-local      -D 5 -R 1 -N 0 -L 25 -i S,1,2.00" + "\n"
			+ "   --fast-local           -D 10 -R 2 -N 0 -L 22 -i S,1,1.75" + "\n"
			+ "   --sensitive-local      -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 (default)" + "\n"
			+ "   --very-sensitive-local -D 20 -R 3 -N 0 -L 20 -i S,1,0.50" + "\n"
			+ "\n"
		    + " Alignment:" + "\n"
			+ "  -N <int>           max # mismatches in seed alignment; can be 0 or 1 (0)"+ "\n"
			+ "  -L <int>           length of seed substrings; must be >3, <32 (22)" + "\n"
			+ "  -i <func>          interval between seed substrings w/r/t read len (S,1,1.15)" + "\n"
			+ "  --n-ceil <func>    func for max # non-A/C/G/Ts permitted in aln (L,0,0.15)" + "\n"
			+ "  --dpad <int>       include <int> extra ref chars on sides of DP table (15)" + "\n"
			+ "  --gbar <int>       disallow gaps within <int> nucs of read extremes (4)" + "\n"
			+ "  --ignore-quals     treat all quality values as 30 on Phred scale (off)" + "\n"
		    + "  --nofw             do not align forward (original) version of read (off)" + "\n"
		    + "  --norc             do not align reverse-complement version of read (off)" + "\n"
		    + "  --no-1mm-upfront   do not allow 1 mismatch alignments before attempting to" + "\n"
		    + "                     scan for the optimal seeded alignments"
		    + "\n"
			+ "  --end-to-end       entire read must align; no clipping (on)" + "\n"
			+ "   OR" + "\n"
			+ "  --local            local alignment; ends might be soft clipped (off)" + "\n"
			+ "\n"
		    + " Scoring:" + "\n"
			+ "  --ma <int>         match bonus (0 for --end-to-end, 2 for --local) " + "\n"
			+ "  --mp <int>         max penalty for mismatch; lower qual = lower penalty (6)" + "\n"
			+ "  --np <int>         penalty for non-A/C/G/Ts in read/ref (1)" + "\n"
			+ "  --rdg <int>,<int>  read gap open, extend penalties (5,3)" + "\n"
			+ "  --rfg <int>,<int>  reference gap open, extend penalties (5,3)" + "\n"
			+ "  --score-min <func> min acceptable alignment score w/r/t read length" + "\n"
			+ "                     (G,20,8 for local, L,-0.6,-0.6 for end-to-end)" + "\n"
			+ "\n"
		    + " Reporting:" + "\n"
		    + "  (default)          look for multiple alignments, report best, with MAPQ" + "\n"
			+ "   OR" + "\n"
		    + "  -k <int>           report up to <int> alns per read; MAPQ not meaningful" + "\n"
			+ "   OR" + "\n"
		    + "  -a/--all           report all alignments; very slow, MAPQ not meaningful" + "\n"
		    + "\n"
		    + " Effort:" + "\n"
		    + "  -D <int>           give up extending after <int> failed extends in a row (15)" + "\n"
		    + "  -R <int>           for reads w/ repetitive seeds, try <int> sets of seeds (2)" + "\n"
		    + "\n"
			+ " Paired-end:" + "\n"
		    + "  -I/--minins <int>  minimum fragment length (0)" + "\n"
		    + "  -X/--maxins <int>  maximum fragment length (500)" + "\n"
		    + "  --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (--fr)"+ "\n"
			+ "  --no-mixed         suppress unpaired alignments for paired reads" + "\n"
			+ "  --no-discordant    suppress discordant alignments for paired reads" + "\n"
			+ "  --dovetail         concordant when mates extend past each other" + "\n"
			+ "  --no-contain       not concordant when one mate alignment contains other" + "\n"
			+ "  --no-overlap       not concordant when mates overlap at all" + "\n"
			+ "\n"
		    + " Output:");

		out.write("  -t/--time          print wall-clock time taken by search phases" + "\n");
		if(wrapper == "basic-0") {
		out.write("  --un <path>        write unpaired reads that didn't align to <path>" + "\n"
		    + "  --al <path>        write unpaired reads that aligned at least once to <path>" + "\n"
		    + "  --un-conc <path>   write pairs that didn't align concordantly to <path>" + "\n"
		    + "  --al-conc <path>   write pairs that aligned concordantly at least once to <path>" + "\n"
		    + "    (Note: for --un, --al, --un-conc, or --al-conc, add '-gz' to the option name, e.g." + "\n"
		    + "    --un-gz <path>, to gzip compress output, or add '-bz2' to bzip2 compress output.)" + "\n");
		}
		out.write("  --quiet            print nothing to stderr except serious errors" + "\n"
			+ "  --met-file <path>  send metrics to file at <path> (off)" + "\n"
			+ "  --met-stderr       send metrics to stderr (off)" + "\n"
			+ "  --met <int>        report internal counters & metrics every <int> secs (1)" + "\n"
		// Following is supported in the wrapper instead
		    + "  --no-unal          suppress SAM records for unaligned reads" + "\n"
		    + "  --no-head          suppress header lines, i.e. lines starting with @" + "\n"
		    + "  --no-sq            suppress @SQ header lines" + "\n"
		    + "  --rg-id <text>     set read group id, reflected in @RG line and RG:Z: opt field" + "\n"
		    + "  --rg <text>        add <text> (\"lab:value\") to @RG line of SAM header." + "\n"
		    + "                     Note: @RG line only printed when --rg-id is set." + "\n"
		    + "  --omit-sec-seq     put '*' in SEQ and QUAL fields for secondary alignments." + "\n"
		    + "  --sam-no-qname-trunc Suppress standard behavior of truncating readname at first whitespace " + "\n"
		    + "                      at the expense of generating non-standard SAM." + "\n"
		    + "  --xeq              Use '='/'X', instead of 'M,' to specify matches/mismatches in SAM record." + "\n"
		    + "  --soft-clipped-unmapped-tlen Exclude soft-clipped bases when reporting TLEN" + "\n"
		    + "\n"
		    + " Performance:" + "\n"
		    + "  -p/--threads <int> number of alignment threads to launch (1)" + "\n"
		    + "  --reorder          force SAM output order to match order of input reads"+ "\n"
		    + "\n"
		    + " Other:" + "\n"
			+ "  --qc-filter        filter out reads that are bad according to QSEQ filter" + "\n"
		    + "  --seed <int>       seed for random number generator (0)" + "\n"
		    + "  --non-deterministic seed rand. gen. arbitrarily instead of using read attributes" + "\n"
		    + "  --version          print version information and quit" + "\n"
		    + "  -h/--help          print this usage message" + "\n");
		if(wrapper.empty()) {
			System.err.println("*** Warning ***" + "\n" + "'bowtie2-align' was run directly.  It is recommended that you run the wrapper script 'bowtie2' instead.");
		}
	}
	
	public static int parseInt(int lower, int upper, String errmsg, String arg) {
		long l = Long.parseLong(arg);
		
			if (l < lower || l > upper) {
				System.err.println(errmsg);
				//printUsage();
			}
			return (int)l;
		System.err.println(errmsg);
		//printUsage();
		return -1;
	}
	
	public static int parseInt(int lower, String errmsg, String arg) {
		return parseInt(lower, Integer.MAX_VALUE, errmsg, arg);
	}
	
	public Pair<T, T> parsePair(String s, char delim) {
		EList<String> ss = tokenize(s, delim);
		Pair<T, T> ret;
		ret.first = parse<T>(ss[0]);
		ret.second = parse<T>(ss[1]);
		return ret;
	}
	
	public void parseTuple(String s, char delim, EList<T> ret) {
		EList<String> ss = tokenize(s, delim);
		for(double i = 0; i < ss.size(); i++) {
			ret.push_back(parse<T>(ss[i]));
		}
	}
	
	public static String applyPreset(String sorig, Presets presets) {
		String s = sorig;
		double found = s.find("%LOCAL%");
		if(found != string::npos) {
			s.replace(found, strlen("%LOCAL%"), localAlign ? "-local" : "");
		}
		if(gVerbose) {
			System.err.println("Applying preset: '" + s + "' using preset menu '"
				 + presets.name() + "'" + "\n");
		}
		String pol;
		presets.apply(s, pol, extra_opts);
		return pol;
	}
}