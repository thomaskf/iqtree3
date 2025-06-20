#include "libiqtree_fun.h"

#if defined WIN32 || defined _WIN32 || defined __WIN32__ || defined WIN64
#include <winsock2.h>
#endif

class input_options {
public:
    vector<string> flags;
    vector<string> values;

    void insert(string flag, string value="") {
        flags.push_back(flag);
        values.push_back(value);
    }
    
    // set Params according to the input options from PiQTREE
    // only invoke this function after the default values of parameters are set
    // this function defines which IQ-TREE options are availble for PiQTREE
    void set_params(Params& params);
};

void cleanup(Params& params) {
    if (params.state_freq_set != NULL) {
        delete[] params.state_freq_set;
        params.state_freq_set = NULL;
    }
}

void convertToVectorStr(StringArray& names, StringArray& seqs, vector<string>& names_vec, vector<string>& seqs_vec) {
	names_vec.clear();
	seqs_vec.clear();
	if (names.length != seqs.length)
	    outError("The number of names must equal to the number of sequences");
	for (int i = 0; i < names.length; i++) {
		names_vec.push_back(string(names.strings[i]));
		seqs_vec.push_back(string(seqs.strings[i]));
	}
}

char* build_phylogenetic(StringArray& cnames, StringArray& cseqs, const char* cmodel, const char* cintree,
                          int rand_seed, string& prog, input_options* in_options);

// Calculates the robinson fould distance between two trees
extern "C" IntegerResult robinson_fould(const char* ctree1, const char* ctree2) {
    IntegerResult output;
    output.errorStr = strdup("");

    try {
        string tree1 = string(ctree1);
        string tree2 = string(ctree2);

        MTree first_tree;
        bool is_rooted = false;
        std::vector<double> rfdist;

        // read in the first tree
        first_tree.read_TreeString(tree1, is_rooted);
        
        // second tree
        stringstream second_tree_str;
        second_tree_str << tree2;
        second_tree_str.seekg(0, ios::beg);
        
        // compute the RF distance
        first_tree.computeRFDist(second_tree_str, rfdist);
        
        output.value = (int)rfdist[0];
    } catch (const exception& e) {
    	if (strlen(e.what()) > 0) {
    	    output.errorStr = new char[strlen(e.what())+1];
    	    strcpy(output.errorStr, e.what());
    	}
        // reset the output and error buffers
        funcExit();
    }

    return output;
}

// Generates a set of random phylogenetic trees
// tree_gen_mode allows:"YULE_HARDING", "UNIFORM", "CATERPILLAR", "BALANCED", "BIRTH_DEATH", "STAR_TREE"
extern "C" StringResult random_tree(int num_taxa, const char* tree_gen_mode, int num_trees, int rand_seed) {
    StringResult output;
    output.errorStr = strdup("");
    
    try {
    	ostringstream ostring;
        PhyloTree ptree;
        int seed = rand_seed;
        if (seed == 0)
            seed = make_new_seed();
        cout << "seed: " << seed << endl;
        init_random(seed);
        
        TreeGenType tree_mode;
        if (strcmp(tree_gen_mode, "YULE_HARDING")==0) {
            tree_mode = YULE_HARDING;
        } else if (strcmp(tree_gen_mode, "UNIFORM")==0) {
            tree_mode = UNIFORM;
        } else if (strcmp(tree_gen_mode, "CATERPILLAR")==0) {
            tree_mode = CATERPILLAR;
        } else if (strcmp(tree_gen_mode, "BALANCED")==0) {
            tree_mode = BALANCED;
        } else if (strcmp(tree_gen_mode, "BIRTH_DEATH")==0) {
            tree_mode = BIRTH_DEATH;
        } else if (strcmp(tree_gen_mode, "STAR_TREE")==0) {
            tree_mode = STAR_TREE;
        } else {
            outError("Unknown mode: " + string(tree_gen_mode));
        }
        
        Params params = Params::getInstance();
        params.setDefault();
        params.sub_size = num_taxa;
        params.tree_gen = tree_mode;
        params.repeated_time = num_trees;
        params.ignore_checkpoint = true; // overrid the output file if exists
        params.user_file = (char*) "";
        
        generateRandomTree(params, ostring);
        string s = ostring.str();
        
        if (s.length() > 0) {
            output.value = new char[s.length()+1];
    	    strcpy(output.value, s.c_str());
        }
        
    } catch (const exception& e) {
    	if (strlen(e.what()) > 0) {
    	    output.errorStr = new char[strlen(e.what())+1];
    	    strcpy(output.errorStr, e.what());
    	}
        // reset the output and error buffers
        funcExit();
    }
    return output;
}

// Perform phylogenetic analysis on the input alignment (in string format)
// With estimation of the best topology
// num_thres -- number of cpu threads to be used, default: 1; 0 - auto detection of the optimal number of cpu threads
extern "C" StringResult build_tree(StringArray& names, StringArray& seqs, const char* model, int rand_seed, int bootstrap_rep, int num_thres) {
    StringResult output;
    output.errorStr = strdup("");
    
    try {
    	const char* intree = "";
        input_options* in_options = NULL;
        if (bootstrap_rep > 0 || num_thres >= 0) {
            in_options = new input_options();
            if (bootstrap_rep > 0)
                in_options->insert("-bb", convertIntToString(bootstrap_rep));
            if (num_thres >= 0)
                in_options->insert("-nt", convertIntToString(num_thres));
        }
        string prog = "build_tree";
        output.value = build_phylogenetic(names, seqs, model, intree, rand_seed, prog, in_options);
        if (in_options != NULL)
            delete in_options;
    } catch (const exception& e) {
    	if (strlen(e.what()) > 0) {
    	    output.errorStr = new char[strlen(e.what())+1];
    	    strcpy(output.errorStr, e.what());
    	}
        // reset the output and error buffers
        funcExit();
    }
    return output;
}

// Perform phylogenetic analysis on the input alignment (in string format)
// With restriction to the input toplogy
// num_thres -- number of cpu threads to be used, default: 1; 0 - auto detection of the optimal number of cpu threads
extern "C" StringResult fit_tree(StringArray& names, StringArray& seqs, const char* model, const char* intree, int rand_seed, int num_thres) {
    StringResult output;
    output.errorStr = strdup("");
    
    try {
        input_options* in_options = NULL;
        if (num_thres >= 0) {
            in_options = new input_options();
            in_options->insert("-nt", convertIntToString(num_thres));
        }
        string prog = "fit_tree";
        output.value = build_phylogenetic(names, seqs, model, intree, rand_seed, prog, in_options);
        if (in_options != NULL)
            delete in_options;
    } catch (const exception& e) {
    	if (strlen(e.what()) > 0) {
    	    output.errorStr = new char[strlen(e.what())+1];
    	    strcpy(output.errorStr, e.what());
    	}
        // reset the output and error buffers
        funcExit();
    }
    return output;
}

// Perform phylogenetic analysis with ModelFinder
// on the input alignment (in string format)
// model_set -- a set of models to consider
// freq_set -- a set of frequency types
// rate_set -- a set of RHAS models
// num_thres -- number of cpu threads to be used, default: 1; 0 - auto detection of the optimal number of cpu threads
extern "C" StringResult modelfinder(StringArray& names, StringArray& seqs, int rand_seed, const char* model_set, const char* freq_set, const char* rate_set, int num_thres) {
	StringResult output;
    output.errorStr = strdup("");
	
    try {
		input_options* in_options = NULL;
		const char* intree = "";
		const char* model = "MF"; // modelfinder
		int i;
        in_options = new input_options();
        // handle model_set, freq_set, rate_set
        if (strlen(model_set) > 0)
            in_options->insert("-mset", string(model_set));
        if (strlen(freq_set) > 0)
            in_options->insert("-mfreq", string(freq_set));
        if (strlen(rate_set) > 0)
            in_options->insert("-mrate", string(rate_set));
        if (num_thres >= 0)
            in_options->insert("-nt", convertIntToString(num_thres));
		string prog = "modelfinder";
        output.value = build_phylogenetic(names, seqs, model, intree, rand_seed, prog, in_options);
        
        delete in_options;
    } catch (const exception& e) {
    	if (strlen(e.what()) > 0) {
    	    output.errorStr = new char[strlen(e.what())+1];
    	    strcpy(output.errorStr, e.what());
    	}
        // reset the output and error buffers
        funcExit();
    }
    return output;
}

// Build pairwise JC distance matrix
// output: set of distances
// (n * i + j)-th element of the list represents the distance between i-th and j-th sequence,
// where n is the number of sequences
// num_thres -- number of cpu threads to be used, default: 1; 0 - use all available cpu threads on the machine
extern "C" DoubleArrayResult build_distmatrix(StringArray& cnames, StringArray& cseqs, int num_thres) {
    DoubleArrayResult output;
    output.errorStr = strdup("");
    
    try {
		string prog = "build_matrix";
		output.length = 0;
		output.value = NULL;
    	vector<string> names;
    	vector<string> seqs;
    	convertToVectorStr(cnames, cseqs, names, seqs);
        int n = names.size();
        int n_sq = n * n;
        if (n_sq >= 1) {
            output.length = n_sq;
            output.value = new double[n_sq];
        }
        if (n == 1) {
            output.value[0] = 0.0;
        } else if (n > 1) {
            // extern VerboseMode verbose_mode;
            progress_display::setProgressDisplay(false);
            verbose_mode = VB_QUIET; // (quiet mode)
            Params params = Params::getInstance();
            params.setDefault();

            int rand_seed = make_new_seed();
            string out_prefix_str = prog + "_" + convertIntToString(rand_seed);
            _log_file = out_prefix_str + ".log";
            bool append_log = false;
            startLogFile(append_log);

        #ifdef _OPENMP
            int max_procs = countPhysicalCPUCores();
            if (num_thres > max_procs)
                num_thres = max_procs;
            if (num_thres > 0) {
                Params::getInstance().num_threads = num_thres;
                omp_set_num_threads(num_thres);
            } else if (num_thres == 0) {
                Params::getInstance().num_threads = max_procs;
                omp_set_num_threads(max_procs);
            }
         #endif
            
            PhyloTree ptree;
            ptree.aln = new Alignment(names, seqs, params.sequence_type, params.model_name);
            
            
            // compute the matrix
            #ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic, 1)
            #endif
            for (int i = 0; i < n; i++) {
                double* dmat = &output.value[i * n];
                dmat[i] = 0.0;
                for (int j = i+1; j < n; j++) {
                    dmat[j] = ptree.aln->computeJCDist(i, j);
                    output.value[j * n + i] = dmat[j];
                }
            }
            
            delete ptree.aln;
            funcExit();
        }
    } catch (const exception& e) {
    	if (strlen(e.what()) > 0) {
    	    output.errorStr = new char[strlen(e.what())+1];
    	    strcpy(output.errorStr, e.what());
    	}
        // reset the output and error buffers
        funcExit();
    }
    return output;
}

// Using Rapid-NJ to build tree from a distance matrix
extern "C" StringResult build_njtree(StringArray& cnames, DoubleArray& distances) {
    StringResult output;
    output.errorStr = strdup("");

    try {
    	// convert to vector
    	vector<string> names;
		for (int i = 0; i < cnames.length; i++) {
			names.push_back(string(cnames.strings[i]));
		}
        // check the size of names and distances
        if (names.size() < 3)
            outError("The size of names must be at least 3");
        size_t n = names.size();
        size_t sq_n = n * n;
        if (distances.length != sq_n)
            outError("The size of distances must equal to the square of the size of names");
        
        string prog = "build_njtree";
        // extern VerboseMode verbose_mode;
        progress_display::setProgressDisplay(false);
        verbose_mode = VB_QUIET; // (quiet mode)
        Params params = Params::getInstance();
        params.setDefault();
        
        int rand_seed = make_new_seed();
        string out_prefix_str = prog + "_" + convertIntToString(rand_seed);
        _log_file = out_prefix_str + ".log";
        bool append_log = false;
        startLogFile(append_log);
        
        string algn_name = "NJ-R"; // Rapid NJ
        StartTree::BuilderInterface* algorithm = StartTree::Factory::getTreeBuilderByName(algn_name);
        stringstream stree;
        if (!algorithm->constructTreeInMemory2(names, distances.doubles, stree)) {
            outError("Tree construction failed.");
        }
        string s = stree.str();
        if (s.length() > 0) {
            output.value = new char[s.length()+1];
    	    strcpy(output.value, s.c_str());
        }
        funcExit();
    } catch (const exception& e) {
    	if (strlen(e.what()) > 0) {
    	    output.errorStr = new char[strlen(e.what())+1];
    	    strcpy(output.errorStr, e.what());
    	}
        // reset the output and error buffers
        funcExit();
    }
    return output;
}

/*
 * Compute a consensus tree
 * trees -- a set of input trees
 * minsup -- a threshold to build the majority consensus, default is 0.0
 * output: the consensus tree of the set of input trees
 */
extern "C" StringResult consensus_tree(StringArray& trees, double minsup) {

    StringResult output;
    output.errorStr = strdup("");

    try {

       	// convert to vector
       	vector<string> trees_vec;
		for (int i = 0; i < trees.length; i++) {
			trees_vec.push_back(string(trees.strings[i]));
		}
		if (trees.length < 1) {
		    outError("Error! The number of input trees < 1");
		}

        Params params = Params::getInstance();
        params.setDefault();
        params.split_threshold = minsup;
        bool rooted = false;

        // // get a set of taxa names
        vector<string> taxa_names;
        MTree mtree;
        mtree.read_TreeString(trees_vec[0], rooted);
        mtree.getTaxaName(taxa_names);

        // read the trees
        MTreeSet boot_trees;
        boot_trees.init(trees_vec, taxa_names, rooted);

        SplitGraph sg;
        SplitIntMap hash_ss;
        double scale = 100.0;
        if (params.scaling_factor > 0)
            scale = params.scaling_factor;
        boot_trees.convertSplits(sg, params.split_threshold, SW_COUNT, params.split_weight_threshold);
        scale /= boot_trees.sumTreeWeights();
        cout << sg.size() << " splits found" << endl;

        if (params.scaling_factor < 0)
            sg.scaleWeight(scale, true);
        else {
            sg.scaleWeight(scale, false, params.numeric_precision);
        }

        MTree mytree;
        SplitGraph maxsg;
        sg.findMaxCompatibleSplits(maxsg);
        mytree.convertToTree(maxsg);
        if (!mytree.rooted) {
            string taxname;
            if (params.root)
                taxname = params.root;
            else
                taxname = sg.getTaxa()->GetTaxonLabel(0);
            Node *node = mytree.findLeafName(taxname);
            if (node)
                mytree.root = node;
        }

       	ostringstream ostring;
        mytree.printTree(ostring, WT_BR_CLADE);
        string s = ostring.str();
        if (s.length() > 0) {
            output.value = new char[s.length()+1];
    	    strcpy(output.value, s.c_str());
        }

    } catch (const exception& e) {
    	if (strlen(e.what()) > 0) {
    	    output.errorStr = new char[strlen(e.what())+1];
    	    strcpy(output.errorStr, e.what());
    	}
        // reset the output and error buffers
        funcExit();
    }
    return output;
}

// verion number
extern "C" StringResult version() {
    StringResult output;
    output.errorStr = strdup("");
    try {
		stringstream ss;
		ss << iqtree_VERSION_MAJOR << "." << iqtree_VERSION_MINOR << iqtree_VERSION_PATCH;
		string s = ss.str();
		if (s.length() > 0) {
			output.value = new char[s.length()+1];
			strcpy(output.value, s.c_str());
		}
	} catch (const exception& e) {
    	if (strlen(e.what()) > 0) {
    	    output.errorStr = new char[strlen(e.what())+1];
    	    strcpy(output.errorStr, e.what());
    	}
    }
    return output;
}

// ----------------------------------------------
// function for performing plylogenetic analysis
// ----------------------------------------------

// Perform phylogenetic analysis on the input alignment (in string format)
// if intree exists, then the topology will be restricted to the intree
char* build_phylogenetic(StringArray& cnames, StringArray& cseqs, const char* cmodel, const char* cintree,
                          int rand_seed, string& prog, input_options* in_options) {
    // perform phylogenetic analysis on the input sequences
    // all sequences have to be the same length

    int instruction_set;
    
    vector<string> names, seqs;
    
    convertToVectorStr(cnames, cseqs, names, seqs);
    string model = string(cmodel);
    string intree = string(cintree);
    
    // checking whether all seqs are in the same length
    if (seqs.size() > 0) {
        int slen = seqs[0].length();
        for (int i=1; i<seqs.size(); i++) {
            if (seqs[i].length() != slen) {
                outError("The input sequences are not in the same length");
            }
        }
    }

    extern VerboseMode verbose_mode;
    progress_display::setProgressDisplay(false);
    // verbose_mode = VB_MIN;
    verbose_mode = VB_QUIET; // (quiet mode)
    Params::getInstance().setDefault();
    Params::getInstance().num_threads = 1; // default
    Params::getInstance().aln_file = (char*) "";
    Params::getInstance().model_name = model;
    Params::getInstance().ignore_identical_seqs = false; // keep the identical seqs
    
    if (intree != "") {
        // tree exists, then the resulting phylogenetic tree will be restricted to the input topology
        Params::getInstance().min_iterations = 0;
        Params::getInstance().stop_condition = SC_FIXED_ITERATION;
        Params::getInstance().start_tree = STT_USER_TREE;
        Params::getInstance().intree_str = intree;
    }

    if (in_options != NULL) {
        // assign the input options to Params
        in_options->set_params(Params::getInstance());
    }

    if (rand_seed == 0)
        rand_seed = make_new_seed();
    Params::getInstance().ran_seed = rand_seed;
    // cout << "Seed: " << Params::getInstance().ran_seed << endl << flush;
    init_random(Params::getInstance().ran_seed);

    string out_prefix_str = prog + "_" + convertIntToString(rand_seed);
    Params::getInstance().out_prefix = (char *) out_prefix_str.c_str();

    Checkpoint *checkpoint = new Checkpoint;
    string filename = (string)Params::getInstance().out_prefix +".ckp.gz";
    checkpoint->setFileName(filename);

    bool append_log = false;

    if (!Params::getInstance().ignore_checkpoint && fileExists(filename)) {
        checkpoint->load();
        if (checkpoint->hasKey("finished")) {
            if (checkpoint->getBool("finished")) {
                if (Params::getInstance().force_unfinished) {
                    if (MPIHelper::getInstance().isMaster())
                        cout << "NOTE: Continue analysis although a previous run already finished" << endl;
                } else {
                    delete checkpoint;
                    if (MPIHelper::getInstance().isMaster())
                        outError("Checkpoint (" + filename + ") indicates that a previous run successfully finished\n" +
                            "Use `-redo` option if you really want to redo the analysis and overwrite all output files.\n" +
                            "Use `--redo-tree` option if you want to restore ModelFinder and only redo tree search.\n" +
                            "Use `--undo` option if you want to continue previous run when changing/adding options."
                        );
                    else
                        exit(EXIT_SUCCESS);
                    exit(EXIT_FAILURE);
                }
            } else {
                append_log = true;
            }
        } else {
            if (MPIHelper::getInstance().isMaster())
                outWarning("Ignore invalid checkpoint file " + filename);
            checkpoint->clear();
        }
    }

    if (MPIHelper::getInstance().isWorker())
        checkpoint->setFileName("");

    _log_file = Params::getInstance().out_prefix;
    _log_file += ".log";
    startLogFile(append_log);
    time_t start_time;

    if (append_log) {
        cout << endl << "******************************************************"
             << endl << "CHECKPOINT: Resuming analysis from " << filename << endl << endl;
    }

    MPIHelper::getInstance().syncRandomSeed();

    signal(SIGABRT, &funcAbort);
    signal(SIGFPE, &funcAbort);
    signal(SIGILL, &funcAbort);
    signal(SIGSEGV, &funcAbort);
#if !defined WIN32 && !defined _WIN32 && !defined __WIN32__ && !defined WIN64
    signal(SIGBUS, &funcAbort);
#endif
    printCopyright(cout);

    char hostname[100];
#if defined WIN32 || defined _WIN32 || defined __WIN32__ || defined WIN64
    WSADATA wsaData;
    WSAStartup(MAKEWORD(2, 2), &wsaData);
    gethostname(hostname, sizeof(hostname));
    WSACleanup();
#else
    gethostname(hostname, sizeof(hostname));
#endif

    instruction_set = instrset_detect();
#if defined(BINARY32) || defined(__NOAVX__)
    instruction_set = min(instruction_set, (int)LK_SSE42);
#endif
    if (instruction_set < LK_SSE2) outError("Your CPU does not support SSE2!");
    bool has_fma3 = (instruction_set >= LK_AVX) && hasFMA3();

#ifdef __FMA__
    bool has_fma =  has_fma3;
    if (!has_fma) {
        outError("Your CPU does not support FMA instruction, quiting now...");
    }
#endif

    cout << "Host:    " << hostname << " (";
    switch (instruction_set) {
    case 0: cout << "x86, "; break;
    case 1: cout << "SSE, "; break;
    case 2: cout << "SSE2, "; break;
    case 3: cout << "SSE3, "; break;
    case 4: cout << "SSSE3, "; break;
    case 5: cout << "SSE4.1, "; break;
    case 6: cout << "SSE4.2, "; break;
    case 7: cout << "AVX, "; break;
    case 8: cout << "AVX2, "; break;
    default: cout << "AVX512, "; break;
    }
    if (has_fma3) cout << "FMA3, ";
    cout << (int)(((getMemorySize()/1024.0)/1024)/1024) << " GB RAM)" << endl;

    time(&start_time);
    cout << "Time:    " << ctime(&start_time);

    // increase instruction set level with FMA
    if (has_fma3 && instruction_set < LK_AVX_FMA)
        instruction_set = LK_AVX_FMA;

    Params::getInstance().SSE = min(Params::getInstance().SSE, (LikelihoodKernel)instruction_set);

    cout << "Kernel:  ";

    if (Params::getInstance().lk_safe_scaling) {
        cout << "Safe ";
    }

    if (Params::getInstance().pll) {
#ifdef __AVX__
        cout << "PLL-AVX";
#else
        cout << "PLL-SSE3";
#endif
    } else {
        if (Params::getInstance().SSE >= LK_AVX512)
            cout << "AVX-512";
        else if (Params::getInstance().SSE >= LK_AVX_FMA) {
            cout << "AVX+FMA";
        } else if (Params::getInstance().SSE >= LK_AVX) {
            cout << "AVX";
        } else if (Params::getInstance().SSE >= LK_SSE2){
            cout << "SSE2";
        } else
            cout << "x86";
    }

#ifdef _OPENMP
    if (Params::getInstance().num_threads >= 1) {
        omp_set_num_threads(Params::getInstance().num_threads);
        Params::getInstance().num_threads = omp_get_max_threads();
    }
//    int max_threads = omp_get_max_threads();
    int max_procs = countPhysicalCPUCores();
    cout << " - ";
    if (Params::getInstance().num_threads > 0)
        cout << Params::getInstance().num_threads  << " threads";
    else
        cout << "auto-detect threads";
    cout << " (" << max_procs << " CPU cores detected)";
    if (Params::getInstance().num_threads  > max_procs) {
        cout << endl;
        outError("You have specified more threads than CPU cores available");
    }
    // omp_set_nested(false); // don't allow nested OpenMP parallelism
    omp_set_max_active_levels(1);
#else
    if (Params::getInstance().num_threads != 1) {
        cout << endl << endl;
        outError("Number of threads must be 1 for sequential version.");
    }
#endif

    int num_procs = countPhysicalCPUCores();

    //cout << "sizeof(int)=" << sizeof(int) << endl;
    cout << endl << endl;
    
    // show msgs which are delayed to show
    cout << Params::getInstance().delay_msgs;

    cout.precision(3);
    cout.setf(ios::fixed);
    
    // checkpoint general run information
    checkpoint->startStruct("iqtree");
    int seed = Params::getInstance().ran_seed;
    CKP_SAVE(seed);
    CKP_SAVE(start_time);

    // check for incompatible version
    string version;
    stringstream sversion;
    sversion << iqtree_VERSION_MAJOR << "." << iqtree_VERSION_MINOR << iqtree_VERSION_PATCH;
    version = sversion.str();
    CKP_SAVE(version);
    checkpoint->endStruct();
    
    Params params = Params::getInstance();
    IQTree *tree;
    Alignment *alignment = new Alignment(names, seqs, params.sequence_type, params.model_name);
    bool align_is_given = true;
    ModelCheckpoint* model_info = NULL;
    if (model == "MF") {
        model_info = new ModelCheckpoint;
    }

    runPhyloAnalysis(params, checkpoint, tree, alignment, align_is_given, model_info);
    
    stringstream ss;
    if (model_info != NULL) {
        // output the modelfinder results in YAML format
        model_info->dump(ss);
    } else {
        // output the checkpoint in YAML format
        checkpoint->dump(ss);
    }
    
    alignment = tree->aln;
    delete tree;
    delete alignment;
    cleanup(params);

    time(&start_time);
    cout << "Date and Time: " << ctime(&start_time);
    try{
        delete checkpoint;
        if (model_info != NULL)
            delete model_info;
    }catch(int err_num){}

    finish_random();
    funcExit();

	string output = ss.str();
    if (output.length() > 0) {
    	char* out_cstr = new char[output.length()+1];
    	strcpy(out_cstr, output.c_str());
    	return out_cstr;
    } else {
        static char empty[] = "";
        return empty;
    }
}

// --------------------------------------------------
// Handle the input options of PiQTREE
// --------------------------------------------------

void input_options::set_params(Params& params) {
    ASSERT(flags.size() == values.size());
    int n = flags.size();
    for (int i = 0; i < n; i++) {
        if (flags[i] == "-keep-indent") {
            params.ignore_identical_seqs = false;
            cout << "params.ignore_identical_seqs = " << params.ignore_identical_seqs << endl;
        }
        else if (flags[i] == "-mset") {
            params.model_set = values[i];
            cout << "params.model_set = " << params.model_set << endl;
        }
        else if (flags[i] == "-mfreq") {
            int clen = values[i].length();
            if (clen > 0) {
                params.state_freq_set = new char[clen + 1];
                strcpy(params.state_freq_set, values[i].c_str());
            }
        }
        else if (flags[i] == "-mrate") {
            params.ratehet_set = values[i];
            cout << "params.ratehet_set = " << params.ratehet_set << endl;
        }
        else if (flags[i] == "-bb") {
            params.gbo_replicates = atoi(values[i].c_str());
            if (params.gbo_replicates < 1000)
                outError("#replicates must be >= 1000");
            params.consensus_type = CT_CONSENSUS_TREE;
            params.stop_condition = SC_BOOTSTRAP_CORRELATION;
        }
        else if (flags[0] == "-nt") {
            params.num_threads = atoi(values[i].c_str());
        }
    }
}

