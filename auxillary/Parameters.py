# Joseph S. Wirth
# September 2023
from __future__ import annotations
from Bio import SeqIO
import getopt, os, sys

class Parameters():
    """class for storing and accessing parameters
    """
    # global constants
    _ARIBA_DB = "ariba_db"
    _ARIBA_DIR = "ariba_results"
    _BLAST_DB = "blast_db"
    _BLAST_DIR = "blast_results"
    __DEFAULT_CLEANUP = True
    __DEFAULT_HELP = False
    __DEFAULT_THREADS = 1
    __ESPW_FNA = "espW_alleles.fna"
    
    def __init__(self, fna:str, read1:str, read2:str, threads:int, cleanup:bool, helpRequested:bool) -> Parameters:
        """constructor for the Parameters class

        Args:
            fna (str): a nucleotide fasta filename and path
            read1 (str): a forward reads filename and path
            read2 (str): a reverse reads filename and path
            threads (int): the number of threads to use for parallel processes
            cleanup (bool): indicates whether intermediate files should be removed
            helpRequested (bool): indicates if help was requested

        Returns:
            Parameters: a Parameters object
        """
        # store inputs
        self.fna:str = os.path.abspath(fna)
        self.reads:tuple[str] = tuple(map(os.path.abspath, (read1, read2)))
        self.threads:int = threads
        self.cleanup:bool = cleanup
        self.helpRequested:bool = helpRequested
        
        # determine the `data` directory
        dataDir = os.path.join(os.path.dirname(__file__), "..", "data")
        
        # import values from the params file
        self._aribaDb:str = os.path.join(os.curdir, Parameters._ARIBA_DB)
        self._aribaResultsDir:str = os.path.join(os.curdir, Parameters._ARIBA_DIR)
        self._blastResultsDir:str = os.path.join(os.curdir, Parameters._BLAST_DIR)
        self._blastDbDir:str = os.path.join(os.curdir, Parameters._BLAST_DB)
        self._espwFna:str = os.path.join(dataDir, Parameters.__ESPW_FNA)
        
        # initialize a few other variables
        self._blastDb:str = None
        self._blastFn:str = None
        self._aribaDir:str = None

    
    def __isFasta(fn:str) -> bool:
        """checks if an ipnut file is in fact a fasta

        Args:
            fn (str): the file to check

        Returns:
            bool: indicates if the file is formatted properly
        """
        # true if one or more records could be parsed; otherwise false
        for rec in SeqIO.parse(fn, 'fasta'):
            return True
        return False
    
    
    def parseArgs() -> Parameters:
        """parses command line arguments and builds a Parameters object from them

        Raises:
            RuntimeError: reads were not specified properly
            FileNotFoundError: read 1 was not found
            FileNotFoundError: read 2 was not found
            FileNotFoundError: fasta was not found
            ValueError: input fasta is improperly formatted
            ValueError: invalid input for threads
            BaseException: one or more required arguments omitted

        Returns:
            Parameters: a Parameters object
        """
        # flags
        READS_FLAGS = ("-r", "--reads")
        FNA_FLAGS = ("-f", "--fasta")
        THREADS_FLAGS = ("-n", "--num_threads")
        CLEANUP_FLAGS = ("-s", "--save_files")
        HELP_FLAGS = ("-h", "--help")
        SHORT_OPTS = READS_FLAGS[0][-1] + ":" + \
                     FNA_FLAGS[0][-1] + ":" + \
                     THREADS_FLAGS[0][-1] + ":" + \
                     CLEANUP_FLAGS[0][-1] + \
                     HELP_FLAGS[0][-1]
        LONG_OPTS = [READS_FLAGS[1][2:] + "=",
                     FNA_FLAGS[1][2:] + "=",
                     THREADS_FLAGS[1][2:] + "=",
                     CLEANUP_FLAGS[1][2:],
                     HELP_FLAGS[1][2:]]
        
        # messages
        ERR_MSG_1  = "must specify exactly two reads separated by a comma"
        ERR_MSG_2  = "file is not a fasta: "
        ERR_MSG_3A = "invalid value for '"
        ERR_MSG_3B = "':   "
        ERR_MSG_4  = "must provide all required arguments"
        IGNORE_MSG = "ignoring unused argument: "
        
        # helper function to print help
        def printHelp():
            """helper function for printing the help message
            """
            GAP = 4*" "
            EOL = "\n"
            SEP = ", "
            MSG = EOL + "determines the allele of espW from genomic sequence data" + EOL + \
                   GAP + "Joseph S. Wirth, 2023" + EOL*2 + \
                   "usage:" + EOL + \
                   GAP + "espwAlleleCaller [-" + \
                         READS_FLAGS[0][-1] + \
                         FNA_FLAGS[0][-1] + \
                         THREADS_FLAGS[0][-1] + \
                         CLEANUP_FLAGS[0][-1] + \
                         HELP_FLAGS[0][-1] + "]" + EOL*2 + \
                   "required arguments:" + EOL + \
                   GAP + f"{READS_FLAGS[0] + SEP + READS_FLAGS[1]:<21}{'input reads as a pair of files separated by a comma (read1,read2)'}" + EOL + \
                   GAP + f"{FNA_FLAGS[0] + SEP + FNA_FLAGS[1]:<21}{'input nucleotide fasta file'}" + EOL*2 + \
                   "optional arguments:" + EOL + \
                   GAP + f"{THREADS_FLAGS[0] + SEP + THREADS_FLAGS[1]:<21}{'the number of threads to use (default: 1)'}" + EOL + \
                   GAP + f"{CLEANUP_FLAGS[0] + SEP + CLEANUP_FLAGS[1]:<21}{'save intermediate files generated by this program (default: delete files)'}" + EOL + \
                   GAP + f"{HELP_FLAGS[0] + SEP + HELP_FLAGS[1]:<21}{'print this message'}" + EOL
            
            print(MSG)
                    
        # set default values
        fna = None
        read1 = None
        read2 = None
        cleanup = Parameters.__DEFAULT_CLEANUP
        threads = Parameters.__DEFAULT_THREADS
        helpRequested = Parameters.__DEFAULT_HELP

        # print help message if requested; update values
        if len(sys.argv) == 1 or HELP_FLAGS[0] in sys.argv or HELP_FLAGS[1] in sys.argv:
            printHelp()
            fna = ''
            read1 = ''
            read2 = ''
            helpRequested = True
        
        # otherwise parse command line arguments
        else:
            opts,args = getopt.getopt(sys.argv[1:], SHORT_OPTS, LONG_OPTS)
            for opt,arg in opts:
                # parse read flags
                if opt in READS_FLAGS:
                    # make sure exactly two reads were provided
                    try:
                        read1,read2 = arg.split(",")
                    except:
                        raise RuntimeError(ERR_MSG_1)
                    
                    # check that both files exist
                    if not os.path.exists(read1):
                        raise FileNotFoundError(read1)
                    if not os.path.exists(read2):
                        raise FileNotFoundError(read2)
                
                # parse fasta flags
                elif opt in FNA_FLAGS:
                    # make sure the file exists
                    if not os.path.exists(arg):
                        raise FileNotFoundError(arg)
                    
                    # make sure the file is a fasta
                    if not Parameters.__isFasta(arg):
                        raise ValueError(ERR_MSG_2 + arg)
                    
                    # store value
                    fna = arg

                # parse threads flags
                elif opt in THREADS_FLAGS:
                    # make sure the number is an integer
                    try:
                        threads = int(arg)
                    except:
                        raise ValueError(ERR_MSG_3A + opt + ERR_MSG_3B + arg)
                
                # parse cleanup flags
                elif opt in CLEANUP_FLAGS:
                    cleanup = False
                
                # ignore all other arguments
                else:
                    print(IGNORE_MSG + opt + " " + arg)
        
        # make sure all required arguments were specified
        if None in {read1, read2, fna}:
            raise BaseException(ERR_MSG_4)
        
        return Parameters(fna, read1, read2, threads, cleanup, helpRequested)
