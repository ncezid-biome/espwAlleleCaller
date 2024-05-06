#!/usr/bin/env python3

import getopt, glob, multiprocessing, os, shutil, subprocess, sys
from miscDirectory import MISC_DIR
sys.path.append(MISC_DIR)
from auxillary.blastn import _blastn
from auxillary.Parameters import Parameters
from downloadSRA import _runner as _downloadSrrs
from auxillary.ariba import _ariba, _buildAribaDb
from downloadAssemblies import _runner as _downloadAssemblies


__author__ = "Joseph S. Wirth"
__version__ = "2.0.4"


def _parseFile(fn:str) -> dict[str,tuple[str,str]]:
    """parses input file

    Args:
        fn (str): the filename to parse

    Returns:
        dict[str,tuple[str,str]]: key=genome key; val=(accession, srr id)
    """
    # constant
    SEP = "\t"
    
    # initialize output
    out = dict()
    
    # go through each line of the file
    with open(fn, 'r') as fh:
        for line in fh:
            # parse the line
            row = line.rstrip().split(SEP)
            
            # skip empty lines
            if row != ['']:
                # extract the data from the line
                key = row[0]
                acn = row[1]
                srr = row[2]
                
                # save the data
                out[key] = (acn,srr)
    
    return out


def _removeExistingAccessions(accessions:list[str], seqdir:str) -> None:
    """removes accession numbers for any fastas that already exist
    does not return; modifies input list

    Args:
        accessions (list[str]): a list of accession numbers
        seqdir (str): the directory where sequences are saved
    """
    # constant
    EXT = ".fna"
    
    # initialize variables
    remove = list()
    
    # get a list of all the filenames (no path)
    files = set(map(os.path.basename, glob.glob(os.path.join(seqdir, "*" + EXT))))
    
    # for each index (in reverse order for on-the-fly popping)
    for idx in range(len(accessions)-1, -1, -1):
        # remove any accession for which a file was already downloaded
        if (accessions[idx] + EXT) in files:
            remove.append(accessions.pop(idx))


def _findMissingAccessions(accessions:list[str], seqdir:str) -> set[str]:
    """identifies any accession numbers that were not downloaded

    Args:
        accessions (list[str]): a list of accession numbers
        seqdir (str): the directory containing the sequences

    Returns:
        set[str]: a collection of accessions that are missing
    """
    # constant
    PATTERN = "*fna"
    
    # get a list of all the existing accession numbers (based on downloaded files)
    exists = glob.glob(os.path.join(seqdir, PATTERN))
    exists = set(map(lambda x: os.path.basename(os.path.splitext(x)[0]), exists))
    
    # return all the accession numbers that do not have associated genomes
    return {x for x in accessions if x not in exists}


def _runOneBlastn(key:str, params:Parameters) -> tuple[str,str]:
    """runs a single blastn command

    Args:
        key (str): the genome key
        params (Parameters): a Parameters object

    Returns:
        tuple[str,str]: (key, allele)
    """
    # run the blast then link the result to the key
    allele,contig = _blastn(params)
    return (key, allele)


def _runAllBlasts(parsed:dict[str,tuple[str,str]], missing:set[str], seqdir:str, cpus:int) -> dict[str,str]:
    """runs all the blastn commands in parallel

    Args:
        parsed (dict[str,tuple[str,str]]): the dictionary produced by _parseFile
        missing (set[str]): the set produced by _findMissingAccessions
        seqdir (str): the directory where sequence files are saved
        cpus (int): the number of cpus for parallel processing

    Returns:
        dict[str,str]: key=genome key; val=allele
    """
    # constant
    EXT = ".fna"
    
    # initialize argument list
    args = list()

    # for each key
    for key in parsed.keys():
        # get the accession number
        accn = parsed[key][0]
        
        # add the associated fna file to the argument list if it isn't missing
        if accn not in missing:
            # get the Parameters object
            params = Parameters(os.path.join(seqdir, accn + EXT), '', '', 1, True, False)
            
            # create the blast directory if it doesn't already exist
            if not os.path.isdir(params._blastResultsDir):
                os.mkdir(params._blastResultsDir)
            
            # add the arguments to the list
            args.append((key, params))
    
    # run the blasts in parallel
    pool = multiprocessing.Pool(cpus)
    results = pool.starmap(_runOneBlastn, args)
    pool.close()
    pool.join()
    
    # convert the list of tuples to a dictionary before returning
    return dict(results)


def _runOneAriba(key:str, params:Parameters) -> tuple[str,str]:
    """runs a single ariba command

    Args:
        key (str): the genome key
        params (Parameters): a Parameters object

    Returns:
        tuple[str,str]: (genome key, allele)
    """
    # run one ariba then link the key to the allele
    allele = _ariba(params)
    return (key, allele)


def _runAllAribas(srrs:dict[str,str], seqdir:str, cpus:int) -> dict[str,str]:
    """runs all ariba commands in parallel

    Args:
        srrs (dict[str,str]): key=genome key; val=srr id
        seqdir (str): the directory containing sequence files
        cpus (int): the number of cpus for parallel processing

    Returns:
        dict[str,str]: key=genome key; val=allele
    """
    # constants
    EXT_1 = "_1.fastq"
    EXT_2 = "_2.fastq"
    
    # initialize variables
    args = list()
    built = False
    
    # for each key/srr pair
    for key,srr in srrs.items():
        # get the read file names and create a params object
        r1 = os.path.join(seqdir, srr + EXT_1)
        r2 = os.path.join(seqdir, srr + EXT_2)
        params = Parameters(srr, r1, r2, 1, True, False)
        
        # build the ariba database only once
        if not built:
            _buildAribaDb(params)
            built = True
        
        # create the ariba results directory if it does not already exist
        if not os.path.isdir(params._aribaResultsDir):
            os.mkdir(params._aribaResultsDir)
        
        # add the key and parameters to the argument list
        args.append((key, params))
    
    # run ariba in parallel
    pool = multiprocessing.Pool(cpus)
    results = pool.starmap(_runOneAriba, args)
    pool.close()
    pool.join()
    
    # convert the list of tuples to a dictionary before returning
    return dict(results)


def _writeResults(results:dict[str,str], fn:str) -> None:
    """writes the results to file

    Args:
        results (dict[str,str]): key=genome key; val=allele
        fn (str): output filename
    """
    #constants
    SEP = "\t"
    EOL = "\n"
    
    # open the file
    with open(fn, 'w') as fh:
        # write each key and its allele to the file
        for key,allele in results.items():
            fh.write(f'{key}{SEP}{allele}{EOL}')


def _cleanup(delete:bool, seqdir:str) -> None:
    """deletes intermediate files

    Args:
        delete (bool): indicates if the sequence directory should be deleted
        seqdir (str): the sequence directory
    """
    # delete the directory containing the blast results
    shutil.rmtree(os.path.join(os.curdir, Parameters._BLAST_DIR))
    
    # delete the directory containing the blast databases
    shutil.rmtree(os.path.join(os.curdir, Parameters._BLAST_DB))
    
    # delete the directory containing the ariba database
    shutil.rmtree(os.path.join(os.curdir, Parameters._ARIBA_DB))
    
    # delete the direcotry containing the ariba results
    shutil.rmtree(os.path.join(os.curdir, Parameters._ARIBA_DIR))
    
    # delete the sequence files if requested
    if delete:
        shutil.rmtree(seqdir)


def _runner(infn:str, email:str, seqdir:str, outfn:str, cpus:int, delete:bool) -> None:
    """main runner function

    Args:
        infn (str): input filename
        email (str): email address
        seqdir (str): directory to save sequence files to
        outfn (str): output filename
        cpus (int): number of cpus for parallel processing
        delete (bool): indicates if sequence files should be deleted
    """
    # constants
    FORMAT = "fasta"
    RENAME = True
    ABSENT = "absent"
    
    # parse the file
    parsed = _parseFile(infn)
    
    # get a list of the accession numbers
    accessions = [x for x,y in parsed.values()]
    
    # create the sequence directory if it does not exist
    if not os.path.isdir(seqdir):
        os.mkdir(seqdir)
    
    # remove any accessions from genomes that have already been downloaded
    _removeExistingAccessions(accessions, seqdir)
    
    # download the assemblies and identify those that failed to download
    _downloadAssemblies(email, accessions, FORMAT, seqdir, RENAME)
    missing = _findMissingAccessions(accessions, seqdir)
    
    # run all the blasts
    print('running blastn ... ', end='', flush=True)
    results = _runAllBlasts(parsed, missing, seqdir, cpus)
    print('done')
    
    # get a list of srrs that need to be downloaded for ariba
    srrs = {x:parsed[x][1] for x,y in results.items() if y == ABSENT}
    missingSrr = {x:parsed[x][1] for x in parsed.keys() if parsed[x][0] in missing}
    srrs.update(missingSrr)
    
    # only process reads if required
    if len(srrs) > 0:
        # download srrs
        _downloadSrrs(list(srrs.values()), seqdir, False, cpus, 1)

        # run ariba and update the results
        print('running ariba ... ', end='', flush=True)
        results.update(_runAllAribas(srrs, seqdir, cpus))
        print('done')
    
    # write the results to file
    _writeResults(results, outfn)
    
    # remove intermediate files
    print('removing files ... ', end='', flush=True)
    _cleanup(delete, seqdir)
    print('done')


def _parseArgs() -> tuple[str,str,str,str,int,bool,bool]:
    """parses command line arguments

    Raises:
        EnvironmentError: incompatible python version
        EnvironmentError: incompatible biopython version
        EnvironmentError: biopython not installed
        EnvironmentError: ncbi-blast+ not found
        EnvironmentError: incompatible ariba version
        EnvironmentError: ariba not found
        EnvironmentError: sratoolkit not found
        FileNotFoundError: input file does not exist
        ValueError: output file not writeable
        ValueError: invalid value for num_threads
        BaseException: missing required arguments

    Returns:
        tuple[str,str,str,str,int,bool]: input filename, email, output filename, sequence directory, number threads, delete seqs, help requested
    """
    # flags
    IN_FLAGS = ("-i", "--in")
    OUT_FLAGS = ("-o", "--out")
    HELP_FLAGS = ("-h", "--help")
    EMAIL_FLAGS = ("-e", "--email")
    CHECK_FLAGS = ("-c", "--check_env")
    SEQDIR_FLAGS = ("-s", "--seq_dir")
    DELETE_FLAGS = ("-d", "--delete")
    THREADS_FLAGS = ("-n", "--num_threads")
    VERSION_FLAGS = ("-v", "--version")
    SHORT_OPTS = IN_FLAGS[0][-1] + ":" + \
                 EMAIL_FLAGS[0][-1] + ":" + \
                 OUT_FLAGS[0][-1] + ":" + \
                 SEQDIR_FLAGS[0][-1] + ":" + \
                 THREADS_FLAGS[0][-1] + ":" + \
                 DELETE_FLAGS[0][-1] + \
                 HELP_FLAGS[0][-1] + \
                 VERSION_FLAGS[0][-1] + \
                 CHECK_FLAGS[0][-1]
    LONG_OPTS = (IN_FLAGS[1][2:] + "=",
                 OUT_FLAGS[1][2:] + "=",
                 HELP_FLAGS[1][2:],
                 EMAIL_FLAGS[1][2:] + "=",
                 CHECK_FLAGS[1][2:],
                 DELETE_FLAGS[1][2:],
                 SEQDIR_FLAGS[1][2:] + "=",
                 THREADS_FLAGS[1][2:] + "=",
                 VERSION_FLAGS[1][2:])
    
    # default values
    DEF_OUT = os.path.join(os.curdir, "espw_alleles.tsv")
    DEF_THREADS = 1
    DEF_HELP = False
    DEF_SEQDIR = os.path.join(os.curdir, "seqs")
    DEF_DELETE = False
    
    # messages
    ERR_MSG_1 = "input file does not exist"
    ERR_MSG_2 = "cannot write output file"
    ERR_MSG_3 = "invalid value for number of threads"
    ERR_MSG_4 = "must provide all required arguments"
    
    def printHelp():
        GAP = " "*4
        EOL = "\n"
        SEP = ", "
        WIDTH = 21
        DEFAULT = " (default: "
        HELP_MSG = f"{EOL}Determines the espW allele in an E. coli genome{EOL}" + \
                   f"{GAP}{__author__}, 2024{EOL*2}" + \
                   f"usage:{EOL}" + \
                   f"{GAP}{os.path.basename(__file__)} [-{SHORT_OPTS.replace(':', '')}]{EOL*2}" + \
                   f"required arguments:{EOL}" + \
                   f"{GAP}{IN_FLAGS[0] + SEP + IN_FLAGS[1]:<{WIDTH}}[file] filename of a tab-separated file with three columns and no headers: key, ncbi accession, srr id{EOL}" + \
                   f"{GAP}{EMAIL_FLAGS[0] + SEP + EMAIL_FLAGS[1]:<{WIDTH}}[str] email address (used to query NCBI){EOL*2}" + \
                   f"optional arguments:{EOL}" + \
                   f"{GAP}{OUT_FLAGS[0] + SEP + OUT_FLAGS[1]:<{WIDTH}}[file] filename to write the output{DEFAULT}'{DEF_OUT}'){EOL}" + \
                   f"{GAP}{SEQDIR_FLAGS[0] + SEP + SEQDIR_FLAGS[1]:<{WIDTH}}[directory] the directory where sequence files will be downloaded {DEFAULT}'{DEF_SEQDIR}'){EOL}" + \
                   f"{GAP}{THREADS_FLAGS[0] + SEP + THREADS_FLAGS[1]:<{WIDTH}}[int] the number of threads to use for parallel processing{DEFAULT}{DEF_THREADS}){EOL}" + \
                   f"{GAP}{DELETE_FLAGS[0] + SEP + DELETE_FLAGS[1]:<{WIDTH}}delete sequence files after finishing{DEFAULT}{DEF_DELETE}){EOL*2}" + \
                   f"troubleshooting:{EOL}" + \
                   f"{GAP}{HELP_FLAGS[0] + SEP + HELP_FLAGS[1]:<{WIDTH}}print this help message{EOL}" + \
                   f"{GAP}{VERSION_FLAGS[0] + SEP + VERSION_FLAGS[1]:<{WIDTH}}print the version{EOL}" + \
                   f"{GAP}{CHECK_FLAGS[0] + SEP + CHECK_FLAGS[1]:<{WIDTH}}check that all dependencies are installed{EOL}"
        
        print(HELP_MSG)
    
    def checkEnvironment() -> None:
        # commands
        MAKEDB_CMD = ["makeblastdb", "-help"]
        BLASTN_CMD = ["blastn", "-help"]
        ARIBA_CMD = ["ariba", "version"]
        SRA_CMD = ["fasterq-dump", "-h"]
        ARIBA_VER = (2, 12)
        BIO_VER = (1, 81)
        PY_VER = (3, 9)
        
        # messages
        BAD_BIO = "biopython is not installed"
        BAD_VER = ' version is incompatible (requires '
        BAD_SRA = 'sratoolkit is not installed or not in the path'
        BAD_NCBI = "ncbi-blast+ is not installed or not in the path"
        BAD_ARIBA = "ariba is not installed or not in the path"
        SUCCESS = "\nenvironment is suitable\n"
        
        # check python version
        if sys.version_info.major != PY_VER[0] or sys.version_info.minor < PY_VER[1]:
            raise EnvironmentError(f"python {BAD_VER}{'.'.join(map(str, PY_VER))} or above)")
        
        # check biopython install
        try:
            import Bio
            
            # check biopython version
            vers = tuple(map(int, Bio.__version__.split('.')))
            if vers[0] < BIO_VER[0] or (vers[0] == BIO_VER[0] and vers[1] < BIO_VER[1]):
                raise EnvironmentError(f"'Bio'{BAD_VER}{'.'.join(map(str,BIO_VER))} or above)")
        except Exception as e:
            if not isinstance(e, EnvironmentError):
                raise EnvironmentError(BAD_BIO)
        
        # check ncbi-blast+
        try:
            subprocess.run(MAKEDB_CMD, capture_output=True, check=True)
            subprocess.run(BLASTN_CMD, capture_output=True, check=True)
        except:
            raise EnvironmentError(BAD_NCBI)
        
        # check ariba
        try:
            result = subprocess.run(ARIBA_CMD, capture_output=True, check=True)
            
            # check ariba version
            vers = tuple(map(int, result.stdout.decode().split('\n')[0][len('ARIBA version: '):].split('.')))
            if vers[0] < ARIBA_VER[0] or (vers[0]==ARIBA_VER[0] and vers[1] < ARIBA_VER[1]):
                raise EnvironmentError(f"'ariba'{BAD_VER}{'.'.join(map(str,ARIBA_VER))})")
        
        except Exception as e:
            if not isinstance(e, EnvironmentError):
                raise EnvironmentError(BAD_ARIBA)

        # check sratoolkit
        try:
            subprocess.run(SRA_CMD, capture_output=True, check=True)
        except:
            raise EnvironmentError(BAD_SRA)
        
        # print success message
        print(SUCCESS)
    
    # set default values
    infn = None
    email = None
    outfn = DEF_OUT
    seqdir = DEF_SEQDIR
    cpus = DEF_THREADS
    delete = DEF_DELETE
    helpRequest = DEF_HELP
    
    # give help if requested
    if HELP_FLAGS[0] in sys.argv or HELP_FLAGS[1] in sys.argv or len(sys.argv) == 1:
        helpRequest = True
        printHelp()
    
    # print the version if requested
    elif VERSION_FLAGS[0] in sys.argv or VERSION_FLAGS[1] in sys.argv:
        helpRequest = True
        print(f"v{__version__}")
    
    # check installation if requested
    elif CHECK_FLAGS[0] in sys.argv or CHECK_FLAGS[1] in sys.argv:
        helpRequest = True
        checkEnvironment()
    
    # parse command line arguments
    else:
        opts,args = getopt.getopt(sys.argv[1:], SHORT_OPTS, LONG_OPTS)
        for opt,arg in opts:
            # get input filename
            if opt in IN_FLAGS:
                if not os.path.exists(arg):
                    raise FileNotFoundError(ERR_MSG_1)
                infn = arg
            
            # get email address
            elif opt in EMAIL_FLAGS:
                email = arg
            
            # get output filename
            elif opt in OUT_FLAGS:
                try:
                    fh = open(arg, 'w')
                    fh.close()
                except:
                    raise ValueError(ERR_MSG_2)
                outfn = arg
            
            # get sequence directory
            elif opt in SEQDIR_FLAGS:
                seqdir = arg
            
            # get number of threads
            elif opt in THREADS_FLAGS:
                try:
                    cpus = int(arg)
                except ValueError:
                    raise ValueError(ERR_MSG_3)
            
            # determine if sequence files should be deleted
            elif opt in DELETE_FLAGS:
                delete = True                
    
        # make sure all required arguments were specified
        if None in (infn, email):
            raise BaseException(ERR_MSG_4)
        
    return infn, email, outfn, seqdir, cpus, delete, helpRequest


def _main():
    """entry point to the program
    """
    # parse command line arguments
    infn, email, outfn, seqdir, cpus, delete, helpRequest = _parseArgs()
    
    # only run the program if help was not requested
    if not helpRequest:
        _runner(infn, email, seqdir, outfn, cpus, delete)


if __name__ == "__main__":
    _main()
