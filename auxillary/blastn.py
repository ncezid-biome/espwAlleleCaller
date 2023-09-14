# Joseph S. Wirth
# September 2023
from auxillary.Parameters import Parameters
from auxillary.BlastHit import BlastHit
import os, subprocess


def _blastn(params:Parameters) -> str:
    """runs BLASTn to detect espW alleles using the values specified in a Parameters object

    Args:
        params (Parameters): a Parameters object

    Returns:
        str: the espW allele detected
    """
    __makeBlastDb(params)
    blastFn = __runBlast(params)
    allele = __parseBlastTable(blastFn)
    return allele


def __makeBlastDb(params:Parameters) -> None:
    """makes a blastn database for the input genome using the values specified in a Parameters object

    Args:
        params (Parameters): a Parameters object
    """
    # constants
    CMD_A = "makeblastdb"
    CMD_B = ("-dbtype", "nucl", "-out")
    DB_DIR = os.path.join(os.path.abspath(os.curdir), "blastdb")
    
    # get the blast db name
    params._blastdb = os.path.join(DB_DIR, os.path.splitext(os.path.basename(params.fna))[0])
    
    # build the command
    cmd = [os.path.join(params._blastExeDir, CMD_A)]
    cmd.extend(CMD_B)
    cmd.extend([params._blastdb, "-in", params.fna])
    
    # make the blastdb
    subprocess.run(cmd, check=True, capture_output=True)


def __runBlast(params:Parameters) -> None:
    """runs blastn using the values specified in a Parameters object

    Args:
        params (Parameters): a Parameters object
    """
    # constants
    FN_SUFFIX = "_vs_espW.blastn"
    BLASTN = "blastn"
    CMD = ("-ungapped", "-perc_identity", "90", "-qcov_hsp_perc", "35", "-max_target_seqs",
           "10000", "-outfmt", '6 qseqid sseqid length qstart qend sstart send qcovhsp pident')
    
    # determine the output file
    params._blastFn = os.path.join(os.curdir, os.path.splitext(os.path.basename(params.fna))[0] + FN_SUFFIX)
    
    # build the blastn command
    cmd = [os.path.join(params._blastExeDir, BLASTN)]
    cmd.extend(CMD)
    cmd.extend(["-query", params._espwFna])
    cmd.extend(['-db', params._blastdb])
    cmd.extend(['-num_threads', str(params.threads)])
    cmd.extend(['-out', params._blastFn])
    
    # run blastn
    subprocess.run(cmd, check=True, capture_output=True)


def __parseBlastTable(fn:str) -> str:
    """parses a blastn table and determines which espW allele is present

    Args:
        fn (str): a blastn table

    Returns:
        str: the espW allele detected
    """
    # constant
    NULL = 0

    # initialize a list of hits
    allHitsL:list[BlastHit] = list()
    
    # parse each line and add it to the list of hits
    with open(fn, 'r') as fh:
        for line in fh:
            allHitsL.append(BlastHit.parseLine(line))
    
    # no hits means espW is absent
    if allHitsL == []:
        allele = 'absent'
    
    # otherwise, choose the best hit's allele
    else:
        allele = max(allHitsL).query
    
    return allele
