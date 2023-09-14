# Joseph S. Wirth
# September 2023

from auxillary.Parameters import Parameters
from Bio.SeqRecord import SeqRecord
import gzip, os, subprocess
from Bio.Seq import Seq
from Bio import SeqIO


def _ariba(params:Parameters) -> str:
    """runs ARIBA and processes the output

    Args:
        params (Parameters): a Parameters object

    Returns:
        str: the detected espW allele
    """
    __buildAribaDb(params)
    aribaFn = __runAriba(params)
    allele = __processAriba(aribaFn)
    return allele


def __buildAribaDb(params:Parameters) -> None:
    """builds the ariba database using the values specified in a Parameters object

    Args:
        params (Parameters): a Parameters object
    """
    # constant
    CMD = ("prepareref", "--all_coding", "no", "-f")
    
    # build the command
    cmd = [params._aribaExe]
    cmd.extend(CMD)
    cmd.extend([params._espwFna, params._aribaDb])

    # run the command
    subprocess.run(cmd, check=True, capture_output=True)


def __runAriba(params:Parameters) -> str:
    """runs ARIBA using the values specified in a Parameters object

    Args:
        params (Parameters): a Parameters object

    Returns:
        str: the file to be processed
    """
    # constants
    CMD = ("run", "--threads")
    OUT_FN = "assembled_seqs.fa.gz"
    
    # determine the output directory
    params._aribaDir = os.path.join(os.curdir, os.path.splitext(os.path.basename(params.fna))[0])
  
    # build the command
    cmd = [params._aribaExe]
    cmd.extend(CMD)
    cmd.extend([str(params.threads), params._aribaDb, params.reads[0], params.reads[1], params._aribaDir])
    
    # run the command
    subprocess.run(cmd, check=True, capture_output=True)

    # return the file to be processed    
    return os.path.join(params._aribaDir, OUT_FN)
    

def __processAriba(aribaFn:str) -> str:
    """processes the output of an ARIBA run

    Args:
        aribaFn (str): the ARIBA output filename

    Raises:
        RuntimeError: multiple alleles detected
        RuntimeError: multiple alleles detected

    Returns:
        str: the espW allele detected
    """
    # constants
    FORMAT = 'fasta'
    INSERTION = Seq("gaaaaaaaaag")
    FULL =      Seq("gaaaaaaaag")
    DELETION =  Seq("gaaaaaaag")
    INS = "insertion"
    FUL = "full length"
    DEL = "deletion"
    AMB = "ambiguous"
    ABS = "absent"
    ERR_MSG = "multiple alleles detected:\n"
    
    # initialize looping variables
    alleles = set()
    rec:SeqRecord
        
    # go through the records in the file
    with gzip.open(aribaFn, 'rt') as fh:
        for rec in SeqIO.parse(fh, FORMAT):
            # check for each allele; add it to the set
            if INSERTION in rec.seq.lower():
                alleles.add(INS)
            elif FULL in rec.seq.lower():
                alleles.add(FUL)
            elif DELETION in rec.seq.lower():
                alleles.add(DEL)
            else:
                alleles.add(AMB)
    
    # if the set is empty then no allele was found
    if alleles == set():
        return ABS
    
    # if only allele is present then return it
    elif len(alleles) == 1:
        return alleles.pop()

    # ignore ambiguous if multiple alleles present
    elif AMB in alleles:
        alleles.remove(AMB)
        
        # if one other allele, then return it
        if len(alleles) == 1:
            return alleles.pop()
        
        # raise error if multiple alleles
        else:
            raise RuntimeError(ERR_MSG + "\n".join(alleles))
    
    # raise error if multiple alleles
    else:
        raise RuntimeError(ERR_MSG + "\n".join(alleles))
