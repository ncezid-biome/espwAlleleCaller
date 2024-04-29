# Joseph S. Wirth
# September 2023

__author__ = "Joseph S. Wirth"
__version__ = "1.0.0"

from auxillary.Parameters import Parameters
from auxillary.blastn import _blastn
from auxillary.ariba import _ariba
import os, shutil


def __cleanup(params:Parameters) -> None:
    """delete any files created by this program using values specified in a Parameters object

    Args:
        params (Parameters): a Parameters object
    """
    # remove blastdb
    if params._blastDb is not None:
        if os.path.exists(os.path.dirname(params._blastDb)):
            shutil.rmtree(os.path.dirname(params._blastDb))
    
    # remove blastn result
    if params._blastFn is not None:
        if os.path.exists(params._blastFn):
            os.remove(params._blastFn)
    
    # remove ariba results
    if params._aribaDir is not None:
        if os.path.exists(params._aribaDir):
            shutil.rmtree(params._aribaDir)
    
    # remove ariba db
    if params._aribaDb is not None:
        if os.path.exists(params._aribaDb):
            shutil.rmtree(params._aribaDb)


def __main() -> None:
    """main runner function
         * searches for espW using blastn
         * if that fails, searches using ariba
         * prints the allele to stdout
    """
    # constant
    ABSENT = "absent"
    SEP = "\t"
    
    # parse command line arguments
    params = Parameters.parseArgs()
    
    # only do work if help wasn't requested
    if not params.helpRequested:
        # attempt to determine the allele using blastn
        allele,contig = _blastn(params)
    
        # if an allele wasn't found, attempt to determine with ariba
        if allele == ABSENT:
            allele = _ariba(params)

        # delete files if requested
        if params.cleanup:
            __cleanup(params)
        
        # print the allele; print the contig if blastn worked
        if contig is not None:
            print(allele + SEP + contig)
        else:
            print(allele)


if __name__ == "__main__":
    __main()
