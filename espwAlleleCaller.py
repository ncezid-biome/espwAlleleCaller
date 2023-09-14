# Joseph S. Wirth
# September 2023

from auxillary.Parameters import Parameters
from auxillary.blastn import _blastn
from auxillary.ariba import _ariba
import os

def __main() -> None:
    """main runner function
         * searches for espW using blastn
         * if that fails, searches using ariba
         * prints the allele to stdout
    """
    # constant
    ABSENT = "absent"
    
    # parse command line arguments
    params = Parameters.parseArgs()
    
    # only do work if help wasn't requested
    if not params.helpRequested:
        # attempt to determine the allele using blastn
        allele = _blastn(params)
    
        # if an allele wasn't found, attempt to determine with ariba
        if allele == ABSENT:
            allele = _ariba(params)

        # delete files if requested
        if params.cleanup:
            os.rmdir(os.path.dirname(params._blastdb))
            os.remove(params._blastFn)
            os.rmdir(params._aribaDir)
        
        print(allele)


if __name__ == "__main__":
    allele = __main()
