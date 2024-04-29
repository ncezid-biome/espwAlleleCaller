# Joseph S. Wirth
# September 2023

from __future__ import annotations

class BlastHit():
    """class for storing blast hit data
    """
    def __init__(self, query:str, subject:str, length:int, coverage:float, pid:float, qstart:int, qend:int, sstart:int, send:int) -> BlastHit:
        """constructor

        Args:
            query (str): query
            subject (str): blast subject
            length (int): hitlen
            coverage (float): qcovhsp
            pid (float): pident
            qstart (int): query start
            qend(int): query end
            sstart (int): subject (hit) start
            send (int): subject (hit) end
        """
        # store inputs
        self.query = query
        self.subject = subject
        self.length = length
        self.coverage = coverage
        self.pid = pid
        self.query_coord = (qstart,qend)
        self.hit_coord = (sstart,send)
    
    def __str__(self) -> str:
        """convert BlastHit to string

        Returns:
            str: string representation
        """
        # constants
        GAP = " "*2
        EOL = "\n"
        
        return self.query + " vs " + self.subject + EOL + \
               GAP + "len:" + GAP + str(self.length) + EOL + \
               GAP + "cov:" + GAP + str(self.coverage) + EOL + \
               GAP + "pid:" + GAP + str(self.pid)
    
    def __repr__(self) -> str:
        """representation of BlastHit

        Returns:
            str: string representation
        """
        return str(self)

    def __gt__(self, other:BlastHit) -> bool:
        """ > operator; needed for sorting

        Args:
            other (BlastHit): a second BlastHit object

        Returns:
            bool: self > other
        """
        # prioritize coverage first
        if self.coverage > other.coverage:
            return True
        elif self.coverage < other.coverage:
            return False
        
        # prioritize percent identity next
        elif self.pid > other.pid:
            return True
        elif self.pid < other.pid:
            return False
        return False
    
    def __lt__(self, other:BlastHit) -> bool:
        """ < operator; needed for sorting

        Args:
            other (BlastHit): a second BlastHit object

        Returns:
            bool: self < other
        """
        # prioritize coverage first
        if self.coverage < other.coverage:
            return True
        elif self.coverage > other.coverage:
            return False
    
        # prioritize percent identity next
        elif self.pid < other.pid:
            return True
        elif self.pid > other.pid:
            return False
        return False

    def parseLine(line:str) -> BlastHit:
        """parses a single line from a blastn table into a BlastHit object

        Args:
            line (str): a line from a blastn table

        Returns:
            BlastHit: a BlastHit object
        """
        # constants
        EOL = "\n"
        SEP = "\t"
        QID = 0
        SID = 1
        LEN = 2
        QST = 3
        QED = 4
        SST = 5
        SED = 6
        QCH = 7
        PID = 8
    
        # chomp newline characters
        if line[-1] == EOL:
            line = line[:-1]
        
        # split the line into a row of columns
        row = line.split(SEP)
    
        # import data into a BlastHit object
        return BlastHit(row[QID],
                        row[SID],
                        int(row[LEN]),
                        float(row[QCH]),
                        float(row[PID]),
                        int(row[QST]),
                        int(row[QED]),
                        int(row[SST]),
                        int(row[SED]))
