import numpy as np
import pysam
import re

def extract_error_probs(fq_path):
    """
    Extract the error probabilities from
    a .fastq file
    
    params
        fq_path : str
            Path to the fastq file.
    
    returns
        ps : list of arrays
            List of arrays containing error
            probabilities for each read
    
    """
    
    with open(fq_path, "r") as fastq:
        ps = []
        for line in fastq:
            if line[0] == "@":
                j = 0
            if j == 3:
                ascii_score = line.rstrip()
                q_score = np.array([ord(c) for c in ascii_score]) - 33
                p_error = 10 ** (q_score / -10)
                ps.append(p_error)
            j += 1
    
    return(ps)

def sam_positional_error(sam_path):
    """
    Extract error probabilities and their positions
    from a SAM file
    
    params
        sam_path : str
            Path to SAM file.
            
    returns
        positions : list of ndarray, int, shape(n_reads, )
            List of numpy arrays. Each array
            encodes positions of aligned bases
            for a read. Note that -998 indicates clipped,
            -999 indicates deletion.
        error_probs : list of ndarray, float, shape(n_reads, )
            List of numpy arrays. Each array
            encodes error probabilities aligned bases
            for a read.
    
    """
    # Open SAM
    with open(sam_path, "r") as sam:

        # Prepare storage
        positions = []
        error_probs = []

        for alignment in sam:
            # Parse alignment & get relevant fields
            fields = alignment.split("\t")
            start = int(fields[3])
            cigar = fields[5]
            quals = fields[10]

            # Compute
            alignment_error_probs = calc_error_probabilities(quals)
            alignment_positions = cigar_positions(start, cigar)

            # Store
            positions.append(alignment_positions)
            error_probs.append(alignment_error_probs)
            
    return positions, error_probs

def calc_error_probabilities(quals):
    """
    Calculate error probabilities from
    a string of ASCII characters `quals` that
    encode error probabilities
    
    params
        quals : str
            Error probabilities encoded in ASCII.
    
    returns
        error_probs : ndarray, float, shape(read_length, )
            Error probabilities in a
            numpy array.
    
    """
    qs = np.array([ord(c) for c in quals]) - 33
    error_probs = 10 ** (qs / -10)
    return error_probs

def cigar_positions(start, cigar):
    """
    Get alignment positions for each base in a read
    given a `start` position and a `cigar` string
    
    params
        start : int
            The start position of the aligment.
        cigar : str
            The CIGAR string for the alignment.
    
    returns
        pos : ndarray, int, shape (read_length,)
            A genomic position for each base in the
            read.
    
    """
    pos = []
    i = start
    tags = re.findall("\d*[MIDSH]", cigar)  # parse the cigar string into 'tags'
    
    for tag in tags:
        n = int(tag[:-1]) # this is the integer value of the tag
        m = tag[-1] # this is the 'm'ethod: M, D, I, S, H ...
        if m == 'M':  # match
            pos.extend(np.arange(i, i + n))
        elif m == "I":  # insertion
            pos.extend(np.repeat(-999, n))
        elif m == "D":  # deletion
            pass
        elif m == "S":  # clip
            pos.extend(np.repeat(-998, n))
        elif m == "H":
            pass
        else:
            Print("Tag %s not identified." % m)
        i += n  # we have moved forward n positions
    
    return np.array(pos)
