import numpy as np

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
