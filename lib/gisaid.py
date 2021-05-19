# ------------------------------------------------------------------------------
# Functions to compute statistics required for GISAID subsmission
# JHendry, 2021/01/01
# ------------------------------------------------------------------------------

import os
import sys
import subprocess
import pandas as pd
import numpy as np


def load_coverage_files(sample_dir):
    """
    Load coverage files from ARCTIC output
    
    Parameters
        sample_dir : str
            Path to the sample directory. Within this
            directory there should be coverage files
            generated from the ARCTIC pipeline which
            have a '.depths' suffix.
            
    Returns
        df : DataFrame, shape (2, n_sites)
            A data frame of two columns containing
            genomic position 'pos' and coverage 'cov'
            data. Each position is representd only once;
            the coverage across multiple primer pools is summed.
    
    """
    
    # Identify coverage files
    cov_files = [fn for fn in os.listdir(sample_dir) if fn.endswith(".depths")]
    
    # Prepare storage
    dt = {
        "pos": [],
        "cov": []
    }

    # Iterate over pools
    for cov_fn in cov_files:

        # Open file
        fn = os.path.join(sample_dir, cov_fn)
        with open(fn, "r") as f:

            # Iterate over lines
            for line in f:

                # Parse line
                _, _, pos, cov = line.strip().split("\t")

                # Store results
                dt["pos"].append(int(pos))  # note data type coercion
                dt["cov"].append(int(cov))  # ""
    
    # Convert to data frame
    df = pd.DataFrame(dt)
    
    # Get total coverage at each position
    df = (df
          .groupby("pos")
          .sum()
          .reset_index()
         )
    
    return df


def load_sample_list(fn):
    """
    Load Sample List into Python
    as a dictionary which maps sample
    names to barcodes
    
    Parameters
        fn : str
            Path to SampleList file.
    
    Returns
        dt : dict
            Maps sample name to barcode.
    
    """
    
    with open(fn, "r") as f:
        dt = {}
        for line in f:
            if not line.isspace():
                barcode, sample = line.rstrip().split("\t")
                dt[sample] = int(barcode)
    
    return dt


def get_contig_lengths(cov):
    """
    Given a coverage profile, return an array
    giving the list of all the contigs
    
    In this context a contig is defined as a contiguous
    region of sequence with a coverage > 0X.
    
    Parameters
        cov : array, int, shape (n_sites)
            Coverage along the genome; number of reads
            covering each consecutive nucleotide.
            
    Returns
        ll : ndarray, int, shape(n_contigs)
            An array containing the length of each
            of the contigs
            
    """
    
    contigs = []
    l = 0
    for covered in cov > 0:
        if covered:
            l += 1
        elif l > 0:
            contigs.append(l)
            l = 0
    
    if l > 0: 
        contigs.append(l)  # append last segment

    return np.array(contigs)


def calc_contig_n50(cov):
    """
    Compute the N50 contig length
    
    Parameters
        cov : array, int, shape (n_sites)
            Coverage along the genome; number of reads
            covering each consecutive nucleotide.
            
    Returns
        _ : int
            The N50 contig length; 50% of the total
            contig length is contained in contigs
            of this length or smaller.

    """
    
    contigs = get_contig_lengths(cov)
    
    if len(contigs) > 0:
        contigs.sort()
        total_contig = contigs.sum()
        cum_contig_frac = (contigs / total_contig).cumsum()
        ix = np.argmax(cum_contig_frac > 0.5)
        
    else:
        return np.nan
    
    return contigs[ix]


def calc_gisaid_stats(df, depth_threshold=10):
    """
    Calculate bioinformatics statistics required for 
    GISAD submission
    
    Parameters
        df : DataFrame, shape (2, n_sites)
            Pandas DataFrame containing two columns,
            one for genomic position "pos" and one for
            coverage "cov". Rows are consecutive
            sites in the genome.
        depth_threshold : int, optional
            Used to calculate breadth.
            
    Returns
        dt : dict
            Dictionary containing all GISAID QC
            statistics, excluding 'Ns per 100 kbp'
            
    """
    
    # Reference genome length
    genome_length = df["pos"].max() - df["pos"].min()
    
    # Consensus genome length
    covered = df["cov"] > 0
    consensus_genome_length = df[covered]["pos"].max() - df[covered]["pos"].min() + 1
    
    # Total base pairs sequenced
    # NB: this is post-mapping and depth truncation
    number_base_pairs = df["cov"].sum()
    
    # Sequence depth (average)
    sequencing_depth_avg = number_base_pairs / genome_length
    
    # Coverage depth (average)
    coverage_depth = number_base_pairs / consensus_genome_length
    
    # Assembly coverage breadth (as percentage)
    coverage_breadth = 100 * (df["cov"] >= depth_threshold).sum() / genome_length
    
    # Mean contig length
    contigs = get_contig_lengths(df["cov"])
    mean_contig_length = contigs.mean()
    
    # N50 contig length
    n50_contig_length = calc_contig_n50(df["cov"])
    
    
    # Package key statistics into a dictionary
    dt = {
        "sequencing_depth_avg": sequencing_depth_avg,
        "coverage_breadth": coverage_breadth,
        "assembly_coverage_depth": coverage_depth,
        "number_base_pairs": number_base_pairs,
        "consensus_genome_length": consensus_genome_length,
        "ref_genome_length": genome_length,
        "mean_contig_length": mean_contig_length,
        "N50": n50_contig_length
    }
    
    return dt


def calc_ns_per_100kbp(consensus_fasta, verbose=True):
    """
    Calculate the number of Ns per 100kbp from
    a consensus FASTA
    
    Parameters
        consensus_fasta : str
            Path to the consensus fasta file.
        verbose : bool
            Switch for print diagnostics.
    
    Returns
        ns_per_100kbp : float
            The average number of Ns per 100kbp.
    
    """
    
    # Open consensus fasta
    with open(consensus_fasta, "r") as f:
    
        # Parse header
        header = f.readline().strip()

        # Parse sequence
        seq = ""
        for line in f:
            seq += line.strip()

        # Calculate statistics
        seq_length = len(seq)

    # Number of Ns
    num_ns = sum([1 for c in seq if c == "N"])

    # Ns per 100 kbp
    ns_per_100kbp = num_ns / seq_length * 100 * 10 **3

    if verbose:
        print("FASTA name: %s" % consensus_fasta)
        print("Header: %s" % header)
        print("Sequence length (bp): %d" % seq_length)
        print("Number of Ns: %d" % num_ns)
        print("Ns per 100kbp: %.02f" % ns_per_100kbp)
        
    return ns_per_100kbp


def calc_avg_seq_depth(fastq_fn, genome_length):
    """
    Calculate average sequencing depth, given
    a fastq file named `fastq_fn` and a
    `genome_length`
    
    Parameters
        fastq_fn : str
            Path to fastq file.
        genome_length : int
            Length of the genome from which the
            reads in the fastq were generated.
            
    Returns
        _ : int
            Average sequencing depth given all
            the data within the fastq.
    
    """
    
    # Code golf style
    cmd = "cat %s | paste - - - - | cut -f2 | tr -d '\n' | wc -c" % fastq_fn
    total_bp = int(subprocess.check_output(cmd, shell=True))
    return total_bp / genome_length

def calc_fastq_total_reads(fastq_fn):
    """
    Calculate number of reads in a given fastq file 
    named `fastq_fn` 
    
    Parameters
        fastq_fn : str
            Path to fastq file.
            
    Returns
        _ : int
            Number of reads associated with that barcode
            from the fastq.
    
    """
        
    c = 0
    with open(fastq_fn, "r") as f: # Open fastq_fn
        for line in f:
            c += 0.25 #There are four lines per read
    return int(c)

