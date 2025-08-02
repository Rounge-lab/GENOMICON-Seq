#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss.util import load, rev_comp

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqRecord import SeqRecord
from shutil import copyfileobj

import os
import sys
import random
import logging
import numpy as np
import logging

logging.basicConfig(level=logging.DEBUG)


def reads(record, ErrorModel, n_pairs, cpu_number, seed, gc_bias=False, mode="default"):
    logger = logging.getLogger(__name__)

    # Load the record from disk if mode is memmap
    if mode == "memmap":
        record_mmap = load(record)
        record = record_mmap

    if seed is not None:
        random.seed(seed + cpu_number)
        np.random.seed(seed + cpu_number)

    read_tuple_list = []
    i = 0
    while i < n_pairs:
        try:
            forward, reverse = simulate_read(record, ErrorModel, i, cpu_number)
            if gc_bias:
                stitched_seq = forward.seq + reverse.seq
                gc_content = gc_fraction(stitched_seq)
                if not (40 < gc_content < 60 or np.random.rand() < 0.90):
                    continue
            read_tuple_list.append((forward, reverse))
            i += 1
        except AssertionError as e:
            logger.warning(f'{record.id} shorter than read length for this ErrorModel: {e}')
            logger.warning(f'Skipping {record.id}. You will have fewer reads than specified')
            break

    return read_tuple_list


def simulate_read(record, ErrorModel, i, cpu_number):
    """From a read pair from one genome (or sequence) according to an
    ErrorModel

    Each read is a SeqRecord object
    returns a tuple containing the forward and reverse read.

    Args:
        record (SeqRecord): sequence or genome of reference
        ErrorModel (ErrorModel): an ErrorModel class
        i (int): a number identifying the read
        cpu_number (int): cpu number. Is added to the read id.

    Returns:
        tuple: tuple containg a forward read and a reverse read
    """
    logger = logging.getLogger(__name__)
    sequence = record.seq
    header = record.id

    read_length = ErrorModel.read_length
    insert_size = ErrorModel.random_insert_size()

    # Modification starts here
    # Always start sequencing from the beginning of the fragment for R1
    forward_start = 0
    forward_end = min(read_length, len(record.seq))  # Ensure we don't go beyond the sequence length

    # Always start sequencing from the end of the fragment for R2
    reverse_end = len(record.seq)
    reverse_start = max(0, reverse_end - read_length)  # Ensure we don't go beyond the sequence length
    # Modification ends here

    bounds_forward = (forward_start, forward_end)
    bounds_reverse = (reverse_start, reverse_end)

    # create a perfect read for R1
    forward = SeqRecord(
        Seq(str(sequence[forward_start:forward_end])),
        id='%s_%s_%s/1' % (header, i, cpu_number),
        description=''
    )
    # add the indels, the qual scores and modify the record accordingly
    forward.seq = ErrorModel.introduce_indels(
        forward, 'forward', sequence, bounds_forward)
    forward = ErrorModel.introduce_error_scores(forward, 'forward')
    forward.seq = ErrorModel.mut_sequence(forward, 'forward')

    # create a perfect read f   or R2
    reverse = SeqRecord(
        Seq(rev_comp(str(sequence[reverse_start:reverse_end]))),
        id='%s_%s_%s/2' % (header, i, cpu_number),
        description=''
    )
    # add the indels, the qual scores and modify the record accordingly
    reverse.seq = ErrorModel.introduce_indels(
        reverse, 'reverse', sequence, bounds_reverse)
    reverse = ErrorModel.introduce_error_scores(reverse, 'reverse')
    reverse.seq = ErrorModel.mut_sequence(reverse, 'reverse')

    return (forward, reverse)

def to_fastq(read_tuple_list, output_R1, output_R2):
    """Write reads to a fastq file

    Take a list containing read pairs (tuples) and write them
    in two fastq files: output_R1.fastq and output_R2.fastq

    Args:
        read_tuple_list (list): a list of tuples, each containing a pair of reads
        output_R1 (str): the output file for R1 reads
        output_R2 (str): the output file for R2 reads
    """
    logger = logging.getLogger(__name__)

    # Open the files in append mode since this function can be called multiple times
    with open(output_R1, 'a') as f_r1, open(output_R2, 'a') as f_r2:
        for forward_read, reverse_read in read_tuple_list:
            SeqIO.write(forward_read, f_r1, 'fastq-sanger')
            SeqIO.write(reverse_read, f_r2, 'fastq-sanger')


