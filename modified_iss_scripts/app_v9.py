#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import bam
from iss import util
from iss import download
from iss import abundance
from iss import generator
from iss.version import __version__

from Bio import SeqIO
from joblib import Parallel, delayed, load, dump

import gc
import os
import sys
import pickle
import random
import logging
import argparse
import numpy as np
import glob
from multiprocessing import Pool
from decimal import Decimal
from collections import OrderedDict
import csv
import io
import subprocess
import re
import os
import logging


def read_fragment_lengths(file_path):
    fragment_lengths = {}
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter=';')
        for row in reader:
            genome_id, length = row
            fragment_lengths[genome_id] = int(length)
    return fragment_lengths


def process_genome_chunk(args_tuple):
    # Unpack all 8 arguments, including buffer_size
    fasta_file, csv_file, err_mod, output, seed, gc_bias, cpu_number, buffer_size = args_tuple
    output_files_for_core = []

    total_reads_per_genome = {}
    with open(csv_file, 'r') as file:
        reader = csv.reader(file, delimiter=';')
        for row in reader:
            genome_id, n_pairs = row
            total_reads_per_genome[genome_id] = int(n_pairs)

    # Generate unique filenames for this core's output
    temp_file_name_R1 = f"{output}_core_{cpu_number}_R1.fastq"
    temp_file_name_R2 = f"{output}_core_{cpu_number}_R2.fastq"

    read_buffer = []
    buffer_size = 10000

    for record in SeqIO.parse(fasta_file, 'fasta'):
        n_pairs = total_reads_per_genome.get(record.id, 0)
        if n_pairs < 1:
            continue

        mode = "memmap" if sys.getsizeof(str(record.seq)) >= 2**31 - 1 else "default"
        if mode == "memmap":
            record_mmap = f"{output}.memmap"
            if os.path.exists(record_mmap):
                os.unlink(record_mmap)
            util.dump(record, record_mmap)
            record = load(record_mmap)

        temp_reads = generator.reads(record, err_mod, n_pairs, cpu_number, seed, gc_bias, mode)
        read_buffer.extend(temp_reads)

        if len(read_buffer) >= buffer_size:
            generator.to_fastq(read_buffer, temp_file_name_R1, temp_file_name_R2)
            read_buffer.clear()

    if read_buffer:
        generator.to_fastq(read_buffer, temp_file_name_R1, temp_file_name_R2)

    return temp_file_name_R1, temp_file_name_R2



def concatenate_files(file_list, output_file):
    with open(output_file, 'wb') as outfile:
        for filename in file_list:
            with open(filename, 'rb') as readfile:
                outfile.write(readfile.read())

def generate_reads(args):
    """Main function for the `iss generate` submodule.

    This submodule generates reads from an ErrorModel and writes them to
    args.output + _R(1|2).fastq.

    Args:
        args (object): The command-line arguments from argparse.
    """
    logger = logging.getLogger(__name__)
    logger.debug('iss version %s' % __version__)
    logger.debug('Using verbose logger')

    try:
        logger.info('Starting iss generate')
        logger.info('Using %s ErrorModel' % args.mode)
        if args.seed:
            logger.info('Setting random seed to %i' % args.seed)
            random.seed(args.seed)
            np.random.seed(args.seed)
        if args.mode == 'kde':
            from iss.error_models import kde
            if args.model is None:
                logger.error('--model is required in --mode kde')
                sys.exit(1)
            elif args.model.lower() == 'hiseq':
                npz = os.path.join(
                    os.path.dirname(__file__),
                    'profiles/HiSeq')
            elif args.model.lower() == 'novaseq':
                npz = os.path.join(
                    os.path.dirname(__file__),
                    'profiles/NovaSeq')
            elif args.model.lower() == 'miseq':
                npz = os.path.join(
                    os.path.dirname(__file__),
                    'profiles/MiSeq')
            else:
                npz = args.model
            err_mod = kde.KDErrorModel(npz)
        elif args.mode == 'basic':
            if args.model is not None:
                logger.warning('--model %s will be ignored in --mode %s' % (
                    args.model, args.mode))
            from iss.error_models import basic
            err_mod = basic.BasicErrorModel()
        elif args.mode == 'perfect':
            if args.model is not None:
                logger.warning('--model %s will be ignored in --mode %s' % (
                    args.model, args.mode))
            from iss.error_models import perfect
            err_mod = perfect.PerfectErrorModel()
    except ImportError as e:
        logger.error('Failed to import ErrorModel module: %s' % e)
        sys.exit(1)

    try:
        # Read and process the fasta file list
        with open(args.fasta_list, 'r') as fasta_file:
            fasta_files = [line.strip() for line in fasta_file if line.strip()]

        # Read and process the CSV file list
        with open(args.csv_reads_list, 'r') as csv_file:
            csv_files = [line.strip() for line in csv_file if line.strip()]

        # Create a mapping of indices to filenames for fasta and CSV files
        fasta_mapping = {int(re.search(r'\d+', os.path.basename(f)).group()): f for f in fasta_files}
        csv_mapping = {int(re.search(r'\d+', os.path.basename(f)).group()): f for f in csv_files}

        # Check for equal number of files and matching indices
        if len(fasta_mapping) != len(csv_mapping) or not all(idx in csv_mapping for idx in fasta_mapping):
            logger.error('Mismatch in the number of fasta and CSV files or their indices')
            sys.exit(1)

        cpus = args.cpus
        logger.info('Using %s cpus for read generation' % cpus)
        buffer_size = 1000     
        process_args = []
        for idx in range(1, cpus+1):
            if idx in fasta_mapping and idx in csv_mapping:
                process_args.append((fasta_mapping[idx], csv_mapping[idx], err_mod, args.output, args.seed, args.gc_bias, idx-1,buffer_size))
            else:
                logger.error(f'File pair for index {idx} not found')
                sys.exit(1)

        with Pool(processes=cpus) as pool:
            core_output_files = pool.map(process_genome_chunk, process_args)

    except KeyboardInterrupt as e:
        logger.error('iss generate interrupted: %s' % e)
        # Cleanup all matching temporary files
        for temp_file in glob.glob(f"{args.output}_core_*"):
            os.remove(temp_file)
        sys.exit(1)
    else:
        # Collect R1 and R2 files using glob
        r1_files = sorted(glob.glob(f"{args.output}_core_*_R1.fastq"))
        r2_files = sorted(glob.glob(f"{args.output}_core_*_R2.fastq"))

        # Concatenate R1 files
        concatenate_files(r1_files, f"{args.output}_R1.fastq")

        # Concatenate R2 files
        concatenate_files(r2_files, f"{args.output}_R2.fastq")

        # Cleanup the temporary R1 and R2 files
        for temp_file in r1_files + r2_files:
            os.remove(temp_file)

        if args.compress:
            util.compress(args.output + '_R1.fastq')
            util.compress(args.output + '_R2.fastq')

        logger.info('Read generation complete')

def model_from_bam(args):
    """Main function for the `iss model` submodule

    This submodule write all variables necessary for building an ErrorModel
    to args.output + .npz

    Args:
        args (object): the command-line arguments from argparse
    """
    logger = logging.getLogger(__name__)
    logger.debug('iss version %s' % __version__)
    logger.debug('Using verbose logger')

    try:  # try to import bam module and write model data to file
        logger.info('Starting iss model')
        from iss import bam
    except ImportError as e:
        logger.error('Failed to import bam module: %s' % e)
        sys.exit(1)
    else:
        logger.info('Using KDE ErrorModel')
        bam.to_model(args.bam, args.output)
        logger.info('Model generation complete')


def main():
    parser = argparse.ArgumentParser(
        prog='InSilicoSeq',
        usage='iss [subcommand] [options]',
        description='InSilicoSeq: A sequencing simulator'
    )
    parser.add_argument(
        '-v',
        '--version',
        action='store_true',
        default=False,
        help='print software version and exit'
    )
    subparsers = parser.add_subparsers(
        title='available subcommands',
        metavar=''
    )

    parser_mod = subparsers.add_parser(
        'model',
        prog='iss model',
        description='generate an error model from a bam file',
        help='generate an error model from a bam file'
    )
    parser_gen = subparsers.add_parser(
        'generate',
        prog='iss generate',
        description='simulate reads from an error model',
        help='simulate reads from an error model'
    )

    # arguments form the read generator module
    param_logging = parser_gen.add_mutually_exclusive_group()
    # --genomes and --ncbi should not be exclusive anymore
    # input_genomes = parser_gen.add_mutually_exclusive_group()
    input_abundance = parser_gen.add_mutually_exclusive_group()
    param_logging.add_argument(
        '--quiet',
        '-q',
        action='store_true',
        default=False,
        help='Disable info logging. (default: %(default)s).'
    )
    param_logging.add_argument(
        '--debug',
        '-d',
        action='store_true',
        default=False,
        help='Enable debug logging. (default: %(default)s).'
    )
    parser_gen.add_argument(
        '--seed',
        type=int,
        metavar='<int>',
        help='Seed all the random number generators',
        default=None
    )
    parser_gen.add_argument(
        '--cpus',
        '-p',
        default=2,
        type=int,
        metavar='<int>',
        help='number of cpus to use. (default: %(default)s).'
    )
    parser_gen.add_argument(
        '--genomes',
        '-g',
        metavar='<genomes.fasta>',
        nargs="+",
        help='Input genome(s) from where the reads will originate'
    )
    parser_gen.add_argument(
        '--draft',
        metavar='<draft.fasta>',
        nargs="+",
        help='Input draft genome(s) from where the reads will originate'
    )
    parser_gen.add_argument(
        '--n_genomes',
        '-u',
        type=int,
        metavar='<int>',
        help='How many genomes will be used for the simulation. is set with \
            --genomes/-g or/and --draft to take random genomes from the \
            input multifasta'
    )
    parser_gen.add_argument(
        '--ncbi',
        '-k',
        choices=['bacteria', 'viruses', 'archaea'],
        action='append',
        nargs='*',
        metavar='<str>',
        help='Download input genomes from NCBI. Requires --n_genomes/-u\
            option. Can be bacteria, viruses, archaea or a combination of the\
            three (space-separated)'
    )
    parser_gen.add_argument(
        '--n_genomes_ncbi',
        '-U',
        type=int,
        action='append',
        metavar='<int>',
        nargs='*',
        help='How many genomes will be downloaded from NCBI. Required if\
            --ncbi/-k is set. If more than one kingdom is set with --ncbi,\
            multiple values are necessary (space-separated).'
    )
    input_abundance.add_argument(
        '--abundance',
        '-a',
        choices=['uniform', 'halfnormal',
                 'exponential', 'lognormal', 'zero_inflated_lognormal'],
        metavar='<str>',
        default='lognormal',
        help='abundance distribution (default: %(default)s). Can be uniform,\
            halfnormal, exponential, lognormal or zero-inflated-lognormal.'
    )
    input_abundance.add_argument(
        '--abundance_file',
        '-b',
        metavar='<abundance.txt>',
        help='abundance file for coverage calculations (default: %(default)s).'
    )
    input_abundance.add_argument(
        '--coverage',
        '-C',
        choices=['uniform', 'halfnormal',
                 'exponential', 'lognormal', 'zero_inflated_lognormal'],
        metavar='<str>',
        help='coverage distribution. Can be uniform,\
            halfnormal, exponential, lognormal or zero-inflated-lognormal.'
    )
    input_abundance.add_argument(
        '--coverage_file',
        '-D',
        metavar='<coverage.txt>',
        help='file containing coverage information (default: %(default)s).'
    )
    parser_gen.add_argument(
        '--n_reads',
        '-n',
        metavar='<int>',
        default='1000000',
        help='Number of reads to generate (default: %(default)s). Allows \
        suffixes k, K, m, M, g and G (ex 0.5M for 500000).'
    )
    parser_gen.add_argument(
        '--mode',
        '-e',
        metavar='<str>',
        choices=['kde', 'basic', 'perfect'],
        default='kde',
        help='Error model. If not specified, using kernel density estimation \
        (default: %(default)s). Can be kde, basic or perfect'
    )
    parser_gen.add_argument(
        '--model',
        '-m',
        metavar='<npz>',
        default=None,
        help='Error model file. (default: %(default)s). Use HiSeq, NovaSeq or \
        MiSeq for a pre-computed error model provided with the software, or a \
        file generated with iss model. If you do not wish to use a model, use \
        --mode basic or --mode perfect. The name of the built-in models are  \
        case insensitive.'
    )

    parser_gen.add_argument(
        '--gc_bias',
        '-c',
        action='store_true',
        default=False,
        help='If set, may fail to sequence reads with abnormal GC content. \
        (default: %(default)s)'
    )
    parser_gen.add_argument(
        '--compress',
        '-z',
        action='store_true',
        default=False,
        help='Compress the output in gzip format (default: %(default)s).'
    )
    parser_gen.add_argument(
        '--output',
        '-o',
        metavar='<fastq>',
        help='Output file path and prefix (Required)',
        required=True
    )
    parser_gen._optionals.title = 'arguments'
    parser_gen.set_defaults(func=generate_reads)

    # arguments for the error model module
    parser_mod.add_argument(
        '--quiet',
        '-q',
        action='store_true',
        default=False,
        help='Disable info logging. (default: %(default)s).'
    )
    parser_mod.add_argument(
        '--debug',
        '-d',
        action='store_true',
        default=False,
        help='Enable debug logging. (default: %(default)s).'
    )
    parser_mod.add_argument(
        '--bam',
        '-b',
        metavar='<bam>',
        help='aligned reads from which the model will be inferred (Required)',
        required=True
    )
    parser_gen.add_argument(
        '--fasta_list',
        required=True,
        help='File containing the list of fasta files'
    )
    parser_gen.add_argument(
        '--csv_reads_list',
        required=True,
        help='File containing the list of CSV files'
    )  

    parser_mod.add_argument(
        '--output',
        '-o',
        metavar='<npz>',
        help='Output file path and prefix (Required)',
        required=True
    )
    parser_mod._optionals.title = 'arguments'
    parser_mod.set_defaults(func=model_from_bam)
    args = parser.parse_args()

    # set logger and display version if args.version
    try:
        if args.version:
            print('iss version %s' % __version__)
            sys.exit(0)
        elif args.quiet:
            logging.basicConfig(level=logging.ERROR)
        elif args.debug:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)

        args.func(args)
        logging.shutdown()
    except AttributeError as e:
        logger = logging.getLogger(__name__)
        logger.debug(e)
        parser.print_help()
        # raise  # extra traceback to uncomment if all hell breaks lose

