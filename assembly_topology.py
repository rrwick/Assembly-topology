#!/usr/bin/env python3
"""
Count types of contig-end alignments from a PAF file to help determine the topology of contigs
in a long-read assembly.

Inputs (positional):
1) Assembly FASTA file.
2) PAF-format read alignments (from minimap2 with -c).

Output to stdout:
- TSV with one row per assembly sequence and these columns:
    name, circular_reads, hairpin_reads, terminating_reads, clipping_reads, ambiguous_reads.

Definitions (within a tolerance of bases, set by --tolerance):
- circular: The read crosses a contig boundary on the same strand — two alignments to the same
    target and same strand, one near the left end and the other near the right end of the
    reference, and the two query segments are adjacent on the read.
- hairpin: The read folds back at the same contig end — two alignments to the same target but
    opposite strands, both near the same contig end (left or right), and their query segments
    are adjacent.
- terminating: A single alignment ends at a contig end and the read also ends.
- clipping: Any read that hits a contig end but is neither circular, hairpin nor terminating.
- ambiguous: Any read that qualifies for more than one of the above categories. For example, 
    when a linear sequenceas has a terminal inverted repeat and the read is contained within this
    repeat, its alignments can be compatible with both circular and hairpin structures. Also,
    any read with excessively an short alignment (as set by --min_align_len) will go in this
    category.

Usage example:
    minimap2 -c -t 16 assembly.fasta reads.fastq.gz > alignments.paf
    ./assembly_topology.py assembly.fasta alignments.paf
"""

import argparse
import collections
import gzip
import sys


def get_arguments():
    parser = argparse.ArgumentParser(description='Assembly topology from long-read alignments',
                                     add_help=False)

    required_args = parser.add_argument_group('Positional arguments')
    required_args.add_argument('assembly_fasta', type=str,
                               help='Assembly sequence in FASTA format')
    required_args.add_argument('alignments_paf', type=str,
                               help='Long-read alignments in PAF format')

    setting_args = parser.add_argument_group('Settings')
    setting_args.add_argument('--tolerance', type=int, default=10,
                              help='Allowed gap size when comparing alignment positions '
                                   '(default: 10)')
    setting_args.add_argument('--min_align_len', type=int, default=500,
                              help='Alignments shorter than this lead to an ambiguous read '
                                   '(default: 500)')

    other_args = parser.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help',
                            help='Show this help message and exit')

    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    target_names = get_fasta_names(args.assembly_fasta)
    alignments = load_alignments(args.alignments_paf)
    print_header()
    counts = count_sequence_end_reads(target_names, alignments, args.tolerance, args.min_align_len)
    for name in target_names:
        circular, hairpin, terminating, clipped, ambiguous = counts[name]
        print(f'{name}\t{circular}\t{hairpin}\t{terminating}\t{clipped}\t{ambiguous}')


def count_sequence_end_reads(target_names, alignments, tolerance, min_align_len):
    counts = {}
    for name in target_names:
        counts[name] = [0, 0, 0, 0, 0]

    alignments_by_read = collections.defaultdict(list)
    for a in alignments:
        alignments_by_read[a.query_name].append(a)

    for read_alignments in alignments_by_read.values():
        if alignments_have_multiple_targets(read_alignments):
            continue
        if not alignments_hit_sequence_end(read_alignments, tolerance):
            continue
        ref_name = read_alignments[0].target_name
        circ = alignments_are_circular(read_alignments, tolerance)
        hairpin = alignments_are_hairpin(read_alignments, tolerance)
        term = alignments_are_terminating(read_alignments, tolerance)
        too_short = alignments_are_too_short(read_alignments, min_align_len)
        if sum([circ, hairpin, term]) > 1 or too_short:  # ambiguous
            counts[ref_name][4] += 1
        elif circ:
            counts[ref_name][0] += 1
        elif hairpin:
            counts[ref_name][1] += 1
        elif term:
            counts[ref_name][2] += 1
        else:  # clipped
            counts[ref_name][3] += 1

    return counts


def alignments_have_multiple_targets(alignments):
    return len({a.target_name for a in alignments}) > 1


def alignments_hit_sequence_end(alignments, tolerance):
    # Return True if any alignment touches either end of the sequence (within tolerance).
    for a in alignments:
        if hits_left_end(a, tolerance) or hits_right_end(a, tolerance):
            return True
    return False


def alignments_are_circular(alignments, tolerance):
    # Return True if any pair of alignments indicates a circular join: one hits the left end and
    # the other the right end of the same target, on the same strand, and their query segments are
    # adjacent (within tolerance).
    n = len(alignments)
    for i in range(n):
        ai = alignments[i]
        for j in range(i + 1, n):
            aj = alignments[j]
            if ai.target_name != aj.target_name:
                continue
            if ai.strand != aj.strand:
                continue
            opposite_ends = (hits_left_end(ai, tolerance) and hits_right_end(aj, tolerance)) or \
                            (hits_left_end(aj, tolerance) and hits_right_end(ai, tolerance))
            if not opposite_ends:
                continue
            if segments_adjacent_on_read(ai, aj, tolerance):
                return True
    return False


def alignments_are_hairpin(alignments, tolerance):
    # Return True if any pair of alignments indicates a hairpin: both hit the same end of the same
    # target, on opposite strands, and their query segments are adjacent (within tolerance).
    n = len(alignments)
    for i in range(n):
        ai = alignments[i]
        for j in range(i + 1, n):
            aj = alignments[j]
            if ai.target_name != aj.target_name:
                continue
            if ai.strand == aj.strand:
                continue
            same_left = hits_left_end(ai, tolerance) and hits_left_end(aj, tolerance)
            same_right = hits_right_end(ai, tolerance) and hits_right_end(aj, tolerance)
            if not (same_left or same_right):
                continue
            if segments_adjacent_on_read(ai, aj, tolerance):
                return True
    return False


def alignments_are_terminating(alignments, tolerance):
    # Return True if any single alignment ends at the sequence end and the read also ends (within
    # tolerance) on the corresponding read side based on strand and end.
    for a in alignments:
        if hits_left_end(a, tolerance):
            if (a.strand == '+' and read_hits_start(a, tolerance)) or (a.strand == '-' and read_hits_end(a, tolerance)):
                return True
        if hits_right_end(a, tolerance):
            if (a.strand == '+' and read_hits_end(a, tolerance)) or (a.strand == '-' and read_hits_start(a, tolerance)):
                return True
    return False


def alignments_are_too_short(alignments, min_align_len):
    for a in alignments:
        if a.query_end - a.query_start < min_align_len:
            return True
    return False


def hits_left_end(a: "Alignment", tolerance) -> bool:
    return a.target_start <= tolerance


def hits_right_end(a: "Alignment", tolerance) -> bool:
    return a.target_end >= (a.target_length - tolerance)


def read_hits_start(a: "Alignment", tolerance) -> bool:
    return a.query_start <= tolerance


def read_hits_end(a: "Alignment", tolerance) -> bool:
    return a.query_end >= (a.query_length - tolerance)


def segments_adjacent_on_read(a: "Alignment", b: "Alignment", tolerance) -> bool:
    # Are two alignments adjacent along the read within tolerance?
    return (abs(a.query_end - b.query_start) <= tolerance or
            abs(b.query_end - a.query_start) <= tolerance)


def print_header():
    print('name\tcircular_reads\thairpin_reads\tterminating_reads\tclipping_reads\tambiguous_reads')


def get_fasta_names(filename):
    names = []
    with get_open_func(filename)(filename, 'rt') as fasta_file:
        for line in fasta_file:
            if line.startswith('>'):
                names.append(line[1:].strip().split()[0])
    return names


def load_alignments(filename):
    alignments = []
    with get_open_func(filename)(filename, 'rt') as paf_file:
        for line in paf_file:
            alignments.append(Alignment(line))
    return alignments


def get_compression_type(filename):
    # Attempts to guess the compression (if any) on a file using the first few bytes.
    # https://stackoverflow.com/questions/13044562
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)
    unknown_file = open(str(filename), 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('\nError: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('\nError: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


class Alignment(object):

    def __init__(self, paf_line):
        line_parts = paf_line.strip().split('\t')
        if len(line_parts) < 11:
            sys.exit('Error: alignment file does not seem to be in PAF format')

        self.query_name = line_parts[0]
        self.query_length = int(line_parts[1])
        self.query_start = int(line_parts[2])
        self.query_end = int(line_parts[3])
        self.strand = line_parts[4]

        self.target_name = line_parts[5]
        self.target_length = int(line_parts[6])
        self.target_start = int(line_parts[7])
        self.target_end = int(line_parts[8])

        self.cigar = None
        for part in line_parts:
            if part.startswith('cg:Z:'):
                self.cigar = part[5:]

        if self.cigar is None:            
            sys.exit('Error: no CIGAR string found in PAF line - did you use the -c option?')


if __name__ == '__main__':
    main()
