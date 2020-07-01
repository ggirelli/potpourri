#!/usr/bin/env python3

# -----------------------------------------------------------------------------
#
# Author: ...
# Email: ...
# Version: X.X.X
# Date: YYYYMMDD
# Project: ...
# Description: ...
#
# -----------------------------------------------------------------------------


import argparse
from Bio import SeqIO
import gzip
import os
import re
from tqdm import tqdm
import sys


version = "0.0.1"
parser = argparse.ArgumentParser(description='''
Description...
''', formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument(
    'input', type=str,
    help='Path to fasta or fastq file.')
parser.add_argument(
    'trim_len', type=int,
    help='Number of characters to be trimmed from the beginning of each read')

parser.add_argument(
    '-p', type=str, metavar='pattern',
    help="""Optional pattern at the beginning of untrimmed part.""")
parser.add_argument(
    '-a', type=int, metavar='allowance',
    help="""Number of characters allowed as shift. Default: 0""", default=0)
parser.add_argument(
    '-o', type=str, metavar='output',
    help="""Path to output file. Defaults to input with .trimmed suffix.""")

parser.add_argument(
    '--no-skip', action='store_const', dest='do_skip',
    const=False, default=True,
    help='Use this option to avoid skipping reads with large shifts.')
parser.add_argument(
    '--overwrite', action='store_const', dest='overwrite',
    const=True, default=False,
    help='Overwrite output if found instead of stopping.')

parser.add_argument('--version', action='version',
                    version='%s v%s' % (sys.argv[0], version,))

args = parser.parse_args()

assert os.path.isfile(args.input)

input_is_gzipped = False
possible_input_types = {".fastq": "fastq", ".fasta": "fasta"}

input_basename, input_ext = os.path.splitext(args.input)
if ".gz" == input_ext:
    input_is_gzipped = True
    input_ext = os.path.splitext(input_basename)[-1]
input_type = possible_input_types[input_ext]

if args.o is not None:
    output_path = args.o
else:
    if input_is_gzipped:
        input_basename = os.path.splitext(input_basename)[0]
    output_path = f"{input_basename}.trimmed.{input_type}"
    if input_is_gzipped:
        output_path += ".gz"

print(f"Input: {args.input}")
print(f"Type: {input_type}")
print(f"Gzipped: {input_is_gzipped}")
print(f"Output: {output_path}")

if args.overwrite:
    if os.path.isfile(output_path):
        print("Warning: overwriting output.")
if not args:
    assert not os.path.isfile(output_path)

if input_is_gzipped:
    IH = gzip.open(args.input, "rt")
    OH = gzip.open(output_path, "wt")
else:
    IH = open(args.input, "r")
    OH = open(output_path, "w+")


def trim_record(record, trim_len, pattern=None, allowance=0, skip=True) -> int:
    if pattern is not None:
        match = re.search(pattern, str(record.seq))
        if match is not None:
            match_start = match.start()
            shift = abs(match_start-trim_len)
            if shift <= allowance:
                return record[match_start:]
            else:
                if skip:
                    print(f"{shift} characters shift in {record.id}, skipped.")
                    return None
                # else:
                #     print(f"default trim for {record.id}, shift: {shift}")
    return record[trim_len:]


for record in tqdm(SeqIO.parse(IH, input_type)):
    trimmed_record = trim_record(record, args.trim_len,
                                 args.p, args.a, args.do_skip)
    if trimmed_record is not None:
        OH.write(trimmed_record.format(input_type))

IH.close()
OH.close()
print("done")
