# Copyright 2013 Graham McVicker and Bryce van de Geijn
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Modification - OS:
#  - Either remove chromosomes with different Null hypothesis (chrX, chrY and chrM) or fit coefficients only for specified chromosomes
#  - header was removed twice in the original code (in the function for reading input files and in main())


import sys
import os
import math
import gzip
import argparse
from scipy.optimize import *
from scipy.special import betaln
import scipy.stats
import numpy as np
import util
import datetime

def print_with_time(message):
    """Prints a message with the current timestamp."""
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{current_time}] {message}")


def parse_options():
    parser = argparse.ArgumentParser(description="This script estimates the "
                                     "overdispersion parameter for the allele-specific "
                                     "(Beta-Binomial) half of the combined halotype test."
                                     " A single overdispersion parameter is estimated for "
                                     "each individual (across sites), under the assumption "
                                     "that all sites come from the null hypothesis (no "
                                     "genetic association)");


    parser.add_argument("--read_error_rate", "-e", action='store', dest='read_error_rate',
                        help="sequence read error rate (default=0.005)",
                        type=float, default=0.005)

    parser.add_argument('--remove_chromosomes', "-r", help="Remove some chromosomes? Can't be used together with --keep-chromosomes flag",
                        action="store_true", dest="remove_chromosomes", default=False)

    parser.add_argument('--keep_chromosomes', "-k", help="Keep only specified chromosomes? Can't be used together with --remove-chromosomes flag",
                        action="store_true", dest="keep_chromosomes", default=False)

    parser.add_argument('--chrom_list', nargs='*', dest='chrom_list',
                        help="List of chromosomes to remove (with -r/--remove_chromosomes flag) or to keep (with -k/--keep_chromosomes flag). Space separated lise, e.g. chrX chrM",
                        default=None)

    parser.add_argument('--aggregate_per_region', "-a", help="Aggregate all counts per region (default=False). If used, heterozygous probabilities are assumed to be 1 for all heterozygous sites",
                        action="store_true", dest="aggregate_per_region", default=False)

    parser.add_argument("infile_list",
                        help="Path to file containing list of CHT input "
                        "files (one for each individual)",
                        action='store', default=None)

    parser.add_argument("out_file",
                        help="File to write overdispersion parameter estimates to",
                        action='store', default=None)

    return parser.parse_args()




def open_input_files(in_filename):
    if not os.path.exists(in_filename) or not os.path.isfile(in_filename):
        sys.stderr.write("input file %s does not exist or is not a regular file\n" %
                         in_filename)
        exit(2)

    # read file that contains list of input files
    in_file = open(in_filename, "rt")

    infiles = []
    for line in in_file:
        # open each input file and read first line
        filename = line.rstrip()
        if not filename or not os.path.exists(filename) or not os.path.isfile(filename):
            sys.stderr.write("input file '%s' does not exist or is not a regular file\n"
                             % line)
            exit(2)
        if util.is_gzipped(filename):
            f = gzip.open(filename, "rt")
        else:
            f = open(filename, "rt")

        # skip header
        f.readline()

        infiles.append(f)
    in_file.close()

    if len(infiles) == 0:
        sys.stderr.write("no input files specified in file '%s'\n" % options.infile_list)
        exit(2)

    return infiles



def main():
    options = parse_options()

    if options.remove_chromosomes and options.chrom_list==None:
        sys.stderr.write("Warning: no chromosomes to remove specified, skipping\n")

    sys.stderr.write("reading input filenames from %s\n" % options.infile_list)
    infiles = open_input_files(options.infile_list)
    outfile = open(options.out_file, "wt")
    dup_snp_warn = True

    sys.stderr.write("estimating allele-specific dispersion coefficients\n")
    # read input data and estimate dispersion coefficient
    # for one individual at a time
    for i in range(len(infiles)):
        cur_file = infiles[i]

        AS_ref = []
        AS_alt = []
        hetps =[]

        # header was already removed in the function 'open_input_files'
        #header = cur_file.readline()

        # combine allele-specific read counts into one large
        # array for this individual
        for line in cur_file:
            snpinfo = line.strip().split()
            
            # skip homozygous positions
            if snpinfo[5] != '1.00':
                continue

            ## remove some chromosomes, if --remove_chromosomes option is activated and chrom_list is specified
            if options.remove_chromosomes and options.chrom_list!=None:
                if snpinfo[0] in options.chrom_list:
                    continue
            ## keep only specified chromosomes, if --keep_chromosomes option is activated and chrom_list is specified
            if options.keep_chromosomes and options.chrom_list!=None:
                if not snpinfo[0] in options.chrom_list:
                    continue

            if snpinfo[12] != "NA":
                # read information aout target SNPs
                snp_locs = np.array([int(y.strip()) for y in snpinfo[10].split(';')],
                                    dtype=np.int32)
                snp_as_ref = np.array([int(y) for y in snpinfo[12].split(';')],
                                      dtype=np.int32)
                snp_as_alt = np.array([int(y) for y in snpinfo[13].split(';')],
                                      dtype=np.int32)
                snp_hetps = np.array([float(y.strip()) for y in snpinfo[11].split(';')],
                                     dtype=np.float64)

                # same SNP should not be provided multiple times, this
                # can create problems with combined test. Warn and filter
                # duplicate SNPs
                uniq_loc, uniq_idx = np.unique(snp_locs, return_index=True)

                if dup_snp_warn and uniq_loc.shape[0] != snp_locs.shape[0]:
                    sys.stderr.write("WARNING: discarding SNPs that are repeated "
                                     "multiple times in same line\n")
                    dup_snp_warn = False
                
                if options.aggregate_per_region:
                    AS_ref.extend([np.sum(snp_as_ref[uniq_idx])])
                    AS_alt.extend([np.sum(snp_as_alt[uniq_idx])])
                    hetps.extend([1]) # all counts come from heterozygous sites
                else:
                    AS_ref.extend(snp_as_ref[uniq_idx])
                    AS_alt.extend(snp_as_alt[uniq_idx])
                    hetps.extend(snp_hetps[uniq_idx])

        AS_ref = np.array(AS_ref)
        AS_alt = np.array(AS_alt)
        hetps = np.array(hetps)

        ms = f"{i}: {len(AS_ref)} heterozygous variants"
        print_with_time(ms)

        # find maximu likelihood estimate for overdispersion parameter
        res = minimize_scalar(likelihood, bounds=(0.01, 0.99),
                              args=(AS_ref, AS_alt, hetps, options.read_error_rate),
                              options = {'xatol' : 0.001},
                              method="Bounded")
        LL = res.fun
        dispersion = res.x

        sys.stderr.write("Sample %d:\n" % i)
        sys.stderr.write("  dispersion estimate: %g\n" % dispersion)
        sys.stderr.write("  LL: -%g\n" % LL)

        outfile.write(str(dispersion)+"\n")
        outfile.flush()

    sys.stderr.write("output written to file: %s\n" % options.out_file)
    sys.stderr.write("done!\n")
        


def likelihood(dispersion, AS_ref, AS_alt, hetps, error):
    cur_like = 0
    for i in range(len(AS_ref)):
        # calculate likelihood for each heterozygous site, under
        # assumption that true reference proportion is 50%
        cur_like +=  AS_betabinom_loglike([math.log(0.5), math.log(0.5)],
                                          dispersion, AS_ref[i],
                                          AS_alt[i], hetps[i], error)
    return -cur_like




def addlogs(loga, logb):
    """Helper function: perform numerically-stable addition in log space"""
    return max(loga, logb) + math.log(1 + math.exp(-abs(loga - logb)))



#Given parameters, returns log likelihood.  Note that some parts have been cancelled out
def AS_betabinom_loglike(logps, sigma, AS1, AS2, hetp, error):
    if sigma >= 1.0 or sigma <= 0.0:
        return -99999999999.0

    a = math.exp(logps[0] + math.log(1/sigma**2 - 1))
    b = math.exp(logps[1] + math.log(1/sigma**2 - 1))

    part1 = 0
    part1 += betaln(AS1 + a, AS2 + b)
    part1 -= betaln(a, b)

    if hetp == 1:
        return part1

    e1 = math.log(error) * AS1 + math.log(1 - error) * AS2
    e2 = math.log(error) * AS2 + math.log(1 - error) * AS1
    if hetp == 0:
        return addlogs(e1, e2)

    return addlogs(math.log(hetp)+part1, math.log(1-hetp) + addlogs(e1, e2))


if __name__ == "__main__":
    main()
