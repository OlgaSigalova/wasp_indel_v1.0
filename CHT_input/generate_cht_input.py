import argparse
import collections
import counts
import gzip


class InputReader(object):

    def __init__(
        self, path
    ):
        # Store path and open file
        self.path = path
        if self.path.endswith('.gz'):
            self.infile = gzip.open(self.path, 'rt')
        else:
            self.infile = open(self.path, 'rt')
        # Create named tuple to store test variants
        self.test = collections.namedtuple(
            'test', [
                'chrom', 'test_start', 'test_end', 'ref', 'alt',
                'region_starts', 'region_ends'
            ]
        )

    def get_tests(
        self
    ):
        # Loop through lines in infile and get data
        for line in self.infile:
            # Extract line data
            line_list = line.strip().split('\t')
            chrom = line_list[0]
            test_start = int(line_list[1]) - 1
            ref, alt = line_list[2:4]
            test_end = test_start + len(ref)
            region_starts = [int(rs) - 1 for rs in line_list[4].split(';')]
            region_ends = [int(re) for re in line_list[5].split(';')]
            # Create tuple and return
            test = self.test(
                chrom=chrom, test_start=test_start, test_end=test_end, ref=ref,
                alt=alt, region_starts=region_starts, region_ends=region_ends
            )
            yield(test)

    def close(
        self
    ):
        self.infile.close()


class OutputWriter(object):

    def __init__(
        self, path
    ):
        # Store input arguments
        self.path = path
        # Open outfile
        if self.path.endswith('.gz'):
            self.outfile = gzip.open(self.path, 'wt')
        else:
            self.outfile = open(self.path, 'wt')
        # Open files and write header to each
        self.header = (
            "CHROM TEST.SNP.POS TEST.SNP.ID TEST.SNP.REF.ALLELE "
            "TEST.SNP.ALT.ALLELE TEST.SNP.GENOTYPE TEST.SNP.HAPLOTYPE "
            "REGION.START REGION.END REGION.SNP.POS REGION.SNP.HET.PROB "
            "REGION.SNP.REF.HAP.COUNT REGION.SNP.ALT.HAP.COUNT REGION.SNP.OTHER.HAP.COUNT "
            "REGION.READ.COUNT GENOMEWIDE.READ.COUNT\n"
        )
        self.outfile.write(self.header)

    def close(
        self
    ):
        self.outfile.close()

    def write(
        self, test_str, region_str, read_str
    ):
        # Create blank output line
        line_list = [test_str, region_str, read_str]
        line_str = ' '.join(line_list) + '\n'
        self.outfile.write(line_str)


if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(
        description="This program generates target regions for subsequent "
        "analysis using the CHT test. The input count files should be "
        "generated using the 'get_counts.py' script and each count file "
        "should contain data for an identical set of variants."
    )
    parser.add_argument(
        "--variants", required=True, help="Input variant count file."
    )
    parser.add_argument(
        "--bam", required=True, help="Input variant BAM file."
    )
    parser.add_argument(
        "--regions", required=True, help="Input region file."
    )
    parser.add_argument(
        "--outfile", required=True, help="Output file."
    )
    parser.add_argument(
        "--adjust_hetp", default=None, choices=['as_counts', 'total_counts'],
        help="Adjust heterozygous probabilities."
    )
    args = parser.parse_args()
    # Open input file and output files
    variants = counts.CountTree(args.variants, args.adjust_hetp)
    bam = counts.BamCounts(args.bam)
    regions = InputReader(args.regions)
    outfile = OutputWriter(args.outfile)
    current_chrom = None
    # Loop through tests in region file
    for test in regions.get_tests():
        # Read in data for new chromosome
        if current_chrom != test.chrom:
            variants.read_counts(test.chrom)
            current_chrom = test.chrom
        # Get test variant and its haplotype
        test_str = variants.get_test_string(
            start=test.test_start, end=test.test_end, ref=test.ref,
            alt=test.alt
        )
        test_haplotype = test_str.split(' ')[-2]
        test_phasing_block = test_str.split(' ')[-1]
        # Get region variants
        region_str = variants.get_region_string(
            test_haplotype=test_haplotype, test_phasing_block = test_phasing_block, 
            starts=test.region_starts, ends=test.region_ends
        )
        # Get counts
        read_str = bam.get_read_string(
            chrom=test.chrom, starts=test.region_starts,
            ends=test.region_ends
        )
        # Write output file
        outfile.write(
            test_str=test_str, region_str=region_str, read_str=read_str
        )
    # Close files
    bam.close()
    regions.close()
    outfile.close()
