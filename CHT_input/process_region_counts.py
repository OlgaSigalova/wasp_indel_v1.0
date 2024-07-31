import polars as pl
import re
import os
import glob
import argparse


def extract_info(filename):
    pattern = r"(?P<cell_type>.+?)_(?P<donorID>ASA_[0-9]+)_(?P<n_cells>[0-9]+)"
    match = re.search(pattern, filename)
    if match:
        return match.group("cell_type"), match.group("donorID"), match.group("n_cells")
    else:
        return None, None, None
    

if __name__ == "__main__":
    # Create argument parser
    parser = argparse.ArgumentParser(
        description="This script parses output of generate_cht_input.py script and write a .txt file"
                    "with region AS and total counts for all donors and cell types"
    )
    parser.add_argument(
        "--counts", required=True, nargs='+', help="Files with region counts (strictly output of )."
    )
    parser.add_argument(
        "--outfile", required=True, help="Output file."
    )

    
    args = parser.parse_args()
    outfile = args.outfile
    count_files = args.counts

    # write header
    header_list = ["chr", "test_snp_pos", "test_snp_id", "REF", "ALT", "region_start", "region_end", "phasing_block" ,
                  "region_total_count", "genomewide_total_count", "region_ref_count", "region_alt_count", "region_as_count", 
                   "cell_type", "donor_id", "n_cells"]
    header = "\t".join(header_list)
    with open(outfile, "w") as file:
        file.write(f"{header}\n")

    # process all counts files
    for path in count_files:

        base_path = os.path.basename(path)
        cell_type, donor_id, n_cells = extract_info(base_path)
        #print(f"Cell Type: {cell_type}, Donor ID: {donorID}, Number of Cells: {n_cells}")

        df = pl.read_csv(path, separator = " ")

        df = (df.filter(pl.col("PHASING.BLOCK") != "None")
        .with_columns([
            pl.lit(cell_type).alias("cell_type"),
            pl.lit(donor_id).alias("donor_id"),
            pl.lit(n_cells).cast(pl.Int32).alias("n_cells"),
            pl.col("REGION.SNP.REF.HAP.COUNT").str.split(";").list.eval(pl.element().cast(pl.Int32)).list.sum().alias("region_ref_count"),
            pl.col("REGION.SNP.ALT.HAP.COUNT").str.split(";").list.eval(pl.element().cast(pl.Int32)).list.sum().alias("region_alt_count")
        ])
         .with_columns([(pl.col("region_ref_count") + pl.col("region_alt_count")).alias("region_as_count"),
                         pl.col("REGION.READ.COUNT").alias("region_total_count"),
                         pl.col("GENOMEWIDE.READ.COUNT").alias("genomewide_total_count")
                       ])
         .select(["CHROM", "TEST.SNP.POS", "TEST.SNP.ID", "TEST.SNP.REF.ALLELE", "TEST.SNP.ALT.ALLELE", "REGION.START", "REGION.END", "PHASING.BLOCK" ,
                  "region_total_count", "genomewide_total_count", "region_ref_count", "region_alt_count", "region_as_count", "cell_type", "donor_id", "n_cells"])
        )

         # Convert to pandas to append to file
        df = df.to_pandas()
        df.to_csv(outfile, mode='a', index=False, header=False, sep = "\t")

    