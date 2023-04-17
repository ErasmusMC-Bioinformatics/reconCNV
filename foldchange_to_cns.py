import argparse
import csv
import math

from collections import defaultdict

def get_gene_pos_dict(norm_bin_file):
    gene_pos_dict = defaultdict(lambda: dict(start=2**40, end=0))
    with open(norm_bin_file) as norm_bin_handle:
        norm_bin_reader = csv.DictReader(
            norm_bin_handle,
            delimiter="\t",
        )

        for row in norm_bin_reader:
            chrom = row["Chr"]
            start = int(row["Start"])
            end = int(row["End"])
            gene = row["Gene"]

            gene_pos = gene_pos_dict[gene]
            
            gene_pos["chr"] = chrom
            gene_pos["start"] = min([gene_pos["start"], start])
            gene_pos["end"] = max([gene_pos["end"], end])
    
    return gene_pos_dict

def create_recon_cns_output(fold_change_file, gene_pos_dict, cns_output):
        with open(fold_change_file) as fold_change_handle, open(cns_output, 'w', newline='', encoding='utf-8') as cns_handle:
            next(fold_change_handle) # gender line
            next(fold_change_handle) # gender score line
            input_reader = csv.DictReader(
                fold_change_handle,
                delimiter="\t",
            )
            cns_writer = csv.DictWriter(
                cns_handle,
                fieldnames=["chromosome", "start", "end", "gene", "log2", "depth", "weight", "probes", "mcn", "lcn", "cf"],
                delimiter="\t",
            )
            cns_writer.writeheader()
            for row in input_reader:
                gene = row["GeneName"]
                if gene not in gene_pos_dict:
                    print(f"{gene} not in gene_pos_dict")
                    continue
                
                gene_ref = gene_pos_dict[gene]
                chrom = gene_ref["chr"]
                start = gene_ref["start"]
                end = gene_ref["end"]
                
                fold_change = row["FC"]
                try:
                    fold_change = float(fold_change)
                except ValueError as e:
                    print(f"Could not convert {fold_change} to float")
                    continue
                log2_fc = math.log2(fold_change)
                
                qScore = row["qScore"]
                
                cns_writer.writerow({
                    "chromosome": chrom[3:],
                    "start": start,
                    "end": end,
                    "gene": gene,
                    "log2": log2_fc,
                    # "depth": "0.1",
                    # "weight": "0.4",
                    # "probes": "500",
                    "mcn": "1",
                    "lcn": "1",
                    # "cf": "0.8"
                })
    
def main(normbin_file, fold_change_file, cns_file):
    print("NormalizedBin:", normbin_file)
    print("FoldChange:", fold_change_file)
    print("CNS Output:", cns_file)

    gene_pos_dict = get_gene_pos_dict(normbin_file)

    create_recon_cns_output(fold_change_file, gene_pos_dict, cns_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="Normbin to CNR",
        description="A script to convert the TSO500 <sample>_normalizedBinCount.tsv to a reconCNV .cnr file",
    )

    parser.add_argument(
        "--normbin", 
        help="The input <sample>_normalizedBinCount.tsv",
        required=True
    )
    parser.add_argument(
        "--fold-change", 
        help="The input <sample>_normalizedBinCount.tsv",
        required=True
    )
    parser.add_argument(
        "--cns",
        help="The output .cns ready to be used by reconCNV",
        required=True
    )

    args = parser.parse_args()

    normbin_file = args.normbin
    fold_change_file = args.fold_change
    cns_file = args.cns

    main(normbin_file, fold_change_file, cns_file)