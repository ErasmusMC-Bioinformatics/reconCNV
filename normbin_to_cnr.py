import argparse
import csv
import math

def get_manifest_id_dict(manifest_file):
    id_dict = {}
    with open(manifest_file) as manifest_handle:
        manifest_reader = csv.DictReader(
            manifest_handle,
            fieldnames=["chr", "start", "end", "id"],
            delimiter="\t"
        )
        for row in manifest_reader:
            chrom = row["chr"]
            chrom = chrom[3:]
            start = int(row["start"])
            end = row["end"]
            manifest_id = row["id"]

            id_dict[(chrom, str(start+1), end)] = manifest_id

    return id_dict

def create_recon_cnr_output(norm_bin_file, cnr_output, manifest_id_dict):
    with open(norm_bin_file) as norm_bin_handle, open(cnr_output, 'w', newline='', encoding='utf-8') as cnr_output_handle:
        norm_bin_reader = csv.DictReader(
            norm_bin_handle,
            delimiter="\t"
        )
        cnr_output_writer = csv.DictWriter(
            cnr_output_handle,
            fieldnames=["chromosome", "start", "end", "gene", "log2", "depth", "weight"],
            delimiter="\t"
        )
        cnr_output_writer.writeheader()

        for row in norm_bin_reader:
            chrom = row["Chr"]
            chrom = chrom.replace("chr", "")

            start = row["Start"]
            end = row["End"]

            gene = row["Gene"]

            # get the full <gene>_<exon>_<transcript>
            position = (chrom, start, end) 
            if position in manifest_id_dict:
                gene = manifest_id_dict[position]

            norm_count = row["NormalizedCount"]
            try:
                norm_count = float(norm_count)
                log2 = math.log2(norm_count)
            except ValueError as e:
                log2 = 0

            depth = "0.1"
            weight = "0.4"
            
            cnr_output_writer.writerow({
                "chromosome": chrom,
                "start": start,
                "end": end,
                "gene": gene,
                "log2": log2,
                "depth": depth,
                "weight": weight,
            })
    
def main(normbin_file, manifest_file, cnr_file):
    print("NormalizedBin:", normbin_file)
    print("Manifest:", manifest_file)
    print("CNR Output:", cnr_file)

    manifest_id_dict = get_manifest_id_dict(manifest_file)

    create_recon_cnr_output(normbin_file, cnr_file, manifest_id_dict)

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
        "--manifest", 
        help="The Illumina TST500C_manifest.bed", 
        default="/staging/TSO500_RUO_LocalApp/resources/TST500C_manifest.bed",
    )
    parser.add_argument(
        "--cnr",
        help="The output .cnr ready to be used by reconCNV",
        required=True
    )

    args = parser.parse_args()

    normbin_file = args.normbin
    manifest_file = args.manifest
    cnr_file = args.cnr

    main(normbin_file, manifest_file, cnr_file)