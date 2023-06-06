import json
import argparse
import csv

def main(logratio_file, cnr_file, rs_gene):
    print("Logratio:", logratio_file)
    print("CNR Output:", cnr_file)

    if rs_gene:
        with open(rs_gene) as f:
            rs_gene = json.loads(f.read())
    else:
        rs_gene = {}

    # so we can sort everything at the end
    lines = {}

    with open(logratio_file) as logratio_handle, open(cnr_file, 'w') as cnr_handle:
        logratio_reader = csv.DictReader(
            logratio_handle,
            fieldnames=["chrom", "pos", "rs", "log2"],
            delimiter="\t"
        )
        
        for row in logratio_reader:
            chrom = row["chrom"]
            chrom = chrom.replace("chr", "")
            
            pos = row["pos"]
            rs = row["rs"]
            if not rs.startswith("rs"):
                continue
            gene = rs_gene.get(rs, rs)
            log2 = row["log2"]

            if chrom not in lines:
                lines[chrom] = []

            lines[chrom].append({
                "chromosome": chrom,
                "start": pos,
                "end": pos,
                "gene": gene,
                "log2": log2,
                "depth": "0.1",
                "weight": "0.4"
            })

        # sort only by chr, within chromosomes everything should still be sorted
        chroms = [c.replace("chr", "") for c in lines.keys()]
        chroms = sorted(chroms, key=lambda c: int(c) if c.isnumeric() else ord(c))

        cnr_writer = csv.DictWriter(
            cnr_handle,
            fieldnames=["chromosome", "start", "end", "gene", "log2", "depth", "weight"],
            delimiter="\t"
        )
        cnr_writer.writeheader()
        for chrom in chroms:
            for line in lines[chrom]:
                cnr_writer.writerow(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="Illumina Gis raw file to CNR",
        description="A script to convert the TSO500 Logs_Intermediates/Gis/ to a reconCNV .cnr file",
    )

    parser.add_argument(
        "--logratio", 
        help="The input <sample>_logRatio.csv from the Gis directory",
        required=True
    )
    parser.add_argument(
        "--cnr",
        help="The output .cnr ready to be used by reconCNV",
        required=True
    )
    parser.add_argument(
        "--rs-gene",
        help="A json file with {rs#: gene}",
        default=None
    )

    args = parser.parse_args()

    logratio_file = args.logratio
    cnr_file = args.cnr
    rs_gene = args.rs_gene

    main(logratio_file, cnr_file, rs_gene)