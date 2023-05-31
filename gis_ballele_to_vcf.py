import os
import json
import argparse
import csv


def main(ballele_file, vcf_file, rs_gene):
    print("B-Allele:", ballele_file)
    print("VCF Output:", vcf_file)

    if rs_gene:
        with open(rs_gene) as f:
            rs_gene = json.loads(f.read())
    else:
        rs_gene = {}

    with open(ballele_file) as ballele_handle, open(vcf_file, 'w') as vcf_handle:
        ballele_reader = csv.DictReader(
            ballele_handle,
            fieldnames=["chrom", "pos", "rs", "freq"],
            delimiter="\t"
        )
        vcf_handle.write("##fileformat=VCFv4.2\n")
        vcf_handle.write("""##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (counting only informative reads out of the total reads) for the ref and alt alleles in the order listed">\n""")
        vcf_handle.write("""##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions for alt alleles in the order listed">\n""")
        vcf_handle.write("""##FORMAT=<ID=VF,Number=A,Type=Float,Description="Allele fractions for alt alleles in the order listed">\n""")
        vcf_handle.write("""##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">\n""")
        vcf_handle.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample\n")
        vcf_writer = csv.DictWriter(
            vcf_handle,
            fieldnames=["chrom","pos", "id", "ref", "alt", "qual", "filter", "info", "format", "sample"],
            delimiter="\t"
        )
        for row in ballele_reader:
            chrom = row["chrom"]
            chrom = chrom.replace("chr", "")
            
            pos = row["pos"]
            
            rs = row["rs"]
            if not rs.startswith("rs"):
                continue
            gene = rs_gene.get(rs, rs)

            freq = row["freq"]
            try:
                freq = float(freq)
            except Exception as e:
                print(f"Could float({freq}):", e)
                continue

            ref_freq = 1 - freq
            ref_dp = int(100 * ref_freq)
            alt_dp = int(100 * freq)
            dp = ref_dp + alt_dp

            ref = "N"
            alt = "N"
            qual = "."
            fltr = "PASS"
            info = ";".join([
                f"DP={dp}",
                f"AF={freq}",
                f"VF={freq}",
            ])

            frmt = "AD:VF:AF"

            sample = ":".join([
                f"{ref_dp},{alt_dp}",
                f"{freq}",
                f"{freq}",
            ])

            vcf_writer.writerow({
                "chrom": chrom,
                "pos": pos,
                "id": gene,
                "ref": ref,
                "alt": alt,
                "qual": qual,
                "filter": fltr,
                "info": info,
                "format": frmt,
                "sample": sample
            })

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="Illumina Gis ballele file to a VCF that can be read by reconCNV",
        description="A script to convert the TSO500 Logs_Intermediates/Gis/ to a reconCNV .VCF file",
    )

    parser.add_argument(
        "--ballele", 
        help="The input <sample>_bAllele.csv",
        required=True
    )
    parser.add_argument(
        "--vcf",
        help="The output .vcf ready to be used by reconCNV",
        required=True
    )
    parser.add_argument(
        "--rs-gene",
        help="A json file with {rs#: gene}",
        default=None
    )

    args = parser.parse_args()

    ballele_file = args.ballele
    vcf_file = args.vcf
    rs_gene = args.rs_gene

    main(ballele_file, vcf_file, rs_gene)