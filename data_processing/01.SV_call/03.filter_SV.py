# add ALT to INFO, edit ALT to SVTYPE, filter vcf
import pysam
import argparse
import subprocess

# parse arguments
parser = argparse.ArgumentParser(description="SV hardfilter")
parser.add_argument("-i", "--invcf", metavar="vcf", type=str, required=True,
                    help="input vcf")
parser.add_argument("-o", "--outvcf", metavar="vcf", type=str, required=True,
                    help="input vcf")
parser.add_argument("--mode", choices=["Sniffles", "SVIM", "cuteSV"], required=True,
                    help="Select mode for Sniffles | SVIM | cuteSV")
parser.add_argument("--prefix", metavar="Sample.aligner", required=True,
                    help="Prefix name to add to SV ID")
parser.add_argument("--maxlen", metavar="n", type=int, required=True, default=1000000,
					help="Maximum SVLEN")
parser.add_argument("--qual", metavar="n", type=int, required=False,
                    help="Minimum QUAL score for svim")
parser.add_argument("--chr-prefix", type=str, required=False, default='chr',
                    help="Prefix for chromosome names")
args = parser.parse_args()

# read input file
myvcf = pysam.VariantFile(args.invcf, "r")
id_prefix = args.prefix

# set keep chrs
keep_chr = set([f"{args.chr_prefix}{x}" for x in range(1,23)]).union([f"{args.chr_prefix}X"])

# write new vcf
myvcf.header.info.add("SEQ_raw", ".", "String", "Raw insertion sequence")
outvcf = pysam.VariantFile(args.outvcf, "w", header=myvcf.header)

# process Sniffles output
if args.mode == "Sniffles":
    for variant in myvcf:
        if variant.chrom not in keep_chr:
            continue
        if variant.info["SVTYPE"] == "BND":
            continue

        SVLEN = abs(variant.info["SVLEN"])
        SVTYPE = variant.info["SVTYPE"]
        
        if SVLEN > args.maxlen:
            continue

        if SVTYPE == "INS":
            if variant.alts[0] != "<INS>":
                variant.info["SEQ_raw"] = variant.alts[0]

        variant.id = f"{id_prefix}.{variant.id}"
        variant.alts = (f"<{SVTYPE}>",)

        outvcf.write(variant)

elif args.mode == "SVIM":
    for variant in myvcf:
        if variant.chrom not in keep_chr:
            continue
        if variant.qual < args.qual or variant.info["SVTYPE"] == "BND":
            continue
        
        SVTYPE = variant.info["SVTYPE"]

        if SVTYPE == "INV":
            variant.info["SVLEN"] = len(variant.alts[0])

        SVLEN = abs(variant.info["SVLEN"])
        
        if SVLEN > args.maxlen:
            continue
        
        if SVTYPE == "INS":
            if variant.alts[0] != "<INS>":
                variant.info["SEQ_raw"] = variant.alts[0]

        # convert DUP:INT DUP:TANDAM to DUP
        if "DUP" in SVTYPE:
            SVTYPE = "DUP"
            variant.info["SVTYPE"] = SVTYPE

        variant.id = f"{id_prefix}.{variant.id}"
        variant.alts = (f"<{SVTYPE}>",)

        outvcf.write(variant)
        
elif args.mode == "cuteSV":
    for variant in myvcf:
        if variant.chrom not in keep_chr:
            continue
        if "q5" in variant.filter.keys() or variant.info["SVTYPE"] == "BND":
            continue

        SVLEN = abs(variant.info["SVLEN"])
        SVTYPE = variant.info["SVTYPE"]
        
        if SVLEN > args.maxlen:
            continue

        if SVTYPE == "INS":
            if variant.alts[0] != "<INS>":
                variant.info["SEQ_raw"] = variant.alts[0]

        variant.id = f"{id_prefix}.{variant.id}"
        variant.alts = (f"<{SVTYPE}>",)

        outvcf.write(variant)
        
outvcf.close()

# sort and index output vcf
subprocess.run(["bcftools", "sort",
                "-O", "z",
                "-o", args.outvcf + ".tmp", args.outvcf])
subprocess.run(["mv", args.outvcf + ".tmp", args.outvcf])
subprocess.run(["tabix", args.outvcf])
