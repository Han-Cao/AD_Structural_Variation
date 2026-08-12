# Add VCF dosage for LD analysis

# For LD analysis, we set dosage as:
# by length: DS = x - min / (max - min) * 2
# by cutoff: allele_ds = 1 if nCCG > cutoff, else = 0

import argparse

import pandas as pd
import pysam


def refine_tr(var: pysam.VariantRecord, ds_dict: dict) -> pysam.VariantRecord:
    """ Refine TR dosage for FAM193B-CCG """

    new_var = var.copy()
    for x in var.samples:
        # if fail filter, set to missing
        if x not in ds_dict:
            new_var.samples[x]["DS"] = None
            new_var.samples[x]["GT"] = (None, None)
            continue

        # set dosage
        ds = ds_dict[x]
        new_var.samples[x]["DS"] = ds

        # set pseudo GT, this is useful for LD by SMR, which does not support dosage input
        if ds < 1/3:
            new_var.samples[x]["GT"] = (0, 0)
        elif ds >= 2/3:
            new_var.samples[x]["GT"] = (1, 1)
        else:
            new_var.samples[x]["GT"] = (0, 1)
    
    return new_var


def get_cutoff_dosage(df_dosage: pd.DataFrame, cutoff: int) -> dict:
    """ Generate dosage dict by cutoff """

    df_new = df_dosage[["eid", "a1_nCCG", "a2_nCCG"]].copy()
    # if allele_nCCG >= cutoff, 1
    # else 0
    def get_cutoff_ds(x):
        if x >= cutoff:
            return 1
        else:
            return 0

    df_new["DS"] = df_new["a1_nCCG"].apply(get_cutoff_ds) + df_new["a2_nCCG"].apply(get_cutoff_ds)

    # drop DS with missing
    df_new = df_new.dropna()

    return df_new.set_index("eid")["DS"].to_dict()


def get_new_var(temp_var: pysam.VariantRecord, id: str, pos: int) -> pysam.VariantRecord:
    new_var = temp_var.copy()
    new_var.id = id
    new_var.pos = pos
    new_var.ref = "G"
    new_var.alts = ("GCCG",)
    return new_var


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add dosage")
    parser.add_argument("-i", "--invcf", metavar="vcf", type=str, required=True,
                        help="input vcf")
    parser.add_argument('-d', '--dosage', metavar="txt", type=str, required=True,
                        help="dosage table file, require columns names IID and DS")
    parser.add_argument("-o", "--outvcf", metavar="vcf", type=str, required=True,
                        help="input vcf")
    parser.add_argument('--max-cap', metavar="int", type=int, required=False, default=None,
                        help="max cap for dosage")
    parser.add_argument('--min-cap', metavar="int", type=int, required=False, default=None,
                        help="min cap for dosage")

    args = parser.parse_args()

    # read input
    invcf = pysam.VariantFile(args.invcf, "r")
    df_dosage = pd.read_csv(args.dosage, sep="\t", dtype={"eid": str})

    # prepare dosage dict (IID -> DS)
    # if no cap, use the same dosage for regression analysis
    if args.max_cap is None:
        ds_dict_add = df_dosage[["eid", "DS_add"]].set_index("eid").to_dict()['DS_add']
    # if cap, re-calculate dosage
    else:
        # convert repeat count by cap
        df_dosage["total_nCCG_cap"] = df_dosage["total_nCCG"].apply(lambda x: min(x, args.max_cap))
        df_dosage["total_nCCG_cap"] = df_dosage["total_nCCG_cap"].apply(lambda x: max(x, args.min_cap))
        # min-max conversion to 0-2
        min_total_nCCG = df_dosage["total_nCCG_cap"].min()
        max_total_nCCG = df_dosage["total_nCCG_cap"].max()
        df_dosage["DS_cap"] = (df_dosage["total_nCCG_cap"] - min_total_nCCG) / (max_total_nCCG - min_total_nCCG) * 2
        ds_dict_add = df_dosage[["eid", "DS_cap"]].set_index("eid").to_dict()['DS_cap']


    # setup output VCF
    header = invcf.header
    header.add_line('##FORMAT=<ID=DS,Number=1,Type=Float,Description="Normalized repeat length dosage">')
    outvcf = pysam.VariantFile(args.outvcf, "w", header=header)

    temp_var = next(invcf)

    # generate dosage records
    pos_orig = 177554497 # true position of FAM193B-CCG, we generate pseudo pos to avoid duplication
    var_dosage_add = get_new_var(temp_var, "FAM193B-CCG-ADD", pos_orig)
    var_dosage_add = refine_tr(var_dosage_add, ds_dict_add)
    outvcf.write(var_dosage_add)

    # generate GT by cutoff
    add_pos = 1
    for x in range(12,36):
        add_pos += 1
        var_cutoff = get_new_var(temp_var, f"FAM193B-CCG-CUTOFF-{x}", pos_orig + add_pos)
        var_cutoff = refine_tr(var_cutoff, get_cutoff_dosage(df_dosage, x))

        outvcf.write(var_cutoff)

    invcf.close()
    outvcf.close()


    
