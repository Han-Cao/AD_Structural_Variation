# Add VCF dosage for regenie analysis

# For regression analysis, we set dosage as:
# if max_nCCG >= cutoff, 1
# elif max_nCCG <= 11, 0
# else, missing

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

        # set pseudo GT, this will not be used for analysis
        if ds <= 0.5:
            new_var.samples[x]["GT"] = (0, 0)
        elif ds >= 1.5:
            new_var.samples[x]["GT"] = (1, 1)
        else:
            new_var.samples[x]["GT"] = (0, 1)
    
    return new_var


def get_cutoff_dosage(df_dosage: pd.DataFrame, cutoff: int) -> dict:
    """ Generate dosage dict by cutoff """

    df_new = df_dosage[["eid", "max_nCCG"]].copy()
    # if max_nCCG >= cutoff, 1
    # elif max_nCCG <= 11, 0
    # else, missing
    def get_cutoff_ds(x):
        if x >= cutoff:
            return 1
        elif x <= 11:
            return 0
        else:
            return None

    df_new["DS"] = df_new["max_nCCG"].apply(get_cutoff_ds)

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
    args = parser.parse_args()

    # read input
    invcf = pysam.VariantFile(args.invcf, "r")
    df_dosage = pd.read_csv(args.dosage, sep="\t", dtype={"eid": str})

    # prepare dosage dict (IID -> DS)
    ds_dict_add = df_dosage[["eid", "DS_add"]].set_index("eid").to_dict()['DS_add']
    ds_dict_dom = df_dosage[["eid", "DS_dom"]].set_index("eid").to_dict()['DS_dom']

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

    var_dosage_dom = get_new_var(temp_var, "FAM193B-CCG-DOM", pos_orig + 1)
    var_dosage_dom = refine_tr(var_dosage_dom, ds_dict_dom)
    outvcf.write(var_dosage_dom)

    # generate GT by cutoff
    add_pos = 1
    for x in range(12,36):
        add_pos += 1
        var_cutoff = get_new_var(temp_var, f"FAM193B-CCG-CUTOFF-{x}", pos_orig + add_pos)
        var_cutoff = refine_tr(var_cutoff, get_cutoff_dosage(df_dosage, x))

        outvcf.write(var_cutoff)

    invcf.close()
    outvcf.close()


    
