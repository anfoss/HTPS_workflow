# !/usr/bin/env python3

import pandas as pd
import sys
import math as m
import numpy as np
import glob
import seaborn as sns
import matplotlib.pyplot as plt
import pylab
from scipy.spatial import distance_matrix


def explode(df, lst_cols, fill_value=""):
    # make sure `lst_cols` is a list
    if lst_cols and not isinstance(lst_cols, list):
        lst_cols = [lst_cols]
    # all columns except `lst_cols`
    idx_cols = df.columns.difference(lst_cols)

    # calculate lengths of lists
    lens = df[lst_cols[0]].str.len()

    if (lens > 0).all():
        # ALL lists in cells aren't empty
        return (
            pd.DataFrame(
                {
                    col: np.repeat(df[col].values, df[lst_cols[0]].str.len())
                    for col in idx_cols
                }
            )
            .assign(**{col: np.concatenate(df[col].values) for col in lst_cols})
            .loc[:, df.columns]
        )
    else:
        # at least one list in cells is empty
        return (
            pd.DataFrame(
                {
                    col: np.repeat(df[col].values, df[lst_cols[0]].str.len())
                    for col in idx_cols
                }
            )
            .assign(**{col: np.concatenate(df[col].values) for col in lst_cols})
            .append(df.loc[lens == 0, idx_cols])
            .fillna(fill_value)
            .loc[:, df.columns]
        )


def merge_cleavage_windows(row):
    """Receive a row from MQ and returns the correct sequence of N-C window

    Args:
    row: pd apply series from a peptides.txt file

    Returns:
    N + C term window without overlapping sequences
    """
    nterm = list(row["N-term cleavage window"])
    cterm = list(row["C-term cleavage window"])
    # inefficient but get the job done
    i = 0
    for idx, x in enumerate(cterm):
        idx += 1
        # assumption is that the end of N is the beginning of C
        if nterm[:-idx] == x:
            i += 1
        else:
            break
    seq = nterm + cterm[:i]
    return "".join(seq)


def convert_to_cleavage(filename, search="MQ"):
    """Receive a search file and returns a formatted cleavage matrix

    Args:
    filename: txt file generated from a search engine
    search: flag signalling the search engine used

    Returns:
    formatted cleavage matrix
    """
    df = pd.read_table(filename, sep="\t")
    df = df[df["Intensity"] > 0]
    # df = df[df['Intensity 1']>0]
    if search == "MQ":
        df = df[(df["Score"] > 40) & (df["PEP"] <= 0.05)]
        df = df[["N-term cleavage window", "C-term cleavage window"]]
        df = df[df["N-term cleavage window"] != "________________"]
        df = df[df["C-term cleavage window"] != "________________"]
        df.dropna(inplace=True)
        # C-term is in N-term window only thing to do is to remove
        # df = df['N-term cleavage window'] + df['C-term cleavage window']
        df = df.apply(merge_cleavage_windows, axis=1)
        df = df.drop_duplicates()
        df = pd.DataFrame(list(df.str.split("").dropna()))
        df.drop([0, 17], axis=1, inplace=True)
    elif search == "MSFragger":
        pass
    elif search == "custom":
        pass
    w = len(list(df)) // 2
    # convert it to the form p-1 and p 1' for p prime
    df.columns = [-x for x in range(1, w + 1)][::-1] + [x for x in range(1, w + 1)]
    return df


def expand_counts(df):
    """Receive a cleavage matrix and add missing aa with 0 counts per position

    Args:
    df: cleavage matrix

    Returns:
    filled cleavage matrix
    """

    aa = [
        "A",
        "R",
        "N",
        "D",
        "C",
        "Q",
        "E",
        "G",
        "H",
        "I",
        "L",
        "K",
        "M",
        "F",
        "P",
        "S",
        "T",
        "W",
        "Y",
        "V",
    ]
    # rmeove - nan or no aa
    df = df.loc[df.index.isin(aa)]
    for k in aa:
        if k not in df.index:
            row = dict(zip(list(df.columns), [0] * len(list(df.columns))))
            # df = df.append(pd.Series(row, index=df.columns, name=k))
            df = df.append(row, ignore_index=True)
    return df


def entropy(arr):
    """Calculates shannon entropy per position

    Args:
    arr: a position (column) from a cleavage matrix

    Returns:
    one entropy value per position
    """
    import math as m

    d = [m.log(x, 20) * x for x in list(arr.values)]
    entr = -1 * np.nansum(np.array(d))
    return entr


def subpocket_e(entr_df):

    v = entr_df.apply(lambda col: entropy(col), axis=0)
    v = pd.Series([(x - min(v)) / (max(v) - min(v)) for x in v])
    return v


def blocks_e(entropy_arr):
    """Calculates subpocket entropy per position according to
    https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003007
    local entropy by summing block wise the entropies
    i.e p4 to p1 available
    4 blocks
    b4 = entropy(p4,p3,p2,p1)
    b3 = entropy(p3,p2,p1)
    b2 = entropy (p2,p1)
    b1 = entropy (p1)

    Args:
    entropy_arr: entropy series

    Returns:
    block entropy per position
    """
    tmp = list(entropy_arr)
    w = len(tmp) // 2
    left = [sum(tmp[x:w]) for x in range(0, w)]
    right = [sum(tmp[w:x]) for x in range(w + 1, len(tmp) + 1)]
    return pd.Series(left + right)


def get_htps_db(htps_db):
    from Bio import SeqIO

    sq = []
    for record in SeqIO.parse(htps_db, "fasta"):
        sq.append(str(record.seq))
    return list("".join(sq))


def sequence_decoy(htps_db, df, n=50):
    """generate a cleavage matrix by random sampling peptides in htps DB
    Args:
    htps db: list of aa from string concatenated htps db
    df: cleavage matrix (peptide  level)
    n: number of iterations

    Returns:
    normalized frequency matrix
    """

    # generate random cleavage matrix N times
    mtrx = []
    nrow, ncol = df.shape
    seed = 111
    for x in range(1, n):
        np.random.seed(seed)
        seed += 1
        arr = np.empty(df.shape, dtype="object")
        for i in range(0, ncol):
            arr[:, i] = np.random.choice(htps_db, nrow)
        dec_df = pd.DataFrame(arr)
        dec_df = dec_df.apply(pd.value_counts, result_type="expand")
        dec_df = expand_counts(dec_df)
        mtrx.append(dec_df)
    # average frequency across all matrixes
    base_df = mtrx.pop()
    for m1 in mtrx:
        base_df = np.add(base_df, m1)
    base_df = base_df / (len(mtrx) + 1)
    return base_df


def normalize_aa_abundance(df, aa):
    """normalize a cleavage frequency matrix by the frequency of aa in htps db fasta
    Args:
    df: cleavage frequency matrix

    Returns:
    normalized frequency matrix
    """
    df = df.apply((lambda x: x / aa["value"]), axis=0)
    df = df.apply((lambda x: x / np.sum(x)), axis=0)
    return df


def convert_cleavage_matrx(cleav_mtrx, pref, write=True):
    """process a single cleavage matrix into its specificity, entropy and block entropy
    Args:
    cleav_mtrx: cleavage frequency matrix (n aa, n sites)
    pref: string prefix to be added to the specificity, cleavage and block entropy df
    write: bool, write output to file

    Returns:
    spec_df = cleavage specificity dataframe per position (aa, n sites)
    entr_df = entropy per position (n sites)
    block_df = block entropy per position (n sites)
    """
    # sort to ensure same order between all matrixes
    cleav_mtrx.sort_index(inplace=True)
    spec_df = cleav_mtrx.apply((lambda x: x / np.sum(x)), axis=0)
    aa = pd.read_csv("test/distribution_AA.csv").set_index("AA")
    cleav_mtrx = normalize_aa_abundance(cleav_mtrx, aa)
    entr_mtrx = cleav_mtrx.apply(entropy, axis=0)
    bl_entr_mtrx = blocks_e(entr_mtrx)
    entr_mtrx.name='{} Entropy'.format(pref)
    bl_entr_mtrx.name='{} Block entropy'.format(pref)
    if write:
        cleav_mtrx.to_csv("{}_cleavage_prob.csv".format(pref))
        entr_mtrx.to_csv("{}_entropy.csv".format(pref))
        bl_entr_mtrx.to_csv("{}_block_entropy.csv".format(pref))
        spec_df.to_csv("{}_specificity.csv".format(pref))
    return spec_df, cleav_mtrx, entr_mtrx, bl_entr_mtrx


def main():
    """
    process a single file
    """
    cleav_mtrx = convert_to_cleavage("test/peptides.txt", "MQ")
    db = get_htps_db("test/HTPS_db.fasta")
    decoys = sequence_decoy(db, cleav_mtrx, 50)

    cleav_mtrx = cleav_mtrx.apply(pd.value_counts, result_type="expand")
    cleav_mtrx = expand_counts(cleav_mtrx)
    spec_df, cleav_mtrx, entr_mtrx, bl_entr_mtrx = convert_cleavage_matrx(
        cleav_mtrx, "Tryp", True
    )
    spec_dec, cleav_dec, entr_dec, bl_entr_dec = convert_cleavage_matrx(
        decoys, "", False
    )
    dff = spec_df.values - spec_dec.values
    dff = pd.DataFrame(dff, index=spec_df.index)


if __name__ == "__main__":
    main()
