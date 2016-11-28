import itertools
import logging
import os
import re
import subprocess

import numpy as np
import pandas as pd
import scipy.stats as stats
from Bio import SeqIO

logger = logging.getLogger(__name__)


def run(parser, args):
    """
    Generate a network from an MCL file of clustered proteins by
    guilt-by-association, as described in Lima-Mendez
    :param parser:
    :param args:
    :return:
    """
    s = parse_mcl_dump_file(args.mcl_dump)
    t = load_gene_map(args.gene_map)
    u = add_pc_labels(s, t)
    v = create_composition_mtx(u, presence_absence=True)
    w = calculate_shared_matrix(v)
    x = calculate_hypergeometric_survival(u, w)
    x.to_csv(args.output)
    return x




def extract_fastas_from_list(fasta_file, id_list, list_file=True):
    """
    Given a list of ids, and a fasta file (can be zipped)
    this will return a list of SeqRecord objects
    that match that id
    :param fasta_file: The fasta file to be parsed
    :param id_list: The filename containing a list of ids to be parsed
    :param list_file: True, if the list is in a file.
    :return: A list of SeqRecords
    """

    file_dir = os.path.dirname(fasta_file)
    zipped = False

    m = re.search('(.*)\.gz$', fasta_file)
    if m:
        zipped = True
        file_to_read = m.group(1)
        subprocess.call(["gunzip", file_dir + "/" + os.path.basename(fasta_file)])
    else:
        file_to_read = fasta_file

    reads = SeqIO.index(file_to_read, 'fasta')

    rtn_value = []

    ids = []

    if list_file:
        with open(id_list, 'r') as handle:
            for line in handle.readlines():
                ids.append(line.strip())
    else:
        ids = id_list

    for i in ids:
        try:
            rtn_value.append(reads[i])
        except KeyError:
            logger.error('Could not find read %s in %s' % (i, fasta_file))

    logger.info('Found %i reads in %s' % (len(rtn_value), fasta_file))

    if zipped:
        subprocess.call(["gzip", file_dir + "/" + os.path.basename(file_to_read)])

    return rtn_value


def parse_mcl_dump_file(dump_file):
    """
    Reads in a dump file from MCL and converts it into a dictionary of
    protein clusters
    :param dump_file: The file to be read in
    :return: a panda dataframe
    """
    counter = 1
    rows = []
    with open(dump_file, 'r') as handle:
        for line in handle.readlines():
            bits = line.strip().split()
            for i in bits:
                rows.append(('PC_%05d' % counter, i))
            counter += 1
    return pd.DataFrame.from_records(rows, columns=['PCid', 'node'])


def load_gene_map(map_file):
    """
    Reads in a file mapping genes to contigs
    :param map_file: The mapping file with gene ids and contigs comma-separated
    :return: a pandas dataframe
    """
    return pd.read_csv(map_file)


def add_pc_labels(pc_dict, gm_df):
    """
    Merge the PC dictionary with the gene map data frame
    :param pc_dict: The dictionary of PCs created with inverted mcl parse
    :param gm_df: The dataframe gene map
    :return: a pandas dataframe
    """
    return pd.merge(gm_df, pc_dict, how='left', on='node')


def create_composition_mtx(pc_df, presence_absence=False):
    """
    Create a composition matrix where M_ij if contig i encodes at least 1 PC j
    :param pc_df: a pandas dataframe of contigs, nodes and PC.ids
    :param presence_absence: If True, rescales values >1 to 1
    :return: an i x j matrix (rows are contigs, columns are PCs
    """
    v = pc_df.groupby(['contig', 'PCid']).size().reset_index()
    v.columns = ['contig', 'PCid', 'count']
    rtnValue = v.pivot(index='contig', columns='PCid', values='count').fillna(0)
    if presence_absence:
        rtnValue[rtnValue > 1] = 1

    return rtnValue


def calculate_shared_content(taxon1, taxon2, pc_df):
    """
    Calculate how many PCs two taxa share
    :param taxon1: The first taxon
    :param taxon2: The second taxon
    :param pc_df: The PC pandas dataframe
    :return: The number of shared PCs
    """
    t_pc = pc_df.transpose()
    t_pc = t_pc[[taxon1, taxon2]]
    rowsums = t_pc.sum(axis=1, numeric_only=True)
    try:
        rtnValue = rowsums.value_counts()[2.0]
    except KeyError:
        rtnValue = 0
    return rtnValue


def calculate_shared_matrix(pc_df):
    """
    Calculate a matrix where each contig is compared to each other contig and
    the number of shared PCs is stored as an upper triangle matrix
    :param pc_df: The pandas dataframe containing the contig / PC data
    :return: an n x n symmetrical identity matrix where n is the number of contigs
    """
    contigs = list(pc_df.index)
    rtnValue = pd.DataFrame(0, index=contigs, columns=contigs)
    for i in itertools.combinations(contigs, 2):
        rtnValue.loc[i[0], i[1]] = calculate_shared_content(i[0], i[1], pc_df)

    # Make it symmetrical
    rtnValue = rtnValue.as_matrix() + rtnValue.as_matrix().T

    np.fill_diagonal(rtnValue, 1)
    rtnValue = pd.DataFrame(rtnValue, index=contigs, columns=contigs)

    return rtnValue


def calculate_hypergeometric_survival(pc_df, shared_mtx):
    """
    Creates a symmetrical matrix containing the hypergeometric survival function (1-cdf)
    for each contig pair
    :param pc_df: The dataframe of PCs - to work out the number of PCs per contig
    :param shared_mtx: The matrix with the number of shared PCs between contigs
    :return: a symmetrical matrix of probabilities
    """
    pc_counts = pc_df.groupby(['contig']).PCid.nunique().reset_index()
    pc_counts.columns = ['contig', 'count']
    pc_counts.index = pc_counts['contig']
    contigs = list(shared_mtx.index)

    rtnValue = pd.DataFrame(0, index=contigs, columns=contigs)

    for i in itertools.combinations(contigs, 2):
        number_common_pcs = shared_mtx.loc[i[0], i[1]]
        a_pc_count = pc_counts.at[i[0], 'count']
        b_pc_count = pc_counts.at[i[1], 'count']
        a, b = sorted((a_pc_count, b_pc_count))
        total_pcs = pc_df.PCid.nunique()
        T = 0.5 * total_pcs * (total_pcs - 1)
        logT = np.log10(T)

        # The -1 is needed here to calculate the inverse cumlutative
        # density function
        pval = stats.hypergeom.sf(number_common_pcs - 1, total_pcs, a, b)
        sig = min(300, np.nan_to_num(-np.log10(pval) - logT))

        # If sig > 1, it is considered a match
        if sig <= 1:
            sig = 0
        rtnValue.at[i[0], i[1]] = sig

    return rtnValue
