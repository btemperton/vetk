import argparse
import logging

import ViralEcologyToolkit


def main(args):
    logging.basicConfig()
    logger = logging.getLogger('extract.signal.peptides.py.logger')

    s = ViralEcologyToolkit.parse_mcl_dump_file(args.mcl_dump)
    s.to_csv('~/Downloads/s.csv')
    t = ViralEcologyToolkit.load_gene_map(args.gene_map)
    t.to_csv('~/Downloads/t.csv')
    u = ViralEcologyToolkit.add_pc_labels(s, t)
    u.to_csv('~/Downloads/u.csv')
    v = ViralEcologyToolkit.create_composition_mtx(u, presence_absence=True)
    v.to_csv('~/Downloads/v.csv')
    w = ViralEcologyToolkit.calculate_shared_matrix(v)
    w.to_csv('~/Downloads/w.csv')
    x = ViralEcologyToolkit.calculate_hypergeometric_survival(u, w)
    x.to_csv('~/Downloads/x.csv')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mcl_dump', type=str, required=True)
    parser.add_argument('--gene_map', type=str, required=True)
    args = parser.parse_args()
    main(args)
