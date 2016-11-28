import argparse

import logging

logger = logging.getLogger('vetk')


def run_subtool(parser, args):
    if args.command == 'reticulator':
        import ViralEcologyToolkit.reticulator as submodule
    elif args.command == 'coverage':
        import ViralEcologyToolkit.coverage as submodule

    submodule.run(parser, args)


class ArgumentParserWithDefaults(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
        self.add_argument("-q", "--quiet", help="Do not output warnings to stderr",
                          action="store_true",
                          dest="quiet")


def main():
    logging.basicConfig()
    parser = argparse.ArgumentParser(prog='vetk', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command', parser_class=ArgumentParserWithDefaults)

    parser_coverage = subparsers.add_parser('coverage')
    parser_coverage.add_argument('--bam_file', dest='bam_file', metavar='STRING', required=True)
    parser_coverage.add_argument('--out_file', dest='out_file', metavar='STRING', required=True)
    parser_coverage.add_argument('-pct_id', dest='pct_id', type=float, default=90.0)

    parser_coverage.set_defaults(func=run_subtool)

    parser_reticulator = subparsers.add_parser('reticulator')
    parser_reticulator.add_argument('--mcl_dump', dest='mcl_dump', metavar='STRING', required=True)
    parser_reticulator.add_argument('--gene_map', dest='gene_map', metavar='STRING', required=True)
    parser_reticulator.add_argument('--output', dest='output', metavar='STRING', required=True)
    parser_reticulator.set_defaults(func=run_subtool)

    args = parser.parse_args()

    if args.quiet:
        logger.setLevel(logging.ERROR)

    try:
        args.func(parser, args)
    except IOError as e:
        if e.errno != 32:  # ignore SIGPIPE
            raise


if __name__ == '__main__':
    main()
