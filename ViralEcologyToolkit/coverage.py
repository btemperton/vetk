import logging
import re
import subprocess

logger = logging.getLogger('VETK')


def run(parser, args):
    logger.info('Running coverage')

    bam_file = args.bam
    pct_cutoff = args.min_pct_cov


def extract_good_reads(bamfile, output, pct_id=90):
    bam_contents = subprocess.check_output(['samtools', 'view', bamfile])
    out_handle = open(output, 'w')
    for line in bam_contents:
        if is_good_enough_read(line, pct_id=pct_id):
            out_handle.write(line)
    out_handle.close()


def calculate_pct_identity(NM_str, MD_str):
    matches = 0.0
    mismatches = 0.0
    for m in re.findall('\d+', MD_str):
        matches += int(m)
    m = re.search('NM\:i\:(\d+)', NM_str)
    if m:
        mismatches = int(m.group(1))

    return (matches / (matches + mismatches)) * 100


def is_good_enough_read(line, pct_id):
    m = re.search('(NM\:i\:\d+)', line)
    if m:
        NM_str = m.group(1)
    else:
        return False

    m = re.search('(MD\:Z\:[A-Z0-9^]+)', line)
    if m:
        MD_str = m.group(1)
    else:
        return False

    return calculate_pct_identity(NM_str, MD_str) >= pct_id
