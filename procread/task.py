#!/usr/bin/env python

import logging
import os
from itertools import product
from .util import sh_c


def do_qc_checks(config, work_dir, cpus):
    logging.info('Start quality control checks for fastq files')

    qc_dir = os.path.join(work_dir, 'qc')
    os.makedirs(qc_dir, exist_ok=True)
    logging.debug('qc_dir: {}'.format(qc_dir))
    qc_log_txt = os.path.join(qc_dir, 'fastqc_log.txt')

    for t, r in product(['foreground', 'background'], ['read1', 'read2']):
        cmd = (
            'fastqc --nogroup --threads {0} --outdir {1} {2}'
            ' 2>&1 | tee -a {3}'
        ).format(
            cpus, qc_dir, config['fastq'][t][r], qc_log_txt
        )
        with open(qc_log_txt, 'a') as f:
            f.write('$ {0}{1}'.format(cmd, os.linesep))
        sh_c(cmd)


def trim_adapters(config, work_dir):
    pass


def map_reads(config, work_dir):
    pass


def call_variants(config, work_dir):
    pass
