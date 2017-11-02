#!/usr/bin/env python

import logging
import os
import re
from itertools import product
from .util import sh_c


def do_qc_checks(config, work_dir, cpus):
    logging.info('Start quality control checks for fastq files')
    qc_dir = os.path.join(work_dir, 'qc')
    logging.debug('qc_dir: {}'.format(qc_dir))
    os.makedirs(qc_dir, exist_ok=True)
    for t, r in product(['foreground', 'background'], ['read1', 'read2']):
        fq = config['fastq'][t][r]
        logging.debug('fq: {}'.format(fq))
        qc_fq_dir = os.path.join(
            qc_dir, re.sub(r'\.fastq\.gz$', '', os.path.basename(fq))
        )
        logging.debug('qc_fq_dir: {}'.format(qc_fq_dir))
        os.makedirs(qc_fq_dir, exist_ok=True)
        sh_c('fastqc --threads {0} --nogroup --outdir {1} {2} > {3}'.format(
            cpus, qc_fq_dir, fq, os.path.join(qc_fq_dir, 'fastqc.log')
        ))


def trim_adapters(config, work_dir):
    pass


def map_reads(config, work_dir):
    pass


def call_variants(config, work_dir):
    pass
