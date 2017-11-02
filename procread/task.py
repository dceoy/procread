#!/usr/bin/env python

import os
from itertools import product
from .util import bash


def do_qc_checks(config, work_dir, cpus):
    qc_dir = os.path.join(work_dir, 'qc')
    os.makedirs(qc_dir, exist_ok=True)
    for t, r in product(['foreground', 'background'], ['read1', 'read2']):
        fq = config['fastq'][t][r]
        qc_fq_dir = os.path.join(
            qc_dir, os.path.splitext(os.path.splitext(fq)[0])[0]
        )
        bash('fastqc --threads {0} --nogroup --outdir {1} {2} > {3}'.format(
            cpus, qc_fq_dir, fq, os.path.join(qc_fq_dir, 'fastqc.log')
        ))


def trim_adapters(config, work_dir):
    pass


def map_reads(config, work_dir):
    pass


def call_variants(config, work_dir):
    pass
