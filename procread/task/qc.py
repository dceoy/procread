#!/usr/bin/env python

from itertools import product
import logging
import os
import re
from ..util.shell_operator import ShellOperator


def do_qc_checks(cf, cpus):
    logger = logging.getLogger(__name__)
    output_files = [
        os.path.join(
            cf['paths']['dir']['qc'],
            re.sub(r'\.fastq\.gz$', s, os.path.basename(d[r]))
        ) for d, r, s in product(
            cf['paths']['fastq'], ['read1', 'read2'],
            ['_fastqc.html', '_fastqc.zip']
        )
    ]
    logger.debug('output_files: {}'.format(output_files))
    if all([os.path.isfile(p) for p in output_files]):
        logger.debug('Skip quality control checks')
        return
    else:
        os.makedirs(cf['paths']['dir']['qc'], exist_ok=True)

    logger.info('Do quality control checks for reads')
    sh = ShellOperator(
        log_txt=os.path.join(cf['paths']['dir']['qc'], 'qc_log.txt')
    )
    sh.run('fastqc --version')

    for d, r in product(cf['paths']['fastq'], ['read1', 'read2']):
        sh.run(
            'fastqc {0} --threads {1} --outdir {2} {3}'.format(
                cf['cmd_args']['fastqc'], cpus, cf['paths']['dir']['qc'], d[r]
            ),
            cwd=cf['paths']['dir']['qc'],
        )
