#!/usr/bin/env python

import logging
import os
from ..util.shell_operator import ShellOperator


def make_ref_index(cf):
    logger = logging.getLogger(__name__)
    output_files = [
        cf['paths']['ref']['faidx'] + s
        for s in ['', '.amb', '.ann', '.pac', '.bwt', '.sa']
    ]
    logger.debug('output_files: {}'.format(output_files))
    if all([os.path.isfile(p) for p in output_files]):
        logger.debug('Skip making reference index files')
        return
    logger.info('Make reference index files')
    sh = ShellOperator(
        log_txt=os.path.join(
            cf['paths']['dir']['input'], 'ref_log.txt'
        )
    )
    sh.run(['samtools --version', 'bwa'], check=False)
    sh.run([
        'samtools faidx {}'.format(cf['paths']['ref']['fasta']),
        'bwa index -p {0} {1}'.format(
            cf['paths']['ref']['faidx'], cf['paths']['ref']['fasta']
        )
    ])
