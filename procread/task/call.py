#!/usr/bin/env python

from itertools import product
import logging
import os
from ..util.shell_operator import ShellOperator


def call_variants(cf):
    logger = logging.getLogger(__name__)
    output_files = [
        os.path.join(cf['paths']['dir']['map'], d['id'], f)
        for d, f in product(cf['paths']['fastq'], ['raw.bcf', 'flt.vcf'])
    ]
    logger.debug('output_files: {}'.format(output_files))
    if all([os.path.isfile(p) for p in output_files]):
        logger.debug('Skip SNVs/indels calling')
        return
    else:
        os.makedirs(cf['paths']['dir']['call'], exist_ok=True)
        for d in cf['paths']['fastq']:
            os.makedirs(
                os.path.join(cf['paths']['dir']['call'], d['id']),
                exist_ok=True
            )

    logger.info('Call SNVs/indels')
    sh = ShellOperator(
        log_txt=os.path.join(
            cf['paths']['dir']['call'], 'call_log.txt'
        )
    )

    sh.run_parallel([
        'samtools mpileup -uf {0} {1} {2} | bcftools view -bvcg - > {3}'
        ' && bcftools view {3} | vcfutils.pl varFilter -D100 > {4}'.format(
            cf['paths']['ref']['fasta'],
            os.path.join(
                cf['paths']['dir']['map'], d['id'],
                'foreground_markdup.bam'
            ),
            os.path.join(
                cf['paths']['dir']['map'], d['id'],
                'background_markdup.bam'
            ),
            os.path.join(cf['paths']['dir']['call'], d['id'], 'raw.bcf'),
            os.path.join(cf['paths']['dir']['call'], d['id'], 'flt.vcf')
        )
        for d in cf['paths']['fastq']
    ])
