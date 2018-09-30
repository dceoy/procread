#!/usr/bin/env python

from itertools import product
import logging
import os
from ..util.shell_operator import ShellOperator


def prepare_paths(cf, cpus):
    logger = logging.getLogger(__name__)
    output_files = [
        d[r]
        for d, r
        in product(cf['paths']['fastq'], ['read1', 'read2'])
    ] + [cf['paths']['ref']['fasta']]
    logger.debug('output_files: {}'.format(output_files))
    if all([os.path.isfile(p) for p in output_files]):
        logger.debug('Skip file and directory preparation')
        return
    else:
        for d in ['work', 'input']:
            os.makedirs(cf['paths']['dir'][d], exist_ok=True)

    logger.info('Prepare paths to directories and files')
    sh = ShellOperator(
        log_txt=os.path.join(cf['paths']['dir']['input'], 'prep_log.txt')
    )
    for c in ['pigz', 'pbzip2']:
        sh.run('{} --version'.format(c))

    for src, dst in zip(cf['yml']['path']['fastq'], cf['paths']['fastq']):
        for r in ['read1', 'read2']:
            fq_dst = dst[r]
            if not os.path.isfile(fq_dst):
                fq_src = os.path.abspath(src[r])
                fq_src_ext = os.path.splitext(fq_src)[1]
                if fq_src_ext in ['.fastq', '.fq']:
                    sh.run(
                        'pigz -p {0} -c {1} > {2}'.format(cpus, fq_src, fq_dst)
                    )
                elif fq_src_ext == '.gz':
                    os.symlink(fq_src, fq_dst)
                elif fq_src_ext == '.bz2':
                    sh.run(
                        'pbzip2 -p {0} -dc {1} '
                        '| pigz -p {0} -c - > {2}'.format(cpus, fq_src, fq_dst)
                    )
                else:
                    raise RuntimeError(
                        'Supported extension for input fastq: '
                        '.fastq, .fastq.gz, .fastq.*.gz, '
                        '.fastq.bz2, .fastq.*.bz2, '
                        '.fq, .fq.gz, .fq.*.gz, .fq.bz2, .fq.*.bz2'
                    )

    ref_dst = cf['paths']['ref']['fasta']
    ref_src = os.path.abspath(cf['yml']['path']['fasta']['reference'])
    ref_src_ext = os.path.splitext(ref_src)[1]
    if not os.path.isfile(ref_dst):
        if ref_src_ext in ['.fasta', '.fa']:
            os.symlink(ref_src, ref_dst)
        elif ref_src_ext == '.gz':
            sh.run('pigz -p {0} -dc {1} > {2}'.format(cpus, ref_src, ref_dst))
        elif ref_src_ext == '.bz2':
            sh.run(
                'pbzip2 -p {0} -dc {1} > {2}'.format(cpus, ref_src, ref_dst)
            )
        else:
            raise RuntimeError(
                'Supported extension for input reference fasta: '
                '.fa, .fa.gz, .fa.bz2, .fasta, .fasta.gz, .fasta.bz'
            )
