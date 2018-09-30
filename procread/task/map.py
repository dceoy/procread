#!/usr/bin/env python

from itertools import product
import logging
import os
import re
from ..util.shell_operator import ShellOperator


def map_reads(cf, cpus):
    logger = logging.getLogger(__name__)
    output_files = [
        os.path.join(cf['paths']['dir']['map'], d['id'], f)
        for d, f
        in product(
            cf['paths']['fastq'],
            [
                '{0}.bam{1}'.format(n, e)
                for n, e
                in product(
                    ['sort', 'fixmate', 'markdup'],
                    ['', '.bai', '.idxstats.txt', '.flagstat.txt',
                     '.stats.txt']
                )
            ]
        )
    ]
    logger.debug('output_files: {}'.format(output_files))
    if all([os.path.isfile(p) for p in output_files]):
        logger.debug('Skip mapping read')
        return
    else:
        os.makedirs(cf['paths']['dir']['map'], exist_ok=True)
        for d in cf['paths']['fastq']:
            os.makedirs(
                os.path.join(cf['paths']['dir']['map'], d['id']),
                exist_ok=True
            )

    logger.info('Map reads to a reference')
    sh = ShellOperator(
        log_txt=os.path.join(
            cf['paths']['dir']['map'], 'map_log.txt'
        )
    )

    io = [
        {
            'id': d['id'],
            'fq': (
                {
                    k: os.path.join(
                        cf['paths']['dir']['trim'],
                        re.sub(
                            r'(\.r[12]\.fastq\.gz)$', r'.trimmed\1',
                            os.path.basename(v)
                        )
                    )
                    for k, v in d.items()
                }
                if os.path.isdir(cf['paths']['dir']['trim']) else d
            ),
            'bam': {
                '{0}.bam{1}'.format(n, e):
                os.path.join(
                    cf['paths']['dir']['map'], d['id'],
                    '{0}.bam{1}'.format(n, e)
                )
                for n, e
                in product(
                    ['sort', 'fixmate', 'markdup'],
                    ['', '.bai', '.idxstats.txt', '.flagstat.txt',
                     '.stats.txt']
                )
            }
        }
        for d in cf['paths']['fastq']
    ]

    for d in io:
        if not os.path.isfile(d['bam']['sort.bam']):
            sh.run([
                'bwa mem -t {0} {1} {2} {3} '
                '| samtools view -@ {0} -bS - '
                '| samtools sort -@ {0} -o {4} -'.format(
                    cpus, cf['paths']['ref']['faidx'], d['fq']['read1'],
                    d['fq']['read2'], d['bam']['sort.bam']
                )
            ])
        if not os.path.isfile(d['bam']['fixmate.bam']):
            sh.run([
                'samtools sort -n -@ {0} {1} '
                '| samtools fixmate -cm - {2}'.format(
                    cpus, d['bam']['sort.bam'], d['bam']['fixmate.bam']
                )
            ])
        if not os.path.isfile(d['bam']['markdup.bam']):
            sh.run([
                'samtools sort -@ {0} {1} '
                '| samtools markdup -rs - {2}'.format(
                    cpus, d['bam']['fixmate.bam'], d['bam']['markdup.bam']
                )
            ])

    sh.run([
        'samtools index -@ {0} {1}'.format(
            cpus, d['bam']['{}.bam'.format(b)]
        )
        for d, b
        in product(io, ['sort', 'markdup'])
        if not os.path.isfile(d['bam']['{}.bam.bai'.format(b)])
    ])

    sh.run_parallel([
        'samtools flagstat {0} > {1}'.format(
            d['bam']['{}.bam'.format(b)],
            d['bam']['{0}.bam.{1}.txt'.format(b, c)]
        )
        for d, b, c
        in product(
            io, ['sort', 'fixmate', 'markdup'],
            ['idxstats', 'flagstat', 'stats']
        )
        if not os.path.isfile(d['bam']['{0}.bam.{1}.txt'.format(b, c)])
    ])
