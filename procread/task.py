#!/usr/bin/env python

from itertools import product
import logging
import os
import re
from .util import ProcreadRuntimeError, Shell


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
    sh = Shell(
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
                    raise ProcreadRuntimeError(
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
            raise ProcreadRuntimeError(
                'Supported extension for input reference fasta: '
                '.fa, .fa.gz, .fa.bz2, .fasta, .fasta.gz, .fasta.bz'
            )


def do_qc_checks(cf, cpus):
    logger = logging.getLogger(__name__)
    output_files = [
        os.path.join(
            cf['paths']['dir']['qc'],
            re.sub(r'\.fastq\.gz$', s, os.path.basename(d[r]))
        )
        for d, r, s
        in product(
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
    sh = Shell(log_txt=os.path.join(cf['paths']['dir']['qc'], 'qc_log.txt'))
    sh.run('fastqc --version')

    for d, r in product(cf['paths']['fastq'], ['read1', 'read2']):
        sh.run(
            'fastqc {0} --threads {1} --outdir {2} {3}'.format(
                cf['cmd_args']['fastqc'], cpus, cf['paths']['dir']['qc'], d[r]
            ),
            cwd=cf['paths']['dir']['qc'],
        )


def trim_adapters(cf, cpus):
    logger = logging.getLogger(__name__)
    output_files = [
        os.path.join(
            cf['paths']['dir']['trim'],
            re.sub(
                r'(\.r[12]\.fastq\.gz)$', r'.trimmed\1',
                os.path.basename(d['{0}_{1}'.format(t, r)])
            )
        )
        for d, t, r in product(cf['paths']['fastq'], ['read1', 'read2'])
    ]
    logger.debug('output_files: {}'.format(output_files))
    if all([os.path.isfile(p) for p in output_files]):
        logger.debug('Skip adapter trimming')
        return
    else:
        os.makedirs(cf['paths']['dir']['trim'], exist_ok=True)

    logger.info('Trim adapter sequences in reads')
    sh = Shell(log_txt=os.path.join(
       cf['paths']['dir']['trim'], 'trim_log.txt'
    ))
    sh.run('cutadapt --version')

    io = [
        {
            'out': {
                r: os.path.join(
                    cf['paths']['dir']['trim'],
                    re.sub(
                        r'(\.r[12]\.fastq\.gz)$', r'.trimmed\1',
                        os.path.basename(d[r])
                    )
                )
                for r in ['read1', 'read2']
            },
            'in': {
                r: d[r] for r in ['read1', 'read2']
            }
        }
        for d in cf['paths']['fastq']
    ]

    sh.run([
        'cutadapt {0} -a {1} -A {2} -o {3} -p {4} {5} {6}'.format(
            cf['cmd_args']['cutadapt'], cpus, cf['yml']['adapter']['forward'],
            cf['yml']['adapter']['reverse'], f['out']['read1'],
            f['out']['read2'], f['in']['read1'], f['in']['read2']
        )
        for f in io
        if not all([os.path.isfile(p) for p in f['out'].values()])
    ])


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
    sh = Shell(log_txt=os.path.join(
        cf['paths']['dir']['input'], 'ref_log.txt'
    ))

    sh.run(['samtools --version', 'bwa'], check=False)
    sh.run([
        'samtools faidx {}'.format(cf['paths']['ref']['fasta']),
        'bwa index -p {0} {1}'.format(
            cf['paths']['ref']['faidx'], cf['paths']['ref']['fasta']
        )
    ])


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
    sh = Shell(log_txt=os.path.join(
        cf['paths']['dir']['map'], 'map_log.txt'
    ))

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
    sh = Shell(log_txt=os.path.join(
        cf['paths']['dir']['call'], 'call_log.txt'
    ))

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
