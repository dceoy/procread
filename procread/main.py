#!/usr/bin/env python
"""
Run read-to-variant pipelines for DNA-seq analyses

Usage:
    procread init [--debug] [--file=<yaml>] [--work=<dir>]
    procread qc [--debug] [--file=<yaml>] [--cpus=<int>] [--work=<dir>]
    procread trim [--debug] [--file=<yaml>] [--cpus=<int>] [--work=<dir>]
    procread map [--debug] [--file=<yaml>] [--cpus=<int>] [--work=<dir>]
    procread call [--debug] [--file=<yaml>] [--cpus=<int>] [--work=<dir>]
    procread run [--debug] [--file=<yaml>] [--cpus=<int>] [--work=<dir>]
    procread -h|--help
    procread -v|--version

Options:
    -h, --help      Print help and exit
    -v, --version   Print version and exit
    --debug         Execute a command with debug messages
    --file=<yaml>   Set a path to a configuration YAML [default: ./pread.yml]
    --cpus=<int>    Limit CPU cores for use
    --work=<dir>    Set a working directory [default: .]

Commands:
    init            Generate a YAML template for configuration
    qc              Do quality control checks for reads
    trim            Trim adapter sequences in reads
    map             Map reads to a reference
    call            Call SNVs/indels
    run             Run a variant calling pipeline
"""

from concurrent.futures import ProcessPoolExecutor
import logging
from multiprocessing import cpu_count
import os
from docopt import docopt
from . import __version__
from .util import dump_yaml, generate_config, set_log_config, write_config_yml
from .task import call_variants, do_qc_checks, map_reads, make_ref_index, \
                  prepare_paths, trim_adapters


def main():
    args = docopt(__doc__, version='procread {}'.format(__version__))
    set_log_config(debug=args['--debug'])
    logging.debug('args:{0}{1}'.format(os.linesep, args))

    if args['init']:
        logging.debug('Initiation')
        write_config_yml(path=args['--file'])
    else:
        cf = generate_config(yml_path=args['--file'], work_dir=args['--work'])
        logging.debug('cf:{0}{1}'.format(os.linesep, dump_yaml(cf)))
        cpus = int(args['--cpus']) if args['--cpus'] else cpu_count()

        prepare_paths(config=cf, cpus=cpus)
        if args['qc']:
            do_qc_checks(config=cf, cpus=cpus)
        elif args['trim']:
            trim_adapters(config=cf)
        elif args['map']:
            make_ref_index(config=cf)
            map_reads(config=cf, cpus=cpus)
        elif args['call']:
            call_variants(config=cf)
        elif args['run']:
            with ProcessPoolExecutor(max_workers=cpus) as ppe:
                ppe.submit(
                    do_qc_checks, config=cf, cpus=(cpus - 2 if cpus > 2 else 1)
                )
                ppe.submit(trim_adapters, config=cf)
                ppe.submit(make_ref_index, config=cf)
            map_reads(config=cf, cpus=cpus)
            call_variants(config=cf)
