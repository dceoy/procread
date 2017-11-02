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
    qc              Do quality control checks
    trim            Trim adapter sequences in reads
    map             Map reads to a reference
    call            Call SNVs/indels
    run             Run a variant calling pipeline
"""

import logging
from multiprocessing import cpu_count
import os
from docopt import docopt
from . import __version__
from .util import dump_yaml, read_yaml, set_log_config, write_config_yml
from .task import do_qc_checks, trim_adapters, map_reads, call_variants


def main():
    args = docopt(__doc__, version='procread {}'.format(__version__))
    set_log_config(debug=args['--debug'])
    logging.debug('args:{0}{1}'.format(os.linesep, args))

    if args['init']:
        logging.debug('Initiation')
        write_config_yml(path=args['--file'])
    else:
        config = read_yaml(path=args['--file'])
        logging.debug('config:{0}{1}'.format(os.linesep, dump_yaml(config)))

        wd = args['--work']
        cpus = int(args['--cpus']) if args['--cpus'] else cpu_count()
        logging.debug('working dir: {}'.format(wd))
        os.makedirs(wd, exist_ok=True)

        if args['qc']:
            do_qc_checks(config=config, work_dir=wd, cpus=cpus)
        elif args['trim']:
            trim_adapters(config=config, work_dir=wd, cpus=cpus)
        elif args['map']:
            map_reads(config=config, work_dir=wd, cpus=cpus)
        elif args['call']:
            call_variants(config=config, work_dir=wd, cpus=cpus)
        elif args['run']:
            trim_adapters(config=config, work_dir=wd, cpus=cpus)
            map_reads(config=config, work_dir=wd, cpus=cpus)
            call_variants(config=config, work_dir=wd, cpus=cpus)
