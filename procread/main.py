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

import logging
from multiprocessing import cpu_count
import os
from docopt import docopt
from . import __version__
from .util import dump_yaml, read_yaml, set_log_config, write_config_yml
from .task import call_variants, do_qc_checks, map_reads, prepare_paths, \
                  trim_adapters


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

        cpus = int(args['--cpus']) if args['--cpus'] else cpu_count()
        paths = prepare_paths(
            config=config, work_dir=args['--work'], cpus=cpus
        )

        if args['qc']:
            do_qc_checks(config=config, paths=paths, cpus=cpus)
        elif args['trim']:
            trim_adapters(config=config, paths=paths, cpus=cpus)
        elif args['map']:
            map_reads(config=config, paths=paths, cpus=cpus)
        elif args['call']:
            call_variants(config=config, paths=paths, cpus=cpus)
        elif args['run']:
            trim_adapters(config=config, paths=paths, cpus=cpus)
            map_reads(config=config, paths=paths, cpus=cpus)
            call_variants(config=config, paths=paths, cpus=cpus)
