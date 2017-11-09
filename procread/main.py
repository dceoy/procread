#!/usr/bin/env python
"""
Run read-to-variant pipelines for DNA-seq analyses

Usage:
    procread init [--debug] [--file=<yaml>] [--work=<dir>]
    procread qc [--debug] [--file=<yaml>] [--cpus=<int>] [--work=<dir>]
    procread trim [--debug] [--file=<yaml>] [--cpus=<int>] [--work=<dir>]
    procread ref [--debug] [--file=<yaml>] [--cpus=<int>] [--work=<dir>]
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
    ref             Prepare reference index files
    map             Map reads to a reference
    call            Call SNVs/indels
    run             Run a variant calling pipeline
"""

import logging
from multiprocessing import cpu_count
import os
from docopt import docopt
import yaml
from . import __version__
from .util import generate_param_config, set_log_config, write_config_yml
from .task import call_variants, do_qc_checks, map_reads, make_ref_index, \
                  prepare_paths, trim_adapters


def main():
    args = docopt(__doc__, version='procread {}'.format(__version__))
    set_log_config(debug=args['--debug'])
    logger = logging.getLogger(__name__)
    logger.debug('args:{0}{1}'.format(os.linesep, args))

    if args['init']:
        logger.debug('Initiation')
        write_config_yml(path=args['--file'])
    else:
        cf = generate_param_config(
            yml_path=args['--file'], work_dir=args['--work']
        )
        logger.debug('cf:{0}{1}'.format(
            os.linesep, yaml.dump(cf, default_flow_style=False)
        ))
        cpus = int(args['--cpus']) if args['--cpus'] else cpu_count()

        prepare_paths(cf=cf, cpus=cpus)
        if args['qc']:
            do_qc_checks(cf=cf, cpus=cpus)
        elif args['trim']:
            trim_adapters(cf=cf)
        elif args['ref']:
            make_ref_index(cf=cf)
        elif args['map']:
            map_reads(cf=cf, cpus=cpus)
        elif args['call']:
            call_variants(cf=cf)
        elif args['run']:
            do_qc_checks(cf=cf, cpus=cpus)
            trim_adapters(cf=cf, cpus=cpus)
            make_ref_index(cf=cf)
            map_reads(cf=cf, cpus=cpus)
            call_variants(cf=cf)
