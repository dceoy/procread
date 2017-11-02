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
    --file=<yaml>   Set a path to a YAML for configurations [$PROCREAD_YML]
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
import os
from docopt import docopt
from . import __version__
from .util import set_log_config, set_config_yml, write_config_yml, read_yaml
from .task import do_qc_checks, trim_adapters, map_reads, call_variants


def main():
    args = docopt(__doc__, version='procread {}'.format(__version__))
    set_log_config(debug=args['--debug'])
    logging.debug('args:{0}{1}'.format(os.linesep, args))
    config_yml = set_config_yml(path=args['--file'])

    if args['init']:
        logging.debug('Initiation')
        write_config_yml(path=config_yml)
    else:
        logging.debug('config_yml: {}'.format(config_yml))
        config = read_yaml(path=config_yml)
        logging.debug('config:{0}{1}'.format(os.linesep, config))
        if args['qc']:
            do_qc_checks(config=config)
        elif args['trim']:
            trim_adapters(config=config)
        elif args['map']:
            map_reads(config=config)
        elif args['call']:
            call_variants(config=config)
        elif args['run']:
            trim_adapters(config=config)
            map_reads(config=config)
            call_variants(config=config)
