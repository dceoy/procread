#!/usr/bin/env python
"""
Run read-to-variant pipelines for DNA-seq analyses

Usage:
    procread init [--debug|--info] [--file=<yml>] [--work=<dir>]
    procread qc [--debug|--info] [--file=<yml>] [--cpus=<int>] [--work=<dir>]
    procread trim [--debug|--info] [--file=<yml>] [--cpus=<int>] [--work=<dir>]
    procread ref [--debug|--info] [--file=<yml>] [--cpus=<int>] [--work=<dir>]
    procread map [--debug|--info] [--file=<yml>] [--cpus=<int>] [--work=<dir>]
    procread call [--debug|--info] [--file=<yml>] [--cpus=<int>] [--work=<dir>]
    procread run [--debug|--info] [--file=<yml>] [--cpus=<int>] [--work=<dir>]
    procread -h|--help
    procread -v|--version

Options:
    -h, --help          Print help and exit
    -v, --version       Print version and exit
    --debug, --info     Execute a command with debug messages
    --file=<yml>        Set a configuration YAML path [default: ./pread.yml]
    --cpus=<int>        Limit CPU cores for use
    --work=<dir>        Set a working directory [default: .]

Commands:
    init                Generate a YAML template for configuration
    qc                  Do quality control checks for reads
    trim                Trim adapter sequences in reads
    ref                 Prepare reference index files
    map                 Map reads to a reference
    call                Call SNVs/indels
    run                 Run a variant calling pipeline
"""

import logging
from multiprocessing import cpu_count
import os
from pprint import pformat
from docopt import docopt
from .. import __version__
from ..util.config import load_param_config, set_log_config, write_config_yml
from ..util.prep import prepare_paths
from ..task.call import call_variants
from ..task.map import map_reads
from ..task.qc import do_qc_checks
from ..task.ref import make_ref_index
from ..task.trim import trim_adapters


def main():
    args = docopt(__doc__, version='procread {}'.format(__version__))
    set_log_config(debug=args['--debug'])
    logger = logging.getLogger(__name__)
    logger.debug('args:' + os.linesep + pformat(args))
    if args['init']:
        logger.debug('Initiation')
        write_config_yml(path=args['--file'])
    else:
        cf = load_param_config(
            yml_path=args['--file'], work_dir=args['--work']
        )
        logger.debug('cf:' + os.linesep + pformat(cf))
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
