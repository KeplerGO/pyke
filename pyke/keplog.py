import logging

_logger = logging.getLogger('pyke')
_logger.setLevel(logging.INFO)

formatter = logging.Formatter('%(asctime)s %(levelname)s - %(message)s',
                              '%Y-%m-%d %H:%M:%S')

console_handler = logging.StreamHandler()
console_handler.setFormatter(formatter)

_logger.addHandler(console_handler)


def debug(*args, **kwargs):
    _logger.debug(*args, **kwargs)


def info(*args, **kwargs):
    _logger.info(*args, **kwargs)


def warn(*args, **kwargs):
    _logger.warn(*args, **kwargs)


def error(*args, **kwargs):
    _logger.error(*args, **kwargs)


def exception(*args, **kwargs):
    _logger.exception(*args, **kwargs)


def set_level(level):
    try:
        _logger.setLevel(getattr(logging, level))
    except TypeError:
        _logger.setLevel(level)


def set_logfile(log_filepath):
    fh = logging.FileHandler(log_filepath)
    fh.setFormatter(formatter)
    _logger.addHandler(fh)


def setup_args(parser):
    g = parser.add_mutually_exclusive_group()
    g.add_argument('--silent',
                   action='store_true',
                   help='Set log level to ERROR or above')
    g.add_argument('--debug',
                   action='store_true',
                   help='Set log level to DEBUG or above')


def handle_args(args):
    if args.silent:
        _logger.setLevel(logging.ERROR)
    if args.debug:
        _logger.setLevel(logging.DEBUG)

    if args.logfile is not None:
        set_logfile(args.logfile)
