from argparse import HelpFormatter, SUPPRESS, OPTIONAL, ZERO_OR_MORE

class PyKEArgumentHelpFormatter(HelpFormatter):
    """Help message formatter which adds default values to argument help,
    except for boolean arguments.
    """

    def _get_help_string(self, action):
        help = action.help
        if '%(default)' not in action.help:
            if (action.default is not SUPPRESS
                and not isinstance(action.default, bool)):
                defaulting_nargs = [OPTIONAL, ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help += ' (default: %(default)s)'
        return help
