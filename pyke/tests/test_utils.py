import pytest
import argparse

from ..utils import PyKEArgumentHelpFormatter

def test_PyKEArgumentHelpFormatter():
    parser = argparse.ArgumentParser(
                description=('Test PyKEArgumentHelpFormatter'),
                formatter_class=PyKEArgumentHelpFormatter)
    parser.add_argument('--str_arg', help='string type argument',
                        default='oi')
    parser.add_argument('--bool_arg', help='bool type argument',
                        action='store_true')

    formatter = parser._get_formatter()
    for actions in parser._action_groups:
        formatter.add_arguments(actions._group_actions)

    ans = ["-h, --help         show this help message and exit",
           "--str_arg STR_ARG  string type argument (default: oi)",
           "--bool_arg         bool type argument", '']
    assert ans == formatter.format_help().split('\n')

