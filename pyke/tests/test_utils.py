import pytest
import argparse

from ..utils import PyKEArgumentHelpFormatter
from ..utils import module_output_to_channel, channel_to_module_output
from ..utils import KeplerQualityFlags


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


def test_channel_to_module_output():
    assert channel_to_module_output(1) == (2, 1)
    assert channel_to_module_output(42) == (13, 2)
    assert channel_to_module_output(84) == (24, 4)
    assert channel_to_module_output(33) == (11, 1)


def test_module_output_to_channel():
    assert module_output_to_channel(2, 1) == 1
    assert module_output_to_channel(13, 2) == 42
    assert module_output_to_channel(24, 4) == 84
    assert module_output_to_channel(11, 1) == 33


def test_parse_kepler_quality_flags():
    flags = list(KeplerQualityFlags.flags.items())
    # Can we recover each individual quality flag?
    for key, value in flags:
        assert KeplerQualityFlags.parse(key)[0] == value
    # Can we recover combinations of flags?
    assert KeplerQualityFlags.parse(flags[5][0] + flags[7][0]) == [flags[5][1], flags[7][1]]
    assert KeplerQualityFlags.parse(flags[3][0] + flags[4][0] + flags[5][0]) == [flags[3][1], flags[4][1], flags[5][1]]
