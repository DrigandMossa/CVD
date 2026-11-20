"""Command line: argument parsing, the text wizard and the entry point."""

from .args import build_parser, config_from_args, print_catalog, read_seq_arg
from .run import main
from .wizard import wizard

__all__ = ["build_parser", "config_from_args", "main", "print_catalog", "read_seq_arg",
           "wizard"]
