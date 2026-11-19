"""Reusable argparse type converters with range and finiteness checks."""

import argparse
from math import isfinite


def positive_int(value: str) -> int:
    """Parse an integer greater than or equal to one."""
    try:
        number = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("muss eine ganze Zahl sein") from error
    if number < 1:
        raise argparse.ArgumentTypeError("muss mindestens 1 sein")
    return number


def non_negative_int(value: str) -> int:
    """Parse an integer greater than or equal to zero."""
    try:
        number = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("muss eine ganze Zahl sein") from error
    if number < 0:
        raise argparse.ArgumentTypeError("muss mindestens 0 sein")
    return number


def positive_float(value: str) -> float:
    """Parse a finite float greater than zero."""
    try:
        number = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("muss eine positive Zahl sein") from error
    if not isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("muss endlich und größer als 0 sein")
    return number


def non_negative_float(value: str) -> float:
    """Parse a finite float greater than or equal to zero."""
    try:
        number = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("muss eine nicht negative Zahl sein") from error
    if not isfinite(number) or number < 0:
        raise argparse.ArgumentTypeError("muss endlich und mindestens 0 sein")
    return number