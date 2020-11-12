#!/usr/bin/env python3

import sys
import argparse

def get_molecules(bedfile):
    """Read molecule lengths from a bed file."""
    mol_len = []
    with open(bedfile) as bed:

