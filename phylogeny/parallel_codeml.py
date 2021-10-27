#!/usr/bin/env python3
"""Run codeml in parallel."""

import argparse
import os
import re
import sys
from multiprocessing import Process, cpu_count
from subprocess import Popen, PIPE
from Bio.Phylo.PAML import codeml

def generate_ctl(alignment, tree, model, outfile):
    """
    Generates a Codeml object with options for the given model. Accepts "site" 
    model (model=0, NSsites="0 1 2"), "branch_site_M2a" (model=2, NSsites=2), and
    "branch_site_null" (model=2, NSsites=2, fix_omega=1, omega=1).
    """
    working_dir = re.search("(.+).paml", os.path.basename(alignment))[1]
    out_file = "mlc" if outfile is None else outfile
    
    # Initialize Codeml object
    cml = codeml.Codeml(
        alignment=alignment,
        tree=tree,
        out_file=f"{working_dir}/{out_file}",
        working_dir=f"./{working_dir}"
        )

    # Set constant params
    cml.set_options(
        noisy=9,
        verbose=0,
        runmode=0,
        seqtype=1,
        CodonFreq=2,
        clock=0,
        aaDist=0,
        fix_alpha=1,
        alpha=0.0,
        Malpha=0,
        aaRatefile="dat/jones.dat",
        icode=0,
        fix_kappa=0,
        kappa=2,
        getSE=0,
        ncatG=8,
        RateAncestor=1,
        Small_Diff=.5e-6,
        cleandata=1,
        fix_blength=None,
        ndata=None
        )

    # Set model-specific params
    if model == "site":
        cml.set_options(
            model=0,
            NSsites=[0, 1, 2]
        )
    elif model == "branch_site_null":
        cml.set_options(
            model=2,
            NSsites=[2],
            fix_omega=1,
            omega=1
        )
    elif model == "branch_site_M2a":
        cml.set_options(
            model=2,
            NSsites=[2]
        )
    else:
        raise ValueError("model must be 'site', 'branch_site_null' or 'branch_site_M2a'")
    
    return cml


def run_codeml(alignment_files, tree, model, outfile):
    """Run codeml sequentially for all alignment_files."""
    for alignment in alignment_files:
        cml = generate_ctl(alignment, tree, model, outfile=outfile)
        cml.run(verbose=True, parse=False)


def launch_run_codeml(num_processes, partitioned_alignments, tree, model, outfile):
    """Launch num_processes parallel processes for codeml."""
    processes = []
    
    for i in range(num_processes):
        p = Process(target=run_codeml, args=(partitioned_alignments[i], tree, model, outfile))
        processes.append(p)
        p.start()
    
    for p in processes:
        p.join()
    

def get_alignment_files(filename, num_processes):
    """Returns a partitioned list of alignment files. Length of list = num_processes."""
    partitioned_alignments = []
    
    for i in range(num_processes):
        partitioned_alignments.append([])
    
    with open(filename, "r") as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            partition = i % num_processes
            partitioned_alignments[partition].append(line)

    return partitioned_alignments    


def parse_args():
    parser = argparse.ArgumentParser(description="Run codeml in in parallel")
    parser.add_argument("alignments",
                        type=str,
                        help="Line-separated list of alignment files")
    parser.add_argument("tree",
                        type=str,
                        help="Tree file for codeml")
    parser.add_argument("model",
                        choices=["site", "branch_site_null", "branch_site_M2a"],
                        help="Codeml model")
    parser.add_argument("-n", "--num_processes",
                        type=int,
                        default=8,
                        help="Number of parallel processes [8]")
    parser.add_argument("-o", "--outfile",
                        type=str,
                        default=None,
                        help="Output file name for codeml runs (each run will produce an output file")
    
    return parser.parse_args()


def main():
    args = parse_args()
    num_cpus = cpu_count()
    if args.num_processes > num_cpus:
        print(f"Not enough CPUs available. Scaling down to {num_cpus} processes")
        args.num_processes = num_cpus
    
    p = Popen("which codeml", shell=True, stdout=PIPE)
    path = p.communicate()[0].decode("utf-8")
    if len(path) == 0:
        print("Error: Cannot locate codeml on PATH")
        sys.exit(1)
    
    print(f"Found codeml executable: {path}")

    partitioned_alignments = get_alignment_files(args.alignments, args.num_processes)
    launch_run_codeml(args.num_processes, partitioned_alignments, args.tree, args.model, args.outfile)


if __name__ == "__main__":
    main()