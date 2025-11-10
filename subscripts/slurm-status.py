#!/usr/bin/env python3
# slurm-status.py
# Snakemake cluster-status script for SLURM
#
# Reports one of: "success", "failed", "running", "pending"

import sys
import subprocess

jobid = sys.argv[1]

def run_cmd(cmd):
    try:
        res = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        return res.stdout.decode().strip()
    except subprocess.CalledProcessError:
        return ""

# First try sacct (more reliable for finished jobs)
sacct_out = run_cmd(f"sacct -P -b -j {jobid} --format=State")
if sacct_out:
    # sacct output looks like:
    # State
    # COMPLETED
    lines = sacct_out.splitlines()
    if len(lines) > 1:
        state = lines[1].split()[0]  # second line, first word
        if state in ["COMPLETED"]:
            print("success")
            sys.exit(0)
        elif state in ["FAILED", "CANCELLED", "TIMEOUT", "NODE_FAIL", "OUT_OF_MEMORY"]:
            print("failed")
            sys.exit(0)
        elif state in ["PENDING"]:
            print("pending")
            sys.exit(0)
        else:
            print("running")
            sys.exit(0)

# If sacct didn’t return info, fall back to squeue (for queued/running jobs)
squeue_out = run_cmd(f"squeue -j {jobid} -h -o %T")
if squeue_out:
    state = squeue_out.strip()
    if state in ["PENDING", "CONFIGURING"]:
        print("running")
    elif state in ["RUNNING", "COMPLETING"]:
        print("running")
    else:
        print("failed")
else:
    # If job is not found in sacct or squeue, assume finished successfully
    # (depending on site policy, you may want "failed" instead)
    print("success")