#!/usr/bin/env python3
"""chrM haplogroup calling for the KZ-219 cohort.

BAM -> chrM reads -> remap to rCRS (NC_012920.1) -> mutserve -> haplogrep3,
aggregated into a single table consumed by 05_chrM_analysis.ipynb
(KZ_CHRM_DIR/kz_219_mtdna_haplogroups.tsv).

Versions used for the article: samtools v1.21, BWA-MEM v0.7.12-r1039,
mutserve v2.0.1, HaploGrep3 v3.2.2 with the phylotree-fu-rcrs@1.2 tree.
Resumable: completed samples are recorded in <output>.done and skipped.

    python 04_mtdna_haplogroups.py --work-dir DIR --bam-dir DIR [--jobs 8] [--threads 12]
"""
import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

TREE = 'phylotree-fu-rcrs@1.2'


def run(cmd, desc):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return (True, '') if r.returncode == 0 else (False, f'{desc}: {r.stderr.strip()[:300]}')


def process_sample(sid, cfg):
    chrm_bam = f'{cfg.work_dir}/mt_fastq/{sid}_chrM.bam'
    fastq = f'{cfg.work_dir}/mt_fastq/{sid}_chrM.fastq'
    bam_rcrs = f'{cfg.work_dir}/mt_bam_rCRS/{sid}_rCRS.bam'
    vcf = f'{cfg.work_dir}/mt_vcf/{sid}_remap.vcf.gz'
    hg = f'{cfg.work_dir}/mutserve/{sid}_hg.tsv'

    steps = [
        (f'{cfg.samtools} view -b {cfg.bam_dir}/{sid}/{sid}.bam chrM -o {chrm_bam}',
         'samtools view'),
        (f'{cfg.samtools} fastq {chrm_bam} > {fastq}', 'samtools fastq'),
    ]
    for cmd, desc in steps:
        ok, msg = run(cmd, desc)
        if not ok:
            return False, msg
    if not os.path.exists(fastq) or os.path.getsize(fastq) == 0:
        return False, 'no chrM reads extracted'

    steps = [
        (f'{cfg.bwa} mem -t {cfg.threads} {cfg.rcrs} {fastq} | '
         f'{cfg.samtools} sort -o {bam_rcrs}', 'bwa mem'),
        (f'{cfg.samtools} index {bam_rcrs}', 'samtools index'),
        (f'{cfg.mutserve} call --reference {cfg.rcrs} --output {vcf} '
         f'--threads {cfg.threads} --mapQ 20 {bam_rcrs}', 'mutserve'),
        (f'{cfg.haplogrep} classify --in {vcf} --tree {TREE} --out {hg} '
         f'--extend-report', 'haplogrep3'),
    ]
    for cmd, desc in steps:
        ok, msg = run(cmd, desc)
        if not ok:
            return False, msg

    for f in (fastq, chrm_bam):
        if os.path.exists(f):
            os.remove(f)
    return True, ''


def main(cfg):
    for d in ('mt_fastq', 'mt_bam_rCRS', 'mt_vcf', 'mutserve'):
        os.makedirs(f'{cfg.work_dir}/{d}', exist_ok=True)

    samples = sorted(s for s in os.listdir(cfg.bam_dir)
                     if os.path.exists(f'{cfg.bam_dir}/{s}/{s}.bam'))
    done_log = cfg.out + '.done'
    done = set()
    if os.path.exists(done_log):
        done = {l.strip() for l in open(done_log) if l.strip()}
    pending = [s for s in samples if s not in done]
    print(f'{len(samples)} samples, {len(pending)} pending, '
          f'{cfg.jobs} workers x {cfg.threads} threads')

    failed = []
    if pending:
        with open(done_log, 'a') as fh, ThreadPoolExecutor(cfg.jobs) as pool:
            futures = {pool.submit(process_sample, s, cfg): s for s in pending}
            for i, fut in enumerate(as_completed(futures), 1):
                sid = futures[fut]
                ok, msg = fut.result()
                if ok:
                    fh.write(sid + '\n')
                    fh.flush()
                else:
                    failed.append((sid, msg))
                    print(f'  FAILED {sid}: {msg}', file=sys.stderr)
                print(f'\r  {i}/{len(pending)}', end='', file=sys.stderr)
        print(file=sys.stderr)

    frames, missing = [], []
    for sid in samples:
        p = f'{cfg.work_dir}/mutserve/{sid}_hg.tsv'
        if not os.path.exists(p):
            missing.append(sid)
            continue
        d = pd.read_csv(p, sep='\t')
        if len(d):
            d['SampleID_clean'] = sid
            frames.append(d)
    if not frames:
        raise RuntimeError('no haplogrep3 results to aggregate')

    pd.concat(frames, ignore_index=True).to_csv(cfg.out, sep='\t', index=False)
    print(f'aggregated {len(frames)} samples -> {cfg.out}')
    if missing:
        print(f'missing: {len(missing)}')
    if failed:
        print(f'failed: {len(failed)}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--work-dir', required=True)
    ap.add_argument('--bam-dir', required=True)
    ap.add_argument('--out', default=None)
    ap.add_argument('--rcrs', default=None)
    ap.add_argument('--mutserve', default=None)
    ap.add_argument('--haplogrep', default=None)
    ap.add_argument('--samtools', default='samtools')
    ap.add_argument('--bwa', default='bwa')
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--threads', type=int, default=12)
    cfg = ap.parse_args()
    cfg.out = cfg.out or f'{cfg.work_dir}/kz_219_mtdna_haplogroups.tsv'
    cfg.rcrs = cfg.rcrs or f'{cfg.work_dir}/rCRS.fasta'
    cfg.mutserve = cfg.mutserve or f'{cfg.work_dir}/mutserve'
    cfg.haplogrep = cfg.haplogrep or f'{cfg.work_dir}/haplogrep3'
    main(cfg)
