# -*- coding: utf-8 -*-
"""
Estimate SNP density (D), heterozygosity (H), and sample size (N) from TSV.
Assumptions:
  - TSV columns: CHROM, POS, SAMPLE1, SAMPLE2, ...
  - Heterozygotes encoded among: '0/1','1/0','AB','BA','H'
"""
from __future__ import annotations
import csv
import collections

HET_CODES = set(['0/1', '1/0', 'AB', 'BA', 'H'])


def median(lst):
    if not lst:
        return 0.0
    s = sorted(lst)
    m = len(s) // 2
    return float(s[m]) if len(s) % 2 else 0.5 * (s[m - 1] + s[m])


def estimate_stats(tsv_path, window_bp=200000, sample_limit=50):
    snp_per_chr = collections.defaultdict(list)
    het_count = 0
    var_count = 0
    n_samples = None

    with open(tsv_path, newline='') as f:
        rdr = csv.reader(f, delimiter='\t')
        header = next(rdr)
        samples = header[2:] if len(header) > 2 else []
        if sample_limit and len(samples) > sample_limit:
            samples = samples[:sample_limit]
        n_samples = len(samples)

        cur_chr = None
        window_start = None
        window_snps = 0

        for row in rdr:
            if not row:
                continue
            chrom = row[0]
            try:
                pos = int(row[1])
            except Exception:
                continue
            gts = row[2:2 + n_samples]

            if cur_chr != chrom:
                if cur_chr is not None and window_start is not None:
                    snp_per_chr[cur_chr].append(window_snps)
                cur_chr, window_start, window_snps = chrom, pos, 0

            window_snps += 1
            hets_here = sum(1 for g in gts if g in HET_CODES)
            het_count += hets_here
            var_count += len(gts)

            if pos - window_start >= window_bp:
                snp_per_chr[chrom].append(window_snps)
                window_start, window_snps = pos, 0

        if cur_chr is not None and window_start is not None:
            snp_per_chr[cur_chr].append(window_snps)

    window_mb = window_bp / 1_000_000.0
    dens_list = []
    for counts in snp_per_chr.values():
        dens_list.extend([(c / window_mb) for c in counts])

    D = median(dens_list) if dens_list else 0.0
    H = (het_count / float(var_count)) if var_count else 0.0
    N = n_samples or 0
    return {'D': D, 'H': H, 'N': N}
