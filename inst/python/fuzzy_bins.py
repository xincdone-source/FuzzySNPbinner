#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fuzzy_bins.py
---------------

"""
import os
import csv
import math
import numpy as np
from collections import defaultdict

try:
    from vendor_snpbinner import bins as _bins
except ImportError:
    _bins = None


# =============== 辅助函数 ===============
def _estimate_stats_from_crosp(crosp_path: str):
    """
  
    """
    import re

    positions = []
    snp_count = 0
    sample_count = 0

    with open(crosp_path, newline='', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip().strip(',')
            if not line:
                continue
            parts = [p.strip() for p in re.split(r'[,\t]+', line) if p.strip() != '']
            if len(parts) < 3:
                continue
            sample_count += 1

            # 格式②: R01, 0, h, 6726336, b, 46081603, ...
            # 提取偶数位（1,3,5...）的数字作为 position
            pos_candidates = []
            for i, val in enumerate(parts[1:], start=1):
                if i % 2 == 1:  # 偶数序列（位置列）
                    try:
                        pos = int(float(val))
                        pos_candidates.append(pos)
                    except ValueError:
                        continue
            if pos_candidates:
                positions.extend(pos_candidates)
                snp_count += len(pos_candidates)

    chrom_len = max(positions) if positions else 1
    if chrom_len <= 0:
        chrom_len = 1_000_000

    if snp_count == 0:
        snp_count = 100

    D = snp_count / (chrom_len / 1e6)
    N = sample_count if sample_count > 0 else 200

    return {'D': D, 'N': N, 'chrom_len': chrom_len}


# =============== 主函数入口 ===============
def bins_fuzzy(input_path: str, output_path: str, min_bin_size: int):
    """
    对 crosspoints 输出进行分箱 (手动指定 bin_size)
    """
    print(f"📘 [Fuzzy] Running binning on {os.path.basename(input_path)}")

    if not os.path.exists(input_path):
        raise FileNotFoundError(f"❌ 找不到输入文件：{input_path}")

    stats = _estimate_stats_from_crosp(input_path)
    bin_size = int(min_bin_size) 

    print(f"[STATS] Estimated D={stats['D']:.2f}, N={stats['N']}, chrom_len={stats['chrom_len']}")
    print(f"[MANUAL] Using min_bin_size={bin_size}")

    # 运行分箱算法
    if _bins is not None:
        try:
            _bins(
                input_path=input_path,
                output_path=output_path,
                min_bin_size=bin_size
            )
            print(f"✅ [Fuzzy] bins 完成: {output_path}")
        except Exception as e:
            print(f"❌ [Fuzzy] bins 出错: {e}")
    else:
        print("⚠️ [Fuzzy] vendor_snpbinner 未加载，尝试生成简单占位输出。")
        # 占位输出（防止中断）
        with open(output_path, "w", newline='') as out:
            writer = csv.writer(out)
            writer.writerow(["chrom", "start", "end", "bin_id"])
            writer.writerow(["NC_unknown", 1, 1000000, "BIN001"])
        print(f"✅ [Fuzzy] 占位输出完成: {output_path}")
