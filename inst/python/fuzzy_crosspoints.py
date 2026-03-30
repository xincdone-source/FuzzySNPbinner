# -*- coding: utf-8 -*-
"""
Fuzzy SNPbinner crosspoints module (for TSV genotype matrix)
------------------------------------------------------------
- 输入：原始 SNP 基因型矩阵（制表符分隔）
  表头示例：marker <TAB> position(bp) <TAB> R01 ... R202 <TAB> p1 <TAB> p2
  行示例：  NC_063452.1_3624 <TAB> 3624 <TAB> h b - ... <TAB> a <TAB> b
- 目标：从 TSV 估算 D/H/N -> 模糊推理 HMM 参数；min_length/min_bin_size 由用户手动给定
- 依赖：vendor_snpbinner.crosspoints（Python 3 兼容版），fuzzy_module（含 infer_hmm_base）
"""

from __future__ import annotations
import os, sys, csv

# 确保当前目录在路径中
sys.path.append(os.path.dirname(__file__))

# 优先使用 vendored 内核
try:
    import vendor_snpbinner.crosspoints as vsc
    _crosspoints = getattr(vsc, "crosspoints", vsc)
except Exception as e:
    print("⚠️ FuzzySNPbinner: vendor_snpbinner.crosspoints 未能导入，尝试系统 snpbinner：", e)
    from snpbinner.crosspoints import crosspoints as _crosspoints

# 模糊模块
from fuzzy_module import calibrate_memberships, infer_hmm_overrides, infer_hmm_base


HET_CODES = {'h', '0/1', '1/0', 'ab', 'ba', 'H', 'AB', 'BA'}

def _estimate_stats_from_tsv(tsv_path: str):
    """
    稳健读取TSV基因型矩阵，自动处理：
      - 真正的 tab (\t) 或者 字面“\\t”
      - 多个空格分隔
      - UTF-8/GBK/CRLF/BOM
    计算：
      D = 行数 / (max(position)/1e6)
      H = 非缺失基因型中 'h'(含0/1,ab等) 比例
      N = 样本列数（header-2）
    """
    import io, re

    def split_smart(line: str):
        # 优先真tab；否则把字面 \t 转真tab；再尝试按≥2空格分割；最后退化为任何空白
        if '\t' in line:
            return line.rstrip('\r\n').split('\t')
        if '\\t' in line:
            line = line.replace('\\t', '\t')
            return line.rstrip('\r\n').split('\t')
        parts = re.split(r'\s{2,}', line.strip())
        if len(parts) > 1:
            return parts
        return line.strip().split()  # 任意空白

    # 试探多种编码
    raw = None
    for enc in ('utf-8-sig', 'utf-8', 'gbk', 'latin1'):
        try:
            with open(tsv_path, 'r', encoding=enc, errors='strict') as f:
                raw = f.read()
                break
        except Exception:
            continue
    if raw is None:  # 最后兜底
        with open(tsv_path, 'r', encoding='utf-8', errors='ignore') as f:
            raw = f.read()

    lines = raw.splitlines()
    if not lines:
        raise ValueError("空文件：%s" % tsv_path)

    header = split_smart(lines[0])
    if len(header) < 3:
        raise ValueError("表头至少需要3列：marker, position(bp), samples... 实际=%d" % len(header))
    n_samples = max(0, len(header) - 2)

    snp_rows = 0
    het_cells = 0
    non_missing_cells = 0
    max_pos = 0

    # 常见杂合编码
    HET_CODES = {'h','H','0/1','1/0','ab','ba','AB','BA'}

    for line in lines[1:]:
        parts = split_smart(line)
        if len(parts) < 3:
            continue
        # 位置列（第二列）
        try:
            pos = int(parts[1])
        except Exception:
            continue
        max_pos = max(max_pos, pos)
        snp_rows += 1

        for g in parts[2:]:
            gs = (g or '').strip()
            if gs == '-' or gs == '':
                continue
            non_missing_cells += 1
            if gs in HET_CODES or gs.lower() == 'h':
                het_cells += 1

    # 统计量
    if max_pos <= 0 or snp_rows <= 0:
        print("⚠️ 未读到有效SNP位点，使用默认值")
        chrom_len = 64950979  # 你这条染色体长度
        D = 100.0
        H = 0.2
        N = n_samples if n_samples > 0 else 1
        return {'D': D, 'H': H, 'N': N, 'chrom_len': chrom_len}

    chrom_len = max_pos
    D = snp_rows / (chrom_len / 1_000_000.0)
    H = (het_cells / float(non_missing_cells)) if non_missing_cells > 0 else 0.0
    N = n_samples if n_samples > 0 else 1

    # 防止0触发log
    if D <= 0: D = 1e-6
    if N <= 0: N = 1

    return {'D': D, 'H': H, 'N': N, 'chrom_len': chrom_len}



def crosspoints_fuzzy(input_path: str,
                      output_path: str,
                      min_length: int = 5000,
                      min_bin_size: int = 500000,
                      seq_depth: float = 10.0): # <--- 1. 新增参数
    """
    执行基于模糊数学的 crosspoints 识别流程（输入为原始TSV基因型矩阵）
    新增 seq_depth 参数，默认为 10X
    """
    print(f"📘 [Fuzzy] Running fuzzy crosspoints on {os.path.basename(input_path)}")
    # <--- log 中显示深度
    print(f"[FUZZY] Settings: min_len={min_length}, min_bin={min_bin_size}, Depth={seq_depth}X") 

    # 1) 从TSV估算统计量
    stats = _estimate_stats_from_tsv(input_path)
    
    # <--- 2. 将深度注入 stats 字典
    stats['Depth'] = float(seq_depth) 
    
    print("[FUZZY] Estimated stats:", stats)

    # 2) 模糊隶属度 & HMM参数
    mu = calibrate_memberships(stats)
    print("[FUZZY] Memberships:", mu)

    base = infer_hmm_base(mu, stats)   
    overrides = infer_hmm_overrides(mu, stats) 
    print("[FUZZY] HMM base params:", base)
    print("[FUZZY] HMM overrides:", overrides)

    # 将 overrides 通过环境变量传给内核
    if 'intra_p' in overrides:
        os.environ['FSNP_INTRA_P'] = str(overrides['intra_p'])
    if 'trans_mult' in overrides:
        os.environ['FSNP_TRANS_MULT'] = str(overrides['trans_mult'])

    # 3) 灵敏度保护
    pred_cross = max(100, int(base.get('predicted_cross_count', 300)))
    pred_homo  = float(base.get('predicted_homogeneity', 0.80))
    min_len = int(min_length) 
    print(f"[FUZZY] Adjusted HMM sensitivity: cross_count={pred_cross}, min_length={min_len}, homo={pred_homo:.3f}")

    # 4) 调用原始 crosspoints 内核
    try:
        _crosspoints(
            input_path=input_path,
            output_path=output_path,
            predicted_homogeneity=pred_homo,
            predicted_cross_count=pred_cross,
            chrom_len=int(stats['chrom_len']),
            min_state_length=min_len,
            min_state_ratio=None 
        )
        print(f"✅ [Fuzzy] crosspoints 完成: {output_path}")
    except Exception as e:
        print("❌ [Fuzzy] crosspoints 出错：", e)
    finally:
        for k in ('FSNP_INTRA_P', 'FSNP_TRANS_MULT'):
            if k in os.environ:
                del os.environ[k]
