# -*- coding: utf-8 -*-
# kontrol kachestva readov
# Петрова М.С. 2024
# FINAL2 !!! рабочая версия, не трогать
# 02.04 добавил графики
# 05.04 поменял порог качества (Саша сказал 25)

import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import time
# import pandas as pd
# from Bio import SeqIO
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = BASE_DIR / "data" / "raw" / "reads_2024_04_02.fastq"
DEFAULT_RESULTS = BASE_DIR / "results"

ADAPTER = "AGATCGGAAGAGC"

STATS = []
SKIPPED = []



def read_fastq(path):
    f = open(path, encoding="utf-8")
    lines = f.readlines()
    f.close()
    print("строк в файле", len(lines))
    reads = []
    i = 0
    while i < len(lines):
        h = lines[i].replace("\n", "")
        s = lines[i + 1].replace("\n", "")
        q = lines[i + 3].replace("\n", "")
        if len(s) != len(q):
            SKIPPED.append(h)
            i = i + 4
            continue
        reads.append([h, s, q])
        i = i + 4
    return reads


def gc(s):
    g = 0
    c = 0
    for x in s:
        if x == "G":
            g = g + 1
        if x == "C":
            c = c + 1
    return (g + c) / len(s) * 100


def meanq(q):
    t = 0
    for x in q:
        t = t + (ord(x) - 33)
    return t / len(q)

def main():

    t0 = time.time()

    input_path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_INPUT

    results_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_RESULTS

    results_dir.mkdir(parents=True, exist_ok=True)

    r = read_fastq(input_path)
    print("прочтений", len(r))

    good = []
    bad = []
    for x in r:
        try:
            z = meanq(x[2])
        except:
            pass
        s = x[1]
        if s[0:13] == ADAPTER:
            s = s[13:]
            x[2] = x[2][13:]
        if z >= 25 and len(s) >= 100:
            good.append([x[0], s, x[2]])
        else:
            bad.append(x[0])
        STATS.append([x[0], len(s), gc(s), z])

    print("хороших", len(good))
    print("плохих", len(bad))
    print("отброшено при чтении", len(SKIPPED))

    # ---- GC по хорошим ----
    g1 = []
    for x in good:
        g1.append(gc(x[1]))
    m1 = sum(g1) / len(g1)
    print("средний GC хороших", m1)

    # ---- GC по всем ----
    g2 = []
    for x in r:
        g2.append(gc(x[1]))
    m2 = sum(g2) / len(g2)
    print("средний GC всех", m2)

    qq = []
    for x in STATS:
        qq.append(x[3])
    print("среднее качество по всем прочтениям", sum(qq) / len(qq))

    nn = 0
    for x in r:
        if "N" in x[1]:
            nn = nn + 1
    print("прочтений с N", nn)
    print("доля N-содержащих", nn / len(r))

    # for x in STATS:
    #     if x[3] < 15:
    #         print("очень плохое", x[0])

    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1)
    plt.hist(g2, bins=40)
    plt.axvline(m2, color="r")
    plt.xlabel("GC, %")
    plt.ylabel("число прочтений")
    plt.title("GC-состав")
    plt.subplot(1, 2, 2)
    ll = []
    for x in STATS:
        ll.append(x[3])
    plt.hist(ll, bins=40)
    plt.axvline(25, color="r")
    plt.xlabel("среднее качество")
    plt.ylabel("число прочтений")
    plt.title("Качество прочтений")
    plt.tight_layout()
    plt.savefig(results_dir/"qc.png", dpi=200)
    print("график готов")

    f2 = open(results_dir / "clean_reads.fasta", "w", encoding="utf-8")
    for x in good:
        f2.write(">" + x[0][1:] + "\n")
        f2.write(x[1] + "\n")
    f2.close()

    f3 = open(results_dir / "qc_summary.tsv", "w", encoding="utf-8")
    f3.write("read\tlength\tgc\tmeanq\n")
    for x in STATS:
        f3.write(x[0][1:] + "\t" + str(x[1]) + "\t" + str(round(x[2], 2)) + "\t" + str(round(x[3], 2)) + "\n")
    f3.close()

    f4 = open(results_dir/"report.txt", "w", encoding="utf-8")
    f4.write("QC ОТЧЁТ\n")
    f4.write("файл: " + str(DEFAULT_INPUT) + "\n")
    f4.write("всего прочтений: " + str(len(r)) + "\n")
    f4.write("прошло фильтр: " + str(len(good)) + "\n")
    f4.write("средний GC: " + str(round(m2, 2)) + "\n")
    f4.write("прочтений с N: " + str(nn) + "\n")
    f4.close()

    print("отчёт записан в", DEFAULT_RESULTS/"report.txt")
    print("ГОТОВО")

if __name__ == "__main__":

    main()