# -*- coding: utf-8 -*-
# kontrol kachestva readov
# Петрова М.С. 2024
# FINAL2 !!! рабочая версия, не трогать
# 02.04 добавил графики
# 05.04 поменял порог качества (Саша сказал 25)

import os
import sys
import math
import random
import numpy as np
import matplotlib.pyplot as plt
# import pandas as pd
# from Bio import SeqIO

IN = "D:\\Работа\\секвенирование\\апрель 2024\\reads_2024_04_02.fastq"
OUTDIR = "D:\\Работа\\секвенирование\\апрель 2024\\результаты\\"

ADAPTER = "AGATCGGAAGAGC"

STATS = []
SKIPPED = []


def read_fastq():
    f = open(IN)
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


r = read_fastq()
print("прочтений", len(r))
print("!!! проверка", r[0][0])

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

# случайная подвыборка для оценки
sub = random.sample(STATS, 500)
qq = []
for x in sub:
    qq.append(x[3])
print("среднее качество по подвыборке", sum(qq) / 500)

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
plt.savefig(OUTDIR + "qc.png", dpi=200)
print("график готов")

f2 = open(OUTDIR + "clean_reads.fasta", "w")
for x in good:
    f2.write(">" + x[0][1:] + "\n")
    f2.write(x[1] + "\n")
f2.close()

f3 = open(OUTDIR + "qc_summary.tsv", "w")
f3.write("read\tlength\tgc\tmeanq\n")
for x in STATS:
    f3.write(x[0][1:] + "\t" + str(x[1]) + "\t" + str(round(x[2], 2)) + "\t" + str(round(x[3], 2)) + "\n")
f3.close()

f4 = open(OUTDIR + "report.txt", "w")
f4.write("QC ОТЧЁТ\n")
f4.write("файл: " + IN + "\n")
f4.write("всего прочтений: " + str(len(r)) + "\n")
f4.write("прошло фильтр: " + str(len(good)) + "\n")
f4.write("средний GC: " + str(round(m2, 2)) + "\n")
f4.write("прочтений с N: " + str(nn) + "\n")
f4.close()

print("отчёт записан в", OUTDIR + "report.txt")
print("ГОТОВО")
