import os
import re
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

genes = []
with open("genes_HLA_full.txt", "r") as f:
    genes = [gene.strip() for gene in f.readlines()]

modes = ["0", "3", "5"]
recs = ["0", "1", "2"]

def parse_logs(directory, rec):
    times = {}
    memory = {}
    edit_distance = {}
    for gene in genes:
        for mode in modes:
            dir_path = f"{directory}/{gene}/{rec}/reads_{mode}_split/"
            count = 0
            count_gaf = 0
            for filename in os.listdir(dir_path):
                if filename.endswith(".log"):
                    count += 1
                    with open(os.path.join(dir_path, filename), "r") as f:
                        for line in f:
                            if "Elapsed (wall clock) time (h:mm:ss or m:ss):" in line:
                                pattern = r"(?:(\d+):)?(\d{1,2}):(\d{2})(?:\.(\d+))?"

                                # Cerca il tempo nella stringa
                                match = re.search(pattern, line)

                                if match:
                                    # Estrai ore, minuti, secondi e frazioni di secondo
                                    ore = int(match.group(1)) if match.group(1) else 0
                                    minuti = int(match.group(2)) if match.group(2) else 0
                                    secondi = int(match.group(3)) if match.group(3) else 0
                                    frazione_di_secondo = int(match.group(4)) if match.group(4) else 0
                                    
                                    # Converte tutto in millisecondi
                                    millisecondi_totali = (
                                        (ore * 3600 * 1000) +  # ore in millisecondi
                                        (minuti * 60 * 1000) +  # minuti in millisecondi
                                        (secondi * 1000) +  # secondi in millisecondi
                                        (frazione_di_secondo * (10 ** (3 - len(match.group(4)))) if match.group(4) else 0)  # frazioni di secondo
                                    )
                                    times[(gene, mode, count)] = millisecondi_totali
                            if "Maximum resident set size (kbytes)" in line:
                                memory[(gene, mode, count)] = int(line.split()[-1])
                if filename.endswith(".gaf"):
                    count_gaf += 1
                    if "ga" in directory:
                        with open(os.path.join(dir_path, filename), "r") as f:
                            for line in f.readlines():
                                cigar = line.strip().split("\t")[-1].split(",")[0].split(":")[-1]
                                edit_score = edit_distance_from_cigar(cigar)
                                if "recombination" in line:
                                    edit_score += 4
                            edit_distance[(gene, mode, count_gaf)] = edit_score
                    else:
                        with open(os.path.join(dir_path, filename), "r") as f:
                            for line in f.readlines():
                                if line.startswith("@CO"):
                                    edit_distance[(gene, mode, count_gaf)] = int(line.split("\t")[1])
                    
    return times, memory, edit_distance

def parse_cigar(cigar_string):
    
    operations = []
    current_length = ""
    for char in cigar_string:
        if char.isdigit():
            current_length += char
        else:
            operations.append((int(current_length), char))
            current_length = ""
    return operations

def edit_distance_from_cigar(cigar_string):
    
    operations = parse_cigar(cigar_string)
    edit_distance = 0
    for length, op in operations:
        if op in 'IDX':
            edit_distance += length
    return edit_distance

old_times, old_memory, old_edit = parse_logs("output/HLA/ra_f", 0)
new_times, new_memory, new_edit = parse_logs("output/HLA/ra_s", 0)


data_time = []
data_memory = []

for gene in genes:
    for mode in modes:
        for count in range(1, max(len(old_times), len(new_times)) + 1):
            if (gene, mode, count) in old_times and (gene, mode, count) in new_times:
                data_time.append((gene, mode, old_times[(gene, mode, count)], "Fast"))
                data_time.append((gene, mode, new_times[(gene, mode, count)], "Seed"))

            if (gene, mode, count) in old_memory and (gene, mode, count) in new_memory:
                data_memory.append((gene, mode, old_memory[(gene, mode, count)], "Fast"))
                data_memory.append((gene, mode, new_memory[(gene, mode, count)], "Seed"))



data_edit = []

for gene, mode, count in old_edit:
    data_edit.append((gene, mode, old_edit[(gene, mode, count)], "Fast"))
    data_edit.append((gene, mode, new_edit[(gene, mode, count)], "Seed"))

df_edit = pd.DataFrame(data_edit, columns=["Gene", "Mode", "Edit", "Version"])

import math
                    
genes_part1 = genes[:10]
genes_part2 = genes[10:20]

# Numero di colonne e righe
ncols = 2
nrows = 5

fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15, 4.5 * nrows))
axes = axes.flatten()

for i, gene in enumerate(genes_part1):
    df_edit_gene = df_edit[df_edit["Gene"] == gene]
    sns.boxplot(ax=axes[i], x="Mode", y="Edit", hue="Version", data=df_edit_gene)
    axes[i].set_title(f"ED for Gene: {gene}")
    axes[i].set_xlabel("Mode")
    axes[i].set_ylabel("Edit distance")
    axes[i].legend(title='Version')

# Rimuove gli assi non necessari se il numero di geni è dispari
for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.subplots_adjust(hspace=0.4)  # Aumenta lo spazio verticale

plt.savefig("/home/dcmonti/Scaricati/ga-vs-ra_ed_1.png")
plt.show()

fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15, 4.5 * nrows))
axes = axes.flatten()

for i, gene in enumerate(genes_part2):
    df_edit_gene = df_edit[df_edit["Gene"] == gene]
    sns.boxplot(ax=axes[i], x="Mode", y="Edit", hue="Version", data=df_edit_gene)
    axes[i].set_title(f"ED for Gene: {gene}")
    axes[i].set_xlabel("Mode")
    axes[i].set_ylabel("Edit distance")
    axes[i].legend(title='Version')

# Rimuove gli assi non necessari se il numero di geni è dispari
for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.subplots_adjust(hspace=0.4)  # Aumenta lo spazio verticale

plt.savefig("/home/dcmonti/Scaricati/ga-vs-ra_ed_2.png")
plt.show()

df_time = pd.DataFrame(data_time, columns=["Gene", "Mode", "Time", "Version"])

fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15, 4.5 * nrows))
axes = axes.flatten()

for i, gene in enumerate(genes_part1):
    df_time_gene = df_time[df_time["Gene"] == gene]
    sns.boxplot(ax=axes[i], x="Mode", y="Time", hue="Version", data=df_time_gene)
    axes[i].set_title(f"Time for Gene: {gene}")
    axes[i].set_yscale("log")
    axes[i].set_xlabel("Mode")
    axes[i].set_ylabel("Time (ms)")
    axes[i].legend(title='Version')

for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.subplots_adjust(hspace=0.4)  # Aumenta lo spazio verticale

plt.savefig("/home/dcmonti/Scaricati/ga-vs-ra_time_1.png")
plt.show()

fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15, 4.5 * nrows))
axes = axes.flatten()

for i, gene in enumerate(genes_part2):
    df_time_gene = df_time[df_time["Gene"] == gene]
    sns.boxplot(ax=axes[i], x="Mode", y="Time", hue="Version", data=df_time_gene)
    axes[i].set_title(f"Time for Gene: {gene}")
    axes[i].set_yscale("log")
    axes[i].set_xlabel("Mode")
    axes[i].set_ylabel("Time (ms)")
    axes[i].legend(title='Version')

for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.subplots_adjust(hspace=0.4)  # Aumenta lo spazio verticale

plt.savefig("/home/dcmonti/Scaricati/ga-vs-ra_time_2.png")
plt.show()

df_memory = pd.DataFrame(data_memory, columns=["Gene", "Mode", "Memory", "Version"])
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15, 4.5 * nrows))
axes = axes.flatten()

for i, gene in enumerate(genes_part1):
    df_memory_gene = df_memory[df_memory["Gene"] == gene]
    sns.boxplot(ax=axes[i], x="Mode", y="Memory", hue="Version", data=df_memory_gene)
    axes[i].set_title(f"Memory for Gene: {gene}")
    axes[i].set_yscale("log")
    axes[i].set_xlabel("Mode")
    axes[i].set_ylabel("Memory used (KB)")
    axes[i].legend(title='Version')

for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.subplots_adjust(hspace=0.4)  # Aumenta lo spazio verticale

plt.savefig("/home/dcmonti/Scaricati/ga-vs-ra_memory_1.png")
plt.show()

fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(15, 4.5 * nrows))
axes = axes.flatten()

for i, gene in enumerate(genes_part2):
    df_memory_gene = df_memory[df_memory["Gene"] == gene]
    sns.boxplot(ax=axes[i], x="Mode", y="Memory", hue="Version", data=df_memory_gene)
    axes[i].set_title(f"Memory for Gene: {gene}")
    axes[i].set_yscale("log")
    axes[i].set_xlabel("Mode")
    axes[i].set_ylabel("Memory used (KB)")
    axes[i].legend(title='Version')

for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.subplots_adjust(hspace=0.4)  # Aumenta lo spazio verticale

plt.savefig("/home/dcmonti/Scaricati/ga-vs-ra_memory_2.png")
plt.show()
