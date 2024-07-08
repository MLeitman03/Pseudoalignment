import pandas as pd
import matplotlib.pyplot as plt


file = open("/Users/madelaineleitman/Downloads/chr11_transcriptome.fasta")
sequences = []
classes = []

i = 0
toAppend = ""
for line in file:
    if line[0] == ">":
        sequences.append(toAppend)
        classes.append(line[1:].strip())
        toAppend = ""
    else:
        toAppend = toAppend + line.strip()
sequences.append(toAppend)

sequences = sequences[1:]


class Node:
    def __init__(self, sequence, eqClass):
        self.seq = sequence
        self.eqClass = eqClass
        self.distToBreakPoint = -1
        self.next = []
        self.back = []
        self.visited = False


def addToGraph(dbGraph, len_kmer, sequence, eqClass):
    priorNode = None
    for i in range(len(sequence) - len_kmer + 1):
        kmer = sequence[i:i + len_kmer]
        if kmer in dbGraph:
            if eqClass:
                dbGraph[kmer].eqClass.append(eqClass)
            if priorNode != None:
                dbGraph[priorNode].next.append(dbGraph[kmer])
                dbGraph[kmer].back.append(dbGraph[priorNode])
        else:
            dbGraph[kmer] = Node(kmer, [eqClass])
            if priorNode != None:
                dbGraph[priorNode].next.append(dbGraph[kmer])
                dbGraph[kmer].back.append(dbGraph[priorNode])
        priorNode = kmer


dbGraph = dict()
len_kmer = 20

i = 0
for seq in sequences:
    addToGraph(dbGraph, len_kmer, seq, classes[i])
    i += 1

for key in dbGraph.keys():
    dbGraph[key].next = list(set(dbGraph[key].next))
    dbGraph[key].back = list(set(dbGraph[key].back))
    dbGraph[key].eqClass = list(set(dbGraph[key].eqClass))


def count_to_breakpoint(dbGraph, key):
    if len(list(dbGraph[key].back)) == 0 or dbGraph[key].distToBreakPoint == -1:
        return
    stack = []
    stack.append(dbGraph[key])
    while len(stack) > 0:
        node = stack.pop()
        node.visited = True
        for temp_node in node.back:
            if temp_node.visited == False:
                if temp_node.distToBreakPoint != 0:
                    temp_node.distToBreakPoint = node.distToBreakPoint + 1
                stack.append(temp_node)


for key in dbGraph.keys():
    dbGraph[key].distToBreakPoint = -1  # Reset
    if len(dbGraph[key].next) != 1:  # Check if there are multiple edges
        dbGraph[key].distToBreakPoint = 0
    else:
        if dbGraph[key].next[0].seq == key:
            dbGraph[key].distToBreakPoint = 0
    dbGraph[key].visited = False  # Reset

for key in dbGraph.keys():
    if dbGraph[key].visited == False:
        count_to_breakpoint(dbGraph, key)

file_reads = open("/Users/madelaineleitman/Downloads/reads.fasta")
reads = []
metadata = []

toAppend = ""
for line in file_reads:
    if line[0] == ">":
        reads.append(toAppend)
        metadata.append(line[1:].strip())
        toAppend = ""
    else:
        toAppend = toAppend + line.strip()

reads.append(toAppend)

reads = reads[1:]

from functools import reduce


def intersection(lists):
    if not lists:
        return []
    return list(reduce(set.intersection, map(set, lists)))


def findEqClass(read, len_kmer, dbGraph):
    eqClass = []
    i = 0
    while i < (len(read) - len_kmer + 1):
        kmer = read[i:i + len_kmer]
        if kmer not in dbGraph:
            return None
        else:
            node = dbGraph[kmer]
            eqClass.append(node.eqClass)
            if node.distToBreakPoint < 0:
                return None
            elif i + node.distToBreakPoint < 0:
                break
            else:
                i += node.distToBreakPoint + 1
    return intersection(eqClass)


def reverseComplement(sequence):
    output = ""
    sequence = sequence[::-1]
    for char in sequence:
        if char == "A":
            output = output + "T"
        if char == "T":
            output = output + "A"
        if char == "C":
            output = output + "G"
        if char == "G":
            output = output + "C"
    return output


read_library = dict()
read_library["None"] = 1
toAdd = ""
m = 1

for read in reads:
    forward = findEqClass(read, 20, dbGraph)
    backwards = findEqClass(reverseComplement(read), 20, dbGraph)
    if forward == None and backwards == None:
        toAdd = None
    elif forward == None:
        toAdd = ';'.join(sorted(backwards))
    elif backwards == None:
        toAdd = ';'.join(sorted(list(set(forward))))
    if toAdd == None:
        read_library["None"] += 1
    elif toAdd in read_library:
        read_library[toAdd] += 1
    else:
        read_library[toAdd] = 1
df = pd.DataFrame(columns=['counts', 'number of items in equivalence class', 'isoforms in equivalence class'])

t = 0
for read in read_library:
    temp_row = pd.DataFrame({'counts': read_library[read],
                             'number of items in equivalence class': (read.count(';') + 1),
                             'isoforms in equivalence class': str(read)}, index=[0])
    df = pd.concat([df, temp_row], ignore_index=True)
df.iloc[0, 1] = 0

total_equivalence_classes = df.shape[0]
mean_counts = df['counts'].mean()
median_counts = df['counts'].median()
mode_counts = df['counts'].mode()[0]
std_counts = df['counts'].std()
range_counts = df['counts'].max() - df['counts'].min()

mean_isoforms = df['number of items in equivalence class'].mean()
median_isoforms = df['number of items in equivalence class'].median()
mode_isoforms = df['number of items in equivalence class'].mode()[0]
std_isoforms = df['number of items in equivalence class'].std()
range_isoforms = df['number of items in equivalence class'].max() - df['number of items in equivalence class'].min()

stats_data = {
    'Statistic': ['Total Equivalence Classes', 'Mean Counts', 'Median Counts', 'Mode Counts',
                  'Standard Deviation of Counts', 'Range of Counts',
                  'Mean Number of Isoforms per Equivalence Class', 'Median Number of Isoforms per Equivalence Class',
                  'Mode Number of Isoforms per Equivalence Class',
                  'Standard Deviation of Isoforms per Equivalence Class', 'Range of Isoforms per Equivalence Class'],
    'Value': [total_equivalence_classes, mean_counts, median_counts, mode_counts, std_counts, range_counts,
              mean_isoforms, median_isoforms, mode_isoforms, std_isoforms, range_isoforms]
}

stats_df = pd.DataFrame(stats_data)

number_of_bins = 50
plt.figure(figsize=(10, 6))
plt.hist(df['counts'], bins=number_of_bins, alpha=0.75)
plt.yscale('log')
plt.title('Distribution of Counts per Equivalence Class')
plt.xlabel('Counts')
plt.ylabel('Frequency (log scale)')
plt.xlim(0, 20000)
plt.grid(True)

plt.savefig('/Users/madelaineleitman/Downloads/CS121/project_histogram.png')

stats_df.to_csv('/Users/madelaineleitman/Downloads/CS121/statistics.csv', index=False)

df.to_csv('/Users/madelaineleitman/Downloads/CS121/equiv_classes.csv')

