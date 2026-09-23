#!/usr/bin/env python3
"""Split the CDS contigs into N balanced shards for countreads.

Contigs are assigned largest-first to the currently-emptiest shard, so the
slowest shard is bounded by the largest single contig (chr1, ~10% of a mammalian
transcriptome). Sharding by contig is exact rather than approximate - see the
comment in CountingFullSeq_Apos.pl.

Usage: make_shards.py <cds.fa> <n_shards> <out_dir>
"""
import collections, os, sys

cds, n_shards, out_dir = sys.argv[1], int(sys.argv[2]), sys.argv[3]
os.makedirs(out_dir, exist_ok=True)

# same acceptance test as load_cds() in CountingFullSeq_Apos.pl
per_contig = collections.Counter()
name, seq = None, []
def take(name, seq):
    if not name:
        return
    a = name.split()
    if len(a) > 5 and 'protein_coding' in a[5] and len(''.join(seq)) > 120:
        per_contig[a[2].split(':')[2]] += 1
with open(cds) as fh:
    for line in fh:
        if line.startswith('>'):
            take(name, seq); name = line[1:].rstrip(); seq = []
        else:
            seq.append(line.strip())
take(name, seq)

buckets = [[] for _ in range(n_shards)]
load = [0] * n_shards
for contig, n in per_contig.most_common():           # largest contig first
    i = min(range(n_shards), key=lambda j: load[j])  # into the emptiest shard
    buckets[i].append(contig); load[i] += n

for i, b in enumerate(buckets):
    with open(os.path.join(out_dir, "%d.txt" % i), 'w') as fh:
        fh.write(''.join(c + '\n' for c in b))

total = sum(per_contig.values())
print("%d transcripts on %d contigs -> %d shards" % (total, len(per_contig), n_shards))
print("largest shard: %d transcripts (%.1f%%) => ceiling %.1fx"
      % (max(load), 100.0 * max(load) / total, total / max(load)))
