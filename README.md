# Assembly topology from long-read alignments

The repo contains a Python script for checking the topology of contigs in a complete long-read bacterial genome assembly. By examining reads which align to contig ends, it can spot evidence for contigs that look: circular (reads loop to other end of the contig), hairpin (reads loop to the same end of the contig but opposite strand) or open (reads terminate at the contig end).

Run it like this:
```bash
minimap2 -c -t 16 assembly.fasta reads.fastq.gz > alignments.paf
./assembly_topology.py assembly.fasta alignments.paf
```

And it will produce a tab-delimited output like this:
```
name          circular_reads  hairpin_reads  terminating_reads  clipping_reads  ambiguous_reads
chromosome_1  248             0              1                  2               8
chromosome_2  0               88             0                  5               0
plasmid_1     219             0              0                  2               9
plasmid_2     167             0              0                  2               4
plasmid_3     174             0              0                  6               3
plasmid_4     0               319            3                  21              3
plasmid_5     316             0              1                  16              4
```

Note: This script is experimental. It assumes a finished, high-quality assembly and is intended to provide evidence for each sequence's topology. On draft or lower-quality assemblies (e.g. fragmented, poorly trimmed or misassembled), the results may be ambiguous or misleading.


## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
