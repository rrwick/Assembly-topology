# Assembly topology from long-read alignments

The repo contains a Python script for checking the topology of contigs in a complete long-read bacterial genome assembly. By examining reads which align to contig ends, it can spot evidence for contigs that look: circular (reads loop to other end of the contig), hairpin (reads loop to the same end of the contig but opposite strand) or open (reads terminate at the contig end).

Run it like this:
```bash
minimap2 -c -t 16 assembly.fasta reads.fastq.gz > alignments.paf
./assembly_topology.py assembly.fasta alignments.paf
```

And it will produce a tab-delimited output like this:
```
name          circular  hairpin  terminating  clipping  ambiguous
chromosome_1  248       0        1            2         8
chromosome_2  0         88       0            5         0
plasmid_1     219       0        0            2         9
plasmid_2     167       0        0            2         4
plasmid_3     174       0        0            6         3
plasmid_4     0         319      3            21        3
plasmid_5     316       0        1            16        4
```

Note: This script is experimental. It assumes a finished, high-quality assembly and is intended to provide evidence for each sequence's topology. On draft or lower-quality assemblies (e.g. fragmented, poorly trimmed or misassembled), the results may be ambiguous or misleading.



## Full usage

```
usage: assembly_topology.py [--tolerance TOLERANCE] [-h] assembly_fasta alignments_paf

Assembly topology from long-read alignments

Positional arguments:
  assembly_fasta        Assembly sequence in FASTA format
  alignments_paf        Long-read alignments in PAF format

Settings:
  --tolerance TOLERANCE
                        Allowed gap size when comparing alignment positions (default: 10)

Other:
  -h, --help            Show this help message and exit
```



## Details

This script processes each read in the alignments using the following steps:
* Discards the read if it has alignments to multiple different assembly contigs.
* Discards the read if its alignments do not reach the end of an assembly contig.
* Checks to see which (if any) of the following apply to a read:
  * Has circular alignments: The read crosses a contig boundary on the same strand – two alignments to the same target and same strand, one near the left end and the other near the right end of the contig, and the two alignments are adjacent on the read.
  * Has hairpin alignments: The read folds back at the same contig end – two alignments to the same target but opposite strands, both near the same contig end (left or right), and the two alignments are adjacent on the read.
  * Has terminating alignments: A single alignment ends at a contig end and the read also ends, consistent with a linear sequence with an open (non-hairpin) end.
* If only one of the above applies to a read, then count that read as circular, hairpin or terminating.
* If none of the above apply to a read, then count that read as clipping.
* If more than one of the above apply to a read, then count that read as ambiguous.
  * An example would be when a linear replicon has a terminal inverted repeat (TIR) with hairpin ends. If a read is entirely contained within the TIR, it may be consistent with both circular and hairpin topologies.



## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
