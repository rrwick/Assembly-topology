# Assembly topology from long-read alignments

The repo contains a Python script for checking the topology of contigs in a complete long-read bacterial genome assembly. By examining reads which align to contig ends, it can spot evidence for contigs that look: circular (reads loop to other end of the contig), hairpin (reads loop to the same end of the contig but opposite strand) or open (reads terminate at the contig end).



## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
