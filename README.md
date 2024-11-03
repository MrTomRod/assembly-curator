# Assembly Curator: Automated Comparison and Selection Tool for Bacterial Genome Assemblies

This tool streamlines the process of comparing and selecting the best bacterial genome assembly from multiple assemblers
such as Flye, hifiasm, PacBio SMRTtools, and LJA. By supporting a flexible plugin system, it can handle various output
formats—including FASTA, GFA, and non-standard types—and is adaptable to future assemblers and file formats.

# weird stuff
spooky lja contigs:
  [FAM27700-lja-8.fasta](../../Downloads/FAM27700-lja-8.fasta)
  [FAM25886-lja-12.fasta](../../Downloads/FAM25886-lja-12.fasta)
  [FAM25181-lja-8.fasta](../../Downloads/FAM25181-lja-8.fasta)
  [FAM25180-lja-6.fasta](../../Downloads/FAM25180-lja-6.fasta)
  [FAM25180-lja-4.fasta](../../Downloads/FAM25180-lja-4.fasta)
  [FAM23202-lja-2.fasta](../../Downloads/FAM23202-lja-2.fasta)
  [FAM22009-lja-6.fasta](../../Downloads/FAM22009-lja-6.fasta)
  [FAM21963-lja-1.fasta](../../Downloads/FAM21963-lja-1.fasta)
  [FAM21963-lja-7.fasta](../../Downloads/FAM21963-lja-7.fasta)
  [FAM17935-lja-6.fasta](../../Downloads/FAM17935-lja-6.fasta)
  [FAM17993-lja-10.fasta](../../Downloads/FAM17993-lja-10.fasta)

# Todo

- [ ] Button: blast random 500 bp on NCBI? Or at least copy to clipboard?
  - https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html
- [ ] Dotplot interpretation
  tutorial: https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/650/Examples_interpretations_dot_plots.html
- [ ] Noam: _L. helveticus_ might be a more interesting example because it's more complex
- [ ] assemblies.html: link to folder, reset sample button, show files per assembly
- [ ] Dotplots: same order as ANI-matrix
- [ ] overview.html: click on sample updates icon to loading icon
- [ ] Hover on dotplots to show conting names and metadata like ANI matrix
- [ ] Context menus for contigs: show sequence
- [ ] Expected GC-content: if below 25 or above 65 -> warning
- [ ] Expected GC-content: if different from chromosome by more than 10 % -> warning
- [x] small pie chart ATGC
- [ ] gtdb-tk matches?
- [ ] PlasmidFinder 2.1 -> taxonomy of plasmid?
- [x] Curation: manual fixing of headers
- [ ] Curation: write output files
    - [x] FASTA
    - [ ] GFA?
- [x] Coverage
- [x] gfapy issue: https://github.com/ggonnella/gfapy/issues/31
    - [x] replace gfapy
- [ ] Autosolve easy cases
- [ ] Dotplot: ani-matrix -> dotplot; whole assembly dotplot ("I forgot to mention also the duplicated/chimeric
  pseudo-plasmids/contigs. They will clearly show up in dotplot. Actually with gepard I was doing all-against-all
  comparisons, not just a self")
- [ ] DnaA centering
    - https://github.com/sanger-pathogens/circlator/wiki/Task:-fixstart
    - Dnaapler https://joss.theoj.org/papers/10.21105/joss.05968
        - Way too slow! Can be automated for circular contigs. Should be part of the post-assembly process?
- [x] Hatice: "I will send you some examples of dotplots that highlight the artefactual circularized contigs looking like plasmids"
- [x] Interactive web server: more lightweight than Django?
    - https://www.pyweb.io/
    - https://bottlepy.org/docs/dev/
- [x] Dotplot: What I have is nice, but maybe overkill?
    - Maybe I should create a more static multiplot per cluster
      like [trycycler](https://github.com/rrwick/Trycycler/blob/77d1be8cf786969a3821f3c04086df60b39ff781/trycycler/dotplot.py#L38)
    - Ignore for now?

See
also: [Snakemake workflows for long-read bacterial genome assembly and evaluation](https://gigabytejournal.com/articles/116)

# Workflow ideas

Run the tool (assembly-curator load) once to create htmls for all assemblies.

Then, run the tool (assembly-curator serve) to start a web server:

- Interactive webserver
- User selects contigs
- Click on export creates curated assemblies

# Assembler Tools

## Problem

Today's long-read HiFi sequencing technologies usually enable us to create near-perfect genome assemblies for bacterial
genomes: chromosomes and plasmids can often be reconstructed into circular sequences. However, different assemblers
often give different results: one assembly is lacking a plasmid, another failed to close the chromosome, etc.

## Solution

My strategy is to use 3-4 different assemblers on each dataset and then manually decide which one is best. However, this
is a time-consuming process. I would like to automate this process by creating a tool that shows nice summaries of all
assemblies next to each other to facilitate the decision-making process.

Default assemblers:

- Flye
- hifiasm
- PacBio SMRTtools
- LJA

## Features

1) **Read in multiple assemblies**

Each assember produces **different output**: a fasta file, a gfa file, files in non-standard formats. Moreover, there
may be new assemblers in the future that produce new file formats. The tool should be able to read in all of these by
implementing a simple **plugin system**.

Bioinformaticians may want to run **additional tools** on the assemblies before deciding which one is best, for example
BUSCO to check for completeness or some other tools to check for contamination.

Error handling: Sometimes the assemblers fail to produce an output. The tool should be able to handle this and display
an error message instead of crashing, and propagate warnings.

2) **Overview of all samples**

An overview table that shows all samples, maybe some metadata, the number of successful assemblies, whether the assembly
has already been curated.

3) **Show a summary of each assembly**

Since the output may consist of different media types (text, images, tables, etc.), the tool should be able to display
all of these. Moreover, they should be displayed in a way that is easy to compare. For example, the assembly graphs
should be displayed next to each other.

Moreover, it would be nice to have a phylogenetic tree of all contigs, which could be used to identify equivalent
contigs. (GenDisCal is fast enough to be run on the fly.)

This necessitates building a dynamic **user interface**, which should be web-based to make it easy to use on any
platform.

4) **Allow the user to select the best assembly**

The user should be able to select the best assembly, or create a hybrid assembly by combining select contigs of multiple
assemblies. This could be done by drag-and-drop or by selecting via click...

Originally, I thought that this means the tool needs a **backend** that can handle such requests, and that static html
files are not enough. However, it might be possible to implement clicking and selection of contigs in a static html file
which then generates json to copy-paste back into the tool. This would simplify the implementation and reduce the number
of dependencies.

It might be nice to generate an interactive dotplot (https://github.com/MrTomRod/dot) on the fly. This would require a
backend, though...

This also includes curation of the selected contigs, e.g. classifying them as chromosome, plasmid, etc.

5) **Output the best assembly**

Write the selected contigs to a new fasta file:

- Compatible with NCBI PGAP: https://github.com/ncbi/pgap/wiki/Input-Files
    - Be less than 50 characters long
    - Only `^[A-Za-z0-9\-_\.:\*#]+$`
    - Be unique within a genome
- Include relevant metadata in the header if available, e.g.:
    - `[topology=circular] [completeness=complete]` (PGAP standard)
    - `[location=chromosome]` or `[location=plasmid] [plasmid-name=pORGANISM_1]` (PGAP standard)
    - `[length=1234567]` or `[coverage=123x]` (custom, but nice to have)
    - `[assembler=flye]` (custom, especially useful for hybrid assemblies)
- A log file and/or a markdown file for OGB

Example:

```
>ORGANISM_scf1 [topology=circular] [completeness=complete] [location=chromosome] [length=1234567] [coverage=123x] [assembler=flye]
ATGC...
>ORGANISM_scf2 [topology=circular] [completeness=complete] [location=plasmid] [plasmid-name=pORGANISM_1] [length=85642] [coverage=354x] [assembler=hifiasm]
ATGC...
```
