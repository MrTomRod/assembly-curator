<p align="center">
    <img src="media/icon-assembly-curator.svg" alt="Logo" width=170>
    <h1 align="center">Assembly Curator</h1>
</p>
<p align="center">
    <b>Quickly generate consensus assemblies for bacterial genomes</b>
</p>
<p align="center">
    <img src="media/imagined-assembly-curator.webp" width=500>
</p>
<p align="center">
    <sub><em>This is what an assembly curator looks like, at least according to OpenAI's DALL·E.</em></sub>
</p>

## 🚀 Overview

**Assembly Curator** is a semi-automated genome assembly curation tool designed to simplify and accelerate the creation
of high-quality bacterial genome assemblies. It bridges the gap between raw assemblies and annotated genomes, empowering
bioinformaticians to efficiently curate hundreds of genome assemblies with ease.

## 🎯 Motivation

Modern long-read sequencing technologies like PacBio HiFi generate near-complete bacterial genomes, yet different
assemblers often yield inconsistent results—missing plasmids, incomplete chromosomes, or misassemblies. Tools like
Trycycler offer a solution but require complex workflows and extensive processing time.

**Assembly Curator** offers a faster, more user-friendly alternative:

- Achieve **90% of the quality** in **10% of the time** compared to Trycycler.
- Designed to handle large datasets (50+ samples) **HiFi data** from PacBio Revio reads.
- Easily integrate **any assembler** through a flexible plugin system.

## ⚙️ Key Features

- **Interactive Web Interface:** Curate assemblies through a visual, browser-based UI.
- **Contig Clustering & Visualization:** ANI matrix clustering and interactive dotplots for contig comparison.
- **Plugin Architecture:** Seamless integration of new or custom assemblers.
- **High Performance:** Leveraging **pyskani** for clustering and **minimap2** for alignment visualization.
- **Parallel Processing:** Optimize curation time with background parallelization.
- **Direct Export:** Generate ready-to-annotate assemblies in [PGAP](https://github.com/ncbi/pgap/wiki/Input-Files)-ready FASTA format.

## 🛠️ Quick Start

See [QuickStart.md](QuickStart.md) for instructions on downloading test data.

## 📊 Core Workflow

1. **Preprocessing:** Generate multiple genome assemblies using any assembler.
2. **Curation:** Use the interactive app to cluster, compare, and select contigs.
3. **Export:** Save curated assemblies ready for genome annotation.

## 🔌 Plugin System

> [!CAUTION]
> Danger: work in progress, will be refactored

Add support for custom assemblers via Python-based plugins:

```python
from assembly_curator.AssemblyImporter import AssemblyImporter
from assembly_curator.Assembly import Assembly


class FlyeImporter1(AssemblyImporter):
    assembler = 'flye1'
    assembly_dir = 'flye1'
    assembly = 'assembly.fasta'
    gfa = 'assembly_graph.gfa'

    def load_assembly(self) -> Assembly:
        pass  # Override this method
```

## 📦 Output

- Interactive ANI matrix and dotplots for visual curation.
- Exported hybrid assemblies (`hybrid.fasta`) ready for annotation.

## 📚 Inspiration and Dependencies

- [Trycycler Inspiration](https://github.com/rrwick/Trycycler)
- [gfaviz](https://github.com/ggonnella/gfaviz)
- [Minimap2](https://github.com/lh3/minimap2)
- [Pyskani](https://github.com/althonos/pyskani)
- [Plotly](https://plotly.com/)
