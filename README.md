# 🧬 Assembly Curator

**Quickly generate consensus assemblies for bacterial genomes**

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

### 1. Download Test Data

```bash
wget https://.../assembly-curator-test-data.tar.xz  # Todo
tar -xvf assembly-curator-test-data.tar.xz
```

### 2. Pull Docker Image

```bash
docker pull troder/assembler-curator
```

### 3. Run Assembly Curator

```bash
docker run -it --rm \
-v ./data-share:/data:Z \
-v ./plugins-share:/plugins:Z \
-p 8080:8080 \
--name assembly-curator \
troder/assembly-curator
```

### 4. Access the Interface

Open your browser and navigate to:  
`http://localhost:8080`

## 📊 Core Workflow

1. **Preprocessing:** Generate multiple genome assemblies using any assembler.
2. **Curation:** Use the interactive app to cluster, compare, and select contigs.
3. **Export:** Save curated assemblies ready for genome annotation.

## 🔌 Plugin System

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

## 📚 Documentation

- [Trycycler Inspiration](https://github.com/rrwick/Trycycler)
- [gfaviz](https://github.com/ggonnella/gfaviz)
- [Minimap2](https://github.com/lh3/minimap2)
- [Pyskani](https://github.com/althonos/pyskani)
- [Plotly](https://plotly.com/)

## ✨ Future Development

- Enhanced scalability for large datasets.
- Advanced visualization options.
- Community-driven plugin contributions.
