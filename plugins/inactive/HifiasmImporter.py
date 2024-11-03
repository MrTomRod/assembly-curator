from assembly_curator.AssemblyImporter import AssemblyImporter
from assembly_curator.Assembly import Assembly


class HifiasmImporter(AssemblyImporter):
    def __init__(self):
        pass

    def load_assembly(self, sample_dir: str) -> Assembly:
        return "HifiasmImporter"
