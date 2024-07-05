from assembler_tools.AssemblyImporter import AssemblyImporter
from assembler_tools.Assembly import Assembly


class PacbioImporter(AssemblyImporter):
    def __init__(self):
        pass

    def load_assembly(self, sample_dir: str) -> Assembly:
        return "PacbioImporter"
