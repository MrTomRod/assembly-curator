# Plugin system inspired by https://gist.github.com/dorneanu/cce1cd6711969d581873a88e0257e312
import os
import logging
from typing import List, Type
from importlib import util
from math import log, floor
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .AssemblyImporter import AssemblyImporter


class AssemblyFailedException(Exception):
    """
    Raised when the assembler failed to produce an assembly
    """

    def __init__(self, message, severity: str = 'danger'):
        super().__init__(message)
        self.severity = severity


class AssemblyImportError(Exception):
    """
    Raised when the assembly import process fails, i.e. a bug
    """
    pass


class MinorAssemblyException(Exception):
    """
    Raised when minor problems arise that should not prevent the assembly from being processed
    """
    pass


def load_module(path):
    name = os.path.split(path)[-1]
    spec = util.spec_from_file_location(name, path)
    module = util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_plugins(dir_path: str) -> List[Type['AssemblyImporter']]:
    plugins = []
    for file_name in os.listdir(dir_path):
        if file_name.endswith('.py') and not file_name.startswith('.') and not file_name.startswith('__'):
            logging.info(f"Loading plugin: {dir_path}/{file_name}")
            module_name = file_name.rstrip('.py')
            module = load_module(os.path.join(dir_path, file_name))
            assert hasattr(module, module_name), f"Plugin {module_name} does not have a class with the same name"
            plugins.append(getattr(module, module_name))
    return plugins


def human_bp(bp: int, decimals: int = 1, zero_val='0bp') -> str:
    if bp == 0:
        return zero_val
    magnitude = floor(log(bp, 1000))
    shortened = bp / 1000 ** magnitude
    unit = ['', 'kbp', 'mbp', 'gbp', 'tbp', 'pbp'][magnitude]
    # return f'{shortened:.1f}{unit}'
    return f'{shortened:.{decimals}f}{unit}'


def get_relative_path(overview_html: str, sample_dir: str) -> str:
    # Get the directory of the overview.html file
    overview_dir = os.path.dirname(overview_html)
    # Calculate the relative path from the overview directory to the folder
    relative_path = os.path.relpath(sample_dir, start=overview_dir)
    return relative_path
