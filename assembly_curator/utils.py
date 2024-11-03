# Plugin system inspired by https://gist.github.com/dorneanu/cce1cd6711969d581873a88e0257e312
import os
import logging
import subprocess
import multiprocessing
from email.policy import default
from typing import List, Type, Callable
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


_importers = None
_get_custom_html = None


def load_importers(plugins_dir: str) -> List[Type['AssemblyImporter']]:
    global _importers

    if _importers is not None:
        return _importers

    _importers = []
    for file_name in os.listdir(plugins_dir):
        if file_name.endswith('Importer.py') and not file_name.startswith('.') and not file_name.startswith('__'):
            logging.info(f"Loading plugin: {plugins_dir}/{file_name}")
            module_name = file_name.rstrip('.py')
            module = load_module(os.path.join(plugins_dir, file_name))
            if not hasattr(module, module_name):
                raise ModuleNotFoundError(f"Plugin {module_name} does not have a class with the same name")
            _importers.append(getattr(module, module_name))
    return _importers


def load_get_custom_html(plugins_dir: str) -> Callable:
    global _get_custom_html

    if _get_custom_html is not None:
        return _get_custom_html

    # Construct the path to the get_custom_html.py file
    file_path = os.path.join(plugins_dir, 'get_custom_html.py')

    if not os.path.isfile(file_path):
        logging.info(f"No get_custom_html.py found in {plugins_dir}, skipping.")
        # default: simply write return name of the sample
        return lambda sample: sample

    try:
        # Load the function get_custom_html from the file
        module = load_module(file_path)
        _get_custom_html = module.get_custom_html
    except Exception as e:
        raise ModuleNotFoundError(f"Error loading get_custom_html from {file_path}: {e}")

    return _get_custom_html


def human_bp(bp: int, decimals: int = 1, zero_val='0bp') -> str:
    if bp == 0:
        return zero_val
    magnitude = floor(log(bp, 1000))
    shortened = bp / 1000 ** magnitude
    unit = ['', 'kbp', 'mbp', 'gbp', 'tbp', 'pbp'][magnitude]
    # return f'{shortened:.1f}{unit}'
    return f'{shortened:.{decimals}f}{unit}'


def get_relative_path(overview_html: str, sample_dir: str) -> str:
    # Get the directory of the index.html.jinja2 file
    overview_dir = os.path.dirname(overview_html)
    # Calculate the relative path from the overview directory to the folder
    relative_path = os.path.relpath(sample_dir, start=overview_dir)
    return relative_path


def rgb_array_to_css(rgb_array):
    return f"rgb({', '.join(str(int(value * 255)) for value in rgb_array)})"


def run_command(cmd: str, **kwargs):
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs)
    if result.returncode != 0:
        print(f"Command failed with return code {result.returncode}")
        print(f"stdout:\n{result.stdout.decode()}")
        print(f"stderr:\n{result.stderr.decode()}")
    return result.returncode


def css_escape(s, to_escape='#@+.'):
    for char in to_escape:
        s = s.replace(char, f'\\{char}')
    return s


def detach_process(target, *args, **kwargs):
    process = multiprocessing.Process(target=target, args=args, kwargs=kwargs)
    process.start()
    return process


def invariant_atgc_count(atgc_count: dict):
    """
    GC content is invariant to the orientation of contigs during assembly, but ATGC counts are not.
    If the contig is read in the forward direction, the counts are "inverted" compared to the reverse direction.
    This function ensures that the ATGC counts are consistent regardless of the orientation.

    Args:
        atgc_count (dict): A dictionary with keys 'A', 'T', 'G', 'C' and their respective counts as values.

    Returns:
        dict: The dictionary with the lower counts of 'A' and 'G' when compared to its inverse.
    """
    inverse_atgc = {'A': atgc_count['T'], 'T': atgc_count['A'], 'G': atgc_count['C'], 'C': atgc_count['G']}
    return min(atgc_count, inverse_atgc, key=lambda x: (x['A'], x['G']))
