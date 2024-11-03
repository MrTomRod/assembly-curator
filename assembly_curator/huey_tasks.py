import os
import dill

from assembly_curator.huey_config import get_huey
from assembly_curator.main_base import process_sample

huey = get_huey()

_importers = None


def load_importers(plugins_dir: str = None):
    global _importers

    if plugins_dir is None:
        plugins_dir = os.environ['PLUGIN_DIR']

    if _importers is None:
        from assembly_curator.utils import load_importers
        _importers = load_importers(plugins_dir)
    return _importers


@huey.task()
def process_assembly(sample, sample_dir):
    print(f'>>>>>>>>>>>>>>> processing assembly {sample}')
    importers = load_importers()
    assemblies = process_sample(sample, sample_dir, importers)

    if not assemblies:
        # touch /assembly-curator/fail
        with open(f"{sample_dir}/assembly-curator/failed", 'w') as f:
            f.write('Assembly failed.')
        return

    dill.dump(assemblies, open(f"{sample_dir}/assembly-curator/assemblies.pkl", 'wb'))
    # delete processingif it exists
    if os.path.isfile(f"{sample_dir}/processing"):
        os.remove(f"{sample_dir}/processing")
