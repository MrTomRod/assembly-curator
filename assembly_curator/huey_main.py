import os
import logging

from huey.constants import WORKER_PROCESS

logging.basicConfig(level=logging.INFO)


def run_huey(
        samples_dir: str = './data-pb',
        plugin_dir: str = './plugins-pb',
        n_workers: int = None
):
    os.environ['HUEY_DB_PATH'] = os.path.join(samples_dir, 'huey.db')
    os.environ['PLUGIN_DIR'] = plugin_dir
    os.environ['MULTIPROCESSING_DOTPLOTS'] = 'FALSE'

    if n_workers is None:
        n_workers = max(os.cpu_count() - 1, 1)

    print('****************************************************')
    print(f"Starting Huey with {n_workers=} for {samples_dir=}")
    print('****************************************************')

    from assembly_curator.huey_tasks import huey
    n_workers = max(os.cpu_count() - 1, 1)
    huey_consumer = huey.create_consumer(
        workers=n_workers,
        worker_type=WORKER_PROCESS,  # Has to be separate processes. Threading is not supported by matplotlib
        periodic=False
    )
    huey_consumer.run()


if __name__ == '__main__':
    from fire import Fire

    Fire(run_huey)
