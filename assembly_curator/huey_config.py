import os
from huey import SqliteHuey


def get_huey(db_path: str = None):
    if not db_path:
        assert 'HUEY_DB_PATH' in os.environ, 'HUEY_DB_PATH environment variable must be set'
        db_path = os.environ['HUEY_DB_PATH']
    return SqliteHuey('preprocessor', filename=db_path)
