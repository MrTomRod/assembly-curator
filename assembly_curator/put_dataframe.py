import pandas as pd
from typing import Union, List, Callable, Any, Sequence, Mapping, Tuple
from pywebio.output import put_datatable, JSFunction, output_register_callback
from pywebio.utils import random_str


def editable_datatable_inplace_update_callback(df: pd.DataFrame):
    """
    Return a callback function that can be used in `put_editable_datatable()` to update the DataFrame inplace.
    All new data will be stored as string.
    """

    def callback(row_id, column, new_value):
        if len(column) == 1:
            df.at[row_id, column[0]] = new_value
        else:
            raise ValueError("Nested fields are not supported with DataFrame")

        print(f"Updated row {row_id} field {column} to {new_value}")
        print(df)

    return callback


def put_dataframe(
        df: pd.DataFrame,
        edit_callback: Callable[[Union[str, int], List[str], 'Any'], None],
        editable_columns: List[str] = None,
        **kwargs
):
    """
    Put an editable datatable.

    :param df: the DataFrame to display in the table
    :param edit_callback: callback function to be called when cell of the table is edited.
        The callback function should accept three parameters: `row_id`, `field_path`, `new_value`,

        - `row_id`: the index of the row being edited in the DataFrame.
        - `field_path`: the field being edited, since the table may be nested, this is a list of field names.
        - `new_value`: the new value of the cell.

    :param editable_columns: list of column names that can be edited. e.g. ['name', 'info', 'age'].
        Leave it to None to make all columns editable
    :param kwargs: other arguments passed to `put_datatable()`
    """
    assert kwargs.get('id_field') is None, \
        'id_field is not allowed in put_editable_datatable()'

    def on_edit(info):
        row_id = int(info['row'])
        field_path = info['path']
        new_value = info['value']
        edit_callback(row_id, field_path, new_value)

    instance_id = kwargs.get('instance_id')
    if instance_id is None:
        instance_id = random_str(10)
        kwargs['instance_id'] = instance_id

    callback_id = output_register_callback(on_edit)
    edit_value_parser = JSFunction('params', """
        if (params.newValue === params.oldValue)
            return params.newValue;
        WebIO.pushData({
            row: params.node.id,
            path: ag_grid_%(instance_id)s.field2path(params.colDef.field),
            value: params.newValue
        }, %(callback_id)r);
        return params.newValue;
    """ % dict(callback_id=callback_id, instance_id=instance_id))

    if editable_columns is None:
        grid_args = kwargs.get('grid_args') or {}
        grid_args.setdefault('defaultColDef', {})
        grid_args['defaultColDef']['editable'] = True
        grid_args['defaultColDef']['valueParser'] = edit_value_parser
        kwargs['grid_args'] = grid_args
    else:
        column_args = kwargs.get('column_args') or {}
        for column in editable_columns:
            column_args.setdefault(column, {})
            column_args[column]['editable'] = True
            column_args[column]['valueParser'] = edit_value_parser
        kwargs['column_args'] = column_args
    return put_datatable(df.to_dict(orient='records'), **kwargs)


def main():
    import urllib.request
    import json

    from pywebio.output import put_button, put_markdown, popup, put_code

    def flatten_json(data):
        def flatten(item, parent_key='', sep='_'):
            flattened = {}
            for k, v in item.items():
                new_key = f"{parent_key}{sep}{k}" if parent_key else k
                if isinstance(v, dict):
                    flattened.update(flatten(v, new_key, sep=sep))
                else:
                    flattened[new_key] = v
            return flattened

        return [flatten(item) for item in data]

    put_markdown("""
    # Editable Datatable Demo
    The table below is editable. The edit will be synced to the original data.
    """)

    with urllib.request.urlopen('https://fakerapi.it/api/v1/persons?_quantity=10') as f:
        data = flatten_json(json.load(f)['data'])

    df = pd.DataFrame(data)
    del data

    put_button("Show Table Data",
               lambda: popup("data of table", put_code(df.to_json(orient='records', indent=2)), size='large'))

    put_dataframe(
        df,
        editable_columns=['firstname', 'lastname', 'email'],
        edit_callback=editable_datatable_inplace_update_callback(df),
    )


if __name__ == '__main__':
    from pywebio import start_server

    start_server(main, port=8080, cdn=False)
