import os
from collections.abc import Iterable

import numpy as np

from sklearn.preprocessing import MinMaxScaler


def scale(
        x: Iterable[float],
        min_x: float,
        max_x: float
) -> list[float]:
    """
    Scale the input array `x` between `min_x` and `max_x`.

    :param x: Iterable array to scale
    :param min_x: The left boundary of scaled array
    :param max_x: The right boundary of scaled array

    :return: The scaled array `x`
    """

    return list(
        MinMaxScaler(
            feature_range=(min_x, max_x)
        ).fit_transform(
            np.expand_dims(np.array(x), axis=1)
        )[:, 0]
    )


def scale_int(
        x: Iterable[float],
        min_x: float,
        max_x: float
) -> list[float]:
    """
    Scale the input array `x` between `min_x` and `max_x` with additional rounding of output values.

    :param x: Iterable array to scale
    :param min_x: The left boundary of scaled array
    :param max_x: The right boundary of scaled array

    :return: The scaled array `x` with only integer values
    """

    return [
        int(el) for el in list(
            MinMaxScaler(
                feature_range=(min_x, max_x)
            ).fit_transform(
                np.expand_dims(np.array(x), axis=1)
            )[:, 0]
        )
    ]


def is_non_empty(fn: str) -> bool:
    """
    Check if non-empty file exists.

    :param fn: The filepath to the file

    :return: True of non-empty file exists, False otherwise
    """

    return os.path.exists(fn) and (os.stat(fn).st_size != 0)


def style_bool_df(df):
    """
    Style (color) the dataframe consisting of True, False and NaN values.
    """
    return df.style.apply(lambda x: [
        "background-color: green" if v else "background-color: red" if not v else 'background: white' for v in
        x], axis=1)


def save_pickle(f, fn):
    """
    Save object as a pickle file (usually used for dicts).

    f: file object
    fn: file name
    """
    import pickle

    with open(fn, 'wb') as fo:
        pickle.dump(f, fo)


def load_pickle(fn):
    """
    Load object from pickle file (usually used for dicts).

    fn: file name

    return: loaded object
    """
    import pickle

    with open(fn, 'rb') as fo:
        f = pickle.load(fo)

    return f
