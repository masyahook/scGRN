import os

import numpy as np

from sklearn.preprocessing import MinMaxScaler


scale = lambda x, min_y, max_y: list(MinMaxScaler(feature_range=(min_y, max_y)).fit_transform(np.expand_dims(np.array(x), axis=1))[:, 0])
scale_int = lambda x, min_y, max_y: [int(el) for el in list(MinMaxScaler(feature_range=(min_y, max_y)).fit_transform(np.expand_dims(np.array(x), axis=1))[:, 0])]


def is_non_empty(fn):
    return os.path.exists(fn) and (os.stat(fn).st_size != 0)


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