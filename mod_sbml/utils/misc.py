from collections import defaultdict

__author__ = 'anna'


def remove_from_map(mapping, key, value):
    """
    Removes a value value associated with the key key from a dict {key: values} where values is
    a collection of values associated with key.
    :param key: key of interest
    :param value: value to be removed
    :param mapping: dict {key: values}
    :return: void
    """
    if key in mapping and value in mapping[key]:
        mapping[key].remove(value)
        if not mapping[key]:
            del mapping[key]


def invert_map(key2value, factory=set):
    """
    Inverts a dict {key: value} into a dict {value: keys} where keys is a collection
    of type specified in the parameter factory, and contains all keys that mapped to this value in the input dictionary.
    :param key2value: dict {key: value}
    :param factory: type of collection to be used to store keys that map to the same value
    :return: dict {value: keys}
    """
    value2keys = defaultdict(factory)
    for key, value in key2value.items():
        if isinstance(value, list) or isinstance(value, set):
            for v in value:
                value2keys[v].add(key)
        else:
            value2keys[value].add(key)
    return value2keys
