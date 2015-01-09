from collections import defaultdict

__author__ = 'anna'


def remove_from_map(mp, key, value):
    if key in mp:
        mp[key] -= {value}
        if not mp[key]:
            del mp[key]


def invert_map(key2value):
    value2keys = defaultdict(set)
    for key, value in key2value.iteritems():
        if isinstance(value, list) or isinstance(value, set):
            for v in value:
                value2keys[v].add(key)
        else:
            value2keys[value].add(key)
    return value2keys