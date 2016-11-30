import os
import shutil

__author__ = 'anna'


def create_dirs(path, clean=True):
    """
    Removes the path if it exists and clean is set to True, and then recreates it.
    :param path: A path to create
    :param clean: Whether to first remove the path if it existed.
    :return: void
    """
    if clean:
        shutil.rmtree(path, True)
    if not os.path.exists(path):
        os.makedirs(path)
