import os
import shutil

__author__ = 'anna'


def create_dirs(path, clean=True):
    if clean:
        shutil.rmtree(path, True)
    if not os.path.exists(path):
        os.makedirs(path)