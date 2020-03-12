import json
import os


def get_from_config(name):
    with open("config.json") as f:
        val = json.load(f)[name]
    return val
