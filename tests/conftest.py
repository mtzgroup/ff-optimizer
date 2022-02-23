# content of conftest.py

import pytest


def pytest_addoption(parser):
    parser.addoption(
            "--runexternal", action="store_true", default=False, help="run slow tests"
            )

def pytest_configure(config):
    config.addinivalue_line("markers", "external: test calls external libraries/code and is slow to run")

def pytest_collection_modifyitems(config, items):
    if config.getoption("--runexternal"):
        return
    skip_slow = pytest.mark.skip(reason="need --runexternal option to run")
    for item in items:
        if "external" in item.keywords:
            item.add_marker(skip_slow)
