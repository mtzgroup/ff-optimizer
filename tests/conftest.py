# content of conftest.py

import pytest


def pytest_addoption(parser):
    parser.addoption(
            "--runexternal", action="store_true", default=False, help="run slow tests"
            )
    parser.addoption(
            "--runtccloud", action="store_true", default=False, help="run additional tccloud tests"
            )

def pytest_configure(config):
    config.addinivalue_line("markers", "external: test calls external libraries/code and is slow to run")
    config.addinivalue_line("markers", "tccloud: test calls tccloud and may be slow to run")

def pytest_collection_modifyitems(config, items):
    print(config.getoption("--runexternal"))
    skip_slow = pytest.mark.skip(reason="need --runexternal option to run")
    for item in items:
        if "external" in item.keywords and not config.getoption("--runexternal"):
            if not "tccloud" in item.keywords or not config.getoption("--runtccloud"):
                item.add_marker(skip_slow)
