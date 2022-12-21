# content of conftest.py

import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--amber",
        action="store_true",
        default=False,
        help="run tests which call amber code",
    )
    parser.addoption(
        "--chemcloud",
        action="store_true",
        default=False,
        help="run additional chemcloud tests",
    )
    parser.addoption(
        "--full",
        action="store_true",
        default=False,
        help="run full optimization cycle tests",
    )
    parser.addoption(
        "--gpu", action="store_true", default=False, help="run tests which require gpus"
    )
    parser.addoption(
        "--queue",
        action="store_true",
        default=False,
        help="run tests which use the queue",
    )
    parser.addoption("--all", action="store_true", default=False, help="run all tests")


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "amber: test calls amber libraries/code and may be slow to run"
    )
    config.addinivalue_line(
        "markers", "chemcloud: test calls chemcloud and may be slow to run"
    )
    config.addinivalue_line(
        "markers", "full: test does a full optimization cycle and is very slow to run"
    )
    config.addinivalue_line("markers", "gpu: test requires a gpu to run")
    config.addinivalue_line(
        "markers", "queue: test runs on queue and may take long to run"
    )
    config.addinivalue_line("markers", "debug: run only these tests")

def pytest_collection_modifyitems(config, items):
    skips = {}
    skips["amber"] = pytest.mark.skip(reason="need --amber option to run")
    skips["chemcloud"] = pytest.mark.skip(reason="need --chemcloud option to run")
    skips["full"] = pytest.mark.skip(reason="need --full option to run")
    skips["gpu"] = pytest.mark.skip(reason="need --gpu option to run")
    skips["queue"] = pytest.mark.skip(reason="need --queue option to run")
    debug = pytest.mark.skip(reason="skipping unmarked tests in debug mode")
    items2 = []
    for item in items:
        for keyword in skips.keys():
            if (
                not config.getoption(f"--{keyword}")
                and keyword in item.keywords
                and not config.getoption("--all")
            ):
                item.add_marker(skips[keyword])
        if "debug" in item.keywords:
            items2.append(item)
    if len(items2) > 0:
        for item in items:
            if item not in items2:
                item.add_marker(debug)
