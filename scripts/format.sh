#!/bin/sh -e

set -x

autoflake --remove-all-unused-imports --recursive --remove-unused-variables --in-place scripts ff_optimizer tests --exclude=__init__.py
black scripts ff_optimizer tests
isort scripts ff_optimizer tests
