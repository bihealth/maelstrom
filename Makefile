.PHONY: default black flake8 test test-v test-vv

default: black flake8

black:
	snakefmt -l 120 .
	black -l 120 .

black-check:
	snakefmt -l 120 --check .
	black -l 120 --check .

flake8:
	flake8 .

test:
	pytest

test-v:
	pytest -v

test-vv:
	pytest -vv
