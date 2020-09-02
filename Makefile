.PHONY: default
default: black pylama mypy

.PHONY: black
black:
	snakefmt -l 100 .
	black -l 100 .

.PHONY: black-check
black-check:
	snakefmt -l 100 --check .
	black -l 100 --check .

.PHONY: pylama
flake8:
	pylama

.PHONY: test
test:
	pytest

.PHONY: test-v
test-v:
	pytest -v

.PHONY: test-vv
test-vv:
	pytest -vv

.PHONY: mypy
mypy:
	mypy maelstrom
