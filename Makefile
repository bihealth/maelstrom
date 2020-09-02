.PHONY: default
default: black flake8 maelstrom

.PHONY: black
black:
	snakefmt -l 100 .
	black -l 100 .

.PHONY: black-check
black-check:
	snakefmt -l 100 --check .
	black -l 100 --check .

.PHONY: flake8
flake8:
	flake8 .

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
