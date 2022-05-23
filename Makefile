.DEFAULT_GOAL := help

VENV_DIR ?= venv

TESTS_DIR=./tests
NOTEBOOKS_DIR=./notebooks
NOTEBOOKS_SANITIZE_FILE=$(TESTS_DIR)/notebook-tests.cfg

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

.PHONY: help
help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

.PHONY: checks
checks: $(VENV_DIR)  ## run all the checks
	@echo "=== bandit ==="; $(VENV_DIR)/bin/bandit -c .bandit.yml -r src || echo "--- bandit failed ---" >&2; \
		echo "\n\n=== black ==="; $(VENV_DIR)/bin/black --check src tests setup.py  || echo "--- black failed ---" >&2; \
		echo "\n\n=== flake8 ==="; $(VENV_DIR)/bin/flake8 src tests setup.py || echo "--- flake8 failed ---" >&2; \
		echo "\n\n=== isort ==="; $(VENV_DIR)/bin/isort --check-only --quiet src tests setup.py || echo "--- isort failed ---" >&2; \
		echo "\n\n=== pydocstyle ==="; $(VENV_DIR)/bin/pydocstyle src || echo "--- pydocstyle failed ---" >&2; \
		echo "\n\n=== pylint ==="; $(VENV_DIR)/bin/pylint src || echo "--- pylint failed ---" >&2; \
		echo "\n\n=== notebook tests ==="; $(VENV_DIR)/bin/pytest notebooks -r a --nbval --sanitize-with $(NOTEBOOKS_SANITIZE_FILE) || echo "--- notebook tests failed ---" >&2; \
		echo "\n\n=== tests ==="; $(VENV_DIR)/bin/pytest tests -r a --cov=ciceroscm --cov-report='' \
			&& $(VENV_DIR)/bin/coverage report --fail-under=95 || echo "--- tests failed ---" >&2; \
		echo

format-checks: $(VENV_DIR)  ## run all the checks
	@echo "=== bandit ==="; $(VENV_DIR)/bin/bandit -c .bandit.yml -r src || echo "--- bandit failed ---" >&2; \
		echo "\n\n=== black ==="; $(VENV_DIR)/bin/black --check src tests setup.py  || echo "--- black failed ---" >&2; \
		echo "\n\n=== flake8 ==="; $(VENV_DIR)/bin/flake8 src tests setup.py || echo "--- flake8 failed ---" >&2; \
		echo "\n\n=== isort ==="; $(VENV_DIR)/bin/isort --check-only --quiet src tests setup.py || echo "--- isort failed ---" >&2; \
		echo "\n\n=== pydocstyle ==="; $(VENV_DIR)/bin/pydocstyle src || echo "--- pydocstyle failed ---" >&2; \
		echo "\n\n=== pylint ==="; $(VENV_DIR)/bin/pylint src || echo "--- pylint failed ---" >&2; \
		echo


.PHONY: format
format:  ## re-format files
	make isort
	make black

.PHONY: format-notebooks
format-notebooks: $(VENV_DIR)  ## format the notebooks
	@status=$$(git status --porcelain $(NOTEBOOKS_DIR)); \
	if test ${FORCE} || test "x$${status}" = x; then \
		$(VENV_DIR)/bin/black-nb $(NOTEBOOKS_DIR); \
	else \
		echo Not trying any formatting. Working directory is dirty ... >&2; \
	fi;

black: $(VENV_DIR)  ## apply black formatter to source and tests
	@status=$$(git status --porcelain src tests docs scripts); \
	if test "x$${status}" = x; then \
		$(VENV_DIR)/bin/black --exclude _version.py setup.py src tests docs/source/conf.py scripts/*.py; \
	else \
		echo Not trying any formatting. Working directory is dirty ... >&2; \
	fi;

isort: $(VENV_DIR)  ## format the code
	@status=$$(git status --porcelain src tests); \
	if test "x$${status}" = x; then \
		$(VENV_DIR)/bin/isort src tests setup.py; \
	else \
		echo Not trying any formatting. Working directory is dirty ... >&2; \
	fi;

docs: $(VENV_DIR)  ## build the docs
	$(VENV_DIR)/bin/sphinx-build -M html docs/source docs/build

test:  $(VENV_DIR) ## run the full testsuite
	$(VENV_DIR)/bin/pytest tests --cov -rfsxEX --cov-report term-missing

test-pypi-install: $(VENV_DIR)  ## test whether installing from PyPI works
	$(eval TEMPVENV := $(shell mktemp -d))
	python3 -m venv $(TEMPVENV)
	$(TEMPVENV)/bin/pip install pip wheel --upgrade
	$(TEMPVENV)/bin/pip install ciceroscm --pre
	$(TEMPVENV)/bin/python scripts/test_install.py

test-install: $(VENV_DIR)  ## test installing works
	$(eval TEMPVENV := $(shell mktemp -d))
	python3 -m venv $(TEMPVENV)
	$(TEMPVENV)/bin/pip install pip wheel --upgrade
	$(TEMPVENV)/bin/pip install .
	$(TEMPVENV)/bin/python scripts/test_install.py

virtual-environment: $(VENV_DIR)  ## update venv, create a new venv if it doesn't exist make
	echo "If you want this to be rerun, run make clean first"
$(VENV_DIR): setup.py setup.cfg
	[ -d $(VENV_DIR) ] || python3 -m venv $(VENV_DIR)
	$(VENV_DIR)/bin/pip install --upgrade pip wheel
	$(VENV_DIR)/bin/pip install -e .[dev]
	$(VENV_DIR)/bin/jupyter nbextension enable --py widgetsnbextension


	touch $(VENV_DIR)
clean: $(VENV_DIR)
	touch setup.py

first-venv: ## create a new virtual environment for the very first repo setup
	python3 -m venv $(VENV_DIR)

	$(VENV_DIR)/bin/pip install --upgrade pip
	$(VENV_DIR)/bin/pip install versioneer
	# don't touch here as we don't want this venv to persist anyway


