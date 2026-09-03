NAME = metacellsgraphspy

MAX_LINE_LENGTH = 120

ALL_SOURCE_FILES = $(shell git ls-files)

PY_SOURCE_FILES = $(filter %.py, $(ALL_SOURCE_FILES))

RST_SOURCE_FILES = $(filter %.rst, $(ALL_SOURCE_FILES))

DOCS_SOURCE_FILES = $(filter-out docs/v0.1.0/%, $(filter docs/%, $(ALL_SOURCE_FILES)))

.DEFAULT_GOAL := help

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help.replace('TODO-', 'TODO')))
endef
export PRINT_HELP_PYSCRIPT

help:
	@python3 -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

clean: clean-make clean-build clean-pyc clean-test clean-docs  ## remove all build, test, coverage and Python artifacts

clean-make:
	rm -fr .make.*

clean-build:
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc:
	find . -name .mypy_cache -exec rm -fr {} +
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test:
	rm -f .coverage
	rm -fr htmlcov/
	rm -fr .pytest_cache

clean-docs:
	rm -fr docs/v0.1.0

TODO = todo$()x

pc: $(TODO) ci staged  ## check everything before commit

ci: format smells pytest docs dist  ## check everything in a CI server

staged:  ## check everything is staged for git commit
	@if git status . | grep -q 'Changes not staged\|Untracked files'; \
	then \
	    git status; \
	    echo 'There are unstaged changes (run `git add .`).'; \
	    false; \
	else true; \
	fi

format: trailingspaces linebreaks backticks fstrings isort black flake8 ## check code format

trailingspaces: .make.trailingspaces  ## check for trailing spaces

REAL_SOURCE_FILES = \
    $(filter-out %.png, \
    $(filter-out %.svg, \
    $(filter-out docs/v0.1.0/%, \
    $(ALL_SOURCE_FILES))))

.make.trailingspaces: $(REAL_SOURCE_FILES)
	@echo "trailingspaces"
	@if grep -Hn '\s$$' $(REAL_SOURCE_FILES); \
	then \
	    echo 'Files contain trailing spaces (run `make reformat` or `make stripspaces`).'; \
	    false; \
	else true; \
	fi
	touch $@

linebreaks: .make.linebreaks  ## check line breaks in Python code

.make.linebreaks: $(PY_SOURCE_FILES)
	@echo "linebreaks"
	@if grep -Hn "[^=*][^][/<>\"'a-zA-Z0-9_,:()#}{.?!\\=\`+-]$$" $(PY_SOURCE_FILES) | grep -v -- '--$$\|import \*$$'; \
	then \
	    echo 'Files wrap lines after instead of before an operator (fix manually).'; \
	    false; \
	fi
	touch $@

backticks: .make.backticks  ## check usage of backticks in documentation

.make.backticks: $(PY_SOURCE_FILES) $(RST_SOURCE_FILES)
	@echo "backticks"
	@OK=true; \
	for FILE in $(PY_SOURCE_FILES) $(RST_SOURCE_FILES); \
	do \
	    if ( sed 's/`genindex`/genindex/;s/``\([^`]*\)``/\1/g;s/:[^ :]*:`[^`]*`/_/g;s/`\([^`]*\)`_/\1_/g' "$$FILE" \
               | grep --label "$$FILE" -n -H '`' \
               ) \
	    then OK=false; \
	    fi; \
	done; \
	if $$OK; \
	then true; \
	else \
	    echo 'Documentation contains invalid ` markers (fix manually).'; \
	    false; \
	fi
	touch $@

fstrings: .make.fstrings  ## check f-strings in Python code

.make.fstrings: $(PY_SOURCE_FILES)
	@echo "fstrings"
	@if grep -Hn '^[^"]*\("\([^"]\|\\"\)*"[^"]*\)*[^f]"\([^"]\|\\"\)*{' $(PY_SOURCE_FILES) | grep -v 'NOT F-STRING'; \
	then \
	    echo 'Strings appear to be f-strings, but are not (fix manually).'; \
	    false; \
	fi
	touch $@

isort: .make.isort  ## check imports with isort

.make.isort: $(PY_SOURCE_FILES)
	isort --line-length $(MAX_LINE_LENGTH) --force-single-line-imports --check $(NAME) tests
	touch $@

$(TODO): .make.$(TODO)  ## check there are no leftover TODO-X

.make.$(TODO): $(REAL_SOURCE_FILES)
	@echo 'grep -n -i $(TODO) `git ls-files`'
	@if grep -n -i $(TODO) `git ls-files`; \
	then \
	    echo "Files contain $(TODO) markers (fix manually)."; \
	    false; \
	else true; \
	fi
	touch $@

black: .make.black  ## check format with black

.make.black: $(PY_SOURCE_FILES)
	black --line-length $(MAX_LINE_LENGTH) --check $(NAME) tests
	touch $@

flake8: .make.flake8  ## check format with flake8

.make.flake8: $(PY_SOURCE_FILES)
	flake8 --max-line-length $(MAX_LINE_LENGTH) --ignore F401,E402,F403,W503,E704,E722,E501 $(NAME) tests
	touch $@

reformat: stripspaces isortify blackify  ## reformat code

stripspaces:  # strip trailing spaces
	@echo stripspaces
	@for FILE in $$(grep -l '\s$$' $$(git ls-files | grep -v docs/v)); \
	do sed -i -s 's/\s\s*$$//' $$FILE; \
	done

isortify:  ## sort imports with isort
	isort --line-length $(MAX_LINE_LENGTH) --force-single-line-imports $(NAME) tests

blackify:  ## reformat with black
	black --line-length $(MAX_LINE_LENGTH) $(NAME) tests

smells: mypy pylint  ## check for code smells

pylint: .make.pylint  ## check code with pylint

.make.pylint: $(PY_SOURCE_FILES)
	pylint --good-names=X --disable=fixme,protected-access,wrong-import-position,no-name-in-module,too-many-arguments,too-many-locals,too-many-public-methods,too-few-public-methods,bare-except,line-too-long $(NAME) tests
	touch $@

mypy: .make.mypy  ## check code with mypy

.make.mypy: $(PY_SOURCE_FILES)
	mypy $(NAME) tests
	touch $@

pytest: .make.pytest  ## run tests on the active Python with pytest

# Julia runs a single thread here because juliacall segfaults at interpreter exit under pytest when Julia has several
# threads. This is unrelated to the tests themselves, which all pass before the crash; see SG_PLAN.md for the analysis.
.make.pytest: $(PY_SOURCE_FILES)
	PYTHON_JULIACALL_THREADS=1 `which pytest` -vv -s --cov=$(NAME) --cov-report=html --cov-report=term --no-cov-on-fail --maxfail=1 tests
	touch $@

.PHONY: docs
docs: .make.docs  ## generate HTML documentation

.make.docs: $(DOCS_SOURCE_FILES) $(PY_SOURCE_FILES) $(RST_SOURCE_FILES)
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	sed -i 's/</\n</g' docs/v0.1.0/html/*.html
	@echo "Results in docs/v0.1.0/html/index.html"
	touch $@

dist: .make.dist  ## builds the release distribution package

.make.dist: staged $(ALL_SOURCE_FILES)
	rm -rf dist/
	python3 setup.py sdist
	twine check dist/*
	touch $@

tags: $(PY_SOURCE_FILES)  ## generate a tags file for vi
	ctags $(PY_SOURCE_FILES)
