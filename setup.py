import os

from setuptools import find_packages, setup
from setuptools.command.test import test as TestCommand

import versioneer

PACKAGE_NAME = "ciceroscm"
AUTHORS = [
    ("Marit Sandstad", "marit.sandstad@cicero.oslo.no"),
    ("Ragnhild Bieltvedt Skeie", "r.b.skeie@cicero.oslo.no"),
    ("Ane Nordlie Johansen", "ane.nordlie.johansen@cicero.oslo.no"),
    ("Benjamin Sanderson", "benjamin.sanderson@cicero.oslo.no"),
]
URL = "https://github.com/ciceroOslo/ciceroscm"

DESCRIPTION = "Python version of the CICERO-SCM simple climate model"
README = "README.rst"

SOURCE_DIR = "src"

REQUIREMENTS = [
    "click",
    "openscm-units>=0.5.0",
    "pyam-iamc",
    "python-dotenv",
    "scmdata>=0.7.4",
    "tqdm",
    "matplotlib>=3.4",
    "openscm_runner",
]

REQUIREMENTS_NOTEBOOKS = [
    "ipywidgets",
    "notebook",
    "seaborn",
]
REQUIREMENTS_TESTS = [
    "codecov",
    "coverage",
    "nbval",
    "pytest-cov",
    "pytest>=4.0",
    "xlrd",
]
REQUIREMENTS_DOCS = ["sphinx>=1.4", "sphinx_rtd_theme", "sphinx-click"]
REQUIREMENTS_DEPLOY = ["twine>=1.11.0", "setuptools>=41.2", "wheel>=0.31.0"]

REQUIREMENTS_DEV = [
    *[
        "bandit",
        "black>=22.3.0",
        "black-nb",
        "flake8",
        "isort>5",
        "mypy",
        "nbdime",
        "pydocstyle",
        "pylint>=2.4.4",
    ],
    *REQUIREMENTS_DEPLOY,
    *REQUIREMENTS_DOCS,
    *REQUIREMENTS_NOTEBOOKS,
    *REQUIREMENTS_TESTS,
]

REQUIREMENTS_EXTRAS = {
    "deploy": REQUIREMENTS_DEPLOY,
    "dev": REQUIREMENTS_DEV,
    "docs": REQUIREMENTS_DOCS,
    "notebooks": REQUIREMENTS_NOTEBOOKS,
    "tests": REQUIREMENTS_TESTS,
}

# no tests/docs in `src` so don't need exclude
PACKAGES = find_packages(SOURCE_DIR)
PACKAGE_DIR = {"": SOURCE_DIR}
PACKAGE_DATA = {"ciceroscm": [os.path.join("default_data", "*.txt")]}

# Get the long description from the README file
with open(README, "r") as f:
    README_LINES = ["ciceroscm", "==============", ""]
    add_line = False
    for line in f:
        if line.strip() == ".. sec-begin-long-description":
            add_line = True
        elif line.strip() == ".. sec-end-long-description":
            break
        elif add_line:
            README_LINES.append(line.strip())

if len(README_LINES) < 3:
    raise RuntimeError("Insufficient description given")


class ciceroscm(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest

        pytest.main(self.test_args)


cmdclass = versioneer.get_cmdclass()
cmdclass.update({"test": ciceroscm})

setup(
    name=PACKAGE_NAME,
    version=versioneer.get_version(),
    description=DESCRIPTION,
    long_description="\n".join(README_LINES),
    long_description_content_type="text/x-rst",
    author=", ".join([author[0] for author in AUTHORS]),
    author_email=", ".join([author[1] for author in AUTHORS]),
    url=URL,
    license="Apache 2.0",
    classifiers=[  # full list at https://pypi.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: Apache Software License",
        "Intended Audience :: Developers",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    keywords=["cicero-scm", "cicero", "python", "repo", "simple", "climate", "model"],
    packages=PACKAGES,
    package_dir=PACKAGE_DIR,
    package_data=PACKAGE_DATA,
    include_package_data=True,
    install_requires=REQUIREMENTS,
    extras_require=REQUIREMENTS_EXTRAS,
    cmdclass=cmdclass,
)
