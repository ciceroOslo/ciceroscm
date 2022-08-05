"""
Test that all of our modules can be imported
Thanks https://stackoverflow.com/a/25562415/10473080
and openscm-runner
"""
import importlib
import os.path
import pkgutil

import ciceroscm

def import_submodules(package_name):
    package = importlib.import_module(package_name)

    for _, name, is_pkg in pkgutil.walk_packages(package.__path__):
        full_name = package.__name__ + "." + name
        print(full_name)
        importlib.import_module(full_name)
        if is_pkg:
            import_submodules(full_name)


import_submodules("ciceroscm")

# make sure input data etc. are included
ciceroscm_root = os.path.dirname(ciceroscm.__file__)
assert os.path.isfile(
    os.path.join(
        ciceroscm_root, "default_data/IPCC_LUCalbedo.txt"
    )
)
assert os.path.isfile(
    os.path.join(
        ciceroscm_root, "default_data/solar_IPCC.txt"
    )
)
assert os.path.isfile(
    os.path.join(
        ciceroscm_root, "default_data/meanVOLCmnd_ipcc_NH.txt"
    )
)

print(ciceroscm.__version__)
