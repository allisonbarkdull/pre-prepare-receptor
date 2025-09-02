from setuptools import setup

setup(
    name="pre_prepare_receptor",
    version="0.1.0",
    py_modules=["pre_prepare_receptor"],
    install_requires=[
        "numpy",
        "requests",
        "prody",
        "rdkit-pypi",
        "openmm",
        "openff-toolkit",
        "pdbfixer",
        "parmed",
        "autopath",
        "mdanalysis",
        "deeptime",
        "pyemma"
    ],
    entry_points={
        "console_scripts": [
            "pre_prepare_receptor=pre_prepare_receptor:main",
        ],
    },
    python_requires=">=3.9",
)

