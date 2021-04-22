# Bioinformatics Scripts

This is a collection of scripts for analyzing different forms of bioinformatics data.

These scripts require Python version >= 3.6, as well as some package dependences. To set up a virtual environment and install the dependencies, run:

```bash
python3 -m venv env
source env/bin/activate
python3 -m pip install -r requirements.txt
```

The virtual environment must be activated before running the Python scripts, as well.

Additionally, make sure to set the environment variable `PYTHONPATH` to the root of this directory:

```bash
export PYTHONPATH=/path/to/scripts
```

This is necessary for importing modules between scripts.