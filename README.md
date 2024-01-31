# EffectorO-1.0

Note: If you have any trouble running this, please email collaborators at:

- Kelsey Wood: **<kjwood@ucdavis.edu>**
- Munir Nur: **<mjnur@ucdavis.edu>**
- Nikko Sanchez: **<niksanchez@ucdavis.edu>**

## About

Please see the [EffectorO paper](https://www.biorxiv.org/content/10.1101/2021.03.19.436227v2) for more information.

## Getting Started (CLI)

Here are the required dependencies for EffectorO:

- Python (versions <= 3.7)
- scikit-learn (versions <= 0.22)
- joblib
- Biopython
- pandas
- NumPy

There are multiple ways to install these dependencies (listed below).

### Using conda

First you have to install conda. If you do not have conda installed on your machine, use this [guide](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) to aid your installation process.

After conda is successfully installed on your machine, create an conda environment using the following command:

```shell
conda env create -f environment_setup/cli-conda-environment.yml
```

Once the environment is created, the `(base)` next to your username should now be `(effectoro-cli)`.

<!-- TODO: see if this works
### Using `requirements.txt`

In order to use requirements.txt, you must have version 3.5-3.8 of Python on your machine. If you do, you can install the dependencies using this command:

```bash
pip install -r environment_setup/requirements.txt
```
-->

## Running EffectorO

After dependencies are installed either on your machine or provided by an environment, you can run EffectorO using this command:

```bash
python3 effectorO.py -i PATH/TO/INPUT_FILE.fasta
```

Once EffectorO is done running, your results can be found in the `effectoro_results` directory.
