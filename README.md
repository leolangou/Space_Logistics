# Space Logistics

Research code for GEO satellite analysis, historical orbit reconstruction, and topological data analysis.

## Environment Setup

This project uses a dedicated Conda environment with Python 3.11.

### 1. Clone the repository

```bash (PowerShell)
git clone <YOUR-REPOSITORY-URL>
cd Space_Logistics
```

### 2. Create the Conda environment

```bash
conda create -n space-logistics python=3.11 -y
```

Activate it:

```bash
conda activate space-logistics
```

### 3. Install dependencies

```bash
python -m pip install -r requirements.txt
```

### 4. Register the Jupyter kernel

```bash
python -m ipykernel install --user --name space-logistics --display-name "Python (space-logistics)"
```

### 5. Select the environment in VS Code

For Python files:

1. Press `Ctrl+Shift+P`
2. Select `Python: Select Interpreter`
3. Choose `space-logistics (Python 3.11.x)`

For Jupyter notebooks, also select the kernel in the upper-right corner of the notebook:

```text
Python (space-logistics)
```

## Verify the Active Python Environment

Run:

```python
import sys

print(sys.executable)
print(sys.version)
```

The executable path should point to something similar to:

```text
...\miniconda3\envs\space-logistics\python.exe
```

To verify that `pip` is installing into the same environment:

```python
import sys
!{sys.executable} -m pip --version
```

## Installing New Packages

Always activate the project environment first:

```bash
conda activate space-logistics
```

Then install packages using:

```bash
python -m pip install <package-name>
```

Using `python -m pip` helps ensure that the package is installed into the Python interpreter currently being used.

After adding a dependency, update `requirements.txt` as appropriate.

## Space-Track Credentials

Do **not** place Space-Track usernames or passwords directly in source code or commit them to GitHub.

In PowerShell, credentials can be stored temporarily as environment variables:

```powershell
$env:SPACETRACK_USER="your_email"
$env:SPACETRACK_PASSWORD="your_password"
```

Files containing secrets, such as `.env`, should remain excluded through `.gitignore`.

## Recommended Repository Structure

```text
Space_Logistics/
├── .gitignore
├── README.md
├── requirements.txt
├── us_geo_timeseries.py
├── TDA_test.ipynb
├── data/
└── src/
```

## Git / Environment Rules

Do not commit Python environments to Git. The repository should contain the code and the files needed to recreate the environment, not the environment itself.

Recommended rule:

```text
One Git repository
    -> one named Conda environment
    -> one Python version
    -> one requirements.txt
```

For this project:

```text
Repository:   Space_Logistics
Environment:  space-logistics
Python:       3.11
```
