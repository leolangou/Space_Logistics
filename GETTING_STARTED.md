# Getting Started

## 1. Clone the Repository

Open VS Code, then open a PowerShell terminal:

```powershell
git clone <REPOSITORY-URL>
cd Space_Logistics
```

Then open the repository folder in VS Code:

```powershell
code .
```

## 2. Create the Python Environment

Create the project Conda environment:

```powershell
conda create -n space-logistics python=3.11 -y
conda activate space-logistics
```

Install the project requirements:

```powershell
python -m pip install -r requirements.txt
```

## 3. Select the Environment in VS Code

Press:

```text
Ctrl + Shift + P
```

Select:

```text
Python: Select Interpreter
```

and choose the `space-logistics` Python 3.11 environment.

For Jupyter notebooks, also select `space-logistics` as the notebook kernel in the upper-right corner.

## 4. Run the Code

Once the environment is selected and the requirements are installed, run the notebooks/scripts normally in VS Code.

For the satellite time-series workflow, use:

```text
us_geo_timeseries.py
```

For the current TDA analysis, use:

```text
TDA_test.ipynb
```

## Updating the Repository

Before starting work:

```powershell
git pull
```

After making changes:

```powershell
git add .
git commit -m "Describe your changes"
git push
```

Do not commit credentials, virtual environments, or large raw datasets.
