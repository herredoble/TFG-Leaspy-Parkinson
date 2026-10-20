# AI4Health 2023 Leaspy Workshop

Longitudinal data consist of the repeated observations of subjects or objects over time. They are ubiquitous in biology and medicine as they inform about the progression of a biological phenomenon such as growth or the progression of a chronic disease.

Analysing longitudinal data requires a careful attention because repeated data of the same subjects are not independent. Linear mixed-effect models have long been a piece of choice to address this problem. They model of the progression of the underlying phenomenon and how it manifests itself in variable forms across subjects.

Recent developments made it possible to address some limitations of these methods. They account for the non-linear dynamics of progression. They do not use age a regressor, so that they can compare subjects data even if the subjects differ in their age at onset or pace of progression.

In this workshop, you will practice using several data longitudinal sets from patients developing neurodegenerative diseases: Alzheimer and Parkinson disease. You will learn how disease course mapping allows you to characterise the variability of the progression profiles across subjects, impute missing data, resample data sets at intermediate time-points, predict the future progression of new patients, and even simulate cohorts of virtual patients.

## Installation

### Get the repository

On the terminal of your computer:
```
git clone https://gitlab.com/icm-institute/aramislab/disease-course-mapping-solutions.git
```

You may also download the files directly from Gitlab.

### Have a Python interpreter ready

If not already done: install a Python interpreter **(version >= 3.7)** for your OS

### Install required packages

Open a terminal on your computer and head into the folder where you just cloned the workshop repository,
and then **head into the `challenges` sub-folder**.

The virtualenv package is required to create virtual environments. You can install it with pip:
```
pip install virtualenv
```

Note: if you are familiar with other tools to create virtual environments you can also use them
([Conda](https://docs.conda.io/en/latest/miniconda.html), [PyEnv](https://github.com/pyenv/pyenv#simple-python-version-management-pyenv), ...)

Create the environment
```
virtualenv leaspy_tutorial --python=python3.8
```

Activate the virtual environment :
_Mac OS / Linux_
```
source leaspy_tutorial/bin/activate
```

_Windows_
```
leaspy_tutorial\Scripts\activate
```

Any python commands you use will now work with your virtual environment.

Install all the needed packages
```
pip install -r requirements.txt
```

Now open `jupyter notebook`:

```
jupyter notebook
```
Or (according to your preferences)
```
jupyter lab
```

(note that you might need to run `deactivate` and re-run `source leaspy_tutorial/bin/activate` after running `pip install -r requirements.txt` if you happen to have a system-wide installation of `jupyter` as well; you can check which installation is used by running `which jupyter`, which should point to `[...]/leaspy_tutorial/bin/jupyter`)

You are ready to start!
You can find information on the content of each part in the following section.

## Workshop Content

The workshop is split into three parts :
 - TP1_LMM [1h30]
 - TP2_leaspy_beginner [1h]
 - TP3_advanced_leaspy [1h30]

Times are indicative and we advice you to make them in the order. You will find below a quick description of each.
You can find the *challenge* version of each part in the `challenges` folder.
We recommand that you start

### TP1_LMM

This TP should give you a better idea of medical data, especially longitudinal ones. You will understand mixed-effects models and get a taste of state-of-the-art techniques.

### TP2_leaspy_beginner

The second TP is here to teach you how to use Leaspy, the package developed by the Aramis team. Do not hesitate to have a look to [its official documentation](https://leaspy.readthedocs.io).

### TP3_advanced_leaspy

The last TP is made for you to test if you well understand Leaspy package and get you into real life problems.
