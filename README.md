# Welcome to Clock Constraints

This project is designed to generate simulated data based on the characteristics of your atomic clocks and use Markov Chain Monte Carlo (MCMC) methods to estimate unknown parameters of various physics models. The results are then used to establish constraints on the proposed theories.

## Project Goals
- Take atomic clock characteristics e.g. as stability and deadtime as input.
- Generate simulated data for atomic clocks.
- Inject signals and perform MCMC analysis.
- Estimate unknown parameters and derive constraints for the theories described below.

# How It Works

## Data Generation
We generate simulated data for atomic clocks based on input characteristics such as stability and deadtime.

## Signal Injection
Once the data is generated, we inject the desired signal based on the theory model at hand.

## MCMC Analysis
Using the injected data, we apply Markov Chain Monte Carlo (MCMC) methods to estimate unknown parameters of various physics models. MCMC generates a posterior distribution of the unknown parameters.

## Theory Constraints
The estimated parameters are then used to derive constraints for the proposed physics theories. 


# Installation

To install, clone this repository.  It is recommended to also set up
a python virtual environment to ensure package compatibility.
To this end, a requirements.txt file is included which lists the
various packages and versions.

For example, from the top of the directory do:
```
> python -m virtualenv venv
> python -m pip install -r requirements.txt
```
This generates the folder venv/ for the virtual environment.  The virtual environment can then be used in the usual way: before running any of the other code, enter
```
> source venv/bin/activate
```
and this will make python use the virtual environment.

# Usage

A makefile is included to facilitate generating plots, as well as a few
other common useful tasks.  To generate the plots, run the command
```
> make plots
```
The plot-generation scripts can also be run individually, if desired, e.g.
```
> python DM-A-vs-omega.py
```

---

```
@misc{sfs,
  author       = {sdv},
  title        = {sdf},
  month        = August,
  year         = 2024,
  publisher    = {sdfs},
  version      = {gdr},
  doi          = {gdfg},
  howpublished = {\url{dfgvd}}
}
```

# Authors


# Contact

For more information or to collaborate, please contact:

- **Email**: [your.email@example.com](mailto:your.email@example.com)
- **GitHub**: [elizabeth-pa](https://github.com/elizabeth-pa)
- **GitHub**: [bencelder](https://github.com/bencelder)
