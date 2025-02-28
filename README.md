# Clock constraints on fundamental physics

This project is designed to generate simulated data based on the characteristics of your atomic clocks and use Markov Chain Monte Carlo (MCMC) methods to estimate unknown parameters of various physics models. The results are then used to establish constraints on the proposed theories.

## Project Goals
- Take atomic clock characteristics e.g. as stability and systematic uncertainty as input.
- Generate simulated data for different atomic clocks.
- Inject signals and perform MCMC analysis.
- Estimate unknown parameters and derive constraints for the theories described below.

# How It Works

## Data Generation
We generate simulated data for atomic clocks based on input characteristics such as stability (white noise) and systematic uncertainty (pink noise). 

## Signal Injection
Once the data is generated, we inject the desired signal based on the theory model at hand.

## MCMC Analysis
Using the injected data, we apply the Markov Chain Monte Carlo (MCMC) method to estimate unknown parameters of various physics models. MCMC generates a posterior distribution of the unknown parameters.

## Theory Constraints
The estimated parameters are then used to derive constraints for the proposed physics theories. 

# a few words about the models used?

# Installation

To install, clone this repository.  It is recommended to also set up
a python virtual environment to ensure package compatibility. The first step is to initialize a virtual environment:
```
> python3 -m venv venv
```
This creates a venv/ directory.  Next we can activate the environment
and populate it with the packages from the requirements list:
```
> source venv/bin/activate
> pip install -r requirements.txt
```
Optionally, one could instead use a makefile recipe automating these steps which has been included for convenience:
```
> make init
```

# Usage

Begin by activating the virtual environment:
```
> source venv/bin/activate
```
The software is now ready to be used.

A makefile recipe is included to facilitate generating plots, as well as a few other common useful tasks.  To generate the plots, run the command
```
> make plots -j3
```
The `j3` flag is optional; it makes the plotting scripts run in parallel.
The plot-generation scripts can also be run individually, if desired, e.g.
```
> python plot_MG.py
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
