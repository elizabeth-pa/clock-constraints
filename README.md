# Clock constraints on fundamental physics

This project is designed to generate simulated data based on the characteristics of your atomic clocks and use Markov Chain Monte Carlo (MCMC) methods to estimate unknown parameters of various physics models. The results are then used to establish constraints on the proposed theories.

## Project Goals
- Take atomic clock characteristics e.g. as stability and systematic uncertainty as input.
- Generate simulated data for different atomic clocks.
- Inject signals and perform MCMC analysis.
- Estimate unknown parameters and derive constraints for the theories described below.

# How It Works

## Data Generation
We generate simulated data for atomic clocks based on input characteristics such as stability (white noise) and systematic uncertainty (pink noise). The values used by the customers should be put in table stats/clock_pars.csv in the form of h0, and h-1 coefficients, according to the conversion from the Allen variance explained in https://en.wikipedia.org/wiki/Allan_variance. Then pdate the couples of clocks used in the table stats/clock_pars.csv and you will be ready to run the data generation and analysis!
Firstly, you generate in frequency domain the simulated noise, given the two clocks noise contributions considered.

## Signal Injection
Once the data is generated, we inject the desired signal based on the phenomenological model at hand. We consider three different phenomenological models:
1. A Modified-Gravity-inspired sinusoidal time modulation of the relative shift of clocks datastream, with known phase and frequency (that corresponds to the Earth's motion) and unknown amplitude
2. A Dark-Energy-inspired linear time drift of the relative shift of clocks datastream
3. A Dark-Matter-inspired sinusoidal time modulation of the relative shift of clocks datastream, with unknown phase, frequency and amplitude

All these functions are present in stats/data_to_model.py and stats/utils.py

## Parameter estimation techniques
Using the injected data, we apply the Markov Chain Monte Carlo (MCMC) and Fisher methods to forecast the reconstruction of the posterior distribution for the parameters of the phenomenological models considered.
The code is implemented to run for the three aforementioned phenomenological models, but an external user can provide a different model, provided that they add its functional form in stats/data_to_model.py and stats/utils.py as indicated

The output of both the MCMC and Fisher methods are a fair set of samples from the posterior distributions of the phenomenological parameters.

In particular, in the table stats/sigmas.csv there are the forecast width of the posterior distribution in the direction of the amplitude of the signal, for each of the three phenomenological models considered, assuming flat prior in the parameters (and after marginalizing in the frequency and phase parameters for the DM-inspired model)

## Theory Constraints
The estimated parameters are then used to derive constraints for the proposed physics theories. 

# a few words about the models used?

# Installation

To install, clone this repository.  It is recommended to also set up
a python virtual environment to ensure package compatibility. This has been automated with a makefile recipe: simply run
```
> make init
```
and the environment is set up in venv/ and is ready to use.

If a manual setup is preferred, the first step is to initialize a new virtual environment:
```
> python3 -m venv venv
```
This creates a venv/ directory.  Next we can activate the environment
and populate it with the packages from the requirements list:
```
> source venv/bin/activate
> pip install -r requirements.txt
```

# Usage

Begin by activating the virtual environment:
```
> source venv/bin/activate
```
The software is now ready to be used.

A makefile recipe is included to facilitate generating plots, as well as a few other common useful tasks.  To generate the plots, run the command
```
> make -j3
```
The `-j3` flag is optional; it makes the plotting scripts run in parallel.  This first generates a simulated signal using the supplied clock parameters, and then undertakes a signal extraction for the generalized signals associated with dark energy, dark matter, and modified gravity and stores the result in sigmas.csv.  Finally it generates constraint plots for specific theories.  The statistical analysis need only be run once for a given set of clock parameters, thereafter that step can be skipped by using the command
```
> make plots -j3
```
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

# Requirements
- Python v3.13, with the packages pip and venv
- GNU Make v4.4+ (optional, used to automate some tasks)
- All other dependencies are added during the installation process.

# Authors

# Contact

For more information or to collaborate, please contact:

- **Email**: [your.email@example.com](mailto:your.email@example.com)
- **GitHub**: [elizabeth-pa](https://github.com/elizabeth-pa)
- **GitHub**: [bencelder](https://github.com/bencelder)
