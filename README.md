# Clock constraints on fundamental physics

This software generates simulated data based on the characteristics of atomic clocks and uses Markov Chain Monte Carlo (MCMC) methods to estimate unknown parameters of various physics models. The results are then used to establish constraints on the proposed theories.

## Project Goals
- Take atomic clock characteristics e.g. as stability and systematic uncertainty as input.
- Generate simulated data for different atomic clocks.
- Inject signals and perform MCMC analysis.
- Estimate unknown parameters and derive constraints for the theories described below.

# How It Works

## Data Generation
We generate simulated data for atomic clocks based on input characteristics such as stability (white noise) and systematic uncertainty (pink noise). The values are supplied via the table stats/clocks_parameters.csv in the form of h0, and h-1 coefficients, according to the conversion from the Allan variance explained in https://en.wikipedia.org/wiki/Allan_variance.
The pairs of clocks that are to be considered are listed in stats/clocks_pairs.csv and you will be ready to run the data generation and analysis!
Firstly, you generate in frequency domain the simulated noise, given the two clocks noise contributions considered.

## Signal Injection
Once the data is generated, we inject the desired signal based on the phenomenological model at hand. We consider three different phenomenological models:
1. A Modified-Gravity-inspired sinusoidal time modulation of the relative shift of clocks datastream, with known phase and frequency (that corresponds to the Earth's motion) and unknown amplitude
2. A Dark-Energy-inspired linear time drift of the relative shift of clocks datastream
3. A Dark-Matter-inspired sinusoidal time modulation of the relative shift of clocks datastream, with unknown phase, frequency and amplitude

All these functions are present in stats/data_to_model.py and stats/utils.py

## Parameter estimation techniques
Using the injected data, we apply the Markov Chain Monte Carlo (MCMC) and Fisher methods to forecast the reconstruction of the posterior distribution for the parameters of the phenomenological models considered.
The code is implemented to run for the three aforementioned phenomenological models, but an external user can provide a different model, provided that they add its functional form in stats/data_to_model.py and stats/utils.py as indicated.

The output of both the MCMC and Fisher methods are a fair set of samples from the posterior distributions of the phenomenological parameters.

In particular, in the table stats/sigmas.csv there are the forecast width of the posterior distribution in the direction of the amplitude of the signal, for each of the three phenomenological models considered, assuming flat prior in the parameters (and after marginalizing in the frequency and phase parameters for the DM-inspired model)

## Theory Constraints
The estimated parameters are then used to derive constraints for the proposed physics theories.  The basic signals that are searched for in the data are 1) linear drift, corresponding to dark energy 2) yearly oscillations, corresponding to modified gravity via the Sun-clock interaction and 3) oscillations with arbitrary frequency, corresponding to dark matter, where the oscillation frequency is set by the mass of the dark matter particle.

Those signal amplitudes are computed and stored in stats/sigma_A_table.csv.  For ease of use, there is also an API to simplify programmatic queries of these signals.  This is supplied in mu_constraints.py and is detailed there, with one API call for each of the above three genreral signals.

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
The `-j3` flag is optional; it allows the plotting scripts to run in parallel.

This command runs the full data analyisis forecast which consists of the following steps:
1) A simulated realization of the noise datastream is generated using the clock parameters (supplied in `stats/clocks_pars.csv`) for each clock pair specified in `stats/clocks_couples.csv`.
2) A realization of the forecast signal from the generalized theoretical models associated with dark energy, dark matter, and modified gravity is added to the datastream. The models considered are reported in `stats/utils.py` and if the user wishes to forecast a different one, the structure of the forecast allows them to do so.
The data analysis runs a Markov Chain Monte Carlo (MCMC) Bayesian analysis on the parameters of the models considered, producing samples distributed as the posterior distribution of the parameters.
3) The statistics of the posterior distribution is obtained, outputting the forecast posterior uncertainty of the overall amplitude of the simulated signal (which is one of the theoretical parameters employed). The analysis is run for each of the three phenomenological models and each of the provided clock pairs. The numerical values of the forecast uncertainties are stored in `stats/sigmas.csv`.
4) Constraint plots for specific dark energy, dark matter, and modified gravity theories are generated and stored in the `plots/` directory.

The statistical analysis need only be run once for a given set of clock parameters.  Thereafter it can be skipped by using the command
```
> make plots -j3
```
The plot-generation scripts can also be run individually if desired via e.g.
```
> python plot_MG.py
```

# Adding a clock pair
TO DO: walk through an example of adding a new clock pair.  Which files to edit, what to run, what is returned.

# Clock pair API
TO DO: describe the interface between the stats and theory sections.

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
