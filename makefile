# Implements some helpful procedures for generating and cleaning up
# files.
# Requires GNU Make v4.4 or greater

.PHONY: init plots requirements clean plot_MG plot_DE plot_DM stats


# The name of the virtual environment folder,
# and the virtual python and pip executables
VENV = venv
PYTHON = $(VENV)/bin/python
PIP = $(VENV)/bin/pip


# The .WAIT command forces stats to run before plots, even
# when running in parallel with the -j flag.
# This requires make version 4.4 or greater.
all: stats .WAIT plots

init:
	python3.13 -m venv $(VENV)
	$(PIP) install -r requirements.txt

stats:
	@echo "Starting stats..."
	$(PYTHON) stats/data_to_model.py
	@echo "Done with stats"

plot_MG:
	$(PYTHON) plot_MG.py

plot_DE:
	$(PYTHON) plot_DE.py

plot_DM:
	$(PYTHON) plot_DM.py

# Run with the -j3 flag to run these in parallel
plots: plot_MG plot_DE plot_DM
	@echo "Done with plots"

requirements:
	$(PIP) freeze > requirements.txt


clean:
	rm plots/*.png
	#rm stats/sigmas.csv

