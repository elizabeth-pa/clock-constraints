# Implements some helpful procedures for generating and cleaning up
# files.

.PHONY: init plots requirements clean plot_MG plot_DE plot_DM

# The name of the virtual environment folder,
# and the virtual python and pip executables
VENV = venv
PYTHON = $(VENV)/bin/python
PIP = $(VENV)/bin/pip

init:
	python3 -m venv $(VENV)
	$(VENV)/bin/pip install -r requirements.txt

plot_MG:
	$(PYTHON) plot_MG.py

plot_DE:
	$(PYTHON) plot_DE.py

plot_DM:
	$(PYTHON) plot_DM.py

# Run with the -j3 flag to run these in parallel
plots: plot_MG plot_DE plot_DM
	@echo "Done"

requirements:
	$(PIP) freeze > requirements.txt

clean:
	rm plots/*.png

