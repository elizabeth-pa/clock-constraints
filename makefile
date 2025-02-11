# Implements some helpful procedures for generating and cleaning up
# files.

.PHONY: plots
plots:
	python make_plots.py
	python DM-A-vs-omega.py

.PHONY: clean
clean:
	rm plots/*.png

.PHONY: requirements
requirements:
	pip freeze > requirements.txt
