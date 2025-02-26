# Implements some helpful procedures for generating and cleaning up
# files.

.PHONY: plots
plots:
	# The & makes them run in parallel
	python MG-plots.py
	python DE-plots.py
	python DM-plots.py

.PHONY: clean
clean:
	rm plots/*.png

.PHONY: requirements
requirements:
	pip freeze > requirements.txt
