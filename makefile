# Implements some helpful procedures for generating and cleaning up
# files.

.PHONY: plots
plots:
	python plot_MG.py
	python plot_DE.py
	python plot_DM.py

.PHONY: clean
clean:
	rm plots/*.png

.PHONY: requirements
requirements:
	pip freeze > requirements.txt
