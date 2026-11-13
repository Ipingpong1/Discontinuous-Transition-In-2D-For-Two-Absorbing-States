TECTONIC ?= tectonic
.PHONY: all pdf
all: pdf
pdf:
	$(TECTONIC) --keep-logs --keep-intermediates main.tex
	$(TECTONIC) --keep-logs --keep-intermediates supplement.tex
	$(TECTONIC) --keep-logs --keep-intermediates main.tex
	$(TECTONIC) --keep-logs --keep-intermediates supplement.tex
