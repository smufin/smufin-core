DOT = $(wildcard *.dot)
PNG = $(DOT:.dot=.png)
PDF = $(DOT:.dot=.pdf)

all: $(PNG) $(PDF)

%.png: %.dot
	dot -Tpng -Gsize=9,15\! -Gdpi=100 < $< > $@

%.pdf: %.dot
	dot -Tpdf < $< > $@

clean:
	rm -f $(PNG) $(PDF)
