.PHONY: all clean submission

ALLFIGS = $(shell grep "^[^%]*includegr" isolation_by_coalescence.tex| sed -e 's/.*{//'| sed -e 's/}.*//')
FIGS = $(shell for x in $(ALLFIGS); do if [ -f $$x.pdf ]; then echo $$x.pdf; else echo $$x.png; fi; done)

all: isolation_by_coalescence.pdf 
	
submission: ibc_main.pdf ibc_supmat.pdf # cover_letter.pdf review-responses.pdf

latex_bundle.tar.gz : 
	tar -cvzhf $@ isolation_by_coalescence.tex references.bib review-response-commands.tex $(FIGS)

isolation_by_coalescence.pdf : references.bib figs/conceptn.pdf figs/all_flat_divergences.fonts.pdf figs/barrier_sample_locations_pretty.fonts.pdf figs/fancy_watershed_assignments.fonts.pdf

ibc_main.pdf : isolation_by_coalescence.pdf 
	pdfjam --outfile $@ $< 1-27

ibc_supmat.pdf : isolation_by_coalescence.pdf 
	pdfjam --outfile $@ $< 28-44

cover_letter.pdf : isolation_by_coalescence.pdf 
	pdfjam --outfile $@ $< 45

review-responses.pdf : isolation_by_coalescence.pdf 
	pdfjam --outfile $@ $< 46-

clean: 
	-rm *.aux *.log *.lof *.lot *.fff *.ttt *.out *.bbl *.blg

%.pdf : %.tex %.bbl
	while ( pdflatex $<;  grep -q "Rerun to get" $*.log ) do true ; done

%.aux : %.tex
	-pdflatex $<

%.bbl : %.aux
	-bibtex $<

%.html : %.md
	Rscript -e "templater::render_template(md.file='$<', output='$@')"

%.svg : %.pdf
	inkscape $< --export-plain-svg=$@

%.png : %.pdf
	convert -density 300 $< -flatten $@

%.pdf : %.ink.svg
	inkscape $< --export-pdf=$@

%.eps : %.pdf
	inkscape --without-gui --export-eps=$@ $<

%.fonts.pdf : %.pdf
	# this will embed fonts properly
	inkscape $< --export-pdf=$@
