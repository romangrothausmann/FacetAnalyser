### create distribution plot of inter-planar angles

.PHONY: all

all : 8faced-rhombic-dodecahedron_twinned_010_deformed_wss_FA.svg
all : 8faced-rhombic-dodecahedron_twinned_simp+def_01_04c_l01.png

# %_FA.vtp %_FA-hull.vtp %_FA-edges.vtp : %.vtp
# 	FacetAnalyserCLI $< 51 15 0.2 0.001 2 $*_FA.vtp $*_FA-hull.vtp $*_FA-edges.vtp # do not use $@ as this varies dep. on what Make demands first

%_FA.vtp : ../../../demos/%_FA.vtp # use output of basicExampleTest01
	ln -s $< $@

%_FA.csv : %_FA.vtp
	/opt/VTK-CLIs/bin/vtpFD2csv $< | datamash -t, --no-strict --filler="" transpose > $@ # vtpFD2csv from VTK-CLIs

%.svg : %.csv
	/opt/gnuplot/bin/gnuplot -e "datafile='$<'; outfile='$@'; bin=10; sigma=5; xmin=0; xmax=180; col='interplanarAngles'; sep=','; xlabel='interplanar angle (IPA)'" ~/gnuplot/hist-ar+kd.gp # gp-5.x allows referencing columns by header-name ;-)
	@test -f $@ # gnuplot' exit dose not allow specifying a return value

define BASE_TEX
\documentclass{article}
\usepackage{svg}
\pagenumbering{gobble}

\begin{document}
\includesvg[width=\textwidth, pretex=\large]{SVGFILE}
\end{document}
endef
export BASE_TEX

%_.tex : %.svg # _ to avoid collision with PDF created by inkscape
	$(file >$@,$(subst SVGFILE,$(basename $<),${BASE_TEX}))

%_.pdf : %_.tex
	pdflatex --shell-escape $<

%.png : %_.pdf
	convert -define png:exclude-chunks=date,time +set date:create +set date:modify  -density 150 $< -trim $@
