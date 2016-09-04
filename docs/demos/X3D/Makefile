
### convert X3D to XHTML (X3DOM/WebGL) with XSLT (also extract X3D from WebGL):
## http://www.web3d.org/wiki/index.php/STEP-NC_WebGL_Demo
## X3dToX3dom.xslt needs XSL2, so cannot be run with xsltproc (XSL1 only), but e.g. xalan is fine


.PHONY: all

all : 8faced-rhombic-dodecahedron_twinned_simp+def.xhtml

%.xhtml : %.x3d X3dToX3dom.xslt
	xalan -in $< -xsl $(lastword $^) -out $@


### let browser do it (does not work?):
## http://stackoverflow.com/questions/14324930/how-can-i-compile-an-xslt-file-to-bytecode-and-run-it-in-a-browser?rq=1
# %.xhtml : %.x3d
# 	cat - $< > $@ << EOF
# 		<?xml-stylesheet type="text/xsl" href="X3dToX3dom.xslt" ?>
# 	EOF
# 	sed -i '1N;s/\(.*\)\n\(.*\)/\2\n\1/' $@ # swap first two lines
# 	firefox $@
