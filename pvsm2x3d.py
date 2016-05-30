#!pvpython

## API: http://www.paraview.org/ParaView/Doc/Nightly/www/py-doc/index.html


import paraview.simple as pvs
import os.path


def main():
    import sys       # to get command line args
    import argparse  # to parse options for us and print a nice help message

    argv = sys.argv

    if __file__ in argv:
        argv = argv[argv.index(__file__) + 1:]  # get all args after __file__


    usage_text =  "Run ParaView in background mode with this script:"
    usage_text += "  pvpython " + __file__ + " [options]"

    parser = argparse.ArgumentParser(description=usage_text)

    parser.add_argument("-i", "--input", dest="input", metavar='FILE', required=True, help="Input path contained in a text file.")
    parser.add_argument("-p", "--plugin", dest="plugin", metavar='FILE', required=True, help="Path to FA-plugin (libFacetAnalyser.so).")
    parser.add_argument("-o", "--output", dest="output", metavar='FILE', required=True, help="Output name to export WebGL")

    args = parser.parse_args(argv)

    if not argv:
        parser.print_help()
        return

    if not args.input:
       print('Need an input file')
       parser.print_help()
       sys.exit(1)

    if not args.plugin:
       print('Need a plugin file')
       parser.print_help()
       sys.exit(1)

    if not args.output:
       print('Need an output file')
       parser.print_help()
       sys.exit(1)

    pvs.LoadPlugin(args.plugin, remote='False')

    ## read pvsm
    pvs.LoadState(args.input)

    rv= pvs.GetRenderView()

    exporters= pvs.servermanager.createModule("exporters")
    # dir(exporters) # lists export modules: CSVExporter, CinemaExporter, GL2PSContextViewExporterBase, GL2PSContextViewExporterEPS, GL2PSContextViewExporterPDF, GL2PSContextViewExporterPS, GL2PSContextViewExporterSVG, GL2PSExporterBase, GL2PSRenderViewExporterBase, GL2PSRenderViewExporterEPS, GL2PSRenderViewExporterPDF, GL2PSRenderViewExporterPS, GL2PSRenderViewExporterSVG, POVExporter, VRMLExporter, WebGLExporter, X3DExporter, X3DExporterBinary
    
    exporter= exporters.WebGLExporter(FileName=args.output)
    exporter.SetView(rv)
    exporter.Write()

    ## output for ctest for regex-check because script itself fails with:
    ## Inconsistency detected by ld.so: dl-close.c: _dl_close: Assertion `map->l_init_called' failed!
    ## this happens even without qt-at-spi: http://www.cfd-online.com/Forums/openfoam-paraview/128851-pvpython-ubuntu-deb-package.html
    if os.path.isfile(os.path.splitext(args.output)[0]+".html"):
        print("Export succeded")
    else:
        print("Export failed")

if __name__ == "__main__":
    main()

