#!/bin/env python3

from ZZAnalysis.AnalysisStep.hadd import haddChunks


if __name__ == '__main__':

    import os
    import sys
    from optparse import OptionParser

    parser = OptionParser()
    parser.usage = """
    %prog <dir>
    Find chunks in dir, and run recursive hadd to group all chunks.
    For example: 
    DYJets_Chunk0/, DYJets_Chunk1/ ... -> hadd -> DYJets/
    WJets_Chunk0/, WJets_Chunk1/ ... -> hadd -> WJets/
    """
    parser.add_option("-r","--remove", dest="remove",
                      default=False,action="store_true",
                      help="remove existing destination directories.")
    parser.add_option("-c","--clean", dest="clean",
                      default=False,action="store_true",
                      help="move chunks to Chunks/ after processing.")
    parser.add_option("-d","--destdir", dest="destdir",
                      default=None,
                      help="write to a destination directpory")

    (options,args) = parser.parse_args()

    if options.destdir != None:
        import click
        if not click.confirm('Output will be written to '+options.destdir+'\nDo you want to continue?', default=True) :
            sys.exit(0)

    if len(args)!=1:
        print('provide exactly one directory in argument.')
        sys.exit(1)

    dir = args[0]




    haddChunks(dir, options.remove, options.clean, options.destdir)

