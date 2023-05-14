# Imported from: https://raw.githubusercontent.com/CERN-PH-CMG/cmg-cmssw/CMGTools-from-CMSSW_7_2_3/CMGTools/Production/python/hadd.py
import os
import pprint
import pickle
import shutil

# Check if the specified file is a nanoAOD file.
def checkNano(file):
    from ROOT import TFile
    rf = TFile.Open(file)
    keys = rf.GetListOfKeys()
    if (keys.FindObject("Events") != None and keys.FindObject("Runs") != None and  keys.FindObject("LuminosityBlocks") != None) :
        isNano = True
    else :
        isNano= False
    rf.Close()
    return isNano



def haddPck(file, odir, idirs):
    '''add pck files in directories idirs to a directory outdir.
    All dirs in idirs must have the same subdirectory structure.
    Each pickle file will be opened, and the corresponding objects added to a destination pickle in odir.
    '''
    sum = None
    for dir in idirs:
        fileName = file.replace( idirs[0], dir )
        pckfile = open(fileName)
        obj = pickle.load(pckfile)
        if sum is None:
            sum = obj
        else:
            try:
                sum += obj
            except TypeError:
                # += not implemented, nevermind
                pass
                
    oFileName = file.replace( idirs[0], odir )
    pckfile = open(oFileName, 'w')
    pickle.dump(sum, pckfile)
    txtFileName = oFileName.replace('.pck','.txt')
    txtFile = open(txtFileName, 'w')
    txtFile.write( str(sum) )
    txtFile.write( '\n' )
    txtFile.close()
    

def hadd(file, odir, idirs):
    if file.endswith('.pck'):
        try:
            haddPck( file, odir, idirs)
        except ImportError:
            pass
        return
    elif not file.endswith('.root'):
        return
    ofile =  file.replace( idirs[0], odir ) 
    destdir =os.path.dirname(ofile)
    if not os.path.exists(destdir) :
        os.mkdir(destdir)
    if len(idirs) == 1 :
        haddCmd = ['cp ',file,ofile]
    else:
        if (checkNano(file.replace( idirs[0], idirs[1]))) : # Check if the first file to be hadded is a nanoAOD file
            haddCmd = ['haddnano.py']
        else:
            haddCmd = ['hadd -ff']
        haddCmd.append(ofile)
        for dir in idirs:
            haddCmd.append( file.replace( idirs[0], dir ) )
    # import pdb; pdb.set_trace()
    cmd = ' '.join(haddCmd)
    print(cmd)
    exc = os.system(cmd)
    if exc != 0 :
        print('---> ABORTING <---')
        exit(1)


def haddRec(odir, idirs):
    print('adding', idirs)
    print('to', odir)

    try:
        os.makedirs( odir )
    except OSError:
        print() 
        print('ERROR: directory in the way. Maybe you ran hadd already in this directory? Remove it and try again, or run with -r')
        print()
        raise
    for root,dirs,files in os.walk( idirs[0] ):
        for file in files:
            hadd('/'.join([root, file]), odir, idirs)

def haddChunks(idir, removeDestDir, cleanUp=False, destdir=None ):
    if destdir == None : destdir = idir
        
    chunks = {}
    for file in sorted(os.listdir(idir)):
        filepath = '/'.join( [idir, file] )
        # print filepath
        if os.path.isdir(filepath):
            compdir = file
            try:
                prefix,num = compdir.split('_Chunk')
            except ValueError:
                # ok, not a chunk
                continue
            # print prefix, num
            chunks.setdefault( prefix, list() ).append(filepath)
    if len(chunks)==0:
        print('warning: no chunk found.')
        return
    for i, (comp, cchunks) in enumerate(chunks.iteritems(), start=1):
        odir = '/'.join( [destdir, comp] )
        print()
        print("======================")
        print("hadding folder", i, "/", len(chunks))
        print("======================")
        print(odir, cchunks)
        if removeDestDir:
            if os.path.isdir( odir ):
                shutil.rmtree(odir)
        haddRec(odir, cchunks)
    if cleanUp:
        chunkDir = 'Chunks'
        if os.path.isdir('Chunks'):
            shutil.rmtree(chunkDir)
        os.mkdir(chunkDir)
        print(chunks)
        for comp, chunks in chunks.iteritems():
            for chunk in chunks:
                shutil.move(chunk, chunkDir)
        
if __name__ == '__main__':
    import sys
    args = sys.argv
    # odir = args[1]
    # idirs = args[2:]
    # haddRec(odir, idirs)
    haddChunks(sys.argv[1])
