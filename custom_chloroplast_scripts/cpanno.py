#!/usr/bin/python3
import subprocess
import os.path
import sys
import re
import tempfile
import argparse

version = "0.9.2"

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--ref", "--reference", dest="reference",
help="specify the reference file (mandatory)")
parser.add_argument("-i", "-s", "--input", "--sequence", dest="sequence",
help="specify the input sequence (mandatory)")
parser.add_argument("-q", "--quiet", action="store_true", dest="quiet",
help="don't print status messages to stdout")
parser.add_argument("-v", "--verbose", action="count", dest="verbosity",
default=1,       help="do print lots of status messages to stdout")
parser.add_argument("-t", "--threads", dest="threads", default="1", help="number
of threads to use (optional, for MAFFT)")
parser.add_argument("--tmp", "--tmpdir", dest="tmpdir", default=None,
help="optionally specify a temp directory")
parser.add_argument("--ali", dest="ali", default=None, help="optionally specify
a ready-made alignment FASTA file")
parser.add_argument("-o", "--out", "--outfile", default="out.gb",
dest="outfile", help="specify the output filename, defaults to './out.gb'")
parser.add_argument("-g", "--graph", dest="graphout", default=None,
help="specify PNG output file if you want your outfile graphed.")
parser.add_argument("--vers", "--version", action="store_true", default=False,
dest="vers", help="display version number and exit")
args = parser.parse_args()

if args.vers:
        print("Version " + version + ".")
        print("Bye.")
        sys.exit()
if args.quiet :
        args.verbosity = 0

def log(mess, lev=1):
        """Prints out mess if verbosity is at least lev"""
        if args.verbosity >= lev:
                print(mess)
                return True

log(args, 2)
log("cpanno starting up.")

if(not args.graphout is None):
        log("Trying to load dna_feature_viewer", 2)
        from dna_features_viewer import BiopythonTranslator

if args.reference is None or not os.path.isfile(args.reference):
        print("\nReference sequence missing or not readable.\n")
        sys.exit()
if args.sequence is None or not os.path.isfile(args.sequence):
        print("\nInput sequence missing or not readable.\n")
        sys.exit()

log('Handling alignment')
if args.tmpdir is None:
        tmpdir = tempfile.mkdtemp(prefix = 'cpanno_')
else:
        if os.path.isdir(args.tmpdir):
                tmpdir = args.tmpdir
        else:
                print("\nTemp dir not found.\n")
                sys.exit()

seqret_ref = ["seqret", "-auto", "-sequence", args.reference, "-snucleotide", "-
osformat", "fasta", "-feature", "-offormat", "gff", "-ofname", tmpdir +
"/ref.gff" ,"-outseq", tmpdir + "/ref.fa"]
seqret_inp = ["seqret", "-auto", "-sequence", args.sequence, "-snucleotide", "-
osformat", "fasta", "-outseq", tmpdir + "/inp.fa"]
log("Calling: " + ' '.join(seqret_ref), 2)
subprocess.call(seqret_ref)
log("Calling: " + ' '.join(seqret_inp), 2)
subprocess.call(seqret_inp)
log("Concatenating " + tmpdir + "/ref.fa and " + tmpdir + "/inp.fa", 2)
#subprocess.call(["cat", tmpdir + "/ref.fa", tmpdir + "/inp.fa", ">",  tmpdir +
"/mafinp.fa"], shell=True)
with open(tmpdir + '/mafinp.fa', 'w') as outfile:
    for infilename in [tmpdir + "/ref.fa", tmpdir + "/inp.fa"]:
        with open(infilename) as infile:
            outfile.write(infile.read())
if args.ali is None:
        mafft = ["mafft", "--nuc", "--memsave", "--fft", "--thread",
args.threads, "--inputorder", "--maxiterate", "2"]
        if args.verbosity < 2:
                mafft.append("--quiet")
        mafft.append(tmpdir + "/mafinp.fa")
        log("Calling: " + ' '.join(mafft), 2)
        mafout = open(tmpdir + "/ali.fa", 'w')
        subprocess.call(mafft, stdout=mafout)
        mafout.close()
        alif = open(tmpdir + "/ali.fa", 'r')
else:
        alif = open(args.ali, 'r')
log("Extracting sequences from alignment.", 2)
ali = alif.read()
alif.close()
ali = ali.split("\n")
gfff = open(tmpdir + '/ref.gff', 'r')
gff = gfff.read()
gfff.close()
refali = ali[0:int(len(ali)/2)]
inpali = ali[int(len(ali)/2):len(ali)]
refhead = refali.pop(0)
inphead = inpali.pop(0)
refali = ''.join(refali)
inpali = ''.join(inpali)
inpseq = inpali.replace('-', '')
inplen = len(inpseq)
log("Matching positions.")
i=0
j=0
k=0
rpos = []
spos = []
for br in refali:
        bs = inpali[j]
        if br != '-':
                i +=1
        rpos.append(i)
        if bs != '-':
                k +=1
        j +=1
        spos.append(k)
gff = gff.split("\n")
print(rpos)
print(spos)
i=0
for gfl in gff:
        gff[i] = gfl.split("\t")
        i +=1
gff2 = []
for i in range(len(gff)):
        gf = gff[i]
        gf2 = gf.copy()
        start2 = 0
        end2 =0
        if len(gf) > 5:
                start = gf[3]
                end = gf[4]
                print(start)
                start2 = spos[rpos.index(int(start))]
                end2 = spos[rpos.index(int(end))]
                if start != start2:
                        gf2[3] = str(start2)
                        gf2[4] = str(end2)
        if start2 < inplen and end2 < inplen: # skip if feature is after end of
sample sequence
                gff2.append(gf2)

        i +=1

gffout = []
for i in range(len(gff2)):
        g = '\t'.join(gff2[i])
        gffout.append(g)
log("Writing new features file.", 2)
gfnew = open(tmpdir + '/seq.gff', 'w')
gfnew.write('\n'.join(gffout))
gfnew.close()
seqret_out = ["seqret", "-auto", "-sequence", tmpdir + "/inp.fa", "-
snucleotide", "-scircular", "-feature", "-fopenfile", tmpdir + "/seq.gff", "-
fformat", "gff", "-osformat", "genbank", "-outseq", args.outfile]
log("Calling: " + ' '.join(seqret_out), 2)
subprocess.call(seqret_out)
log("Output written to '" + args.outfile + "'.")
if not args.graphout is None:
        graphic_record = BiopythonTranslator().translate_record(args.outfile)
        ax, _ = graphic_record.plot(figure_width=10,
strand_in_label_threshold=7)
        ax.figure.savefig(args.graphout)
