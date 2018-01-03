#!/usr/bin/env python

"""GluExons - tool that extracts entire reference transcriptome from genome and GTF-formatted (GENCODE) annotation"""

__author__ = "Yury Barbitoff"
__copyright__ = "Copyright 2015, Yury Barbitoff"
__license__ = "CC"
__version__ = "1.0"
__email__ = "barbitoff@bk.ru"
__status__ = "Release"


import sys
import getopt

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'N': 'N', 'T': 'A'}

def revCompl(string):
    """Make a DNA pattern's reverse complement"""
    complementary = ''
    for letter in string:
        complementary += complement[letter]
    return complementary[::-1]


class transcript:
    """Transcript representation - holds
    ID (ENSTXXXXXXXXXXX.X),
    place (chromosome/patch name),
    locus (start/end offsets),
    strand (True/False),
    list of exons (each represented as a tuple (number in transcript, start offset, end offset)"""
    def __init__(self, no, place, startpos, endpos, strand, exons=[]):
        self.id = no
        self.place = place
        self.locus = (startpos, endpos)
        self.strand = strand
        self.exons = exons

    def build(self, piece):
        """Builds transcript sequence by exon concatenation. If no exons found, uses the transcript's locus"""
        seq = ''
        exons = sorted(self.exons)
        if self.strand:
            if len(exons) == 0:
                seq = piece[self.locus[0] - 1:self.locus[1]]
            else:
                for _, startpos, endpos in exons:
                    seq += piece[startpos - 1:endpos]
        else:
            if len(exons) == 0:
                seq = revCompl(piece[self.locus[0] - 1:self.locus[1]])
            else:
                for _, startpos, endpos in exons:
                    seq += revCompl(piece[startpos - 1:endpos])
        return seq


def optRecognize(argv):
    """Processing given command line options"""
    try:
        opts, args = getopt.getopt(argv, "ha:g:o:", ["help", "annot=", "genome=", "ofile="])
    except:
        print('\nInvalid options. To see usage info, type gluexons.py -h\n')
        sys.exit()
    annotation, genome, outfile = None, None, None
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('''\n\nTo run transcriptome construction, use gluexons.py with options:
-a (--annot) <file> - GTF/GFF file containing basic annotation of the genome
-g (--genome) <file> - FASTA formatted file with the reference genome sequence
-o (--ofile) <file> - filename (FASTA) to put output transcriptome in\n
=================================
-h (--help) To show the help message again\n''')
            sys.exit()
        elif opt in ("-a", "--annot"):
            annotation = arg
        elif opt in ("-g", "--genome"):
            genome = arg
        elif opt in ("-o", "--ofile"):
            outfile = arg
    return annotation, genome, outfile


def constructTranscriptome(annotation):
    """Extracts all transcripts included into annotation and puts them to a dictionary"""
    transcriptome = {}
    for line in annotation:
        content = line.split()
        if len(content) < 15:
            continue
        if content[2] == 'transcript':
            this_transcript = transcript(content[11][1:len(content[11])-2], content[0], content[3], content[4], content[6] == '+')
            transcriptome[this_transcript.id] = this_transcript
    return transcriptome


def mapExons(annotation, transcriptome):
    """Maps exons to transcripts"""
    for line in annotation:
        content = line.split()
        if len(content) < 15:
            continue
        if content[2] == 'exon':
            tr_id = content[11][1:len(content[11])-2]
            transcriptome[tr_id].exons = transcriptome[tr_id].exons + [(int(content[25][0]), int(content[3]), int(content[4]))]


def constructPiece(transcriptome, piece, place):
    """Iterate over transcripts, build sequences originating from chromosome/patch"""
    for transcript in transcriptome.values():
        if transcript.place != place:
            continue
        out.write('>' + transcript.id + '\n' + transcript.build(piece) + '\n')


def extractSeq(transcriptome, genome):
    """Builds and writes transcript sequences"""
    piece = ''
    for line in genome:
        if line.startswith('>'):
            if piece != '':
                constructPiece(transcriptome, piece, place)
            place = line.split()[0][1:]
            piece = ''
            continue
        piece += line.rstrip()
    constructPiece(transcriptome, piece, place)


# Get input options
if __name__ == "__main__":
    annotationfile, genomefile, outfile = optRecognize(sys.argv[1:])

# Try/catch files' opening
try:
    print '\n[progress] Opening input files...'
    annotation = open(annotationfile, 'r')
    genome = open(genomefile, 'r')
    out = open(outfile, 'w')
except:
    print '\nInvalid file names! Check specified paths! For usage info, type gluexons.py -h\n'
    sys.exit()

# Build process
print '[progress] Constructing transcriptome...'
transcriptome = constructTranscriptome(annotation)
annotation.close()

print '[progress] Indexing all exons to a transcriptome...'
with open(annotationfile, 'r') as annotation:
    mapExons(annotation, transcriptome)

print '[progress] Extracting and writing sequences...'
extractSeq(transcriptome, genome)

print '[progress] Done constructing reference transcriptome, removing trash and exiting...\n'

genome.close()
out.close()

sys.exit()
