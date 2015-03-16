import os
import sys
import re
import os.path

class IncorrectSequenceLetter(ValueError):

        def __init__(self, letter, classname):
                self.message = "The sequence item %s is not found in the alphabet of class %s\n" %(letter, classname)

class AlignmentLengthError(ValueError):
        
        def __init__(self,classname):
                self.message = "The sequences of {} are not of equal length (=unaligned)".format(classname)


def FASTA_iterator(fasta_filename,SequenceType):

        fd = open(fasta_filename,"r")
        sequence = ""
        for line in fd:
                if line[0]==">":
                        if len(sequence)>0:
                                try:
                                        yield SequenceType(identifier, sequence)
                                except IncorrectSequenceLetter as e:
                                        sys.stderr.write(e.message)
                        identifier = line[1:].strip()
                        sequence = ""
                else:
                        sequence+=line.strip()
        fd.close()

        if len(sequence)>0:
                try:
                        yield SequenceType(identifier, sequence)
                except IncorrectSequenceLetter as e:
                        sys.stderr.write(e.message)

#Sequence Objects

class Sequence(object):

        alphabet = set()
        mw = {}

        def __init__(self, identifier, sequence):

                self.identifier = identifier

                for letter in sequence:
                        if letter not in self.alphabet:
                                raise IncorrectSequenceLetter( letter,  self.__class__.__name__ )

                self.sequence = sequence

                self.__mw = None

        def __eq__(self,other):
                return self.sequence==other.sequence

        def __ne__(self, other):
                return self.sequence!=other.sequence

        def __len__(self):
                return len(self.sequence)

        def get_fasta_string(self):
                splitted_seq = []
                for x in xrange(0,len(self),80):
                        splitted_seq.append(self.sequence[x:x+80])
                return ">%s\n%s\n" %(self.identifier,"\n".join(splitted_seq))
        
        def __cmp__(self, other):
                if len(self)==len(other):       return 0
                elif len(self)<len(other):      return -1
                elif len(self)>len(other):      return 1

        def __getitem__(self, k):
                return self.sequence[k]

        def __add__(self, other):
                return self.__class__( identifier = "%s_%s" %(self.identifier, other.identifier),
                                       sequence = self.sequence + other.sequence )


        def get_identifier(self):
                return self.identifier

        def get_sequence(self):
                return self.sequence

        def get_mw(self):
                if self.__mw is None:
                        self.__mw = sum( self.mw[letter] for letter in self.sequence )
                return self.__mw

        def has_subsequence( sequence_obj ):
                return sequence_obj.get_sequence() in self.sequence

        
class ProteinSequence(Sequence):
        
        #First, override all specific Class attributes for ProteinSequence
        #alphabet = set('ACDEFGHIKLMNPQRSTVWY')
        alphabet = 'ACDEFGHIKLMNPQRSTVWY'
        mw = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, \
              'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, \
              'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, \
              'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}
        reverse_RNA = {'A': 'GCU', 'C': 'UGU', None: 'UAA', 'E': 'GAG', \
                       'D': 'GAU', 'G': 'GGU', 'F': 'UUU', 'I': 'AUU', \
                       'H': 'CAU', 'K': 'AAG', 'M': 'AUG', 'L': 'UUG', \
                       'N': 'AAU', 'Q': 'CAG', 'P': 'CCU', 'S': 'UCU', \
                       'R': 'CGU', 'T': 'ACU', 'W': 'UGG', 'V': 'GUU', \
                       'Y': 'UAU'}
        reverse_DNA = {'A': 'GCT', 'C': 'TGT', None: 'TAA', 'E': 'GAG', \
                       'D': 'GAT', 'G': 'GGT', 'F': 'TTT', 'I': 'ATT', \
                       'H': 'CAT', 'K': 'AAG', 'M': 'ATG', 'L': 'TTG', \
                       'N': 'AAT', 'Q': 'CAG', 'P': 'CCT', 'S': 'TCT', \
                       'R': 'CGT', 'T': 'ACT', 'W': 'TGG', 'V': 'GTT', \
                       'Y': 'TAT'}

        def reverse_translate_to_RNA(self):
                return RNASequence(     identifier = self.identifier+"_reverse_translated",
                                        sequence = "".join( [ self.reverse_RNA[aa] for aa in self.sequence ] ) )

        def reverse_translate_to_DNA(self):
                return DNASequence(     identifier = self.identifier+"_reverse_translated",
                                        sequence = "".join( [ self.reverse_DNA[aa] for aa in self.sequence ] ) )


class NucleotideSequence(Sequence):
        complement = {}
        mw = {}
        start_codons = set()
        end_codons = set()
        
        def _get_complement(self):
                return "".join([ self.complement[letter] for letter in self.sequence ])

        def get_complement(self):
                return self.__class__( self.identifier+"_complement", \
                                      sequence="".join([ self.complement[letter] for letter in self.sequence ]) )

        def translate(self):
                started = False
                translated_sequence = ""
                for i in range(0,len(self.sequence),3):
                        codon = self.sequence[i:i+3]
                        if started is False:
                                if codon in self.start_codons:
                                        started = True
                        elif codon in self.stop_codons:
                                break
                        else:
                                translated_sequence+=self.translation_dict[codon]
                return ProteinSequence( identifier = self.identifier+"_translated", 
                                        sequence = translated_sequence )


class DNASequence(NucleotideSequence):
        alphabet = set('ACTGNactgn')
        mw = {'A': 347.0, 'C': 323.0, 'T': 322.0, 'G': 363.0} 
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'\
                      , 'a': 't', 'c' : 'g', 'g': 'c','t': 'a','n':'n'}
        
        stop_codons = set(['TAA', 'TAG', 'TGA'])
        start_codons = set(['TTG', 'CTG', 'ATG'])

        translation_dict = {'CTT': 'L', 'ATG': 'M', 'ACA': 'T', 'ACG': 'T', \
                            'ATC': 'I', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', \
                            'AGC': 'S', 'AGA': 'R', 'ATT': 'I', 'CTG': 'L', \
                            'CTA': 'L', 'ACT': 'T', 'CCG': 'P', 'AGT': 'S', \
                            'CCA': 'P', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', \
                            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'GGG': 'G', \
                            'GGA': 'G', 'GGC': 'G', 'TAC': 'Y', 'CGT': 'R', \
                            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GAG': 'E', \
                            'GTT': 'V', 'GAC': 'D', 'GAA': 'E', 'AAG': 'K', \
                            'AAA': 'K', 'AAC': 'N', 'CTC': 'L', 'CAT': 'H', \
                            'AAT': 'N', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', \
                            'TGT': 'C', 'TCT': 'S', 'GAT': 'D', 'TTT': 'F', \
                            'TGC': 'C', 'TGG': 'W', 'TTC': 'F', 'TCG': 'S', \
                            'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', \
                            'TCA': 'S', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', \
                            'GCT': 'A'}
        
        def get_complement(self):
                return DNASequence(identifier = self.identifier+"_complement",\
                                   sequence = self._get_complement() )[::-1]

        def transcribe(self):
                return RNASequence( identifier = self.identifier+"_transcribed",\
                                   sequence = self.sequence.replace("T","U") )


class DNASequenceAligned(DNASequence):
        """
        This Class is identical to DNASequence, except it alows for gaps ("-")
        in the alphabet. These will be translated to "-". Furthermore it has
        start and stop coordinates to alow for local alignments
        """
        
        alphabet = set('ACTGNactgn-')
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', \
                      'a': 't', 'c' : 'g', 'g': 'c','t': 'a','n':'n','-':'-'}        
        
        def __init__(self, identifier, sequence, start = None, stop = None, strand = None):
                self.sequence = sequence
                self.identifier = identifier
                #start and stop coordinates for potential local alignments
                
                if start is not None:
                        self.start = int(start)
                else:
                        self.start = None
                if stop is not None:
                        
                        self.stop = int(stop)
                else:
                        self.stop = None
                
                self.strand = strand

        def get_coordinates(self):
                
                return (self.start, self.stop)





class RNASequence(NucleotideSequence):
        alphabet = set('ACUG')
        mw = {'A': 363.0, 'C': 339.0, 'U': 340.0, 'G': 379.0}   
        complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}

        stop_codons = set(['UAA', 'UAG', 'UGA'])
        start_codons = set(['UUG', 'CUG', 'AUG'])

        translation_dict = {'GUC': 'V', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', \
                            'GUU': 'V', 'AAC': 'N', 'AGG': 'R', 'UGG': 'W', \
                            'AGC': 'S', 'AUC': 'I', 'AGA': 'R', 'AAU': 'N', \
                            'ACU': 'T', 'CAC': 'H', 'GUG': 'V', 'CCG': 'P', \
                            'CCA': 'P', 'AGU': 'S', 'CCC': 'P', 'GGU': 'G', \
                            'UCU': 'S', 'GCG': 'A', 'CGA': 'R', 'CAG': 'Q', \
                            'CGC': 'R', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', \
                            'CCU': 'P', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', \
                            'GAG': 'E', 'UCC': 'S', 'UAC': 'Y', 'CGU': 'R', \
                            'GAA': 'E', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', \
                            'UCA': 'S', 'AUG': 'M', 'CUG': 'L', 'AUU': 'I', \
                            'CAU': 'H', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', \
                            'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'GAC': 'D', \
                            'GUA': 'V', 'UGC': 'C', 'GCU': 'A', 'UGU': 'C', \
                            'CUC': 'L', 'UUG': 'L', 'UUA': 'L', 'GAU': 'D', \
                            'UUC': 'F'}

        def get_complement(self):
                return RNASequence(identifier = self.identifier+"_complement", \
                                   sequence = self._get_complement() )

        def reverse_transcribe(self):
                return DNASequence( identifier = self.identifier+"_transcribed",\
                                   sequence = self.sequence.replace("U","T") )


#Alignemnt Object:
class PairwiseAlignment(object):
        
        def __init__(self, SequenceA, SequenceB):
                """
                SequenceA and B are instances of DNASequenceAligned
                """
                self.SequenceA = SequenceA
                self.SequenceB = SequenceB
                
                if len(self.SequenceA.get_sequence()) != len(self.SequenceB.get_sequence()):
                        raise AlignmentLengthError(self.__class__.__name__)
                
        def get_sequenceA(self):
                return self.SequenceA
        
        def get_sequenceB(self):
                return self.SequenceB
        
def AXT_iterator(AXTfile,sampleA = "sampleA", sampleB = "sampleB"):
        """
        Iterates through axt file, returns alignment object for 2 sequences
        sampleA/B arguments are prefixes for the identiefiers of the sequence
        objects
         identifier, sequence, start = None, stop = None, strand = None)
        """
        
        with open(AXTfile) as af:
                for line in af:
                        if re.match(r'^\d', line):
                                alnN, chrTarget, startTarget, stopTarget, \
                                chrQuery, startQuery, stopQuery, \
                                strangQuery, blastZscore = line.rstrip().split()
                                
                                seqA = af.next().rstrip()
                                seqB = af.next().rstrip()
                                
                                SequenceA = DNASequenceAligned(sampleA+"_"+chrTarget+"_"+alnN, \
                                                               seqA, startTarget,stopTarget,"+")
                                SequenceB = DNASequenceAligned(sampleB+"_"+chrQuery+"_"+alnN, \
                                                               seqB, startQuery,stopQuery, strangQuery)
                                yield PairwiseAlignment(SequenceA,SequenceB)

 
 
class BedContainer(object):
    """
    A container for BedFile entries
    """
        
    def __init__(self, bedfile):
               
        self.chromosomes = {}
    
        self.bedfile = bedfile
        with open(self.bedfile) as bf:
            for line in bf:
                placeHolders = [None] * 12
                line = line.rstrip().split("\t")
                placeHolders[0:len(line)] = line
                chrom, start, stop, name, score, strand,\
                thickStart, thickEnd, itemRgb, blockCount,\
                blockSizes, blockStart = placeHolders
                
                start = int(start)
                stop = int(stop)


                if chrom in self.chromosomes:
                    self.chromosomes[chrom].append((\
                                        start, stop, name, score, strand,\
                                        thickStart, thickEnd, itemRgb, \
                                        blockCount, blockSizes, blockStart ))
                else:
                    self.chromosomes[chrom] = [(\
                                        start, stop, name, score, strand,\
                                        thickStart, thickEnd, itemRgb, \
                                        blockCount, blockSizes, blockStart )]

    def get_all_chromosomes(self):
        
        return keys(self.chromosomes)
                    
    def get_all_entries(self):
        """
        returns all entries from a given bed file as a set with
        chr=>entryTuples
        """
        return self.chromosomes

    def get_entry_for_chrom(self,chrom):
        """
        returns a list of tuples with all entries for a certain chromos
        """
            
        if chrom in self.chromosomes:
            return sorted(self.chromosomes[chrom])
        else:
            return None