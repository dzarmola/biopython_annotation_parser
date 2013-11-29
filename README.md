biopython_annotation_parser
===========================

Creates SeqRecord objects from automatic annotation files

Obviously both Seq and SeqRecord modules from BioPython are required.

Use examples:
1.
  genes={}
  genes=in_ncbi("ncbi/NC_000913.gff",genes)
  records=to_fasta(genes,open("../ecoli/sequence.fasta"),"ecoli as annotated in ncbi")
2. 
  records=workflow("../ecoli/sequence.fasta",["ncbi/NC_000913.gff"],[],[],[],[],"ecoli as annotated in ncbi")

Included functions:

read_file(plik,name,n,x,y,z,genes={}):

takes the annotation file name, annotation description,  and important data locations as specified in sugested routines
if no genes dictionary specified creates a new one. It is preferable to create one dictionary that is passed to all input routines when working on multiple files
returns modified dictionary of probable genes locations, with following format:
dictionary[(start_position,stop_position)]=["programs","that","found","this","CDS",True/False],where True represents "+" strand and False "-"

in_ncbi(plik,genes={}): 

reads .gff file as downloaded from NCBI


in_mga(plik,genes={}): 

reads .out file from MetaGeneAnnotator


in_mgm(plik,genes={}): 

reads output file from MetaGeneMark


in_gli(plik,genes={}): 

reads .predict file from Glimmer3


in_pro(plik,genes={}): 

reads .gff file from Prodigal


convert(seq): 

creates complementary strand sequence 


seq_write(out,seq):

when writing to file breaks the sequence into multiple lines, takes output file handle and the sequence


to_fasta(genes,plik,desc="",outfile=""): 

requires properly formated annotation dictionary and genome sequence file handle
it is strongly sugested to specify a protocol description which will appear in SeqRecord objects. When given an outfile name a fasta file containing all the sequences will be created


workflow(genome, ncbi, pro, gli, mga, mgm, desc="", output=""): 

suggested workflow, takes genome file name, names of files containing annotation in indicated order (lists if multiple files from one program) and protocol decription and output file name as decribed in to_fasta function
