import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#takes the annotation file name, annotation description,  and important data locations as specified in sugested routines
def read_file(plik,name,n,x,y,z,genes={}): #if no genes dictionary specified creates a new one. It is preferable to create one dictionary that is passed to all input routines when working on multiple files
	plik=open(plik)  
	for i in xrange(n):
		a=plik.readline()
	while a and a[0]!="#":
			b=a.split()
			c=(min(int(b[x]),int(b[y])),max(int(b[x]),int(b[y])))
			if b[z][0] == "+":
				if c in genes.keys():
					genes[c]=genes[c][:-1]+[name,True]
				else:
					genes[c]=[name,True]
			elif b[z][0] == "-":
				if c in genes.keys():
					genes[c]=genes[c][:-1]+[name,False]
				else:
					genes[c]=[name,False]
			a=plik.readline()
	return genes # returns modified dictionary of probable genes locations, with following format:
	#dictionary[(start_position,stop_position)]=["programs","that","found","this","CDS",True/False],where True represents "+" strand and False "-"

def in_ncbi(plik,genes={}): #reads .gff file as downloaded from NCBI
		if type(plik)==type([]):
				for i in plik:
						genes=read_file(i,"ncbi",7,3,4,6,genes)
		else:
				genes=read_file(plik,"ncbi",7,3,4,6,genes)
		return genes

def in_mga(plik,genes={}): #reads .out file from MetaGeneAnnotator
		if type(plik)==type([]):
				for i in plik:
						genes=read_file(i,"mga",4,1,2,3,genes)
		else:
				genes=read_file(plik,"mga",4,1,2,3,genes)
		return genes

def in_mgm(plik,genes={}): #reads output file from MetaGeneMark
		if type(plik)==type([]):
				for i in plik:
						genes=read_file(i,"mgm",12,3,2,1,genes)
		else:
				genes=read_file(plik,"mgm",12,3,2,1,genes)
		return genes

def in_gli(plik,genes={}): # reads .predict file from Glimmer3
		if type(plik)==type([]):
			for i in plik:
					genes=read_file(i,"gli",2,1,2,3,genes)
		else:
				genes=read_file(plik,"gli",2,1,2,3,genes)
		return genes

def in_pro(plik,genes={}): # reads .gff file from Prodigal
		if type(plik)==type([]):
				for i in plik:
						genes=read_file(i,"pro",4,4,3,6,genes)
		else:
				genes=read_file(plik,"pro",4,4,3,6,genes)
		return genes

def convert(seq): # creates seqnece from the complementary strand 
		dict={'A':'T','G':'C','C':'G','T':'A'}
		new=''
		for i in seq[::-1]:
				new+=dict[i]
		return new

def seq_write(out,seq): #when writing to file breaks the sequence into multiple lines
		while len(seq)>=80:
				out.write(seq[:80]+"\n")
				seq=seq[80:]
		if seq: out.write(seq+"\n")

def to_fasta(genes,plik,desc="",outfile=""): #requires properly formated annotation dictionary and genome sequence file handle
# it is strongly sugested to specify a protocol description which will appear in SeqRecord objects. When given an outfile name a fasta file containing all the sequences will be created
		if outfile: outfile=open(outfile,"w",0)
		records=[]
		log=0
		positions=sorted(genes.keys(), key=lambda x: int(x[0]))
		set_pos=[]
		seqs={}
		for poz in positions:
				seqs[poz]=''
				set_pos.append(poz)
		cnt=0
		pos_len=len(positions)
		while True:
				linia=plik.readline()
				if not linia: break
				elif linia[0]==">": pass
				linia=linia.strip()
				xlen=len(linia)
				for index in xrange(pos_len):
						positionsindex=positions[index]
						if positionsindex[1]>0:
								if positionsindex[1]<=xlen:
										seqs[set_pos[index]]+=linia[:positionsindex[1]]
										positions[index]=(-1,-1)
										reversed="+"
										if not genes[set_pos[index]][-1]:
												seqs[set_pos[index]]=convert(seqs[set_pos[index]])
												reversed="-"
										gene_id=str(genes[set_pos[index]][:1])+"_"+str(log)+"_"+reversed+"_"+str(set_pos[index])
										record = SeqRecord(Seq(seqs[set_pos[index]]), id=gene_id, description=desc, annotations={"loc":set_pos[index], "source":genes[set_pos[index]]})#, all_similar=genes[set_pos[index]])
										records.append(record)
										if outfile:
											outfile.write('>'+str(log+1)+"_"+all[set_pos[index]][0]+"_"+str(len(seqs[set_pos[index]]))+'\t'+str(set_pos[index])+'\n')
											seq_write(outfile,seqs[set_pos[index]])
										seqs[set_pos[index]]=''
										log+=1
										print log,"/",pos_len
								elif positionsindex[0]<xlen:
										seqs[set_pos[index]]+=linia[positionsindex[0]:]
								positions[index]=(positionsindex[0]-xlen,positionsindex[1]-xlen)
				cnt+=1
		plik.close()
		return records

def workflow(genome, ncbi, pro, gli, mga, mgm, desc="", output=""): #suggested workflow, takes genome file name, names of files containing annotation in indicated order (lists if multiple files from one program),
#and protocol decription and output file name as decribed in to_fasta function
		genes={}
		if ncbi: in_ncbi(ncbi,genes)
		if pro: in_pro(pro,genes)
		if gli: in_gli(gli,genes)
		if mga: in_mga(mga,genes)
		if mgm: in_mgm(mgm,genes)
		records=to_fasta(all,open(genome),desc,output)
		return records
