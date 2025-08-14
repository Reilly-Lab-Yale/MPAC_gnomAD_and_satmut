#best run in pypy3 for speed. 

import sys

def main():
	for line in sys.stdin.buffer:
		line=line.decode()
		if line.startswith('#'):
			continue#skip commented lines
		#split on tabs
		line=line.rstrip('\n').split('\t')

		#split out info fields
		info=line[-1].split(';')
		#remove labels. 'ref, alt, skew' is the order
		info=[i.split('=')[1] for i in info]
		
		
		#remove old info representation
		line=line[:-1]
		#add processed info back to the line
		line=line+info
		
		#fields of inbound line are, in order; 0: CHROM, 1: POS, 2: ID, 3: REF, 4: ALT, 5: QUAL, 6: FILTER, 7: ref, 8: alt, 9: skew
		
		
		emvar=abs(float(line[9]))>0.5 and max(float(line[7]),float(line[8]))>1
		#variant is an emvar if the absolute skew is greater than 0.5 and the max of ref or alt is greater than 1.
		color=None
		if emvar:
			color='255,0,0'
		else:
			color='0,0,255'
		emvar=int(emvar)#convert to int to fit in an 8 bit char

		#VCF files (the input) are 1-based, BED files (the output) are 0-based.
		
		spit=[]
		spit.append(line[0])#chrom
		spit.append(str(int(line[1])-1))#start
		spit.append(line[1])#end
		spit.append(line[2])#rsid
		spit.append('0')#score, not used
		spit.append('.')#strand, not used
		spit.append(str(int(line[1])-1))#thick start
		spit.append(line[1])#thick end
		spit.append(color)#color
		spit.append(line[7])#ref score
		spit.append(line[8])#alt score
		spit.append(line[9])#skew
		spit.append(line[3])#ref
		spit.append(line[4])#alt
		

		#spit.append(color)#color
		#spit.append(str(emvar))#emvar
		print('\t'.join(spit))
		#print('\t'.join([line[0],str(int(line[1])-1),line[1]]+line[2:5]+line[7:]+[color,str(emvar)]))


if __name__=="__main__":
	main()
