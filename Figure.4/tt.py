import sys

fin = open(sys.argv[1])
linenum = 0 
for line in fin:
	linenum += 1
	if linenum > 2:
		tmp = line.replace("\n", "").split("\t")
		tax = tmp[-1].split("; ")
		tax[-1] = tmp[0]
		tmp[-1] = "; ".join(tax)
		print "\t".join(tmp)
	else:
		print line

fin.close()
