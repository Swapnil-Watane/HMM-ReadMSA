import math
# take msa file as input
fname = raw_input("Please enter input MSA filename:")
print 'You entered:', fname 
stck = []
stck1 = []
# dictionary for mapping array location to amino acids
dictry = {'A':0, 'Q':1, 'L':2, 'S':3, 'R':4, 'E':5, 'K':6, 'T':7, 'N':8, 'G':9, 'M':10, 'W':11, 'D':12, 'H':13, 'F':14, 'Y':15, 'C':16, 'I':17, 'P':18, 'V':19}
bagprob = {'A':8.26, 'Q':3.39, 'L':9.66, 'S':6.58, 'R':5.53, 'E':6.67, 'K':5.83, 'T':5.34, 'N':4.06, 'G':7.08, 'M':2.41, 'W':1.09, 'D':5.46, 'H':2.27, 'F':3.86, 'Y':2.92, 'C':1.37, 'I':5.94, 'P':4.71, 'V':6.87}

# parse msa
ifile = open(fname, 'r')
tmp = ''
while 1:
	line = ifile.readline()
	if not line:
		stck.append(tmp)
		tmp = ''
		break
	if ">" in line and tmp != '':
		stck.append(tmp)
		tmp = ''
	elif ">" not in line:
		b = line.rstrip()
		tmp = tmp + b
ifile.close()
nmseq = len(stck)
if nmseq % 2 ==0:
	cutoff = nmseq/2
else:
	cutoff = (nmseq/2) + 1
print 'Threshold to find conserved gaps: ', cutoff
print 'Number of sequences: ', nmseq
lnmsa = len(stck[0])
print 'Length of MSA: ', lnmsa
cgaps = 0
insertst = []
deletest = []
mainst = []
for i in range(0, lnmsa):
	a = ''
	for j in range(0, nmseq):
		a = a + stck[j][i]
	stck1.insert(i, a)
	for c in a:
		ct = a.count(c)
		if ct >= cutoff and c == '-':
			insertst.insert(cgaps, i)
			cgaps = cgaps + 1
			break
lnmod = lnmsa - cgaps
print "Length of model: ", lnmod
print "Conserved gaps: ",cgaps, ": ", insertst
insst = len(insertst)
for i in range(0, len(insertst)-1):
	if insertst[i] == (insertst[i+1]-1):
		insst = insst - 1
print "Insert states: ", insst
dictinc = {}
for i in range(insst):
	t = 'I' + str(i+1)
	dictinc[t] = lnmod+ 1 + i

#print dictinc
emmision = [[ 0 for i in range(lnmod + insst)] for j in range(20)]

mainct = 0
insertct = lnmod
flag = 0
ct = 0
for i in range(0, lnmsa):
	if i not in insertst:
		mn = i
	mainst.insert(i, mn)
	if flag == 1:
		flag = 0
		continue
	
	if i not in insertst:
		tmp = stck1[i]
		if tmp.count('-') > 0:
			deletest.insert(ct, i)
			ct = ct + 1
		ttl = len(tmp.replace('-',''))
		for j in ("".join(set(tmp))).replace('-',''):
			lval = (float(tmp.count(j))+float(20*bagprob[j]))/float(20+float(ttl))
			emmision[dictry[j]][mainct] = math.log(lval)
		mainct = mainct + 1
	if i in insertst:
		if i+1 in insertst:
			tmp = stck1[i] + stck1[i+1] 
			ttl = len(tmp.replace('-',''))
			for j in ("".join(set(tmp))).replace('-',''):
				lval = (float(tmp.count(j))+float(20*bagprob[j]))/float(20+float(ttl)) 
				emmision[dictry[j]][insertct] = math.log(lval) 
				flag =1
		else:
			tmp = stck1[i]
			ttl = len(tmp.replace('-',''))
			for j in ("".join(set(tmp))).replace('-',''):
				lval = (float(tmp.count(j))+float(20*bagprob[j]))/float(20+float(ttl))
				emmision[dictry[j]][insertct] = math.log(lval)
		insertct = insertct + 1


emmi = {}
ct = 1
for i in range(lnmod + insst):
	a = ''
	if i < lnmod:
		Key = i + 1
		Val = ''
		for amino, pos in dictry.iteritems():
			if emmision[pos][i] != 0:
				b = amino + ': '+ str(emmision[pos][i]) + ' '
				a = a + b
		Val = a
		emmi[Key] = Val
	else:
		Key1 = 'I' + str(ct)
		Val1 = ''
		for amino, pos in dictry.iteritems():
			if emmision[pos][i] != 0:
				b = amino + ': ' + str(emmision[pos][i]) + ' '
				a = a + b
		Val1 = a
		emmi[Key1] = Val1
		ct = ct + 1

with open('Value.txt', 'a') as f:
	print >> f,'############################# Emission Below #################################'
	print >> f, ''
f.close

for key, val in emmi.iteritems():
	with open('Value.txt', 'a') as f:
		print >> f, key, '-> ', val
f.close()


mainst = list(set(mainst))
mainst.insert(len(mainst), lnmsa)
print 'Delete States:', len(deletest)
#print 'Main states:', mainst
delst = len(deletest)
dictdel = {}
for i in range(delst):
	t = 'D' + str(i+1)
	dictdel[t] = len(mainst) + insst + i
#print dictdel
transition = [[ 0 for i in range(len(mainst)+insst+delst)] for j in range(len(mainst)+insst+delst)]
delct = 1
incct = 1
for i in range(len(mainst) - 1): 
	if mainst[i] == (mainst[i+1] - 1):
		if mainst[i+1] not in deletest:
			transition[i][i+1] = math.log(1 + 1)
		elif mainst[i+1] in deletest:
			g = stck1[mainst[i+1]].count('-')
			transition[i][i+1] = math.log(float(1) + float(1) - (float(g)/float(nmseq)))
			transition[i][dictdel['D'+str(delct)]] = math.log(float(1) + (float(g)/float(nmseq)))
			delct = delct + 1
	
	if mainst[i] == (mainst[i+1] - 2):
		if mainst[i+1] not in deletest:
			g = stck1[mainst[i+1]-1].count('-')
			g = nmseq - g
			transition[i][i+1] = math.log(float(1) + float(1) - (float(g)/float(nmseq)))
			transition[i][dictinc['I'+str(incct)]] = math.log(float(1) + (float(g)/float(nmseq)))
			transition[dictinc['I' + str(incct)]][i+1] = math.log(float(1) + float(1))
			incct = incct + 1
		elif mainst[i+1] in deletest:
			g = stck1[mainst[i+1]-1].count('-')
			g = nmseq - g
			transition[i][dictinc['I'+str(incct)]] = math.log(float(1) + (float(g)/float(nmseq)))
			transition[dictinc['I' + str(incct)]][i+1] = math.log(float(1) + float(1))
			incct = incct + 1
			d = stck1[mainst[i+1]].count('-')
			transition[i][dictdel['D'+str(delct)]] = math.log(float(1) + float(d)/float(nmseq))
			delct = delct + 1
			transition[i][i+1] = math.log(float(1) + float(1) - (float(g)/float(nmseq)) - (float(d)/float(nmseq)))
	if mainst[i] == (mainst[i+1] - 3):
		if mainst[i+1] not in deletest:
				a = ("".join(set(stck1[mainst[i+1]-2]))).replace('-','')
				b = ("".join(set(stck1[mainst[i+1]-1]))).replace('-','')
				c = len(a) + len(b)
				ct = 0
				seq1 = stck1[mainst[i+1]-2]
				seq2 = stck1[mainst[i+1]-1]
				# transition from insert to insert
				bothamino = 0
				# transition from insert to main
				gapamino = 0
				for k in range(nmseq):
					if seq1[k] == seq2[k] == '-':
						ct = ct + 1
					else:
						if seq1[k] != '-' and seq2[k] != '-':
							bothamino = bothamino + 1
						elif seq1[k] != '_' and seq2[k] == '-':
							gapamino = gapamino + 1
						else:
							gapamino = gapamino + 1
				ct = (bothamino*2) + gapamino 
				transition[i][i+1] = math.log(float(1) + float(1) - (float(ct)/float(nmseq)))
				transition[i][dictinc['I'+ str(incct)]] = math.log(float(1) + float(ct)/float(nmseq))
				transition[dictinc['I' + str(incct)]][i+1] = math.log(float(1) + float(gapamino+bothamino)/float(ct))
				transition[dictinc['I' + str(incct)]][dictinc['I' + str(incct)]] = -2 #float(1) + float(bothamino)/float(ct)
				incct = incct + 1
delct = 1
for i in range(len(deletest)):
	st = deletest[i]
	if i < len(deletest)-1:
		if st == deletest[i+1]-1:
			transition[dictdel['D'+str(delct)]][dictdel['D'+str(delct+1)]] = -2
			delct = delct + 1
		else:
			a = next(x for x in enumerate(mainst) if x[1] > st)[0]
			transition[dictdel['D'+str(delct)]][a] = math.log(1+1)
			delct = delct + 1
	else:
		a = next(x for x in enumerate(mainst) if x[1] > st)[0]
		#import pdb; pdb.set_trace()
		transition[dictdel['D'+str(delct)]][a] = math.log(1+1)
		delct = delct + 1


dictmn = {}
for i in range(len(mainst)):
	dictmn['M' + str(i+1)] = i

#import pdb; pdb.set_trace()
dictmn.update(dictinc)
dictmn.update(dictdel)
#print dictmn

trans = {}

m = 1
ii = 1
d = 1
for i in range(len(mainst)+ insst + delst):
	a= ''
	if i < len(mainst):
		key = i 
		for st, pos in dictmn.iteritems():
			if transition[i][pos] != 0:
				b = 'M' + str(m) + '->'+ st + ': ' + str(transition[i][pos])+ ' '
				a = a + b
		m = m + 1
		if i == len(mainst)-1:
			a = 'This is end state'
		trans[i] = a
	elif i >= len(mainst) and i < len(mainst)+insst:
		key = i
		for st, pos in dictmn.iteritems():
			if transition[i][pos] != 0:
				b = 'I' + str(ii) + '->'+ st + ': ' + str(transition[i][pos])+ ' '
				a = a + b
		ii = ii + 1
		trans[i] = a
	else:
		key = i 
		for st, pos in dictmn.iteritems():
			if transition[i][pos] != 0:
				b = 'D' + str(d) + '->'+ st + ': ' + str(transition[i][pos])+ ' '
				a = a + b
		d = d + 1
		trans[i] = a
	
with open('Value.txt', 'a') as f:
	print >> f, ''
	print >> f,'############################# Transition Below #################################'
	print >> f, ''
f.close

for key, val in trans.iteritems():
	with open('Value.txt', 'a') as f:
		print >> f, key+1, ':', val
f.close()

