import os

cesna = open('paris.cesna.100.coms.txt','w')
walktrap = open('paris.walktrap.100.coms.txt','w')
#paris = open('paris.paris.100.coms.txt','w')
greedy = open('paris.greedy.100.coms.txt','w')

for f in os.listdir('SubComsBySize'):
	write_to = None
	if 'paris.cesna' in f:
		write_to = cesna
	elif 'paris.walktrap' in f:
		write_to = walktrap
	elif 'paris.String' in f:
		write_to = greedy
	for line in open('SubComsBySize/' + f):
		write_to.write(line)
cesna.close()
walktrap.close()
greedy.close()
