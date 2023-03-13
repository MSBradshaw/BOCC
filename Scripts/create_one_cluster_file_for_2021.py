import os
with open('all_paris_coms.2021.txt','w') as outfile:
	for file in os.listdir('SubComs/2021/'):
		if 'paris' not in file:
			continue
		print(file)
		prefix = file.split('.coms.')[0]
		for line in open('SubComs/2021/'+file,'r'):
			outfile.write(prefix + ':' + line)
