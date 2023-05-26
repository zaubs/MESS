import numpy as np
import glob
import os

print('starting')

for path, dirs, files in os.walk('/srv/meteor/klingon/evcorr'):
	for filename in files:
		fullpath = os.path.join(path,filename)

		if filename == 'event.txt':
			# file_path = os.path.split(event_file_name[0])[0]
			# file_name = os.path.split(event_file_name[0])[1]
			# event_name = file_path.split('/')[-1]

			with open(fullpath) as f:
				lines = f.readlines()

				# print(lines[2].split())
				# print(lines[3].split())
				# print(lines[4].split())
				
				try:
					if lines[3].split()[13] == "'KTJ'" or lines[3].split()[15] == "'KTJ'" or lines[4].split()[15] == "'KTJ'":
						# print(fullpath)
						print(path.split('/')[-1])
						spectralfile = path.split('/')[-1]
						spectraldate = spectralfile.split('_')[0]
						if not os.path.isdir('./SpectralEvents/%s' % spectraldate):
							os.makedirs('./SpectralEvents/%s' % spectraldate)

						os.symlink('/srv/meteor/klingon/evcorr/%s/%s/event.txt' % (spectraldate, spectralfile), './SpectralEvents/%s/%s.txt' % (spectraldate,spectralfile))
				except:
					pass
