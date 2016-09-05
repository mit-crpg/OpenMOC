import os
from fractions import gcd

comp = True
submit = True

directory = 'pin-cell'
case = 'pin-cell-3d'

ppn = 12
wall = '72:00:00'
mem = 6 #GB

def get_unique_variants(num, max_decomp=10**8):
  x = set()
  for i in range(1, min(num, max_decomp) + 1):
    for j in range(1, min(num/i, max_decomp) + 1):
      for k in range(1, min(num/(i*j), max_decomp) + 1):
        if (i*j*k == num):
          x.add((i,j,k))
  return list(x)

def lcm(a, b):
  return a * b / gcd(a, b)

os.system('mkdir ' + directory)

dom = list()
for i in [1, 2, 4, 8]:
  dom += get_unique_variants(i, 2)
dom += [(2,2,2)]
print dom

nd = [1,1,1]
for i, d in enumerate(dom):
  for j in range(3):
    nd[j] = lcm(nd[j], d[j])
print nd

if comp:
  for d in dom:
    d_str = str(d[0]) + '-' + str(d[1]) + '-' + str(d[2])
    print d_str
    with open('base_makefile','r') as fh:
      lines = fh.readlines()
      with open(directory + '/Makefile', 'w') as fhw:
        for line in lines:
          line = line.replace('xxxDIRxxx', directory)
          line = line.replace('xxxCASExxx', case)
          line = line.replace('xxxDOMxxx', d_str)
          fhw.write(line)
    fname = '../models/' + directory + '/' + case + '.cpp'
    with open(fname, 'r') as fh:
      lines = fh.readlines()
      fname = directory + '/' + case + '-' + str(d[0]) + '-' + str(d[1]) \
          + '-' + str(d[2]) + '.cpp'
      with open(fname, 'w') as fhw:
        for line in lines:
          if (line.find('setNumDomainModules') == -1 and
              line.find('setDomainDecomposition') == -1):
            fhw.write(line)

          if (line.find('setRootUniverse') != -1):
            words = line.split()
            words = words[0].split('->')
            words = words[0].split('.')
            geometry = words[0]
            pieces = line.split(geometry)
            white_space = pieces[0]
            breaker = pieces[1].split('setRootUniverse')[0]
            write_line = white_space + geometry + breaker + \
                'setDomainDecomposition(' + \
                str(d[0]) + ', ' + str(d[1]) + ', ' + str(d[2]) + ');\n'
            fhw.write(write_line)
            write_line = white_space + geometry + breaker + \
                'setNumDomainModules(' + \
                str(nd[0]/d[0]) + ', ' + str(nd[1]/d[1]) + ', ' + \
                str(nd[2]/d[2]) + ');\n'
            fhw.write(write_line)

    os.chdir(directory)
    os.system('ls')
    os.system('make')
    os.chdir('..')


if submit:
  for d in dom:
    append = str(d[0]) + '-' + str(d[1]) + '-' + str(d[2])
    with open('strip_file.py', 'r') as fh:
      lines = fh.readlines()
      with open(directory + '/strip_file.py', 'w') as fhw:
        for line in lines:
          line = line.replace('xxxMACHINExxx', "'" + append + "'")
          fhw.write(line)

    with open('submit_base.pbs','r') as fh:
      lines = fh.readlines()
      with open(directory + '/submit.pbs', 'w') as fhw:
        name = case + '-' + append
        nodes = d[0] * d[1] * d[2]
        machine = 'machine-' + append
        executable = name
        for line in lines:
          line = line.replace('xxxNAMExxx', name)
          line = line.replace('xxxNODESxxx', str(nodes))
          line = line.replace('xxxPPNxxx', str(ppn))
          line = line.replace('xxxWTxxx', wall)
          line = line.replace('xxxMEMxxx', str(mem))
          line = line.replace('xxxMACHINExxx', machine)
          line = line.replace('xxxEXECUTExxx', executable)
          fhw.write(line)

    os.chdir(directory)
    os.system('qsub sumbmit.pbs')
    os.chdir('..')
