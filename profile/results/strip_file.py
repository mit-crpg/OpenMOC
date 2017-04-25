import os

x = set()
with open(os.environ['PBS_NODEFILE'], 'r') as fh:
  lines = fh.readlines()
  for line in lines:
    x.add(line)

name = 'machine-' + xxxMACHINExxx
with open(name, 'w') as fh:
  for xi in x:
    fh.write(xi)
