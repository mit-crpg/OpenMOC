import os
import time
from fractions import gcd

comp = True
submit = True

directory = 'Takeda'
case = 'Takeda-unrodded'

ppn = 12
wall = '72:00:00'
mem = 6 #GB
qlim = 4
node_lim = 8
task_lim = 96
wait_time = 3

def get_unique_variants(num, max_decomp=10**8):
    x = set()
    for i in range(1, min(num, max_decomp) + 1):
        for j in range(1, min(num//i, max_decomp) + 1):
            for k in range(1, min(num//(i*j), max_decomp) + 1):
                if (i*j*k == num):
                    x.add((i,j,k))
    return list(x)

def lcm(a, b):
    return a * b / gcd(a, b)

def checkQueue(req_nodes, req_tasks):
    curr_queue = qlim + 1
    curr_nodes = node_lim + 1
    curr_tasks = task_lim + 1
    user = os.environ['USER']
    while curr_queue > qlim or curr_nodes > node_lim or curr_tasks > task_lim:
        os.system('qstat -u {} > temp_queue'.format(user))
        q_num = 0
        task_num = 0
        node_num = 0
        with open('temp_queue','r') as fh:
            lines = fh.readlines()
            for line in lines:
                if line.find(user) != -1:
                    q_num += 1
                    words = line.split()
                    curr_tasks += int(words[-5])
                    curr_nodes += int(words[-6])
        os.system('rm temp_queue')
        curr_queue = q_num + 1
        curr_tasks = task_num + req_tasks
        curr_nodes = node_num + req_nodes
        if curr_queue > qlim or curr_nodes > node_lim or curr_tasks > task_lim:
            time.sleep(wait_time)


os.system('mkdir ' + directory)

dom = list()
for i in [1, 2, 4]:
    dom += get_unique_variants(i, 2)
print(dom)

nd = [1,1,1]
for i, d in enumerate(dom):
    for j in range(3):
        nd[j] = lcm(nd[j], d[j])
print(nd)

if comp:
    for d in dom:
        d_str = str(d[0]) + '-' + str(d[1]) + '-' + str(d[2])
        print(d_str)
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
                                str(d[0]) + ', ' + str(d[1]) + ', ' + \
                                str(d[2]) + ', MPI_COMM_WORLD);\n'
                        fhw.write(write_line)
                        write_line = white_space + geometry + breaker + \
                                'setNumDomainModules(' + \
                                str(nd[0]/d[0]) + ', ' + \
                                str(nd[1]/d[1]) + ', ' + \
                                str(nd[2]/d[2]) + ');\n'
                        fhw.write(write_line)

        os.chdir(directory)
        os.system('ls')
        os.system('make -j')
        os.chdir('..')
    os.chdir(directory)
    os.system('mv *.o obj')
    os.chdir('..')

if submit:
    for d in dom:
        nodes = d[0] * d[1] * d[2]
        append = str(d[0]) + '-' + str(d[1]) + '-' + str(d[2])
        strip_file = 'strip_file-' + append + '.py'
        with open('strip_file.py', 'r') as fh:
            lines = fh.readlines()
            with open(directory + '/' + strip_file, 'w') as fhw:
                for line in lines:
                    line = line.replace('xxxMACHINExxx', "'" + append + "'")
                    fhw.write(line)

        with open('submit_base.pbs','r') as fh:
            lines = fh.readlines()
            with open(directory + '/submit.pbs', 'w') as fhw:
                name = case + '-' + append
                machine = 'machine-' + append
                executable = name
                for line in lines:
                    line = line.replace('xxxNAMExxx', name)
                    line = line.replace('xxxNODESxxx', str(nodes))
                    line = line.replace('xxxPPNxxx', str(ppn))
                    line = line.replace('xxxWTxxx', wall)
                    line = line.replace('xxxMEMxxx', str(mem))
                    line = line.replace('xxxSTRIPxxx', strip_file)
                    line = line.replace('xxxMACHINExxx', machine)
                    line = line.replace('xxxEXECUTExxx', executable)
                    fhw.write(line)

        checkQueue(nodes, nodes*ppn)
        os.chdir(directory)
        os.system('qsub submit.pbs')
        os.chdir('..')

if False:
    os.chdir(directory)
    os.system('mkdir temp')
    os.system('mv * temp')
    os.system('mv *.o* ..')
    os.system('mv *.e* ..')
    os.system('../../')
