import subprocess

# This test should only be collected if OpenMOC was installed with an MPI
# wrapped compiler
output = subprocess.call(["mpirun", "-n", "12", "--oversubscribe", "python", "3D_lattice.py"])

if output != 0:
    raise RuntimeError("Test failed")
