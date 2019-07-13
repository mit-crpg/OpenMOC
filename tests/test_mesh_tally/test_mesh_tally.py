#!/usr/bin/env python

import os
import sys
import hashlib
from collections import OrderedDict
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import SimpleLatticeInput
import openmoc.process as process
import numpy as np
import openmoc


class MeshTallyTestHarness(TestHarness):
    """An eigenvalue calculation with a mesh tally of the rates
    using the C++ Mesh class."""

    def __init__(self):
        super(MeshTallyTestHarness, self).__init__()
        self.input_set = SimpleLatticeInput(3)

        # Change spacing to avoid having rays start on lattice planes
        # Those rays are problematic because they cross through fuel pins
        # parallelly to sector planes.
        self.spacing = 0.12

    def _run_openmoc(self):
        """Run an OpenMOC eigenvalue calculation."""
        super(MeshTallyTestHarness, self)._run_openmoc()

    def _get_results(self, num_iters=True, keff=True, fluxes=True,
                     num_fsrs=True, num_tracks=True, num_segments=True,
                     hash_output=True):
        """Digest info from the mesh tallies and return as a string."""
        #NOTE Indexing in mesh tally is x,y,z and the mesh is a tuple

        # Create OpenMOC Mesh on which to tally fission rates
        mesh = openmoc.Mesh(self.solver)
        mesh.createLattice(4, 4)

        # List reaction rates wanted
        rxn_types = OrderedDict(
            {'flux'      : openmoc.FLUX_RX,
             'total'     : openmoc.TOTAL_RX,
             'fission'   : openmoc.FISSION_RX,
             'nu-fission': openmoc.NUFISSION_RX,
             'absorption': openmoc.ABSORPTION_RX})

        # Tally volume-integrated OpenMOC rates on the Mesh
        tallies = OrderedDict()
        for rxn, rxn_type in rxn_types.items():
            tallies[rxn] = mesh.getReactionRates(rxn_types[rxn])

        # Test defining a lattice from Python
        lattice = openmoc.Lattice()
        lattice.setNumX(4)
        lattice.setNumY(4)
        lattice.setNumZ(1)
        lattice.setWidth(10, 10, 1)
        lattice.computeSizes()
        mesh.setLattice(lattice)

        # Tally volume-averaged OpenMOC rates on the Mesh
        for rxn, rxn_type in rxn_types.items():
            tallies[rxn+"-ave"] = mesh.getReactionRates(rxn_types[rxn], True)

        # Append reaction rates to the output string
        outstr = ''
        for rxn, rates in tallies.items():
            outstr += rxn.title() + ' Rate Mesh Tally\n'
            rates = ['{0:12.6E}'.format(rate) for rate in rates]
            outstr += '\n'.join(rates) + '\n'

        return outstr


if __name__ == '__main__':
    harness = MeshTallyTestHarness()
    harness.main()
