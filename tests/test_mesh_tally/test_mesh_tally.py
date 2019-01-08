#!/usr/bin/env python

import os
import sys
import hashlib
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import SimpleLatticeInput
import openmoc.process as process
import numpy as np


class MeshTallyTestHarness(TestHarness):
    """An eigenvalue calculation with a mesh tally of the fission rates
    using the openmoc.process module."""

    def __init__(self):
        super(MeshTallyTestHarness, self).__init__()
        self.input_set = SimpleLatticeInput()

    def _run_openmoc(self):
        """Run an OpenMOC eigenvalue calculation."""
        super(MeshTallyTestHarness, self)._run_openmoc()

    def _get_results(self, num_iters=True, keff=True, fluxes=True,
                     num_fsrs=True, num_tracks=True, num_segments=True,
                     hash_output=True):
        """Digest info from the mesh tallies and return as a string."""

        # Create OpenMOC Mesh on which to tally fission rates
        mesh = process.Mesh()
        mesh.dimension = [4, 4]
        mesh.lower_left = [-2., -2.]
        mesh.upper_right = [2., 2.]
        mesh.width = [1., 1.]

        # Tally OpenMOC fission rates on the Mesh
        fiss_rates = mesh.tally_fission_rates(self.solver)

        # Append fission rates to the output string
        outstr = 'Fission Rate Mesh Tally\n'
        rates = ['{0:12.6E}'.format(rate) for rate in fiss_rates.ravel()]
        outstr += '\n'.join(rates) + '\n'

        # Retrieve the Materials and number of groups from the geometry
        materials = self.input_set.geometry.getAllMaterials()
        num_groups = self.input_set.geometry.getNumEnergyGroups()

        # Aggregate the total cross sections for each Material
        # into a dictionary to pass to the mesh tally
        domains_to_coeffs = {}
        for material_id in materials:
            domains_to_coeffs[material_id] = np.zeros(num_groups)
            for group in range(num_groups):
                domains_to_coeffs[material_id][group] = \
                    materials[material_id].getSigmaTByGroup(group+1)

        # Tally volume-averaged OpenMOC total rates on the Mesh
        tot_rates = mesh.tally_on_mesh(self.solver, domains_to_coeffs, 
                                       domain_type='material', 
                                       volume='integrated')

        # Append total rates to the output string
        outstr += 'Total Rate Mesh Tally\n'
        rates = ['{0:12.6E}'.format(rate) for rate in tot_rates.ravel()]
        outstr += '\n'.join(rates) + '\n'

        return outstr


if __name__ == '__main__':
    harness = MeshTallyTestHarness()
    harness.main()
