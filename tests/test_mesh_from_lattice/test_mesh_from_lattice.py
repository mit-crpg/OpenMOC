#!/usr/bin/env python

import os
import sys
sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import OblongLatticeGridInput
import openmoc.process as process
import openmoc


def get_cell_by_name(cells, name):
    for c in cells:
        if c.getName() == name:
            return c


class MeshFromLatticeTestHarness(TestHarness):
    """An eigenvalue calculation with a mesh tally of the rates
    using the C++ Mesh class."""

    def __init__(self):
        super(MeshFromLatticeTestHarness, self).__init__()
        self.input_set = OblongLatticeGridInput()

        # Change spacing to avoid having rays start on lattice planes
        # Those rays are problematic because they cross through fuel pins
        # parallelly to sector planes.
        self.spacing = 0.12

    def _run_openmoc(self):
        """Run an OpenMOC eigenvalue calculation."""
        super(MeshFromLatticeTestHarness, self)._run_openmoc()

    def _get_results(self, num_iters=True, keff=True, fluxes=True,
                     num_fsrs=True, num_tracks=True, num_segments=True,
                     hash_output=True):
        """Create Meshes from a lattice and return their shapes as a string"""
        
        # Find the test oblong lattice
        all_cells = harness.input_set.geometry.getAllCells().values()
        root_cell = get_cell_by_name(all_cells, "root cell")
        ulat = root_cell.getFillUniverse()
        lat = openmoc.castUniverseToLattice(ulat)
        
        
        meshes = [None]*2
        meshes[0] = process.Mesh.from_lattice(lat)
        meshes[1] = process.Mesh.from_lattice(lat, division=2)

        # Append mesh properties to the output string
        fmt_str = """\
Mesh {mesh_num}:
{dimension}
{width}
{lleft}
{uright}"""
        outstr = ''
        for i, mesh in enumerate(meshes):
            outstr += fmt_str.format(
                mesh_num=i+1, dimension=mesh.dimension, width=mesh.width,
                lleft=mesh.lower_left, uright=mesh.upper_right)
            outstr += "\n"
        
        return outstr
    
    def _execute_test(self):
        """Build geometry, ray trace, run calculation, and verify results."""
        try:
            self._create_geometry()
            results = self._get_results()
            self._write_results(results)
            self._compare_results()
        finally:
            self._cleanup()
    

if __name__ == '__main__':
    harness = MeshFromLatticeTestHarness()
    harness.main()
