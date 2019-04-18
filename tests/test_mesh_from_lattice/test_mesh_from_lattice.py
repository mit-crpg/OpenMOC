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
    """Create and subdivide an openmoc.process.Mesh from an openmoc.Lattice"""

    def __init__(self):
        super(MeshFromLatticeTestHarness, self).__init__()
        self.input_set = OblongLatticeGridInput()

        # Change spacing to avoid having rays start on lattice planes
        # Those rays are problematic because they cross through fuel pins
        # parallelly to sector planes.
        self.spacing = 0.12

    def _get_results(self, num_iters=False, keff=False, fluxes=False,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=False):
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
        fmt_str = \
            "Mesh {mesh_num}:\n" \
            "{dimension}\n" \
            "{width}\n" \
            "{lleft}\n" \
            "{uright}\n"
        outstr = ""
        for i, mesh in enumerate(meshes):
            outstr += fmt_str.format(
                mesh_num=i+1, dimension=mesh.dimension, width=mesh.width,
                lleft=mesh.lower_left, uright=mesh.upper_right)
        
        return outstr
    
    def _execute_test(self):
        """Build geometry and verify results."""
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
