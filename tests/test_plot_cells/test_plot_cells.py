#!/usr/bin/env python

import os
import sys
import shutil
import glob
import hashlib
from PIL import Image


import matplotlib
import matplotlib.pyplot as plt


sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import TestHarness
from input_set import SimpleLatticeInput
import openmoc
from openmoc.plotter import plot_cells


class PlotCellsTestHarness(TestHarness):
    """Test cell plotting with a 4x4 lattice."""

    def __init__(self):
        super(PlotCellsTestHarness, self).__init__()
        self.input_set = SimpleLatticeInput()
        self.figures = []

    def _run_openmoc(self):
        """Plot the cells in the geometry."""

        # Create a series of Matplotlib Figures / PIL Images for different
        # plotting parameters and append to figures list
        geometry = self.input_set.geometry
        self.figures.append(
            plot_cells(geometry, gridsize=100, get_figure=True))
        self.figures.append(
            plot_cells(geometry, gridsize=100, zcoord=10., get_figure=True))
        self.figures.append(
            plot_cells(geometry, gridsize=100, get_figure=True, 
                       xlim=(0., 2.), ylim=(0., 2.)))
        self.figures.append(
            plot_cells(geometry, gridsize=100, get_figure=True, library='pil'))

    def _get_results(self, num_iters=True, keff=True, fluxes=True,
                     num_fsrs=False, num_tracks=False, num_segments=False,
                     hash_output=True):

        outstr = ''

        # Loop over each Matplotlib figure / PIL Image and hash it
        for i, fig in enumerate(self.figures):
            plot_filename = 'cells-{0}.png'.format(i)

            # Save the figure to a file
            if isinstance(fig, matplotlib.figure.Figure):
                fig.savefig(plot_filename, bbox_inches='tight')
                plt.close(fig)
            else:
                fig.save(plot_filename)
            
            # Open the image file in PIL
            img = Image.open(plot_filename)

            # Hash the image and append to output string
            plot_hash = hashlib.md5(img.tobytes()).hexdigest()
            outstr += '{}\n'.format(plot_hash)

        return outstr

    def _cleanup(self):
        """Delete track, plot, etc. directories and test files."""

        # Find all plot files
        outputs = glob.glob(os.path.join(os.getcwd(), '*.png'))

        # Remove each plot file if it exists
        for output in outputs:
            if os.path.isfile(output):
                os.remove(output)
            elif os.path.isdir(output):
                shutil.rmtree(output)

        super(PlotCellsTestHarness, self)._cleanup()


if __name__ == '__main__':
    harness = PlotCellsTestHarness()
    harness.main()
