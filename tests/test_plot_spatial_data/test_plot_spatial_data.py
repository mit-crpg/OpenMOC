#!/usr/bin/env python

import os
import sys

import pandas as pd
import numpy as np

sys.path.insert(0, os.pardir)
sys.path.insert(0, os.path.join(os.pardir, 'openmoc'))
from testing_harness import PlottingTestHarness
from input_set import SimpleLatticeInput
from openmoc.plotter import plot_spatial_data, PlotParams


class PlotSpatialDataTestHarness(PlottingTestHarness):
    """Test general spatial data plotting with a 4x4 lattice."""

    def __init__(self):
        super(PlotSpatialDataTestHarness, self).__init__()
        self.input_set = SimpleLatticeInput()

    def _run_openmoc(self):
        """Plot spatial data from Pandas DataFrames."""

        # Run an eigenvalue calculation to setup FSR centroids
        super(PlotSpatialDataTestHarness, self)._run_openmoc()

        # Seed NumPy's random number generator to ensure reproducible results
        np.random.seed(1)

        # Initialize a Pandas DataFrame with normally distributed random data
        num_fsrs = self.input_set.geometry.getNumFSRs()
        df = pd.DataFrame(np.random.randn(num_fsrs,3), columns=list('ABC'))

        # Initialize a PlotParams object
        plot_params = PlotParams()
        plot_params.geometry = self.input_set.geometry
        plot_params.suptitle = 'Pandas DataFrame'
        plot_params.filename = 'pandas-df'
        plot_params.colorbar = True

        # Enforce consistent color scheme across figures
        plot_params.vmin = df.values.min()
        plot_params.vmax = df.values.max()
        
        # Create a series Matplotlib Figures / PIL Images for different
        # plotting parameters and append to figures list
        self.figures.extend(plot_spatial_data(df, plot_params, get_figure=True))

        plot_params.library = 'pil'
        self.figures.extend(plot_spatial_data(df, plot_params, get_figure=True))

if __name__ == '__main__':
    harness = PlotSpatialDataTestHarness()
    harness.main()
