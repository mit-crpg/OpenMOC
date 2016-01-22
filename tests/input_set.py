import sys
from abc import ABCMeta, abstractmethod

sys.path.insert(0, 'openmoc')
import openmoc


class InputSet(object):
    """An abstract class for defining OpenMOC test geometries."""

    __metaclass__ = ABCMeta

    def __init__(self):
        self.materials = {}
        self.geometry = None

    @abstractmethod
    def create_materials(self):
        """Instantiate Materials."""
        return

    @abstractmethod
    def create_geometry(self):
        """Instantiate a Geometry."""
        return


class PinCellInput(InputSet):
    """A pin cell problem from sample-input/pin-cell."""

    def create_materials(self):
        """Instantiate C5G7 Materials."""
        self.materials = \
            openmoc.materialize.load_from_hdf5(filename='c5g7-mgxs.h5',
                                               directory='../../sample-input/')

    def create_geometry(self):
        """Instantiate a pin cell Geometry."""

        zcylinder = openmoc.ZCylinder(x=0.0, y=0.0, radius=1.0, name='pin')
        xmin = openmoc.XPlane(x=-2.0, name='xmin')
        xmax = openmoc.XPlane(x=+2.0, name='xmax')
        ymin = openmoc.YPlane(y=-2.0, name='ymin')
        ymax = openmoc.YPlane(y=+2.0, name='ymax')

        xmin.setBoundaryType(openmoc.REFLECTIVE)
        xmax.setBoundaryType(openmoc.REFLECTIVE)
        ymin.setBoundaryType(openmoc.REFLECTIVE)
        ymax.setBoundaryType(openmoc.REFLECTIVE)

        fuel = openmoc.Cell(name='fuel')
        fuel.setFill(self.materials['UO2'])
        fuel.addSurface(halfspace=-1, surface=zcylinder)

        moderator = openmoc.Cell(name='moderator')
        moderator.setFill(self.materials['Water'])
        moderator.addSurface(halfspace=+1, surface=zcylinder)
        moderator.addSurface(halfspace=+1, surface=xmin)
        moderator.addSurface(halfspace=-1, surface=xmax)
        moderator.addSurface(halfspace=+1, surface=ymin)
        moderator.addSurface(halfspace=-1, surface=ymax)

        root_universe = openmoc.Universe(name='root universe')
        root_universe.addCell(fuel)
        root_universe.addCell(moderator)

        self.geometry = openmoc.Geometry()
        self.geometry.setRootUniverse(root_universe)

        super(PinCellInput, self).create_geometry()
