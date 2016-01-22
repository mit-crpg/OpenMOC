from abc import ABCMeta, abstractmethod

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
        self.geometry.initializeFlatSourceRegions()
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
        left = openmoc.XPlane(x=-2.0, name='left')
        right = openmoc.XPlane(x=2.0, name='right')
        top = openmoc.YPlane(y=2.0, name='top')
        bottom = openmoc.YPlane(y=-2.0, name='bottom')

        left.setBoundaryType(openmoc.REFLECTIVE)
        right.setBoundaryType(openmoc.REFLECTIVE)
        top.setBoundaryType(openmoc.REFLECTIVE)
        bottom.setBoundaryType(openmoc.REFLECTIVE)

        fuel = openmoc.Cell(name='fuel')
        fuel.setFill(self.materials['UO2'])
        fuel.addSurface(halfspace=-1, surface=zcylinder)

        moderator = openmoc.Cell(name='moderator')
        moderator.setFill(self.materials['Water'])
        moderator.addSurface(halfspace=+1, surface=zcylinder)
        moderator.addSurface(halfspace=+1, surface=left)
        moderator.addSurface(halfspace=-1, surface=right)
        moderator.addSurface(halfspace=+1, surface=bottom)
        moderator.addSurface(halfspace=-1, surface=top)

        root_universe = openmoc.Universe(name='root universe')
        root_universe.addCell(fuel)
        root_universe.addCell(moderator)

        self.geometry = openmoc.Geometry()
        self.geometry.setRootUniverse(root_universe)

        super(PinCellInput, self).create_geometry()
