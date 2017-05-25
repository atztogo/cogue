import sys
from cogue.crystal.builder import CellBuilder
from cogue.crystal.cell import Cell

class PointDefect(Cell):
    def __init__(self, cell):
        self._original_cell = cell.copy()
        self.set_cell(cell)

    def set_cell(self, cell):
        Cell.__init__(self,
                      lattice=cell.lattice,
                      magmoms=cell.get_magnetic_moments(),
                      masses=cell.get_masses(),
                      numbers=cell.numbers,
                      points=cell.get_points())
        
    def set_point_vacancy(self, index):
        n = len(self._symbols)
        if index < 0 or index > n - 1:
            sys.stderr.write(
                "At least a pair of point and symbol or "
                "a pair of point and number have to be set.\n")
        else:
            builder = CellBuilder(self._original_cell)
            builder.pop(index)
            self.set_cell(builder.get_cell())
            
            
