"""Test VASP-io."""

import io
import unittest

from cogue.interface.vasp_io import (
    Vasprunxml,
    VasprunxmlExpat,
    read_poscar_yaml,
    write_poscar,
)


class TestVASPIO(unittest.TestCase):
    """Test VASP-io."""

    def setUp(self):
        """Set up."""
        pass

    def tearDown(self):
        """Tear down."""
        pass

    def test_read_poscar_yaml(self):
        """Test reading POSCAR.yaml."""
        filename = "POSCAR.yaml"
        cell, poscar_order = read_poscar_yaml(filename)
        for s, p in zip(cell.get_symbols(), cell.get_points().T):
            print(s, p)
        print(write_poscar(cell))
        print(poscar_order)

    def test_VasprunxmlExpat(self):
        """Test VasprunxmlExpat."""
        with io.open("vasprun-stropt.xml", "rb") as f:
            vxml = VasprunxmlExpat(f)
        vxml.parse()

        print("VasprunxmlExpat")
        print("Forces:")
        print(vxml.get_forces())
        print("Stress:")
        print(vxml.get_stress())
        print("Atomic points:")
        print(vxml.get_points())
        print("Lattice:")
        print(vxml.get_lattice())
        print("Energy:")
        print(vxml.get_energies())

    def test_Vasprunxml(self):
        """Test Vasprunxml."""
        vxml = Vasprunxml("vasprun-energy.xml")
        vxml.parse_calculation()
        vxml.parse_eigenvalues()
        print("Forces:")
        print(vxml.get_forces())
        print("Stress:")
        print(vxml.get_stress())
        print("Lattice:")
        print(vxml.get_lattice())
        print("Atomic points:")
        print(vxml.get_points())
        print("Energy:")
        print(vxml.get_energies())
        print("Eigenvalues:")
        print(vxml.get_eigenvalues())
        print("Kpoints:")
        print(vxml.get_kpoints())
        print("Occupancies:")
        print(vxml.get_occupancies())
        print("Born charges")
        print(vxml.get_born_charges())
        print("Epsilon")
        print(vxml.get_epsilon())


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestVASPIO)
    unittest.TextTestRunner(verbosity=2).run(suite)
    # unittest.main()
