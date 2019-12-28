import unittest
import eqlib as eq
import numpy as np
from numpy.testing import assert_almost_equal, assert_equal


class TestRowMajorSparseStructure(unittest.TestCase):
    def test_rectangular(self):
        ia = [0, 2, 4, 6]
        ja = [0, 2, 2, 3, 1, 3]

        structure = eq.RowMajorSparseStructure(3, 4, ia, ja)

        assert_equal(structure.ia, ia)
        assert_equal(structure.ja, ja)

        for idx, (i, j) in enumerate([(0, 0), (0, 2), (1, 2), (1, 3), (2, 1), (2, 3)]):
            assert_equal(structure.get_index(i, j), idx)

        assert_equal(structure.get_index(1, 0), -1)

        assert_equal(structure.density, 6 / 12)

    def test_square(self):
        ia = [0, 3, 5, 6, 7]
        ja = [0, 1, 2, 1, 2, 2, 3]

        structure = eq.RowMajorSparseStructure(4, 4, ia, ja)

        assert_equal(structure.ia, ia)
        assert_equal(structure.ja, ja)

        for idx, (i, j) in enumerate([(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2), (3, 3)]):
            assert_equal(structure.get_index(i, j), idx)

        assert_equal(structure.get_index(1, 0), -1)

        assert_equal(structure.density, 7 / 16)

    def test_convert_from(self):
        ia = [0, 2, 4, 6]
        ja = [0, 2, 2, 3, 1, 3]
        a = np.array([0, 1, 2, 3, 4, 5], float)

        structure = eq.RowMajorSparseStructure(3, 4, ia, ja)

        structure = eq.ColMajorSparseStructure.convert_from(structure, a)

        assert_equal(structure.rows, 3)
        assert_equal(structure.cols, 4)
        assert_equal(structure.ia, [0, 1, 2, 4, 6])
        assert_equal(structure.ja, [0, 2, 0, 1, 1, 2])
        assert_equal(a, [0, 4, 1, 2, 3, 5])


class TestColMajorSparseStructure(unittest.TestCase):
    def test_rectangular(self):
        ia = [0, 2, 4, 6]
        ja = [0, 2, 2, 3, 1, 3]

        structure = eq.ColMajorSparseStructure(4, 3, ia, ja)

        assert_equal(structure.ia, ia)
        assert_equal(structure.ja, ja)

        for idx, (i, j) in enumerate([(0, 0), (2, 0), (2, 1), (3, 1), (1, 2), (3, 2)]):
            assert_equal(structure.get_index(i, j), idx)

        assert_equal(structure.get_index(0, 1), -1)

        assert_equal(structure.density, 6 / 12)

    def test_square(self):
        ia = [0, 3, 5, 6, 7]
        ja = [0, 1, 2, 1, 2, 2, 3]

        structure = eq.ColMajorSparseStructure(4, 4, ia, ja)

        assert_equal(structure.ia, ia)
        assert_equal(structure.ja, ja)

        for idx, (i, j) in enumerate([(0, 0), (1, 0), (2, 0), (1, 1), (2, 1), (2, 2), (3, 3)]):
            assert_equal(structure.get_index(i, j), idx)

        assert_equal(structure.get_index(0, 1), -1)

        assert_equal(structure.density, 7 / 16)

    def test_convert_from(self):
        ia = [0, 2, 4, 6]
        ja = [0, 2, 2, 3, 1, 3]
        a = np.array([0, 1, 2, 3, 4, 5], float)

        structure = eq.ColMajorSparseStructure(4, 3, ia, ja)

        structure = eq.RowMajorSparseStructure.convert_from(structure, a)

        assert_equal(structure.rows, 4)
        assert_equal(structure.cols, 3)
        assert_equal(structure.ia, [0, 1, 2, 4, 6])
        assert_equal(structure.ja, [0, 2, 0, 1, 1, 2])
        assert_equal(a, [0, 4, 1, 2, 3, 5])


if __name__ == '__main__':
    unittest.main()