import eqlib as eq

import numpy as np
import pytest

from numpy.testing import assert_equal

if __name__ == '__main__':
    import sys
    import os
    print(f'pid: {os.getpid()}')
    pytest.main(sys.argv)


# csr


@pytest.fixture
def csr_rectangular():
    ia = [0, 2, 4, 6]
    ja = [0, 2, 2, 3, 1, 3]

    return eq.RowMajorSparseStructure(3, 4, ia, ja)


@pytest.fixture
def csr_square():
    ia = [0, 3, 5, 6, 7]
    ja = [0, 1, 2, 1, 2, 2, 3]

    return eq.RowMajorSparseStructure(4, 4, ia, ja)


def test_csr_rectangular(csr_rectangular):
    structure = csr_rectangular

    assert_equal(structure.ia, [0, 2, 4, 6])
    assert_equal(structure.ja, [0, 2, 2, 3, 1, 3])

    for idx, (i, j) in enumerate([(0, 0), (0, 2), (1, 2), (1, 3), (2, 1), (2, 3)]):
        assert_equal(structure.get_index(i, j), idx)

    assert_equal(structure.get_index(1, 0), -1)

    assert_equal(structure.density, 6 / 12)


def test_csr_square(csr_square):
    structure = csr_square

    assert_equal(structure.ia, [0, 3, 5, 6, 7])
    assert_equal(structure.ja, [0, 1, 2, 1, 2, 2, 3])

    for idx, (i, j) in enumerate([(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2), (3, 3)]):
        assert_equal(structure.get_index(i, j), idx)

    assert_equal(structure.get_index(1, 0), -1)

    assert_equal(structure.density, 7 / 16)


def test_csr_convert_from(csr_rectangular):
    a = np.array([0, 1, 2, 3, 4, 5], float)

    structure = eq.ColMajorSparseStructure.convert_from(csr_rectangular, a)

    assert_equal(structure.rows, 3)
    assert_equal(structure.cols, 4)
    assert_equal(structure.ia, [0, 1, 2, 4, 6])
    assert_equal(structure.ja, [0, 2, 0, 1, 1, 2])
    assert_equal(a, [0, 4, 1, 2, 3, 5])


def test_csr_from_pattern():
    structure = eq.RowMajorSparseStructure.from_pattern(3, 4, [{0, 2}, {2, 3}, {1, 3}])

    assert_equal(structure.rows, 3)
    assert_equal(structure.cols, 4)
    assert_equal(structure.ia, [0, 2, 4, 6])
    assert_equal(structure.ja, [0, 2, 2, 3, 1, 3])


def test_csr_for_each(csr_rectangular):
    results = []

    csr_rectangular.for_each(lambda idx, row, col: results.append((idx, row, col)))

    assert_equal(results, [(0, 0, 0), (0, 2, 1), (1, 2, 2), (1, 3, 3), (2, 1, 4), (2, 3, 5)])


# csc


@pytest.fixture
def csc_rectangular():
    ia = [0, 2, 4, 6]
    ja = [0, 2, 2, 3, 1, 3]

    return eq.ColMajorSparseStructure(4, 3, ia, ja)


@pytest.fixture
def csc_square():
    ia = [0, 3, 5, 6, 7]
    ja = [0, 1, 2, 1, 2, 2, 3]

    return eq.ColMajorSparseStructure(4, 4, ia, ja)


def test_csc_rectangular(csc_rectangular):
    structure = csc_rectangular

    assert_equal(structure.ia, [0, 2, 4, 6])
    assert_equal(structure.ja, [0, 2, 2, 3, 1, 3])

    for idx, (i, j) in enumerate([(0, 0), (2, 0), (2, 1), (3, 1), (1, 2), (3, 2)]):
        assert_equal(structure.get_index(i, j), idx)

    assert_equal(structure.get_index(0, 1), -1)

    assert_equal(structure.density, 6 / 12)


def test_csc_square(csc_square):
    structure = csc_square

    assert_equal(structure.ia, [0, 3, 5, 6, 7])
    assert_equal(structure.ja, [0, 1, 2, 1, 2, 2, 3])

    for idx, (i, j) in enumerate([(0, 0), (1, 0), (2, 0), (1, 1), (2, 1), (2, 2), (3, 3)]):
        assert_equal(structure.get_index(i, j), idx)

    assert_equal(structure.get_index(0, 1), -1)

    assert_equal(structure.density, 7 / 16)


def test_csc_convert_from(csc_rectangular):
    a = np.array([0, 1, 2, 3, 4, 5], float)

    structure = eq.RowMajorSparseStructure.convert_from(csc_rectangular, a)

    assert_equal(structure.rows, 4)
    assert_equal(structure.cols, 3)
    assert_equal(structure.ia, [0, 1, 2, 4, 6])
    assert_equal(structure.ja, [0, 2, 0, 1, 1, 2])
    assert_equal(a, [0, 4, 1, 2, 3, 5])


def test_csc_from_pattern(csc_rectangular):
    structure = eq.ColMajorSparseStructure.from_pattern(4, 3, [{0, 2}, {2, 3}, {1, 3}])

    assert_equal(structure.rows, 4)
    assert_equal(structure.cols, 3)
    assert_equal(structure.ia, [0, 2, 4, 6])
    assert_equal(structure.ja, [0, 2, 2, 3, 1, 3])


def test_csc_for_each(csc_rectangular):
    results = []

    csc_rectangular.for_each(lambda idx, row, col: results.append((idx, row, col)))

    assert_equal(results, [(0, 0, 0), (2, 0, 1), (2, 1, 2), (3, 1, 3), (1, 2, 4), (3, 2, 5)])
