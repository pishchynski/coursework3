import numpy as np
import scipy.linalg as la
import copy


def kron(A, B):
    """
    Just np.linalg.kron() function.

    :param A: np.array
    :param B: np.array
    :return: np.array as a result of kronecker multiplication of A and B
    """
    return la.kron(A, B)


def kronsum(A, B):
    """
    Kronecker sum function.

    :param A: np.array
    :param B: np.array
    :return: np.array as a result of kronecker sum of A and B
    """
    if A.shape[0] != A.shape[1]:
        raise ValueError('A is not square')

    if B.shape[0] != B.shape[1]:
        raise ValueError('B is not square')

    L = kron(A, np.eye(B.shape[0]))
    R = kron(np.eye(A.shape[0]), B)

    return L + R


def system_solve(matr):
    """
    Solves system of type (vect * matr = 0) & (vect * e = 1).

    :param matr: np.array
    :return: np.array with vect
    """
    matr_a = np.array(matr)

    for i in range(matr_a.shape[0]):
        matr_a[i][0] = 1

    matr_b = np.zeros((matr_a.shape[0], 1))
    matr_b[0][0] = 1
    matr_a = np.transpose(matr_a)

    result = np.transpose(la.solve(matr_a, matr_b))

    return result[0]


def r_multiply_e(matr):
    """
    Multiplies matrix matr on unitary vector of matching shape from the right: matr * e.

    :param matr: np.array
    :return: any type with the result of multiplication
    """
    if len(matr.shape) > 1:
        vect_e = np.array([[1.] for _ in range(matr.shape[1])])
    else:
        vect_e = np.array([[1.] for _ in range(matr.shape[0])])
    result = np.dot(matr, vect_e)

    return result


def e_col(dim):
    """
    Generates unitary vector-column of given dimension.

    :param dim: int dimension
    :return: np.array with unitary vector-column
    """
    return np.array([[1.] for _ in range(dim)])
