from sys import stderr
from time import sleep

import numpy as np
import scipy.linalg as la


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
    elif len(matr.shape) > 0:
        vect_e = np.array([[1.] for _ in range(matr.shape[0])])
    else:
        vect_e = 1.
    result = np.dot(matr, vect_e)

    return result


def e_col(dim):
    """
    Generates unitary vector-column of given dimension.

    :param dim: int dimension
    :return: np.array with unitary vector-column
    """
    return np.array([[1.] for _ in range(dim)])


def matr_print(matr):
    """
    Pretty prints the given matrix.

    :param matr: iterable with matrix to print
    :return: None
    """

    s = [[str(e) for e in row] for row in matr]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    print('\n'.join(table))


def linux_check_cpu_temperature():
    try:
        sleep_flag = False
        while True:
            with open('/sys/class/thermal/thermal_zone0/temp', mode='r') as temp_file:
                cpu_temp = float(temp_file.readline()) / 1000
            if cpu_temp >= 90.:
                print('Your CPU is too hot to proceed!', file=stderr)
                print('Please wait for 60 seconds to cool CPU...', file=stderr)
                sleep_flag = True
                sleep(60)
            else:
                if sleep_flag:
                    print('CPU temp is Ok. Proceeding...', file=stderr)
                return
    except:
        return
