from utils import *
from math import sqrt
import pandas as pd

class MAPStream:
    """
    MAP stream class.
    Contains two transition matrices, stream intensity,
    sum of transition matrices, variation coefficient and correlation coefficient.
    """

    def print_characteristics(self, matrix_name):
        """
        Prints characteristics of MAP stream:
        Average intensity
        Variation coefficient
        Correlation coefficient
        :return: None
        """

        for i, matr in enumerate(self.transition_matrices):
            print(matrix_name + '_' + str(i), ':')
            matr_print(matr)

        print('Average intensity:', self.avg_intensity)
        print('Variation coefficient:', self.c_var)
        print('Correlation coefficient:', self.c_cor)
        print('=======END=======', '\n')

    def __init__(self, transition_matr0, transition_matr1):
        """
        Constructor for MAPStream.

        :param transition_matr0: np.array or list with transition matrix 0
        :param transition_matr1: np.array or list with transition matrix 1
        """

        self.transition_matrices = [np.array(transition_matr0)]
        self.transition_matrices.append(np.array(transition_matr1))
        self.transition_matrices_sum = self.transition_matrices[0] + self.transition_matrices[1]
        gamma = system_solve(self.transition_matrices_sum)
        self.avg_intensity = r_multiply_e(np.dot(gamma, self.transition_matrices[1]))[0]
        self.dim_ = self.transition_matrices[1].shape[0]
        self.dim = self.dim_ - 1
        c_var2 = 2 * self.avg_intensity * r_multiply_e(np.dot(gamma,
                                                              la.inv(-self.transition_matrices[0])))[0] - 1
        self.c_var = sqrt(c_var2)
        self.c_cor = self.avg_intensity * (r_multiply_e(np.dot(np.dot(np.dot(gamma,
                                                                             la.inv(-self.transition_matrices[0])),
                                                                      self.transition_matrices[1]),
                                                               la.inv(-self.transition_matrices[0])))[0] - 1) / c_var2


class BMAPStream:
    """
    BMAP stream class.
    Contains list of transition matrices, stream average intensity, stream batches intensity,
    variation coefficient and correlation coefficient.
    """

    def print_characteristics(self, matrix_name):
        """
        Prints characteristics of BMAP stream:
        Matrices
        Average intensity
        Average batch intensity
        Variation coefficient
        Correlation coefficient
        :return: None
        """

        for i, matr in enumerate(self.transition_matrices):
            print(matrix_name + '_' + str(i), ':')
            matr_print(matr)

        print('Average intensity:', self.avg_intensity)
        print('Average batch intensity:', self.batch_intensity)
        print('Variation coefficient:', self.c_var)
        print('Correlation coefficient:', self.c_cor)
        print('=======END=======', '\n')

    def __init__(self, matrD_0, matrD, q, n):
        """
        Constructor for BMAPStream.

        :param matrD_0: np.array or list with matrix D_0
        :param matrD: np.array or list with matrix that will be used to generate other matrices
        :param q: float coefficient for generating other matrices
        :param n: int number of matrices to be generated (excluding matrix D_0)
        """

        self.transition_matrices = [np.array(matrD_0)]
        matrD = np.array(matrD)
        for k in range(1, n + 1):
            self.transition_matrices.append(matrD * (q ** (k - 1)) * (1 - q) / (1 - q ** 3))

        matrD_1_ = np.zeros(self.transition_matrices[0].shape)
        for matr in self.transition_matrices:
            matrD_1_ += matr
        theta = system_solve(matrD_1_)
        matrdD_1_ = np.array(copy.deepcopy(self.transition_matrices[1]))
        for i in range(2, n + 1):
            matrdD_1_ += self.transition_matrices[i] * i
        self.avg_intensity = r_multiply_e(np.dot(theta, matrdD_1_))[0]
        self.batch_intensity = r_multiply_e(np.dot(theta, -self.transition_matrices[0]))[0]
        self.dim_ = self.transition_matrices[0].shape[0]
        self.dim = self.dim_ - 1
        c_var2 = 2 * self.batch_intensity * r_multiply_e(np.dot(theta,
                                                                la.inv(-self.transition_matrices[0])))[0] - 1
        self.c_var = sqrt(c_var2)
        self.c_cor = self.batch_intensity * (r_multiply_e(np.dot(np.dot(np.dot(theta,
                                                                               la.inv(-self.transition_matrices[0])),
                                                                        matrD_1_ - self.transition_matrices[0]),
                                                                 la.inv(-self.transition_matrices[0])))[0] - 1) / c_var2


class PHStream:
    """
    PH stream class.
    Contains representation vector, representation matrix, representation matrix_0,
    stream control Markov chain dimensions, stream intensity,
    variation coefficient and correlation coefficient.
    """

    def print_characteristics(self, matrix_name, vector_name):
        """
        Prints characteristics of PH stream:
        Matrix
        Vector
        Average intensity
        Variation coefficient
        Correlation coefficient
        :return: None
        """

        print(matrix_name, ':')
        matr_print(self.repres_matr)
        print(vector_name, ':')
        print(self.repres_vect[0])

        print('Average intensity:', self.avg_intensity)
        print('Variation coefficient:', self.c_var)
        print('=======END=======', '\n')

    def __init__(self, repres_vect, repres_matr):
        """
        Constructor for PHStream

        :param repres_vect: np.array or list with representation vector
        :param repres_matr: np.array or list with representation matrix
        """

        self.repres_vect = np.array(repres_vect)
        self.repres_matr = np.array(repres_matr)
        self.repres_matr_0 = -r_multiply_e(self.repres_matr)
        self.avg_intensity = -la.inv(r_multiply_e(np.dot(self.repres_vect,
                                                         la.inv(self.repres_matr))))[0, 0]
        self.dim = self.repres_matr.shape[0]
        self.dim_ = self.dim + 1
        b1 = r_multiply_e(np.dot(self.repres_vect,
                                 la.inv(-self.repres_matr)))[0]
        b2 = 2 * r_multiply_e(np.dot(self.repres_vect,
                                 np.linalg.matrix_power(-self.repres_matr, -2)))[0]
        c_var2 = (b2 - b1 ** 2) / b1 ** 2
        self.c_var = sqrt(c_var2)
