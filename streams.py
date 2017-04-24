from utils import *


class MAPStream:
    """
    MAP stream class.
    Contains two transition matrices, stream intensity, variance coefficient and correlation coefficient.
    """
    def __init__(self, transition_matr0, transition_matr1):
        """
        Constructor for MAPStream.

        :param transition_matr0: np.array or list with transition matrix 0
        :param transition_matr1: np.array or list with transition matrix 1
        """
        self.transition_matrices = [np.array(transition_matr0)]
        self.transition_matrices.append(np.array(transition_matr1))
        gamma = system_solve(self.transition_matrices[0] + self.transition_matrices[1])
        self.intensity = r_multiply_e(np.dot(gamma, self.transition_matrices[1]))
        self.dim_ = self.transition_matrices[1].shape[0]
        self.dim = self.dim_ - 1
        self.transition_matrices_sum = self.transition_matrices[0] + self.transition_matrices[1]


class BMAPStream:
    """
    BMAP stream class.
    Contains list of transition matrices, stream average intensity, stream batches intensity,
    variance coefficient and correlation coefficient.
    """
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
        for k in range(1, n+1):
            self.transition_matrices.append(matrD * (q ** (k-1)) * (1 - q) / (1 - q ** 3))

        matrD_1_ = np.zeros(self.transition_matrices[0].shape)
        for matr in self.transition_matrices:
            matrD_1_ += matr
        theta = system_solve(matrD_1_)
        matrdD_1_ = np.array(copy.deepcopy(self.transition_matrices[1]))
        for i in range(2, n + 1):
            matrdD_1_ += self.transition_matrices[i] * i
        self.avg_intensity = r_multiply_e(np.dot(theta, matrdD_1_))
        self.batch_intensity = r_multiply_e(np.dot(theta, -self.transition_matrices[0]))
        self.dim_ = self.transition_matrices[0].shape[0]
        self.dim = self.dim_ - 1


class PHStream:
    """
    PH stream class.
    Contains representation vector and representation matrix, stream intensity.
    """
    def __init__(self, repres_vect, repres_matr):
        """
        Constructor for PHStream

        :param repres_vect: np.array or list with representation vector
        :param repres_matr: np.array or list with representation matrix
        """
        self.repres_vect = np.array(repres_vect)
        self.repres_matr = np.array(repres_matr)
        self.repres_matr_0 = -r_multiply_e(self.repres_matr)
        self.intensity = -la.inv(r_multiply_e(np.dot(self.repres_vect, la.inv(self.repres_matr))))
        self.dim = self.repres_matr.shape[0]
        self.dim_ = self.dim + 1
