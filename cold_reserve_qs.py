import numpy as np
import sys

from streams import *


np.set_printoptions(threshold=np.inf, suppress=True, formatter={'float': '{: 0.8f}'.format}, linewidth=75)


class ColdReserveQueueingSystem:
    """
    Class describing Unreliable Queueing System with Infinite Buffer and
    Reserve Service Unit with Cold Reserve.
    In Kendall designation (informal): BMAP|PH1|PH2|Switch12PH|Switch21PH|BreakMAP|RecoverPH

    Default values are for stationary Poisson case.
    """

    def __init__(self):
        self._p_num = 100
        self._eps_G = 10 ** (-6)
        self._eps_Phi = 10 ** (-6)
        self._eps_p_i = 10 ** (-6)
        self.queries_stream = MAPStream([[-19]], [[19]])
        self.break_stream = MAPStream([[-0.00001]], [[0.00001]])
        self.serv_unit1_stream = PHStream([[1]], [[-20]])
        self.serv_unit2_stream = PHStream([[1]], [[-5]])
        self.switch1_2_stream = PHStream([[0.05, 0.95]], [[-1.86075, 0.], [0., -146.9994]])
        self.switch2_1_stream = PHStream([[0.05, 0.95]], [[-5.58225, 0.], [0., -440.9982]])
        self.recover_stream = PHStream([[1]], [[-10000]])
        self.a = self.queries_stream.dim_ * self.break_stream.dim_
        self.n = 1

    def set_BMAP_queries_stream(self, matrD_0, matrD, q, n):
        self.queries_stream = BMAPStream(matrD_0, matrD, q, n)
        self.a = self.queries_stream.dim_ * self.break_stream.dim_
        self.n = n

    def set_MAP_queries_stream(self, matr0, matr1):
        self.queries_stream = MAPStream(matr0, matr1)
        self.a = self.queries_stream.dim_ * self.break_stream.dim_

    def set_MAP_break_stream(self, matr0, matr1):
        self.break_stream = MAPStream(matr0, matr1)
        self.a = self.queries_stream.dim_ * self.break_stream.dim_

    def set_PH_serv_unit1_stream(self, vect, matr):
        self.serv_unit1_stream = PHStream(vect, matr)

    def set_PH_serv_unit2_stream(self, vect, matr):
        self.serv_unit2_stream = PHStream(vect, matr)

    def set_PH_switch1_2_stream(self, vect, matr):
        self.switch1_2_stream = PHStream(vect, matr)

    def set_PH_switch2_1_stream(self, vect, matr):
        self.switch2_1_stream = PHStream(vect, matr)

    def set_PH_recover_stream(self, vect, matr):
        self.recover_stream = PHStream(vect, matr)

    def _calc_Qw_0(self):
        block00 = kronsum(self.queries_stream.transition_matrices[0], self.break_stream.transition_matrices[0])
        block03 = kron(kron(kron(np.eye(self.queries_stream.dim_),
                                 self.break_stream.transition_matrices[1]),
                            self.recover_stream.repres_vect),
                       self.switch1_2_stream.repres_vect)
        block10 = kron(np.eye(self.a), self.switch2_1_stream.repres_matr_0)
        block11 = kronsum(kronsum(self.queries_stream.transition_matrices[0],
                                  self.break_stream.transition_matrices[0]),
                          self.switch2_1_stream.repres_matr)
        block12 = kron(kron(kron(np.eye(self.queries_stream.dim_),
                                 self.break_stream.transition_matrices[1]),
                            self.recover_stream.repres_vect), e_col(self.switch2_1_stream.dim))
        block21 = kron(kron(np.eye(self.a),
                            self.recover_stream.repres_matr_0),
                       self.switch2_1_stream.repres_vect)
        block22 = kronsum(kronsum(self.queries_stream.transition_matrices[0],
                                  self.break_stream.transition_matrices_sum),
                          self.recover_stream.repres_matr)
        block30 = kron(kron(np.eye(self.a),
                            self.recover_stream.repres_matr_0),
                       e_col(self.switch1_2_stream.dim))
        block32 = kron(kron(np.eye(self.a),
                            np.eye(self.recover_stream.dim)),
                       self.switch1_2_stream.repres_matr_0)
        block33 = kronsum(kronsum(kronsum(self.queries_stream.transition_matrices[0],
                                          self.break_stream.transition_matrices_sum),
                                  self.recover_stream.repres_matr),
                          self.switch1_2_stream.repres_matr)
        block01 = np.zeros((block00.shape[0], block11.shape[1]))
        block02 = np.zeros((block00.shape[0], block12.shape[1]))
        block13 = np.zeros((block10.shape[0], block03.shape[1]))
        block20 = np.zeros((block21.shape[0], block10.shape[1]))
        block23 = np.zeros((block21.shape[0], block03.shape[1]))
        block31 = np.zeros((block30.shape[0], block11.shape[1]))

        print(block30.shape)

        matrQw_0 = np.bmat([[block00, block01, block02, block03],
                            [block10, block11, block12, block13],
                            [block20, block21, block22, block23],
                            [block30, block31, block32, block33]])
        return matrQw_0

    def _calc_Qw_k(self):
        matrQw_k = [self._calc_Qw_0()]

        for i in range(1, self.n + 1):
            block0 = kron(kron(self.queries_stream.transition_matrices[i],
                               np.eye(self.break_stream.dim_)),
                          self.serv_unit1_stream.repres_vect)
            block1 = kron(kron(kron(self.queries_stream.transition_matrices[i],
                                    np.eye(self.break_stream.dim_)),
                               self.serv_unit2_stream.repres_vect),
                          e_col(self.switch2_1_stream.dim))
            block2 = kron(kron(kron(self.queries_stream.transition_matrices[i],
                                    np.eye(self.break_stream.dim_)),
                               self.serv_unit2_stream.repres_vect),
                          np.eye(self.recover_stream.dim))
            block3 = kron(kron(self.queries_stream.transition_matrices[i],
                               np.eye(self.break_stream.dim_)),
                          np.eye(self.recover_stream.dim * self.switch1_2_stream.dim))
            matr_temp = la.block_diag(block0, block1, block2, block3)
            matrQw_k.append(matr_temp)

        # print(matrQw_k[0].shape)
        # print(matrQw_k[1].shape)

        for i in range(matrQw_k[0].shape[0]):
            sum = 0
            for matr in matrQw_k:
                sum += np.sum(matr[i])
            print('matrQw_k[' + str(i) + '] = ', sum)

        return matrQw_k

    def _calc_Qv_0(self):
        block0 = kron(np.eye(self.a),
                      self.serv_unit1_stream.repres_matr_0)
        block1 = kron(kron(np.eye(self.a), self.serv_unit2_stream.repres_matr_0),
                      np.eye(self.switch2_1_stream.dim))
        block2 = kron(kron(np.eye(self.a),
                           self.serv_unit2_stream.repres_matr_0),
                      np.eye(self.recover_stream.dim))
        block3 = np.zeros((self.a * self.recover_stream.dim * self.serv_unit1_stream.dim,
                           self.a * self.recover_stream.dim * self.serv_unit1_stream.dim))
        matrQv_0 = la.block_diag(block0, block1, block2, block3)

        return matrQv_0

    def _calc_Q_0(self):
        block0 = kron(np.eye(self.a),
                      np.dot(self.serv_unit1_stream.repres_matr_0,
                             self.serv_unit1_stream.repres_vect))
        block1 = kron(kron(np.eye(self.a),
                           np.dot(self.serv_unit2_stream.repres_matr_0,
                                  self.serv_unit2_stream.repres_vect)),
                      np.eye(self.switch2_1_stream.dim))
        block2 = kron(kron(np.eye(self.a), np.dot(self.serv_unit2_stream.repres_matr_0,
                                                  self.serv_unit2_stream.repres_vect)),
                      np.eye(self.recover_stream.dim))
        block3 = np.zeros((self.a * self.recover_stream.dim * self.switch1_2_stream.dim,
                           self.a * self.recover_stream.dim * self.switch1_2_stream.dim))
        matrQ_0 = la.block_diag(block0, block1, block2, block3)

        return matrQ_0

    def _calc_Q_1(self):
        block00 = kronsum(kronsum(self.queries_stream.transition_matrices[0],
                                  self.break_stream.transition_matrices[0]),
                          self.serv_unit1_stream.repres_matr)
        block03 = kron(kron(kron(kron(np.eye(self.queries_stream.dim_),
                                      self.break_stream.transition_matrices[1]),
                                 e_col(self.serv_unit1_stream.dim)),
                            self.recover_stream.repres_vect),
                       self.switch1_2_stream.repres_vect)
        block10 = kron(kron(kron(np.eye(self.a),
                                 e_col(self.serv_unit2_stream.dim)),
                            self.serv_unit1_stream.repres_vect),
                       self.switch2_1_stream.repres_matr_0)
        block11 = kronsum(kronsum(kronsum(self.queries_stream.transition_matrices[0],
                                          self.break_stream.transition_matrices[0]),
                                  self.serv_unit2_stream.repres_matr),
                          self.switch2_1_stream.repres_matr)
        block12 = kron(kron(kron(kron(np.eye(self.queries_stream.dim_),
                                      self.break_stream.transition_matrices[1]),
                                 np.eye(self.serv_unit2_stream.dim)),
                            self.recover_stream.repres_vect),
                       e_col(self.switch2_1_stream.dim))
        block21 = kron(kron(kron(np.eye(self.a),
                                 np.eye(self.serv_unit2_stream.dim)),
                            self.recover_stream.repres_matr_0),
                       self.switch1_2_stream.repres_vect)
        block22 = kronsum(kronsum(kronsum(self.queries_stream.transition_matrices[0],
                                          self.break_stream.transition_matrices_sum),
                                  self.serv_unit2_stream.repres_matr),
                          self.recover_stream.repres_matr)
        block30 = kron(kron(kron(np.eye(self.a),
                                 self.serv_unit1_stream.repres_vect),
                            self.recover_stream.repres_matr_0),
                       e_col(self.switch1_2_stream.dim))
        block32 = kron(kron(kron(np.eye(self.a),
                                 self.serv_unit2_stream.repres_vect),
                            np.eye(self.recover_stream.dim)),
                       self.switch1_2_stream.repres_matr_0)
        block33 = kronsum(kronsum(kronsum(self.queries_stream.transition_matrices[0],
                                          self.break_stream.transition_matrices_sum),
                                  self.recover_stream.repres_matr),
                          self.switch1_2_stream.repres_matr)
        block01 = np.zeros((block00.shape[0], block11.shape[1]))
        block02 = np.zeros((block00.shape[0], block12.shape[1]))
        block13 = np.zeros((block10.shape[0], block03.shape[1]))
        block20 = np.zeros((block21.shape[0], block10.shape[1]))
        block23 = np.zeros((block21.shape[0], block03.shape[1]))
        block31 = np.zeros((block30.shape[0], block11.shape[1]))

        matrQ_1 = np.bmat([[block00, block01, block02, block03],
                           [block10, block11, block12, block13],
                           [block20, block21, block22, block23],
                           [block30, block31, block32, block33]])

        return matrQ_1

    def _calc_Q_k(self):
        matrQ_k = [self._calc_Q_0(), self._calc_Q_1()]
        for k in range(2, self.n + 2):
            block0 = kron(self.queries_stream.transition_matrices[k - 1],
                          np.eye(self.break_stream.dim_ * self.serv_unit1_stream.dim))
            block1 = kron(self.queries_stream.transition_matrices[k - 1],
                          np.eye(self.break_stream.dim_ * self.serv_unit2_stream.dim * self.switch2_1_stream.dim))
            block2 = kron(self.queries_stream.transition_matrices[k - 1],
                          np.eye(self.break_stream.dim_ * self.serv_unit2_stream.dim * self.recover_stream.dim))
            block3 = kron(self.queries_stream.transition_matrices[k - 1],
                          np.eye(self.break_stream.dim_ * self.recover_stream.dim * self.switch1_2_stream.dim))
            matr_temp = la.block_diag(block0, block1, block2, block3)
            matrQ_k.append(matr_temp)
        return matrQ_k

    def ergodicity_check(self, matrQ_k):
        """
        Checks ergodicity condition of given system

        :param matrQ_k: np.array
        :return: bool, True if ergodicity condition is fulfilled
        """
        matr_Q_1_ = copy.deepcopy(matrQ_k[0])
        for l in range(1, self.n + 2):
            matr_Q_1_ += matrQ_k[l]

        matr_dQ_1_ = copy.deepcopy(matrQ_k[1])
        for l in range(2, self.n + 2):
            matr_dQ_1_ += l * matrQ_k[l]

        vect_y = system_solve(matr_Q_1_)

        vect_e = np.array([[1.] for _ in range(matr_dQ_1_.shape[1])])
        ergodicity = np.dot(np.dot(vect_y, matr_dQ_1_), vect_e)[0, 0]

        return ergodicity < 0, ergodicity, vect_y

    def _iter_G(self, matrG_prev, matrQ_1_neg_inv, matrQ_k):
        """
        Calculates matrix G according to algorithm

        :param matrG_prev: np.array with matrix G calculated on previous iteration
        :return: np.array with matrix G after iteration
        """
        temp_sum = np.array(copy.deepcopy(matrQ_k[0]))
        for k in range(2, self.n + 2):
            temp_sum += np.dot(matrQ_k[k], np.linalg.matrix_power(matrG_prev, k))
        matrG_new = np.dot(matrQ_1_neg_inv, temp_sum)
        return matrG_new

    def _calc_G(self):
        matrG_old = np.eye(self.queries_stream.transition_matrices[1].shape[0])
        matrG = self._iter_G(matrG_old)

        i = 1
        while la.norm(matrG - matrG_old, ord=np.inf) >= self.eps_G:
            matrG_old = matrG
            matrG = self._iter_G(matrG_old)
            i += 1
        return matrG

    def _calc_G_0(self, matrG, matrQ_k, matrQv_0):
        temp_sum = np.array(copy.deepcopy(matrQ_k[1]))
        for k in range(2, self.n + 2):
            temp_sum += np.dot(matrQ_k[k], np.linalg.matrix_power(matrG, k - 1))
        matrG_0 = la.inv(temp_sum)
        matrG_0 = -np.dot(matrG_0, matrQv_0)
        return matrG_0

    def _calc_Q_il(self, matrQw_k, matrQ_k, matrG, matrG_0):
        matrQ_il = []
        for i in range(0, self._p_num):
            matrQ_il.append([])
            if i == 0:
                for l in range(0, self.n + 1):
                    # здесь до n, т.к. нет больше матриц Q_k
                    temp_matr = np.array(copy.deepcopy(matrQw_k[l]))
                    for k in range(l + 1, self.n + 1):
                        mult_matr = np.array(copy.deepcopy(matrQw_k[k]))
                        for kk in range(k - 1, l - 1, -1):
                            if kk == 0:
                                mult_matr = np.dot(mult_matr, matrG_0)
                            else:
                                mult_matr = np.dot(mult_matr, matrG)
                        # print("mult matr 0 shape : ", mult_matr.shape)
                        temp_matr += mult_matr
                    matrQ_il[i].append(temp_matr)
                for l in range(self.n + 1, self._p_num):
                    matrQ_il[i].append(np.zeros(matrQw_k[1].shape))
            else:
                for l in range(0, self._p_num):
                    if l >= i and (l - i) <= (self.n + 1):
                        if (l - i + 1) <= (self.n + 1):
                            temp_matr = np.array(copy.deepcopy(matrQ_k[l - i + 1]))
                        else:
                            temp_matr = np.zeros(matrQ_k[0].shape)

                        for k in range(l + 1, self._p_num):  # sum from l+1 to inf
                            if (k - i + 1) <= (self.n + 1):
                                mult_matr = np.array(copy.deepcopy(matrQ_k[k - i + 1]))
                                for kk in range(l, k):
                                    mult_matr = np.dot(mult_matr, matrG)

                                temp_matr += mult_matr
                        matrQ_il[i].append(temp_matr)
                    else:
                        matrQ_il[i].append(np.zeros(matrQ_k[0].shape))
        return matrQ_il

    def _calc_Phi_l(self, matrQ_il):
        matrPhi_0 = np.eye(matrQ_il[0][0].shape[0])
        matrPhi_l = [matrPhi_0]
        for l in range(1, self._p_num):
            temp_matr = np.dot(np.dot(matrPhi_l[0], matrQ_il[0][l]), la.inv(-matrQ_il[l][l]))
            for i in range(1, l):
                # print(matrPhi_l[i].dot(matrQ_il[i][l]).dot(la.inv(-matrQ_il[l][l])).shape)
                temp_matr += np.dot(np.dot(matrPhi_l[i], matrQ_il[i][l]), la.inv(-matrQ_il[l][l]))
            matrPhi_l.append(temp_matr)
        return matrPhi_l

    def _calc_p_0(self, matrQ_il, matrPhi_l):
        # Вычисление p_0
        matr_a = np.array(- copy.deepcopy(matrQ_il[0][0]))
        vect_eaR = e_col(matrPhi_l[0].shape[1])
        for i in range(1, self._p_num):
            vect_e = e_col(matrPhi_l[i].shape[1])
            vect_eaR += np.dot(matrPhi_l[i], vect_e)

        for i in range(matr_a.shape[0]):
            matr_a[i][0] = vect_eaR[i][0]

        matr_b = np.zeros((matr_a.shape[0], 1))
        matr_b[0][0] = 1.
        matr_a = np.transpose(matr_a)
        p0 = np.transpose(la.solve(matr_a, matr_b))

        return p0

    def _calc_p_l(self, matrQ_il, matrPhi_l):
        p0 = self._calc_p_0(matrQ_il, matrPhi_l)
        vect_p_l = [p0]
        p_sums = [np.sum(p0)]
        # print('p0 = ', vect_p_l[0][0])
        print('sum0 = ', p_sums[0])
        for l in range(1, self._p_num):
            vect_p_l.append(np.dot(vect_p_l[0], matrPhi_l[l]))
            p_sums.append(np.sum(vect_p_l[l]))
            # print('p' + str(l) + ' = ', vect_p_l[l][0])
            print('sum' + str(l) + ' = ', p_sums[l])
        print('sum = ', np.sum(p_sums))

        return vect_p_l

    def calc_system_load(self, matrQ_0, vect_y):
        vect_e = e_col(matrQ_0.shape[1])
        denom = np.dot(np.dot(vect_y, matrQ_0), vect_e)[0]
        rho = self.queries_stream.intensity / denom
        return rho

    def calc_system_capacity(self):
        block00 = kronsum(self.break_stream.transition_matrices[0],
                          self.serv_unit1_stream.repres_matr) + kron(np.eye(self.break_stream.dim_),
                                                                     np.dot(self.serv_unit1_stream.repres_matr_0,
                                                                            self.serv_unit1_stream.repres_vect)
                                                                     )
        block03 = kron(kron(kron(self.break_stream.transition_matrices[1],
                                 e_col(self.serv_unit1_stream.dim)),
                            self.recover_stream.repres_vect),
                       self.switch1_2_stream.repres_vect)
        block10 = kron(kron(np.eye(self.break_stream.dim_),
                            np.dot(e_col(self.serv_unit2_stream.dim),
                                   self.serv_unit1_stream.repres_vect)),
                       self.switch2_1_stream.repres_matr_0)
        block11 = kronsum(kronsum(self.break_stream.transition_matrices[0],
                                  self.serv_unit2_stream.repres_matr),
                          self.switch2_1_stream.repres_matr) + kron(kron(np.eye(self.break_stream.dim_),
                                                                         np.dot(self.serv_unit2_stream.repres_matr_0,
                                                                                self.serv_unit2_stream.repres_vect)),
                                                                    np.eye(self.switch2_1_stream.dim))
        block12 = kron(kron(kron(self.break_stream.transition_matrices[1],
                                 np.eye(self.serv_unit2_stream.dim)),
                            self.recover_stream.repres_vect),
                       e_col(self.switch2_1_stream.dim))
        block21 = kron(kron(kron(np.eye(self.break_stream.dim_),
                                 np.eye(self.serv_unit2_stream.dim)),
                            self.recover_stream.repres_matr_0),
                       self.switch1_2_stream.repres_vect)
        block22 = kronsum(kronsum(self.break_stream.transition_matrices_sum,
                                  self.serv_unit2_stream.repres_matr),
                          self.recover_stream.repres_matr) + kron(kron(np.eye(self.break_stream.dim_),
                                                                       np.dot(self.serv_unit2_stream.repres_matr_0,
                                                                              self.serv_unit2_stream.repres_vect)),
                                                                  np.eye(self.recover_stream.dim))
        block30 = kron(kron(kron(np.eye(self.break_stream.dim_),
                                 self.serv_unit1_stream.repres_vect),
                            self.recover_stream.repres_matr_0),
                       e_col(self.switch1_2_stream.dim))
        block32 = kron(kron(kron(np.eye(self.break_stream.dim_),
                                 self.serv_unit2_stream.repres_vect),
                            np.eye(self.recover_stream.dim)),
                       self.switch1_2_stream.repres_matr_0)
        block33 = kronsum(kronsum(self.break_stream.transition_matrices_sum,
                                  self.recover_stream.repres_matr),
                          self.switch1_2_stream.repres_matr)
        block01 = np.zeros((block00.shape[0], block11.shape[1]))
        block02 = np.zeros((block00.shape[0], block12.shape[1]))
        block13 = np.zeros((block10.shape[0], block03.shape[1]))
        block20 = np.zeros((block21.shape[0], block00.shape[1]))
        block23 = np.zeros((block21.shape[0], block03.shape[1]))
        block31 = np.zeros((block30.shape[0], block11.shape[1]))
        matrGamma = np.bmat([[block00, block01, block02, block03],
                             [block10, block11, block12, block13],
                             [block20, block21, block22, block23],
                             [block30, block31, block32, block33]])

        x = system_solve(matrGamma)

        # print('x = ', x)
        x1 = x[0:self.break_stream.dim_ * self.serv_unit1_stream.dim]
        x2 = x[self.break_stream.dim_ * self.serv_unit1_stream.dim:self.break_stream.dim_ * self.serv_unit1_stream.dim + self.break_stream.dim_ * self.serv_unit2_stream.dim * self.switch2_1_stream.dim]
        x3 = x[self.break_stream.dim_ * self.serv_unit1_stream.dim + self.break_stream.dim_ * self.serv_unit2_stream.dim * self.switch2_1_stream.dim:self.break_stream.dim_ * self.serv_unit1_stream.dim + self.break_stream.dim_ * self.serv_unit2_stream.dim * (self.switch2_1_stream.dim + self.recover_stream.dim)]

        # print('x1 = ', x1)
        # print('x2 = ', x2)
        # print('x3 = ', x3)

        e_V_ = e_col(self.break_stream.dim)
        e_R = e_col(self.recover_stream.dim)
        pi1 = x1.dot(kron(e_V_, np.eye(self.serv_unit1_stream.dim)))
        pi2 = x2.dot(kron(kron(e_V_, np.eye(self.serv_unit2_stream.dim)), e_col(self.switch2_1_stream.dim)))
        pi3 = x3.dot(kron(kron(e_V_, np.eye(self.serv_unit2_stream.dim)), e_R))

        varrho = np.dot(pi1, self.serv_unit1_stream.repres_matr_0) + np.dot((pi2 + pi3), self.serv_unit2_stream.repres_matr_0)

        return varrho[0]

    def get_prod_func(self, vectors):
        """
        Generates production function P(1)
        :param vectors: iterable of np.arrays
        :return: np.array with prod function P(1)
        """
        vect_1_ = vectors[1]
        for l in range(2, self._p_num):
            vect_1_ += vectors[l]
        return vect_1_

    def get_prod_func_deriv(self, vectors):
        vect_dP_1_ = vectors[1]
        for l in range(2, self._p_num):
            vect_dP_1_ += l * vectors[l]
        return vect_dP_1_

    def get_prod_func_2deriv(self, vectors):
        vect_ddP_1_ = copy.deepcopy(vectors[2]) * 2
        for l in range(3, self._p_num):
            vect_ddP_1_ += l * (l - 1) * vectors[l]
        return vect_ddP_1_

    def calc_avg_queries_num(self, vect_dP_1_):
        L = r_multiply_e(vect_dP_1_)[0]
        return L

    def calc_queries_num_dispersion(self, vect_ddP_1_, avg_queries_num):
        dispV = r_multiply_e(vect_ddP_1_)[0] + avg_queries_num - avg_queries_num ** 2
        return dispV[0]

    def calc_prob_1_work_serves(self, vect_P_1_):
        temp_matr = np.dot(vect_P_1_,
                           la.block_diag(np.eye(self.a * self.serv_unit1_stream.dim),
                                         np.zeros((self.a * (self.serv_unit2_stream.dim * self.switch2_1_stream.dim + self.serv_unit2_stream.dim * self.recover_stream.dim + self.recover_stream.dim * self.switch1_2_stream.dim),
                                                   self.a * (self.serv_unit2_stream.dim * self.switch2_1_stream.dim + self.serv_unit2_stream.dim * self.recover_stream.dim + self.recover_stream.dim * self.switch1_2_stream.dim)))))
        return r_multiply_e(temp_matr)[0, 0]

    def calc_prob_1_broken_2_serves(self, vect_P_1_):
        temp_matr = np.dot(vect_P_1_,
                           la.block_diag(np.zeros((self.a * self.serv_unit1_stream.dim,
                                                   self.a * self.serv_unit1_stream.dim)),
                                         np.eye(self.a * self.serv_unit2_stream.dim * (self.switch2_1_stream.dim + self.recover_stream.dim)),
                                                    np.zeros((self.a * self.recover_stream.dim * self.switch1_2_stream.dim,
                                                              self.a * self.recover_stream.dim * self.switch1_2_stream.dim))))
        return r_multiply_e(temp_matr)[0, 0]

    def calc_prob_1_broken_switch_1_2(self, vect_P_1_):
        temp_matr = np.dot(vect_P_1_,
                           la.block_diag(np.zeros((self.a * (self.serv_unit1_stream.dim + self.serv_unit2_stream.dim * (self.switch2_1_stream.dim + self.recover_stream.dim)),
                                                   self.a * (self.serv_unit1_stream.dim + self.serv_unit2_stream.dim * (self.switch2_1_stream.dim + self.recover_stream.dim)))),
                                                    np.eye(self.a * self.recover_stream.dim * self.switch1_2_stream.dim)))
        return r_multiply_e(temp_matr)

    def calc_prob_1_work_switch_2_1(self, vect_P_1_):
        temp_matr = np.dot(vect_P_1_, la.block_diag(np.zeros((self.a * self.serv_unit1_stream.dim, self.a * self.serv_unit1_stream.dim)),
                                                    np.eye(self.a * self.serv_unit2_stream.dim * self.switch2_1_stream.dim),
                                                    np.zeros((self.a * (self.recover_stream.dim * self.switch1_2_stream.dim + self.serv_unit2_stream.dim * self.recover_stream.dim),
                                                              self.a * (self.recover_stream.dim * self.switch1_2_stream.dim + self.serv_unit2_stream.dim * self.recover_stream.dim)))))
        return r_multiply_e(temp_matr)[0, 0]

    def calc_prob_1_available(self, vect_p_l, vect_P_1_):
        temp_matr1 = np.dot(vect_P_1_, la.block_diag(np.eye(self.a * self.serv_unit1_stream.dim),
                                                     np.zeros((self.a * (self.serv_unit2_stream.dim * self.switch2_1_stream.dim + self.serv_unit2_stream.dim * self.recover_stream.dim + self.recover_stream.dim * self.switch1_2_stream.dim),
                                                               self.a * (self.serv_unit2_stream.dim * self.switch2_1_stream.dim + self.serv_unit2_stream.dim * self.recover_stream.dim + self.recover_stream.dim * self.switch1_2_stream.dim)))))
        temp_matr = r_multiply_e(temp_matr1)

        temp_matr2 = np.dot(vect_p_l[0],
                            la.block_diag(np.eye(self.a),
                                          np.zeros((self.a * (self.switch2_1_stream.dim + self.recover_stream.dim + self.recover_stream.dim * self.switch1_2_stream.dim),
                                                    self.a * (self.switch2_1_stream.dim + self.recover_stream.dim + self.recover_stream.dim * self.switch1_2_stream.dim)))))
        temp_matr += r_multiply_e(temp_matr2)
        return temp_matr[0, 0]

    def calc_prob_1_unavail_2_avail(self, vect_p_l, vect_P_1_):
        temp_matr1 = np.dot(vect_P_1_, la.block_diag(np.zeros((self.a * self.serv_unit1_stream.dim, self.a * self.serv_unit1_stream.dim)),
                                                     np.eye(self.a * self.serv_unit2_stream.dim * (self.switch2_1_stream.dim + self.recover_stream.dim)),
                                                     np.zeros((self.a * self.recover_stream.dim * self.switch1_2_stream.dim,
                                                               self.a * self.recover_stream.dim * self.switch1_2_stream.dim))))
        temp_matr = r_multiply_e(temp_matr1)

        temp_matr2 = np.dot(vect_p_l[0],
                            la.block_diag(np.zeros((self.a,
                                                    self.a)),
                                          np.eye(self.a * (self.switch2_1_stream.dim + self.recover_stream.dim)),
                                          np.zeros((self.a * self.recover_stream.dim * self.switch1_2_stream.dim,
                                                    self.a * self.recover_stream.dim * self.switch1_2_stream.dim))))
        temp_matr += r_multiply_e(temp_matr2)

        return temp_matr[0, 0]

    def calc_prob_1_2_unavail(self, vect_p_l, vect_P_1_):
        temp_matr1 = np.dot(vect_P_1_, la.block_diag(np.zeros((self.a * (self.serv_unit1_stream.dim + self.serv_unit2_stream.dim * (self.switch2_1_stream.dim + self.recover_stream.dim)),
                                                               self.a * (self.serv_unit1_stream.dim + self.serv_unit2_stream.dim * (self.switch2_1_stream.dim + self.recover_stream.dim)))),
                                                     np.eye(self.a * self.recover_stream.dim * self.switch1_2_stream.dim)))
        temp_matr = r_multiply_e(temp_matr1)

        temp_matr2 = np.dot(vect_p_l[0],
                            la.block_diag(np.zeros((self.a * (1 + self.switch2_1_stream.dim + self.recover_stream.dim),
                                                    self.a * (1 + self.switch2_1_stream.dim + self.recover_stream.dim))),
                                          np.eye(self.a * self.recover_stream.dim * self.switch1_2_stream.dim)))
        temp_matr += r_multiply_e(temp_matr2)

        return temp_matr[0, 0]

    def calc_avg_switch_1_2_num(self, vect_p_l, vect_P_1_):
        temp_matr1 = np.dot(vect_P_1_,
                            la.block_diag(kron(kron(np.eye(self.queries_stream.dim_),
                                                    self.break_stream.transition_matrices[1]),
                                               np.eye(self.serv_unit1_stream.dim)),
                                          np.zeros((self.a * (self.serv_unit2_stream.dim * (self.switch2_1_stream.dim + self.recover_stream.dim) + self.recover_stream.dim * self.switch1_2_stream.dim),
                                                    self.a * (self.serv_unit2_stream.dim * (self.switch2_1_stream.dim + self.recover_stream.dim) + self.recover_stream.dim * self.switch1_2_stream.dim)))))
        temp_matr = r_multiply_e(temp_matr1)

        temp_matr2 = np.dot(vect_p_l[0],
                            la.block_diag(kron(np.eye(self.queries_stream.dim_),
                                               self.break_stream.transition_matrices[1]),
                                          np.zeros((self.a * (self.switch2_1_stream.dim + self.recover_stream.dim + self.recover_stream.dim * self.switch1_2_stream.dim),
                                                    self.a * (self.switch2_1_stream.dim + self.recover_stream.dim + self.recover_stream.dim * self.switch1_2_stream.dim)))))
        temp_matr += r_multiply_e(temp_matr2)

        return temp_matr[0, 0]

    def calc_avg_switch_2_1_num(self, vect_p_l, vect_P_1_):
        temp_matr1 = np.dot(vect_P_1_,
                            la.block_diag(np.zeros((self.a * (self.serv_unit1_stream.dim + self.serv_unit2_stream.dim * self.switch2_1_stream.dim),
                                                    self.a * (self.serv_unit1_stream.dim + self.serv_unit2_stream.dim * self.switch2_1_stream.dim))),
                                          kron(np.eye(self.a * self.serv_unit2_stream.dim),
                                               self.recover_stream.repres_matr),
                                          np.zeros((self.a * self.recover_stream.dim * self.switch1_2_stream.dim,
                                                    self.a * self.recover_stream.dim * self.switch1_2_stream.dim))))
        temp_matr = r_multiply_e(temp_matr1)
        temp_matr2 = np.dot(vect_p_l[0],
                            la.block_diag(np.zeros((self.a * (1 + self.switch2_1_stream.dim),
                                                    self.a * (1 + self.switch2_1_stream.dim))),
                                          kron(np.eye(self.a),
                                               self.recover_stream.repres_matr),
                                          np.zeros((self.a * self.recover_stream.dim * self.switch1_2_stream.dim,
                                                    self.a * self.recover_stream.dim * self.switch1_2_stream.dim))))
        temp_matr += r_multiply_e(temp_matr2)
        return temp_matr[0, 0]

    def calc_characteristics(self):
        matrQw_k = self._calc_Qw_k()
        matrQ_k = self._calc_Q_k()
        matrQv_0 = self._calc_Qv_0()

        # Check ergodicity condition
        is_ergodicic, ergodicity, vect_y = self.ergodicity_check(matrQ_k)
        if not is_ergodicic:
            print('Условие эргодичности не выполняется!', file=sys.stderr)
            print(ergodicity, '>= 0', file=sys.stderr)
            print('Выполнение программы остановлено!', file=sys.stderr)
            exit()  # Replace with return or smth else to detect iteration
        else:
            print('Условие эргодичности выполняется!')
            print(ergodicity, '< 0')

        matrG = self._calc_G()
        matrG_0 = self._calc_G_0(matrG, matrQ_k, matrQv_0)

        matrQ_il = self._calc_Q_il(matrQw_k, matrQ_k, matrG, matrG_0)

        matrPhi_l = self._calc_Phi_l(matrQ_il)

        vect_p_l = self._calc_p_l(matrQ_il, matrPhi_l)

        system_load = self.calc_system_load(matrQ_k[0], vect_y)

        system_capacity = self.calc_system_capacity()

        vect_P_1_ = self.get_prod_func(vect_p_l)
        vect_dP_1_ = self.get_prod_func_deriv(vect_p_l)
        vect_ddP_1_ = self.get_prod_func_2deriv(vect_p_l)

        avg_queries_num = self.calc_avg_queries_num(vect_dP_1_)

        queries_num_dispersion = self.calc_queries_num_dispersion(vect_ddP_1_, avg_queries_num)

        prob_1_work_serves = self.calc_prob_1_work_serves(vect_P_1_)

        prob_1_broken_2_serves = self.calc_prob_1_broken_2_serves(vect_P_1_)

        prob_1_broken_switch_1_2 = self.calc_prob_1_broken_switch_1_2(vect_P_1_)

        prob_1_work_switch_2_1 = self.calc_prob_1_work_switch_2_1(vect_P_1_)

        prob_1_available = self.calc_prob_1_available(vect_p_l, vect_P_1_)

        prob_1_unavail_2_avail = self.calc_prob_1_unavail_2_avail(vect_p_l, vect_P_1_)

        prob_1_2_unavail = self.calc_prob_1_2_unavail(vect_p_l, vect_P_1_)

        avg_switch_1_2_num = self.calc_avg_switch_1_2_num(vect_p_l, vect_P_1_)

        avg_switch_2_1_num = self.calc_avg_switch_2_1_num(vect_p_l, vect_P_1_)

if __name__ == '__main__':
    system = ColdReserveQueueingSystem()
    print(system.a)
    system.set_BMAP_queries_stream(np.array([[-86., 0.01], [0.02, -2.76]]),
                                   np.array([[85, 0.99], [0.2, 2.54]]),
                                   n=3, q=0.8)
    print(system.a)
