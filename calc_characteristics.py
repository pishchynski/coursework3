import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la
import copy

import sys

np.set_printoptions(threshold=np.inf, suppress=True, formatter={'float': '{: 0.8f}'.format}, linewidth=75)


def kron(A, B):
    return la.kron(A, B)


def kronsum(A, B):
    if A.shape[0] != A.shape[1]:
        raise ValueError('A is not square')

    if B.shape[0] != B.shape[1]:
        raise ValueError('B is not square')

    L = kron(A, np.eye(B.shape[0]))
    R = kron(np.eye(A.shape[0]), B)

    return L + R


def calc_G(matrG_prev):
    temp_sum = np.array(copy.deepcopy(matrQ_k[0]))
    for k in range(2, n + 2):
        temp_sum += np.dot(matrQ_k[k], np.linalg.matrix_power(matrG_prev, k))
    matrG_new = np.dot(matrQ_1_neg_inv, temp_sum)
    return matrG_new


# Consts to be in .ini file

n = 3  # Количество матриц D_k
p_num = 150  # Количество векторов стац. вероятн. p
eps_G = 10 ** (-6)  # Точность нахождения матр. G
eps_Phi = 10 ** (-6)  # Точность нахождения матр. Phi
eps_p_i = 10 ** (-6)  # Точность нахождения векторов p
#


kappa_inv_list_list = []
P_neg_list_list = []
kappa_2_inv_list = []
khi1_khi2_list_list = []

for j_experimental in (1,):
    kappa_inv_list = list()
    P_neg_list = list()
    khi1_khi2_list = list()
    for i_experimental in (1,):

        # Входной BMAP
        matrD_0 = np.array([[-86., 0.01], [0.02, -2.76]])
        matrD = np.array([[85., 0.99], [0.2, 2.54]])
        matrD_k = [matrD_0]

        W_ = matrD_0.shape[0]
        W = W_ - 1

        q = 0.8
        for k in range(1, n + 1):
            matrD_k.append(matrD * (q ** (k - 1)) * (1 - q) / (1 - q ** 3))

        for i, matr in enumerate(matrD_k):
            None
            # print('D_{} ='.format(i), matr)
        #


        # Характеристики входного BMAP

        matrD_1_ = np.zeros(matrD_k[0].shape)  # D(1)
        for matr in matrD_k:
            matrD_1_ += matr

        # print('D(1) =', matrD_1_)

        matr_a = np.array(copy.deepcopy(matrD_1_))
        for i in range(matr_a.shape[0]):
            matr_a[i][0] = 1

        matr_b = np.zeros((matr_a.shape[0], 1))
        matr_b[0][0] = 1

        matr_a = np.transpose(matr_a)

        theta = la.solve(matr_a, matr_b).reshape(-1)  # Алгоритм проверен.

        matrdD_1_ = np.array(copy.deepcopy(matrD_k[1]))  # D'(1)

        # print('\\theta =', theta)

        for i in range(2, n + 1):
            matrdD_1_ += matrD_k[i] * i
        vect_e = np.array([[1.] for i in range(0, matrD_1_.shape[1])])

        lamD = np.dot(np.dot(theta, matrdD_1_), vect_e)[0]  # Средняя интенсивность прихода заявок
        # print('\\lambda =', lamD)

        lamDb = np.dot(np.dot(theta, -matrD_0), vect_e)[0]  # Средняя интенсивность прихода групп
        # print('\\lambda_b =', lamDb)

        c2var = np.sum(2 * lamDb * (np.dot(theta, la.inv(-matrD_0)))) - 1  # c_{var}^2
        # print('c_{var}^2 =', c2var)

        vect_e = [[1.] for _ in range(matrD_1_.shape[1])]

        # c_cor = np.sum(
        #     lamDb * np.dot(theta, la.inv(-matrD_0)), matrD_1_ - np.dot(matrD_0, la.inv(-matrD_0))) - 1 / c2var
        # print(c_cor)
        #


        # Поток поломок MAP
        matrH0 = np.array([[-8., 1.], [2., -12.]])
        matrH1 = np.array([[2., 5.], [4., 6.]])
        V_ = matrH1.shape[0]
        V = V_ - 1
        matrH = matrH0 + matrH1
        matr_a = copy.deepcopy(matrH)
        for i in range(matr_a.shape[0]):
            # print(matr_a)
            matr_a[i][0] = 1

        matr_b = np.zeros((matr_a.shape[0], 1))
        matr_b[0][0] = 1

        matr_a = np.transpose(matr_a)

        # print('H_0 =', matrH0)
        # print('H_1 =', matrH1)

        gamma = la.solve(matr_a, matr_b).reshape(-1)

        vect_e = np.array([[1.] for i in range(0, matrD_1_.shape[1])])
        h = np.dot(np.dot(gamma, matrH1), vect_e)[0]
        # print('h =', h)

        # Поток обслуживания PH1
        beta1 = np.array([[0.2, 0.8]])
        matrS1 = np.array([[-170., 15.], [40., -210.]])
        M1 = matrS1.shape[0]
        M1_ = M1 + 1
        M1_e = np.array([[1] for _ in range(matrS1.shape[1])])
        matrS1_0 = np.dot(- matrS1, M1_e)
        vect_e = np.array([[1.] for i in range(0, matrS1.shape[1])])
        # print(np.dot(beta1, la.inv(matrS1)))
        mu_1 = -la.inv(np.dot(np.dot(beta1, la.inv(matrS1)), vect_e))[0, 0]
        # print('\\mu_1 =', mu_1)
        #


        # Поток обслуживания PH2
        beta2 = np.array([[0.9, 0.1]])
        matrS2 = np.array([[-110., 80.], [10., -150.]])
        M2 = matrS2.shape[0]
        M2_ = M2 + 1
        M2_e = np.array([[1] for _ in range(matrS2.shape[1])])
        matrS2_0 = np.dot(- matrS2, M2_e)

        vect_e = np.array([[1.] for i in range(0, matrS2.shape[1])])
        mu_2 = -la.inv(np.dot(np.dot(beta2, la.inv(matrS2)), vect_e))[0, 0]
        # print('\\mu_2 =', mu_2)
        #

        matrS_w = kron(np.dot(matrS1_0, beta1),
                       np.dot(M2_e, beta2)) + kron(np.dot(M1_e, beta1), np.dot(matrS2_0, beta2))

        # Поток переключения с прибора-1 на прибор-2
        alpha1 = np.array([[0.9, 0.1]])
        matrA1 = np.array([[-220., 160.], [20., -300.]])
        L1 = matrA1.shape[0]
        L1_ = L1 + 1
        L1_e = np.array([[1] for _ in range(matrA1.shape[1])])
        matrA1_0 = - np.dot(matrA1, L1_e)

        vect_e = np.array([[1.] for i in range(0, matrA1.shape[1])])
        kappa_1 = -la.inv(np.dot(np.dot(alpha1, la.inv(matrA1)), vect_e))[0, 0]
        print('\\kappa_1 =', kappa_1)
        kappa_inv_list.append(1 / copy.deepcopy(kappa_1))
        #

        # Поток переключения с прибора-2 на прибор-1
        alpha2 = np.array([[0.05, 0.95]])
        matrA2 = np.array([[-5.58225, 0.], [0., -440.9982]]) * j_experimental
        L2 = matrA2.shape[0]
        L2_ = L2 + 1
        L2_e = np.array([[1] for _ in range(matrA2.shape[1])])
        matrA2_0 = - np.dot(matrA2, L2_e)

        vect_e = np.array([[1.] for i in range(0, matrA2.shape[1])])
        kappa_2 = -la.inv(np.dot(np.dot(alpha2, la.inv(matrA2)), vect_e))[0, 0]
        print('\\kappa_2 =', kappa_2)
        #

        # Поток ремонта PH
        tau = np.array([[0.2, 0.8]])
        matrT = np.array([[-17., 1.5], [4., -21.]])
        T_e = np.array([[1] for _ in range(matrT.shape[1])])
        matrT0 = - np.dot(matrT, T_e)

        R = matrT.shape[0]
        R_ = R + 1
        vect_e = np.array([[1.] for i in range(0, matrT.shape[1])])
        phi = -la.inv(np.dot(np.dot(tau, la.inv(matrT)), vect_e))[0, 0]
        # print('\\phi =', phi)
        #

        a = W_ * V_
        # print('a =', a)

        # Q~0
        block00 = kronsum(matrD_k[0], matrH0)
        block03 = kron(kron(kron(np.eye(W_), matrH1), tau), alpha1)
        block10 = kron(np.eye(a), matrA2_0)
        block11 = kronsum(kronsum(matrD_k[0], matrH0), matrA2)
        block12 = kron(kron(kron(np.eye(W_), matrH1), tau), L2_e)
        block21 = kron(kron(np.eye(a), matrT0), alpha2)
        block22 = kronsum(kronsum(matrD_k[0], matrH), matrT)
        block30 = kron(kron(np.eye(a), matrT0), L1_e)
        block32 = kron(kron(np.eye(a), np.eye(R)), matrA1_0)
        block33 = kronsum(kronsum(kronsum(matrD_k[0], matrH), matrT), matrA1)
        block01 = np.zeros((block00.shape[0], block11.shape[1]))
        block02 = np.zeros((block00.shape[0], block12.shape[1]))
        block13 = np.zeros((block10.shape[0], block03.shape[1]))
        block20 = np.zeros((block21.shape[0], block10.shape[1]))
        block23 = np.zeros((block21.shape[0], block03.shape[1]))
        block31 = np.zeros((block30.shape[0], block11.shape[1]))

        # print(block30.shape)

        matrQw_0 = np.bmat([[block00, block01, block02, block03],
                            [block10, block11, block12, block13],
                            [block20, block21, block22, block23],
                            [block30, block31, block32, block33]])
        # print(matrQw_0.shape)
        #

        # Q~k
        matrQw_k = [matrQw_0]
        for i in range(1, n + 1):
            block0 = kron(kron(matrD_k[i], np.eye(V_)), beta1)
            block1 = kron(kron(kron(matrD_k[i], np.eye(V_)), beta2), np.eye(L2))
            block2 = kron(kron(kron(matrD_k[i], np.eye(V_)), beta2), np.eye(R))
            block3 = kron(kron(matrD_k[i], np.eye(V_)), np.eye(R * L1))
            matr_temp = la.block_diag(block0, block1, block2, block3)
            matrQw_k.append(matr_temp)

        # print(matrQw_k[0].shape)
        # print(matrQw_k[1].shape)

        for i in range(matrQw_k[0].shape[0]):
            sum = 0
            for matr in matrQw_k:
                sum += np.sum(matr[i])
                # print('matrQw_k[' + str(i) + '] = ', sum)
        #

        # Q^0
        block0 = kron(np.eye(a), matrS1_0)
        block1 = kron(kron(np.eye(a), matrS2_0), np.eye(L2))
        block2 = kron(kron(np.eye(a), matrS2_0), np.eye(R))
        block3 = np.zeros((a * R * L1, a * R * L1))
        matrQv_0 = la.block_diag(block0, block1, block2, block3)

        # print(matrQv_0.shape)
        #

        # Q_0
        block0 = kron(np.eye(a), np.dot(matrS1_0, beta1))
        block1 = kron(kron(np.eye(a), np.dot(matrS2_0, beta2)), np.eye(L2))
        block2 = kron(kron(np.eye(a), np.dot(matrS2_0, beta2)), np.eye(R))
        block3 = np.zeros((a * R * L1, a * R * L1))
        matrQ_0 = la.block_diag(block0, block1, block2, block3)

        # print(matrQ_0.shape)
        #

        # Q_1
        block00 = kronsum(kronsum(matrD_k[0], matrH0), matrS1)
        block03 = kron(kron(kron(kron(np.eye(W_), matrH1), M1_e), tau), alpha1)
        block10 = kron(kron(kron(np.eye(a), M2_e), beta1), matrA2_0)
        block11 = kronsum(kronsum(kronsum(matrD_k[0], matrH0), matrS2), matrA2)
        block12 = kron(kron(kron(kron(np.eye(W_), matrH1), np.eye(M2)), tau), L2_e)
        block21 = kron(kron(kron(np.eye(a), np.eye(M2)), matrT0), alpha1)
        block22 = kronsum(kronsum(kronsum(matrD_k[0], matrH), matrS2), matrT)
        block30 = kron(kron(kron(np.eye(a), beta1), matrT0), L1_e)
        block32 = kron(kron(kron(np.eye(a), beta2), np.eye(R)), matrA1_0)
        block33 = kronsum(kronsum(kronsum(matrD_k[0], matrH), matrT), matrA1)
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

        # print(matrQ_1.shape)
        #

        # Q_k
        matrQ_k = [matrQ_0, matrQ_1]
        for k in range(2, n + 2):
            block0 = kron(matrD_k[k - 1], np.eye(V_ * M1))
            block1 = kron(matrD_k[k - 1], np.eye(V_ * M2 * L2))
            block2 = kron(matrD_k[k - 1], np.eye(V_ * M2 * R))
            block3 = kron(matrD_k[k - 1], np.eye(V_ * R * L1))
            matr_temp = la.block_diag(block0, block1, block2, block3)
            matrQ_k.append(matr_temp)
        #


        # Проверка условия эргодичности
        matr_Q_1_ = copy.deepcopy(matrQ_k[0])
        for l in range(1, n + 2):
            matr_Q_1_ += matrQ_k[l]

        matr_dQ_1_ = copy.deepcopy(matrQ_k[1])
        for l in range(2, n + 2):
            matr_dQ_1_ += l * matrQ_k[l]

        matr_a = copy.deepcopy(np.array(matr_Q_1_))
        for i in range(matr_a.shape[0]):
            matr_a[i][0] = 1
        matr_b = np.zeros((matr_a.shape[0], 1))
        matr_b[0][0] = 1
        matr_a = np.transpose(matr_a)

        vect_y = la.solve(matr_a, matr_b).reshape(-1)

        vect_e = np.array([[1.] for _ in range(matr_dQ_1_.shape[1])])
        ergodicity = np.dot(np.dot(vect_y, matr_dQ_1_), vect_e)[0, 0]
        if ergodicity >= 0:
            print('Условие эргодичности не выполняется на итерации', i_experimental, file=sys.stderr)
            print(ergodicity, '>= 0', file=sys.stderr)
            exit()
        else:
            print('iteration:', i_experimental)
            print(ergodicity, '< 0')
        #

        # ---------- Алгоритм вычисления стационарного распределения ---------------

        # 1. Вычисление матрицы G

        matrQ_1_neg_inv = np.linalg.inv(-matrQ_k[1])

        matrG_old = np.eye(matrQ_k[1].shape[0])
        matrG = calc_G(matrG_old)

        i = 1
        while la.norm(matrG - matrG_old, ord=np.inf) >= eps_G:
            matrG_old = matrG
            matrG = calc_G(matrG_old)
            i += 1
        # print(i)
        # print(la.norm(matrG, ord=np.inf))
        # print(matrG.shape)
        #

        # 2. Вычисление G_0
        temp_sum = np.array(copy.deepcopy(matrQ_k[1]))
        for k in range(2, n + 2):
            temp_sum += np.dot(matrQ_k[k], np.linalg.matrix_power(matrG, k - 1))
        matrG_0 = la.inv(temp_sum)
        matrG_0 = -np.dot(matrG_0, matrQv_0)
        #

        # 3. Вычисление матриц Q_il
        matrQ_il = []
        for i in range(0, p_num):
            matrQ_il.append([])
            if i == 0:
                for l in range(0, n + 1):
                    # здесь до n, т.к. нет больше матриц Q_k
                    temp_matr = np.array(copy.deepcopy(matrQw_k[l]))
                    for k in range(l + 1, n + 1):
                        mult_matr = np.array(copy.deepcopy(matrQw_k[k]))
                        for kk in range(k - 1, l - 1, -1):
                            if kk == 0:
                                mult_matr = np.dot(mult_matr, matrG_0)
                            else:
                                mult_matr = np.dot(mult_matr, matrG)
                        temp_matr += mult_matr
                    matrQ_il[i].append(temp_matr)
                for l in range(n + 1, p_num):
                    matrQ_il[i].append(np.zeros(matrQw_k[1].shape))
            else:
                for l in range(0, p_num):
                    if l >= i and (l - i) <= (n + 1):
                        if (l - i + 1) <= (n + 1):
                            temp_matr = np.array(copy.deepcopy(matrQ_k[l - i + 1]))
                        else:
                            temp_matr = np.zeros(matrQ_k[0].shape)

                        for k in range(l + 1, p_num):  # sum from l+1 to inf
                            if (k - i + 1) <= (n + 1):
                                mult_matr = np.array(copy.deepcopy(matrQ_k[k - i + 1]))
                                for kk in range(l, k):
                                    mult_matr = np.dot(mult_matr, matrG)

                                temp_matr += mult_matr
                        matrQ_il[i].append(temp_matr)
                    else:
                        matrQ_il[i].append(np.zeros(matrQ_k[0].shape))
        #

        # 4. Вычисление матриц Phi_l
        matrPhi_0 = np.eye(matrQ_il[0][0].shape[0])
        matrPhi_l = [matrPhi_0]
        for l in range(1, p_num):
            temp_matr = np.dot(np.dot(matrPhi_l[0], matrQ_il[0][l]), la.inv(-matrQ_il[l][l]))
            for i in range(1, l):
                # print(matrPhi_l[i].dot(matrQ_il[i][l]).dot(la.inv(-matrQ_il[l][l])).shape)
                temp_matr += np.dot(np.dot(matrPhi_l[i], matrQ_il[i][l]), la.inv(-matrQ_il[l][l]))
            matrPhi_l.append(temp_matr)
        #

        # 5. Вычисление p_0
        matr_a = np.array(- copy.deepcopy(matrQ_il[0][0]))
        vect_eaR = np.array([[1.] for _ in range(matrPhi_l[0].shape[1])])
        for i in range(1, p_num):
            vect_e = np.array([[1.] for _ in range(matrPhi_l[i].shape[1])])
            vect_eaR += np.dot(matrPhi_l[i], vect_e)

        for i in range(matr_a.shape[0]):
            matr_a[i][0] = vect_eaR[i][0]

        matr_b = np.zeros((matr_a.shape[0], 1))
        matr_b[0][0] = 1.
        matr_a = np.transpose(matr_a)
        p0 = np.transpose(la.solve(matr_a, matr_b))
        #

        # Вычисление векторов cтационарных вероятностей p_l
        vect_p_l = [p0]
        p_sums = [np.sum(p0)]
        # print('p0 = ', vect_p_l[0][0])
        # print('sum0 = ', p_sums[0])
        for l in range(1, p_num):
            vect_p_l.append(np.dot(vect_p_l[0], matrPhi_l[l]))
            p_sums.append(np.sum(vect_p_l[l]))
            # print('p' + str(l) + ' = ', vect_p_l[l][0])
            # print('sum' + str(l) + ' = ', p_sums[l])
        # print('sum = ', np.sum(p_sums))
        #

        # Построение графика вероятностей
        # plt.scatter(x=[_ for _ in range(150)], y=p_sums)
        # plt.show()
        # print(p_sums)
        #

        # Коэффициент загруженности системы
        vect_e = [[1] for _ in range(matrQ_0.shape[1])]
        denom = np.dot(np.dot(vect_y, matrQ_0), vect_e)[0]
        rho = lamD / denom
        # print('\\rho =', rho)
        #

        # Вычисление пропускной способности
        block00 = kronsum(matrH0, matrS1) + kron(np.eye(V_), np.dot(matrS1_0, beta1))
        block03 = kron(kron(kron(matrH1, M1_e), tau), alpha1)
        block10 = kron(kron(np.eye(V_), np.dot(M2_e, beta1)), matrA2_0)
        block11 = kronsum(kronsum(matrH0, matrS2), matrA2) + kron(kron(np.eye(V_), np.dot(matrS2_0, beta2)), np.eye(L2))
        block12 = kron(kron(kron(matrH1, np.eye(M2)), tau), L2_e)
        block21 = kron(kron(kron(np.eye(V_), np.eye(M2)), matrT0), alpha1)
        block22 = kronsum(kronsum(matrH, matrS2), matrT) + kron(kron(np.eye(V_), np.dot(matrS2_0, beta2)), np.eye(R))
        block30 = kron(kron(kron(np.eye(V_), beta1), matrT0), L1_e)
        block32 = kron(kron(kron(np.eye(V_), beta2), np.eye(R)), matrA1_0)
        block33 = kronsum(kronsum(matrH, matrT), matrA1)
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

        matr_a = copy.deepcopy(np.array(matrGamma))
        for i in range(matr_a.shape[0]):
            matr_a[i][0] = 1
        matr_b = np.zeros((matr_a.shape[0], 1))
        matr_b[0][0] = 1
        matr_a = np.transpose(matr_a)

        x = la.solve(matr_a, matr_b).reshape(-1)

        # print('x = ', x)
        x1 = x[0:V_ * M1]
        x2 = x[V_ * M1:V_ * M1 + V_ * M2 * L2]
        x3 = x[V_ * M1 + V_ * M2 * L2: V_ * M1 + V_ * M2 * L2 + V_ * M2 * R]

        # print('x1 = ', x1)
        # print('x2 = ', x2)
        # print('x3 = ', x3)

        e_V_ = np.array([[1.] for i in range(0, V_)])
        e_R = np.array([[1.] for i in range(0, R)])
        pi1 = x1.dot(kron(e_V_, np.eye(M1)))
        pi2 = x2.dot(kron(kron(e_V_, np.eye(M2)), L2_e))
        pi3 = x3.dot(kron(kron(e_V_, np.eye(M2)), e_R))
        # print('pi1 = ', pi1)
        # print('pi2 = ', pi2)
        # print('pi3 = ', pi3)
        #

        # Пропускная способность системы
        varrho = np.dot(pi1, matrS1_0) + np.dot((pi2 + pi3), matrS2_0)
        varrho = varrho[0]
        # print('\\varrho = ', varrho)
        #

        # Векторная производящая функция
        # P(1)
        vect_P_1_ = copy.deepcopy(vect_p_l[1])
        for l in range(2, p_num):
            vect_P_1_ += vect_p_l[l]
        # print('P(1) =', vect_P_1_[0])

        # P'(1)
        vect_dP_1_ = copy.deepcopy(vect_p_l[1])
        for l in range(2, p_num):
            vect_dP_1_ += l * vect_p_l[l]

        # P''(1)
        vect_ddP_1_ = copy.deepcopy(vect_p_l[2]) * 2
        for l in range(3, p_num):
            vect_ddP_1_ += l * (l - 1) * vect_p_l[l]
        #

        # Среднее число запросов в системе
        vect_e = np.array([[1.] for _ in range(vect_dP_1_.shape[1])])
        L = np.dot(vect_dP_1_, vect_e)[0, 0]
        # print('L =', L)
        #

        # Дисперсия числа запросов в системе
        dispV = np.dot(vect_ddP_1_, vect_e)[0] + L - L ** 2
        dispV = dispV[0]
        # print('V =', dispV)
        #

        # Вероятность того, что прибор 1 исправен и обслуживает запрос
        temp_matr = np.dot(vect_P_1_, la.block_diag(np.eye(a * M1), np.zeros(
            (a * (M2 * L2 + M2 * R + R * L1), a * (M2 * L2 + M2 * R + R * L1)))))
        vect_e = np.array([[1.] for i in range(0, temp_matr.shape[1])])
        prob1work = np.dot(temp_matr, vect_e)[0, 0]
        # print('P^{(1,0)} =', prob1work)
        #

        # Вероятность того, что прибор-1 в неисправном состоянии и прибор-2 обслуживает запрос
        temp_matr = np.dot(vect_P_1_, la.block_diag(np.zeros((a * M1, a * M1)), np.eye(a * M2 * (L2 + R)),
                                                    np.zeros((a * R * L1, a * R * L1))))
        vect_e = np.array([[1.] for i in range(0, temp_matr.shape[1])])
        prob1notwork = np.dot(temp_matr, vect_e)[0, 0]
        # print('P^{(2,0),(1,2)} =', prob1notwork)
        #

        # Вероятность того, что в системе есть запросы, прибор 1 в неисправном состоянии и идет переключение с этого прибора на прибор 2 (при этом оба  прибора не обслуживают заявки)
        temp_matr = np.dot(vect_P_1_, la.block_diag(np.zeros((a * (M1 + M2 * (L2 + R)), a * (M1 + M2 * (L2 + R)))),
                                                    np.eye(a * R * L1)))
        vect_e = np.array([[1.] for i in range(0, temp_matr.shape[1])])
        prob1notworkswitch2 = np.dot(temp_matr, vect_e)[0, 0]
        # print('P^{(2,1)} =', prob1notworkswitch2)
        #

        # Вероятность того, что в системе есть запросы, прибор 1 в исправном состоянии и идет переключение с  прибора 2 на прибор 1 (при этом прибор 2 продолжает обслуживать запросы)
        temp_matr = np.dot(vect_P_1_, la.block_diag(np.zeros((a * M1, a * M1)), np.eye(a * M2 * L2),
                                                    np.zeros((a * (R * L1 + M2 * R), a * (R * L1 + M2 * R)))))
        vect_e = np.array([[1.] for i in range(0, temp_matr.shape[1])])
        prob1workswitch21 = np.dot(temp_matr, vect_e)[0, 0]
        # print('P^{(2,2)} =', prob1workswitch21)
        #

        # Вероятность того, что прибор 1 доступен (средняя доля времени, в течение которого прибор 1 доступен)
        temp_matr1 = np.dot(vect_P_1_, la.block_diag(np.eye(a * M1), np.zeros(
            (a * (M2 * L2 + M2 * R + R * L1), a * (M2 * L2 + M2 * R + R * L1)))))
        vect_e = np.array([[1.] for i in range(0, temp_matr1.shape[1])])
        temp_matr = np.dot(temp_matr1, vect_e)

        temp_matr2 = np.dot(vect_p_l[0],
                            la.block_diag(np.eye(a), np.zeros((a * (L2 + R + R * L1), a * (L2 + R + R * L1)))))
        vect_e = np.array([[1.] for i in range(0, temp_matr2.shape[1])])
        temp_matr += np.dot(temp_matr2, vect_e)
        prob1avail = temp_matr[0, 0]
        # print('P_1+ =', prob1avail)
        #

        # Вероятность того, что прибор 1 недоступен (средняя доля времени, в течение которого прибор 1 недоступен)
        temp_matr1 = np.dot(vect_P_1_,
                            la.block_diag(np.zeros((a * M1, a * M1)), np.eye(a * (M2 * L2 + M2 * R + R * L1))))
        vect_e = np.array([[1.] for i in range(0, temp_matr1.shape[1])])
        temp_matr = np.dot(temp_matr1, vect_e)

        temp_matr2 = np.dot(vect_p_l[0], la.block_diag(np.zeros((a, a)), np.eye(a * (L2 + R + R * L1))))
        vect_e = np.array([[1.] for i in range(0, temp_matr2.shape[1])])
        temp_matr += np.dot(temp_matr2, vect_e)
        prob1notavail = temp_matr[0, 0]
        # print('P_2 =', prob1notavail)
        #

        # 10.1 Вероятность того, что прибор 1 недоступен, а прибор 2 доступен (средняя доля времени, в течение которого прибор 2 доступен)
        temp_matr1 = np.dot(vect_P_1_, la.block_diag(np.zeros((a * M1, a * M1)), np.eye(a * M2 * (L2 + R)),
                                                     np.zeros((a * R * L1, a * R * L1))))
        vect_e = np.array([[1.] for i in range(0, temp_matr1.shape[1])])
        temp_matr = np.dot(temp_matr1, vect_e)

        temp_matr2 = np.dot(vect_p_l[0],
                            la.block_diag(np.zeros((a, a)), np.eye(a * (L2 + R)), np.zeros((a * R * L1, a * R * L1))))
        vect_e = np.array([[1.] for i in range(0, temp_matr2.shape[1])])
        temp_matr += np.dot(temp_matr2, vect_e)

        prob2_avail = temp_matr[0, 0]
        # print('P_2 =', prob2_avail)
        #

        # 11. Вероятность того, что оба прибора недоступны, т.е. идет переключение с прибора 1 на прибор 2 (средняя доля времени, в течение которого оба прибора недоступны)
        temp_matr1 = np.dot(vect_P_1_, la.block_diag(np.zeros((a * (M1 + M2 * (L2 + R)), a * (M1 + M2 * (L2 + R)))),
                                                     np.eye(a * R * L1)))
        vect_e = np.array([[1.] for i in range(0, temp_matr1.shape[1])])
        temp_matr = np.dot(temp_matr1, vect_e)

        temp_matr2 = np.dot(vect_p_l[0],
                            la.block_diag(np.zeros((a * (1 + L2 + R), a * (1 + L2 + R))), np.eye(a * R * L1)))
        vect_e = np.array([[1.] for i in range(0, temp_matr2.shape[1])])
        temp_matr += np.dot(temp_matr2, vect_e)

        prob_both_not_avail = temp_matr[0, 0]
        # print('P- =', prob_both_not_avail)
        P_neg_list.append(copy.deepcopy(prob_both_not_avail))
        #

        # 12. Среднее число переключений с прибора 1 на прибор 2 в единицу времени
        temp_matr1 = np.dot(vect_P_1_, la.block_diag(kron(kron(np.eye(W_), matrH1), np.eye(M1)), np.zeros(
            (a * (M2 * (L2 + R) + R * L1), a * (M2 * (L2 + R) + R * L1)))))
        vect_e = np.array([[1.] for i in range(0, temp_matr1.shape[1])])
        temp_matr = np.dot(temp_matr1, vect_e)

        temp_matr2 = np.dot(vect_p_l[0], la.block_diag(kron(np.eye(W_), matrH1),
                                                       np.zeros((a * (L2 + R + R * L1), a * (L2 + R + R * L1)))))
        vect_e = np.array([[1.] for i in range(0, temp_matr2.shape[1])])
        temp_matr += np.dot(temp_matr2, vect_e)

        switches12_num = temp_matr[0, 0]
        # print('Khi_1_2 =', switches12_num)
        #

        # 13. Среднее число переключений с прибора 2 на прибор 1 в единицу времени
        temp_matr1 = np.dot(vect_P_1_, la.block_diag(np.zeros((a * (M1 + M2 * L2), a * (M1 + M2 * L2))),
                                                     kron(np.eye(a * M2), matrT),
                                                     np.zeros((a * R * L1, a * R * L1))))
        vect_e = np.array([[1.] for i in range(0, temp_matr1.shape[1])])
        temp_matr = np.dot(temp_matr1, vect_e)
        temp_matr2 = np.dot(vect_p_l[0], la.block_diag(np.zeros((a * (1 + L2), a * (1 + L2))), kron(np.eye(a), matrT),
                                                       np.zeros((a * R * L1, a * R * L1))))
        vect_e = np.array([[1.] for i in range(0, temp_matr2.shape[1])])
        temp_matr += np.dot(temp_matr2, vect_e)
        switches21_num = temp_matr[0, 0]
        # print('Khi_2_1 =', -switches21_num)
        #
        khi1_khi2_list.append(switches12_num - switches21_num)

    kappa_2_inv_list.append(1 / copy.deepcopy(kappa_2))
    kappa_inv_list_list.append(copy.deepcopy(kappa_inv_list))
    P_neg_list_list.append(copy.deepcopy(P_neg_list))
    khi1_khi2_list_list.append(copy.deepcopy(khi1_khi2_list))

with open("experiment1_2_2.txt", mode='w') as file:
    file.write('lambda = ' + str(lamD) + '\n')
    file.write('h = ' + str(h) + '\n')
    file.write('mu_1 = ' + str(mu_1) + '\n')
    file.write('mu_2 = ' + str(mu_2) + '\n')
    file.write('phi = ' + str(phi) + '\n')
    for i in range(1):
        file.write('1 / kappa_2 = ' + str(kappa_2_inv_list[i]) + '\n')
        file.write('1 / kappa_1 = ' + ', '.join(map(str, kappa_inv_list_list[i])) + '\n')
        file.write('P- = ' + ', '.join(map(str, P_neg_list_list[i])) + '\n')
        file.write('khi1 + khi2 = ' + ', '.join(map(str, khi1_khi2_list_list[i])) + '\n')
