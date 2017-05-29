import ast
import sys
import traceback
from itertools import cycle

import matplotlib.pyplot as plt
from beautifultable import BeautifulTable
from tqdm import tqdm

sys.path.append("../")
from src.cold_reserve_qs import *


def build_plot(experiment_result_list, experiment_name, x_label, y_label, leg_label, file_name='experiment_plot', file_type='png', loc=1, display_title=False):
    """
    Builds experiment plots and saves them to file.

    :param experiment_result_list: list containing experiment results in two dimensions
    :param experiment_name: str with experiment name to be plot's title
    :param x_label: str with x-axis label to be displayed on plot
    :param y_label: str with y-axis label to be displayed on plot
    :param file_type: str with file type ('eps', 'png'). Default is 'png'
    :param loc: int with legend location (1 is right-top, then anti-clockwise) 1 by default
    :param display_title: boolean True for displaying plot title on the top. False by default
    :return: None
    """

    plt.clf()

    fig = plt.figure(1)

    if display_title:
        fig.suptitle(experiment_name)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    lines = ["-", "--", "-.", ":"]
    linecycler = cycle(lines)

    for num, experiments in enumerate(experiment_result_list):
        exp_param = experiments[0]
        x_list = [x[0] for x in experiments[1]]
        y_list = [x[1] for x in experiments[1]]
        plt.plot(x_list, y_list, next(linecycler), label=leg_label + ' = ' + str(round(exp_param, 3)))

    plt.legend(loc=loc)
    plot_filename = '../experiment_plots/' + file_name + '.' + file_type
    fig.savefig(filename=plot_filename, dpi=400, format=file_type)


def experiment_1(queueing_system: ColdReserveQueueingSystem, read_file=False):
    """
    Зависимость \bar{v} от \lambda при различных коэффициентах корреляции c_{cor} во входящем потоке

    :param queueing_system: ColdReserveQueueingSystem
    :param read_file: boolean, if True, read data for plotting from qsr file. False is default
    :return:
    """

    linux_check_cpu_temperature()

    print('Experiment 1 launched!')

    experiment_1_result_list = []

    if not read_file:
        try:
            local_queueing_system = copy.deepcopy(queueing_system)

            queries_matrices = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices)

            for cor_coef in range(3):
                linux_check_cpu_temperature()

                queries_matrices_0 = copy.deepcopy(queries_matrices)

                if cor_coef == 1:
                    matrD_0_g = np.array([[-6.3408, 1.87977 * (10 ** (-6))], [1.87977 * (10 ** (-6)), -0.13888]]) / 17.4
                    matrD_g = np.array([[6.3214, 0.01939], [0.10822, 0.03066]]) / 17.4

                    local_queueing_system.set_BMAP_queries_stream(matrD_0_g, matrD_g,
                                                                  q=local_queueing_system.queries_stream.q,
                                                                  n=local_queueing_system.n)
                    queries_matrices_0 = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices)
                elif cor_coef == 2:
                    # MAP queries stream
                    matrD_0_map = np.array([[-0.575]])
                    matrD_1_map = np.array([[0.575]])
                    local_queueing_system.set_MAP_queries_stream(matrD_0_map, matrD_1_map)
                    queries_matrices_0 = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices)

                characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                filename = '../experiment_results/' + 'experiment_1_' + queueing_system.name + '.qsc'
                local_queueing_system.print_characteristics(filename)
                with open(filename, mode="a") as file:
                    for i, vect in enumerate(vect_p_l):
                        print("P_{} = ".format(str(i)), np.sum(vect), file=file)

                    for i, charact in enumerate(characteristics):
                        print("{}. ".format(str(i)), characteristics_loc[i], ':', charact, file=file)
                    print("============== END SYSTEM =================\n", file=file)

                experiment_1_sublist = [local_queueing_system.queries_stream.c_cor, []]

                output_table = BeautifulTable()
                output_table.column_headers = ['\\rho', '\\lambda', '\\bar{v}']

                for queries_coef in tqdm([i / 12 for i in range(1, 25)]):
                    linux_check_cpu_temperature(notify=False)

                    queries_matrices_1 = [matr * queries_coef for matr in queries_matrices_0]

                    if cor_coef != 2:
                        local_queueing_system.queries_stream.set_transition_matrices(queries_matrices_1)
                    else:
                        local_queueing_system.set_MAP_queries_stream(queries_matrices_1[0], queries_matrices_1[1])

                    characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                    experiment_1_sublist[1].append([local_queueing_system.queries_stream.avg_intensity,
                                                     characteristics[13]])
                    output_table.append_row([characteristics[0], local_queueing_system.queries_stream.avg_intensity, characteristics[13]])

                with open(filename, mode="a") as file:
                    print("c_{cor} = ", local_queueing_system.queries_stream.c_cor, file=file)
                    print("", file=file)
                    print(output_table, file=file)
                    print('\n', file=file)

                experiment_1_result_list.append(copy.deepcopy(experiment_1_sublist))

            file_name = 'experiment_1_' + local_queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_1_result_list))
        except ValueError as e:
            print(str(e))
            traceback.print_exc(file=sys.stderr)
            file_name = 'experiment_1_except' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_1_result_list))
        except Exception:
            traceback.print_exc(file=sys.stderr)
            file_name = 'experiment_1_except' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_1_result_list))
    else:
        file_name = 'experiment_1_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_1_result_list = ast.literal_eval(res_line)

    build_plot(experiment_1_result_list,
               r'Зависимость $\bar v$ от $\lambda$ при различных' + '\nкоэффициентах корреляции $c_{cor}$ во входящем потоке',
               r'$\lambda$',
               r'$\bar v$',
               '$c_{cor}$',
               'experiment_1_gserv')


def experiment_1_1(queueing_system: ColdReserveQueueingSystem, read_file=False):
    """
    Зависимость L от \lambda при различных коэффициентах корреляции c_{cor} во входящем потоке

    :param queueing_system: ColdReserveQueueingSystem
    :param read_file: boolean, if True, read data for plotting from qsr file. False is default
    :return:
    """

    linux_check_cpu_temperature()

    print('Experiment 1_1 launched!')

    experiment_1_1_result_list = []

    if not read_file:
        try:
            local_queueing_system = copy.deepcopy(queueing_system)

            queries_matrices = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices)

            for cor_coef in range(3):
                linux_check_cpu_temperature()

                queries_matrices_0 = copy.deepcopy(queries_matrices)

                if cor_coef == 1:
                    matrD_0_g = np.array([[-6.3408, 1.87977 * (10 ** (-6))], [1.87977 * (10 ** (-6)), -0.13888]]) / 17.4
                    matrD_g = np.array([[6.3214, 0.01939], [0.10822, 0.03066]]) / 17.4

                    local_queueing_system.set_BMAP_queries_stream(matrD_0_g, matrD_g,
                                                                  q=local_queueing_system.queries_stream.q,
                                                                  n=local_queueing_system.n)
                    queries_matrices_0 = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices)
                elif cor_coef == 2:
                    # MAP queries stream
                    matrD_0_map = np.array([[-0.575]])
                    matrD_1_map = np.array([[0.575]])
                    local_queueing_system.set_MAP_queries_stream(matrD_0_map, matrD_1_map)
                    queries_matrices_0 = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices)

                characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                filename = '../experiment_results/' + 'experiment_1_1_' + queueing_system.name + '.qsc'
                local_queueing_system.print_characteristics(filename)
                with open(filename, mode="a") as file:
                    for i, vect in enumerate(vect_p_l):
                        print("P_{} = ".format(str(i)), np.sum(vect), file=file)

                    for i, charact in enumerate(characteristics):
                        print("{}. ".format(str(i)), characteristics_loc[i], ':', charact, file=file)
                    print("============== END SYSTEM =================\n", file=file)

                experiment_1_1_sublist = [local_queueing_system.queries_stream.c_cor, []]

                output_table = BeautifulTable()
                output_table.column_headers = ['\\rho', '\\lambda', 'L']

                for queries_coef in tqdm([i / 12 for i in range(1, 25)]):
                    linux_check_cpu_temperature(notify=False)

                    queries_matrices_1 = [matr * queries_coef for matr in queries_matrices_0]

                    if cor_coef != 2:
                        local_queueing_system.queries_stream.set_transition_matrices(queries_matrices_1)
                    else:
                        local_queueing_system.set_MAP_queries_stream(queries_matrices_1[0], queries_matrices_1[1])

                    characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                    experiment_1_1_sublist[1].append([local_queueing_system.queries_stream.avg_intensity,
                                                     characteristics[2]])
                    output_table.append_row([characteristics[0], local_queueing_system.queries_stream.avg_intensity, characteristics[2]])

                with open(filename, mode="a") as file:
                    print("c_{cor} = ", local_queueing_system.queries_stream.c_cor, file=file)
                    print("", file=file)
                    print(output_table, file=file)
                    print('\n', file=file)

                experiment_1_1_result_list.append(copy.deepcopy(experiment_1_1_sublist))

            file_name = 'experiment_1_1_' + local_queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_1_1_result_list))
        except ValueError as e:
            print(str(e))
            traceback.print_exc(file=sys.stderr)
            file_name = 'experiment_1_1_except_' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_1_1_result_list))
        except Exception:
            traceback.print_exc(file=sys.stderr)
            file_name = 'experiment_1_1_except_' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_1_1_result_list))
    else:
        file_name = 'experiment_1_1_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_1_1_result_list = ast.literal_eval(res_line)

    build_plot(experiment_1_1_result_list,
               r'Зависимость $L$ от $\lambda$ при различных' + '\nкоэффициентах корреляции $c_{cor}$ во входящем потоке',
               r'$\lambda$',
               r'$L$',
               '$c_{cor}$',
               'experiment_1_1')


def experiment_2(queueing_system: ColdReserveQueueingSystem, read_file=False):
    """
    Зависимость P_1^{+} от h при различных коэффициентах корреляции c_{cor} в потоке поломок

    :param queueing_system: ColdReserveQueueingSystem
    :param read_file: boolean, if True, read data for plotting from qsr file. False is default
    :return:
    """

    linux_check_cpu_temperature()

    print('Experiment 2_1 launched!')

    experiment_2_1_result_list = []
    experiment_2_2_result_list = []
    experiment_2_3_result_list = []

    if not read_file:
        try:
            local_queueing_system = copy.deepcopy(queueing_system)

            for cor_coef in range(3):
                linux_check_cpu_temperature()

                break_matrices_0 = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)

                if cor_coef == 1:
                    matrH_0 = np.array([[-6.3408, 1.87977 * (10 ** (-6))], [1.87977 * (10 ** (-6)), -0.13888]])
                    matrH_1 = np.array([[6.3214, 0.01939], [0.10822, 0.03066]])

                    local_queueing_system.set_MAP_break_stream(matrH_0, matrH_1)
                    break_matrices_0 = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)
                elif cor_coef == 2:
                    # MAP queries stream
                    matrH_0 = np.array([[-4.5]])
                    matrH_1 = np.array([[4.5]])
                    local_queueing_system.set_MAP_break_stream(matrH_0, matrH_1)
                    break_matrices_0 = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)

                characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                filename = '../experiment_results/' + 'experiment_2' + queueing_system.name + '.qsc'
                local_queueing_system.print_characteristics(filename)
                with open(filename, mode="a") as file:
                    for i, vect in enumerate(vect_p_l):
                        print("P_{} = ".format(str(i)), np.sum(vect), file=file)

                    for i, charact in enumerate(characteristics):
                        print("{}. ".format(str(i)), characteristics_loc[i], ':', charact, file=file)
                    print("============== END SYSTEM =================\n", file=file)

                experiment_2_1_sublist = [local_queueing_system.break_stream.c_cor, []]
                experiment_2_2_sublist = [local_queueing_system.break_stream.c_cor, []]
                experiment_2_3_sublist = [local_queueing_system.break_stream.c_cor, []]

                output_table = BeautifulTable()
                output_table.column_headers = ['\\rho', 'h', 'P_1^{+}', 'P_2', 'P^{-}']

                for queries_coef in tqdm([i / 15 for i in range(1, 37)]):
                    linux_check_cpu_temperature(notify=False)

                    break_matrices_1 = [matr * queries_coef for matr in break_matrices_0]

                    local_queueing_system.set_MAP_break_stream(break_matrices_1[0], break_matrices_1[1])

                    characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                    experiment_2_1_sublist[1].append([local_queueing_system.break_stream.avg_intensity,
                                                     characteristics[4]])
                    experiment_2_2_sublist[1].append([local_queueing_system.break_stream.avg_intensity,
                                                      characteristics[9]])
                    experiment_2_3_sublist[1].append([local_queueing_system.break_stream.avg_intensity,
                                                      characteristics[10]])

                    output_table.append_row([characteristics[0],
                                             local_queueing_system.break_stream.avg_intensity,
                                             characteristics[4],
                                             characteristics[9],
                                             characteristics[10]])

                with open(filename, mode="a") as file:
                    print("c_{cor} = ", local_queueing_system.break_stream.c_cor, file=file)
                    print("", file=file)
                    print(output_table, file=file)
                    print('\n', file=file)

                experiment_2_1_result_list.append(copy.deepcopy(experiment_2_1_sublist))
                experiment_2_2_result_list.append(copy.deepcopy(experiment_2_2_sublist))
                experiment_2_3_result_list.append(copy.deepcopy(experiment_2_3_sublist))

            file_name = 'experiment_2_1_' + local_queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_2_1_result_list))

            file_name = 'experiment_2_2_' + local_queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_2_2_result_list))

            file_name = 'experiment_2_3_' + local_queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_2_3_result_list))

        except ValueError as e:
            print(str(e))
            traceback.print_exc(file=sys.stderr)
            file_name = 'experiment_2_1_except_' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_2_1_result_list))

            file_name = 'experiment_2_2_except_' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_2_2_result_list))

            file_name = 'experiment_2_3_except_' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_2_3_result_list))
        except Exception:
            traceback.print_exc(file=sys.stderr)
            file_name = 'experiment_2_1_except_' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_2_1_result_list))

            file_name = 'experiment_2_2_except_' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_2_2_result_list))

            file_name = 'experiment_2_3_except_' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_2_3_result_list))
    else:
        file_name = 'experiment_2_1_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_2_1_result_list = ast.literal_eval(res_line)

        file_name = 'experiment_2_2_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_2_2_result_list = ast.literal_eval(res_line)

        file_name = 'experiment_2_3_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_2_3_result_list = ast.literal_eval(res_line)

    build_plot(experiment_2_1_result_list,
               r'Зависимость $P_1^{+}$ от $h$ при различных' + '\nкоэффициентах корреляции $c_{cor}$ в потоке поломок',
               r'h',
               r'$P_1^{+}$',
               '$c_{cor}$',
               'experiment_2_1',
               loc=1)
    build_plot(experiment_2_2_result_list,
               r'Зависимость $P_2$ от $h$ при различных' + '\nкоэффициентах корреляции $c_{cor}$ в потоке поломок',
               r'h',
               r'$P_2$',
               '$c_{cor}$',
               'experiment_2_2',
               loc=4)
    build_plot(experiment_2_3_result_list,
               r'Зависимость $P^{-}$ от $h$ при различных' + '\nкоэффициентах корреляции $c_{cor}$ в потоке поломок',
               r'h',
               r'$P^{-}$',
               '$c_{cor}$',
               'experiment_2_3',
               loc=4)


def experiment_3(queueing_system: ColdReserveQueueingSystem, read_file=False):
    """
    Зависимость \bar{v} от h при различных коэффициентах корреляции c_{cor} в потоке ремонта

    :param queueing_system: ColdReserveQueueingSystem
    :param read_file: boolean, if True, read data for plotting from qsr file. False is default
    :return:
    """

    linux_check_cpu_temperature()

    print('Experiment 3 launched!')

    experiment_3_result_list = []

    if not read_file:
        try:
            local_queueing_system = copy.deepcopy(queueing_system)

            break_matrices = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)

            for cor_coef in range(3):
                linux_check_cpu_temperature()

                # repair_vect = copy.deepcopy(local_queueing_system.repair_stream.repres_vect)
                # repair_matr = copy.deepcopy(local_queueing_system.repair_stream.repres_matr)

                break_matrices_0 = copy.deepcopy(break_matrices)

                local_queueing_system.set_MAP_break_stream(break_matrices_0[0], break_matrices_0[1])

                break_coefficients = [i / 33 for i in range(1, 67, 2)] + [2]

                if cor_coef == 1:
                    repair_vect = np.array([[1., 0.]])
                    repair_matr = np.array([[-1., 1.], [0., -1.]]) / 5

                    local_queueing_system.set_PH_repair_stream(repair_vect, repair_matr)

                    break_coefficients = [i / 15 for i in range(1, 130, 4)]

                elif cor_coef == 2:
                    repair_vect = np.array([[1.]])
                    repair_matr = np.array([[-0.1]])

                    local_queueing_system.set_PH_repair_stream(repair_vect, repair_matr)

                    break_coefficients = [i / 15 for i in range(1, 130, 4)]

                characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                filename = '../experiment_results/' + 'experiment_3_' + queueing_system.name + '.qsc'
                local_queueing_system.print_characteristics(filename)
                with open(filename, mode="a") as file:
                    for i, vect in enumerate(vect_p_l):
                        print("P_{} = ".format(str(i)), np.sum(vect), file=file)

                    for i, charact in enumerate(characteristics):
                        print("{}. ".format(str(i)), characteristics_loc[i], ':', charact, file=file)
                    print("============== END SYSTEM =================\n", file=file)

                experiment_3_sublist = [local_queueing_system.repair_stream.c_var, []]

                output_table = BeautifulTable()
                output_table.column_headers = ['\\rho', 'h', '\\bar{v}']

                for break_coef in tqdm(break_coefficients):
                    linux_check_cpu_temperature(notify=False)

                    break_matrices_1 = [matr * break_coef for matr in break_matrices_0]

                    local_queueing_system.set_MAP_break_stream(break_matrices_1[0], break_matrices_1[1])

                    characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                    experiment_3_sublist[1].append([local_queueing_system.break_stream.avg_intensity,
                                                     characteristics[13]])
                    output_table.append_row([characteristics[0], local_queueing_system.break_stream.avg_intensity, characteristics[13]])

                with open(filename, mode="a") as file:
                    print("c_{var} = ", local_queueing_system.repair_stream.c_var, file=file)
                    print("", file=file)
                    print(output_table, file=file)
                    print('\n', file=file)

                experiment_3_result_list.append(copy.deepcopy(experiment_3_sublist))

            file_name = 'experiment_3_' + local_queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_3_result_list))
        except ValueError as e:
            print(str(e))
            traceback.print_exc(file=sys.stderr)
            file_name = 'experiment_3_except_' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_3_result_list))
        except Exception:
            traceback.print_exc(file=sys.stderr)
            file_name = 'experiment_3_except_' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_3_result_list))
    else:
        file_name = 'experiment_3_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_3_result_list = ast.literal_eval(res_line)

    build_plot(experiment_3_result_list,
               r'Зависимость $\bar{v}$ от $h$ при различных' + '\nкоэффициентах вариации $c_{var}$ в потоке восстановления',
               r'$h$',
               r'$\bar{v}$',
               '$c_{var}$',
               'experiment_3',
               loc=4)

def experiment_4(queueing_system: ColdReserveQueueingSystem, read_file=False):
    """
    Зависимость \bar{v} от h при различных коэффициентах корреляции c_{cor} в потоке поломок

    :param queueing_system: ColdReserveQueueingSystem
    :param read_file: boolean, if True, read data for plotting from qsr file. False is default
    :return:
    """

    linux_check_cpu_temperature()

    print('Experiment 4 launched!')

    experiment_4_result_list = []

    if not read_file:
        try:
            local_queueing_system = copy.deepcopy(queueing_system)

            for cor_coef in range(3):
                linux_check_cpu_temperature()

                break_matrices_0 = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)

                if cor_coef == 1:
                    matrH_0 = np.array([[-6.3408, 1.87977 * (10 ** (-6))], [1.87977 * (10 ** (-6)), -0.13888]])
                    matrH_1 = np.array([[6.3214, 0.01939], [0.10822, 0.03066]])

                    local_queueing_system.set_MAP_break_stream(matrH_0, matrH_1)
                    break_matrices_0 = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)
                elif cor_coef == 2:
                    matrH_0 = np.array([[-4.5]])
                    matrH_1 = np.array([[4.5]])
                    local_queueing_system.set_MAP_break_stream(matrH_0, matrH_1)
                    break_matrices_0 = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)

                characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                filename = '../experiment_results/' + 'experiment_4' + queueing_system.name + '.qsc'
                local_queueing_system.print_characteristics(filename)
                with open(filename, mode="a") as file:
                    for i, vect in enumerate(vect_p_l):
                        print("P_{} = ".format(str(i)), np.sum(vect), file=file)

                    for i, charact in enumerate(characteristics):
                        print("{}. ".format(str(i)), characteristics_loc[i], ':', charact, file=file)
                    print("============== END SYSTEM =================\n", file=file)

                experiment_4_sublist = [local_queueing_system.break_stream.c_cor, []]

                output_table = BeautifulTable()
                output_table.column_headers = ['\\rho', 'h', '\\bar{v}']

                for queries_coef in tqdm([i / 15 for i in range(1, 37)]):
                    linux_check_cpu_temperature(notify=False)

                    break_matrices_1 = [matr * queries_coef for matr in break_matrices_0]

                    local_queueing_system.set_MAP_break_stream(break_matrices_1[0], break_matrices_1[1])

                    characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                    experiment_4_sublist[1].append([local_queueing_system.break_stream.avg_intensity,
                                                    characteristics[13]])

                    output_table.append_row([characteristics[0],
                                             local_queueing_system.break_stream.avg_intensity,
                                             characteristics[13]])

                with open(filename, mode="a") as file:
                    print("c_{cor} = ", local_queueing_system.break_stream.c_cor, file=file)
                    print("", file=file)
                    print(output_table, file=file)
                    print('\n', file=file)

                experiment_4_result_list.append(copy.deepcopy(experiment_4_sublist))

            file_name = 'experiment_4_' + local_queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_4_result_list))

        except ValueError as e:
            print(str(e))
            traceback.print_exc(file=sys.stderr)
            file_name = 'experiment_4_except_' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_4_result_list))

        except Exception:
            traceback.print_exc(file=sys.stderr)
            file_name = 'experiment_4_except_' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_4_result_list))

    else:
        file_name = 'experiment_4_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_4_result_list = ast.literal_eval(res_line)

    build_plot(experiment_4_result_list,
               r'Зависимость $\bar{v}$ от $h$ при различных' + '\nкоэффициентах корреляции $c_{cor}$ в потоке поломок',
               r'h',
               r'$\bar{v}$',
               '$c_{cor}$',
               'experiment_4',
               loc=4)


def experiment_5(queueing_system: ColdReserveQueueingSystem, read_file=False):
    """
    Зависимость \bar{v} от \lambda при различных интенсивностях h потока поломок

    :param queueing_system: ColdReserveQueueingSystem
    :param read_file: boolean, if True, read data for plotting from qsr file. False is default
    :return:
    """

    linux_check_cpu_temperature()

    print('Experiment 5 launched!')

    experiment_5_result_list = []

    if not read_file:
        try:
            local_queueing_system = copy.deepcopy(queueing_system)

            break_matrices = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)
            queries_matrices = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices)

            for break_coef in [1/4, 1, 2.4]:
                linux_check_cpu_temperature()

                break_matrices_0 = copy.deepcopy(break_matrices)

                matrH_0 = break_matrices_0[0] * break_coef
                matrH_1 = break_matrices_0[1] * break_coef

                local_queueing_system.set_MAP_break_stream(matrH_0, matrH_1)

                characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                filename = '../experiment_results/' + 'experiment_5' + queueing_system.name + '.qsc'
                local_queueing_system.print_characteristics(filename)
                with open(filename, mode="a") as file:
                    for i, vect in enumerate(vect_p_l):
                        print("P_{} = ".format(str(i)), np.sum(vect), file=file)

                    for i, charact in enumerate(characteristics):
                        print("{}. ".format(str(i)), characteristics_loc[i], ':', charact, file=file)
                    print("============== END SYSTEM =================\n", file=file)

                experiment_5_sublist = [local_queueing_system.break_stream.avg_intensity, []]

                output_table = BeautifulTable()
                output_table.column_headers = ['\\rho', '\\lambda', '\\bar{v}']

                queries_matrices_0 = copy.deepcopy(queries_matrices)

                for queries_coef in tqdm([i / 12 for i in range(1, 25)]):
                    linux_check_cpu_temperature(notify=False)

                    queries_matrices_1 = [matr * queries_coef for matr in queries_matrices_0]

                    local_queueing_system.queries_stream.set_transition_matrices(queries_matrices_1)

                    characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                    experiment_5_sublist[1].append([local_queueing_system.queries_stream.avg_intensity,
                                                    characteristics[13]])

                    output_table.append_row([characteristics[0],
                                             local_queueing_system.queries_stream.avg_intensity,
                                             characteristics[13]])

                with open(filename, mode="a") as file:
                    print("h = ", local_queueing_system.break_stream.avg_intensity, file=file)
                    print("", file=file)
                    print(output_table, file=file)
                    print('\n', file=file)

                experiment_5_result_list.append(copy.deepcopy(experiment_5_sublist))

            file_name = 'experiment_5_' + local_queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_5_result_list))

        except ValueError as e:
            print(str(e))
            traceback.print_exc(file=sys.stderr)
            file_name = 'experiment_5_except_' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_5_result_list))

        except Exception:
            traceback.print_exc(file=sys.stderr)
            file_name = 'experiment_5_except_' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_5_result_list))

    else:
        file_name = 'experiment_5_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_5_result_list = ast.literal_eval(res_line)

    build_plot(experiment_5_result_list,
               r'Зависимость $\bar{v}$ от $\lambda$ при различных' + '\nинтенсивностях $h$ потока поломок',
               r'$\lambda$',
               r'$\bar{v}$',
               '$h$',
               'experiment_5',
               loc=2,
               display_title=True)

def experiment_6(queueing_system: ColdReserveQueueingSystem, read_file=False):
    """
    Зависимость \chi_{1,2} от среднего времени переключения с прибора-1 на прибор-2 a_1 при различных интенсивностях h потока поломок

    :param queueing_system: ColdReserveQueueingSystem
    :param read_file: boolean, if True, read data for plotting from qsr file. False is default
    :return:
    """

    linux_check_cpu_temperature()

    print('Experiment 6 launched!')

    experiment_6_result_list = []

    if not read_file:
        try:
            local_queueing_system = copy.deepcopy(queueing_system)

            break_matrices = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)
            switch1_2_matr = copy.deepcopy(local_queueing_system.switch1_2_stream.repres_matr)
            switch1_2_vect = copy.deepcopy(local_queueing_system.switch1_2_stream.repres_vect)

            for iter, break_coef in enumerate([1/4, 1, 2.4]):
                linux_check_cpu_temperature()

                break_matrices_0 = copy.deepcopy(break_matrices)

                matrH_0 = break_matrices_0[0] * break_coef
                matrH_1 = break_matrices_0[1] * break_coef

                local_queueing_system.set_MAP_break_stream(matrH_0, matrH_1)

                characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                filename = '../experiment_results/' + 'experiment_6' + queueing_system.name + '.qsc'
                local_queueing_system.print_characteristics(filename)
                with open(filename, mode="a") as file:
                    for i, vect in enumerate(vect_p_l):
                        print("P_{} = ".format(str(i)), np.sum(vect), file=file)

                    for i, charact in enumerate(characteristics):
                        print("{}. ".format(str(i)), characteristics_loc[i], ':', charact, file=file)
                    print("============== END SYSTEM =================\n", file=file)

                experiment_6_sublist = [local_queueing_system.break_stream.avg_intensity, []]

                output_table = BeautifulTable()
                output_table.column_headers = ['\\rho', 'a_1', '\\chi_{1,2}']

                switch1_2_matr_0 = copy.deepcopy(switch1_2_matr)

                s = 1.13547

                for queries_coef in tqdm([0.000013 * (s ** i) for i in range(0, 50)]):
                    linux_check_cpu_temperature(notify=False)

                    switch1_2_matr_1 = switch1_2_matr_0 * queries_coef

                    local_queueing_system.set_PH_switch1_2_stream(switch1_2_vect, switch1_2_matr_1)

                    characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                    experiment_6_sublist[1].append([1 / local_queueing_system.switch1_2_stream.avg_intensity,
                                                    characteristics[11]])

                    output_table.append_row([characteristics[0],
                                             1 / local_queueing_system.switch1_2_stream.avg_intensity,
                                             characteristics[11]])

                with open(filename, mode="a") as file:
                    print("h = ", local_queueing_system.break_stream.avg_intensity, file=file)
                    print("", file=file)
                    print(output_table, file=file)
                    print('\n', file=file)

                experiment_6_result_list.append(copy.deepcopy(experiment_6_sublist))

            file_name = 'experiment_6_' + local_queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_6_result_list))

        except ValueError as e:
            print(str(e))
            traceback.print_exc(file=sys.stderr)
            file_name = 'experiment_6_except_' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_6_result_list))

        except Exception:
            traceback.print_exc(file=sys.stderr)
            file_name = 'experiment_6_except_' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_6_result_list))

    else:
        file_name = 'experiment_6_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_6_result_list = ast.literal_eval(res_line)

    build_plot(experiment_6_result_list,
               r'Зависимость $\chi_{1,2}$ от $a_1$ при различных' + '\nинтенсивностях $h$ потока поломок',
               r'$a_1$',
               r'$\chi_{1,2}$',
               '$h$',
               'experiment_6',
               loc=2)


def experiment_7(queueing_system: ColdReserveQueueingSystem, read_file=False):
    """
    Зависимость \chi_{2,1} от среднего времени переключения с прибора-1 на прибор-2 a_2 при различных интенсивностях h потока поломок

    :param queueing_system: ColdReserveQueueingSystem
    :param read_file: boolean, if True, read data for plotting from qsr file. False is default
    :return:
    """

    linux_check_cpu_temperature()

    print('Experiment 7 launched!')

    experiment_7_result_list = []

    if not read_file:
        try:
            local_queueing_system = copy.deepcopy(queueing_system)

            break_matrices = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)
            switch2_1_matr = copy.deepcopy(local_queueing_system.switch2_1_stream.repres_matr)
            switch2_1_vect = copy.deepcopy(local_queueing_system.switch2_1_stream.repres_vect)

            for iter, break_coef in enumerate([1/4, 1, 2.4]):
                linux_check_cpu_temperature()

                break_matrices_0 = copy.deepcopy(break_matrices)

                matrH_0 = break_matrices_0[0] * break_coef
                matrH_1 = break_matrices_0[1] * break_coef

                local_queueing_system.set_MAP_break_stream(matrH_0, matrH_1)

                characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                filename = '../experiment_results/' + 'experiment_7' + queueing_system.name + '.qsc'
                local_queueing_system.print_characteristics(filename)
                with open(filename, mode="a") as file:
                    for i, vect in enumerate(vect_p_l):
                        print("P_{} = ".format(str(i)), np.sum(vect), file=file)

                    for i, charact in enumerate(characteristics):
                        print("{}. ".format(str(i)), characteristics_loc[i], ':', charact, file=file)
                    print("============== END SYSTEM =================\n", file=file)

                experiment_7_sublist = [local_queueing_system.break_stream.avg_intensity, []]

                output_table = BeautifulTable()
                output_table.column_headers = ['\\rho', 'a_2', '\\chi_{2,1}']

                switch2_1_matr_0 = copy.deepcopy(switch2_1_matr)

                s = 1.13547

                for queries_coef in tqdm([0.0000013 * (s ** i) for i in range(0, 50)]):
                    linux_check_cpu_temperature(notify=False)

                    switch2_1_matr_1 = switch2_1_matr_0 * queries_coef

                    local_queueing_system.set_PH_switch2_1_stream(switch2_1_vect, switch2_1_matr_1)

                    characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                    experiment_7_sublist[1].append([1 / local_queueing_system.switch2_1_stream.avg_intensity,
                                                    characteristics[12]])

                    output_table.append_row([characteristics[0],
                                             1 / local_queueing_system.switch2_1_stream.avg_intensity,
                                             characteristics[12]])

                with open(filename, mode="a") as file:
                    print("h = ", local_queueing_system.break_stream.avg_intensity, file=file)
                    print("", file=file)
                    print(output_table, file=file)
                    print('\n', file=file)

                experiment_7_result_list.append(copy.deepcopy(experiment_7_sublist))

            file_name = 'experiment_7_' + local_queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_7_result_list))

        except ValueError as e:
            print(str(e))
            traceback.print_exc(file=sys.stderr)
            file_name = 'experiment_7_except_' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_7_result_list))

        except Exception:
            traceback.print_exc(file=sys.stderr)
            file_name = 'experiment_7_except_' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_7_result_list))

    else:
        file_name = 'experiment_7_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_7_result_list = ast.literal_eval(res_line)

    build_plot(experiment_7_result_list,
               r'Зависимость $\chi_{2,1}$ от $a_2$ при различных' + '\nинтенсивностях $h$ потока поломок',
               r'$a_2$',
               r'$\chi_{2,1}$',
               '$h$',
               'experiment_7',
               loc=4)