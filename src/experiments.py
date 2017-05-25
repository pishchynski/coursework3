import ast
import sys
import traceback
from itertools import cycle

import matplotlib.pyplot as plt
from tqdm import tqdm

sys.path.append("../")
from src.cold_reserve_qs import *


# todo Check ALL paths in Linux!


def build_plot(experiment_result_list, experiment_name, x_label, y_label, leg_label, file_name='experiment_plot', file_type='png'):
    """
    Builds experiment plots and saves them to file.

    :param experiment_result_list: list containing experiment results in two dimensions
    :param experiment_name: str with experiment name to be plot's title
    :param x_label: str with x-axis label to be displayed on plot
    :param y_label: str with y-axis label to be displayed on plot
    :param file_type: str with file type ('eps', 'png'). Default is 'eps'
    :return: None
    """

    plt.clf()

    fig = plt.figure(1)

    fig.suptitle(experiment_name)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    lines = ["-", "--", "-.", ":"]
    linecycler = cycle(lines)

    for num, experiments in enumerate(experiment_result_list):
        exp_param = experiments[0]
        x_list = [x[0] for x in experiments[1]]
        y_list = [x[1] for x in experiments[1]]
        plt.plot(x_list, y_list, next(linecycler), label=leg_label + ' = ' + str(exp_param))

    plt.legend(loc=2)
    plot_filename = '../experiment_plots/' + file_name + '.' + file_type
    fig.savefig(filename=plot_filename, dpi=1000, format=file_type)


def experiment_1(queueing_system: ColdReserveQueueingSystem, read_file=False):
    """
    Зависимость вероятности P- того, что оба прибора недоступны, от среднего времени 1/kappa_1
    переключения с прибора-1 на прибор-2 при разных средних временах 1/kappa_2 переключения с прибора-2 на прибор-1.

    :param queueing_system: ColdReserveQueueingSystem
    :return: None
    """

    linux_check_cpu_temperature()

    print('Experiment 1 launched!')



    experiment_1_result_list = []

    if not read_file:

        try:
            local_queueing_system = copy.deepcopy(queueing_system)  # todo Check how impacts memory consumption!

            switch2_1_matr = copy.deepcopy(local_queueing_system.switch2_1_stream.repres_matr)
            switch2_1_vect = copy.deepcopy(local_queueing_system.switch2_1_stream.repres_vect)

            switch1_2_matr = copy.deepcopy(local_queueing_system.switch1_2_stream.repres_matr)
            switch1_2_vect = copy.deepcopy(local_queueing_system.switch1_2_stream.repres_vect)

            s = 1.0868994010299229698

            for switch2_1_coef in (0.000001, 0.00001, 0.0001):
                linux_check_cpu_temperature()

                switch2_1_matr_1 = copy.deepcopy(switch2_1_matr) * switch2_1_coef
                local_queueing_system.set_PH_switch2_1_stream(switch2_1_vect, switch2_1_matr_1)

                characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                filename = '../experiment_results/' + 'experiment_1_' + queueing_system.name + '.qsc'
                local_queueing_system.print_characteristics(filename)
                with open(filename, mode="a") as file:
                    for i, vect in enumerate(vect_p_l):
                        print("P_{} = ".format(str(i)), np.sum(vect), file=file)

                    for i, charact in enumerate(characteristics):
                        print("{}. ".format(str(i)), characteristics_loc[i], ':', charact, file=file)
                    print("============== END SYSTEM =================", file=file)

                experiment_1_sublist = [1 / local_queueing_system.switch2_1_stream.avg_intensity, []]
                for switch1_2_coef in tqdm([0.000001 * (s ** i) for i in range(0, 200)]):

                    switch1_2_matr_1 = copy.deepcopy(switch1_2_matr) * switch1_2_coef

                    local_queueing_system.set_PH_switch1_2_stream(switch1_2_vect, switch1_2_matr_1)
                    characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                    experiment_1_sublist[1].append(
                        [1 / local_queueing_system.switch1_2_stream.avg_intensity, characteristics[10]])
                experiment_1_result_list.append(copy.deepcopy(experiment_1_sublist))

            file_name = 'experiment_1_' + local_queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_1_result_list))
        except ValueError as e:
            print(str(e))
            file_name = 'experiment_1_except' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_1_result_list))
        except Exception:
            traceback.print_exc(file=sys.stderr)
    else:
        file_name = 'experiment_1_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_1_result_list = ast.literal_eval(res_line)

    build_plot(experiment_1_result_list,
               'P- от 1/kappa_1',
               '1/kappa_1',
               'P-',
               'a_2',
               'experiment_1')


def experiment_2(queueing_system: ColdReserveQueueingSystem, read_file=False):
    """
     Зависимость числа переключений режимов 1_2 + 2_1 от среднего времени 1/kappa_1 переключения с прибора-1 на прибор-2
     при различных средних временах 1/kappa_2 переключений с прибора-2 на прибор-1

    :param queueing_system: ColdReserveQueueingSystem
    :return: None
    """

    linux_check_cpu_temperature()

    print('Experiment 2 launched!')

    experiment_2_result_list = []

    if not read_file:
        try:
            local_queueing_system = copy.deepcopy(queueing_system)

            switch2_1_matr = copy.deepcopy(local_queueing_system.switch2_1_stream.repres_matr)
            switch2_1_vect = copy.deepcopy(local_queueing_system.switch2_1_stream.repres_vect)

            switch1_2_matr = copy.deepcopy(local_queueing_system.switch1_2_stream.repres_matr)
            switch1_2_vect = copy.deepcopy(local_queueing_system.switch1_2_stream.repres_vect)

            for switch2_1_coef in (0.1, 1, 10):
                linux_check_cpu_temperature()

                switch2_1_matr_1 = copy.deepcopy(switch2_1_matr) * switch2_1_coef  # todo Check if copy necessary
                local_queueing_system.set_PH_switch2_1_stream(switch2_1_vect, switch2_1_matr_1)

                filename = '../experiment_results/' + 'experiment_2_' + queueing_system.name + '.qsc'
                local_queueing_system.print_characteristics(filename)

                experiment_2_sublist = [local_queueing_system.switch2_1_stream.avg_intensity, []]
                for switch1_2_coef in tqdm([i / 50 for i in range(1, 1001)]):
                    switch1_2_matr_1 = copy.deepcopy(switch1_2_matr) * switch1_2_coef
                    local_queueing_system.set_PH_switch1_2_stream(switch1_2_vect, switch1_2_matr_1)
                    characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                    experiment_2_sublist[1].append([1 / local_queueing_system.switch1_2_stream.avg_intensity,
                                                   characteristics[14]])
                experiment_2_result_list.append(copy.deepcopy(experiment_2_sublist))

                file_name = 'experiment_2_' + local_queueing_system.name + '.qsr'
                with open('../experiment_results/' + file_name, mode='w') as res_file:
                    res_file.write(str(experiment_2_result_list))
        except ValueError as e:
            print(str(e))
            file_name = 'experiment_2_except' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_2_result_list))
        except Exception:
            traceback.print_exc(file=sys.stderr)
    else:
        file_name = 'experiment_2_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_2_result_list = ast.literal_eval(res_line)

    build_plot(experiment_2_result_list,
               'Число переключений от 1/kappa_1',
               '1/kappa_1',
               'switch_num',
               'experiment_2')


def experiment_3(queueing_system: ColdReserveQueueingSystem, read_file=False):
    """
    Зависимость L от lambda при различных значениях h.

    :param queueing_system: ColdReserveQueueingSystem
    :return: None
    """

    linux_check_cpu_temperature()

    print('Experiment 3 launched!')

    experiment_3_result_list = []

    if not read_file:
        try:
            local_queueing_system = copy.deepcopy(queueing_system)

            queries_matrices = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices)
            break_matrices = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)

            for break_coef in (0.25, 1, 4):
                linux_check_cpu_temperature()

                break_matrices_1 = [matr * break_coef for matr in break_matrices]  # todo Check if copy necessary
                local_queueing_system.set_MAP_break_stream(break_matrices_1[0], break_matrices_1[1])

                filename = '../experiment_results/' + 'experiment_3_' + queueing_system.name + '.qsc'
                local_queueing_system.print_characteristics(filename)

                experiment_3_sublist = [local_queueing_system.break_stream.avg_intensity, []]
                for queries_coef in tqdm([i / 500 for i in range(1, 1001)]):
                    queries_matrices_1 = [matr * queries_coef for matr in queries_matrices]
                    local_queueing_system.queries_stream.set_transition_matrices(queries_matrices_1)
                    characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                    experiment_3_sublist[1].append([local_queueing_system.queries_stream.avg_intensity,
                                                    characteristics[2]])
                experiment_3_result_list.append(copy.deepcopy(experiment_3_sublist))

            file_name = 'experiment_3_' + local_queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_3_result_list))
        except ValueError as e:
            print(str(e))
            file_name = 'experiment_3_except' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_3_result_list))
        except Exception:
            traceback.print_exc(file=sys.stderr)
    else:
        file_name = 'experiment_3_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_3_result_list = ast.literal_eval(res_line)

    build_plot(experiment_3_result_list,
               'L от \\lambda',
               '\\lambda',
               'L',
               'experiment_3')


def experiment_4(queueing_system: ColdReserveQueueingSystem, read_file=False):
    """
    Зависимость L от lambda при различных коэффициентах корреляции c_cor во входном потоке.
    :param queueing_system: ColdReserveQueueingSystem
    :return: None
    """

    linux_check_cpu_temperature()

    print('Experiment 4 launched!')

    experiment_4_result_list = []

    if not read_file:
        try:
            local_queueing_system = copy.deepcopy(queueing_system)

            queries_matrices = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices)
            matrD = copy.deepcopy(local_queueing_system.queries_stream.matrD)

            for cor_coef in (0.01, 1, 100):
                linux_check_cpu_temperature()

                matrD_10 = copy.deepcopy(matrD)
                matrD_10[0] *= cor_coef
                local_queueing_system.set_BMAP_queries_stream(queries_matrices[0], matrD_10, q=local_queueing_system.queries_stream.q, n=local_queueing_system.n)

                filename = '../experiment_results/' + 'experiment_4_' + queueing_system.name + '.qsc'
                local_queueing_system.print_characteristics(filename)

                experiment_4_sublist = [local_queueing_system.queries_stream.c_cor, []]
                for queries_coef in tqdm([i / 500 for i in range(1, 1001)]):
                    matrD_0_1 = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices[0]) * queries_coef
                    matrD_1 = copy.deepcopy(local_queueing_system.queries_stream.matrD) * queries_coef
                    local_queueing_system.set_BMAP_queries_stream(matrD_0_1, matrD_1, q=local_queueing_system.queries_stream.q, n=local_queueing_system.n)
                    characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                    experiment_4_sublist[1].append([local_queueing_system.queries_stream.avg_intensity,
                                                 characteristics[2]])
                experiment_4_result_list.append(copy.deepcopy(experiment_4_sublist))

            file_name = 'experiment_4_' + local_queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_4_result_list))
        except ValueError as e:
            print(str(e))
            file_name = 'experiment_4_except' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_4_result_list))
        except Exception:
            traceback.print_exc(file=sys.stderr)
    else:
        file_name = 'experiment_4_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_4_result_list = ast.literal_eval(res_line)

    build_plot(experiment_4_result_list,
               'L от \\lambda',
               '\\lambda',
               'L',
               'experiment_4')


def experiment_5(queueing_system: ColdReserveQueueingSystem, read_file=False):
    """
    Зависимость V от lambda при различных коэффициентах корреляции c_cor во входном потоке.

    :param queueing_system: ColdReserveQueueingSystem
    :return: None
    """

    linux_check_cpu_temperature()

    print('Experiment 5 launched!')

    experiment_5_result_list = []

    if not read_file:
        try:
            local_queueing_system = copy.deepcopy(queueing_system)

            queries_matrices = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices)
            matrD = copy.deepcopy(local_queueing_system.queries_stream.matrD)

            for cor_coef in (0.01, 1, 100):
                linux_check_cpu_temperature()

                matrD_10 = copy.deepcopy(matrD)
                matrD_10[0] *= cor_coef
                local_queueing_system.set_BMAP_queries_stream(queries_matrices[0], matrD_10, q=local_queueing_system.queries_stream.q, n=local_queueing_system.n)

                filename = '../experiment_results/' + 'experiment_5_' + queueing_system.name + '.qsc'
                local_queueing_system.print_characteristics(filename)

                experiment_5_sublist = [local_queueing_system.queries_stream.c_cor, []]
                for queries_coef in tqdm([i / 100 for i in range(1, 1001)]):
                    matrD_0_1 = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices[0]) * queries_coef
                    matrD_1 = copy.deepcopy(local_queueing_system.queries_stream.matrD) * queries_coef
                    local_queueing_system.set_BMAP_queries_stream(matrD_0_1, matrD_1, q=local_queueing_system.queries_stream.q, n=local_queueing_system.n)
                    characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                    experiment_5_sublist[1].append([local_queueing_system.queries_stream.avg_intensity,
                                                 list(characteristics.items())[3][1]])
                experiment_5_result_list.append(copy.deepcopy(experiment_5_sublist))

            file_name = 'experiment_5_' + local_queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_5_result_list))
        except ValueError as e:
            print(str(e))
            file_name = 'experiment_5_except' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_5_result_list))
        except Exception:
            traceback.print_exc(file=sys.stderr)
    else:
        file_name = 'experiment_5_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_5_result_list = ast.literal_eval(res_line)

    build_plot(experiment_5_result_list,
               'V от \\lambda',
               '\\lambda',
               'V',
               'experiment_5')


def experiment_6(queueing_system: ColdReserveQueueingSystem, read_file=False):
    """
    Зависимость вероятности доступности P_1+ доступности прибора-1 от h при различных коэффициентах корреляции c_cor
    в потоке поломок
    :param queueing_system: ColdReserveQueueingSystem
    :return: None
    """

    linux_check_cpu_temperature()

    print('Experiment 6 launched!')

    experiment_6_result_list = []

    if not read_file:
        try:
            local_queueing_system = copy.deepcopy(queueing_system)

            break_matrices = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)

            for cor_coef in (0.01, 1, 100):
                linux_check_cpu_temperature()

                break_matr0 = copy.deepcopy(break_matrices[0])
                break_matr0[0] *= cor_coef
                break_matr1 = copy.deepcopy(break_matrices[1])
                break_matr1[0] *= cor_coef

                break_matrices_0 = [break_matr0, break_matr1]

                local_queueing_system.set_MAP_break_stream(break_matrices_0[0], break_matrices_0[1])

                filename = '../experiment_results/' + 'experiment_6_' + queueing_system.name + '.qsc'
                local_queueing_system.print_characteristics(filename)

                experiment_6_sublist = [local_queueing_system.break_stream.c_cor, []]
                for break_coef in tqdm([i / 50 for i in range(1, 1001)]):
                    break_matrices_1 = copy.deepcopy(break_matrices_0)
                    break_matrices_1[0] *= break_coef
                    break_matrices_1[1] *= break_coef
                    local_queueing_system.set_MAP_break_stream(break_matrices_1[0], break_matrices_1[1])
                    characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                    experiment_6_sublist[1].append([local_queueing_system.break_stream.avg_intensity,
                                                 list(characteristics.items())[8][1]])
                experiment_6_result_list.append(copy.deepcopy(experiment_6_sublist))

            file_name = 'experiment_6_' + local_queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_6_result_list))
        except ValueError as e:
            print(str(e))
            file_name = 'experiment_6_except' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_6_result_list))
        except Exception:
            traceback.print_exc(file=sys.stderr)
    else:
        file_name = 'experiment_6_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_6_result_list = ast.literal_eval(res_line)

    build_plot(experiment_6_result_list,
               'P_1+ от h',
               'h',
               'P_1+',
               'experiment_6')


def experiment_7(queueing_system: ColdReserveQueueingSystem, read_file=False):
    """
    Зависимость вероятности доступности P_1- недоступности прибора-1 от h при различных коэффициентах корреляции c_cor
    в потоке поломок
    :param queueing_system: ColdReserveQueueingSystem
    :return: None
    """

    linux_check_cpu_temperature()

    print('Experiment 7 launched!')

    experiment_7_result_list = []

    if not read_file:
        try:
            local_queueing_system = copy.deepcopy(queueing_system)

            break_matrices = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)

            for cor_coef in (0.01, 1, 100):
                linux_check_cpu_temperature()

                break_matr0 = copy.deepcopy(break_matrices[0])
                break_matr0[0] *= cor_coef
                break_matr1 = copy.deepcopy(break_matrices[1])
                break_matr1[0] *= cor_coef

                break_matrices_0 = [break_matr0, break_matr1]

                local_queueing_system.set_MAP_break_stream(break_matrices_0[0], break_matrices_0[1])

                filename = '../experiment_results/' + 'experiment_7_' + queueing_system.name + '.qsc'
                local_queueing_system.print_characteristics(filename)

                experiment_7_sublist = [local_queueing_system.break_stream.c_cor, []]
                for break_coef in tqdm([i / 50 for i in range(1, 1001)]):
                    break_matrices_1 = copy.deepcopy(break_matrices_0)
                    break_matrices_1[0] *= break_coef
                    break_matrices_1[1] *= break_coef
                    local_queueing_system.set_MAP_break_stream(break_matrices_1[0], break_matrices_1[1])
                    characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                    experiment_7_sublist[1].append([local_queueing_system.break_stream.avg_intensity,
                                                 1 - list(characteristics.items())[8][1]])
                experiment_7_result_list.append(copy.deepcopy(experiment_7_sublist))

            file_name = 'experiment_7_' + local_queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_7_result_list))
        except ValueError as e:
            print(str(e))
            file_name = 'experiment_7_except' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_7_result_list))
        except Exception:
            traceback.print_exc(file=sys.stderr)
    else:
        file_name = 'experiment_7_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_7_result_list = ast.literal_eval(res_line)

    build_plot(experiment_7_result_list,
               'P_1- от h',
               'h',
               'P_1-',
               'experiment_7')


def experiment_8(queueing_system: ColdReserveQueueingSystem, read_file=False):
    """
    Зависимость L от h при различных коэффициентах вариации c_var времени ремонта

    :param queueing_system: ColdReserveQueueingSystem
    :return: None
    """

    linux_check_cpu_temperature()

    print('Experiment 8 launched!')

    experiment_8_result_list = []

    if not read_file:
        try:
            local_queueing_system = copy.deepcopy(queueing_system)

            break_matrices = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)

            for repair_vect_elem in (0.1, 0.5, 0.9):
                linux_check_cpu_temperature()

                local_queueing_system.set_PH_repair_stream(np.array([[repair_vect_elem, 1 - repair_vect_elem]]),
                                                           local_queueing_system.repair_stream.repres_matr)

                filename = '../experiment_results/' + 'experiment_8_' + queueing_system.name + '.qsc'
                local_queueing_system.print_characteristics(filename)

                experiment_8_sublist = [local_queueing_system.repair_stream.c_var, []]
                for break_coef in tqdm([i / 50 for i in range(1, 1001)]):
                    break_matrices_1 = copy.deepcopy(break_matrices)
                    break_matrices_1[0] *= break_coef
                    break_matrices_1[1] *= break_coef
                    local_queueing_system.set_MAP_break_stream(break_matrices_1[0], break_matrices_1[1])
                    characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                    experiment_8_sublist[1].append([local_queueing_system.break_stream.avg_intensity,
                                                 list(characteristics.items())[2][1]])
                experiment_8_result_list.append(copy.deepcopy(experiment_8_sublist))

            file_name = 'experiment_8_' + local_queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_8_result_list))
        except ValueError as e:
            print(str(e))
            file_name = 'experiment_8_except' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_8_result_list))
        except Exception:
            traceback.print_exc(file=sys.stderr)
    else:
        file_name = 'experiment_8_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_8_result_list = ast.literal_eval(res_line)

    build_plot(experiment_8_result_list,
               'L от h',
               'h',
               'L',
               'experiment_8')


def experiment_9(queueing_system: ColdReserveQueueingSystem, read_file=False):
    """
    Зависимость пропускной способности от lambda при различных значениях h.

    :param queueing_system: ColdReserveQueueingSystem
    :return: None
    """

    linux_check_cpu_temperature()

    print('Experiment 9 launched!')

    experiment_9_result_list = []

    if not read_file:
        try:
            local_queueing_system = copy.deepcopy(queueing_system)

            queries_matrices = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices)
            break_matrices = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)

            for break_coef in (0.1, 1, 10):
                linux_check_cpu_temperature()

                break_matrices_1 = [matr * break_coef for matr in break_matrices]  # todo Check if copy necessary
                local_queueing_system.set_MAP_break_stream(break_matrices_1[0], break_matrices_1[1])

                filename = '../experiment_results/' + 'experiment_9_' + queueing_system.name + '.qsc'
                local_queueing_system.print_characteristics(filename)

                experiment_9_sublist = [local_queueing_system.break_stream.avg_intensity, []]
                for queries_coef in tqdm([i / 50 for i in range(1, 1001)]):
                    queries_matrices_1 = [matr * queries_coef for matr in queries_matrices]
                    local_queueing_system.queries_stream.set_transition_matrices(queries_matrices_1)
                    characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                    experiment_9_sublist[1].append([local_queueing_system.queries_stream.avg_intensity,
                                                 list(characteristics.items())[1][1]])
                experiment_9_result_list.append(copy.deepcopy(experiment_9_sublist))

            file_name = 'experiment_9_' + local_queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_9_result_list))
        except ValueError as e:
            print(str(e))
            file_name = 'experiment_9_except' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_9_result_list))
        except Exception:
            traceback.print_exc(file=sys.stderr)
    else:
        file_name = 'experiment_9_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_9_result_list = ast.literal_eval(res_line)

    build_plot(experiment_9_result_list,
               '\\rho от \\lambda',
               '\\lambda',
               '\\rho',
               'experiment_9')


def experiment_10(queueing_system: ColdReserveQueueingSystem, read_file=False):
    """
    Зависимость u_ от \lambda при различных коэффициентах корреляции во входящем потоке

    :param queueing_system:
    :param read_file:
    :return:
    """

    linux_check_cpu_temperature()

    print('Experiment 10 launched!')

    experiment_10_result_list = []

    if not read_file:
        try:
            local_queueing_system = copy.deepcopy(queueing_system)

            queries_matrices = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices)
            matrD = copy.deepcopy(local_queueing_system.queries_stream.matrD)

            for cor_coef in range(2):
                linux_check_cpu_temperature()

                queries_matrices_0 = copy.deepcopy(queries_matrices)

                if cor_coef == 1:
                    matrD_0_g = np.array([[-6., 0.01], [0.02, -2.76]]) * 5
                    matrD_g = np.array([[5., 0.99], [0.2, 2.54]]) * 5
                
                    local_queueing_system.set_BMAP_queries_stream(matrD_0_g, matrD_g,
                                                                  q=local_queueing_system.queries_stream.q,
                                                                  n=local_queueing_system.n)
                    queries_matrices_0 = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices)
                characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                filename = '../experiment_results/' + 'experiment_10_' + queueing_system.name + '.qsc'
                local_queueing_system.print_characteristics(filename)
                with open(filename, mode="a") as file:
                    for i, vect in enumerate(vect_p_l):
                        print("P_{} = ".format(str(i)), np.sum(vect), file=file)

                    for i, charact in enumerate(characteristics):
                        print("{}. ".format(str(i)), characteristics_loc[i], ':', charact, file=file)
                    print("============== END SYSTEM =================", file=file)
                experiment_10_sublist = [local_queueing_system.queries_stream.c_cor, []]

                for queries_coef in tqdm([i / 25 for i in range(1, 51)]):
                    queries_matrices_1 = [matr * queries_coef for matr in queries_matrices_0]

                    local_queueing_system.queries_stream.set_transition_matrices(queries_matrices_1)
                    
                    characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                    experiment_10_sublist[1].append([local_queueing_system.queries_stream.avg_intensity,
                                                    characteristics[13]])
                experiment_10_result_list.append(copy.deepcopy(experiment_10_sublist))

            file_name = 'experiment_10_' + local_queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_10_result_list))
        except ValueError as e:
            print(str(e))
            file_name = 'experiment_10_except' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_10_result_list))
        except Exception:
            traceback.print_exc(file=sys.stderr)
            file_name = 'experiment_10_except' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_10_result_list))
    else:
        file_name = 'experiment_10_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_10_result_list = ast.literal_eval(res_line)

    build_plot(experiment_10_result_list,
               'u_ от lambda',
               'lambda',
               'u_',
               'c_cor_q',
               'experiment_10')


def experiment_11(queueing_system: ColdReserveQueueingSystem, read_file=False):
    """
    Зависимость u_ от h при различных коэффициентах вариации c_var времени ремонта

    :param queueing_system: ColdReserveQueueingSystem
    :return: None
    """

    linux_check_cpu_temperature()

    print('Experiment 11 launched!')

    experiment_11_result_list = []

    if not read_file:
        try:
            local_queueing_system = copy.deepcopy(queueing_system)

            break_matrices = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)

            for repair_vect_elem in (0.1, 0.5, 0.9):
                linux_check_cpu_temperature()

                local_queueing_system.set_PH_repair_stream(np.array([[repair_vect_elem, 1 - repair_vect_elem]]),
                                                           local_queueing_system.repair_stream.repres_matr)

                filename = '../experiment_results/' + 'experiment_11_' + queueing_system.name + '.qsc'
                local_queueing_system.print_characteristics(filename)

                experiment_11_sublist = [local_queueing_system.repair_stream.c_var, []]
                for break_coef in tqdm([i / 5 for i in range(1, 101)]):
                    break_matrices_1 = copy.deepcopy(break_matrices)
                    break_matrices_1[0] *= break_coef
                    break_matrices_1[1] *= break_coef
                    local_queueing_system.set_MAP_break_stream(break_matrices_1[0], break_matrices_1[1])
                    characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                    experiment_11_sublist[1].append([local_queueing_system.break_stream.avg_intensity,
                                                    list(characteristics.items())[13][1]])
                experiment_11_result_list.append(copy.deepcopy(experiment_11_sublist))

            file_name = 'experiment_11_' + local_queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_11_result_list))
        except ValueError as e:
            print(str(e))
            file_name = 'experiment_11_except' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_11_result_list))
        except Exception:
            traceback.print_exc(file=sys.stderr)
    else:
        file_name = 'experiment_11_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_11_result_list = ast.literal_eval(res_line)

    build_plot(experiment_11_result_list,
               'u_ от h',
               'h',
               'u_',
               'experiment_11')


def experiment_12(queueing_system: ColdReserveQueueingSystem, read_file=False):
    """
    Зависимость u_ от h при различных коэффициентах корреляции c_cor времени ремонта

    :param queueing_system: ColdReserveQueueingSystem
    :return: None
    """

    linux_check_cpu_temperature()

    print('Experiment 12 launched!')

    experiment_12_result_list = []

    if not read_file:
        try:
            local_queueing_system = copy.deepcopy(queueing_system)

            break_matrices = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)

            for repair_vect_elem in (0.1, 0.5, 0.9):
                linux_check_cpu_temperature()

                local_queueing_system.set_PH_repair_stream(np.array([[repair_vect_elem, 1 - repair_vect_elem]]),
                                                           local_queueing_system.repair_stream.repres_matr)

                filename = '../experiment_results/' + 'experiment_12_' + queueing_system.name + '.qsc'
                local_queueing_system.print_characteristics(filename)

                experiment_12_sublist = [local_queueing_system.repair_stream.c_cor, []]
                for break_coef in tqdm([i / 5 for i in range(1, 101)]):
                    break_matrices_1 = copy.deepcopy(break_matrices)
                    break_matrices_1[0] *= break_coef
                    break_matrices_1[1] *= break_coef
                    local_queueing_system.set_MAP_break_stream(break_matrices_1[0], break_matrices_1[1])
                    characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                    experiment_12_sublist[1].append([local_queueing_system.break_stream.avg_intensity,
                                                    list(characteristics.items())[13][1]])
                experiment_12_result_list.append(copy.deepcopy(experiment_12_sublist))

            file_name = 'experiment_12_' + local_queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_12_result_list))
        except ValueError as e:
            print(str(e))
            file_name = 'experiment_12_except' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_12_result_list))
        except Exception:
            traceback.print_exc(file=sys.stderr)
    else:
        file_name = 'experiment_12_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_12_result_list = ast.literal_eval(res_line)

    build_plot(experiment_12_result_list,
               'u_ от h',
               'h',
               'u_',
               'experiment_12')


def experiment_13(queueing_system: ColdReserveQueueingSystem, read_file=False):
    """
    Зависимость u_ от lambda при различных значениях h.

    :param queueing_system: ColdReserveQueueingSystem
    :return: None
    """

    linux_check_cpu_temperature()

    print('Experiment 13 launched!')

    experiment_13_result_list = []

    if not read_file:
        try:
            local_queueing_system = copy.deepcopy(queueing_system)

            queries_matrices = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices)
            break_matrices = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)

            for break_coef in (0.1, 1, 10):
                linux_check_cpu_temperature()

                break_matrices_1 = [matr * break_coef for matr in break_matrices]  # todo Check if copy necessary
                local_queueing_system.set_MAP_break_stream(break_matrices_1[0], break_matrices_1[1])

                filename = '../experiment_results/' + 'experiment_13_' + queueing_system.name + '.qsc'
                local_queueing_system.print_characteristics(filename)

                experiment_13_sublist = [local_queueing_system.break_stream.avg_intensity, []]
                for queries_coef in tqdm([i / 500 for i in range(1, 1001)]):
                    queries_matrices_1 = [matr * queries_coef for matr in queries_matrices]
                    local_queueing_system.queries_stream.set_transition_matrices(queries_matrices_1)
                    characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                    experiment_13_sublist[1].append([local_queueing_system.queries_stream.avg_intensity,
                                                 list(characteristics.items())[2][1]])
                experiment_13_result_list.append(copy.deepcopy(experiment_13_sublist))

            file_name = 'experiment_13_' + local_queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_13_result_list))
        except ValueError as e:
            print(str(e))
            file_name = 'experiment_13_except' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_13_result_list))
        except Exception:
            traceback.print_exc(file=sys.stderr)
    else:
        file_name = 'experiment_13_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_13_result_list = ast.literal_eval(res_line)

    build_plot(experiment_13_result_list,
               'u_ от \\lambda',
               '\\lambda',
               'u_',
               'experiment_13')


def experiment_14(queueing_system: ColdReserveQueueingSystem, read_file=False):
    """
    Зависимость количества переключений с прибора-1 на прибор-2, \khi_{1,2} от 1/kappa_1 при различных значениях h.

    :param queueing_system: ColdReserveQueueingSystem
    :return: None
    """

    linux_check_cpu_temperature()

    print('Experiment 14 launched!')

    experiment_14_result_list = []

    if not read_file:
        try:
            local_queueing_system = copy.deepcopy(queueing_system)

            queries_matrices = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices)
            break_matrices = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)

            for break_coef in (0.1, 1, 10):
                linux_check_cpu_temperature()

                break_matrices_1 = [matr * break_coef for matr in break_matrices]  # todo Check if copy necessary
                local_queueing_system.set_MAP_break_stream(break_matrices_1[0], break_matrices_1[1])

                filename = '../experiment_results/' + 'experiment_14_' + queueing_system.name + '.qsc'
                local_queueing_system.print_characteristics(filename)

                experiment_14_sublist = [local_queueing_system.break_stream.avg_intensity, []]
                for queries_coef in tqdm([i / 50 for i in range(1, 1001)]):
                    queries_matrices_1 = [matr * queries_coef for matr in queries_matrices]
                    local_queueing_system.queries_stream.set_transition_matrices(queries_matrices_1)
                    characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                    experiment_14_sublist[1].append([1 / local_queueing_system.switch1_2_stream.avg_intensity,
                                                 list(characteristics.items())[11][1]])
                experiment_14_result_list.append(copy.deepcopy(experiment_14_sublist))

            file_name = 'experiment_14_' + local_queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_14_result_list))
        except ValueError as e:
            print(str(e))
            file_name = 'experiment_14_except' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_14_result_list))
        except Exception:
            traceback.print_exc(file=sys.stderr)
    else:
        file_name = 'experiment_14_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_14_result_list = ast.literal_eval(res_line)

    build_plot(experiment_14_result_list,
               '\\khi_{1,2} от h',
               'h',
               '\\khi_{1,2}',
               'experiment_14')


def experiment_15(queueing_system: ColdReserveQueueingSystem, read_file=False):
    """
    Зависимость количества переключений с прибора-2 на прибор-1, \khi_{2,1} от 1/kappa_2 при различных значениях h.

    :param queueing_system: ColdReserveQueueingSystem
    :return: None
    """

    linux_check_cpu_temperature()

    print('Experiment 15 launched!')

    experiment_15_result_list = []

    if not read_file:
        try:
            local_queueing_system = copy.deepcopy(queueing_system)

            queries_matrices = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices)
            break_matrices = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)

            for break_coef in (0.1, 1, 10):
                linux_check_cpu_temperature()

                break_matrices_1 = [matr * break_coef for matr in break_matrices]  # todo Check if copy necessary
                local_queueing_system.set_MAP_break_stream(break_matrices_1[0], break_matrices_1[1])

                filename = '../experiment_results/' + 'experiment_15_' + queueing_system.name + '.qsc'
                local_queueing_system.print_characteristics(filename)

                experiment_15_sublist = [local_queueing_system.break_stream.avg_intensity, []]
                for queries_coef in tqdm([i / 50 for i in range(1, 1001)]):
                    queries_matrices_1 = [matr * queries_coef for matr in queries_matrices]
                    local_queueing_system.queries_stream.set_transition_matrices(queries_matrices_1)
                    characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                    experiment_15_sublist[1].append([1 / local_queueing_system.switch2_1_stream.avg_intensity,
                                                 list(characteristics.items())[12][1]])
                experiment_15_result_list.append(copy.deepcopy(experiment_15_sublist))

            file_name = 'experiment_15_' + local_queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_15_result_list))
        except ValueError as e:
            print(str(e))
            file_name = 'experiment_15_except' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_15_result_list))
        except Exception:
            traceback.print_exc(file=sys.stderr)
    else:
        file_name = 'experiment_15_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_15_result_list = ast.literal_eval(res_line)

    build_plot(experiment_15_result_list,
               '\\khi_{2,1} от h',
               'h',
               '\\khi_{2,1}',
               'experiment_15')
