import ast
import sys
import traceback
from itertools import cycle

import matplotlib.pyplot as plt
from tqdm import tqdm

sys.path.append("../")
from src.cold_reserve_qs import *

def build_plot(experiment_result_list, experiment_name, x_label, y_label, leg_label, file_name='experiment_plot', file_type='png'):
    """
    Builds experiment plots and saves them to file.

    :param experiment_result_list: list containing experiment results in two dimensions
    :param experiment_name: str with experiment name to be plot's title
    :param x_label: str with x-axis label to be displayed on plot
    :param y_label: str with y-axis label to be displayed on plot
    :param file_type: str with file type ('eps', 'png'). Default is 'png'
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
                experiment_1_sublist = [local_queueing_system.queries_stream.c_cor, []]

                for queries_coef in tqdm([i / 12.9 for i in range(1, 52)]):
                    linux_check_cpu_temperature(notify=False)

                    queries_matrices_1 = [matr * queries_coef for matr in queries_matrices_0]

                    local_queueing_system.queries_stream.set_transition_matrices(queries_matrices_1)
                    # local_queueing_system.set_MAP_queries_stream(queries_matrices_1[0], queries_matrices_1[1])

                    characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

                    experiment_1_sublist[1].append([local_queueing_system.queries_stream.avg_intensity,
                                                     characteristics[13]])
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
            file_name = 'experiment_1_except' + queueing_system.name + '.qsr'
            with open('../experiment_results/' + file_name, mode='w') as res_file:
                res_file.write(str(experiment_1_result_list))
    else:
        file_name = 'experiment_1_' + queueing_system.name + '.qsr'
        with open('../experiment_results/' + file_name, mode='r') as res_file:
            res_line = res_file.readline()
        experiment_1_result_list = ast.literal_eval(res_line)

    build_plot(experiment_1_result_list,
               'u_ от lambda',
               'lambda',
               'u_',
               'c_cor_q',
               'experiment_1')