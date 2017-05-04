from email._encoded_words import _QByteMap

from cold_reserve_qs import *
from tqdm import tqdm


def build_graph():
    None


def experiment_1(queueing_system: ColdReserveQueueingSystem):
    """
    Зависимость вероятности P- того, что оба прибора недоступны, от среднего времени 1/kappa_1
    переключения с прибора-1 на прибор-2 при разных средних временах 1/kappa_2 переключения с прибора-2 на прибор-1.

    :param queueing_system: ColdReserveQueueingSystem
    :return: None
    """

    linux_check_cpu_temperature()

    print('Experiment 1 launched!')

    local_queueing_system = copy.deepcopy(queueing_system)  # todo Check how impacts memory consumption!

    switch2_1_matr = copy.deepcopy(local_queueing_system.switch2_1_stream.repres_matr)
    switch2_1_vect = copy.deepcopy(local_queueing_system.switch2_1_stream.repres_vect)

    switch1_2_matr = copy.deepcopy(local_queueing_system.switch1_2_stream.repres_matr)
    switch1_2_vect = copy.deepcopy(local_queueing_system.switch1_2_stream.repres_vect)

    experiment_1_result_list = []

    for switch2_1_coef in (0.01, 1, 100):
        switch2_1_matr_1 = copy.deepcopy(switch2_1_matr) * switch2_1_coef
        local_queueing_system.set_PH_switch2_1_stream(switch2_1_vect, switch2_1_matr_1)
        experiment_1_sublist = [local_queueing_system.switch2_1_stream.avg_intensity]
        for switch1_2_coef in tqdm([i / 1000 if i < 10 else i / 100 if i < 50 else i / 10 for i in range(1, 10)]):
            switch1_2_matr_1 = copy.deepcopy(switch1_2_matr) * switch1_2_coef
            local_queueing_system.set_PH_switch1_2_stream(switch1_2_vect, switch1_2_matr_1)
            characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

            experiment_1_sublist.append(
                [1 / local_queueing_system.switch1_2_stream.avg_intensity, list(characteristics.items())[10][1]])
        experiment_1_result_list.append(copy.deepcopy(experiment_1_sublist))

    file_name = 'experiment_1_' + local_queueing_system.name + '.qsr'
    with open('experiment_results/' + file_name, mode='w') as res_file:
        res_file.write(str(experiment_1_result_list))


def experiment_2(queueing_system: ColdReserveQueueingSystem):
    """
     Зависимость числа переключений режимов 1_2 + 2_1 от среднего времени 1/kappa_1 переключения с прибора-1 на прибор-2
     при различных средних временах 1/kappa_2 переключений с прибора-2 на прибор-1

    :param queueing_system: ColdReserveQueueingSystem
    :return: None
    """

    linux_check_cpu_temperature()

    print('Experiment 2 launched!')

    local_queueing_system = copy.deepcopy(queueing_system)

    switch2_1_matr = copy.deepcopy(local_queueing_system.switch2_1_stream.repres_matr)
    switch2_1_vect = copy.deepcopy(local_queueing_system.switch2_1_stream.repres_vect)

    switch1_2_matr = copy.deepcopy(local_queueing_system.switch1_2_stream.repres_matr)
    switch1_2_vect = copy.deepcopy(local_queueing_system.switch1_2_stream.repres_vect)

    experiment_2_result_list = []

    for switch2_1_coef in (0.01, 1, 100):
        switch2_1_matr_1 = copy.deepcopy(switch2_1_matr) * switch2_1_coef  # todo Check if copy necessary
        local_queueing_system.set_PH_switch2_1_stream(switch2_1_vect, switch2_1_matr_1)
        experiment_2_sublist = [local_queueing_system.switch2_1_stream.avg_intensity]
        for switch1_2_coef in tqdm([i / 1000 if i < 10 else i / 100 if i < 50 else i / 10 for i in range(1, 10)]):
            switch1_2_matr_1 = copy.deepcopy(switch1_2_matr) * switch1_2_coef
            local_queueing_system.set_PH_switch1_2_stream(switch1_2_vect, switch1_2_matr_1)
            characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

            experiment_2_sublist.append([1 / local_queueing_system.switch1_2_stream.avg_intensity,
                                         list(characteristics.items())[11][1] + list(characteristics.items())[12][1]])
        experiment_2_result_list.append(copy.deepcopy(experiment_2_sublist))

        file_name = 'experiment_2_' + local_queueing_system.name + '.qsr'
        with open('experiment_results/' + file_name, mode='w') as res_file:
            res_file.write(str(experiment_2_result_list))


def experiment_3(queueing_system: ColdReserveQueueingSystem):
    linux_check_cpu_temperature()

    print('Experiment 3 launched!')

    local_queueing_system = copy.deepcopy(queueing_system)

    queries_matrices = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices)
    break_matrices = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)

    experiment_3_result_list = []

    for break_coef in (0.01, 1, 100):
        break_matrices_1 = [matr * break_coef for matr in break_matrices]  # todo Check if copy necessary
        local_queueing_system.set_MAP_break_stream(break_matrices_1[0], break_matrices_1[1])
        experiment_3_sublist = [local_queueing_system.break_stream.avg_intensity]
        for queries_coef in tqdm([i / 1000 if i < 10 else i / 100 if i < 50 else i / 10 for i in range(1, 10)]):
            queries_matrices_1 = [matr * queries_coef for matr in queries_matrices]
            local_queueing_system.queries_stream.set_transition_matrices(queries_matrices_1)
            characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

            experiment_3_sublist.append([local_queueing_system.queries_stream.avg_intensity,
                                         list(characteristics.items())[2][1]])
        experiment_3_result_list.append(copy.deepcopy(experiment_3_sublist))

        file_name = 'experiment_3_' + local_queueing_system.name + '.qsr'
        with open('experiment_results/' + file_name, mode='w') as res_file:
            res_file.write(str(experiment_3_result_list))


def experiment_4(queueing_system: ColdReserveQueueingSystem):
    linux_check_cpu_temperature()

    print('Experiment 4 launched!')

    local_queueing_system = copy.deepcopy(queueing_system)

    queries_matrices = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices)
    matrD = copy.deepcopy(local_queueing_system.queries_stream.matr_D)

    experiment_4_result_list = []

    for q_value in (0.1, 0.5, 0.8):
        local_queueing_system.set_BMAP_queries_stream(queries_matrices[0], matrD)
        experiment_4_sublist = [local_queueing_system.queries_stream.c_cor]
        for queries_coef in tqdm([i / 1000 if i < 10 else i / 100 if i < 50 else i / 10 for i in range(1, 10)]):
            matrD_0_1 = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices[0]) * queries_coef
            matrD_1 = copy.deepcopy(local_queueing_system.queries_stream.matrD) * queries_coef
            local_queueing_system.set_BMAP_queries_stream(matrD_0_1, matrD_1, q=q_value, n=local_queueing_system.n)
            characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

            experiment_4_sublist.append([local_queueing_system.queries_stream.avg_intensity,
                                         list(characteristics.items())[2][1]])
        experiment_4_result_list.append(copy.deepcopy(experiment_4_sublist))

        file_name = 'experiment_4_' + local_queueing_system.name + '.qsr'
        with open('experiment_results/' + file_name, mode='w') as res_file:
            res_file.write(str(experiment_4_result_list))


def experiment_5(queueing_system: ColdReserveQueueingSystem):
    linux_check_cpu_temperature()

    print('Experiment 5 launched!')

    local_queueing_system = copy.deepcopy(queueing_system)

    queries_matrices = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices)
    matrD = copy.deepcopy(local_queueing_system.queries_stream.matr_D)

    experiment_5_result_list = []

    for q_value in (0.1, 0.5, 0.8):
        local_queueing_system.set_BMAP_queries_stream(queries_matrices[0], matrD)
        experiment_5_sublist = [local_queueing_system.queries_stream.c_cor]
        for queries_coef in tqdm([i / 1000 if i < 10 else i / 100 if i < 50 else i / 10 for i in range(1, 10)]):
            matrD_0_1 = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices[0]) * queries_coef
            matrD_1 = copy.deepcopy(local_queueing_system.queries_stream.matrD) * queries_coef
            local_queueing_system.set_BMAP_queries_stream(matrD_0_1, matrD_1, q=q_value, n=local_queueing_system.n)
            characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

            experiment_5_sublist.append([local_queueing_system.queries_stream.avg_intensity,
                                         list(characteristics.items())[3][1]])
        experiment_5_result_list.append(copy.deepcopy(experiment_5_sublist))

        file_name = 'experiment_5_' + local_queueing_system.name + '.qsr'
        with open('experiment_results/' + file_name, mode='w') as res_file:
            res_file.write(str(experiment_5_result_list))