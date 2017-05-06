from tqdm import tqdm

from cold_reserve_qs import *


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

    for switch2_1_coef in (0.1, 1, 10):
        linux_check_cpu_temperature()

        switch2_1_matr_1 = copy.deepcopy(switch2_1_matr) * switch2_1_coef
        local_queueing_system.set_PH_switch2_1_stream(switch2_1_vect, switch2_1_matr_1)
        experiment_1_sublist = [local_queueing_system.switch2_1_stream.avg_intensity]
        for switch1_2_coef in tqdm([i / 1000 if i < 10 else i / 100 if i < 50 else i / 10 for i in range(1, 101)]):
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

    for switch2_1_coef in (0.1, 1, 10):
        linux_check_cpu_temperature()

        switch2_1_matr_1 = copy.deepcopy(switch2_1_matr) * switch2_1_coef  # todo Check if copy necessary
        local_queueing_system.set_PH_switch2_1_stream(switch2_1_vect, switch2_1_matr_1)
        experiment_2_sublist = [local_queueing_system.switch2_1_stream.avg_intensity]
        for switch1_2_coef in tqdm([i / 1000 if i < 10 else i / 100 if i < 50 else i / 10 for i in range(1, 101)]):
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
    """
    Зависимость L от lambda при различных значениях h.

    :param queueing_system: ColdReserveQueueingSystem
    :return: None
    """

    linux_check_cpu_temperature()

    print('Experiment 3 launched!')

    local_queueing_system = copy.deepcopy(queueing_system)

    queries_matrices = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices)
    break_matrices = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)

    experiment_3_result_list = []

    for break_coef in (0.1, 1, 10):
        linux_check_cpu_temperature()

        break_matrices_1 = [matr * break_coef for matr in break_matrices]  # todo Check if copy necessary
        local_queueing_system.set_MAP_break_stream(break_matrices_1[0], break_matrices_1[1])
        experiment_3_sublist = [local_queueing_system.break_stream.avg_intensity]
        for queries_coef in tqdm([i / 1000 if i < 10 else i / 100 if i < 50 else i / 50 for i in range(1, 101)]):
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
    """
    Зависимость L от lambda при различных коэффициентах корреляции c_cor во входном потоке.
    :param queueing_system: ColdReserveQueueingSystem
    :return: None
    """

    linux_check_cpu_temperature()

    print('Experiment 4 launched!')

    local_queueing_system = copy.deepcopy(queueing_system)

    queries_matrices = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices)
    matrD = copy.deepcopy(local_queueing_system.queries_stream.matrD)

    experiment_4_result_list = []

    for cor_coef in (0.1, 1, 10):
        linux_check_cpu_temperature()

        matrD_10 = copy.deepcopy(matrD)
        matrD_10[0] *= cor_coef
        local_queueing_system.set_BMAP_queries_stream(queries_matrices[0], matrD_10, q=local_queueing_system.queries_stream.q, n=local_queueing_system.n)
        experiment_4_sublist = [local_queueing_system.queries_stream.c_cor]
        for queries_coef in tqdm([i / 1000 if i < 10 else i / 100 if i < 50 else i / 50 for i in range(1, 101)]):
            matrD_0_1 = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices[0]) * queries_coef
            matrD_1 = copy.deepcopy(local_queueing_system.queries_stream.matrD) * queries_coef
            local_queueing_system.set_BMAP_queries_stream(matrD_0_1, matrD_1, q=local_queueing_system.queries_stream.q, n=local_queueing_system.n)
            characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

            experiment_4_sublist.append([local_queueing_system.queries_stream.avg_intensity,
                                         list(characteristics.items())[2][1]])
        experiment_4_result_list.append(copy.deepcopy(experiment_4_sublist))

    file_name = 'experiment_4_' + local_queueing_system.name + '.qsr'
    with open('experiment_results/' + file_name, mode='w') as res_file:
        res_file.write(str(experiment_4_result_list))


def experiment_5(queueing_system: ColdReserveQueueingSystem):
    """
    Зависимость V от lambda при различных коэффициентах корреляции c_cor во входном потоке.

    :param queueing_system: ColdReserveQueueingSystem
    :return: None
    """

    linux_check_cpu_temperature()

    print('Experiment 5 launched!')

    local_queueing_system = copy.deepcopy(queueing_system)

    queries_matrices = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices)
    matrD = copy.deepcopy(local_queueing_system.queries_stream.matrD)

    experiment_5_result_list = []

    for cor_coef in (0.1, 1, 10):
        linux_check_cpu_temperature()

        matrD_10 = copy.deepcopy(matrD)
        matrD_10[0] *= cor_coef
        local_queueing_system.set_BMAP_queries_stream(queries_matrices[0], matrD_10, q=local_queueing_system.queries_stream.q, n=local_queueing_system.n)
        experiment_5_sublist = [local_queueing_system.queries_stream.c_cor]
        for queries_coef in tqdm([i / 1000 if i < 10 else i / 100 if i < 50 else i / 50 for i in range(1, 101)]):
            matrD_0_1 = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices[0]) * queries_coef
            matrD_1 = copy.deepcopy(local_queueing_system.queries_stream.matrD) * queries_coef
            local_queueing_system.set_BMAP_queries_stream(matrD_0_1, matrD_1, q=local_queueing_system.queries_stream.q, n=local_queueing_system.n)
            characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

            experiment_5_sublist.append([local_queueing_system.queries_stream.avg_intensity,
                                         list(characteristics.items())[3][1]])
        experiment_5_result_list.append(copy.deepcopy(experiment_5_sublist))

    file_name = 'experiment_5_' + local_queueing_system.name + '.qsr'
    with open('experiment_results/' + file_name, mode='w') as res_file:
        res_file.write(str(experiment_5_result_list))


def experiment_6(queueing_system: ColdReserveQueueingSystem):
    """
    Зависимость вероятности доступности P_1+ доступности прибора-1 от h при различных коэффициентах корреляции c_cor
    в потоке поломок
    :param queueing_system: ColdReserveQueueingSystem
    :return: None
    """

    linux_check_cpu_temperature()

    print('Experiment 6 launched!')

    local_queueing_system = copy.deepcopy(queueing_system)

    break_matrices = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)

    experiment_6_result_list = []

    for cor_coef in (0.1, 1, 10):
        linux_check_cpu_temperature()

        break_matr0 = copy.deepcopy(break_matrices[0])
        break_matr0[0] *= cor_coef
        break_matr1 = copy.deepcopy(break_matrices[1])
        break_matr1[0] *= cor_coef

        break_matrices_0 = [break_matr0, break_matr1]

        local_queueing_system.set_MAP_break_stream(break_matrices_0[0], break_matrices_0[1])
        experiment_6_sublist = [local_queueing_system.break_stream.c_cor]
        for break_coef in tqdm([i / 1000 if i < 10 else i / 100 if i < 50 else i / 10 for i in range(1, 101)]):
            break_matrices_1 = copy.deepcopy(break_matrices_0)
            break_matrices_1[0] *= break_coef
            break_matrices_1[1] *= break_coef
            local_queueing_system.set_MAP_break_stream(break_matrices_1[0], break_matrices_1[1])
            characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

            experiment_6_sublist.append([local_queueing_system.break_stream.avg_intensity,
                                         list(characteristics.items())[8][1]])
        experiment_6_result_list.append(copy.deepcopy(experiment_6_sublist))

    file_name = 'experiment_6_' + local_queueing_system.name + '.qsr'
    with open('experiment_results/' + file_name, mode='w') as res_file:
        res_file.write(str(experiment_6_result_list))


def experiment_7(queueing_system: ColdReserveQueueingSystem):
    """
    Зависимость вероятности доступности P_1- недоступности прибора-1 от h при различных коэффициентах корреляции c_cor
    в потоке поломок
    :param queueing_system: ColdReserveQueueingSystem
    :return: None
    """

    linux_check_cpu_temperature()

    print('Experiment 7 launched!')

    local_queueing_system = copy.deepcopy(queueing_system)

    break_matrices = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)

    experiment_7_result_list = []

    for cor_coef in (0.1, 1, 10):
        linux_check_cpu_temperature()

        break_matr0 = copy.deepcopy(break_matrices[0])
        break_matr0[0] *= cor_coef
        break_matr1 = copy.deepcopy(break_matrices[1])
        break_matr1[0] *= cor_coef

        break_matrices_0 = [break_matr0, break_matr1]

        local_queueing_system.set_MAP_break_stream(break_matrices_0[0], break_matrices_0[1])
        experiment_7_sublist = [local_queueing_system.break_stream.c_cor]
        for break_coef in tqdm([i / 1000 if i < 10 else i / 100 if i < 50 else i / 10 for i in range(1, 101)]):
            break_matrices_1 = copy.deepcopy(break_matrices_0)
            break_matrices_1[0] *= break_coef
            break_matrices_1[1] *= break_coef
            local_queueing_system.set_MAP_break_stream(break_matrices_1[0], break_matrices_1[1])
            characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

            experiment_7_sublist.append([local_queueing_system.break_stream.avg_intensity,
                                         1 - list(characteristics.items())[8][1]])
        experiment_7_result_list.append(copy.deepcopy(experiment_7_sublist))

    file_name = 'experiment_7_' + local_queueing_system.name + '.qsr'
    with open('experiment_results/' + file_name, mode='w') as res_file:
        res_file.write(str(experiment_7_result_list))


def experiment_8(queueing_system: ColdReserveQueueingSystem):
    """
    Зависимость L от h при различных коэффициентах вариации c_var времени ремонта

    :param queueing_system: ColdReserveQueueingSystem
    :return: None
    """

    linux_check_cpu_temperature()

    print('Experiment 8 launched!')

    local_queueing_system = copy.deepcopy(queueing_system)

    break_matrices = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)

    experiment_8_result_list = []

    for repair_vect_elem in (0.1, 0.5, 0.9):
        linux_check_cpu_temperature()

        local_queueing_system.set_PH_repair_stream(np.array([[repair_vect_elem, 1 - repair_vect_elem]]),
                                                   local_queueing_system.repair_stream.repres_matr)
        experiment_8_sublist = [local_queueing_system.repair_stream.c_var]
        for break_coef in tqdm([i / 1000 if i < 10 else i / 100 if i < 50 else i / 10 for i in range(1, 101)]):
            break_matrices_1 = copy.deepcopy(break_matrices)
            break_matrices_1[0] *= break_coef
            break_matrices_1[1] *= break_coef
            local_queueing_system.set_MAP_break_stream(break_matrices_1[0], break_matrices_1[1])
            characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

            experiment_8_sublist.append([local_queueing_system.break_stream.avg_intensity,
                                         list(characteristics.items())[2][1]])
        experiment_8_result_list.append(copy.deepcopy(experiment_8_sublist))

    file_name = 'experiment_8_' + local_queueing_system.name + '.qsr'
    with open('experiment_results/' + file_name, mode='w') as res_file:
        res_file.write(str(experiment_8_result_list))

def experiment_9(queueing_system: ColdReserveQueueingSystem):
    """
    Зависимость пропускной способности от lambda при различных значениях h.

    :param queueing_system: ColdReserveQueueingSystem
    :return: None
    """

    linux_check_cpu_temperature()

    print('Experiment 9 launched!')

    local_queueing_system = copy.deepcopy(queueing_system)

    queries_matrices = copy.deepcopy(local_queueing_system.queries_stream.transition_matrices)
    break_matrices = copy.deepcopy(local_queueing_system.break_stream.transition_matrices)

    experiment_9_result_list = []

    for break_coef in (0.1, 1, 10):
        linux_check_cpu_temperature()

        break_matrices_1 = [matr * break_coef for matr in break_matrices]  # todo Check if copy necessary
        local_queueing_system.set_MAP_break_stream(break_matrices_1[0], break_matrices_1[1])
        experiment_9_sublist = [local_queueing_system.break_stream.avg_intensity]
        for queries_coef in tqdm([i / 1000 if i < 10 else i / 100 if i < 50 else i / 50 for i in range(1, 101)]):
            queries_matrices_1 = [matr * queries_coef for matr in queries_matrices]
            local_queueing_system.queries_stream.set_transition_matrices(queries_matrices_1)
            characteristics, vect_p_l = local_queueing_system.calc_characteristics(verbose=False)

            experiment_9_sublist.append([local_queueing_system.queries_stream.avg_intensity,
                                         list(characteristics.items())[1][1]])
        experiment_9_result_list.append(copy.deepcopy(experiment_9_sublist))

    file_name = 'experiment_9_' + local_queueing_system.name + '.qsr'
    with open('experiment_results/' + file_name, mode='w') as res_file:
        res_file.write(str(experiment_9_result_list))
