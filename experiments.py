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

    switch2_1_matr = copy.deepcopy(queueing_system.switch2_1_stream.repres_matr)
    switch2_1_vect = copy.deepcopy(queueing_system.switch2_1_stream.repres_vect)

    switch1_2_matr = copy.deepcopy(queueing_system.switch1_2_stream.repres_matr)
    switch1_2_vect = copy.deepcopy(queueing_system.switch1_2_stream.repres_vect)

    experiment_1_list = []

    for switch1_2_coef in tqdm(np.arange(0.05, 10., 0.05)):
        switch1_2_matr_1 = copy.deepcopy(switch1_2_matr) * switch1_2_coef
        queueing_system.set_PH_switch1_2_stream(switch1_2_vect, switch1_2_matr_1)
        characteristics, vect_p_l = queueing_system.calc_characteristics(verbose=False)

        experiment_1_list.append([1 / queueing_system.switch1_2_stream.avg_intensity, list(characteristics.items())[10]])

    print(experiment_1_list)


def experiment_2(queueing_system: ColdReserveQueueingSystem):
    """
     Зависимость числа переключений режимов 1_2 + 2_1 от среднего времени 1/kappa_1 переключения с прибора-1 на прибор-2
     при различных средних временах 1/kappa_2 переключений с прибора-2 на прибор-1

    :param queueing_system: ColdReserveQueueingSystem
    :return: None
    """

    switch2_1_matr = copy.deepcopy(queueing_system.switch2_1_stream.repres_matr)
    switch2_1_vect = copy.deepcopy(queueing_system.switch2_1_stream.repres_vect)

    switch1_2_matr = copy.deepcopy(queueing_system.switch1_2_stream.repres_matr)
    switch1_2_vect = copy.deepcopy(queueing_system.switch1_2_stream.repres_vect)
