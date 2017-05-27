import sys
sys.path.append("../")
from src.experiments import *
from src.utils import characteristics_loc


def run_test(verbose=False):
    queueing_system = ColdReserveQueueingSystem(name="Simple System", p_num=200)
    queueing_system.set_MAP_queries_stream(np.array([[-33]]), np.array([[33]]))
    queueing_system.set_MAP_break_stream(np.array([[-33]]), np.array([[33]]))
    queueing_system.set_PH_serv_unit1_stream(np.array([[1]]), np.array([[-140]]))
    queueing_system.set_PH_serv_unit2_stream(np.array([[1]]), np.array([[-70]]))
    queueing_system.set_PH_switch1_2_stream(np.array([[1]]), np.array([[-10]]))
    queueing_system.set_PH_switch2_1_stream(np.array([[1]]), np.array([[-100]]))
    queueing_system.set_PH_repair_stream(np.array([[1]]), np.array([[-18]]))
    characteristics, vect_p_l = queueing_system.calc_characteristics(verbose=verbose)

    for i, vect in enumerate(vect_p_l):
        print("P_{} = ".format(str(i)), np.sum(vect))

    for i, charact in enumerate(characteristics):
        print("{}. ".format(str(i)), characteristics_loc[i], ':', charact)

    return queueing_system


if __name__ == '__main__':
    q_system = run_test(verbose=False)
    read_file = False
    experiment_1(q_system, read_file=read_file)
    # experiment_2(q_system, read_file=read_file)
    # experiment_3(q_system, read_file=read_file)
    # experiment_4(q_system, read_file=read_file)
    # experiment_5(q_system, read_file=read_file)
    # experiment_6(q_system, read_file=read_file)
    # experiment_7(q_system, read_file=read_file)
    # experiment_8(q_system, read_file=read_file)
    # experiment_9(q_system, read_file=read_file)
    # experiment_10(q_system, read_file=read_file)
    # experiment_11(q_system, read_file=read_file)
    # experiment_12(q_system, read_file=read_file)
