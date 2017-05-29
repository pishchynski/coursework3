import sys
sys.path.append("../")
from src.experiments import *
from src.utils import characteristics_loc


def run_test_1(verbose=False):
    queueing_system = ColdReserveQueueingSystem(name="System_Main_Final_1", p_num=3)
    queueing_system.set_BMAP_queries_stream(matrD_0=np.array([[-86., 0.01], [0.02, -2.76]]) / 57.2,
                                            matrD=np.array([[85., 0.99], [0.2, 2.54]]) / 57.2,
                                            n=3,
                                            q=0.8)
    # queueing_system.set_BMAP_queries_stream(np.array([[-6.3408, 1.87977 * (10 ** (-6))], [1.87977 * (10 ** (-6)), -0.13888]]) / 8.7,
    #                                         np.array([[6.3214, 0.01939], [0.10822, 0.03066]]) / 8.7,
    #                                         q=0.8,
    #                                         n=3)
    # queueing_system.set_BMAP_queries_stream(np.array([[-6., 0.01], [0.02, -2.76]]) * 5.32,
    #                                         np.array([[5., 0.99], [0.2, 2.54]]) * 5.32,
    #                                         q=0.8,
    #                                         n=3)
    # queueing_system.set_MAP_queries_stream(np.array([[-0.575]]), np.array([[0.575]]))
    queueing_system.set_MAP_break_stream(np.array([[-8.110725, 0.], [0., -0.26325]]),
                                         np.array([[8.0568, 0.053925], [0.146625, 0.116625]]))
    # queueing_system.set_MAP_break_stream(np.array([[-6.3408, 1.87977 / (10**6)], [1.87977 / (10**6), -0.13888]]),
    #                                      np.array([[6.3214, 0.01939], [0.10822, 0.03066]]))
    # queueing_system.set_MAP_break_stream(np.array([[-8.5]]), np.array([[8.5]]))
    queueing_system.set_PH_serv_unit1_stream(np.array([[1., 0.]]), np.array([[-20., 20.], [0., -20.]]))
    queueing_system.set_PH_serv_unit2_stream(np.array([[1., 0.]]), np.array([[-2., 2.], [0., -2.]]))
    queueing_system.set_PH_switch1_2_stream(np.array([[0.05, 0.95]]), np.array([[-1.86075, 0.], [0., -146.9994]]))
    queueing_system.set_PH_switch2_1_stream(np.array([[0.05, 0.95]]), np.array([[-18.6075, 0.], [0., -1469.994]]))
    queueing_system.set_PH_repair_stream(np.array([[0.98, 0.02]]), np.array([[-100., 0.], [0., -0.002]]))
    characteristics, vect_p_l = queueing_system.calc_characteristics(verbose=verbose)

    for i, vect in enumerate(vect_p_l):
        print("P_{} = ".format(str(i)), np.sum(vect))

    for i, charact in enumerate(characteristics):
        print("{}. ".format(str(i)), characteristics_loc[i], ':', charact)

    return queueing_system


def run_test_2(verbose=False):
    queueing_system = ColdReserveQueueingSystem(name="System_Main_Final_2", p_num=3)

    queueing_system.set_BMAP_queries_stream(np.array([[-8.110725, 0.], [0., -0.26325]]) / 21,
                                            np.array([[8.0568, 0.053925], [0.146625, 0.116625]]) / 21,
                                            q=0.8,
                                            n=3)

    # queueing_system.set_MAP_break_stream(np.array([[-6.3408, 1.87977 / (10**6)], [1.87977 / (10**6), -0.13888]]),
    #                                      np.array([[6.3214, 0.01939], [0.10822, 0.03066]]))
    queueing_system.set_MAP_break_stream(np.array([[-86., 0.01], [0.02, -2.76]]) * 1.57 / 2.4,
                                         np.array([[85., 0.99], [0.2, 2.54]]) * 1.57 / 2.4)
    # queueing_system.set_MAP_break_stream(np.array([[-4.5]]), np.array([[4.5]]))
    queueing_system.set_PH_serv_unit1_stream(np.array([[1., 0.]]), np.array([[-20., 20.], [0., -20.]]))
    queueing_system.set_PH_serv_unit2_stream(np.array([[1., 0.]]), np.array([[-2., 2.], [0., -2.]]))
    queueing_system.set_PH_switch1_2_stream(np.array([[0.05, 0.95]]), np.array([[-1.86075, 0.], [0., -146.9994]]))
    queueing_system.set_PH_switch2_1_stream(np.array([[0.05, 0.95]]), np.array([[-18.6075, 0.], [0., -1469.994]]))
    queueing_system.set_PH_repair_stream(np.array([[0.98, 0.02]]), np.array([[-100., 0.], [0., -0.002]]))
    characteristics, vect_p_l = queueing_system.calc_characteristics(verbose=verbose)

    for i, vect in enumerate(vect_p_l):
        print("P_{} = ".format(str(i)), np.sum(vect))

    for i, charact in enumerate(characteristics):
        print("{}. ".format(str(i)), characteristics_loc[i], ':', charact)

    return queueing_system


def run_test_3(verbose=False):
    queueing_system = ColdReserveQueueingSystem(name="System_Main_Final_3", p_num=3)

    queueing_system.set_BMAP_queries_stream(matrD_0=np.array([[-86., 0.01], [0.02, -2.76]]) / 30.7,
                                            matrD=np.array([[85., 0.99], [0.2, 2.54]]) / 30.7,
                                            n=3,
                                            q=0.8)

    queueing_system.set_MAP_break_stream(np.array([[-8.110725, 0.], [0., -0.26325]]),
                                         np.array([[8.0568, 0.053925], [0.146625, 0.116625]]))

    queueing_system.set_PH_serv_unit1_stream(np.array([[1., 0.]]), np.array([[-20., 20.], [0., -20.]]))
    queueing_system.set_PH_serv_unit2_stream(np.array([[1., 0.]]), np.array([[-2., 2.], [0., -2.]]))
    queueing_system.set_PH_switch1_2_stream(np.array([[0.05, 0.95]]), np.array([[-1.86075, 0.], [0., -146.9994]]))
    queueing_system.set_PH_switch2_1_stream(np.array([[0.05, 0.95]]), np.array([[-18.6075, 0.], [0., -1469.994]]))

    queueing_system.set_PH_repair_stream(np.array([[0.98, 0.02]]), np.array([[-100., 0.], [0., -0.002]]))
    # queueing_system.set_PH_repair_stream(np.array([[1., 0.]]), np.array([[-1., 1.], [0., -1.]]) / 5)
    # queueing_system.set_PH_repair_stream(np.array([[1.]]), np.array([[-0.1]]))

    characteristics, vect_p_l = queueing_system.calc_characteristics(verbose=verbose)

    for i, vect in enumerate(vect_p_l):
        print("P_{} = ".format(str(i)), np.sum(vect))

    for i, charact in enumerate(characteristics):
        print("{}. ".format(str(i)), characteristics_loc[i], ':', charact)

    return queueing_system


if __name__ == '__main__':
    read_file = True

    # q_system = run_test_1(verbose=True)
    # experiment_1(q_system, read_file=read_file)
    # experiment_1_1(q_system, read_file=read_file)

    q_system = run_test_2(verbose=True)
    # experiment_2(q_system, read_file=read_file)
    # experiment_4(q_system, read_file=read_file)
    experiment_5(q_system, read_file=read_file)
    # experiment_6(q_system, read_file=read_file)
    # experiment_7(q_system, read_file=read_file)

    # q_system = run_test_3(verbose=True)
    # experiment_3(q_system, read_file=read_file)
