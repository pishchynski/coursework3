from cold_reserve_qs import *


def run_test():
    queueing_system = ColdReserveQueueingSystem(p_num=150)
    queueing_system.set_MAP_queries_stream(np.array([[-19]]), np.array([[19]]))
    queueing_system.set_MAP_break_stream(np.array([[-0.00001]]), np.array([[0.00001]]))
    queueing_system.set_PH_serv_unit1_stream(np.array([[1]]), np.array([[-20]]))
    queueing_system.set_PH_serv_unit2_stream(np.array([[1]]), np.array([[-5]]))
    queueing_system.set_PH_switch1_2_stream(np.array([[0.05, 0.95]]), np.array([[-1.86075, 0.], [0., -146.9994]]))
    queueing_system.set_PH_switch2_1_stream(np.array([[0.05, 0.95]]), np.array([[-5.58225, 0.], [0., -440.9982]]))
    queueing_system.set_PH_repair_stream(np.array([[1]]), np.array([[-10000]]))
    characteristics, vect_p_l = queueing_system.calc_characteristics()

    for i, vect in enumerate(vect_p_l):
        print("P_{} = ".format(str(i)), np.sum(vect))

    for i, charact in enumerate(characteristics):
        print("{}. ".format(str(i)), charact)


if __name__ == '__main__':
    run_test()
