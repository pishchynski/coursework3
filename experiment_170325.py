from cold_reserve_qs import *

def run_test():
    queueing_system = ColdReserveQueueingSystem(p_num=150)
    queueing_system.set_BMAP_queries_stream(matrD_0=np.array([[-86., 0.01], [0.02, -2.76]]),
                                            matrD=np.array([[85., 0.99], [0.2, 2.54]]),
                                            n=3,
                                            q=0.8)
    queueing_system.set_MAP_break_stream(np.array([[-8., 1.], [2., -12.]]), np.array([[2., 5.], [4., 6.]]))
    queueing_system.set_PH_serv_unit1_stream(np.array([[0.2, 0.8]]), np.array([[-170., 15.], [40., -210.]]))
    queueing_system.set_PH_serv_unit2_stream(np.array([[0.9, 0.1]]), np.array([[-110., 80.], [10., -150.]]))
    queueing_system.set_PH_switch1_2_stream(np.array([[0.9, 0.1]]), np.array([[-220., 160.], [20., -300.]]))
    queueing_system.set_PH_switch2_1_stream(np.array([[0.9999, 0.0001]]), np.array([[-100000., 0.], [0., -230.]]))
    queueing_system.set_PH_recover_stream(np.array([[0.2, 0.8]]), np.array([[-17., 1.5], [4., -21.]]))
    characteristics, vect_p_l= queueing_system.calc_characteristics()

    for i, vect in enumerate(vect_p_l):
        print("P_{} = ".format(str(i)), np.sum(vect))

    for i, charact in enumerate(characteristics):
        print("{}. ".format(str(i)), charact)


if __name__ == '__main__':
    run_test()