import asymmetric_immiscibility_simulations.SlabEquilibrate as equilibrate


def test_nothing():
    args = equilibrate.parse_args('-T 100 '.split())
    equilibrate.main(args)
