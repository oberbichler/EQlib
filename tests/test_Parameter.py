import unittest
import EQlib as eq

class TestParameter(unittest.TestCase):

    def test_constructor(self):
        parameter = eq.Parameter(ref_value=1, act_value=2, target=3, result=4)

        self.assertEqual(parameter.ref_value, 1.0)
        self.assertEqual(parameter.act_value, 2.0)
        self.assertEqual(parameter.target, 3.0)
        self.assertEqual(parameter.result, 4.0)
        self.assertEqual(parameter.lower_bound, -float('inf'))
        self.assertEqual(parameter.upper_bound, float('inf'))
        self.assertEqual(parameter.isfixed, False)
        self.assertEqual(parameter.name, '')

    def test_getters_and_setters(self):
        parameter = eq.Parameter(ref_value=1, act_value=2, target=3, result=4)

        parameter.ref_value = 5.1
        self.assertEqual(parameter.ref_value, 5.1)

        parameter.act_value = 6.2
        self.assertEqual(parameter.act_value, 6.2)

        parameter.target = 7.3
        self.assertEqual(parameter.target, 7.3)

        parameter.result = 8.4
        self.assertEqual(parameter.result, 8.4)

        parameter.lower_bound = 9.1
        self.assertEqual(parameter.lower_bound, 9.1)

        parameter.upper_bound = 2.3
        self.assertEqual(parameter.upper_bound, 2.3)

        parameter.isfixed = True
        self.assertEqual(parameter.isfixed, True)

        parameter.name = 'Example'
        self.assertEqual(parameter.name, 'Example')

    def test_delta(self):
        parameter = eq.Parameter(ref_value=1, act_value=3, target=4, result=7)

        self.assertEqual(parameter.delta, 2.0)

        parameter.delta = 5

        self.assertEqual(parameter.delta, 5.0)

        self.assertEqual(parameter.ref_value, 1.0)
        self.assertEqual(parameter.act_value, 6.0)
        self.assertEqual(parameter.target, 4.0)
        self.assertEqual(parameter.result, 7.0)

    def test_residual(self):
        parameter = eq.Parameter(ref_value=1, act_value=3, target=4, result=7)

        self.assertEqual(parameter.residual, -3.0)

        parameter.residual = -5

        self.assertEqual(parameter.residual, -5.0)

        self.assertEqual(parameter.ref_value, 1.0)
        self.assertEqual(parameter.act_value, 3.0)
        self.assertEqual(parameter.target, 4.0)
        self.assertEqual(parameter.result, 9.0)

    def test_pickle_and_copy(self):
        import pickle
        from copy import copy, deepcopy

        for op in [copy, deepcopy, lambda o: pickle.loads(pickle.dumps(o))]:
            parameter_a = eq.Parameter(ref_value=1, act_value=2, target=3,
                                       result=4)

            parameter_b = op(parameter_a)

            self.assertEqual(parameter_a.ref_value, parameter_b.ref_value)
            self.assertEqual(parameter_a.act_value, parameter_b.act_value)
            self.assertEqual(parameter_a.target, parameter_b.target)
            self.assertEqual(parameter_a.result, parameter_b.result)

            parameter_a.ref_value = 5.1
            parameter_a.act_value = 6.2
            parameter_a.target = 7.3
            parameter_a.result = 8.4

            self.assertEqual(parameter_a.ref_value, 5.1)
            self.assertEqual(parameter_a.act_value, 6.2)
            self.assertEqual(parameter_a.target, 7.3)
            self.assertEqual(parameter_a.result, 8.4)

            self.assertEqual(parameter_b.ref_value, 1.0)
            self.assertEqual(parameter_b.act_value, 2.0)
            self.assertEqual(parameter_b.target, 3.0)
            self.assertEqual(parameter_b.result, 4.0)

    def test_equality(self):
        parameter_a = eq.Parameter(ref_value=1, act_value=2, target=3, result=4)
        parameter_b = eq.Parameter(ref_value=1, act_value=2, target=3, result=4)

        assert(parameter_a == parameter_a)
        assert(parameter_a != parameter_b)
        assert(parameter_b == parameter_b)

if __name__ == '__main__':
    unittest.main()