import unittest
import EQlib as eq

class TestParameter(unittest.TestCase):

    def test_parameter_constructor(self):
        parameter = eq.Parameter(ref_value=1, act_value=2, target=3, result=4)

        self.assertEqual(parameter.ref_value, 1.0)
        self.assertEqual(parameter.act_value, 2.0)
        self.assertEqual(parameter.target, 3.0)
        self.assertEqual(parameter.result, 4.0)

    def test_parameter_getters_and_setters(self):
        parameter = eq.Parameter(ref_value=1, act_value=2, target=3, result=4)

        parameter.ref_value = 5.1
        self.assertEqual(parameter.ref_value, 5.1)

        parameter.act_value = 6.2
        self.assertEqual(parameter.act_value, 6.2)

        parameter.target = 7.3
        self.assertEqual(parameter.target, 7.3)

        parameter.result = 8.4
        self.assertEqual(parameter.result, 8.4)

    def test_parameter_dof_delta(self):
        parameter = eq.Parameter(ref_value=1, act_value=3, target=4, result=7)

        dof = parameter.dof

        self.assertEqual(dof.delta, 2.0)

        dof.delta = 5

        self.assertEqual(dof.delta, 5.0)

        self.assertEqual(parameter.ref_value, 1.0)
        self.assertEqual(parameter.act_value, 6.0)
        self.assertEqual(parameter.target, 4.0)
        self.assertEqual(parameter.result, 7.0)

    def test_parameter_dof_residual(self):
        parameter = eq.Parameter(ref_value=1, act_value=3, target=4, result=7)

        dof = parameter.dof

        self.assertEqual(dof.residual, 3.0)

        dof.residual = 5

        self.assertEqual(dof.residual, 5.0)

        self.assertEqual(parameter.ref_value, 1.0)
        self.assertEqual(parameter.act_value, 3.0)
        self.assertEqual(parameter.target, 4.0)
        self.assertEqual(parameter.result, 9.0)

    def test_parameter_pickle_and_copy(self):
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

    def test_parameter_dof_equality(self):
        parameter_a = eq.Parameter(ref_value=1, act_value=2, target=3, result=4)
        parameter_b = eq.Parameter(ref_value=1, act_value=2, target=3, result=4)

        dof_a_1 = parameter_a.dof
        dof_a_2 = parameter_a.dof

        dof_b_1 = parameter_b.dof
        dof_b_2 = parameter_b.dof

        assert(dof_a_1 == dof_a_2)
        assert(dof_a_1 != dof_b_1)
        assert(dof_b_1 == dof_b_2)

if __name__ == '__main__':
    unittest.main()