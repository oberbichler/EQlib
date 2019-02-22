import unittest
import EQLib as eq

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


if __name__ == '__main__':
    unittest.main()