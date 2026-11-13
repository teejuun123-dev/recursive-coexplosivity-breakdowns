"""Small smoke tests of function imports and simple mathematical helpers.

These are not tests of empirical results or econometric validity.
Run manually: python -m unittest discover -s tests
"""
import unittest
import numpy as np
from analysis import weight_one, sPSY_df, recursive_residual

class SmokeTests(unittest.TestCase):
    def test_weight_one_shape(self):
        w = weight_one(np.array([1.0, 2.0, 3.0, 4.0, 5.0]), 3)
        self.assertEqual(w.shape, (5, 5))

    def test_spsy_df_scalar(self):
        result = sPSY_df(np.arange(1.0, 20.0))
        self.assertTrue(np.isscalar(result) or np.isnan(result))

    def test_recursive_residual_output(self):
        y = np.arange(1.0, 20.0)
        x = y * 1.1
        yy, xx, res = recursive_residual(y, x, 0)
        self.assertEqual(xx.shape[1], 2)
        self.assertEqual(res.shape[1], 1)

if __name__ == '__main__':
    unittest.main()
