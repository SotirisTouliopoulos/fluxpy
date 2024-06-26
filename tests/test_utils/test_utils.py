# fluxpy: a Python toolkit for metabolic modelling analysis utils
# Authors: Haris Zafeiropoulos (2024)
# Licensed under GNU LGPL.3, see LICENCE file

import os
import pandas as pd
import unittest
import cobra
from fluxpy.utils import utils

tests_dir = os.path.dirname(__file__)
root_dir = os.path.dirname(tests_dir)
models_dir = os.path.join(root_dir, "ext_data")
ecoli = os.path.join(models_dir, "e_coli_core.xml")

class TestUtils(unittest.TestCase):

    # def compare_models(self):

    #     input_model_list = [os.path.join(models_dir, filename) for filename in os.listdir(models_dir)]
    #     comp = utils.CompareModels(input_model_list)

    #     set.assertTrue(comp.shared_metabolites() == 47)

    def test_non_zero_fluxes(self):

        # model = cobra.io.read_sbml_model(ecoli)
        # solution = model.optimize()

        test_series = pd.Series([0,0,0,4,0,0,6])
        assumed_non_zero = utils.nonzero_fluxes(test_series)
        self.assertTrue( len(assumed_non_zero) == 2)


if __name__ == "__main__":
    unittest.main()


