# fluxpy: a Python toolkit for metabolic modelling analysis utils
# Authors: Haris Zafeiropoulos (2024)
# Licensed under GNU LGPL.3, see LICENCE file

import os
import unittest
from fluxpy import utils

tests_dir = os.path.dirname(__file__)
root_dir = os.path.dirname(tests_dir)
models_dir = os.path.join(root_dir, "ext_data")

class TestUtils(unittest.TestCase):

    def compare_models(self):

        input_model_list = [os.path.join(models_dir, filename) for filename in os.listdir(models_dir)]
        comp = utils.CompareModels(input_model_list)

        set.assertTrue(comp.shared_metabolites() == 47)



utils.CompareModels

if __name__ == "__main__":
    unittest.main()


