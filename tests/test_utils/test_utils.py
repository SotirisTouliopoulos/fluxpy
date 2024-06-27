# fluxpy: a Python toolkit for metabolic modelling analysis utils
# Authors: Haris Zafeiropoulos (2024)
# Licensed under GNU LGPL.3, see LICENCE file

import os
import pandas as pd
import unittest
import cobra
from fluxpy.utils import utils

tests_dir = os.path.dirname(__file__)
root_dir = os.path.dirname(os.path.dirname(tests_dir))
models_dir = os.path.join(root_dir, "ext_data")
ecoli = os.path.join(models_dir, "e_coli_core.xml")

class TestUtils(unittest.TestCase):

    # def compare_models(self):

    #     input_model_list = [os.path.join(models_dir, filename) for filename in os.listdir(models_dir)]
    #     comp = utils.CompareModels(input_model_list)

    #     set.assertTrue(comp.shared_metabolites() == 47)

    def test_non_zero_fluxes(self):
        test_series = pd.Series([0,0,0,4,0,0,6])
        assumed_non_zero = utils.nonzero_fluxes(test_series)
        self.assertTrue( len(assumed_non_zero) == 2)


    def test_get_rxns_producing_consuming_met(self):
        model = cobra.io.read_sbml_model(ecoli)
        fba_cons, fba_prod = utils.get_rxns_producing_consuming_met(model.metabolites.accoa_c.id, model= model)
        pfba_sol = cobra.flux_analysis.pfba(model)
        pfba_cons, pfba_prod = utils.get_rxns_producing_consuming_met(model.metabolites.accoa_c.id, model= model,
                                                                      flux_vector = pfba_sol.fluxes)
        self.assertTrue(len(fba_cons) == 2 and len(fba_prod) == 1)
        self.assertTrue(len(pfba_cons) == 2 and len(pfba_prod) == 1)




if __name__ == "__main__":
    unittest.main()


