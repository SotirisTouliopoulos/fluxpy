import os
import unittest
from fluxpy.utils.analysis import producing_or_consuming_a_met
import cobra


tests_dir = os.path.dirname(__file__)
root_dir = os.path.dirname(os.path.dirname(tests_dir))
models_dir = os.path.join(root_dir, "ext_data/models")
ecoli = os.path.join(models_dir, "e_coli_core.xml")


class TestUtils(unittest.TestCase):

    def test_NameIdConverter(self):
        self.assertTrue(producing_or_consuming_a_met(model=ecoli, reaction_id="CS", metabolite_id="accoa_c") == "consuming")




