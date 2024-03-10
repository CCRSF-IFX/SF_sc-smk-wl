import unittest
import sys
from utils.utils import generate_DME_link 

class TestGenerateDMELink(unittest.TestCase):
    def test_valid_project_id(self):
        project_id = "NyallLondon_CS032345_16scRNASeqTCRLib_081122"
        expected_link = "https://hpcdmeweb.nci.nih.gov/collection?action=view&path=/FNL_SF_Archive/PI_Lab_Nyall_London/Project_NyallLondon_CS032345_16scRNASeqTCRLib_081122&source=browse&init"
        actual_link = generate_DME_link(project_id)
        self.assertEqual(actual_link, expected_link)
    
    #def test_invalid_project_id(self):
    #    project_id = "Invalid_CS032445_16scRNASeqTCRLib"
    #    with self.assertRaises(SystemExit) as context:
    #        generate_DME_link(project_id)
    #    # Check the exit code or message if needed
    #    self.assertIn("{DMELINK}", str(context.exception))

    def test_invalid_project_id(self):
        project_id = "Invalid_CS032445_16scRNASeqTCRLib"
        self.assertEqual(generate_DME_link(project_id), "{DMELINK}")

if __name__ == '__main__':
    unittest.main()
