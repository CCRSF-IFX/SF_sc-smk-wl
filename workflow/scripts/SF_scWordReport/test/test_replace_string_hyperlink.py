import unittest
from unittest.mock import Mock, patch
from utils.utils import replace_string_hyperlink, add_hyperlink 

class TestReplaceStringHyperlink(unittest.TestCase):
    @patch('utils.utils.add_hyperlink')
    def test_replace_string_hyperlink(self, mock_add_hyperlink):
        # Mock the document object
        doc = Mock()

        # Mock the paragraphs and runs
        p1 = Mock()
        p1.text = 'This is a {DMELINK} example.'
        p1.runs = [Mock()]  # Mock the runs attribute

        p2 = Mock()
        p2.text = 'Another {DMELINK} example.'
        p2.runs = [Mock()]  # Mock the runs attribute

        doc.paragraphs = [p1, p2]

        # Run the function
        result = replace_string_hyperlink(doc, '{DMELINK}', 'http://www.example.com')

        # Verify that the add_hyperlink function was called with the correct arguments
        mock_add_hyperlink.assert_called_once_with(p1, 'http://www.example.com', 'http://www.example.com', '4169e1', False)

        # Verify that the text replacement occurred in the document
        #for inline in [p1.runs, p2.runs]:
        #    for run_mock in inline:
        #        run_mock.text = run_mock.text.replace('{DMELINK}', '')

        # Verify the return value
        self.assertEqual(result, 1)

if __name__ == '__main__':
    unittest.main()


