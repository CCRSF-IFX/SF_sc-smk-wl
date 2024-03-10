import unittest
from unittest.mock import patch, MagicMock

from utils.utils import add_hyperlink


class TestAddHyperlink(unittest.TestCase):

    @patch('docx.document.DocumentPart.relate_to')
    def test_add_hyperlink(self, mock_relate_to):
        """Tests basic functionality of add_hyperlink."""

        # Mock document part and relationship ID
        mock_part = MagicMock()
        mock_relate_to.return_value = 'rel_id'
        mock_paragraph = MagicMock()
        mock_paragraph.part = mock_part

        # Call the function
        url = 'https://www.example.com'
        text = 'Click here'
        color = 'FF0000'  # Red
        underline = True
        hyperlink = add_hyperlink(mock_paragraph, url, text, color, underline)

        # Assertions
        mock_relate_to.assert_called_once_with(url, docx.opc.constants.RELATIONSHIP_TYPE.HYPERLINK, is_external=True)

        # Verify hyperlink element structure (limited for brevity)
        self.assertIsInstance(hyperlink, docx.oxml.shared.OxmlElement)
        self.assertEqual(hyperlink.tag, docx.oxml.shared.qn('w:hyperlink'))

        # Verify run element with text and color (assuming color is set)
        run = hyperlink.element.find('w:r')
        self.assertIsInstance(run, docx.oxml.shared.OxmlElement)
        self.assertEqual(run.text, text)
        color_element = run.find('w:rPr/w:color')
        if color:  # Check color only if provided
            self.assertIsInstance(color_element, docx.oxml.shared.OxmlElement)
            self.assertEqual(color_element.get(docx.oxml.shared.qn('w:val')), color.upper())

        # Verify underline element (if not disabled)
        underline_element = run.find('w:rPr/w:u')
        if underline:
            self.assertIsInstance(underline_element, docx.oxml.shared.OxmlElement)
            self.assertEqual(underline_element.get(docx.oxml.shared.qn('w:val')), 'single')
        else:
            self.assertIsNone(underline_element)

    @patch('docx.document.DocumentPart.relate_to')
    def test_add_hyperlink_no_color(self, mock_relate_to):
        """Tests add_hyperlink without color argument."""

        mock_part = MagicMock()
        mock_relate_to.return_value = 'rel_id'
        mock_paragraph = MagicMock()
        mock_paragraph.part = mock_part

        url = 'https://www.example.com'
        text = 'Click here'
        color = None  # No color
        underline = True
        add_hyperlink(mock_paragraph, url, text, color, underline)

        # Verify no color element is added
        run = mock_paragraph._p.element.find('w:r')
        color_element = run.find('w:rPr/w:color')
        self.assertIsNone(color_element)

    @patch('docx.document.DocumentPart.relate_to')
    def test_add_hyperlink_no_underline(self, mock_relate_to):
        """Tests add_hyperlink with underline set to False."""

        mock_part = MagicMock()
        mock_relate_to.return_value = 'rel_id'
        mock_paragraph = MagicMock()
        mock_paragraph.part = mock_part

        url = 'https://www.example.com'
        text = 'Click here'
        color = 'FF0000'  # Red
        underline = False  # No underline
        add_hyperlink(mock_paragraph, url, text, color, underline)

        # Verify no underline element is added
        run = mock_paragraph._p.element.find('w:r')
        underline_element = run.find('w:rPr/w:u')
        self.assertIsNone(underline_element)
