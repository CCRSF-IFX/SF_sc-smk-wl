# xiess4 on 2022/09/06: add function getDMElink, replace_string_hyperlink and add_hyperlink

import sys
import docx 

def generate_DME_link(project_id):
    """Generate DME link based on project ID
   
    This function should work for the cases in which first name only has one uppcase letter. 
    For first name with two uppercase letters, the DME link will be incorrect. But the it should 
    be rare. 

    Parameters
    ----------
    project_id : str
        The project ID. e.g. NyallLondon_CS032345_16scRNASeqTCRLib_081122

    Returns
    -------
    dmelink : str
        Http link to the DME 
    """
    ele = project_id.split("_")
    PI_name = ele[0]
    
    loc_uc_letter = [idx for idx in range(len(PI_name)) if PI_name[idx].isupper()]
    # First letter 
    if len(loc_uc_letter) < 2:  
        sys.stderr.write("At least two uppercase letters are expected. Please check the project ID")
        return "{DMELINK}"
    else:
        first_name = PI_name[0:loc_uc_letter[1]]
        last_name = PI_name[loc_uc_letter[1]:]
        dmelink = f"https://hpcdmeweb.nci.nih.gov/collection?action=view&path=/FNL_SF_Archive/PI_Lab_{first_name}_{last_name}/Project_{project_id}&source=browse&init"
        return dmelink

def add_hyperlink(paragraph, url, text, color, underline):
    """
    A function that places a hyperlink within a paragraph object.

    :param paragraph: The paragraph we are adding the hyperlink to.
    :param url: A string containing the required url
    :param text: The text displayed for the url
    :return: The hyperlink object
    :reference: https://github.com/python-openxml/python-docx/issues/74#issuecomment-261169410
    """

    # This gets access to the document.xml.rels file and gets a new relation id value
    part = paragraph.part
    r_id = part.relate_to(url, docx.opc.constants.RELATIONSHIP_TYPE.HYPERLINK, is_external=True)

    # Create the w:hyperlink tag and add needed values
    hyperlink = docx.oxml.shared.OxmlElement('w:hyperlink')
    hyperlink.set(docx.oxml.shared.qn('r:id'), r_id, )

    # Create a w:r element
    new_run = docx.oxml.shared.OxmlElement('w:r')

    # Create a new w:rPr element
    rPr = docx.oxml.shared.OxmlElement('w:rPr')

    # Add color if it is given
    if not color is None:
      c = docx.oxml.shared.OxmlElement('w:color')
      c.set(docx.oxml.shared.qn('w:val'), color)
      rPr.append(c)

    # Remove underlining if it is requested
    if not underline:
      u = docx.oxml.shared.OxmlElement('w:u')
      u.set(docx.oxml.shared.qn('w:val'), 'none')
      rPr.append(u)

    # Join all the xml elements together add add the required text to the w:r element
    new_run.append(rPr)
    new_run.text = text
    hyperlink.append(new_run)

    paragraph._p.append(hyperlink)

    return hyperlink

def replace_string_hyperlink(doc, before, after):
    """Replace text (before) with new text (after) with hyperlink
    Parameters
    ----------
    doc : word doc object
        Word document object
    before: str
        The string to be replaced (e.g. '{DMELINK}')
    after: str
        The string should be used (e.g. 'httt://www.google.com')
    Returns
    -------
    return: 1 

    """
    for p in doc.paragraphs:
        if before in p.text:
            inline = p.runs

            for i in range(len(inline)):
                if before in inline[i].text:
                    # remove {DMELINK} 
                    text = inline[i].text.replace(before, "")
                    inline[i].text = text
                    #Add a hyperlink with the normal formatting (royalblue with underline)
                    hyperlink = add_hyperlink(p, after, after, "4169e1", False)
    return 1
