import xml.etree.ElementTree as ET
from xml.dom import minidom
import re

sequencer = {'NEXTSEQ': 'NextSeq 500', 'NEXTSEQ 2000': 'NextSeq 2000', 'HISEQ': 'HiSeq 4000', 'ISEQ': 'iSeq', 'MISEQ': 'MiSeq', 'NOVASEQ': 'NovaSeq 6000'}

#def mergeDict(dict1, dict2):
#    dict3 = {**dict1, **dict2}
#    return dict3

def to_parseable(tree):
    t = ET.tostring(tree)
    t = t.lower()
    return ET.fromstring(t)

def getRTA(root):
    #return next(root.iter('rtaversion')).text
    if len(root.findall('./rtaversion')) == 1:
      return root.findall('./rtaversion')[0].text

def getSequencer(root):
    if len(root.findall('instrumenttype')) == 1:
        return(root.findall('instrumenttype')[0].text)
    elif len(root.findall('application')) == 1:
        return root.findall('application')[0].text.split(' ')[0]
    elif len(root.findall('.//applicationname')) == 1:
        return root.findall('.//applicationname')[0].text.split(' ')[0]
    else:
        return 'No Info'

def getChemistry(root):
    if len(root.findall('chemistry')) == 1:
        return root.findall('chemistry')[0].text
    elif len(root.findall('.//flowcellmode')) == 1:
        if len(re.findall('P\d', root.findall('.//flowcellmode')[0].text.upper())) == 1:
            return(re.findall('P\d', root.findall('.//flowcellmode')[0].text.upper())[0])
        else:
            return root.findall('.//flowcellmode')[0].text.upper()
    elif len(root.findall('./setup/sbs')) == 1:
        return root.findall('./setup/sbs')[0].text.upper()


def getWorkflowType(root):
    if len(root.findall('workflowtype')) >0:
        return root.findall('workflowtype')[0].text
    else:
        return 'No Info'

def getInstrumentID(root):
    if len(root.findall('instrumentid')) == 1:
        return root.findall('instrumentid')[0].text.upper()
    elif len(root.findall('instrumentname')) == 1:
        return root.findall('instrumentname')[0].text.upper()
    elif len(root.findall('scannerid')) == 1:
        return root.findall('scannerid')[0].text.upper()
    elif len(root.findall('.//scannerid')) == 1:
        return root.findall('.//scannerid')[0].text.upper()
    elif len(root.findall('instrumentserialnumber')) == 1:
        return root.findall('instrumentserialnumber')[0].text.upper()


def getFlowcell(root):
    if len(root.findall('./flowcellrfidtag/serialnumber')) == 1:
        return root.findall('./flowcellrfidtag/serialnumber')[0].text.upper()
    elif len(root.findall('.//flowcellserialbarcode')) == 1:
        return root.findall('.//flowcellserialbarcode')[0].text.upper()
    elif len(root.findall('./flowcelleepromtag/serialnumber')) == 1:
        return root.findall('./flowcelleepromtag/serialnumber')[0].text.upper()
    elif len(root.findall('.//barcode')) == 1:
        return root.findall('.//barcode')[0].text.upper()
    elif len(root.findall('./flowcellserialnumber')) == 1:
        return root.findall('./flowcellserialnumber')[0].text.upper()
    else:
        return 'No Info'

def getReadLength(root):
    reads = list()
    indices = list()
    if len(root.findall('./setup/read1')) == 1:
      reads.append(root.findall('./setup/read1')[0].text)
      reads.append(root.findall('./setup/read2')[0].text)
    elif len(root.findall('read1numberofcycles')) == 1:
      reads.append(root.findall('read1numberofcycles')[0].text)
      reads.append(root.findall('read2numberofcycles')[0].text)
    elif len(root.findall('reads')) == 1:
      for read in root.find('reads'):
        if read.attrib['isindexedread'] == 'n':
          reads.append(read.attrib['numcycles'])
        else:
          indices.append(read.attrib['numcycles'])
    elif len(root.findall('./plannedcycles/read1')) == 1:
        reads.append(root.findall('./plannedcycles/read1')[0].text)
        if len(root.findall('./plannedcycles/read2')) == 1:
            reads.append(root.findall('./plannedcycles/read2')[0].text)
    elif len(root.findall('./plannedreads')) == 1:
      for read in root.findall('./plannedreads/read'):
        if 'read' in read.attrib['readname']:
          reads.append(read.attrib['cycles'])
        else:
          indices.append(read.attrib['cycles'])
    else:
      return "Format Not Expected"
    if len(root.findall('./setup/indexread1')) == 1:
        indices.append(root.findall('./setup/indexread1')[0].text)
        if len(root.findall('./setup/indexread2')) == 1:
            indices.append(root.findall('./setup/indexread2')[0].text)
    elif len(root.findall('./setup/index1read')) == 1:
        indices.append(root.findall('./setup/index1read')[0].text)
        if len(root.findall('./setup/index2read')) == 1:
            indices.append(root.findall('./setup/index2read')[0].text)
    elif len(root.findall('indexread1numberofcycles')) == 1:
      indices.append(root.findall('indexread1numberofcycles')[0].text)
      if len(root.findall('indexread2numberofcycles')) == 1:
          indices.append(root.findall('indexread2numberofcycles')[0].text)
    elif len(root.findall('./plannedcycles/index1')) == 1:
        indices.append(root.findall('./plannedcycles/index1')[0].text)
        if len(root.findall('./plannedcycles/index2')) == 1:
            indices.append(root.findall('./plannedcycles/index2')[0].text)

    if len(reads) > 2 or len(indices) > 2:
        return "Format Not Expected\nReads: " + str(reads) + "\nIndices: " + str(indices)
    else:
        output = dict()
        for i in range(2):
            if len(reads) < i+1:
                output['read'+str(i+1)] = str(0)
            else:
                output['read'+str(i+1)] = reads[i]
            if len(indices) < i+1:
                output['idx'+str(i+1)] = str(0)
            else:
                output['idx'+str(i+1)] = indices[i]
        return(output)
    #if len(reads) == 1:
    #  return (reads[0] + ' (1x' + reads[0] + ' cycles)')
    #elif len(reads) > 2:
    #  return "Format Not Expected: " + str(reads)
    #else:
    #  if reads[0] == reads[1]:
    #    return(reads[0] + ' (2x' + reads[0]+ ' cycles)'
    #  else:
    #    return 'R1: ' + reads[0] + ', R2: ' + reads[1] + ' cycles'


#import glob
#files = glob.glob('runParameters/*xml')
#doc = '/is2/projects/CCR-SF/scratch/illumina/RawData_MiSeq/200309_M01595_0119_000000000-CV4CT/runParameters.xml'
#tree = ET.parse(doc)
#root = tree.getroot()

#root = to_parseable(root)
#getRTA(root)
#getSequencer(root)
#getReadLength(root)
#getChemistry(root)

def parseXML(xmlfile):
    # create element tree object
    tree = ET.parse(xmlfile)
    root = tree.getroot()
    root = to_parseable(root)
    info = getReadLength(root)
    info['rta'] = getRTA(root)
    if info['rta'] == None:
        info['rta'] = "No Info"
    info['chemistry'] = getChemistry(root)
    if info['chemistry'] == None:
        info['chemistry'] = "No Info"
    info['platform'] = sequencer.get(getSequencer(root).upper(), getSequencer(root))
    if 'novaseq' in info['platform'].lower():
        info['workflow'] = getWorkflowType(root)
    else:
        info['workflow'] = 'No Info'
    info['instrumentid'] = getInstrumentID(root)
    info['flowcell'] = getFlowcell(root)
    return(info)

#for doc in files:
#    print(parseXML(doc))
