#!/usr/bin/env python3
#2021/04/21 shent2: added attributes of data_curator and data_owner in createPILabJson
#2021/04/21 shent2: retrieve read lengths from RanParameters.xml and the reference from count_path/config.py
#2021/05/07 chenv3: changed to use only first run name when trying to grab run parameters"
#2021/07/22 shent2: 1. added {os.path.dirname(fastqPath)}/input_samplesheet.csv to run_name_supplement.tar
#                   2. added retrieveReadLength() to retrieve read lengths from fastq.gz when RunParameters.xml/RunInfo.xml does not exist
#2022/02/23 shent2: added the regular expression to parse the format of the PI name as "LastName, FirstName"
#2022/03/13 shent2: added "meta2json_md5.py -p 4 -l fastqFileList.txt/analysisFileList.txt --rewrite" to generate the checksum and add to metadata.json for uploading
#2022/04/04 shent2: modified to parse library.csv and tar the count results based on column name in libraries.csv
import re,os,sys,warnings,argparse,subprocess,glob,gzip,datetime
from xml.dom import minidom
from pathlib import PurePath

parser = argparse.ArgumentParser(description="""Read metadata and gnerate json files
                                                for collections and objects.""")
parser.add_argument("-m", "--metadata", metavar="metadata.txt", dest="fmetadata",
                    action="store", type=str, required=False,
                    help="input metadata.txt (optional)")
parser.add_argument("-c", "--count_path", metavar="count_path",
                    dest="count_path", action="store", type=str,
                    required=False, help="count folder (optional)")
parser.add_argument("-a", "--aggr", metavar="aggr_path",
                    dest="aggr_path", action="store", type=str,
                    required=False, help="aggregate folder (optional, comma-separated if multiple)")
parser.add_argument("-p", "--project", metavar="projectName",
                    dest="fprojectName", action="store", type=str,
                    required=False, help="project name (optional)")
parser.add_argument("-r", "--run_name", metavar="run_name", dest="run_name", 
                    action="store", type=str, required=True, 
                    help="run name (required, comma-separated if multiple), e.g., 180405_J00170_0093_AHV2CKBBXX")
parser.add_argument("-f", "--fastq_path", metavar="fastq_path", dest="fastq_path",
                    action="store", type=str, required=True, 
                    help="fastq_path folder (required, comma-separated if multiple), e.g., outs/fastq_path")
parser.add_argument("-o", "--output", metavar="file_list", dest="foutput",
                    action="store", default="file_list.txt", type=str, 
                    help="file list to archive (default: %(default)s)")
args = parser.parse_args()

class classSample:
    'sample object containing metadata'
    def __init__(self, attributes, values):
        self.attribute2value = dict() 
        for i in range(0, len(attributes)):
            if values[i] == "":
                values[i] = "Unspecified"
            self.attribute2value[attributes[i]] = values[i]
        self.flag = 1
    def addAttribute(self, attribute, value):
        self.attribute2value[attribute] = value

def retrieveReadLength(fileName):
    header = ""
    readSeq = ""
    with gzip.open(fileName, 'rb') as IN:
        lineNum = 1
        for line in IN:
            if lineNum == 1:
                header = line[0:-1]
            elif lineNum == 2:
                readSeq = line[0:-1]
            else:
                break
            lineNum += 1
    return(len(readSeq))

def parseXml(xmlRunParametersPath, xmlRunInfoPath, instrument, readLengths, indexLengths):
    RTAVersion = "Unknown"
    flowcell = "Unknown"
    workFlowTypeValue = "Unknown"
    flowcellModeValue = "Unknown"
    chemistryValue = "Unknown"
    chemistryVersionValue = "Unknown"
    xmldoc = minidom.parse(xmlRunParametersPath)
    tagRTA = ""
    if instrument == "NovaSeq 6000" or instrument == "NextSeq 2000" or instrument == "iSeq":
        tagRTA = "RtaVersion"
    elif instrument == "HiSeq" or  instrument == "NextSeq 550" or instrument == "MiSeq":
        tagRTA = "RTAVersion"
    elif instrument == "NovaSeq Xplus":
        RTAVersion = "RTA4"
    if tagRTA:
        cyc = xmldoc.getElementsByTagName(tagRTA)
    else:
        cyc = []
    if len(cyc) > 0:
        RTAVersion = cyc[0].childNodes[0].nodeValue
    workFlowType = xmldoc.getElementsByTagName('WorkflowType')
    if len(workFlowType) > 0:
        workFlowTypeValue = workFlowType[0].childNodes[0].nodeValue
    flowcellMode = xmldoc.getElementsByTagName('FlowCellMode')
    if len(flowcellMode) > 0:
        flowcellModeValue = flowcellMode[0].childNodes[0].nodeValue
    chemistry = xmldoc.getElementsByTagName('Chemistry')
    if len(chemistry) > 0:
        chemistryValue = chemistry[0].childNodes[0].nodeValue
    chemistryVersion = xmldoc.getElementsByTagName('ChemistryVersion')
    if len(chemistryVersion) > 0:
        chemistryVersionValue = chemistryVersion[0].childNodes[0].nodeValue
    application = xmldoc.getElementsByTagName('Application')
    if len(application) > 0:
        applicaitonValue = application[0].childNodes[0].nodeValue
    consumableInfo = xmldoc.getElementsByTagName('ConsumableInfo')
    if len(consumableInfo) > 0:
        for consumableInfoIndex in range(len(consumableInfo)):
            if consumableInfo[consumableInfoIndex].getElementsByTagName('Type')[0].childNodes[0].nodeValue == 'FlowCell':
                if len(consumableInfo[consumableInfoIndex].getElementsByTagName('Name')) > 0:
                    flowcellModeValue = consumableInfo[consumableInfoIndex].getElementsByTagName('Name')[0].childNodes[0].nodeValue
                else:
                    flowcellModeValue = consumableInfo[consumableInfoIndex].getElementsByTagName('Mode')[0].childNodes[0].nodeValue
            else:
                continue

    xmldoc = minidom.parse(xmlRunInfoPath)
    runInfoReads = xmldoc.getElementsByTagName('Read')
    if len(runInfoReads) > 0:
        for runInfoRead in runInfoReads:
            if runInfoRead.getAttribute("IsIndexedRead") == "N":
                readLengths.append(runInfoRead.getAttribute("NumCycles"))
            elif runInfoRead.getAttribute("IsIndexedRead") == "Y":
                indexLengths.append(runInfoRead.getAttribute("NumCycles"))
    runInfoFlowcell = xmldoc.getElementsByTagName('Flowcell')
    if len(runInfoFlowcell) > 0:
        flowcell = runInfoFlowcell[0].childNodes[0].nodeValue
    return (RTAVersion, flowcell, workFlowTypeValue, flowcellModeValue, chemistryValue, chemistryVersionValue)

def genome2annotation(sampleName, **samples):
    if "hg19" in samples[sampleName].attribute2value["ReferenceGenome"]:
        return("ENSEMBL_v70")
    elif "hg38" in samples[sampleName].attribute2value["ReferenceGenome"]:
        return("ENSEMBL_GRCh38_v79")
    elif "mm9" in samples[sampleName].attribute2value["ReferenceGenome"]:
        return("ENSEMBL_NCBI37_mm9")
    elif "mm10" in samples[sampleName].attribute2value["ReferenceGenome"]:
        return("ENSEMBL_GRCm38_mm10")
    else:
        return("Unknown")

def runCommand(command):
    try:
        retcode = subprocess.call(command, shell=True)
        if retcode < 0:
            print >>sys.stderr, command + " was terminated by signal", -retcode
        #else:
        #    print >>sys.stderr, command + " returned", retcode
    except OSError as e:
        print >>sys.stderr, "Execution failed:", e

def createPILabJson(PIPath, PIName):
    OUT = open(PIPath + ".metadata.json", "w")
    OUT.write("{\"metadataEntries\":[\n" +
              "    {\"attribute\":\"collection_type\",\"value\":\"PI_Lab\"},\n" +
              #"    {\"attribute\":\"pi_name\",\"value\":\"" + PIName + "\"},\n" +
              "    {\"attribute\":\"data_owner_designee\",\"value\":\"" + PIName + "\"},\n" +
              "    {\"attribute\":\"data_generator\",\"value\":\"Zhao, Yongmei\"},\n" +
              "    {\"attribute\":\"data_owner\",\"value\":\"" + PIName + "\"}\n" +
              "    ]\n}"
    )
    OUT.close()

def createProjectJson(projectPath, sampleName, **samples):
    OUT = open(projectPath + ".metadata.json", "w")
    OUT.write("{\"metadataEntries\":[\n" +
              "    {\"attribute\":\"collection_type\",\"value\":\"Project\"},\n" +
              "    {\"attribute\":\"project_id_CSAS_NAS\",\"value\":\"" + samples[sampleName].attribute2value["CSAS_NAS"] + "\"},\n" +
              "    {\"attribute\":\"project_name\",\"value\":\"" + samples[sampleName].attribute2value["ProjectName"] + "\"},\n" +
              "    {\"attribute\":\"contact_name\",\"value\":\"" + samples[sampleName].attribute2value["LabContact"] + "\"},\n" +
              "    {\"attribute\":\"bioinformatics_contact\",\"value\":\"" + samples[sampleName].attribute2value["BioinformaticsContact"] + "\"},\n" +
              "    {\"attribute\":\"project_start_date\",\"value\":\"" + samples[sampleName].attribute2value["ProjectOpenDate"] + "\"},\n" +
              "    {\"attribute\":\"project_completed_date\",\"value\":\"" + datetime.date.today().isoformat() + "\"},\n" +
              "    {\"attribute\":\"grant_funding_agent\",\"value\":\"NIH\"}\n" +
              "    ]\n}"
    )
    OUT.close()

def createFlowcellJson(flowcellPath, flowcellID, runName, sampleName, **samples):
    runDate = ""
    oRunDate = re.compile("(\d{6,})_(.+)")
    sRunDate = oRunDate.search(runName)
    if sRunDate:
        year = 2000 + int(sRunDate.group(1)[0:2])
        month = sRunDate.group(1)[2:4]
        day = sRunDate.group(1)[4:6]
        runDate = str(year) + "-" + month + "-" + day
    else:
        sys.stderr.write("Could not identify the run date from the convention" + 
                         " of the run name, " + runName + ".\n")
        sys.exit(1)

    OUT = open(flowcellPath + ".metadata.json", "w")
    OUT.write("{\"metadataEntries\":[\n" +
              "    {\"attribute\":\"collection_type\",\"value\":\"Flowcell\"},\n" +
              "    {\"attribute\":\"flowcell_id\",\"value\":\"" + flowcellID + "\"},\n" +
              "    {\"attribute\":\"run_name\",\"value\":\"" + runName + "\"},\n" +
              "    {\"attribute\":\"run_date\",\"value\":\"" + runDate + "\"},\n" +
              "    {\"attribute\":\"sequencing_platform\",\"value\":\"" + samples[sampleName].attribute2value["MachineType"] + "\"},\n" +
              "    {\"attribute\":\"sequencing_application_type\",\"value\":\"" + samples[sampleName].attribute2value["Application"] + "\"},\n" +
              "    {\"attribute\":\"read_length\",\"value\":\"" + samples[sampleName].attribute2value["ReadLength"] + "\"},\n" +
              "    {\"attribute\":\"pooling\",\"value\":\"" + samples[sampleName].attribute2value["Pooling"] + "\"}\n" +
              "    ]\n}"
    )
    OUT.close()

def createSampleJson(samplePath, sampleName, **samples):
    OUT = open(samplePath + ".metadata.json", "w")
    OUT.write("{\"metadataEntries\":[\n" +
              "    {\"attribute\":\"collection_type\",\"value\":\"Sample\"},\n" +
              "    {\"attribute\":\"sample_id\",\"value\":\"" + samples[sampleName].attribute2value["SampleID"]  + "\"},\n" +
              "    {\"attribute\":\"sample_name\",\"value\":\"" + samples[sampleName].attribute2value["SampleName"]  + "\"},\n" +
              "    {\"attribute\":\"initial_sample_concentration_ngul\",\"value\":\"" + samples[sampleName].attribute2value["SubmittedSampleConc"]  + "\"},\n" +
              "    {\"attribute\":\"initial_sample_volume_ul\",\"value\":\"" + samples[sampleName].attribute2value["SubmittedSampleVolumeUL"]  + "\"},\n" +
              "    {\"attribute\":\"library_id\",\"value\":\"Unspecified\"},\n" +
              "    {\"attribute\":\"library_lot\",\"value\":\"Unspecified\"},\n" +
              "    {\"attribute\":\"sfqc_library_concentration_nM\",\"value\":\"" + samples[sampleName].attribute2value["LibraryConc"]  + "\"},\n" +
              "    {\"attribute\":\"sfqc_library_size\",\"value\":\"" + samples[sampleName].attribute2value["ConsensusLibrarySize"]  + "\"},\n" +
              "    {\"attribute\":\"source_id\",\"value\":\"Unspecified\"},\n" +
              "    {\"attribute\":\"source_name\",\"value\":\"Unspecified\"},\n" +
              "    {\"attribute\":\"source_organism\",\"value\":\"" + samples[sampleName].attribute2value["ReferenceGenome"]  + "\"},\n" +
              "    {\"attribute\":\"source_provider\",\"value\":\"Unspecified\"}\n    ]\n}"
    )
    OUT.close()

def createObjectJson(samplePath, fileName, refAnnotation, sampleName, **samples):
    fileType = "Unspecified"
    dataCompressionStatus = "Not Compressed"
    if fileName[-3:] == "bam":
        fileType = "BAM"
        dataCompressionStatus = "Compressed"
    elif fileName[-3:] == "bai":
        fileType = "INDEX"
        dataCompressionStatus = "Not Compressed"
    elif fileName[-3:] == "tab":
        fileType = "TAB"
        dataCompressionStatus = "Not Compressed"
    elif fileName[-3:] == "csv":
        fileType = "CSV"
        dataCompressionStatus = "Not Compressed"
    elif fileName[-3:] == "log":
        fileType = "LOG"
        dataCompressionStatus = "Not Compressed"
    elif fileName[-3:] == "tbi":
        fileType = "TBI"
        dataCompressionStatus = "Not Compressed"
    elif "vcf.gz" in fileName:
        fileType = "VCF"
        dataCompressionStatus = "Compressed"
    elif "fastq.gz" in fileName:
        fileType = "FASTQ"
        dataCompressionStatus = "Compressed"
    OUT = open(samplePath + "/" + fileName + ".metadata.json", "w")
    OUT.write("{\"metadataEntries\":[\n" +
              "    {\"attribute\":\"object_name\",\"value\":\"" + fileName + "\"},\n" +
              "    {\"attribute\":\"file_type\",\"value\":\"" + fileType + "\"},\n" +
              "    {\"attribute\":\"reference_genome\",\"value\":\"" + samples[sampleName].attribute2value["ReferenceGenome"]  + "\"},\n" +
              "    {\"attribute\":\"reference_annotation\",\"value\":\"" + refAnnotation + "\"},\n" +
              "    {\"attribute\":\"software_tool\",\"value\":\"" + samples[sampleName].attribute2value["Software"]  + "\"},\n" +
              "    {\"attribute\":\"data_compression_status\",\"value\":\"" + dataCompressionStatus + "\"}\n    ]\n}"
    )
    OUT.close()

def configWorkingDirectory(flowcellID, sampleName, **samples):
    pathHPCDMECLU='/mnt/ccrsf-ifx/Software/tools/HPC_DME_APIs'
    pathPythonScripts='/mnt/ccrsf-ifx/Software/scripts/bin'
    homePath = os.getcwd()
    PIName = samples[sampleName].attribute2value["PrincipalInvestigator"]
    PINamePath = ''
    oPIName = re.compile('(.+)(\s{1})([A-Za-z\-]+)')
    if PIName == 'CCRSF':
        PIName = PIName
    elif "," not in PIName:
        sPIName = oPIName.search(PIName)
        if sPIName:
            firstName = re.sub("\s", "", sPIName.group(1))
            PIName = f'{sPIName.group(3)}, {firstName}'
            PINamePath = firstName + "_" + sPIName.group(3)
        else:
            sys.stderr.write("The PI name, " + samples[sampleName].attribute2value["PrincipalInvestigator"] +
                             ", in metadata.txt is not correct.\n")
            sys.exit(1)
    else:
        oListedPIName = re.compile('^([A-Za-z\-\s]+), ([A-Za-z\-\s\.]+)$')
        sListedPIName = oListedPIName.search(PIName)
        if sListedPIName:
            PIName = sListedPIName.group(0)
            PINamePath = f'{re.sub(" ", "", sListedPIName.group(2))}_{re.sub(" ", "", sListedPIName.group(1))}'
        else:
            sys.stdout.write(f"There are mulitple PI names in metadata.txt. The first one is used in the archive path.\n")
            names = PIName.split(',')
            sPIName = oPIName.search(names[0])
            if sPIName:
                firstName = re.sub("\s", "", sPIName.group(1))
                PIName = f'{sPIName.group(3)}, {firstName}'
                PINamePath = firstName + "_" + sPIName.group(3)
            else:
                sys.stderr.write("The PI name, " + samples[sampleName].attribute2value["PrincipalInvestigator"] +
                                 ", in metadata.txt is not correct.\n")
                sys.exit(1)
    path = "PI_Lab_" + PINamePath
    if not os.path.exists(path):
        os.mkdir(path)
    createPILabJson(path, PIName)

    OUT = open(args.foutput, "w")
    refAnnotation = "10X Genomics"
    path += "/Project_" + samples[sampleName].attribute2value["ProjectName"]
    if not os.path.exists(path):
        os.mkdir(path)
    createProjectJson(path, sampleName, **samples)

    fastqPaths = args.fastq_path.split(",")
    runNames = args.run_name.split(",")
    flowcellIDs = []
    i = 0
    tarCmd = ""
    for fastqPath in fastqPaths:
        #180405_J00170_0093_AHV2CKBBXX; NextSeq 2000: 200903_VH00271_2_AAAGNJCM5
        if fastqPath[-1] == '/':
            fastqPath = fastqPath[0:-1]
        elements = runNames[i].split("_")
        if elements[1].startswith('VH'):
            flowcellID = elements[-1] 
        else:    
            if '-' in elements[-1]:
                flowcellID = elements[-1].split('-')[-1]
            else:
                flowcellID = elements[-1][1:] 
        flowcellIDs.append(flowcellID)
        flowcellPath = path + "/Flowcell_" + flowcellID
        if not os.path.exists(flowcellPath):
            os.mkdir(flowcellPath)
        createFlowcellJson(flowcellPath, flowcellID, runNames[i], sampleName, **samples)
        if os.path.exists(fastqPath + "/" + flowcellID):
            for folderName in os.listdir(fastqPath + "/" + flowcellID):
                modifiedFolderName = folderName
                if "__" in folderName:
                    modifiedFolderName = folderName.replace("__", "_")
                else:
                    modifiedFolderName = folderName
                if os.path.isfile(f'{fastqPath}/{flowcellID}/{folderName}'):
                    oFastqFileName = re.compile(r'(.+)_S(\d+)_L(\d+)_[IR](\d)_001.fastq.gz')
                    sFastqFileName = oFastqFileName.search(modifiedFolderName)
                    if sFastqFileName:
                        if not os.path.exists(f'{flowcellPath}/Sample_{sFastqFileName.group(1)}'):
                            os.mkdir(flowcellPath + "/Sample_" + sFastqFileName.group(1))
                            createSampleJson(flowcellPath + "/Sample_" + sFastqFileName.group(1), sFastqFileName.group(1), **samples)
                        cmd = f'ln -s {fastqPath}/{flowcellID}/{modifiedFolderName} {flowcellPath}/Sample_{sFastqFileName.group(1)}/{modifiedFolderName}'
                        runCommand(cmd)
                        createObjectJson(f'{flowcellPath}/Sample_{sFastqFileName.group(1)}', modifiedFolderName, refAnnotation, sFastqFileName.group(1), **samples)
                        OUT.write(f'{flowcellPath}/Sample_{sFastqFileName.group(1)}/{modifiedFolderName}\n')
                else:
                    if not os.path.exists(flowcellPath + "/Sample_" + modifiedFolderName):
                        os.mkdir(flowcellPath + "/Sample_" + modifiedFolderName)
                        createSampleJson(flowcellPath + "/Sample_" + modifiedFolderName, modifiedFolderName, **samples)
                    for fileName in os.listdir(fastqPath + "/" + flowcellID + "/" + folderName):
                        if "__" in fileName:
                            modifiedFileName = fileName.replace("__", "_")
                        else:
                            modifiedFileName = fileName
                        cmd = ("ln -s " + fastqPath + "/" + flowcellID + "/" + folderName +
                               "/" + fileName + " " + flowcellPath + "/Sample_" +
                               modifiedFolderName + "/" + modifiedFileName)
                        runCommand(cmd)
                        createObjectJson(flowcellPath + "/Sample_" + modifiedFolderName, modifiedFileName, refAnnotation, modifiedFolderName, **samples)
                        OUT.write(flowcellPath + "/Sample_" + modifiedFolderName + "/" + modifiedFileName + "\n")
        elif os.path.exists(fastqPath + "/" + samples[sampleName].attribute2value["ProjectName"]):
            for folderName in os.listdir(fastqPath + "/" + samples[sampleName].attribute2value["ProjectName"]):
                if "__" in folderName:
                    modifiedFolderName = folderName.replace("__", "_")
                else:
                    modifiedFolderName = folderName
                if not os.path.exists(flowcellPath + "/" + modifiedFolderName):
                    os.mkdir(flowcellPath + "/" + modifiedFolderName)
                    createSampleJson(flowcellPath + "/" + modifiedFolderName, modifiedFolderName.split('Sample_')[1], **samples)
                for fileName in os.listdir(fastqPath + "/" + samples[sampleName].attribute2value["ProjectName"] + "/" + folderName):
                    if "__" in fileName:
                        modifiedFileName = fileName.replace("__", "_")
                    else:
                        modifiedFileName = fileName
                    cmd = ("ln -s " + fastqPath + "/" + samples[sampleName].attribute2value["ProjectName"] + "/" + folderName +
                           "/" + fileName + " " + flowcellPath + "/" +
                           modifiedFolderName + "/" + modifiedFileName)
                    runCommand(cmd)
                    createObjectJson(flowcellPath + "/" + modifiedFolderName, modifiedFileName, refAnnotation, modifiedFolderName.split('Sample_')[1], **samples)
                    OUT.write(flowcellPath + "/" + modifiedFolderName + "/" + modifiedFileName + "\n")
        else:
            sys.stdout.write(f'No Fastq files are identified in {fastqPath}\n')

        #write metadata.json for undetermined fastq files
        if not os.path.exists(flowcellPath + "/Undetermined"):
            os.mkdir(flowcellPath + "/Undetermined")
        UndeterminedJSON = open(flowcellPath + "/Undetermined.metadata.json", "w")
        UndeterminedJSON.write("{\"metadataEntries\":[\n" +
                               "    {\"attribute\":\"collection_type\",\"value\":\"Sample\"},\n" + 
                               "    {\"attribute\":\"sample_id\",\"value\":\"Unspecified\"},\n" +
                               "    {\"attribute\":\"sample_name\",\"value\":\"Undetermined\"},\n" +
                               "    {\"attribute\":\"source_organism\",\"value\":\"Unspecified\"}\n    ]\n}"
        )
        UndeterminedJSON.close()
        for fileName in os.listdir(fastqPath):
            if "Undetermined" in fileName:
                cmd = ("ln -s " + fastqPath + "/" + fileName + " " + flowcellPath + "/Undetermined/" + fileName)
                runCommand(cmd)
                UndeterminedFastqJSON = open(flowcellPath + "/Undetermined/" + fileName + ".metadata.json", "w")
                UndeterminedFastqJSON.write("{\"metadataEntries\":[\n" +
                              "    {\"attribute\":\"object_name\",\"value\":\"" + fileName + "\"},\n" +
                              "    {\"attribute\":\"file_type\",\"value\":\"FASTQ\"},\n" +
                              "    {\"attribute\":\"software_tool\",\"value\":\"cellranger\"},\n" +
                              "    {\"attribute\":\"data_compression_status\",\"value\":\"Compressed\"}\n    ]\n}"
                )
                UndeterminedFastqJSON.close()
                OUT.write(flowcellPath + "/Undetermined/" + fileName + "\n")

        #archive flowcellID.mri.tgz generated by Cellranger and write metadata.json for flowcellID.mri.tgz
        if os.path.isfile(f'{PurePath(fastqPath).parents[1]}/{flowcellID}.mri.tgz'):
            cmd = (f'ln -s {PurePath(fastqPath).parents[1]}/{flowcellID}.mri.tgz {flowcellPath}/{flowcellID}.mri.tgz')
            runCommand(cmd)
            TARJSON = open(f'{flowcellPath}/{flowcellID}.mri.tgz.metadata.json', 'w')
            TARJSON.write(f'{{"metadataEntries":[\n'
                          f'    {{"attribute":"object_name","value":"{flowcellID}.mri.tgz"}},\n'
                          f'    {{"attribute":"file_type","value":"TGZ"}},\n'
                          f'    {{"attribute":"reference_genome","value":"{samples[sampleName].attribute2value["ReferenceGenome"]}"}},\n'
                          f'    {{"attribute":"reference_annotation","value":"{refAnnotation}"}},\n'
                          f'    {{"attribute":"software_tool","value":"Cellranger"}},\n'
                          f'    {{"attribute":"data_compression_status","value":"Compressed"}}\n    ]\n}}'
            )
            TARJSON.close()
            OUT.write(f'{flowcellPath}/{flowcellID}.mri.tgz\n')
        else:
            sys.stdout.write(f'No {flowcellID}.mri.tgz in {PurePath(fastqPath).parents[1]} is identified to be archived.\n')
            #create run_name_supplement.tar and write metadata.json for run_name_supplement.tar
            sampleSheet = ""
            if os.path.isfile(os.path.dirname(fastqPath) + '/input_samplesheet.csv'):
                sampleSheet = os.path.dirname(fastqPath) + '/input_samplesheet.csv'
            elif os.path.isfile(f'{os.path.dirname(fastqPath)}/{flowcellID}_mod.csv'):
                sampleSheet = f'{os.path.dirname(fastqPath)}/{flowcellID}_mod.csv'
            else:
                sys.stdout.write(f'No sample sheet in {fastqPath} is identified to be archived.\n')
            tarCmd += (f'tar -cvhf {homePath}/{flowcellPath}/{runNames[i]}_supplement.tar {sampleSheet} {fastqPath}/Reports {fastqPath}/Stats --transform="s,{fastqPath[1:]}/,," --transform="s,{os.path.dirname(fastqPath[1:])}/,,"\n')
            TARJSON = open(flowcellPath + "/" + runNames[i] + "_supplement.tar.metadata.json", "w")
            TARJSON.write("{\"metadataEntries\":[\n" +
                          "    {\"attribute\":\"object_name\",\"value\":\"" + runNames[i] + "_supplement.tar\"},\n" +
                          "    {\"attribute\":\"file_type\",\"value\":\"TAR\"},\n" +
                          "    {\"attribute\":\"reference_genome\",\"value\":\"" + samples[sampleName].attribute2value["ReferenceGenome"]  + "\"},\n" +
                          "    {\"attribute\":\"reference_annotation\",\"value\":\"" + refAnnotation + "\"},\n" +
                          "    {\"attribute\":\"software_tool\",\"value\":\"cellranger\"},\n" +
                          "    {\"attribute\":\"data_compression_status\",\"value\":\"Compressed\"}\n    ]\n}"
            )
            TARJSON.close()
            OUT.write(flowcellPath + "/" + runNames[i] + "_supplement.tar\n")
        i += 1

    SLURMOUT = open("dm_register_directory.sh", "w")
    SLURMOUT.write(
        f'#!/bin/bash\n'
        f'#SBATCH --partition=norm\n'
        f'#SBATCH --ntasks=4\n'
        f'#SBATCH --mem=8g\n'
        f'#SBATCH --nodes=1\n'
        f'#SBATCH --time=48:00:00\n'
        f'#SBATCH --no-requeue\n'
        f'{tarCmd}\n'
    )
    if args.count_path:
        if os.path.exists(f'{args.count_path}/libraries.csv'):
            libraryName2cmd = dict()
            with open (f'{args.count_path}/libraries.csv', 'r') as LIBRARIES:
                header = ''
                for line in LIBRARIES:
                    if 'Name,Flowcell,Sample,Type' in line:
                        header = line
                    else:
                        columns = line[:-1].split(',')
                        if columns[0] not in libraryName2cmd:
                            libraryName2cmd[columns[0]] = (f'tar -cvhf {homePath}/{path}/{columns[0]}_count.tar -C {args.count_path} {columns[0]}/outs\n')
                            with open(path + "/" + columns[0] + "_count.tar.metadata.json", "w") as TARJSON:
                                TARJSON.write("{\"metadataEntries\":[\n" +
                                              "    {\"attribute\":\"object_name\",\"value\":\"" + columns[0] + "_count.tar\"},\n" +
                                              "    {\"attribute\":\"file_type\",\"value\":\"TAR\"},\n" +
                                              "    {\"attribute\":\"reference_genome\",\"value\":\"" + samples[sampleName].attribute2value["ReferenceGenome"]  + "\"},\n" +
                                              "    {\"attribute\":\"reference_annotation\",\"value\":\"" + refAnnotation + "\"},\n" +
                                              "    {\"attribute\":\"software_tool\",\"value\":\"cellranger\"},\n" +
                                              "    {\"attribute\":\"data_compression_status\",\"value\":\"Compressed\"}\n    ]\n}"
                                )
                            OUT.write(path + "/" + columns[0] + "_count.tar\n")
                            SLURMOUT.write(f'{libraryName2cmd[columns[0]]}\n')
                        else:
                            continue
        else:
            for entryName in os.listdir(args.count_path):
                if "__" in entryName:
                    modifiedEntryName = entryName.replace("__", "_")
                else:
                    modifiedEntryName = entryName
                if modifiedEntryName in samples:
                    # 10x data
                    if os.path.exists(f'{args.count_path}/{entryName}/outs') and os.path.isdir(f'{args.count_path}/{entryName}/outs'):
                        cmd = (f'tar -cvhf {homePath}/{path}/{modifiedEntryName}_count.tar -C {args.count_path} {entryName}/outs')
                        #runCommand(cmd)
                        SLURMOUT.write(f'{cmd}\n')
                        os.chdir(homePath)
                        TARJSON = open(path + "/" + modifiedEntryName + "_count.tar.metadata.json", "w")
                        TARJSON.write("{\"metadataEntries\":[\n" +
                                      "    {\"attribute\":\"object_name\",\"value\":\"" + modifiedEntryName + "_count.tar\"},\n" +
                                      "    {\"attribute\":\"file_type\",\"value\":\"TAR\"},\n" +
                                      "    {\"attribute\":\"reference_genome\",\"value\":\"" + samples[sampleName].attribute2value["ReferenceGenome"]  + "\"},\n" +
                                      "    {\"attribute\":\"reference_annotation\",\"value\":\"" + refAnnotation + "\"},\n" +
                                      "    {\"attribute\":\"software_tool\",\"value\":\"cellranger\"},\n" +
                                      "    {\"attribute\":\"data_compression_status\",\"value\":\"Compressed\"}\n    ]\n}"
                        )
                        TARJSON.close()
                        OUT.write(path + "/" + modifiedEntryName + "_count.tar\n")
                    # PIPseq data 
                    if os.path.exists(f'{args.count_path}/{entryName}/barcodes') and os.path.isdir(f'{args.count_path}/{entryName}/barcodes'):
                        cmd = (f'tar -cvhf {homePath}/{path}/{modifiedEntryName}_pipseeker.tar -C {args.count_path} {entryName}/')
                        #runCommand(cmd)
                        SLURMOUT.write(f'{cmd}\n')
                        os.chdir(homePath)
                        TARJSON = open(path + "/" + modifiedEntryName + "_pipseeker.tar.metadata.json", "w")
                        TARJSON.write("{\"metadataEntries\":[\n" +
                                      "    {\"attribute\":\"object_name\",\"value\":\"" + modifiedEntryName + "_pipseeker.tar\"},\n" +
                                      "    {\"attribute\":\"file_type\",\"value\":\"TAR\"},\n" +
                                      "    {\"attribute\":\"reference_genome\",\"value\":\"" + samples[sampleName].attribute2value["ReferenceGenome"]  + "\"},\n" +
                                      "    {\"attribute\":\"reference_annotation\",\"value\":\"" + refAnnotation + "\"},\n" +
                                      "    {\"attribute\":\"software_tool\",\"value\":\"cellranger\"},\n" +
                                      "    {\"attribute\":\"data_compression_status\",\"value\":\"Compressed\"}\n    ]\n}"
                        )
                        TARJSON.close()
                        OUT.write(path + "/" + modifiedEntryName + "_pipseeker.tar\n")
            
    if args.aggr_path:
        aggrPaths = args.aggr_path.split(",")
        for aggrPath in aggrPaths:
            aggrPath = aggrPath.rstrip("/")
            elements = aggrPath.split("/")
            #os.chdir(aggrPath + "/..")
            #cmd = (f'tar -cvf {homePath}/{path}/{elements[-1]}.tar {args.aggr_path}/outs')
            cmd = (f'tar -cvhf {homePath}/{path}/{elements[-1]}.tar -C {"/".join(elements[:-1])} {elements[-1]}/outs')
            #runCommand(cmd)
            SLURMOUT.write(f'{cmd}\n')
            #os.chdir(homePath)
            TARJSON = open(path + "/" + elements[-1] + ".tar.metadata.json", "w")
            TARJSON.write("{\"metadataEntries\":[\n" +
                          "    {\"attribute\":\"object_name\",\"value\":\"" + elements[-1] + ".tar\"},\n" +
                          "    {\"attribute\":\"file_type\",\"value\":\"TAR\"},\n" +
                          "    {\"attribute\":\"reference_genome\",\"value\":\"" + samples[sampleName].attribute2value["ReferenceGenome"]  + "\"},\n" +
                          "    {\"attribute\":\"reference_annotation\",\"value\":\"" + refAnnotation + "\"},\n" +
                          "    {\"attribute\":\"software_tool\",\"value\":\"cellranger\"},\n" +
                          "    {\"attribute\":\"data_compression_status\",\"value\":\"Compressed\"}\n    ]\n}"
            )
            TARJSON.close()
            OUT.write(path + "/" + elements[-1] + ".tar\n")
    OUT.close()

    SLURMOUT.write(f'module load java/1.8.0\n'
                   f'export HPC_DM_UTILS={pathHPCDMECLU}/utils\n'
                   f'source $HPC_DM_UTILS/functions\n'
                   f'dm_register_collection PI_Lab_{PINamePath}.metadata.json /FNL_SF_Archive/PI_Lab_{PINamePath}\n'
                   f'dm_register_collection PI_Lab_{PINamePath}/Project_{samples[sampleName].attribute2value["ProjectName"]}.metadata.json /FNL_SF_Archive/PI_Lab_{PINamePath}/Project_{samples[sampleName].attribute2value["ProjectName"]}\n')

    for flowcellID in flowcellIDs:
        SLURMOUT.write(f'dm_register_collection PI_Lab_{PINamePath}/Project_{samples[sampleName].attribute2value["ProjectName"]}/Flowcell_{flowcellID}.metadata.json /FNL_SF_Archive/PI_Lab_{PINamePath}/Project_{samples[sampleName].attribute2value["ProjectName"]}/Flowcell_{flowcellID}\n')

    SLURMOUT.write(f'python {pathPythonScripts}/meta2json_md5.py -p 4 -l {args.foutput} --rewrite\n')
    SLURMOUT.write(f'dm_register_directory -s -l {args.foutput} . /FNL_SF_Archive 1> dm_register_directory.log 2> dm_register_directory.err\n')
    fout_dme_chk = re.sub(".txt", "_dme_check.txt", args.foutput)
    SLURMOUT.write(f'dme_check.py {args.foutput} > {fout_dme_chk}\n')
    SLURMOUT.close()
    sys.stdout.write("sbatch dm_register_directory.sh\n")

def checkMetadata(flowcellID, runName, readLengths, instrument, fastqPath):
    samples = dict()
    sampleName = ""
    projectName = ""

    if args.fmetadata:
        with open(args.fmetadata, "r") as IN:
            line = IN.readline()
            headers = line[0:-1].split("\t")
            for line in IN:
                columns = line[0:-1].split("\t")
                columns[6] = columns[6].replace("__", "_")
                samples[columns[6]] = classSample(headers, columns)
            #sampleName = columns[6] #the last sample name in metadata.txt
            for sampleName in samples:
                if len(readLengths) == 1:
                    samples[sampleName].attribute2value["ReadLength"] = f'1x{readLengths[0]}'
                elif len(readLengths) == 2:
                    if readLengths[0] == readLengths[1]:
                        samples[sampleName].attribute2value["ReadLength"] = f'2x{readLengths[0]}'
                    else:
                        samples[sampleName].attribute2value["ReadLength"] = f'R1: {readLengths[0]}, R2: {readLengths[1]}'
                else:
                    samples[sampleName].attribute2value["ReadLength"] = "Customized" 
            configAttr2Value = {}
            if args.count_path and os.path.exists(f'{args.count_path}/config.py'):
                with open(f'{args.count_path}/config.py', 'r') as CONFIG:
                    line = CONFIG.readline()
                    for line in CONFIG:
                        columns = line[0:-1].split("=")
                        configAttr2Value[columns[0]] = columns[1].replace('"', '')
                    for sampleName in samples:
                        samples[sampleName].attribute2value["ReferenceGenome"] = configAttr2Value['ref']
            else:
                configAttr2Value['ref'] = samples[sampleName].attribute2value["ReferenceGenome"]
                sys.stdout.write(f'{samples[sampleName].attribute2value["ProjectName"]}/config.py does not exist, or --count_path is not assgined.\n'
                                 f'The reference of {configAttr2Value["ref"]} will be retrieved from LabQC metadata...\n')
        configWorkingDirectory(flowcellID, sampleName, **samples)
    else:
        sys.stderr.write(f'No metadata is provided. Trying to generate a template metadata table...\n')
        if args.fprojectName:
            projectName = args.fprojectName
        else:
            oSampleSheetName = re.compile(r'(\d{8})_([-A-Z0-9]+).csv')
            demultiplexPath = f'/scratch/ccrsf_scratch/scratch/Illumina_Demultiplex/{instrument}/{runName}'
            for fileName in os.listdir(demultiplexPath):
                sSampleSheetName = oSampleSheetName.search(fileName)
                if sSampleSheetName:
                    sys.stdout.write(f'No project name is assigned. Trying to retrieve from {demultiplexPath}/{fileName}\n')
                    with open (f'{demultiplexPath}/{fileName}', 'r') as fSampleSheet:
                        for line in fSampleSheet:
                            if line.startswith('[Data]'):
                                break
                        header = fSampleSheet.readline()
                        columns = header.split(',')
                        index = 0
                        columnIndexProject = 0
                        for column in columns:
                            if 'Project' in column:
                                columnIndexProject = index
                            index += 1
                        line = fSampleSheet.readline()
                        projectName = line[:-1].split(',')[columnIndexProject]
                    sys.stdout.write(f'{projectName} is used for uploading\n')
        if projectName:
            with open(f'{projectName}_{flowcellID}_Metadata.txt', 'w') as OUT:
                OUT.write(
                    f'SampleID\tPrincipalInvestigator\tLabContact\tProjectName\t'
                    f'BioinformaticsContact\tProjectOpenDate\tSampleName\t'
                    f'SampleType\tApplication\tReferenceGenome\tCSAS_NAS\t'
                    f'MachineType\tReadLength\tPooling\tLibraryKit\t'
                    f'SubmittedSampleVolumeUL\tSubmittedSampleConc\t'
                    f'SubmittedSampleConcUnit\tDateReceived\tLabQCSampleConc\t'
                    f'LabQCSampleConcUnit\tRIN\t28s18s\tRnaArea\tSamplePeakSize\t'
                    f'SamplePeakSizeUnit\tSampleRegionSize\tSampleRegionSizeUnit\t'
                    f'LibraryConc\tLibraryConcUnit\tConsensusLibrarySize\t'
                    f'ConsensusLibrarySizeUnit\tLibraryPeakSize\t'
                    f'LibraryPeakSizeUnit\tLibraryRegionSize\t'
                    f'LibraryRegionSizeUnit\tSoftware\n'
                )
                sys.stderr.write(f'Checking sample names in {fastqPath}...\n')
                if fastqPath:
                    oSampleName = re.compile('Sample_(.+)')
                    oFastqName = re.compile("(.+)_R1_001.fastq.gz")
                    for dirName in os.listdir(fastqPath):
                        if os.path.isdir(f'{fastqPath}/{dirName}') and dirName == projectName:
                            for dirSampleName in os.listdir(f'{fastqPath}/{dirName}'):
                                sSampleName = oSampleName.search(dirSampleName)
                                sFastqName = oFastqName.search(dirSampleName)
                                if sSampleName:
                                    sampleName = sSampleName.group(1)
                                elif sFastqName:
                                    oSampleName = re.compile('(.+)_S(\d+)_R1_001.fastq.gz')
                                    sSampleName = oSampleName.search(dirSampleName)
                                    if sSampleName:
                                        sampleName = sSampleName.group(1)
                                    else:
                                        sampleName = dirSampleName.replace("__", "_").split('_R1_001.fastq.gz')[0]
                                else:
                                    sampleName = dirName
                                tabs = '\t' * 30
                                OUT.write(f'\tCCRSF\t\t{projectName}\t\t\t{sampleName}{tabs}\n')
                        elif os.path.isdir(f'{fastqPath}/{dirName}') and dirName != projectName:
                            sys.stderr.write(f'{fastqPath}/{dirName} is not consistent with {projectName}. Skipped...\n')
                            continue
                        else:
                            sys.stderr.write(f'{fastqPath}/{dirName} is Skipped...\n')
                            continue
                    sys.stderr.write(
                        f'A template metadata table, {projectName}_{flowcellID}_Metadata.txt'
                        f', has been genreated. Pleaase edit it accordingly and include it '
                        f'in the command line, e.g.,\n'
                        f'meta2json.py -m {projectName}_{flowcellID}_Metadata.txt -uf {args.fastq_path} -r {args.run_name}\n'
                    ) 
                else:   
                    sys.stderr.write(
                        f'The unaligned folder is not defined and no metadata is '
                        f'avaiable. No archive will be done.\n'
                    )
                    sys.exit(1)
        else:
            sys.stderr.write(f'No project name is identified for uploading.\nPlease assign the project name using -p in the command line.\n')
            sys.exit(1)

def getRunPath(sRunName):
    instrument = ''
    runPath = ''
    rawDataPathIslon2 = "/is2/projects/CCR-SF/scratch/illumina"
    rawDataPathQumulo = "/mnt/ccrsf-raw/illumina"
    if sRunName.group(2)[0] == "A":
        instrument = "NovaSeq 6000"
        if os.path.isdir(os.path.join(rawDataPathIslon2, "RawData", sRunName.group(0))):
            runPath = os.path.join(rawDataPathIslon2, "RawData")
        elif os.path.isdir(os.path.join(rawDataPathQumulo, "RawData_NovaSeq", sRunName.group(0))):
            runPath = os.path.join(rawDataPathQumulo, "RawData_NovaSeq")
    elif sRunName.group(2)[0] == "L":
        instrument = "NovaSeq Xplus"
        if os.path.isdir(os.path.join('/is2/projects/CCR-SF/scratch/illumina', "RawData", sRunName.group(0))):
            runPath = os.path.join('/is2/projects/CCR-SF/scratch/illumina', "RawData")
        elif os.path.isdir(os.path.join('/mnt/ccrsf-raw/illumina', "RawData_Xplus", sRunName.group(0))):
            runPath = os.path.join('/mnt/ccrsf-raw/illumina', "RawData_Xplus")
    elif sRunName.group(2)[0] == "J" or sRunName.group(2)[0] == "D":
        instrument = "HiSeq"
        if os.path.isdir(os.path.join(rawDataPathIslon2, "RawData", sRunName.group(0))):
            runPath = os.path.join(rawDataPathIslon2, "RawData")
        elif os.path.isdir(os.path.join(rawDataPathQumulo, "RawData_HiSeq", sRunName.group(0))):
            runPath = os.path.join(rawDataPathQumulo, "RawData_HiSeq")
    elif sRunName.group(2)[0] == "N":
        instrument = "NextSeq 550"
        if os.path.isdir(os.path.join(rawDataPathIslon2, "RawData_NextSeq", sRunName.group(0))):
            runPath = os.path.join(rawDataPathIslon2, "RawData_NextSeq")
        elif os.path.isdir(os.path.join(rawDataPathQumulo, "RawData_NextSeq", sRunName.group(0))):
            runPath = os.path.join(rawDataPathQumulo, "RawData_NextSeq")
    elif sRunName.group(2)[0] == "V":
        instrument = "NextSeq 2000"
        if os.path.isdir(os.path.join(rawDataPathIslon2, "RawData_NextSeq", sRunName.group(0))):
            runPath = os.path.join(rawDataPathIslon2, "RawData_NextSeq")
        elif os.path.isdir(os.path.join(rawDataPathQumulo, "RawData_NextSeq", sRunName.group(0))):
            runPath = os.path.join(rawDataPathQumulo, "RawData_NextSeq")
    elif sRunName.group(2)[0] == "M":
        instrument = "MiSeq"
        if os.path.isdir(os.path.join(rawDataPathIslon2, "RawData_MiSeq", sRunName.group(0))):
            runPath = os.path.join(rawDataPathIslon2, "RawData_MiSeq")
        elif os.path.isdir(os.path.join(rawDataPathQumulo, "RawData_MiSeq", sRunName.group(0))):
            runPath = os.path.join(rawDataPathQumulo, "RawData_MiSeq")
    elif sRunName.group(2)[0] == "F":
        instrument = "iSeq"
        if os.path.isdir(os.path.join(rawDataPathIslon2, "RawData_iSeq", sRunName.group(0))):
            runPath = os.path.join(rawDataPathIslon2, "RawData_iSeq")
        elif os.path.isdir(os.path.join(rawDataPathQumulo, "RawData_iSeq", sRunName.group(0))):
            runPath = os.path.join(rawDataPathQumulo, "RawData_iSeq")
    return(instrument, runPath)
    
def main(argv):
    samples = dict()
    path = os.getcwd()
    paths = path.split("/")
    if args.run_name is None:
        runName = paths[-1]
    else:
        runName = args.run_name.split(',')[0]
        fastqPath = args.fastq_path.split(',')[0]
    instrument = "Unknown"
    runPath = ''
    #rawDataPath = "/mnt/ccrsf-raw/illumina" 
    xmlRunParametersPath = ""
    xmlRunInfoPath = ""
    readLengths = [] 
    indexLengths = []
    flowcell = "Unknown"
    workFlowType = "Unknown"
    flowcellMode = "Unknown"
    chemistry = "Unknown"
    chemistryVersion = "Unknown"
    oRunName = re.compile("(\d{6,})_([A-Z0-9]+)_(\d+)_([-A-Z0-9]+)")
    sRunName = oRunName.search(runName)
    if sRunName:
        (instrument, runPath) = getRunPath(sRunName)
        if os.path.isdir(os.path.join(runPath, runName)):
            xmlRunParametersPath = glob.glob(os.path.join(runPath, runName, "[rR]unParameters.xml" ))[0]
            xmlRunInfoPath = glob.glob(os.path.join(runPath, runName, "[rR]unInfo.xml" ))[0]
            (RTAVersion, flowcell, workFlowType, flowcellMode, chemistry, chemistryVersion) = parseXml(xmlRunParametersPath, xmlRunInfoPath, instrument, readLengths, indexLengths)
            if flowcellMode == 'NextSeq 2000 P3 Flow Cell Cartridge':
                flowcellMode = 'P3'
            elif flowcellMode == 'NextSeq 1000/2000 P2 Flow Cell Cartridge':
                flowcellMode = 'P2'
        else:
            sys.stdout.write(f'{xmlRunParametersPath}/{xmlRunInfoPath} does not exist.\nRetriving read lengths from fastq.gz...\n')
            for fileName in os.listdir(fastqPath):
                if fileName == 'Undetermined_S0_L001_R1_001.fastq.gz':
                    readLengths.append(retrieveReadLength(f'{fastqPath}/{fileName}'))
                elif fileName == 'Undetermined_S0_L001_R2_001.fastq.gz':
                    readLengths.append(retrieveReadLength(f'{fastqPath}/{fileName}'))
        checkMetadata(flowcell, runName, readLengths, instrument, fastqPath)
    else:
        sys.stdout.write(
            f'The run name is not recognized from the current path.\n'
            f'Please assign the run name using,\n'
            f'meta2json_single_cell.py -m metadata.txt -r run_name\n'
        )
        sys.exit(1)

if __name__ == "__main__":
    main(sys.argv[1:])
