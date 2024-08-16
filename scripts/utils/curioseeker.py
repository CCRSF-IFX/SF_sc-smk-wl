import pandas as pd

def meta2json(args, SLURMOUT, OUT, analysis_folder, path, analysis_subfolder):
    libraryName2cmd = {}
    df = pd.read_csv(args.library_file)
    for index, row in df.iterrows():
        if row['sample'] not in libraryName2cmd:
            libraryName2cmd[row['sample']] = (f'tar -cvhf {analysis_folder}/{row['sample']}_curioseeker.tar -C {args.count_path} {row['sample']}\n')
            with open(path + f"/{analysis_subfolder}/" + row['sample'] + "_curioseeker.tar.metadata.json", "w") as TARJSON:
                TARJSON.write("{\"metadataEntries\":[\n" +
                              "    {\"attribute\":\"object_name\",\"value\":\"" + row['sample'] + "_count.tar\"},\n" +
                              "    {\"attribute\":\"file_type\",\"value\":\"TAR\"},\n" +
                              "    {\"attribute\":\"reference_genome\",\"value\":\"" + row['genome'] + "\"},\n" +
                              "    {\"attribute\":\"software_tool\",\"value\":\"cellranger\"},\n" +
                              "    {\"attribute\":\"data_compression_status\",\"value\":\"Compressed\"}\n    ]\n}"
                )
            OUT.write(path + f"/{analysis_subfolder}/" + row['sample'] + "_curioseeker.tar\n")
            SLURMOUT.write(f'{libraryName2cmd[row['sample']]}\n') 
