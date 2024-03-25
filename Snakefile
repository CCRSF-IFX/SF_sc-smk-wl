import logging as sflog
sflog.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=sflog.DEBUG)

def get_config_val(config_file):
    pipeline, fullanalysis = 'rna', False
    with open(config_file) as fconfig: 
        for line in fconfig:
            if line.startswith('pipeline='):
                pipeline = line.rstrip().split('=')[-1].strip('"')
            if line.startswith('fullanalysis='):
                fullanalysis = line.rstrip().split('=')[-1].strip('"')
    return pipeline, fullanalysis 

pipeline, fullanalysis = get_config_val("config.py")
sflog.debug(pipeline)
sflog.debug(fullanalysis)

if pipeline == "rna":
    if fullanalysis: 
        include: "workflow/Snakefile_rna_fullanalysis"
    else:
        sflog.info("Running rna pipeline...") 
        include: "workflow/Snakefile_rna"
