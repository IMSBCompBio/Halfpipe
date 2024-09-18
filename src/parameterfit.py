import os

def ParameterFit(summaryfile, output, timepoints, threads, model):

    task = f"Rscript src/parameterfit.R -i {summaryfile} -o {output} -tp {' '.join(map(str, timepoints))} -c {threads} -m {model}"
    os.system(task)
    return None