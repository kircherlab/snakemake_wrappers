import os

def get_wrapper_path(name):
    return os.path.join(os.path.dirname(workflow.main_snakefile), name)



include: "evaluate/evaluate.smk"
