from snakemake.workflow import srcdir

WRAPPER_DIR = srcdir("./")


def get_wrapper_path(name):
    return "%s/%s" % (WRAPPER_DIR, name)


include: "evaluate/evaluate.smk"
