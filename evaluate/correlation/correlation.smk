rule correlation:
    conda:
        "environment.yaml"
    input:
        script=get_wrapper_path("evaluate/corelation/correlate.py"),
        files=["a.tsv.gz", "b.tsv.gz"],
    output:
        "output.tsv.gz",
    params:
        value_a="ID",
        value_b="ID",
    log:
        "correlation.log",
    shell:
        """
        fileB=""
        files=({input.files})
        if [ "${{#files[@]}}" -eq 2 ]; then
            fileB="--fileB ${{files[1]}}"
        fi
        python  {input.script} \
        --values {params.value_a} {params.value_b} \
        --fileA {input.files[0]}  `echo "${{fileB}}"`  \
        --output {output} &> {log}
        """
