
rule pyclone:
    input:
        lambda wildcards: input_files.loc[wildcards.sample]["input_file"]
    params:
        c=config["pyclone_c"],
        d=config["pyclone_d"],
        g=config["pyclone_g"],
        r=config["pyclone_r"]
    output:
        h5="results/pyclone-vi_h5/{sample}.pyclone-vi.h5",
        tsv="results/pyclone-vi_tsv/{sample}.pyclone-vi.tsv"
    log: 
        stdout="logs/pyclone_fit.{sample}.stdout", 
        stderr="logs/pyclone_fit.{sample}.stderr"
    conda:
        "../envs/pyclone-vi.yml"
    shell:
        """
        pyclone-vi fit -i {input} -o {output.h5} -c {params.c} -d {params.d} -g {params.g} -r {params.r} \
        > {log.stdout} 2> {log.stderr}
        pyclone-vi write-results-file -i {output.h5} -o {output.tsv}\
        >> {log.stdout} 2>> {log.stderr}
        """
        