rule clonevol:
    input: 
        "results/pyclone-vi_filtered/{sample}.pyclone-vi.filt.tsv"
    params:
        model = config["clonevol_model"]
    output:
        ccf_plot = "results/clonevol/plots/{sample}/{sample}.ccf_jitterplot.pdf",
        model="results/clonevol/{sample}.pyclone-vi.filt.clonevol.rda"
    log: 
        stdout="logs/clonevol.{sample}.stdout", 
        stderr="logs/clonevol.{sample}.stderr" 
    conda:
        "../envs/clonevol.yml"
    shell: 
        """
        Rscript workflow/scripts/clonevol.R {input} {params.model} {output.ccf_plot} {output.model}\
        >> {log.stdout} 2> {log.stderr}
        """