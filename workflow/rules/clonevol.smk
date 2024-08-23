rule clonevol:
    input: 
        "results/pyclone-vi_filtered/{sample}.pyclone-vi.filt.tsv"
    output:
        ccf_plot = "results/clonevol/plots/{sample}.ccf_jitterplot.pdf",
        monoclonal="results/clonevol/{sample}.pyclone-vi.filt.monoclonal.clonevol.rda",
        polyclonal="results/clonevol/{sample}.pyclone-vi.filt.polyclonal.clonevol.rda"
    log: 
        stdout="logs/clonevol.{sample}.stdout", 
        stderr="logs/clonevol.{sample}.stderr" 
    conda:
        "../envs/clonevol.yml"
    shell: 
        """
        Rscript workflow/scripts/clonevol.R {input} {output.ccf_plot} {output.monoclonal} {output.polyclonal}\
        >> {log.stdout} 2> {log.stderr}
        """