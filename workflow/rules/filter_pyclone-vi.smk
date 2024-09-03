rule filter_pyclone:
    input:
        "results/pyclone-vi_tsv/{sample}.pyclone-vi.tsv"
    params:
        min_cluster_size = config["min_cluster_size"],
        min_founder_size = config["min_founder_size"],
        adjust_ccf = config["adjust_ccf"]
    output:
        "results/pyclone-vi_filtered/{sample}.pyclone-vi.filt.tsv"
    log: 
        stdout="logs/filter_pyclone.{sample}.stdout", 
        stderr="logs/filter_pyclone.{sample}.stderr" 
    conda:
        "../envs/filter_pyclone-vi.yml"
    shell:
        """
        python3 workflow/scripts/filter_pyclone-vi.py\
        --input {input}\
        --min_cluster_size {params.min_cluster_size}\
        --min_founder_size {params.min_founder_size}\
        --adjust_ccf {params.adjust_ccf}\
        --out {output}\
        >> {log.stdout} 2> {log.stderr}
        """