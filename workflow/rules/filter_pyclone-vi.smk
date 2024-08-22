rule filter_pyclone:
    input:
        "results/pyclone-vi_tsv/{sample}.pyclone-vi.tsv"
    params:
        min_cluster_size = config["min_cluster_size"],
        min_founder_size = config["min_founder_size"]
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
        --out {output}\
        >> {log.stdout} 2> {log.stderr}
        """
        
        python3 workflow/scripts/filter_pyclone-vi.py\
        --input results/pyclone-vi_tsv/T-3001.pyclone-vi.lala.tsv \
        --min_cluster_size 100 \
        --min_founder_size 0.1 \
        --out testo.tsv \
        >> testo.out 2> testo.err