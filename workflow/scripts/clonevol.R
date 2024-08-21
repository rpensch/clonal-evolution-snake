library(clonevol, lib.loc = "workflow/envs/R_packages/")
library(trees, lib.loc = "workflow/envs/R_packages/")

args <- commandArgs(trailingOnly = TRUE)

pyclone = read.table(args[1], header = TRUE)
prefix = sub(".pyclone-vi.filt.tsv","",sub("results/pyclone-vi_filtered/", "",args[1]))

# Prep data for clonevol 
pyclone$cellular_prevalence = pyclone$cellular_prevalence*100
pyclone$cluster = pyclone$cluster_id

samples = unique(pyclone$sample_id)
if (length(samples) == 1) {

    # If only one sample is present:
    # Because of a bug, Clonevol does not run with only one sample 
    # -> create an artificial second sample with the same values
    merged = merge(x = pyclone[pyclone$sample_id == samples[1],], 
                   y = pyclone[pyclone$sample_id == samples[1],c('mutation_id','cellular_prevalence')], 
                   by = "mutation_id", all = TRUE, 
                   suffixes = c(paste("_", sample, sep = ''),
                                paste("_art_", sample, sep = '')))

    vaf.col.names = c(paste("cellular_prevalence_", sample, sep = ''),
                      paste("cellular_prevalence_art_", sample, sep = ''))

} else {

    # If multiple samples:
    # Reformat the df so that the CCF of every sample is in a unique column
    merged = pyclone[pyclone$sample_id == samples[1],]
    vaf.col.names =c(paste('cellular_prevalence', samples[1], sep = '_'))

    for (sample in samples[2:length(samples)]) {

        merged = merge(x = merged, 
                   y = pyclone[pyclone$sample_id == sample,c('mutation_id','cellular_prevalence')], 
                   by = "mutation_id", all = TRUE, suffixes = c("",paste("_", sample, sep = '')))
        vaf.col.names = append(vaf.col.names, paste('cellular_prevalence', sample, sep = '_'))

    }

    names(merged)[4] = paste('cellular_prevalence', samples[1], sep = '_')
    
}

# More formatting prep for clonevol
sample.names = unique(pyclone$sample_id)
merged[, sample.names] <- merged[, vaf.col.names]
vaf.col.names <- sample.names
sample.groups <- sample.names
names(sample.groups) <- vaf.col.names

# Setup the order of clusters to display in various plots (later)
merged <- merged[order(merged$cluster),]

# Plot CCF
pdf(args[3], width = 5, height = 3, useDingbats = FALSE, title='')
pp <- plot.variant.clusters(merged,
                            cluster.col.name = 'cluster',
                            show.cluster.size = FALSE,
                            cluster.size.text.color = 'blue',
                            vaf.col.names = vaf.col.names,
                            vaf.limits = 100,
                            violin = FALSE,
                            box = FALSE,
                            jitter = TRUE,
                            jitter.shape = 1,
                            jitter.color = NULL,
                            jitter.size = 2,
                            jitter.alpha = 1,
                            jitter.center.method = 'mean',
                            jitter.center.size = 1,
                            jitter.center.color = 'darkgray',
                            jitter.center.display.value = 'median',
                            highlight = 'is.driver',
                            highlight.shape = 21,
                            highlight.color = 'blue',
                            highlight.fill.color = 'green',
                            highlight.note.col.name = 'gene',
                            highlight.note.size = 2,
                            order.by.total.vaf = TRUE, 
                            show.cluster.axis.label = FALSE,
                            base_size = 12)
dev.off()

# Run clonevol
y = infer.clonal.models(variants = merged,
                        cluster.col.name = 'cluster',
                        sample.groups = sample.groups,
                        ccf.col.names = vaf.col.names,
                        cancer.initiation.model=args[2],
                        subclonal.test = 'bootstrap',
                        subclonal.test.model = 'non-parametric',
                        num.boots = 1000,
                        founding.cluster = 1,
                        cluster.center = 'mean',
                        ignore.clusters = NULL,
                        clone.colors = NULL,
                        min.cluster.vaf = 0.01,
                        # min probability that CCF(clone) is non-negative
                        sum.p = 0.01,
                        # alpha level in confidence interval estimate for CCF(clone)
                        alpha = 0.01, 
                        vaf.in.percent= TRUE)

y <- convert.consensus.tree.clone.to.branch(y, branch.scale = 'sqrt')

save(y, file = args[4])

plot.clonal.models(y,
                   # box plot parameters
                   box.plot = FALSE,
                   # bell plot parameters
                   clone.shape = 'bell',
                   bell.event = TRUE,
                   bell.event.label.color = 'blue',
                   bell.event.label.angle = 60,
                   clone.time.step.scale = 1,
                   bell.curve.step = 2,
                   # node-based consensus tree parameters
                   merged.tree.plot = FALSE,
                   # branch-based consensus tree parameters
                   merged.tree.clone.as.branch = TRUE,
                   mtcab.event.sep.char = ',',
                   mtcab.branch.text.size = 1,
                   mtcab.branch.width = 1,
                   mtcab.node.size = 3,
                   mtcab.node.label.size = 1,
                   mtcab.node.text.size = 1,
                   #meta-parameters
                   scale.monoclonal.cell.frac = TRUE,
                   show.score = FALSE,
                   cell.frac.ci = FALSE,
                   disable.cell.frac = FALSE,
                   # output figure parameters
                   out.dir = paste('results/clonevol/plots/',prefix,sep=''),
                   out.format = 'pdf',
                   overwrite.output = TRUE,
                   width = 5,
                   height = 4,
                   # vector of width scales for each panel from left to right
                   panel.widths = c(3,1),
                   max.num.models.to.plot = 10)

#}