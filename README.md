# Lab_GeneExpression

## General Setup

First, hop onto Poseidon and clone this repo.

Request some designated space on the HPC:\
`srun --time=02:00:00 --mem=10gb -n 1 -p compute --pty bash`

Now, let's set up a conda environment to play in using the provided .yml file. It can be complicated to install R in a conda environment, so I'm giving you a file that should let you work with R v4 in a jupyter notebook (the final line helps with that).

```
mamba env create -f R_DE.yml
mamba activate r_diffex
R -e 'IRkernel::installspec()'
```

We're going to install the packages we need directly in R. Open R interactively by typing:\
`R`

Your prompt should change to a greater-than sign:\
`>`

First, we need to install some packages we'll want to work with: `DESeq2` (which is used for differential expression), and `ggplot2`, which is used for plotting. The standard approach to installing R packages uses the R `install.packages()` command. However, some R packages are available on a more specialized sciencey corner of R called Bioconductor (similar to biopython), and these need to be installed somewhat differently. If `install.packages()` says that a package doesn't exist, search for it online - it is likely installed via Bioconductor instead.

Install `DESeq2` through Bioconductor. You'll need to first install `BiocManager` (the Bioconductor installer). For good measure, and so we can play with data display, load `ggplot2` and `gridExtra`.

```
install.packages("BiocManager")
install.packages("ggplot2")
install.packages("gridExtra")
BiocManager::install("DESeq2")
```

Note: In theory, you should be able to run all of this in a Jupyter notebook in an R kernel. In practice, sometimes complex R packages can be a bit difficult to install in Jupyter environments. Today we'll work directly in R on the HPC. Later, we'll cover how to make images in R on the screen-less HPC interface.

## Differential Expression

Nearly all gene expression analysis programs work in R. Much like python, a lot of the good stuff in R is done through "add-on" modules (packages) for more specialized tasks.

*BONUS R factoid for baseball lovers: you can install the entire Sean Lahman Baseball Database in R if you want to try your hand at sabermetrics. The package is called `Lahman`, and it contains statistics from 1871-present.*

Many folks prefer to run R locally (not on the HPC) because it can be finicky to install R packages and deal with images on the HPC. If you want a fuller-featured interface, you can use R on the HPC via Jupyter (should be possible with this class environment), or even RStudio with some legwork. Arianna provides some guidance for setting up RStudio on the HPC [here](https://alexanderlabwhoi.github.io/post/2021-03-17-xquartz/). For simplicity, today we're just going to work directly in R on Poseidon.

Once packages are installed, you need to load them in your environment (or script) in order to use them. R packages are loaded with this syntax: `library(PACKAGE_NAME)`. Like so many things in UNIX, package names are case-sensitive. Now, load the DESeq2 library that you previously installed:

```
library(DESeq2)
library(ggplot2)
library(gridExtra)
```

Now, we'll run through a simple DE analysis, following along with [the very detailed DESeq2 vignette](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html). I've provided you with two files: `Rh_NC_sal_counts_no-contam.tsv` contains raw read counts (number of reads from each sample mapping to each contig); and `Rh_NC_sal_sample-info.csv` contains metadata on samples, including important information on experimental exposures.

First we need to load the samples:

```
cts = as.matrix(read.csv("Rh_NC_sal_counts_no-contam.tsv", sep="\t", row.names = "Contig"))
coldata = read.csv("Rh_NC_sal_sample-info.csv", row.names="Sample")
```

Let's take a look at them:

```
dim(cts)
head(cts)
dim(coldata)
head(coldata)
```

There's a LOT of info in the `coldata` file, and we don't need all of it. Let's cut that down to just the factors we're interested in - Source and Acclimation.

```
coldata = coldata[,c("Source","Acclimation")]
coldata$Source = factor(coldata$Source)
coldata$Acclimation = factor(coldata$Acclimation)
```

Take another look at our new coldata file:

```
dim(coldata)
head(coldata)
```

Make sure the samples are listed in the **same order** in the cts and coldata files. This is important! DESeq2 assumes this is true, and won't run right if sample order differs across the two files.

```
cts = cts[, rownames(coldata)]
```

Like the coldata file, the cts file includes some extra samples and conditions we aren't going to test today. (Samples acclimated at "Fresh" conditions - 0.8 PSU - because this treatment wasn't replicated across all sources.) Remove them, so both data sets include only samples in the "Low" and "Moderate" conditions. 

```
levels(coldata$Acclimation)
cts_trim = cts[,!colnames(cts) %in% c('P00_07','P00_08','P00_01','P00_02')]
keep = coldata$Acclimation != "Fresh"
coldata_trim = coldata[keep,]
```

Check that we have indeed removed the "Fresh" samples.

```
coldata_trim$Acclimation = factor(coldata_trim$Acclimation)
levels(coldata_trim$Acclimation)
```

A little more prep - we need to create a DESeq-formatted data set from our matrix. This will bring together the cts data, the coldata information, and we'll specify the experimental design to be used.

```
dds = DESeqDataSetFromMatrix(countData = cts_trim, colData = coldata_trim, design = ~ Source + Acclimation)
dds
```

Below a certain point, expression is too low for some genes to usefully compared. To ease computation, we'll remove very low-expression genes from the data set now. In this case, we'll specify that if there are <10 reads *in total* across all samples, the gene will be dropped.

```
keep = rowSums(counts(dds)) >= 10
dds = dds[keep,]
dds
```

Now, the main event....we run DESeq2! Lots of stuff happens under the hood, so you'll see several updates as it proceeds (it will take a couple minutes).

```
dds = DESeq(dds)
```

We've actually generated a ton of comparison data, but we need to parse it further to get something usable out of it. A big part of this is multiple test correction - keep in mind that we are testing >137K genes, so we definitely need to compensate for that. Set the alpha level (significance to correct to) to a standard 0.05.

```
res = results(dds, alpha = 0.05)
```

For downstream analysis, save the results as a comma-separated file:

```
write.csv(res, file = "DE_results.csv")
```

Take a look at your results, and use R to count the number of significantly-different genes.

```
summary(res)
resOrder = res[order(res$pvalue),]
head(resOrder)
sum(res$padj < 0.05, na.rm=TRUE)
```

How about the genes with the biggest (not necessarily most significant!) changes?

```
resOrder = res[order(res$log2FoldChange),]
head(resOrder)
resOrder = res[order(-res$log2FoldChange),]
head(resOrder)
```

Now let's think about plotting some of these top counts!

IMPORTANT NOTE: Do NOT directly enter any R commands that output figures on Poseidon (in Jupyter is fine, on the HPC directly is not). Always redirect the output to a file. By default, R prints to screen...and  the HPC doesn’t have a screen. (It will not print to your local computer, since it’s not running locally.) If you do directly enter a command that outputs a figure, you will get nothing more than a sad, empty file named `Rplots.pdf`. If you're using a more complex command with multiple outputs, there are ways to get around this issue, but for now, please just redirect anything you want to look at to pdf using this general approach:

```
pdf(file = "FILENAME.pdf")
COMMAND(S) THAT MAKE(S) THE FIGURE
dev.off()
```

Plot counts for the genes with the biggest adjusted p-values:

```
resOrder = res[order(res$padj),]

pdf(file = "DE_top6.pdf")

plots = list()
for (i in 1:6) {
 name = rownames(resOrder[i,])
 d = plotCounts(dds, gene=name, intgroup=c("Source", "Acclimation"), returnData=TRUE)
 plots[[i]] = ggplot(d, aes(x=Acclimation, y=count, color = Source)) + geom_point(position=position_jitter(w=0.1,h=0), size = 3) + theme_bw() + scale_y_log10() + ggtitle(name)
}
do.call(grid.arrange,plots)

dev.off()
```

Let's save our session info, in case we want to revisit our results without rerunning everything:\
`save.image(file="DE_lab.RData")`

To pick up where we left off:\
`load(file="DE_lab.RData")`
