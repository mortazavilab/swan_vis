# Understanding Swan visualizations

Swan utilizes a representation of transcript structure and alternative splicing that most people aren't familiar with. The goal of this guide is to help the user understand and interpret the graphical output format from Swan, in a step-by-step manner.

## Table of contents

* [Basics](understanding_swan_vis.md#basics)
* [Alternative splice site usage](understanding_swan_vis.md#alternative-splice-site-usage)
* [Exon skipping and intron retention](understanding_swan_vis.md#exon-skipping-and-intron-retention)

## Basics

### 1. Swan pulls transcript models from a GTF or TALON database

![](../.gitbook/assets/teaching_1.png)

### 2. Swan assigns a node to each unique splice site. Blue nodes are transcription start sites, orange nodes are transcription end sites, and yellow nodes are internal splice sites.

![](../.gitbook/assets/teaching_3.png)

### 3. Splice junction pairs that span an exon are connected by exonic \(green\) edges. Splice junction pairs that span an intron are connected by intronic \(pink\) edges.

![](../.gitbook/assets/teaching_8.png)

### 4. Datasets containing additional transcript models can be added.

![](../.gitbook/assets/teaching_9.png)

### 5. New nodes correspond to splice sites not seen in the transcript models already in the SwanGraph.

![](../.gitbook/assets/teaching_10.png)

### 6. Finally, new edges are added to connect the new nodes to the preexisting SwanGraph.

![](../.gitbook/assets/teaching_11.png)

## Alternative splice site usage

Phenomena such as alternative 5'/3', and TSS/TES usage can be visualized from the SwanGraph.

Alternative 5' splice site usage for an exon can be seen when there are multiple incoming exonic \(green\) edges into a splice site that represents the start of an intron.

![](../.gitbook/assets/alt_5.png)

Alternative 3' splice site usage for an exon can be seen when there are multiple outgoing exonic \(green\) edges from a splice site that represents the end of an intron.

![](../.gitbook/assets/alt_3.png)

## Exon skipping and intron retention

You can also visualize exon skipping and intron retention events in Swan graphs.

Exon skipping in a SwanGraph consists of an intronic \(pink\) edge that completely spans an exonic \(green\) edge. This means that an exonic region of a one transcript model has been completely skipped over in another transcript model.

![](../.gitbook/assets/exon_skipping.png)

In a SwanGraph, intron retention is the opposite of exon skipping. It is seen when an exonic \(green\) edge completely spans an intronic \(pink\) edge. This means that an intronic region from one transcript model has been included in a different transcript model.

![](../.gitbook/assets/intron_retention.png)

