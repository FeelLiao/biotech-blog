---
title: 差异分析三巨头 DESeq2、edgeR、limma
date: 2024-06-01 17:25:03
categories:
- 转录组分析
tags:
- 生物信息
- 转录组分析
- 生信软件使用
- R语言
---

## 差异分析三巨头 DESeq2、edgeR、limma

DESeq2、edgeR、limma是转录组差异分析的金标准，大多数转录组的文章和公司都是使用这三个R包做转录组的差异基因分析。

做差异基因表达分析需要的数据有：**表达矩阵、分组信息**

1. **表达矩阵**：即上游分析得到的每个基因在每个样本中的reads数，这里所使用的差异分析包都要求原始reads数的格式输入。
2. **分组信息**：即一个实验设计矩阵，包含每一个测序样本属于哪一个实验组或者对照组。


> 如果想要实现**每个组之间的差异比较**（相当于比较两个组之间哪些基因的表达存在差异，并求得这种差异的统计显著性），则还需要设置比较矩阵，具体设计方法可以参考每个包的说明文档。一般来说如果不设置的话，默认输出中不包含每个组的多重比较结果。

## 包简介

对于有生物学重复的样品，使用DESeq2进行样品组间的差异表达分析，获得两个生物学条件之间的差异表达基因集；对于无生物学重复我们使用edgeR。这里只是一般的使用状态，没有绝对性，比如，edgeR在有生物学重复的情况下也能使用，只是在识别差异基因的能力和出现假阳性结果等方面两个包存在不同而已。上文提到的三个R包做差异分析都要求输入基因的未经过标准化的 reads 计数数据，而不能是 RPKM、FPKM 等经过标准化的数据。

差异分析之后，还需要用 Benjamini-Hochberg 方法对假设检验概率（P value）进行多重假设检验校正，得到错误发现率（False Discovery Rate，FDR）。差异基因的筛选条件为 |log2Fold Change| >= 1（建议值，可根据实际情况进行调整），且 FDR < 0.05


### DESeq2

DESeq2 基于负二项分布模型，考虑了基因表达数据的离散性和变异性，以及库大小差异对差异分析的影响。DESeq2 通过正态化转换和归一化来减少样本间的技术变异，然后估计基因表达的离散性。它使用负二项分布模型来鉴定差异表达基因，并校正多重检验问题。DESeq2 适用于小样本 RNA 测序数据，特别是在样本数较少的情况下表现较好，能够有效地处理样本间的差异、技术性噪音和批次效应。

### edgeR

edgeR 与 DESeq2类似，edgeR 也考虑了数据的离散性。它使用负二项分布模型和不同的归一化方法，如 TMM（Trimmed Mean of M values），来处理样本之间的技术变异。edgeR 使用一个假设检验框架来鉴定差异表达基因，并采用了类似于 Benjamini-Hochberg 方法的多重检验校正。edgeR 在样本较多的情况下表现较好，适用于中等规模的 RNA 测序数据，具有较高的灵敏度和精确度。使用 edgeR 时注意选择合适的分析算法，**官方建议 bulk RNA-seq 选择 quasi-likelihood(QL) F-test tests，scRNA-seq 或是没有重复样品的数据选用 likelihood ratio test**。

### limma

**DESeq2和edgeR都由limma二次开发而来**

limma 最初是针对基因芯片数据开发的，但后来也被应用于 RNA 测序数据。limma 基于线性模型，使用加权最小二乘法来估计基因表达的差异，并通过贝叶斯方法来校正多重检验问题。limma 在处理大规模数据时表现出色，适用于高通量数据分析，如芯片和大规模RNA测序数据，能够很好地控制假阳性率。

limma 适用于各种类型的高通量数据，包括芯片数据和 RNA-seq。它要求每个基因的表达值，可以是原始计数也可以是已经归一化的表达值。limma 进行差异分析有两种方法 ：limma-trend 和 voom。在样本测序深度相差不大时两种方法差距不大，而测序深度相差大时 voom 更有优势，因此我们一般都选择 voom 的方法进行差异分析

### 总结

1. limma 包做差异分析要求数据满足正态分布或近似正态分布，如基因芯片、TPM 格式的高通量测序数据。
2. 通常认为 Counts 数据不符合正态分布而服从泊松分布。所以对于 Counts 数据来说，用 limma 包进行差异分析，误差较大。
3. edgeR 差异分析速度快，得到的基因数目比较多，假阳性高（实际不差异，结果差异）；DESeq2 差异分析速度慢，得到的基因数目比较少，假阴性高（实际差异，结果不差异）。
4. DESeq2 更关注在小样本条件下的差异分析，提供了一些特有的统计模型和方法。edgeR 强调对技术性噪音和样本之间的变异性进行建模，特别适用于样本数量较大的数据集。limma 则是一种通用的差异分析方法，适用于各种高通量表达数据，不仅适用于 RNA-seq 数据，也适用于其他类型的表达数据，如芯片数据。

## 设计对比矩阵

**DESeq2**

在使用`DESeq()`函数计算差异基因过后，可以使用`resultsNames()`函数查看默认对比组别，这里默认使用第一个组来做对照组，组别之间的差异基因都是跟对照组来比的。如果需要自定义各个组别之间的对比，或者组别为时间序列等其他特殊情况的话，则需要自己设计对比矩阵，可以参考[tavareshugo/tutorial_DESeq2_contrasts (github.com)](https://github.com/tavareshugo/tutorial_DESeq2_contrasts)这篇文章。注意，**DESeq2对直接的所有组别两两对比并不友好，如果要使用，建议用edgeR包。**

**edgeR**

在使用`glmQLFTest()`函数计算差异基因时，可以指定`contrast`参数，该参数接受一个`makeContrasts()`，该函数有两个参数`contrasts`和`levels`，`contrasts`描述了组别的对比关系，`levels`描述了不同的组别。

```R
#A, B, C 三组分别两两比较
makeContrasts(B-A,C-B,C-A,levels=c("A","B","C"))

#A, B, C 三组满足更复杂的函数关系式
makeContrasts(contrasts="A(B+C)/2",levels=c("A","B","C"))
```

更详细的对比设计参考[A guide to creating design matrices for gene expression experiments (bioconductor.org)](https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/designmatrices.html)


## 参考

- [看完还不会来揍我 | 差异分析三巨头 —— DESeq2、edgeR 和 limma 包 | 附完整代码 + 注释 - 知乎 (zhihu.com)](https://zhuanlan.zhihu.com/p/653841949)
- [使用limma、Glimma和edgeR，RNA-seq数据分析易如反掌 (bioconductor.org)](https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow_CHN.html)
- [From reads to genes to pathways: differential... | F1000Research](https://f1000research.com/articles/5-1438/v2)
- [Analyzing RNA-seq data with DESeq2 (bioconductor.org)](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#contrasts)
- [RNA-seq workflow: gene-level exploratory analysis and differential expression (bioconductor.org)](https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#differential-expression-analysis)
- [Exaggerated false positives by popular differential expression methods when analyzing human population samples | Genome Biology | Full Text (biomedcentral.com)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02648-4)