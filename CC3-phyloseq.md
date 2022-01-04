Analyse des données avec Phyloseq
================
Mélody Lebrun
3 janvier 2022

-   [Préparation de l’environnement](#préparation-de-lenvironnement)
-   [Construire un arbre phylogénétique permettant de relié les variants
    de
    séquences](#construire-un-arbre-phylogénétique-permettant-de-relié-les-variants-de-séquences)
-   [Combiner des données dans un objet
    phyloseq](#combiner-des-données-dans-un-objet-phyloseq)
-   [Création de plusieurs objects](#création-de-plusieurs-objects)
    -   [Utilisation de phyloseq](#utilisation-de-phyloseq)
-   [Filtration “non supervisées” reposant sur les données de
    l’article](#filtration-non-supervisées-reposant-sur-les-données-de-larticle)
-   [Taxons agglomérés](#taxons-agglomérés)
-   [Transformation de la valeur de
    l’abondance](#transformation-de-la-valeur-de-labondance)
-   [Analyse complémentaire](#analyse-complémentaire)

# Préparation de l’environnement

``` r
library(phyloseq)
library(ggplot2)
library(readr)
```

``` r
library("knitr")
library("BiocStyle")
.cran_packages <-c ("ggplot2", "gridExtra", "devtools")
install.packages(.cran_packages) 
```

    ## Installing packages into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)

``` r
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")
BiocManager::install(.bioc_packages)# Load packages into session, and print package version
```

    ## 'getOption("repos")' replaces Bioconductor standard repositories, see
    ## '?repositories' for details
    ## 
    ## replacement repositories:
    ##     CRAN: https://packagemanager.rstudio.com/all/__linux__/focal/latest

    ## Bioconductor version 3.14 (BiocManager 1.30.16), R 4.1.2 (2021-11-01)

    ## Warning: package(s) not installed when version(s) same as current; use `force = TRUE` to
    ##   re-install: 'dada2' 'phyloseq' 'DECIPHER' 'phangorn'

    ## Installation paths not writeable, unable to update packages
    ##   path: /usr/local/lib/R/library
    ##   packages:
    ##     Matrix

    ## Old packages: 'broom', 'gert'

``` r
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
```

    ## Loading required package: gridExtra

    ## Loading required package: devtools

    ## Loading required package: usethis

    ## Loading required package: dada2

    ## Loading required package: Rcpp

    ## Loading required package: DECIPHER

    ## Loading required package: Biostrings

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     distance

    ## Loading required package: XVector

    ## Loading required package: GenomeInfoDb

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## Loading required package: RSQLite

    ## Loading required package: parallel

    ## Loading required package: phangorn

    ## Loading required package: ape

    ## 
    ## Attaching package: 'ape'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     complement

    ##   ggplot2 gridExtra  devtools     dada2  phyloseq  DECIPHER  phangorn 
    ##      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE      TRUE

``` r
load("02_Data-analysis-with-DADA2_FinalEnv")
```

# Construire un arbre phylogénétique permettant de relié les variants de séquences

``` r
seqs <- getSequences(seqtabNoC)
names(seqs) <- seqs 
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
```

``` r
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm)
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
        rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)
```

# Combiner des données dans un objet phyloseq

``` r
Metadonnees<-read.csv(file="SraRunTable(2).txt")
```

``` r
Metadonnees$Sample.Name <- paste0(gsub("00", "", Metadonnees$Sample.Name), "D")
Metadonnees <- Metadonnees[!duplicated(Metadonnees$Sample.Name),] # Remove dupicate entries for reverse reads
rownames(seqtabAll) <- gsub("124", "125", rownames(seqtabAll)) # Fix discrepancy
all(rownames(seqtabAll) %in% Metadonnees$Sample.Name) # TRUE
```

    ## [1] FALSE

``` r
rownames(Metadonnees) <- Metadonnees$Sample.Name
keep.cols <- c("Run","env_medium","host_issue_sampled","lat_lon") 
```

######## pb

Metadonnees \<- Metadonnees\[rownames(seqtabAll), keep.cols\]

# Création de plusieurs objects

L’idée est de stratifier les données de tel sortes à simplifier le
tableau et de séparer les types d’échantillon, entre tissu intestinal,
contenu intestinal et sédiment et entre station.

``` r
SediCH=(Metadonnees$env_medium=="shallow marine sediment [ENVO:00000428]" | Metadonnees$lat_lon=="62.21879 S 58.95786 W")
SediA=(Metadonnees$env_medium=="shallow marine sediment [ENVO:00000428]" | Metadonnees$lat_lon=="62.209556 S 58.92900 W")
TissuCH=(Metadonnees$host_tissue_sampled=="Gut membrane" | Metadonnees$lat_lon=="62.21879 S 58.95786 W")
TissuA=(Metadonnees$host_tissue_sampled=="Gut membrane" | Metadonnees$lat_lon=="62.209556 S 58.92900 W")
ContenuCH=(Metadonnees$host_tissue_sampled=="Gut content" | Metadonnees$lat_lon=="62.21879 S 58.95786 W")
ContenuA=(Metadonnees$host_tissue_sampled=="Gut content" | Metadonnees$lat_lon=="62.209556 S 58.92900 W")
```

##########pb##########

library(phyloseq) ps \<- phyloseq(otu_table(seqtabNoC,
taxa_are_rows=FALSE), sample_data(Metadonnees),
tax_table(taxTab),phy_tree(fitGTR$tree)) ps \<-
prune_samples(sample_names(ps) != “Mock”, ps) # Remove mock sample ps

## Utilisation de phyloseq

``` r
rank_names(ps) #rang de la base de donnée
```

    ## [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"

``` r
table(tax_table(ps)[, "Phylum"], exclude = NULL) # creation d'une table montrant les phylums avec une caractéristique observée
```

    ## 
    ##              Actinobacteria               Bacteroidetes 
    ##                          13                          23 
    ## Candidatus_Saccharibacteria   Cyanobacteria/Chloroplast 
    ##                           1                           4 
    ##         Deinococcus-Thermus                  Firmicutes 
    ##                           1                         327 
    ##                Fusobacteria              Proteobacteria 
    ##                           1                          11 
    ##                 Tenericutes             Verrucomicrobia 
    ##                           1                           1

``` r
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized")) # supression des phylums non caractérisés
```

``` r
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))
# prévalence des caractéristiques dans l'ensemble des données cad le nb de fois que le taxon apparait dans l'échantillon
```

``` r
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
```

    ##                         Phylum         1     2
    ## 1               Actinobacteria 120.15385  1562
    ## 2                Bacteroidetes 265.52174  6107
    ## 3  Candidatus_Saccharibacteria 280.00000   280
    ## 4    Cyanobacteria/Chloroplast  64.25000   257
    ## 5          Deinococcus-Thermus  52.00000    52
    ## 6                   Firmicutes 179.24771 58614
    ## 7                 Fusobacteria   2.00000     2
    ## 8               Proteobacteria  59.09091   650
    ## 9                  Tenericutes 234.00000   234
    ## 10             Verrucomicrobia 104.00000   104

``` r
# calcule des prévalences totales et moyennes des caractéristiques dans chaque embranchement
```

``` r
filterPhyla = c("Fusobacteria", "Deinococcus-Thermus")
# Filtrations "supervisées" des phylums dans la base de donnée de référence taxonomique
ps1 = subset_taxa(ps, !Phylum %in% filterPhyla)
ps1
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 381 taxa and 360 samples ]
    ## sample_data() Sample Data:       [ 360 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 381 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 381 tips and 379 internal nodes ]

# Filtration “non supervisées” reposant sur les données de l’article

``` r
#lien entre la prévalence et le nombre total de lectures pour chaque caractéristique.

prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +

  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
```

![](CC3-phyloseq_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
# Le seuil de prévalence est définie à 5% pour l'ensemble des échantillons
prevalenceThreshold = 0.05 * nsamples(ps)
prevalenceThreshold
```

    ## [1] 18

``` r
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps)
```

# Taxons agglomérés

``` r
length(get_taxa_unique(ps2, taxonomic.rank = "Genus"))
```

    ## [1] 49

``` r
ps3 = tax_glom(ps2, "Genus", NArm = TRUE)
h1 = 0.4
ps4 = tip_glom(ps2, h = h1)
```

``` r
multiPlotTitleTextSize = 15
p2tree = plot_tree(ps2, method = "treeonly",
                   ladderize = "left",
                   title = "Before Agglomeration") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p3tree = plot_tree(ps3, method = "treeonly",
                   ladderize = "left", title = "By Genus") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p4tree = plot_tree(ps4, method = "treeonly",
                   ladderize = "left", title = "By Height") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
```

# Transformation de la valeur de l’abondance

Transforme les données afin de tenir compte des différences de taille,
de variance, d’échelle pour les abondances.

``` r
plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "sex",y = "Abundance",
                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}
```

Abondance relative

``` r
ps3ra = transform_sample_counts(ps3, function(x){x / sum(x)})
```

Abondance avant et après transformation

``` r
plotBefore = plot_abundance(ps3,"")
plotAfter = plot_abundance(ps3ra,"")
grid.arrange(nrow = 2,  plotBefore, plotAfter)
```

![](CC3-phyloseq_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

# Analyse complémentaire

``` r
.cran_packages <- c( "shiny","miniUI", "caret", "pls", "e1071", "ggplot2", "randomForest", "dplyr", "ggrepel", "nlme", "devtools",
                  "reshape2", "PMA", "structSSI", "ade4",
                  "ggnetwork", "intergraph", "scales")
.github_packages <- c("jfukuyama/phyloseqGraphTest")
.bioc_packages <- c("genefilter", "impute")
# Install CRAN packages (if not already installed)
.inst <- .cran_packages %in% installed.packages()
if (any(!.inst)){
  install.packages(.cran_packages[!.inst],repos = "http://cran.rstudio.com/")
}
```

    ## Installing package into '/usr/local/lib/R/site-library'
    ## (as 'lib' is unspecified)

    ## Warning: package 'structSSI' is not available for this version of R
    ## 
    ## A version of this package for your version of R might be available elsewhere,
    ## see the ideas at
    ## https://cran.r-project.org/doc/manuals/r-patched/R-admin.html#Installing-packages

``` r
.inst <- .github_packages %in% installed.packages()
if (any(!.inst)){
  devtools::install_github(.github_packages[!.inst])
}
```

    ## Skipping install of 'phyloseqGraphTest' from a github remote, the SHA1 (3fb6c274) has not changed since last install.
    ##   Use `force = TRUE` to force installation

``` r
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)){
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst])
}
```
