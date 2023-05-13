Tree building with ape and phangorn
ape

>install.packages("adegenet", dep=TRUE)
>install.packages("phangorn", dep=TRUE)

>library(ape)
>library(adegenet)
>library(phangorn)

>dna <- fasta2DNAbin(file="http://adegenet.r-forge.r-project.org/files/usflu.fasta")

>D <- dist.dna(dna, model="TN93")

>tre <- nj(D)

>plot(tre, cex=.6)
>title("A simple NJ tree")

nj- neighbor joinng algorithm

phangorn
>dna2 <- as.phyDat(dna)

>tre.ini <- nj(dist.dna(dna,model="raw"))
>parsimony(tre.ini, dna2)

(returns 422)

> tre.pars <- optim.parsimony(tre.ini, dna2)

(returns Final p-score 420 after  2 nni operations)

>plot(tre.pars, cex=0.6)

For my data:
ape try 1

>library(ape)
>library(adegenet)
>library(phangorn)

>aa <- fasta2DNAbin("/Users/brook/Desktop/bot563/data/hh_sequences1.fasta")

>A <- dist.aa(aa) #aa is type DNAbin, not AAbin but it still works, does this affect out come?

>tre <- nj(A)

>plot(tre, cex=.6)

phangorn try 1
>aa2 <- as.phyDat.AAbin(aa)

Warning message:
In phyDat.AA(data, return.index = return.index, ...) :
  Found unknown characters (not supplied in levels). Deleted sites with unknown states.
> aa2

32 sequences with 212 character and 135 different site patterns.
The states are A R N D C Q E G H I L K M F P S T W Y V 

> tre.ini <- nj(dist.aa(aa,model="raw"))

Error in dist.aa(aa, model = "raw") : unused argument (model = "raw")

> parsimony(tre.ini, aa2)

Error in is.binary(tree) : object 'tre.ini' not found

ape try 2
> aa <- read.FASTA("/Users/brook/Desktop/bot563/data/hh_sequences1.fasta", type="AA")

>A <- dist.aa(aa)

Error in numeric(n * (n - 1)/2) : invalid 'length' argument

#I need aligned file to be read

ape try 3 (with cds seq)

 dna <- fasta2DNAbin("/Users/brook/Desktop/bot563/data/hh_cds2.fasta")
 D <- dist.dna(dna, model="TN93")
 tre <- nj(D)
 tre <- ladderize(tre)
 plot(tre, cex=.6)

phangorn try 2 (cds)
dna2 <- as.phyDat(dna)
tre.ini <- nj(dist.dna(dna,model="raw"))
parsimony(tre.ini, dna2)
tre.pars <- optim.parsimony(tre.ini, dna2)
plot(tre.pars, cex=0.6)

#check if nj is good method
tre2 <- root(tre, out=1)
tre2 <- ladderize(tre2)
x <- as.vector(D)
y <- as.vector(as.dist(cophenetic(tre2)))
plot(x, y, xlab="original pairwise distances", ylab= "pairwise distances on the tree", main= "Is NJ appropriate?", pch=20, col=transp("black",.1), cex=3)
abline(lm(y~x), col="red")
cor(x,y)^2
#correlation of .9975