### R code from vignette source 'hypred.Rnw'

###################################################
### code chunk number 1: setOps
###################################################
op <- options(); utils::str(op) 
options(SweaveSyntax="SweaveSyntaxNoweb")
options(continue=" ")
options(prompt = "R> ")


###################################################
### code chunk number 2: load
###################################################
library(hypred)


###################################################
### code chunk number 3: define.genome
###################################################
genome <- hypredGenome(num.chr = 2, 
                       len.chr = c(1.0, 1.0), 
                       num.snp = 100)


###################################################
### code chunk number 4: summary
###################################################
summary(genome)


###################################################
### code chunk number 5: show.random.map
###################################################
slot(genome, "pos.snp")[1:5] ## pos. of first 5 


###################################################
### code chunk number 6: new.map
###################################################
map <- rep(seq(0, 0.99, by = 0.01 ), 2)
genome <- hypredNewMap(genome, 
                       new.map = map)
                       


###################################################
### code chunk number 7: orderly.map
###################################################
slot(genome, "pos.snp")[1:5] ## pos. of first 5 


###################################################
### code chunk number 8: snp.dist
###################################################
hypredSNPdist(genome, 
              chromosome = 1, 
              SNP1 = 1, 
              SNP2 = 3)              


###################################################
### code chunk number 9: qtl.vect
###################################################
qtl.ids <- c(1, 20, 40, 60, 80,  ## chr. 1
             101, 120, 140, 160, 180)  ## chr. 2


###################################################
### code chunk number 10: qtl.vect2
###################################################
qtl.dom.ids <- c(40,  ## chr. 1
                 140)  ## chr. 2


###################################################
### code chunk number 11: qtl.vect3
###################################################
per.mar.ids <- c(20,  ## chr. 1
                 120)  ## chr. 2


###################################################
### code chunk number 12: qtl
###################################################
genome <- hypredNewQTL(genome, 
                       new.id.add = qtl.ids, 
                       new.id.dom = qtl.dom.ids, 
                       new.id.per.mar = per.mar.ids, 
                       new.eff.add =  rep(1, 10), 
                       new.eff.dom = c(0.5, 0.5)
                       )


###################################################
### code chunk number 13: summary2
###################################################
summary(genome)


###################################################
### code chunk number 14: founder
###################################################
founder1 <- hypredFounder(genome, 
                          prob.snp = 1)[1, ]
founder2 <- hypredFounder(genome, 
                         prob.snp = 1)[2, ]


###################################################
### code chunk number 15: gen.founder
###################################################
founder1[1:5] ## first 5 loci
founder2[1:5] ## first 5 loci


###################################################
### code chunk number 16: intro.hypred
###################################################
set.seed(114)## to see a recomb. event
new.gamete <- hypredRecombine(genome, 
                              genomeA = founder1, 
                              genomeB = founder2, 
                              mutate = TRUE, 
                              mutation.rate.snp = 2.5 * 10^-5, 
                              mutation.rate.qtl = 2.5 * 10^-5, 
                              block = FALSE)


###################################################
### code chunk number 17: show.gam
###################################################
new.gamete[1:20] ## first 20 loci


###################################################
### code chunk number 18: pre.all.mat.F2
###################################################
N <- 250 ## number of individuals
F2 <- matrix(nrow = N * 2, 
             ncol = 200)


###################################################
### code chunk number 19: f2
###################################################
for(i in 1:(N*2))
  {
    F2[i, ] <- hypredRecombine(genome, 
                               genomeA = founder1, 
                               genomeB = founder2, 
                               mutate = TRUE, 
                               mutation.rate.snp = 2.5 * 10^-5, 
                               mutation.rate.qtl = 2.5 * 10^-5, 
                               block = FALSE)
  }


###################################################
### code chunk number 20: random.mate.var
###################################################
G <- 50 ## number of generations

random.mate.pop <- F2 ## identical to the F2 at start
random.mate.pop.temp <- matrix(nrow = N*2, ncol = 200)


###################################################
### code chunk number 21: random.mate
###################################################

for(generation in 1 : G) ## loop over generations
  {
    gameteIndex1 <- 1 ## indexing variables
    gameteIndex2 <- 2
    for(indiv in 1 : N) ## loop over individuals
      {
        ## gamete 1
        random.mate.pop.temp[gameteIndex1,] <-         
          hypredRecombine(genome,
                          genomeA = random.mate.pop[gameteIndex1,],
                          genomeB = random.mate.pop[gameteIndex2,],
                          mutate = TRUE,
                          mutation.rate.snp = 2.5 * 10^(-5),
                          mutation.rate.qtl = 2.5 * 10^(-5),
                          block = FALSE)
        ## gamete 2
        random.mate.pop.temp[gameteIndex2,] <-         
          hypredRecombine(genome,
                          genomeA = random.mate.pop[gameteIndex1,],
                          genomeB = random.mate.pop[gameteIndex2,],
                          mutate = TRUE,
                          mutation.rate.snp = 2.5 * 10^(-5),
                          mutation.rate.qtl = 2.5 * 10^(-5),
                          block = FALSE)

        ## increment to next individual
        gameteIndex1 <- gameteIndex1 + 2 
        gameteIndex2 <- gameteIndex2 + 2
      } ## end for N
    
    ## permutate
    random.mate.pop <- random.mate.pop.temp[sample(1 : (N * 2)), ]
  } ## end for G


###################################################
### code chunk number 22: g.value
###################################################
g.value <- hypredTruePerformance(genome, 
                                 random.mate.pop, 
                                 DH = FALSE)
summary(g.value)


###################################################
### code chunk number 23: varAndmean
###################################################

##-------------------------
## compute expected variance and mean

## additive

## purely additive loci
p.add <- (colSums(random.mate.pop)/nrow(random.mate.pop))[qtl.ids[!qtl.ids%in%qtl.dom.ids]]
q.add <- 1 - p.add

var.add <- sum((2*p.add*q.add*1^2))
mean.add <- sum(2*p.add*1)

## additive + dominance loci
p.add.dom <-  (colSums(random.mate.pop)/nrow(random.mate.pop))[qtl.dom.ids]
q.add.dom <-  1 - (colSums(random.mate.pop)/nrow(random.mate.pop))[qtl.dom.ids]

var.add.dom <- sum(2 * p.add.dom * q.add.dom * ((1 + 0.5*(q.add.dom - p.add.dom))^2))
mean.add.dom <- sum(2*p.add.dom*(1 + (1 - p.add.dom ) * 0.5))

## dominance variance
var.dom <- sum((2 * p.add.dom * q.add.dom * 0.5)^2)

## expectations

meanExp <- mean.add + mean.add.dom
varExp <- var.add + var.add.dom + var.dom


###################################################
### code chunk number 24: pheno
###################################################
var.env <- var(g.value)
phen.value <- g.value + rnorm(250, 0,  sqrt(var.env))


###################################################
### code chunk number 25: select
###################################################
index.individual <- which(phen.value == max(phen.value)) 
index.row1 <- index.individual * 2 - 1
index.row2 <- index.individual * 2 
selected <- random.mate.pop[c(index.row1, index.row2), ]


###################################################
### code chunk number 26: pre.all.mat.DH
###################################################
DH <- matrix(nrow = N, 
             ncol = 200)


###################################################
### code chunk number 27: DH
###################################################
for(i in 1:N)
  {
    DH[i, ] <- hypredRecombine(genome, 
                               genomeA = selected[1, ], 
                               genomeB = selected[2, ], 
                               mutate = TRUE, 
                               mutation.rate.snp = 2.5 * 10^-5, 
                               mutation.rate.qtl = 2.5 * 10^-5, 
                               block = FALSE)
  }


###################################################
### code chunk number 28: precalc.mean.DH
###################################################
mean.DH <- mean(hypredTruePerformance(genome, 
                                      DH, 
                                      DH = TRUE))


###################################################
### code chunk number 29: QTL.geno
###################################################
selected[, qtl.ids] ## all QTL
selected[, qtl.dom.ids] ## QTL with dominance

## effect a = 1 times the number of 1 alleles
sum(selected[, qtl.ids] == 1) +
  ## 0.5 if the dominance QTL have genotype 01 or 10 (hence sum = 1)
  (sum(selected[, qtl.dom.ids[1]]) == 1)*0.5 +
  (sum(selected[, qtl.dom.ids[2]]) == 1)*0.5 


###################################################
### code chunk number 30: pheno.DH
###################################################
g.value.DH <- hypredTruePerformance(genome, 
                                    DH, 
                                    DH = TRUE)

phen.value.DH <- g.value.DH + rnorm(250, 0,  sqrt(var.env))


###################################################
### code chunk number 31: desg
###################################################
design.DH <- hypredCode(genome, 
                        genotypes = DH, 
                        DH = TRUE, 
                        type = "012")


###################################################
### code chunk number 32: not.per.mar
###################################################
design.DH[1, 1:5]


###################################################
### code chunk number 33: per.mar
###################################################
design.DH[1, 17:21]


