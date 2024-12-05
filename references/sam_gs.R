library(qvalue)

# SAMGS : original method from Dinu et al - BMC Bioinformatics - 2007


sd.na <- function(z){sd(z,na.rm=TRUE)}
mean.na <- function(z){mean(z,na.rm=TRUE)}



rowMeansVars <- function(d,margin=1){
  if(margin==2) d <- t(d)
  m <- rowMeans(d,na.rm=T)
  dif <- d - m
  ssd <- rowSums(dif^2,na.rm=T)
  list("means"=m,
       "sumsquaredif" =ssd,
       "vars" =ssd/(ncol(d)-1),
       "centered rows"=dif )

}

sam.TlikeStat <- function(DATA,           # <=> SAM variant of t-test with equal variances hypothesis
                          cl=NULL,        # as implemented in SAMGS R code furnished by Irina Dinu
                          C1=NULL,
                          C2=NULL,
                          s0=NULL,
                          s0.param=list(nb.Groups=100,mad.Const=.64),
                          alternative=c("two.sided", "greater", "less")[1],
                          conf.level = 0.95 ){
 # DATA : expression data
 #     -> dataframe with  rows=genes,
 #                        columns=samples,
 # C1 : index of samples in class 1    (refering to DATA columns order)
 # C2 : index of samples in class 2    (refering to DATA columns order)
 # cl : factor defining a bipartition of the samples
 #      IN THE SAME ORDER AS IN DATA
 # NB : use     'C1' + 'C2'     XOR   'cl'  alone

    # print("Running sam.TlikeStat")

     if(!is.null(cl)){
         cl <- as.factor(as.character(cl))
         C1 <- which(as.numeric(cl)==1)
         C2 <- which(as.numeric(cl)==2)
     }
     if(is.null(C1) | is.null(C2))stop("Error -  sam.TlikeStat : classes 1 and 2 are undefined.")

     nb.Genes    <- nrow(DATA)
     nb.Samples  <- ncol(DATA)
     C1.size     <- length(C1)
     C2.size     <- length(C2)

     stat.C1            <- rowMeansVars(DATA[,C1])
     stat.C2            <- rowMeansVars(DATA[,C2])
     
    #  print(paste0("typeof(DATA):", typeof(DATA)))
    #  print(paste0("DATA:", DATA))
    #  print(paste0("stat.C1:", stat.C1))
    #  print(paste0("stat.C1.sumsquaredif:", stat.C1$sumsquaredif))
    #  print(paste0("typeof(stat.C1):", typeof(stat.C1)))
    #  print(paste0("typeof(stat.C1[1]):", typeof(stat.C1[1])))
    #  print(paste0("typeof(stat.C1[1][1]):", typeof(stat.C1[1][1])))
    #  print(typeof(stat.C1[1][1][1][1][1][1][1][1]))
    #  print(length(stat.C1))
    #  print(paste0("C1.size:", C1.size))
    #  print(paste0("C2.size:", C2.size))
    #  print(paste0("nb.Samples:", nb.Samples))
    #   print(paste0("nb.Genes:", nb.Genes))

     diffmean.C1C2      <- stat.C1$means - stat.C2$means
     pooledSqrtVar.C1C2 <- sqrt((1/C1.size+1/C2.size)*(stat.C1$sumsquaredif+stat.C2$sumsquaredif)/(nb.Samples-2))

     if(is.null(s0)){
           nb.Groups <- s0.param$nb.Groups
           mad.Const <- s0.param$mad.Const

           tmp <- as.data.frame(cbind(pooledSqrtVar.C1C2,diffmean.C1C2))
           tmp <- tmp[order(tmp[,1]),]

           group.Size     <- as.integer(nb.Genes/nb.Groups)
           percentiles    <- seq(0,1,.05)
           nb.Percentiles <- length(percentiles)
           s0.quantiles   <- quantile(pooledSqrtVar.C1C2,percentiles)

           tt <- matrix(NA,nb.Groups,nb.Percentiles)
           coeffvar <- as.data.frame(cbind(s0.quantiles,rep(NA,nb.Percentiles)))
           for(j in 1:nb.Percentiles){
              x <- matrix(tmp[1:(group.Size*nb.Groups),1]/(tmp[1:(group.Size*nb.Groups),2]+ s0.quantiles[j]),group.Size,nb.Groups)
              tt[,j]=apply(x,2,mad,constant=mad.Const,"na.rm" =TRUE)
              coeffvar[j,2] <- sd.na(tt[,j])/mean.na(tt[,j])
           }

           s0 <-  min(s0.quantiles[coeffvar[,2]==min(coeffvar[,2])])
     }

     tstat <- diffmean.C1C2/(pooledSqrtVar.C1C2+s0)
     df <- nb.Samples-3 # <- s0 is a supplemental parameter to estimate

   	 if (alternative == "less") {
    			pval <- pt(tstat, df)
    			cint <- cbind(rep(-Inf,nb.Genes), tstat + qt(conf.level, df) )
     }else if (alternative == "greater") {
    			pval <- pt(tstat, df, lower = FALSE)
    			cint <- cbind(tstat - qt(conf.level, df), rep(Inf,nb.Genes))
 		 }else {
          pval <- 2 * pt(-abs(tstat), df)
          alpha <- 1 - conf.level
          cint <- qt(1 - alpha/2, df)
          cint <- cbind(tstat -cint,tstat + cint)
 		 }

     list(s0            = s0,
          diffmean      = diffmean.C1C2,
          pooledSqrtVar = pooledSqrtVar.C1C2,
          TlikeStat     = tstat,
          "p.values (using Student law)"=pval,
          gm.C1 = stat.C1$means,
          gm.C2 = stat.C2$means,
          confidence.intervals=cint)
}


GS.format.dataframe.to.list <- function(GS){
     if(is.data.frame(GS)){
         genes <- rownames(GS)
         L <- NULL
         for(ags in names(GS)){
            w <- which(GS[,ags]==1)
            if(length(w)>0)  {
               L <- c(L,list(genes[w]))
               names(L)[length(L)] <- ags
            }
         }
         L
     }else{
         GS
     }
}

SAMGS <- function(GS,
                  DATA,
                  cl,
                  nbPermutations=1, # Changed to 1 because we're not using permutations in this R code
                  silent=F){

 # GS : gene sets
 #      -> a dataframe with rows=genes,
 #                        columns= gene sets,
 #                        GS[i,j]=1 if gene i in gene set j
 #                        GS[i,j]=0 otherwise
 #         OR
 #         a list with each element corresponding to a gene set = a vector of strings (genes identifiers)
 #
 #
 # DATA : expression data
 #     -> a dataframe with  rows=genes,
 #                        columns=samples
 #
 # cl : a factor defining a bipartition of the samples IN THE SAME ORDER AS IN DATA
 #
    # print("Running SAMGS")

     GS <-  GS.format.dataframe.to.list(GS)   
     if(!silent) print("GS dataframe-to-list : done.")

    # genes <- unique(unlist(GS))
     genes <- unique(intersect(dimnames(DATA)[[1]],unlist(GS)))
     genes <- genes[!is.na(genes)]  #####only the row numbers with gene in at least 1 set

     nb.Samples  <- ncol(DATA)
     C1.size <- table(cl)[1]  # nb of samples in class 1

     # finding constant s0 for SAM-like test
     tmp     <- sam.TlikeStat(DATA,cl=cl)
     s0      <- tmp$s0
     if(!silent) print(paste0("s0 estimation : ", s0, " done."))

     DATA <- DATA[genes,]
     GS <-  lapply(GS,function(z) which(genes %in% z))
	 nb.GeneSets <- length(GS)   # nb of gene sets
     GeneSets.sizes <- sapply(GS,length) # size of each gene set

     # stats obtained on 'true' data
	 samT.ok <-as.data.frame(tmp$TlikeStat)[genes,]	# SAM T-like statistic for each gene

     sam.sumsquareT.ok <- sapply(GS,function(z) sum(samT.ok[z]^2))  # SAMGS statitic for each gene set

     # stats obtained on 'permuted' data
     permut.C1 <- matrix(NA,nbPermutations,C1.size)
     sam.sumsquareT.permut <- matrix(NA,nbPermutations,nb.GeneSets)
     for(i in 1:nbPermutations) {
         C1.permut <-  permut.C1[i,] <- sample(nb.Samples,C1.size)
         C2.permut <-  (1:nb.Samples)[-C1.permut]
         samT.permut <- sam.TlikeStat(DATA,C1=C1.permut,C2=C2.permut,s0=s0)$TlikeStat
         sam.sumsquareT.permut[i,] <- sapply(GS,function(z) sum(samT.permut[z]^2)) # SAMGS statitic for each gene set  - for current permutation
         if(!silent & i%%50 == 0)print(paste(i," permutations done."))
     }

     GeneSets.pval <- apply(t(sam.sumsquareT.permut) >= sam.sumsquareT.ok ,1,sum)/nbPermutations
     qobj <- NULL
     try(qobj <- qvalue(GeneSets.pval))
     GeneSets.qval <- rep(NA,nb.GeneSets)
     PI0 <- NA

     if(!is.null(attr(qobj,"class"))){
         GeneSets.qval <- qobj$qvalues
         PI0 <- qobj$pi0
     }

     res <- as.data.frame(cbind("GS size"              = GeneSets.sizes,
                                "GS p-value (0 <=> < 1/nb permutations)" = GeneSets.pval,
                                "GS q-value"           = GeneSets.qval ))

     res <- cbind(res,"GS name"= names(GS))[c(4,1:3)]


     L <- list("GS stats"=res,
               "pi0"=PI0,
               "genes stats"= as.data.frame(cbind(tmp[[4]],tmp[[5]])))
     names(L[[3]])<- c("SAM T like stat","p-value (Student law)")
     attr(L$pi0,"info") <- "estimated proportion of Gene Sets non significantly associated to the phenotype of interest"
     L
}





# You just need to copy the above function to R once, then you can run SAM-GS by calling the function SAMGS.
# Please see the example below.


#---------------------------------------------------------------------------------------
###############examples can be downloaded from our website ################
# library(qvalue)

# #---------------------------------------------------------------------------------------
# #--------------------------- how to use? -----------------------------------------------
# #---------------------------------------------------------------------------------------


# #---------------------------------------------------------------------------------------------------------#
# #-- if you already have the gene set file in the the binary format, like c2part1.csv on SAM-GS website#
# #---------------------------------------------------------------------------------------------------------#


# #-- read the set definitions

# try1=read.csv(file="c2part1.csv",header = TRUE, sep = ",",fill = T)
# try2=read.csv(file="c2part2.csv",header = TRUE, sep = ",",fill = T)
# gene.sets=cbind(try1,try2)

# # read the gene expressions

# exprs=read.csv("p53.csv",sep=",",header=T)
# exprs==as.data.frame(exprs)
# dim(exprs)
# p53=exprs[,-c(1,2)]


# # define the number of samples in each group
# n1=33
# n2=17
# cl=c(rep(1,n1),rep(0,n2))

# # Run SAM-GS

# out=SAMGS(GS=gene.sets,DATA=p53,cl=cl,nbPermutations=100) 
# out1=as.data.frame(out[[1]])
# out1[,1] <- as.character(out1[,1])
# out1[1:5,]



# #----------------------------------------------------------------------------------------------
# #---------   if you have the set file in the GSEA .gmt files.--- --- --- --- --- --- --- ---  #
# #####		no need to make the 0/1 binary sets file 		####
# #----------------------------------------------------------------------------------------------

# C2 <- readLines( "c2.v2.symbols.gmt" )
# C2 <- strsplit(C2,"\t")
# nam <- sapply(C2,function(z)z[1])
# C2 <- lapply(C2,function(z)z[-(1:2)])
# names(C2) <- nam   #####a list

# #read gene expression data

# exprs=read.csv("p53.csv",sep=",",header=T)
# exprs=as.data.frame(exprs)
# dim(exprs)
# p53=as.matrix(log2(exprs[,-c(1,2)]))
# dimnames(p53)[[1]] <- as.character(exprs[,1])



# # define the number of samples in each group
# n1=33
# n2=17
# cl=c(rep(1,n1),rep(0,n2))

# # Run SAM-GS
# set.seed(1981)
# out=SAMGS(GS=C2,DATA=p53,cl=cl,nbPermutations=100)
# out1=as.data.frame(out[[1]])
# out1[,1] <- as.character(out1[,1])
# out1[1:5,]
