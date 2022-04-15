## Article title: The effects of species abundance, spatial distribution, and phylogeny on a plant-ectomycorrhizal fungal network
## Authors: Chunchao Zhu, Zihui, Wang, David C. Deane, Wenqi Luo, Yongfa Chen, Yongjun Cao, Yumiao Lin, Minhua Zhang 
## E-mail: zhuchunchao11@163.com; 641246135@qq.com


##--------------------------------------------------------------------------------------------
## Data:
## wEMspmat: a plant-EM fungal interaction matrix including 862 EM fungal species and 43 plant species
##--------------------------------------------------------------------------------------------

##---------------------------------------------------------------------------------------------
#(1)some R codes to generate proprability matrix from abundance, spatial overlap, null and observed networks

#   Np.web.em: a propability matrix generated from relative abundances of EM and host plants
#   S.overlap.EM.PL.ses: a propability matrix generated from spatial overlap for EM fungi and plants
#   NullEM.pweb: a null probability matrix
#   EM.NandS.pweb.ses: a propability matrix generated from abundance and spatital information
#   obsEM.pweb: an obs probability matrix generated from observed plant-EM matrix

##---------------------------------------------------------------------------------------------


##--------------------------------------------------------------------------------------------
#  relative abundance of host plants
##--------------------------------------------------------------------------------------------
#  abund213df     # dataframe of abundance of 213 plant species in the HSD forest plot
    id =  match(  rownames( wEMspmat) , rownames ( abund213df) )   
    Host43abund <-   as.data.frame(abund213df[id,])
    rownames( Host43abund) <- rownames(abund213df)[id]
    colnames( Host43abund) <- "plantabund"
    rel.abund.host <-  as.matrix( Host43abund/sum(Host43abund[,1]))

##---------------------------------------------------------------------------------------------
##   relative abundance of EM   
##---------------------------------------------------------------------------------------------
    rel.abund.em <-  as.matrix(apply( wEMspmat,2, sum)/sum(wEMspmat))

##--------------------------------------------------------------------------------------------
##  a propability matrix generated from relative abundances of EM and host plants
##  Np.web.em
##--------------------------------------------------------------------------------------------
##  EM and hosts
    t.rel.abund.em  <- t(rel.abund.em)
    rel.abund.em.new <- matrix(0,43,1)
    rel.abund.em.new[,1] <- rel.abund.host[,1]
    Np.web.em <- rel.abund.host[,1] %*%  t.rel.abund.em

##---------------------------------------------------------------------------------------------
##  a propability matrix generated from spatial overlap for EM fungi and plants
##  S.overlap.EM.PL
##---------------------------------------------------------------------------------------------
## PLquad.mat.bin: a matrix indicating plant presence/obsence on quadrats
## EMquad.mat.bin: a matrix indicating EM presence/obsence on quadrats

   S.overlap.EM.PL <- t(PLquad.mat.bin) %*% (EMquad.mat.bin)

## standardized spatital overlaped likelihood matrix

   S.overlap.EM.PL.ses  <-  S.overlap.EM.PL/sum( S.overlap.EM.PL)

##--------------------------------------------------------------------------------------------
##  a null probability matrix
##  NullEM.pweb
##---------------------------------------------------------------------------------------------

    nr = nrow(wEMspmat);
    nc = ncol(wEMspmat); 
    NullEM.pweb <-  matrix(1/( nr* nc), nr, nc)

##---------------------------------------------------------------------------------------------
##  an obs probability matrix generated from observed plant-EM matrix
##  obsEM.pweb
##---------------------------------------------------------------------------------------------
    obsEM.pweb <-   wEMspmat/ sum(wEMspmat)

##---------------------------------------------------------------------------------------------

##---------------------------------------------------------------------------------------------          
##  a Absp propability matrix generated from abundance and spatital information
## EM.NandS.pweb.ses
##---------------------------------------------------------------------------------------------
   EM.NandS.pweb <-  Np.web.em * S.overlap.EM.PL
   EM.NandS.pweb.ses <-  EM.NandS.pweb/sum(EM.NandS.pweb)   # standardized 
 
##---------------------------------------------------------------------------------------------



##---------------------------------------------------------------------------------------------- 
## (2) some R functions to estimate network metircs and
##----------------------------------------------------------------------------------------------
##  These R functions are downloaded from V芍zquez's published article (2009)
##  https://figshare.com/collections/Evaluating_multiple_determinants_of_the_structure_of_plant_animal_mutualistic_networks/3301145/1
##----------------------------------------------------------------------------------------------

##----------------------------------------------------------------------------------------------
##     "quant2bin" function
##----------------------------------------------------------------------------------------------
## Transformation of quantitative matrix into binary

 quant2bin<-function(matr)
 {
  cn=colnames(matr)
  rn=rownames(matr)
  ij<-which(matr!=0,arr.ind=T)
  matrb<-matrix(0,dim(matr)[1],dim(matr)[2])
  matrb[ij]<-1
  colnames(matrb)=cn
  rownames(matrb)=rn
  return(matrb)
 }       
##----------------------------------------------------------------------------------------------


##----------------------------------------------------------------------------------------------
##      "confint" function  
##----------------------------------------------------------------------------------------------
 confint<-function(x,alpha=0.05){
   if (is.na(ncol(x))==TRUE){
    xs<-sort(x)
    n<-length(xs)
    lo<-xs[round(n*alpha/2)]
    up<-xs[round(n*(1-alpha/2))]
    res<-data.frame(lo=lo,up=up)
   }
  if (is.na(ncol(x))==FALSE){
    res<-data.frame(matrix(0,ncol(x),2))
    colnames(res)=c("lo","up")
    rownames(res)=colnames(x)
    for (i in 1:ncol(x)){
      xs<-sort(x[,i])
      n<-length(xs)
      if(round(n*alpha/2)==0) ii=1
      if(round(n*alpha/2)>0) ii=round(n*alpha/2)
      lo<-xs[ii]
      up<-xs[round(n*(1-alpha/2))]
      res[i,]<-c(lo,up)
     }
   }
   return(res)
  }

##----------------------------------------------------------------------------------------------
##       "intereven" function 
##----------------------------------------------------------------------------------------------
  intereven<-function(imatr) {
   p_i.mat <- imatr/sum(imatr)
   ie<--sum(p_i.mat*log2(p_i.mat),na.rm=TRUE)/log2(sum(imatr!=0))
   ie
  }

#----------------------------------------------------------------------------------------------
#        "intasymm" function
#----------------------------------------------------------------------------------------------
 intasymm<-function(imatr){
  res<-list()
  nr=nrow(imatr)
  nc=ncol(imatr)
  sumr=matrix(rowSums(imatr),nr,nc)
  sumc=matrix(colSums(imatr),nr,nc,byrow=TRUE)
  s.r<-imatr/sumr
  s.c<-imatr/sumc
  dij<-s.c-s.r
  a.r<-colSums(-dij)/colSums(imatr>0)
  a.c<-rowSums(dij)/rowSums(imatr>0)
  res$s.r<-s.c
  res$s.c<-s.r
  res$a.r<-a.r
  res$a.c<-a.c
  res
 }

#-----------------------------------------------------------------------------------------------
#         "mgen" function
#-----------------------------------------------------------------------------------------------
 mgen<-function(m,n,e=TRUE,zs=FALSE){
  #m<-r%*%t(c) #Interaction probability matrix
  mac<-matrix(cumsum(m),nrow(m),ncol(m)) #Cumulative probability matrix
  mint<-matrix(0,nrow(m),ncol(m)) #Interaction matrix

  if (zs==FALSE){
    for (i in 1:nrow(m)){
      c1<-sample(ncol(m),replace=T,prob=colSums(m)); c1<-c1[1]
      mint[i,c1]<-1
    }
    for (i in 1:ncol(m)){
      r1<-sample(nrow(m),replace=T,prob=rowSums(m)); r1<-r1[1]
      mint[r1,i]<-1
    }
  }
  while(sum(mint)<n){
    rand<-runif(1,0,1)
    ri<-min(which(mac>=rand))
    mint[ri]<-mint[ri]+1
  }
  #if((e==T)&(zs=T)) mint<-mred(mint)
  return(mint)
 }
##-----------------------------------------------------------------------------------------------
##         mlik function
##-----------------------------------------------------------------------------------------------
 mlik<-function(o,m.p,par){
    lik<--dmultinom(o,prob=m.p,log=T) #Negative log likelihood
    aic=2*lik+2*par
    res=data.frame(lik,aic)
    return(res)
 }
##-----------------------------------------------------------------------------------------------

##-----------------------------------------------------------------------------------------------
##        "netstats" function
##-----------------------------------------------------------------------------------------------

 netstats  <-  function(imatr,randomize=TRUE,iter= 100,pmatr=NULL){  #iter=1000
    conn = sum(quant2bin(imatr))/(nrow(imatr)*ncol(imatr)) # Connectance
    nest = nested( imatr, method= "weighted NODF" ,rescale=FALSE, normalised=TRUE) # Nestedness
    mod  = computeModules(imatr, method="Beckett", deep = FALSE, deleteOriginalFiles = TRUE,
           steps = 1000000, tolerance = 1e-10, experimental = FALSE, forceLPA=FALSE)@likelihood # modularity
    even = intereven(imatr) # Interaction eveness
    H2   = H2fun(imatr, H2_integer=TRUE)[[1]]  # specificity
    intas = intasymm(imatr) # Interaction asymmetry
    asymm.c = mean(intas$a.c) # Mean asymmetry for fungi   
    netstat.o = c(conn, nest ,mod, even, H2, asymm.r, asymm.c)
    names(netstat.o) = c("conn","nest", "mod", "even","H2","asymm.r","asymm.c")

  imatrsum = sum(imatr)
    netstat.r = matrix (0, iter, 7) # 7, Matrix to store simulation results for each model
    colnames (netstat.r) = c("conn","nest", "mod","even","H2", "asymm.pla","asymm.pol")
    for (it in 1:iter){
      mit  = mgen(pmatr,imatrsum, zs = F) # Interaction matrix generated by null model
      conn = sum(quant2bin(mit))/(nrow(mit)*ncol(mit)) # Connectance
      nest = nested( mit, method="weighted NODF" ,rescale=FALSE, normalised=TRUE) # Nestedness
      mod  = computeModules(mit,method="Beckett", deep = FALSE, deleteOriginalFiles = TRUE,
             steps = 1000000, tolerance = 1e-10, experimental = FALSE, forceLPA=FALSE)@likelihood 
      H2   = H2fun( mit , H2_integer=TRUE)[[1]]
      even = intereven(mit) # Interaction eveness
      intas = intasymm(mit) # Interaction asymmetry
      asymm.c = mean(intas$a.c) # Mean asymmetry for fungi
      netstat.r[it,] = c(conn, nest, mod, even, H2, asymm.r, asymm.c)
          }
      netstat.ci = confint(netstat.r)
      netstat = data.frame(netstat.o, netstat.ci[,1], netstat.ci[,2])
      names(netstat) = c("observed","ci.lo","ci.up")
      netstat
   }

###-------------------------------------------------------------------------------------------------
##  predict structural metrices of plant-EM network based on propability matrix information
##--------------------------------------------------------------------------------------------------

   netstats.EM.Np.H2 <-  netstats( wEMspmat, randomize=TRUE,iter= 100,pmatr = Np.web.em)
   netstats.EM.Sp.H2 <-  netstats( wEMspmat, randomize=TRUE,iter= 100,pmatr = S.overlap.EM.PL.ses)
   netstats.EM.NSp.H2 <-  netstats( wEMspmat, randomize=TRUE,iter= 100,pmatr = EM.NandS.pweb.ses)
   netstats.EM.null.H2 <-  netstats( wEMspmat, randomize=TRUE,iter= 100,pmatr = NullEM.pweb)
   netstats.EM.obs.H2 <-  netstats( wEMspmat, randomize=TRUE,iter= 100, pmatr = obsEM.pweb)

##--------------------------------------------------------------------------------------------------


##--------------------------------------------------------------------------------------------------
##       The details about these R functions 
##--------------------------------------------------------------------------------------------------
## Function descriptions
## confint -- Confidence intervals of a vector or matrix of simulated values. Used by netstats to calculate 95% confidence intervals of simulated aggregate network statistics.

## intasymm -- Interaction strength asymmetry, calculated following V芍zquez et al. (2007).

## intereven -- Interaction evenness, calculated following Tylianakis et al. (2007).

## mgen -- Matrix generating algorithm used by netstats to generate simulated interaction matrices according to a given probability matrix.

## mlik -- Calculation of multinomial likelihood and AIC according to a given probability matrix. Usage: mlik(imatr, m.p, par), where imatr is the observed interaction matrix, m.p is the probability matrix, and par is the number of parameters used to calculate AIC.

## netstats -- Aggregate network statistics: connectance, nestedness, interaction evenness and mean interaction strength asymmetry for each group in the network. Usage: netstats(imatr, randomize=TRUE, iter=1000, pmatr=NULL), where imatr is the observed interaction matrix and pmatr is the probability matrix used for generating predicted matrices. 
## The R package bipartite (Dormann et al. 2008) is required by "netstats" function to calculate nestedness,modularity and specificity.

## quant2bin -- Transformation of quantitative matrix into binary. Used by netstats to calculate connectance.

##------------------------------------------------------------------------------------------------------------
References

Dormann, C. F., B. Gruber, and J. Fr邦nd. 2008. Introducing the bipartite package: analysing ecological networks. R News 8/2:8每11.
Tylianakis, J. M., T. Tscharntke, and O. T. Lewis. 2007. Habitat modification alters the structure of tropical host-parasitoid food webs. Nature 445:202每205.
V芍zquez, D. P., C. J. Meli芍n, N. M. Williams, N. Bl邦thgen, B. R. Krasnov, and R. Poulin. 2007. Species abundance and asymmetric interaction strength in ecologicalnetworks. Oikos 116:1120每1127.
V芍zquez, D. P., Chaocff, N. P., and Cagnolo, L. (2009). Evaluating multiple determinants of the structure of plant每animal mutualistic networks. Ecology 90, 2039每2046.
##-------------------------------------------------------------------------------------------------------------














