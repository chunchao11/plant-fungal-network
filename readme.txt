## Article title: The effects of species abundance, spatial distribution, and phylogeny on a plant-ectomycorrhizal fungal network
## Authors: Chunchao Zhu, Zihui, Wang, David C. Deane, Wenqi Luo, Yongfa Chen, Yongjun Cao, Yumiao Lin, Minhua Zhang 
## E-mail: zhuchunchao11@163.com; 641246135@qq.com

##--------------------------------------------------------------------------------------------
## Data:
## wEMspmat: a plant-EM fungal interaction matrix including 862 EM fungal species and 43 plant species
##--------------------------------------------------------------------------------------------

##---------------------------------------------------------------------------------------------
#(1)some R codes to generate proprability matrix from abundance, spatial overlap, null and obs

#   Np.web.em: a propability matrix generated from relative abundances of EM and host plants
#   S.overlap.EM.PL.ses: a propability matrix generated from spatial overlap for EM fungi and plants
#   NullEM.pweb: a null probability matrix
#   EM.NandS.pweb.ses: a propability matrix generated from abundance and spatital information
#   obsEM.pweb: an obs probability matrix generated from observed plant-EM matrix

##---------------------------------------------------------------------------------------------- 
## (2) some R functions 
##----------------------------------------------------------------------------------------------
##  These R functions are downloaded from Vázquez's published article (2009)
##  https://figshare.com/collections/Evaluating_multiple_determinants_of_the_structure_of_plant_animal_mutualistic_networks/3301145/1
##----------------------------------------------------------------------------------------------


##--------------------------------------------------------------------------------------------------
##       The details about these R functions 
##--------------------------------------------------------------------------------------------------
## Function descriptions
## confint -- Confidence intervals of a vector or matrix of simulated values. Used by netstats to calculate 95% confidence intervals of simulated aggregate network statistics.

## intasymm -- Interaction strength asymmetry, calculated following Vázquez et al. (2007).

## intereven -- Interaction evenness, calculated following Tylianakis et al. (2007).

## mgen -- Matrix generating algorithm used by netstats to generate simulated interaction matrices according to a given probability matrix.

## mlik -- Calculation of multinomial likelihood and AIC according to a given probability matrix. Usage: mlik(imatr, m.p, par), where imatr is the observed interaction matrix, m.p is the probability matrix, and par is the number of parameters used to calculate AIC.

## netstats -- Aggregate network statistics: connectance, nestedness, interaction evenness and mean interaction strength asymmetry for each group in the network. Usage: netstats(imatr, randomize=TRUE, iter=1000, pmatr=NULL), where imatr is the observed interaction matrix and pmatr is the probability matrix used for generating predicted matrices. 
## The R package bipartite (Dormann et al. 2008) is required by "netstats" function to calculate nestedness,modularity and specificity.

## quant2bin -- Transformation of quantitative matrix into binary. Used by netstats to calculate connectance.

##------------------------------------------------------------------------------------------------------------
References

Dormann, C. F., B. Gruber, and J. Fründ. 2008. Introducing the bipartite package: analysing ecological networks. R News 8/2:8–11.
Tylianakis, J. M., T. Tscharntke, and O. T. Lewis. 2007. Habitat modification alters the structure of tropical host-parasitoid food webs. Nature 445:202–205.
Vázquez, D. P., C. J. Melián, N. M. Williams, N. Blüthgen, B. R. Krasnov, and R. Poulin. 2007. Species abundance and asymmetric interaction strength in ecologicalnetworks. Oikos 116:1120–1127.
Vázquez, D. P., Chaocff, N. P., and Cagnolo, L. (2009). Evaluating multiple determinants of the structure of plant–animal mutualistic networks. Ecology 90, 2039–2046.
##-------------------------------------------------------------------------------------------------------------



