library(ggplot2)
library(data.table)
library(p.exact)
per_chrom <- function(x){ # function to find per chromosome position
  if (x <= 34964571){
    as.numeric(x)
  }
  else if(x <= 57002136){
    as.numeric(x - 34964571)
  }
  else if(x <= 82501170){
    as.numeric(x - 57002136)
  }
  else if(x <= 103363881){
    as.numeric(x - 82501170)
  }
  else if(x > 103363881){
    as.numeric(x - 103363881)
  }
}

data(arab) # load data
pheno <- phdata(arab) # phenotype data
geno <- gtdata(arab) # genotype data
map <- geno@map # map positions
chrom <- geno@chromosome # chromosome positions (1-5)
snp_names <- geno@snpnames # snp names by index
per_chrom_map <- lapply(map, per_chrom) # map position by chromosome
#snp_geno_df <- as.genotype.snp.data(geno) # convert genotype data to dataframe
write.csv(pheno, file="Arabidopsis_pheno.csv") # save data files
write.csv(snp_geno_df, file="Arabidopsis_geno.txt")
sum_vec <- 0
ptm1 <- proc.time()
convert <- function(x) as.numeric(factor(x, levels = names(sort(-table(x)))))

dat <- fread("Arabidopsis_geno.txt", sep=",")
dat <- as.data.frame(lapply(dat, FUN = convert))
res <- 0
pheno <- fread("Arabidopsis_pheno.txt", sep=",")[["X47_0W_GH_FT"]] #dat[,length(dat[1,])]
for (i in 2:(length(dat[1,]))){
su <- 0
geno <- dat[,i]
dim(geno)
ave_0 <- mean( pheno[geno == 1], na.rm=TRUE) #changed to 1 from 0
ave_1 <- mean( pheno[geno == 2], na.rm=TRUE) #changed to 2 from 1
ord <- c(ave_0, ave_1)
ord <- rank(ord)
ref <- c(1,2) #changed from (0,1)
who_len <- length(geno)
min_len <- length( which( geno == ref[ which(ord == 1) ]) )
max_len <- length( which( geno == ref[ which(ord == 2) ]) )
i_list <- sort(pheno,decreasing=F)
top_g0 <- i_list[1 : min_len]
top_g1 <- i_list[(min_len + 1) : (min_len + max_len) ]
r_list <- sort( pheno,decreasing=T)
bot_g0 <- r_list[1 : min_len]
bot_g1 <- r_list[(min_len + 1) : (min_len + max_len) ]
top_s0 <- sum(top_g0, na.rm=TRUE)
top_s1 <- sum(top_g1, na.rm=TRUE)
bot_s0 <- sum(bot_g0, na.rm=TRUE)
bot_s1 <- sum(bot_g1, na.rm=TRUE)
obs_s0 <- sum( pheno[which( geno == ref[ which(ord == 1) ])], na.rm=TRUE)
obs_s1 <- sum( pheno[which( geno == ref[ which(ord == 2) ])], na.rm=TRUE)
su <- obs_s0 + obs_s1
x2 <- (su*log(obs_s1) - obs_s1) - (su*log(bot_s1) - bot_s1)
x1 <- (su*log(top_s1) - top_s1) - (su*log(bot_s1) - bot_s1)
res_1 <- x2/x1
if(is.na(res_1) | res_1 < 0 | res_1 == Inf){
  res[i-1] <- 0
}
else {
  res[i-1] <- res_1
}
}

sum_vec <- sum_vec + res

ptm2 <- proc.time()
cat(ptm2-ptm1)

HA_data <- as.data.frame(cbind(snp_names, chrom, map, per_chrom_map, res))
plot(map, res, col=chrom, xlab="Map Position", ylab="HA-coefficient", main="HA Coefficient GWAS of Phenotype 0W_GH_FT")

colnames(HA_data) <- c("ID", "Chromosome", "Map","Per_chrom Map", "Res")
write.csv(HA_data, "HA_res_SD.csv")
HA_data <- fread("HA_res_SD.csv", sep=",")
HA_data <- as.data.frame(lapply(HA_data, FUN = as.numeric))
HA_no_zeros <- HA_data[(HA_data$Res > 0), ]

ord_snps <- order(HA_data$Res, decreasing=TRUE)
top_hit_indexes <- ord_snps[1:10]
top_hits <- as.data.frame(cbind(as.numeric(as.matrix(HA_data$Chromosome[top_hit_indexes],ncol=1)), as.numeric(as.matrix(per_chrom_map[top_hit_indexes], ncol=1))))
colnames(top_hits) <- c("Chrom", "Pos")
range <- 500
top_hits$search_string <- rep("",10)
for (i in 1:10){
  top_hits$search_string[i] <- paste0("Chr", top_hits$Chrom[i], ":", round(top_hits$Pos[i]) - range, "..", round(top_hits$Pos[i]) + range)
}
top_hits
write.csv(top_hits, "Top_genes_HA_coed_0W_GH_FT.csv")