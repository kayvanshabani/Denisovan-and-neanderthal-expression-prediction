#####Creating plink files to run faster by choosing proper range
for (i in c(22:1)) {
  if(!dir.exists(paste("plinks/plink_", i, sep = "")))
    dir.create(paste("plinks/plink_", i, sep = ""))
  argplink = paste("plink/plink --vcf ", "ALL.chr", i , ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --make-bed --out ", "plinks/plink_", i, "/ --chr ", i, " --keep ", "finalselected", sep = "")
  print(argplink)
  system(argplink , ignore.stdout=T,ignore.stderr=T)
}

#####filterng and training

fs = fread("finalselected")
View(fs)
nrow(isosamp)
nrow(filter(isosamp, expressionmean > 0.05))
iss = fread("/PROJECTS/Kayvan/modified_isotrans.csv")
iss = na.omit(iss)
mean(iss$expressionmean)
mean(iss$expressionvar)
isosamp = fread("/PROJECTS/Kayvan/modified_isotrans.csv")
isosamp = na.omit(isosamp)
isosamp = filter(isosamp, zeropercent < 0.1)
isosamp = filter(isosamp, expressionmean > 0.05)
isosamp = filter(isosamp, expressionmean != 1)
for(i in c(1:nrow(isosamp))){
  a = Sys.time()
  x = as.data.frame(isosamp[i,])
  if(!dir.exists(paste("transcripts/", x[1,7], sep = ""))){
    dir.create(paste("transcripts/", x[1,7], sep = ""))
    if (!file.exists(paste("transcripts/", x[1,7], "/plinkextract_", x[1,7], ".txt/", sep = ""))) {
      extract = data.frame(
        n1 = substr(x[1,10], 4, nchar(x[1,10])),
        n2 = as.integer(x[1,11])-500000,
        n3 = x[1,12]+500000,
        n4 = 1
      )
      fwrite(extract, paste("transcripts/", x[1,7], "/plinkextract_", x[1,7], ".txt", sep = ""), col.names = FALSE, sep = " ")
    }
    argplink = paste("plink/plink --bfile plinks/plink_", substr(x[1,10], 4, nchar(x[1,10])), "/", " --extract range transcripts/",x[1,7] ,"/plinkextract_", x[1,7], ".txt --make-bed --out ", "transcripts/", x[1,7], "/", sep = "")
    print(argplink)
    system(argplink , ignore.stdout=T,ignore.stderr=T)
    print("plink done")
    famdf = data.frame(
      n1 = fs[,1],
      n2 = fs[,2],
      n3 = 0,
      n4 = 0,
      n5 = 0,
      n6 = 0
    )
    famdf$n6 = as.character(x[1, 13:ncol(x)])
    argcleanfam = paste("rm -f", " transcripts/", x[1,7], "/.fam", sep = "")
    system(argcleanfam , ignore.stdout=T,ignore.stderr=T)
    print("cleanfam done")
    fwrite(famdf, paste("transcripts/", x[1,7], "/.fam", sep = ""), col.names = FALSE, sep = " ")
    argfusion = paste("Rscript FUSION.compute_weights.R --bfile ","transcripts/", x[1,7] ,"/ --tmp tmp/", " --out ", "transcripts/", x[1,7], "/ --models enet --PATH_gcta gcta64 --PATH_plink plink/plink --verbose 2 --noclean TRUE --hsq_p 1", sep = "")
    system(argfusion , ignore.stdout=T,ignore.stderr=T)
    print("fusion done")
  }
  print(i)
  print(Sys.time()-a)
}
