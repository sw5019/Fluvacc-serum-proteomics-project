
##### Function required to do normalization
BANDWIDTH=5
MAD.THRES=5

gauss.kernel.smooth = function(xx, yy, k.width, outlier=TRUE) {

  yest = yy
  n = length(yy)
  medy = median(yy, na.rm=TRUE)
  mady = mad(yy, na.rm=TRUE)
  oid = abs(yy - medy) / mady > MAD.THRES
  oid[is.na(oid)] = FALSE
  yy[oid] = NA
  for(i in 1:n) {
    wt = (xx - xx[i]) / k.width
    wt = dnorm(wt, 0, 1)
    wt[is.na(yy)] = NA
    yest[i] = sum(wt * yy, na.rm=TRUE) / sum(wt, na.rm=TRUE)
  }
  yest
}

sdata = read.delim("Samples_run time_and_order_infomation.txt", header=T, as.is=T, check.names=F)
ord = sdata$`Actual Run Order`
sdata = sdata[ord, ]

edata = read.delim("Spectronaut_custom_report.xls", header=T, as.is=T, check.names=F)
edata = edata[edata$F.NormalizedPeakArea > 1, ]  ### There are a lot of "1"s in the entries


##### First make sure all sample names are in the quant data
cid = sdata$FileName %in% edata$R.FileName
sdata = sdata[cid, ]

##### Re-organize quant data into a spreadsheet
id = paste(edata$PG.ProteinAccessions, edata$EG.PrecursorId, edata$F.FrgIon)
uid = unique(id)
nid = length(uid)
rid = match(id, uid)

usp = unique(edata$R.FileName)
nsp = length(usp)
cid = match(edata$R.FileName, usp)

mdata = matrix(NA, nid, nsp)

for(i in 1:nrow(edata)) {
  if(i %% 100000 == 0) print(i)
  x = rid[i]
  y = cid[i]
  mdata[x,y] = edata$F.NormalizedPeakArea[i]
}
colnames(mdata) = usp

prot = pep = frag = rep(NA, nid)
for(i in 1:nid) {
  tmp = strsplit(uid[i], " ")[[1]]
  prot[i] = tmp[1]
  pep[i] = tmp[2]
  frag[i] = tmp[3]
}

mdata = data.frame(Protein=prot, 
                   Peptide=pep,
                   Fragment=frag,
                   mdata, stringsAsFactors=FALSE, check.names=FALSE)

cc = apply(mdata, 2, function(x) sum(is.na(x)))
mdata = mdata[ , cc<=40000]
sdata = sdata[sdata$FileName %in% colnames(mdata), ]

mm = match(sdata$FileName, colnames(mdata))
mdata = mdata[,c(1:3,mm)]

dir.create("cleanData")
write.table(mdata, "cleanData/matrix_data.txt", sep="\t", quote=F, row.names=F)
write.table(sdata, "cleanData/meta_data.txt", sep="\t", quote=F, row.names=F)

rm(edata)


################# Missing data criteria

ct = apply(mdata[,-c(1:3)], 1, function(x) mean(is.na(x)))
mdata = mdata[ct <= 0.2, ]
mdata[mdata == 0] = NA

########################### Time trend and batch correction
tmp = mdata[,-c(1:3)]
tmp = log2(tmp)
med = apply(tmp, 2, function(x) median(x, na.rm=TRUE))
meanmed = mean(med, na.rm=TRUE)
med = med - meanmed
tmp = sweep(tmp, 2, med)
pdf("boxplot.pdf", height=6, width=18)
boxplot(tmp, las=2, cex=0.1, cex.axis=0.3)
dev.off()
mdata[,-c(1:3)] = 2^tmp
rm(tmp)


############################ Time trend and Batch correction
tmpdat = as.matrix(mdata[,-c(1:3)])
tmpdat.clean = tmpdat
nmol = nrow(tmpdat)
ubatch = unique(sdata$Batch)
nbatch = length(ubatch)

for(j in 1:nmol) {
  
  ### Smoothing of time trends within each batch
  ### by subtracting off Gaussian kernel regression (of bandwidth covering X samples)
  if(j %% 1000 == 0) print(j)
  
  b.med = rep(NA, nbatch)
  for(b in 1:nbatch) {
    wid = which(sdata$Batch == ubatch[b])
    x = sdata$`Actual Run Order`[wid]
    tmp.y = as.numeric(tmpdat[j,wid])
    tmp.y[tmp.y <= 0] = NA
    y = log2(tmp.y)
    yest = gauss.kernel.smooth(x, y, BANDWIDTH)  ### Kernel Width = 5 samples
    yest = yest - median(yest)
    y.corr = y - yest
    tmpdat.clean[j,wid] = 2^y.corr
    b.med[b] = 2^median(y.corr, na.rm=TRUE)
  }
  
  ### Take out the batch-wise median shifts
  tmp = as.numeric(tmpdat.clean[j,])
  tmp[tmp <= 0] = NA
  tmp = log2(tmp)
  
  ### First identify non-outlier data points 
  use.obs = rep(TRUE, length(tmp))
  for(b in 1:nbatch) {
    bid = sdata$Batch == b 
    tmp.x = tmp[bid]
    tmp.use = rep(TRUE, sum(bid))
    if(any(!is.na(tmp.x))) {
      medy = median(tmp.x, na.rm=TRUE)
      mady = mad(tmp.x, na.rm=TRUE)
      oid = abs(tmp.x - medy) / mady > MAD.THRES
      tmp.use[oid] = FALSE
      use.obs[bid] = tmp.use
    }
  }
  
  ### Execute this only if there are >1 batches
  if(nbatch > 1) {
    tmp.fit = lm(tmp ~ factor(sdata$Batch), subset=use.obs)
    var.names = names(coef(tmp.fit))
    gg = grep("Batch", var.names)
    nn = length(gg)
    for(k in 1:nn) {
      batch.num = as.numeric(strsplit(var.names[gg[k]], ")")[[1]][2])
      batch.coef = coef(tmp.fit)[gg[k]]
      if(!is.na(batch.coef)) {
        bid = sdata$Batch == batch.num
        tmp[bid] = tmp[bid] - batch.coef
      }
    }
  }
  
  ### Re-exponentiate to original peak area scale
  #tmp.min = min(tmp[tmp > 0], na.rm=TRUE)
  #tmp[tmp <= 0] = tmp.min / 2
  tmp = 2^tmp
  tmpdat.clean[j,] = tmp
}


par(mfrow=c(3,3))
for(k in 1:nbatch) {
  ks = sample(1:nmol, 1)
  plot(as.numeric(mdata[ks,-c(1:3)]), as.numeric(tmpdat.clean[ks,]), cex=.2, log="xy", main=ks)
}

### Exponentiate back
tmpdat.clean = round(tmpdat.clean, 2)
mdata[,-c(1:3)] = tmpdat.clean

any(mdata <= 0, na.rm=TRUE)

write.table(mdata, "normalized_data.txt", sep="\t", quote=F, row.names=F)

dir.create("mapDIA")
write.table(mdata, "mapDIA/normalized_data.txt", sep="\t", quote=F, row.names=F)


#pdf("boxplots.pdf", height=12, width=12)
#boxplot(log2(mdata[,-c(1:3)]), las=2)
#dev.off()

#length(unique(d1u$Protein))
length(unique(mdata$Protein))


