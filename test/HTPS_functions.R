# read the file peptide.txt
preparing_input<-function (input, debugs_plot=FALSE) {
  tmp<-read.csv(input, header = T, stringsAsFactors = F, sep ="\t")
  remove_con_decoy<- function(x) {grep("CON\\_\\_|REV\\_\\_", x[,35], invert=TRUE) }
  # score >40
  tmp<-subset(tmp, tmp$Score>40)
  # remove PEP <0.05
  tmp<-subset(tmp, tmp$PEP < 0.05)
  if (debugs_plot == TRUE) {
    plot(sort(tmp$Score))
    plot(sort(tmp$PEP))
  }
  # remove decoy and reverse
  tmp<-tmp[remove_con_decoy(tmp),]
  return(tmp)}

# split
split_cleavage<-function (input, intename, debugs_plot=FALSE) {
  tog = paste(c('N.term.cleavage.window', 'C.term.cleavage.window', intename), collapse="|")
  tmp<-input[,grepl(tog, colnames(input))]
  tmp<-data.frame(subset(tmp,tmp[intename]>0))
  if (debugs_plot == TRUE) {
    plot(sort(unlist(log10(tmp[intename]))))
  }
  tmp<-data.frame(peptide =c(tmp[,"N.term.cleavage.window"], tmp[,"C.term.cleavage.window"]))
  tmp<-data.frame(unique(tmp$peptide))
  tmp<-subset(tmp, tmp$unique.tmp.peptide.!="________________")
  return (tmp)
}


cbind.fill <- function(...){
  nm <- list(...)
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow))
  do.call(cbind, lapply(nm, function (x)
    rbind(x, matrix(, n-nrow(x), ncol(x)))))
}

# change location X.csv file
start_fun<-function (input, aa_pos) {
  tmp_a<-read.csv(aa_pos, header = F, stringsAsFactors = F, sep =",")
  colnames(tmp_a)<-c('cleavage')
  tmp1<-data.frame(input[,1], stringsAsFactors = F)
  tmp1<-data.frame(tmp1[!apply(tmp1 == "", 1, all),], stringsAsFactors = F)
  colnames(tmp1)<-c('cleavage')
  tmp2<-data.frame(input[,2], stringsAsFactors = F)
  tmp2<-data.frame(tmp2[!apply(tmp2 == "", 1, all),], stringsAsFactors = F)
  colnames(tmp2)<-c('cleavage')
  tmp3<-data.frame(input[,3], stringsAsFactors = F)
  tmp3<-data.frame(tmp3[!apply(tmp3 == "", 1, all),], stringsAsFactors = F)
  colnames(tmp3)<-c('cleavage')
  tmp1<-rbind(tmp1,tmp_a)
  tmp2<-rbind(tmp2,tmp_a)
  tmp3<-rbind(tmp3,tmp_a)
  tmp<-list(tmp1, tmp2, tmp3)
  return (tmp)
}

separate_fun<- function (input) {
  tmp<-sapply(input, as.character)
  tmp<-data.frame(strsplit(tmp, ""))
  tmp<-data.frame(lapply(tmp, function (x) gsub("_", "X", x)))
  tmp<-data.frame(t(tmp))
  return(tmp)
}

specificity_function<-function (input) {
  tmp<-apply(input, 2, function (x) {x/(sum(as.matrix(input)))})
  tmp<-tmp[-which(rownames(tmp) %in% "X"),]
  require(reshape)
  tmp<-melt(tmp)
  colnames(tmp)<-c("AA", "position", "value")
  tmp$ID<-paste(tmp$AA, tmp$position, sep ="_")
  return(tmp)
}


create_decoy<- function(input, db, seed) {
  set.seed(seed)
  tmp<-data.frame(matrix(sample(as.character(unlist(strsplit(db, ""))),nrow(input)*ncol(input)),nrow(input),ncol(input)))
  return (tmp)}


specificity_function_decoy<- function(input){
  tmp<-apply(input, 2, function (x) {x/(sum(as.matrix(input)))})
  tmp<-melt(tmp)
  colnames(tmp)<-c("AA", "position", "value")
  tmp$ID<-paste(tmp$AA, tmp$position, sep ="_")
  return(tmp)
}

applyBy <- function(x, by, fun, ...)
{
  if (length(by) == 1)
  {
    nc <- ncol(x)
    split.index <- rep(1:ceiling(nc / by), each = by, length.out = nc)
  } else
  {
    nc <- length(by)
    split.index <- by
  }
  index.list <- split(seq(from = 1, to = nc), split.index)


  sapply(index.list, function(i)
  {
    do.call(fun, list(x[, i], ...))
  })
}

calc_stats<-function (input, debugs_plot=FALSE) {
  tmp_2<-data.frame(applyBy(input[,c(2:7)], 3, apply, 1, mean))
  tmp<-cbind(input, tmp_2)
  tmp$FC<-tmp$X1-tmp$X2
  t.res<-apply(tmp[,c(2:7)], 1, function (x) t.test(x[1:3],x[4:6],paired=FALSE, var.equal = TRUE, alternative =  "two.sided"))
  tmp$pval<-unlist(lapply(t.res, function (x) x$p.val))
  tmp$pval_adj<-p.adjust(p= tmp$pval, method = "BH")
  tmp$FC_cor<-ifelse(tmp$pval_adj < 0.01, tmp$FC, 0)
  tmp$FC_cor_pos<-ifelse(tmp$FC_cor > 0, tmp$FC_cor, 0)
  if (debugs_plot==TRUE) {
    plot(sort(tmp$FC_cor_pos))
    plot(x= tmp$FC, y= -log10(tmp$pval_adj))
  }
  return(tmp)
}

# change location distribution_AA.csv file
entropy_function<-function (input) {
  tmp<-input
  tmp<-tmp-1
  dist_AA<-read.csv("distribution_AA.csv")
  tmp<-data.frame(apply(tmp, 2, function (x) {x/dist_AA$value}))
  tmp<-data.frame(apply(tmp, 2, function (x) {(x)/sum(x)}))
  tmp<-(log(tmp,20)*tmp)
  tmp<-data.frame(apply(tmp, 2, function(x) {sum(x, na.rm = T)*(-1)}))
  return(tmp)}

block_entropy_function<-function(input, debugs_plot=FALSE) {
  tmp<-input
  tmp_1<-data.frame(0)
  tmp_1$B1<-tmp[8,7]
  tmp_1$B2<-tmp[7,7]+tmp[8,7]
  tmp_1$B3<-tmp[6,7]+tmp[7,7]+tmp[8,7]
  tmp_1$B3<-tmp[5,7]+tmp[6,7]+tmp[7,7]+tmp[8,7]
  tmp_1$B4<-tmp[4,7]+tmp[5,7]+tmp[6,7]+tmp[7,7]+tmp[8,7]
  tmp_1$B1_1<-tmp[9,7]
  tmp_1$B2_1<-tmp[9,7]+tmp[10,7]
  tmp_1$B3_1<-tmp[9,7]+tmp[10,7]+tmp[11,7]
  tmp_1$B3_1<-tmp[9,7]+tmp[10,7]+tmp[11,7]+tmp[12,7]
  tmp_1$B4_1<-tmp[9,7]+tmp[10,7]+tmp[11,7]+tmp[12,7]+tmp[13,7]
  tmp_1<-tmp_1[,c(5,4,3,2,6,7,8,9)]

  tmp_2<-data.frame(0)
  tmp_2$B1<-tmp[8,8]
  tmp_2$B2<-tmp[7,8]+tmp[8,8]
  tmp_2$B3<-tmp[6,8]+tmp[7,8]+tmp[8,8]
  tmp_2$B3<-tmp[5,8]+tmp[6,8]+tmp[7,8]+tmp[8,8]
  tmp_2$B4<-tmp[4,8]+tmp[5,8]+tmp[6,8]+tmp[7,8]+tmp[8,8]
  tmp_2$B1_1<-tmp[9,8]
  tmp_2$B2_1<-tmp[9,8]+tmp[10,8]
  tmp_2$B3_1<-tmp[9,8]+tmp[10,8]+tmp[11,8]
  tmp_2$B3_1<-tmp[9,8]+tmp[10,8]+tmp[11,8]+tmp[12,8]
  tmp_2$B4_1<-tmp[9,8]+tmp[10,8]+tmp[11,8]+tmp[12,8]+tmp[13,8]
  tmp_2<-tmp_2[,c(5,4,3,2,6,7,8,9)]

  tmp_3<-rbind(tmp_1,tmp_2)
  tmp_4<-tmp_3[2,]-tmp_3[1,]
  if (debugs_plot==TRUE) {
    plot(t(tmp_1))
    plot(t(tmp_2))
    plot(t(tmp_4))
  }
  tmp_3<-rbind(tmp_3, tmp_4)
  rownames(tmp_3)<-c("ctrl", "protease", "DELTA")
  return(tmp_3)

}

specificity_entropy_tda<-function(cols, HTPS_DB) {
  df <- cols[!is.na(cols)]
  table_cleavages<-separate_fun(df)
  table_cleavages_table<-data.frame(apply(table_cleavages,2, table))
  table_cleavages_spec<-specificity_function(table_cleavages_table)
  # decoy
  decoy_cleavages<- create_decoy(table_cleavages, HTPS_DB,123)
  decoy_cleavages_table<-data.frame(apply(decoy_cleavages,2, table))
  decoy_cleavages_spec<-specificity_function_decoy(decoy_cleavages_table)

  table_cleavages_table<-table_cleavages_table[-which(rownames(table_cleavages_table) %in% "X"),]
  entropy_1<-entropy_function(table_cleavages_table)
  entropy_decoy1<-entropy_function(decoy_cleavages_table)

  comb = list("targ_spec" = table_cleavages_spec, "dec_spec" = decoy_cleavages_spec, 'targ_entr'=entropy_1, 'dec_entr'=entropy_decoy1)
  return (comb)
}
