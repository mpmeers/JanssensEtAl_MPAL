## Read in library size files for normalization ## 

frame1_lib<-read.table("allcelltypes.libsizes1")
frame2_lib<-read.table("allcelltypes.libsizes2")
frame3_lib<-read.table("allcelltypes.libsizes3")

## Read in bed files representing merged peaks from all cell types with mapped fragments from each cell line on peaks ##

frame1<-read.table("allcelltypes.counts1.bed", header=TRUE)
frame2<-read.table("allcelltypes.counts2.bed", header=TRUE)
frame3<-read.table("allcelltypes.counts3.bed", header=TRUE)

## Generate normalized counts matrix ##

frame1<-as.data.frame(t(t(frame1[,c(4:ncol(frame1))])*(1000000/frame1_lib$V1)))
frame2<-as.data.frame(t(t(frame2[,c(4:ncol(frame2))])*(1000000/frame2_lib$V1)))
frame3<-as.data.frame(t(t(frame3[,c(4:ncol(frame3))])*(1000000/frame3_lib$V1)))
counts<-data.frame(H1_N=frame2$JFS1339+frame2$JFS1340, H1_C=frame2$JFS1341+frame2$JFS1342, K562_N=frame3$N1_K5+frame3$N2_K5, K562_C=frame3$C1_K5+frame3$C2_K5, CD34_N=frame3$N1_CD34+frame3$N2_CD34, CD34_C=frame3$C1_CD34+frame3$C2_CD34, SEM_N=frame1$JFS1292+frame1$JFS1293, SEM_C=frame1$JFS1294+frame1$JFS1295, RS411_N=frame1$JFS1306+frame1$JFS1307, RS411_C=frame1$JFS1308+frame1$JFS1309, KOPN_N=frame1$JFS1313+frame1$JFS1314, KOPN_C=frame1$JFS1315+frame1$JFS1316, A4_N=frame2$JFS1345+frame2$JFS1346, A4_C=frame2$JFS1347+frame2$JFS1348, A5_N=frame2$JFS1353+frame2$JFS1354, A5_C=frame2$JFS1355+frame2$JFS1356, A6_N=frame2$JFS1361+frame2$JFS1362, A6_C=frame2$JFS1363+frame2$JFS1364, TB11_N=frame2$JFS1369+frame2$JFS1370, TB11_C=frame2$JFS1371+frame2$JFS1372, TB13_N=frame2$JFS1377+frame2$JFS1378, TB13_C=frame2$JFS1379+frame2$JFS1380, A107C_N=frame3$N1_A107C+frame3$N2_A107C, A107C_C=frame3$C1_A107C+frame3$C2_A107C,  A109C_N=frame3$N1_A109C+frame3$N2_A109C, A109C_C=frame3$C1_A109C+frame3$C2_A109C, A384C_N=frame3$N1_A384C+frame3$N2_A384C, A384C_C=frame3$C1_A384C+frame3$C2_A384C)

## Generate log10-transformed N-C scores for counts from peaks called in specific cell types ##

SEM<-read.table("SEM.counts.bed", header=TRUE)
SEM_counts<-SEM[,c(4:ncol(SEM))]
colnames(SEM_counts)<-c("N1", "N2", "C1", "C2")
SEM_counts_norm<-as.data.frame(t(t(SEM_counts)*(1000000/frame1_lib$V1[1:4])))
SEM_NoverC<-data.frame(N=SEM_counts_norm$N1+SEM_counts_norm$N2, C=SEM_counts_norm$C1+SEM_counts_norm$C2)
SEM_NoverC$mean<-(SEM_NoverC$N+SEM_NoverC$C)/2
SEM_NoverC$ratio<-log10((SEM_NoverC$N+min(SEM_NoverC$N[SEM_NoverC$N > 0]))/(SEM_NoverC$C+min(SEM_NoverC$C[SEM_NoverC$C > 0])))*ecdf(abs(log((counts$SEM_N+counts$SEM_C)/2)))((abs(log(SEM_NoverC$N+SEM_NoverC$C)/2)))
SEM_NoverC$length<-SEM$end-SEM$start
RS411<-read.table("RS411.counts.bed", header=TRUE)
RS411_counts<-RS411[,c(4:ncol(RS411))]
colnames(RS411_counts)<-c("N1", "N2", "C1", "C2")
RS411_counts_norm<-as.data.frame(t(t(RS411_counts)*(1000000/frame1_lib$V1[9:12])))
RS411_NoverC<-data.frame(N=RS411_counts_norm$N1+RS411_counts_norm$N2, C=RS411_counts_norm$C1+RS411_counts_norm$C2)
RS411_NoverC$mean<-(RS411_NoverC$N+RS411_NoverC$C)/2
RS411_NoverC$ratio<-log10((RS411_NoverC$N+min(RS411_NoverC$N[RS411_NoverC$N > 0]))/(RS411_NoverC$C+min(RS411_NoverC$C[RS411_NoverC$C > 0])))*ecdf(abs(log((counts$RS411_N+counts$RS411_C)/2)))((abs(log(RS411_NoverC$N+RS411_NoverC$C)/2)))
RS411_NoverC$length<-RS411$end-RS411$start
KOPN<-read.table("KOPN-8.counts.bed", header=TRUE)
KOPN_counts<-KOPN[,c(4:ncol(KOPN))]
colnames(KOPN_counts)<-c("N1", "N2", "C1", "C2")
KOPN_counts_norm<-as.data.frame(t(t(KOPN_counts)*(1000000/frame1_lib$V1[13:16])))
KOPN_NoverC<-data.frame(N=KOPN_counts_norm$N1+KOPN_counts_norm$N2, C=KOPN_counts_norm$C1+KOPN_counts_norm$C2)
KOPN_NoverC$mean<-(KOPN_NoverC$N+KOPN_NoverC$C)/2
KOPN_NoverC$ratio<-log10((KOPN_NoverC$N+min(KOPN_NoverC$N[KOPN_NoverC$N > 0]))/(KOPN_NoverC$C+min(KOPN_NoverC$C[KOPN_NoverC$C > 0])))*ecdf(abs(log((counts$KOPN_N+counts$KOPN_C)/2)))((abs(log(KOPN_NoverC$N+KOPN_NoverC$C)/2)))
KOPN_NoverC$length<-KOPN$end-KOPN$start
H1<-read.table("H1.counts.bed", header=TRUE)
H1_counts<-H1[,c(4:ncol(H1))]
colnames(H1_counts)<-c("N1", "N2", "C1", "C2")
H1_counts_norm<-as.data.frame(t(t(H1_counts)*(1000000/frame2_lib$V1[1:4])))
H1_NoverC<-data.frame(N=H1_counts_norm$N1+H1_counts_norm$N2, C=H1_counts_norm$C1+H1_counts_norm$C2)
H1_NoverC$mean<-(H1_NoverC$N+H1_NoverC$C)/2
H1_NoverC$ratio<-log10((H1_NoverC$N+min(H1_NoverC$N[H1_NoverC$N > 0]))/(H1_NoverC$C+min(H1_NoverC$C[H1_NoverC$C > 0])))*ecdf(abs(log((counts$H1_N+counts$H1_C)/2)))((abs(log(H1_NoverC$N+H1_NoverC$C)/2)))
H1_NoverC$length<-H1$end-H1$start
A4<-read.table("AML-1.counts.bed", header=TRUE)
A4_counts<-A4[,c(4:ncol(A4))]
colnames(A4_counts)<-c("N1", "N2", "C1", "C2")
A4_counts_norm<-as.data.frame(t(t(A4_counts)*(1000000/frame2_lib$V1[5:8])))
A4_NoverC<-data.frame(N=A4_counts_norm$N1+A4_counts_norm$N2, C=A4_counts_norm$C1+A4_counts_norm$C2)
A4_NoverC$mean<-(A4_NoverC$N+A4_NoverC$C)/2
A4_NoverC$ratio<-log10((A4_NoverC$N+min(A4_NoverC$N[A4_NoverC$N > 0]))/(A4_NoverC$C+min(A4_NoverC$C[A4_NoverC$C > 0])))*ecdf(abs(log((counts$A4_N+counts$A4_C)/2)))((abs(log(A4_NoverC$N+A4_NoverC$C)/2)))
A4_NoverC$length<-A4$end-A4$start
A5<-read.table("MPAL-2.counts.bed", header=TRUE)
A5_counts<-A5[,c(4:ncol(A5))]
colnames(A5_counts)<-c("N1", "N2", "C1", "C2")
A5_counts_norm<-as.data.frame(t(t(A5_counts)*(1000000/frame2_lib$V1[9:12])))
A5_NoverC<-data.frame(N=A5_counts_norm$N1+A5_counts_norm$N2, C=A5_counts_norm$C1+A5_counts_norm$C2)
A5_NoverC$mean<-(A5_NoverC$N+A5_NoverC$C)/2
A5_NoverC$ratio<-log10((A5_NoverC$N+min(A5_NoverC$N[A5_NoverC$N > 0]))/(A5_NoverC$C+min(A5_NoverC$C[A5_NoverC$C > 0])))*ecdf(abs(log((counts$A5_N+counts$A5_C)/2)))((abs(log(A5_NoverC$N+A5_NoverC$C)/2)))
A5_NoverC$length<-A5$end-A5$start
A6<-read.table("AML-2.counts.bed", header=TRUE)
A6_counts<-A6[,c(4:ncol(A6))]
colnames(A6_counts)<-c("N1", "N2", "C1", "C2")
A6_counts_norm<-as.data.frame(t(t(A6_counts)*(1000000/frame2_lib$V1[13:16])))
A6_NoverC<-data.frame(N=A6_counts_norm$N1+A6_counts_norm$N2, C=A6_counts_norm$C1+A6_counts_norm$C2)
A6_NoverC$mean<-(A6_NoverC$N+A6_NoverC$C)/2
A6_NoverC$ratio<-log10((A6_NoverC$N+min(A6_NoverC$N[A6_NoverC$N > 0]))/(A6_NoverC$C+min(A6_NoverC$C[A6_NoverC$C > 0])))*ecdf(abs(log((counts$A6_N+counts$A6_C)/2)))((abs(log(A6_NoverC$N+A6_NoverC$C)/2)))
A6_NoverC$length<-A6$end-A6$start
TB11<-read.table("TB11_C1C2N1N2_merged.relaxed.bed.counts.bed", header=TRUE)
TB11_counts<-TB11[,c(4:ncol(TB11))]
colnames(TB11_counts)<-c("N1", "N2", "C1", "C2")
TB11_counts_norm<-as.data.frame(t(t(TB11_counts)*(1000000/frame2_lib$V1[17:20])))
TB11_NoverC<-data.frame(N=TB11_counts_norm$N1+TB11_counts_norm$N2, C=TB11_counts_norm$C1+TB11_counts_norm$C2)
TB11_NoverC$mean<-(TB11_NoverC$N+TB11_NoverC$C)/2
TB11_NoverC$ratio<-log10((TB11_NoverC$N+min(TB11_NoverC$N[TB11_NoverC$N > 0]))/(TB11_NoverC$C+min(TB11_NoverC$C[TB11_NoverC$C > 0])))*ecdf(abs(log((counts$TB11_N+counts$TB11_C)/2)))((abs(log(TB11_NoverC$N+TB11_NoverC$C)/2)))
TB11_NoverC$length<-TB11$end-TB11$start
TB13<-read.table("ALL-1.counts.bed", header=TRUE)
TB13_counts<-TB13[,c(4:ncol(TB13))]
colnames(TB13_counts)<-c("N1", "N2", "C1", "C2")
TB13_counts_norm<-as.data.frame(t(t(TB13_counts)*(1000000/frame2_lib$V1[21:24])))
TB13_NoverC<-data.frame(N=TB13_counts_norm$N1+TB13_counts_norm$N2, C=TB13_counts_norm$C1+TB13_counts_norm$C2)
TB13_NoverC$mean<-(TB13_NoverC$N+TB13_NoverC$C)/2
TB13_NoverC$ratio<-log10((TB13_NoverC$N+min(TB13_NoverC$N[TB13_NoverC$N > 0]))/(TB13_NoverC$C+min(TB13_NoverC$C[TB13_NoverC$C > 0])))*ecdf(abs(log((counts$TB13_N+counts$TB13_C)/2)))((abs(log(TB13_NoverC$N+TB13_NoverC$C)/2)))
TB13_NoverC$length<-TB13$end-TB13$start
A107C<-read.table("AML-3.counts.bed", header=TRUE)
A107C_counts<-A107C[,c(4:ncol(A107C))]
colnames(A107C_counts)<-c("C1", "C2", "N1", "N2")
A107C_counts_norm<-as.data.frame(t(t(A107C_counts)*(1000000/frame3_lib$V1[c(1,6,11,16)])))
A107C_NoverC<-data.frame(N=A107C_counts_norm$N1+A107C_counts_norm$N2, C=A107C_counts_norm$C1+A107C_counts_norm$C2)
A107C_NoverC$mean<-(A107C_NoverC$N+A107C_NoverC$C)/2
A107C_NoverC$ratio<-log10((A107C_NoverC$N+min(A107C_NoverC$N[A107C_NoverC$N > 0]))/(A107C_NoverC$C+min(A107C_NoverC$C[A107C_NoverC$C > 0])))*ecdf(abs(log((counts$A107C_N+counts$A107C_C)/2)))((abs(log(A107C_NoverC$N+A107C_NoverC$C)/2)))
A107C_NoverC$length<-A107C$end-A107C$start
A109C<-read.table("AML-5.counts.bed", header=TRUE)
A109C_counts<-A109C[,c(4:ncol(A109C))]
colnames(A109C_counts)<-c("C1", "C2", "N1", "N2")
A109C_counts_norm<-as.data.frame(t(t(A109C_counts)*(1000000/frame3_lib$V1[c(2,7,12,17)])))
A109C_NoverC<-data.frame(N=A109C_counts_norm$N1+A109C_counts_norm$N2, C=A109C_counts_norm$C1+A109C_counts_norm$C2)
A109C_NoverC$mean<-(A109C_NoverC$N+A109C_NoverC$C)/2
A109C_NoverC$ratio<-log10((A109C_NoverC$N+min(A109C_NoverC$N[A109C_NoverC$N > 0]))/(A109C_NoverC$C+min(A109C_NoverC$C[A109C_NoverC$C > 0])))*ecdf(abs(log((counts$A109C_N+counts$A109C_C)/2)))((abs(log(A109C_NoverC$N+A109C_NoverC$C)/2)))
A109C_NoverC$length<-A109C$end-A109C$start
A384C<-read.table("AML-4.counts.bed", header=TRUE)
A384C_counts<-A384C[,c(4:ncol(A384C))]
colnames(A384C_counts)<-c("C1", "C2", "N1", "N2")
A384C_counts_norm<-as.data.frame(t(t(A384C_counts)*(1000000/frame3_lib$V1[c(3,8,13,18)])))
A384C_NoverC<-data.frame(N=A384C_counts_norm$N1+A384C_counts_norm$N2, C=A384C_counts_norm$C1+A384C_counts_norm$C2)
A384C_NoverC$mean<-(A384C_NoverC$N+A384C_NoverC$C)/2
A384C_NoverC$ratio<-log10((A384C_NoverC$N+min(A384C_NoverC$N[A384C_NoverC$N > 0]))/(A384C_NoverC$C+min(A384C_NoverC$C[A384C_NoverC$C > 0])))*ecdf(abs(log((counts$A384C_N+counts$A384C_C)/2)))((abs(log(A384C_NoverC$N+A384C_NoverC$C)/2)))
A384C_NoverC$length<-A384C$end-A384C$start
CD34<-read.table("CD34.counts.bed", header=TRUE)
CD34_counts<-CD34[,c(4:ncol(CD34))]
colnames(CD34_counts)<-c("C1", "C2", "N1", "N2")
CD34_counts_norm<-as.data.frame(t(t(CD34_counts)*(1000000/frame3_lib$V1[c(4,9,14,19)])))
CD34_NoverC<-data.frame(N=CD34_counts_norm$N1+CD34_counts_norm$N2, C=CD34_counts_norm$C1+CD34_counts_norm$C2)
CD34_NoverC$mean<-(CD34_NoverC$N+CD34_NoverC$C)/2
CD34_NoverC$ratio<-log10((CD34_NoverC$N+min(CD34_NoverC$N[CD34_NoverC$N > 0]))/(CD34_NoverC$C+min(CD34_NoverC$C[CD34_NoverC$C > 0])))*ecdf(abs(log((counts$CD34_N+counts$CD34_C)/2)))((abs(log(CD34_NoverC$N+CD34_NoverC$C)/2)))
CD34_NoverC$length<-CD34$end-CD34$start
K5<-read.table("K562.counts.bed", header=TRUE)
K5_counts<-K5[,c(4:ncol(K5))]
colnames(K5_counts)<-c("C1", "C2", "N1", "N2")
K5_counts_norm<-as.data.frame(t(t(K5_counts)*(1000000/frame3_lib$V1[c(5,10,15,20)])))
K5_NoverC<-data.frame(N=K5_counts_norm$N1+K5_counts_norm$N2, C=K5_counts_norm$C1+K5_counts_norm$C2)
K5_NoverC$mean<-(K5_NoverC$N+K5_NoverC$C)/2
K5_NoverC$ratio<-log10((K5_NoverC$N+min(K5_NoverC$N[K5_NoverC$N > 0]))/(K5_NoverC$C+min(K5_NoverC$C[K5_NoverC$C > 0])))*ecdf(abs(log((counts$K562_N+counts$K562_C)/2)))((abs(log(K5_NoverC$N+K5_NoverC$C)/2)))
K5_NoverC$length<-K5$end-K5$start

## Generate mixture models for cell type-specific N-C score distributions and set threshold for shared vs. N-specific ##

library(mixtools)
SEM_mix<-normalmixEM(SEM_NoverC$ratio, k=2)
SEM_comp<-data.frame(thresh=SEM_NoverC$ratio, SEM_comp1=dnorm(SEM_NoverC$ratio, mean=as.numeric(min(SEM_mix$mu)), sd=as.numeric(SEM_mix$sigma[SEM_mix$mu==min(SEM_mix$mu)]))*SEM_mix$lambda[SEM_mix$mu==min(SEM_mix$mu)], SEM_comp2=dnorm(SEM_NoverC$ratio, mean=as.numeric(max(SEM_mix$mu)), sd=as.numeric(SEM_mix$sigma[SEM_mix$mu==max(SEM_mix$mu)]))*SEM_mix$lambda[SEM_mix$mu==max(SEM_mix$mu)])
SEM_bed<-SEM[,c(1:3)]
SEM_bed$ratio<-SEM_NoverC$ratio
SEM_bed$Nspec<-0
SEM_bed$Nspec[SEM_NoverC$ratio > min(na.omit(apply(SEM_comp, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(SEM_mix$mu))][1])))]<-1
RS411_mix<-normalmixEM(RS411_NoverC$ratio, k=2)
RS411_comp<-data.frame(thresh=RS411_NoverC$ratio, RS411_comp1=dnorm(RS411_NoverC$ratio, mean=as.numeric(min(RS411_mix$mu)), sd=as.numeric(RS411_mix$sigma[RS411_mix$mu==min(RS411_mix$mu)]))*RS411_mix$lambda[RS411_mix$mu==min(RS411_mix$mu)], RS411_comp2=dnorm(RS411_NoverC$ratio, mean=as.numeric(max(RS411_mix$mu)), sd=as.numeric(RS411_mix$sigma[RS411_mix$mu==max(RS411_mix$mu)]))*RS411_mix$lambda[RS411_mix$mu==max(RS411_mix$mu)])
RS411_bed<-RS411[,c(1:3)]
RS411_bed$ratio<-RS411_NoverC$ratio
RS411_bed$Nspec<-0
RS411_bed$Nspec[RS411_NoverC$ratio > min(na.omit(apply(RS411_comp, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(RS411_mix$mu))][1])))]<-1
KOPN_mix<-normalmixEM(KOPN_NoverC$ratio, k=2)
KOPN_comp<-data.frame(thresh=KOPN_NoverC$ratio, KOPN_comp1=dnorm(KOPN_NoverC$ratio, mean=as.numeric(min(KOPN_mix$mu)), sd=as.numeric(KOPN_mix$sigma[KOPN_mix$mu==min(KOPN_mix$mu)]))*KOPN_mix$lambda[KOPN_mix$mu==min(KOPN_mix$mu)], KOPN_comp2=dnorm(KOPN_NoverC$ratio, mean=as.numeric(max(KOPN_mix$mu)), sd=as.numeric(KOPN_mix$sigma[KOPN_mix$mu==max(KOPN_mix$mu)]))*KOPN_mix$lambda[KOPN_mix$mu==max(KOPN_mix$mu)])
KOPN_bed<-KOPN[,c(1:3)]
KOPN_bed$ratio<-KOPN_NoverC$ratio
KOPN_bed$Nspec<-0
KOPN_bed$Nspec[KOPN_NoverC$ratio > min(na.omit(apply(KOPN_comp, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(KOPN_mix$mu))][1])))]<-1
H1_mix<-normalmixEM(H1_NoverC$ratio, k=2)
H1_comp<-data.frame(thresh=H1_NoverC$ratio, H1_comp1=dnorm(H1_NoverC$ratio, mean=as.numeric(min(H1_mix$mu)), sd=as.numeric(H1_mix$sigma[H1_mix$mu==min(H1_mix$mu)]))*H1_mix$lambda[H1_mix$mu==min(H1_mix$mu)], H1_comp2=dnorm(H1_NoverC$ratio, mean=as.numeric(max(H1_mix$mu)), sd=as.numeric(H1_mix$sigma[H1_mix$mu==max(H1_mix$mu)]))*H1_mix$lambda[H1_mix$mu==max(H1_mix$mu)])
H1_bed<-H1[,c(1:3)]
H1_bed$ratio<-H1_NoverC$ratio
H1_bed$Nspec<-0
H1_bed$Nspec[H1_NoverC$ratio > min(na.omit(apply(H1_comp, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(H1_mix$mu))][1])))]<-1
A4_mix<-normalmixEM(A4_NoverC$ratio, k=2)
A4_comp<-data.frame(thresh=A4_NoverC$ratio, A4_comp1=dnorm(A4_NoverC$ratio, mean=as.numeric(min(A4_mix$mu)), sd=as.numeric(A4_mix$sigma[A4_mix$mu==min(A4_mix$mu)]))*A4_mix$lambda[A4_mix$mu==min(A4_mix$mu)], A4_comp2=dnorm(A4_NoverC$ratio, mean=as.numeric(max(A4_mix$mu)), sd=as.numeric(A4_mix$sigma[A4_mix$mu==max(A4_mix$mu)]))*A4_mix$lambda[A4_mix$mu==max(A4_mix$mu)])
A4_bed<-A4[,c(1:3)]
A4_bed$ratio<-A4_NoverC$ratio
A4_bed$Nspec<-0
A4_bed$Nspec[A4_NoverC$ratio > min(na.omit(apply(A4_comp, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(A4_mix$mu))][1])))]<-1
A5_mix<-normalmixEM(A5_NoverC$ratio, k=2)
A5_comp<-data.frame(thresh=A5_NoverC$ratio, A5_comp1=dnorm(A5_NoverC$ratio, mean=as.numeric(min(A5_mix$mu)), sd=as.numeric(A5_mix$sigma[A5_mix$mu==min(A5_mix$mu)]))*A5_mix$lambda[A5_mix$mu==min(A5_mix$mu)], A5_comp2=dnorm(A5_NoverC$ratio, mean=as.numeric(max(A5_mix$mu)), sd=as.numeric(A5_mix$sigma[A5_mix$mu==max(A5_mix$mu)]))*A5_mix$lambda[A5_mix$mu==max(A5_mix$mu)])
A5_bed<-A5[,c(1:3)]
A5_bed$ratio<-A5_NoverC$ratio
A5_bed$Nspec<-0
A5_bed$Nspec[A5_NoverC$ratio > min(na.omit(apply(A5_comp, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(A5_mix$mu))][1])))]<-1
A6_mix<-normalmixEM(A6_NoverC$ratio, k=2)
A6_comp<-data.frame(thresh=A6_NoverC$ratio, A6_comp1=dnorm(A6_NoverC$ratio, mean=as.numeric(min(A6_mix$mu)), sd=as.numeric(A6_mix$sigma[A6_mix$mu==min(A6_mix$mu)]))*A6_mix$lambda[A6_mix$mu==min(A6_mix$mu)], A6_comp2=dnorm(A6_NoverC$ratio, mean=as.numeric(max(A6_mix$mu)), sd=as.numeric(A6_mix$sigma[A6_mix$mu==max(A6_mix$mu)]))*A6_mix$lambda[A6_mix$mu==max(A6_mix$mu)])
A6_bed<-A6[,c(1:3)]
A6_bed$ratio<-A6_NoverC$ratio
A6_bed$Nspec<-0
A6_bed$Nspec[A6_NoverC$ratio > min(na.omit(apply(A6_comp, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(A6_mix$mu))][1])))]<-1
TB11_mix<-normalmixEM(TB11_NoverC$ratio, k=2)
TB11_comp<-data.frame(thresh=TB11_NoverC$ratio, TB11_comp1=dnorm(TB11_NoverC$ratio, mean=as.numeric(min(TB11_mix$mu)), sd=as.numeric(TB11_mix$sigma[TB11_mix$mu==min(TB11_mix$mu)]))*TB11_mix$lambda[TB11_mix$mu==min(TB11_mix$mu)], TB11_comp2=dnorm(TB11_NoverC$ratio, mean=as.numeric(max(TB11_mix$mu)), sd=as.numeric(TB11_mix$sigma[TB11_mix$mu==max(TB11_mix$mu)]))*TB11_mix$lambda[TB11_mix$mu==max(TB11_mix$mu)])
TB11_bed<-TB11[,c(1:3)]
TB11_bed$ratio<-TB11_NoverC$ratio
TB11_bed$Nspec<-0
TB11_bed$Nspec[TB11_NoverC$ratio > min(na.omit(apply(TB11_comp, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(TB11_mix$mu))][1])))]<-1
TB13_mix<-normalmixEM(TB13_NoverC$ratio, k=2)
TB13_comp<-data.frame(thresh=TB13_NoverC$ratio, TB13_comp1=dnorm(TB13_NoverC$ratio, mean=as.numeric(min(TB13_mix$mu)), sd=as.numeric(TB13_mix$sigma[TB13_mix$mu==min(TB13_mix$mu)]))*TB13_mix$lambda[TB13_mix$mu==min(TB13_mix$mu)], TB13_comp2=dnorm(TB13_NoverC$ratio, mean=as.numeric(max(TB13_mix$mu)), sd=as.numeric(TB13_mix$sigma[TB13_mix$mu==max(TB13_mix$mu)]))*TB13_mix$lambda[TB13_mix$mu==max(TB13_mix$mu)])
TB13_bed<-TB13[,c(1:3)]
TB13_bed$ratio<-TB13_NoverC$ratio
TB13_bed$Nspec<-0
TB13_bed$Nspec[TB13_NoverC$ratio > min(na.omit(apply(TB13_comp, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(TB13_mix$mu))][1])))]<-1
A107C_mix<-normalmixEM(A107C_NoverC$ratio, k=2)
A107C_comp<-data.frame(thresh=A107C_NoverC$ratio, A107C_comp1=dnorm(A107C_NoverC$ratio, mean=as.numeric(min(A107C_mix$mu)), sd=as.numeric(A107C_mix$sigma[A107C_mix$mu==min(A107C_mix$mu)]))*A107C_mix$lambda[A107C_mix$mu==min(A107C_mix$mu)], A107C_comp2=dnorm(A107C_NoverC$ratio, mean=as.numeric(max(A107C_mix$mu)), sd=as.numeric(A107C_mix$sigma[A107C_mix$mu==max(A107C_mix$mu)]))*A107C_mix$lambda[A107C_mix$mu==max(A107C_mix$mu)])
A107C_bed<-A107C[,c(1:3)]
A107C_bed$ratio<-A107C_NoverC$ratio
A107C_bed$Nspec<-0
A107C_bed$Nspec[A107C_NoverC$ratio > min(na.omit(apply(A107C_comp, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(A107C_mix$mu))][1])))]<-1
A109C_mix<-normalmixEM(A109C_NoverC$ratio, k=2)
A109C_comp<-data.frame(thresh=A109C_NoverC$ratio, A109C_comp1=dnorm(A109C_NoverC$ratio, mean=as.numeric(min(A109C_mix$mu)), sd=as.numeric(A109C_mix$sigma[A109C_mix$mu==min(A109C_mix$mu)]))*A109C_mix$lambda[A109C_mix$mu==min(A109C_mix$mu)], A109C_comp2=dnorm(A109C_NoverC$ratio, mean=as.numeric(max(A109C_mix$mu)), sd=as.numeric(A109C_mix$sigma[A109C_mix$mu==max(A109C_mix$mu)]))*A109C_mix$lambda[A109C_mix$mu==max(A109C_mix$mu)])
A109C_bed<-A109C[,c(1:3)]
A109C_bed$ratio<-A109C_NoverC$ratio
A109C_bed$Nspec<-0
A109C_bed$Nspec[A109C_NoverC$ratio > min(na.omit(apply(A109C_comp, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(A109C_mix$mu))][1])))]<-1
A384C_mix<-normalmixEM(A384C_NoverC$ratio, k=2)
A384C_comp<-data.frame(thresh=A384C_NoverC$ratio, A384C_comp1=dnorm(A384C_NoverC$ratio, mean=as.numeric(min(A384C_mix$mu)), sd=as.numeric(A384C_mix$sigma[A384C_mix$mu==min(A384C_mix$mu)]))*A384C_mix$lambda[A384C_mix$mu==min(A384C_mix$mu)], A384C_comp2=dnorm(A384C_NoverC$ratio, mean=as.numeric(max(A384C_mix$mu)), sd=as.numeric(A384C_mix$sigma[A384C_mix$mu==max(A384C_mix$mu)]))*A384C_mix$lambda[A384C_mix$mu==max(A384C_mix$mu)])
A384C_bed<-A384C[,c(1:3)]
A384C_bed$ratio<-A384C_NoverC$ratio
A384C_bed$Nspec<-0
A384C_bed$Nspec[A384C_NoverC$ratio > min(na.omit(apply(A384C_comp, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(A384C_mix$mu))][1])))]<-1
CD34_mix<-normalmixEM(CD34_NoverC$ratio, k=2)
CD34_comp<-data.frame(thresh=CD34_NoverC$ratio, CD34_comp1=dnorm(CD34_NoverC$ratio, mean=as.numeric(min(CD34_mix$mu)), sd=as.numeric(CD34_mix$sigma[CD34_mix$mu==min(CD34_mix$mu)]))*CD34_mix$lambda[CD34_mix$mu==min(CD34_mix$mu)], CD34_comp2=dnorm(CD34_NoverC$ratio, mean=as.numeric(max(CD34_mix$mu)), sd=as.numeric(CD34_mix$sigma[CD34_mix$mu==max(CD34_mix$mu)]))*CD34_mix$lambda[CD34_mix$mu==max(CD34_mix$mu)])
CD34_bed<-CD34[,c(1:3)]
CD34_bed$ratio<-CD34_NoverC$ratio
CD34_bed$Nspec<-0
CD34_bed$Nspec[CD34_NoverC$ratio > min(na.omit(apply(CD34_comp, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(CD34_mix$mu))][1])))]<-1
K5_mix<-normalmixEM(K5_NoverC$ratio, k=2)
K5_comp<-data.frame(thresh=K5_NoverC$ratio, K5_comp1=dnorm(K5_NoverC$ratio, mean=as.numeric(min(K5_mix$mu)), sd=as.numeric(K5_mix$sigma[K5_mix$mu==min(K5_mix$mu)]))*K5_mix$lambda[K5_mix$mu==min(K5_mix$mu)], K5_comp2=dnorm(K5_NoverC$ratio, mean=as.numeric(max(K5_mix$mu)), sd=as.numeric(K5_mix$sigma[K5_mix$mu==max(K5_mix$mu)]))*K5_mix$lambda[K5_mix$mu==max(K5_mix$mu)])
K5_bed<-K5[,c(1:3)]
K5_bed$ratio<-K5_NoverC$ratio
K5_bed$Nspec<-0
K5_bed$Nspec[K5_NoverC$ratio > min(na.omit(apply(K5_comp, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(K5_mix$mu))][1])))]<-1

## Generate mixture models for cell type-specific peak length distributions and set threshold for short vs. long peak ##

SEM_mix2<-normalmixEM(log10(SEM_NoverC$length), k=2)
SEM_comp2<-data.frame(thresh=log10(SEM_NoverC$length), SEM_comp1=dnorm(log10(SEM_NoverC$length), mean=as.numeric(min(SEM_mix2$mu)), sd=as.numeric(SEM_mix2$sigma[SEM_mix2$mu==min(SEM_mix2$mu)]))*SEM_mix2$lambda[SEM_mix2$mu==min(SEM_mix2$mu)], SEM_comp2=dnorm(log10(SEM_NoverC$length), mean=as.numeric(max(SEM_mix2$mu)), sd=as.numeric(SEM_mix2$sigma[SEM_mix2$mu==max(SEM_mix2$mu)]))*SEM_mix2$lambda[SEM_mix2$mu==max(SEM_mix2$mu)])
SEM_bed$Nspec2<-0
SEM_bed$Nspec2[log10(SEM_NoverC$length) > min(na.omit(apply(SEM_comp2, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(SEM_mix2$mu))][1])))]<-1
RS411_mix2<-normalmixEM(log10(RS411_NoverC$length), k=2)
RS411_comp2<-data.frame(thresh=log10(RS411_NoverC$length), RS411_comp1=dnorm(log10(RS411_NoverC$length), mean=as.numeric(min(RS411_mix2$mu)), sd=as.numeric(RS411_mix2$sigma[RS411_mix2$mu==min(RS411_mix2$mu)]))*RS411_mix2$lambda[RS411_mix2$mu==min(RS411_mix2$mu)], RS411_comp2=dnorm(log10(RS411_NoverC$length), mean=as.numeric(max(RS411_mix2$mu)), sd=as.numeric(RS411_mix2$sigma[RS411_mix2$mu==max(RS411_mix2$mu)]))*RS411_mix2$lambda[RS411_mix2$mu==max(RS411_mix2$mu)])
RS411_bed$Nspec2<-0
RS411_bed$Nspec2[log10(RS411_NoverC$length) > min(na.omit(apply(RS411_comp2, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(RS411_mix2$mu))][1])))]<-1
KOPN_mix2<-normalmixEM(log10(KOPN_NoverC$length), k=2)
KOPN_comp2<-data.frame(thresh=log10(KOPN_NoverC$length), KOPN_comp1=dnorm(log10(KOPN_NoverC$length), mean=as.numeric(min(KOPN_mix2$mu)), sd=as.numeric(KOPN_mix2$sigma[KOPN_mix2$mu==min(KOPN_mix2$mu)]))*KOPN_mix2$lambda[KOPN_mix2$mu==min(KOPN_mix2$mu)], KOPN_comp2=dnorm(log10(KOPN_NoverC$length), mean=as.numeric(max(KOPN_mix2$mu)), sd=as.numeric(KOPN_mix2$sigma[KOPN_mix2$mu==max(KOPN_mix2$mu)]))*KOPN_mix2$lambda[KOPN_mix2$mu==max(KOPN_mix2$mu)])
KOPN_bed$Nspec2<-0
KOPN_bed$Nspec2[log10(KOPN_NoverC$length) > min(na.omit(apply(KOPN_comp2, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(KOPN_mix2$mu))][1])))]<-1
H1_mix2<-normalmixEM(log10(H1_NoverC$length), k=2)
H1_comp2<-data.frame(thresh=log10(H1_NoverC$length), H1_comp1=dnorm(log10(H1_NoverC$length), mean=as.numeric(min(H1_mix2$mu)), sd=as.numeric(H1_mix2$sigma[H1_mix2$mu==min(H1_mix2$mu)]))*H1_mix2$lambda[H1_mix2$mu==min(H1_mix2$mu)], H1_comp2=dnorm(log10(H1_NoverC$length), mean=as.numeric(max(H1_mix2$mu)), sd=as.numeric(H1_mix2$sigma[H1_mix2$mu==max(H1_mix2$mu)]))*H1_mix2$lambda[H1_mix2$mu==max(H1_mix2$mu)])
H1_bed$Nspec2<-0
H1_bed$Nspec2[log10(H1_NoverC$length) > min(na.omit(apply(H1_comp2, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(H1_mix2$mu))][1])))]<-1
A4_mix2<-normalmixEM(log10(A4_NoverC$length), k=2)
A4_comp2<-data.frame(thresh=log10(A4_NoverC$length), A4_comp1=dnorm(log10(A4_NoverC$length), mean=as.numeric(min(A4_mix2$mu)), sd=as.numeric(A4_mix2$sigma[A4_mix2$mu==min(A4_mix2$mu)]))*A4_mix2$lambda[A4_mix2$mu==min(A4_mix2$mu)], A4_comp2=dnorm(log10(A4_NoverC$length), mean=as.numeric(max(A4_mix2$mu)), sd=as.numeric(A4_mix2$sigma[A4_mix2$mu==max(A4_mix2$mu)]))*A4_mix2$lambda[A4_mix2$mu==max(A4_mix2$mu)])
A4_bed$Nspec2<-0
A4_bed$Nspec2[log10(A4_NoverC$length) > min(na.omit(apply(A4_comp2, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(A4_mix2$mu))][1])))]<-1
A5_mix2<-normalmixEM(log10(A5_NoverC$length), k=2)
A5_comp2<-data.frame(thresh=log10(A5_NoverC$length), A5_comp1=dnorm(log10(A5_NoverC$length), mean=as.numeric(min(A5_mix2$mu)), sd=as.numeric(A5_mix2$sigma[A5_mix2$mu==min(A5_mix2$mu)]))*A5_mix2$lambda[A5_mix2$mu==min(A5_mix2$mu)], A5_comp2=dnorm(log10(A5_NoverC$length), mean=as.numeric(max(A5_mix2$mu)), sd=as.numeric(A5_mix2$sigma[A5_mix2$mu==max(A5_mix2$mu)]))*A5_mix2$lambda[A5_mix2$mu==max(A5_mix2$mu)])
A5_bed$Nspec2<-0
A5_bed$Nspec2[log10(A5_NoverC$length) > min(na.omit(apply(A5_comp2, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(A5_mix2$mu))][1])))]<-1
A6_mix2<-normalmixEM(log10(A6_NoverC$length), k=2)
A6_comp2<-data.frame(thresh=log10(A6_NoverC$length), A6_comp1=dnorm(log10(A6_NoverC$length), mean=as.numeric(min(A6_mix2$mu)), sd=as.numeric(A6_mix2$sigma[A6_mix2$mu==min(A6_mix2$mu)]))*A6_mix2$lambda[A6_mix2$mu==min(A6_mix2$mu)], A6_comp2=dnorm(log10(A6_NoverC$length), mean=as.numeric(max(A6_mix2$mu)), sd=as.numeric(A6_mix2$sigma[A6_mix2$mu==max(A6_mix2$mu)]))*A6_mix2$lambda[A6_mix2$mu==max(A6_mix2$mu)])
A6_bed$Nspec2<-0
A6_bed$Nspec2[log10(A6_NoverC$length) > min(na.omit(apply(A6_comp2, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(A6_mix2$mu))][1])))]<-1
TB11_mix2<-normalmixEM(log10(TB11_NoverC$length), k=3)
TB11_comp2<-data.frame(thresh=log10(TB11_NoverC$length), TB11_comp1=dnorm(log10(TB11_NoverC$length), mean=as.numeric(TB11_mix2$mu[2]), sd=as.numeric(TB11_mix2$sigma[2]))*TB11_mix2$lambda[2], TB11_comp2=dnorm(log10(TB11_NoverC$length), mean=as.numeric(max(TB11_mix2$mu)), sd=as.numeric(TB11_mix2$sigma[TB11_mix2$mu==max(TB11_mix2$mu)]))*TB11_mix2$lambda[TB11_mix2$mu==max(TB11_mix2$mu)])
TB11_bed$Nspec2<-0
TB11_bed$Nspec2[log10(TB11_NoverC$length) > min(na.omit(apply(TB11_comp2, 1, function(x) x[x[3]>x[2] & x[1] > as.numeric(min(TB11_mix2$mu))][1])))]<-1
TB13_mix2<-normalmixEM(log10(TB13_NoverC$length), k=2)
TB13_comp2<-data.frame(thresh=log10(TB13_NoverC$length), TB13_comp1=dnorm(log10(TB13_NoverC$length), mean=as.numeric(min(TB13_mix2$mu)), sd=as.numeric(TB13_mix2$sigma[TB13_mix2$mu==min(TB13_mix2$mu)]))*TB13_mix2$lambda[TB13_mix2$mu==min(TB13_mix2$mu)], TB13_comp2=dnorm(log10(TB13_NoverC$length), mean=as.numeric(max(TB13_mix2$mu)), sd=as.numeric(TB13_mix2$sigma[TB13_mix2$mu==max(TB13_mix2$mu)]))*TB13_mix2$lambda[TB13_mix2$mu==max(TB13_mix2$mu)])
TB13_bed$Nspec2<-0
TB13_bed$Nspec2[log10(TB13_NoverC$length) > min(na.omit(apply(TB13_comp2, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(TB13_mix2$mu))][1])))]<-1
A107C_mix2<-normalmixEM(log10(A107C_NoverC$length), k=2)
A107C_comp2<-data.frame(thresh=log10(A107C_NoverC$length), A107C_comp1=dnorm(log10(A107C_NoverC$length), mean=as.numeric(min(A107C_mix2$mu)), sd=as.numeric(A107C_mix2$sigma[A107C_mix2$mu==min(A107C_mix2$mu)]))*A107C_mix2$lambda[A107C_mix2$mu==min(A107C_mix2$mu)], A107C_comp2=dnorm(log10(A107C_NoverC$length), mean=as.numeric(max(A107C_mix2$mu)), sd=as.numeric(A107C_mix2$sigma[A107C_mix2$mu==max(A107C_mix2$mu)]))*A107C_mix2$lambda[A107C_mix2$mu==max(A107C_mix2$mu)])
A107C_bed$Nspec2<-0
A107C_bed$Nspec2[log10(A107C_NoverC$length) > min(na.omit(apply(A107C_comp2, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(A107C_mix2$mu))][1])))]<-1
A109C_mix2<-normalmixEM(log10(A109C_NoverC$length), k=2)
A109C_comp2<-data.frame(thresh=log10(A109C_NoverC$length), A109C_comp1=dnorm(log10(A109C_NoverC$length), mean=as.numeric(min(A109C_mix2$mu)), sd=as.numeric(A109C_mix2$sigma[A109C_mix2$mu==min(A109C_mix2$mu)]))*A109C_mix2$lambda[A109C_mix2$mu==min(A109C_mix2$mu)], A109C_comp2=dnorm(log10(A109C_NoverC$length), mean=as.numeric(max(A109C_mix2$mu)), sd=as.numeric(A109C_mix2$sigma[A109C_mix2$mu==max(A109C_mix2$mu)]))*A109C_mix2$lambda[A109C_mix2$mu==max(A109C_mix2$mu)])
A109C_bed$Nspec2<-0
A109C_bed$Nspec2[log10(A109C_NoverC$length) > min(na.omit(apply(A109C_comp2, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(A109C_mix2$mu))][1])))]<-1
A384C_mix2<-normalmixEM(log10(A384C_NoverC$length), k=2)
A384C_comp2<-data.frame(thresh=log10(A384C_NoverC$length), A384C_comp1=dnorm(log10(A384C_NoverC$length), mean=as.numeric(min(A384C_mix2$mu)), sd=as.numeric(A384C_mix2$sigma[A384C_mix2$mu==min(A384C_mix2$mu)]))*A384C_mix2$lambda[A384C_mix2$mu==min(A384C_mix2$mu)], A384C_comp2=dnorm(log10(A384C_NoverC$length), mean=as.numeric(max(A384C_mix2$mu)), sd=as.numeric(A384C_mix2$sigma[A384C_mix2$mu==max(A384C_mix2$mu)]))*A384C_mix2$lambda[A384C_mix2$mu==max(A384C_mix2$mu)])
A384C_bed$Nspec2<-0
A384C_bed$Nspec2[log10(A384C_NoverC$length) > min(na.omit(apply(A384C_comp2, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(A384C_mix2$mu))][1])))]<-1
CD34_mix2<-normalmixEM(log10(CD34_NoverC$length), k=2)
CD34_comp2<-data.frame(thresh=log10(CD34_NoverC$length), CD34_comp1=dnorm(log10(CD34_NoverC$length), mean=as.numeric(min(CD34_mix2$mu)), sd=as.numeric(CD34_mix2$sigma[CD34_mix2$mu==min(CD34_mix2$mu)]))*CD34_mix2$lambda[CD34_mix2$mu==min(CD34_mix2$mu)], CD34_comp2=dnorm(log10(CD34_NoverC$length), mean=as.numeric(max(CD34_mix2$mu)), sd=as.numeric(CD34_mix2$sigma[CD34_mix2$mu==max(CD34_mix2$mu)]))*CD34_mix2$lambda[CD34_mix2$mu==max(CD34_mix2$mu)])
CD34_bed$Nspec2<-0
CD34_bed$Nspec2[log10(CD34_NoverC$length) > min(na.omit(apply(CD34_comp2, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(CD34_mix2$mu))][1])))]<-1
K5_mix2<-normalmixEM(log10(K5_NoverC$length), k=2)
K5_comp2<-data.frame(thresh=log10(K5_NoverC$length), K5_comp1=dnorm(log10(K5_NoverC$length), mean=as.numeric(min(K5_mix2$mu)), sd=as.numeric(K5_mix2$sigma[K5_mix2$mu==min(K5_mix2$mu)]))*K5_mix2$lambda[K5_mix2$mu==min(K5_mix2$mu)], K5_comp2=dnorm(log10(K5_NoverC$length), mean=as.numeric(max(K5_mix2$mu)), sd=as.numeric(K5_mix2$sigma[K5_mix2$mu==max(K5_mix2$mu)]))*K5_mix2$lambda[K5_mix2$mu==max(K5_mix2$mu)])
K5_bed$Nspec2<-0
K5_bed$Nspec2[log10(K5_NoverC$length) > min(na.omit(apply(K5_comp2, 1, function(x) x[x[3]>=x[2] & x[1] > as.numeric(min(K5_mix2$mu))][1])))]<-1

## Write threshold bed files ##

write.table(SEM_bed, "SEM_C1C2N1N2_merged.relaxed.NoverCvalues_210209.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(RS411_bed, "RS411_C1C2N1N2_merged.relaxed.NoverCvalues_210209.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(KOPN_bed, "KOPN_C1C2N1N2_merged.relaxed.NoverCvalues_210209.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(H1_bed, "H1_C1C2N1N2_merged.relaxed.NoverCvalues_210209.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(A4_bed, "A4_CcombNcomb_merged.relaxed.NoverCvalues_210209.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(A5_bed, "A5_C1C2N1N2_merged.relaxed.NoverCvalues_210209.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(A6_bed, "A6_C1C2N1N2_merged.relaxed.NoverCvalues_210209.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(TB11_bed, "TB11_C1C2N1N2_merged.relaxed.NoverCvalues_210209.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(TB13_bed, "TB13_C1C2N1N2_merged.relaxed.NoverCvalues_210209.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(A107C_bed, "A107C_C1C2N1N2_merged.relaxed.NoverCvalues_210209.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(A109C_bed, "A109C_C1C2N1N2_merged.relaxed.NoverCvalues_210209.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(A384C_bed, "A384C_C1C2N1N2_merged.relaxed.NoverCvalues_210209.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(CD34_bed, "CD34_C1C2N1N2_merged.relaxed.NoverCvalues_210209.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(K5_bed, "K5_C1C2N1N2_merged.relaxed.NoverCvalues_210209.bed", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

## Assign KMT2A fusion identity ##

SEM_bed$Nspecific<-"No"
SEM_bed$Nspecific[SEM_bed$Nspec==1 & SEM_bed$Nspec2==1]<-"Yes"
RS411_bed$Nspecific<-"No"
RS411_bed$Nspecific[RS411_bed$Nspec==1 & RS411_bed$Nspec2==1]<-"Yes"
KOPN_bed$Nspecific<-"No"
KOPN_bed$Nspecific[KOPN_bed$Nspec==1 & KOPN_bed$Nspec2==1]<-"Yes"
A4_bed$Nspecific<-"No"
A4_bed$Nspecific[A4_bed$Nspec==1 & A4_bed$Nspec2==1]<-"Yes"
A5_bed$Nspecific<-"No"
A5_bed$Nspecific[A5_bed$Nspec==1 & A5_bed$Nspec2==1]<-"Yes"
A6_bed$Nspecific<-"No"
A6_bed$Nspecific[A6_bed$Nspec==1 & A6_bed$Nspec2==1]<-"Yes"
TB11_bed$Nspecific<-"No"
TB11_bed$Nspecific[TB11_bed$Nspec==1 & TB11_bed$Nspec2==1]<-"Yes"
TB13_bed$Nspecific<-"No"
TB13_bed$Nspecific[TB13_bed$Nspec==1 & TB13_bed$Nspec2==1]<-"Yes"
A107C_bed$Nspecific<-"No"
A107C_bed$Nspecific[A107C_bed$Nspec==1 & A107C_bed$Nspec2==1]<-"Yes"
A109C_bed$Nspecific<-"No"
A109C_bed$Nspecific[A109C_bed$Nspec==1 & A109C_bed$Nspec2==1]<-"Yes"
A384C_bed$Nspecific<-"No"
A384C_bed$Nspecific[A384C_bed$Nspec==1 & A384C_bed$Nspec2==1]<-"Yes"
SEM_bed$length<-SEM_bed$end-SEM_bed$start
RS411_bed$length<-RS411_bed$end-RS411_bed$start
KOPN_bed$length<-KOPN_bed$end-KOPN_bed$start
A4_bed$length<-A4_bed$end-A4_bed$start
A5_bed$length<-A5_bed$end-A5_bed$start
A6_bed$length<-A6_bed$end-A6_bed$start
TB11_bed$length<-TB11_bed$end-TB11_bed$start
TB13_bed$length<-TB13_bed$end-TB13_bed$start
A107C_bed$length<-A107C_bed$end-A107C_bed$start
A109C_bed$length<-A109C_bed$end-A109C_bed$start
A384C_bed$length<-A384C_bed$end-A384C_bed$start

## Plot N-C score vs. length for each cell type and save ##

SEMplot<-ggplot(SEM_bed, aes(x=length, y=ratio, color=Nspecific)) + geom_point(size=0.75, pch=16) + scale_x_log10(limits=c(50,150000)) + scale_y_continuous(limits=c(-1.3,2.6)) + scale_color_manual(values=c("4E70EA", "#DF1C1C")) + theme_light() + xlab("KMT2A peak length") + ylab("log10(KMT2A N-terminal/C-terminal ratio)") + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text=element_text(size=14, face="bold"), axis.title=element_text(size=16), legend.text=element_text(size=12), legend.title=element_text(size=14), axis.line=element_line(size=2))
RS411plot<-ggplot(RS411_bed, aes(x=length, y=ratio, color=Nspecific)) + geom_point(size=0.75, pch=16) + scale_x_log10(limits=c(50,150000)) + scale_y_continuous(limits=c(-1.3,2.6)) + scale_color_manual(values=c("4E70EA", "#DF1C1C")) + theme_light() + xlab("KMT2A peak length") + ylab("log10(KMT2A N-terminal/C-terminal ratio)") + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text=element_text(size=14, face="bold"), axis.title=element_text(size=16), legend.text=element_text(size=12), legend.title=element_text(size=14), axis.line=element_line(size=2))
KOPNplot<-ggplot(KOPN_bed, aes(x=length, y=ratio, color=Nspecific)) + geom_point(size=0.75, pch=16) + scale_x_log10(limits=c(50,150000)) + scale_y_continuous(limits=c(-1.3,2.6)) + scale_color_manual(values=c("4E70EA", "#DF1C1C")) + theme_light() + xlab("KMT2A peak length") + ylab("log10(KMT2A N-terminal/C-terminal ratio)") + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text=element_text(size=14, face="bold"), axis.title=element_text(size=16), legend.text=element_text(size=12), legend.title=element_text(size=14), axis.line=element_line(size=2))
A4plot<-ggplot(A4_bed, aes(x=length, y=ratio, color=Nspecific)) + geom_point(size=0.75, pch=16) + scale_x_log10(limits=c(50,150000)) + scale_y_continuous(limits=c(-1.3,2.6)) + scale_color_manual(values=c("4E70EA", "#DF1C1C")) + theme_light() + xlab("KMT2A peak length") + ylab("log10(KMT2A N-terminal/C-terminal ratio)") + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text=element_text(size=14, face="bold"), axis.title=element_text(size=16), legend.text=element_text(size=12), legend.title=element_text(size=14), axis.line=element_line(size=2))
A5plot<-ggplot(A5_bed, aes(x=length, y=ratio, color=Nspecific)) + geom_point(size=0.75, pch=16) + scale_x_log10(limits=c(50,150000)) + scale_y_continuous(limits=c(-1.3,2.6)) + scale_color_manual(values=c("4E70EA", "#DF1C1C")) + theme_light() + xlab("KMT2A peak length") + ylab("log10(KMT2A N-terminal/C-terminal ratio)") + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text=element_text(size=14, face="bold"), axis.title=element_text(size=16), legend.text=element_text(size=12), legend.title=element_text(size=14), axis.line=element_line(size=2))
A6plot<-ggplot(A6_bed, aes(x=length, y=ratio, color=Nspecific)) + geom_point(size=0.75, pch=16) + scale_x_log10(limits=c(50,150000)) + scale_y_continuous(limits=c(-1.3,2.6)) + scale_color_manual(values=c("4E70EA", "#DF1C1C")) + theme_light() + xlab("KMT2A peak length") + ylab("log10(KMT2A N-terminal/C-terminal ratio)") + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text=element_text(size=14, face="bold"), axis.title=element_text(size=16), legend.text=element_text(size=12), legend.title=element_text(size=14), axis.line=element_line(size=2))
TB11plot<-ggplot(TB11_bed, aes(x=length, y=ratio, color=Nspecific)) + geom_point(size=0.75, pch=16) + scale_x_log10(limits=c(50,150000)) + scale_y_continuous(limits=c(-1.3,2.6)) + scale_color_manual(values=c("4E70EA", "#DF1C1C")) + theme_light() + xlab("KMT2A peak length") + ylab("log10(KMT2A N-terminal/C-terminal ratio)") + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text=element_text(size=14, face="bold"), axis.title=element_text(size=16), legend.text=element_text(size=12), legend.title=element_text(size=14), axis.line=element_line(size=2))
TB13plot<-ggplot(TB13_bed, aes(x=length, y=ratio, color=Nspecific)) + geom_point(size=0.75, pch=16) + scale_x_log10(limits=c(50,150000)) + scale_y_continuous(limits=c(-1.3,2.6)) + scale_color_manual(values=c("4E70EA", "#DF1C1C")) + theme_light() + xlab("KMT2A peak length") + ylab("log10(KMT2A N-terminal/C-terminal ratio)") + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text=element_text(size=14, face="bold"), axis.title=element_text(size=16), legend.text=element_text(size=12), legend.title=element_text(size=14), axis.line=element_line(size=2))
A107Cplot<-ggplot(A107C_bed, aes(x=length, y=ratio, color=Nspecific)) + geom_point(size=0.75, pch=16) + scale_x_log10(limits=c(50,150000)) + scale_y_continuous(limits=c(-1.3,2.6)) + scale_color_manual(values=c("4E70EA", "#DF1C1C")) + theme_light() + xlab("KMT2A peak length") + ylab("log10(KMT2A N-terminal/C-terminal ratio)") + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text=element_text(size=14, face="bold"), axis.title=element_text(size=16), legend.text=element_text(size=12), legend.title=element_text(size=14), axis.line=element_line(size=2))
A109Cplot<-ggplot(A109C_bed, aes(x=length, y=ratio, color=Nspecific)) + geom_point(size=0.75, pch=16) + scale_x_log10(limits=c(50,150000)) + scale_y_continuous(limits=c(-1.3,2.6)) + scale_color_manual(values=c("4E70EA", "#DF1C1C")) + theme_light() + xlab("KMT2A peak length") + ylab("log10(KMT2A N-terminal/C-terminal ratio)") + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text=element_text(size=14, face="bold"), axis.title=element_text(size=16), legend.text=element_text(size=12), legend.title=element_text(size=14), axis.line=element_line(size=2))
A384Cplot<-ggplot(A384C_bed, aes(x=length, y=ratio, color=Nspecific)) + geom_point(size=0.75, pch=16) + scale_x_log10(limits=c(50,150000)) + scale_y_continuous(limits=c(-1.3,2.6)) + scale_color_manual(values=c("4E70EA", "#DF1C1C")) + theme_light() + xlab("KMT2A peak length") + ylab("log10(KMT2A N-terminal/C-terminal ratio)") + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text=element_text(size=14, face="bold"), axis.title=element_text(size=16), legend.text=element_text(size=12), legend.title=element_text(size=14), axis.line=element_line(size=2))

ggsave("plots/SEM_lengthvsNoverC_210607.pdf", SEMplot, width=8, height=7)
ggsave("plots/RS411_lengthvsNoverC_210607.pdf", RS411plot, width=8, height=7)
ggsave("plots/KOPN_lengthvsNoverC_210607.pdf", KOPNplot, width=8, height=7)
ggsave("plots/A4_lengthvsNoverC_210607.pdf", A4plot, width=8, height=7)
ggsave("plots/A5_lengthvsNoverC_210607.pdf", A5plot, width=8, height=7)
ggsave("plots/A6_lengthvsNoverC_210607.pdf", A6plot, width=8, height=7)
ggsave("plots/TB11_lengthvsNoverC_210607.pdf", TB11plot, width=8, height=7)
ggsave("plots/TB13_lengthvsNoverC_210607.pdf", TB13plot, width=8, height=7)
ggsave("plots/A107C_lengthvsNoverC_210607.pdf", A107Cplot, width=8, height=7)
ggsave("plots/A109C_lengthvsNoverC_210607.pdf", A109Cplot, width=8, height=7)
ggsave("plots/A384C_lengthvsNoverC_210607.pdf", A384Cplot, width=8, height=7)

