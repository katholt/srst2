# input = ST calls (including ?, *, NF)
# strip of ? and * and return (to aid aggregation in barplots)
simplifyST<-function(x) {
x<-gsub("[*]", "", x)
x<-gsub("[?]", "", x)
x
}

# calculate distance matrix indicating number of non-shared alleles between STs
# input = data matrix with samples in rows, MLST allele numbers in columns
stDist<-function(m) {
d_m<-matrix(0,ncol=nrow(m),nrow=nrow(m))
for (i in 1:(nrow(m)-1)){
for (j in (i+1):nrow(m)) {
dij<-sum(m[i,]!=m[j,])
d_m[i,j]<-dij
d_m[j,i]<-dij
}
}
colnames(d_m)<-rownames(m)
rownames(d_m)<-rownames(m)
d_m
}

# cluster data matrix on the basis of MLST data
# input = data matrix with samples in rows, MLST allele numbers in columns
clusterByST<-function(m) {
hclust(as.dist(stDist(m)),"single")
}

# extract ST profiles
# input = matrix of mlst data; samples in rows, ST in column 1, alleles in other columns
# return = matrix of STs encountered; one row per ST, ST labels in rownames, alleles in columns, sorted by ST
stProfiles<-function(m, suppressSNPs=F, suppressUncertainty=F) {
if (suppressSNPs) {
for (j in 1:ncol(m)) { m[,j]<-gsub("[*]", "", m[,j]) }
}
if (suppressUncertainty) {
for (j in 1:ncol(m)) { m[,j]<-gsub("[?]", "", m[,j]) }
}
st_profiles<-unique(data.frame(m))
st_list <- character()
for (i in 1:nrow(st_profiles)){
if (st_profiles[i,1] != "NF" & st_profiles[i,1] != "NF*" & st_profiles[i,1] != "NF*?" & st_profiles[i,1] != "NF?*") {
st<-paste0("ST",st_profiles[i,1],collapse=NULL)
}
else {st <- st_profiles[i,1]}
profile<-paste(st_profiles[i,-1],collapse="_")
rownames(st_profiles)[i] <- paste(st,profile)
st_list[i] <- st
}
st_labels <- character()
for (i in 1:nrow(m)){
if (m[i,1] != "NF" & m[i,1] != "NF*" & m[i,1] != "NF*?" & m[i,1] != "NF?*") {
st<-paste0("ST",m[i,1],collapse=NULL)
}
else {st <- as.character(m[i,1])}
profile<-paste(m[i,-1],collapse="_")
st_labels[i] <- paste(st,profile)
}
return(list(profiles=st_profiles[order(st_profiles[,1]),-1],st=st_list[order(st_profiles[,1])],sts=st_labels))
}


# generate binary matrix of gene presence/absence data
binaryMatrix<-function(m){
bm<-as.matrix(m)
b<-replace(bm, bm =="-", 0.0)
b<-replace(b, bm!="-", 1.0)
b<-replace(b, bm=="?", 0.0)
storage.mode(b)<-'integer'
rownames(b)<-rownames(m)
b
}

# calculate the relative frequency (percentage present) for each gene in each ST
# input matrix is the binary gene matrix, samples in rows
# also requires vector indicating the ST of each sample, ordered the same as rows of the matrix
geneRateMatrix<-function(m,st){
n<-length(unique(st))
gene_rates_ST<-apply(m, 2, function(x) aggregate(x~st,FUN=mean))
gene_rates_ST_matrix<-matrix(unlist(gene_rates_ST),nrow=n*2) [(n+1):(n*2),]
colnames(gene_rates_ST_matrix)<-colnames(m)
rownames(gene_rates_ST_matrix)<-gene_rates_ST[[1]][,1]
storage.mode(gene_rates_ST_matrix)<-'double'
gene_rates_ST_matrix
}

# count number of times each gene observed in each ST
# input matrix is the binary gene matrix, samples in rows
# also requires vector indicating the ST of each sample, ordered the same as rows of the matrix
geneCountMatrix<-function(m,st){
apply(m, 2, function(x) xtabs(x~st))
}


# generate ST vs gene plot (rate or count)
geneSTplot<-function(m,mlst_columns=c(2:9),gene_columns=NULL,plot_type="rate",suppressSNPs=F, suppressUncertainty=F,cluster=F, dataWidth=20, matrix.colours=colorRampPalette(c("white","yellow","blue"),space="rgb")(100)) {
if (is.null(gene_columns)) {gene_columns=13:ncol(m)} # guess
st_analysis<-stProfiles(m[,mlst_columns],suppressSNPs=suppressSNPs, suppressUncertainty=suppressUncertainty)
h<-clusterByST(st_analysis$profiles) # names are pasted profiles
b<-binaryMatrix(m[,gene_columns])
if (ncol(b)==1){
colnames(b)<-colnames(m)[gene_columns]
}
# get gene matrix
if (plot_type=="rate"){
gm<-geneRateMatrix(b,st_analysis$sts)
}
else if (plot_type=="count") {
gm<-geneCountMatrix(b,st_analysis$sts)
}
# plot
plotTreeData(as.phylo(h),matrixFile=gm,infoFile=cbind(st_analysis$st, st_analysis$profiles), freq=table(st_analysis$sts), tip.labels=F, infoCex=1, dataWidth=dataWidth, cluster=cluster, labelHeight=20, colLabelCex=1, matrix.colours=matrix.colours, freqCol=1, freqWidth=5)
}

# generate sample vs gene content plot
geneContentPlot<-function(m,mlst_columns,gene_columns,strain_names=0,suppressSNPs=F, suppressUncertainty=F,cluster=F, dataWidth=20, infoWidth=10, treeWidth=10, matrix.colours=colorRampPalette(c("white","yellow","blue"),space="rgb")(100), labelHeight=20, infoCex = 1) {
h<-clusterByST(m[,mlst_columns])
b<-binaryMatrix(m[,gene_columns])
plotTreeData(as.phylo(h),matrixFile=b,infoFile=m[,c(strain_names,mlst_columns)], tip.labels=F, infoCex=infoCex, dataWidth=dataWidth, infoWidth=infoWidth, treeWidth=treeWidth, cluster=cluster, labelHeight=labelHeight, colLabelCex=1, matrix.colours=matrix.colours, freqCol=1, freqWidth=1)
}

# generate ST barplot, split by some gene
geneSTbarPlot<-function(m,mlst_columns,gene_column,suppressSNPs=F, suppressUncertainty=F,horiz=T,las=NULL,col=1:2){
sts<-m[,mlst_columns[1]]
if (suppressSNPs) {
sts<-gsub("[*]", "", sts)
}
if (suppressUncertainty) {
sts<-gsub("[?]", "", sts)
}
b<-binaryMatrix(m[,gene_column])
colnames(b)<-colnames(m)[gene_column]
if (is.null(las)) {
if (horiz) {las<-1}
else {las<-2}
}
barplot(table(b,sts),col=col,horiz=horiz,las=las)
}

plotTreeData<-function(treeFile,matrixFile=NULL,infoFile=NULL,locFile=NULL,outputPDF=NULL,outputPNG=NULL,w,h,matrix.colours=rev(gray(seq(0,1,0.1))),matrix.legend=F,tip.labels=F,tipLabelSize=1,offset=0,tip.colour.cex=0.5,legend=T,legend.pos="bottomleft",ancestral.reconstruction=F,boundaries=c(0.5,0.75),cluster=F,locColours=NULL,lwd=1.5,axis=T,axisPos=3,edge.color="black",infoCex=0.5,colLabelCex=0.5,treeWidth=10,infoWidth=10,dataWidth=30,edgeWidth=1,labelHeight=10,mainHeight=100,freq=NULL,freqWidth=10,freqCol=2) {
require(ape)
# ladderize tree and extract tip order
if (is.character(treeFile)){
t<-read.tree(treeFile)
}
else t<-treeFile
tl<-ladderize(t)
tips<-tl$edge[,2]
tip.order<-tips[tips<=length(tl$tip.label)]
tip.label.order<-tl$tip.label[tip.order]

# prepare heatmap matrix
if (!is.null(matrixFile)) {
if (is.matrix(matrixFile)) {
x = data.frame(matrixFile,check.names=F)
}
else if (is.data.frame(matrixFile,check.names=F)) {
x = matrixFile
}
else {
x<-read.csv(matrixFile,row.names=1,check.names=F)
}
y.ordered<-x[tip.label.order,]
if (cluster) {
h<-hclust(dist(t(na.omit(y.ordered))),"ward")
y.ordered<-y.ordered[,h$order]
}
}

# prepare frequency plot
if (!is.null(freq)) {
freq<-freq[tip.label.order]
}

# prepare coloured labels for tree leaves
if (!is.null(locFile)) { 
loc<-read.csv(locFile,row.names=1)
loc1<-as.matrix(loc)[row.names(loc) %in% tl$tip.label,] #vector
tipLabelSet <- character(length(loc1))
names(tipLabelSet) <- names(loc1)
groups<-table(loc1)
n<-length(groups)
groupNames<-names(groups)
if (is.null(locColours)){
colours<-rainbow(n)
}
else{
colours<-locColours
}
for (i in 1:n) {
g<-groupNames[i]
tipLabelSet[loc1==g]<-colours[i]
}
# ancestral reconstruction
if (ancestral.reconstruction) {
ancestral<-ace(loc1,tl,type="discrete")
}
else{
ancestral=NULL
}
}
else{
ancestral=NULL
}

# order additional info
if (!is.null(infoFile)) {
if (is.character(infoFile)) { ids<-read.csv(infoFile,row.names=1) }
else {ids <- infoFile}
ids.ordered<-ids[rev(tip.label.order),]
}
else {ids.ordered=NULL}

# open PDF for drawing
if (!is.null(outputPDF)) {
pdf(width=w,height=h,file=outputPDF)
}
# open PNG for drawing
if (!is.null(outputPNG)) {
png(width=w,height=h,file=outputPNG)
}

# set up layout
if (!is.null(infoFile) & !is.null(matrixFile)) {
# tree, info, data heatmap, barplot
layout(matrix(c(0,0,0,4,6,0,0,1,2,3,5,0,0,0,0,0,0,0), nrow=3, byrow=TRUE), width=c(edgeWidth,treeWidth,infoWidth,dataWidth,freqWidth,edgeWidth),height=c(labelHeight,mainHeight,edgeWidth))
}
else if (!is.null(matrixFile)) {
# tree, data heatmap, barplot
layout(matrix(c(0,0,3,5,0,0,1,2,4,0,0,0,0,0,0), nrow=3, byrow=TRUE), width=c(edgeWidth,treeWidth,dataWidth,freqWidth,edgeWidth),height=c(labelHeight,mainHeight,edgeWidth))
}
else {
# only have tree, info
layout(matrix(c(0,0,0,0,0,1,2,0,0,0,0,0), nrow=3, byrow=TRUE), width=c(edgeWidth,treeWidth,infoWidth,edgeWidth),height=c(labelHeight,mainHeight,edgeWidth))
}

# plot tree
par(mar=rep(0,4))
tlp<-plot.phylo(tl,no.margin=T,show.tip.label=tip.labels,label.offset=offset,edge.width=lwd,edge.color=edge.color,xaxs="i", yaxs="i", y.lim=c(0.5,length(tl$tip)+0.5),cex=tipLabelSize)
if (!is.null(locFile)) { 
tiplabels(col= tipLabelSet[tl$tip.label],pch=16,cex=tip.colour.cex) 
if (ancestral.reconstruction) {
nodelabels(pie=ancestral$lik.anc, cex=0.5, piecol=colours)
}
if (axis) {
axisPhylo(axisPos)
}
}
if (matrix.legend && ncol(y.ordered)<20) { text(labels=colnames(y.ordered),x=rep(tlp$x.lim[2]/2,ncol(y.ordered)),y=c(ncol(y.ordered):1)*10) }
if (legend && !is.null(locFile)) {
legend(legend.pos,legend=groupNames,fill=colours)
}

# plot info
if (!is.null(infoFile)) {
par(mar=rep(0,4))
plot(NA,axes=F,pch="",xlim=c(0,ncol(ids.ordered)+1.5),ylim=c(0.5,length(tl$tip)+0.5),xaxs="i",yaxs="i")
# print ST itself with more room
text(x=rep(1,nrow(ids.ordered)),y=c((nrow(ids.ordered)):1),ids.ordered[,1],cex=infoCex)
for (i in 2:ncol(ids.ordered)) {
# print allele profiles
text(x=rep(i+1,nrow(ids.ordered)+1),y=c((nrow(ids.ordered)):1),ids.ordered[,i],cex=infoCex)
# text(rep(i,nrow(ids.ordered)),c((nrow(ids.ordered)):1),ids.ordered[,i],cex=infoCex)
}
}
#else{plot(NA,ylim=c(0,1),xlim=c(0,1),axes=F)}

# plot heatmap
if (!is.null(matrixFile)) {

par(mar=rep(0,4), xpd=TRUE)
image((1:ncol(y.ordered))-0.5, (1:nrow(y.ordered))-0.5, as.matrix(t(y.ordered)),col=matrix.colours,axes=F,xaxs="i", yaxs="i", xlab="",ylab="")

# data labels for heatmap
par(mar=rep(0,4))
plot(NA, axes=F, xaxs="i", yaxs="i", ylim=c(0,2), xlim=c(0.5,ncol(y.ordered)+0.5))
text(1:ncol(y.ordered)-0.5,rep(0,ncol(x)),colnames(y.ordered), srt=90, cex=colLabelCex, pos=4)
}

# frequency barplot
if (!is.null(freq)) {
par(mar=rep(0,4))
barplot(freq, horiz=T, axes=F, xaxs="i", yaxs="i", xlab="", ylab="", ylim=c(0.25,length(freq)+0.25),xlim=c(-1,max(freq,na.rm=T)),col=freqCol,border=0,width=0.5,space=1,names.arg=NA)

# scale for freq plot
par(mar=c(2,0,0,0))
plot(NA, yaxt="n", xaxs="i", yaxs="i", xlab="", ylab="", ylim=c(0,2), xlim=c(-1,max(freq,na.rm=T)),frame.plot=F)
}

# if no freq plot
else{
plot(NA,axes=F)
plot(NA,axes=F)
}

# close drawing device
if (!is.null(outputPDF) | !is.null(outputPNG)) {
dev.off()
}

# return ordered info and ancestral reconstruction object
return(list(id=ids.ordered,anc=ancestral,mat=as.matrix(t(y.ordered))))
}
