#############
# Libraries #
#############

library(ggplot2)
library(ggpubr)
library(data.table)
library(scales)
library(Seurat)
library(plotly)
library(htmlwidgets)
library(SoupX)
library(tidyverse)
library(DropletUtils)
library(vcfR)
library(remotes)
library(DoubletFinder)
library(RColorBrewer)
library(edgeR)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(gplots)
library(scater)
library(monocle)
library(Matrix)
library(dplyr)
library(SingleCellExperiment)
library(ggrepel)
library(colorRamps)
library(enrichplot)
library(ggnewscale)
library(gridExtra)
library(rstatix)

#############
# Functions #
#############

GeneNames<-function(data,symbol=":",what="Names",remove.dot=TRUE){
    ss<-strsplit(as.character(data), symbol)
    k<-ifelse(what=="Names",2,1)
    res<-rep(0,length(ss))
    for(i in 1:length(ss)){
        res[i]<-ss[[i]][min(k,length(ss[[i]]))]
        if(k==1 & remove.dot==TRUE){
            res[i]<-unlist(strsplit(as.character(res[i]),".",fixed=TRUE))[1]
        }
    }
    return(res)
}


tukey_fun<-function(data,quantile,fit,pcut,disease_state=NULL,cell_type=NULL){
    
    if(fit=="~Cell_Type"){
        if(is.null(disease_state)){
            stop("Please provide a disease state")
        }
        x<-data[data$Disease_State==disease_state,]
        mod<-aov(x[,which(colnames(x)==quantile)]~x$Cell_Type,data=x)
        tuk<-TukeyHSD(mod, conf.level=.95)
        res<-tuk$`x$Cell_Type`[tuk$`x$Cell_Type`[,4]<=pcut,]
    }

    if(fit=="~Disease_State"){
        if(is.null(cell_type)){
            stop("Please provide a cell type")
        }
        x<-data[data$Cell_Type==cell_type,]
        mod<-aov(x[,which(colnames(x)==quantile)]~x$Disease_State,data=x)
        tuk<-TukeyHSD(mod, conf.level=.95)
        res<-tuk$`x$Disease_State`[tuk$`x$Disease_State`[,4]<=pcut,]
    }
 return(res)
}


pseudobulk_conversion<-function(data){
    
    data.ND<-data[,data$Condition=="ND"]
    data.PRE<-data[,data$Condition=="PD"]
    data.T2D<-data[,data$Condition=="T2D"]
    donors.ND<-as.character(data.ND$Islet)
    donors.PRE<-as.character(data.PRE$Islet)
    donors.T2D<-as.character(data.T2D$Islet)

    nn<-table(donors.ND)
    uu<-names(nn)
    NDcounts<-matrix(0,nrow(data.ND),length(uu))
    NDdesign<-matrix(0,length(uu),8)
    for(i in 1:length(uu)){
           w<-which(donors.ND==uu[i])
           if(length(w)>0){
               mat<-matrix(data.ND@assays$RNA@counts[,w],nrow=nrow(data.ND@assays$RNA@counts))
               NDcounts[,i]<-rowSums(mat)
               des<-as.matrix(data.ND@meta.data[w,])
               NDdesign[i,]<-c(des[1,8:14],nn[i])
           }
    }
    NDcounts<-data.frame(NDcounts)
    NDdesign<-data.frame(NDdesign)
    colnames(NDcounts)<-uu
    rownames(NDcounts)<-rownames(data.ND@assays$RNA@counts)
    rownames(NDdesign)<-uu
    colnames(NDdesign)<-c("Chemistry","Condition","Sex","BMI","Race","HbA1c","Age","Ncells")
    NDcounts<-round(NDcounts,0)
    keep<-apply(NDcounts,2,sum)>0
    NDcounts<-NDcounts[,keep]
    NDdesign<-NDdesign[keep,]
    
    
    nn<-table(donors.PRE)
    uu<-names(nn)
    PREcounts<-matrix(0,nrow(data.PRE),length(uu))
    PREdesign<-matrix(0,length(uu),8)
    for(i in 1:length(uu)){
            w<-which(donors.PRE==uu[i])
            if(length(w)>0){
                mat<-matrix(data.PRE@assays$RNA@counts[,w],nrow=nrow(data.PRE@assays$RNA@counts))
                PREcounts[,i]<-rowSums(mat)
                des<-as.matrix(data.PRE@meta.data[w,])
                PREdesign[i,]<-c(des[1,8:14],nn[i])
            }
    }
    PREcounts<-data.frame(PREcounts)
    PREdesign<-data.frame(PREdesign)
    colnames(PREcounts)<-uu
    rownames(PREcounts)<-rownames(data.PRE@assays$RNA@counts)
    rownames(PREdesign)<-uu
    colnames(PREdesign)<-c("Chemistry","Condition","Sex","BMI","Race","HbA1c","Age","Ncells")
    PREcounts<-round(PREcounts,0)
    keep<-apply(PREcounts,2,sum)>0
    PREcounts<-PREcounts[,keep]
    PREdesign<-PREdesign[keep,]

    
    nn<-table(donors.T2D)
    uu<-names(nn)
    T2Dcounts<-matrix(0,nrow(data.T2D),length(uu))
    T2Ddesign<-matrix(0,length(uu),8)
    for(i in 1:length(uu)){
            w<-which(donors.T2D==uu[i])
            if(length(w)>0){
                mat<-matrix(data.T2D@assays$RNA@counts[,w],nrow=nrow(data.T2D@assays$RNA@counts))
                T2Dcounts[,i]<-rowSums(mat)
                des<-as.matrix(data.T2D@meta.data[w,])
                T2Ddesign[i,]<-c(des[1,8:14],nn[i])
            }
    }
    T2Dcounts<-data.frame(T2Dcounts)
    T2Ddesign<-data.frame(T2Ddesign)
    colnames(T2Dcounts)<-uu
    rownames(T2Dcounts)<-rownames(data.T2D@assays$RNA@counts)
    rownames(T2Ddesign)<-uu
    colnames(T2Ddesign)<-c("Chemistry","Condition","Sex","BMI","Race","HbA1c","Age","Ncells")
    T2Dcounts<-round(T2Dcounts,0)
    keep<-apply(T2Dcounts,2,sum)>0
    T2Dcounts<-T2Dcounts[,keep]
    T2Ddesign<-T2Ddesign[keep,]

    
    counts<-cbind(NDcounts,PREcounts,T2Dcounts)
    design<-rbind(NDdesign,PREdesign,T2Ddesign)

 return(list(Counts=counts,Design=design))
}


MDS_fun<-function(data,top=1000){
    
    y<-data$Counts
    desi<-data$Design
    
    group<-factor(desi$Condition)
    y<-DGEList(counts=y,group=group)
    rND<-rowSums(y$counts[,y$samples$group=="ND"]>0)
    rPRE<-rowSums(y$counts[,y$samples$group=="PD"]>0)
    rT2D<-rowSums(y$counts[,y$samples$group=="T2D"]>0)
    keepND<-which(rND>0.3*ncol(y$counts[,y$samples$group=="ND"]))
    keepPRE<-which(rPRE>0.3*ncol(y$counts[,y$samples$group=="PD"]))
    keepT2D<-which(rT2D>0.3*ncol(y$counts[,y$samples$group=="T2D"]))
    keep<-unique(c(keepND,keepPRE,keepT2D))
    yy<-y[keep,,keep.lib.sizes=FALSE]
    yy <- calcNormFactors(yy)

    f1<-factor(desi$Condition)
    f2<-factor(desi$Chemistry)
    f3<-factor(desi$Sex)
    f4<-factor(desi$Race)
    
    if(length(which(colnames(desi)=="CellType"))>0){
        f5<-factor(desi$CellType)
    }
    
    mds<-plotMDS.DGEList(yy,top=top,plot = FALSE)
    if(length(which(names(mds)=="cmdscale.out"))>0){
        x<-data.frame(mds$cmdscale.out)
        colnames(x)<-c("MDS_1","MDS_2")
    } else {
        x<-data.frame(MDS_1=mds$x,MDS_2=mds$y)
    }
    x$Condition<-f1
    x$Sex<-f3
    x$Chemistry<-f2
    x$Race<-f4
    if(length(which(colnames(desi)=="CellType"))>0){
        x$CellType<-f5
    }

    p1<-ggplot(x,aes(x=MDS_1,y=MDS_2,color=Condition))+geom_point()
    p2<-ggplot(x,aes(x=MDS_1,y=MDS_2,color=Sex))+geom_point()
    p3<-ggplot(x,aes(x=MDS_1,y=MDS_2,color=Chemistry))+geom_point()
    p4<-ggplot(x,aes(x=MDS_1,y=MDS_2,color=Race))+geom_point()
    p<-ggarrange(p1,p2,p3,p4,ncol=2,nrow=2)
    
    if(length(which(colnames(desi)=="CellType"))>0){
        pp<-ggplot(x,aes(x=MDS_1,y=MDS_2,color=CellType))+geom_point()
        p<-list(p,pp)
    }
    
 return(list(x,p))
}


edgeR_fun<-function(data,cond.1,cond.2,include.chemistry=TRUE){
    
    y<-data$Counts
    desi<-data$Design
    group<-factor(desi$Condition)
    y<-DGEList(counts=y,group=group)
   

    if(include.chemistry){
        ww<-which(desi$Condition==cond.1 | desi$Condition==cond.2)
    } else {
        ww<-which(desi$Condition==cond.1 & desi$Chemistry=="V3" | desi$Condition==cond.2 & desi$Chemistry=="V3")
    }
    y1<-y[,ww]
    desi1<-desi[ww,]
    r1<-rowSums(y1$counts[,y1$samples$group==cond.1]>0)
    r2<-rowSums(y1$counts[,y1$samples$group==cond.2]>0)
    keep1<-which(r1>0.3*ncol(y1$counts[,y1$samples$group==cond.1]))
    keep2<-which(r2>0.3*ncol(y1$counts[,y1$samples$group==cond.2]))
    keep<-unique(c(keep1,keep2))
    y1<-y1[keep,,keep.lib.sizes=FALSE]
    f1<-factor(desi1$Condition)
    f2<-factor(desi1$Chemistry)
    f3<-factor(desi1$Sex)
    f4<-factor(desi1$Race)
    f5<-scale(as.numeric(as.character(desi1$BMI)),scale=T)
    f6<-as.numeric(as.character(desi1$HbA1c))
    f7<-scale(as.numeric(as.character(desi1$Age)),scale=T)

    conds<-sort(c(cond.1,cond.2))
    if(include.chemistry){
        designC<-model.matrix(~0+f1+f2+f3+f4+f5+f7)
        colnames(designC)<-c(conds[1],conds[2],"V3","Male","Hispanic","White","BMI","Age")
    } else {
        designC<-model.matrix(~0+f1+f3+f4+f5+f7)
        colnames(designC)<-c(conds[1],conds[2],"Male","Hispanic","White","BMI","Age")
    }
    y1 <- calcNormFactors(y1)
    y1<-estimateDisp(y1, designC)
    fit <- glmFit(y1, designC)
    fit1 <- glmLRT(fit, contrast=c(-1,1,rep(0,(ncol(designC)-2))))
    res1<-data.frame(topTags(fit1,n=nrow(fit1$counts)),Comp=paste(conds[2],"-",conds[1]))
    res1.edger<-cbind(Gene=rownames(res1),res1)

 return(res1.edger)
}


Pseudo_marker_celltype<-function(data,condition,celltype){
    y <- data$Counts
    desi <- data$Design
    w<-which(match(desi$Condition,condition,nomatch=0)>0)
    y<-y[,w]
    desi<-desi[w,]
    desi$cluster_id<-"ZZZ"
    desi$cluster_id[desi$CellType==celltype]<-celltype
    group<-factor(desi$cluster_id)
    y<-DGEList(counts=y,group=group)
    
    rC1 <-rowSums(y$counts[,y$samples$group== paste0(celltype)]>0)
    rOther<-rowSums(y$counts[,y$samples$group=="ZZZ"]>0)
    keepC1<-which(rC1>0.3*ncol(y$counts[,y$samples$group==paste0(celltype)]))
    keepOther<-which(rOther>0.3*ncol(y$counts[,y$samples$group=="ZZZ"]))
    keep<-unique(c(keepC1,keepOther))
    y<-y[keep,,keep.lib.sizes=FALSE]
    f1<-factor(desi$cluster_id)
    f2<-factor(desi$Sex)
    f3<-factor(desi$Race)
    f4<-factor(desi$Chemistry)
    f5<-factor(desi$Condition)
    f6<-as.numeric(as.character(desi$HbA1c))
    f7<-as.numeric(as.character(desi$Age))
    f8 <- as.numeric(as.character(desi$n_cells))
    designC<-model.matrix(~0+f1+f2+f3+f4+f7)
    colnames(designC)<-c("Acinar", "Other", "Male", "Hispanic", "White", "V3", "Age")

    y <- calcNormFactors(y)
    y <-estimateDisp(y, designC)
    fit <- glmFit(y, designC)
    fit <- glmLRT(fit, contrast=c(1,-1,0,0,0,0,0))
    res<-data.frame(topTags(fit,n=nrow(fit$counts)),Comp=(paste0(celltype,"-Other")))
    res<-cbind(Gene=rownames(res),res)

 return(res)
}



pairwise_analysis<-function(data,celltypes,condition,include.chemistry=TRUE){
    
    for(i in 1:(length(celltypes)-1)){
      
      w<-which(data$Design$CellType==celltypes[i] &
               data$Design$Condition==condition)
      my.counts.1<-data$Counts[,w]
      my.design.1<-data$Design[w,]

      for(j in (i+1):length(celltypes)){
        w<-which(data$Design$CellType==celltypes[j] &
                 data$Design$Condition==condition)
        my.counts.2<-data$Counts[,w]
        my.design.2<-data$Design[w,]

        y<-cbind(my.counts.1,my.counts.2)
        desi<-rbind(my.design.1,my.design.2)
        
        group<-factor(desi$CellType)
        y<-DGEList(counts=y,group=group)
       
        if(include.chemistry){
            ww<-1:nrow(desi)
        } else {
            ww<-which(desi$Chemistry=="V3")
        }
        y1<-y[,ww]
        desi1<-desi[ww,]
        r1<-rowSums(y1$counts[,y1$samples$group==celltypes[i]]>0)
        r2<-rowSums(y1$counts[,y1$samples$group==celltypes[j]]>0)
        keep1<-which(r1>0.3*ncol(y1$counts[,y1$samples$group==celltypes[i]]))
        keep2<-which(r2>0.3*ncol(y1$counts[,y1$samples$group==celltypes[j]]))
        keep<-unique(c(keep1,keep2))
        y1<-y1[keep,,keep.lib.sizes=FALSE]
        f1<-factor(desi1$CellType)
        f2<-factor(desi1$Chemistry)
        f3<-factor(desi1$Sex)
        f4<-factor(desi1$Race)
        f5<-scale(as.numeric(as.character(desi1$BMI)),scale=T)
        f6<-as.numeric(as.character(desi1$HbA1c))
        f7<-scale(as.numeric(as.character(desi1$Age)),scale=T)

        if(include.chemistry){
            designC<-model.matrix(~0+f1+f2+f3+f4+f5+f7)
            colnames(designC)<-c(celltypes[i],celltypes[j],"V3","Male","Hispanic","White","BMI","Age")
        } else {
            designC<-model.matrix(~0+f1+f3+f4+f5+f7)
            colnames(designC)<-c(celltypes[i],celltypes[j],"Male","Hispanic","White","BMI","Age")
        }
        y1 <- calcNormFactors(y1)
        y1<-estimateDisp(y1, designC)
        fit <- glmFit(y1, designC)
        fit1 <- glmLRT(fit, contrast=c(1,-1,rep(0,(ncol(designC)-2))))
        res1<-data.frame(topTags(fit1,n=nrow(fit1$counts)),Comp=paste(celltypes[i],"-",celltypes[j]))
        res1.edger<-cbind(Gene=rownames(res1),res1)

        if(i==1 & j==2){
          pairwise_results<-res1.edger
        } else {
          pairwise_results<-rbind(pairwise_results,res1.edger)
        }
        
      }
    }

 return(pairwise_results)
}


pairwise_unique<-function(data,celltypes,logfc.cut,fdr.cut){
    
    genes<-as.list(rep(0,length(celltypes)))
    names(genes)<-celltypes
    for(i in 1:length(celltypes)){
        g<-grep(celltypes[i],data$Comp)
        dd<-data[g,]
        dd<-dd[abs(as.numeric(as.character(dd$logFC)))>=logfc.cut & dd$FDR<=fdr.cut,]
        uu<-names(table(dd$Comp))
        genes1<-NULL
        for(j in 1:length(uu)){
            g<-grep(paste(celltypes[i]," -",sep=""),uu[j])
            if(j==1){
                if(length(g)>0){
                    ww<-which(dd$logFC>0 & dd$Comp==uu[j])
                } else {
                    ww<-which(dd$logFC<0 & dd$Comp==uu[j])
                }
                genes1<-dd[ww,c(1,ncol(dd))]
            } else {
                if(length(g)>0){
                    ww<-which(dd$logFC>0 & dd$Comp==uu[j])
                } else {
                    ww<-which(dd$logFC<0 & dd$Comp==uu[j])
                }
                genes1<-rbind(genes1,dd[ww,c(1,ncol(dd))])
            }
        }
        genes1<-table(genes1[,1])
        genes[[i]]<-names(genes1)[genes1==3]
    }
    v<-venn(genes,show.plot=T)
    v<-attributes(v)$intersections
    v<-cbind(unlist(v),rep(names(v),lapply(v,length)))
    
 return(v)
}



pseudoQC<-function(data){
    
    cou<-data$Counts
    desi<-data$Design

    oo<-order(desi$Condition,desi$Chemistry)
    cou<-cou[,oo]
    desi<-desi[oo,]

    # 1
    lib<-apply(cou,2,sum)
    adjlib<-lib/desi$Ncells
    
    # 2.1
    mito<-grep("MT-",rownames(cou))
    ins<-which(rownames(cou)=="INS")
    gcg<-which(rownames(cou)=="GCG")
    ppy<-which(rownames(cou)=="PPY")
    sst<-which(rownames(cou)=="SST")
    out<-c(mito,ins,gcg,ppy,sst)
    lib.nomito<-apply(cou[-out,],2,sum)
    adjlib.nomito<-lib.nomito/desi$Ncells

    # 2.2
    out<-mito
    lib.nomito1<-apply(cou[-out,],2,sum)
    adjlib.nomito1<-lib.nomito1/desi$Ncells
    
    # 2.3
    out<-c(ins,gcg,ppy,sst)
    lib.nomito2<-apply(cou[-out,],2,sum)
    adjlib.nomito2<-lib.nomito2/desi$Ncells

    # 3
    mitoperc<-cou[mito,]
    mitoperc<-100 * apply(mitoperc, 2, sum) / lib

    # 4
    detg <- apply(cou, 2, function(x) length(which(x > 0)))

    # 5
    sce <- SingleCellExperiment(list(counts=as.matrix(cou)))
    pp<-plotHighestExprs(sce, exprs_values = "counts",as_percentage = FALSE)
    pp<-as.character(pp$data$Tag[1:50])
    mm<-match(pp,rownames(cou),nomatch=0)
    tops<-t(cou[mm,])

    dat<-data.frame(ID=rownames(desi),
                    Chemistry=desi$Chemistry,
                    Condition=desi$Condition,
                    LibSize=lib/1e+6,
                    adj.LibSize=adjlib,
                    LibSize.nomito<-lib.nomito/1e+6,
                    adj.LibSize.nomito=adjlib.nomito,
                    LibSize.nomito1<-lib.nomito1/1e+6,
                    adj.LibSize.nomito1=adjlib.nomito1,
                    LibSize.nomito2<-lib.nomito2/1e+6,
                    adj.LibSize.nomito2=adjlib.nomito2,
                    MTperc=mitoperc,
                    DetGenes=detg)
    dat2<-cbind(tops,dat[,c(1:3)])
    dat2<-reshape2::melt(dat2)
    dat2$value<-log(as.numeric(as.character(dat2$value))+1,2)
    
    dat$ID<-factor(dat$ID,levels=unique(dat$ID))
    dat2$ID<-factor(dat2$ID,levels=unique(dat2$ID))
    p1<-ggplot(dat,aes(x=ID,y=LibSize,fill=Condition))+geom_bar(stat="identity") + xlab("") + ylab("Library Size (x 1M)") + ylim(0,100) +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    p2<-ggplot(dat,aes(x=ID,y=LibSize,fill=Chemistry))+geom_bar(stat="identity") + xlab("") + ylab("Library Size (x 1M)") + ylim(0,100) +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    p3<-ggplot(dat,aes(x=ID,y=adj.LibSize,fill=Condition))+geom_bar(stat="identity") + xlab("") + ylab("Adj. Library Size (by number of cells)") + ylim(0,40000) +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    p4<-ggplot(dat,aes(x=ID,y=adj.LibSize,fill=Chemistry))+geom_bar(stat="identity") + xlab("") + ylab("Adj. Library Size (by number of cells)") + ylim(0,40000) +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    p5<-ggplot(dat,aes(x=ID,y=LibSize.nomito,fill=Condition))+geom_bar(stat="identity") + xlab("") + ylab("Library Size w/o MT and Endo markers (x 1M)") + ylim(0,100) +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    p6<-ggplot(dat,aes(x=ID,y=LibSize.nomito,fill=Chemistry))+geom_bar(stat="identity") + xlab("") + ylab("Library Size w/o MT and Endo markers (x 1M)") + ylim(0,100) +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    p7<-ggplot(dat,aes(x=ID,y=adj.LibSize.nomito,fill=Condition))+geom_bar(stat="identity") + xlab("") + ylab("Adj. Library Size w/o MT and Endo markers (by number of cells)") + ylim(0,40000) +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    p8<-ggplot(dat,aes(x=ID,y=adj.LibSize.nomito,fill=Chemistry))+geom_bar(stat="identity") + xlab("") + ylab("Adj. Library Size w/o MT and Endo markers (by number of cells)") + ylim(0,40000) +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    p9<-ggplot(dat,aes(x=ID,y=MTperc,fill=Condition))+geom_bar(stat="identity") + xlab("") + ylab("% MT") + ylim(0,25) +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
        
    p5.1<-ggplot(dat,aes(x=ID,y=LibSize.nomito1,fill=Condition))+geom_bar(stat="identity") + xlab("") + ylab("Library Size w/o MT (x 1M)") + ylim(0,100) +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    p6.1<-ggplot(dat,aes(x=ID,y=LibSize.nomito1,fill=Chemistry))+geom_bar(stat="identity") + xlab("") + ylab("Library Size w/o MT (x 1M)") + ylim(0,100) +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    p7.1<-ggplot(dat,aes(x=ID,y=adj.LibSize.nomito1,fill=Condition))+geom_bar(stat="identity") + xlab("") + ylab("Adj. Library Size w/o MT (by number of cells)") + ylim(0,40000) +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    p8.1<-ggplot(dat,aes(x=ID,y=adj.LibSize.nomito1,fill=Chemistry))+geom_bar(stat="identity") + xlab("") + ylab("Adj. Library Size w/o MT (by number of cells)") + ylim(0,40000) +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
        
    p5.2<-ggplot(dat,aes(x=ID,y=LibSize.nomito2,fill=Condition))+geom_bar(stat="identity") + xlab("") + ylab("Library Size w/o Endo markers (x 1M)") + ylim(0,100) +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    p6.2<-ggplot(dat,aes(x=ID,y=LibSize.nomito2,fill=Chemistry))+geom_bar(stat="identity") + xlab("") + ylab("Library Size w/o Endo markers (x 1M)") + ylim(0,100) +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    p7.2<-ggplot(dat,aes(x=ID,y=adj.LibSize.nomito2,fill=Condition))+geom_bar(stat="identity") + xlab("") + ylab("Adj. Library Size w/o Endo markers (by number of cells)") + ylim(0,40000) +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    p8.2<-ggplot(dat,aes(x=ID,y=adj.LibSize.nomito2,fill=Chemistry))+geom_bar(stat="identity") + xlab("") + ylab("Adj. Library Size w/o Endo markers (by number of cells)") + ylim(0,40000) +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
        
    p10<-ggplot(dat,aes(x=ID,y=MTperc,fill=Chemistry))+geom_bar(stat="identity") + xlab("") + ylab("% MT") + ylim(0,25) +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    p11<-ggplot(dat,aes(x=ID,y=DetGenes,fill=Condition))+geom_bar(stat="identity") + xlab("") + ylab("Detected Genes") + ylim(0,25000) +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    p12<-ggplot(dat,aes(x=ID,y=DetGenes,fill=Chemistry))+geom_bar(stat="identity") + xlab("") + ylab("Detected Genes") + ylim(0,25000) +
        theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    p13<-ggplot(dat2,aes(x=variable,y=value,fill=Condition))+geom_boxplot() + xlab("") + ylab("log2 Expression") + ylim(0,25) +
        theme(text = element_text(size=10),axis.text.x = element_text(angle=90, hjust=1))
    p14<-ggplot(dat2,aes(x=variable,y=value,fill=Chemistry))+geom_boxplot() + xlab("") + ylab("log2 Expression") + ylim(0,25) +
        theme(text = element_text(size=10),axis.text.x = element_text(angle=90, hjust=1))

    pp<-list(p1,p2,p3,p4,p5,p6,p7,p8,
                         p5.1,p6.1,p7.1,p8.1,
                         p5.2,p6.2,p7.2,p8.2,
             p9,p10,p11,p12,p13,p14)
        
 return(list(dat,dat2,pp))
}



#################################
# Functions for HbA1c analysis  #
#################################

#' Generates a Lorgnette object from the raw counts and metadata matrices.
#'
#' Generates a Lorgnette object from the raw counts and metadata matrices. It is essentially based on
#'   the construction of a Monocle2 object. It also requires a gene description table with
#'   the Ensembl IDs and gene symbols.
#' @param counts data frame. The raw counts with the genes in rows and samples in columns.
#' @param metadata data frame. The metadata table of the samples.
#' @param gene.attr data frame. A table with the gene information. it should always contain
#'   the columns 'Ensembl' and 'Symbol' with the EnsemblIDs and the gene symbols respectively.
#'   It provides the information to convert the rownames(counts) into the format 'EnsemblID:Symbol'.
#' @param gene.id character. If NULL (default), it assumes that rownames(counts) are already in the
#'   'EnsemblID:Symbol' format, so no conversion is needed. If its value is 'Symbol' or 'EnsemblID',
#'   it assums that rownames(counts) are gene symbols or ensembl IDs respectively and it will convert
#'   them to EnsemblID:Symbol format.
#' @keywords Lorgnette_from_matrix
#' @return A Lorgnette object with the size factors and dispersions
#'
Lorgnette_from_matrix<-function(counts,metadata,gene.attr,gene.id=NULL){
    
    mm<-match(c("Symbol","Ensembl"),colnames(gene.attr),nomatch=0)
    if(length(mm[mm==0])>0){
        stop("gene.attr should always contain the columns Ensembl and Symbol. Either they are missing or the colnames need to be editted.")
    }
    
    if(!is.null(gene.id)){
        mm<-match(gene.id,c("Symbol","Ensembl"),nomatch=0)
        if(mm==0){
            stop("gene.id should take only the values Ensembl, Symbol or NULL!")
        }
        
        if(gene.id=="Symbol"){
            
            print("Assuming that rownames(counts) are gene symbols. Converting them to EnsemblID:Symbol!")
            
            mm<-match(rownames(counts),gene.attr$Symbol,nomatch=0)
            gene.attr1<-gene.attr[mm,]
            if(length(mm[mm==0])>0){
                if(ncol(gene.attr)>2){
                    gene.attr2<-cbind(Ensembl=paste("ENSG",1:length(mm[mm==0]),sep=""),
                                      Symbol=rownames(counts)[mm==0],
                                      matrix("---",length(mm[mm==0]),(ncol(gene.attr)-2)))
                 } else {
                     gene.attr2<-cbind(Ensembl=paste("ENSG",1:length(mm[mm==0]),sep=""),
                                       Symbol=rownames(counts)[mm==0])
                 }
                 gene.attr2<-data.frame(gene.attr2)
                 colnames(gene.attr2)<-colnames(gene.attr)
                 gene.attr1<-rbind(gene.attr1,gene.attr2)
            }
            gene.attr1<-cbind(Gene=apply(gene.attr1[,1:2],1,paste,collapse=":"),gene.attr1)
            counts1<-counts[match(gene.attr1$Symbol,rownames(counts)),]
            rownames(counts1)<-gene.attr1$Gene

        }
        
        if(gene.id=="Ensembl"){
            
            print("Assuming that rownames(counts) are Ensembl IDs. Converting them to EnsemblID:Symbol!")
            
            mm<-match(rownames(counts),gene.attr$Ensembl,nomatch=0)
            gene.attr1<-gene.attr[mm,]
            if(length(mm[mm==0])>0){
                if(ncol(gene.attr)>2){
                    gene.attr2<-cbind(Ensembl=rownames(counts)[mm==0],
                                      Symbol=paste("Gene",1:length(mm[mm==0]),sep=""),
                                      matrix("---",length(mm[mm==0]),(ncol(gene.attr)-2)))
                 } else {
                     gene.attr2<-cbind(Ensembl=rownames(counts)[mm==0],
                                       Symbol=paste("Gene",1:length(mm[mm==0]),sep=""))
                 }
                 gene.attr2<-data.frame(gene.attr2)
                 colnames(gene.attr2)<-colnames(gene.attr)
                 gene.attr1<-rbind(gene.attr1,gene.attr2)
            }
            gene.attr1<-cbind(Gene=apply(gene.attr1[,1:2],1,paste,collapse=":"),gene.attr1)
            mm<-match(gene.attr1$Ensembl,rownames(counts),nomatch=0)
            counts1<-counts[mm,]
            rownames(counts1)<-gene.attr1$Gene
            gene.attr1<-gene.attr1[mm>0,]

        }
        
    }
    
    
    if(is.null(gene.id)){
        
        print("Assuming that rownames(counts) are EnsemblID:Symbol.")
        
        mm<-match(GeneNames(rownames(counts),what="Ensg"),gene.attr$Ensembl,nomatch=0)
        gene.attr1<-gene.attr[mm,]
        if(length(mm[mm==0])>0){
            if(ncol(gene.attr)>2){
                gene.attr2<-cbind(Ensembl=rownames(counts)[mm==0],
                                  Symbol=paste("Gene",1:length(mm[mm==0]),sep=""),
                                  matrix("---",length(mm[mm==0]),(ncol(gene.attr)-2)))
             } else {
                 gene.attr2<-cbind(Ensembl=rownames(counts)[mm==0],
                                   Symbol=paste("Gene",1:length(mm[mm==0]),sep=""))
             }
             gene.attr2<-data.frame(gene.attr2)
             colnames(gene.attr2)<-colnames(gene.attr)
             gene.attr1<-rbind(gene.attr1,gene.attr2)
        }
        gene.attr1<-cbind(Gene=apply(gene.attr1[,1:2],1,paste,collapse=":"),gene.attr1)
        mm<-match(gene.attr1$Ensembl,GeneNames(rownames(counts),what="Ensg"),nomatch=0)
        counts1<-counts[mm,]
        rownames(counts1)<-gene.attr1$Gene
        gene.attr1<-gene.attr1[mm>0,]
        
    }

    gene.attr1$gene_short_name<-gene.attr1$Symbol
    rownames(gene.attr1)<-gene.attr1$Gene
    
    mm<-match(colnames(counts1),rownames(metadata),nomatch=0)
    if(length(mm[mm==0])>0){
        stop("The colnames(counts) and rownames(metadata) do not match!")
    }
    metadata<-metadata[mm,]
    
    
    pd <- new("AnnotatedDataFrame", data = metadata)
    fd <- new("AnnotatedDataFrame", data = gene.attr1)
    obj <- monocle::newCellDataSet(as(counts1, "sparseMatrix"),
                                   phenoData = pd,
                                   featureData = fd,
                                   expressionFamily=VGAM::negbinomial.size())

 return(obj)
}





#' Performs filtering on the metadata table.
#'
#' Performs filtering on the metadata table using a filters expression.
#' @param data data frame. The metadata table.
#' @param filters character. The filters in the format used by the filter() function
#'   of dplyr R package.
#' @keywords filterData
#' @return A filtered metadata table.
#'
filterData<-function(data,filters){
    
    uu<-unlist(strsplit(filters,"[&|]"))
    o<-which(uu=="")
    if(length(o)>0){
        uu<-uu[-o]
    }
    for(i in 1:length(uu)){
        uu1<-unlist(strsplit(uu[i],"[==><!=]"))
        uuStart<-uu1[1]
        uuStart<-unlist(strsplit(uuStart,""))
        o<-which(uuStart=="'" | uuStart==" ")
        if(length(o)>0){
            uuStart<-uuStart[-o]
        }
        uuStart<-paste(uuStart,collapse="")
        
        uuEnd<-uu1[length(uu1)]
        uuEnd<-unlist(strsplit(uuEnd,""))
        o<-which(uuEnd=="'" | uuEnd==" ")
        if(length(o)>0){
            uuEnd<-uuEnd[-o]
        }
        uuEnd<-paste(uuEnd,collapse="")
        
        if(suppressWarnings(is.na(as.numeric(uuEnd)))){
            uuEnd<-unlist(strsplit(uuEnd,""))
            o<-which(uuEnd=="'" | uuEnd==" ")
            if(length(o)>0){
                uuEnd<-uuEnd[-o]
            }
            uuEnd<-paste(uuEnd,collapse="")
        }
        
        w<-which(colnames(data)==uuStart)
        if(length(w)==0){
            stop(paste("Filters variable ",uuStart," does not exist in the pData(obj)!",sep=""))
        }
        ll<-length(which(data[,w]==uuEnd))
        if(ll==0){
            stop(paste("Factor level ",uuEnd," of filters variable ",uuStart," does not exist in the pData(obj)!",sep=""))
        }
        
    }
    
    filters<-as.expression(parse(text=filters))
    res<-tryCatch({
        data %>% filter(eval(filters[1]))},
        error=function(x){
      return(NULL)
    })
  
    if(is.null(res)){
        stop("One or more filter variables does not exist in the pData(obj)!",sep="")
  }

 return(res)
}




#' Constructs the Lorgnette dataset for DE analysis.
#'
#' Constructs the Lorgnette dataset for DE analysis using data filtering, a pseudotime variable and
#'   a state variable. The user has to define which characteristics serves as pseudotime and which
#'   as state. The samples will be annotated accordingly for downstream analysis.
#' @param obj object. A Lorgnette object.
#' @param filters character. The filters in the format used by the filter() function
#'   of dplyr R package. Default is NULL where no filtering is performed
#' @param pseudotime.by character. The name of the variable of pData(obj) to serve as the pseudotime
#'   variable.
#' @param state.by character. The name of the variable of pData(obj) to serve as the state variable.
#' @keywords LorgnetteSet
#' @return An updated Lorgnette object.
#'
LorgnetteSet<-function(obj,pseudotime.by,state.by,filters=NULL){
    
    mm<-match("Size_Factor",colnames(pData(obj)),nomatch=0)
    if(mm==0){
        pData(obj)$Size_Factor<-NA
        print("Estimating size factors and dispersions...")
            obj <- estimateSizeFactors(obj)
            obj <- estimateDispersions(obj)
    }
    
    mm<-match("Size_Factor",colnames(pData(obj)),nomatch=0)
    if(mm>0){
        if(is.na(pData(obj)$Size_Factor[1])){
            print("Estimating size factors and dispersions...")
                obj <- estimateSizeFactors(obj)
                obj <- estimateDispersions(obj)
        }
    }
    
    if(is.null(filters)){
        print("Variable filters was not defined. All data will be used for analysis!")
    } else {
        ff<-filterData(data=pData(obj),filters)
        mm<-match(rownames(ff),rownames(pData(obj)))
        obj<-obj[,mm]
    }
    
    if(is.null(pseudotime.by)){
        w<-which(colnames(pData(obj))=="Pseudotime")
        if(length(w)==0){
            stop("The Pseudotime column does not exist in the pData(obj). Please define it via pseudotime.by!")
        }
    } else {
        w<-which(colnames(pData(obj))==pseudotime.by)
        if(length(w)==0){
            stop("The value of the pseudotime.by parameter does not exist in the pData(obj).")
        }
        colnames(pData(obj))[w]<-"Pseudotime"
    }
    print(paste("Setting ",pseudotime.by," as pseudotime.",sep=""))
    
    
    if(is.null(state.by)){
        w<-which(colnames(pData(obj))=="State")
        if(length(w)==0){
            stop("The State column does not exist in the pData(obj). Please define it via state.by!")
        }
    } else {
        w<-which(colnames(pData(obj))==state.by)
        if(length(w)==0){
            stop("The value of the state.by parameter does not exist in the pData(obj).")
        }
        colnames(pData(obj))[w]<-"State"
    }
    
    
    br<-which(colnames(pData(obj))=="Branch")
    if(length(br)>0){
       colnames(pData(obj))[br]<-"Branch.1"
    } else {
        br<-rep(NA,nrow(pData(obj)))
        ss<-sort(unique(pData(obj)$State))
        for(i in 1:length(ss)){
            br[pData(obj)$State==ss[i]]<- i
        }
        pData(obj)$Branch<-br
    }
    
    ll<-table(pData(obj)$Pseudotime,pData(obj)$Branch)
    ll1<-apply(ll,2,function(x) length(which(x>0)))
    for(i in 1:length(ll1)){
        print(paste(ll1[i]," unique pseudotimes (out of ",sum(ll[,i]),") in branch ",names(ll1)[i],sep=""))
    }
    
    
    obj<-obj[,order(pData(obj)$Branch,pData(obj)$Pseudotime)]
    
 return(obj)
}



#' Selects the data to be analyzed downstream of the LorgnetteDE() analysis.
#'
#' Selects the data to be analyzed downstream od the LorgnetteDE() analysis. It is
#'   useful in case where parameter branch of LorgnetteDE() took a vector of values.
#'   This function wil select the results of interest for downstream analysis.
#' @param obj object. A Lorgnette object. Typically, the output of LorgnetteSet().
#' @param de list. The output of LorgnetteDE() or NULL. If NULL only the obj data
#'   wil be subsetted.
#' @param branch numeric. The branch to be extracted for downstream analysis.
#'   Default is 1. It can never take NULL value.
#' @param branch2 numeric. The second (set of) branch(es) that may need to be extracted for
#'   downstream analysis (the case of multi-branch analysis). If NULL (default) or
#'   branch=branch2, the data of a single branch will only be selected.
#' @param joined logical. If TRUE, the joined data of all branches are selected. Thus, the branch
#'   parameter is not considered. This is different from specifying all data from the branch parameter
#'   because it will select the differentially expressed genes branch-wise. Default is FALSE.
#' @keywords selectData
#' @return A list with the selected data of the Lorgnette object and the differential
#'   expression analysis variable (if any).
#'
selectData<-function(obj,de,branch=1,branch2=NULL,joined=FALSE){
    
    if(joined){
    if(match("Joined_Branches",names(de),nomatch=0)==0){
       stop("There is no differential expressin analysis for all joined branches. Run MonocleDE() with joinedDE = TRUE!")
    }
    }

    if(joined & !is.null(branch2)){
    stop("Selecting all data ignores branch specification. Set branch2 = NULL or joinedDE = FALSE")
    }

    if(is.null(branch) & !joined){
        stop("Parameter branch must be defined to select the appropriate single-branch data (branch2 should be NULL).")
    }

    if(!is.null(branch)){
        mm<-match(branch,pData(obj)$Branch,nomatch=0)
        if(sum(mm)==0){
       stop("The value of branch parameter is not in the pData(obj)$Branch!")
        }
    }

    if(!is.null(branch2)){
    if(any(branch==branch2)){
        print("Parameters branch and branch2 cannot take the same value. Setting branch2 to a different value or NULL!")
        branch2<-branch2[which(branch!=branch2)]
        if(length(branch2)==0){
            branch2<-NULL
        }
    }
    }
    
    if(is.null(branch2)){

    if(!is.null(branch)){
       print(paste("Only the data of branch ",branch," are extracted",sep=""))
        }
    if(is.null(branch) & joined){
           print(paste("All data are extracted",sep=""))
       branch<-sort(unique(pData(obj)$Branch))
        }

        mm<-match(branch,pData(obj)$Branch,nomatch=0)
        mm.in<-branch[mm>0]
        mm.out<-branch[mm==0]
        if(length(mm.in)==0){
            stop("None of the specified branches is in the pData(obj)$Branch!")
        }
        if(length(mm.out)>0){
            print(paste("Branch ",paste(mm.out,collapse=", ")," cannot be found in the pData(obj)$Branch!",sep=""))
        }
        obj<-obj[,match(pData(obj)$Branch,branch,nomatch=0)>0]
            
        if(!is.null(de)){
        if(!joined){
                de<-de[match(paste("Branch",branch,sep=""),names(de))]
        } else {
        de<-de[match("Joined_Branches",names(de))]
        }
        }
    } else {
            
    branch2<-c(branch,branch2)
    mm<-match(branch2,pData(obj)$Branch,nomatch=0)
        mm.in<-branch2[mm>0]
        mm.out<-branch2[mm==0]
        if(length(mm.in)==0){
            stop("None of the specified branches is in the pData(obj)$Branch!")
        }
        if(length(mm.out)>0){
        stop("One of branch or branch2 does not exist in the pData(obj)$Branch!")
        }
    
    print(paste("The data of branches ",paste(branch2,collapse=", ")," are extracted",sep=""))
        obj<-obj[,match(pData(obj)$Branch,branch2,nomatch=0)>0]
            
        if(!is.null(de)){
            cc<-apply(t(combn(branch2,2)),1,paste,collapse="&")
            de<-de[match(paste("Branch",cc,sep=""),names(de))]
        }
            
    }
   return(list(Obj=obj,DE=de))
}


#' Selects the data to be analyzed downstream od the simpleHeat() analysis.
#'
#' Selects the data to be analyzed downstream od the LorgnetteDE() analysis. It is
#'   useful in case where parameter branch of LorgnetteDE() took a vector of values.
#'   This function wil select the results of interest for downstream analysis.
#' @param obj object. A Lorgnette object. Typically, the output of simpleHeat().
#' @param clusters numeric. If simpleHeat() is run for a series of specified clusters,
#'   this parameter specifies the number of clusters selected. If NULL (default), the
#'   minimum number of clusters from simpleHeat() will be selected.
#' @keywords selectClusters
#' @return A Lorgnette object with the selected clusters from simpleHeat().
#'
selectClusters<-function(obj,clusters=NULL){

    if(!is.null(clusters)){
        cl<-paste("Clusters",clusters,sep="")
    } else {
    cl<-names(obj$Obj)[1]
    }
    mm<-match(cl,names(obj$Obj),nomatch=0)
    mm.in<-cl[mm>0]
    mm.out<-cl[mm==0]
    if(length(mm.in)==0){
        stop("None of the specified clustering is in the data!")
    }
    if(length(mm.out)>0){
        print(paste("Clustering ",paste(mm.out,collapse=", ")," cannot be found in the data!",sep=""))
    }
    if(length(mm.in)>1){
        print(paste("Only the first clustering ",mm.in[1]," will be selected!",sep=""))
    }
    for(i in 1:length(obj)){
        if(i!=2){
            obj[[i]]<-obj[[i]][[match(mm.in[1],names(obj[[i]]))]]
        }
    }

 return(obj)
}




#' Runs Lorgnette's DE analysis for a predefined pseudotime.
#'
#' Runs Lorgnette's DE analysis for a predefined pseudotime, branch (state) and full / reduced models.
#' @param obj object. A Lorgnette object.
#' @param num.exprs.samples numeric. The minimum number of expressed samples (in bulk) or cells (in single-cell)
#'   a gene should have to be included in the analysis. Default is 10.
#' @param branch numeric. The branch (state) to be analyzed. This is defined in LorgnetteSet(). Default is 1.
#'   If it is a vector of values, then all the specified branches will be anayzed and stored as different
#'   componnts of a list.
#' @param full.model character. The full model of the DE analysis. Default is ~sm.ns(Pseudotime,df=3).
#'   More terms can be added to serve as control variables.
#' @param reduced.model character. The reduced model against which the full model will be contrasted.
#'   Default is ~1. Adding extra terms should be treated with care as it depends on the reserch questions
#'   and the full model specification.
#' @param joinedDE logical. If TRUE, the states are joined into a single one and the analysis is done on all
#'   data. Default is FALSE.
#' @keywords LorgnetteDE
#' @return A list of data frames with the DE estimates for each branch.
#'
LorgnetteDE<-function(obj,num.exprs.samples=10,branch=1,branch2=NULL,
                    full.model="~sm.ns(Pseudotime,df=3)",reduced.model="~1",joinedDE=FALSE){
 
    #pData(obj)$Branch<-paste(pData(obj)$State,"_",pData(obj)$Branch,sep="")
    if(joinedDE & !is.null(branch2)){
    stop("Running differential expression to all data ignores branch specification. Set branch2 = NULL or joinedDE = FALSE")
    }

    if(is.null(branch) & !joinedDE){
        stop("When joinedDE = FALSE, Parameter branch must be defined to select the appropriate single-branch data (branch2 should be NULL).")
    }

    if(!is.null(branch)){
        mm<-match(branch,pData(obj)$Branch,nomatch=0)
        if(sum(mm)==0){
           stop("The value of branch parameter is not in the pData(obj)$Branch!")
        }
    }

    if(!is.null(branch2)){
        if(any(branch==branch2)){
            print("Parameters branch and branch2 cannot take the same value. Setting branch2 to a different value or NULL!")
            branch2<-branch2[which(branch!=branch2)]
            if(length(branch2)==0){
                branch2<-NULL
            }
        }
    }
   
    if(!is.null(branch2)){
    bb<-c(branch,branch2)
    } else {
    bb<-branch
    }
    mm<-match(bb,pData(obj)$Branch,nomatch=0)
    mm.in<-bb[mm>0]
    mm.out<-bb[mm==0]
    if(length(mm.in)==0 & !joinedDE){
        stop("None of the specified branches is in the pData(obj)$Branch!")
    }
    if(length(mm.out)>0 & !joinedDE){
        print(paste("Branch ",paste(mm.out,collapse=", ")," cannot be found in the pData(obj)$Branch!",sep=""))
    }
    
    if(!is.null(branch2)){
    branch<-branch[match(mm.in,branch,nomatch=0)]
    branch2<-branch2[match(mm.in,branch2,nomatch=0)]
    mm.in<-matrix(0,1,2)
    for(i in 1:length(branch)){
      for(j in 1:length(branch2)){
        mm.in<-rbind(mm.in,matrix(c(branch[i],branch2[j]),nrow=1))
      }
    }
    mm.in<-matrix(mm.in[-1,],ncol=2)
    }
    DE<-NULL
    if(length(mm.in)>0){
        if(!is.matrix(mm.in)){
           DE<-as.list(rep(0,length(mm.in)))
           names(DE)<-paste("Branch",mm.in,sep="")
           for(i in 1:length(mm.in)){
        
               obj1 <- detectGenes(obj[,pData(obj)$Branch==mm.in[i]], min_expr = 0.1)
                expressed_genes <-  row.names(subset(fData(obj1),num_cells_expressed >= num.exprs.samples))
        
                print(paste("Now running MonocleDE for branch ",mm.in[i],"...",sep=""))
                de<-differentialGeneTest(obj1[expressed_genes,],
                                     fullModelFormulaStr = full.model,
                                     reducedModelFormulaStr=reduced.model)
                de<-de[,c(5:ncol(de),3:4)]
                DE[[i]]<-de
           }
       } else {
      DE<-as.list(rep(0,nrow(mm.in)))
          names(DE)<-paste("Branch",apply(mm.in,1,paste,collapse="&"),sep="")
          for(i in 1:nrow(mm.in)){

                obj1 <- detectGenes(obj[,match(pData(obj)$Branch,mm.in[i,],nomatch=0)>0], min_expr = 0.1)
                expressed_genes <-  row.names(subset(fData(obj1),num_cells_expressed >= num.exprs.samples))

                print(paste("Now running LorgnetteDE for branch ",mm.in[i,1]," vs ",mm.in[i,2],"...",sep=""))
                de<-differentialGeneTest(obj1[expressed_genes,],
                                     fullModelFormulaStr = full.model,
                                     reducedModelFormulaStr=reduced.model)
                de<-de[,c(5:ncol(de),3:4)]
                DE[[i]]<-de
          }
       }
    }

    if(joinedDE){
    obj1<-detectGenes(obj, min_expr = 0.1)
    expressed_genes <-  row.names(subset(fData(obj1),num_cells_expressed >= num.exprs.samples))
    print(paste("Now running LorgnetteDE for all data (joined branches ",paste(mm.in,collapse=","),")...",sep=""))
            de<-differentialGeneTest(obj1[expressed_genes,],
                                     fullModelFormulaStr = full.model,
                                     reducedModelFormulaStr=reduced.model)
            de<-de[,c(5:ncol(de),3:4)]
    DE<-c(DE,list(de))
    names(DE)[length(DE)]<-"Joined_Branches"
    }

 return(DE)
}





#' Converts a DESeqDataSetFromMatrix object into a Monocle2 object
#'
#' Converts a DESeqDataSetFromMatrix object into a Monocle2 object which is utilized by Lorgnette
#'   for the differential expression analysis via LorgnetteDE().
#' @param obj object. A DESeqDataSetFromMatrix object
#' @keywords DESeq2Monocle
#' @return A Monocle2 SingleCellExperiment object
#'
DESeq2Monocle<-function(obj){

   pd<-data.frame(colData(obj))
   fd<-data.frame(rowData(obj)[,which(colnames(rowData(obj))=="Gene"):ncol(rowData(obj))])
   rownames(fd)<-fd[,1]
   pd <- new("AnnotatedDataFrame", data = pd)
   fd <- new("AnnotatedDataFrame", data = fd)
   obj <- newCellDataSet(as.matrix(counts(obj)),
                         phenoData = pd,
                         featureData = fd,
                         expressionFamily=negbinomial.size())

   obj <- estimateSizeFactors(obj)
   obj <- estimateDispersions(obj)

 return(obj)

}




#' Generates a heatmap of single-branch Monocle's DE genes and their cluster membership.
#'
#' Generates a heatmap of single-branch Monocle's DE genes and their cluster membership.
#' @param obj object. A Monocle object.
#' @param de data frame. The DE statistics generated by MonocleDE().
#' @param id character. An ID to be used in the filename of the generated heatmap.
#' @param heat.clusters numeric. The number of gene clusters to show on the heatmap.
#' @param outdir character. A folder to store the heatmap.
#' @param trend_formula character. A formula specifying the model formula used in fitting the spline curve
#'   for each gene/feature. Default is ~sm.ns(Pseudotime,df=3).
#' @param hmcols character. The color scheme for drawing the heatmap. Default is NULL that uses the blue2green2red
#'   scheme of the R package colorRamps.
#' @param cluster_rows logical. Whether to cluster the rows of the heatmap. Default is TRUE.
#' @param show_rownames logical. Whether to show the names and the membership for each gene on the heatmap.
#'   Default is FALSE.
#' @param clustering_method character. A gene / sample clustering method. Accepts the same values as hclust.
#'   Default is ward.D2.
#' @param add_annotation_row data frame.  Additional annotations to show for each row in the heatmap. Must be
#'   a data frame with one row for each row in the fData table of obj with matching IDs. Default is NULL.
#' @param add_annotation_col data frame.  Additional annotations to show for each row in the heatmap. Must be
#'   a data frame with one row for each row in the pData table of obj with matching IDs. Default is NULL.
#' @param scale.range numeric. The minimum and the maximum value (in standard deviations) to show in the
#'   heatmap. Values smaller / larger than this are set to the min / max. Default is c(-3,3).
#' @param plot.width numeric. The width of the graphics region in inches. Default is 7.
#' @param plot.height numeric. The height of the graphics region in inches. Default is 7.
#' @param screen.it logical. If TRUE, the plot will be shown on screen and not stored. Default is FALSE.
#' @keywords simpleHeat
#' @return A list with components a Monocle object with the DE genes only, a data frame of normalized and smoothed
#'   gene expression profiles (all data and DE only), the list of DE genes with the top cluster memberships and the
#'   full memberships matrix.
#'
simpleHeat<-function(obj,de,id,heat.clusters,outdir,
                     trend_formula="~sm.ns(Pseudotime,df=3)",
                     scale.range=c(-3,3),hmcols=NULL,cluster_rows=TRUE,
                     show_rownames=FALSE,clustering_method="ward.D2",
                     add_annotation_row=NULL,add_annotation_col=NULL,
                     plot.width=7,plot.height=7,screen.it=FALSE){
    
    if(!dir.exists(outdir)){
        dir.create(outdir)
        print(paste("Folder ",outdir," has been created to store the simpleHeat() plots.",sep=""))
    }
    normAll<-smoother(obj,trend_formula=trend_formula,scale.range=scale.range)
    ii<-intersect(as.character(rownames(normAll)),as.character(de[,1]))
    norm<-normAll[match(ii,as.character(rownames(normAll))),]
    de<-de[match(ii,as.character(de[,1])),]
    obj<-obj[match(de$Gene,rownames(obj)),]

    tt<-table(GeneNames(rownames(obj)))
    tt<-which(tt>1)
    if(length(tt)>0){
        out<-c()
        for(i in 1:length(tt)){
            g<-grep(paste(":",names(tt[i]),"$",sep=""),rownames(obj))
            d1<- obj@assayData$exprs[g,]
            a<-apply(d1,1,mean)
            out<-c(out,g[which.min(a)])
        }
        obj<-obj[-out,]
    }
    
    heat.clusters<-sort(heat.clusters)
    list.de<-list.normAll<-list.norm<-list.obj<-list.membership<-as.list(rep(0,length(heat.clusters)))
    for(k in 1:length(heat.clusters)){
    print(paste("Generating heatmap and ",heat.clusters[k]," clusters for ",nrow(obj)," genes...",sep=""))
        p<-plot_pseudotime_heatmap(obj,num_clusters = heat.clusters[k],return_heatmap=T)
        dev.off()
        clus <- as.data.frame(cutree(p$tree_row, k=heat.clusters[k]))
        colnames(clus) <- "Cluster"
        clus$Gene <- rownames(clus)
        clus<-clus[,2:1]
        mm<-match(de[,1],clus[,1],nomatch=0)
        list.de[[k]]<-cbind(de[mm>0,],Cluster=as.numeric(as.character(clus[mm,2])),k=heat.clusters[k])
        
        norm1<-norm[match(list.de[[k]]$Gene,rownames(norm)),]
        
        cc<-matrix(0,heat.clusters[k],ncol(norm1))
        for(i in 1:heat.clusters[k]){
           nn<-matrix(norm1[list.de[[k]]$Cluster==i,],ncol=ncol(norm1))
           cc[i,]<-apply(nn,2,mean)
        }
        rownames(cc)<-1:heat.clusters[k]
    
        dm <- sapply(seq_len(nrow(norm1)),
                     function(i) apply(cc, 1, function(v) sqrt(sum((norm1[i, ]-v)^2))))

        m <- 2
        ms <- t(apply(dm, 2,function(x) {
                                     tmp <- 1/((x/sum(x))^(2/(m-1)))  # formula above
                                     tmp/sum(tmp)  # normalization
                                   }))
        rownames(ms)<-list.de[[k]][,1]
    
        memb<-matrix(0,nrow(ms),3)
        for(i in 1:nrow(list.de[[k]])){
            a<-ms[i,]
            memb[i,]<-c(a[as.numeric(as.character(list.de[[k]]$Cluster[i]))],
                        as.numeric(as.character(names(sort(a[-as.numeric(as.character(list.de[[k]]$Cluster[i]))],decreasing=T)[1]))),
                    sort(a[-as.numeric(as.character(list.de[[k]]$Cluster[i]))],decreasing=T)[1])
        }
        memb<-data.frame(memb)
        colnames(memb)<-c("Membership","Cluster2","Membership2")
        rownames(memb)<-list.de[[k]][,1]
        list.de[[k]]<-cbind(list.de[[k]],memb)
    
        if(!screen.it){
           pdf(paste(outdir,"DE_Heatmap_",id,"_",heat.clusters[k],"clusters.pdf",sep=""),height=plot.height,width=plot.width)
              print(plot_pseudotime_heatmap(obj,cluster_rows=cluster_rows,hclust_method=clustering_method,
                                            hmcols = hmcols,add_annotation_row=add_annotation_row,add_annotation_col=add_annotation_col,
                                            show_rownames=show_rownames,use_gene_short_name=TRUE,norm_method="log",
                                            scale_max=scale.range[2],scale_min=scale.range[1],trend_formula=trend_formula,
                                            num_clusters = heat.clusters[k],return_heatmap=FALSE))
           dev.off()
    
           # message
           print(paste("Download plot: ",outdir,"DE_Heatmap_",id,"_",heat.clusters[k],"clusters.pdf",sep=""))
        } else {

       print(plot_pseudotime_heatmap(obj,cluster_rows=cluster_rows,hclust_method=clustering_method,
                                            hmcols = hmcols,add_annotation_row=add_annotation_row,add_annotation_col=add_annotation_col,
                                            show_rownames=show_rownames,use_gene_short_name=TRUE,norm_method="log",
                                            scale_max=scale.range[2],scale_min=scale.range[1],trend_formula=trend_formula,
                                            num_clusters = heat.clusters[k],return_heatmap=FALSE))
    
    }
        mm<-match(rownames(obj)[p$tree_row$order],rownames(obj))
        list.obj[[k]]<-obj[mm,]
        mm<-match(rownames(list.obj[[k]]),rownames(norm1),nomatch=0)
    list.obj[[k]]<-list.obj[[k]][mm>0,]
        list.norm[[k]]<-norm1[mm,]
        list.de[[k]]<-list.de[[k]][mm,]
        list.membership[[k]]<-memb[mm,]
 
    }
    names(list.obj)<-names(list.norm)<-names(list.de)<-names(list.membership)<-paste("Clusters",heat.clusters,sep="")


  return(list(Obj=list.obj,NormDataAll=normAll,NormData=list.norm,DE=list.de,Membership=list.membership))
}





#' Generates the smoothed data for a set of single branch Monocle DE genes.
#'
#' Generates the smoothed data for a set of single branch Monocle DE genes.
#' @param data object. A Monocle object.
#' @param trend_formula character. A formula specifying the model formula used in fitting the spline curve
#'   for each gene/feature.
#' @param scale.range numeric. The minimum and the maximum value (in standard deviations) to show in the
#'   heatmap. Values smaller / larger than this are set to the min / max.
#' @keywords smoother
#' @return A data frame of smoothed and normalized values for heatmap plotting and membership estimation.
#'
smoother<-function(data,trend_formula,scale.range){

    newdata_prep <- data.frame(Pseudotime = seq(min(pData(data)$Pseudotime),
                                                max(pData(data)$Pseudotime), length.out = nrow(pData(data))),
                               Pseudotime_old = pData(data)$Pseudotime)
    uu<-unique(newdata_prep$Pseudotime_old)
    newdata<-c()
    for(i in 1:length(uu)){
        d<-newdata_prep[newdata_prep$Pseudotime_old==uu[i],]
        newdata<-c(newdata,rep(mean(d[,1]),nrow(d)))
    }
    newdata<-data.frame(Pseudotime = newdata)
    m <- genSmoothCurves(data, cores = 1, trend_formula = trend_formula, relative_expr = T, new_data = newdata)
    m = m[!apply(m, 1, sum) == 0, ]
    m = log2(m + 1)
    m = m[!apply(m, 1, sd) == 0, ]
    m = Matrix::t(scale(Matrix::t(m), center = TRUE))
    m = m[is.na(row.names(m)) == FALSE, ]
    m[is.nan(m)] = 0
    m[m > scale.range[2]] = scale.range[2]
    m[m < scale.range[1]] = scale.range[1]

  return(m)
}



#' Generates the smoothed data for two-branch Monocle DE genes.
#'
#' Generates the smoothed data for two-branch Monocle DE genes.
#' @param data object. A Monocle object.
#' @param trend_formula character. A formula specifying the model formula used in fitting the spline curve
#'   for each gene/feature.
#' @keywords branchSmoother
#' @return A data frame of smoothed and normalized values for heatmap plotting and membership estimation.
#'
branchSmoother<-function(obj,trend_formula){

    ww<-which(colnames(pData(obj))=="Branch")
    if(length(ww)==0){
        stop("The Branch variable is missing from pData(obj). Run MonocleBranchedSet() first!")
    }
    br<-table(pData(obj)$Branch)
    #if(length(br)>2){
    #    stop("Only two branches can be compared")
    #}
    if(length(grep("*Branch", trend_formula))==0){
        stop("For branch comparison the trend_formula should include the branch interaction!")
    }


    newdata<-as.list(rep(0,length(br)))
    for(i in 1:length(br)){
    data<-obj[,pData(obj)$Branch==as.numeric(names(br)[i])]
        newdata_prep <- data.frame(Pseudotime = seq(min(pData(data)$Pseudotime),
                            max(pData(data)$Pseudotime), length.out = nrow(pData(data))),
                                  Pseudotime_old = pData(data)$Pseudotime)
        uu<-unique(newdata_prep$Pseudotime_old)
        newdataA<-c()
        for(j in 1:length(uu)){
            d<-newdata_prep[newdata_prep$Pseudotime_old==uu[j],]
            newdataA<-c(newdataA,rep(mean(d[,1]),nrow(d)))
        }
        newdata[[i]]<-data.frame(Pseudotime = as.numeric(as.character(newdataA)),Branch = as.factor(rep(names(br)[i],br[i])))
    }

    if(length(br)==2){
        Branch_exprs <- genSmoothCurves(obj, cores = 1, trend_formula = trend_formula, relative_expr = T, new_data = do.call(rbind,newdata))
        BranchA_exprs <- Branch_exprs[, 1:br[1]]
        BranchB_exprs <- Branch_exprs[, (br[1]+1):sum(br)]
        BranchA_exprs <- log2(BranchA_exprs + 1)    # ---> normalized data A
        BranchB_exprs <- log2(BranchB_exprs + 1)    # ---> normalized data B
        Branch_exprs<-list(A=BranchA_exprs,B=BranchB_exprs,BranchInfo=br)
    } else {
       if(length(grep("*Branch", trend_formula))>0){
        trend_formula<-gsub("\\*Branch","",trend_formula)
       }
      
       start<-1
       end<-br[1]
       Branch_exprs<-as.list(rep(0,length(br)))
       for(i in 1:length(br)){
           Be<-genSmoothCurves(obj[,start:end], cores = 1, trend_formula = trend_formula, relative_expr = T, new_data = newdata[[i]])
           Be<-log2(Be + 1)
           Branch_exprs[[i]]<-Be
           start<-br[i]+1
           end<-sum(br[1:(i+1)])
        }
       Branch_exprs<-c(Branch_exprs,list(br))
       names(Branch_exprs)<-c(LETTERS[1:length(br)],"BranchInfo")
    }
       
 return(Branch_exprs)
}



#' Generates a heatmap of two-branch Monocle's DE genes and their cluster membership.
#'
#' Generates a heatmap of two-branch Monocle's DE genes and their cluster membership.
#' @param obj object. A Monocle object.
#' @param de data frame. The DE statistics generated by MonocleDE().
#' @param id character. An ID to be used in the filename of the generated heatmap.
#' @param heat.clusters numeric. The number of gene clusters to show on the heatmap.
#' @param outdir character. A folder to store the heatmap.
#' @param branch_labels character. Labels for the two branches.
#' @param branch_colors character. Colors for the two branches. Default is c("#F05662","#7990C8").
#' @param trend_formula character. A formula specifying the model formula used in fitting the spline curve
#'   for each gene/feature. Default is ~sm.ns(Pseudotime,df=3).
#' @param reverse.time logical. If TRUE (default), the times of the first branch are reversed to mimic
#'   the default Monocle branch analysis. It is only used when teo branches are compared. For the comparison
#'   of more than two branches the times are not reversed.
#' @param hmcols character. The color scheme for drawing the heatmap. Default is NULL that uses the blue2green2red
#'   scheme of the R package colorRamps.
#' @param cluster_rows logical. Whether to cluster the rows of the heatmap. Default is TRUE.
#' @param show_rownames logical. Whether to show the names and the membership for each gene on the heatmap.
#'   Default is FALSE.
#' @param clustering_method character. A gene / sample clustering method. Accepts the same values as hclust.
#'   Default is ward.D2.
#' @param add_annotation_row data frame.  Additional annotations to show for each row in the heatmap. Must be
#'   a data frame with one row for each row in the fData table of obj with matching IDs. Default is NULL.
#' @param add_annotation_col data frame.  Additional annotations to show for each row in the heatmap. Must be
#'   a data frame with one row for each row in the pData table of obj with matching IDs. Default is NULL.
#' @param scale.range numeric. The minimum and the maximum value (in standard deviations) to show in the
#'   heatmap. Values smaller / larger than this are set to the min / max. Default is c(-3,3).
#' @param plot.width numeric. The width of the graphics region in inches. Default is 7.
#' @param plot.height numeric. The height of the graphics region in inches. Default is 7.
#' @param screen.it logical. If TRUE, the plot will be shown on screen and not stored. Default is FALSE.
#' @keywords branchHeat
#' @return A list with components a Monocle object with the DE genes only, a data frame of normalized and smoothed
#'   gene expression profiles (all data and DE only), the list of DE genes with the top cluster memberships and the
#'   full memberships matrix.
#'
branchHeat<-function(obj,de=de,id,heat.clusters,outdir,branch_labels,
                     branch_colors=c("#F05662","#7990C8"),scale.range=c(-3,3),
                     trend_formula="~sm.ns(Pseudotime,df=3)*Branch",reverse.time=TRUE,
                     cluster_rows=TRUE,show_rownames=FALSE,clustering_method="ward.D2",
                     add_annotation_row=NULL,add_annotation_col=NULL,hmcols=NULL,
                     plot.width=7,plot.height=7,screen.it=FALSE){

    dat<-branchSmoother(obj=obj,trend_formula=trend_formula)
    if(is.data.frame(de)){
    stop("Parameter de must be a list!")
    }
    

    Branch_exprs<-as.list(rep(0,(length(dat)-1)))
    names(Branch_exprs)<-names(dat)[1:(length(dat)-1)]
    if(length(de)==1){

        for(i in 1:length(Branch_exprs)){
            Branch_exprs[[i]]<-dat[[i]][match(de[[1]][,1],rownames(dat[[i]]),nomatch=0),]
        }
      obj<-obj[match(unique(de[[1]][,1]),rownames(obj),nomatch=0),]
  
    } else {
    
    gg<-unique(do.call(rbind,de)$Gene)
    for(i in 1:length(Branch_exprs)){
            Branch_exprs[[i]]<-dat[[i]][match(gg,rownames(dat[[i]]),nomatch=0),]
        }
    obj<-obj[match(gg,rownames(obj),nomatch=0),]

    }

    if(length(branch_colors)!=length(Branch_exprs)){
    stop("Use a branch color for each branch!")
    }
    if(length(branch_labels)!=length(Branch_exprs)){
        stop("Use a branch label for each branch!")
    }


    tt<-table(GeneNames(rownames(obj)))
    tt<-which(tt>1)
    if(length(tt)>0){
        out<-c()
        for(i in 1:length(tt)){
            g<-grep(paste(":",names(tt[i]),"$",sep=""),rownames(obj))
            d1<- obj@assayData$exprs[g,]
            a<-apply(d1,1,mean)
            out<-c(out,g[which.min(a)])
        }
        obj<-obj[-out,]
    }

    if(length(Branch_exprs)==2){
    if(reverse.time){
            print("Reversing pseudotimes for the first branch!")
            Branch_exprs[[1]]<-Branch_exprs[[1]][,ncol(Branch_exprs[[1]]):1]
        }
    }
    
    heatmap_matrix <- do.call(cbind,Branch_exprs)
    heatmap_matrix <- heatmap_matrix[!apply(heatmap_matrix, 1, sd) == 0, ]
    heatmap_matrix <- Matrix::t(scale(Matrix::t(heatmap_matrix),center = TRUE))
    heatmap_matrix <- heatmap_matrix[is.na(row.names(heatmap_matrix)) == FALSE, ]
    heatmap_matrix[is.nan(heatmap_matrix)] <- 0
    heatmap_matrix[heatmap_matrix > scale.range[2]] <- scale.range[2]
    heatmap_matrix[heatmap_matrix < scale.range[1]] <- scale.range[1]
    heatmap_matrix_ori <- heatmap_matrix
    heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[,1]) & is.finite(heatmap_matrix[, ncol(heatmap_matrix)]), ]

    row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
    row_dist[is.na(row_dist)] <- 1
    exp_rng <- range(heatmap_matrix)
    bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by = 0.1)
    if(is.null(hmcols)){
        hmcols <- blue2green2red(length(bks) - 1)
    }
    
    heat.clusters<-sort(heat.clusters)
    list.de<-list.normAll<-list.norm<-list.obj<-list.membership<-as.list(rep(0,length(heat.clusters)))
    for(k in 1:length(heat.clusters)){
        print(paste("Generating heatmap and ",heat.clusters[k]," clusters for ",nrow(heatmap_matrix)," genes...",sep=""))
        
        ph <- pheatmap(heatmap_matrix,useRaster=TRUE,cluster_cols=FALSE,
                       cluster_rows=cluster_rows,show_rownames=show_rownames,show_colnames=FALSE,
                       clustering_method=clustering_method,clustering_distance_rows = row_dist,
                       cutree_rows = heat.clusters, silent = TRUE, filename = NA,
                       breaks = bks, color = hmcols)

        # membership estimation!
        clus <- as.data.frame(cutree(ph$tree_row, k=heat.clusters[k]))
        colnames(clus) <- "Cluster"
        clus$Gene <- rownames(clus)
        clus<-clus[,2:1]
    if(length(de)==1){
        mm<-match(unique(de[[1]][,1]),clus[,1],nomatch=0)
            list.de[[k]]<-cbind(de[[1]][mm>0,],Cluster=as.numeric(as.character(clus[mm,2])),k=heat.clusters[k])
    } else {
        mm<-match(unique(de[[1]][,1]),clus[,1],nomatch=0)
            list.de[[k]]<-cbind(de[[1]][mm>0,],Cluster=as.numeric(as.character(clus[mm,2])),k=heat.clusters[k],Branch=names(de)[1])
        for(i in 2:length(de)){
        mm<-match(unique(de[[i]][,1]),clus[,1],nomatch=0)
        list.de[[k]]<-rbind(list.de[[k]],cbind(de[[i]][mm>0,],Cluster=as.numeric(as.character(clus[mm,2])),k=heat.clusters[k],Branch=names(de)[i]))
        }
    }
        
        heatmap_matrix1<-heatmap_matrix[match(unique(list.de[[k]]$Gene),rownames(heatmap_matrix)),]
        
        cc<-matrix(0,heat.clusters[k],ncol(heatmap_matrix1))
    xx<-unique(cbind(as.character(list.de[[k]]$Gene),as.character(list.de[[k]]$Cluster)))
        xx<-xx[match(rownames(heatmap_matrix1),xx[,1]),]
        for(i in 1:heat.clusters[k]){
            nn<-matrix(heatmap_matrix1[as.numeric(xx[,2])==i,],ncol=ncol(heatmap_matrix1))
            cc[i,]<-apply(nn,2,mean)
        }
        rownames(cc)<-1:heat.clusters[k]

        dm <- sapply(seq_len(nrow(heatmap_matrix1)),
                    function(i) apply(cc, 1, function(v) sqrt(sum((heatmap_matrix1[i, ]-v)^2))))

        m <- 2
        ms <- t(apply(dm, 2,function(x) {
                                tmp <- 1/((x/sum(x))^(2/(m-1)))  # formula above
                                tmp/sum(tmp)  # normalization
                                    }))
        rownames(ms)<-xx[,1]

        memb<-matrix(0,nrow(ms),3)
        for(i in 1:nrow(xx)){
            a<-ms[i,]
            memb[i,]<-c(a[as.numeric(as.character(xx[i,2]))],
                          as.numeric(as.character(names(sort(a[-as.numeric(as.character(xx[i,2]))],decreasing=T)[1]))),
                          sort(a[-as.numeric(as.character(xx[i,2]))],decreasing=T)[1])
        }
        memb<-data.frame(memb)
        colnames(memb)<-c("Membership","Cluster2","Membership2")
        rownames(memb)<-xx[,1]
    xx<-cbind(xx,memb)
    xx<-xx[match(list.de[[k]][,1],xx[,1]),]
        list.de[[k]]<-cbind(list.de[[k]],xx[,3:ncol(xx)])

        annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row,heat.clusters[k])))
        if (!is.null(add_annotation_row)) {
            annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row),])
        }
        colnames(heatmap_matrix1) <- c(1:ncol(heatmap_matrix1))
        annotation_col <- data.frame(row.names = c(1:ncol(heatmap_matrix1)),Condition = rep(branch_labels, as.numeric(dat$BranchInfo)))

        colnames(annotation_col) <- "Condition"
        if(!is.null(add_annotation_col)) {
            annotation_col <- cbind(annotation_col, add_annotation_col[fData(obj[row.names(annotation_col),])$gene_short_name, 1])
        }

        names(branch_colors) <- branch_labels
        annotation_colors = list(Condition = branch_colors)
        names(annotation_colors$Condition) = branch_labels
        if(is.null(fData(obj)$gene_short_name) == FALSE) {
            feature_label <- as.character(fData(obj)[row.names(heatmap_matrix1),"gene_short_name"])
            feature_label[is.na(feature_label)] <- row.names(heatmap_matrix1)
            row_ann_labels <- as.character(fData(obj)[row.names(annotation_row),"gene_short_name"])
            row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
        } else {
            feature_label <- row.names(heatmap_matrix1)
            row_ann_labels <- row.names(annotation_row)
        }

        row.names(heatmap_matrix1) <- feature_label
        row.names(annotation_row) <- row_ann_labels
        if(!screen.it){
            pdf(paste(outdir,"DE_Heatmap_",id,"_",heat.clusters[k],"clusters.pdf",sep=""),width=plot.width,height=plot.height)
                ph_res <- pheatmap(heatmap_matrix1[, ],useRaster=TRUE,cluster_cols=FALSE,
                                   cluster_rows=cluster_rows,show_rownames=show_rownames,show_colnames=FALSE,
                                   clustering_distance_rows = row_dist, clustering_method=clustering_method,
                                   cutree_rows=heat.clusters[k],annotation_row=annotation_row,
                                   annotation_col=annotation_col,gaps_col = cumsum(dat$BranchInfo[1:(length(dat$BranchInfo)-1)]),breaks = bks,
                                   color = hmcols)
                grid::grid.rect(gp = grid::gpar("fill", col = NA))
                grid::grid.draw(ph_res$gtable)
                print(ph_res)
            dev.off()

         # message
         print(paste("Download plot: ",outdir,"DE_Heatmap_",id,"_",heat.clusters[k],"clusters.pdf",sep=""))
         
        } else {
         ph_res <- pheatmap(heatmap_matrix1[, ],useRaster=TRUE,cluster_cols=FALSE,
                            cluster_rows=cluster_rows,show_rownames=show_rownames,show_colnames=FALSE,
                            clustering_distance_rows = row_dist, clustering_method=clustering_method,
                            cutree_rows=heat.clusters[k],annotation_row=annotation_row,
                            annotation_col=annotation_col,gaps_col = cumsum(dat$BranchInfo[1:(length(dat$BranchInfo)-1)]),breaks = bks,
                            color = hmcols)
         grid::grid.rect(gp = grid::gpar("fill", col = NA))
         grid::grid.draw(ph_res$gtable)
         print(ph_res)

        }

        mm<-match(rownames(obj)[ph$tree_row$order],rownames(obj))
        list.obj[[k]]<-obj[mm,]
        mm<-match(GeneNames(rownames(list.obj[[k]])),rownames(heatmap_matrix1),nomatch=0)
        list.obj[[k]]<-list.obj[[k]][mm>0,]
        list.norm[[k]]<-heatmap_matrix1[mm,]
    ii<-intersect(list.de[[k]][,1],list.de[[k]][,1])
        list.de[[k]]<-list.de[[k]][match(list.de[[k]][,1],ii,nomatch=0)>0,]
        list.membership[[k]]<-memb[mm,]

    }
    names(list.obj)<-names(list.norm)<-names(list.de)<-names(list.membership)<-paste("Clusters",heat.clusters,sep="")
         
 return(list(Obj=list.obj,NormDataAll=dat,NormData=list.norm,DE=list.de,Membership=list.membership))
}




#' Prepare the dataset for plotting of selected genes and coloring based on gene names.
#'
#' Prepare the dataset for plotting of selected genes and coloring based on gene names.
#' @param Dat data frame. A data frame with the normalized, scaled and smoothed expression data. The genes
#'   are in columns. The last column contains the time variable.
#' @param relexpr data frame. A data frame with the normalized data. The genes are in columns. The last
#'   column contains the time variable.
#' @param de data frame. The differential expression results.
#' @param genes.to.plot character. The genes to be plotted. They must be included at the short_gene_name
#'   slot of the DE data frame.
#' @param colors character. A vector of colors.
#' @param timeVar character. A name for the x variable of the regression model.
#' @param highlightVar character. The data of the variable that highlights the points of the plot.
#' @keywords fun1
#' @return A list with components the Dat, de, genes.to.plot and the normalized expressions with the genes
#'   to be plotted.
#'
fun1<-function(Dat,relexpr,de,genes.to.plot,colors,timeVar,highlightVar){
    
    mm<-match(unique(genes.to.plot),colnames(Dat),nomatch=0)
    if(sum(mm)==0){
        stop("None of these genes is found in the data")
    }
    ww<-which(mm==0)
    if(length(ww)>0){
        print(paste("Genes ",paste(genes.to.plot[mm[ww]],sep=",")," do not exist in the data",sep=""))
    }
    
    de<-de[match(unique(genes.to.plot),de$gene_short_name,nomatch=0),]
    genes.to.plot<-de$gene_short_name
    names(genes.to.plot)<-colors[1:length(genes.to.plot)]

    dat<-Dat[,match(c(genes.to.plot,highlightVar,timeVar),colnames(Dat))]
    relexpr<-relexpr[,match(c(genes.to.plot,highlightVar,timeVar),colnames(relexpr))]
   
 return(list(dat,de,genes.to.plot,relexpr))
}


#' Prepare the dataset for plotting of selected genes and coloring based on membership score.
#'
#' Prepare the dataset for plotting of selected genes and coloring based on membership score.
#'   The colors are automatically specified.
#' @param Dat data frame. A data frame with the normalized, scaled and smoothed expression data. The genes
#'   are in columns. The last column contains the time variable.
#' @param relexpr data frame. A data frame with the normalized data. The genes are in columns. The last
#'   column contains the time variable.
#' @param de data frame. The differential expression results.
#' @param genes.to.plot character. The genes to be plotted. They must be included at the short_gene_name
#'   slot of the DE data frame.
#' @param memb.cut numeric. A cutoff above while the plotted genes are highlighted in #D54C4C (red).
#' @param timeVar character. A name for the x variable of the regression model.
#' @param highlightVar character. The data of the variable that highlights the points of the plot.
#' @keywords fun2
#' @return A list with components the Dat, de, genes.to.plot and the normalized expressions with the genes
#'   to be plotted.
#'
fun2<-function(Dat,relexpr,de,genes.to.plot,memb.cut,timeVar,highlightVar){
    
    mm<-match(unique(genes.to.plot),colnames(Dat),nomatch=0)
    if(sum(mm)==0){
        stop("None of these genes is found in the data")
    }
    ww<-which(mm==0)
    if(length(ww)>0){
        print(paste("Genes ",paste(genes.to.plot[mm[ww]],sep=",")," do not exist in the data",sep=""))
    }


    cols<-rep("#D6D2C4",length(unique(genes.to.plot)))
    de<-de[match(unique(genes.to.plot),de$gene_short_name,nomatch=0),]
    cols[as.numeric(as.character(de$Membership))>=memb.cut]<-"#D54C4C"
    genes.to.plot<-de$gene_short_name
    names(genes.to.plot)<-cols
    
    dat<-Dat[,match(c(genes.to.plot,highlightVar,timeVar),colnames(Dat))]
    relexpr<-relexpr[,match(c(genes.to.plot,highlightVar,timeVar),colnames(relexpr))]
  
 return(list(dat,de,genes.to.plot,relexpr))
}


#' Prepare the dataset for plotting of selected genes and coloring based on cluster number.
#'
#' Prepare the dataset for plotting of selected genes and coloring based on cluster number.
#' @param Dat data frame. A data frame with the normalized, scaled and smoothed expression data. The genes
#'   are in columns. The last column contains the time variable.
#' @param relexpr data frame. A data frame with the normalized data. The genes are in columns. The last
#'   column contains the time variable.
#' @param de data frame. The differential expression results.
#' @param genes.to.plot character. The genes to be plotted. They must be included at the short_gene_name
#'   slot of the DE data frame.
#' @param colors character. A vector of colors with names(colors) being the cluster numbers.
#' @param timeVar character. A name for the x variable of the regression model.
#' @param highlightVar character. The data of the variable that highlights the points of the plot.
#' @keywords fun3
#' @return A list with components the Dat, de, genes.to.plot and the normalized expressions with the genes
#'   to be plotted.
#'
fun3<-function(Dat,relexpr,de,genes.to.plot,colors,timeVar,highlightVar){
 
    mm<-match(unique(genes.to.plot),colnames(Dat),nomatch=0)
    if(sum(mm)==0){
        stop("None of these genes is found in the data")
    }
    ww<-which(mm==0)
    if(length(ww)>0){
        print(paste("Genes ",paste(genes.to.plot[mm[ww]],sep=",")," do not exist in the data",sep=""))
    }

    de<-de[match(unique(genes.to.plot),de$gene_short_name,nomatch=0),]
    cols<-colors[match(as.numeric(as.character(de$Cluster)),as.numeric(as.character(names(colors))))]
    genes.to.plot<-de$gene_short_name
    names(genes.to.plot)<-cols
    
    dat<-Dat[,match(c(genes.to.plot,highlightVar,timeVar),colnames(Dat))]
    relexpr<-relexpr[,match(c(genes.to.plot,highlightVar,timeVar),colnames(relexpr))]

 return(list(dat,de,genes.to.plot,relexpr))
}


#' Prepare the dataset for plotting of all cluster-specific genes and coloring based on membership
#'   cutoff score.
#'
#' Prepare the dataset for plotting of all cluster-specific genes and coloring based on membership
#'   cutoff score. The colors are automatically specified.
#' @param Dat data frame. A data frame with the normalized, scaled and smoothed expression data. The genes
#'   are in columns. The last column contains the time variable.
#' @param relexpr data frame. A data frame with the normalized data. The genes are in columns. The last
#'   column contains the time variable.
#' @param de data frame. The differential expression results.
#' @param cluster.to.plot. The cluster whose genes are to be plotted. They must be included at the Cluster
#'   slot of the DE data frame.
#' @param memb.cut numeric. A cutoff above while the plotted genes are highlighted in #D54C4C (red).
#' @param timeVar character. A name for the x variable of the regression model.
#' @param highlightVar character. The data of the variable that highlights the points of the plot.
#' @keywords fun4
#' @return A list with components the Dat, de, genes.to.plot and the normalized expressions with the genes
#'   to be plotted.
#'
fun4<-function(Dat,relexpr,de,cluster.to.plot,memb.cut,timeVar,highlightVar){
    
    de<-de[match(as.numeric(as.character(de$Cluster)),cluster.to.plot,nomatch=0)>0,]
    cols<-rep("#D6D2C4",nrow(de))
    cols[as.numeric(as.character(de$Membership))>=memb.cut]<-"#D54C4C"
    genes.to.plot<-de$gene_short_name
    names(genes.to.plot)<-cols
    genes.to.plot<-genes.to.plot[!duplicated(genes.to.plot)]

    mm<-match(genes.to.plot,colnames(Dat),nomatch=0)
    if(sum(mm)==0){
        stop("None of these genes is found in the data")
    }
    ww<-which(mm==0)
    if(length(ww)>0){
        print(paste("Genes ",paste(genes.to.plot[mm[ww]],sep=",")," do not exist in the data",sep=""))
    }
    
    dat<-Dat[,match(c(genes.to.plot,highlightVar,timeVar),colnames(Dat))]
    relexpr<-relexpr[,match(c(genes.to.plot,highlightVar,timeVar),colnames(relexpr))]

 return(list(dat,de,genes.to.plot,relexpr))
}


#' Prepare the dataset for plotting of all cluster-specific genes and coloring based on membership top x
#'   score.
#'
#' Prepare the dataset for plotting of all cluster-specific genes and coloring based on membership top x
#'   score. Activated only when membership.cut = NULL. The colors are automatically specified.
#' @param Dat data frame. A data frame with the normalized, scaled and smoothed expression data. The genes
#'   are in columns. The last column contains the time variable.
#' @param relexpr data frame. A data frame with the normalized data. The genes are in columns. The last
#'   column contains the time variable.
#' @param de data frame. The differential expression results.
#' @param cluster.to.plot. The cluster whose genes are to be plotted. They must be included at the Cluster
#'   slot of the DE data frame.
#' @param top numeric. The number of genes with the highest membership score to be highlighted in #D54C4C (red).
#' @param timeVar character. A name for the x variable of the regression model.
#' @param highlightVar character. The data of the variable that highlights the points of the plot.
#' @keywords fun5
#' @return A list with components the Dat, de, genes.to.plot and the normalized expressions with the genes
#'   to be plotted.
#'
fun5<-function(Dat,relexpr,de,cluster.to.plot,top,timeVar,highlightVar){
    
    de<-de[match(as.numeric(as.character(de$Cluster)),cluster.to.plot,nomatch=0)>0,]
    
    if(nrow(de)<top){
        top<-nrow(de)
    }
    
    cols<-rep("#D6D2C4",nrow(de))
    cols[sort.list(as.numeric(as.character(de$Membership)),decreasing=TRUE)[1:top]]<-"#D54C4C"
    genes.to.plot<-de$gene_short_name
    names(genes.to.plot)<-cols
    genes.to.plot<-genes.to.plot[!duplicated(genes.to.plot)]
   
    mm<-match(genes.to.plot,colnames(Dat),nomatch=0)
    if(sum(mm)==0){
        stop("None of these genes is found in the data")
    }
    ww<-which(mm==0)
    if(length(ww)>0){
        print(paste("Genes ",paste(genes.to.plot[mm[ww]],sep=",")," do not exist in the data",sep=""))
    }
    
    dat<-Dat[,match(c(genes.to.plot,highlightVar,timeVar),colnames(Dat))]
    relexpr<-relexpr[,match(c(genes.to.plot,highlightVar,timeVar),colnames(relexpr))]

 return(list(dat,de,genes.to.plot,relexpr))
}



#' Prepare the dataset for plotting of cluster-specific genes with the top x qval and coloring based on
#'   gene names.
#'
#' Prepare the dataset for plotting of cluster-specific genes with the top x qval and coloring based on
#'   gene names. Activated only when qval.cut = NULL.
#' @param Dat data frame. A data frame with the normalized, scaled and smoothed expression data. The genes
#'   are in columns. The last column contains the time variable.
#' @param relexpr data frame. A data frame with the normalized data. The genes are in columns. The last
#'   column contains the time variable.
#' @param de data frame. The differential expression results.
#' @param cluster.to.plot. The cluster whose genes are to be plotted. They must be included at the Cluster
#'   slot of the DE data frame.
#' @param top numeric. The number of genes with the highest qval to be plotted.
#' @param colors character. A vector of colors.
#' @param timeVar character. A name for the x variable of the regression model.
#' @param highlightVar character. The data of the variable that highlights the points of the plot.
#' @keywords fun6
#' @return A list with components the Dat, de, genes.to.plot and the normalized expressions with the genes
#'   to be plotted.
#'
fun6<-function(Dat,relexpr,de,cluster.to.plot,colors,top,timeVar,highlightVar){
    
    de<-de[match(as.numeric(as.character(de$Cluster)),cluster.to.plot,nomatch=0)>0,]
    genes.to.plot<-de$gene_short_name[sort.list(as.numeric(as.character(de$qval)))[1:top]]
    de<-de[match(genes.to.plot,de$gene_short_name),]
    if(nrow(de)<top){
        top<-nrow(de)
    }
    names(genes.to.plot)<-colors[1:top]
    genes.to.plot<-genes.to.plot[!duplicated(genes.to.plot)]

    mm<-match(genes.to.plot,colnames(Dat),nomatch=0)
    if(sum(mm)==0){
        stop("None of these genes is found in the data")
    }
    ww<-which(mm==0)
    if(length(ww)>0){
        print(paste("Genes ",paste(genes.to.plot[mm[ww]],sep=",")," do not exist in the data",sep=""))
    }
    
    dat<-Dat[,match(c(genes.to.plot,highlightVar,timeVar),colnames(Dat))]
    relexpr<-relexpr[,match(c(genes.to.plot,highlightVar,timeVar),colnames(relexpr))]
        
 return(list(dat,de,genes.to.plot,relexpr))
}


#' Prepare the dataset for plotting of cluster-specific genes with qval<cutoff and coloring based on
#'   gene names.
#'
#' Prepare the dataset for plotting of cluster-specific genes with qval<cutoff and coloring based on
#'   gene names.
#' @param Dat data frame. A data frame with the normalized, scaled and smoothed expression data. The genes
#'   are in columns. The last column contains the time variable.
#' @param relexpr data frame. A data frame with the normalized data. The genes are in columns. The last
#'   column contains the time variable.
#' @param de data frame. The differential expression results.
#' @param cluster.to.plot. The cluster whose genes are to be plotted. They must be included at the Cluster
#'   slot of the DE data frame.
#' @param qcut numeric. The qvalue cutoff.
#' @param colors character. A vector of colors.
#' @param timeVar character. A name for the x variable of the regression model.
#' @param highlightVar character. The data of the variable that highlights the points of the plot.
#' @keywords fun7
#' @return A list with components the Dat, de, genes.to.plot and the normalized expressions with the genes
#'   to be plotted.
#'
fun7<-function(Dat,relexpr,de,cluster.to.plot,colors,qcut,timeVar,highlightVar){
    
    de<-de[match(as.numeric(as.character(de$Cluster)),cluster.to.plot,nomatch=0)>0,]
    genes.to.plot<-de$gene_short_name[as.numeric(as.character(de$qval))<=qcut]
    de<-de[match(genes.to.plot,de$gene_short_name),]
    names(genes.to.plot)<-colors[1:nrow(de)]
    genes.to.plot<-genes.to.plot[!duplicated(genes.to.plot)]

    mm<-match(genes.to.plot,colnames(Dat),nomatch=0)
    if(sum(mm)==0){
        stop("None of these genes is found in the data")
    }
    ww<-which(mm==0)
    if(length(ww)>0){
        print(paste("Genes ",paste(genes.to.plot[mm[ww]],sep=",")," do not exist in the data",sep=""))
    }
    
    dat<-Dat[,match(c(genes.to.plot,highlightVar,timeVar),colnames(Dat))]
    relexpr<-relexpr[,match(c(genes.to.plot,highlightVar,timeVar),colnames(relexpr))]
        
 return(list(dat,de,genes.to.plot,relexpr))
}


#' Prepare the dataset for plotting of cluster-specific genes with qval<cutoff and coloring based on
#'   membership scores.
#'
#' Prepare the dataset for plotting of cluster-specific genes with qval<cutoff and coloring based on
#'   membership scores.
#' @param Dat data frame. A data frame with the normalized, scaled and smoothed expression data. The genes
#'   are in columns. The last column contains the time variable.
#' @param relexpr data frame. A data frame with the normalized data. The genes are in columns. The last
#'   column contains the time variable.
#' @param de data frame. The differential expression results.
#' @param cluster.to.plot. The cluster whose genes are to be plotted. They must be included at the Cluster
#'   slot of the DE data frame.
#' @param qcut numeric. The qvalue cutoff.
#' @param colors character. A vector of colors.
#' @param memb.cut numeric. A membership score cutoff. The genes that pass it will be highlighted in #D54C4C (red).
#' @param timeVar character. A name for the x variable of the regression model.
#' @param highlightVar character. The data of the variable that highlights the points of the plot.
#' @keywords fun8
#' @return A list with components the Dat, de, genes.to.plot and the normalized expressions with the genes
#'   to be plotted.
#'
fun8<-function(Dat,relexpr,de,cluster.to.plot,colors,qcut,memb.cut,timeVar,highlightVar){
    
    de<-de[match(as.numeric(as.character(de$Cluster)),cluster.to.plot,nomatch=0)>0,]
    genes.to.plot<-de$gene_short_name[as.numeric(as.character(de$qval))<=qcut]
    de<-de[match(genes.to.plot,de$gene_short_name),]
    
    cols<-rep("#D6D2C4",nrow(de))
    cols[as.numeric(as.character(de$Membership))>=memb.cut]<-"#D54C4C"
    names(genes.to.plot)<-cols
    genes.to.plot<-genes.to.plot[!duplicated(genes.to.plot)]
   
    mm<-match(genes.to.plot,colnames(Dat),nomatch=0)
    if(sum(mm)==0){
        stop("None of these genes is found in the data")
    }
    ww<-which(mm==0)
    if(length(ww)>0){
        print(paste("Genes ",paste(genes.to.plot[mm[ww]],sep=",")," do not exist in the data",sep=""))
    }
    
    dat<-Dat[,match(c(genes.to.plot,highlightVar,timeVar),colnames(Dat))]
    relexpr<-relexpr[,match(c(genes.to.plot,highlightVar,timeVar),colnames(relexpr))]
        
 return(list(dat,de,genes.to.plot,relexpr))
}


#' Prepare the dataset for plotting of cluster-specific genes with the top x qval and coloring based on
#'   membership scores.
#'
#' Prepare the dataset for plotting of cluster-specific genes with the top x qval and coloring based on
#'   membership scores. Activated only when qval.cut = NULL.
#' @param Dat data frame. A data frame with the normalized, scaled and smoothed expression data. The genes
#'   are in columns. The last column contains the time variable.
#' @param relexpr data frame. A data frame with the normalized data. The genes are in columns. The last
#'   column contains the time variable.
#' @param de data frame. The differential expression results.
#' @param cluster.to.plot. The cluster whose genes are to be plotted. They must be included at the Cluster
#'   slot of the DE data frame.
#' @param top numeric. The number of genes with the highest qval to be plotted.
#' @param colors character. A vector of colors.
#' @param memb.cut numeric. A membership score cutoff. The genes that pass it will be highlighted in #D54C4C (red).
#' @param timeVar character. A name for the x variable of the regression model.
#' @param highlightVar character. The data of the variable that highlights the points of the plot.
#' @keywords fun9
#' @return A list with components the Dat, de, genes.to.plot and the normalized expressions with the genes
#'   to be plotted.
#'
fun9<-function(Dat,relexpr,de,cluster.to.plot,colors,top,memb.cut,timeVar,highlightVar){
    
    de<-de[match(as.numeric(as.character(de$Cluster)),cluster.to.plot,nomatch=0)>0,]
    genes.to.plot<-de$gene_short_name[sort.list(as.numeric(as.character(de$qval)))[1:top]]
    de<-de[match(genes.to.plot,de$gene_short_name),]
    if(nrow(de)<top){
        top<-nrow(de)
    }
    cols<-rep("#D6D2C4",nrow(de))
    cols[as.numeric(as.character(de$Membership))>=memb.cut]<-"#D54C4C"
    names(genes.to.plot)<-cols
    genes.to.plot<-genes.to.plot[!duplicated(genes.to.plot)]

    mm<-match(genes.to.plot,colnames(Dat),nomatch=0)
    if(sum(mm)==0){
        stop("None of these genes is found in the data")
    }
    ww<-which(mm==0)
    if(length(ww)>0){
        print(paste("Genes ",paste(genes.to.plot[mm[ww]],sep=",")," do not exist in the data",sep=""))
    }
    
    dat<-Dat[,match(c(genes.to.plot,highlightVar,timeVar),colnames(Dat))]
    relexpr<-relexpr[,match(c(genes.to.plot,highlightVar,timeVar),colnames(relexpr))]
        
 return(list(dat,de,genes.to.plot,relexpr))
}



#' Converts the plotly formatted data to ggplot formatted data.
#'
#' Converts the plotly formatted data to ggplot formatted data.
#' @param Dat data frame. A data frame with the expression data of the genes to be plotted.
#'   The genes are in columns.
#' @param genes character. The genes to be plotted. The vector names should be the colors of interest.
#' @param timeVar character. A name for the x variable of the regression model.
#' @param highlightVar character. The data of the variable that highlights the points of the plot.
#' @keywords ggdata
#' @return A data frame with data formatted for a ggplot type of plot.
#'
ggdata<-function(Dat,genes,timeVar,highlightVar){
    
    for(i in 1:(ncol(Dat)-2)){
        if(i==1){
            res<-cbind(Dat[,i],Dat[,(ncol(Dat)-1)],Dat[,ncol(Dat)],rep(colnames(Dat)[i],nrow(Dat)),rep(names(genes)[which(genes==colnames(Dat)[i])]),rep(c("Norm","Raw"),each=0.5*nrow(Dat)))
        } else {
            res<-rbind(res,cbind(Dat[,i],Dat[,(ncol(Dat)-1)],Dat[,ncol(Dat)],rep(colnames(Dat)[i],nrow(Dat)),rep(names(genes)[which(genes==colnames(Dat)[i])]),rep(c("Norm","Raw"),each=0.5*nrow(Dat))))
        }
    }
    res<-data.frame(res)
    colnames(res)<-c("Expr",highlightVar,timeVar,"Gene","Color","Type")
    res$Expr<-as.numeric(as.character(res$Expr))
    res[,which(colnames(res)==timeVar)]<-as.numeric(as.character(res[,which(colnames(res)==timeVar)]))
    res$Gene<-factor(res$Gene)
    
 return(res)
}


#' Converts the plotly formatted data to ggplot formatted data for two-branch analysis.
#'
#' Converts the plotly formatted data to ggplot formatted data for two-branch analysis.
#' @param Dat data frame. A data frame with the expression data of the genes to be plotted.
#'   The genes are in columns.
#' @param genes character. The genes to be plotted. The vector names should be the colors of interest.
#' @param ind list. Each component contains an index, internally calculated in plotBranch(), that
#'   separates the branch-specific data.
#' @param timeVar character. A name for the x variable of the regression model.
#' @param highlightVar character. The data of the variable that highlights the points of the plot.
#' @keywords ggdata
#' @return A data frame with data formatted for a ggplot type of plot.
#'
ggdata2<-function(Dat,genes,ind,timeVar,highlightVar){

    for(i in 1:(ncol(Dat)-2)){
        if(i==1){
            res<-cbind(Dat[,i],Dat[,(ncol(Dat)-1)],Dat[,ncol(Dat)],rep(colnames(Dat)[i],nrow(Dat)),rep(names(genes)[which(genes==colnames(Dat)[i])]),rep(names(ind),unlist(lapply(ind,length))),rep(c("Norm","Raw"),each=sum(unlist(lapply(ind,length)))))
        } else {
            res<-rbind(res,cbind(Dat[,i],Dat[,(ncol(Dat)-1)],Dat[,ncol(Dat)],rep(colnames(Dat)[i],nrow(Dat)),rep(names(genes)[which(genes==colnames(Dat)[i])]),rep(names(ind),unlist(lapply(ind,length))),rep(c("Norm","Raw"),each=sum(unlist(lapply(ind,length))))))
        }
    }
    res<-data.frame(res)
    colnames(res)<-c("Expr",highlightVar,timeVar,"Gene","Color","Ind","Type")
    res$Expr<-as.numeric(as.character(res$Expr))
    res[,2]<-as.numeric(as.character(res[,2]))
    res$Gene<-factor(res$Gene)

 return(res)
}


#' It plots the data with plotly or ggplot.
#'
#' It plots the data with plotly or ggplot.
#' @param dat data frame. A data frame with the normalized, scaled and smoothed expression
#'   data of the genes to be plotted. The genes are in columns.
#' @param relexpr data frame. A data frame with the normalized expressin data of the genes
#'   to be plotted. The genes are in columns.
#' @param de data frame. The differential expression results for the genes to be plotted.
#' @param genes.to.plot character. The genes to be plotted.
#' @param timeVar character. A name for the x variable of the regression model.
#' @param highlightVar character. The data of the variable that highlights the points of
#'   the plot.
#' @param id character. An id to be used in the plot's filename.
#' @param m list. The figure margins.
#' @param figure.dimensions list. The figure dimensions.
#' @param line.size numeric. The line width in ggplot.
#' @param title.size numeric. The font size for the plot title.
#' @param plot.type character. One of interactive (plotly) or static (ggplot).
#' @param screen.it logical. If TRUE, the plot will be shown on screen and not stored.
#' @param outdir character. A folder to store the plot.
#' @keywords plotFun
#' @return A plotly or a ggplot object.
#'
plotFun<-function(dat,relexpr,de,genes.to.plot,timeVar,highlightVar,id,m,figure.dimensions,line.size,title.size,plot.type,outdir,screen.it){
    
    dat<-rbind(dat,relexpr)

    if(plot.type=="interactive"){
        p1<-plotly::plot_ly(dat[1:(0.5*nrow(dat)),],
                            x=as.numeric(as.character(dat[1:(0.5*nrow(dat)),ncol(dat)])),width = figure.dimensions$width, height = figure.dimensions$height) %>%
                layout(title=list(text=paste("Modeling of Expression ~ ",timeVar,sep=""),font=list(color="black",size=title.size)),
                xaxis=list(title=paste("<b> ",timeVar," </b>",sep="")),
                yaxis=list(title=paste("<b> Normalized Expression </b>",sep=""),
                range=range(as.numeric(as.character(dat[,1:(ncol(dat)-2)])))+c(-0.5,0.5)),
                autosize = F, margin = m)

        for(i in 1:(ncol(dat)-2)){
            
            p1<-p1 %>% add_trace(y=as.numeric(as.character(dat[1:(0.5*nrow(dat)),i])),
                                 mode="lines",type="scatter",
                                 line=list(color=names(genes.to.plot)[i]),
                                 name=genes.to.plot[i],
                                 text=paste(genes.to.plot[i]," ~ ",timeVar,
                                           "<br> qval = ",round(as.numeric(as.character(de$qval[i])),3),
                                           "<br> Cluster: ",as.numeric(as.character(de$Cluster[i])),
                                           "<br> Membership: ",round(as.numeric(as.character(de$Membership[i])),3),sep=""))

            if(length(genes.to.plot)==1){
                if(length(unique(dat[,(ncol(dat)-1)]))==1){
                    p1<-p1 %>% add_trace(y=as.numeric(as.character(dat[(0.5*nrow(dat)+1):nrow(dat),i])),
                                     mode="markers",type="scatter",
                                  symbols="circle",color=I("gray60"),
                                  marker = list(size = 5))
                } else {
                    p1<-p1 %>% add_trace(y=as.numeric(as.character(dat[(0.5*nrow(dat)+1):nrow(dat),i])),
                                         mode="markers",type="scatter",
                                         symbols="circle",color=dat[(0.5*nrow(dat)+1):nrow(dat),which(colnames(dat)==highlightVar)],
                                         marker = list(size = 5))
                }
            }

        }
        
        if(!screen.it){
            # plot and message
            saveWidget(p1,paste(outdir,"simplePlotlyLines_",id,".html",sep=""),selfcontained = FALSE)
            print(paste("Download plot file: ",outdir,"simplePlotlyLines_",id,".html & plot folder: ",outdir,"simplePlotlyLines_",id,"_files",sep=""))
        } else {
        p1
        }

    } else {
        
        dat<-ggdata(Dat=dat,genes=genes.to.plot,timeVar=timeVar,highlightVar=highlightVar)
        vn<-as.character(names(genes.to.plot))
        vn<-vn[sort.list(genes.to.plot)]
    
    if(!combine){
      if(length(genes.to.plot)==1){
               if(length(unique(dat[,which(colnames(dat)==highlightVar)]))==1){
               p1<-ggplot(dat[1:(0.5*nrow(dat)),],aes(x=Pseudotime,y=Expr,color=Gene)) + geom_line(size=line.size) +
                       scale_color_manual(values=vn) +
                       xlab(timeVar) + ylab("Normalized Expression") + ggtitle(paste("Modeling of Expression ~ ",timeVar,sep="")) +
                    theme_bw() + geom_hline(yintercept=0) + theme(plot.title = element_text(size=title.size))
                p1<-p1 + geom_point(data=dat[(0.5*nrow(dat)+1):nrow(dat),],color="grey60")
            } else {
                if(is.character(dat[,which(colnames(dat)==highlightVar)]) | is.factor(dat[,which(colnames(dat)==highlightVar)])){
                    feature<-dat[(0.5*nrow(dat)+1):nrow(dat),which(colnames(dat)==highlightVar)]
                    p1<-ggplot(dat[(0.5*nrow(dat)+1):nrow(dat),],aes(x=Pseudotime,y=Expr,color=feature)) + geom_point() +
                        geom_line(data=dat[1:(0.5*nrow(dat)),],aes(x=Pseudotime,y=Expr,color=Gene),size=line.size) +
                        xlab(timeVar) + ylab("Normalized Expression") + ggtitle(paste("Modeling of Expression ~ ",timeVar,sep="")) +
                        theme_bw() + geom_hline(yintercept=0) + theme(plot.title = element_text(size=title.size))
                } else {
                    feature<-dat[(0.5*nrow(dat)+1):nrow(dat),which(colnames(dat)==highlightVar)]
                    p1<-ggplot(dat[(0.5*nrow(dat)+1):nrow(dat),],aes(x=Pseudotime,y=Expr)) + geom_point(aes(color=feature)) + scale_color_gradient(low="blue", high="red") +
                        geom_line(data=dat[1:(0.5*nrow(dat)),],aes(x=Pseudotime,y=Expr),size=line.size,color=unique(dat$Color)) +
                        xlab(timeVar) + ylab("Normalized Expression") + ggtitle(paste("Modeling of Expression ~ ",timeVar,sep="")) +
                        theme_bw() + geom_hline(yintercept=0) + theme(plot.title = element_text(size=title.size))
            }
            }
    } else {
        dat1<-dat[dat$Type=="Norm",]
        dat1$Expr<-as.numeric(as.character(dat1$Expr))
            p1<-ggplot(dat1,aes(x=Pseudotime,y=Expr,color=Gene)) + geom_line(size=line.size) +
                scale_color_manual(values=vn) +
                xlab(timeVar) + ylab("Normalized Expression") + ggtitle(paste("Modeling of Expression ~ ",timeVar,sep="")) +
                theme_bw() + geom_hline(yintercept=0) + theme(plot.title = element_text(size=title.size))
        }
        if(ncol(dat)>21){
            p1<-p1+theme(legend.position = "none")
        }
    }


    if(combine){
    if(is.null(genes.to.plot)){
       stop("Combine = TRUE is used when genes.to.plot is speficied!")
    } else {
       p1<-as.list(rep(0,length(genes.to.plot)))
       for(i in 1:length(genes.to.plot)){
        dat.comb<-dat[dat$Gene==genes.to.plot[i],]
        if(length(unique(dat.comb[,which(colnames(dat.comb)==highlightVar)]))==1){
                p1[[i]]<-ggplot(dat.comb[1:(0.5*nrow(dat.comb)),],aes(x=Pseudotime,y=Expr,color=Gene)) + geom_line(size=line.size) +
                    scale_color_manual(values=vn) +
                    xlab(timeVar) + ylab("Normalized Expression") + ggtitle(paste("Modeling of Expression ~ ",timeVar,sep="")) +
                    theme_bw() + geom_hline(yintercept=0) + theme(plot.title = element_text(size=title.size))
                p1[[i]]<-p1[[i]] + geom_point(data=dat.comb[(0.5*nrow(dat.comb)+1):nrow(dat.comb),],color="grey60")
            } else {
                if(is.character(dat.comb[,which(colnames(dat.comb)==highlightVar)]) | is.factor(dat.comb[,which(colnames(dat.comb)==highlightVar)])){
                    feature<-dat.comb[(0.5*nrow(dat.comb)+1):nrow(dat.comb),which(colnames(dat.comb)==highlightVar)]
                    p1[[i]]<-ggplot(dat.comb[(0.5*nrow(dat.comb)+1):nrow(dat.comb),],aes(x=Pseudotime,y=Expr,color=feature)) + geom_point() +
                        geom_line(data=dat.comb[1:(0.5*nrow(dat.comb)),],aes(x=Pseudotime,y=Expr,color=Gene),size=line.size) +
                        xlab(timeVar) + ylab("Normalized Expression") + ggtitle(paste("Modeling of Expression ~ ",timeVar,sep="")) +
                        theme_bw() + geom_hline(yintercept=0) + theme(plot.title = element_text(size=title.size))
                } else {
                    feature<-dat.comb[(0.5*nrow(dat.comb)+1):nrow(dat.comb),which(colnames(dat.comb)==highlightVar)]
                    p1[[i]]<-ggplot(dat[(0.5*nrow(dat.comb)+1):nrow(dat.comb),],aes(x=Pseudotime,y=Expr)) + geom_point(aes(color=feature)) + scale_color_gradient(low="blue", high="red") +
                        geom_line(data=dat.comb[1:(0.5*nrow(dat.comb)),],aes(x=Pseudotime,y=Expr),size=line.size,color=unique(dat.comb$Color)) +
                        xlab(timeVar) + ylab("Normalized Expression") + ggtitle(paste("Modeling of Expression ~ ",timeVar,sep="")) +
                        theme_bw() + geom_hline(yintercept=0) + theme(plot.title = element_text(size=title.size))
                }
            }
    }
      }
   }




        if(!screen.it){
            # plot and message
            pdf(paste(outdir,"simpleggplotLines_",id,".pdf",sep=""),width=figure.dimensions$width,height=figure.dimensions$height)
               print(p1)
            dev.off()
            print(paste("Download plot: ",outdir,"simpleggplotLines_",id,".pdf",sep=""))
        } else {
            print(p1)
        }
    }
    
 return(dat)
}


#' Generates the line plot of the selected features coming from the single-branch heatmap data.
#'
#' Generates the line plot of the selected features coming from the single-branch heatmap data.
#' @param obj list. A list with components a Monocle object with the DE genes only, a data frame
#'   of normalized and smoothed gene expression profiles, the list of DE genes with the top cluster
#'   memberships and the full memberships matrix. Typically the output of simpleHeat().
#' @param timeVar character. A name for the x variable of the regression model.
#' @param highlightVar character. The data of the variable that highlights the points of the plot. If
#'   NULL (default), the data ponts are not highlighted (colored gray).
#' @param outdir character. A folder to store the plot.
#' @param id character. An id to be used in the plot's filename.
#' @param plotting.feature character. One of gene, cluster or cluster:gene. It specifies what to plot.
#'   If gene, it expects that the names of genes to be plotted are given in parameter genes.to.plot.
#'   If cluster, it expects that the cluster number to be plotted is given in parameter cluster.to.plot.
#'   If cluster:gene, it expects that the cluster number to be plotted is given in parameter cluster.to.plot
#'   as well as an entry in one of top.hits, qval.cut or membership.cut parameters. It will plot those genes
#'   of the provided cluster which satisfy one of the above cutoffs / criteria. Default is cluster.
#' @param cluster.to.plot numeric. The cluster number to plot. Default is 1.
#' @param color.by character. One of gene, cluster and membership (default). If gene, it will give a different
#    color to each plotted gene. It is only activated when plotting.feature is gene or cluster:gene. If cluster,
#'   it will give a different color to each cluster. It is only activated when plotting.feature is gene. If
#'   membership the plotted genes will be colored using the different membership scores. Default is membership.
#' @param genes.to.plot character. The genes to be plotted. They must be included at the short_gene_name.
#'   Default is NULL.
#'   slot of the DE data frame.
#' @param colors character. A vector of colors. If color.by is cluster, the vector names should be the cluster numbers.
#'   Default is NULL where 432 colors will be automatically generated.
#' @param plot.type character. One of interactive (plotly) or static (ggplot). Default is interactive.
#' @param top.hits numeric. The number of genes with the x highest qval or the highest memberships to be plotted (depending
#'   on the values of plotting.feature and color.by parameters). Default is 10.
#' @param qval.cut numeric. The qvalue cutoff to select the top genes to be plotted. Default is 0.01.
#' @param membership.cut numeric. The membership cutoff to select the top genes to be plotted. Default is 0.5.
#' @param figure.margins list. The plotly figure margins. Default is list(l = 100,r = 100,b = 100,t = 100,pad = 8).
#' @param figure.dimensions list. The plotly and pdf (ggplot) figure dimensions. Default is list(width=700,height=700).
#' @param line.size numeric. The line width in ggplot. Default is 1
#' @param title.size numeric. The font size for the plot title. Default is 20.
#' @param screen.it logical. If TRUE, the plot will be shown on screen and not stored. Default is FALSE.
#' @param seed numeric. A seed number to randomize the vector of automatically generated colors if color.palette = NULL.
#'   Default is 999.
#' @keywords plotSimple
#' @return A list with components the Monocle object, the Dat, the de and the genes.to.plot with the plotted genes.
#'
plotSimple<-function(obj,timeVar,highlightVar,outdir,id,
                     plotting.feature="cluster",
                     cluster.to.plot=1,
                     color.by="membership",
                     genes.to.plot=NULL,
                     color.palette=NULL,
                     plot.type="interactive",
                     top.hits=10,
                     qval.cut=0.01,
                     membership.cut=0.5,
                     figure.margins=list(l = 100,r = 100,b = 100,t = 100,pad = 8),
                     figure.dimensions=list(width=700,height=700),
                     title.size=20,
             line.size=1,
             screen.it=FALSE,
                     seed=999){

    # generate the color palette if it is NULL
    if(is.null(color.palette)){
        set.seed(seed)
        gencol = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
        gencol<-gencol[-1]
        gencol<-sample(gencol)
        names(gencol)<-1:length(gencol)
    color.palette<-gencol
    }
    
    # make the baseline data
    if(length(which(names(obj)=="DE"))==0 | length(which(names(obj)=="NormData"))==0){
        stop("obj does not contain the NormData and / or the DE slots. Run simpleHeat() first!")
    }

    relexpr<-t(t(obj$Obj@assayData$exprs) /  pData(obj$Obj)[, 'Size_Factor'])
    relexpr<-log(relexpr+1,2)
    relexpr<-t(scale(t(relexpr),scale=T))

    if(is.null(highlightVar)){
    hv<-rep(1,nrow(pData(obj$Obj)))
    } else {
        hv<-which(colnames(pData(obj$Obj))==highlightVar)
    if(length(hv)==0){
        stop("The highlightVar is not present in the pData(obj)!")
    } else {
        hv<-pData(obj$Obj)[,hv]
    }
    }

    DE<-obj$DE
    Dat<-t(obj$NormData)
    Dat<-cbind(Dat,hv,pData(obj$Obj)$Pseudotime)
    Dat<-data.frame(Dat,check.names=F)
    colnames(Dat)<-c(GeneNames(colnames(Dat)[1:(ncol(Dat)-2)]),highlightVar,timeVar)

    relexpr<-t(relexpr[match(rownames(obj$NormData),rownames(relexpr)),])
    relexpr<-cbind(relexpr,hv,pData(obj$Obj)$Pseudotime)
    relexpr<-data.frame(relexpr,check.names=F)
    colnames(relexpr)<-c(GeneNames(colnames(relexpr)[1:(ncol(relexpr)-2)]),highlightVar,timeVar)

    obj<-obj$Obj[match(DE[,1],rownames(obj$Obj)),]
     
       
    marg<-match.arg(plotting.feature,c("gene","cluster","cluster:gene"))

    # 1. when selecting genes and color based on gene names
    if(plotting.feature=="gene" & color.by=="gene"){
        if(is.null(genes.to.plot[1]) | is.na(genes.to.plot[1]) | length(genes.to.plot)==0){
            stop("Specify the genes to be plotted!")
        }
        if(length(genes.to.plot)>432){
            stop("The maximum number of genes.to.plot is 432! To plot more, use a custom color palette with more colors.")
        }
        if(length(genes.to.plot)>length(color.palette)){
            stop("The number of genes are more than the number of colors. Use more colors in color.palette!")
        }
        print("Plotting selected genes and coloring by gene.")
        dat<-fun1(Dat=Dat,relexpr=relexpr,de=DE,genes.to.plot=genes.to.plot,colors=gencol,timeVar=timeVar,highlightVar=highlightVar)
        if(nrow(dat[[1]])==0){
            stop("No data to plot.")
        }
        p1<-plotFun(dat=dat[[1]],de=dat[[2]],genes.to.plot=dat[[3]],relexpr=dat[[4]],timeVar=timeVar,highlightVar=highlightVar,id=id,outdir=outdir,
            m=figure.margins,figure.dimensions=figure.dimensions,line.size=line.size,title.size=title.size,plot.type=plot.type,screen.it=screen.it)
    }

    # 2. when selecting genes and color based on memberships
    if(plotting.feature=="gene" & color.by=="membership"){
        if(is.null(genes.to.plot[1]) | is.na(genes.to.plot[1]) | length(genes.to.plot)==0){
            stop("Specify the genes to be plotted!")
        }
        if(length(genes.to.plot)>432){
            stop("The maximum number of genes.to.plot is 432! To plot more, use a custom color palette with more colors.")
        }
        if(length(genes.to.plot)>length(color.palette)){
            stop("The number of genes are more than the number of colors. Use more colors in color.palette!")
        }
        if(membership.cut<0 | membership.cut>1){
            stop("membership.cut must be a number in [0,1].")
        }
        print("Plotting selected genes and coloring by membership cutoff.")
        dat<-fun2(Dat=Dat,relexpr=relexpr,de=DE,genes.to.plot=genes.to.plot,memb.cut=membership.cut,timeVar=timeVar,highlightVar=highlightVar)
        if(nrow(dat[[1]])==0){
            stop("No data to plot.")
        }
        p1<-plotFun(dat=dat[[1]],de=dat[[2]],genes.to.plot=dat[[3]],relexpr=dat[[4]],timeVar=timeVar,highlightVar=highlightVar,id=id,outdir=outdir,
            m=figure.margins,figure.dimensions=figure.dimensions,line.size=line.size,title.size=title.size,plot.type=plot.type,screen.it=screen.it)
    }

    # 3. when selecting genes and color based on cluster
    if(plotting.feature=="gene" & color.by=="cluster"){
        if(is.null(genes.to.plot[1]) | is.na(genes.to.plot[1]) | length(genes.to.plot)==0){
            stop("Specify the genes to be plotted!")
        }
        if(length(genes.to.plot)>432){
            stop("The maximum number of genes.to.plot is 432! To plot more, use a custom color palette with more colors.")
        }
        if(length(genes.to.plot)>length(color.palette)){
            stop("The number of genes are more than the number of colors. Use more colors in color.palette!")
        }
        print("Plotting selected genes and coloring by cluster.")
        dat<-fun3(Dat=Dat,relexpr=relexpr,de=DE,genes.to.plot=genes.to.plot,colors=gencol,timeVar=timeVar,highlightVar=highlightVar)
        if(nrow(dat[[1]])==0){
            stop("No data to plot.")
        }
        p1<-plotFun(dat=dat[[1]],de=dat[[2]],genes.to.plot=dat[[3]],relexpr=dat[[4]],timeVar=timeVar,highlightVar=highlightVar,id=id,outdir=outdir,
            m=figure.margins,figure.dimensions=figure.dimensions,line.size=line.size,title.size=title.size,plot.type=plot.type,screen.it=screen.it)
    }

    # 4. when selecting cluster and color based on membership cutoff
    if(plotting.feature=="cluster" & color.by=="membership"){
        
        if(is.null(cluster.to.plot[1]) | is.na(cluster.to.plot[1])){
            stop("Specify the cluster to be plotted!")
        }
        if(length(cluster.to.plot)>1){
            print(paste("Only one cluster can be plotted at each time. Only cluster ",cluster.to.plot[1]," will be plotted.",sep=""))
        }
        if(match(cluster.to.plot,as.numeric(as.character(DE$Cluster)),nomatch=0)==0){
            stop("This cluster does not exist in the data.")
        }
        if(membership.cut<0 | membership.cut>1){
            stop("membership.cut must be a number in [0,1].")
        }
        if(color.by!="membership"){
            color.by<-"membership"
            print("The value of parameter color.by is changed to membership")
        }
        print("Plotting selected cluster and coloring by membership cutoff. To color by top membership scores set membership.cut = NULL.")
        dat<-fun4(Dat=Dat,relexpr=relexpr,de=DE,cluster.to.plot=cluster.to.plot,memb.cut=membership.cut,timeVar=timeVar,highlightVar=highlightVar)
        if(nrow(dat[[1]])==0){
            stop("No data to plot.")
        }
        p1<-plotFun(dat=dat[[1]],de=dat[[2]],genes.to.plot=dat[[3]],relexpr=dat[[4]],timeVar=timeVar,highlightVar=highlightVar,id=id,outdir=outdir,
            m=figure.margins,figure.dimensions=figure.dimensions,line.size=line.size,title.size=title.size,plot.type=plot.type,screen.it=screen.it)
    }

    # 5. when selecting cluster and color based on membership top
    if(plotting.feature=="cluster" & color.by=="membership" & is.null(membership.cut)){
        
        if(match(cluster.to.plot,as.numeric(as.character(DE$Cluster)),nomatch=0)==0){
            stop("This cluster does not exist in the data.")
        }
        if(is.null(top[1]) | is.na(top[1]) | top==0){
            stop("Specify top parameter. It should be a positive integer.")
        }
        if(color.by!="membership"){
            color.by<-"membership"
            print("The value of parameter color.by is changed to membership")
        }
        print("Plotting selected cluster and coloring by top membership scores.")
        dat<-fun5(Dat=Dat,relexpr=relexpr,de=DE,cluster.to.plot=cluster.to.plot,top=top.hits,timeVar=timeVar,highlightVar=highlightVar)
        if(nrow(dat[[1]])==0){
            stop("No data to plot.")
        }
        p1<-plotFun(dat=dat[[1]],de=dat[[2]],genes.to.plot=dat[[3]],relexpr=dat[[4]],timeVar=timeVar,highlightVar=highlightVar,id=id,outdir=outdir,
            m=figure.margins,figure.dimensions=figure.dimensions,line.size=line.size,title.size=title.size,plot.type=plot.type,screen.it=screen.it)
    }

    # 6. when selecting cluster and top x genes (by top qval) and color based on gene
    if(plotting.feature=="cluster:gene" & color.by=="gene" & is.null(qval.cut)){
        if(match(cluster.to.plot,as.numeric(as.character(DE$Cluster)),nomatch=0)==0){
            stop("This cluster does not exist in the data.")
        }
        if(is.null(top[1]) | is.na(top[1]) | top==0){
            stop("Specify top parameter. It should be a positive integer.")
        }
        print("Plotting selected cluster-specific genes and coloring by gene. Only the top genes in terms of qval are used.")
        dat<-fun6(Dat=Dat,relexpr=relexpr,de=DE,cluster.to.plot=cluster.to.plot,colors=gencol,top=top.hits,timeVar=timeVar,highlightVar=highlightVar)
        if(nrow(dat[[1]])==0){
            stop("No data to plot.")
        }
        p1<-plotFun(dat=dat[[1]],de=dat[[2]],genes.to.plot=dat[[3]],relexpr=dat[[4]],timeVar=timeVar,highlightVar=highlightVar,id=id,outdir=outdir,
            m=figure.margins,figure.dimensions=figure.dimensions,line.size=line.size,title.size=title.size,plot.type=plot.type,screen.it=screen.it)
    }

    # 7. when selecting cluster and top x genes (by cutoff qval) and color based on gene
    if(plotting.feature=="cluster:gene" & color.by=="gene"){
        if(match(cluster.to.plot,as.numeric(as.character(DE$Cluster)),nomatch=0)==0){
            stop("This cluster does not exist in the data.")
        }
        if(qval.cut<=0 | qval.cut>1){
            stop("qval.cut must be a number in (0,1].")
        }
        print("Plotting selected cluster-specific genes and coloring by gene. Only the genes that pass the qval citerion are used. To use the top criterion set qval.cut = NULL.")
        dat<-fun7(Dat=Dat,relexpr=relexpr,de=DE,cluster.to.plot=cluster.to.plot,colors=gencol,qcut=qval.cut,timeVar=timeVar,highlightVar=highlightVar)
        if(nrow(dat[[1]])==0){
            stop("No data to plot.")
        }
        p1<-plotFun(dat=dat[[1]],de=dat[[2]],genes.to.plot=dat[[3]],relexpr=dat[[4]],timeVar=timeVar,highlightVar=highlightVar,id=id,outdir=outdir,
            m=figure.margins,figure.dimensions=figure.dimensions,line.size=line.size,title.size=title.size,plot.type=plot.type,screen.it=screen.it)
    }

    # 8. when selecting cluster and top x genes (by cutoff qval) and color based on membership
    if(plotting.feature=="cluster:gene" & color.by=="membership"){
        if(match(cluster.to.plot,as.numeric(as.character(DE$Cluster)),nomatch=0)==0){
            stop("This cluster does not exist in the data.")
        }
        if(qval.cut<=0 | qval.cut>1){
            stop("qval.cut must be a number in (0,1].")
        }
        if(membership.cut<0 | membership.cut>1){
            stop("membership.cut must be a number in [0,1].")
        }
        print("Plotting selected cluster-specific genes and coloring by membership. Only the genes that pass the qval citerion are used. To use the top criterion set qval.cut = NULL.")
        dat<-fun8(Dat=Dat,relexpr=relexpr,de=DE,cluster.to.plot=cluster.to.plot,colors=gencol,qcut=qval.cut,memb.cut=membership.cut,timeVar=timeVar,highlightVar=highlightVar)
        if(nrow(dat[[1]])==0){
            stop("No data to plot.")
        }
        p1<-plotFun(dat=dat[[1]],de=dat[[2]],genes.to.plot=dat[[3]],relexpr=dat[[4]],timeVar=timeVar,highlightVar=highlightVar,id=id,outdir=outdir,
            m=figure.margins,figure.dimensions=figure.dimensions,line.size=line.size,title.size=title.size,plot.type=plot.type,screen.it=screen.it)
    }


    # 9. when selecting cluster and top x genes (by top qval) and color based on membership
    if(plotting.feature=="cluster:gene" & color.by=="gene" & is.null(qval.cut)){
        if(match(cluster.to.plot,as.numeric(as.character(DE$Cluster)),nomatch=0)==0){
            stop("This cluster does not exist in the data.")
        }
        if(membership.cut<0 | membership.cut>1){
            stop("membership.cut must be a number in [0,1].")
        }
        if(is.null(top[1]) | is.na(top[1]) | top==0){
            stop("Specify top parameter. It should be a positive integer.")
        }
        print("Plotting selected cluster-specific genes and coloring by membership. Only the top genes in terms of qval are used.")
        dat<-fun9(Dat=Dat,relexpr=relexpr,de=DE,cluster.to.plot=cluster.to.plot,colors=gencol,top=top.hits,memb.cut=membership.cut,timeVar=timeVar,highlightVar=highlightVar)
        if(nrow(dat[[1]])==0){
            stop("No data to plot.")
        }
        p1<-plotFun(dat=dat[[1]],de=dat[[2]],genes.to.plot=dat[[3]],relexpr=dat[[4]],timeVar=timeVar,highlightVar=highlightVar,id=id,outdir=outdir,
            m=figure.margins,figure.dimensions=figure.dimensions,line.size=line.size,title.size=title.size,plot.type=plot.type,screen.it=screen.it)
    }
    
 return(list(Obj=obj,Data=dat))
}
   



#' Extracts the requested component from the rownames (gene IDs) of the Seurat object.
#'
#' Extracts the requested component from the rownames (gene IDs) of the Seurat object. The gene IDs of the Seurat object
#'   have been formatted as EnsemblID:GeneName. The function will extract either the EnsemblID or the GeneName.
#' @param data character. A vector of gene IDs.
#' @param symbol character. The symbol that separates the EnsemblID from the GeneName. Defaut is ':'.
#' @param what character. If Names (default), the gene names will be extracted, otherwise the Ensembl IDs will be extracted.
#' @param remove.dot logical. If TRUE (default) and what != Names, the dot and what follows it from the Ensembl IDs will be
#'   removed.
#' @keywords GeneNames
#' @return A vector with the extracted gene names or Ensembl IDs.
#'
GeneNames<-function(data,symbol=":",what="Names",remove.dot=TRUE){
    ss<-strsplit(as.character(data), symbol)
    k<-ifelse(what=="Names",2,1)
    res<-rep(0,length(ss))
    for(i in 1:length(ss)){
        res[i]<-ss[[i]][min(k,length(ss[[i]]))]
        if(k==1 & remove.dot==TRUE){
            res[i]<-unlist(strsplit(as.character(res[i]),".",fixed=TRUE))[1]
        }
    }
    return(res)
}





#' A helper funtion that converts the GSEA downloaded .gmt pathway files in the
#'   format accepted for my_fgsea() function.
#'
#' A helper funtion that converts the GSEA downloaded .gmt pathway files in the
#'   format accepted for my_fgsea() function.
#' @param gmtFie character. The filename of the .gmt file to be processed.
#' @param type character. One of BP, CC, MF, KEGG or REACTOME.
#' @keywords readGSEAdata
#' @return The pathways in the my_fgsea() format.
#'
readGSEAdata<-function(gmtFile,type){
    
    x<-scan(gmtFile,what="character")
    pat<-NULL
    if(type=="BP"){
        pat<-"GOBP_"
    }
    if(type=="MF"){
        pat<-"GOMF_"
    }
    if(type=="CC"){
        pat<-"GOCC_"
    }
    if(type=="KEGG"){
        pat<-"KEGG_"
    }
    if(type=="REACTOME"){
        pat<-"REACTOME_"
    }
    if(is.null(pat)){
        stop("Parameter type should be one of BP, MF, CC, KEGG or REACTOME")
    }
    
    ww<-grep(pat,x)
    if(length(ww)==0){
        stop("The gmtFile does not contain data of that type!")
    }
    out<-ww[seq(2,length(ww),2)]
    x<-x[-out]
    ww<-grep(pat,x)
    ww<-c(ww,length(x))
    res<-as.list(rep(0,(length(ww)-1)))
    names(res)<-x[ww[-length(ww)]]
    for(i in 1:(length(ww)-1)){
        rr<-x[(ww[i]+1):(ww[(i+1)]-1)]
        res[[i]]<-rr
    }
 return(res)
    
}



#' It runs the fgsea analysis.
#'
#' It runs the fgsea analysis for a predifined set of pathways.
#' @param de data frame. A data frame of DE results where the different clusters and qval
#'   must be present.
#' @param pathwayFile character. The filename with the pathways of interest in the form of .gmt.
#' @param type character. The type of pathways to be interrogated. It must be one of BP (default),
#'   CC, MF, KEGG or REACTOME. Only pathways downloaded from https://www.gsea-msigdb.org/gsea/downloads.jsp
#'   are accepted.
#' @param cluster numeric. The DE cluster whose genes will be interrogated. Default is 1.
#' @param size.range numeric. Only the pathways whose sample size are within the specified range
#'   will be analyzed.
#' @param qcut numeric. Only the DE genes of the predefined cluster will be analyzed. In Monocle analysis
#'   logFC is not meaningfull (we study cluster patterns), so the up- and down-regulated gene input that
#'   fgsea typically expected cannot be used. For this reason we utilize the -log10 qvalue of the estimation
#'   and we estiate enrichment from these genes. Default is 0.05.
#' @keywords my_fgsea
#' @return A data frame with the fgsea results.
#'
my_fgsea<-function(de,pathwayFile,type="BP",cluster=1,size.range=c(10,500),qcut=0.05){
    
    pathways<-readGSEAdata(gmtFile=pathwayFile,type=type)
    de<-de[de$Cluster==cluster,]
    de<-de[de$qval<=qcut,]
    nn<-GeneNames(de$Gene)
    pp<- -log(de$qval,10)
    
    tt<-table(nn)
    tt<-which(tt>1)
    if(length(tt)>0){
         out<-c()
         for(i in 1:length(tt)){
             ww<-which(nn==names(tt)[i])
             pp1<- pp[ww]
             a<-which.min(pp1)
             out<-c(out,ww[a])
         }
         nn<-nn[-out,]
         pp<-pp[-out]
    }
    stats<-pp
    names(stats)<-nn
     
    fgseaRes <- fgsea(pathways = pathways,
                      stats    = stats,
                      minSize  = size.range[1],
                      maxSize  = size.range[2])
    fgsea<-cbind(fgsea,Type=rep(type,nrow(fgsea)),Cluster=rep(cluster,nrow(fgsea)))
    
 return(fgsea)
    
}




#' Generates the data for GSEA analysis using the correlation (increasing profile) method.
#'
#' Generates the data for GSEA analysis using the correlation (increasing profile) method.
#' @param obj list. A list with components a Monocle object with the DE genes only, a data frame of normalized and smoothed
#'   gene expression profiles (all data and DE only), the list of DE genes with the top cluster memberships and the
#'   full memberships matrix. Typically the output of simpleHeat().
#' @param id character. An id to be used in the GSEA data filenames.
#' @param outdir character. A folder to store the GSEA data.
#' @param qval.cut numeric. The qvalue cutoff to select the top genes to be plotted. Default is 0.01.
#' @param membership.cut numeric. The membership cutoff to select the top genes to be plotted. Default is 0.5.
#' @keywords GSEAdataset
#' @return A list with the normalized data of all genes and the cluser-specific cls files. The data are also
#'   stored in the specified folder.
#'
GSEAdataset<-function(obj,id,outdir,qval.cut=0.01,membership.cut=0.5){
    
    norm<-obj$NormDataAll
    #norm<-cbind(norm$A,norm$B)
    de<-obj$DE
    de<-de[de$qval<=qval.cut,]
    cl<-de$Cluster
    uu<-sort(unique(cl))
    GSEAdata<-as.list(rep(0,length(uu)))
    
    if(!dir.exists(outdir)){
        dir.create(outdir)
        print(paste("Folder ",outdir," has been created to store the GSEAdataset() output files.",sep=""))
    }
    
    for(i in 1:length(uu)){
        
        dem<-de[de$Membership>=membership.cut & de$Cluster==uu[i],]
        if(nrow(dem)<5){
            dem<-de[sort.list(de$Membership,decreasing=T)[1:5],]
        }
        mm<-match(dem[,1],rownames(norm))
        normm<-norm[mm,]
        ts<-apply(normm,2,mean)
        GSEAdata[[i]]<-matrix(c("#numeric","#IncreasingProfile",paste(round(ts,3),collapse=" ")),nrow=3)
        write.table(GSEAdata[[i]],paste(outdir,"GSEA_cluster",i,"_",id,".cls",sep=""),sep="\t",row.names=F,col.names=F,quote=F)
    }
    ss<-paste("Sample",1:ncol(norm),sep="")
    norm<-cbind(GeneNames(rownames(norm)),rownames(norm),norm)
    colnames(norm)<-c("NAME","DESCRIPTION",ss)
    write.table(norm,paste(outdir,"GSEA_data_",id,".txt",sep=""),row.names=F,quote=F)
    GSEAdata<-list(Data=norm,cls=GSEAdata)

    # message
    print(paste("Download data: ",outdir,"GSEA_cluster* and ",outdir,"GSEA_data*.",sep=""))


 return(GSEAdata)
}



#' Performs the clusterProfiler enrichGO analysis for each cluster of DE genes.
#'
#' Performs the clusterProfiler enrichGO analysis for each cluster of DE genes.
#' @param obj object. A list with components a Monocle object with the DE genes only, a data frame of normalized and smoothed
#'   gene expression profiles (all data and DE only), the list of DE genes with the top cluster memberships and the
#'   full memberships matrix. Typically the output of simpleHeat(). Only the DE slot is required.
#' @param id character. An id to be used in the clusterProfiler plot filename.
#' @param outdir character. A folder to store the clusterProfiler plot.
#' @param cluster numeric. The cluster(s) whose genes will be tested in the enrichment analysis model. Default is 1.
#' @param organism character. One of human (default) or mouse.
#' @param keyType character. One of SYMBOL (default) or ENSEMBL. It specifies the gene identifiers.
#' @param ontology character. one of "BP" (default), "CC", "MF" or "ALL" for all three gene ontologies.
#' @param pvalue.cut numeric. The adjusted pvalue cutoff on enrichment tests to report. Default is 0.05.
#' @param qvalue.cut numeric. The qvalue cutoff on enrichment tests to report as significant. Tests must pass i) pvalueCutoff on unadjusted
#'   pvalues, ii) pvalueCutoff on adjusted pvalues and iii) qvalueCutoff on qvalues to be reported.
#' @param adjMethod character. One of "holm", "hochberg", "hommel", "bonferroni", "BH" (default) , "BY", "fdr" or "none" for pvalue adjustment.
#' @param universe character. The background genes. If NULL (default), the all genes listed in the database will be used as background.
#' @param range.size numeric. The minimal and maximal size of genes annotated by Ontology term for testing. Default is c(10,500).
#' @param gene.freq numeric. A cutoff that specifies the minimum number of genes a significant GO term should contain in order to
#'   be considered for appearance in the plot (see also showCategory parameter for the final list of terms to be displayed). Default
#'   is 5. If the significant terms are less than 20, the cutoff automatically drops to 2.
#' @param showCategory numeric. The maximum number of enriched terms to appear in the plot. Default is 50.
#' @param cex_gene numeric. A number indicating the amount by which plotting gene nodes should be scaled relative to the default. Default
#'   is 1.
#' @param cex_label_gene numeric. The scale of gene node label size. Default is 1.
#' @param cex_category numeric. A number indicating the amount by which plotting category nodes should be scaled relative to the default.
#'   Default is 1.
#' @param cex_label_category numeric. The scale of category node label size. Default is 1.
#' @keywords clusterProfilerFun
#' @return A data frame with the enrichGO results for the genes of the selected cluster(s).
#'
clusterProfilerFun<-function(obj,id,outdir,cluster=1,
                            organism="human",keyType="SYMBOL",ontology="BP",
                            pvalue.cut=0.05,qvalue.cut=0.2,adjMethod="BH",universe=NULL,
                size.range=c(10,500),gene.freq=5,showCategory=50,
                            cex_gene=1,cex_label_gene=1,cex_category=1,cex_label_category=1){
                                
    keyType<-match.arg(keyType,c("SYMBOL","ENSEMBL"))
    de<-obj$DE[match(obj$DE$Cluster,cluster,nomatch=0)>0,]
    if(keyType=="SYMBOL"){
        de[,1]<-GeneNames(de[,1])
    }
    if(keyType=="ENSEMBL"){
        de[,1]<-GeneNames(de[,1],what="Ensg")
    }
        
    genes<-as.numeric(as.character(de$qval))
    names(genes)<-de[,1]
    genes<-sort(genes)
    if(organism=="human"){
        ego <- enrichGO(gene     = names(genes),
                    OrgDb        = org.Hs.eg.db,
                    ont          = ontology,
                    pvalueCutoff = pvalue.cut,
            qvalueCutoff = qvalue.cut,
            pAdjustMethod = adjMethod,
            universe=universe,
                    minGSSize = size.range[1],
                    maxGSSize = size.range[2],
                    keyType = keyType)
    }
        
    if(organism=="mouse"){
        ego <- enrichGO(gene     = names(genes),
                    OrgDb        = org.Mm.eg.db,
                    ont          = ontology,
                    pvalueCutoff = pvalue.cut,
            qvalueCutoff = qvalue.cut,
             pAdjustMethod = adjMethod,
            universe=universe,
                    minGSSize = size.range[1],
                    maxGSSize = size.range[2],
                    keyType = keyType)
    }

    ego <- simplify(ego)
    egop<-ego
    freq<-ifelse(nrow(egop@result)>20,gene.freq,2)
    egop@result<-egop@result[egop@result$Count>freq,]

    pdf(paste(outdir,"clusterProfiler_",id,"_",ontology,".pdf",sep=""))
        print(enrichplot::cnetplot(egop,showCategory = showCategory,foldChange= -log(genes,10),
                   cex_gene=cex_gene,cex_category=cex_category,
                   cex_label_category=cex_label_category,cex_label_gene=cex_label_gene))
    dev.off()

    # message
    print(paste("Download plot: ", outdir,"clusterProfiler_",id,"_",ontology,".pdf",sep=""))
        
    ego<-cbind(ego@result,Cluster=paste(cluster,collapse="&"))
    
 return(ego)
}


#' Generates the line plot of the selected features coming from the two-branch heatmap data.
#'
#' Generates the line plot of the selected features coming from the two-branch heatmap data.
#' @param obj list. A list with components a Monocle object with the DE genes only, a data frame
#'   of normalized and smoothed gene expression profiles, the list of DE genes with the top cluster
#'   memberships and the full memberships matrix. Typically the output of branchHeat().
#' @param timeVar character. A name for the x variable of the regression model.
#' @param outdir character. A folder to store the plot.
#' @param id character. An id to be used in the plot's filename.
#' @param plotting.feature character. One of gene, cluster or cluster:gene. It specifies what to plot.
#'   If gene, it expects that the names of genes to be plotted are given in parameter genes.to.plot.
#'   If cluster, it expects that the cluster number to be plotted is given in parameter cluster.to.plot.
#'   If cluster:gene, it expects that the cluster number to be plotted is given in parameter cluster.to.plot
#'   as well as an entry in one of top.hits, qval.cut or membership.cut parameters. It will plot those genes
#'   of the provided cluster which satisfy one of the above cutoffs / criteria. Default is cluster.
#' @param cluster.to.plot numeric. The cluster number to plot. Default is 1.
#' @param color.by character. One of gene, cluster and membership (default). If gene, it will give a different
#    color to each plotted gene. It is only activated when plotting.feature is gene or cluster:gene. If cluster,
#'   it will give a different color to each cluster. It is only activated when plotting.feature is gene. If
#'   membership the plotted genes will be colored using the different membership scores. Default is membership.
#' @param genes.to.plot character. The genes to be plotted. They must be included at the short_gene_name.
#'   Default is NULL.
#'   slot of the DE data frame.
#' @param colors character. A vector of colors. If color.by is cluster, the vector names should be the cluster numbers.
#'   Default is NULL where 432 colors will be automatically generated.
#' @param plot.type character. One of interactive (plotly) or static (ggplot). Default is interactive.
#' @param top.hits numeric. The number of genes with the x highest qval or the highest memberships to be plotted (depending
#'   on the values of plotting.feature and color.by parameters). Default is 10.
#' @param qval.cut numeric. The qvalue cutoff to select the top genes to be plotted. Default is 0.01.
#' @param membership.cut numeric. The membership cutoff to select the top genes to be plotted. Default is 0.5.
#' @param figure.margins list. The plotly figure margins. Default is list(l = 100,r = 100,b = 100,t = 100,pad = 8).
#' @param figure.dimensions list. The plotly and pdf (ggplot) figure dimensions. Default is list(width=700,height=700).
#' @param line.size numeric. The line width in ggplot. Default is 1
#' @param title.size numeric. The font size for the plot title. Default is 20.
#' @param seed numeric. A seed number to randomize the vector of automatically generated colors if color.palette = NULL.
#'   Default is 999.
#' @param screen.it logical. If TRUE, the plot will be shown on screen and not stored. Default is FALSE.
#' @keywords plotBranch
#' @return A list with components the Monocle object, the Dat, the de and the genes.to.plot with the plotted genes.
#'
plotBranch<-function(obj,timeVar,outdir,id,
                     plotting.feature="cluster",
                     cluster.to.plot=1,
                     color.by="membership",
                     genes.to.plot=NULL,
                     color.palette=NULL,
                     plot.type="interactive",
                     top.hits=10,
                     qval.cut=0.01,
                     membership.cut=0.5,
                     figure.margins=list(l = 100,r = 100,b = 100,t = 100,pad = 8),
                     figure.dimensions=list(width=700,height=700),
                     title.size=20,
                     line.size=1,
                     seed=999,
             screen.it=FALSE){

    # generate the color palette if it is NULL
    if(is.null(color.palette)){
        set.seed(seed)
        gencol = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
        gencol<-gencol[-1]
        gencol<-sample(gencol)
        names(gencol)<-1:length(gencol)
    }
    
    # make the baseline data
    if(length(which(names(obj)=="DE"))==0 | length(which(names(obj)=="NormData"))==0){
        stop("obj does not contain the NormData and / or the DE slots. Run branchHeat() first!")
    }


    relexpr<-t(t(exprs(obj$Obj)) /  pData(obj$Obj)[, 'Size_Factor'])
    relexpr<-log(relexpr+1,2)
    relexpr<-t(scale(t(relexpr),scale=T))

    DE<-obj$DE
    Dat<-t(obj$NormData)
    pp<-pData(obj$Obj)
    uu<-sort(unique(pp$Branch))
    pseudo<-ind<-as.list(rep(0,length(uu)))
    start<-1
    nn<-c()
    for(i in 1:length(uu)){
        pseudo[[i]]<-pp$Pseudotime[pp$Branch==uu[i]]
        end<-(start+length(pseudo[[i]])-1)
        ind[[i]]<-start:end
        start<-(end+1)
        nn<-c(nn,unique(pp$State[pp$Branch==uu[i]]))
    }
    names(ind)<-nn
    pseudo<-unlist(pseudo)
    Dat<-cbind(Dat,as.numeric(as.character(pseudo)))
    Dat<-data.frame(Dat,check.names=F)
    colnames(Dat)<-c(GeneNames(colnames(Dat)[-ncol(Dat)]),timeVar)
    Dat[,ncol(Dat)]<-as.numeric(as.character(Dat[,ncol(Dat)]))

    relexpr<-t(relexpr[match(rownames(obj$NormData),GeneNames(rownames(relexpr))),])
    relexpr<-cbind(relexpr,pData(obj$Obj)$Pseudotime)
    relexpr<-data.frame(relexpr,check.names=F)
    colnames(relexpr)<-c(GeneNames(colnames(relexpr)[-ncol(relexpr)]),timeVar)

    obj<-obj$Obj[match(DE[,1],rownames(obj$Obj)),]

    marg<-match.arg(plotting.feature,c("gene","cluster","cluster:gene"))

    # 1. when selecting genes and color based on gene names
    if(plotting.feature=="gene" & color.by=="gene"){
        if(is.null(genes.to.plot[1]) | is.na(genes.to.plot[1]) | genes.to.plot==0){
            stop("Specify the genes to be plotted!")
        }
        if(length(genes.to.plot)>432){
            stop("The maximum number of genes.to.plot is 432! To plot more, use a custom color palette with more colors.")
        }
        if(length(genes.to.plot)>length(color.palette)){
            stop("The number of genes are more than the number of colors. Use more colors in color.palette!")
        }
        print("Plotting selected genes and coloring by gene.")
        dat<-fun1(Dat=Dat,relexpr=relexpr,de=DE,genes.to.plot=genes.to.plot,colors=gencol,timeVar=timeVar)
        if(nrow(dat[[1]])==0){
            stop("No data to plot.")
        }
        p1<-plotFun2(dat=dat[[1]],de=dat[[2]],genes.to.plot=dat[[3]],relexpr=dat[[4]],timeVar=timeVar,id=id,
             ind=ind,m=figure.margins,figure.dimensions=figure.dimensions,line.size=line.size,title.size=title.size,plot.type=plot.type,screen.it=screen.it)
    }

    # 2. when selecting genes and color based on memberships
    if(plotting.feature=="gene" & color.by=="membership"){
        if(is.null(genes.to.plot[1]) | is.na(genes.to.plot[1]) | genes.to.plot==0){
            stop("Specify the genes to be plotted!")
        }
        if(length(genes.to.plot)>432){
            stop("The maximum number of genes.to.plot is 432! To plot more, use a custom color palette with more colors.")
        }
        if(length(genes.to.plot)>length(color.palette)){
            stop("The number of genes are more than the number of colors. Use more colors in color.palette!")
        }
        if(membership.cut<0 | membership.cut>1){
            stop("membership.cut must be a number in [0,1].")
        }
        print("Plotting selected genes and coloring by membership cutoff.")
        dat<-fun2(Dat=Dat,relexpr=relexpr,de=DE,genes.to.plot=genes.to.plot,memb.cut=membership.cut,timeVar=timeVar)
        if(nrow(dat[[1]])==0){
            stop("No data to plot.")
        }
        p1<-plotFun2(dat=dat[[1]],de=dat[[2]],genes.to.plot=dat[[3]],relexpr=dat[[4]],timeVar=timeVar,id=id,
             ind=ind,m=figure.margins,figure.dimensions=figure.dimensions,line.size=line.size,title.size=title.size,plot.type=plot.type,screen.it=screen.it)
    }

    # 3. when selecting genes and color based on cluster
    if(plotting.feature=="gene" & color.by=="cluster"){
        if(is.null(genes.to.plot[1]) | is.na(genes.to.plot[1]) | genes.to.plot==0){
            stop("Specify the genes to be plotted!")
        }
        if(length(genes.to.plot)>432){
            stop("The maximum number of genes.to.plot is 432! To plot more, use a custom color palette with more colors.")
        }
        if(length(genes.to.plot)>length(color.palette)){
            stop("The number of genes are more than the number of colors. Use more colors in color.palette!")
        }
        print("Plotting selected genes and coloring by cluster.")
        dat<-fun3(Dat=Dat,relexpr=relexpr,de=DE,genes.to.plot=genes.to.plot,colors=gencol,timeVar=timeVar)
        if(nrow(dat[[1]])==0){
            stop("No data to plot.")
        }
        p1<-plotFun2(dat=dat[[1]],de=dat[[2]],genes.to.plot=dat[[3]],relexpr=dat[[4]],timeVar=timeVar,id=id,
             ind=ind,m=figure.margins,figure.dimensions=figure.dimensions,line.size=line.size,title.size=title.size,plot.type=plot.type,screen.it=screen.it)
    }

    # 4. when selecting cluster and color based on membership cutoff
    if(plotting.feature=="cluster" & color.by=="membership"){
        
        if(is.null(cluster.to.plot[1]) | is.na(cluster.to.plot[1])){
            stop("Specify the cluster to be plotted!")
        }
        if(length(cluster.to.plot)>1){
            print(paste("Only one cluster can be plotted at each time. Only cluster ",cluster.to.plot[1]," will be plotted.",sep=""))
        }
        if(match(cluster.to.plot,as.numeric(as.character(DE$Cluster)),nomatch=0)==0){
            stop("This cluster does not exist in the data.")
        }
        if(membership.cut<0 | membership.cut>1){
            stop("membership.cut must be a number in [0,1].")
        }
        if(color.by!="membership"){
            color.by<-"membership"
            print("The value of parameter color.by is changed to membership")
        }
        print("Plotting selected cluster and coloring by membership cutoff. To color by top membership scores set membership.cut = NULL.")
        dat<-fun4(Dat=Dat,relexpr=relexpr,de=DE,cluster.to.plot=cluster.to.plot,memb.cut=membership.cut,timeVar=timeVar)
        if(nrow(dat[[1]])==0){
            stop("No data to plot.")
        }
        p1<-plotFun2(dat=dat[[1]],de=dat[[2]],genes.to.plot=dat[[3]],relexpr=dat[[4]],timeVar=timeVar,id=id,
             ind=ind,m=figure.margins,figure.dimensions=figure.dimensions,line.size=line.size,title.size=title.size,plot.type=plot.type,screen.it=screen.it)
    }

    # 5. when selecting cluster and color based on membership top
    if(plotting.feature=="cluster" & color.by=="membership" & is.null(membership.cut)){
        
        if(match(cluster.to.plot,as.numeric(as.character(DE$Cluster)),nomatch=0)==0){
            stop("This cluster does not exist in the data.")
        }
        if(is.null(top[1]) | is.na(top[1]) | top==0){
            stop("Specify top parameter. It should be a positive integer.")
        }
        if(color.by!="membership"){
            color.by<-"membership"
            print("The value of parameter color.by is changed to membership")
        }
        print("Plotting selected cluster and coloring by top membership scores.")
        dat<-fun5(Dat=Dat,relexpr=relexpr,de=DE,cluster.to.plot=cluster.to.plot,top=top.hits,timeVar=timeVar)
        if(nrow(dat[[1]])==0){
            stop("No data to plot.")
        }
        p1<-plotFun2(dat=dat[[1]],de=dat[[2]],genes.to.plot=dat[[3]],relexpr=dat[[4]],timeVar=timeVar,id=id,
             ind=ind,m=figure.margins,figure.dimensions=figure.dimensions,line.size=line.size,title.size=title.size,plot.type=plot.type,screen.it=screen.it)
    }

    # 6. when selecting cluster and top x genes (by top qval) and color based on gene
    if(plotting.feature=="cluster:gene" & color.by=="gene" & is.null(qval.cut)){
        if(match(cluster.to.plot,as.numeric(as.character(DE$Cluster)),nomatch=0)==0){
            stop("This cluster does not exist in the data.")
        }
        if(is.null(top[1]) | is.na(top[1]) | top==0){
            stop("Specify top parameter. It should be a positive integer.")
        }
        print("Plotting selected cluster-specific genes and coloring by gene. Only the top genes in terms of qval are used.")
        dat<-fun6(Dat=Dat,relexpr=relexpr,de=DE,cluster.to.plot=cluster.to.plot,colors=gencol,top=top.hits,timeVar=timeVar)
        if(nrow(dat[[1]])==0){
            stop("No data to plot.")
        }
        p1<-plotFun2(dat=dat[[1]],de=dat[[2]],genes.to.plot=dat[[3]],relexpr=dat[[4]],timeVar=timeVar,id=id,
             ind=ind,m=figure.margins,figure.dimensions=figure.dimensions,line.size=line.size,title.size=title.size,plot.type=plot.type,screen.it=screen.it)
    }

    # 7. when selecting cluster and top x genes (by cutoff qval) and color based on gene
    if(plotting.feature=="cluster:gene" & color.by=="gene"){
        if(match(cluster.to.plot,as.numeric(as.character(DE$Cluster)),nomatch=0)==0){
            stop("This cluster does not exist in the data.")
        }
        if(qval.cut<=0 | qval.cut>1){
            stop("qval.cut must be a number in (0,1].")
        }
        print("Plotting selected cluster-specific genes and coloring by gene. Only the genes that pass the qval citerion are used. To use the top criterion set qval.cut = NULL.")
        dat<-fun7(Dat=Dat,relexpr=relexpr,de=DE,cluster.to.plot=cluster.to.plot,colors=gencol,qcut=qval.cut,timeVar=timeVar)
        if(nrow(dat[[1]])==0){
            stop("No data to plot.")
        }
        p1<-plotFun2(dat=dat[[1]],de=dat[[2]],genes.to.plot=dat[[3]],relexpr=dat[[4]],timeVar=timeVar,id=id,
             ind=ind,m=figure.margins,figure.dimensions=figure.dimensions,line.size=line.size,title.size=title.size,plot.type=plot.type,screen.it=screen.it)
    }

    # 8. when selecting cluster and top x genes (by cutoff qval) and color based on membership
    if(plotting.feature=="cluster:gene" & color.by=="membership"){
        if(match(cluster.to.plot,as.numeric(as.character(DE$Cluster)),nomatch=0)==0){
            stop("This cluster does not exist in the data.")
        }
        if(qval.cut<=0 | qval.cut>1){
            stop("qval.cut must be a number in (0,1].")
        }
        if(membership.cut<0 | membership.cut>1){
            stop("membership.cut must be a number in [0,1].")
        }
        print("Plotting selected cluster-specific genes and coloring by membership. Only the genes that pass the qval citerion are used. To use the top criterion set qval.cut = NULL.")
        dat<-fun8(Dat=Dat,relexpr=relexpr,de=DE,cluster.to.plot=cluster.to.plot,colors=gencol,qcut=qval.cut,memb.cut=membership.cut,timeVar=timeVar)
        if(nrow(dat[[1]])==0){
            stop("No data to plot.")
        }
        p1<-plotFun2(dat=dat[[1]],de=dat[[2]],genes.to.plot=dat[[3]],relexpr=dat[[4]],timeVar=timeVar,id=id,
             ind=ind,m=figure.margins,figure.dimensions=figure.dimensions,line.size=line.size,title.size=title.size,plot.type=plot.type,screen.it=screen.it)
    }


    # 9. when selecting cluster and top x genes (by top qval) and color based on membership
    if(plotting.feature=="cluster:gene" & color.by=="gene" & is.null(qval.cut)){
        if(match(cluster.to.plot,as.numeric(as.character(DE$Cluster)),nomatch=0)==0){
            stop("This cluster does not exist in the data.")
        }
        if(membership.cut<0 | membership.cut>1){
            stop("membership.cut must be a number in [0,1].")
        }
        if(is.null(top[1]) | is.na(top[1]) | top==0){
            stop("Specify top parameter. It should be a positive integer.")
        }
        print("Plotting selected cluster-specific genes and coloring by membership. Only the top genes in terms of qval are used.")
        dat<-fun9(Dat=Dat,relexpr=relexpr,de=DE,cluster.to.plot=cluster.to.plot,colors=gencol,top=top.hits,memb.cut=membership.cut,timeVar=timeVar)
        if(nrow(dat[[1]])==0){
            stop("No data to plot.")
        }
        p1<-plotFun2(dat=dat[[1]],de=dat[[2]],genes.to.plot=dat[[3]],relexpr=dat[[4]],timeVar=timeVar,id=id,
             ind=ind,m=figure.margins,figure.dimensions=figure.dimensions,line.size=line.size,title.size=title.size,plot.type=plot.type,screen.it=screen.it)
    }
    
 return(list(Obj=obj,Data=dat))
}




#' It plots the two-branch data with plotly or ggplot.
#'
#' It plots the two-branch data with plotly or ggplot.
#' @param dat data frame. A data frame with the expression data of the genes to be plotted.
#'   The genes are in columns.
#' @param de data frame. The differential expression results for the genes to be plotted.
#' @param genes.to.plot character. The genes to be plotted.
#' @param timeVar character. A name for the x variable of the regression model.
#' @param id character. An id to be used in the plot's filename.
#' @param ind list. Each component contains an index, internally calculated in plotBranch(),
#'   that separates the branch-specific data.
#' @param m list. The figure margins.
#' @param figure.dimensions list. The figure dimensions.
#' @param line.size numeric. The line width in ggplot.
#' @param title.size numeric. The font size for the plot title.
#' @param plot.type character. One of interactive (plotly) or static (ggplot).
#' @param screen.it logical. If TRUE, the plot will be shown on screen and not stored.
#'   Default is FALSE.
#' @keywords plotFun2
#' @return A plotly or a ggplot object.
#'
plotFun2<-function(dat,relexpr,de,genes.to.plot,timeVar,id,ind,m,figure.dimensions,line.size,title.size,plot.type,screen.it){
    
    dat<-rbind(dat,relexpr)
    if(plot.type=="interactive"){
        plo<-as.list(rep(0,length(ind)))
        for(k in 1:length(plo)){
            plo[[k]]<-plotly::plot_ly(dat[ind[[k]],],
                                      x=as.numeric(as.character(dat[ind[[k]],ncol(dat)])),
                                      width = figure.dimensions$width, height = figure.dimensions$height) %>%
                                      layout(title=list(text=paste("Modeling of Expression ~ ",timeVar," for States ",paste(names(ind),collapse=" & "),sep=""),
                                      font=list(color="black",size=title.size)),
                                      xaxis=list(title=paste("<b> ",timeVar," </b>",sep="")),
                                      yaxis=list(title=paste("<b> Normalized Expression </b>",sep=""),
                                      range=range(as.numeric(as.character(dat[,1:(ncol(dat)-2)])))+c(-0.5,0.5)),
                                      autosize = F, margin = m)
        
            for(i in 1:(ncol(dat)-2)){
                
                dd<-de[GeneNames(de$Gene)==colnames(dat)[i],]
                
                plo[[k]]<-plo[[k]] %>% add_trace(y=as.numeric(as.character(dat[ind[[k]],i])),
                                     mode="lines",type="scatter",
                                     line=list(color=names(genes.to.plot)[i]),
                                     name=genes.to.plot[i],
                                     text=paste(genes.to.plot[i]," ~ ",timeVar,
                                           "<br> qval = ",paste(dd$qval," (",dd$Branch,")",collapse=";"),
                                           "<br> Cluster: ",as.numeric(as.character(dd$Cluster[1])),
                                           "<br> Membership: ",round(as.numeric(as.character(dd$Membership[1])),3),
                                           "<br> State: ",names(ind)[k],sep=""))

                if(length(genes.to.plot)==1){
                     plo[[k]]<-plo[[k]] %>% add_trace(y=as.numeric(as.character(dat[ind[[k]]+max(ind[[length(ind)]]),i])),
                     mode="markers",type="scatter",
                     symbols="circle",color=I("gray60"),
                     marker = list(size = 5))

          
                }
            }
        }
       
        pall<-subplot(plo)
        
        if(!screen.it){
            # plot and message
            saveWidget(pall,paste(outdir,"branchPlotlyLines_",id,".html",sep=""),selfcontained = FALSE)
            print(paste("Download plot file: ",outdir,"branchPlotlyLines_",id,".html & plot folder: ",outdir,"branchPlotlyLines_",id,"_files",sep=""))
        } else {
            pall
        }
   
    } else {
        
        dat<-ggdata2(Dat=dat,genes=genes.to.plot,timeVar=timeVar,ind=ind)
        vn<-as.character(names(genes.to.plot))
        vn<-vn[sort.list(genes.to.plot)]
        uu<-names(ind)
        pall<-as.list(rep(0,length(uu)))
        
        for(i in 1:length(uu)){
            d1<-dat[dat$Ind==uu[i],]
            d1.1<-d1[d1$Type=="Norm",]
            d1.2<-d1[d1$Type=="Raw",]
    
            if(length(genes.to.plot)==1){
                pall[[i]]<-ggplot(d1.1,aes(x=Age,y=Expr,color=Gene)) + geom_line(size=line.size) + geom_point(data=d1.2[,2:1],color="black") +
                           scale_color_manual(values=vn) +
                           xlab(timeVar) + ylab("Normalized Expression") + ggtitle(paste("Modeling of Expression ~ ",timeVar," for State ",uu[i],sep="")) +
                           theme_bw() + geom_hline(yintercept=0) + theme(plot.title = element_text(size=title.size))
                if(length(unique(d1.1$Gene))>20){
                    pall[[i]]<-pall[[i]]+theme(legend.position = "none")
                }

            } else {
                pall[[i]]<-ggplot(d1.1,aes(x=Age,y=Expr,color=Gene)) + geom_line(size=line.size) +
                           scale_color_manual(values=vn) +
                           xlab(timeVar) + ylab("Normalized Expression") + ggtitle(paste("Modeling of Expression ~ ",timeVar," for State ",names(ind)[1],sep="")) +
                           theme_bw() + geom_hline(yintercept=0) + theme(plot.title = element_text(size=title.size))
                if(length(unique(d1.1$Gene))>20){
                    pall[[i]]<-pall[[i]]+theme(legend.position = "none")
                }
            }
        }
        
        if(!screen.it){
            # plot and message
            pdf(paste(outdir,"branchggplotLines_",id,".pdf",sep=""),width=figure.dimensions$width,height=figure.dimensions$height)
                do.call("grid.arrange", c(pall, nrow=1))
            dev.off()
            print(paste("Download plot: ",outdir,"branchggplotLines_",id,".pdf",sep=""))
        } else {
            do.call("grid.arrange", c(pall, nrow=1))
        }
    }
    
 return(dat)
}
    



#' Generates a heatmap of single-branch Monocle's DE genes with updated clusters using the membership scores
#'   or a manually derived clustering scheme.
#'
#' Generates a heatmap of single-branch Monocle's DE genes with updated clusters using the membership scores
#'   or a manually derived clustering scheme.
#' @param obj object. A Monocle object.
#' @param de data frame. The DE statistics generated by MonocleDE().
#' @param id character. An ID to be used in the filename of the generated heatmap.
#' @param outdir character. A folder to store the heatmap.
#' @param trend_formula character. A formula specifying the model formula used in fitting the spline curve
#'   for each gene/feature. Default is ~sm.ns(Pseudotime,df=3).
#' @param hmcols character. The color scheme for drawing the heatmap. Default is NULL that uses the blue2green2red
#'   scheme of the R package colorRamps.
#' @param show_rownames logical. Whether to show the names and the membership for each gene on the heatmap.
#'   Default is FALSE.
#' @param add_annotation_row data frame. The manual clustering scheme can be defined there. It should be a data
#'   frame whose rownames are the rownames(obj) and a single column with the updated clusters. Default is NULL
#'   under which the membership scores are used for refining the clustering (see membership.cut below).
#' @param add_annotation_col data frame.  Additional annotations to show for each row in the heatmap. Must be
#'   a data frame with one row for each row in the pData table of obj with matching IDs. Default is NULL.
#' @param scale.range numeric. The minimum and the maximum value (in standard deviations) to show in the
#'   heatmap. Values smaller / larger than this are set to the min / max. Default is c(-3,3).
#' @param plot.width numeric. The width of the graphics region in inches. Default is 7.
#' @param plot.height numeric. The height of the graphics region in inches. Default is 7.
#' @param membership.cut list. Each component may contain a single value or many sorted values depending on the
#'   number of cutoffs (and ultimately subclusters) that we want to define. The genes with membership scores above,
#'   in-between and below the respective value(s) are separated into subclusters. If one component is present, the
#'   cutoff(s) will be applied on all pre-existing clusters. Otherwise, it should contain as many components as the
#'   number of the pre-existing clusters (in which case the component names should be the cluster numbers). Default
#'   is 1 that will not update any clusters.
#' @param membership2 logical. If TRUE (default), it will compare the top two membership scores and it will assign the cluster
#'   and the associated membership of the maximum membership. It preceeds the membership.cut criterion, so if both are
#'   used, the function will update the clusters based on membership2 first and then will run a second update using
#'   membership.cut.
#' @keywords simpleReHeat
#' @return A list with components a Monocle object with the DE genes only, a data frame of normalized and smoothed
#'   gene expression profiles (all and DE only) and the list of DE genes and all memberships.
#'
simpleReHeat<-function(obj,id,outdir,trend_formula="~sm.ns(Pseudotime,df=3)",
                     scale.range=c(-3,3),hmcols=NULL,
                     show_rownames=FALSE,membership.cut=1,
                     membership2=TRUE,
                     add_annotation_row=NULL,add_annotation_col=NULL,
                     plot.width=7,plot.height=7){
    
    de<-obj$DE
    all<-obj$NormDataAll
    
    if(membership2){
        
        x<-ifelse(de$Membership<de$Membership2,TRUE,FALSE)
        tt<-table(c(x,TRUE,FALSE))
        print(paste("Updating clusters using the maximum membership score. The clusters of ",as.numeric(as.character(tt[names(tt)==TRUE]))-1," genes have been updated!",sep=""))
        x1<-ifelse(as.numeric(as.character(de$Membership))<as.numeric(as.character(de$Membership2)),de$Cluster2,de$Cluster)
        x2<-ifelse(as.numeric(as.character(de$Membership))<as.numeric(as.character(de$Membership2)),de$Membership2,de$Membership)
        de$Cluster<-x1
        de$Membership<-x2
    
    }
    
    if(is.null(add_annotation_row[1])){
        
        if(length(membership.cut)==1 & length(unique(de$Cluster))>1){
            print(paste("Setting membership.cut = ",membership.cut," for all clusters!",sep=""))
            membership.cut<-as.list(rep(membership.cut,length(unique(de$Cluster))))
            names(membership.cut)<-sort(unique(de$Cluster))
        }
        if(length(membership.cut)>1 & length(membership.cut)!=length(unique(de$Cluster))){
            stop("Parameter membership.cut should either be 1 value (common to all clusters) or should contain as many values at the number of clusters.")
        }
        if(is.null(names(membership.cut))){
            stop("Provide cluster names for each membership.cut value!")
        }

        de1<-matrix(0,1,ncol(de))
        for(i in 1:length(membership.cut)){
            x<-de[de$Cluster==as.numeric(names(membership.cut)[i]),]
        
            mcut<-membership.cut[[which(as.numeric(names(membership.cut))==unique(x$Cluster))]]
            mcut<-sort(unique(mcut[mcut>0 | mcut<=1]))
            if(length(mcut)==0){
            stop("membership.cut takes values i (0,1]!")
            }
            if(any(mcut<1)){
            print(paste("Splitting cluster ",as.numeric(names(membership.cut)[i])," using the membership scores!",sep=""))
            mcut<-mcut[mcut<1]
                x$Cluster<-x$Cluster+findInterval(x$Membership,mcut)/100
            }
            de1<-rbind(de1,as.matrix(x))
        
        }
        de<-data.frame(de1[-1,])
        colnames(de)<-colnames(obj$DE)
        
    } else {
        print("Updating clusters using the user defined cluster information!")
        mm<-match(de$Gene,rownames(add_annotation_row))
        de$Cluster<-add_annotation_row$Cluster
    }
    de<-de[sort.list(de$Cluster),]
    norm<-obj$NormData[match(de$Gene,rownames(obj$NormData)),]
    all<-all[match(de$Gene,rownames(all)),]
    memb<-obj$Membership[match(de$Gene,rownames(obj$Membership)),]
    obj<-obj$Obj
    obj<-obj[match(de$Gene,rownames(obj)),]
   
    annot_row<-data.frame(Cluster=factor(de$Cluster))
    rownames(annot_row)<-fData(obj)$Gene
    
    if (is.null(hmcols)) {
            bks <- seq(-3.1, 3.1, by = 0.1)
            hmcols <- blue2green2red(length(bks) - 1)
    }
    
   gap_row <- match(unique(annot_row$Cluster), annot_row$Cluster)
   pdf(paste(outdir,"DE_ReHeatmap_",id,".pdf",sep=""),height=plot.height,width=plot.width)
        print(annotated_plot_pseudotime_heatmap(obj,cluster_rows=FALSE,
                                                hmcols = hmcols,add_annotation_row=annot_row,
                                                add_annotation_col=add_annotation_col,
                                                show_rownames=show_rownames,use_gene_short_name=TRUE,norm_method="log",
                                                scale_max=scale.range[2],scale_min=scale.range[1],trend_formula=trend_formula,
                                                num_clusters = length(unique(de$Cluster)),return_heatmap=FALSE,gaps_row=(gap_row[-1]-1)))
    dev.off()

    # message
    print(paste("Download plot: ",outdir,"DE_ReHeatmap_",id,".pdf",sep=""))

    oo<-list(Obj=obj,NormDataAll=all,NormData=norm,DE=de,Membership=memb)
     
  return(oo)
}
     


#' A helper function to generate the plot of plotMembership().
#'
#' A helper function to generate the plot of plotMembership().
#' @param data data frame. The DE output of simpleHeat().
#' @param varname character. The variable for plotting.
#' @param groupnames character. The grouping variable, i.e. the estimated clusters.
#' @keywords data_summary
#' @return A data frame for ggplot plotting.
#'
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}



#' Plots the cluster memberships of single or branched data.
#'
#' Plots the cluster memberships of single or branched data. It is useful to check the
#'   distribution of memberships in different number of clusters or smoothing parameters.
#' @param obj data frame. The DE output of simpleHeat() or branchHeat().
#' @param id character. An id to be used in the filename of the stored plot.
#' @param outdir character. A folder to store the plot.
#' @param plot.width numeric. The width of the graphics region in inches. Default is 7.
#' @param plot.height numeric. The height of the graphics region in inches. Default is 7.
#' @param screen.it logical. If TRUE, the plot will be shown on screen and not stored.
#'   Default is FALSE.
#' @keywords plotMemberships
#' @return A plot.
#'
plotMembership<-function(obj,id,outdir,plot.width,plot.height,screen.it=FALSE){
    
    if(!dir.exists(outdir)){
        dir.create(outdir)
        print(paste("Folder ",outdir," has been created to store the plotMembership() plots.",sep=""))
    }
    
    if(is.data.frame(obj)){
        obj<-list(obj)
        names(obj)<-paste("Clusters",max(obj[[1]]$Cluster),sep="")
    }
    
    df<-as.list(rep(0,length(obj)))
    for(k in 1:length(obj)){
    xx<-data.frame(unique(cbind(obj[[k]]$Gene,obj[[k]]$Cluster,obj[[k]]$Membership)))
    colnames(xx)<-c("Gene","Cluster","Membership")
    obj[[k]]<-xx
        obj[[k]]$Membership<-as.numeric(as.character(obj[[k]]$Membership))
        obj[[k]]$Cluster<-factor(obj[[k]]$Cluster)
        df[[k]] <- data_summary(obj[[k]], varname="Membership",groupnames="Cluster")
    df[[k]]$sd[is.na(df[[k]]$sd)]<-0
        for(i in 1:nrow(df[[k]])){
            xx<-df[[k]]$Membership[i]+c(-df[[k]]$sd[i],df[[k]]$sd[i])
            if(xx[1]<0){
                df[[k]]$sd[i]<-df[[k]]$Membership[i]
            }
            if(xx[2]>1){
                df[[k]]$sd[i]<-1-df[[k]]$Membership[i]
            }
        }
        df[[k]]<-cbind(df[[k]],ID=names(obj)[k])
    }
    df<-do.call(rbind,df)
    
    df$Membership<-as.numeric(as.character(df$Membership))
    df$Cluster<-factor(df$Cluster)
    ss<-unlist(lapply(strsplit(unique(df$ID),"s"),function(x) x[length(x)]))
    df$ID<-factor(df$ID,levels=paste("Clusters",sort(as.numeric(as.character(ss))),sep=""))
    p1<- ggplot(df, aes(x=ID, y=Membership,fill=Cluster)) +
                geom_bar(stat="identity", color="black",position=position_dodge()) + theme_bw() +
                geom_errorbar(aes(ymin=Membership-sd, ymax=Membership+sd), width=.2,position=position_dodge(.9))+ylim(0,1)
    
    aa<-aggregate(df$Membership,list(df$ID),mean)
    aa[,2]<-as.numeric(as.character(aa[,2]))
    dd<-t(matrix(unlist(strsplit(as.character(aa[,1]),"rs",fixed=T)),nrow=2))
    aa[,1]<-factor(aa[,1],levels=paste("Clusters",sort(as.numeric(as.character(dd[,2]))),sep=""))
    colnames(aa)<-c("Cluster","Membership")
    if(nrow(aa)>1){
        p2<-ggplot(aa,aes(x=Cluster,y=Membership,group=1))+geom_line()+geom_point()+theme_bw()+labs(y="Average Membership")
    }

    if(!screen.it){
        pdf(paste(outdir,"MembershipPlots_",id,".pdf",sep=""),height=plot.height,width=plot.width)
            print(p1)
        if(nrow(aa)>1){
            print(p2)
        }
        dev.off()
    } else {
    print(p1)
    if(nrow(aa)>1){
        x11()
        print(p2)
    }
    }

    # message
    print(paste("Download plot: ",outdir,"MembershipPlots_",id,".pdf",sep=""))

}





#' Edit of the plot_pseudotime_heatmap() of Monocle to accomodate for the annotation_row parameter.
#'
#' Edit of the plot_pseudotime_heatmap() of Monocle to accomodate for the annotation_row parameter.
#' @param cds_subset object. A Monocle object.
#' @param cluster_rows logical. Whether to cluster the rows of the heatmap.
#' @param hclust_method character. The method used by pheatmap to perform hirearchical clustering of
#'   the rows.
#' @param num_clusters numeric. Number of clusters for the heatmap of branch genes.
#' @param hmcols character. The color scheme for drawing the heatmap.
#' @param add_annotation_row data frame. Additional annotations to show for each row in the heatmap.
#'   Must be a dataframe with one row for each row in the fData table of cds_subset, with matching IDs.
#' @param add_annotation_col data frame. Additional annotations to show for each column in the heatmap.
#'   Must be a dataframe with one row for each cell in the pData table of cds_subset, with matching IDs.
#' @param show_rownames logical. Whether to show the names for each row in the table.
#' @param use_gene_short_name logical. Whether to use the short names for each row. If FALSE, uses row IDs
#'   from the fData table.
#' @param norm_method character. Determines how to transform expression values prior to rendering.
#' @param scale_max numeric. The maximum value (in standard deviations) to show in the heatmap. Values larger
#'   than this are set to the max.
#' @param scale_min numeric. The minimum value (in standard deviations) to show in the heatmap. Values smaller
#'   than this are set to the min.
#' @param trend_formula character. A formula string specifying the model used in fitting the spline curve for
#'   each gene/feature.
#' @param return_heatmap logical. Whether to return the pheatmap object to the user.
#' @param cores numeric. Number of cores to use when smoothing the expression curves shown in the heatmap.
#' @keywords annotated_plot_pseudotime_heatmap
#' @return A heatmap object or a message.
#'
annotated_plot_pseudotime_heatmap<-function (cds_subset, cluster_rows = TRUE, hclust_method = "ward.D2",
    num_clusters = 6, hmcols = NULL, add_annotation_row = NULL,
    add_annotation_col = NULL, show_rownames = FALSE, use_gene_short_name = TRUE,
    norm_method = c("log", "vstExprs"), scale_max = 3, scale_min = -3,
    trend_formula = "~sm.ns(Pseudotime, df=3)", return_heatmap = FALSE,
    cores = 1,gaps_row=NULL)
{
    num_clusters <- min(num_clusters, nrow(cds_subset))
    pseudocount <- 1
    newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime),
        max(pData(cds_subset)$Pseudotime), length.out = 100))
    m <- genSmoothCurves(cds_subset, cores = cores, trend_formula = trend_formula,
        relative_expr = T, new_data = newdata)
    m = m[!apply(m, 1, sum) == 0, ]
    norm_method <- match.arg(norm_method)
    if (norm_method == "vstExprs" && is.null(cds_subset@dispFitInfo[["blind"]]$disp_func) ==
        FALSE) {
        m = vstExprs(cds_subset, expr_matrix = m)
    }
    else if (norm_method == "log") {
        m = log2(m + pseudocount)
    }
    m = m[!apply(m, 1, sd) == 0, ]
    m = Matrix::t(scale(Matrix::t(m), center = TRUE))
    m = m[is.na(row.names(m)) == FALSE, ]
    m[is.nan(m)] = 0
    m[m > scale_max] = scale_max
    m[m < scale_min] = scale_min
    heatmap_matrix <- m
    row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
    row_dist[is.na(row_dist)] <- 1
    if (is.null(hmcols)) {
        bks <- seq(-3.1, 3.1, by = 0.1)
        hmcols <- blue2green2red(length(bks) - 1)
    }
    else {
        bks <- seq(-3.1, 3.1, length.out = length(hmcols))
    }
    ph <- pheatmap(heatmap_matrix, useRaster = T, cluster_cols = FALSE,
        cluster_rows = cluster_rows, show_rownames = F, show_colnames = F,
        clustering_distance_rows = row_dist, clustering_method = hclust_method,
        cutree_rows = num_clusters, silent = TRUE, filename = NA,
        breaks = bks, border_color = NA, color = hmcols,gaps_row=gaps_row)
    
    annotation_row<-add_annotation_row
    
    
    if (!is.null(add_annotation_col)) {
        if (nrow(add_annotation_col) != 100) {
            stop("add_annotation_col should have only 100 rows (check genSmoothCurves before you supply the annotation data)!")
        }
        annotation_col <- add_annotation_col
    }
    else {
        annotation_col <- NA
    }
    if (use_gene_short_name == TRUE) {
        if (is.null(fData(cds_subset)$gene_short_name) == FALSE) {
            feature_label <- as.character(fData(cds_subset)[row.names(heatmap_matrix),
                "gene_short_name"])
            feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
            row_ann_labels <- as.character(fData(cds_subset)[row.names(annotation_row),
                "gene_short_name"])
            row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
        }
        else {
            feature_label <- row.names(heatmap_matrix)
            row_ann_labels <- row.names(annotation_row)
        }
    }
    else {
        feature_label <- row.names(heatmap_matrix)
        if (!is.null(annotation_row))
            row_ann_labels <- row.names(annotation_row)
    }
    row.names(heatmap_matrix) <- feature_label
    if (!is.null(annotation_row))
        row.names(annotation_row) <- row_ann_labels
    colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
    ph_res <- pheatmap(heatmap_matrix[, ], useRaster = T, cluster_cols = FALSE,
        cluster_rows = cluster_rows, show_rownames = show_rownames,
        show_colnames = F, clustering_distance_rows = row_dist,
        clustering_method = hclust_method, cutree_rows = num_clusters,
        annotation_row = annotation_row, annotation_col = annotation_col,
        treeheight_row = 20, breaks = bks, fontsize = 6, color = hmcols,
        border_color = NA, silent = TRUE, filename = NA,gaps_row=gaps_row)
    grid::grid.rect(gp = grid::gpar("fill", col = NA))
    grid::grid.draw(ph_res$gtable)
    if (return_heatmap) {
        return(ph_res)
    } else {
        return("The clusters have been updated")
    }
    
}




#' Generates a heatmap of two-branch Monocle's DE genes with updated clusters using the membership scores
#'   or a manually derived clustering scheme.
#'
#' Generates a heatmap of two-branch Monocle's DE genes with updated clusters using the membership scores
#'   or a manually derived clustering scheme.
#' @param obj object. A Monocle object.
#' @param de data frame. The DE statistics generated by MonocleDE().
#' @param id character. An ID to be used in the filename of the generated heatmap.
#' @param outdir character. A folder to store the heatmap.
#' @param hmcols character. The color scheme for drawing the heatmap. Default is NULL that uses the blue2green2red
#'   scheme of the R package colorRamps.
#' @param show_rownames logical. Whether to show the names and the membership for each gene on the heatmap.
#'   Default is FALSE.
#' @param add_annotation_row data frame. The manual clustering scheme can be defined there. It should be a data
#'   frame whose rownames are the rownames(obj) and a single column with the updated clusters. Default is NULL
#'   under which the membership scores are used for refining the clustering (see membership.cut below).
#' @param add_annotation_col data frame.  Additional annotations to show for each row in the heatmap. Must be
#'   a data frame with one row for each row in the pData table of obj with matching IDs. Default is NULL.
#' @param scale.range numeric. The minimum and the maximum value (in standard deviations) to show in the
#'   heatmap. Values smaller / larger than this are set to the min / max. Default is c(-3,3).
#' @param plot.width numeric. The width of the graphics region in inches. Default is 7.
#' @param plot.height numeric. The height of the graphics region in inches. Default is 7.
#' @param membership.cut list. Each component may contain a single value or many sorted values depending on the
#'   number of cutoffs (and ultimately subclusters) that we want to define. The genes with membership scores above,
#'   in-between and below the respective value(s) are separated into subclusters. If one component is present, the
#'   cutoff(s) will be applied on all pre-existing clusters. Otherwise, it should contain as many components as the
#'   number of the pre-existing clusters (in which case the component names should be the cluster numbers). Default
#'   is 1 that will not update any clusters.
#' @param membership2 logical. If TRUE (default), it will compare the top two membership scores and it will assign the cluster
#'   and the associated membership of the maximum membership. It preceeds the membership.cut criterion, so if both are
#'   used, the function will update the clusters based on membership2 first and then will run a second update using
#'   membership.cut.
#' @param screen.it logical.  If TRUE, the plot will be shown on screen and not stored. Default is FALSE.
#' @keywords branchReHeat
#' @return A list with components a Monocle object with the DE genes only, a data frame of normalized and smoothed
#'   gene expression profiles (all and DE only) and the list of DE genes and all memberships.
#'
branchReHeat<-function(obj,id,outdir,
                     scale.range=c(-3,3),hmcols=NULL,
                     show_rownames=FALSE,membership.cut=1,
                     membership2=TRUE,
                     add_annotation_row=NULL,add_annotation_col=NULL,
                     plot.width=7,plot.height=7,screen.it=FALSE){

    de<-obj$DE
    all<-obj$NormDataAll

    if(membership2){
        
    xx<-data.frame(Gene=de$Gene,Membership=de$Membership,Membership2=de$Membership2)
    xx<-unique(xx)
        x<-ifelse(xx$Membership<xx$Membership2,TRUE,FALSE)
        tt<-table(c(x,TRUE,FALSE))
        print(paste("Updating clusters using the maximum membership score. The clusters of ",as.numeric(as.character(tt[names(tt)==TRUE]))-1," genes have been updated!",sep=""))
        x1<-ifelse(de$Membership<de$Membership2,de$Cluster2,de$Cluster)
        x2<-ifelse(de$Membership<de$Membership2,de$Membership2,de$Membership)
        de$Cluster<-x1
        de$Membership<-x2

    }

    if(is.null(add_annotation_row[1])){

        if(length(membership.cut)==1 & length(unique(de$Cluster))>1){
            print(paste("Setting membership.cut = ",membership.cut," for all clusters!",sep=""))
            membership.cut<-as.list(rep(membership.cut,length(unique(de$Cluster))))
            names(membership.cut)<-sort(unique(de$Cluster))
        }
        if(length(membership.cut)>1 & length(membership.cut)!=length(unique(de$Cluster))){
            stop("Parameter membership.cut should either be 1 value (common to all clusters) or should contain as many values at the number of clusters.")
        }
        if(is.null(names(membership.cut))){
            stop("Provide cluster names for each membership.cut value!")
        }

        de1<-matrix(0,1,ncol(de))
        for(i in 1:length(membership.cut)){
            x<-de[de$Cluster==as.numeric(names(membership.cut)[i]),]

            mcut<-membership.cut[[which(as.numeric(names(membership.cut))==unique(x$Cluster))]]
            mcut<-sort(unique(mcut[mcut>0 | mcut<=1]))
            if(length(mcut)==0){
                stop("membership.cut takes values i (0,1]!")
            }
            if(any(mcut<1)){
                print(paste("Splitting cluster ",as.numeric(names(membership.cut)[i])," using the membership scores!",sep=""))
                mcut<-mcut[mcut<1]
                x$Cluster<-x$Cluster+findInterval(x$Membership,mcut)/100
            }
            de1<-rbind(de1,as.matrix(x))

        }
        de<-data.frame(de1[-1,])
        colnames(de)<-colnames(obj$DE)

    } else {
        print("Updating clusters using the user defined cluster information!")
        mm<-match(de$Gene,rownames(add_annotation_row))
        de$Cluster<-add_annotation_row$Cluster
     }
     de<-de[sort.list(de$Cluster),]
     xx<-data.frame(Gene=de$Gene,Cluster=de$Cluster,Membership=de$Membership,Cluster2=de$Cluster2,Membership2=de$Membership2)
     xx<-unique(xx)
     norm<-obj$NormData[match(GeneNames(xx$Gene),rownames(obj$NormData)),]
     for(i in 1:(length(all)-1)){
         all[[i]]<-all[[i]][match(xx$Gene,rownames(all[[i]])),]
     }
     info<-all[[length(all)]]
     memb<-obj$Membership[match(xx$Gene,rownames(obj$Membership)),]
     obj<-obj$Obj
     obj<-obj[match(xx$Gene,rownames(obj)),]

     annot_row<-data.frame(Cluster=factor(xx$Cluster))
     rownames(annot_row)<-GeneNames(fData(obj)$Gene)

     if (is.null(hmcols)) {
             bks <- seq(-3.1, 3.1, by = 0.1)
             hmcols <- blue2green2red(length(bks) - 1)
     }

    gap_row <- match(unique(annot_row$Cluster), annot_row$Cluster)
    if(!screen.it){
        pdf(paste(outdir,"DE_ReHeatmap_",id,".pdf",sep=""),height=plot.height,width=plot.width)
            ph_res<-pheatmap(norm[,],useRaster=TRUE,cluster_cols=FALSE,
                             cluster_rows=FALSE,show_rownames=show_rownames,show_colnames=FALSE,
                             cutree_rows=length(unique(de$Cluster)),annotation_row=annot_row,
                             annotation_col=add_annotation_col,gaps_col = cumsum(info[1:(length(info)-1)]),breaks = bks,
                             color = hmcols,gaps_row=(gap_row[-1]-1))
            grid::grid.rect(gp = grid::gpar("fill", col = NA))
            grid::grid.draw(ph_res$gtable)
            print(ph_res)
         dev.off()

         # message
         print(paste("Download plot: ",outdir,"DE_ReHeatmap_",id,".pdf",sep=""))
    } else {
        ph_res<-pheatmap(norm[,],useRaster=TRUE,cluster_cols=FALSE,
                         cluster_rows=FALSE,show_rownames=show_rownames,show_colnames=FALSE,
                         cutree_rows=length(unique(de$Cluster)),annotation_row=annot_row,
                         annotation_col=add_annotation_col,gaps_col = cumsum(info[1:(length(info)-1)]),breaks = bks,
                         color = hmcols,gaps_row=(gap_row[-1]-1))
        grid::grid.rect(gp = grid::gpar("fill", col = NA))
        grid::grid.draw(ph_res$gtable)
        print(ph_res)
    }
     oo<-list(Obj=obj,NormDataAll=all,NormData=norm,DE=de,Membership=memb)

   return(oo)
}


#' Summarizes single cell expression profiles into pseudobulk expression profiles for Seurat objects.
#'
#' Summarizes single cell expression profiles into pseudobulk expression profiles for Seurat objects.
#' @param obj object. A seurat object.
#' @param metaVars character. The column names from obj@meta.data whose information will be stored
#'   in the design of the pseudobulk data.
#' @param replCol character. A column name from obj@meta.data whose information is used to split and
#'   generate the pseudobulk. If NULL, all cells are summarized, otherwise only the cells of each
#'   replCol level will be summarized.
#' @keywords sc_seurat
#' @return A list with the pseudobulk expression profiles, the design and the genes of each seurat object.
#'
sc_seurat<-function(obj,metaVars,replCol){
    
    EXPRS<-DESI<-GENES<-as.list(rep(0,length(obj)))
    names(EXPRS)<-names(DESI)<-names(GENES)<-nn
    for(i in 1:length(obj)){
        
        #print(paste("Now processing object ",nn[i],"...",sep=""))
        
        if(is.list(obj[[i]])){
            obj[[i]]<-obj[[i]][[1]]
        }
        if(!is.null(replCol)){
            w<-which(colnames(obj[[i]]@meta.data)==replCol)
            if(length(w)==0){
                stop(paste("The replCol does not exist in the obj@meta.data of object ",nn[i],sep=""))
            }
            tt<-table(as.character(obj[[i]]@meta.data[,w]))
            mat<-matrix(0,nrow(obj[[i]]@assays$RNA@counts),length(tt))
            for(j in 1:nrow(mat)){
                mat[j,]<-as.numeric(as.character(aggregate(as.numeric(as.character(obj[[i]]@assays$RNA@counts[j,])),list(as.character(obj[[i]]@meta.data[,w])),sum)[,2]))
            }
            mat<-data.frame(mat)
            rownames(mat)<-rownames(obj[[i]]@assays$RNA@counts)
            colnames(mat)<-paste(nn[i],"_",names(tt),sep="")
            desi<-matrix(0,1,length(metaVars))
            for(j in 1:length(metaVars)){
                w<-which(colnames(obj[[i]]@meta.data)==metaVars[j])
                if(length(w)==0){
                    desi[j]<-NA
                    print(paste("metaVar ",metaVars[j]," does not exist in the obj@meta.data of object ",nn[i],sep=""))
                } else {
                    tt1<-table(as.character(obj[[i]]@meta.data[,w]))
                    if(length(tt1)>1){
                        print(paste("metaVar ",metaVars[j]," contains two unique values in the obj@meta.data of object ",nn[i],". The data will be summarized!",sep=""))
                        desi[j]<-paste(as.character(names(tt1)),collapse="_")
                    } else {
                        desi[j]<-as.character(names(tt1))
                    }
                }
            }
            desi<-data.frame(desi)
            colnames(desi)<-metaVars
            rownames(desi)<-paste(nn[i],"_",names(tt),sep="")
        } else {
            mat<-matrix(rowSums(obj[[i]]@assays$RNA@counts),ncol=1)
            rownames(mat)<-rownames(obj[[i]]@assays$RNA@counts)
            colnames(mat)<-nn[i]
            desi<-matrix(0,1,length(metaVars))
            for(j in 1:length(metaVars)){
                w<-which(colnames(obj[[i]]@meta.data)==metaVars[j])
                if(length(w)==0){
                    desi[j]<-NA
                    print(paste("metaVar ",metaVars[j]," does not exist in the obj@meta.data of object ",nn[i],sep=""))
                } else {
                    tt1<-table(as.character(obj[[i]]@meta.data[,w]))
                    if(length(tt1)>1){
                        print(paste("metaVar ",metaVars[j]," contains two unique values in the obj@meta.data of object ",nn[i],". The data will be summarized!",sep=""))
                        desi[j]<-paste(as.character(names(tt1)),collapse="_")
                    } else {
                        desi[j]<-as.character(names(tt1))
                    }
                }
            }
            desi<-data.frame(desi)
            colnames(desi)<-metaVars
            rownames(desi)<-nn[i]
        }
        EXPRS[[i]]<-mat
        DESI[[i]]<-desi
        GENES[[i]]<-rownames(mat)
    }
    
 return(list(EXPRS,DESI,GENES))
}



#' Summarizes single cell expression profiles into pseudobulk expression profiles for SingleCellExperiment
#'   objects.
#'
#' Summarizes single cell expression profiles into pseudobulk expression profiles for SingleCellExperiment
#'   objects.
#' @param obj object. A SingleCellExperiment object.
#' @param metaVars character. The column names from pData(obj) whose information will be stored
#'   in the design of the pseudobulk data.
#' @param replCol character. A column name from pData(obj) whose information is used to split and generate
#'   the pseudobulk. If NULL, all cells are summarized, otherwise only the cells of each replCol level
#'   will be summarized.
#' @keywords sc_sce
#' @return A list with the pseudobulk expression profiles, the design and the genes of each SingleCellExperiment
#'   object.
#'
sc_sce<-function(obj,metaVars,replCol){
    EXPRS<-DESI<-GENES<-as.list(rep(0,length(obj)))
    names(EXPRS)<-names(DESI)<-names(GENES)<-nn
    for(i in 1:length(obj)){
        
        #print(paste("Now processing object ",nn[i],"...",sep=""))
        
        if(is.list(obj[[i]])){
            obj[[i]]<-obj[[i]][[1]]
        }
        if(!is.null(replCol)){
            w<-which(colnames(pData(obj[[i]]))==replCol)
            if(length(w)==0){
                stop(paste("The replCol does not exist in the pData(obj) of object ",nn[i],sep=""))
            }
            tt<-table(as.character(pData(obj[[i]])[,w]))
            mat<-matrix(0,nrow(exprs(obj[[i]])),length(tt))
            for(j in 1:nrow(mat)){
                mat[j,]<-as.numeric(as.character(aggregate(as.numeric(as.character(exprs(obj[[i]])[j,])),list(as.character(pData(obj[[i]])[,w])),sum)[,2]))
            }
            mat<-data.frame(mat)
            rownames(mat)<-rownames(exprs(obj[[i]]))
            colnames(mat)<-paste(nn[i],"_",names(tt),sep="")
            desi<-matrix(0,1,length(metaVars))
            for(j in 1:length(metaVars)){
                w<-which(colnames(pData(obj[[i]]))==metaVars[j])
                if(length(w)==0){
                    desi[j]<-NA
                    print(paste("metaVar ",metaVars[j]," does not exist in the pData(obj) of object ",nn[i],sep=""))
                } else {
                    tt1<-table(as.character(pData(obj[[i]])[,w]))
                    if(length(tt1)>1){
                        print(paste("metaVar ",metaVars[j]," contains two unique values in the pData(obj) of object ",nn[i],". The data will be summarized!",sep=""))
                        desi[j]<-paste(as.character(names(tt1)),collapse="_")
                    } else {
                        desi[j]<-as.character(names(tt1))
                    }
                }
            }
            desi<-data.frame(desi)
            colnames(desi)<-metaVars
            rownames(desi)<-paste(nn[i],"_",names(tt),sep="")
        } else {
            mat<-matrix(rowSums(exprs(obj[[i]])),ncol=1)
            rownames(mat)<-rownames(exprs(obj[[i]]))
            colnames(mat)<-nn[i]
            desi<-matrix(0,1,length(metaVars))
            for(j in 1:length(metaVars)){
                w<-which(colnames(pData(obj[[i]]))==metaVars[j])
                if(length(w)==0){
                    desi[j]<-NA
                    print(paste("metaVar ",metaVars[j]," does not exist in the pData(obj) of object ",nn[i],sep=""))
                } else {
                    tt1<-table(as.character(pData(obj[[i]])[,w]))
                    if(length(tt1)>1){
                        print(paste("metaVar ",metaVars[j]," contains two unique values in the pData(obj) of object ",nn[i],". The data will be summarized!",sep=""))
                        desi[j]<-paste(as.character(names(tt1)),collapse="_")
                    } else {
                        desi[j]<-as.character(names(tt1))
                    }
                }
            }
            desi<-data.frame(desi)
            colnames(desi)<-metaVars
            rownames(desi)<-nn[i]
        }
        EXPRS[[i]]<-mat
        DESI[[i]]<-desi
        GENES[[i]]<-rownames(mat)
    }
    
 return(list(EXPRS,DESI,GENES))
}



#' Summarizes single cell expression profiles into pseudobulk expression profiles for Seurat or SingleCellExperiment
#'   objects.
#'
#' Summarizes single cell expression profiles into pseudobulk expression profiles for Seurat or SingleCellExperiment
#'   objects. This is the main function that calls either sc_seurat or sc_sce.
#' @param obj object. A Seurat or SingleCellExperiment object.
#' @param metaVars character. The column names from obj@meta.data or pData(obj) whose information will be stored
#'   in the design of the pseudobulk data.
#' @param replCol character. A column name from obj@meta.data or pData(obj) whose information is used to split and
#'   generate the pseudobulk. If NULL (default), all cells are summarized, otherwise only the cells of each replCol
#'   level will be summarized.
#' @param obj.type character. One of Seurat (default) or SCE.
#' @param common.genes logical. If TRUE, only the common genes across all libraries will be stored and analyzed. If FALSE
#'   (default) all genes will be stored. In libraries where the genes are not present their counts will be 0.
#' @keywords sc2pseudoBulk
#' @return A list with the pseudobulk expression profiles, the design and the genes of each Seurat or SingleCellExperiment
#'   object.
#'

sc2pseudoBulk<-function(obj,metaVars,replCol=NULL,obj.type="Seurat",common.genes=FALSE){
    
    nn<-names(obj)
    testname<-grep(":",nn)
    if(length(testname)>0){
        print("Character : is not permitted in names(obj). It will be replaced with _ !")
        nn<-gsub(":","_",nn)
    }
    
    
    if(is.null(nn)){
        stop("Provide the obj list component names")
    }
    if(length(!is.na(nn))!=length(nn)){
        stop("One or more obj list components do(es) not have a name!")
    }
    if(length(unique(nn))!=length(nn)){
        stop("The obj list components should have identical names!")
    }
    obj.type<-match.arg(obj.type,c("Seurat","SCE"))
    
    
    # seurat objects
    if(obj.type=="Seurat"){
        res<-sc_seurat(obj=obj,metaVars=metaVars,replCol=replCol)
    }
    if(obj.type=="SCE"){
        res<-sc_sce(obj=obj,metaVars=metaVars,replCol=replCol)
    }
    
    EXPRS<-res[[1]]
    DESI<-res[[2]]
    GENES<-res[[3]]
    V<-table(unlist(GENES))
    v<-names(V)[V==length(EXPRS)]
    
    if(common.genes){
        
        if(length(v)==0){
            print("There are no common genes across all objects. Setting common.genes = FALSE and summarizing...")
            v1<-names(V)
            for(i in 1:length(EXPRS)){
                mm<-match(v1,rownames(EXPRS[[i]]),nomatch=0)
                EXPRS1<-data.frame(EXPRS[[i]][mm,])
                rownames(EXPRS1)<-rownames(EXPRS[[i]])[mm]
                colnames(EXPRS1)<-colnames(EXPRS[[i]])
                if(length(mm[mm==0])>0){
                    EXPRS2<-matrix(0,length(mm[mm==0]),ncol(EXPRS[[i]]))
                    EXPRS2<-data.frame(EXPRS2)
                    rownames(EXPRS2)<-v1[mm==0]
                    colnames(EXPRS2)<-colnames(EXPRS[[i]])
                } else {
                    EXPRS2<-NULL
                }
                if(!is.null(EXPRS2)){
                    EXPRS[[i]]<-rbind(EXPRS1,EXPRS2)
                } else {
                    EXPRS[[i]]<-EXPRS1
                }
                EXPRS[[i]]<-EXPRS[[i]][sort.list(rownames(EXPRS[[i]])),]
            }
            
            
        } else {
            
            for(i in 1:length(EXPRS)){
                mm<-match(v,rownames(EXPRS[[i]]),nomatch=0)
                EXPRS[[i]]<-data.frame(EXPRS[[i]][mm,])
                rownames(EXPRS[[i]])<-v
                colnames(EXPRS[[i]])<-colnames(res[[1]][[i]])
            }
            
        }

    } else {
        
        v1<-names(V)
        for(i in 1:length(EXPRS)){
            mm<-match(v1,rownames(EXPRS[[i]]),nomatch=0)
            EXPRS1<-data.frame(EXPRS[[i]][mm,])
            rownames(EXPRS1)<-rownames(EXPRS[[i]])[mm]
            colnames(EXPRS1)<-colnames(EXPRS[[i]])
            if(length(mm[mm==0])>0){
                EXPRS2<-matrix(0,length(mm[mm==0]),ncol(EXPRS[[i]]))
                EXPRS2<-data.frame(EXPRS2)
                rownames(EXPRS2)<-v1[mm==0]
                colnames(EXPRS2)<-colnames(EXPRS[[i]])
            } else {
                EXPRS2<-NULL
            }
            if(!is.null(EXPRS2)){
                EXPRS[[i]]<-rbind(EXPRS1,EXPRS2)
            } else {
                EXPRS[[i]]<-EXPRS1
            }
            EXPRS[[i]]<-EXPRS[[i]][sort.list(rownames(EXPRS[[i]])),]
        }
        
    }
    
    exprs<-do.call(cbind,EXPRS)
    exprs<-data.frame(Gene=rownames(exprs),exprs)
    desi<-do.call(rbind,DESI)
    desi<-data.frame(sampleID=rownames(desi),desi)

 return(list(Exprs=exprs,Design=desi))
}





hb_fun<-function(data,genes.attr,conds,include.chemistry=TRUE,model.df=3){
    
    y<-data$Counts
    desi<-data$Design
    group<-factor(desi$Condition)
    y<-DGEList(counts=y,group=group)
   
    if(include.chemistry){
        ww<-1:nrow(desi)
    } else {
        ww<-which(desi$Chemistry=="V3")
    }
    y1<-y[,ww]
    desi1<-desi[ww,]
    r1<-rowSums(y1$counts[,y1$samples$group==conds[1]]>0)
    r2<-rowSums(y1$counts[,y1$samples$group==conds[2]]>0)
    keep1<-which(r1>0.3*ncol(y1$counts[,y1$samples$group==conds[1]]))
    keep2<-which(r2>0.3*ncol(y1$counts[,y1$samples$group==conds[2]]))
    keep<-unique(c(keep1,keep2))
    y1<-y1[keep,,keep.lib.sizes=FALSE]
    
    counts <- y1$counts
    meta <- desi1
    colnames(counts)<-rownames(meta)
    obj<-Lorgnette_from_matrix(counts=counts,metadata=meta,gene.attr=genes.attr,gene.id="Symbol")
    mm<-match(c("ND","PD","T2D"),conds,nomatch=0)
    if(length(mm[mm==0])>1){
        stop("Two conditions need to be specified")
    }
    if(length(mm[mm==0])==0){
        stop("Two conditions need to be specified")
    }
    if(mm[3]==0){
        obj<-LorgnetteSet(obj=obj,filters="Condition!='T2D'",pseudotime.by="HbA1c",state.by="Condition")
        full.mod<-paste("~sm.ns(Pseudotime,df=",model.df,") + Sex + BMI + Race + Age",sep="")
        red.mod<-"~Sex + BMI + Race + Age"
    }
    if(mm[2]==0){
        obj<-LorgnetteSet(obj=obj,filters="Condition!='PD'",pseudotime.by="HbA1c",state.by="Condition")
        full.mod<-paste("~sm.ns(Pseudotime,df=",model.df,") + Sex + BMI + Race + Age + Chemistry",sep="")
        red.mod<-"~Sex + BMI + Race + Age + Chemistry"
    }
    if(mm[1]==0){
        obj<-LorgnetteSet(obj=obj,filters="Condition!='ND'",pseudotime.by="HbA1c",state.by="Condition")
        full.mod<-paste("~sm.ns(Pseudotime,df=",model.df,") + Sex + BMI + Race + Age",sep="")
        red.mod<-"~Sex + BMI + Race + Age"
    }
    pData(obj)$BMI<-scale(as.numeric(as.character(pData(obj)$BMI)),scale=T)
    pData(obj)$Age<-scale(as.numeric(as.character(pData(obj)$Age)),scale=T)
    DE<-LorgnetteDE(obj=obj,num.exprs.samples=3,branch=NULL,branch2=NULL,
                    full.model=full.mod,
                    reduced.model=red.mod,
                    joinedDE=TRUE)

    real<-DE[[1]]
    real<-real[,-5]
    colnames(real)[5]<-"num_samples_expressed"
    
    f1<-factor(pData(obj)$State)
    f2<-factor(pData(obj)$Chemistry)
    f3<-factor(pData(obj)$Sex)
    f4<-factor(pData(obj)$Race)
    f5<-pData(obj)$BMI
    f6<-pData(obj)$Age
    y1<-obj@assayData$exprs
    y1<-DGEList(counts=y1,group=f1)
    
    conds<-sort(conds)
    if(include.chemistry){
        designC<-model.matrix(~0+f1+f2+f3+f4+f5+f6)
        colnames(designC)<-c(conds[1],conds[2],"V3","Male","Hispanic","White","BMI","Age")
    } else {
        designC<-model.matrix(~0+f1+f3+f4+f5+f6)
        colnames(designC)<-c(conds[1],conds[2],"Male","Hispanic","White","BMI","Age")
    }
    yy <- calcNormFactors(y1)
    yy<-estimateDisp(yy, designC)
    fit <- glmFit(yy, designC)
    fit1 <- glmLRT(fit, contrast=c(1,-1,rep(0,(ncol(designC)-2))))
    res1<-data.frame(topTags(fit1,n=nrow(fit1$counts)),Comp=paste(conds[1],"-",conds[2]))

    ii<-intersect(rownames(res1),real$Gene)
    res1<-res1[match(ii,rownames(res1)),]
    real<-real[match(ii,real$Gene),]
    real<-cbind(real,logFC=res1$logFC,logCPM=res1$logCPM,edgeR_PV=res1$PValue,edgeR_FDR=res1$FDR,Comparison=res1$Comp)
    
    dat<-obj@assayData$exprs
    relexpr<-t(t(dat) /  pData(obj)[, 'Size_Factor'])
    relexpr<-log(relexpr+1,10)
    
 return(list(DE=real,NormData=relexpr))
}



plot_hb<-function(data,genes.attr,conds,genes2plot,normalize.data,plot.formula="y~x"){
    
    mm<-match(c("ND","PD","T2D"),conds,nomatch=0)
    
    if(!normalize.data){
        relexpr<-data
        
    } else {
        
        y<-data$Counts
        desi<-data$Design
        group<-factor(desi$Condition)
        y<-DGEList(counts=y,group=group)
        keep<-as.list(rep(0,length(conds)))
        for(i in 1:length(conds)){
            rr<-rowSums(y$counts[,y$samples$group==conds[i]]>0)
            keep[[i]]<-which(rr>0.3*ncol(y$counts[,y$samples$group==conds[i]]))
        }
        keep<-unique(unlist(keep))
        y<-y[keep,,keep.lib.sizes=FALSE]
    
        counts <- y$counts
        meta <- desi
        colnames(counts)<-rownames(meta)
        obj<-Lorgnette_from_matrix(counts=counts,metadata=meta,gene.attr=genes.attr,gene.id="Symbol")
        
        if(length(mm[mm==0])>1){
            stop("At least two conditions need to be specified")
        }
        
        if(length(mm[mm==0])==0){
            obj<-LorgnetteSet(obj=obj,filters=NULL,pseudotime.by="HbA1c",state.by="Condition")
        }
        
        if(mm[3]==0){
            obj<-LorgnetteSet(obj=obj,filters="Condition!='T2D'",pseudotime.by="HbA1c",state.by="Condition")
        }
        if(mm[2]==0){
            obj<-LorgnetteSet(obj=obj,filters="Condition!='PD'",pseudotime.by="HbA1c",state.by="Condition")
        }
        if(mm[1]==0){
            obj<-LorgnetteSet(obj=obj,filters="Condition!='ND'",pseudotime.by="HbA1c",state.by="Condition")
        }
        
        dat<-obj@assayData$exprs
        relexpr<-t(t(dat) /  pData(obj)[, 'Size_Factor'])
        relexpr<-log(relexpr+1,10)
    }
    
    
    if(length(genes2plot)>1){
        relexpr<-relexpr[match(genes2plot,rownames(relexpr),nomatch=0),]
    } else {
        relexpr<-matrix(as.numeric(as.character(relexpr[match(genes2plot,rownames(relexpr),nomatch=0),])),nrow=1)
        rownames(relexpr)<-genes2plot
    }
    ss<-setdiff(genes2plot,rownames(relexpr))
    if(length(ss)>0){
        print(paste("The following genes have been removed from the analysis: ",paste(ss,collapse=", ")))
    }
    genes2plot<-intersect(genes2plot,rownames(relexpr))
    
    for(i in 1:length(genes2plot)){
        w<-which(rownames(relexpr)==genes2plot[i])
        e<-as.numeric(as.character(relexpr[w,]))
        dat<-data.frame(Expr=e,HbA1c=pData(obj)$Pseudotime,Condition=pData(obj)$State)
        p1<-ggplot(dat,aes(x=HbA1c,y=Expr))+geom_point(aes(color=Condition))+geom_smooth(data=dat,aes(x=HbA1c,y=Expr),method = lm, formula = formula(plot.formula)) + ylab(paste(genes2plot[i]," expression",sep=""))
        p2<-ggplot(dat, aes(x=Condition, y=Expr)) + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.2)) + ylab(paste(genes2plot[i]," expression",sep=""))
        grid.arrange(p1, p2, nrow = 2)
    }

  return(relexpr)
}



hb_fun2<-function(data,genes.attr,conds,x,include.chemistry=TRUE){
    
    y<-data$Counts
    desi<-data$Design
    group<-factor(desi$Condition)
    y<-DGEList(counts=y,group=group)
   
    if(include.chemistry){
        ww<-1:nrow(desi)
    } else {
        ww<-which(desi$Chemistry=="V3")
    }
    y1<-y[,ww]
    desi1<-desi[ww,]
    keep<-as.list(rep(0,length(conds)))
    for(i in 1:length(conds)){
        rr<-rowSums(y1$counts[,y1$samples$group==conds[i]]>0)
        keep[[i]]<-which(rr>0.3*ncol(y1$counts[,y1$samples$group==conds[i]]))
    }
    keep<-unique(unlist(keep))
    y1<-y1[keep,,keep.lib.sizes=FALSE]
    
    counts <- y1$counts
    meta <- desi1
    colnames(counts)<-rownames(meta)
    obj<-Lorgnette_from_matrix(counts=counts,metadata=meta,gene.attr=genes.attr,gene.id="Symbol")
    obj<-LorgnetteSet(obj=obj,filters="Condition!='T2D'",pseudotime.by=x,state.by="Condition")
    if(x=="BMI"){
        pData(obj)$Age<-scale(as.numeric(as.character(pData(obj)$Age)),scale=T)
        DE<-LorgnetteDE(obj=obj,num.exprs.samples=3,branch=NULL,branch2=NULL,
                        full.model="~sm.ns(Pseudotime,df=1) + Sex + Age + Race",
                        reduced.model="~Sex + Age + Race",
                        joinedDE=TRUE)
    }
    if(x=="Age"){
        pData(obj)$BMI<-scale(as.numeric(as.character(pData(obj)$BMI)),scale=T)
        DE<-LorgnetteDE(obj=obj,num.exprs.samples=3,branch=NULL,branch2=NULL,
                        full.model="~sm.ns(Pseudotime,df=1) + Sex + BMI + Race + State",
                        reduced.model="~Sex + BMI + Race + State",
                        joinedDE=TRUE)
    }
    
    de<-DE[[1]]
    dat<-obj@assayData$exprs
    relexpr<-t(t(dat) /  pData(obj)[, 'Size_Factor'])
    relexpr<-log(relexpr+1,2)

    cors<-matrix(0,nrow(relexpr),3)
    for(i in 1:nrow(relexpr)){
        rr<-cor.test(relexpr[i,],pData(obj)$Pseudotime)
        cors[i,]<-c(rownames(relexpr)[i],rr[[4]],rr$p.value)
    }
    cors<-cors[match(rownames(relexpr),cors[,1],nomatch=0),]
    de<-cbind(de,CorrCoef=cors[,2],CorrP=cors[,3],CorrFDR=p.adjust(as.numeric(cors[,3]),"BH"))

 return(de)
}


fisher_test<-function(data){
    
    res<-matrix(0,nrow(data),3)
    for(i in 1:nrow(data)){
        x<-matrix(c(data[i,2],data[i,3],sum(data[-i,2]),sum(data[-i,3])),nrow=2)
        ff<-fisher.test(x,alternative="greater")
        res[i,]<-c(data$Cluster[i],as.numeric(as.character(ff[[3]])),ff$p.value)
    }
    res<-data.frame(res)
    res$FDR<-p.adjust(as.numeric(as.character(res[,3])),"BH")
    colnames(res)<-c("Cluster","OR","Fisher.p","FDR")
    
 return(res)
}
