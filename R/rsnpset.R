rsnpset<-function(Y,delta=NULL,G,X=NULL,snp.sets,
                  score=c("cox", "binomial", "gaussian"),B=0, r.method = "permutation",
                  v.method="empirical",v.permute=TRUE,ret.rank=FALSE,
                  pinv.check=FALSE,pinv.method="specdecomp", pinv.tol=7.8e-8) {              
    score <- match.arg(score)
    
    Y    <-as.numeric(Y)
    delta<-as.numeric(delta)
 
    if(any(is.na(Y))) {
        stop("Y should be numeric")
    }
    if(any(is.na(delta))) {
        stop("delta should be numeric")
    }

    if(score=="cox") {
        if(!all(Y>=0)) {
            stop("Y should be greater than or equal to 0")
        }
        if(!all(delta %in% c(0,1))) {
            stop("delta should be 0 or 1")
        }
    }
    else if(score=="binomial") {
        if(length(unique(Y))==2) {
            Y<-as.numeric(c(0,1)[factor(Y)])
        }
        else {
            stop("Y should have two levels.")
        }
    }
    
    allGeneIDs <- colnames(G)
    geneIndexSets<-vector("list",length(snp.sets))
    names(geneIndexSets)<-names(snp.sets)
    for(i in 1:length(snp.sets)) {
       geneIndexSets[[i]]<-na.omit(fmatch(snp.sets[[i]],allGeneIDs))
    }
    geneIndexSets<-geneIndexSets[lapply(geneIndexSets,length)>0]

    if (is.null(getDoParName())) {
        registerDoSEQ()
    }
    result<-foreach(i=0:B, .packages="RSNPset") %dorng%{ 
        .Call("rsnpsetRcpp",Y,delta,G,geneIndexSets,i,score,v.method,pinv.check,pinv.tol,PACKAGE="RSNPset")
    }
    
    attributes(result) <- NULL
    
    names(result)<-paste0("Permutation.",0:(length(result)-1))
    names(result)[1]<-"Observed"
    
    class(result)<-"RSNPset"
    attr(result,"n")<-length(Y)
    attr(result,"KSub")<-length(snp.sets)
    attr(result,"KAna")<-length(geneIndexSets)
    attr(result,"mSub")<-sapply(snp.sets, length)
    attr(result,"mAna")<-sapply(geneIndexSets, length)
    attr(result,"B")<-B
    attr(result,"ret.rank")<-ret.rank
    attr(result,"v.permute")<-v.permute
    attr(result,"pinv.tol")<-pinv.tol
    
    if(pinv.check==FALSE) {
        attr(result,"pinv.check")<-NA
    }
    else {
        attr(result,"pinv.check")<- lapply(result,function(x) x[,3:7])
    }
    
    for(i in 1:(B+1)) {
        rownames(result[[i]])<-names(geneIndexSets)
        
        if(ret.rank==FALSE && i>1) {
            result[[i]]<-result[[i]][,1,drop=FALSE]
        }
        else {
            result[[i]]<-result[[i]][,1:2]
        }
    }
    result[[1]]["m"]<-sapply(geneIndexSets, length)
    return( result )
}
