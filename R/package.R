#################NCP_wrapper
NCP_wrapper1<-function(bam_file,paired=TRUE, result_path, chromosome_path,
                      temp1=NULL,center=TRUE,unique_map_wsize=120, wsize=73,parse_skip=FALSE,
                      chrom_to_include=NULL){

  #result_path: folder where the bam_file is located
  #chromosome_path: folder for fasta file for individual chromosomes
  #temp1 is the single-template file containing one row 8 positions
  #center=TRUE, to compute center-weighted occupancy. otherwise uniform
  #wsize=73: occupancy window size +/-73. sometimes we use 60
  #unique_map_wsize: window size to define unique-map. 120 means +/-120 to be local maximum
  #chrom_to_include only include a list of those chromosomes specified to process parsing and downstream analysis


  #######parse bam file
  setwd(result_path)
  if(parse_skip==FALSE|is.null(parse_skip)==TRUE){
    if(paired==TRUE){
      cat("\n","Parsing bam..","\n")
      parse_bam_paired(bam_file=bam_file, result_path=result_path,chrom_to_include=chrom_to_include)
    }else{
      cat("\n","Parsing bam..","\n")
      parse_bam_single(bam_file=bam_file, result_path=result_path,chrom_to_include=chrom_to_include)
    }
  }

  ########deconvolution
  cat("\n\n","Deconvolution..","\n")

  #### deconvolution
  ### chrom_to_include only include those chromosomes specified to process downstream analysis

  NCP1(result_path=result_path, chromosome_path=chromosome_path,
       temp1=temp1,unique_map_wsize=unique_map_wsize)

  #### occupancy
  cat("\n","Computing occupancy..","\n")
  redundant_map_list=list.files(result_path, pattern= paste("R.T1.*.txt$",sep=""),full.names = TRUE)
  chrom_length=read.table("chrom.sizes",header=FALSE)
  occupancy(redundant_map_list=redundant_map_list, result_path=result_path, T4=FALSE,chrom_length=chrom_length, center=center,wsize=73)
  gc(verbose=FALSE)

}

#################NCP_wrapper
NCP_wrapper4<-function(bam_file,paired=TRUE, result_path, chromosome_path,
                       temp4=NULL,center=TRUE,unique_map_wsize=120, wsize=73, parse_skip=FALSE,
                       chrom_to_include=NULL){

  #result_path: folder where the bam_file is located
  #chromosome_path: folder for fasta file for individual chromosomes
  #temp1 is the single-template file containing one row 8 positions
  #center=TRUE, to compute center-weighted occupancy. otherwise uniform
  #wsize=73: occupancy window size +/-73. sometimes we use 60
  #unique_map_wsize: window size to define unique-map. 120 means +/-120 to be local maximum

  #######parse bam file
  setwd(result_path)
  if(parse_skip==FALSE|is.null(parse_skip)==TRUE){
    if(paired==TRUE){
      cat("\n","Parsing bam..","\n")
      parse_bam_paired(bam_file=bam_file, result_path=result_path,chrom_to_include)
    }else{
    cat("\n","Parsing bam..","\n")
    parse_bam_single(bam_file=bam_file, result_path=result_path,chrom_to_include)
    }
  }
  ########deconvolution
  cat("\n\n","Deconvolution..","\n")
  #### deconvolution
  NCP4(result_path=result_path, chromosome_path=chromosome_path,
         temp4=temp4,unique_map_wsize=unique_map_wsize)
  #### occupancy
  cat("\n","Computing occupancy T4..","\n")
  redundant_map_list=list.files(result_path, pattern= paste("R.T4.*.txt$",sep=""),full.names = TRUE)
  chrom_length=read.table("chrom.sizes",header=FALSE)
  occupancy(redundant_map_list=redundant_map_list, result_path=result_path, T4=TRUE,chrom_length=chrom_length,
            center=center,wsize=73)
  gc(verbose=FALSE)
}

######parse_bam_paired
parse_bam_paired=function(bam_file,result_path,chrom_to_include){
  setwd(result_path)
  bf = BamFile(bam_file, asMates = TRUE, qnameSuffixStart = ".")
  gr = as(seqinfo(bf), "GRanges")
  chrom_length=seqlengths(gr)
  chrom=as.vector(unlist(runValue(seqnames(gr))))
  write.table(data.frame(chr=chrom,length=chrom_length),file="chrom.sizes",col.names=FALSE,
              row.names=FALSE,quote=FALSE)

  if(!file.exists(paste(bam_file,".bai",sep=""))){
    cat("Index", bam_file,"\n")
    indexBam(bam_file)
  }
  what=c("rname","pos","mpos","isize")

  ##one can specificy which chromosomes to parse
  if(is.null(chrom_to_include)==FALSE){
    chrom=chrom_to_include
  }

  for (i in 1:length(chrom)){
    cat("\t", chrom[i])
    ch=chrom[i]
    if (grepl("M", ch)==FALSE) { ##skip mitochondria

      #which <-gr[i]

      which=gr[seqnames(gr)==ch]
      param = ScanBamParam(flag=scanBamFlag(isPaired=TRUE),
                           what=what,which=which,tag=c("NH","HI"))
      bam <- scanBam(bam_file, param=param)

      hits=bam[[1]]$tag$NH
      bamout=data.table(chr=bam[[1]]$rname, pos1=bam[[1]]$pos, mate_pos=bam[[1]]$mpos,
                        frag_size=bam[[1]]$isize, hits=bam[[1]]$tag$NH)
      rm(bam)
      #bamout1=bamout[abs(frag_size)>=140 & abs(frag_size)<=500,]
      bamout1=bamout
      rm(bamout)
      if(nrow(bamout1)>0){
      bamout2=data.table(chr=bamout1$chr,pos=(bamout1$frag_size>0)*bamout1$pos1-1*(bamout1$frag_size<0)*(bamout1$mate_pos-bamout1$frag_size-1),
                         weight=1/bamout1$hits)
      rm(bamout1)


      #bam_file_prefix=strsplit(bam_file, split=".",fixed=TRUE)[[1]][1]
      #fwrite(bamout2,file=paste(bam_file_prefix,".paired.parsed.txt",sep=""),sep="\t")

      ###########
      ### convert parsed bam to cleavage

      colnames(bamout2)=c("chr","start","weight")
      bamout2=bamout2[weight>0.1,]
      W_start=bamout2[start>0,start]
      W_weight=bamout2[start>0,weight]
      C_end=bamout2[start<0,start]
      C_weight=bamout2[start<0,weight]
      C_end=C_end*(-1)
      rm(bamout2)
      W=IRanges(start=W_start,width=rep(1,length(W_start)))
      C=IRanges(start=C_end,width=rep(1,length(C_end)))
      n_W=max(W_start)
      n_C=max(C_end)
      cov_W=round(coverage(W, weight=W_weight),2)
      cov_C=round(coverage(C, weight=C_weight),2)

      watson=data.table(cbind(seq_len(n_W),as.vector(cov_W)))
      crick=data.table(cbind(seq_len(n_C),as.vector(cov_C)))
      colnames(watson)=c("pos","cleavage")
      colnames(crick)=c("pos","cleavage")

      #watson[,chr:=ch]
      #crick[,chr:=ch]

      rm(list=c("W","C","cov_W","cov_C","W_start","C_end"))
      #setcolorder(watson,c("chr","pos","cleavage"))
      #setcolorder(crick,c("chr","pos","cleavage"))

      setcolorder(watson,c("pos","cleavage"))
      setcolorder(crick, c("pos","cleavage"))

      watson=watson[cleavage>0.1,]
      crick=crick[cleavage>0.1,]
      setwd(result_path)
      fwrite(watson,file=paste(ch,".W.txt",sep=""),sep="\t")
      fwrite(crick,file=paste(ch,".C.txt",sep=""),sep="\t")
      }
    }
  }
}

#####################################3
###parse_bam_single file
parse_bam_single=function(bam_file,result_path,chrom_to_include){
  setwd(result_path)
  bf = BamFile(bam_file, asMates = TRUE, qnameSuffixStart = ".")
  gr = as(seqinfo(bf), "GRanges")
  chrom_length=seqlengths(gr)
  chrom=as.vector(unlist(runValue(seqnames(gr))))
  write.table(data.frame(chr=chrom,length=chrom_length),file="chrom.sizes",col.names=FALSE,row.names=FALSE,quote=FALSE)

  ##
  if(!file.exists(paste(bam_file,".bai",sep=""))){
    cat("Index", bam_file,"\n")
    indexBam(bam_file)
  }
  what=c("flag","rname","pos","cigar")

  ##one can specificy which chromosomes to parse
  if(is.null(chrom_to_include)==FALSE){
    chrom=chrom_to_include
  }

  for (i in 1:length(chrom)){
    cat("\t", chrom[i])
    ch=chrom[i]
    if (grepl("M", ch)==FALSE) { ##skip mitochondria

      #which <-gr[i]
      which=gr[seqnames(gr)==ch]
      #param = ScanBamParam(flag=scanBamFlag(isPaired=FALSE, isSecondaryAlignment=FALSE),
      #                 what=what,mapqFilter=255,which=which) ##255: unique map
      param = ScanBamParam(flag=scanBamFlag(isPaired=FALSE),
                           what=what,which=which,tag=c("NH","HI"))
      bam <- scanBam(bam_file, param=param)
      bamout=data.table(flag= bam[[1]]$flag,chr=bam[[1]]$rname, pos1=bam[[1]]$pos,
                        cigar=bam[[1]]$cigar,hits=bam[[1]]$tag$NH)
      if(nrow(bamout)>0){
      off=strsplit(bamout$cigar,"M|S|D|H|I|H|=|X|N")
      off_last=sapply(off,tail,1) #flag=16 reverse
      #if(sum(is.na(as.numeric(off_last)))>1){
      #  cat("i=",i,"\n")
      #}
      offset=(bamout$flag==16|bamout$flag==272)*(as.numeric(off_last)-1) #flag=16 reverse, position will be negative =0: forward
      pos=-(bamout$flag==16|bamout$flag==272)*(bamout$pos+offset)+
        (bamout$flag==0|bamout$flag==256)*(bamout$pos+offset)
      bamout2=data.table(chr=bamout$chr,pos=pos,weight=1/bamout$hits)

      #bam_file_prefix=strsplit(bam_file,".",fixed=TRUE)[[1]][1]
      #fwrite(bamout2,file=paste(bam_file_prefix,".single.parsed.txt",sep=""),sep="\t")

      rm(list=c("bam","bamout"))
      #cat("\t\t","summarizing cleavages..","\n")
      bamout2=bamout2[weight>0.1,]
      colnames(bamout2)=c("chr","start","weight")
      W_start=bamout2[start>0,start]
      W_weight=bamout2[start>0,weight]
      C_end=bamout2[start<0,start]
      C_end=C_end*(-1)
      C_weight=bamout2[start<0,weight]

      rm(bamout2)
      W=IRanges(start=W_start,width=rep(1,length(W_start)))
      C=IRanges(start=C_end,width=rep(1,length(C_end)))
      n_W=max(W_start)
      n_C=max(C_end)
      cov_W=round(coverage(W,weight=W_weight),2)
      cov_C=round(coverage(C,weight=C_weight),2)
      watson=data.table(cbind(seq_len(n_W),as.vector(cov_W)))
      crick=data.table(cbind(seq_len(n_C),as.vector(cov_C)))
      colnames(watson)=c("pos","cleavage")
      colnames(crick)=c("pos","cleavage")

      #watson[,chr:=ch]
      #crick[,chr:=ch]
      rm(list=c("W","C","cov_W","cov_C","W_start","C_end"))

      #setcolorder(watson,c("chr","pos","cleavage"))
      #setcolorder(crick,c("chr","pos","cleavage"))

      setcolorder(watson,c("pos","cleavage"))
      setcolorder(crick, c("pos","cleavage"))

      watson=watson[cleavage>0.1,]
      crick=crick[cleavage>0.1,]
      fwrite(watson,file=paste(ch,".W.txt",sep=""),sep="\t")
      fwrite(crick,file=paste(ch,".C.txt",sep=""),sep="\t")
      }
    }
  }
  gc(verbose=FALSE)
}

##single-template NCP calling and

NCP1<-function(result_path, chromosome_path, temp1=NULL,unique_map_wsize=120){
  setwd(result_path)
  if(is.null(temp1)){
    tempdata=single/max(single) #use stroed single-templated trained from yeast
  }else{
    tempdata=scan(temp1) #temp1  should a file containing one row of 8 numbers
    tempdata=tempdata/sum(tempdata)
  }

  ##read in cleavages files
  w_files=list.files(result_path, pattern= paste(".W.txt$",sep=""),full.names=TRUE)
  c_files=list.files(result_path, pattern= paste(".C.txt$",sep=""),full.names=TRUE)
  #bam_file_prefix=strsplit(w_files[1],".",fixed=TRUE)[[1]][1]

  for (j in 1:length(w_files)){
      #bam_file_prefix=strsplit(w_files[j],".",fixed=TRUE)[[1]][1]

      watson=fread(w_files[j])
      crick=fread(c_files[j])
      #ch=as.character(unique(watson[,1]))
      ch=tail(strsplit(w_files[j],split="/",fixed=TRUE)[[1]],1)
      ch=strsplit(ch,split=".",fixed=TRUE)[[1]][1]

      cat("\t", ch)

      ##deconvolution
      max_n=max(watson[,1],crick[,1])
      w_vector=rep(0,max_n)
      c_vector=rep(0,max_n)
      w_vector[watson$pos]=watson$cleavage
      c_vector[crick$pos]=crick$cleavage
      output=NCP_cal_poisson(w_vector, c_vector, tempdata, 20)
      k1=round(output$k1,2)
      k2=round(output$k2,2)
      noise_w=output$noise_w
      noise_c=output$noise_c
      delta=output$delta
      k=round(k1+k2,2)
      NCPratio=round(k/(noise_w+noise_c),2)
      outlier_thresh=2*sd(delta) ##if k1 and k2 ratio exceeds 2 standard dev. twice the max
      #outlier_thresh=sd(abs(delta)) ##if k1 and k2 ratio exceeds 2 standard dev. twice the max

      cNCPscore=2*pmax(k1,k2)*(abs(delta)>outlier_thresh)+k*(abs(delta)<=outlier_thresh)
      NCP=data.table(cbind(seq(1:max_n),k1,k2,k,
                           NCPratio,round(cNCPscore,2)))
      colnames(NCP)=c("Position","k1","k2","NCPscore","Ratio","cNCPscore")

      #NCP[,chr:=ch]

      #setcolorder(NCP,c(7,1,2,3,4,5,6))
      fwrite(NCP,file=paste("NCP.T1.",ch,".txt",sep=""),sep="\t")

      ##unique map >10% of k
      seg=segment(k,k,gap_size=500,thresh=mean(k))
      #seg=segment(w_vector,c_vector,gap_size=500,thresh=0.01)

      seg_start=as.vector(seg[[1]])
      seg_end=as.vector(seg[[2]])
      if(length(seg_start)==0){
        umap=unique_M_one_seg_long(k,120,0)
      }else{
        umap=unique_M2(k,seg_start,seg_end,wsize=unique_map_wsize)
      }
      threshold_10=quantile(k[umap],0.1)
      umap=umap[k[umap]>threshold_10]
      ##segment chromosome into pieces gap should contain at least 200 zeros.

      #NCP_umap=NCP[umap,1:7]
      NCP_umap=NCP[umap,1:6]

      rid=seq(1:length(k))[k>threshold_10]
      NCP_rmap=NCP[rid,1:6]

      fwrite(NCP_umap,file=paste("U.T1.",ch,".txt",sep=""),sep="\t")
      fwrite(NCP_rmap,file=paste("R.T1.",ch,".txt",sep=""),sep="\t")
      ###plot AATT by chromosomes
      pdf(file=paste("AATT.",ch,".T1.pdf",sep=""))
        umap_file=paste("U.T1.",ch,".txt",sep="")
        genfile=list.files(chromosome_path, pattern= paste("^",ch,"\\.fa$",sep=""),full.names=TRUE)
        plotAATT(chromosome_path=chromosome_path,result_path=result_path, center_list=umap_file,NCP_thresh=0, wsize=73)
      dev.off()
    }
}

###calculate gaussian weighted and uniformly weighted occupancy.
occupancy=function(redundant_map_list,result_path, chrom_length=NULL, T4=FALSE, center=TRUE,wsize=73){
  ##redundent_map_list is a list of redumant map files
  ## chrom_length should be length of chromosome corresponding to teach redundant map
  ## if chrom_length is missing, then the maximum coordinate +60 bp from each redundance map will be used as the chromosome length
  ## center=TRUE will calculate center weighted and center=false will calcualte uniform weighted
  ##wsize is 73 as default, we often use 60 in the paper.
  setwd(result_path)
  if(is.null(center)|center==TRUE){
    weights=dnorm(c(-wsize:wsize),0,20)
    weights=weights/dnorm(0,0,20)
  }else{
    weights=rep(1,2*wsize+1)
  }
  for (i in 1:length(redundant_map_list)){
    rmap=fread(redundant_map_list[i])
    #cat(redundant_map_list[i], "\n")
    #chr=unlist(rmap[1,1])
    #chr=strsplit(w_files[j],split=".",fixed=TRUE)[[1]][3]

    chr=tail(strsplit(redundant_map_list[i],split="/",fixed=TRUE)[[1]],1)
    chr=strsplit(chr,split=".",fixed=TRUE)[[1]][3]

    #pos=rmap[,2]
    pos=rmap[,1]

    #score=rmap[,7] #use cNCPscore corrected center NCP score
    score=rmap[,6] #use cNCPscore corrected center NCP score

    if(is.null(chrom_length)){
      n=max(pos)+wsize
    }else{
      ###########
      #n=chrom_length[i]
      n=chrom_length[chrom_length[,1]==chr,2]
    }
    occu=as.numeric(occu_cal(unlist(pos),unlist(score), weights,n))
    occu=round(occu,2)

    #chrom=as.character(unique(rmap[,1]))
    occu_out=data.table("pos"=c(1:n),"occu"=occu)

    ##define bed and wig formats
    bed=data.frame(chrom=chr,
                   start=c(0:(n-1)),
                   end=c(1:(n)),
                   score=as.numeric(occu))
    bed=bed[bed$score>0,]
    seqin=Seqinfo(seqnames=chr, seqlengths=n)
    bed2=makeGRangesFromDataFrame(bed,keep.extra.columns=TRUE, ignore.strand=TRUE,
                                  seqinfo=seqin,starts.in.df.are.0based=TRUE)
    ####output

    if(is.null(center)|center==TRUE){
      if(T4==FALSE|is.null(T4)==TRUE){
        fwrite(occu_out,file=paste("Occu.center.T1.",chr,".txt",sep=""),sep="\t")
        export.wig(bed2,con=paste("Occu.center.T1.",chr,".wig",sep=""))
      }else{
        fwrite(occu_out,file=paste("Occu.center.T4.",chr,".txt",sep=""),sep="\t")
        export.wig(bed2,con=paste("Occu.center.T4.",chr,".wig",sep=""))
      }
    }else{
      if(T4==FALSE|is.null(T4)==TRUE){
        fwrite(occu_out,file=paste("Occu.unif.T1.",chr,".txt",sep=""),sep="\t")
        export.wig(bed2,con=paste("Occu.unif.T1.",chr,".wig",sep=""))

      }else{
        fwrite(occu_out,file=paste("Occu.unif.T4.",chr,".txt",sep=""),sep="\t")
        export.wig(bed2,con=paste("Occu.unif.T4.",chr,".wig",sep=""))
      }
    }
  }
}

##4-template deconvolution
NCP4 <-function(result_path, chromosome_path, temp4=NULL,unique_map_wsize=120){
  setwd(result_path)
  if(is.null(temp4)){
    tempdata=four/max(four) #use stroed single-templated trained from yeast
  }else{
    tempdata=as.matrix(read.table(temp4)) #temp4 should contain 4 rows of 8 columsn
    tempdata=tempdata/max(tempdata)
  }

  ##read in cleavages files
  w_files=list.files(result_path, pattern= paste(".W.txt$",sep=""),full.names=TRUE)
  c_files=list.files(result_path, pattern= paste(".C.txt$",sep=""),full.names=TRUE)
  #bam_file_prefix=strsplit(w_files[1],".",fixed=TRUE)[[1]][1]

  for (j in 1:length(w_files)){

    watson=fread(w_files[j])
    crick=fread(c_files[j])
    #ch=as.character(unique(watson[,1]))
    #ch=strsplit(w_files[j],split=".",fixed=TRUE)[[1]][1]

    ch=tail(strsplit(w_files[j],split="/",fixed=TRUE)[[1]],1)
    ch=strsplit(ch,split=".",fixed=TRUE)[[1]][1]


    cat("\t", ch)

    chr=list.files(chromosome_path, pattern= paste("^",ch,"\\.fa$",sep=""),full.names = TRUE)
    chr_fasta=read.fasta(chr)

    seq=chr_fasta[[1]]
    seq[seq=="a"]=1
    seq[seq=="c"]=2
    seq[seq=="g"]=3
    seq[seq=="t"]=4
    seq[seq=="n"]=0
    seq1=as.numeric(seq)

    max_n=max(watson[,1],crick[,1])
    w_vector=rep(0,max_n)
    c_vector=rep(0,max_n)

    w_vector[watson$pos]=watson$cleavage
    c_vector[crick$pos]=crick$cleavage

    ##call Rcpp codes
    tempdata=as.matrix(tempdata)
    output=NCP_cal_poisson_4(seq1,w_vector, c_vector, tempdata, 20)

    k1=round(output$k1,2)
    k2=round(output$k2,2)
    noise_w=output$noise_w
    noise_c=output$noise_c
    delta=output$delta
    k=k1+k2
    NCPratio=round(k/(noise_w+noise_c),2)

    outlier_thresh=2*sd(abs(delta)) ##if k1 and k2 ratio exceeds 2 standard dev. twice the max
    cNCPscore=2*pmax(k1,k2)*(abs(delta)>outlier_thresh)+k*(abs(delta)<=outlier_thresh)
    NCP=data.table(cbind(seq(1:max_n),k1, k2,k,NCPratio),round(cNCPscore,2))
    colnames(NCP)=c("Position","k1","k2","NCPscore","Ratio","cNCPscore")

    #NCP[,chr:=ch]
    #setcolorder(NCP,c(7,1,2,3,4,5,6))

    fwrite(NCP,file=paste("NCP.T4.",ch,".txt",sep=""),sep="\t")

    ##unique map >10% of k
    #seg=segment(w_vector,c_vector,gap_size=500,thresh=.001)
    seg=segment(k,k,gap_size=500,thresh=mean(k))

    seg_start=as.vector(seg[[1]])
    seg_end=as.vector(seg[[2]])
    if(length(seg_start)==0){
      umap=unique_M_one_seg_long(k,120,0)
    }else{
      umap=unique_M2(k,seg_start,seg_end,wsize=unique_map_wsize)
    }
    threshold_10=quantile(k[umap],0.1)
    umap=umap[k[umap]>threshold_10]
    ##segment chromosome into pieces gap should contain at least 200 zeros.

    NCP_umap=NCP[umap,1:6]
    rid=seq(1:length(k))[k>threshold_10]
    NCP_rmap=NCP[rid,1:6]
    fwrite(NCP_umap,file=paste("U.T4.",ch,".txt",sep=""),sep="\t")
    fwrite(NCP_rmap,file=paste("R.T4.",ch,".txt",sep=""),sep="\t")

    ###plot AATT by chromosomes
    pdf(file=paste("AATT.",ch,".T4.pdf",sep=""))
    umap_file=paste("U.T4.",ch,".txt",sep="")
    genfile=list.files(chromosome_path, pattern=paste("^",ch,"\\.fa$",sep=""),full.names=TRUE)
    plotAATT(chromosome_path=chromosome_path,result_path=result_path,center_list=umap_file,NCP_thresh=0, wsize=73)
    dev.off()
  }
}

## train single-template profile
trainTEMP1=function(w_files,c_files,result_path){
  setwd(result_path)
  P0=c(0.125,0.25,0.125,0,0,0,0,0.125,0.25,0.125)

  #connect w_vector and c_vector to single, add 200 0 in beginning and end for each
  total_w=numeric() #contain concatenated w_vector
  total_c=numeric() #contain concatenated c_vector
  for(i in 1:length(w_files)){
    #cat("i=",i,"\t")
    watson=fread(w_files[i])
    crick=fread(c_files[i])

    #maxw=max(watson[,2])
    #maxc=max(crick[,2])

    maxw=max(watson[,1])
    maxc=max(crick[,1])

    n=max(maxw,maxc)
    cutw=numeric(n+200)
    cutc=numeric(n+200)
    #pos=unlist(watson[,2])
    #cut=unlist(watson[,3])

    pos=unlist(watson[,1])
    cut=unlist(watson[,2])

    cutw[pos]=cut
    #pos=unlist(crick[,2])
    #cut=unlist(crick[,3])

    pos=unlist(crick[,1])
    cut=unlist(crick[,2])

    cutc[pos]=cut
    if(i==1){
      total_w=unlist(cutw)
      total_c=unlist(cutc)
    }else{
      total_w=c(total_w,cutw)
      total_c=c(total_c,cutc)
    }
  }
  temp1=temp_update(P0, total_w, total_c)

  temp1=as.vector(temp1[-c(5:6)])
  temp1=temp1/sum(temp1)
  write.table(t(temp1),file="temp1.txt",row.names=FALSE,col.names=FALSE)
return(temp1)
}

## train single-template profile based on masked genome
trainTEMP1_masked=function(w_files,c_files,result_path,chromosome_masked_path){
  setwd(result_path)
  P0=c(0.125,0.125,0.125,0.125,0,0,0.125,0.125,0.125,0.125)
  for (i in 1:length(w_files)){
  watson=fread(w_files[i])
  ch=unique(watson[,1])
  crick=fread(c_files[i])
  chr=list.files(chromosome_masked_path, pattern= paste("^",ch,"\\.fa$",sep=""),full.names = TRUE)
  chr_fasta=read.fasta(chr,forceDNAtolower=FALSE)
  seq=chr_fasta[[1]]
  seq[seq=="a"]=0
  seq[seq=="c"]=0
  seq[seq=="g"]=0
  seq[seq=="t"]=0
  seq[seq=="n"]=0
  seq[seq=="A"]=1
  seq[seq=="C"]=2
  seq[seq=="G"]=3
  seq[seq=="T"]=4
  seq[seq=="N"]=0
  seq=as.numeric(seq)
  n=length(seq)
  cutw=numeric(n)
  cutc=numeric(n)
  pos=unlist(watson[,1])
  cut=unlist(watson[,2])
  cutw[pos]=cut
  cutw[1:n]=cutw[1:n]*(seq!=0)
  pos=unlist(crick[,1])
  cut=unlist(crick[,2])
  cutc[pos]=cut
  cutc[1:n]=cutc[1:n]*(seq!=0)
  rm(list=c("seq","pos","cut","chr_fasta"))
  gc()
  if(i==1){
    total_w=unlist(cutw)
    total_c=unlist(cutc)
  }else{
    total_w=c(total_w,cutw)
    total_c=c(total_c,cutc)
  }
}
  temp1=temp_update(P0, total_w, total_c)
  temp1=as.vector(temp1[-c(5:6)])
  temp1=temp1/sum(temp1)
  write.table(t(temp1),file="temp1.txt",row.names=FALSE,col.names=FALSE)
  return(temp1)
}

##AATT plot
plotAATT = function(chromosome_path,center_list,result_path,NCP_thresh=0, wsize=73){
  setwd(result_path)
  AATT=as.vector(numeric(2*wsize))
  total=0
  for (center in center_list){
    umap=fread(center)
    umap=umap[umap$NCPscore >NCP_thresh,]

    #chr=strsplit(center,split=".",fixed=TRUE)[[1]][3]
    ch=strsplit(center,split=".",fixed=TRUE)[[1]][3]

    #ch=unique(umap[,1])
    pos=unlist(umap[,1])

    total=total+length(pos)
    chr=list.files(chromosome_path, pattern= paste("^",ch,"\\.fa$",sep=""),full.names = TRUE)
    chr_fasta=read.fasta(chr)
    seq=chr_fasta[[1]]
    seq[seq=="a"]=1
    seq[seq=="c"]=0
    seq[seq=="g"]=0
    seq[seq=="t"]=1
    seq[seq=="n"]=0
    seq=as.numeric(seq)
    AATT=AATT+AATT_cal(seq,pos,wsize)
  }
  AATT=AATT/total
  #write.table(AATT,file="AATT.freq.txt",row.names=FALSE,col.names=FALSE)
  plot(c(-wsize:(wsize-1)),AATT,type="l",xlab="distance to nucleosome center", ylab="frequency",main="AA/TT/TA/AT frequency")
  return(AATT)
}

## occupancy at TSS
TSS_occu=function(TSS_file, occu_list, wsize=1000){
  ## TSS_file must tab or space_delimited file and contain three columns, chr, strand and start (start is TSS start), 1 for watson and 2 for crick strand
  ## chr must match the chromosome names in Occupancy file
  ## occu_file must have same format from output of occupancy function, first column is chr, and last column is occupancy
  TSS=read.table(TSS_file,header=TRUE)
  chr=TSS[,1]
  strand=TSS[,2]
  start=TSS[,3]

  TSS_occu=numeric(2*wsize+1)
  for ( i in 1:length(occu_list)){
    cat(i,"\t")
    occu=fread(occu_list[i])

    #ch=as.character(unique(occu[,1]))
    #ch=strsplit(occu_list[i],split=".",fixed=TRUE)[[1]][4]

    ch=tail(strsplit(occu_list[i],split="/",fixed=TRUE)[[1]],1)
    ch=strsplit(ch,split=".",fixed=TRUE)[[1]][4]

    score=as.numeric(unlist(occu[,2]))

    strand_ch=strand[which(chr==ch)]
    start_ch=start[which(chr==ch)]
    TSS_occu=TSS_occu+occu_TSS_cal(score, start_ch, strand_ch, 1000)
  }
  plot(TSS_occu/mean(TSS_occu),type="l",ylab="normalized occupancy", xlab="distance to TSS")
  return(TSS_occu)
}

###MNase map
######parse_bam_paired_Mnase
###MNase map
######parse_bam_paired_Mnase

###MNase map
######parse_bam_paired_Mnase
parse_bam_paired_Mnase=function(bam_file,result_path,chrom_to_include=NULL){
  setwd(result_path)
  bf = BamFile(bam_file, asMates = TRUE, qnameSuffixStart = ".")
  gr = as(seqinfo(bf), "GRanges")
  chrom_length=seqlengths(gr)
  chrom=as.vector(unlist(runValue(seqnames(gr))))
  #cat(chrom,"\n")
  #cat(chrom_length)
  write.table(data.frame(chr=chrom,length=chrom_length),file="chrom.sizes",col.names=FALSE,row.names=FALSE,quote=FALSE)

  if(!file.exists(paste(bam_file,".bai",sep=""))){
    cat("Index", bam_file,"\n")
    indexBam(bam_file)
  }
  what=c("rname","pos","mpos","isize")

  if(is.null(chrom_to_include)==FALSE){
    chrom=chrom_to_include
  }

  for (i in 1:length(chrom)){
    cat("\t", chrom[i])
    ch=chrom[i]
    if (grepl("M", ch)==FALSE) { ##skip mitochondria
      which <-gr[i]
      #param = ScanBamParam(flag=scanBamFlag(isPaired=TRUE, isSecondaryAlignment=FALSE),
      #                     what=what,mapqFilter=255,which=which)
      param = ScanBamParam(flag=scanBamFlag(isPaired=TRUE),
                           what=what,which=which,tag=c("NH","HI"))
      bam <- scanBam(bam_file, param=param)

      bamout=data.table(chr=bam[[1]]$rname, pos1=bam[[1]]$pos, mate_pos=bam[[1]]$mpos,
                        frag_size=bam[[1]]$isize, hits=bam[[1]]$tag$NH)
      rm(bam)
      bamout1=bamout[bamout$frag_size>0,]

      if(nrow(bamout1)>0){

      #bamout2=data.table(chr=bamout1$chr,start=bamout1$pos1, end=bamout1$pos1+bamout1$frag_size-1,
      #                   weight=1/bamout1$hits)
      bamout2=data.table(start=bamout1$pos1, end=bamout1$pos1+bamout1$frag_size-1,
                         weight=1/bamout1$hits)
      bamout2$weight=round(bamout2$weight,2)

      rm(bamout1)

      ###########
      ### convert parsed bam to cleavage
      fwrite(bamout2,file=paste(ch,".reads.txt",sep=""),sep="\t")
      }
    }
  }
}

#################
##occupancy_Mnase
occupancy_Mnase_paired=function(reads_file_list, result_path, chrom_length=NULL, center=TRUE,
                                wsize=73,insert_length_low,insert_length_high,
                                maxhits=Inf,nameExt=NULL){
  ##reads_file is a list of paired-end reads file
  ## chrom_length should be length of chromosome corresponding to teach redundant map
  ## if chrom_length is missing, then the maximum coordinate +60 bp from each redundance map will be used as the chromosome length
  ## center=TRUE will calculate center weighted and center=false will calcualte uniform weighted
  ##wsize is 73 as default, we often use 60 in the paper.

  setwd(result_path)
  if(is.null(center)|center==TRUE){
    weights=dnorm(c(-wsize:wsize),0,20)
    weights=weights/dnorm(0,0,20)
  }

  for (i in 1:length(reads_file_list)){
    cat(reads_file_list[i],"\n")
    rmap=fread(reads_file_list[i])
    rmap=rmap[rmap$weight>=1/maxhits,]

    start=rmap$start
    end=rmap$end
    score=rmap$weight
    insert_size=abs(end-start)
    keeper=(insert_size>insert_length_low&insert_size<insert_length_high)
    start=start[keeper]
    end=end[keeper]
    score=score[keeper]

    #chr=strsplit(reads_file_list[i],"\\.reads.txt")[[1]][1]

    chr=tail(strsplit(reads_file_list[i],split="/",fixed=TRUE)[[1]],1)
    chr=strsplit(chr,split=".",fixed=TRUE)[[1]][1]


    if(is.null(chrom_length)){
      n=max(end)+wsize
    }else{
      #chr=strsplit(reads_file_list[i],"\\.reads.txt")[[1]][1]
      n=chrom_length[chrom_length[,1]==chr,2]
    }

    if(is.null(center)|center==TRUE){
      occu=occu_cal_M_paired_center(start=unlist(start),end=unlist(end),weights=weights,score=unlist(score),chrom_len=n)
      occu=round(occu,2)

      #chrom=as.character(unique(rmap[,1]))
      #occu_out=data.table(chr=rep(chrom,n),pos=c(1:n),occu=as.numeric(occu))
      occu_out=data.table(pos=c(1:n),occu=as.numeric(occu))

      #fwrite(occu_out,file=paste("Moccu.center.",chrom,".",nameExt,".txt",sep=""),sep="\t")

      fwrite(occu_out,file=paste("Moccu.center.",chr,".",nameExt,".txt",sep=""),sep="\t")

      bed=data.frame(chrom=chr,
                     start=c(0:(n-1)),
                     end=c(1:(n)),
                     score=as.numeric(occu))
      bed=bed[bed$score>0,]
      seqin=Seqinfo(seqnames=chrom, seqlengths=n)
      bed2=makeGRangesFromDataFrame(bed,keep.extra.columns=TRUE, ignore.strand=TRUE,
                                   seqinfo=seqin,starts.in.df.are.0based=TRUE)
      export.bw(bed2,con=paste("Moccu.center.",chr,".",nameExt,".bw",sep=""))
    }else{
      occu=occu_cal_M_paired_unif(start=unlist(start),end=unlist(end),chrom_len=n)
      occu=round(occu,2)

      chrom=as.character(unique(rmap[,1]))
      occu_out=data.table(chr=rep(chrom,n),pos=c(1:n),occu=as.numeric(occu))
      fwrite(occu_out,file=paste("Moccu.unif.",chr,".",nameExt,".txt",sep=""),sep="\t")
      bed=data.frame(chrom=chrom,
                     start=c(0:(n-1)),
                     end=c(1:(n)),
                     score=as.numeric(occu))
      bed=bed[bed$score>0,]
      seqin=Seqinfo(seqnames=chrom, seqlengths=n)
      bed2=makeGRangesFromDataFrame(bed,keep.extra.columns=TRUE, ignore.strand=TRUE,
                                    seqinfo=seqin,starts.in.df.are.0based=TRUE)
      export.bw(bed2,con=paste("Moccu.unif.",chr,".",nameExt,".bw",sep=""))
    }
  }
}

#################################
##occupancy_Mnase multiple core
#################################
occupancy_Mnase_paired_par=function(reads_file_list, result_path, chrom_length=NULL,
                                center=TRUE, wsize=73,insert_length_low,insert_length_high,
                                cores=NULL, maxhits=Inf, nameExt=NULL){
  ##reads_file is a list of paired-end reads file
  ## chrom_length should be length of chromosome corresponding to teach redundant map
  ## if chrom_length is missing, then the maximum coordinate +60 bp from each redundance map will be used as the chromosome length
  ## center=TRUE will calculate center weighted and center=false will calcualte uniform weighted
  ##wsize is 73 as default, we often use 60 in the paper.

  if(is.null(cores)==TRUE) cores=2

  max.cores=detectCores(logical = FALSE)
  if(max.cores-1<cores) cores <- max.cores-1

  cl <- makeCluster(cores)
  registerDoParallel(cl)

  setwd(result_path)
  if(is.null(center)|center==TRUE){
    weights=dnorm(c(-wsize:wsize),0,20)
    weights=weights/dnorm(0,0,20)
  }

  foreach(i = 1:length(reads_file_list),
          .export=c("occu_cal_M_paired_center","occu_cal_M_paired_unif"),
          .packages = c("data.table","GenomicRanges","rtracklayer","GenomeInfoDb"),
          .inorder=TRUE,.multicombine=FALSE) %dopar%
    {
      cat(reads_file_list[i],"\n")
      rmap=fread(reads_file_list[i])
      rmap=rmap[rmap$weight>=1/maxhits,]
      start=rmap$start
      end=rmap$end
      score=rmap$weight
      insert_size=abs(end-start)
      keeper=(insert_size>=insert_length_low & insert_size<=insert_length_high)
      start=start[keeper]
      end=end[keeper]
      score=score[keeper]

      chr=strsplit(reads_file_list[i],"\\.reads.txt")[[1]][1]

      if(is.null(chrom_length)){
        n=max(end)+wsize
      }else{
        #chr=strsplit(reads_file_list[i],"\\.reads.txt")[[1]][1]
        n=chrom_length[chrom_length[,1]==chr,2]
      }
      if(is.null(center)|center==TRUE){
        occu=occu_cal_M_paired_center(start=unlist(start),end=unlist(end),weights=weights,
                                      score=unlist(score),chrom_len=n)
        occu=round(occu,2)
        chrom=chr

        ##### added
        #chrom=paste("chr",chrom,sep="")
        ####

        occu_out=data.table(pos=c(1:n),occu=as.numeric(occu))
        fwrite(occu_out,file=paste("Moccu.center.",chrom,".",nameExt,".txt",sep=""),sep="\t")

        bed=data.frame(chrom=chrom,
                       start=c(0:(n-1)),
                       end=c(1:(n)),
                       score=as.numeric(occu))
        bed=bed[bed$score>0,]
        seqin=Seqinfo(seqnames=chrom, seqlengths=n)
        bed2=makeGRangesFromDataFrame(bed,keep.extra.columns=TRUE, ignore.strand=TRUE,
                                      seqinfo=seqin,starts.in.df.are.0based=TRUE)
        rm(bed)
        export.bw(bed2,con=paste("Moccu.center.",chrom,".",nameExt,".bw",sep=""))

      }else{
        occu=occu_cal_M_paired_unif(start=unlist(start),end=unlist(end),chrom_len=n)
        occu=round(occu,2)
        chrom=as.character(unique(rmap[,1]))

        ##### added
        #chrom=paste("chr",chrom,sep="")
        ####

        occu_out=data.table(pos=c(1:n),occu=as.numeric(occu))
        fwrite(occu_out,file=paste("Moccu.unif.",chrom,".",nameExt,".txt",sep=""),sep="\t")
        bed=data.frame(chrom=chrom,
                       start=c(0:(n-1)),
                       end=c(1:(n)),
                       score=as.numeric(occu))
        bed=bed[bed$score>0,]
        seqin=Seqinfo(seqnames=chrom, seqlengths=n)
        bed2=makeGRangesFromDataFrame(bed,keep.extra.columns=TRUE, ignore.strand=TRUE,
                                      seqinfo=seqin,starts.in.df.are.0based=TRUE)
        export.bw(bed2,con=paste("Moccu.unif.",chrom,".",nameExt,".bw",sep=""))
      }
    }
  stopImplicitCluster()
  stopCluster(cl)
}

