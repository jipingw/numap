//#include <Rcpp.h>
//#include <algorithm>
//#include <cmath>

#include <RcppArmadillo.h>
#include <algorithm>
//[[Rcpp::depends("RcppArmadillo")]]
using namespace Rcpp;
using namespace arma;

////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List segment(arma::vec w_vector, arma::vec c_vector, int gap_size=500,
             double thresh=0.01){
  int i=0;
  //scan for blocks with both watson and crick have 0's
  //gap_size is number of 0's in the region to define as segment point.default is 200
  //segment coordinates is 0 indexed.
  int n=w_vector.size();
  int zero_len=0;
  int j_s=0;
  int j_e=0;
  //start and for zero string>gap_size
  uvec start(1);
  uvec end(1);

  start(0)=0;
  end(0)=w_vector.size()-1;
  while (i<n){
    //Rcout<<i<<endl;
    if( w_vector(i) < thresh & c_vector(i)<thresh){
      //Rcout <<"i="<<i<<"\t";
      zero_len = zero_len+1;
      //Rcout<<"zero_len="<<zero_len<<endl;
      if(i==n-1 & zero_len>=gap_size){
        start.resize(j_s+1);
        end.resize(j_e+1);
        start(j_s)=i-zero_len+1;
        end(j_e)=i;
        //Rcout<<"i="<<i<<endl;
        //Rcout<<"start"<<start<<endl;
        //Rcout<<"end"<<end<<endl;
      }
    }else{
      if (zero_len>=gap_size){
        //Rcout<<i<<endl;
        //Rcout<<j_e<<endl;
        if(j_e==0){
          start(0)=i-zero_len;
          j_s=j_s+1;
          j_e=j_e+1;
          end(0)=i-1;
        //Rcout<<"start"<<start<<endl;
        //Rcout<<"end"<<end<<endl;
        }else{
          start.resize(j_s+1);
          end.resize(j_e+1);
          start(j_s)=i-zero_len;
          end(j_e)=i-1;
          j_s=j_s+1;
          j_e=j_e+1;
          //Rcout<<"start"<<start<<endl;
          //Rcout<<"end"<<end<<endl;
        }
      }
      zero_len=0;
    }
    i=i+1;
  }

  //Rcout<<"start"<<start<<endl;
  //Rcout<<"end"<<end<<endl;

  uvec seg_start=zeros<uvec>(start.size()+1);
  uvec seg_end=zeros<uvec>(start.size()+1);

  if(!(start(0)==0&&end(0)==w_vector.size()-1)){
    //Rcout<<"here"<<endl;

    if(start(0)==0){
      if(end(end.size()-1)==n-1){
        for(i=0;i<end.size()-1;i++){
          seg_start(i)=end(i)+1;
          seg_end(i)=start(i+1)-1;
        }
      }else{
        for(i=0;i<end.size()-1;i++){
          seg_start(i)=end(i)+1;
          seg_end(i)=start(i+1)-1;
        }
        seg_start(i)=end(i)+1;
        seg_end(i)=n-1;
      }
    }else{
      seg_start(0)=0;
      seg_end(0)=start(0)-1;
      if(end(end.size()-1)==n-1){
        for(i=0;i<end.size()-1;i++){
          seg_start(i+1)=end(i)+1;
          seg_end(i+1)=start(i+1)-1;
        }
        //Rcout<<"seg_start"<<seg_start<<endl;
        //Rcout<<"seg_end"<<seg_end<<endl;

      }else{
        for(i=0;i<end.size()-1;i++){
          seg_start(i+1)=end(i)+1;
          seg_end(i+1)=start(i+1)-1;
        }
        seg_start(i+1)=end(i)+1;
        seg_end(i+1)=n-1;
      }
    }
  }

  for(i=0;i<seg_start.size();i++){
    if(seg_end(i)-seg_start(i)<250){
     seg_start(i)=0;
      seg_end(i)=0;
    }
  }
  seg_start=seg_start(find(seg_end>0));
  seg_end=seg_end(find(seg_end>0));
  //Rcout<<"seg_start"<<seg_start<<endl;
  //Rcout<<"seg_end"<<seg_end<<endl;
  List seg;
  seg("seg_start")=seg_start;
  seg("seg_end")=seg_end;
  return (seg);
}


////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::uvec unique_M_one_seg(arma::vec a, int wsize, int seg_start){

  //seg_start is 0-based index
  //this function return 1-based nucleosome map
  //wsize is window size to define unique, currently set as 120 bp, that means +/-120

  int i;
  int map_size=round(a.size()/wsize)+1;
  uvec final=zeros<uvec>(map_size);
  //Rcout<<final.size()<<endl;
  //int from,end;
  //Rcout<<a<<endl;
  arma::uvec b=sort_index(a,"descend");
  // b is the order of index corresponding to small to larger
  // its not the rank of element in b
  //Rcout<<"b"<<b<<endl;
  i=0;
  //Rcout<<max(a(b))<<endl;
  //Rcout<<b.size()<<endl;

  while(b.size()>0){
    if(a(b(0))<0.1) break;
    final(i)=b(0)+1+seg_start;
    if(b(0) < wsize){
      b=b(find(b>b(0)+wsize));
    }else{
      b=b(find(b>b(0)+wsize||b<b(0)-wsize));
    }
    //Rcout<<max(a(b))<<endl;
    //Rcout<<b<<endl;
    i=i+1;
  }

  final=final(find(final>0));
  final=sort(final);
  return(final);
}


////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::uvec unique_M_one_seg_long(arma::vec a, int wsize, int seg_start){

  //seg_start is 0-based index
  //this function return 1-based nucleosome map
  //wsize is window size to define unique, currently set as 120 bp, that means +/-120
  int i;
  int map_size=round(a.size()/wsize)+1;
  int seg_size=round(a.size()/wsize/100)+1;
  uvec final=zeros<uvec>(map_size);
  arma::uvec b=sort_index(a,"descend");
  // b is the order of index corresponding to small to larger
  // its not the rank of element in b
  //Rcout<<"b="<<b.subvec(0,9)<<endl;
  i=0;
  int current_size=0;

  while(b.size()>0){
    if(i > seg_size) break;
    final(i)=b(0)+1+seg_start;
    //Rcout<<"final(i)="<<final(i)<<endl;
    int from=0;
    //Rcout<<"b(0)="<<b(0))-wsize<<endl;
    if (b(0)>wsize){
      from=b(0)-wsize;
    }
    int end=b(0)+wsize;
    if(end>a.size()-1){
      end=a.size()-1;
    }
    //Rcout<<from<<"+"<<end<<endl;
    if(b(0) < wsize){
      b=b(find(b>b(0)+wsize));
    }else{
      b=b(find(b>b(0)+wsize||b<b(0)-wsize));
    }
    //Rcout<<"b(0)="<<b(0)<<endl;
    a.subvec(from,end)=zeros<vec>(end-from+1);
    i=i+1;
    //Rcout<<"i"<<i<<endl;
    current_size=current_size+1;
  }
  //segment it
  List sub_seg=segment(a,a,241,.1);
  uvec sub_start=sub_seg[0];
  uvec sub_end=sub_seg[1];
  //Rcout<<"sub_seg size"<<sub_seg.size()<<endl;
  for(uword j=0;j<=sub_start.size()-1;j++){
    vec sub_NCP=a.subvec(sub_start(j),sub_end(j));
    uvec current_seg=unique_M_one_seg(sub_NCP, wsize, seg_start+sub_start(j));
    if(current_seg.size()>0){
      final.subvec(current_size,current_size+current_seg.size()-1)=current_seg;
      current_size=current_size+current_seg.size();
    }
  }
  final=final(find(final>0));
  final=sort(final);
  return(final);
}


////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::uvec unique_M2(arma::vec a,arma::uvec seg_start,arma::uvec seg_end, int wsize=120){

  //seg is zero_based,and start with coordinate 0, end with a.size()-1
  //wsize is window size to define unique, currently set as +/-120 bp
  int i;
  int map_size=round(a.size()/wsize);
  map_size=map_size+2*seg_start.size();
  //Rcout<<map_size<<endl;
  arma::uvec final=zeros<uvec>(map_size);
  int current_size=0;
  for(i=0;i<seg_start.size();i++){
    //Rcout<<"i="<<i<<seg_start(i)<<"\t"<<seg_end(i)<<endl;
    vec current_NCP=a.subvec(seg_start(i),seg_end(i));
    if(seg_end(i)-seg_start(i)<10000 ){
      uvec current_seg=unique_M_one_seg(current_NCP, wsize, seg_start(i));
      //Rcout<<"current_size="<<current_seg.size()<<endl;
      if (current_seg.size()>0){
       final.subvec(current_size,current_size+current_seg.size()-1)=current_seg;
       current_size=current_size+current_seg.size();
      }
    }else{
      uvec current_seg=unique_M_one_seg_long(current_NCP, wsize, seg_start(i));
      //uvec current_seg=unique_M_one_seg(current_NCP, wsize, seg_start(i));
      //Rcout<<"current_size="<<current_seg.size()<<endl;
      if (current_seg.size()>0){
       final.subvec(current_size,current_size+current_seg.size()-1)=current_seg;
       current_size=current_size+current_seg.size();
      }
    }
  }
  final=final(find(final>0));
  return(final);
}

/////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List NCP_cal_poisson(arma::vec w_vector,arma::vec c_vector, arma::vec temp1, int ite_num){

  int i,m,j,ite;
  double current_k1,current_k2;
  int n=w_vector.size();
  List output;
  vec delta=zeros<vec>(n);
  //delta is log ratio of total cleavages in +/-73 of two strands
  //Rcout<<"n="<<n<<std::endl;
  w_vector=w_vector+0.1;
  c_vector=c_vector+0.1;
  //noise term
  vec noise_w=0.5*ones<vec>(n);
  vec noise_c=0.5*ones<vec>(n);
  uvec pos(147);
  double median_w,median_c;
  vec nuc(147);

  for(i=73;i<n-73;i++){
    nuc=w_vector.subvec(i-73,i+73);
    delta(i)=log(sum(nuc)+1);
    //delta(i)=log(sum(nuc.subvec(73-12,73+12))+1);

    median_w=median(nuc);
    noise_w(i)=mean(nuc(find(nuc<= median_w)));
    nuc=c_vector.subvec(i-73,i+73);
    median_c=median(nuc);
    noise_c(i)=mean(nuc(find(nuc<= median_c)));
    delta(i)=delta(i)-log(sum(nuc)+1);
    //delta(i)=delta(i)-log(sum(nuc.subvec(73-12,73+12))+1);
  }

  //Rcout<<"noise done"<<endl;
  vec k1=zeros<vec>(n);
  vec k2=zeros<vec>(n);
  k1.subvec(73,(n-74))=sum(w_vector)/n*ones<vec>(n-146);
  k2.subvec(73,(n-74))=sum(c_vector)/n*ones<vec>(n-146);

  //Rcout<<k1<<endl;

  vec t1=zeros<vec>(n);
  vec t2=zeros<vec>(n);


  //template index
  arma::uvec D(8);
  //D << -2 << -1 << 0 << 1 << 4 << 5 << 6 <<7;
  D << 0 << 1 << 2 << 3 << 6 << 7 << 8 <<9;
  //Rcout<<D<<endl;

  //ite is big iteration for algorithm convergence

  for (ite=0;ite<ite_num;ite++){
    t1=noise_w;
    t2=noise_c;

    for(i=73;i<n-73;i++){
      if(i<81|i>n-82){
        for(j=0;j<8;j++){
          if(i+D(j)-2>=73 & i+D(j)-2<=n-74){
            t2(i)=t2(i)+k2(i+D(j)-2)*temp1(j);
          }
          if(i-D(j)+2>=73 & i-D(j)+2<=n-74){
            t1(i)=t1(i)+k1(i-D(j)+2)*temp1(j);
          }
        }
      }else{
        for(j=0;j<8;j++){
          t1(i)=t1(i)+k1(i-D(j)+2)*temp1(j);
          t2(i)=t2(i)+k2(i+D(j)-2)*temp1(j);
        }
      }
    }

    //Rcout<<"t done"<<endl;

    for(m=73;m<n-73;m++){
      //Rcout<<m<<endl;
      current_k1=0.;
      current_k2=0.;
      if(m<81|m>n-82){
        for(j=0;j<8;j++){
          //Rcout<<j<<endl;
          //D(j) should be D(j)+2
          if(m+D(j)-2>=73 & m+D(j)-2<=n-74 & t1(m+D(j)-2)>0){
            current_k1=current_k1+w_vector(m+D(j)-2)*temp1(j)/t1(m+D(j)-2);
          }
          if(m-D(j)+2>=73 & m-D(j)+2<=n-74 & t2(m-D(j)+2)>0){
            current_k2=current_k2+c_vector(m-D(j)+2)*temp1(j)/t2(m-D(j)+2);
          }
        }
      }else{
        for(j=0;j<8;j++){
          if(t1(m+D(j)-2)>0){
            current_k1=current_k1+w_vector(m+D(j)-2)*temp1(j)/t1(m+D(j)-2);
          }
          if(t2(m-D(j)+2)>0){
            current_k2=current_k2+c_vector(m-D(j)+2)*temp1(j)/t2(m-D(j)+2);
          }
        }
      }
      k1(m)=k1(m)*current_k1/sum(temp1);
      k2(m)=k2(m)*current_k2/sum(temp1);
    }
    //Rcout<<"here"<<endl;
  }
  output("k1")=k1;
  output("k2")=k2;
  output("noise_w")=noise_w;
  output("noise_c")=noise_c;
  output("delta")=delta;
  return (output);
}


/////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List NCP_cal_poisson_4(arma::vec seq1,arma::vec w_vector,arma::vec c_vector,
                       arma::mat temp4, int ite_num){
  //four-template NCP score calculation
  int i,m,j,ite;
  double current_k1,current_k2;
  int n=w_vector.size();
  List output;
  vec delta=zeros<vec>(n);
  uvec w_temp=zeros<uvec>(seq1.size());
  uvec c_temp=zeros<uvec>(seq1.size());

  //determine template
  for(i=73;i<seq1.size()-73;i++){
    uword index3= i-3;
    uword index4= i+3;
    if (seq1(index3)==1& seq1(index4)==4){
      w_temp(i)=0;
      c_temp(i)=0;
    }else if(seq1(index3)==1& seq1(index4)!=4){
      w_temp(i)=1;
      c_temp(i)=2;
    }else if(seq1(index3)!=1& seq1(index4)==4){
      w_temp(i)=2;
      c_temp(i)=1;
    }else{
      w_temp(i)=3;
      c_temp(i)=3;
    }
  }
  //Rcout<<"template data done!"<<endl;

  //Rcout<<"n="<<n<<std::endl;
  w_vector=w_vector+0.1;
  c_vector=c_vector+0.1;
  //noise term
  vec noise_w=0.5*ones<vec>(n);
  vec noise_c=0.5*ones<vec>(n);
  uvec pos(147);
  double median_w,median_c;
  vec nuc(147);

  for(i=73;i<n-73;i++){
    nuc=w_vector.subvec(i-73,i+73);
    delta(i)=log(sum(nuc)+1);
    //delta(i)=log(sum(nuc.subvec(73-12,73+12))+1);

    median_w=median(nuc);
    noise_w(i)=mean(nuc(find(nuc<= median_w)));
    nuc=c_vector.subvec(i-73,i+73);
    median_c=median(nuc);
    noise_c(i)=mean(nuc(find(nuc<= median_c)));
    delta(i)=delta(i)-log(sum(nuc)+1);
    //delta(i)=delta(i)-log(sum(nuc.subvec(73-12,73+12))+1);
  }

  //Rcout<<"noise done!"<<endl;

  vec k1=zeros<vec>(n);
  vec k2=zeros<vec>(n);
  k1.subvec(73,(n-74))=sum(w_vector)/n*ones<vec>(n-146);
  k2.subvec(73,(n-74))=sum(c_vector)/n*ones<vec>(n-146);

  //Rcout<<k1<<endl;

  vec t1=zeros<vec>(n);
  vec t2=zeros<vec>(n);

  //template index
  arma::uvec D(8);
  //D << -2 << -1 << 0 << 1 << 4 << 5 << 6 <<7;
  D << 0 << 1 << 2 << 3 << 6 << 7 << 8 <<9;
  //Rcout<<D<<endl;

  //ite is big iteration for algorithm convergence

  for (ite=0;ite<ite_num;ite++){
    t1=noise_w;
    t2=noise_c;

    for(i=73;i<n-73;i++){
      if(i<81|i>n-82){
        for(j=0;j<8;j++){
          if(i+D(j)-2>=73 & i+D(j)-2<=n-74){
            t2(i)=t2(i)+k2(i+D(j)-2)*temp4(c_temp(i+D(j)-2),j);
          }
          if(i-D(j)+2>=73 & i-D(j)+2<=n-74){
            t1(i)=t1(i)+k1(i-D(j)+2)*temp4(w_temp(i-D(j)+2),j);
          }
        }
      }else{
        for(j=0;j<8;j++){
          t1(i)=t1(i)+k1(i-D(j)+2)*temp4(w_temp(i-D(j)+2),j);
          t2(i)=t2(i)+k2(i+D(j)-2)*temp4(c_temp(i+D(j)-2),j);
        }
      }
    }

    //Rcout<<"t done"<<endl;

    for(m=73;m<n-73;m++){
      //Rcout<<m<<endl;
      current_k1=0.;
      current_k2=0.;
      if(m<81|m>n-82){
        for(j=0;j<8;j++){
          //Rcout<<j<<endl;
          //D(j) should be D(j)+2
          if(m+D(j)-2>=73 & m+D(j)-2<=n-74 & t1(m+D(j)-2)>0){
            current_k1=current_k1+w_vector(m+D(j)-2)*temp4(w_temp(m+D(j)-2),j)/t1(m+D(j)-2);
          }
          if(m-D(j)+2>=73 & m-D(j)+2<=n-74 & t2(m-D(j)+2)>0){
            current_k2=current_k2+c_vector(m-D(j)+2)*temp4(c_temp(m-D(j)+2),j)/t2(m-D(j)+2);
          }
        }
      }else{
        for(j=0;j<8;j++){
          if(t1(m+D(j)-2)>0){
            current_k1=current_k1+w_vector(m+D(j)-2)*temp4(w_temp(m+D(j)-2),j)/t1(m+D(j)-2);
          }
          if(t2(m-D(j)+2)>0){
            current_k2=current_k2+c_vector(m-D(j)+2)*temp4(c_temp(m-D(j)+2),j)/t2(m-D(j)+2);
          }
        }
      }
      k1(m)=k1(m)*current_k1/sum(temp4.row(w_temp(m)));
      k2(m)=k2(m)*current_k2/sum(temp4.row(c_temp(m)));
    }
    //Rcout<<"here"<<endl;
  }
  output("k1")=k1;
  output("k2")=k2;
  output("noise_w")=noise_w;
  output("noise_c")=noise_c;
  output("delta")=delta;
  return (output);
}



//train 4 template model


/////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
arma::mat train4(arma::vec seq1, arma::vec w_vector, arma::vec c_vector, arma::uvec umap){
  arma::mat temp4=zeros<mat>(4,41);
  rowvec w_current=zeros<rowvec>(41);
  rowvec c_current=zeros<rowvec>(41);
  //uvec D(8);
  //D << 0 << 1 << 2 << 3 << 6 << 7 << 8 <<9;

  uvec D=regspace<arma::uvec>(0,40);
  int i=0;
  int j=0;
  int n1=0;
  int n2=0;
  int n3=0;
  int n4=0;

  //Rcout<<"here"<<endl;
  //umap is 1 based, however seq is 0 based.

  umap=umap-ones<uvec>(umap.size());
  //Rcout<<umap(0)<<endl;
  //fours rows are: +A+T; +A-T;-A+T,-A-T
  for (i=0;i<umap.size();i++){
    //Rcout<<"i="<<i<<endl;
    //Rcout<<"umpa(i)="<<umap(i)<<endl;

    //for (j = 0;j<8;j++){
    for (j = 0;j<41;j++){

      //Rcout<<"j="<<j<<endl;
      //Rcout<<D(j)<<endl;
      //uword index1=umap(i)+D(j)-2;
      //uword index2=umap(i)-D(j)+2;
      uword index1=umap(i)+D(j)-20;
      uword index2=umap(i)-D(j)+20;

      //Rcout<<index1<<" "<<index<<<<endl;

      w_current(j)=w_vector(index1);
      c_current(j)=c_vector(index2);

    }

    //Rcout<<"here"<<endl;
    //Rcout<<seq1(umap(i)-3)<<endl;
    uword index3= umap(i)-3;
    uword index4= umap(i)+3;

    if (seq1(index3)==1){ //-3 is A
      if(seq1(index4)==4){ //+3 is T
        temp4.row(0)=temp4.row(0)+w_current+c_current;
        n1=n1+2;
      }else{//+3 is not T
        temp4.row(1)=temp4.row(1)+w_current; //A-T
        temp4.row(2)=temp4.row(2)+c_current; //-A+T
        n2=n2+1;
        n3=n3+1;
      }
    }else{
      if (seq1(index4)==4){ //-A+T
        temp4.row(2)=temp4.row(2)+w_current; //-A+T
        temp4.row(1)=temp4.row(1)+c_current; //-A-T
        n3=n3+1;
        n2=n2+1;
      }else{//-A-T
        temp4.row(3)=temp4.row(3)+w_current+c_current; //-A+T
        n4=n4+2;
      }
    }
  }
  temp4.row(0)=temp4.row(0)/n1;
  temp4.row(1)=temp4.row(1)/n2;
  temp4.row(2)=temp4.row(2)/n3;
  temp4.row(3)=temp4.row(3)/n4;
  return(temp4);
}

// [[Rcpp::export]]
arma::mat train42(arma::vec seq1, arma::vec w_vector, arma::vec c_vector, arma::uvec umap){
  arma::mat temp4=zeros<mat>(4,8);
  rowvec w_current=zeros<rowvec>(8);
  rowvec c_current=zeros<rowvec>(8);
  uvec D(8);
  int n1=0;
  int n2=0;
  int n3=0;
  int n4=0;

  D << 0 << 1 << 2 << 3 << 6 << 7 << 8 <<9;
  int i=0;
  int j=0;
  //Rcout<<"here"<<endl;
  //umap is 1 based, however seq is 0 based.

  umap=umap-ones<uvec>(umap.size());
  //shift one as umap is 1-based.
  //Rcout<<umap(0)<<endl;
  //fours rows are: +A+T; +A-T;-A+T,-A-T
  for (i=0;i<umap.size();i++){
    //Rcout<<"i="<<i<<endl;
    //Rcout<<"umpa(i)="<<umap(i)<<endl;

    for (j = 0;j<8;j++){
      //Rcout<<"j="<<j<<endl;
      //Rcout<<D(j)<<endl;
      uword index1=umap(i)+D(j)-2;
      uword index2=umap(i)-D(j)+2;
      //Rcout<<index1<<" "<<index<<<<endl;

      w_current(j)=w_vector(index1);
      c_current(j)=c_vector(index2);

    }

    //Rcout<<"here"<<endl;
    //Rcout<<seq1(umap(i)-3)<<endl;
    uword index3= umap(i)-3;
    uword index4= umap(i)+3;

    if (seq1(index3)==1){ //-3 is A
      if(seq1(index4)==4){ //+3 is T
        n1=n1+2;
        temp4.row(0)=temp4.row(0)+w_current+c_current;

      }else{//+3 is not T
        temp4.row(1)=temp4.row(1)+w_current; //A-T
        n2=n2+1;
        temp4.row(2)=temp4.row(2)+c_current; //-A+T
        n3=n3+1;
      }
    }else{
      if (seq1(index4)==4){ //-A+T
        temp4.row(2)=temp4.row(2)+w_current; //-A+T
        n3=n3+1;
        temp4.row(1)=temp4.row(1)+c_current; //-A-T
        n2=n2+1;
      }else{//-A-T
        temp4.row(3)=temp4.row(3)+w_current+c_current; //-A-T
        n4=n4+2;
      }
    }
  }
  temp4.row(0)=temp4.row(0)/n1;
  temp4.row(1)=temp4.row(1)/n2;
  temp4.row(2)=temp4.row(2)/n3;
  temp4.row(3)=temp4.row(3)/n4;

  return(temp4);
}


// [[Rcpp::export]]
arma::rowvec linker_cal(arma::uvec umap){
  rowvec linker=zeros<rowvec>(500);
  uword i;
  for (i=1;i<umap.size();i++){
    //Rcout<<i<<endl;
    //Rcout<<umap(i)<<"\t"<<umap(i-1)<<endl;
    if(umap(i)-umap(i-1)-147<500 & umap(i)-umap(i-1)-147>0){
      //uword index1=round(umap(i));
      //Rcout<<index1<<endl;
      //uword index2=round(umap(i-1));
      //Rcout<<index2<<endl;
      uword index=umap(i)-umap(i-1)-147;
      //Rcout<<index<<endl;
      linker(index-1)= linker(index-1)+1;
    }
  }
  return(linker);
}

// [[Rcpp::export]]
arma::rowvec auto_corr(arma::uvec rmap, int max_lag){
  rowvec cor=zeros<rowvec>(max_lag);
  uword i;
  uword lag;
  for (lag=1;lag<=max_lag;lag++){
    for (i=0;i<rmap.size()-lag;i++){
      if(rmap(i+lag)==1 & rmap(i)==1){
        cor(lag-1)=cor(lag-1)+1;
      }
    }
  }
  return(cor);
}

// [[Rcpp::export]]
arma::vec occu_cal(arma::uvec pos, arma::vec score, arma::vec weights, int chrom_len){
  vec occu=zeros<vec>(chrom_len);
  uword i;
  //weights should a vector of odd number length
  int margin=(weights.size()-1)/2;
  for (i=0;i<pos.size();i++){
    //Rcout<<i<<endl;
    occu.subvec(pos(i)-margin-1,pos(i)+margin-1)=
    occu.subvec(pos(i)-margin-1,pos(i)+margin-1)+
    score(i)*weights;
  }
  return(occu);
}

// [[Rcpp::export]]
arma::vec occu_cal_M_paired_center(arma::uvec start, arma::vec end, arma::vec weights,
                                   arma::vec score, int chrom_len){
  vec occu=zeros<vec>(chrom_len);
  uword i;
  //weights should a vector of odd number length
  int margin=(weights.size()-1)/2;
  for (i=0;i<start.size();i++){
    //Rcout<<i<<endl;
    int center=floor((start(i)+end(i))/2);
    if(center<chrom_len-73 & center>73){
      occu.subvec(center-margin-1,center+margin-1)=
      occu.subvec(center-margin-1,center+margin-1)+
      score(i)*weights/2;
    center=ceil((start(i)+end(i))/2);
    occu.subvec(center-margin-1,center+margin-1)=
      occu.subvec(center-margin-1,center+margin-1)+
      score(i)*weights/2;
    }
  }
  return(occu);
}

// [[Rcpp::export]]
arma::vec occu_cal_M_paired_unif(arma::uvec start, arma::vec end, int chrom_len){
  vec occu=zeros<vec>(chrom_len);
  uword i;
  //weights should a vector of odd number length
  for (i=0;i<start.size();i++){
    //Rcout<<i<<endl;
    occu.subvec(start(i)-1,end(i)-1)=
      occu.subvec(start(i)-1,end(i)-1)+1;
  }
  return(occu);
}

// [[Rcpp::export]]
arma::vec temp_update(arma::vec temp1, arma::vec w_vector, arma::vec c_vector){
  uword i,j,k;
  //weights should a vector of odd number length
  //weights at pos -2,-1,0,1,2,3,4,5,6,7, pos 2,3 put 0
  vec score=zeros<vec>(w_vector.size());
  vec temp_old=temp1;

  for (k=1;k<20;k++){
    //Rcout<<"k="<<k<<endl;
    temp_old=temp1;
    for (i=7;i<score.size()-7;i++){
      //Rcout<<i<<endl;
      score(i)=sum(dot(temp_old,w_vector.subvec(i-2,i+7)))+sum(dot(temp_old,reverse(c_vector.subvec(i-7,i+2))));
    }
      //call peaks
      double ave=mean(score);
      Rcout<<ave<<endl;
      int num=0;
      List sub_seg=segment(score,score,200,ave);
      //List sub_seg=segment(score,score,200,5);

      //Rcout<<"segment done"<<endl;
      uvec sub_start=sub_seg[0];
      uvec sub_end=sub_seg[1];
      //Rcout<<"segment size"<<sub_start.size()<<endl;

      uvec peak= unique_M2(score,sub_start,sub_end,120);
      //Rcout<<"unique done"<<endl;

      ave=mean(score(peak));
      for(j=0;j<=peak.size()-1;j++){
        //if(score(peak(j))>ave){
        if(score(peak(j))>ave){
          peak(j)=0;
        }
      }
      peak=peak(find(peak>0));
      //Rcout<<"peak done, peak size"<<peak.size()<<endl;

      temp1=0.*temp1;
      for(j=0;j<=peak.size()-1;j++){
        //peak is 1-based, w,c is 0 based
          //Rcout<<peak(j)<<"\t";
          temp1=temp1+w_vector.subvec(peak(j)-3,peak(j)+6)+reverse(c_vector.subvec(peak(j)-8,peak(j)+1));
      }
      temp1=temp1/2/peak.size();
      temp1(4)=0.;
      temp1(5)=0.;
      //Rcout<<"temp1"<<temp1<<endl;
      if(sum(abs(temp1-temp_old))<0.01*sum(temp1)){
        break;
      }
    }
    //update temp1
  return(temp1);
}


// [[Rcpp::export]]
arma::vec occu_TSS_cal(arma::vec occu, arma::uvec pos, arma::uvec strand, int region){
  vec occu_TSS=zeros<vec>(region*2+1);
  uword i;
  //plot TSS+/-500bp
  for (i=0;i<pos.size();i++){
    if(strand(i)==1){
//Rcout<<"strand +"<<endl;
       occu_TSS=occu_TSS+occu.subvec(pos(i)-region-1,pos(i)+region-1);
    }else{
        occu_TSS=occu_TSS+reverse(occu.subvec(pos(i)-region-1,pos(i)+region-1));
    }
  }
  return(occu_TSS);
}

// [[Rcpp::export]]
arma::vec AATT_cal(arma::vec seq1, arma::uvec pos, int region){
  vec AATT=zeros<vec>(region*2);
  uword i;
  for (i=0;i<pos.size();i++){
    if(pos(i)>region & pos(i)+region+1<seq1.size()){
      AATT=AATT+seq1.subvec(pos(i)-region-1,pos(i)+region-2)%seq1.subvec(pos(i)-region,pos(i)+region-1);
    }
  }
  return(AATT);
}
