#include <Rcpp.h>
using namespace Rcpp;
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

// [[Rcpp::export]]
void printProgress (double percentage)
{
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush (stdout);
}

// [[Rcpp::export]]
int one_cat_zero_begin(NumericVector probs){
  int k = probs.size();
  int res = 0;
  double z=sum(probs);
  if(z!=1){
    for(int i=0;i<k;++i){
      probs[i]=probs[i]/z;
    }
  }
  IntegerVector ans(k);
  rmultinom(1, probs.begin(), k, ans.begin());
  for(int j=0;j<k;++j){
    if(ans[j]==1){
      res=j;
      break;
    }
  }
  return res;
}

input: facts(vector) claims fcidx(vector) ctsc beta alpha0 alpha1

for( f in fid){
    
    cf <- claims[claims$fid==f,]
    ## get startidx and endidx of f in claim
    
    for(c in cf){
        calculate p with ctsc
    }

    update facts
    for(c in cf){
        update ctsc
    }
}

// [[Rcpp::export]]
List truthfinding(IntegerVector facts,IntegerVector fcidx, IntegerMatrix claims, IntegerMatrix ctsc, NumericMatrix beta, NumericMatrix alpha0, NumericMatrix alpha1){
  int nfacts=facts.size();
  int startidx=0,endidx=0;	// claim start and end index for each fact
  NumericVector probs(2);
  double conditional_claim=1;
  
  // main gibbs loop
  // for each fact f...
  for(int f = 0; f < nfacts; ++f){
    startidx=fcidx[f];
    endidx=fcidx[f+1];
    // generate sample probabilities
    for(int t = 0; t < 2; ++t){
      // reset conditional claim probability
      conditional_claim=1;
      // for each claim c in fact f...
      for(int c = startidx; c < endidx; ++c){
	conditional_claim = conditional_claim*
      }
      probs[t] = beta(f,t)*conditional_claim
    }
    // sample and update facts
    one_cat_zero_begin(probs);

    
    // update ctsc
    for(int c = startidx; c < endidx; ++c){
      
    }
    
  }
  
  
}


// [[Rcpp::export]]
List truthfinding(IntegerVector facts, IntegerMatrix claims, IntegerMatrix ctsc, int maxit, int sumL, int K,int V, IntegerVector L,double a, double g){
  int doc=0,word=0,topic=0,idx=0;
  NumericVector probs(K);
  int pwt_idx=0;		     // indicating P(w|t) storage index
  int pwt_rows=(int)maxit/100;
  int pwt_columns=(int)(K+K*V);
  int VK=(int)V*K;
  IntegerMatrix pwt_raw(pwt_rows,pwt_columns);
  
  // main gibbs loop
  for(int j = 0; j < maxit; ++j){
    // sub-steps
    for(int i = 0; i < sumL; ++i){
      doc=dwt(i,0);
      word=dwt(i,1);
      topic=dwt(i,2);

      // std::cout << doc << " " << word << " " << topic << std::endl;
      // std::cout << i << std::endl;
      
      // cvk[topic,word] <- cvk[topic,word]-1
      idx=(topic-1)*V+word-1;
      cvk(idx,2)=cvk(idx,2)-1;
      // cik[doc,topic] <- cik[doc,topic]-1
      idx=(doc-1)*K+topic-1;
      cik(idx,2)=cik(idx,2)-1;
      // std::cout <<  cik(idx,2) << std::endl;
      // ck[topic] <- ck[topic]-1 #ck must be ordered!
      ck[topic-1]=ck[topic-1]-1;
      // probs <- (cvk[,word]+g)/(ck+V*g) * (cik[doc,]+a)/(L[doc]+K*a)
      for(int k = 0; k < K; ++k){
	probs[k]= ((double)cvk(k*V+word-1,2)+g)/
	  ((double)ck[k]+(double)V*g)*
	  ((double)cik((doc-1)*K+k,2)+a)/
	  ((double)L[doc-1]+(double)K*a);
      }
      // topic <- Cat(probs)
      topic=one_cat(probs);
      // std::cout << probs[0] << std::endl;
      // dwt$topic[i] <- topic               #replace the old topic
      dwt(i,2)=topic;
      // cvk[topic,word] <- cvk[topic,word]+1
      idx=(topic-1)*V+word-1;
      cvk(idx,2)=cvk(idx,2)+1;
      // cik[doc,topic] <- cik[doc,topic]+1
      idx=(doc-1)*K+topic-1;
      cik(idx,2)=cik(idx,2)+1;
      // ck[topic] <- ck[topic]+1
      ck[topic-1]=ck[topic-1]+1;
      // if(i>=4665500){
      // 	std::cout << i << std::endl;
      // }
    }

    // store raw data in pwt_raw to calculate log(P(w|t))
    if(j%(int)100 == 99){
      pwt_idx=(int)j/(int)100;
      for(int l = 0; l < VK; ++l){
    	pwt_raw(pwt_idx,l)=cvk(l,2);
      }
      for(int l = VK; l < pwt_columns; ++l){
    	pwt_raw(pwt_idx,l)=ck[l-VK];
      }
    }
    // std::cout << j << "th iteration done" << std::endl;
    printProgress((double)j/maxit);
  }

  // return Rcpp::List::create(Rcpp::Named("cvk") = cvk,
  // 			    Rcpp::Named("dwt") = dwt,
  // 			    Rcpp::Named("cik") = cik,
  // 			    Rcpp::Named("ck") = ck);
    return Rcpp::List::create(Rcpp::Named("cvk") = cvk,
  			    Rcpp::Named("dwt") = dwt,
  			    Rcpp::Named("cik") = cik,
  			    Rcpp::Named("ck") = ck,
  			    Rcpp::Named("pwt_raw") = pwt_raw);
}

