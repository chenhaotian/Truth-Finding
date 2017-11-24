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

// mn: each entity has MULTIPLE true attributes, attributes are NOT sharing among entities.
// [[Rcpp::export]]
List truthfinding_mn(IntegerVector facts,IntegerVector fcidx, IntegerMatrix claims, IntegerMatrix ctsc, NumericVector beta, NumericMatrix alpha0, NumericMatrix alpha1,int burnin, int maxit, int sample_step){

  // claims: fid sid o
  // ctsc:   t   sid o count

  int nfacts=facts.size();
  int nsources=alpha0.nrow();
  int expand_source_claim=nsources*2;
  int idx=0;
  int startidx=0,endidx=0;	// claim start and end index for each fact
  int sid=0,o=0,t=0;		// tmp variable in locating claims in ctsc
  int fact_pre=0;
  NumericVector probs(2);
  double conditional_claim0=1, conditional_claim1=1;

  // outputs
  int sample_size = maxit/sample_step - burnin/sample_step;
  if(sample_size <=0){
    std::cout << "sample_size = maxit/sample_step - burnin/sample_step <=0!" << std::endl;
    std::exit(-1);
  }
  NumericVector facts_out(nfacts,(double)0);
  NumericVector ctsc_out(ctsc.nrow(),(double)0);
  NumericVector recall(nsources,(double)0);
  NumericVector specificity(nsources,(double)0);
  NumericVector precision(nsources,(double)0);
  
  // main gibbs loop
  int it = 0;
  while(it < maxit){
    // 1. for each fact f...
    for(int f = 0; f < nfacts; ++f){
      startidx=fcidx[f];		// fact->claim index
      endidx=fcidx[f+1];
      fact_pre=facts[f];
      
      // 2. for each claim c in fact f...
      for(int c = startidx; c < endidx; ++c){
	sid=claims(c,1);
	o=claims(c,2);
	
	// 2.0 generate sample probability for t=0
	// alpha0: prior for specificity
	t=0;			// when f is false
	conditional_claim0=1;	//reset conditional claim probability
	idx=sid*2;
	conditional_claim0 = conditional_claim0*
	  ((double)ctsc(idx+o,3)-(double)1+alpha0(sid,o))/ //equation (2)
	  ((double)ctsc(idx,3)+(double)ctsc(idx+1,3)-(double)1+alpha0(sid,0)+alpha0(sid,1));
	
	// 2.1 generate sample probability for t=1
	// alpha1: prior for recall
	t=1;			// when f is true
	conditional_claim1=1;	// reset conditional claim probability
	idx= expand_source_claim + idx;
	conditional_claim1=conditional_claim1*
	  ((double)ctsc(idx+o,3)-(double)1+alpha1(sid,1-o))/ //equation (2)
	  ((double)ctsc(idx,3)+(double)ctsc(idx+1,3)-(double)1+alpha1(sid,0)+alpha1(sid,1));
      }
      
      probs[0] = beta[1]*conditional_claim0;
      probs[1] = beta[0]*conditional_claim1;
      
      // 3. sample and update facts
      facts[f]=one_cat_zero_begin(probs);
      
      // 4. update ctsc
      if(facts[f] != fact_pre){
	for(int c = startidx; c < endidx; ++c){
	  sid=claims(c,1);
	  o=claims(c,2);
	  
	  // pre - 1
	  idx=fact_pre*expand_source_claim+sid*2+o;
	  ctsc(idx,3)=ctsc(idx,3)-1;
	  
	  // aft + 1
	  idx=facts[f]*expand_source_claim+sid*2+o;
	  ctsc(idx,3)=ctsc(idx,3)+1;
	}
      }
    }

    // sample output
    it=it+1;
    if((it > burnin) & (it % sample_step == 0)){
      for(int l = 0; l < nfacts; ++l){
	facts_out[l] = facts_out[l] + facts[l]/(double)sample_size;
      }
      for(int l = 0; l < ctsc.nrow(); ++l){
	ctsc_out[l] = ctsc_out[l] + (double)ctsc(l,3)/(double)sample_size;
      }
    }

    printProgress((double)it/maxit);
    
  }


  for(int s = 0; s < nsources; ++s){
    // recall
    idx=expand_source_claim+s*(int)2+1;
    recall[s] = ((double)ctsc_out[idx]+alpha1(s,0))/
      ((double)ctsc_out[idx]+(double)ctsc_out[idx-1]+alpha1(s,0)+alpha1(s,1));
    
    // precision
    precision[s] = ((double)ctsc_out[idx]+alpha1(s,1))/
      ((double)ctsc_out[idx] + (double)ctsc_out[s*(int)2+(int)1] + alpha0(s,0) + alpha1(s,1));
    
    // specificity
    idx=s*(int)2;
    specificity[s] = ((double)ctsc_out[idx]+alpha0(s,0))/
      ((double)ctsc_out[idx]+(double)ctsc_out[idx+1]+alpha0(s,0)+ alpha0(s,1));
  }
  
  return Rcpp::List::create(Rcpp::Named("facts_out") = facts_out,
			    Rcpp::Named("ctsc_out") = ctsc_out,
			    Rcpp::Named("sample_size") = sample_size,
			    Rcpp::Named("recall") = recall,
			    Rcpp::Named("precision") = precision,
			    Rcpp::Named("specificity") = specificity);
}

// sn: each entity has SINGLE true attribute, attributes are NOT sharing among entities.
// [[Rcpp::export]]
List truthfinding_sn(IntegerVector ecidx, IntegerVector e_truths, IntegerVector e_n_attr, IntegerVector s_n_right, IntegerVector s_n_claims, IntegerMatrix rawdb, double beta, NumericMatrix alpha1,int max_nattributes, int burnin, int maxit, int sample_step){

  // ctsc:   t   sid o count
  int nentities=e_truths.size();
  int nsources=alpha1.nrow();
  int nattributes=0;
  int startidx=0,endidx=0;	// claim start and end index for each fact
  int sid=0;   		// tmp variable in locating claims
  int truth_pre=0;
  NumericVector probs(max_nattributes+1); // initialize probability vector

  // outputs
  int sample_size = maxit/sample_step - burnin/sample_step;
  if(sample_size <=0){
    std::cout << "sample_size = maxit/sample_step - burnin/sample_step <=0!" << std::endl;
    std::exit(-1);
  }
  NumericVector s_n_right_out(nsources,(double)0);
  NumericVector precision(nsources,(double)0);

  // rawdb:
  // entity attribute source
  
  // main gibbs loop
  int it = 0;
  while(it < maxit){
    // 1. for each entity ...
    for(int eid = 0; eid < nentities; ++eid){
      // initialize indexes
      startidx=ecidx[eid];
      endidx=ecidx[eid+1];
      truth_pre=e_truths[eid];
      nattributes=e_n_attr[eid];

      // 2. for each attribute, calculate it's probability of being right
      // note that a<=nattributes, the nattributes-th iteration is to calculate the probability that none of the claims are right
      for(int a = 0; a < nattributes; ++a){
      // for(int a = 0; a <= nattributes; ++a){
	// initialize & prior
	probs[a]=beta;
	for(int i = startidx; i < endidx; ++i){
	  sid=rawdb(i,2);
	  if(rawdb(i,1)==a){	       // p(a|claim)
	    if(rawdb(i,1)==truth_pre){
	      probs[a] = probs[a]*
		((double)s_n_right[sid]- (double)1 + alpha1(sid,0))/
		((double)s_n_claims[sid]- (double)1 + alpha1(sid,0) + alpha1(sid,0));
	    }else{
	      probs[a] = probs[a]*
		((double)s_n_right[sid] + alpha1(sid,0))/
		((double)s_n_claims[sid]- (double)1 + alpha1(sid,0) + alpha1(sid,0));
	    }
	  }else{		       // p(hat(a)|claim)
	    if(rawdb(i,1)==truth_pre){
	      probs[a] = probs[a]*
		((double)1 - 
		 ((double)s_n_right[sid]- (double)1 + alpha1(sid,0))/
		 ((double)s_n_claims[sid]- (double)1 + alpha1(sid,0) + alpha1(sid,0)));
	    }else{
	      probs[a] = probs[a]*
		((double)1 - 
		 ((double)s_n_right[sid] + alpha1(sid,0))/
		 ((double)s_n_claims[sid]- (double)1 + alpha1(sid,0) + alpha1(sid,0)));
	    }
	  }
	}
      }

      e_truths[eid]=one_cat_zero_begin(probs[seq(0,nattributes-1)]);
      // e_truths[eid]=one_cat_zero_begin(probs[seq(0,nattributes)]);
      
      // update
      if(truth_pre != e_truths[eid]){
	for(int i = startidx; i < endidx; ++i){
	  sid=rawdb(i,2);
	  if(rawdb(i,1)==truth_pre){ // minus
	    s_n_right[sid] = s_n_right[sid] - 1;
	  }else if(rawdb(i,1)==e_truths[eid]){ // plus
	    s_n_right[sid] = s_n_right[sid] + 1;
	  }
	}
      }
      
    }
    
    // sample output
    it=it+1;
    if((it > burnin) & (it % sample_step == 0)){
      // for(int l = 0; l < nentities; ++l){
      // 	e_truths_out[l]
      // 	entities_out[l] = entities_out[l] + entities[l]/(double)sample_size;
      // }
      for(int s = 0; s < nsources; ++s){
	s_n_right_out[s] = s_n_right_out[s] + (double)s_n_right[s] / (double)sample_size;
      }
    }

    printProgress((double)it/maxit);
    
  }
  
  // esitmate source quality
  for(int s = 0; s < nsources; ++s){
    // precision
    precision[s] = ((double)s_n_right_out[s]+alpha1(s,0))/
      ((double)s_n_claims[s] + alpha1(s,0) + alpha1(s,1));
  }
  
  return Rcpp::List::create(Rcpp::Named("e_truths") = e_truths,
			    Rcpp::Named("s_n_right_out") = s_n_right_out,
			    Rcpp::Named("s_n_right") = s_n_right,
			    Rcpp::Named("sample_size") = sample_size,
			    Rcpp::Named("precision") = precision);
}

// ss: each entity has SINGLE true attribute, attributes are SHARING among entities.
// [[Rcpp::export]]
List truthfinding_ss(IntegerVector ecidx, IntegerVector e_truths, IntegerVector s_n_right, IntegerVector s_n_claims, IntegerMatrix rawdb, NumericVector beta, NumericVector pi, NumericMatrix alpha1,int nattributes, int nsources, int nentities, int burnin, int maxit, int sample_step){

  int startidx=0,endidx=0;	// claim start and end index for each fact
  int sid=0;
  int truth_pre=0;
  NumericVector probs(nattributes); // initialize probability vector
  
  if(beta.size()!=nattributes){
    std::cout << "length(beta) not equal to number of unique attributes!" << std::endl;
    std::exit(-1);
  }
  if(alpha1.nrow()!=nsources){
    std::cout << "nrow(alpha1) not equal to number_of_unique_sources!" << std::endl;
    std::exit(-1);
  }

  // outputs
  int sample_size = maxit/sample_step - burnin/sample_step;
  if(sample_size <=0){
    std::cout << "sample_size = maxit/sample_step - burnin/sample_step <=0!" << std::endl;
    std::exit(-1);
  }
  NumericVector s_n_right_out(nsources,(double)0);
  NumericVector precision(nsources,(double)0);

  // rawdb:
  // entity attribute source
  
  // main gibbs loop
  int it = 0;
  while(it < maxit){
    // 1. for each entity ...
    for(int eid = 0; eid < nentities; ++eid){
      // initialize indexes
      startidx=ecidx[eid];
      endidx=ecidx[eid+1];
      truth_pre=e_truths[eid];

      // 2. for each attribute, calculate it's probability of being right
      for(int a = 0; a < nattributes; ++a){
	// initialize & prior
	// probs[a] = beta[a];
	if(truth_pre==a){
	  probs[a] = (double)(pi[a]-1)+beta[a];
	  // probs[a] = 1;
	}else{
	  probs[a] = (double)pi[a]+beta[a];
	  // probs[a] = 1;
	}
	for(int i = startidx; i < endidx; ++i){
	  sid=rawdb(i,2);
	  if(rawdb(i,1)==a){	       // p(a|claim)
	    if(rawdb(i,1)==truth_pre){
	      probs[a] = probs[a]*
		((double)s_n_right[sid]- (double)1 + alpha1(sid,0))/
		((double)s_n_claims[sid]- (double)1 + alpha1(sid,0) + alpha1(sid,0));
	    }else{
	      probs[a] = probs[a]*
		((double)s_n_right[sid] + alpha1(sid,0))/
		((double)s_n_claims[sid]- (double)1 + alpha1(sid,0) + alpha1(sid,0));
	    }
	  }else{		       // p(hat(a)|claim)
	    if(rawdb(i,1)==truth_pre){
	      probs[a] = probs[a]*
		((double)1 -
		 ((double)s_n_right[sid]- (double)1 + alpha1(sid,0))/
		 ((double)s_n_claims[sid]- (double)1 + alpha1(sid,0) + alpha1(sid,0)));
	    }else{
	      probs[a] = probs[a]*
		((double)1 -
		 ((double)s_n_right[sid] + alpha1(sid,0))/
		 ((double)s_n_claims[sid]- (double)1 + alpha1(sid,0) + alpha1(sid,0)));
	    }
	  }
	}
      }

      e_truths[eid]=one_cat_zero_begin(probs);
      
      // update
      if(truth_pre != e_truths[eid]){
	// update pi
	pi[truth_pre] = pi[truth_pre]-1;
	pi[e_truths[eid]] = pi[e_truths[eid]]+1;
	// update s_n_right
	for(int i = startidx; i < endidx; ++i){
	  sid=rawdb(i,2);
	  if(rawdb(i,1)==truth_pre){ // minus
	    s_n_right[sid] = s_n_right[sid] - 1;
	  }else if(rawdb(i,1)==e_truths[eid]){ // plus
	    s_n_right[sid] = s_n_right[sid] + 1;
	  }
	}
      }
      
    }
    
    // sample output
    it=it+1;
    if((it > burnin) & (it % sample_step == 0)){
      // for(int l = 0; l < nentities; ++l){
      // 	e_truths_out[l]
      // 	entities_out[l] = entities_out[l] + entities[l]/(double)sample_size;
      // }
      for(int s = 0; s < nsources; ++s){
	s_n_right_out[s] = s_n_right_out[s] + (double)s_n_right[s] / (double)sample_size;
      }
    }

    printProgress((double)it/maxit);
    
  }
  
  // esitmate source quality
  for(int s = 0; s < nsources; ++s){
    // precision
    precision[s] = ((double)s_n_right_out[s]+alpha1(s,0))/
      ((double)s_n_claims[s] + alpha1(s,0) + alpha1(s,1));
  }
  
  return Rcpp::List::create(Rcpp::Named("e_truths") = e_truths,
			    Rcpp::Named("s_n_right_out") = s_n_right_out,
			    Rcpp::Named("sample_size") = sample_size,
			    Rcpp::Named("precision") = precision,
			    Rcpp::Named("s_n_right") = s_n_right,
			    Rcpp::Named("pi") = pi);
}

// ss with full parameter
// TO DO:
// 1. change alpha1 from double to numericmatrix, to support online learning
// 2. change beta from double to numericvector of length nattributes
// [[Rcpp::export]]
List truthfinding_ss_fullpar(IntegerVector ecidx, IntegerVector e_truths, IntegerVector s_aa_n_claims, IntegerVector s_a_n_claims, IntegerMatrix rawdb, double beta, NumericVector pi, double alpha1,int nattributes, int nsources, int nentities, int burnin, int maxit, int sample_step){

  int nsourceattributes = nsources*nattributes*nattributes;
  int startidx=0,endidx=0;	// claim start and end index for each fact
  int idx1=0,idx2=0;   		// tmp variable in locating claims
  int truth_pre=0;
  NumericVector probs(nattributes); // initialize probability vector
  
  // if(beta.size()!=nattributes){
  //   std::cout << "length(beta) not equal to number of unique attributes!" << std::endl;
  //   std::exit(-1);
  // }
  
  // outputs
  int sample_size = maxit/sample_step - burnin/sample_step;
  if(sample_size <=0){
    std::cout << "sample_size = maxit/sample_step - burnin/sample_step <=0!" << std::endl;
    std::exit(-1);
  }
  NumericVector s_aa_n_claims_out(nsourceattributes,(double)0);
  NumericVector log_likelihood(sample_size,(double)0);
  NumericMatrix e_truths_out(nentities,nattributes); // matrix with nentities rows and nattributes columns, filled with 0s.

  // rawdb:
  // entity attribute source

  
  // main gibbs loop
  int it = 0;
  while(it < maxit){
    // 1. for each entity ...
    for(int eid = 0; eid < nentities; ++eid){
      // initialize indexes
      startidx=ecidx[eid];
      endidx=ecidx[eid+1];
      truth_pre=e_truths[eid];

      // 2. for each attribute, calculate it's probability
      for(int a = 0; a < nattributes; ++a){
	// initialize & prior
	if(truth_pre==a){
	  // probs[a] = (double)(pi[a]-1)+beta[a];
	  probs[a] = (double)(pi[a]-1)+beta;
	}else{
	  // probs[a] = (double)pi[a]+beta[a];
	  probs[a] = (double)pi[a]+beta;
	}
	for(int i = startidx; i < endidx; ++i){
	  idx1=rawdb(i,2)*nattributes*nattributes+a*nattributes+rawdb(i,1);
	  idx2=rawdb(i,2)*nattributes+a;
	  if(truth_pre==a){
	    probs[a] = probs[a]*
	      ((double)s_aa_n_claims[idx1] - (double)1 + alpha1)/
	      ((double)s_a_n_claims[idx2] - (double)1 + (double)nattributes*alpha1);
	  }else{
	    probs[a] = probs[a]*
	      ((double)s_aa_n_claims[idx1] + alpha1)/
	      ((double)s_a_n_claims[idx2] - (double)1 + (double)nattributes*alpha1);
	  }
	}
      }

      e_truths[eid]=one_cat_zero_begin(probs);
      
      // update
      if(truth_pre != e_truths[eid]){
	// update pi
	pi[truth_pre] = pi[truth_pre]-1;
	pi[e_truths[eid]] = pi[e_truths[eid]]+1;
	// update s_aa_n_claims
	for(int i = startidx; i < endidx; ++i){
	  idx1=rawdb(i,2)*nattributes*nattributes+rawdb(i,1);
	  idx2=rawdb(i,2)*nattributes;
	  s_aa_n_claims[idx1+truth_pre*nattributes] = s_aa_n_claims[idx1+truth_pre*nattributes]-1;
	  s_a_n_claims[idx2+truth_pre] = s_a_n_claims[idx2+truth_pre]-1;
	  s_aa_n_claims[idx1+e_truths[eid]*nattributes] = s_aa_n_claims[idx1+e_truths[eid]*nattributes]+1;
	  s_a_n_claims[idx2+e_truths[eid]] = s_a_n_claims[idx2+e_truths[eid]]+1;
	}
      }
      
    }
    
    // sample output
    it=it+1;
    if((it > burnin) & (it % sample_step == 0)){
      for(int l = 0; l < nentities; ++l){
	e_truths_out(l,e_truths[l]) = e_truths_out(l,e_truths[l])+1;
      }
      for(int l = 0; l < nsourceattributes; ++l){
	s_aa_n_claims_out[l] = s_aa_n_claims_out[l] + (double)s_aa_n_claims[l]/(double)sample_size;
      }
    }
    
    printProgress((double)it/maxit);
    
  }
  
  return Rcpp::List::create(Rcpp::Named("e_truths") = e_truths,
			    Rcpp::Named("e_truths_out") = e_truths_out,
			    Rcpp::Named("sample_size") = sample_size,
			    Rcpp::Named("s_aa_n_claims") = s_aa_n_claims,
			    Rcpp::Named("s_aa_n_claims_out") = s_aa_n_claims_out,
			    Rcpp::Named("pi") = pi);
}

// p(t1|a1,a2) = p(a1|t1,a2)*p(t1|a2) / p(a1|a2)
//             = p(a1|t1)*p(t1|a2) / p(a1)
//             = p(a1|t1)*p(a2|t1)*p(t1)/ p(a2) / p(a1)
//             = p(a1|t1)*p(a2|t1)*p(t1)/ p(a2) / p(a1)
//             ~ p(a1|t1)*p(a2|t1)*p(t1)
// ...
// p(t1|a1,a2,...,an) ~ p(a1|t1)*p(a2|t1)*...*p(an|t1)*p(t1)

