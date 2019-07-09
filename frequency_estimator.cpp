#include "frequency_estimator.h"

frequency_estimator::frequency_estimator(Eigen::MatrixXd* _pEigVecs, double _tol, double _maxLambda) :
  skipIf(false), skipInfo(false), siteOnly(false), nelderMead(false), assumeHWD(false), gtError(0.005) 
{
  if ( _pEigVecs == NULL )
    error("[E:%s:%d %s] Invalid eigenvectors", __FILE__, __LINE__, __PRETTY_FUNCTION__ );

  pEigVecs = _pEigVecs;
  nsamples = (int32_t)pEigVecs->rows();
  ndims = (int32_t)pEigVecs->cols();

  pSVD = new Eigen::BDCSVD<Eigen::MatrixXd>(*_pEigVecs, Eigen::ComputeThinU | Eigen::ComputeThinV);

  // set dimensions;
  tol = _tol;
  maxLambda = _maxLambda;  
  
  hdr = NULL; wdr = NULL; iv = NULL;
  hwe0z = hwe1z = ibc0 = ibc1 = 0;

  pls = NULL; n_pls = 0; ploidies = NULL;
  ifs = new float[nsamples];
  betas = new float[ndims];
  theta = 0;
  pooled_af = -1;
  isaf_computed = false;

  gt = NULL;
  gq = NULL;
}

frequency_estimator::~frequency_estimator() {
  if ( ( pEigVecs ) && ( pSVD ) )
    delete pSVD;
  if (ifs != NULL )
    delete [] ifs;
  if ( gt != NULL )
    free(gt);
  if ( gq != NULL )
    free(gq);
}

frequency_estimator::frequency_estimator(Eigen::BDCSVD<Eigen::MatrixXd>* _pSVD, double _tol, double _maxLambda) :
  skipIf(false), skipInfo(false), siteOnly(false), gtError(0.005)   
{
  // set dimensions;
  pEigVecs = NULL;

  pSVD = _pSVD;
  nsamples = (int32_t)pSVD->matrixU().rows();
  ndims = (int32_t)pSVD->matrixU().cols();

  tol = _tol;
  maxLambda = _maxLambda;    
  
  hdr = NULL; wdr = NULL; iv = NULL;
  hwe0z = hwe1z = ibc0 = ibc1 = 0;

  pls = NULL; n_pls = 0; ploidies = NULL;
  ifs = new float[nsamples];
  betas = new float[ndims];
  theta = 0;
  pooled_af = -1;
  isaf_computed = false;  
}

bool frequency_estimator::set_hdr(bcf_hdr_t* _hdr, bcf_hdr_t* _wdr ) {
  if ( hdr != _hdr ) {
    hdr = _hdr;
    wdr = _wdr;
    char buffer[65535];
    if ( !skipInfo ) {
      if ( bcf_hdr_id2int(_hdr, BCF_DT_ID, "HWE_SLP_P" ) < 0 ) {
	sprintf(buffer,"##INFO=<ID=HWE_SLP_P,Number=1,Type=Float,Description=\"Signed log p-value of HWE test with pooled allele frequencies\">\n");
	bcf_hdr_append(wdr, buffer);
      }
      if ( bcf_hdr_id2int(_hdr, BCF_DT_ID, "FIBC_P" ) < 0 ) {      
	sprintf(buffer,"##INFO=<ID=FIBC_P,Number=1,Type=Float,Description=\"Inbreeding coefficient with pooled allele frequencies\">\n");
	bcf_hdr_append(wdr, buffer);
      }
      if ( bcf_hdr_id2int(_hdr, BCF_DT_ID, "HWE_SLP_I" ) < 0 ) {      
	sprintf(buffer,"##INFO=<ID=HWE_SLP_I,Number=1,Type=Float,Description=\"Signed log p-value of HWE test with individual-sepcific allele frequencies\">\n");
	bcf_hdr_append(wdr, buffer);
      }
      if ( bcf_hdr_id2int(_hdr, BCF_DT_ID, "FIBC_I" ) < 0 ) {
	sprintf(buffer,"##INFO=<ID=FIBC_I,Number=1,Type=Float,Description=\"Inbreeding coefficient with individual-sepcific allele frequencies\">\n");
	bcf_hdr_append(wdr, buffer);
      }
      if ( bcf_hdr_id2int(_hdr, BCF_DT_ID, "MAX_IF" ) < 0 ) {
	sprintf(buffer,"##INFO=<ID=MAX_IF,Number=1,Type=Float,Description=\"Maximum Individual-specific allele frequency\">\n");
	bcf_hdr_append(wdr, buffer);
      }
      if ( bcf_hdr_id2int(_hdr, BCF_DT_ID, "MIN_IF" ) < 0 ) {      
	sprintf(buffer,"##INFO=<ID=MIN_IF,Number=1,Type=Float,Description=\"Minimum Individual-specific allele frequency\">\n");
	bcf_hdr_append(wdr, buffer);
      }
      if ( bcf_hdr_id2int(_hdr, BCF_DT_ID, "BETA_IF" ) < 0 ) {      
	sprintf(buffer,"##INFO=<ID=BETA_IF,Number=%d,Type=Float,Description=\"Coefficients for intercept and each eigenvector to obtain ISAF\">\n", ndims);
	bcf_hdr_append(wdr, buffer);
      }      
      if ( bcf_hdr_id2int(_hdr, BCF_DT_ID, "LLK0" ) < 0 ) {      
	sprintf(buffer,"##INFO=<ID=LLK0,Number=1,Type=Float,Description=\"Null likelihood - for debug purpose\">\n");
	bcf_hdr_append(wdr, buffer);
      }
    }
    if ( ( !skipIf ) && ( !siteOnly ) ) {
      sprintf(buffer,"##FORMAT=<ID=IF,Number=1,Type=Float,Description=\"Individual-specific allele frequencies\">\n");
      bcf_hdr_append(wdr, buffer);
    }
    bcf_hdr_sync(wdr);
    return true;
  }
  else return false;
}

bool frequency_estimator::set_variant(bcf1_t* _iv, int8_t* _ploidies, int32_t* _pl) { //, std::vector<int32_t>* p_icols) {
  iv = _iv;
  ploidies = _ploidies;
  if ( iv->n_sample != nsamples )
    error("[E:%s:%d %s] nsamples %d != %d in the EigenVector", __FILE__, __LINE__, __PRETTY_FUNCTION__, iv->n_sample, nsamples);

  //error("%d %d %d",bcf_hdr_nsamples(hdr),iv->n_sample,nsamples);      
    
  // parse PL fields
  bcf_unpack(iv, BCF_UN_ALL);
    
  if ( _pl != NULL ) { pls = _pl; n_pls = 3*nsamples; }
  else {
    bool plfound = false;
    
    if ( field.empty() || ( field.compare("PL") == 0 ) ) {
      if ( bcf_get_format_int32(hdr, iv, "PL", &pls, &n_pls) < 0 ) {
	if ( field.compare("PL") == 0 ) 
	  error("[E:%s:%d %s] Cannot parse PL field", __FILE__, __LINE__, __PRETTY_FUNCTION__);
      }
      else
	plfound = true;
    }
    
    if ( ( (!plfound) && field.empty() ) || field.compare("GL") == 0 ){
      float* gls = NULL;
      int32_t n_gls = 0;
      if ( bcf_get_format_float(hdr, iv, "GL", &gls, &n_gls) < 0 ) {
	error("[E:%s:%d %s] Cannot parse GL field", __FILE__, __LINE__, __PRETTY_FUNCTION__);
      }
      else {
	if ( pls == NULL ) pls = new int32_t[n_gls];      
	for(int32_t i=0; i < nsamples; ++i) {
	  float maxgl = gls[3*i];
	  if ( gls[3*i+1] > maxgl ) maxgl = gls[3*i+1];
	  if ( gls[3*i+2] > maxgl ) maxgl = gls[3*i+2];
	  pls[3*i] = (int32_t)floor(-10*(gls[3*i]-maxgl)+0.5);
	  pls[3*i+1] = (int32_t)floor(-10*(gls[3*i+1]-maxgl)+0.5);
	  pls[3*i+2] = (int32_t)floor(-10*(gls[3*i+2]-maxgl)+0.5);
	  if ( pls[3*i] > 255 ) pls[3*i] = 255;
	  if ( pls[3*i+1] > 255 ) pls[3*i+1] = 255;
	  if ( pls[3*i+2] > 255 ) pls[3*i+2] = 255;
	}
	free(gls);
	n_pls = n_gls;
      }
    }
    else if ( field.compare("GT") == 0 ) {
      int32_t* gts = NULL;
      int32_t n_gts = 0;
      
      double tmp = gtError + gtError*gtError;
      int32_t errM[9] =
	{ 0, (int32_t)floor(-10*log10(gtError/(1-tmp))+0.5), (int32_t)floor(-10*log10(gtError*gtError/(1-tmp))+0.5),
	  (int32_t)floor(-10*log10(tmp/(2-tmp-tmp))+0.5), 0, 0,
	  0, 0, 0 };
      errM[5] = errM[3];
      errM[6] = errM[2];    
      errM[7] = errM[1];
      
      if ( bcf_get_genotypes(hdr, iv, &gts, &n_gts) < 0 ) {
	error("[E:%s:%d %s] Cannot parse GT field", __FILE__, __LINE__, __PRETTY_FUNCTION__);
      }
      else {
	int32_t max_ploidy = n_gts/nsamples;
	if ( max_ploidy != 2 )
	  error("[E:%s:%d %s] Multi-allelic (or Mono-allelic) variants found", __FILE__, __LINE__, __PRETTY_FUNCTION__);	
	if ( pls == NULL ) pls = new int32_t[nsamples*3];      
	for(int32_t i=0; i < nsamples; ++i) {
	  int32_t geno = bcf_gt_allele(gts[2*i])+bcf_gt_allele(gts[2*i+1]);	
	  if ( bcf_gt_is_missing(gts[2*i+1]) ) { // haploid or missing
	    if ( bcf_gt_is_missing(gts[2*i]) ) { // missing
	      pls[3*i] = pls[3*i+1] = pls[3*i+2] = 0;
	      continue;
	    }
	    else { // pretend to be homozygous for haploid
	      geno = bcf_gt_allele(gts[2*i]) + bcf_gt_allele(gts[2*i]);
	    }
	  }
	  pls[3*i] = errM[geno*3];
	  pls[3*i+1] = errM[geno*3+1];
	  pls[3*i+2] = errM[geno*3+2];
	}
	free(gts);
	n_pls = nsamples*3;
      }
    }
    else {
      if ( !plfound ) 
	error("[E:%s:%d %s] Cannot recognize the field [%s]", __FILE__, __LINE__, __PRETTY_FUNCTION__, field.c_str());
    }
  }
  /*
  else if ( bcf_get_format_int32(hdr, iv, "PL", &pls, &n_pls) < 0 ) {
    //float gls[nsamples+1];
    float* gls = NULL;
    int32_t n_gls = 0;
    if ( bcf_get_format_float(hdr, iv, "GL", &gls, &n_gls) < 0 ) {
      error("[E:%s:%d %s] Cannot parse PL or GL field", __FILE__, __LINE__, __PRETTY_FUNCTION__);
    }
    else {
      if ( pls == NULL ) pls = new int32_t[n_gls];      
      for(int32_t i=0; i < nsamples; ++i) {
	float maxgl = gls[3*i];
	if ( gls[3*i+1] > maxgl ) maxgl = gls[3*i+1];
	if ( gls[3*i+2] > maxgl ) maxgl = gls[3*i+2];
	pls[3*i] = (int32_t)floor(-10*(gls[3*i]-maxgl)+0.5);
	pls[3*i+1] = (int32_t)floor(-10*(gls[3*i+1]-maxgl)+0.5);
	pls[3*i+2] = (int32_t)floor(-10*(gls[3*i+2]-maxgl)+0.5);
	if ( pls[3*i] > 255 ) pls[3*i] = 255;
	if ( pls[3*i+1] > 255 ) pls[3*i+1] = 255;
	if ( pls[3*i+2] > 255 ) pls[3*i+2] = 255;
      }
      free(gls);      
    }
  }
  */
  pooled_af = -1;
  isaf_computed = false;
  
  return true;
  //else return false;
}

double frequency_estimator::estimate_pooled_af_em(int32_t maxiter) {
  if ( pooled_af < 0 ) {
    double p = 0.5, q = 0.5;
    double p0, p1, p2, sump, apsum, an;
    for(int32_t i=0; i < maxiter; ++i) {
      apsum = 0;
      an = 0;
      for(int32_t j=0; j < nsamples; ++j) {
	if ( ploidies[j] == 2 ) {
	  p0 = q*q*phredConv.toProb(pls[j*3]);
	  p1 = 2*p*q*phredConv.toProb(pls[j*3+1]);
	  p2 = p*p*phredConv.toProb(pls[j*3+2]);
	  sump = p0+p1+p2;
	  apsum += (p1/sump + p2/sump*2.0);
	  an += 2;
	}
	else if ( ploidies[j] == 1 ) {
	  p0 = q*phredConv.toProb(pls[j*3]);
	  p2 = p*phredConv.toProb(pls[j*3+2]);
	  sump = p0+p2;
	  apsum += (p2/sump);
	  ++an;
	}
      }
      p = apsum/an;
      q = 1.0-p;
    }
    pooled_af = p;
    if ( pooled_af * nsamples < 0.5 ) { pooled_af = 0.5 / ( 1 + 2 *nsamples ); }
    else if ( (1-pooled_af)*nsamples < 0.5 ) { pooled_af = 1 - 0.5/(1 + 2*nsamples); }
  }
  return pooled_af;
}

bool frequency_estimator::score_test_hwe(bool use_isaf) {
  estimate_pooled_af_em();    
    
  double pp0 = (1.-pooled_af)*(1.-pooled_af);
  double pp1 = 2*pooled_af*(1-pooled_af);
  double pp2 = pooled_af*pooled_af;
  double sumU0 = 0, sqU0 = 0, sumU1 = 0, sqU1 = 0;
  double obsH0 = 0, obsH1 = 0, expH1 = 0;
  double l0, l1, l2, sum1, sum0, ph1, ph0, U0, U1;
  int32_t ndiploids = 0;

  // pretend that we have pp0, pp1, pp2 observations of each genotype (pseudocount)
  sumU0 = sumU1 = pp0*(1/pp0) + pp1 * (-2/pp1) + pp2*(1/pp2);
  sqU0 = sqU1 = pp0*1/pp0/pp0 + pp1 * 4/pp1/pp1 + pp2*(1/pp2/pp2);

  assumeHWD = false;

  for(int32_t i=0; i < nsamples; ++i) {
    if ( ploidies[i] != 2 ) continue;
    ++ndiploids;
    l0 = phredConv.toProb(pls[i*3]);
    l1 = phredConv.toProb(pls[i*3+1]);
    l2 = phredConv.toProb(pls[i*3+2]);

    sum0 = l0*pp0 + l1*pp1 + l2*pp2 + 1e-100;
    ph0 = l1*pp1;
    obsH0 += (ph0/sum0);
    U0 = pooled_af*(1-pooled_af)*(l0-2*l1+l2)/sum0;
    sumU0 += U0;
    sqU0 += (U0*U0);
    
    if ( use_isaf ) {
      if ( nelderMead )
	estimate_isaf_simplex();
      else
	estimate_isaf_em();

      sum1 = l0*(1.-ifs[i])*(1.-ifs[i]) + 2*l1*(1.-ifs[i])*ifs[i] + l2*ifs[i]*ifs[i] + 1e-100;
      ph1 = 2*l1*(1.-ifs[i])*ifs[i];
      obsH1 += (ph1/sum1);
      expH1 += (2*(1.-ifs[i])*ifs[i]);
      U1 = (1-ifs[i])*ifs[i]*(l0-2*l1+l2)/sum1;
      sumU1 += U1;
      sqU1 += (U1*U1);
    }
  }

  hwe0z = sumU0/sqrt(sqU0);
  ibc0 = 1.0 - (obsH0+1)/(pp1*ndiploids+1);

  if ( use_isaf ) {
    hwe1z = sumU1/sqrt(sqU1);
    ibc1 = 1.0 - (obsH1+1)/(expH1+1);
  }

  // temporary: calculate llk_null
  Vector v(ndims+1);
  v[0] = pooled_af * 2.0;
  for(int32_t j=0; j < ndims; ++j)
    v[j+1] = betas[j];
  llknull = 0 - Evaluate(v);

  return true;
}

/*
bool frequency_estimator::lr_test_hwe(bool use_isaf) {
  estimate_isaf_em_hwd(); // estimate pooled allele frequency first

  Vector v00(ndims+1); // null hypothesis   : pooled AF with HWE
  Vector v01(ndims+2); // alt1 hypohthesis  : pooled AF with HWD (theta)
  Vector v10(ndims+1); // alt2 hypothesis   : is     AF with HWE
  Vector v11(ndims+2); // alt2 hypothesis   : is     AF with HWD (theta)

  v00[0] = v01[0] = v10[0] = v11[0] = pooled_af * 2.0;

  if ( theta < -1 ) { theta = -0.99999; }
  else if ( theta > 1 ) { theta = 0.99999; }  

  for(int32_t i=0; i < ndims; ++i) {
    v00[i] = v01[i] = 0;
    v10[i] = v11[i] = betas[i];
  }

  double llk00, llk01, llk10, llk11;
  
  assumeHWD = false;
  llk00 = 0-Evaluate(v00);
  llk10 = 0-Evaluate(v10);

  v10[ndims+1] = v11[ndims+1] = 0.5 * (log(1.0+theta) - log(1.0-theta));
  assumeHWD = true;
  llk01 = 0-Evaluate(v01);
  llk11 = 0-Evaluate(v11);  

  ibc0 = theta; 
  hwe0z = (llk01 < llk00) ? 0 : ( ( (ibc0 < 0) ? -1.0 : 1.0 ) * sqrt(2*(llk01 - llk00)) );  

  ibc1 = theta;
  hwe1z = (llk11 < llk10) ? 0 : ( ( (ibc1 < 0) ? -1.0 : 1.0 ) * sqrt(2*(llk11 - llk10)) );

  isaf_computed = true;

  return true;
}
*/


double frequency_estimator::Evaluate(Vector& v) {
  double llk = 0;
  if ( ndims+1+(assumeHWD ? 1 : 0) != v.dim )
    error("[E:%s:%d %s] Dimensions do not match %d vs %d", __FILE__, __LINE__, __PRETTY_FUNCTION__, ndims+1+(assumeHWD ? 1 : 0), v.dim);

  double expGeno, isaf, isafQ, isafRR, isafRA, isafAA, offset, h;

  if ( assumeHWD ) h = tanh(v[ndims+1]);
  else h = 0;

  for(int32_t i=0; i < nsamples; ++i) {
    expGeno = v[0];
    for(int32_t j=0; j < ndims; ++j)
      expGeno += (pEigVecs->operator()(i,j) * v[j+1]);
      
    isaf = expGeno/2.0;
    if ( isaf*(nsamples+nsamples+1) < 0.5 ) isaf = 0.5/(nsamples+nsamples+1);
    if ( (1.0-isaf)*(nsamples+nsamples+1) < 0.5 ) isaf = 1.0-0.5/(nsamples+nsamples+1);
    isafQ = 1.0-isaf;

    if ( assumeHWD ) {
      offset = isaf * isafQ * h;
      isafRR = isafQ * isafQ + offset;
      isafRA = 2 * (1.0-h) * isaf * isafQ;
      isafAA = isaf * isaf + offset;

      if ( isafRR < 0.25/(nsamples+nsamples+1)/(nsamples+nsamples+1) ) {
	isafRR = 0.25/(nsamples+nsamples+1)/(nsamples+nsamples+1);
	offset = isafRR - isafQ * isafQ;
	isafAA = isaf * isaf + offset;
	isafRA = 2 * isaf * isafQ - 2 * offset;
      }
      else if ( isafAA < 0.25/(nsamples+nsamples+1)/(nsamples+nsamples+1) ) {
	isafAA = 0.25/(nsamples+nsamples+1)/(nsamples+nsamples+1);
	offset = isafAA - isaf * isaf;
	isafRR = isafQ * isafQ + offset;
	isafRA = 2 * isaf * isafQ - 2 * offset;	
      }
    }
    else {
      isafRR = isafQ * isafQ;
      isafRA = 2 * isaf * isafQ;
      isafAA = isaf * isaf;
    }

    ifs[i] = isaf;
    if ( ploidies[i] == 2 ) {
      llk += log(isafRR    * phredConv.toProb(pls[i*3]) +
		 isafRA    * phredConv.toProb(pls[i*3+1]) +
		 isafAA    * phredConv.toProb(pls[i*3+2]));
    }
    else if ( ploidies[i] == 1 ) {
      llk += log(isafQ * phredConv.toProb(pls[i*3]) +
		 isaf  * phredConv.toProb(pls[i*3+2]));
    }

  }
  return 0-llk;  
}

void frequency_estimator::estimate_isaf_em(int32_t maxiter) {
  //Eigen::MatrixXd eV = Eigen::MatrixXd::Constant(nsamples,ndims+1,1.0);
  //eV.block(0,1,nsamples,ndims) = *pEigVecs;
  //Eigen::BDCSVD<Eigen::MatrixXd> svd(eV, Eigen::ComputeThinU | Eigen::ComputeThinV);
  if ( !isaf_computed ) {
    estimate_pooled_af_em();     
    Eigen::VectorXd y(nsamples);
    Eigen::VectorXd isaf = Eigen::VectorXd::Constant(nsamples, pooled_af);

    double maf = pooled_af > 0.5 ? 1-pooled_af : pooled_af;

    //Eigen::VectorXd lambda = Eigen::VectorXd::Constant(ndims, maxLambda * (1.-maf) / (maf * nsamples * 2));
    double lambda = maxLambda * (1.-maf) / (maf * nsamples * 2.0);
    double p0, p1, p2;

    for(int32_t i=0; i < maxiter; ++i) { // maxiter = 30
      for(int32_t j=0; j < nsamples; ++j) {
	if ( ploidies[j] == 2 ) {
	  p0 = ( 1.0 - isaf[j] ) * ( 1.0 - isaf[j] ) * phredConv.toProb(pls[3*j]);
	  p1 = 2.0 * isaf[j] * ( 1.0 - isaf[j] ) * phredConv.toProb(pls[3*j+1]);
	  p2 = isaf[j] * isaf[j] * phredConv.toProb(pls[3*j+2]);
	  y[j] = (p1+p2+p2+1e-100)/(p0+p1+p2+1e-100);
	}
	else if ( ploidies[j] == 1 ) {
	  p0 = ( 1.0 - isaf[j] ) * phredConv.toProb(pls[3*j]);
	  p2 = isaf[j] * phredConv.toProb(pls[3*j+2]);
	  y[j] = (p2+p2+1e-100)/(p0+p2+1e-100);
	}
	else {
	  y[j] = isaf[j] * 2;
	}
      }
      // U diag(d_i^2/(d_i^2+lambda)) U'y
      Eigen::VectorXd d2 = pSVD->singularValues();
      Eigen::MatrixXd UD2 = pSVD->matrixU(); //
      for(int32_t j=0; j < nsamples; ++j) {
	for(int32_t k=0; k < ndims; ++k) {
	  //UD2(j,k) *= ( d2[k] / ( d2[k] + lambda ) );
	  UD2(j,k) *= ( d2[k] * d2[k] / ( d2[k] + lambda ) / ( d2[k] + lambda) );
	}
      }
      
      isaf = UD2 * ( pSVD->matrixU().transpose() * y ) / 2.0;
      for(int32_t j=0; j < nsamples; ++j) {
	if ( isaf[j]*(nsamples+nsamples+1) < 0.5 ) isaf[j] = 0.5/(nsamples+nsamples+1.0);
	else if ( (1.0-isaf[j])*(nsamples+nsamples+1) < 0.5 ) isaf[j] = 1.0-0.5/(nsamples+nsamples+1.0);
      }
    }

    for(int32_t j=0; j < nsamples; ++j) {
      ifs[j] = (float)isaf[j];
    }

    Eigen::VectorXd d2 = pSVD->singularValues();    
    Eigen::MatrixXd VD = pSVD->matrixV();
    for(int32_t j=0; j < ndims; ++j) {
      for(int32_t k=0; k < ndims; ++k) {
	VD(j,k) *= ( d2[k] / ( d2[k] + lambda ) / ( d2[k] + lambda) );	
      }
    }
    Eigen::VectorXd vBeta = VD * ( pSVD->matrixU().transpose() * y ) / 2.0;
    for(int32_t k=0; k < ndims; ++k)
      betas[k] = (float)vBeta[k];
    
    isaf_computed = true;
  }
  
  //return 0;
}

void frequency_estimator::estimate_isaf_em_hwd(int32_t maxiter) {
  if ( !isaf_computed ) {
    estimate_pooled_af_em();     
    Eigen::VectorXd yE(nsamples);
    Eigen::VectorXd yD(nsamples);    
    Eigen::VectorXd isafE = Eigen::VectorXd::Constant(nsamples, pooled_af);
    Eigen::VectorXd isafD = Eigen::VectorXd::Constant(nsamples, pooled_af);    

    double maf = pooled_af > 0.5 ? 1-pooled_af : pooled_af;
    double lambda = maxLambda * (1.-maf) / (maf * nsamples * 2.0);
    double p0E, p1E, p2E; // frequencies under HWE - priors
    double p0D, p1D, p2D; // frequencies under HWD - priors
    double x0E, x1E, x2E; // frequencies under HWE - posteriors
    double x0D, x1D, x2D; // frequencies under HWD - posteriors

    //notice("pooled_af = %lf, maf = %lf", pooled_af, maf);    

    double minAF = 0.5/(nsamples+nsamples+1.);
    double minGF = minAF*minAF;

    double thetaDouble = 0;

    // U diag(d_i^2/(d_i^2+lambda)) U'y
    Eigen::VectorXd d2 = pSVD->singularValues();
    Eigen::MatrixXd UD2 = pSVD->matrixU(); //
    for(int32_t j=0; j < nsamples; ++j) {
      for(int32_t k=0; k < ndims; ++k) {
	//UD2(j,k) *= ( d2[k] / ( d2[k] + lambda ) );
	UD2(j,k) *= ( d2[k] * d2[k] / ( d2[k] + lambda ) / ( d2[k] + lambda) );
      }
    }

    Eigen::MatrixXd VD = pSVD->matrixV();
    for(int32_t j=0; j < ndims; ++j) {
      for(int32_t k=0; k < ndims; ++k) {
	VD(j,k) *= ( d2[k] / ( d2[k] + lambda ) / ( d2[k] + lambda) );	
      }
    }

    double l0, l1, l2;    
    for(int32_t i=0; i < maxiter; ++i) { // maxiter = 30
      double thetaNum = 1e-100;
      double thetaDen = 1e-100;	
             
      for(int32_t j=0; j < nsamples; ++j) {
	if ( ploidies[j] == 2 ) {
	  // calculate isaf under HWD
	  p0E = ( 1.0 - isafD[j] ) * ( 1.0 - isafD[j] ); // * phredConv.toProb(pls[3*j]);
	  p1E = 2.0 * isafD[j] * ( 1.0 - isafD[j] ); // * phredConv.toProb(pls[3*j+1]);
	  p2E = isafD[j] * isafD[j]; // * phredConv.toProb(pls[3*j+2]);
	  p0D = p0E + 0.5 * thetaDouble * p1E;
	  p1D = p1E * (1.-thetaDouble);
	  p2D = p2E + 0.5 * thetaDouble * p1E;
	  if ( p0D < minGF ) { // boundary condition
	    p0D = minGF;                  // K = -p0E + minGF
	    p1D = p1E + 2.0*p0E - 2.0*minGF;  // p1D = p1E - 2K = p1E + 2p0E - 2minGF
	    p2D = p2E - p0E + minGF;      // p2D = p2E + K = p2E - p0E + minGF
	  }
	  else if ( p2D < minGF) {
	    p2D = minGF;                  // K = -p0E + minGF
	    p1D = p1E + 2.0*p2E - 2.0*minGF;  // p1D = p1E - 2K = p1E + 2p0E - 2minGF
	    p0D = p0E - p2E + minGF;      // p2D = p2E + K = p2E - p0E + minGF	    
	  }
	  // calculate isaf under HWE	  
	  p0E = ( 1.0 - isafE[j] ) * ( 1.0 - isafE[j] ); // * phredConv.toProb(pls[3*j]);
	  p1E = 2.0 * isafE[j] * ( 1.0 - isafE[j] ); // * phredConv.toProb(pls[3*j+1]);
	  p2E = isafE[j] * isafE[j]; // * phredConv.toProb(pls[3*j+2]);	  

	  l0 = phredConv.toProb(pls[3*j+0]);
	  l1 = phredConv.toProb(pls[3*j+1]);
	  l2 = phredConv.toProb(pls[3*j+2]);	  
	  
	  x0E = p0E * l0;
	  x1E = p1E * l1; 
	  x2E = p2E * l2; 

	  x0D = p0D * l0; 
	  x1D = p1D * l1; 
	  x2D = p2D * l2;

	  yE[j] = (x1E+x2E+x2E+1e-100)/(x0E+x1E+x2E+1e-100);
	  yD[j] = (x1D+x2D+x2D+1e-100)/(x0D+x1D+x2D+1e-100);

	  //double w = fabs(l0 + l2 - l1 - l1); // heuristic for now.
	  
	  thetaNum += (x1D / (x0D + x1D + x2D) );
	  thetaDen += (p1E);
	}
	else if ( ploidies[j] == 1 ) {
	  x0E = ( 1.0 - isafE[j] ) * phredConv.toProb(pls[3*j]);
	  x2E = isafE[j] * phredConv.toProb(pls[3*j+2]);
	  yE[j] = (x2E+x2E+1e-100)/(x0E+x2E+1e-100);
	  
	  x0D = ( 1.0 - isafD[j] ) * phredConv.toProb(pls[3*j]);
	  x2D = isafD[j] * phredConv.toProb(pls[3*j+2]);
	  yD[j] = (x2D+x2D+1e-100)/(x0D+x2D+1e-100);	  
	}
	else { // this might be buggy for multi-allelics
	  yE[j] = isafE[j] * 2.0; // use expectation
	  yD[j] = isafD[j] * 2.0; // use expectation	  
	}
      }
      
      // calculate isaf (which already reflects beta)
      isafD = UD2 * ( pSVD->matrixU().transpose() * yD ) / 2.0;
      isafE = UD2 * ( pSVD->matrixU().transpose() * yE ) / 2.0;

      //notice("pooled_af = %lf, maf = %lf, minMAF = %lf", pooled_af, maf, minAF);
         
      for(int32_t j=0; j < nsamples; ++j) {
	if ( isafD[j] < minAF ) isafD[j] = minAF;
	else if ( 1.0-isafD[j] < minAF ) isafD[j] = 1.0-minAF;
	if ( isafE[j] < minAF ) isafE[j] = minAF;
	else if ( 1.0-isafE[j] < minAF ) isafE[j] = 1.0-minAF;	
      }

      thetaDouble = 1.0 - thetaNum/thetaDen;
    }

    for(int32_t j=0; j < nsamples; ++j) {
      ifs[j] = (float)isafD[j];
    }

    Eigen::VectorXd vBetaD = VD * ( pSVD->matrixU().transpose() * yD ) / 2.0;
    Eigen::VectorXd vBetaE = VD * ( pSVD->matrixU().transpose() * yE ) / 2.0;    
    for(int32_t k=0; k < ndims; ++k) {
      betas[k] = (float)vBetaD[k];
      //betasE[k] = (float)vBetaE[k];      
    }

    theta = (float)thetaDouble;

    // calculate likelihoods
    double llk00 = 0, llk01 = 0, llk10 = 0, llk11 = 0;
    double fE[3] = { (1.0-pooled_af)*(1.0-pooled_af),
		      2*pooled_af*(1.0-pooled_af),
		      pooled_af*pooled_af };
    double fD[3] = { fE[0] + 0.5 * thetaDouble * fE[1],
		     fE[1] * (1.0-thetaDouble),
		     fE[2] + 0.5 * thetaDouble * fE[1] };
    if ( fD[0] < minGF ) {
      fD[0] = minGF;
      fD[1] = fE[1] + 2*fE[0] - 2*minGF;
      fD[2] = fE[2] - fE[0] + minGF;      
    }
    else if ( fD[2] < minGF ) {
      fD[2] = minGF;
      fD[1] = fE[1] + 2*fE[2] - 2*minGF;
      fD[0] = fE[0] - fE[2] + minGF;            
    }

    for(int32_t j=0; j < nsamples; ++j) {
      if ( ploidies[j] == 2 ) {
	p0E = ( 1.0 - isafD[j] ) * ( 1.0 - isafD[j] ); // * phredConv.toProb(pls[3*j]);
	p1E = 2.0 * isafD[j] * ( 1.0 - isafD[j] ); // * phredConv.toProb(pls[3*j+1]);
	p2E = isafD[j] * isafD[j]; // * phredConv.toProb(pls[3*j+2]);
	
	p0D = p0E + 0.5 * thetaDouble * p1E;
	p1D = p1E * (1.0 - thetaDouble);
	p2D = p2E + 0.5 * thetaDouble * p1E;
	if ( p0D < minGF ) { // boundary condition
	  p0D = minGF;                  // K = -p0E + minGF
	  p1D = p1E + 2*p0E - 2*minGF;  // p1D = p1E - 2K = p1E + 2p0E - 2minGF
	  p2D = p2E - p0E + minGF;      // p2D = p2E + K = p2E - p0E + minGF
	}
	else if ( p2D < minGF ) {
	  p2D = minGF;                  // K = -p0E + minGF
	  p1D = p1E + 2*p2E - 2*minGF;  // p1D = p1E - 2K = p1E + 2p0E - 2minGF
	  p0D = p0E - p2E + minGF;      // p2D = p2E + K = p2E - p0E + minGF	    
	}

	p0E = ( 1.0 - isafE[j] ) * ( 1.0 - isafE[j] ); // * phredConv.toProb(pls[3*j]);
	p1E = 2 * isafE[j] * ( 1.0 - isafE[j] ); // * phredConv.toProb(pls[3*j+1]);
	p2E = isafE[j] * isafE[j]; // * phredConv.toProb(pls[3*j+2]);	

	l0 = phredConv.toProb(pls[3*j+0]);
	l1 = phredConv.toProb(pls[3*j+1]);
	l2 = phredConv.toProb(pls[3*j+2]);

	llk00 += log(fE[0] * l0 + fE[1] * l1 + fE[2] * l2);
	llk01 += log(fD[0] * l0 + fD[1] * l1 + fD[2] * l2);
	llk10 += log(p0E   * l0 + p1E   * l1 + p2E   * l2);
	llk11 += log(p0D   * l0 + p1D   * l1 + p2D   * l2);		
      }
      else if ( ploidies[j] == 1 ) {
	p0E = ( 1.0 - isafE[j] );
	p2E = isafE[j];
	p0D = ( 1.0 - isafD[j] );
	p2D = isafD[j];	
	l0 = phredConv.toProb(pls[3*j+0]);
	l2 = phredConv.toProb(pls[3*j+2]);

	llk00 += log( (1.0-pooled_af) * l0 + pooled_af * l2);
	llk01 += log( (1.0-pooled_af) * l0 + pooled_af * l2);
	llk10 += log( p0E * l0 + p2E * l2);
	llk11 += log( p0D * l0 + p2D * l2);
      }
    }

    ibc0 = ibc1 = theta;
    hwe0z = (llk01 < llk00) ? 0 : ( ( (ibc0 < 0) ? -1.0 : 1.0 ) * sqrt(2*(llk01 - llk00)) );
    hwe1z = (llk11 < llk10) ? 0 : ( ( (ibc1 < 0) ? -1.0 : 1.0 ) * sqrt(2*(llk11 - llk10)) );    
    
    isaf_computed = true;
  }
  
  //return 0;
}

void frequency_estimator::estimate_isaf_simplex() {
  if ( isaf_computed ) return;

  assumeHWD = false;
    
  AmoebaMinimizer isafMinimizer;
  Vector startingPoint(ndims+1);
  double emaf = estimate_pooled_af_em();
  
  startingPoint.Zero();
  startingPoint[0] = emaf*2.0;

  isafMinimizer.func = this;

  isafMinimizer.Reset(ndims+1);
  isafMinimizer.point = startingPoint;
  isafMinimizer.Minimize(tol);
  Evaluate(isafMinimizer.point);          

  llknull = 0 - isafMinimizer.fmin;

  //score_test_hwe(true, emaf);

  isaf_computed = true;  

  //iter = isafMinimizer.cycleCount;
  //return 0-isafMinimizer.fmin;
  //return ancMinimizer.fmin;
}

void frequency_estimator::estimate_isaf_lrt() {
  if ( isaf_computed ) return;

  double emaf = estimate_pooled_af_em();    
  double llk0, llk1;
  std::vector<double> p0(ndims+1);
  std::vector<double> p1(ndims+1);  

  // Find MLE assuming HWE
  {
    assumeHWD = false;
    
    AmoebaMinimizer isafMinimizer;
    Vector startingPoint(ndims+1);
    
    startingPoint.Zero();
    startingPoint[0] = emaf*2.0;
    
    isafMinimizer.func = this;
    
    isafMinimizer.Reset(ndims+1);
    isafMinimizer.point = startingPoint;
    isafMinimizer.Minimize(tol);
    Evaluate(isafMinimizer.point);        

    for(int i=0; i < ndims+1; ++i)
      p0[i] = isafMinimizer.point[i];
    llknull = llk0 = 0 - isafMinimizer.fmin;

    //notice("ndims = %d, p0[0] = %.5lg, p0[1] = %.5lg, p0[2] = %.5lg, p0[3] = %.5lg", ndims, p0[0], p0[1], p0[2], p0[3]);
  }

  // Find MLE without assuming HWE  
  {
    assumeHWD = true;
    AmoebaMinimizer isafMinimizer;
    Vector startingPoint(ndims+2);
    
    startingPoint.Zero();

    for(int i=0; i < ndims+1; ++i)
      startingPoint[i] = p0[i];
    
    isafMinimizer.func = this;
    
    isafMinimizer.Reset(ndims+2);
    isafMinimizer.point = startingPoint;
    isafMinimizer.Minimize(tol);
    Evaluate(isafMinimizer.point);    

    for(int i=0; i < ndims+1; ++i)
      p1[i] = isafMinimizer.point[i];
    
    theta = tanh(isafMinimizer.point[ndims+1]);
    llk1 = 0 - isafMinimizer.fmin;

    //notice("ndims = %d, p1[0] = %.5lg, p1[1] = %.5lg, p1[2] = %.5lg, p1[3] = %.5lg", ndims, p1[0], p1[1], p1[2], p1[3]);    
  }


  double pp0 = (1.-emaf)*(1.-emaf);
  double pp1 = 2 * emaf * (1-emaf);
  double pp2 = emaf * emaf;
  double obsH0 = 0, obsH1 = 0, expH1 = 0;
  double l0, l1, l2, sumU0, sqU0, sum1, sum0, ph1, ph0, U0;
  int32_t ndiploids = 0;

  // pretend that we have pp0, pp1, pp2 observations of each genotype (pseudocount)
  sumU0 = pp0*(1/pp0) + pp1 * (-2/pp1) + pp2*(1/pp2);
  sqU0 =  pp0*1/pp0/pp0 + pp1 * 4/pp1/pp1 + pp2*(1/pp2/pp2);
  
  //int32_t ndiploids = 0;
  
  for(int32_t i=0; i < nsamples; ++i) {
    if ( ploidies[i] != 2 ) continue;
    ++ndiploids;
    l0 = phredConv.toProb(pls[i*3]);
    l1 = phredConv.toProb(pls[i*3+1]);
    l2 = phredConv.toProb(pls[i*3+2]);

    sum0 = l0*pp0 + l1*pp1 + l2*pp2 + 1e-100;
    ph0 = l1*pp1;
    obsH0 += (ph0/sum0);
    U0 = emaf*(1-emaf)*(l0-2*l1+l2)/sum0;
    sumU0 += U0;
    sqU0 += (U0*U0);
    
    sum1 = l0*(1.-ifs[i])*(1.-ifs[i]) + 2*l1*(1.-ifs[i])*ifs[i] + l2*ifs[i]*ifs[i] + 1e-100;
    ph1 = 2*l1*(1.-ifs[i])*ifs[i];
    obsH1 += (ph1/sum1);
    expH1 += (2*(1.-ifs[i])*ifs[i]);
  }

  hwe0z = sumU0/sqrt(sqU0);
  ibc0 = 1.0 - (obsH0+1)/(pp1*ndiploids+1);

  //hwe1z = sumU1/sqrt(sqU1);
  ibc1 = 1.0 - (obsH1+1)/(expH1+1);
  hwe1z = (llk1 < llk0) ? 0 : ( ( (ibc1 < 0) ? -1.0 : 1.0 ) * sqrt(2*(llk1 - llk0)) );

  isaf_computed = true;
}

bool frequency_estimator::update_gt_gq(bool update_gq) {
  if ( siteOnly ) return false;
  double gp = 0, gp_sum = 0, max_gp = 0;
  int32_t best_gt = 0;
  int32_t best_a1 = 0, best_a2 = 0;
  int32_t an = 0;
  int32_t acs[2] = {0,0};
  int32_t gcs[3] = {0,0,0};
  float afs[3];
  int32_t max_gq = 0;

  if ( gt == NULL ) 
    gt = (int32_t*) malloc(sizeof(int32_t)*2*nsamples);
  if ( ( update_gq ) && ( gq == NULL ) ) 
    gq = (int32_t*) malloc(sizeof(int32_t)*nsamples);    
  
  for(int32_t i=0; i < nsamples; ++i) {
    int32_t* pli = &pls[ i * 3 ];
    
    if ( ploidies[i] == 1 ) {
      max_gp = gp_sum = gp = ( phredConv.toProb((uint32_t)pli[0]) * (1.0 - ifs[i]) );
      best_gt = 0; best_a1 = 0; best_a2 = 0;
      gp = ( phredConv.toProb((uint32_t)pli[2]) * ifs[i] );
      gp_sum += gp;
      if ( max_gp < gp ) {
	max_gp = gp;
	best_gt = 2; best_a1 = 1; best_a2 = 1;
      }      
    }
    else if ( ploidies[i] == 2 ) {
      max_gp = gp_sum = gp = ( phredConv.toProb((uint32_t)pli[0]) * (1.0-ifs[i]) * (1.0-ifs[i]) );
      best_gt = 0; best_a1 = 0; best_a2 = 0;

      gp = phredConv.toProb((uint32_t)pli[1]) * 2.0 * ifs[i] * (1.0-ifs[i]);
      gp_sum += gp;
      if ( max_gp < gp ) { max_gp = gp; best_gt = 1; best_a1 = 0; best_a2 = 1; }

      gp = phredConv.toProb((uint32_t)pli[2]) * ifs[i] * ifs[i];
      gp_sum += gp;
      if ( max_gp < gp ) { max_gp = gp; best_gt = 2; best_a1 = 1; best_a2 = 1; }      
    }
    else if ( ploidies[i] == 0 ) {
      best_gt = 0;
      max_gp = 0;
      gp_sum = 1e-100;      
    }
    else
      error("[E:%s:%d %s] Unexpected ploidy %d", __FILE__, __LINE__, __PRETTY_FUNCTION__, (int32_t)ploidies[i]);

    if ( update_gq ) {
      double prob = 1.-max_gp/gp_sum;  // to calculate GQ
      if ( prob <= 3.162278e-26 )
	prob = 3.162278e-26;
      if ( prob > 1 )
	prob = 1;
    
      gq[i] = (int32_t)phredConv.err2Phred((double)prob);
      if ( ( best_gt > 0 ) && ( max_gq < gq[i] ) ) {
	max_gq = gq[i];
      }
    }
    
    gt[2*i]   = ((best_a1 + 1) << 1);
    gt[2*i+1] = ((best_a2 + 1) << 1);	    
    an += 2;             // still use diploid representation of chrX for now.
    ++acs[best_a1];
    ++acs[best_a2];
    ++gcs[best_gt];    
  }

  for(size_t i=0; i < 2; ++i) {
    afs[i] = acs[i]/(float)an;
  }

  //notice("Calling bcf_update_format_int32() with nsamples=%d",nsamples);

  bcf_update_format_int32(hdr, iv, "GT", gt, nsamples * 2);
  if ( update_gq )
    bcf_update_format_int32(hdr, iv, "GQ", gq, nsamples );	  

  iv->qual = (float)max_gq;

  bcf_update_info_int32(hdr, iv, "AC", &acs[1], 1);
  bcf_update_info_int32(hdr, iv, "AN", &an, 1);
  bcf_update_info_float(hdr, iv, "AF", &afs[1], 1);
  bcf_update_info_int32(hdr, iv, "GC", gcs, 3);
  bcf_update_info_int32(hdr, iv, "GN", &nsamples, 1);

  return true;
}
 
bool frequency_estimator::update_variant() {
  float hweslp0 = (float)((hwe0z > 0 ? -1 : 1) * log10( erfc(fabs(hwe0z)/sqrt(2.0)) + 1e-100 ));
  float hweslp1 = (float)((hwe1z > 0 ? -1 : 1) * log10( erfc(fabs(hwe1z)/sqrt(2.0)) + 1e-100 ));
  float max_if = 0, min_if = 1;
  for(int32_t j=0; j < nsamples; ++j) {
    if ( ifs[j] > max_if ) max_if = ifs[j];
    if ( ifs[j] < min_if ) min_if = ifs[j];
  }

  if ( siteOnly ) {
    //notice("foo");
    bcf_subset(hdr, iv, 0, 0);
    //notice("goo");    
  }  

  if ( !skipInfo ) {
    float hweaf = (float)pooled_af;
    bcf_update_info_float(wdr, iv, "HWEAF_P", &hweaf, 1);
    bcf_update_info_float(wdr, iv, "FIBC_P", &ibc0, 1);    
    bcf_update_info_float(wdr, iv, "HWE_SLP_P", &hweslp0, 1);
    bcf_update_info_float(wdr, iv, "FIBC_I", &ibc1, 1);    
    bcf_update_info_float(wdr, iv, "HWE_SLP_I", &hweslp1, 1);
    bcf_update_info_float(wdr, iv, "MAX_IF", &max_if, 1);
    bcf_update_info_float(wdr, iv, "MIN_IF", &min_if, 1);
    bcf_update_info_float(wdr, iv, "LLK0", &llknull, 1);    
    bcf_update_info_float(wdr, iv, "BETA_IF", betas, ndims);
  }
  if ( ( !skipIf ) && ( !siteOnly ) ) {
    bcf_update_format_float(wdr, iv, "IF", ifs, nsamples);
  }
  
  return true;
}
