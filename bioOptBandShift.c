/*! Bio-Optical Band Shifting. 
 * Algorithm to shift bands to conform to OCR-VC band set; 
 * 412, 443, 490, 510, 555, 670.
 * Band shifting based on QAA, produces data for bands 
 * missing from the ocr-vc band set.
 * A sensor-dependent "conversion context" structure is developed 
 * in the "setup" function and is referred to throughout the program
 * for the various calculations that need to take place.
 * Assignment of Bricaud and pure water coefficients are built-in to 
 * avoid out-of-program calls within the "line & pixel loops"     
 */

// E. Karakoylu 140710
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"


#define REF_WVL 443
#define QAA_MIN 0
#define QAA_MAX 5
#define WVL_NUM 6

typedef struct _context{
    int *wvl_i,*wvl_o,*wvl_all,wvl_sms;
    float *aw_all,*bbw_all;
    float *aw_i,*bbw_i,*a_i,*b_i;
    float *aw_o,*bbw_o,*a_o,*b_o;
    float aph_ref,adg_ref,bbp_ref;//QAA stuff
    float a_sms,b_sms;
    int nBands;
    int *ibvc_bin;
    float greenBand;
        
} ccstr;

void Setup(l2str*, ccstr*);
void GetBricaudAB(ccstr*);
void GetAwBw(ccstr*);
void QaaCalc(ccstr*,l2str*,int32_t);
float BandShift(ccstr*,l2str*,int32_t, int);
float BandInterp(ccstr*,float*,int,int);
float ToBelowWater(float);
float ToAboveWater(float);
void CheckNULL(void*,char*,int,char*);

void Deallocate(ccstr*); // currently not implemented
static int current_iscan = -1;
static int current_npix = -1;
static float *sh_prod_arr;
static int nband_vc = 7;
static int vc_bands[7] = {412,443,490,510,531,555,670};

void bioOptBandShift(l2str *l2rec, l2prodstr *p, float prod[]){
	
	int32_t ibvc = p->prod_ix;
	int32_t ib,ip,ipb,col;
	static ccstr vcc;
	static int firstRun = 1;
    float temp; // DELETE THIS 
	if (firstRun){
		Setup(l2rec, &vcc);
		firstRun = 0;
	}
    /* sh_prod_arr is a 2D array accessed as a simple pointer of dims npix rows, and nband_vc cols.
       If band shift not necessary, correspoding cols filled w/ Rrs already in l2rec. */
    /* occasional segment below -- happens only when switching to new row 
       but avoid recalculation for multiple products */
	if (l2rec->iscan != current_iscan){
        current_iscan = l2rec->iscan;
	    current_npix = l2rec->npix;
        if (firstRun)
            firstRun=0;
            
        if (sh_prod_arr != NULL)
            free(sh_prod_arr);
	    CheckNULL(sh_prod_arr = (float*)malloc(current_npix * nband_vc* sizeof(float))
                    ,"sh_prod_arr",__LINE__,__FILE__);
        for (ip=0;ip<current_npix;ip++){
		    QaaCalc(&vcc,l2rec,ip); // calculate some reference iops for bandshifting, for each pix.
            for (ib=0;ib<nband_vc;ib++){
                // fill sh_prod_arr, by bandshifting if necessary, or by copying l2rec rrs... 
                if (vcc.ibvc_bin[ib]){
                   *(sh_prod_arr + ip * nband_vc + ib) = BandShift(&vcc,l2rec,ip,vc_bands[ib]);
                }
                else { 
                    
                    ipb = ip * l2rec->nbands + bindex_get(vc_bands[ib]);
                    *(sh_prod_arr +ip * nband_vc + ib) = l2rec->Rrs[ipb];
                }
            }
        }
	}
	for (ip=0;ip<l2rec->npix; ip++){
            prod[ip] = *(sh_prod_arr + ip * nband_vc + ibvc); 
	}
    
}

void Setup(l2str *l2rec, ccstr *vccp){
	int kb,iflag;
	
	int rrs_wvl[4][WVL_NUM] = {
							{412,443,490,510,555,670},//SWF
							{412,443,488,531,547,667},//MODA
							{412,443,490,510,560,665},//MER
							{412,445,488,-1,555,672}  //VIIRS
							};
	
	int lambda_i[4][5] = {
							{510,555,-1,-1,-1},   //SWF
							{488,488,531,547,667},//MODA
							{510,560,560,665,-1}, //MER
							{488,488,555,555,670} //VIIRS
							};
	
	int lambda_o[4][5] = {
							{531,531,-1,-1,-1},	  //SWF
							{490,510,510,555,670},//MODA
							{531,531,555,670,-1}, //MER
							{490,510,510,531,670} //VIIRS
							};
	
	int ibvc_bin[3][7] = { //index of band needing bandshift for each sensor
							{0,0,0,0,1,0,0},	//SWF
							{0,0,1,1,0,1,1},	//MODA
							{0,0,0,0,1,1,1}		//MER
							};	

	switch (l2rec->sensorID){
		case SEAWIFS:{
			vccp->greenBand = 555;
			vccp->nBands = 2; // 2 input bands will be used to compute 531.
			iflag = 0;
			break;
		}

		case HMODISA:{
			vccp->greenBand = 547;
			vccp->nBands = 5; //inp bands: 488,488,531,547,667, resp
			iflag = 1;
			break;
		}

		case MERIS:{
			vccp->greenBand = 560;
			vccp->nBands = 5;
			iflag = 2;
			break;
		}
	}

	CheckNULL(vccp->ibvc_bin = (int*)malloc(nband_vc * sizeof(int)),"vccp->ibvc_bin",__LINE__,__FILE__);
	CheckNULL(vccp->wvl_i = (int*)malloc(vccp->nBands * sizeof(int)),"vccp->wvl_i",__LINE__,__FILE__);
    CheckNULL(vccp->wvl_o = (int*)malloc(vccp->nBands * sizeof(int)),"vccp->wvl_o",__LINE__,__FILE__);
    CheckNULL(vccp->wvl_all = (int*)malloc(WVL_NUM * sizeof(int)),"vccp->wvl_all",__LINE__,__FILE__);
    CheckNULL(vccp->aw_i = (float*)malloc(vccp->nBands * sizeof(float)),"vccp->aw_i",__LINE__,__FILE__);
	CheckNULL(vccp->aw_o = (float*)malloc(vccp->nBands * sizeof(float)),"vccp->aw_o",__LINE__,__FILE__);	
    CheckNULL(vccp->aw_all = (float*)malloc(WVL_NUM * sizeof(float)),"vccp->aw_all",__LINE__,__FILE__);
    CheckNULL(vccp->bbw_i = (float*)malloc(vccp->nBands * sizeof(float)),"vccp->aw_all",__LINE__,__FILE__);
    CheckNULL(vccp->bbw_o = (float*)malloc(vccp->nBands * sizeof(float)),"vccp->bbw_o",__LINE__,__FILE__);
	CheckNULL(vccp->bbw_all = (float*)malloc(WVL_NUM * sizeof(float)),"vccp->bbw_all",__LINE__,__FILE__);
    CheckNULL(vccp->a_i = (float*)malloc(vccp->nBands * sizeof(float)),"vccp->a_i",__LINE__,__FILE__);
    CheckNULL(vccp->a_o = (float*)malloc(vccp->nBands * sizeof(float)),"vccp->a_o",__LINE__,__FILE__);
    CheckNULL(vccp->b_i = (float*)malloc(vccp->nBands * sizeof(float)),"vccp->b_i",__LINE__,__FILE__);
	CheckNULL(vccp->b_o = (float*)malloc(vccp->nBands * sizeof(float)),"vccp->b_o",__LINE__,__FILE__);
	
	for(kb=0;kb<vccp->nBands;kb++){
	    vccp->wvl_i[kb] = lambda_i[iflag][kb];
		vccp->wvl_o[kb] = lambda_o[iflag][kb];
    }
	for(kb=0;kb<WVL_NUM;kb++)
		vccp->wvl_all[kb] = rrs_wvl[iflag][kb];
/*
	for(kb=0;kb<vccp->shWvlNum;kb++)
		vccp->shWvl[kb] = uniq_wvl_o[iflag][kb];
*/	
	for(kb=0;kb<nband_vc;kb++)
		vccp->ibvc_bin[kb] = ibvc_bin[iflag][kb];
	
	GetBricaudAB(vccp);
	GetAwBw(vccp);
	// fill up Bricaud and Pure water coefficient arrays...
}

void GetBricaudAB(ccstr *ccp){
	
	float bricaud_L_A_B[][13] = 
			{
				{412,413,443,488,490,510,531,547,555,560,665,667,670},
				{0.0323,0.032775,0.0394,0.0279,0.0274,0.018,0.0115,0.00845,
									0.007,0.0062,0.0152,0.01685,0.0189},
				{0.286,0.28775,0.3435,0.369,0.361,0.260,0.134,0.0625,
									0.0315,0.016,0.134,0.140,0.149}
			};

	int k_i = 0,k_o = 0;
	int bi,num_wvl = 13;
	
	for (bi=0;bi<num_wvl;bi++){
		
		if (ccp->wvl_i[k_i] < bricaud_L_A_B[0][bi]){ // for repeat wvl_i
			ccp->a_i[k_i] = ccp->a_i[k_i-1];
			ccp->b_i[k_i] = ccp->b_i[k_i-1];
			k_i++;
		}
		
		if (ccp->wvl_i[k_i]  == bricaud_L_A_B[0][bi]){
			ccp->a_i[k_i] = bricaud_L_A_B[1][bi];
			ccp->b_i[k_i] = bricaud_L_A_B[2][bi];
			k_i++;
		}
		
		if (ccp->wvl_o[k_o] < bricaud_L_A_B[0][bi]){// for repeat wvl_o
			ccp->a_o[k_o] = ccp->a_o[k_o - 1];
			ccp->b_o[k_o] = ccp->b_o[k_o - 1];
			k_o++;
		}

		if (ccp->wvl_o[k_o] == bricaud_L_A_B[0][bi]){
			ccp->a_o[k_o] = bricaud_L_A_B[1][bi];
			ccp->b_o[k_o] = bricaud_L_A_B[2][bi];
			k_o++;
		}
		
		if (bricaud_L_A_B[0][bi] == REF_WVL){
			ccp->a_sms = bricaud_L_A_B[1][bi];
			ccp->b_sms = bricaud_L_A_B[2][bi];	
		}
	}
}

void GetAwBw(ccstr *ccp){

	float pwLAwBw[][13] = 
			{
				{412,413,443,488,490,510,531,547,555,560,665,667,670},
				{0.00455056,0.00449607,0.00706914,0.0145167,0.015,0.0325,
					0.0439153,0.0531686,0.0596,0.0619,0.429,0.434888,0.439},
				{0.00579201,0.00573196,0.00424592,0.00281659,0.00276835,0.00233870,
					0.00197385,0.00174280,0.00163999,0.00157958,0.000772104,
					0.000762543,0.000748479}
			};
	int k_i=0,k_o=0,kwl=0;
	int bi,num_wvl=13;

	for(bi=0;bi<num_wvl;bi++){

		if (ccp->wvl_i[k_i] < (int)pwLAwBw[0][bi]){ // for repeating wvl_i...
			ccp->aw_i[k_i] = ccp->aw_i[k_i - 1 ];
			ccp->bbw_i[k_i] = ccp->bbw_i[k_i - 1];
			k_i++;
		}

		if (ccp->wvl_i[k_i] == (int)pwLAwBw[0][bi]){
			ccp->aw_i[k_i] = pwLAwBw[1][bi];
			ccp->bbw_i[k_i] = pwLAwBw[2][bi];
			k_i++;
		}

		if (ccp->wvl_o[k_o] < (int)pwLAwBw[0][bi]){ // for repeating wvl_o...
			ccp->aw_o[k_o] = ccp->aw_o[k_o - 1];
			ccp->bbw_o[k_o] = ccp->bbw_o[k_o - 1];
			k_o++;
		}

		if (ccp->wvl_o[k_o] == (int)pwLAwBw[0][bi]){
			ccp->aw_o[k_o] = pwLAwBw[1][bi];
			ccp->bbw_o[k_o] = pwLAwBw[2][bi];
			k_o++;
		}

		if (ccp->wvl_all[kwl] == (int)pwLAwBw[0][bi]){
			ccp->aw_all[kwl] = pwLAwBw[1][bi];
			ccp->bbw_all[kwl] = pwLAwBw[2][bi];
			kwl++;
		}
	}
}

void QaaCalc(ccstr *vcc,l2str *l2rec,int32_t ip){
	int ib,ipb;
	float up_667,lw_667;
	float ab_surf_rrs[WVL_NUM], blw_surf_rrs[WVL_NUM], u_array[WVL_NUM];
	float c1 = 0.52,c2 = 1.7;
	float g0 = 0.089,g1=0.125;
	float a_ref,bbp_ref,bbp_exp,a_ref_exp;
	float ratio_aph,slope_adg,ratio_adg;
	float bbp[2],a[2];
	float aph_443,adg_443;
	
	for (ib=0;ib<WVL_NUM;ib++){
		ipb = ip * l2rec->nbands + bindex_get(vcc->wvl_all[ib]); 
		*(ab_surf_rrs + ib) = l2rec->Rrs[ipb];
		if (ib == 5){	
			// red band correction
			up_667 = 20 * pow( *(ab_surf_rrs+4), 1.5);
			lw_667 = 0.9 * pow( *(ab_surf_rrs+4), 1.7);
			if ( *(ab_surf_rrs+5) > up_667 || *(ab_surf_rrs+5) < lw_667){
				
				*(ab_surf_rrs+5) = 1.27 * pow( *(ab_surf_rrs+4), 1.47) + 
					0.00018 * pow( *(ab_surf_rrs+2) / *(ab_surf_rrs+4),
									-3.19);
 			
			}
		}
		//*(blw_surf_rrs+ib) = *(ab_surf_rrs+ib)/(c1 + c2 * *(ab_surf_rrs+ib));
		*(blw_surf_rrs+ib) = ToBelowWater(*(ab_surf_rrs+ib));
		// estimate bb = bb / (a+bb)
		*(u_array + ib) = (-g0 + pow(pow(g0,2) + 4 * g1 * 
										*(blw_surf_rrs+ib),0.5)) / (2 * g1);
	}
	// estimate a at ref. wavelength
	a_ref_exp =  log10(
						( *(blw_surf_rrs+1) + *(blw_surf_rrs+2) ) /
					( *(blw_surf_rrs+4) + 5 * ( *(blw_surf_rrs+5) /
													*(blw_surf_rrs+2) ) 
											* *(blw_surf_rrs+5) ) );
	
	/*a_ref = *(vcc->aw+4) + pow(10,(-1.146 + (-1.366 * a_ref_exp) + 
											(-0.469 * pow(a_ref_exp,2)	)	)	);*/
	a_ref = *(vcc->aw_all+4) + pow(10,(-1.146 + (-1.366 * a_ref_exp) + 
											(-0.469 * pow(a_ref_exp,2)	)	)	);
	//bbp_ref = ( *(u_array+4) * a_ref / (1 - *(u_array+4) ) ) - *(vcc->bbw+4);
	bbp_ref = ( *(u_array+4) * a_ref / (1 - *(u_array+4) ) ) - *(vcc->bbw_all+4);
	bbp_exp = 2 * (1 - 1.2 * exp(-0.9 * *(blw_surf_rrs + 1) / 
									*(blw_surf_rrs + 4) ) );
 	
    // estimate ratio of aph443/aph550
	ratio_aph = 0.74 + (0.2 / (0.8 + *(blw_surf_rrs + 1) /
								 *(blw_surf_rrs + 4) ) );
	// estimate ratio of adg443/adg550
	slope_adg = 0.015 + (0.002 / (0.6 + *(blw_surf_rrs + 1) / *(blw_surf_rrs + 4) ) );
	ratio_adg = exp(slope_adg * ( *(vcc->wvl_all+1) - *(vcc->wvl_all) ) );
	for(ib=0;ib<2;ib++){
		*(bbp + ib) = bbp_ref * pow ( ( vcc->greenBand / *(vcc->wvl_all+ib) ),bbp_exp );
		*(a + ib) = ( 1 - *(u_array+ib) ) * ( *(bbp+ib) + *(vcc->bbw_all+ib) ) /
						*(u_array+ib);
	}	
	// estimate adg at 443
	adg_443 = (
				( *(a) - ratio_aph * *(a+1) ) 
				- ( *(vcc->aw_all) - ratio_aph * *(vcc->aw_all+1) )
				) 
				/ (ratio_adg - ratio_aph);	
	//aph_443 =  *(a+1) - *(vcc->aw + 1) - adg_443;
	aph_443 =  *(a+1) - *(vcc->aw_all+1) - adg_443;

	// assign aph(443), adg(443), bbp(443) to the qaa matrix to use as ref.
	// for band shifting...
	vcc->aph_ref = aph_443;
	vcc->adg_ref = adg_443;
	vcc->bbp_ref = bbp[1];
}

float BandInterp(ccstr* vcc, float* rrs_sh, int idx_0,int num_in){
	/* Inverse distance interpolation
     * For example in the case of SeaWiFS data, rrs(510) and rrs(555) are used to calculate
	 * two separate rrs(531). These are used here to calculate a final rrs(531).
     * A limit of 2 input bands is assumed.
	 */ 
    
	float *wts,wts_sum=0;
    float rrs_sh_interp=0;
	int idx,ib;
    CheckNULL(wts = (float*)malloc(num_in*sizeof(float)),"wts",__LINE__,__FILE__);
        
    for (idx=0;idx<num_in;idx++){
        ib = idx_0 + idx;
        wts[idx] = 1.0 / abs(vcc->wvl_o[ib] - vcc->wvl_i[ib]);
        wts_sum += wts[idx];
    }
    
    for (idx = 0;idx<num_in;idx++){
		rrs_sh_interp += (wts[idx] * rrs_sh[idx]) / (wts_sum); 
       
    }
    free(wts);    
    return rrs_sh_interp;
}
float BandShift(ccstr* vcc, l2str* l2rec, int32_t ip,int target_band){
//float  bandShift(ccstr *vcc,l2str *l2rec,int32_t ip,int32_t ibvc){
// compare wave_vc[ibvc] to wvl_o. if not in bandshift output
// then assume it is already in the input, simply copy the record 
// to the output array 
	int ib,ipb,idx;
    int t_start,t_num=0; // target_band stuff
	float s = 0.015,g0 = 0.08945,g1 = 0.1247;
	float rrs_blue = 0,rrs_green = 0,b2g;
	float sdg,yy,chla,ll_i,ll_o;
	float aph_i,adg_i,bbp_i,a_tot_i,bb_tot_i;
	float aph_o,adg_o,bbp_o,a_tot_o,bb_tot_o;
	float qaa_fwd_i,qaa_fwd_o;
	float qaa_rrs_bw_o,qaa_rrs_aw_o,qaa_rrs_aw_i,qaa_rrs_bw_i;
	float correc_fac;
	int blue_idx,green_idx;
	float *rrs_in, *rrs_sh,fin_sh_rrs,estRrs;
    float blueBand=443;
	
    // locating target_band in the input record...
    for (ib=0;ib<vcc->nBands;ib++){
        if (vcc->wvl_o[ib] == target_band) {
            if (t_num==0)
                t_start = ib;
            t_num++;
        }
    }
      
    if (t_num == 0) {
        printf("\n -E- line %d Unable to find target band,%d => t_num = 0 :( \n", __LINE__,target_band);
        exit(FATAL_ERROR);
        }
    
    CheckNULL(rrs_in = (float*)malloc(t_num*sizeof(float)),"rrs_in",__LINE__,__FILE__);
	CheckNULL(rrs_sh = (float*)malloc(t_num*sizeof(float)),"rrs_sh",__LINE__,__FILE__);
    
    for (idx=0;idx<t_num;idx++)
        { // get Rrs input for band shift
        ib = t_start+idx;
		ipb = ip * l2rec->nbands + bindex_get(vcc->wvl_i[ib]);
		//rrs_in[idx] = l2rec->Rrs[ipb];
		rrs_in[idx] = ToBelowWater(l2rec->Rrs[ipb]);
        }
    
    blue_idx = bindex_get(blueBand);
    ipb = ip * l2rec->nbands + blue_idx;
    //rrs_blue = l2rec->Rrs[ipb];
    rrs_blue = ToBelowWater(l2rec->Rrs[ipb]);

    green_idx = bindex_get(vcc->greenBand);
    ipb = ip * l2rec->nbands + green_idx;
    //rrs_green = l2rec->Rrs[ipb];
    rrs_green = ToBelowWater(l2rec->Rrs[ipb]);
	b2g = rrs_blue / rrs_green;
	sdg = s + 0.002 / (0.6 + b2g);
	yy = 2 * (1 - 1.2 * exp(-0.9 * b2g));
	chla = pow( (vcc->aph_ref / vcc->a_sms),(1 / (1 - vcc->b_sms)));
	
	for (idx=0;idx<t_num;idx++){ // compute and apply band_shift factor for each input rrs
        ib  = t_start + idx;
		// compute IOPs for input wvls
		ll_i = vcc->wvl_i[ib] - vcc->wvl_sms;
		aph_i = pow(( vcc->a_i[ib] * chla),(1 - vcc->b_i[ib]));
		adg_i = vcc->adg_ref * exp(-sdg * ll_i);
		bbp_i = vcc->bbp_ref * pow((vcc->wvl_sms / vcc->wvl_i[ib]),yy);

		a_tot_i = aph_i + adg_i + vcc->aw_i[ib];
		bb_tot_i = bbp_i + vcc->bbw_i[ib];
		
		// compute IOPs for output wvls
		ll_o = vcc->wvl_o[ib] - vcc->wvl_sms;
		aph_o = pow(  (vcc->a_o[ib] * chla),(1 - vcc->b_o[ib])  );
		adg_o = vcc->adg_ref * exp(-sdg * ll_o);
		bbp_o = vcc->bbp_ref * pow( (vcc->wvl_sms / vcc->wvl_o[ib]) ,yy );
		
		
		a_tot_o = aph_o + adg_o + vcc->aw_o[ib];
		bb_tot_o = bbp_o + vcc->bbw_o[ib];
		
		// use qaa in fwd mode to get above-water rrs @ input & output wvls.
		qaa_fwd_i = bb_tot_i / (a_tot_i + bb_tot_i);
		qaa_rrs_bw_i = (g0 + g1 * qaa_fwd_i) * qaa_fwd_i;
		qaa_rrs_aw_i = ToAboveWater(qaa_rrs_bw_i);
	
		qaa_fwd_o = bb_tot_o / (a_tot_o + bb_tot_o);
		qaa_rrs_bw_o = (g0 + g1 * qaa_fwd_o) * qaa_fwd_o;
		qaa_rrs_aw_o = ToAboveWater(qaa_rrs_bw_o);
	
		correc_fac = qaa_rrs_aw_o / qaa_rrs_aw_i;
        
		rrs_sh[idx] = ToAboveWater(rrs_in[idx]) * correc_fac;
	}
    if(t_num>1) // Interpolation needed if >1 wavelengths used for bandshift.	
        {
    	fin_sh_rrs = BandInterp(vcc,rrs_sh,t_start,t_num);	
        }
    else
        fin_sh_rrs = rrs_sh[0];  
  
    free(rrs_sh);
    free(rrs_in);
    return fin_sh_rrs; 
}

void CheckNULL(void *ptr,char* ptrName,int lineNum,char* filename)
    {
    if (ptr==NULL){
        printf("\n --E-- pointer %s at line %d is %s not allocated :(\n",ptrName,lineNum,filename);
        exit(FATAL_ERROR);
        }
    }

float ToBelowWater(float rrsAbove)
    {
    float rrsBelow;
    rrsBelow = rrsAbove /(0.52 + 1.7 * rrsAbove);
    return rrsBelow;
    }

float ToAboveWater(float rrsBelow)
    {
    float rrsAbove;
    rrsAbove = 0.52 * rrsBelow / (1 - 1.7 * rrsBelow);
    return rrsAbove;
    }

