#define NMAX 15

#include "TPad.h"
#include "TROOT.h" // gpad


//#include "Riostream.h"
#include "Rtypes.h"
#include "TObject.h"




#ifndef __CINT__
#include "RooGlobalFunc.h"
//#include "TFitResult.h"
//#include "RooFitResult.h"
#endif

#include "RooRealVar.h"
#include "RooDataSet.h"

#include "RooGaussian.h"
#include "RooNovosibirsk.h"
#include "RooChebychev.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"

#include "RooExtendPdf.h"


#include "RooConstVar.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooDataHist.h"

//#include "RooRealVar.h"
//#include "RooDataSet.h"
//#include "RooGaussian.h"
//#include "RooConstVar.h"
//#include "TCanvas.h"
//#include "TAxis.h"
//#include "RooPlot.h"
#include "RooHist.h"

#include "TPRegexp.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TFile.h"


#include "RooCustomizer.h"
#include "RooFitResult.h"

// to get canvas list
#include "TROOT.h"

using namespace std; 
using namespace RooFit;   // THIS LOADS  ROOFIT  RooFit RooFit   RooFit   RooFit  


// version   2014 06 12
//using namespace std; 
//using namespace RooFit;   // THIS LOADS  ROOFIT  RooFit RooFit   RooFit   RooFit  


//  npeak  ,   double *x,  bg
//  n peaks,   positions, pn/p0/p1/p2  etc...
//  v budoucnu by to chtelo 2step - bg per partes; bgfix+ccc 
//===================================================================
RooFitResult* ro_cbs(TH1F *h2,int npeaks,double* xa,const char* bgstring,TPad *pad=NULL){
  int i;
  FILE* log;
  log=fopen("ro_cbsfit.log","a");
  if (log!=NULL){ 
    fprintf( log, "\n##########################################################\n%s", "" );   
    fprintf( log, "%s\n", h2->GetTitle() );
  }

  TVirtualPad *oldpad=NULL;
 printf("oldpad == %lld\n", (Long64_t)oldpad );

  if (gPad==NULL){
    oldpad=NULL;
  printf("oldpad IN GPAD NULL == %lld\n", (Long64_t)oldpad );
 }else{
    oldpad=gPad;
  printf("oldpad OUTOF GPAD NULL == %lld\n", (Long64_t)oldpad );
  }

 printf("oldpad == %lld\n", (Long64_t)oldpad );


 if (pad!=NULL){
   pad->cd();
 }else{
  TCanvas *chi;
  chi=(TCanvas*)gROOT->GetListOfCanvases()->FindObject("c_roofit");
  if (chi!=NULL){
    chi->cd();
  }else{
    chi=new TCanvas("c_roofit","roofit check result window");
    chi->cd();
  }
 }//pad !=NULL


  if (npeaks>NMAX){ printf("TOO MANY PEAKS%s\n",""); return NULL;}

  // in parameters  TH1F *h2;
  double pk[NMAX];
  double si[NMAX];
  double ar[NMAX];
  double al[NMAX];
  double n1;

  //  double min=h2->GetXaxis()->GetFirst(),
  //    max=h2->GetXaxis()->GetLast(),
  /* 1/  min, max -  is x value
   * 2/  however - bin can be whatever
   *
   */
  // there was some struggle here.
  //from verticals- use getxaxis-getxmin
  //from cint - use getfirst??
  //============================== GetXmin is real minimum
  //GetFirst gets 1st bin - 0 if not zoom
  double 
    //    min=h2->GetXaxis()->GetXmin(),
    //    max=h2->GetXaxis()->GetXmax(),

    min=h2->GetBinCenter(h2->GetXaxis()->GetFirst()),
    max=h2->GetBinCenter(h2->GetXaxis()->GetLast()),
    smin=(max-min)/1000,
    smax=(max-min),
    amin=h2->Integral( h2->GetXaxis()->FindBin(min), 
		       h2->GetXaxis()->FindBin(max) ,"width")/1000.,
    
    amax=h2->Integral( h2->GetXaxis()->FindBin(min), 
		       h2->GetXaxis()->FindBin(max) ,"width"),
    almin=0.,
    almax=100.,
    nmin=1.0,
    nmax=100.;


  if (log!=NULL){ 
    fprintf( log, "min=%f\n", min );
    fprintf( log, "max=%f\n", max );

    fprintf( log, "smin=%f\n", smin );
    fprintf( log, "smax=%f\n", smax );

    fprintf( log, "amin=%f\n", amin );
    fprintf( log, "amax=%f\n", amax );
  }

  printf("FROM ====== %f DOWNTO  %f \n" ,max, min );
  if (xa!=NULL){
    for (i=0;i<npeaks;i++){
      pk[i]=xa[i];
    }
  }else{
    for (i=0;i<npeaks;i++){
      pk[i]=min+ (i+1)*(max-min)/(npeaks+1);
    }
  }//--------------------



  for (i=0;i<npeaks;i++){
    //    pk[i]=1.0*(i+1)/(1+npeaks)*(max-min) +min;
    //    si[i]=4.3; 
    //    pk[i]=xa[i];
    si[i]=smax/4  /2; //    /2 EMPIRICAL 1peak
    ar[i]=1.0*(amax+amin)/npeaks/2.  /3.; //     /3 EMPIRICAL 1peak AREA  
    al[i]=1.3;
    if (log!=NULL){ 
      fprintf( log, "pk%02d %9.3f   %9.3f  %9.3f \n", 
	       i+1, 
	       pk[i], ar[i], si[i] );
    }
  }//---------for all npeaks --------
  n1=5.1;

  //  smin=0.5* si[0];
  //  smax=2.*si[0];
 


 printf("oldpad == %lld\n", (Long64_t)oldpad );

  //================ MAIN variable 1D ==============
  RooRealVar       x("x",    "x",   min, max);
  /*
   *
   *  CHEBYSHEV  HERE
   */
  char sbg[10];
  sprintf(sbg,"%s", bgstring );

  RooRealVar   bgarea("bgarea", "bgarea",  1.0*(amin+amax)/2., 0, amax);  

 // Build Chebychev polynomial p.d.f.  
 // RooRealVar a0("a0","a0", 0.) ;
  RooRealVar a0("a0","a0",    0., -10, 10) ;
  RooRealVar a1("a1","a1",    0., -10, 10) ;
  RooRealVar a2("a2","a2",    0., -10, 10) ;
  RooRealVar a3("a3","a3",    0., -10, 10) ;
  RooArgSet setcheb;
  if ( strstr(sbg,"pn")!=NULL ){ setcheb.add(a0);a0.setVal(0.);a0.setConstant(kTRUE);bgarea.setVal(0.);bgarea.setConstant(kTRUE);}
  if ( strstr(sbg,"p0")!=NULL ){ setcheb.add(a0);  a0=0.; a0.setConstant(kTRUE); }
  if ( strstr(sbg,"p1")!=NULL ){ setcheb.add(a0); }
  if ( strstr(sbg,"p2")!=NULL ){ setcheb.add(a1);setcheb.add(a0); }
  if ( strstr(sbg,"p3")!=NULL ){ setcheb.add(a2);setcheb.add(a1); setcheb.add(a0);}
  if ( strstr(sbg,"p4")!=NULL ){ setcheb.add(a3);setcheb.add(a2); setcheb.add(a1); setcheb.add(a0); }
  RooChebychev bkg("bkg","Background",x, setcheb ) ;


  RooRealVar mean1("mean1", "mean_", pk[0], min,max);
  RooRealVar mean2("mean2", "mean_", pk[1], min,max);
  RooRealVar mean3("mean3", "mean_", pk[2], min,max);
  RooRealVar mean4("mean4", "mean_", pk[3], min,max);
  RooRealVar mean5("mean5", "mean_", pk[4], min,max);
  RooRealVar mean6("mean6", "mean_", pk[5], min,max);
  RooRealVar mean7("mean7", "mean_", pk[6], min,max);
  RooRealVar mean8("mean8", "mean_", pk[7], min,max);
  RooRealVar mean9("mean9", "mean_", pk[8], min,max);
  RooRealVar meana("meana", "mean_", pk[9], min,max);
  RooRealVar meanb("meanb", "mean_", pk[10], min,max);
  RooRealVar meanc("meanc", "mean_", pk[11], min,max);
  RooRealVar meand("meand", "mean_", pk[12], min,max);
  RooRealVar meane("meane", "mean_", pk[13], min,max);
  RooRealVar meanf("meanf", "mean_", pk[14], min,max);

  RooRealVar sigm1("sigm1", "sigma", si[0], smin, smax);
  RooRealVar sigm2("sigm2", "sigma", si[1], smin, smax);
  RooRealVar sigm3("sigm3", "sigma", si[2], smin, smax);
  RooRealVar sigm4("sigm4", "sigma", si[3], smin, smax);
  RooRealVar sigm5("sigm5", "sigma", si[4], smin, smax);
  RooRealVar sigm6("sigm6", "sigma", si[5], smin, smax);
  RooRealVar sigm7("sigm7", "sigma", si[6], smin, smax);
  RooRealVar sigm8("sigm8", "sigma", si[7], smin, smax);
  RooRealVar sigm9("sigm9", "sigma", si[8], smin, smax);
  RooRealVar sigma("sigma", "sigma", si[9], smin, smax);
  RooRealVar sigmb("sigmb", "sigma", si[10], smin, smax);
  RooRealVar sigmc("sigmc", "sigma", si[11], smin, smax);
  RooRealVar sigmd("sigmd", "sigma", si[12], smin, smax);
  RooRealVar sigme("sigme", "sigma", si[13], smin, smax);
  RooRealVar sigmf("sigmf", "sigma", si[14], smin, smax);

  RooRealVar    fs2("fs2","fs2", 0.0 ); // factor sigma = 1==free or 0==fix
  RooRealVar    fs3("fs3","fs3", 0.0 ); // factor sigma = 1==free or 0==fix
  RooRealVar    fs4("fs4","fs4", 0.0 ); // factor sigma = 1==free or 0==fix
  RooRealVar    fs5("fs5","fs5", 0.0 ); // factor sigma = 1==free or 0==fix
  RooRealVar    fs6("fs6","fs6", 0.0 ); // factor sigma = 1==free or 0==fix
  RooRealVar    fs7("fs7","fs7", 0.0 ); // factor sigma = 1==free or 0==fix
  RooRealVar    fs8("fs8","fs8", 0.0 ); // factor sigma = 1==free or 0==fix
  RooRealVar    fs9("fs9","fs9", 0.0 ); // factor sigma = 1==free or 0==fix
  RooRealVar    fsa("fsa","fsa", 0.0 ); // factor sigma = 1==free or 0==fix
  RooRealVar    fsb("fsb","fsb", 0.0 ); // factor sigma = 1==free or 0==fix
  RooRealVar    fsc("fsc","fsc", 0.0 ); // factor sigma = 1==free or 0==fix
  RooRealVar    fsd("fsd","fsd", 0.0 ); // factor sigma = 1==free or 0==fix
  RooRealVar    fse("fse","fse", 0.0 ); // factor sigma = 1==free or 0==fix
  RooRealVar    fsf("fsf","fsf", 0.0 ); // factor sigma = 1==free or 0==fix

  RooFormulaVar fsigm2( "fsigm2", " fs2*sigm2 + (1.0-fs2)*sigm1" , RooArgSet(sigm1,sigm2,fs2) );
  RooFormulaVar fsigm3( "fsigm3", " fs3*sigm3 + (1.0-fs3)*sigm1" , RooArgSet(sigm1,sigm3,fs3) );
  RooFormulaVar fsigm4( "fsigm4", " fs4*sigm4 + (1.0-fs4)*sigm1" , RooArgSet(sigm1,sigm4,fs4) );
  RooFormulaVar fsigm5( "fsigm5", " fs5*sigm5 + (1.0-fs5)*sigm1" , RooArgSet(sigm1,sigm5,fs5) );
  RooFormulaVar fsigm6( "fsigm6", " fs6*sigm6 + (1.0-fs6)*sigm1" , RooArgSet(sigm1,sigm6,fs6) );
  RooFormulaVar fsigm7( "fsigm7", " fs7*sigm7 + (1.0-fs7)*sigm1" , RooArgSet(sigm1,sigm7,fs7) );
  RooFormulaVar fsigm8( "fsigm8", " fs8*sigm8 + (1.0-fs8)*sigm1" , RooArgSet(sigm1,sigm8,fs8) );
  RooFormulaVar fsigm9( "fsigm9", " fs9*sigm9 + (1.0-fs9)*sigm1" , RooArgSet(sigm1,sigm9,fs9) );
  RooFormulaVar fsigma( "fsigma", " fsa*sigma + (1.0-fsa)*sigm1" , RooArgSet(sigm1,sigma,fsa) );
  RooFormulaVar fsigmb( "fsigmb", " fsb*sigmb + (1.0-fsb)*sigm1" , RooArgSet(sigm1,sigmb,fsb) );
  RooFormulaVar fsigmc( "fsigmc", " fsc*sigmc + (1.0-fsc)*sigm1" , RooArgSet(sigm1,sigmc,fsc) );
  RooFormulaVar fsigmd( "fsigmd", " fsd*sigmd + (1.0-fsd)*sigm1" , RooArgSet(sigm1,sigmd,fsd) );
  RooFormulaVar fsigme( "fsigme", " fse*sigme + (1.0-fse)*sigm1" , RooArgSet(sigm1,sigme,fse) );
  RooFormulaVar fsigmf( "fsigmf", " fsf*sigmf + (1.0-fsf)*sigm1" , RooArgSet(sigm1,sigmf,fsf) );


  // -------- fsx   0== FIX to sigm1;    1==free sigma
  fs2.setVal( 0.0 );  // 
  fs3.setVal( 1.0 );  
  fs4.setVal( 1.0 );  
  fs5.setVal( 1.0 );  
  fs6.setVal( 1.0 );  
  fs7.setVal( 1.0 );  
  fs8.setVal( 1.0 );  
  fs9.setVal( 1.0 );  
  fsa.setVal( 1.0 );  
  fsb.setVal( 1.0 );  
  fsc.setVal( 1.0 );  
  fsd.setVal( 1.0 );  
  fse.setVal( 1.0 );  
  fsf.setVal( 1.0 );  
  
  fs2.setConstant( kTRUE );  
  fs3.setConstant( kTRUE );  
  fs4.setConstant( kTRUE );  
  fs5.setConstant( kTRUE );  
  fs6.setConstant( kTRUE );  
  fs7.setConstant( kTRUE );  
  fs8.setConstant( kTRUE );  
  fs9.setConstant( kTRUE );  
  fsa.setConstant( kTRUE );  
  fsb.setConstant( kTRUE );  
  fsc.setConstant( kTRUE );  
  fsd.setConstant( kTRUE );  
  fse.setConstant( kTRUE );  
  fsf.setConstant( kTRUE );  

  if (fs2.getVal()==0.0){sigm2.setConstant( kTRUE );}
  if (fs3.getVal()==0.0){sigm3.setConstant( kTRUE );}
  if (fs4.getVal()==0.0){sigm4.setConstant( kTRUE );}
  if (fs5.getVal()==0.0){sigm5.setConstant( kTRUE );}
  if (fs6.getVal()==0.0){sigm6.setConstant( kTRUE );}
  if (fs7.getVal()==0.0){sigm7.setConstant( kTRUE );}
  if (fs8.getVal()==0.0){sigm8.setConstant( kTRUE );}
  if (fs9.getVal()==0.0){sigm9.setConstant( kTRUE );}
  if (fsa.getVal()==0.0){sigma.setConstant( kTRUE );}
  if (fsb.getVal()==0.0){sigmb.setConstant( kTRUE );}
  if (fsc.getVal()==0.0){sigmc.setConstant( kTRUE );}
  if (fsd.getVal()==0.0){sigmd.setConstant( kTRUE );}
  if (fse.getVal()==0.0){sigme.setConstant( kTRUE );}
  if (fsf.getVal()==0.0){sigmf.setConstant( kTRUE );}


  // extra thing that goes to RootSet
  RooRealVar area1("area1", "area_", ar[0], amin, amax);
  RooRealVar area2("area2", "area_", ar[1], amin, amax);
  RooRealVar area3("area3", "area_", ar[2], amin, amax);
  RooRealVar area4("area4", "area_", ar[3], amin, amax);
  RooRealVar area5("area5", "area_", ar[4], amin, amax);
  RooRealVar area6("area6", "area_", ar[5], amin, amax);
  RooRealVar area7("area7", "area_", ar[6], amin, amax);
  RooRealVar area8("area8", "area_", ar[7], amin, amax);
  RooRealVar area9("area9", "area_", ar[8], amin, amax);
  RooRealVar areaa("areaa", "area_", ar[9], amin, amax);
  RooRealVar areab("areab", "area_", ar[10], amin, amax);
  RooRealVar areac("areac", "area_", ar[11], amin, amax);
  RooRealVar aread("aread", "area_", ar[12], amin, amax);
  RooRealVar areae("areae", "area_", ar[13], amin, amax);
  RooRealVar areaf("areaf", "area_", ar[14], amin, amax);



  //=================== only one instance ===================
  RooRealVar alph1("alph1", "alph1",  al[0], almin, almax );
  RooRealVar ncbs1("n1", "n1",  n1, nmin, nmax );

  RooCBShape cb1("cb1","CBS", x,mean1, sigm1,alph1,ncbs1);
  RooCBShape cb2("cb2","CBS", x,mean2,fsigm2,alph1,ncbs1);
  RooCBShape cb3("cb3","CBS", x,mean3,fsigm3,alph1,ncbs1);
  RooCBShape cb4("cb4","CBS", x,mean4,fsigm4,alph1,ncbs1);
  RooCBShape cb5("cb5","CBS", x,mean5,fsigm5,alph1,ncbs1);
  RooCBShape cb6("cb6","CBS", x,mean6,fsigm6,alph1,ncbs1);
  RooCBShape cb7("cb7","CBS", x,mean7,fsigm7,alph1,ncbs1);
  RooCBShape cb8("cb8","CBS", x,mean8,fsigm8,alph1,ncbs1);
  RooCBShape cb9("cb9","CBS", x,mean9,fsigm9,alph1,ncbs1);
  RooCBShape cba("cba","CBS", x,meana,fsigma,alph1,ncbs1);
  RooCBShape cbb("cbb","CBS", x,meanb,fsigmb,alph1,ncbs1);
  RooCBShape cbc("cbc","CBS", x,meanc,fsigmc,alph1,ncbs1);
  RooCBShape cbd("cbd","CBS", x,meand,fsigmd,alph1,ncbs1);
  RooCBShape cbe("cbe","CBS", x,meane,fsigme,alph1,ncbs1);
  RooCBShape cbf("cbf","CBS", x,meanf,fsigmf,alph1,ncbs1);

  RooExtendPdf ecb1("ecb1","ECBS", cb1, area1 );
  RooExtendPdf ecb2("ecb2","ECBS", cb2, area2 );
  RooExtendPdf ecb3("ecb3","ECBS", cb3, area3 );
  RooExtendPdf ecb4("ecb4","ECBS", cb4, area4 );
  RooExtendPdf ecb5("ecb5","ECBS", cb5, area5 );
  RooExtendPdf ecb6("ecb6","ECBS", cb6, area6 );
  RooExtendPdf ecb7("ecb7","ECBS", cb7, area7 );
  RooExtendPdf ecb8("ecb8","ECBS", cb8, area8 );
  RooExtendPdf ecb9("ecb9","ECBS", cb9, area9 );
  RooExtendPdf ecba("ecba","ECBS", cba, areaa );
  RooExtendPdf ecbb("ecbb","ECBS", cbb, areab );
  RooExtendPdf ecbc("ecbc","ECBS", cbc, areac );
  RooExtendPdf ecbd("ecbd","ECBS", cbd, aread );
  RooExtendPdf ecbe("ecbe","ECBS", cbe, areae );
  RooExtendPdf ecbf("ecbf","ECBS", cbf, areaf );

 RooArgList erl;  
 if (npeaks>0)erl.add( ecb1 );  
 if (npeaks>1)erl.add( ecb2 );  
 if (npeaks>2)erl.add( ecb3 );  
 if (npeaks>3)erl.add( ecb4 );  
 if (npeaks>4)erl.add( ecb5 );  
 if (npeaks>5)erl.add( ecb6 );  
 if (npeaks>6)erl.add( ecb7 );  
 if (npeaks>7)erl.add( ecb8 );  
 if (npeaks>8)erl.add( ecb9 );  
 if (npeaks>9)erl.add( ecba );  
 if (npeaks>10)erl.add( ecbb );  
 if (npeaks>11)erl.add( ecbc );  
 if (npeaks>12)erl.add( ecbd );  
 if (npeaks>13)erl.add( ecbe );  
 if (npeaks>14)erl.add( ecbf );  

  RooExtendPdf ebkg("ebkg","EBkg", bkg, bgarea );

  if ( strstr(sbg,"pn")==NULL){ erl.add( ebkg ); }



  //nechci customizer, jenom model
 RooAddPdf emodelV("emodelV","emultiplet", erl  );
 emodelV.Print("v");



 RooDataHist datah("datah","datah with x", x, h2);
 RooPlot* xframe = x.frame();
 datah.plotOn(xframe );
 // xframe.Print("v"); 
 // TFitResult* fitresult;
 RooFitResult* fitresult;

 fitresult = emodelV.fitTo( datah , Save()  ); 
 fitresult->Print("v") ;


 char varm[10];
 char vara[10];
 char vars[10];


 double sigmalog,dsigmalog;
 if (log!=NULL){ 
    fprintf( log, "----------------------------------------------------------\n%s", "" );
    //    fprintf( log, "pk%02d %9.3f   %9.3f  %9.3f\n", i,
    //	     mean1.getVal(),area1.getVal(),sigm1.getVal() );
    for (i=0;i<npeaks;i++){
      sprintf(varm, "mean%x", i+1 );
      sprintf(vara, "area%x", i+1 );
      sprintf(vars, "sigm%x", i+1 );
      RooRealVar* t1=(RooRealVar*)fitresult->floatParsFinal().find( varm ) ;
      RooRealVar* t2=(RooRealVar*)fitresult->floatParsFinal().find( vara ) ;
      RooRealVar* t3=(RooRealVar*)fitresult->floatParsFinal().find( vars ) ;
      if (t3!=NULL){ sigmalog=t3->getVal();dsigmalog=t3->getError();}
    fprintf(log,"pk%02d %9.3f   %9.3f  %9.3f /%9.3f   %9.3f  %9.3f \n", i+1, 
	      t1->getVal() ,  t2->getVal() ,  sigmalog , 
	      t1->getError(), t2->getError(), dsigmalog );
            
    }// for all npeaks
 }// if log..............................................................................


 //NOT  fitresult->plotOn( xframe , DrawOption("F")  );
 emodelV.plotOn( xframe );
 // xframe->Draw();

 emodelV.Print();

 emodelV.plotOn(xframe, LineColor(kRed),   DrawOption("l0z") );

 emodelV.plotOn(xframe, Components(ecb1),LineColor(kGreen),LineStyle(kDashed) );

 emodelV.plotOn(xframe, Components(ebkg),LineColor(kYellow),LineStyle(kSolid)  );


 emodelV.plotOn(xframe, Components(RooArgSet(ecb1,ebkg)),LineColor(kBlue),LineStyle(kDashed) );
 if (npeaks>1){
 emodelV.plotOn(xframe, Components(RooArgSet(ecb2,ebkg)),LineColor(kBlue),LineStyle(kDashed) );
 emodelV.plotOn(xframe, Components(ecb2),LineColor(kGreen),LineStyle(kDashed) );
 }
 if (npeaks>2){
 emodelV.plotOn(xframe, Components(RooArgSet(ecb3,ebkg)),LineColor(kBlue),LineStyle(kDashed) );
 emodelV.plotOn(xframe, Components(ecb3),LineColor(kGreen),LineStyle(kDashed) );
 }
 if (npeaks>3){
 emodelV.plotOn(xframe, Components(RooArgSet(ecb4,ebkg)),LineColor(kBlue),LineStyle(kDashed) );
 emodelV.plotOn(xframe, Components(ecb4),LineColor(kGreen),LineStyle(kDashed) );
 }
 if (npeaks>4){
 emodelV.plotOn(xframe, Components(RooArgSet(ecb5,ebkg)),LineColor(kBlue),LineStyle(kDashed) );
 emodelV.plotOn(xframe, Components(ecb5),LineColor(kGreen),LineStyle(kDashed) );
 }
 if (npeaks>5){
 emodelV.plotOn(xframe, Components(RooArgSet(ecb6,ebkg)),LineColor(kBlue),LineStyle(kDashed) );
 emodelV.plotOn(xframe, Components(ecb6),LineColor(kGreen),LineStyle(kDashed) );
 }
 if (npeaks>6){
 emodelV.plotOn(xframe, Components(RooArgSet(ecb7,ebkg)),LineColor(kBlue),LineStyle(kDashed) );
 emodelV.plotOn(xframe, Components(ecb7),LineColor(kGreen),LineStyle(kDashed) );
 }
 if (npeaks>7){
 emodelV.plotOn(xframe, Components(RooArgSet(ecb8,ebkg)),LineColor(kBlue),LineStyle(kDashed) );
 emodelV.plotOn(xframe, Components(ecb8),LineColor(kGreen),LineStyle(kDashed) );
 }
 if (npeaks>8){
 emodelV.plotOn(xframe, Components(RooArgSet(ecb9,ebkg)),LineColor(kBlue),LineStyle(kDashed) );
 emodelV.plotOn(xframe, Components(ecb9),LineColor(kGreen),LineStyle(kDashed) );
 }
 if (npeaks>9){
 emodelV.plotOn(xframe, Components(RooArgSet(ecba,ebkg)),LineColor(kBlue),LineStyle(kDashed) );
 emodelV.plotOn(xframe, Components(ecba),LineColor(kGreen),LineStyle(kDashed) );
 }
 // ...

 xframe->Draw();


 gPad->Modified();gPad->Update();

 fitresult->floatParsFinal().Print("s") ;//crashes when fixed sigma2??

 printf("sigma  = =  %f chan\n",   sigm1.getVal()    );

 //last things
 if (log!=NULL){ fclose( log ); } 

 // printf("oldpad == %lld\n", (Long64_t)oldpad );

 if (oldpad!=NULL){ 
   //   printf("ODLPAD = %lld\n", (Long64_t) oldpad);
   if (oldpad!=gPad){
   oldpad->cd();
   }
 } 

 // return NULL ;
  return fitresult ;
// return sigm1.getVal() ;

}//==================================ro_cbs ======


