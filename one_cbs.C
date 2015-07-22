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


int one_cbs(double x1,  double sig, double a1, double n1,  
     double limlow, double limhi ){


  FILE* log;
  log=fopen("ro_cbsfit.log","a");
  if (log!=NULL){ 
    fprintf( log, "\n############draw one###################################\n%s", "" );   
  }

  TVirtualPad *oldpad=NULL;
 printf("oldpad == %lld\n", (Long64_t)oldpad );

  if (gPad==NULL){    TCanvas *c=new TCanvas();  }

  if (gPad==NULL){
    oldpad=NULL;
  printf("oldpad IN GPAD NULL == %lld\n", (Long64_t)oldpad );
 }else{
    oldpad=gPad;
  printf("oldpad OUTOF GPAD NULL == %lld\n", (Long64_t)oldpad );
  }

 printf("oldpad == %lld\n", (Long64_t)oldpad );




 // if (pad!=NULL){
 //   pad->cd();
 // }else{
 //  TCanvas *chi;
 //  chi=(TCanvas*)gROOT->GetListOfCanvases()->FindObject("c_roofit");
 //  if (chi!=NULL){
 //    chi->cd();
 //  }else{
 //    chi=new TCanvas("c_roofit","roofit check result window");
 //    chi->cd();
 //  }
 // }//pad !=NULL





  // in parameters  TH1F *h2;
  double pk[1];
  double si[1];
  //  double ar[1];
  double al[1];
  // shadows  double n1;

  double  xmin=0;
  double  xmax=limhi;
  double  smin=sig;
  double  smax=sig;
  //  double  amin=area;
  //  double  amax=area;
  
  double    almin=0.;
  double    almax=100.;
  double    nmin=1.0;
  double    nmax=100.;
  
  pk[0]=x1;






    si[0]=sig;
    //    ar[0]=area;
    al[0]=a1;
    //---    n1=5.1;



    printf("x=%5.3f L(%f..%f)  sigma=%5.3f   alpha=%5.3f   n1= %5.3f\n", 
	   pk[0], xmin,xmax,si[0],  al[0], n1 );
  //================ MAIN variable 1D ==============
  RooRealVar     x("x",    "x",   xmin, xmax);
  RooRealVar mean1("mean1", "mean_", pk[0], xmin,xmax);
  RooRealVar sigm1("sigm1", "sigma", si[0], smin, smax);
  //  RooRealVar area1("area1", "area_", ar[0], amin, amax);
  //=================== only one instance ===================
  RooRealVar alph1("alph1", "alph1",  al[0], almin, almax );
  RooRealVar ncbs1("n1", "n1",  n1, nmin, nmax );


    RooAbsPdf  *cb1   = new  RooCBShape("cb1","CBS", x, mean1, sigm1,alph1,ncbs1);
  //  RooCBShape cb1("cb1","CBS", x, mean1, sigm1,alph1,ncbs1);
  //  RooGaussian cb1("cb1","CBS", x, mean1, sigm1 );
  //  RooExtendPdf ecb1("ecb1","ECBS", cb1, area1 );


  RooPlot* xframe = x.frame(); //  printf("verbose frame:\n","");  xframe->Print("v"); 
  cb1->plotOn(xframe,  LineColor(kBlack) );

 gPad->SetLogy( oldpad->GetLogy() );
 gPad->SetLogy( 0 );
 gPad->Modified();gPad->Update();

RooAbsReal* integ = cb1->createIntegral (x ); 
double IntegralValue = integ->getVal(); 
//cout << IntegralValue << endl; 


//=================range 1 
 x.setRange("fitted", limlow, limhi ) ;  
 cb1->plotOn(xframe,  LineColor(kRed), Range("fitted")   );
 RooAbsReal* igx_sig = cb1->createIntegral(x, Range("fitted") ) ;
 double rangeint=igx_sig->getVal();

//=================range 2
 x.setRange("unseen", 0, limlow ) ;  
 cb1->plotOn(xframe,  LineColor(kGreen), Range("unseen")   );
 RooAbsReal* unseen = cb1->createIntegral(x, Range("unseen") ) ;
 double unseenint=unseen->getVal();



 printf("Total= %f  range= %f  unseen= %f (SUM= %f)  correction= %f\n", 
	IntegralValue, rangeint , unseenint, rangeint+unseenint,
	IntegralValue/rangeint   );
  xframe->Draw();
 // //NOT  fitresult->plotOn( xframe , DrawOption("F")  );
 // emodelV.plotOn( xframe );
 // // xframe->Draw();
 // emodelV.Print();
 // emodelV.plotOn(xframe, LineColor(kRed),   DrawOption("l0z") );
 // // emodelV.plotOn(xframe, Components(ecb1),LineColor(kGreen),LineStyle(kDashed) );
 // emodelV.plotOn(xframe, Components(ebkg),LineColor(kYellow),LineStyle(kSolid)  );
 // emodelV.plotOn(xframe, Components(RooArgSet(ecb1,ebkg)),LineColor(kBlue),LineStyle(kDashed) );
 // if (npeaks>1){
 // emodelV.plotOn(xframe, Components(RooArgSet(ecb2,ebkg)),LineColor(kBlue),LineStyle(kDashed) );
 // // emodelV.plotOn(xframe, Components(ecb2),LineColor(kGreen),LineStyle(kDashed) );
 // }


 if (oldpad!=NULL){ 
   //   printf("ODLPAD = %lld\n", (Long64_t) oldpad);
   if (oldpad!=gPad){
   oldpad->cd();
   }
 } 

 return 0;
}//==================================ro_cbs ======


