void plotDataSet(TString inputFile="mass_life.data") {

  RooRealVar t("t","t",-1.5,5.);
  RooRealVar mass("mass","mass",300.,700.);
  RooDataSet* data = RooDataSet::read(inputFile,RooArgSet(mass,t));

  RooRealVar tSig("tSig","tSig",-1.5,5.);
  RooRealVar massSig("massSig","massSig",300.,700.);
  RooDataSet* dataSig = RooDataSet::read("masslife_signal.data",RooArgSet(massSig,tSig));

  // Plot
  TCanvas * c = new TCanvas();
  c->Divide(2,2);
  c->cd(1);

  RooPlot* tframe = t.frame() ;
  data->plotOn(tframe) ;
  tframe->Draw() ;

  // Plot 
  c->cd(2);
  RooPlot* mframe = mass.frame() ;
  data->plotOn(mframe) ;
  mframe->Draw() ;
  // Plot 
  c->cd(4);
  RooPlot* mSigframe = massSig.frame() ;
  dataSig->plotOn(mSigframe) ;
  mSigframe->Draw() ;
  // Plot 
  c->cd(3);
  RooPlot* tSigframe = tSig.frame() ;
  dataSig->plotOn(tSigframe) ;
  tSigframe->Draw() ;
  c->SaveAs("data.pdf");
}

void signalLifetimeFit(TString inputFile="masslife_signal.data"){
  RooRealVar t("t","t",-2.,5.);
  RooRealVar mass("mass","mass",300.,700.);
  RooDataSet* data = RooDataSet::read(inputFile,RooArgSet(mass,t));
  t.setBins(10000,"cache") ; 
	mass.setBins(10000,"cache") ; 
  RooRealVar sigTau("sigTau","Signal lifetime",1.0,0.,5.);
  RooRealVar sigResolution("sigResolution","Resolution for signal lifetime",0.5,0.,5.0);
  RooRealVar sigResolutionCenter("sigResolutionCenter","Resolution lifetime offset",0.0,-1.,1.);

  RooGaussian resolutionGauss("resolutionGauss","Lifetime Resolution PDF",t,sigResolutionCenter,sigResolution);
  RooGenericPdf signalDecay("signalDecay","Lifetime pdf","0.5*(1+sign(t))/sigTau*exp(-t/sigTau)",RooArgSet(t,sigTau));

  RooFFTConvPdf fittingFunc("Convolution","lifetime (X) gauss",t,signalDecay,resolutionGauss);

  fittingFunc.fitTo(*data,RooFit::NumCPU(8));
  RooPlot* xframe = t.frame();
  data->plotOn(xframe);
  fittingFunc.plotOn(xframe);
  xframe->Draw();
}

void signalMassFit(TString inputFile="masslife_signal.data"){
	RooRealVar t("t","t",-2.,5.);
  RooRealVar mass("mass","mass",450.,480.);
  t.setBins(10000,"cache") ; 
	mass.setBins(10000,"cache") ; 
  RooDataSet* data = RooDataSet::read(inputFile,RooArgSet(mass,t));

  RooRealVar sigMass("sigMass","Signal mass peak position",460.,450.,480.) ;
  RooRealVar sigMassWidth("sigMassSigma","Signal mass peak width",10,0,50);


  RooGaussian fittingFunc("massGauss","Mass gaussian PDF",mass,sigMass,sigMassWidth);
  fittingFunc.fitTo(*data,RooFit::NumCPU(4));
  RooPlot* xframe = mass.frame();
  data->plotOn(xframe);
  fittingFunc.plotOn(xframe);
  xframe->Draw();
}

void lifetimeFit(TString inputFile="mass_life.data"){
 
  RooRealVar t("t","t",-2.,15.);
  RooRealVar mass("mass","mass",300.,700.);
  RooDataSet* data = RooDataSet::read(inputFile,RooArgSet(mass,t));
  t.setBins(10000,"cache") ; 
	mass.setBins(10000,"cache") ; 
  //signal
  RooRealVar sigTau("sigTau","Signal lifetime",0.3,0.,1.);
  RooRealVar sigResolution("sigResolution","Resolution lifetime for signal",0.5,0.,1.0);
  RooRealVar sigResolutionCenter("sigResolutionCenter","Resolution lifetime offset",0.0,0.,1.);

  RooGenericPdf signalDecay("signalDecay","Lifetime pdf","0.5*(1+sign(t))*0.5/sigTau*exp(-t/sigTau)",RooArgSet(t,sigTau));
  RooGaussian resolutionGauss("resolutionGauss","Lifetime Resolution PDF",t,sigResolutionCenter,sigResolution);
  RooFFTConvPdf signalPDF("signalConv","lifetime (X) gauss",t,signalDecay,resolutionGauss);

  //bkg1  gaussian bkg, no lifetime
  RooRealVar bkg1Resolution("bkg1Res","No lifetime resolution",0.1,0.,1.0);
  RooRealVar bkg1Offset("bkg1Offset","No lifetime offset",0.0);
  RooGaussian bkg1("bkg1Gauss","bkg1 PDF",t,bkg1Offset,bkg1Resolution);

  //bkg2 exp * gauss conv
  RooRealVar bkg2Tau("bkg2Tau","Bkg2 lifetime",2.0,1.,10.);
  RooRealVar bkg2Resolution("bkg2Res","Long lifetime resolution",0.1,0.,1.);
  RooRealVar bkg2Offset("bkg1Offset","Long lifetime offset",0.0);
  RooGenericPdf bkg2Decay("bkg2Decay","Long Lifetime pdf","0.5*(1+sign(t))/bkg2Tau*exp(-t/bkg2Tau)",RooArgSet(t,bkg2Tau));
  RooGaussian bkg2ResolutionGauss("bkg2ResolutionGauss","Long lifetime Resolution PDF",t,bkg2Offset,bkg2Resolution);
  RooFFTConvPdf bkg2("bkg2","lifetime (X) gauss bkg",t,bkg2Decay,bkg2ResolutionGauss);
  
  RooRealVar gfrac1("gfrac1","fraction of first bkg in tot bkg",0.5,0.,1.) ;
  RooRealVar gfrac2("gfrac2","fraction of background in data",0.5,0.,1.) ;
  RooAddPdf  totalBkg("totalBkg","Bkg1 + Bkg2",RooArgList(bkg2,bkg1),RooArgList(gfrac1));
  RooAddPdf  fitFunc("fitFunc","Sig + Bkg",RooArgList(totalBkg,signalPDF),RooArgList(gfrac2));

  RooFitResult* fitRes = fitFunc.fitTo(*data,RooFit::Save(),RooFit::NumCPU(20));
  TCanvas * c = new TCanvas();
  c->Divide(2,1);
  c->cd(1);
	gPad-> SetLogy();
  RooPlot* xframe = t.frame();
  data->plotOn(xframe);
  fitFunc.plotOn(xframe);
  xframe->Draw();

  c->cd(2);
  RooAbsReal* nll = fitFunc.createNLL(*data,RooFit::NumCPU(20));
  double minNLL =  fitRes->minNll();
  cout <<"Min Nll "<< minNLL<<"\n";
  RooPlot* Lframe = sigTau.frame(0.6,0.8);
  nll->plotOn(Lframe);  
  Lframe->SetMinimum(minNLL-1.);
  Lframe->SetMaximum(minNLL+5);
  Lframe->Draw();
  c->SaveAs("LifetimeFit.pdf");

}


void massFit(TString inputFile="mass_life.data"){
 
  RooRealVar t("t","t",-1.5,5.);
  RooRealVar mass("mass","mass",400.,600.);
  RooDataSet* data = RooDataSet::read(inputFile,RooArgSet(mass,t));
  t.setBins(10000,"cache") ; 
	mass.setBins(10000,"cache") ; 

  //signal
  RooRealVar sigMass("sigMass","Signal mass peak position",450,400,480) ;
  RooRealVar sigMassWidth("sigMassSigma","Signal mass peak width",15.,10.,30.);
  //background
  RooRealVar x1("x1","first null point",350,340,355) ;
  RooRealVar x2("x2","second null point",655,640,660) ;
  RooRealVar bkgFrac("bkgFrac","Fraction of background",0.8,0.5,1.) ;

  RooGaussian signalFunc("Gaussian","Mass PDF",mass,sigMass,sigMassWidth);
  RooGenericPdf backgroundFunc("Quadratic bkg","Quadratic PDF","-0.1*(x1-mass)*(x2-mass)*(1+sign(mass-x1))*(1+sign(x2-mass))",RooArgSet(x1,x2,mass));

  RooAddPdf  totalPDF("sum","Sig+bkg",RooArgList(backgroundFunc,signalFunc),RooArgList(bkgFrac));

  RooFitResult* fitRes = totalPDF.fitTo(*data,RooFit::Save(),RooFit::NumCPU(20));
  TCanvas * c = new TCanvas();
  c->Divide(2,1);
  c->cd(1);
  RooPlot* xframe = mass.frame();
  data->plotOn(xframe);
  totalPDF.plotOn(xframe,"L") ;
  xframe->Draw();

  c->cd(2);
  RooAbsReal* nll = totalPDF.createNLL(*data);
  double minNLL =  fitRes->minNll();
  cout <<"Min Nll "<< minNLL<<"\n";
  RooPlot* Lframe = sigMass.frame(463,466);
  nll->plotOn(Lframe);  
  Lframe->SetMinimum(minNLL-1.);
  Lframe->SetMaximum(minNLL+5);
  Lframe->Draw();
  c->SaveAs("massFit.pdf");
}

void profiles(TString inputFile="mass_life.data"){
  RooRealVar t("t","t",-2.,15.);
  RooRealVar mass("mass","mass",400.,600.);
  RooDataSet* data = RooDataSet::read(inputFile,RooArgSet(mass,t));
  t.setBins(10000,"cache"); 
	mass.setBins(10000,"cache");

		//signal
  RooRealVar sigMass("sigMass","Signal mass peak position",450,400,480) ;
  RooRealVar sigMassWidth("sigMassSigma","Signal mass peak width",10.4548);
  RooGaussian massSignalPDF("Gaussian","Mass PDF",mass,sigMass,sigMassWidth);

  RooRealVar sigTau("sigTau","Signal lifetime",0.3,0.,1.);
  RooRealVar sigResolution("sigResolution","Resolution lifetime for signal",0.422078);
  RooRealVar sigResolutionCenter("sigResolutionCenter","Resolution lifetime offset",0.0);
  RooGenericPdf signalDecay("signalDecay","Lifetime pdf","0.5*(1+sign(t))*0.5/sigTau*exp(-t/sigTau)",RooArgSet(t,sigTau));
  RooGaussian resolutionGauss("resolutionGauss","Lifetime Resolution PDF",t,sigResolutionCenter,sigResolution);
  RooFFTConvPdf lifetimeSignalPDF("signalConv","lifetime (X) gauss",t,signalDecay,resolutionGauss);

  RooProdPdf signalPDF("signalPDF", "Product of mass and lifetime signal PDFs", massSignalPDF, lifetimeSignalPDF);


  //background 1
  RooRealVar x1("x1","first null point",351.228) ;
  RooRealVar x2("x2","second null point",647.485) ;

  RooGenericPdf backgroundFunc("Quadratic bkg","Quadratic PDF","-0.1*(x1-mass)*(x2-mass)*(1+sign(mass-x1))*(1+sign(x2-mass))",RooArgSet(x1,x2,mass));

  RooRealVar bkg1Resolution("bkg1Res","No lifetime resolution",0.319940);
  RooRealVar bkg1Offset("bkg1Offset","No lifetime offset",0.0);
  RooGaussian bkg1("bkg1Gauss","bkg1 lifetime pdf",t,bkg1Offset,bkg1Resolution);
  RooProdPdf bkg1PDF("bkg1PDF", "Product of bkg1 mass and lifetime PDF", backgroundFunc, bkg1);

  //bkg2 exp * gauss conv
  RooRealVar bkg2Tau("bkg2Tau","Bkg2 lifetime",4.91124);
  RooRealVar bkg2Resolution("bkg2Res","Long lifetime resolution",0.185538);
  RooRealVar bkg2Offset("bkg1Offset","Long lifetime offset",0.0);
  RooGenericPdf bkg2Decay("bkg2Decay","Long Lifetime pdf","0.5*(1+sign(t))/bkg2Tau*exp(-t/bkg2Tau)",RooArgSet(t,bkg2Tau));
  RooGaussian bkg2ResolutionGauss("bkg2ResolutionGauss","Long lifetime Resolution PDF",t,bkg2Offset,bkg2Resolution);
  RooFFTConvPdf bkg2("bkg2","lifetime (X) gauss bkg",t,bkg2Decay,bkg2ResolutionGauss);
  RooProdPdf bkg2PDF("bkg2PDF", "Product of bkg2 mass and lifetime PDF", backgroundFunc, bkg2);
  
  RooRealVar gfrac1("gfrac1","fraction of signal",0.107061) ;
  RooRealVar gfrac2("gfrac2","fraction of background1",0.930948) ;

  RooAddPdf bkgFunc("bkgFunc","Total bkg fitting PDF", RooArgList(bkg1PDF,bkg2PDF),RooArgSet(gfrac2));
  RooAddPdf fitFunc("fitFunc","Total fitting PDF", RooArgList(signalPDF,bkgFunc),RooArgSet(gfrac1));
  RooFitResult* fitRes = fitFunc.fitTo(*data,RooFit::Save(),RooFit::NumCPU(20));

	double minNLL =  fitRes->minNll();
	
	RooAbsReal* nll = fitFunc.createNLL(*data,RooFit::NumCPU(20));

  TCanvas * c3 = new TCanvas();
  c3->Divide(2,1);
  c3->cd(1);
	RooPlot* pTauframe = sigTau.frame(0.6,0.8,20);
	RooAbsReal* pTll = nll->createProfile(sigTau);
	pTll->plotOn(pTauframe) ;
	pTauframe->Draw();

  c3->cd(2);
	RooPlot* pMframe = sigMass.frame(460,470,20);
	RooAbsReal* pMll = nll->createProfile(sigMass);
	pMll->plotOn(pMframe) ;
	pMframe->Draw();
	c3->SaveAs("profiles.pdf");
}


void combinedFit(TString inputFile="mass_life.data"){
  RooRealVar t("t","t",-2.,15.);
  RooRealVar mass("mass","mass",400.,600.);
  RooDataSet* data = RooDataSet::read(inputFile,RooArgSet(mass,t));
  t.setBins(10000,"cache"); 
	mass.setBins(10000,"cache");

		//signal
  RooRealVar sigMass("sigMass","Signal mass peak position",450,400,480) ;
  RooRealVar sigMassWidth("sigMassSigma","Signal mass peak width",15.,10.,30.);
  RooGaussian massSignalPDF("Gaussian","Mass PDF",mass,sigMass,sigMassWidth);

  RooRealVar sigTau("sigTau","Signal lifetime",0.3,0.,1.);
  RooRealVar sigResolution("sigResolution","Resolution lifetime for signal",0.5,0.,1.0);
  RooRealVar sigResolutionCenter("sigResolutionCenter","Resolution lifetime offset",0.0);
  RooGenericPdf signalDecay("signalDecay","Lifetime pdf","0.5*(1+sign(t))*0.5/sigTau*exp(-t/sigTau)",RooArgSet(t,sigTau));
  RooGaussian resolutionGauss("resolutionGauss","Lifetime Resolution PDF",t,sigResolutionCenter,sigResolution);
  RooFFTConvPdf lifetimeSignalPDF("signalConv","lifetime (X) gauss",t,signalDecay,resolutionGauss);

  RooProdPdf signalPDF("signalPDF", "Product of mass and lifetime signal PDFs", massSignalPDF, lifetimeSignalPDF);


  //background 1
  RooRealVar x1("x1","first null point",350,340,355) ;
  RooRealVar x2("x2","second null point",655,640,660) ;

  RooGenericPdf backgroundFunc("Quadratic bkg","Quadratic PDF","-0.1*(x1-mass)*(x2-mass)*(1+sign(mass-x1))*(1+sign(x2-mass))",RooArgSet(x1,x2,mass));

  RooRealVar bkg1Resolution("bkg1Res","No lifetime resolution",0.1,0.,1.0);
  RooRealVar bkg1Offset("bkg1Offset","No lifetime offset",0.0);
  RooGaussian bkg1("bkg1Gauss","bkg1 lifetime pdf",t,bkg1Offset,bkg1Resolution);
  RooProdPdf bkg1PDF("bkg1PDF", "Product of bkg1 mass and lifetime PDF", backgroundFunc, bkg1);

  //bkg2 exp * gauss conv
  RooRealVar bkg2Tau("bkg2Tau","Bkg2 lifetime",2.0,1.,10.);
  RooRealVar bkg2Resolution("bkg2Res","Long lifetime resolution",0.1,0.,1.);
  RooRealVar bkg2Offset("bkg1Offset","Long lifetime offset",0.0);
  RooGenericPdf bkg2Decay("bkg2Decay","Long Lifetime pdf","0.5*(1+sign(t))/bkg2Tau*exp(-t/bkg2Tau)",RooArgSet(t,bkg2Tau));
  RooGaussian bkg2ResolutionGauss("bkg2ResolutionGauss","Long lifetime Resolution PDF",t,bkg2Offset,bkg2Resolution);
  RooFFTConvPdf bkg2("bkg2","lifetime (X) gauss bkg",t,bkg2Decay,bkg2ResolutionGauss);
  RooProdPdf bkg2PDF("bkg2PDF", "Product of bkg2 mass and lifetime PDF", backgroundFunc, bkg2);
  
  RooRealVar gfrac1("gfrac1","fraction of signal",0.1,0.,1.) ;
  RooRealVar gfrac2("gfrac2","fraction of background1",0.5,0.,1.) ;

  RooAddPdf bkgFunc("bkgFunc","Total bkg fitting PDF", RooArgList(bkg1PDF,bkg2PDF),RooArgSet(gfrac2));
  RooAddPdf fitFunc("fitFunc","Total fitting PDF", RooArgList(signalPDF,bkgFunc),RooArgSet(gfrac1));

  RooFitResult* fitRes = fitFunc.fitTo(*data,RooFit::Save(),RooFit::NumCPU(20));


	double minNLL =  fitRes->minNll();

  TCanvas * c = new TCanvas();
  c->Divide(1,1);
	RooPlot* f = new RooPlot(sigTau,sigMass,0.6,0.8,463.,465.);
	fitRes->plotOn(f,sigTau,sigMass,"ME12VHB") ;
	c->cd(1);
	f->Draw() ;
	c->SaveAs("contour.pdf");


  TCanvas * singleFits = new TCanvas();
  singleFits->Divide(2,2);
  singleFits->cd(1);
  RooPlot* frame1 = mass.frame();
  data->plotOn(frame1);
  fitFunc.plotOn(frame1,"L") ;
  frame1->Draw();

  singleFits->cd(2);
  gPad-> SetLogy();
  RooPlot* frame2 = t.frame();
  data->plotOn(frame2);
  fitFunc.plotOn(frame2,"L") ;
  frame2->Draw();

	singleFits->cd(3);
	RooAbsReal* nll = fitFunc.createNLL(*data);
	RooPlot* tauFrame = sigTau.frame(0.62,0.72) ;
	nll->plotOn(tauFrame);
	tauFrame->SetMinimum(minNLL-1.);
	tauFrame->SetMaximum(minNLL+5);
	tauFrame->Draw();

	singleFits->cd(4);
	RooPlot* massFrame = sigMass.frame(463.5,464.5);
	nll->plotOn(massFrame);
	massFrame->SetMinimum(minNLL-1.);
	massFrame->SetMaximum(minNLL+5);
	massFrame->Draw();
  singleFits->SaveAs("fitCurves.pdf");
}