// vim: set sw=4 sts=4 filetype=cpp fdm=marker et: 

// -----------------------------------------------
//       Author: Xin Shi <Xin.Shi@cern.ch> 
//       Created:   [2013-08-15 Thu 14:54] 
// -----------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include <sys/stat.h>
#include <sys/sendfile.h>
#include <math.h>
#include <string.h>
#include <getopt.h> // passing unix-like arguments
//#include <thread> //c++11 feature
//#include <regex> //c++11 feature should be fine using gcc491.

#include <TSystem.h>
#include <TStyle.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1.h>
#include <TH2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TPad.h> 
#include <TCanvas.h> 
#include <TChain.h> 
#include <TChainElement.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TLine.h>
#include <TString.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TLorentzVector.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

#include <RooConstVar.h>
#include <RooRealVar.h>
#include <RooFormulaVar.h>
#include <RooAbsPdf.h>
#include <RooAddPdf.h>
#include <RooGaussian.h>
#include <RooBifurGauss.h>
#include <RooChebychev.h> 
#include <RooGenericPdf.h> 
#include <RooExponential.h>
#include <RooBernstein.h>
#include <RooPolynomial.h>
#include <RooPoisson.h>
#include <RooKeysPdf.h>
#include <RooExtendPdf.h>
#include <RooProdPdf.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooAbsData.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooAddition.h>
#include <RooProduct.h>
#include <RooMinuit.h>
#include <RooWorkspace.h>
#include <RooRandom.h>
#include <RooMultiVarGaussian.h>
//#include <RooMsgService.h>

#include "tools.h" 
#define PI 3.14159265358979

using namespace std; 
using namespace RooFit;

// note : 
//      1~6, 1~19 with peak regions excluded.


// Tags configration
bool redirectStdout = false;
bool is7TeVCheck = false; // Using 2011 efficiency map.
bool gKeepParam = false;
int isCDFcut = 4; // -1 for off, 1 for cdf, 2 for LHCb . 3 for 16Aug reOptimization . 4 for sel_v3p5
bool isToy=false; // For toy, no CDF cut is applied
TChain *ch=new TChain("tree");
TString plotpath="./plots";
TString iwspacepath=".";
TString iCombBkgWspacepath=".";
TString owspacepath=".";
TString idatacardpath=".";
TString odatacardpath=".";
TString ologpath=".";
double  scaleFactor=1.;
//Constants, Fit results for efficiency, etc.. //{{{
const int nSummaryBins = 2;
int summaryBins[nSummaryBins] = {11,12};
const int nQ2Ranges = 14;
char genQ2range[nQ2Ranges][128] = {"genQ2 < 2.00 && genQ2 > 1.00",//0
                                  "genQ2 < 4.30 && genQ2 > 2.00",
                                  "genQ2 < 8.68 && genQ2 > 4.30",
                                  "genQ2 <10.09 && genQ2 > 8.68",//Jpsi
                                  "genQ2 <12.86 && genQ2 >10.09",
                                  "genQ2 <14.18 && genQ2 >12.86",//Psi2S
                                  "genQ2 <16.00 && genQ2 >14.18",
                                  "genQ2 <19.00 && genQ2 >16.00",
                                  "genQ2 < 4.30 && genQ2 > 1.00",// merge 0+1
                                  "genQ2 <19.00 && genQ2 >14.18",// merge 6+7
                                  "genQ2 < 8.68 && genQ2 > 1.00",// merge 0+1+2
                                  "(genQ2 <8.68 && genQ2 > 1.00) || (genQ2 > 10.09 && genQ2 < 12.86) || (genQ2 >14.18 && genQ2 <19.00)",
                                  "genQ2 < 6.00 && genQ2 > 1.00",
                                  "genQ2 < 8.68 && genQ2 > 6.00"};
char q2range[nQ2Ranges][128] = {"Q2 < 2.00 && Q2 > 1.00",
                                "Q2 < 4.30 && Q2 > 2.00",
                                "Q2 < 8.68 && Q2 > 4.30",
                                "Q2 <10.09 && Q2 > 8.68",//3,Jpsi
                                "Q2 <12.86 && Q2 >10.09",
                                "Q2 <14.18 && Q2 >12.86",//5,Psi2S
                                "Q2 <16.00 && Q2 >14.18",
                                "Q2 <19.00 && Q2 >16.00",//7
                                "Q2 < 4.30 && Q2 > 1.00",//0+1
                                "Q2 <19.00 && Q2 >14.18",//6+7
                                "Q2 < 8.68 && Q2 > 1.00",//0+1+2
                                "(Q2 <8.68 && Q2 > 1.00) || (Q2 > 10.09 && Q2 < 12.86) || (Q2 >14.18 && Q2 <19.00)",
                                "Q2 < 6.00 && Q2 > 1.00",
                                "Q2 < 8.68 && Q2 > 6.00"};
char q2rangeLatex[nQ2Ranges][32] = {" 1.00 < q^{2} < 2.00",
                                    " 2.00 < q^{2} < 4.30",
                                    " 4.30 < q^{2} < 8.68",
                                    " 8.68 < q^{2} <10.09",//3,Jpsi
                                    "10.09 < q^{2} <12.86",
                                    "12.86 < q^{2} <14.18",//5,Psi2S
                                    "14.18 < q^{2} <16.00",
                                    "16.00 < q^{2} <19.00",
                                    " 1.00 < q^{2} < 4.30",
                                    "14.18 < q^{2} <19.00",
                                    " 1.00 < q^{2} < 8.68",
                                    " 1.00 < q^{2} <19.00",
                                    " 1.00 < q^{2} < 6.00",
                                    " 6.00 < q^{2} < 8.68"};
double q2rangedn[nQ2Ranges] = {1.00 , 2.00 , 4.30 , 8.68  , 10.09 , 12.86 , 14.18 , 16.00 , 1.00 , 14.18 , 1.00 ,  1.00 , 1.00, 6.00};
double q2rangeup[nQ2Ranges] = {2.00 , 4.30 , 8.68 , 10.09 , 12.86 , 14.18 , 16.00 , 19.00 , 4.30 , 19.00 , 8.68 , 19.00 , 6.00, 8.68};
char nTriggeredPath[3][32] = {"Triggers == 0", "Triggers == 1", "Triggers >= 1"};
char mumuMassWindow[11][512] = {"Mumumass > 0",
                                "(Mumumass > 3.096916+3.5*Mumumasserr || Mumumass < 3.096916-5.5*Mumumasserr) && (Mumumass > 3.686109+3.5*Mumumasserr || Mumumass < 3.686109-3.5*Mumumasserr)",
                                "(Mumumass < 3.096916+3*Mumumasserr && Mumumass > 3.096916-5*Mumumasserr) || (Mumumass < 3.686109+3*Mumumasserr && Mumumass > 3.686109-3*Mumumasserr)",
                                "(Mumumass*Mumumass<8.68 && Bmass-Mumumass>2.182+0.16 && Bmass-Mumumass<2.182-0.16) || (Mumumass*Mumumass>10.09 && Mumumass*Mumumass<12.86 && Bmass-Mumumass>1.593+0.06 && Bmass-Mumumass<1.593-0.06) || (Mumumass*Mumumass>14.18 && Bmass-Mumumass>1.593+0.06 && Bmass-Mumumass<1.593-0.06)",
                                "(Mumumass*Mumumass < 8.68 && ( Bmass-Mumumass < 2.182+0.16 || Bmass-Mumumass > 2.182-0.16)) || (Mumumass*Mumumass>8.68 && Mumumass*Mumumass<10.09) || (Mumumass*Mumumass > 10.09 && Mumumass < 12.86 && ( Bmass-Mumumass < 1.593+0.06 || Bmass-Mumumass > 1.593-0.06)) || (Mumumass*Mumumass>12.86 && Mumumass*Mumumass<14.08) || (Mumumass*Mumumass > 14.08 && (Bmass-Mumumass > 1.593-0.06 || Bmass-Mumumass < 1.593+0.06))",
                                "Mumumass > 0",
                                "Mumumass > 0",
                                "Mumumass > 0",
                                "Mumumass > 0",
                                "abs(Bmass-Mumumass-2.182)>0.09 && abs(Bmass-Mumumass-1.593)>0.03",
                                "abs(Bmass-Mumumass-2.182)<0.09 || abs(Bmass-Mumumass-1.593)<0.03",
};//None, sig, bkg, #Jpsi, #Psi2S, CDF, anti-CDF, LHCb, anti-LHCb, 16Aug reOptimization, anti-16Aug-reOptimization, sel_v3p5, anti-sel_v3p5
char kstarMassWindow[3][128] = {    "Kstarmass > 0",
                                    "Kstarmass > 0.792 && Kstarmass < 0.992",
                                    "(Kstarmass > 0 && Kstarmass < 0.792) || Kstarmass > 0.992"
};//None, sig, sideband
double genAfb   [nQ2Ranges]={-0.160   , -0.066   , 0.182    , 0.317    , 0.374    , 0.412    , 0.421    , 0.376    , -0.098   , 0.398    , 0.069    , 0.231 , -0.036    ,-0.036};
double genAfberr[nQ2Ranges]={0.000434 , 0.000285 , 0.000223 , 0.000383 , 0.000277 , 0.000421 , 0.000395 , 0.000422 , 0.000333 , 0.000146 , 0.000172, 0.000136, 0.000279 ,0.00279};
double genFl    [nQ2Ranges]={0.705    , 0.791    , 0.649    , 0.524    , 0.454    , 0.399    , 0.369    , 0.341    , 0.7620   , 0.355    , 0.694    , 0.544 , 0.748     ,0.748};
double genFlerr [nQ2Ranges]={0.000568 , 0.000409 , 0.000271 , 0.000420 , 0.000287 , 0.000415 , 0.000369 , 0.000361 , 0.000240 , 0.000132 , 0.000258, 0.000144, 0.000206 ,0.000208};
// Lumi = Nreco/(cross section*branch factor*filter efficiency), cross section is 49.59e9 [pb] for 8TeV and 48.44e9 [pb] for 7TeV.
// BF_BuToK*MuMu = 1.07E-6, 1.12E-6(2014)
// BF_BuToK*Jpsi = 1.43E-3, 1.44E-3(2014)
// BF_BuToK*Psi2S = 6.7E-4, 6.7 E-4(2014)
// BF_JpsToMuMu  = 5.96E-2
// BF_Psi2sToMuMu= 6.70E-4
// BF_K*ToK0Pi  = 2/3  (K* decays to Kpi)
// BF_K0ToKs  = 1/2
// BF_KsToPiPi = 2/3
double datasetLumi[5] = {19.98,37378.629,295.761,218.472,9.81};//data, BuToKstarMuMu(16281.440+21097.189), BuToKstarJpsi(118.201+177.560), BuToKstarPsi2S(63.103,155.369), JpsiX
//}}}

double readParam(int iBin, const char parName[], int iColumn, double defVal=0., double forceReturn=999.)
{//{{{
    // Remark: first value is at iColumn=0.
    if (forceReturn != 999.) return forceReturn;

    std::vector<double> output;
    char lineBuff[1024];
    char *valBuff;
    memset(lineBuff,' ',1024*sizeof(char));
    FILE *fp = fopen(TString::Format("%s/fitParameters%d.txt",idatacardpath.Data(),iBin),"r");
    if (!fp){
        printf("WARNING: readParam, missing parameter files, by default return %f.\n",defVal);
        return defVal;
    }
    while(fgets(lineBuff,1024,fp) != NULL ){
        valBuff = strtok(lineBuff," ");
        if ( strcmp(valBuff,parName) == 0 ){
            printf("INFO: readParam, matched %s!\n",valBuff);
            valBuff = strtok(NULL," ");
            while(valBuff != NULL){
                //output.push_back(stof(valBuff));//stof is a c++11 function, use other function
                if (strcmp(valBuff,"nan") == 0 || strcmp(valBuff,"inf") == 0 ){
                    output.push_back(defVal);
                }else{
                    output.push_back(std::atof(valBuff));
                }
                valBuff = strtok(NULL," ");
            }
            break;
        }
        memset(lineBuff,' ',1024*sizeof(char));
    }
    fclose(fp);
    
    if (iColumn < output.size() ){
        printf("INFO: readParam, get %s[%d]=%e\n",parName,iColumn,output.at(iColumn));
        return output.at(iColumn);
    }else{
        printf("WARNING: readParam, empty column! Return %s[%d]=defVal=%f.\n",parName,iColumn,defVal);
        return defVal;
    }
}//}}}
RooRealVar* readParam(const char parName[], const char wspacePathName[])
{//{{{
    TFile *f_wspace = new TFile(wspacePathName);
    RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
    if (wspace){
        RooRealVar *var = (RooRealVar*)wspace->var(parName);
        if (var != 0){
            return var;
        }else{
            printf("ERROR\t\t: %s cannot be found in %s\n", parName, wspacePathName);
            return 0;
        }
    }else{
        printf("ERROR\t\t: wspace cannot be found in %s\n", wspacePathName);
        return 0;
    }
}//}}}
std::string readParam(int iBin, const char parName[], string defVal="", string forceReturn="defaultForceReturn")
{//{{{
    // Remark: first value is at iColumn=0.
    if (forceReturn != "defaultForceReturn") return forceReturn;

    string output;
    char lineBuff[1024];
    char *valBuff;
    memset(lineBuff,' ',1024*sizeof(char));
    FILE *fp = fopen(TString::Format("%s/fitParameters%d.txt",idatacardpath.Data(),iBin),"r");
    if (!fp){
        printf("WARNING: readParam, missing parameter files, by default return %s.",defVal.c_str());
        return defVal;
    }
    while(fgets(lineBuff,1024,fp) != NULL ){
        valBuff = strtok(lineBuff," ");
        if ( strcmp(valBuff,parName) == 0 ){
            printf("INFO: readParam, matched %s!\n",valBuff);
            valBuff = strtok(NULL,"\n");
            output=string(valBuff);
            break;
        }
        memset(lineBuff,' ',1024*sizeof(char));
    }
    fclose(fp);
    
    if (output != ""){
        printf("INFO: readParam, get %s=%s\n",parName,output.c_str());
        return output;
    }else{
        printf("WARNING: readParam, empty item! Return %s=defVal=%s.\n",parName, defVal.c_str());
        return defVal;
    }
}//}}}
void writeParam(int iBin, const char parName[], double *val, int nVal=2, bool overwrite=true)
{//{{{
    if ( !overwrite ) return;

    struct stat fiBuff;
    FILE *fi = 0;
    if (stat(TString::Format("%s/fitParameters%d.txt",odatacardpath.Data(),iBin),&fiBuff) == 0){
        rename(TString::Format("%s/fitParameters%d.txt",odatacardpath.Data(),iBin),TString::Format("%s/fitParameters%d.txt.temp",odatacardpath.Data(),iBin));
        fi = fopen(TString::Format("%s/fitParameters%d.txt.temp",odatacardpath.Data(),iBin),"r");
    }else{
        fi = fopen(TString::Format("%s/fitParameters%d.txt.temp",odatacardpath.Data(),iBin),"w");
    }
    
    bool parExist = false;
    char lineBuff[1024];
    char *valBuff = 0;
    memset(lineBuff,' ',1024*sizeof(char));
    FILE *fp = fopen(TString::Format("%s/fitParameters%d.txt",odatacardpath.Data(),iBin),"w");
    while(fgets(lineBuff,1024,fi) != NULL ){
        valBuff = strtok(lineBuff," ");
        if ( strcmp(valBuff,parName) == 0 ){
            fprintf(fp,"%s",parName);
            int iVal = 0;
            while(iVal < nVal){
                fprintf(fp," %e",val[iVal]);
                iVal++;
            }
            fprintf(fp,"\n");
            parExist = true;
        }else{
            fprintf(fp,"%s",lineBuff);
            valBuff = strtok(NULL," ");
            while( valBuff != NULL ){
                fprintf(fp," %s",valBuff);
                valBuff = strtok(NULL," ");
            }
        }
        memset(lineBuff,' ',512*sizeof(char));
    }
    if (parExist == false){
        fprintf(fp,"%s",parName);
        int iVal = 0;
        while(iVal < nVal){
            fprintf(fp," %e",val[iVal]);
            iVal++;
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
    fclose(fi);
    remove(TString::Format("%s/fitParameters%d.txt.temp",odatacardpath.Data(),iBin));
}//}}}
void writeParam(int iBin, const char parName[], string instring, bool overwrite=true)
{//{{{
    if ( !overwrite ) return;

    struct stat fiBuff;
    FILE *fi = 0;
    if (stat(TString::Format("%s/fitParameters%d.txt",odatacardpath.Data(),iBin),&fiBuff) == 0){
        rename(TString::Format("%s/fitParameters%d.txt",odatacardpath.Data(),iBin),TString::Format("%s/fitParameters%d.txt.temp",odatacardpath.Data(),iBin));
        fi = fopen(TString::Format("%s/fitParameters%d.txt.temp",odatacardpath.Data(),iBin),"r");
    }else{
        fi = fopen(TString::Format("%s/fitParameters%d.txt.temp",odatacardpath.Data(),iBin),"w");
    }
    
    bool parExist = false;
    char lineBuff[1024];
    char *valBuff = 0;
    memset(lineBuff,' ',1024*sizeof(char));
    FILE *fp = fopen(TString::Format("%s/fitParameters%d.txt",odatacardpath.Data(),iBin),"w");
    while(fgets(lineBuff,1024,fi) != NULL ){
        valBuff = strtok(lineBuff," ");
        if ( strcmp(valBuff,parName) == 0 ){
            fprintf(fp,"%s %s\n", parName, instring.c_str());
            parExist = true;
        }else{
            fprintf(fp,"%s",lineBuff);
            valBuff = strtok(NULL," ");
            while( valBuff != NULL ){
                fprintf(fp," %s",valBuff);
                valBuff = strtok(NULL," ");
            }
        }
        memset(lineBuff,' ',512*sizeof(char));
    }
    if (parExist == false){
        fprintf(fp,"%s %s\n", parName, instring.c_str());
    }
    fclose(fp);
    fclose(fi);
    remove(TString::Format("%s/fitParameters%d.txt.temp",odatacardpath.Data(),iBin));
    return;
}//}}}
void switchRedirectStdio(const char outfile[]="_stdio", const char mode[]="a", FILE *tty=stdout)
{//{{{
    // This function works ONLY on unix-like system.
    // Or the stdout cannot be restored after running freopen.
    //
    if (!redirectStdout) return;

    struct stat fiBuff;
    if (stat("/dev/tty",&fiBuff) != 0){
        printf("WARNING\t\t: \"/dev/tty\" is NOT available\n");
        return;
    }

    if (strcmp(mode,"a")*strcmp(mode,"w")!=0){
        freopen("/dev/tty","a",stdout);
        freopen("/dev/tty","a",stderr);
        printf("ERROR\t\t: Only mode \"a\" and \"w\" are supported in switchRedirectStdio\n");
        printf("\t\t: Stdout/stderr now restored to screen\n");
    }

    if (strcmp(outfile,"_stdio")==0){
        printf("INFO\t\t: Direct stdout/stderr to screen.\n");
        freopen("/dev/tty","a",stdout);
        freopen("/dev/tty","a",stderr);
    }else{
        printf("INFO\t\t: Redirect stdout/stderr to \"%s\".\n",outfile);
        freopen(outfile,mode,tty);
    }
    return;
}//}}}

// Physically allowed ranges from AN2015_129_v14, p25.
// Transformation rule comes from AN2014_129_v14, p28.
double toUnboundedFl(double fl){
    if ( fabs(fl-0.5) > 0.5 ){
        fl = 0.5 + (fl>0.5-fl<0.5)*0.499;
    }
    double result = TMath::Tan((fl-0.5)*TMath::Pi());
    return result;
}
double toBoundedFl(double fl_ubd){
    return 0.5+TMath::ATan(fl_ubd)/TMath::Pi();
}
double toUnboundedAfb(double afb, double fl){
    double rad = 2./3.*afb/(1.-fl);
    if ( fabs(rad) > 0.5 ){
        printf("DEBUG\t\t: [toUnboundedAfb] Input rad(afb=%f,fl=%f) = %f > 0.5\n",afb,fl,rad);
        rad = ((rad>0)-(rad<0))*0.499;
    }
    return TMath::Tan(rad*TMath::Pi());
}
double toBoundedAfb(double afb_ubd, double fl_ubd){
    return 3./2.*(0.5-TMath::ATan(fl_ubd)/TMath::Pi())*TMath::ATan(afb_ubd)/TMath::Pi();
}
double toTransformedAs(double fs, double fl, double as)
{
    if (fs*(fs-1) > 0 || fl == 0 || as==0) return 0;// No interference
    return as/(2*sqrt(3*fs*(1-fs)*fl)*0.89);
}
double toOriginAs(double as_tr, double fs, double fl_ubd)
{
    double fl=toBoundedFl(fl_ubd);
    return 2*sqrt(3*fs*(1-fs)*fl)*0.89*as_tr;
}

bool scanAfbFlPositivePdf(double afb, double fl, bool fineScan=false)
{//{{{
    if ( fl < 0.95 && fabs(afb)<(0.95-fl)*0.75 ) return true;

    // Create test function to find possible domain for fl/afb. CosThetaL as x, CosThetaK as y.
    TString f2_format = "2.*[0]*y**2*(1.-x**2)+0.5*(1-[0])*(1.-y**2)*(1.+x**2)+4./3*[1]*(1-y**2)*x";
    TF2 *f2_model = new TF2("f2_model", f2_format.Data(),-1.,1.,-1.,1.);
    f2_model->SetTitle("PDF value;cos#theta_{L};cos#theta_{K}");

    f2_model->SetParameter(0,fl);
    f2_model->SetParameter(1,afb);
    int nScanSteps = 1000;
    if (fineScan) nScanSteps *= 10;
    bool isPositivePDF=true;
    if (f2_model->Eval(1.,0.) < 0 || f2_model->Eval(-1,0) < 0){
        isPositivePDF = false;
    }else{
        for (int i = 0; i < nScanSteps; ++i) {//cosThetaL
            if (afb*2./nScanSteps*i-1. > 0) continue;
            for (int j = 0; j < nScanSteps; ++j) {//cosThetaK
                if (f2_model->Eval(2./nScanSteps*i-1.,2./nScanSteps*j-1.) < 0.){
                    isPositivePDF = false;
                    break;
                }
            }
            if (!isPositivePDF) break;
        }
    }

    return isPositivePDF;
}//}}}

void scanAfbFlPositivePdf()
{//{{{
    // Create test function to find possible domain for fl/afb. CosThetaL as x, CosThetaK as y.
    TCanvas *canvas= new TCanvas("canvas");
    TH2F *h2_minPdfValue = new TH2F("h2_minPdfValue","",200,-1,1,200,0,1);
    h2_minPdfValue->SetStats(false);
    h2_minPdfValue->SetXTitle("A_{FB}");
    h2_minPdfValue->SetYTitle("F_{L}");

    for( int xBin = 1; xBin <= h2_minPdfValue->GetNbinsX(); xBin++){//afb
        for( int yBin = 1; yBin <= h2_minPdfValue->GetNbinsY(); yBin++){//fl
            bool isPositivePDF = scanAfbFlPositivePdf((yBin-0.5)/h2_minPdfValue->GetNbinsY(),(2.*xBin-1.)/h2_minPdfValue->GetNbinsX()-1);
            if(isPositivePDF){
                h2_minPdfValue->SetBinContent(xBin,yBin,1);
            }else{
            }
        }
    }

    // Draw contour
    h2_minPdfValue->Draw("COL");
    canvas->Update();
    canvas->SaveSource(TString::Format("%s/scanAfbFlPositivePdf.cc",plotpath.Data()));
    canvas->Print(TString::Format("%s/scanAfbFlPositivePdf.pdf",plotpath.Data()));
    return;
}//}}}

bool printPdf2D(RooAbsPdf &pdf, RooRealVar &cosK, RooRealVar &cosL, const char oName[]="")
{//{{{
    bool isPositivePDF = true;
    TCanvas *canvas = new TCanvas();
    TH2D *h2 = new TH2D("h2","",200,-1,1,200,-1,1);

    for( double iCosThetaK=-1; iCosThetaK < 1; iCosThetaK+=0.01){
        cosK.setVal(iCosThetaK);
        for( double iCosThetaL=-1; iCosThetaL < 1; iCosThetaL+=0.01){
            cosL.setVal(iCosThetaL);
            double val = pdf.getVal(RooArgSet(cosL,cosK));
            if (val < 0) isPositivePDF = false;
            h2->Fill(iCosThetaL,iCosThetaK,val);
        }
    }
    
    h2->SetMinimum(0.);
    h2->Draw("COLZ");
    if (strcmp(oName,"") != 0){
        canvas->Print(TString::Format("%s/%s.pdf",plotpath.Data(),oName));
    }

    return isPositivePDF;
}//}}}

TF2  *f2_fcn = NULL;
double model_2D(double *x, double *par)
{//{{{
    double xx = x[0];
    double yy = x[1];
    for (int i = 0; i < f2_fcn->GetNpar(); i++) f2_fcn->SetParameter(i,par[i]);
    return f2_fcn->Eval(xx,yy);
}//}}}

TH2F *h2_fcn = NULL;
// nParameters, ???, return fcn value, parameter array, strategy
void fcn_binnedChi2_2D(int &npar, double *gin, double &f, double *par, int iflag)
{//{{{
    f=0;
    for (int i = 1; i <= h2_fcn->GetNbinsX(); i++) {
        for (int j = 1; j <= h2_fcn->GetNbinsY(); j++) {
            int gBin = h2_fcn->GetBin(i,j);
            double x[2] = {h2_fcn->GetXaxis()->GetBinCenter(i),h2_fcn->GetYaxis()->GetBinCenter(j)};
            double measure  = h2_fcn->GetBinContent(gBin);
            double error    = h2_fcn->GetBinError(gBin);
            
            //// Naively test using center value
            //double func     = model_2D(x, par);//Take center value
            //double delta    = (measure-func)/error;
            //if (measure != 0) 
            //    f+=delta*delta;
            
            //// Real run using integral
            for (int k = 0; k < f2_fcn->GetNpar(); k++){//nPar MUST be the same value as f2_fcn
                f2_fcn->SetParameter(k,par[k]);
            }
            double xi = h2_fcn->GetXaxis()->GetBinLowEdge(i);
            double xf = h2_fcn->GetXaxis()->GetBinUpEdge(i);
            double yi = h2_fcn->GetYaxis()->GetBinLowEdge(j);
            double yf = h2_fcn->GetYaxis()->GetBinUpEdge(j);
            //f2_fcn->SetRange(xi,xf,yi,yf);
            //double minX, minY;
            //f2_fcn->GetMinimumXY(minX,minY);
            //if (f2_fcn->Eval(minX,minY) < 0){
            //    f += 100;
            //}else{
                f += pow( (f2_fcn->Integral(xi,xf,yi,yf)/(xf-xi)/(yf-yi)-measure)/error,2);
            //}
        }
    }
    //printf("FCN in calls = %f\n",f);
}//}}}

//_________________________________________________________________________________

void bmass(int iBin, const char outfile[] = "bmass") // Somehow DEPRECATED
{//{{{
  bool test = false; 
  
  ch->SetBranchStatus("*",0);
  ch->SetBranchStatus("Bmass"         , 1);
  ch->SetBranchStatus("Mumumass"      , 1);
  ch->SetBranchStatus("Mumumasserr"   , 1);
  ch->SetBranchStatus("Q2"            , 1);
  ch->SetBranchStatus("Triggers"      , 1);
  RooRealVar Bmass("Bmass", "B^{+/-} mass(GeV/c^{2})", 5.1, 5.6) ;
  RooRealVar Q2("Q2","q^{2}",0.5,20.);
  RooRealVar Mumumass("Mumumass","M^{#mu#mu}",0.,10.);
  RooRealVar Mumumasserr("Mumumasserr","Error of M^{#mu#mu}",0.,10.);
  RooRealVar Triggers("Triggers","",0,100);
  int mumuMassWindowBin = 1+2*isCDFcut;
  if (iBin==3 || iBin==5 || isCDFcut < 0) mumuMassWindowBin = 0; // no cut
  RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass, Mumumass, Mumumasserr, Triggers),TString::Format("(%s) && (%s) && (%s)",nTriggeredPath[2], q2range[iBin],mumuMassWindow[mumuMassWindowBin]),0);
  if (data->sumEntries() == 0){
      return;
  }

  data->Print();

  // Create model and dataset
  // -------------------------------------------------------------------------
  // Gaussian signal 
  RooRealVar mean("mean","mean of gaussians", 5.27, 5.23, 5.32) ;
  RooRealVar sigma1("sigma1","width of Gaussian1", 0.0285, 0.01, 0.05) ;
  RooRealVar sigma2("sigma2","width of Gaussian2", 0.08, 0.05, 0.35) ;
  RooRealVar sigM_frac("sigM_frac","fraction of Gaussians",0,0.,1.);
  RooGaussian sigGauss1("siggauss1","Signal component", Bmass, mean, sigma1) ;  
  RooGaussian sigGauss2("siggauss2","Signal component", Bmass, mean, sigma2) ;  
  RooAddPdf sig("sig","sig",RooArgList(sigGauss1,sigGauss2),RooArgList(sigM_frac));

  // Build Chebychev polynomial p.d.f.  
  RooRealVar a0("a0", "constant", 0.5, -1, 1) ;
  RooRealVar a1("a1", "linear", 0.6, -1, 1) ;
  RooRealVar a2("a2", "quadratic", 0.1, -1, 1) ;
  RooChebychev bkg("bkg", "Background", Bmass, RooArgSet(a0, a1, a2)) ;

  // Construct signal+background PDF
  RooRealVar nsig("nsig", "number of signal events", 4648, 0, 1E8); 
  RooRealVar nbkg("nbkg", "number of background events", 21472, 0, 1E8);
  RooAddPdf  model("model", "g+c", RooArgList(bkg, sig), RooArgList(nbkg, nsig)) ;
  
  // Print structure of composite p.d.f.
  model.Print("t") ;

  // Fit model to data, save fitresult 
  // ------------------------------------------------------------------------
  RooFitResult* fitres; 
  if (! test) {
    fitres = model.fitTo(*data, Extended(kTRUE), Minos(kTRUE), Save(kTRUE)) ;
    fitres->Print("v"); 
  }
  
  // Plot model 
  // ---------------------------------------------------------
  TString title = "B^{+/-} mass";
  int nbins = 20; 
  RooPlot* frame = Bmass.frame(Title(title), Bins(nbins));
  data->plotOn(frame) ;
  model.plotOn(frame) ;

  // Overlay the background component of model with a dashed line
  model.plotOn(frame,Components("bkg"), LineStyle(kDashed), LineColor(2)) ;

  // Draw the frame on the canvas
  TCanvas *c = new TCanvas("c"); 
  set_root_style(); 
  c->UseCurrentStyle() ;

  gPad->SetLeftMargin(0.15) ;
  frame->GetYaxis()->SetTitleOffset(1.7) ; 
  frame->Draw();

  TPaveText* paveText = new TPaveText(0.75, 0.82, 1., 1., "NDC"); 
  paveText->SetBorderSize(0);
  paveText->SetFillColor(kWhite);
  paveText->AddText(Form("nsig = %.0f #pm %.0f " , nsig.getVal()  , nsig.getError())); 
  paveText->AddText(Form("nbkg = %.0f #pm %.0f " , nbkg.getVal()  , nbkg.getError())); 
  paveText->AddText(Form("mean = %.3f #pm %.3f " , mean.getVal()  , mean.getError())); 
  paveText->AddText(Form("sigma1 = %.3f #pm %.3f ", sigma1.getVal(), sigma1.getError())); 
  paveText->AddText(Form("sigma2 = %.3f #pm %.3f ", sigma2.getVal(), sigma2.getError())); 
  paveText->AddText(Form("frac = %.3f #pm %.3f ", sigM_frac.getVal(), sigM_frac.getError())); 
  paveText->Draw(); 

  c->SaveSource(TString::Format("%s/%s_bin%d.cc",plotpath.Data(),outfile,iBin));
  c->Print(TString::Format("%s/%s_bin%d.pdf",plotpath.Data(),outfile,iBin));
  
  // Persist fit result in root file 
  // -------------------------------------------------------------
  //TFile resf(TString::Format("%s/%s.root",plotpath.Data(),outfile), "RECREATE") ;
  //gPad->Write("plot"); 
  //if (! test) fitres->Write("fitres") ;
  //resf.Close() ;

  // In a clean ROOT session retrieve the persisted fit result as follows:
  // RooFitResult* r = gDirectory->Get("fitres") ;
  
  delete paveText; 
  delete c;

}//}}}

//_________________________________________________________________________________

std::vector<double> fl_gen_bin(int iBin, const char outfile[] = "fl_gen")
{//{{{
  // From fomula (2) in LHCb 2012 PRL108, 181806(2012)
  // integrated over theta_l and phi: 
  // 
  // 1/Gamma * d^2 Gamma/d cos(theta_K) dq^2 = 3/2 * F_L cos^2(theta_K)
  // + 3/4(1-F_L)(1-cos^2theta_K)
  // 

  bool test = false;

  RooRealVar genCosThetaK("genCosThetaK", "cos#theta_{K}", -1, 1);
  RooRealVar genQ2("genQ2","q^{2}",0.5,20.);
  RooRealVar Q2("Q2","q^{2}",0.5,20.);
  RooRealVar fl("fl", "F_{L}", 0.5, -0.5, 1.5);

  RooGenericPdf f("f", "1.5*fl*genCosThetaK*genCosThetaK+0.75*(1-fl)*(1-genCosThetaK*genCosThetaK)", RooArgSet(genCosThetaK,fl));
  RooDataSet* data;
  
  if (test){
      fl.setVal(0.5);
      data = f.generate(RooArgSet(genCosThetaK,Q2), 10000);
  }else{
      data = new RooDataSet("data","data",ch,RooArgSet(genCosThetaK,genQ2),genQ2range[iBin],0);
  }
  
  f.fitTo(*data,Extended(kTRUE), Minos(kTRUE)); 

  RooPlot* framecosk = genCosThetaK.frame(); 
  data->plotOn(framecosk); 
  f.plotOn(framecosk); 

  // Draw the frame on the canvas
  TCanvas *c = new TCanvas("c"); 
  framecosk->SetTitle("");
  framecosk->Draw();

  TLatex *t1 = new TLatex();
  t1->SetNDC();
  t1->DrawLatex(.40,.85,TString::Format("%s",genQ2range[iBin]));
  t1->DrawLatex(.40,.79,TString::Format("F_{L}=%5.3f#pm%5.3f",fl.getVal(),fl.getError()));

  c->SaveSource(TString::Format("%s/%s_bin%d.cc",plotpath.Data(),outfile,iBin));
  c->Print(TString::Format("%s/%s_bin%d.pdf",plotpath.Data(),outfile,iBin));

  delete c;
  delete t1;
  delete data;

  std::vector<double> outvect;
  outvect.push_back(fl.getVal());
  outvect.push_back(fl.getError());
  return outvect;
}//}}}

void fl_gen(const char outfile[] = "fl_gen")
{//{{{

    TCanvas *c = new TCanvas();
    TH2F *frame = new TH2F("frame","",18,1,19,10,0,1.5);
    frame->SetStats(kFALSE);
    frame->SetXTitle("q^{2} [(GeV)^{2}]");
    frame->SetYTitle("F_{L}");
    frame->Draw();

    double x[8]={1.5,3.15,6.49,9.385,11.475,13.52,15.09,17.5};
    double xerr[8]={0.5,1.15,2.09,0.705,1.385,0.66,0.91,1.5};
    double yfl[8],yerrfl[8];

    std::vector<double> vbin;
    for(int ibin = 0; ibin < 8; ibin++){
        vbin = fl_gen_bin(ibin);
        yfl[ibin]       =vbin.at(0);
        yerrfl[ibin]    =vbin.at(1);
    }
    
    // Check input data
    for(int ibin = 0; ibin < 8; ibin++){
        printf("yfl [%d]=%6.4f +- %6.4f\n",ibin,yfl[ibin],yerrfl[ibin]);
    }

    TGraphAsymmErrors *g_fl  = new TGraphAsymmErrors(8,x,yfl,xerr,xerr,yerrfl,yerrfl);
    g_fl->Draw("P*");
    c->SaveSource(TString::Format("%s/%s.cc",plotpath.Data(),outfile));
    c->Print(TString::Format("%s/%s.pdf",plotpath.Data(),outfile));

    delete g_fl;
    delete frame;
    delete c;
}//}}}

//_________________________________________________________________________________

std::vector<double> angular_gen_bin(int iBin, const char outfile[] = "angular_gen")
{//{{{

    RooRealVar genCosThetaK("genCosThetaK", "cos#theta_{K}", -1., 1.);
    RooRealVar genCosThetaL("genCosThetaL", "cos#theta_{L}", -1., 1.);
    RooRealVar genQ2("genQ2","q^{2}",0.5,20.);
    RooRealVar fl("fl", "F_{L}", genFl[iBin], 0.1, 0.9);
    RooRealVar afb("afb", "A_{FB}", genAfb[iBin], -0.75, 0.75);
    RooRealVar fs("fs","F_{S}",0.,-0.1,0.1); //Very close to 0.
    RooRealVar as("as","A_{S}",0.01,-1,1.);
    //fs.setConstant(kTRUE);
    //as.setConstant(kTRUE);

    RooRealVar nsig("nsig","nsig",1E6,1E2,1E9);
    RooRealVar nbkg("nbkg","nbkg",10,0.1,1E6);
    
    RooGenericPdf f_sig("f_sig", "9/16*((2/3*fs+4/3*as*genCosThetaK)*(1-genCosThetaL*genCosThetaL)+(1-fs)*(2*fl*genCosThetaK*genCosThetaK*(1-genCosThetaL*genCosThetaL)+1/2*(1-fl)*(1-genCosThetaK*genCosThetaK)*(1+genCosThetaL*genCosThetaL)+4/3*afb*(1-genCosThetaK*genCosThetaK)*genCosThetaL))", RooArgSet(genCosThetaK,genCosThetaL,fl,afb,fs,as));
    RooExtendPdf f("f","",f_sig,nsig);
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(genCosThetaK,genCosThetaL,genQ2),genQ2range[iBin],0);

    RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Minos(kFALSE),NumCPU(4));

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");
    RooPlot* framecosk = genCosThetaK.frame(); 
    data->plotOn(framecosk,Binning(100)); 
    f.plotOn(framecosk); 

    framecosk->SetTitle("");
    framecosk->SetMinimum(0);
    framecosk->Draw();

    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = -0.5;
    switch(iBin){
        case 10:
            fixNDC = 0.;
            break;
        default:
            break;
    }
    //t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",genQ2range[iBin]));
    t1->DrawLatex(.35,.80+fixNDC,TString::Format("F_{L}=%5.3f#pm%8.6f",fl.getVal(),fl.getError()));
    t1->DrawLatex(.35,.74+fixNDC,TString::Format("A_{FB}=%5.3f#pm%8.6f",afb.getVal(),afb.getError()));
    c->SaveSource(TString::Format("%s/%s_cosk_bin%d.cc",plotpath.Data(),outfile,iBin));
    c->Print(TString::Format("%s/%s_cosk_bin%d.pdf",plotpath.Data(),outfile,iBin));

    //
    RooPlot* framecosl = genCosThetaL.frame(); 
    data->plotOn(framecosl,Binning(100)); 
    f.plotOn(framecosl); 

    framecosl->SetTitle("");
    framecosl->SetMinimum(0);
    framecosl->Draw();

    fixNDC = -0.5;
    switch(iBin){
        default:
            break;
    }
    //t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",genQ2range[iBin]));
    t1->DrawLatex(.35,.80+fixNDC,TString::Format("F_{L}=%5.3f#pm%8.6f",fl.getVal(),fl.getError()));
    t1->DrawLatex(.35,.74+fixNDC,TString::Format("A_{FB}=%5.3f#pm%8.6f",afb.getVal(),afb.getError()));
    c->Update();
    c->SaveSource(TString::Format("%s/%s_cosl_bin%d.cc",plotpath.Data(),outfile,iBin));
    c->Print(TString::Format("%s/%s_cosl_bin%d.pdf",plotpath.Data(),outfile,iBin));

    // Make 2-D plot
    TH1 *h1 = data->createHistogram("genCosThetaL,genCosThetaK", 100, 100);
    h1->Draw("LEGO2");
    c->Update();
    c->SaveSource(TString::Format("%s/%s_2D_bin%d.cc",plotpath.Data(),outfile,iBin));
    c->Print(TString::Format("%s/%s_2D_bin%d.pdf",plotpath.Data(),outfile,iBin));

    // clear
    delete t1;
    delete c;
    delete data;

    //write output
    double outputp[4] = {0,0,0,0};
    outputp[0] = fl.getVal();
    outputp[1] = fl.getError();
    outputp[2] = fl.getErrorLo();
    outputp[3] = fl.getErrorHi();
    writeParam(iBin,"fl_gen",outputp,4);
    outputp[0] = afb.getVal();
    outputp[1] = afb.getError();
    outputp[2] = afb.getErrorLo();
    outputp[3] = afb.getErrorHi();
    writeParam(iBin,"afb_gen",outputp,4);

    std::vector<double> output;
    output.push_back(fl.getVal());
    output.push_back(fl.getError());
    output.push_back(afb.getVal());
    output.push_back(afb.getError());
    return output;

}//}}}

void angular_gen(const char outfile[] = "angular_gen")
{//{{{
    bool doFit = true; // Turn to true if you want to fit again.
    printf("INFO\t\t: You'll need UNFILTERED MC to be the input data.\n");

    TCanvas *c = new TCanvas();

    int nWorkBins = 3;
    int workBins[] = {10,4,9};
    double x[nWorkBins];
    double xerr[nWorkBins];
    double yafb[nWorkBins],yerrafbLo[nWorkBins],yerrafbHi[nWorkBins],yfl[nWorkBins],yerrflLo[nWorkBins],yerrflHi[nWorkBins];
    for(int iBin = 0; iBin < nWorkBins; iBin++){
        x[iBin] = (q2rangeup[workBins[iBin]]+q2rangedn[workBins[iBin]])/2;
        xerr[iBin] = (q2rangeup[workBins[iBin]]-q2rangedn[workBins[iBin]])/2;
    }

    if (doFit){
        for(int ibin = 0; ibin < nWorkBins; ibin++){
            angular_gen_bin(workBins[ibin],outfile);
        }
    }

    // plotting
    for(int ibin = 0; ibin < nWorkBins; ibin++){
        yfl[ibin] = -100;
        yerrflLo[ibin] = 0;
        yerrflHi[ibin] = 0;
        yafb[ibin] = -100;
        yerrafbLo[ibin] = 0;
        yerrafbHi[ibin] = 0;
        yfl[ibin]           = readParam(workBins[ibin],"fl_gen",0);
        yerrflLo[ibin]      = readParam(workBins[ibin],"fl_gen",2);
        yerrflHi[ibin]      = readParam(workBins[ibin],"fl_gen",3);
        if (yerrflHi[ibin] == -1){
            yerrflHi[ibin] = 0;
        }
        if (yerrflLo[ibin] == -1){
            yerrflLo[ibin] = 0;
        }
        if (yerrflHi[ibin] == yerrflLo[ibin]){
            yerrflHi[ibin] = 0;
            yerrflLo[ibin] = 0;
        }
        yafb[ibin]          = readParam(workBins[ibin],"afb_gen",0);
        yerrafbLo[ibin]     = readParam(workBins[ibin],"afb_gen",2);
        yerrafbHi[ibin]     = readParam(workBins[ibin],"afb_gen",3);
        if (yerrafbHi[ibin] == -1){
            yerrafbHi[ibin] = 0;
        }
        if (yerrafbLo[ibin] == -1){
            yerrafbLo[ibin] = 0;
        }
        if (yerrafbHi[ibin] == yerrafbLo[ibin]){
            yerrafbHi[ibin] = 0;
            yerrafbLo[ibin] = 0;
        }
        printf("yafb[%d]=%6.4f + %6.4f - %6.4f\n",ibin,yafb[ibin],yerrafbHi[ibin],yerrafbLo[ibin]);
        printf("yfl [%d]=%6.4f + %6.4f - %6.4f\n",ibin,yfl[ibin],yerrflHi[ibin],yerrflLo[ibin]);
    }

    TGraphAsymmErrors *g_fl  = new TGraphAsymmErrors(nWorkBins,x,yfl,xerr,xerr,yerrflLo,yerrflHi);
    g_fl->SetTitle("");
    g_fl->GetXaxis()->SetTitle("q^{2} [(GeV)^{2}]");
    g_fl->GetYaxis()->SetTitle("F_{L}");
    g_fl->GetYaxis()->SetRangeUser(0,1);
    g_fl->SetFillColor(2);
    g_fl->SetFillStyle(3001);
    g_fl->Draw("a2");
    g_fl->Draw("P");
    c->SaveSource(TString::Format("%s/%s_fl.cc",plotpath.Data(),outfile));
    c->Print(TString::Format("%s/%s_fl.pdf",plotpath.Data(),outfile));
    c->Clear();

    TGraphAsymmErrors *g_afb = new TGraphAsymmErrors(nWorkBins,x,yafb,xerr,xerr,yerrafbLo,yerrafbHi);
    g_afb->SetTitle("");
    g_afb->GetXaxis()->SetTitle("q^{2} [(GeV)^{2}]");
    g_afb->GetYaxis()->SetTitle("A_{FB}");
    g_afb->GetYaxis()->SetRangeUser(-1,1);
    g_afb->SetFillColor(2);
    g_afb->SetFillStyle(3001);
    g_afb->Draw("a2");
    g_afb->Draw("P");
    c->Update();
    c->SaveSource(TString::Format("%s/%s_afb.cc",plotpath.Data(),outfile));
    c->Print(TString::Format("%s/%s_afb.pdf",plotpath.Data(),outfile));
}//}}}

//_________________________________________________________________________________
//Fit parameters of acceptance and efficiency using TMinuit instead of RooFit.

std::vector<double> acceptance(int iBin) // acceptance. Not used in standard procedure, just for check...
{//{{{

    double accUpperBound = 0.09;
    double gQ2 = 0;
    double gCosThetaK = 0;
    double gCosThetaL = 0;
    double gmuppt = 0;
    double gmupeta= 0;
    double gmumpt = 0;
    double gmumeta= 0;

    ch->SetBranchStatus("*",0);
    ch->SetBranchStatus("genQ2"         , 1);
    ch->SetBranchStatus("genCosTheta*"  , 1);
    ch->SetBranchStatus("genMu*"        , 1);
    ch->SetBranchAddress("genQ2"        , &gQ2);
    ch->SetBranchAddress("genCosThetaK" , &gCosThetaK);
    ch->SetBranchAddress("genCosThetaL" , &gCosThetaL);
    ch->SetBranchAddress("genMupPt"     , &gmuppt);
    ch->SetBranchAddress("genMupEta"    , &gmupeta);
    ch->SetBranchAddress("genMumPt"     , &gmumpt);
    ch->SetBranchAddress("genMumEta"    , &gmumeta);

    // Fill histograms
    int nbinsL = 20;
    int nbinsK = 20;
    TH2F h2_ngen("h2_ngen","h2_ngen",nbinsL,-1.,1,nbinsK,-1,1);
    TH2F h2_nacc("h2_nacc" ,"h2_nacc" ,nbinsL,-1,1,nbinsK,-1,1); 
    // Read data
    for (int entry = 0; entry < ch->GetEntries(); entry++) {
        ch->GetEntry(entry);
        if (gQ2 > q2rangeup[iBin] || gQ2 < q2rangedn[iBin]) continue;
        if (gQ2 > q2rangeup[3] && gQ2 < q2rangedn[3]) continue;//jpsi
        if (gQ2 > q2rangeup[5] && gQ2 < q2rangedn[5]) continue;//psi2s
        h2_ngen.Fill(gCosThetaL,gCosThetaK);
        if ( fabs(gmumeta) < 2.3 && fabs(gmupeta) < 2.3 && gmumpt > 2.8 && gmuppt > 2.8 ) h2_nacc.Fill(gCosThetaL,gCosThetaK);
    }
    
    // Calculate acceptance
    TH2F h2_acc("h2_acc","",nbinsL,-1,1,nbinsK,-1,1);
    h2_acc.SetAxisRange(0.,1.,"Z");
    for (int i = 1; i <= nbinsL; i++) {
        for (int j = 1; j <= nbinsK; j++) {
            // Generate toy model
            //double ii = -1. + i*2./nbinsL - 1./nbinsL;
            //double jj = -1. + j*2./nbinsK - 1./nbinsK;
            //TF2 f2_gen("f2_gen","(0.06+0.02*(3*y**2-1)/2)-(0.03+0.03*(3*y**2-1)/2)*x**2+(0.005+0.003*(3*y**2-1)/2)*x**4",-1.,1.,-1.,1.);
            //h2_ngen.SetBinContent(i,j,400);
            //h2_nacc.SetBinContent(i,j,400*f2_gen.Eval(ii,jj));
            
            // Fill acceptance
            if (h2_ngen.GetBinContent(i,j) == 0) {
                printf("WARNING: Acceptance(%d,%d)=%f/%f\n",i,j,h2_nacc.GetBinContent(i,j),h2_ngen.GetBinContent(i,j));
                h2_acc.SetBinContent(i,j,0.);
                h2_acc.SetBinError(i,j,1.);
            }else{
                h2_acc.SetBinContent(i,j,h2_nacc.GetBinContent(i,j)/h2_ngen.GetBinContent(i,j));
                if (h2_nacc.GetBinContent(i,j) != 0){
                    h2_acc.SetBinError(i,j,sqrt(h2_acc.GetBinContent(i,j)*(1.-h2_acc.GetBinContent(i,j))/h2_ngen.GetBinContent(i,j)));
                }else{
                    h2_acc.SetBinError(i,j,sqrt(0.05/h2_ngen.GetBinContent(i,j)));
                }
                //printf("INFO: Angular bin(%d,%d)= %f +- %f ( %f / %f).\n",i,j,h2_acc.GetBinContent(i,j),h2_acc.GetBinError(i,j),h2_nacc.GetBinContent(i,j),h2_ngen.GetBinContent(i,j));
            }
        }
    }
    printf("INFO: h2_acc built.\n");
    
    // Using pure TMinuit
    int nPar = 20;
    TMinuit *gMinuit = new TMinuit(nPar);
    h2_fcn = &h2_acc;
    gMinuit->SetFCN(fcn_binnedChi2_2D);
    
    TF2 f2_model("f2_model","([0]+[1]*y+[2]*(3*y**2-1)/2+[3]*(5*y**3-3*y)/2)+([4]+[5]*y+[6]*(3*y**2-1)/2+[7]*(5*y**3-3*y)/2)*x**2+([8]+[9]*y+[10]*(3*y**2-1)/2+[11]*(5*y**3-3*y)/2)*x**3+([12]+[13]*y+[14]*(3*y**2-1)/2+[15]*(5*y**3-3*y)/2)*x**4+([16]+[17]*y+[18]*(3*y**2-1)/2+[19]*(5*y**3-3*y)/2)*x**6",-1.,1.,-1.,1.);
    f2_fcn = &f2_model;
    gMinuit->DefineParameter( 0, "k0l0",  .06,  1E-3,    -1E+0, 1E+2);
    gMinuit->DefineParameter( 1, "k1l0",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter( 2, "k2l0",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter( 3, "k3l0",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter( 4, "k0l2",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter( 5, "k1l2",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter( 6, "k2l2",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter( 7, "k3l2",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter( 8, "k0l3",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter( 9, "k1l3",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter(10, "k2l3",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter(11, "k3l3",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter(12, "k0l4",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter(13, "k1l4",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter(14, "k2l4",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter(15, "k3l4",   0.,  1E-3,    -1E+1, 1E+2);
    gMinuit->DefineParameter(16, "k0l6",   0.,  1E-3,    -1E+2, 1E+4);
    gMinuit->DefineParameter(17, "k1l6",   0.,  1E-3,    -1E+2, 1E+4);
    gMinuit->DefineParameter(18, "k2l6",   0.,  1E-3,    -1E+2, 1E+4);
    gMinuit->DefineParameter(19, "k3l6",   0.,  1E-3,    -1E+2, 1E+4);
    if (iBin == 0) {
        gMinuit->Command("SET PARM 9 0");
        gMinuit->Command("SET PARM 10 0");
        gMinuit->Command("SET PARM 11 0");
        gMinuit->Command("SET PARM 12 0");
        gMinuit->Command("FIX 9");
        gMinuit->Command("FIX 10");
        gMinuit->Command("FIX 11");
        gMinuit->Command("FIX 12");
    }else if (iBin == 1) {
        gMinuit->Command("SET PARM 9 0");
        gMinuit->Command("SET PARM 10 0");
        gMinuit->Command("SET PARM 11 0");
        gMinuit->Command("SET PARM 12 0");
        gMinuit->Command("SET PARM 17 0");
        gMinuit->Command("SET PARM 18 0");
        gMinuit->Command("SET PARM 19 0");
        gMinuit->Command("SET PARM 20 0");
        gMinuit->Command("FIX 9");
        gMinuit->Command("FIX 10");
        gMinuit->Command("FIX 11");
        gMinuit->Command("FIX 12");
        gMinuit->Command("FIX 17");
        gMinuit->Command("FIX 18");
        gMinuit->Command("FIX 19");
        gMinuit->Command("FIX 20");
    }else if (iBin > 1 && iBin < 6 ) {
        gMinuit->Command("SET PARM 17 0");
        gMinuit->Command("SET PARM 18 0");
        gMinuit->Command("SET PARM 19 0");
        gMinuit->Command("SET PARM 20 0");
        gMinuit->Command("FIX 17");
        gMinuit->Command("FIX 18");
        gMinuit->Command("FIX 19");
        gMinuit->Command("FIX 20");
    }else{
        gMinuit->Command("SET PARM 13 0");
        gMinuit->Command("SET PARM 14 0");
        gMinuit->Command("SET PARM 15 0");
        gMinuit->Command("SET PARM 16 0");
        gMinuit->Command("SET PARM 17 0");
        gMinuit->Command("SET PARM 18 0");
        gMinuit->Command("SET PARM 19 0");
        gMinuit->Command("SET PARM 20 0");
        gMinuit->Command("FIX 13");
        gMinuit->Command("FIX 14");
        gMinuit->Command("FIX 15");
        gMinuit->Command("FIX 16");
        gMinuit->Command("FIX 17");
        gMinuit->Command("FIX 18");
        gMinuit->Command("FIX 19");
        gMinuit->Command("FIX 20");
    }
    
    gMinuit->Command("MINI");
    gMinuit->Command("MINI");
    gMinuit->Command("IMPROVE");
    gMinuit->Command("MINOS");

    double arrPar[nPar];
    double arrParErr[nPar];
    for (int iPar = 0; iPar < nPar; iPar++) gMinuit->GetParameter(iPar,arrPar[iPar],arrParErr[iPar]);
    for (int iPar = 0; iPar < nPar; iPar++) f2_model.SetParameter(iPar,arrPar[iPar]);

    // Prepare draw
    double chi2Val=0;
    fcn_binnedChi2_2D(nPar, 0, chi2Val, arrPar, 0);
    printf("Chi2=%f \n",chi2Val);
    
    TCanvas canvas("canvas");
    TLatex *latex = new TLatex();
    latex->SetNDC();
    
    // Draw efficiency
    h2_acc.SetStats(0);
    h2_acc.SetMinimum(0.);
    h2_acc.SetMaximum(accUpperBound);
    h2_acc.SetTitleOffset(2,"XY");
    h2_acc.SetXTitle("genCosThetaL");
    h2_acc.SetYTitle("genCosThetaK");
    h2_acc.SetZTitle("Acceptance");
    h2_acc.Draw("LEGO2");
    latex->DrawLatex(0.35,0.95,TString::Format("Acceptance in Bin%d",iBin));
    
    // Draw FitResult
    f2_model.SetTitle("");
    f2_model.SetMaximum(accUpperBound);
    f2_model.SetLineWidth(1);
    //latex->DrawLatex(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
    //latex->DrawLatex(0.01,0.90,TString::Format("DoF = %d",nbinsK*nbinsL-gMinuit->GetNumFreePars()));
    //f2_model.Draw("SURF SAME ");
    canvas.SaveSource(TString::Format("%s/acceptance_2D_bin%d.cc",plotpath.Data(),iBin));
    canvas.Print(TString::Format("%s/acceptance_2D_bin%d.pdf",plotpath.Data(),iBin));

    //// Draw compare
    TH2F h2_compFit("h2_compFit","",nbinsL,-1,1,nbinsK,-1,1);
    for (int i = 1; i <= nbinsL; i++) {//thetaL
        for (int j = 1; j <= nbinsK; j++) {//thetaK
            if (h2_acc.GetBinContent(i,j) != 0){
                h2_compFit.SetBinContent(i,j,f2_model.Eval(h2_acc.GetXaxis()->GetBinCenter(i),h2_acc.GetYaxis()->GetBinCenter(j))/h2_acc.GetBinContent(i,j));
            }else{
                h2_compFit.SetBinContent(i,j,0.);
            }
        }
    }
    h2_compFit.SetMinimum(0.);
    h2_compFit.SetStats(0);
    h2_compFit.SetTitleOffset(2,"XY");
    h2_compFit.SetXTitle("genCosThetaL");
    h2_compFit.SetYTitle("genCosThetaK");
    h2_compFit.Draw("LEGO2");
    latex->DrawLatex(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
    latex->DrawLatex(0.01,0.90,TString::Format("DoF = %d",nbinsK*nbinsL-gMinuit->GetNumFreePars()));
    latex->DrawLatex(0.30,0.95,TString::Format("acceptance_{measured} / acceptance_{fit} in Bin%d",iBin));
    canvas.Update();
    //canvas.SaveSource(TString::Format("%s/acceptance_compFit_2D_bin%d.cc",plotpath.Data(),iBin));
    //canvas.Print(TString::Format("%s/acceptance_compFit_2D_bin%d.pdf",plotpath.Data(),iBin));
    
    // Draw significance
    TH1F h_pull("Deviation/Error","",15,-3.,3.);
    h_pull.SetXTitle("Significance of deviation");
    h_pull.SetYTitle("Angular bins");
    for (int i = 1; i <= nbinsL; i++) {//thetaL
        for (int j = 1; j <= nbinsK; j++) {//thetaK
            double _xlo = h2_acc.GetXaxis()->GetBinLowEdge(i);
            double _xhi = h2_acc.GetXaxis()->GetBinUpEdge(i);
            double _ylo = h2_acc.GetYaxis()->GetBinLowEdge(j);
            double _yhi = h2_acc.GetYaxis()->GetBinUpEdge(j);
            if (h2_nacc.GetBinContent(i,j) != 0){
                h_pull.Fill((f2_model.Integral(_xlo,_xhi,_ylo,_yhi)/(_xhi-_xlo)/(_yhi-_ylo)-h2_acc.GetBinContent(i,j))/h2_acc.GetBinError(i,j));
            }
        }
    }
    h_pull.Draw();
    canvas.Update();
    //canvas.SaveSource(TString::Format("%s/acceptance_sigma_bin%d.cc",plotpath.Data(),iBin));
    //canvas.Print(TString::Format("%s/acceptance_sigma_bin%d.pdf",plotpath.Data(),iBin));
    
    // Clear
    delete latex;
    delete gMinuit;

    //prepare output
    std::vector<double> output;
    for (int iPar = 0; iPar < nPar; iPar++){
        output.push_back(arrPar[iPar]);
        output.push_back(arrParErr[iPar]);
        
        printf("%f,",arrPar[iPar]);
        if (iPar+1 >= nPar) printf("\n");
    }
    for (int i = 0; i < output.size(); i=i+2) {
        printf("%f,",output[i+1]);
        if (i+2 >= output.size()) printf("\n");
    }
    return output;
}//}}}

std::vector<double> recoEff(int iBin) // reconstruction efficiency. Not used in standard procedure, just for check...
{//{{{
    printf("Evaluate reconstruction efficiency for bin#%d\n",iBin);
    double effUpperBound = 0.5;
    double BMass = 0;
    double Mumumass = 0;
    double Mumumasserr = 0;
    double gQ2 = 0;
    double gCosThetaK = 0;
    double gCosThetaL = 0;
    double gmuppt = 0;
    double gmupeta= 0;
    double gmumpt = 0;
    double gmumeta= 0;

    ch->SetBranchStatus("*",0);
    ch->SetBranchStatus("Bmass"         , 1);
    ch->SetBranchStatus("Mumumass"      , 1);
    ch->SetBranchStatus("Mumumasserr"   , 1);
    ch->SetBranchStatus("genQ2"         , 1);
    ch->SetBranchStatus("genCosTheta*"  , 1);
    ch->SetBranchStatus("genMu*"        , 1);
    ch->SetBranchAddress("Bmass"        , &BMass);
    ch->SetBranchAddress("Mumumass"     , &Mumumass);
    ch->SetBranchAddress("Mumumasserr"  , &Mumumasserr);
    ch->SetBranchAddress("genQ2"        , &gQ2);
    ch->SetBranchAddress("genCosThetaK" , &gCosThetaK);
    ch->SetBranchAddress("genCosThetaL" , &gCosThetaL);
    ch->SetBranchAddress("genMupPt"     , &gmuppt);
    ch->SetBranchAddress("genMupEta"    , &gmupeta);
    ch->SetBranchAddress("genMumPt"     , &gmumpt);
    ch->SetBranchAddress("genMumEta"    , &gmumeta);

    // Fill histograms
    float thetaKBins[6]={-1,-0.7,0.,0.4,0.8,1};
    float thetaLBins[7]={-1,-0.7,-0.3,0.,0.3,0.7,1};
    TH2F h2_nacc("h2_nacc" ,"h2_nacc" ,6,thetaLBins,5,thetaKBins); 
    TH2F h2_nreco("h2_nreco","h2_nreco",6,thetaLBins,5,thetaKBins);
    for (int entry = 0; entry < ch->GetEntries(); entry++) {
        ch->GetEntry(entry);
        if (gQ2 > q2rangeup[iBin] || gQ2 < q2rangedn[iBin]) continue;
        if (gQ2 > q2rangeup[3] && gQ2 < q2rangedn[3]) continue;//jpsi
        if (gQ2 > q2rangeup[5] && gQ2 < q2rangedn[5]) continue;//psi2s
        if ( fabs(gmumeta) < 2.3 && fabs(gmupeta) < 2.3 && gmumpt > 2.8 && gmuppt > 2.8 ) h2_nacc.Fill(gCosThetaL,gCosThetaK);
        if (BMass != 0 && ((Mumumass > 3.096916+3.5*Mumumasserr || Mumumass < 3.096916-5.5*Mumumasserr) && (Mumumass > 3.686109+3.5*Mumumasserr || Mumumass < 3.686109-     3.5*Mumumasserr)) ){
            h2_nreco.Fill(gCosThetaL,gCosThetaK);
        }
    }
    
    // Calculate efficiency
    TH2F h2_rec("h2_rec","",6,thetaLBins,5,thetaKBins);
    h2_rec.SetMinimum(0.);
    h2_rec.SetTitleOffset(2,"XY");
    h2_rec.SetXTitle("genCosThetaL");
    h2_rec.SetYTitle("genCosThetaK");
    for (int i = 1; i <= 6; i++) {
        for (int j = 1; j <= 5; j++) {
            // Build from MC samples
            if (h2_nacc.GetBinContent(i,j) == 0 || h2_nreco.GetBinContent(i,j) == 0) {
                printf("WARNING: Efficiency(%d,%d)=0, set error to be 1.\n",i,j);
                h2_rec.SetBinContent(i,j,0.);
                h2_rec.SetBinError(i,j,1.);
            }else{
                h2_rec.SetBinContent(i,j,h2_nreco.GetBinContent(i,j)/h2_nacc.GetBinContent(i,j));
                h2_rec.SetBinError(i,j,sqrt(h2_rec.GetBinContent(i,j)*(1-h2_rec.GetBinContent(i,j))/h2_nacc.GetBinContent(i,j)));
            }

            // Build a toy sample
            //h2_rec.SetBinContent(i,j,fabs(sin((thetaLBins[i]+thetaLBins[i-1])*3.1415926/4)/(20.-2*j))+0.1);
            //h2_rec.SetBinError(i,j,sqrt(h2_rec.GetBinContent(i,j)*(1-h2_rec.GetBinContent(i,j))/20));
        }
    }

    // Using pure TMinuit
    int nPar = 20;
    TMinuit *gMinuit = new TMinuit(nPar);
    h2_fcn = &h2_rec;
    gMinuit->SetFCN(fcn_binnedChi2_2D);
    
    // Use Legendre polynomial for better convergance
    // 1,x,(3x^2-1)/2,(5x^3-3x)/2
    //TF2 f2_model("f2_model","([0]+[1]*y+[2]*y**2+[3]*y**3)+([4]+[5]*y+[6]*y**2+[7]*y**3)*x**2+([8]+[9]*y+[10]*y**2+[11]*y**3)*x**3+([12]+[13]*y+[14]*y**2+[15]*y**3)*x**4+([16]+[17]*y+[18]*y**2+[19]*y**3)*x**6",-1.,1.,-1.,1.);
    TF2 f2_model("f2_model","([0]+[1]*y+[2]*(3*y**2-1)/2+[3]*(5*y**3-3*y)/2)+([4]+[5]*y+[6]*(3*y**2-1)/2+[7]*(5*y**3-3*y)/2)*x**2+([8]+[9]*y+[10]*(3*y**2-1)/2+[11]*(5*y**3-3*y)/2)*x**3+([12]+[13]*y+[14]*(3*y**2-1)/2+[15]*(5*y**3-3*y)/2)*x**4+([16]+[17]*y+[18]*(3*y**2-1)/2+[19]*(5*y**3-3*y)/2)*x**6",-1.,1.,-1.,1.);
    f2_fcn = &f2_model;
    gMinuit->DefineParameter( 0, "k0l0",  .01,  1E-3,    -1E+1, 1E+1);
    gMinuit->DefineParameter( 1, "k1l0",   0.,  1E-3,    -1E+1, 1E+1);
    gMinuit->DefineParameter( 2, "k2l0",   0.,  1E-3,    -1E+1, 1E+1);
    gMinuit->DefineParameter( 3, "k3l0",   0.,  1E-3,    -1E+1, 1E+1);
    gMinuit->DefineParameter( 4, "k0l2", 1E-2,  1E-3,    -1E+1, 1E+1);
    gMinuit->DefineParameter( 5, "k1l2",   0.,  1E-3,    -1E+1, 1E+1);
    gMinuit->DefineParameter( 6, "k2l2",   0.,  1E-3,    -1E+1, 1E+1);
    gMinuit->DefineParameter( 7, "k3l2",   0.,  1E-3,    -1E+1, 1E+1);
    gMinuit->DefineParameter( 8, "k0l3",   0.,  1E-3,    -1E+1, 1E+1);
    gMinuit->DefineParameter( 9, "k1l3",   0.,  1E-3,    -1E+1, 1E+1);
    gMinuit->DefineParameter(10, "k2l3",   0.,  1E-3,    -1E+1, 1E+1);
    gMinuit->DefineParameter(11, "k3l3",   0.,  1E-3,    -1E+1, 1E+1);
    gMinuit->DefineParameter(12, "k0l4",-1E-2,  1E-3,    -1E+1, 1E+1);
    gMinuit->DefineParameter(13, "k1l4",   0.,  1E-3,    -1E+1, 1E+1);
    gMinuit->DefineParameter(14, "k2l4",   0.,  1E-3,    -1E+1, 1E+1);
    gMinuit->DefineParameter(15, "k3l4",   0.,  1E-3,    -1E+1, 1E+1);
    gMinuit->DefineParameter(16, "k0l6", 1E-2,  1E-3,    -1E+2, 1E+3);
    gMinuit->DefineParameter(17, "k1l6",   0.,  1E-3,    -1E+2, 1E+3);
    gMinuit->DefineParameter(18, "k2l6",   0.,  1E-3,    -1E+2, 1E+3);
    gMinuit->DefineParameter(19, "k3l6",   0.,  1E-3,    -1E+2, 1E+3);

    if (iBin == 0) {
        gMinuit->Command("SET PARM 9 0");
        gMinuit->Command("SET PARM 10 0");
        gMinuit->Command("SET PARM 11 0");
        gMinuit->Command("SET PARM 12 0");
        gMinuit->Command("FIX 9");
        gMinuit->Command("FIX 10");
        gMinuit->Command("FIX 11");
        gMinuit->Command("FIX 12");
    }else if (iBin == 1) {
        gMinuit->Command("SET PARM 9 0");
        gMinuit->Command("SET PARM 10 0");
        gMinuit->Command("SET PARM 11 0");
        gMinuit->Command("SET PARM 12 0");
        gMinuit->Command("SET PARM 17 0");
        gMinuit->Command("SET PARM 18 0");
        gMinuit->Command("SET PARM 19 0");
        gMinuit->Command("SET PARM 20 0");
        gMinuit->Command("FIX 9");
        gMinuit->Command("FIX 10");
        gMinuit->Command("FIX 11");
        gMinuit->Command("FIX 12");
        gMinuit->Command("FIX 17");
        gMinuit->Command("FIX 18");
        gMinuit->Command("FIX 19");
        gMinuit->Command("FIX 20");
    }else if (iBin > 1 && iBin < 6 ) {
        gMinuit->Command("SET PARM 17 0");
        gMinuit->Command("SET PARM 18 0");
        gMinuit->Command("SET PARM 19 0");
        gMinuit->Command("SET PARM 20 0");
        gMinuit->Command("FIX 17");
        gMinuit->Command("FIX 18");
        gMinuit->Command("FIX 19");
        gMinuit->Command("FIX 20");
    }else{
        gMinuit->Command("SET PARM 13 0");
        gMinuit->Command("SET PARM 14 0");
        gMinuit->Command("SET PARM 15 0");
        gMinuit->Command("SET PARM 16 0");
        gMinuit->Command("SET PARM 17 0");
        gMinuit->Command("SET PARM 18 0");
        gMinuit->Command("SET PARM 19 0");
        gMinuit->Command("SET PARM 20 0");
        gMinuit->Command("FIX 13");
        gMinuit->Command("FIX 14");
        gMinuit->Command("FIX 15");
        gMinuit->Command("FIX 16");
        gMinuit->Command("FIX 17");
        gMinuit->Command("FIX 18");
        gMinuit->Command("FIX 19");
        gMinuit->Command("FIX 20");
    }
    
    gMinuit->Command("MINI");
    gMinuit->Command("MINI");
    gMinuit->Command("IMPROVE");
    gMinuit->Command("MINOS");

    double arrPar[nPar];
    double arrParErr[nPar];
    for (int iPar = 0; iPar < nPar; iPar++) gMinuit->GetParameter(iPar,arrPar[iPar],arrParErr[iPar]);
    for (int iPar = 0; iPar < nPar; iPar++) f2_model.SetParameter(iPar,arrPar[iPar]);
    
    // Prepare draw
    TCanvas canvas("canvas");
    TLatex *latex = new TLatex();
    latex->SetNDC();
    
    // Draw efficiency
    h2_rec.SetStats(0);
    h2_rec.SetMaximum(effUpperBound);
    h2_rec.Draw("LEGO2");
    latex->DrawLatex(0.35,0.95,TString::Format("#varepsilon_{RECO} in Bin%d",iBin));
    
    // Draw FitResult
    f2_model.SetTitle("");
    f2_model.SetMaximum(effUpperBound);
    f2_model.SetLineWidth(1);
    f2_model.Draw("SURF SAME ");
    canvas.SaveSource(TString::Format("%s/recoEff_2D_bin%d.cc",plotpath.Data(),iBin));
    canvas.Print(TString::Format("%s/recoEff_2D_bin%d.pdf",plotpath.Data(),iBin));

    //// Draw compare
    double chi2Val=0;
    fcn_binnedChi2_2D(nPar, 0, chi2Val, arrPar, 0);
    printf("Chi2(Bin center)=%f \n",chi2Val);
    
    TH2F h2_compFit("h2_compFit","",6,thetaLBins,5,thetaKBins);
    h2_compFit.SetTitleOffset(2,"XY");
    h2_compFit.SetXTitle("genCosThetaL");
    h2_compFit.SetYTitle("genCosThetaK");
    for (int i = 1; i <= 6; i++) {//thetaL
        for (int j = 1; j <= 5; j++) {//thetaK
            if (h2_rec.GetBinContent(i,j) != 0){
                h2_compFit.SetBinContent(i,j,f2_model.Eval(h2_rec.GetXaxis()->GetBinCenter(i),h2_rec.GetYaxis()->GetBinCenter(j))/h2_rec.GetBinContent(i,j));
            }else{
                h2_compFit.SetBinContent(i,j,0.);
            }
        }
    }
    h2_compFit.SetMinimum(0.);
    h2_compFit.SetStats(0);
    h2_compFit.Draw("LEGO2");
    latex->DrawLatex(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
    latex->DrawLatex(0.3,0.95,TString::Format("#varepsilon_{RECO,fit} / #varepsilon_{RECO,measured} in Bin%d",iBin));
    canvas.Update();
    canvas.SaveSource(TString::Format("%s/recoEff_compFit_2D_bin%d.cc",plotpath.Data(),iBin));
    canvas.Print(TString::Format("%s/recoEff_compFit_2D_bin%d.pdf",plotpath.Data(),iBin));

    // Draw significance of deviation
    TH1F h_pull("Deviation/Error","",15,-3.,3.);
    h_pull.SetXTitle("Significance of deviation");
    h_pull.SetYTitle("Angular bins");
    for (int i = 1; i <= 6; i++) {//thetaL
        for (int j = 1; j <= 5; j++) {//thetaK
            double _xlo = h2_rec.GetXaxis()->GetBinLowEdge(i);
            double _xhi = h2_rec.GetXaxis()->GetBinUpEdge(i);
            double _ylo = h2_rec.GetYaxis()->GetBinLowEdge(j);
            double _yhi = h2_rec.GetYaxis()->GetBinUpEdge(j);
            if (h2_rec.GetBinContent(i,j) != 0){
                h_pull.Fill((f2_model.Integral(_xlo,_xhi,_ylo,_yhi)/(_xhi-_xlo)/(_yhi-_ylo)-h2_rec.GetBinContent(i,j))/h2_rec.GetBinError(i,j));
            }
        }
    }
    h_pull.Draw();
    canvas.Update();
    canvas.SaveSource(TString::Format("%s/recoEff_sigma_bin%d.cc",plotpath.Data(),iBin));
    canvas.Print(TString::Format("%s/recoEff_sigma_bin%d.pdf",plotpath.Data(),iBin));

    // Clear
    delete latex;
    delete gMinuit;

    //prepare output
    std::vector<double> output;
    for (int iPar = 0; iPar < nPar; iPar++){
        output.push_back(arrPar[iPar]);
        output.push_back(arrParErr[iPar]);
        
        printf("%f,",arrPar[iPar]);
        if (iPar+1 >= nPar) printf("\n");
    }
    for (int i = 0; i < output.size(); i=i+2) {
        printf("%f,",output[i+1]);
        if (i+2 >= output.size()) printf("\n");
    }
    return output;
}//}}}

void createAcceptanceHist() // create acceptance histogram from UNFILTERED GEN.
{//{{{
    double accUpperBound = 0.09;
    double gQ2 = 0;
    double gCosThetaK = 0;
    double gCosThetaL = 0;
    double gmuppt = 0;
    double gmupeta= 0;
    double gmupphi= 0;
    double gmumpt = 0;
    double gmumeta= 0;
    double gmumphi= 0;

    TChain *treein=new TChain("tree");
    treein->Add("./sel_v3p5/unfilt-mc-genonly/sel_BToKstarMuMu_NoFilt_8TeV_mc_genonly_s*.root"); // fixed input
    if (treein == NULL) gSystem->Exit(0);
    treein->SetBranchStatus("*",0);
    treein->SetBranchStatus("genQ2"         , 1);
    treein->SetBranchStatus("genCosTheta*"  , 1);
    treein->SetBranchStatus("genMu*"        , 1);
    treein->SetBranchAddress("genQ2"        , &gQ2);
    treein->SetBranchAddress("genCosThetaK" , &gCosThetaK);
    treein->SetBranchAddress("genCosThetaL" , &gCosThetaL);
    treein->SetBranchAddress("genMupPt"     , &gmuppt);
    treein->SetBranchAddress("genMupEta"    , &gmupeta);
    treein->SetBranchAddress("genMupPhi"    , &gmupphi);
    treein->SetBranchAddress("genMumPt"     , &gmumpt);
    treein->SetBranchAddress("genMumEta"    , &gmumeta);
    treein->SetBranchAddress("genMumPhi"    , &gmumphi);

    // Create histograms
    TFile *fout = new TFile("acceptance_8TeV.root","RECREATE");
    float thetaKBins[6]={-1,-0.7,0.,0.4,0.8,1};
    float thetaLBins[7]={-1,-0.7,-0.3,0.,0.3,0.7,1};
    TH2F *h2_ngen[nQ2Ranges];
    TH2F *h2_nacc[nQ2Ranges];
    TH2F *h2_acc[nQ2Ranges];
    TH2F *h2_ngen_fine[nQ2Ranges];
    TH2F *h2_nacc_fine[nQ2Ranges];
    TH2F *h2_acc_fine[nQ2Ranges];
    TH1F *h_ngenL_fine[nQ2Ranges];
    TH1F *h_naccL_fine[nQ2Ranges];
    TH1F *h_accL_fine[nQ2Ranges];
    TH1F *h_ngenK_fine[nQ2Ranges];
    TH1F *h_naccK_fine[nQ2Ranges];
    TH1F *h_accK_fine[nQ2Ranges];
    for(int iBin = 0; iBin < nQ2Ranges; iBin++){
        h2_ngen[iBin] = new TH2F(TString::Format("h2_ngen_bin%d",iBin),"h2_ngen",6,thetaLBins,5,thetaKBins);
        h2_nacc[iBin] = new TH2F(TString::Format("h2_nacc_bin%d",iBin) ,"h2_nacc" ,6,thetaLBins,5,thetaKBins); 
        h2_acc [iBin] = new TH2F(TString::Format("h2_acc_bin%d",iBin),"",6,thetaLBins,5,thetaKBins);
        h2_ngen_fine[iBin] = new TH2F(TString::Format("h2_ngen_fine_bin%d",iBin),"h2_ngen",20,-1,1,20,-1,1);
        h2_nacc_fine[iBin] = new TH2F(TString::Format("h2_nacc_fine_bin%d",iBin) ,"h2_nacc" ,20,-1,1,20,-1,1); 
        h2_acc_fine[iBin]  = new TH2F(TString::Format("h2_acc_fine_bin%d",iBin),"",20,-1,1,20,-1,1);
        h_ngenL_fine[iBin] = new TH1F(TString::Format("h_ngenL_fine_bin%d",iBin),"h_ngenL",20,-1,1);
        h_naccL_fine[iBin] = new TH1F(TString::Format("h_naccL_fine_bin%d",iBin) ,"h_naccL" ,20,-1,1); 
        h_accL_fine[iBin]  = new TH1F(TString::Format("h_accL_fine_bin%d",iBin),"",20,-1,1);
        h_ngenK_fine[iBin] = new TH1F(TString::Format("h_ngenK_fine_bin%d",iBin),"h_ngenK",20,-1,1);
        h_naccK_fine[iBin] = new TH1F(TString::Format("h_naccK_fine_bin%d",iBin) ,"h_naccK" ,20,-1,1); 
        h_accK_fine[iBin]  = new TH1F(TString::Format("h_accK_fine_bin%d",iBin),"",20,-1,1);
        h2_ngen[iBin]->SetTitleOffset(2,"XYZ");
        h2_ngen[iBin]->SetXTitle("genCosThetaL");
        h2_ngen[iBin]->SetYTitle("genCosThetaK");
        h2_ngen[iBin]->SetZTitle("Generated events");
        h2_nacc[iBin]->SetTitleOffset(2,"XY");
        h2_nacc[iBin]->SetXTitle("genCosThetaL");
        h2_nacc[iBin]->SetYTitle("genCosThetaK");
        h2_nacc[iBin]->SetZTitle("Events in acceptance");
        h2_acc [iBin]->SetStats(0);
        h2_acc [iBin]->SetMinimum(0.);
        h2_acc [iBin]->SetMaximum(accUpperBound);
        h2_acc [iBin]->SetTitleOffset(2,"XY");
        h2_acc [iBin]->SetXTitle("genCosThetaL");
        h2_acc [iBin]->SetYTitle("genCosThetaK");
        h2_acc [iBin]->SetZTitle("Acceptance");
        h2_ngen_fine[iBin]->SetTitleOffset(2,"XYZ");
        h2_ngen_fine[iBin]->SetXTitle("genCosThetaL");
        h2_ngen_fine[iBin]->SetYTitle("genCosThetaK");
        h2_ngen_fine[iBin]->SetZTitle("Generated events");
        h2_nacc_fine[iBin]->SetTitleOffset(2,"XY");
        h2_nacc_fine[iBin]->SetXTitle("genCosThetaL");
        h2_nacc_fine[iBin]->SetYTitle("genCosThetaK");
        h2_nacc_fine[iBin]->SetZTitle("Events in acceptance");
        h2_acc_fine [iBin]->SetStats(0);
        h2_acc_fine [iBin]->SetMinimum(0.);
        h2_acc_fine [iBin]->SetMaximum(accUpperBound);
        h2_acc_fine [iBin]->SetTitleOffset(2,"XY");
        h2_acc_fine [iBin]->SetXTitle("genCosThetaL");
        h2_acc_fine [iBin]->SetYTitle("genCosThetaK");
        h2_acc_fine [iBin]->SetZTitle("Acceptance");
        h_ngenL_fine[iBin]->SetXTitle("genCosThetaL");
        h_ngenL_fine[iBin]->SetZTitle("Generated events");
        h_naccL_fine[iBin]->SetXTitle("genCosThetaL");
        h_naccL_fine[iBin]->SetZTitle("Events in acceptance");
        h_accL_fine [iBin]->SetStats(0);
        h_accL_fine [iBin]->SetMinimum(0.);
        h_accL_fine [iBin]->SetMaximum(accUpperBound);
        h_accL_fine [iBin]->SetXTitle("genCosThetaL");
        h_accL_fine [iBin]->SetZTitle("Acceptance");
        h_ngenK_fine[iBin]->SetXTitle("genCosThetaK");
        h_ngenK_fine[iBin]->SetZTitle("Generated events");
        h_naccK_fine[iBin]->SetXTitle("genCosThetaK");
        h_naccK_fine[iBin]->SetZTitle("Events in acceptance");
        h_accK_fine [iBin]->SetStats(0);
        h_accK_fine [iBin]->SetMinimum(0.);
        h_accK_fine [iBin]->SetMaximum(accUpperBound);
        h_accK_fine [iBin]->SetXTitle("genCosThetaK");
        h_accK_fine [iBin]->SetZTitle("Acceptance");
    }
    
    // Fill histograms
        // Read data
    for (int entry = 0; entry < treein->GetEntries(); entry++) {
        treein->GetEntry(entry);
        for(int iBin = 0; iBin < nQ2Ranges; iBin++){
            if (gQ2 > q2rangeup[iBin] || gQ2 < q2rangedn[iBin]) continue;
            if (gQ2 > q2rangeup[3] && gQ2 < q2rangedn[3]) continue;//jpsi
            if (gQ2 > q2rangeup[5] && gQ2 < q2rangedn[5]) continue;//psi2s
            h2_ngen[iBin]->Fill(gCosThetaL,gCosThetaK);
            h2_ngen_fine[iBin]->Fill(gCosThetaL,gCosThetaK);
            h_ngenL_fine[iBin]->Fill(gCosThetaL);
            h_ngenK_fine[iBin]->Fill(gCosThetaK);
            if ( fabs(gmupeta) < 2.3 && gmuppt > 2.8 && fabs(gmumeta) < 2.3 && gmumpt > 2.8){
                h2_nacc[iBin]->Fill(gCosThetaL,gCosThetaK);
                h2_nacc_fine[iBin]->Fill(gCosThetaL,gCosThetaK);
                h_naccL_fine[iBin]->Fill(gCosThetaL);
                h_naccK_fine[iBin]->Fill(gCosThetaK);
            }
        }
    }
        
    for(int iBin = 0; iBin < nQ2Ranges; iBin++){
        // Calculate acceptance
        h2_acc[iBin]->SetAxisRange(0.,1.,"Z");
        for (int i = 1; i <= 6; i++) {
            for (int j = 1; j <= 5; j++) {
                // Fill acceptance
                if (h2_ngen[iBin]->GetBinContent(i,j) == 0) {
                    printf("WARNING: Acceptance(%d,%d)=%f/%f\n",i,j,h2_nacc[iBin]->GetBinContent(i,j),h2_ngen[iBin]->GetBinContent(i,j));
                    h2_acc[iBin]->SetBinContent(i,j,0.);
                    h2_acc[iBin]->SetBinError(i,j,1.);
                }else{
                    h2_acc[iBin]->SetBinContent(i,j,h2_nacc[iBin]->GetBinContent(i,j)/h2_ngen[iBin]->GetBinContent(i,j));
                    if (h2_nacc[iBin]->GetBinContent(i,j) != 0){
                        h2_acc[iBin]->SetBinError(i,j,sqrt(h2_acc[iBin]->GetBinContent(i,j)*(1.-h2_acc[iBin]->GetBinContent(i,j))/h2_ngen[iBin]->GetBinContent(i,j)));
                    }else{
                        h2_acc[iBin]->SetBinError(i,j,0.);
                    }
                }
            }
        }
        printf("INFO: h2_acc_bin%d built.\n",iBin);
        
        h2_acc_fine[iBin]->SetAxisRange(0.,1.,"Z");
        for (int i = 1; i <= 20; i++) {//L
            for (int j = 1; j <= 20; j++) {//K
                // Fill acceptance
                if (h2_ngen_fine[iBin]->GetBinContent(i,j) == 0) {
                    h2_acc_fine[iBin]->SetBinContent(i,j,0.);
                    h2_acc_fine[iBin]->SetBinError(i,j,1.);
                }else{
                    h2_acc_fine[iBin]->SetBinContent(i,j,h2_nacc_fine[iBin]->GetBinContent(i,j)/h2_ngen_fine[iBin]->GetBinContent(i,j));
                    if (h2_nacc_fine[iBin]->GetBinContent(i,j) != 0){
                        h2_acc_fine[iBin]->SetBinError(i,j,sqrt(h2_acc_fine[iBin]->GetBinContent(i,j)*(1.-h2_acc_fine[iBin]->GetBinContent(i,j))/h2_ngen_fine[iBin]->GetBinContent(i,j)));
                    }else{
                        h2_acc_fine[iBin]->SetBinError(i,j,0.);
                    }
                }
                
            }
            
            // 1-D
            if (h_ngenL_fine[iBin]->GetBinContent(i) == 0) {
                h_accL_fine[iBin]->SetBinContent(i,0.);
                h_accL_fine[iBin]->SetBinError(i,1.);
            }else{
                h_accL_fine[iBin]->SetBinContent(i,h_naccL_fine[iBin]->GetBinContent(i)/h_ngenL_fine[iBin]->GetBinContent(i));
                if (h_naccL_fine[iBin]->GetBinContent(i) != 0){
                    h_accL_fine[iBin]->SetBinError(i,sqrt(h_accL_fine[iBin]->GetBinContent(i)*(1.-h_accL_fine[iBin]->GetBinContent(i))/h_ngenL_fine[iBin]->GetBinContent(i)));
                }else{
                    h_accL_fine[iBin]->SetBinError(i,0.);
                }
            }
            
        }
        for (int i = 1; i <= 20; i++) {//K
            // 1-D
            if (h_ngenK_fine[iBin]->GetBinContent(i) == 0) {
                h_accK_fine[iBin]->SetBinContent(i,0.);
                h_accK_fine[iBin]->SetBinError(i,1.);
            }else{
                h_accK_fine[iBin]->SetBinContent(i,h_naccK_fine[iBin]->GetBinContent(i)/h_ngenK_fine[iBin]->GetBinContent(i));
                if (h_naccK_fine[iBin]->GetBinContent(i) != 0){
                    h_accK_fine[iBin]->SetBinError(i,sqrt(h_accK_fine[iBin]->GetBinContent(i)*(1.-h_accK_fine[iBin]->GetBinContent(i))/h_ngenK_fine[iBin]->GetBinContent(i)));
                }else{
                    h_accK_fine[iBin]->SetBinError(i,0.);
                }
            }
        }
        printf("INFO: h2_acc_fine_bin%d built.\n",iBin);
    }
    fout->Write();
    fout->Close();
}//}}}

void createRecoEffHist(int iBin) // create reco efficiency histogram from official MC sample.
{//{{{
    printf("Evaluate reconstruction efficiency for bin#%d\n",iBin);
    double effUpperBound = 0.03;
    double BMass = 0;
    double Mumumass = 0;
    double Mumumasserr = 0;
    double gQ2 = 0;
    double gCosThetaK = 0;
    double gCosThetaL = 0;
    double gmuppt = 0;
    double gmupeta= 0;
    double gmumpt = 0;
    double gmumeta= 0;
    double Triggers=0;

    ch->SetBranchStatus("*",0);
    ch->SetBranchStatus("Bmass"         , 1);
    ch->SetBranchStatus("Mumumass"      , 1);
    ch->SetBranchStatus("Mumumasserr"   , 1);
    ch->SetBranchStatus("genQ2"         , 1);
    ch->SetBranchStatus("genCosTheta*"  , 1);
    ch->SetBranchStatus("genMu*"        , 1);
    ch->SetBranchAddress("Bmass"        , &BMass);
    ch->SetBranchAddress("Mumumass"     , &Mumumass);
    ch->SetBranchAddress("Mumumasserr"  , &Mumumasserr);
    ch->SetBranchAddress("genQ2"        , &gQ2);
    ch->SetBranchAddress("genCosThetaK" , &gCosThetaK);
    ch->SetBranchAddress("genCosThetaL" , &gCosThetaL);
    ch->SetBranchAddress("genMupPt"     , &gmuppt);
    ch->SetBranchAddress("genMupEta"    , &gmupeta);
    ch->SetBranchAddress("genMumPt"     , &gmumpt);
    ch->SetBranchAddress("genMumEta"    , &gmumeta);
    ch->SetBranchAddress("Triggers"     , &Triggers);

    // Fill histograms
    const int nKBins = 20;
    const int nLBins = 20;
    //float thetaKBins[nKBins]={-1,-0.7,0.,0.4,0.8,1};//nKBins=5
    //float thetaLBins[nLBins]={-1,-0.7,-0.3,0.,0.3,0.7,1};//nLBins=6
    //TH2F h2_nacc("h2_nacc" ,"h2_nacc" ,6,thetaLBins,5,thetaKBins); 
    //TH2F h2_nreco("h2_nreco","h2_nreco",6,thetaLBins,5,thetaKBins);
    TH2F h2_nacc("h2_nacc","h2_nacc",nLBins,-1,1,nKBins,-1,1);
    TH2F h2_nreco("h2_nreco","h2_nreco",nLBins,-1,1,nKBins,-1,1);
    for (int entry = 0; entry < ch->GetEntries(); entry++) {
        ch->GetEntry(entry);
        if (gQ2 > q2rangeup[iBin] || gQ2 < q2rangedn[iBin]) continue;
        if (Triggers == 0) continue;
        if ( fabs(gmumeta) < 2.3 && fabs(gmupeta) < 2.3 && gmumpt > 2.8 && gmuppt > 2.8 ) h2_nacc.Fill(gCosThetaL,gCosThetaK);
        if (BMass != 0 && ((Mumumass > 3.096916+3.5*Mumumasserr || Mumumass < 3.096916-5.5*Mumumasserr) && (Mumumass > 3.686109+3.5*Mumumasserr || Mumumass < 3.686109-     3.5*Mumumasserr)) ){
            h2_nreco.Fill(gCosThetaL,gCosThetaK);
        }
    }
    
    // Calculate efficiency
    //TH2F h2_rec("h2_rec","",6,thetaLBins,5,thetaKBins);
    TH2F h2_rec("h2_rec","",nLBins,-1,1,nKBins,-1,1);
    for (int i = 1; i <= nLBins; i++) {
        for (int j = 1; j <= nKBins; j++) {
            // Build from MC samples
            if (h2_nacc.GetBinContent(i,j) == 0 || h2_nreco.GetBinContent(i,j) == 0) {
                printf("WARNING: Efficiency(%d,%d)=0, set error to be 1.\n",i,j);
                h2_rec.SetBinContent(i,j,0.);
                h2_rec.SetBinError(i,j,1.);
            }else{
                h2_rec.SetBinContent(i,j,h2_nreco.GetBinContent(i,j)/h2_nacc.GetBinContent(i,j));
                h2_rec.SetBinError(i,j,sqrt(h2_rec.GetBinContent(i,j)*(1-h2_rec.GetBinContent(i,j))/h2_nreco.GetBinContent(i,j)));
                printf("INFO: Efficiency(%d,%d)=%f +- %f.\n",i,j,h2_rec.GetBinContent(i,j),h2_rec.GetBinError(i,j));
            }
        }
    }
    h2_rec.SetTitleOffset(2,"XY");
    h2_rec.SetXTitle("CosThetaL");
    h2_rec.SetYTitle("CosThetaK");
    h2_rec.SetStats(0);
    h2_rec.SetMinimum(0.);
    
    //TH1F h_recL("h_recL","",6,thetaLBins,5,thetaKBins);
    TH1F h_recL("h_recL","",nLBins,-1,1);
    for (int i = 1; i <= nLBins; i++) {
        double nacc = 0;
        double nreco = 0;
        for (int j = 1; j <= nKBins; j++) {
            nacc+= h2_nacc.GetBinContent(i,j);
            nreco+= h2_nreco.GetBinContent(i,j);
        }
        if (nacc !=0 ){
            h_recL.SetBinContent(i,nreco/nacc);
            h_recL.SetBinError(i,sqrt(h_recL.GetBinContent(i)*(1-h_recL.GetBinContent(i))/nacc));
        }else{
            h_recL.SetBinContent(i,0);
            h_recL.SetBinError(i,1);
        }
    }
    h_recL.SetStats(0);
    h_recL.SetMinimum(0.);
    h_recL.SetXTitle("CosThetaL");
    
    //TH1F h_recK("h_recK","",6,thetaLBins,5,thetaKBins);
    TH1F h_recK("h_recK","",nKBins,-1,1);
    for (int i = 1; i <= nKBins; i++) {
        double nacc = 0;
        double nreco = 0;
        for (int j = 1; j <= nLBins; j++) {
            nacc+= h2_nacc.GetBinContent(j,i);
            nreco+= h2_nreco.GetBinContent(j,i);
        }
        if (nacc !=0 ){
            h_recK.SetBinContent(i,nreco/nacc);
            h_recK.SetBinError(i,sqrt(h_recK.GetBinContent(i)*(1-h_recK.GetBinContent(i))/nacc));
        }else{
            h_recK.SetBinContent(i,0);
            h_recK.SetBinError(i,1);
        }
    }
    h_recK.SetStats(0);
    h_recK.SetMinimum(0.);
    h_recK.SetXTitle("CosThetaK");

    // Print
    TCanvas canvas("canvas");
    h2_rec.Draw("LEGO2 TEXT");
    canvas.SaveSource(TString::Format("%s/recoEff_2D_bin%d.cc",plotpath.Data(),iBin));
    canvas.Print(TString::Format("%s/recoEff_2D_bin%d.pdf",plotpath.Data(),iBin));
    h_recL.Draw("TEXT");
    canvas.Update();
    canvas.SaveSource(TString::Format("%s/recoEff_cosl_bin%d.cc",plotpath.Data(),iBin));
    canvas.Print(TString::Format("%s/recoEff_cosl_bin%d.pdf",plotpath.Data(),iBin));
    h_recK.Draw("TEXT");
    canvas.Update();
    canvas.SaveSource(TString::Format("%s/recoEff_cosk_bin%d.cc",plotpath.Data(),iBin));
    canvas.Print(TString::Format("%s/recoEff_cosk_bin%d.pdf",plotpath.Data(),iBin));
}//}}}

std::string accXrecoEff2(int iBin, bool keepParam = false) // acceptance*reconstruction efficiency
{//{{{
    switchRedirectStdio(TString::Format("%s/accXrecoEff2_stdout_bin%d.log",odatacardpath.Data(),iBin).Data(),"w",stdout);
    switchRedirectStdio(TString::Format("%s/accXrecoEff2_stderr_bin%d.log",odatacardpath.Data(),iBin).Data(),"w",stderr);

    printf("Evaluate full efficiency for bin#%d\n",iBin);
        // cut3-8TeV
    double effUpperBound = 0.00015;
    if (iBin == 3 || iBin == 5) effUpperBound = 2e-5;
    double BMass = 0;
    double Mumumass = 0;
    double Mumumasserr = 0;
    double Kstarmass = 0;
    double gQ2 = 0;
    double gCosThetaK = 0;
    double gCosThetaL = 0;
    double CosThetaK = 0;
    double CosThetaL = 0;
    double gmuppt = 0;
    double gmupeta= 0;
    double gmupphi= 0;
    double gmumpt = 0;
    double gmumeta= 0;
    double gmumphi= 0;
    int    triggers=0;

    ch->SetBranchStatus("*",0);
    ch->SetBranchStatus("Bmass"         , 1);
    ch->SetBranchStatus("Mumumass"      , 1);
    ch->SetBranchStatus("Mumumasserr"   , 1);
    ch->SetBranchStatus("genQ2"         , 1);
    ch->SetBranchStatus("genCosTheta*"  , 1);
    ch->SetBranchStatus("genMu*"        , 1);
    ch->SetBranchStatus("Triggers"      , 1);
    ch->SetBranchAddress("Bmass"        , &BMass);
    ch->SetBranchAddress("Mumumass"     , &Mumumass);
    ch->SetBranchAddress("Mumumasserr"  , &Mumumasserr);
    ch->SetBranchAddress("Kstarmass"    , &Kstarmass);
    ch->SetBranchAddress("genQ2"        , &gQ2);
    ch->SetBranchAddress("genCosThetaK" , &gCosThetaK);
    ch->SetBranchAddress("genCosThetaL" , &gCosThetaL);
    ch->SetBranchAddress("CosThetaK"    , &CosThetaK);
    ch->SetBranchAddress("CosThetaL"    , &CosThetaL);
    ch->SetBranchAddress("genMupPt"     , &gmuppt);
    ch->SetBranchAddress("genMupEta"    , &gmupeta);
    ch->SetBranchAddress("genMupPhi"    , &gmupphi);
    ch->SetBranchAddress("genMumPt"     , &gmumpt);
    ch->SetBranchAddress("genMumEta"    , &gmumeta);
    ch->SetBranchAddress("genMumPhi"    , &gmumphi);
    ch->SetBranchAddress("Triggers"     , &triggers);

    // Load acceptance
    TFile f_acc("acceptance_8TeV.root");
    TH2F *h2_acc = (TH2F*)f_acc.Get(TString::Format("h2_acc_bin%d",iBin));
    TH1F *h_accL = (TH1F*)f_acc.Get(TString::Format("h_accL_fine_bin%d",iBin));
    TH1F *h_accK = (TH1F*)f_acc.Get(TString::Format("h_accK_fine_bin%d",iBin));
    TH2F *h2_ngen = (TH2F*)f_acc.Get(TString::Format("h2_ngen_bin%d",iBin));
    TH1F *h_ngenL = (TH1F*)f_acc.Get(TString::Format("h_ngenL_fine_bin%d",iBin));
    TH1F *h_ngenK = (TH1F*)f_acc.Get(TString::Format("h_ngenK_fine_bin%d",iBin));

    // Fill histograms
    float thetaKBins[6]={-1,-0.7,0.,0.4,0.8,1};
    float thetaLBins[7]={-1,-0.7,-0.3,0.,0.3,0.7,1};
    TH2F h2_nacc("h2_nacc" ,"h2_nacc" ,6,thetaLBins,5,thetaKBins); 
    TH2F h2_nreco("h2_nreco","h2_nreco",6,thetaLBins,5,thetaKBins);
    int nLBins = 20;// The same value as h_acc
    int nKBins = 20;
    TH1F h_naccL("h_naccL" ,"h_naccL" ,nLBins,-1,1); 
    TH1F h_nrecoL("h_nrecoL","h_nrecoL",nLBins,-1,1);
    TH1F h_naccK("h_naccK" ,"h_naccK" ,nKBins,-1,1); 
    TH1F h_nrecoK("h_nrecoK","h_nrecoK",nKBins,-1,1);
    h_naccL.SetStats(0);
    h_naccL.SetMinimum(0.);
    h_naccL.SetXTitle("CosThetaL");
    h_naccL.SetYTitle("#Events/0.2");
    h_nrecoL.SetStats(0);
    h_nrecoL.SetMinimum(0.);
    h_nrecoL.SetXTitle("CosThetaL");
    h_nrecoL.SetYTitle("#Events/0.2");
    h_naccK.SetStats(0);
    h_naccK.SetMinimum(0.);
    h_naccK.SetXTitle("CosThetaK");
    h_naccK.SetYTitle("#Events/0.2");
    h_nrecoK.SetStats(0);
    h_nrecoK.SetMinimum(0.);
    h_nrecoK.SetXTitle("CosThetaK");
    h_nrecoK.SetYTitle("#Events/0.2");
    for (int entry = 0; entry < ch->GetEntries(); entry++) {
        ch->GetEntry(entry);
        if (gQ2 > q2rangeup[3] && gQ2 < q2rangedn[3]) continue;//jpsi
        if (gQ2 > q2rangeup[5] && gQ2 < q2rangedn[5]) continue;//psi2s
        if (gQ2 > q2rangeup[iBin] || gQ2 < q2rangedn[iBin]) continue;
        if ( fabs(gmupeta) < 2.3 && gmuppt > 2.8 && fabs(gmumeta) < 2.3 && gmumpt > 2.8){
            h2_nacc.Fill(gCosThetaL,gCosThetaK);
            h_naccL.Fill(gCosThetaL);
            h_naccK.Fill(gCosThetaK);
            if ( triggers > 0 && BMass != 0 && (Kstarmass > 0.792 && Kstarmass < 0.992)  && ((Mumumass > 3.096916+3.5*Mumumasserr || Mumumass < 3.096916-5.5*Mumumasserr) && (Mumumass > 3.686109+3.5*Mumumasserr || Mumumass < 3.686109-3.5*Mumumasserr)) ){
                // isCDFcut: 1 for CDF . 2 for LHCb . 3 for 16Aug reOptimization . 4 for sel_v3p5
                if (isCDFcut == 0){
                    h2_nreco.Fill(CosThetaL,CosThetaK);
                    h_nrecoL.Fill(CosThetaL);
                    h_nrecoK.Fill(CosThetaK);
                }else if (isCDFcut==1){
                    if( (iBin < 3 && fabs(BMass-Mumumass-2.182)>0.16) || (iBin ==4 && fabs(BMass-Mumumass-1.593)>0.06) || (iBin > 5 && fabs(BMass-Mumumass-1.593)>0.06) ){
                        h2_nreco.Fill(CosThetaL,CosThetaK);
                        h_nrecoL.Fill(CosThetaL);
                        h_nrecoK.Fill(CosThetaK);
                    }
                }else if (isCDFcut==2 || isCDFcut==3){
                    h2_nreco.Fill(CosThetaL,CosThetaK);
                    h_nrecoL.Fill(CosThetaL);
                    h_nrecoK.Fill(CosThetaK);
                }else if (isCDFcut==4){
                    if( fabs(BMass-Mumumass-2.182)>0.09 && fabs(BMass-Mumumass-1.593)>0.03 ){
                        h2_nreco.Fill(CosThetaL,CosThetaK);
                        h_nrecoL.Fill(CosThetaL);
                        h_nrecoK.Fill(CosThetaK);
                    }
                }

            }
        }
    }
    
    // Calculate efficiency
    TH1F h_recL("h_recL","",nLBins,-1,1);
    for (int i = 1; i <= nLBins; i++) {
        h_recL.SetBinContent(i,h_nrecoL.GetBinContent(i)/h_naccL.GetBinContent(i));
        h_recL.SetBinError(i,sqrt(h_recL.GetBinContent(i)*(1-h_recL.GetBinContent(i))/h_naccL.GetBinContent(i)));
    }
    TH1F h_recK("h_recK","",nKBins,-1,1);
    for (int i = 1; i <= nKBins; i++) {
        h_recK.SetBinContent(i,h_nrecoK.GetBinContent(i)/h_naccK.GetBinContent(i));
        h_recK.SetBinError(i,sqrt(h_recK.GetBinContent(i)*(1-h_recK.GetBinContent(i))/h_naccK.GetBinContent(i)));
    }

    TH2F h2_eff("h2_eff","",6,thetaLBins,5,thetaKBins);
    for (int i = 1; i <= 6; i++) {//L
        for (int j = 1; j <= 5; j++) {//K
            // Build from MC samples
            if (h2_nacc.GetBinContent(i,j) == 0 || h2_nreco.GetBinContent(i,j) == 0) {
                printf("WARNING: Efficiency(%d,%d)=0, set error to be 1.\n",i,j);
                h2_eff.SetBinContent(i,j,0.);
                h2_eff.SetBinError(i,j,1.);
            }else{
                h2_eff.SetBinContent(i,j,h2_nreco.GetBinContent(i,j)/h2_nacc.GetBinContent(i,j)*h2_acc->GetBinContent(i,j));
                h2_eff.SetBinError(i,j,h2_eff.GetBinContent(i,j)*sqrt(-1./h2_nacc.GetBinContent(i,j)+1./h2_nreco.GetBinContent(i,j)+pow(h2_acc->GetBinError(i,j)/h2_acc->GetBinContent(i,j),2)));
                printf("INFO: Efficiency(%d,%d)=%f +- %f.\n",i,j,h2_eff.GetBinContent(i,j),h2_eff.GetBinError(i,j));
            }
        }
    }

    TH1F h_effL("h_effL","",nLBins,-1,1);
    for (int i = 1; i <= nLBins; i++) {
        if (h_naccL.GetBinContent(i) == 0 || h_nrecoL.GetBinContent(i) == 0) {
            printf("WARNING: EfficiencyL(%d)=0, set error to be 1.\n",i);
            h_effL.SetBinContent(i,0.);
            h_effL.SetBinError(i,1.);
        }else{
            //h_effL.SetBinContent(i,h_nrecoL.GetBinContent(i)/h_naccL.GetBinContent(i)*h_accL->GetBinContent(i));
            //h_effL.SetBinError(i,h_effL.GetBinContent(i)*sqrt(-1./h_naccL.GetBinContent(i)+1./h_nrecoL.GetBinContent(i)+pow(h_accL->GetBinError(i)/h_accL->GetBinContent(i),2)));
            h_effL.SetBinContent(i,h_recL.GetBinContent(i)*h_accL->GetBinContent(i));
            h_effL.SetBinError(i,h_effL.GetBinContent(i)*h_recL.GetBinError(i)/h_recL.GetBinContent(i));
            //printf("INFO: EfficiencyL(%d)=%f +- %f.\n",i,h_effL.GetBinContent(i),h_effL.GetBinError(i));
        }
    }
    h_effL.SetStats(0);
    h_effL.SetMinimum(0.);
    h_effL.SetXTitle("CosThetaL");
    h_effL.SetYTitle("Efficiency");

    TH1F h_effK("h_effK","",nKBins,-1,1);
    for (int i = 1; i <= nKBins; i++) {
        if (h_naccK.GetBinContent(i) == 0 || h_nrecoK.GetBinContent(i) == 0) {
            printf("WARNING: EfficiencyK(%d)=0, set error to be 1.\n",i);
            h_effK.SetBinContent(i,0.);
            h_effK.SetBinError(i,1.);
        }else{
            //h_effK.SetBinContent(i,h_nrecoK.GetBinContent(i)/h_naccK.GetBinContent(i)*h_accK->GetBinContent(i));
            //h_effK.SetBinError(i,h_effK.GetBinContent(i)*sqrt(-1./h_naccK.GetBinContent(i)+1./h_nrecoK.GetBinContent(i)+pow(h_accK->GetBinError(i)/h_accK->GetBinContent(i),2)));
            h_effK.SetBinContent(i,h_recK.GetBinContent(i)*h_accK->GetBinContent(i));
            h_effK.SetBinError(i,h_effK.GetBinContent(i)*h_recK.GetBinError(i)/h_recK.GetBinContent(i));
            //printf("INFO: EfficiencyK(%d)=%f +- %f.\n",i,h_effK.GetBinContent(i),h_effK.GetBinError(i));
        }
    }
    h_effK.SetStats(0);
    h_effK.SetMinimum(0.);
    h_effK.SetXTitle("CosThetaK");
    h_effK.SetYTitle("Efficiency");

    // Quick check order0 (decoupled)
    TF1 *f_effK_ord0 = new TF1("f_effK_ord0","pol6",-1,1);//y, pol6
    TF1 *f_effL_ord0 = 0;
    switch ( iBin ){
        case 0:
        case 1:
        case 8:
        case 12://summaryBin[1]
            f_effL_ord0 = new TF1("f_effL_ord0","[2]*exp(-0.5*((x-[0])/[1])**2)+[5]*exp(-0.5*((x-[3])/[4])**2)+[8]*exp(-0.5*((x-[6])/[7])**2)",-1,1);//x
            f_effL_ord0->SetParameter(1,0.5);// width must be non-zero.
            f_effL_ord0->SetParameter(4,0.5);// width must be non-zero.
            f_effL_ord0->SetParameter(7,0.5);// width must be non-zero.
            break;
        default:
            f_effL_ord0 = new TF1("f_effL_ord0","pol6",-1,1);//x

    }
    h_effL.Fit("f_effL_ord0","S");
    h_effK.Fit("f_effK_ord0","S");

    // Using pure TMinuit for order1+
    int nPar = 21;
    TMinuit *gMinuit = new TMinuit(nPar);
    h2_fcn = &h2_eff;
    gMinuit->SetFCN(fcn_binnedChi2_2D);
    
    // 3 Gaussians or 6th order polynomial as 0th order
    // (0~3th order Legendre poly of CosThetaK)*(0~6th order power poly of CosThetaL) as 1st order
    TString f2_model_format_ord0 = TString::Format("(%e*exp(-0.5*((x-(%e))/%e)**2)%+e*exp(-0.5*((x-(%e))/%e)**2)%+e*exp(-0.5*((x-(%e))/%e)**2))*(%e%+e*y%+e*y**2%+e*y**3%+e*y**4%+e*y**5%+e*y**6)",\
                f_effL_ord0->GetParameter(2),\
                f_effL_ord0->GetParameter(0),\
                f_effL_ord0->GetParameter(1),\
                f_effL_ord0->GetParameter(5),\
                f_effL_ord0->GetParameter(3),\
                f_effL_ord0->GetParameter(4),\
                f_effL_ord0->GetParameter(8),\
                f_effL_ord0->GetParameter(6),\
                f_effL_ord0->GetParameter(7),\
                f_effK_ord0->GetParameter(0),\
                f_effK_ord0->GetParameter(1),\
                f_effK_ord0->GetParameter(2),\
                f_effK_ord0->GetParameter(3),\
                f_effK_ord0->GetParameter(4),\
                f_effK_ord0->GetParameter(5),\
                f_effK_ord0->GetParameter(6));
    switch(iBin){
        case 0:
        case 1:
        case 8:
        case 12:// summaryBin[1]
            break;
        default:
            f2_model_format_ord0 = TString::Format("(%e%+e*x%+e*x**2%+e*x**3%+e*x**4%+e*x**5%+e*x**6)*(%e%+e*y%+e*y**2%+e*y**3%+e*y**4%+e*y**5%+e*y**6)",\
                f_effL_ord0->GetParameter(0),\
                f_effL_ord0->GetParameter(1),\
                f_effL_ord0->GetParameter(2),\
                f_effL_ord0->GetParameter(3),\
                f_effL_ord0->GetParameter(4),\
                f_effL_ord0->GetParameter(5),\
                f_effL_ord0->GetParameter(6),\
                f_effK_ord0->GetParameter(0),\
                f_effK_ord0->GetParameter(1),\
                f_effK_ord0->GetParameter(2),\
                f_effK_ord0->GetParameter(3),\
                f_effK_ord0->GetParameter(4),\
                f_effK_ord0->GetParameter(5),\
                f_effK_ord0->GetParameter(6));
            break;
    }
    printf("DEBUG\t\t: f2_model_format_ord0=%s\n",f2_model_format_ord0.Data());
    TString f2_model_format_ord1 = "([0]+[1]*y+[2]*(3*y**2-1)/2+[3]*(5*y**3-3*y)/2)+([4]+[5]*y+[6]*(3*y**2-1)/2+[7]*(5*y**3-3*y)/2)*x**2+([8]+[9]*y+[10]*(3*y**2-1)/2+[11]*(5*y**3-3*y)/2)*x**3+([12]+[13]*y+[14]*(3*y**2-1)/2+[15]*(5*y**3-3*y)/2)*x**4+([16]+[17]*y+[18]*(3*y**2-1)/2+[19]*(5*y**3-3*y)/2)*x**6";
    TF2 f2_model_ord0("f2_model_ord0",f2_model_format_ord0,-1.,1.,-1.,1.);
    TF2 f2_model_ord1("f2_model_ord1",f2_model_format_ord1,-1.,1.,-1.,1.);
    TF2 f2_model("f2_model",TString::Format("[20]*%s*(1+%s)",f2_model_format_ord0.Data(),f2_model_format_ord1.Data()).Data(),-1.,1.,-1.,1.);
    f2_fcn = &f2_model;
    gMinuit->DefineParameter( 0, "k0l0",   0.,  1E-4,    -1E+1, 1E+1);
    gMinuit->DefineParameter( 1, "k1l0",   0.,  1E-4,    -1E+1, 1E+1);
    gMinuit->DefineParameter( 2, "k2l0",   0.,  1E-4,    -1E+1, 1E+1);
    gMinuit->DefineParameter( 3, "k3l0",   0.,  1E-4,    -1E+1, 1E+1);
    gMinuit->DefineParameter( 4, "k0l2",   0.,  1E-4,    -1E+1, 1E+1);
    gMinuit->DefineParameter( 5, "k1l2",   0.,  1E-4,    -1E+1, 1E+1);
    gMinuit->DefineParameter( 6, "k2l2",   0.,  1E-4,    -1E+1, 1E+1);
    gMinuit->DefineParameter( 7, "k3l2",   0.,  1E-4,    -1E+1, 1E+1);
    gMinuit->DefineParameter( 8, "k0l3",   0.,  1E-4,    -1E+1, 1E+1);
    gMinuit->DefineParameter( 9, "k1l3",   0.,  1E-4,    -1E+1, 1E+1);
    gMinuit->DefineParameter(10, "k2l3",   0.,  1E-4,    -1E+1, 1E+1);
    gMinuit->DefineParameter(11, "k3l3",   0.,  1E-4,    -1E+1, 1E+1);
    gMinuit->DefineParameter(12, "k0l4",   0.,  1E-4,    -1E+1, 1E+1);
    gMinuit->DefineParameter(13, "k1l4",   0.,  1E-4,    -1E+1, 1E+1);
    gMinuit->DefineParameter(14, "k2l4",   0.,  1E-4,    -1E+1, 1E+1);
    gMinuit->DefineParameter(15, "k3l4",   0.,  1E-4,    -1E+1, 1E+1);
    gMinuit->DefineParameter(16, "k0l6",   0.,  1E-4,    -1E+2, 1E+3);
    gMinuit->DefineParameter(17, "k1l6",   0.,  1E-4,    -1E+2, 1E+3);
    gMinuit->DefineParameter(18, "k2l6",   0.,  1E-4,    -1E+2, 1E+3);
    gMinuit->DefineParameter(19, "k3l6",   0.,  1E-4,    -1E+2, 1E+3);
    gMinuit->DefineParameter(20, "ord0",8000.,  1E-2,       10, 1E+6);

    // Calculte normalize factor for 0th order
    for (int i = 1; i < 21; i++) {
        gMinuit->Command(TString::Format("SET PARM %d 0",i));
        gMinuit->Command(TString::Format("FIX %d",i));
    }
    gMinuit->Command("MINI");
    gMinuit->Command("MINI");
    gMinuit->Command("MINOS");
    gMinuit->Command("FIX 21");

    // Start processing 1st order
    for (int i = 1; i < 21; i++) {gMinuit->Command(TString::Format("REL %d",i));}
    if (iBin == 0 ) {
        gMinuit->Command("SET PARM 9 0");
        gMinuit->Command("SET PARM 10 0");
        gMinuit->Command("SET PARM 11 0");
        gMinuit->Command("SET PARM 12 0");
        gMinuit->Command("FIX 9");
        gMinuit->Command("FIX 10");
        gMinuit->Command("FIX 11");
        gMinuit->Command("FIX 12");

        gMinuit->Command("FIX 13");
        gMinuit->Command("FIX 14");
        gMinuit->Command("FIX 15");
        gMinuit->Command("FIX 16");
        gMinuit->Command("FIX 17");
        gMinuit->Command("FIX 18");
        gMinuit->Command("FIX 19");
        gMinuit->Command("FIX 20");
    }else if (iBin == 1) {
        gMinuit->Command("SET PARM 9 0");
        gMinuit->Command("SET PARM 10 0");
        gMinuit->Command("SET PARM 11 0");
        gMinuit->Command("SET PARM 12 0");
        gMinuit->Command("SET PARM 17 0");
        gMinuit->Command("SET PARM 18 0");
        gMinuit->Command("SET PARM 19 0");
        gMinuit->Command("SET PARM 20 0");
        gMinuit->Command("FIX 9");
        gMinuit->Command("FIX 10");
        gMinuit->Command("FIX 11");
        gMinuit->Command("FIX 12");
        gMinuit->Command("FIX 17");
        gMinuit->Command("FIX 18");
        gMinuit->Command("FIX 19");
        gMinuit->Command("FIX 20");
    }else if (iBin == 7 || iBin == 9) {
        gMinuit->Command("SET PARM 13 0");
        gMinuit->Command("SET PARM 14 0");
        gMinuit->Command("SET PARM 15 0");
        gMinuit->Command("SET PARM 16 0");
        gMinuit->Command("SET PARM 17 0");
        gMinuit->Command("SET PARM 18 0");
        gMinuit->Command("SET PARM 19 0");
        gMinuit->Command("SET PARM 20 0");
        gMinuit->Command("FIX 13");
        gMinuit->Command("FIX 14");
        gMinuit->Command("FIX 15");
        gMinuit->Command("FIX 16");
        gMinuit->Command("FIX 17");
        gMinuit->Command("FIX 18");
        gMinuit->Command("FIX 19");
        gMinuit->Command("FIX 20");
    }else if (iBin == 10 || iBin == 13) {
        gMinuit->Command("SET PARM 9 0");
        gMinuit->Command("SET PARM 13 0");
        gMinuit->Command("SET PARM 15 0");
        gMinuit->Command("SET PARM 16 0");
        gMinuit->Command("SET PARM 18 0");
        gMinuit->Command("SET PARM 19 0");
        gMinuit->Command("SET PARM 20 0");
        gMinuit->Command("FIX 9");
        gMinuit->Command("FIX 13");
        gMinuit->Command("FIX 15");
        gMinuit->Command("FIX 16");
        gMinuit->Command("FIX 18");
        gMinuit->Command("FIX 19");
        gMinuit->Command("FIX 20");
    }else if (iBin > 1 || iBin == summaryBins[1]) {
        gMinuit->Command("SET PARM 17 0");
        gMinuit->Command("SET PARM 18 0");
        gMinuit->Command("SET PARM 19 0");
        gMinuit->Command("SET PARM 20 0");
        gMinuit->Command("FIX 17");
        gMinuit->Command("FIX 18");
        gMinuit->Command("FIX 19");
        gMinuit->Command("FIX 20");
    }else{
        gMinuit->Command("SET PARM 13 0");
        gMinuit->Command("SET PARM 14 0");
        gMinuit->Command("SET PARM 15 0");
        gMinuit->Command("SET PARM 16 0");
        gMinuit->Command("SET PARM 17 0");
        gMinuit->Command("SET PARM 18 0");
        gMinuit->Command("SET PARM 19 0");
        gMinuit->Command("SET PARM 20 0");
        gMinuit->Command("FIX 13");
        gMinuit->Command("FIX 14");
        gMinuit->Command("FIX 15");
        gMinuit->Command("FIX 16");
        gMinuit->Command("FIX 17");
        gMinuit->Command("FIX 18");
        gMinuit->Command("FIX 19");
        gMinuit->Command("FIX 20");
    }

    gMinuit->Command("MINI");
    gMinuit->Command("MINI");
    gMinuit->Command("IMPROVE");
    gMinuit->Command("MINOS");

    double arrPar[nPar];
    double arrParErr[nPar];
    for (int iPar = 0; iPar < nPar; iPar++) gMinuit->GetParameter(iPar,arrPar[iPar],arrParErr[iPar]);
    for (int iPar = 0; iPar < nPar; iPar++) f2_model.SetParameter(iPar,arrPar[iPar]);

    double covMatrix[nPar][nPar];
    gMinuit->mnemat(&covMatrix[0][0],nPar);
    TVectorD cenVec(nPar);
    double covMatrix1D[nPar*nPar];
    TMatrixDSym errMtx(0,nPar-1, covMatrix1D);//arrPar[nPar-1] is normalization.
    for (int iPar=0; iPar<nPar; iPar++){
        cenVec[iPar] = arrPar[iPar];
        for (int jPar=0; jPar<nPar; jPar++){
            if (arrPar[iPar] == 0 || arrPar[jPar] == 0){
                errMtx[iPar][jPar] = 1e-30; // Just pick an extermely small number.
            }else{
                errMtx[iPar][jPar] = covMatrix[iPar][jPar];
            }
        }
    }
    cenVec.Print();
    errMtx.Print();

    // Draw and write config
    TCanvas canvas("canvas");
    TLatex *latex = new TLatex();
    latex->SetNDC();
    double chi2Val=0;
    fcn_binnedChi2_2D(nPar, 0, chi2Val, arrPar, 0);
    printf("Chi2(Bin center)=%f \n",chi2Val);
        
    h_effL.Draw();
    if (iBin == 0) h_effL.SetMaximum(0.0003);
    canvas.Update();
    canvas.SaveSource(TString::Format("%s/accXrecoEff2_cosl_order0_bin%d.cc",plotpath.Data(),iBin));
    canvas.Print(TString::Format("%s/accXrecoEff2_cosl_order0_bin%d.pdf",plotpath.Data(),iBin));
    h_effK.Draw();
    canvas.Update();
    canvas.SaveSource(TString::Format("%s/accXrecoEff2_cosk_order0_bin%d.cc",plotpath.Data(),iBin));
    canvas.Print(TString::Format("%s/accXrecoEff2_cosk_order0_bin%d.pdf",plotpath.Data(),iBin));
        
    if (keepParam){
        //// Draw 1-D
        h_naccL.Draw();
        canvas.Update();
        canvas.SaveSource(TString::Format("%s/accXrecoEff2_naccL_bin%d.cc",plotpath.Data(),iBin));
        canvas.Print(TString::Format("%s/accXrecoEff2_naccL_bin%d.pdf",plotpath.Data(),iBin));
        h_naccK.Draw();
        canvas.Update();
        canvas.SaveSource(TString::Format("%s/accXrecoEff2_naccK_bin%d.cc",plotpath.Data(),iBin));
        canvas.Print(TString::Format("%s/accXrecoEff2_naccK_bin%d.pdf",plotpath.Data(),iBin));

        h_nrecoL.Draw();
        canvas.Update();
        canvas.SaveSource(TString::Format("%s/accXrecoEff2_nrecoL_bin%d.cc",plotpath.Data(),iBin));
        canvas.Print(TString::Format("%s/accXrecoEff2_nrecoL_bin%d.pdf",plotpath.Data(),iBin));
        h_nrecoK.Draw();
        canvas.Update();
        canvas.SaveSource(TString::Format("%s/accXrecoEff2_nrecoK_bin%d.cc",plotpath.Data(),iBin));
        canvas.Print(TString::Format("%s/accXrecoEff2_nrecoK_bin%d.pdf",plotpath.Data(),iBin));
        
        h_recL.SetStats(0);
        h_recL.SetMinimum(0.);
        h_recL.SetXTitle("CosThetaL");
        h_recL.Draw("TEXT");
        canvas.Update();
        canvas.SaveSource(TString::Format("%s/accXrecoEff2_recoEffL_bin%d.cc",plotpath.Data(),iBin));
        canvas.Print(TString::Format("%s/accXrecoEff2_recoEffL_bin%d.pdf",plotpath.Data(),iBin));
        
        h_recK.SetStats(0);
        h_recK.SetMinimum(0.);
        h_recK.SetXTitle("CosThetaK");
        h_recK.Draw("TEXT");
        canvas.Update();
        canvas.SaveSource(TString::Format("%s/accXrecoEff2_recoEffK_bin%d.cc",plotpath.Data(),iBin));
        canvas.Print(TString::Format("%s/accXrecoEff2_recoEffK_bin%d.pdf",plotpath.Data(),iBin));
        
        TH1F h_theoK("h_theoK" ,"h_theoK" ,nKBins,-1,1); 
        h_theoK.SetStats(0);
        h_theoK.SetMinimum(0.);
        h_theoK.SetXTitle("CosThetaK");
        h_theoK.SetYTitle("#Events / 0.2");
        for (int kBin = 1; kBin <= nKBins; kBin++) {
            h_theoK.SetBinContent(kBin,h_ngenK->GetBinContent(kBin)*h_effK.GetBinContent(kBin));
            if (h_effK.GetBinContent(kBin) != 0){
                h_theoK.SetBinError(kBin,h_theoK.GetBinContent(kBin)*h_effK.GetBinError(kBin)/h_effK.GetBinContent(kBin));
            }else{
                h_theoK.SetBinError(kBin,sqrt(h_theoK.GetBinContent(kBin)));
            }
        }
        h_theoK.Draw();
        canvas.Update();
        canvas.SaveSource(TString::Format("%s/accXrecoEff2_theoK_bin%d.cc",plotpath.Data(),iBin));
        canvas.Print(TString::Format("%s/accXrecoEff2_theoK_bin%d.pdf",plotpath.Data(),iBin));
        
        TH1F h_theoL("h_theoL" ,"h_theoL" ,nLBins,-1,1); 
        h_theoL.SetStats(0);
        h_theoL.SetMinimum(0.);
        h_theoL.SetXTitle("CosThetaL");
        h_theoL.SetYTitle("#Events / 0.2");
        for (int lBin = 1; lBin <= nLBins; lBin++) {
            h_theoL.SetBinContent(lBin,h_ngenL->GetBinContent(lBin)*h_effL.GetBinContent(lBin));
            if (h_effL.GetBinContent(lBin) != 0){
                h_theoL.SetBinError(lBin,h_theoL.GetBinContent(lBin)*h_effL.GetBinError(lBin)/h_effL.GetBinContent(lBin));
            }else{
                h_theoL.SetBinError(lBin,sqrt(h_theoL.GetBinContent(lBin)));
            }
        }
        h_theoL.Draw();
        canvas.Update();
        canvas.SaveSource(TString::Format("%s/accXrecoEff2_theoL_bin%d.cc",plotpath.Data(),iBin));
        canvas.Print(TString::Format("%s/accXrecoEff2_theoL_bin%d.pdf",plotpath.Data(),iBin));
        
        //// Draw 2-D
        h2_eff.SetMinimum(0.);
        h2_eff.SetTitleOffset(2,"XY");
        h2_eff.SetXTitle("CosThetaL");
        h2_eff.SetYTitle("CosThetaK");
        h2_eff.SetStats(0);
        h2_eff.SetMaximum(effUpperBound);
        h2_eff.Draw("LEGO2");
        latex->DrawLatex(0.35,0.95,TString::Format("#varepsilon in Bin%d",iBin));
        
        // Draw FitResult
        f2_model.SetTitle("");
        f2_model.SetMaximum(effUpperBound);
        f2_model.SetLineWidth(1);
        f2_model.Draw("SURF SAME ");
        canvas.Update();
        canvas.SaveSource(TString::Format("%s/accXrecoEff2_2D_bin%d.cc",plotpath.Data(),iBin));
        canvas.Print(TString::Format("%s/accXrecoEff2_2D_bin%d.pdf",plotpath.Data(),iBin));

        //// Draw compare
        TH2F h2_compFit("h2_compFit","",6,thetaLBins,5,thetaKBins);
        h2_compFit.SetTitleOffset(2,"XY");
        h2_compFit.SetXTitle("CosThetaL");
        h2_compFit.SetYTitle("CosThetaK");
        TH2F h2_pullFit("h2_pullFit","",6,thetaLBins,5,thetaKBins);
        h2_pullFit.SetTitleOffset(1,"XY");
        h2_pullFit.SetXTitle("CosThetaL");
        h2_pullFit.SetYTitle("CosThetaK");
        for (int i = 1; i <= 6; i++) {//thetaL
            for (int j = 1; j <= 5; j++) {//thetaK
                if (h2_eff.GetBinContent(i,j) != 0){
                    h2_compFit.SetBinContent(i,j,f2_model.Eval(h2_eff.GetXaxis()->GetBinCenter(i),h2_eff.GetYaxis()->GetBinCenter(j))/h2_eff.GetBinContent(i,j));
                    double _xlo = h2_eff.GetXaxis()->GetBinLowEdge(i);
                    double _xhi = h2_eff.GetXaxis()->GetBinUpEdge(i);
                    double _ylo = h2_eff.GetYaxis()->GetBinLowEdge(j);
                    double _yhi = h2_eff.GetYaxis()->GetBinUpEdge(j);
                    h2_pullFit.SetBinContent(i,j,(f2_model.Integral(_xlo,_xhi,_ylo,_yhi)/(_xhi-_xlo)/(_yhi-_ylo)-h2_eff.GetBinContent(i,j))/h2_eff.GetBinError(i,j));
                }else{
                    h2_compFit.SetBinContent(i,j,0.);
                    h2_pullFit.SetBinContent(i,j,0.);
                }
            }
        }
        h2_compFit.SetMinimum(0.);
        h2_compFit.SetStats(0);
        h2_compFit.Draw("LEGO2");
        latex->DrawLatex(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
        latex->DrawLatex(0.3,0.95,TString::Format("#varepsilon_{fit} / #varepsilon_{measured} in Bin%d",iBin));
        canvas.Update();
        canvas.SaveSource(TString::Format("%s/accXrecoEff2_compFit_2D_bin%d.cc",plotpath.Data(),iBin));
        canvas.Print(TString::Format("%s/accXrecoEff2_compFit_2D_bin%d.pdf",plotpath.Data(),iBin));
        
        h2_pullFit.SetStats(0);
        h2_pullFit.Draw("COLZ TEXT");
        latex->DrawLatex(0.01,0.95,TString::Format("#chi^{2} = %f",chi2Val));
        latex->DrawLatex(0.3,0.95,TString::Format("(#varepsilon_{fit} - #varepsilon_{measured})/Error in Bin%d",iBin));
        canvas.Update();
        canvas.SaveSource(TString::Format("%s/accXrecoEff2_pullFit_2D_bin%d.cc",plotpath.Data(),iBin));
        canvas.Print(TString::Format("%s/accXrecoEff2_pullFit_2D_bin%d.pdf",plotpath.Data(),iBin));

        // Draw significance of deviation
        TH1F h_pull("Deviation/Error","",15,-3.,3.);
        h_pull.SetXTitle("Significance of deviation");
        for (int i = 1; i <= 6; i++) {//thetaL
            for (int j = 1; j <= 5; j++) {//thetaK
                double _xlo = h2_eff.GetXaxis()->GetBinLowEdge(i);
                double _xhi = h2_eff.GetXaxis()->GetBinUpEdge(i);
                double _ylo = h2_eff.GetYaxis()->GetBinLowEdge(j);
                double _yhi = h2_eff.GetYaxis()->GetBinUpEdge(j);
                if (h2_eff.GetBinContent(i,j) != 0){
                    h_pull.Fill((f2_model.Integral(_xlo,_xhi,_ylo,_yhi)/(_xhi-_xlo)/(_yhi-_ylo)-h2_eff.GetBinContent(i,j))/h2_eff.GetBinError(i,j));
                }
            }
        }
        h_pull.Draw();
        canvas.Update();
        canvas.SaveSource(TString::Format("%s/accXrecoEff2_sigma_bin%d.cc",plotpath.Data(),iBin));
        canvas.Print(TString::Format("%s/accXrecoEff2_sigma_bin%d.pdf",plotpath.Data(),iBin));

        delete latex;
        
    }

    // prepare output
    string out_accXrecoEff2_ord0;
    out_accXrecoEff2_ord0 = TString::Format("%e*(%s)",arrPar[20],f2_model_format_ord0.Data());
    while(out_accXrecoEff2_ord0.find("*y") !=std::string::npos ){
        out_accXrecoEff2_ord0.replace(out_accXrecoEff2_ord0.find("*y"), 2, "*CosThetaK");
    }
    while(out_accXrecoEff2_ord0.find("*x") !=std::string::npos ){
        out_accXrecoEff2_ord0.replace(out_accXrecoEff2_ord0.find("*x"), 2, "*CosThetaL");
    }
    while(out_accXrecoEff2_ord0.find("(x") !=std::string::npos ){
        out_accXrecoEff2_ord0.replace(out_accXrecoEff2_ord0.find("(x"), 2, "(CosThetaL");
    }
    printf("\"%s\",\n",out_accXrecoEff2_ord0.c_str());
    writeParam(iBin,"f_accXrecoEff_ord0",out_accXrecoEff2_ord0);
    writeParam(iBin,"accXrecoEff2",arrPar,20,true);
    writeParam(iBin,"accXrecoEff2Err",arrParErr,20,true);
    // write to wspace as well
    if (keepParam){
        RooRealVar CosThetaK("CosThetaK", "cos#theta_{K}", -1., 1.);
        RooRealVar CosThetaL("CosThetaL", "cos#theta_{L}", -1., 1.);
        //RooRealVar fl("fl", "F_{L}", genFl[iBin], 0., 1.);
        //RooRealVar afb("afb", "A_{FB}", genAfb[iBin], -1., 1.);
        RooRealVar fl("fl", "F_{L}", genFl[iBin], -100, 100.);// unbounded fl
        RooRealVar afb("afb", "A_{FB}", genAfb[iBin], -100., 100.);// unbounded afb
        RooRealVar fs("fs","F_{S}",0.0129254,-0.3,0.3);
        //RooRealVar as("as","A_{S}",-0.0975919,-0.3,0.3);
        RooRealVar as("as","A_{S}",-1.1,1.1);
        RooRealVar recK0L0("recK0L0","recK0L0",arrPar[ 0]);
        RooRealVar recK1L0("recK1L0","recK1L0",arrPar[ 1]);
        RooRealVar recK2L0("recK2L0","recK2L0",arrPar[ 2]);
        RooRealVar recK3L0("recK3L0","recK3L0",arrPar[ 3]);
        RooRealVar recK0L2("recK0L2","recK0L2",arrPar[ 4]);
        RooRealVar recK1L2("recK1L2","recK1L2",arrPar[ 5]);
        RooRealVar recK2L2("recK2L2","recK2L2",arrPar[ 6]);
        RooRealVar recK3L2("recK3L2","recK3L2",arrPar[ 7]);
        RooRealVar recK0L3("recK0L3","recK0L3",arrPar[ 8]);
        RooRealVar recK1L3("recK1L3","recK1L3",arrPar[ 9]);
        RooRealVar recK2L3("recK2L3","recK2L3",arrPar[10]);
        RooRealVar recK3L3("recK3L3","recK3L3",arrPar[11]);
        RooRealVar recK0L4("recK0L4","recK0L4",arrPar[12]);
        RooRealVar recK1L4("recK1L4","recK1L4",arrPar[13]);
        RooRealVar recK2L4("recK2L4","recK2L4",arrPar[14]);
        RooRealVar recK3L4("recK3L4","recK3L4",arrPar[15]);
        RooRealVar recK0L6("recK0L6","recK0L6",arrPar[16]);
        RooRealVar recK1L6("recK1L6","recK1L6",arrPar[17]);
        RooRealVar recK2L6("recK2L6","recK2L6",arrPar[18]);
        RooRealVar recK3L6("recK3L6","recK3L6",arrPar[19]);
        RooRealVar effNorm("effNorm","effNorm",arrPar[20]);
        recK0L0.setError(arrPar[ 0]);
        recK1L0.setError(arrPar[ 1]);
        recK2L0.setError(arrPar[ 2]);
        recK3L0.setError(arrPar[ 3]);
        recK0L2.setError(arrPar[ 4]);
        recK1L2.setError(arrPar[ 5]);
        recK2L2.setError(arrPar[ 6]);
        recK3L2.setError(arrPar[ 7]);
        recK0L3.setError(arrPar[ 8]);
        recK1L3.setError(arrPar[ 9]);
        recK2L3.setError(arrPar[10]);
        recK3L3.setError(arrPar[11]);
        recK0L4.setError(arrPar[12]);
        recK1L4.setError(arrPar[13]);
        recK2L4.setError(arrPar[14]);
        recK3L4.setError(arrPar[15]);
        recK0L6.setError(arrPar[16]);
        recK1L6.setError(arrPar[17]);
        recK2L6.setError(arrPar[18]);
        recK3L6.setError(arrPar[19]);
        effNorm.setError(arrPar[20]);
        RooArgSet f_sigA_argset(CosThetaL,CosThetaK);
        f_sigA_argset.add(RooArgSet(fl,afb,fs,as));
        f_sigA_argset.add(RooArgSet(recK0L0,recK1L0,recK2L0,recK3L0));
        f_sigA_argset.add(RooArgSet(recK0L2,recK1L2,recK2L2,recK3L2));
        RooArgSet f_sigA_MK_argset(CosThetaK);
        f_sigA_MK_argset.add(RooArgSet(fl,fs,as));
        TString f_sigA_format;
        TString f_sigA_MK_format;
        //TString f_ang_format = "9/16*((2/3*fs+4/3*as*CosThetaK)*(1-CosThetaL**2)+(1-fs)*(2*fl*CosThetaK**2*(1-CosThetaL**2)+1/2*(1-fl)*(1-CosThetaK**2)*(1+CosThetaL**2)+4/3*afb*(1-CosThetaK**2)*CosThetaL))";
        //TString f_ang_format = "9/16*((2/3*fs+4/3*as*CosThetaK)*(1-CosThetaL**2)+(1-fs)*(2*(0.5+TMath::ATan(fl)/TMath::Pi())*CosThetaK**2*(1-CosThetaL**2)+0.5*(0.5-TMath::ATan(fl)/TMath::Pi())*(1-CosThetaK**2)*(1+CosThetaL**2)+4/3*(3/2*(1/2-TMath::ATan(fl)/TMath::Pi())*TMath::ATan(afb)/TMath::Pi())*(1-CosThetaK**2)*CosThetaL))";// unbounded fl, afb
        TString f_ang_format = "9/16*((2/3*fs+4/3*as*2*sqrt(3*fs*(1-fs)*(0.5+TMath::ATan(fl)/TMath::Pi()))*CosThetaK)*(1-CosThetaL**2)+(1-fs)*(2*(0.5+TMath::ATan(fl)/TMath::Pi())*CosThetaK**2*(1-CosThetaL**2)+0.5*(0.5-TMath::ATan(fl)/TMath::Pi())*(1-CosThetaK**2)*(1+CosThetaL**2)+4/3*(3/2*(1/2-TMath::ATan(fl)/TMath::Pi())*TMath::ATan(afb)/TMath::Pi())*(1-CosThetaK**2)*CosThetaL))";// unbounded fl, afb, transformed as.

        TString f_accXrecoEff2_ord0 = out_accXrecoEff2_ord0;
        TString f_accXrecoEff2_format, f_accXrecoEff2_L0, f_accXrecoEff2_L2, f_accXrecoEff2_L3, f_accXrecoEff2_L4, f_accXrecoEff2_L6;
        f_accXrecoEff2_L0 = "(recK0L0+recK1L0*CosThetaK+recK2L0*(3*CosThetaK**2-1)/2+recK3L0*(5*CosThetaK**3-3*CosThetaK)/2)";
        f_accXrecoEff2_L2 = "(recK0L2+recK1L2*CosThetaK+recK2L2*(3*CosThetaK**2-1)/2+recK3L2*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**2";
        f_accXrecoEff2_L3 = "(recK0L3+recK1L3*CosThetaK+recK2L3*(3*CosThetaK**2-1)/2+recK3L3*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**3";
        f_accXrecoEff2_L4 = "(recK0L4+recK1L4*CosThetaK+recK2L4*(3*CosThetaK**2-1)/2+recK3L4*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**4";
        f_accXrecoEff2_L6 = "(recK0L6+recK1L6*CosThetaK+recK2L6*(3*CosThetaK**2-1)/2+recK3L6*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**6";

        TString f_ang_MK_format = "3/4*((2/3*fs+4/3*as*CosThetaK)+(1-fs)*(2*(1/2+TMath::ATan(fl)/TMath::Pi())*CosThetaK**2+(1/2-TMath::ATan(fl)/TMath::Pi())*(1-CosThetaK**2)))";// unbounded fl, integrated over cosThetaL
        TString f_accXrecoEff2_MK_format = TString::Format("(%e%+e*CosThetaK%+e*CosThetaK**2%+e*CosThetaK**3%+e*CosThetaK**4%+e*CosThetaK**5%+e*CosThetaK**6)",
            f_effK_ord0->GetParameter(0),\
            f_effK_ord0->GetParameter(1),\
            f_effK_ord0->GetParameter(2),\
            f_effK_ord0->GetParameter(3),\
            f_effK_ord0->GetParameter(4),\
            f_effK_ord0->GetParameter(5),\
            f_effK_ord0->GetParameter(6));
        f_sigA_MK_format = TString::Format("%s*%s", f_accXrecoEff2_MK_format.Data(), f_ang_MK_format.Data());

        if (iBin == 0 || iBin == 10 || iBin == 13 || iBin == summaryBins[1]) {
            f_sigA_argset.add(RooArgSet(recK0L4,recK1L4,recK2L4,recK3L4));
            f_sigA_argset.add(RooArgSet(recK0L6,recK1L6,recK2L6,recK3L6));
            f_sigA_format = TString::Format("%s*(1+%s+%s+%s+%s)*%s",f_accXrecoEff2_ord0.Data(),f_accXrecoEff2_L0.Data(),f_accXrecoEff2_L2.Data(),f_accXrecoEff2_L4.Data(),f_accXrecoEff2_L6.Data(),f_ang_format.Data());
        }else if (iBin == 1) {
            f_sigA_argset.add(RooArgSet(recK0L4,recK1L4,recK2L4,recK3L4));
            f_sigA_format = TString::Format("%s*(1+%s+%s+%s)*%s",f_accXrecoEff2_ord0.Data(),f_accXrecoEff2_L0.Data(),f_accXrecoEff2_L2.Data(),f_accXrecoEff2_L4.Data(),f_ang_format.Data());
        }else if (iBin == 7) {
            f_sigA_argset.add(RooArgSet(recK0L3,recK1L3,recK2L3,recK3L3));
            f_sigA_format = TString::Format("%s*(1+%s+%s+%s)*%s",f_accXrecoEff2_ord0.Data(),f_accXrecoEff2_L0.Data(),f_accXrecoEff2_L2.Data(),f_accXrecoEff2_L3.Data(),f_ang_format.Data());
        }else if ((iBin > 1 && iBin < 6) ){ 
            f_sigA_argset.add(RooArgSet(recK0L3,recK1L3,recK2L3,recK3L3));
            f_sigA_argset.add(RooArgSet(recK0L4,recK1L4,recK2L4,recK3L4));
            f_sigA_format = TString::Format("%s*(1+%s+%s+%s+%s)*%s",f_accXrecoEff2_ord0.Data(),f_accXrecoEff2_L0.Data(),f_accXrecoEff2_L2.Data(),f_accXrecoEff2_L3.Data(),f_accXrecoEff2_L4.Data(),f_ang_format.Data());
        }else{
            f_sigA_argset.add(RooArgSet(recK0L3,recK1L3,recK2L3,recK3L3));
            f_sigA_format = TString::Format("%s*(1+%s+%s+%s)*%s",f_accXrecoEff2_ord0.Data(),f_accXrecoEff2_L0.Data(),f_accXrecoEff2_L2.Data(),f_accXrecoEff2_L3.Data(),f_ang_format.Data());
        }
            // angular map of signal
        RooGenericPdf f_sigA("f_sigA", f_sigA_format, f_sigA_argset);
        RooGenericPdf f_sigA_MK("f_sigA_MK", f_sigA_MK_format, f_sigA_MK_argset);

            // angular map of signal
        RooWorkspace *wspace = new RooWorkspace("wspace","wspace");
        wspace->import(f_sigA);
        wspace->import(f_sigA_MK);
        wspace->import(effNorm);
        wspace->import(errMtx,"errMtx");
        wspace->writeToFile(TString::Format("%s/wspace_sigA_bin%d.root",owspacepath.Data(),iBin),true);
    }

    delete gMinuit;// delete before return.
    switchRedirectStdio("_stdio");
    return out_accXrecoEff2_ord0.c_str();
}//}}}  

//_________________________________________________________________________________
void angular2D_bin(int iBin, const char outfile[] = "angular2D")
{//{{{
    // Remark: You must use RooFit!! It's better in unbinned fit.
    //         Extended ML fit is adopted by Mauro, just follow!!
    if (iBin == 3 || iBin == 5) return;
    switchRedirectStdio(TString::Format("%s/angular2D_bin_stdout_bin%d.log",odatacardpath.Data(),iBin).Data(),"w",stdout);
    switchRedirectStdio(TString::Format("%s/angular2D_bin_stderr_bin%d.log",odatacardpath.Data(),iBin).Data(),"w",stderr);

    // Read data
    ch->SetBranchStatus("*",0);
    ch->SetBranchStatus("Bmass"         , 1);
    ch->SetBranchStatus("Mumumass"      , 1);
    ch->SetBranchStatus("Mumumasserr"   , 1);
    ch->SetBranchStatus("Kstarmass"     , 1);
    ch->SetBranchStatus("CosTheta*"     , 1);
    ch->SetBranchStatus("Q2"            , 1);
    ch->SetBranchStatus("Triggers"      , 1);
    RooRealVar Bmass("Bmass","M_{K^{*}#Mu#Mu}",5.1,5.6);
    RooRealVar CosThetaK("CosThetaK"     , "cos#theta_{K}"       , -1. , 1.   ) ;
    RooRealVar CosThetaL("CosThetaL"     , "cos#theta_{L}"       , -1. , 1.   ) ;
    RooRealVar Mumumass("Mumumass","M^{#mu#mu}",0.,10.);
    RooRealVar Mumumasserr("Mumumasserr","Error of M^{#mu#mu}",0.,10.);
    RooRealVar Kstarmass("Kstarmass","M_{K^{*}}", 0,1.5);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    RooRealVar Triggers("Triggers","",0,100);

    TFile *f_wspace_sigA = new TFile(TString::Format("%s/wspace_sigA_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace_sigA = (RooWorkspace*)f_wspace_sigA->Get("wspace");
    RooGenericPdf *f_sigA = 0;
    RooRealVar *fl = 0;
    RooRealVar *fs = 0;
    RooRealVar *afb = 0;
    RooRealVar *as = 0;
    if (wspace_sigA){
        f_sigA = (RooGenericPdf*)wspace_sigA->pdf("f_sigA");
        fl = (RooRealVar*)wspace_sigA->var("fl");
        fs = (RooRealVar*)wspace_sigA->var("fs");
        afb = (RooRealVar*)wspace_sigA->var("afb");
        as = (RooRealVar*)wspace_sigA->var("as");
        fl      ->setRange(-100,100);// unbounded fl
        afb     ->setRange(-100,100);// unbounded afb
        fl      ->setConstant(kFALSE);
        afb     ->setConstant(kFALSE);
        fs->setVal(0);
        fs->setConstant(kTRUE);
        as->setVal(0);
        as->setConstant(kTRUE);
    }else{
        printf("ERROR\t\t: Please have wspace_sigA_bin?.root prepared.\n");
        return;
    }

    RooRealVar nsig("nsig","nsig",1E4,1E2,1E8);
    RooExtendPdf f_ext("f_ext","f_ext",*f_sigA,nsig);
    
    // Get data and apply unbinned fit
    int mumuMassWindowBin = 1+2*isCDFcut;
    if (iBin==3 || iBin==5 || isCDFcut < 0) mumuMassWindowBin = 0; // no cut
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaK, CosThetaL, Bmass, Q2, Mumumass, Mumumasserr, Kstarmass, Triggers),TString::Format("(%s) && (%s) && (%s) && (%s)",nTriggeredPath[2], q2range[iBin],mumuMassWindow[mumuMassWindowBin],kstarMassWindow[1]),0);
    if (data->sumEntries() == 0){
        return;
    }

    // Fitting procedure in TMinuit
    double isMigradConverge[2] = {-1,0};
    double isMinosValid = -1;
    RooAbsReal *nll = f_ext.createNLL(*data,Extended(kTRUE),Offset(kFALSE),NumCPU(4));
    RooMinuit minuit(*nll);
    printf("INFO\t\t: Start MIGRAD loop\n");
    for(int iLoop = 0; iLoop < 10; iLoop++){
        isMigradConverge[0] = minuit.migrad();
        printf("INFO\t\t: MIGRAD return code=%.0f\n",isMigradConverge[0]);
        if (isMigradConverge[0] == 0) break;
    }
    isMigradConverge[1] = minuit.save()->minNll();
    writeParam(iBin, "migrad2D", isMigradConverge);
    double isHesseValid = minuit.hesse();
    writeParam(iBin, "hesse2D", &isHesseValid, 1);
    minuit.save();
    double val[4]={0,0,0,0};
    val[0] = fl->getVal();val[1] = fl->getError();val[2]=fl->getErrorLo();val[3]=fl->getErrorHi();
    writeParam(iBin, "fl_hesse2D", val, 4);
    val[0] = afb->getVal();val[1] = afb->getError();val[2]=afb->getErrorLo();val[3]=afb->getErrorHi();
    writeParam(iBin, "afb_hesse2D",val, 4);
    val[0] = fs->getVal();val[1] = fs->getError();val[2]=fs->getErrorLo();val[3]=fs->getErrorHi();
    writeParam(iBin, "fs_hesse2D", val, 4);
    val[0] = as->getVal();val[1] = as->getError();val[2]=as->getErrorLo();val[3]=as->getErrorHi();
    writeParam(iBin, "as_hesse2D", val, 4);
    printf("INFO\t\t: Start MINOS loop\n");
    for(int iLoop = 0; iLoop < 5; iLoop++){
        isMinosValid = minuit.minos(RooArgSet(*afb,*fl));
        printf("INFO\t\t: MINOS return code=%.0f\n",isMinosValid);
        if (isMinosValid == 0) break;
    }
    writeParam(iBin, "minos2D", &isMinosValid, 1);

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");
    TLatex *t1 = new TLatex();
    t1->SetNDC();
    
    RooPlot* framecosk = CosThetaK.frame(); 
    data->plotOn(framecosk,Binning(20)); 
    f_ext.plotOn(framecosk); 
    if (false) { // Draw dashed line using generator level values
        double buffFl = fl->getVal();
        double buffAfb = afb->getVal();
        fl->setVal(genFl[iBin]); afb->setVal(genAfb[iBin]);
        f_ext.plotOn(framecosk, LineColor(2),LineWidth(2),LineStyle(2)); 
        fl->setVal(buffFl); afb->setVal(buffAfb);
    }
    framecosk->SetTitle("");
    framecosk->SetMinimum(0);
    framecosk->Draw();
    
    double fixNDC = 0;
    if (iBin > 3) fixNDC = -0.5;
    double fl_bdd[3] = {toBoundedFl(fl->getVal()),toBoundedFl(fl->getVal()+fl->getErrorHi()),toBoundedFl(fl->getVal()+fl->getErrorLo())};
    fl_bdd[1] -= fl_bdd[0];
    fl_bdd[2] -= fl_bdd[0];
    double afb_bdd[3] = {toBoundedAfb(afb->getVal(),fl->getVal()),toBoundedAfb(afb->getVal()+afb->getErrorHi(),fl->getVal()),toBoundedAfb(afb->getVal()+afb->getErrorLo(),fl->getVal())};
    afb_bdd[1] -= afb_bdd[0];
    afb_bdd[2] -= afb_bdd[0];
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    t1->DrawLatex(.35,.80+fixNDC,TString::Format("F_{L,ubd}=%.3f%+5f%+5f",fl_bdd[0],fl_bdd[1],fl_bdd[2]));
    t1->DrawLatex(.35,.74+fixNDC,TString::Format("A_{FB,ubd}=%.3f%+5f%+5f",afb_bdd[0],afb_bdd[1],afb_bdd[2]));
    c->SaveSource(TString::Format("%s/%s_cosk_bin%d.cc",plotpath.Data(),outfile,iBin));
    c->Print(TString::Format("%s/%s_cosk_bin%d.pdf",plotpath.Data(),outfile,iBin));

    // Draw projection to CosThetaL
    RooPlot* framecosl = CosThetaL.frame(); 
    data->plotOn(framecosl,Binning(20)); 
    f_ext.plotOn(framecosl); 
    if (false) { // put generator level curve for comparison
        double buffFl = fl->getVal();
        double buffAfb = afb->getVal();
        fl->setVal(genFl[iBin]); afb->setVal(genAfb[iBin]);
        f_ext.plotOn(framecosl, LineColor(2),LineWidth(2),LineStyle(2)); 
        fl->setVal(buffFl); afb->setVal(buffAfb);
    }
    framecosl->SetTitle("");
    framecosl->SetMinimum(0);
    framecosl->Draw();

    fixNDC = 0.;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    t1->DrawLatex(.35,.80+fixNDC,TString::Format("F_{L,ubd}=%.3f%+5f%+5f",fl_bdd[0],fl_bdd[1],fl_bdd[2]));
    t1->DrawLatex(.35,.74+fixNDC,TString::Format("A_{FB,ubd}=%.3f%+5f%+5f",afb_bdd[0],afb_bdd[1],afb_bdd[2]));
    c->Update();
    c->SaveSource(TString::Format("%s/%s_cosl_bin%d.cc",plotpath.Data(),outfile,iBin));
    c->Print(TString::Format("%s/%s_cosl_bin%d.pdf",plotpath.Data(),outfile,iBin));

    // Make 2-D plot
    TH1 *h1 = data->createHistogram("CosThetaL,CosThetaK", 6, 5);
    h1->SetXTitle("CosThetaL");
    h1->SetYTitle("CosThetaK");
    h1->Draw("LEGO2");
    c->Update();
    c->SaveSource(TString::Format("%s/%s_bin%d.cc",plotpath.Data(),outfile,iBin));
    c->Print(TString::Format("%s/%s_bin%d.pdf",plotpath.Data(),outfile,iBin));
    
    if (true){
        RooWorkspace *wspace = new RooWorkspace("wspace","wspace");
        nsig.setConstant(kTRUE);
        wspace->import(nsig);
        wspace->import(f_ext);
        wspace->writeToFile(TString::Format("%s/wspace_angular2D_bin%d.root",owspacepath.Data(),iBin),true);
    }

    //write output
    double output[4] = {0,0,0,0};
    output[0] = fl->getVal();
    output[1] = fl->getError();
    output[2] = fl->getErrorLo();
    output[3] = fl->getErrorHi();
    writeParam(iBin,"fl_2d",output,4);
    output[0] = afb->getVal();
    output[1] = afb->getError();
    output[2] = afb->getErrorLo();
    output[3] = afb->getErrorHi();
    writeParam(iBin,"afb_2d",output,4);
    output[0] = as->getVal();
    output[1] = as->getError();
    writeParam(iBin,"as_2d",output);
    output[0] = fs->getVal();
    output[1] = fs->getError();
    writeParam(iBin,"fs_2d",output);
    
    // clear
    delete t1;
    delete c;
    delete data;

    switchRedirectStdio("_stdio");

    return;
}//}}}

void angular2D(const char outfile[] = "angular2D", bool doFit=false) // 2D check for signal MC at RECO level
{//{{{
    int nWorkBins = 3;
    int workBins[] = {10,4,9};
    double x[nWorkBins];
    double xerr[nWorkBins];
    double yafb[nWorkBins],yerrafbLo[nWorkBins],yerrafbHi[nWorkBins],yfl[nWorkBins],yerrflLo[nWorkBins],yerrflHi[nWorkBins];
    for(int iBin = 0; iBin < nWorkBins; iBin++){
        x[iBin] = (q2rangeup[workBins[iBin]]+q2rangedn[workBins[iBin]])/2;
        xerr[iBin] = (q2rangeup[workBins[iBin]]-q2rangedn[workBins[iBin]])/2;
    }

    if (doFit){
        for(int ibin = 0; ibin < nWorkBins; ibin++){
            angular2D_bin(workBins[ibin],outfile);
        }
    }


    // Checkout input data
    for(int ibin = 0; ibin < nWorkBins; ibin++){
        yfl[ibin] = -100;
        yerrflLo[ibin] = 0;
        yerrflHi[ibin] = 0;
        yafb[ibin] = -100;
        yerrafbLo[ibin] = 0;
        yerrafbHi[ibin] = 0;
        yfl[ibin]           = toBoundedFl(readParam(workBins[ibin],"fl_2d",0));
        yerrflLo[ibin]      = fabs(toBoundedFl(readParam(workBins[ibin],"fl_2d",0)+readParam(workBins[ibin],"fl_2d",2))-yfl[ibin]);
        yerrflHi[ibin]      = fabs(toBoundedFl(readParam(workBins[ibin],"fl_2d",0)+readParam(workBins[ibin],"fl_2d",3))-yfl[ibin]);
        if (yerrflHi[ibin] == -1){
            yerrflHi[ibin] = 0;
        }
        if (yerrflLo[ibin] == -1){
            yerrflLo[ibin] = 0;
        }
        yafb[ibin]          = toBoundedAfb(readParam(workBins[ibin],"afb_2d",0),readParam(ibin,"fl_2d",0));
        yerrafbLo[ibin]     = fabs(toBoundedAfb(readParam(workBins[ibin],"afb_2d",0)+readParam(workBins[ibin],"afb_2d",2),readParam(workBins[ibin],"fl_2d",0))-yafb[ibin]);
        yerrafbHi[ibin]     = fabs(toBoundedAfb(readParam(workBins[ibin],"afb_2d",0)+readParam(workBins[ibin],"afb_2d",3),readParam(workBins[ibin],"fl_2d",0))-yafb[ibin]);
        if (yerrafbHi[ibin] == -1){
            yerrafbHi[ibin] = 0;
        }
        if (yerrafbLo[ibin] == -1){
            yerrafbLo[ibin] = 0;
        }
        printf("yafb[%d]=%6.4f + %6.4f - %6.4f\n",workBins[ibin],yafb[ibin],yerrafbHi[ibin],yerrafbLo[ibin]);
        printf("yfl [%d]=%6.4f + %6.4f - %6.4f\n",workBins[ibin],yfl[ibin],yerrflHi[ibin],yerrflLo[ibin]);
    }
    
    // Draw
    TCanvas *c = new TCanvas("c");

    TGraphAsymmErrors *g_fl  = new TGraphAsymmErrors(nWorkBins,x,yfl,xerr,xerr,yerrflLo,yerrflHi);
    g_fl->SetTitle("");
    g_fl->GetXaxis()->SetTitle("q^{2} [(GeV)^{2}]");
    g_fl->GetYaxis()->SetTitle("F_{L}");
    g_fl->GetYaxis()->SetRangeUser(0,1);
    g_fl->SetFillColor(2);
    g_fl->SetFillStyle(3001);
    g_fl->Draw("a2");
    g_fl->Draw("P");
    double work_genFl[nWorkBins];
    double work_genFlerr[nWorkBins];
    for(int ibin = 0; ibin < nWorkBins; ibin++){
        work_genFl[ibin] = genFl[workBins[ibin]];
        work_genFlerr[ibin] = genFlerr[workBins[ibin]];
    }
    TGraphAsymmErrors *gen_fl  = new TGraphAsymmErrors(nWorkBins,x,work_genFl,xerr,xerr,work_genFlerr,work_genFlerr);
    gen_fl->SetMarkerStyle(21);
    gen_fl->SetFillColor(4);
    gen_fl->SetFillStyle(3001);
    gen_fl->Draw("P2 same");
    c->SaveSource(TString::Format("%s/%s_fl.cc",plotpath.Data(),outfile));
    c->Print(TString::Format("%s/%s_fl.pdf",plotpath.Data(),outfile));
    c->Clear();

    TGraphAsymmErrors *g_afb = new TGraphAsymmErrors(nWorkBins,x,yafb,xerr,xerr,yerrafbLo,yerrafbHi);
    g_afb->SetTitle("");
    g_afb->GetXaxis()->SetTitle("q^{2} [(GeV)^{2}]");
    g_afb->GetYaxis()->SetTitle("A_{FB}");
    g_afb->GetYaxis()->SetRangeUser(-1,1);
    g_afb->SetFillColor(2);
    g_afb->SetFillStyle(3001);
    g_afb->Draw("a2");
    g_afb->Draw("P");
    double work_genAfb[nWorkBins];
    double work_genAfberr[nWorkBins];
    for(int ibin = 0; ibin < nWorkBins; ibin++){
        work_genAfb[ibin] = genAfb[workBins[ibin]];
        work_genAfberr[ibin] = genAfberr[workBins[ibin]];
    }
    TGraphAsymmErrors *gen_afb = new TGraphAsymmErrors(nWorkBins,x,work_genAfb,xerr,xerr,work_genAfberr,work_genAfberr);
    gen_afb->SetMarkerStyle(21);
    gen_afb->SetFillColor(4);
    gen_afb->SetFillStyle(3001);
    gen_afb->Draw("P2 same");
    c->SaveSource(TString::Format("%s/%s_afb.cc",plotpath.Data(),outfile));
    c->Print(TString::Format("%s/%s_afb.pdf",plotpath.Data(),outfile));
}//}}}

//_________________________________________________________________________________

void angular3D_1a_Sm(int iBin, const char outfile[] = "angular3D_1a_Sm", bool keepParam = false)
{//{{{
    // Fit to signal simulation by YsSm+YcCm to determine Sm
    ch->SetBranchStatus("*",0);
    ch->SetBranchStatus("Bmass"         , 1);
    ch->SetBranchStatus("Mumumass"      , 1);
    ch->SetBranchStatus("Mumumasserr"   , 1);
    ch->SetBranchStatus("Kstarmass"     , 1);
    ch->SetBranchStatus("Q2"            , 1);
    ch->SetBranchStatus("Triggers"      , 1);
    RooRealVar Bmass("Bmass","M_{K^{*}#Mu#Mu}",5.1,5.6);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    RooRealVar Mumumass("Mumumass","M^{#mu#mu}",0.,10.);
    RooRealVar Mumumasserr("Mumumasserr","Error of M^{#mu#mu}",0.,10.);
    RooRealVar Kstarmass("Kstarmass", "M_{K^{*}}",0,1.5);
    RooRealVar Triggers("Triggers","",0,100);
    
    // Create parameters and PDFs
        // Signal double gaussian
    RooRealVar sigGauss_mean("sigGauss_mean","M_{K*#Mu#Mu}",5.28,5.25,5.30);
    RooRealVar sigGauss1_sigma("sigGauss1_sigma","#sigma_{1}",readParam(iBin,"sigGauss1_sigma",0,0.02),.01,.05);
    RooRealVar sigGauss2_sigma("sigGauss2_sigma","#sigma_{2}",readParam(iBin,"sigGauss2_sigma",0,0.08),.05,.40);
    RooRealVar sigM_frac("sigM_frac","sigM_frac",readParam(iBin,"sigM_frac",0,0.5),0.,1.);
    
    // Create signal distribution
        // mass distro of signal
    RooGaussian f_sigMGauss1("f_sigMGauss1","f_sigMGauss1", Bmass, sigGauss_mean, sigGauss1_sigma);//double gaussian with shared mean
    RooGaussian f_sigMGauss2("f_sigMGauss2","f_sigMGauss2", Bmass, sigGauss_mean, sigGauss2_sigma);//double gaussian with shared mean
    RooAddPdf f_sigM("f_sigM","f_sigM", f_sigMGauss1, f_sigMGauss2, sigM_frac);
    
    RooRealVar nsig("nsig_MC","nsig",5E4,2000,1E8);
    RooExtendPdf f("f_sigM_ext","f",f_sigM,nsig);

    // Get data and apply unbinned fit
    int mumuMassWindowBin = 1+2*isCDFcut;
    if (iBin==3 || iBin==5) mumuMassWindowBin = 2+2*isCDFcut;
    if (isCDFcut < 0) mumuMassWindowBin =0;
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass, Mumumass, Mumumasserr, Kstarmass,Triggers),TString::Format("(%s) && (%s) && (%s)", nTriggeredPath[2],q2range[iBin],mumuMassWindow[mumuMassWindowBin], kstarMassWindow[1]),0);
    RooFitResult *f_fitresult = f.fitTo(*data,Extended(kTRUE),Save(kTRUE),Minimizer("Minuit"),Minos(kTRUE));

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");
    RooPlot* frame = Bmass.frame(); 
    data->plotOn(frame,Binning(20)); 
    f.plotOn(frame); 
    f.plotOn(frame,Components(f_sigM),LineColor(2),LineWidth(2));

    frame->SetTitle("");
    frame->SetMinimum(0);
    frame->Draw();

    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = 0;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    t1->DrawLatex(.13,.86,TString::Format("%.2f/fb",datasetLumi[1]));
    t1->DrawLatex(.13,.78,TString::Format("nsig=%5.1f#pm%5.1f",nsig.getVal(),nsig.getError()));
    c->SaveSource(TString::Format("%s/%s_bin%d.cc",plotpath.Data(),outfile,iBin));
    c->Print(TString::Format("%s/%s_bin%d.pdf",plotpath.Data(),outfile,iBin));

    // Create workspace
    if (keepParam){
        RooWorkspace *wspace = new RooWorkspace("wspace","wspace");
        nsig.setConstant(kTRUE);
        sigGauss_mean.setConstant(kTRUE);
        sigGauss1_sigma.setConstant(kTRUE);
        sigGauss2_sigma.setConstant(kTRUE);
        sigM_frac.setConstant(kTRUE);
        wspace->import(nsig);
        wspace->import(f_sigM);
        wspace->import(f);
        wspace->writeToFile(TString::Format("%s/wspace_Sm_bin%d.root",owspacepath.Data(),iBin),true);
    }

    // Prepare datacard
    double val[4]={0,0,0,0};
    writeParam(iBin, "iBin", new double((double)iBin), 1, keepParam);
    val[0]=sigGauss1_sigma.getVal();val[1]=sigGauss1_sigma.getError();val[2]=sigGauss1_sigma.getErrorLo();val[3]=sigGauss1_sigma.getErrorHi();
    writeParam(iBin, "sigGauss1_sigma", val, 4 , keepParam);
    val[0]=sigGauss2_sigma.getVal();val[1]=sigGauss2_sigma.getError();val[2]=sigGauss2_sigma.getErrorLo();val[3]=sigGauss2_sigma.getErrorHi();
    writeParam(iBin, "sigGauss2_sigma", val, 4 , keepParam);
    val[0]=sigM_frac.getVal();val[1]=sigM_frac.getError();val[2]=sigM_frac.getErrorLo();val[3]=sigM_frac.getErrorHi();
    writeParam(iBin, "sigM_frac", val, 4 , keepParam);
    
    // clear
    delete t1;
    delete c;
    delete data;

    return;
}//}}}

void angular3D_1b_YpPm(int iBin, const char outfile[] = "angular3D_1b_YpPm", bool keepParam = false)
{//{{{
    static char decmode[10];
    while(strcmp(decmode,"jpsi")*strcmp(decmode, "psi2s") != 0){
        printf("Please insert background type [ jpsi / psi2s ]:");
        scanf("%9s",decmode);
    }
    printf("Processing Bin#%d, decmode=%s\n",iBin,decmode);

    //if (iBin ==0 || iBin%2 == 1){//0,1,3,5,7
    if (iBin < 2 || iBin == 7 || iBin == 8 ){//0,1,7,{0,1}
        double val[3]={1,0,0};
        //writeParam(iBin , TString::Format("bkg%sGauss1_mean1"      , decmode) , val , 2 , keepParam);
        //writeParam(iBin , TString::Format("bkg%sGauss1_mean2"      , decmode) , val , 2 , keepParam);
        //writeParam(iBin , TString::Format("bkg%sGauss1_sigma1"     , decmode) , val , 2 , keepParam);
        //writeParam(iBin , TString::Format("bkg%sGauss1_sigma2"     , decmode) , val , 2 , keepParam);
        writeParam(iBin , TString::Format("bkg%sBifruGauss_mean"   , decmode) , val , 2 , keepParam);
        writeParam(iBin , TString::Format("bkg%sBifruGauss_Lsigma" , decmode) , val , 2 , keepParam);
        writeParam(iBin , TString::Format("bkg%sBifruGauss_Rsigma" , decmode) , val , 2 , keepParam);
        writeParam(iBin , TString::Format("bkg%sExp_index"         , decmode) , val , 2 , keepParam);
        writeParam(iBin , TString::Format("bkg%sM_frac1"           , decmode) , val , 2 , keepParam);
        writeParam(iBin , TString::Format("bkg%sGauss2_mean1"      , decmode) , val , 2 , keepParam);
        writeParam(iBin , TString::Format("bkg%sGauss2_mean2"      , decmode) , val , 2 , keepParam);
        writeParam(iBin , TString::Format("bkg%sGauss2_sigma1"     , decmode) , val , 2 , keepParam);
        writeParam(iBin , TString::Format("bkg%sGauss2_sigma2"     , decmode) , val , 2 , keepParam);
        writeParam(iBin , TString::Format("bkg%sM_frac2"           , decmode) , val , 2 , keepParam);
        writeParam(iBin , TString::Format("bkg%sM_frac12"          , decmode) , val , 2 , keepParam);
        val[0]=0;val[1]=0.1;
        writeParam(iBin , TString::Format("nbkg%sPeak"        ,decmode) , val , 2 , keepParam);
        return;
    }

    // Fit to control channel simulations by YpPm to determine Yp,Pm.
    ch->SetBranchStatus("*",0);
    ch->SetBranchStatus("Bmass"         , 1);
    ch->SetBranchStatus("Mumumass"      , 1);
    ch->SetBranchStatus("Mumumasserr"   , 1);
    ch->SetBranchStatus("Kstarmass"     , 1);
    ch->SetBranchStatus("Q2"            , 1);
    ch->SetBranchStatus("Triggers"      , 1);
    RooRealVar Bmass("Bmass","M_{K^{*}#Mu#Mu}",5.1,5.6);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    RooRealVar Mumumass("Mumumass","M^{#mu#mu}",0.,10.);
    RooRealVar Mumumasserr("Mumumasserr","Error of M^{#mu#mu}",0.,10.);
    RooRealVar Kstarmass("Kstarmass","M_{K^{*}}",0,1.5);
    RooRealVar Triggers("Triggers","",0,100);
    RooAddition Bmass_offset("Bmass_offset","Bmass_offset",RooArgSet(Bmass,RooConst(-5)));
    
    // Create peak background distribution
    RooRealVar bkgGauss1_mean1(TString::Format("bkg%sPeakGauss1_mean1",decmode)   , "M_{K*#Mu#Mu}" , readParam(iBin , TString::Format("bkg%sGauss1_mean1" ,decmode) , 0 , 5.05) , 4.00 , 5.10);
    //RooRealVar bkgGauss1_mean2(TString::Format("bkg%sPeakGauss1_mean2",decmode)   , "M_{K*#Mu#Mu}" , readParam(iBin , TString::Format("bkg%sGauss1_mean2" ,decmode) , 0 , 5.15) , 4.00 , 5.30);
    RooRealVar bkgGauss1_sigma1(TString::Format("bkg%sPeakGauss1_sigma1",decmode) , "#sigma_{11}"  , readParam(iBin , TString::Format("bkg%sGauss1_sigma1",decmode) , 0 , .03)  , .01  , .80);
    //RooRealVar bkgGauss1_sigma2(TString::Format("bkg%sPeakGauss1_sigma2",decmode) , "#sigma_{12}"  , readParam(iBin , TString::Format("bkg%sGauss1_sigma2",decmode) , 0 , .12)  , .02  , 1.2);
    RooRealVar bkgBifurGauss_mean(TString::Format("bkg%sPeakBifurGauss_mean",decmode), "M_{K*#Mu#Mu}", 5.17);
    RooRealVar bkgBifurGauss_Lsigma(TString::Format("bkg%sPeakBifurGauss_Lsigma",decmode), "M_{K*#Mu#Mu}", -0.001);
    RooRealVar bkgBifurGauss_Rsigma(TString::Format("bkg%sPeakBifurGauss_Rsigma",decmode), "M_{K*#Mu#Mu}", readParam(iBin, TString::Format("bkg%sBifurGauss_Rsigma", decmode), 0 , 0.01), 0.005, 0.04);
    RooRealVar bkgExp_index(TString::Format("bkg%sPeakExp_index",decmode) , "M_{K*#Mu#Mu}", readParam(iBin, TString::Format("bkg%sExp_index",decmode), 0, -20), -50., 0.);
    RooRealVar bkgM_frac1(TString::Format("bkg%sPeakM_frac1",decmode)             , "bkgM_frac1"   , readParam(iBin , TString::Format("bkg%sM_frac1"      ,decmode) , 0 , 1.)   , 0.   , 1.);
    RooRealVar bkgGauss2_mean1(TString::Format("bkg%sPeakGauss2_mean1",decmode)   , "M_{K*#Mu#Mu}" , readParam(iBin , TString::Format("bkg%sGauss2_mean1" ,decmode) , 0 , 5.40) , 5.30 , 5.60);
    RooRealVar bkgGauss2_mean2(TString::Format("bkg%sPeakGauss2_mean2",decmode)   , "M_{K*#Mu#Mu}" , readParam(iBin , TString::Format("bkg%sGauss2_mean2" ,decmode) , 0 , 5.40) , 5.30 , 5.60);
    RooRealVar bkgGauss2_sigma1(TString::Format("bkg%sPeakGauss2_sigma1",decmode) , "#sigma_{21}"  , readParam(iBin , TString::Format("bkg%sGauss2_sigma1",decmode) , 0 , .03)  , .005  , .10);
    RooRealVar bkgGauss2_sigma2(TString::Format("bkg%sPeakGauss2_sigma2",decmode) , "#sigma_{22}"  , readParam(iBin , TString::Format("bkg%sGauss2_sigma2",decmode) , 0 , .12)  , .10  , 1.2);
    RooRealVar bkgM_frac2(TString::Format("bkg%sPeakM_frac2",decmode)             , "bkgM_frac2"   , readParam(iBin , TString::Format("bkg%sM_frac2"      ,decmode) , 0 , 1.)   , 0.   , 1.);
    RooRealVar bkgM_frac12(TString::Format("bkg%sPeakM_frac12",decmode)           , "bkgM_frac12"  , readParam(iBin , TString::Format("bkg%sM_frac12"     ,decmode) , 0 , 1.)   , 0.   , 1.);

    RooExponential f_bkgPeakMExp(TString::Format("f_bkg%sPeakMExp",decmode) , "f_bkgPeakMExp", Bmass_offset, bkgExp_index);
    RooGaussian f_bkgPeakMGauss11(TString::Format("f_bkg%sPeakMGauss11",decmode) ,"f_bkgPeakMGauss11", Bmass, bkgGauss1_mean1, bkgGauss1_sigma1);
    //RooGaussian f_bkgPeakMGauss12(TString::Format("f_bkg%sPeakMGauss12",decmode) ,"f_bkgPeakMGauss12", Bmass, bkgGauss1_mean2, bkgGauss1_sigma2);
    RooBifurGauss f_bkgPeakMBifurGauss(TString::Format("f_bkg%sPeakMBifurGauss",decmode) , "f_bkgPeakMBifurGauss", Bmass, bkgBifurGauss_mean, bkgBifurGauss_Lsigma, bkgBifurGauss_Rsigma);

    RooGaussian f_bkgPeakMGauss21(TString::Format("f_bkg%sPeakMGauss21",decmode) ,"f_bkgPeakMGauss21", Bmass, bkgGauss2_mean1, bkgGauss2_sigma1);
    RooGaussian f_bkgPeakMGauss22(TString::Format("f_bkg%sPeakMGauss22",decmode) ,"f_bkgPeakMGauss22", Bmass, bkgGauss2_mean2, bkgGauss2_sigma2);

    RooAddPdf f_bkgPeakM1(TString::Format("f_bkg%sPeakM1",decmode) , "f_bkgPeakM1", RooArgList(f_bkgPeakMBifurGauss, f_bkgPeakMExp),RooArgList(bkgM_frac1));
    RooAddPdf f_bkgPeakM2(TString::Format("f_bkg%sPeakM2",decmode)   , "f_bkgPeakM2"  , RooArgList(f_bkgPeakMGauss21   , f_bkgPeakMGauss22) , bkgM_frac2);
    RooAddPdf f_bkgPeakM12(TString::Format("f_bkg%sPeakM12",decmode) , "f_bkgPeakM12" , RooArgList(f_bkgPeakMExp       , f_bkgPeakMGauss21) , bkgM_frac12);
    RooRealVar nbkgPeak(TString::Format("nbkg%sPeak_MC",decmode) ,"nbkgPeak",2E1,0,1E7);
    RooExtendPdf *f = 0;
    f = new RooExtendPdf(TString::Format("f_bkg%sPeak_ext",decmode) , "f", f_bkgPeakM12, nbkgPeak);
        // Fit with the bifurcated term, but drop it in the final 3-D fit.
    RooRealVar nbkgPeakMBifur(TString::Format("nbkg%sPeakM_Bifur",decmode) , "nbkgPeakMBifur",  2E1, 0. , 1E7);
    RooExtendPdf f_bkgPeakMBifur(TString::Format("f_bkg%sPeakMBifur",decmode), "f_bkgPeakMBifur", f_bkgPeakMBifurGauss, nbkgPeakMBifur);
    RooAddPdf f_bkgPeakM12B(TString::Format("f_bkg%sPeakM12B",decmode), "f_bkgPeakM12B", RooArgList(f_bkgPeakMBifur,*f));

    RooAddition nbkgPeakFull("nbkgPeakFull" ,"nbkgPeakFull", RooArgList(nbkgPeak, nbkgPeakMBifur));

    bkgM_frac2.setVal(1);
    bkgM_frac2.setConstant(kTRUE);
    bkgGauss2_mean2.setConstant(kTRUE);
    bkgGauss2_sigma2.setConstant(kTRUE);
    nbkgPeakMBifur.setVal(0);
    nbkgPeakMBifur.setConstant(kTRUE);
    bkgBifurGauss_Rsigma.setConstant(kTRUE);
    if (strcmp(decmode,"jpsi")==0){
        switch (iBin) {
            case 11:
            case 10:
                bkgM_frac12.setVal(1);
                bkgM_frac12.setConstant(kTRUE);
                bkgGauss2_mean1.setConstant(kTRUE);
                bkgGauss2_sigma1.setConstant(kTRUE);
                bkgM_frac1.setVal(0.1);
                break;
            case 3:
                bkgM_frac1.setVal(1);
                bkgM_frac1.setConstant(kTRUE);
                bkgExp_index.setConstant(kTRUE);
                break;
            case 4:
                bkgM_frac12.setVal(0.);
                bkgM_frac12.setConstant(kTRUE);
                bkgM_frac1.setConstant(kTRUE);
                bkgExp_index.setConstant(kTRUE);
                bkgBifurGauss_Rsigma.setConstant(kTRUE);
                break;
            case 12:// No events in bin12 
            default:
                return;
        }
    }else{
        switch (iBin) {
            case 4:
            case 11://summaryBin[0]
                bkgM_frac12.setVal(1);
                bkgM_frac12.setConstant(kTRUE);
                bkgGauss2_mean1.setConstant(kTRUE);
                bkgGauss2_sigma1.setConstant(kTRUE);
                break;
            case 5:
                bkgM_frac12.setVal(1);
                bkgM_frac12.setConstant(kTRUE);
                bkgGauss2_mean1.setConstant(kTRUE);
                bkgGauss2_sigma1.setConstant(kTRUE);
                bkgM_frac1.setVal(1);
                bkgM_frac1.setConstant(kTRUE);
                bkgExp_index.setConstant(kTRUE);
                break;
            case 6:
            case 9:
                bkgM_frac12.setVal(0.);
                bkgM_frac12.setConstant(kTRUE);
                bkgM_frac1.setConstant(kTRUE);
                bkgExp_index.setConstant(kTRUE);
                bkgBifurGauss_Rsigma.setConstant(kTRUE);
                break;
            default:
                return;
        }
    }

    // Get data and apply unbinned fit
    int mumuMassWindowBin = 1+2*isCDFcut;
    if (iBin==3 || iBin==5) mumuMassWindowBin = 2+2*isCDFcut;
    if (isCDFcut < 0) mumuMassWindowBin =0;
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, Bmass, Mumumass, Mumumasserr, Kstarmass, Triggers),TString::Format("(%s) && (%s) && (%s) && (%s)", nTriggeredPath[2],q2range[iBin],mumuMassWindow[mumuMassWindowBin], kstarMassWindow[1]),0);
    RooFitResult *f_fitresult = f_bkgPeakM12B.fitTo(*data,Save(kTRUE),Minimizer("Minuit"),Extended(),Minos(kTRUE));

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");
    RooPlot* frame = Bmass.frame(); 
    data->plotOn(frame,Binning(20)); 
    f_bkgPeakM12B.plotOn(frame,LineColor(1)); 
    f_bkgPeakM12B.plotOn(frame,Components(*f),LineColor(4),LineWidth(2),LineStyle(2)); 
    //f_bkgPeakM12B.plotOn(frame,Components(f_bkgPeakMBifur),LineColor(2),LineWidth(2),LineStyle(2)); 

    frame->SetTitle("");
    frame->SetMinimum(0);
    frame->Draw();

    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = 0;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    //t1->DrawLatex(.13,.78,TString::Format("nsig=%5.2f#pm%5.2f",nbkgPeak.getVal(),nbkgPeak.getError()));
    c->SaveSource(TString::Format("%s/%s_%s_bin%d.cc",plotpath.Data(),outfile,decmode,iBin));
    c->Print(TString::Format("%s/%s_%s_bin%d.pdf",plotpath.Data(),outfile,decmode,iBin));

    // Write to workspace
    if (keepParam){
        RooWorkspace *wspace = new RooWorkspace("wspace","wspace");
        bkgBifurGauss_mean.setConstant(kTRUE);
        bkgBifurGauss_Lsigma.setConstant(kTRUE);
        bkgBifurGauss_Rsigma.setConstant(kTRUE);
        bkgExp_index.setConstant(kTRUE);
        bkgM_frac1.setConstant(kTRUE);
        bkgGauss2_mean1.setConstant(kTRUE);
        bkgGauss2_mean2.setConstant(kTRUE);
        bkgGauss2_sigma1.setConstant(kTRUE);
        bkgGauss2_sigma2.setConstant(kTRUE);
        bkgM_frac2.setConstant(kTRUE);
        bkgM_frac12.setConstant(kTRUE);
        nbkgPeakMBifur.setConstant(kTRUE);
        nbkgPeak.setConstant(kTRUE);
        wspace->import(nbkgPeakFull);
        wspace->import(f_bkgPeakM12B);
        wspace->writeToFile(TString::Format("%s/wspace_YpPm_%s_bin%d.root",owspacepath.Data(),decmode,iBin),true);
    }

    // Write datacard
    double val[3]={0,0,0};
    //val[0] = bkgGauss1_mean1.getVal();val[1] = bkgGauss1_mean1.getError();
    //writeParam(iBin , TString::Format("bkg%sGauss1_mean1"  , decmode) , val , 2 , keepParam);
    //val[0] = bkgGauss1_mean2.getVal();val[1] = bkgGauss1_mean2.getError();
    //writeParam(iBin , TString::Format("bkg%sGauss1_mean2"  , decmode) , val , 2 , keepParam);
    //val[0] = bkgGauss1_sigma1.getVal();val[1] = bkgGauss1_sigma1.getError();
    //writeParam(iBin , TString::Format("bkg%sGauss1_sigma1" , decmode) , val , 2 , keepParam);
    //val[0] = bkgGauss1_sigma2.getVal();val[1] = bkgGauss1_sigma2.getError();
    //writeParam(iBin , TString::Format("bkg%sGauss1_sigma2" , decmode) , val , 2 , keepParam);
    val[0] = bkgM_frac1.getVal();val[1] = bkgM_frac1.getError();
    writeParam(iBin , TString::Format("bkg%sM_frac1"       , decmode) , val , 2 , keepParam);
    val[0] = bkgGauss2_mean1.getVal();val[1] = bkgGauss2_mean1.getError();
    writeParam(iBin , TString::Format("bkg%sGauss2_mean1"  , decmode) , val , 2 , keepParam);
    val[0] = bkgGauss2_mean2.getVal();val[1] = bkgGauss2_mean2.getError();
    writeParam(iBin , TString::Format("bkg%sGauss2_mean2"  , decmode) , val , 2 , keepParam);
    val[0] = bkgGauss2_sigma1.getVal();val[1] = bkgGauss2_sigma1.getError();
    writeParam(iBin , TString::Format("bkg%sGauss2_sigma1" , decmode) , val , 2 , keepParam);
    val[0] = bkgGauss2_sigma2.getVal();val[1] = bkgGauss2_sigma2.getError();
    writeParam(iBin , TString::Format("bkg%sGauss2_sigma2" , decmode) , val , 2 , keepParam);
    val[0] = bkgM_frac2.getVal();val[1] = bkgM_frac2.getError();
    writeParam(iBin , TString::Format("bkg%sM_frac2"       , decmode) , val , 2 , keepParam);
    val[0] = bkgM_frac12.getVal();val[1] = bkgM_frac12.getError();
    writeParam(iBin , TString::Format("bkg%sM_frac12"      , decmode) , val , 2 , keepParam);
    val[0] = nbkgPeak.getVal();val[1] = nbkgPeak.getError();
    writeParam(iBin, TString::Format("nbkg%sPeak",decmode), val , 2 , keepParam);
    
    // clear
    delete t1;
    delete c;
    delete data;

    printf("INFO\t\t: End of angular3D_1b_YpPm(%d)\n",iBin);
}//}}}

void angular3D_2a_PkPl(int iBin, const char outfile[] = "angular3D_2a_PkPl", bool keepParam = false)
{//{{{
    static char decmode[20];
    while(strcmp(decmode,"jpsi")*strcmp(decmode, "psi2s") != 0){
        printf("Please insert background type [ jpsi / psi2s ]:");
        scanf("%19s",decmode);
    }

    // Gaussian constraint on yields and mass is needed.
    //if (iBin ==0 || iBin%2 == 1){//0,1,3,5,7
    if (iBin < 2 || iBin == 7 || iBin == 8 ){//0,1,7,{0,1}
        // Pm is flat(and the yield is 0) for bins other than 2,4,6,10
        double val[3]={0,0,0};
        writeParam(iBin, TString::Format("bkg%sPeak_c1",decmode), val , 2 , keepParam);
        writeParam(iBin, TString::Format("bkg%sPeak_c2",decmode), val , 2 , keepParam);
        writeParam(iBin, TString::Format("bkg%sPeak_c3",decmode), val , 2 , keepParam);
        writeParam(iBin, TString::Format("bkg%sPeak_c4",decmode), val , 2 , keepParam);
        writeParam(iBin, TString::Format("bkg%sPeak_c5",decmode), val , 2 , keepParam);
        writeParam(iBin, TString::Format("bkg%sPeak_c6",decmode), val , 2 , keepParam);
        writeParam(iBin, TString::Format("bkg%sPeak_c7",decmode), val , 2 , keepParam);
        writeParam(iBin, TString::Format("bkg%sPeak_c8",decmode), val , 2 , keepParam);
        return;
    }

    ch->SetBranchStatus("*",0);
    ch->SetBranchStatus("Bmass"         , 1);
    ch->SetBranchStatus("Mumumass"      , 1);
    ch->SetBranchStatus("Mumumasserr"   , 1);
    ch->SetBranchStatus("Kstarmass"     , 1);
    ch->SetBranchStatus("CosTheta*"     , 1);
    ch->SetBranchStatus("Q2"            , 1);
    ch->SetBranchStatus("Triggers"      , 1);
    RooRealVar Bmass("Bmass","M_{K^{*}#Mu#Mu}",5.1,5.6);
    RooRealVar CosThetaK("CosThetaK"     , "cos#theta_{K}"       , -1. , 1.   ) ;
    RooRealVar CosThetaL("CosThetaL"     , "cos#theta_{L}"       , -1. , 1.   ) ;
    RooRealVar Mumumass("Mumumass","M^{#mu#mu}",0.,10.);
    RooRealVar Mumumasserr("Mumumasserr","Error of M^{#mu#mu}",0.,10.);
    RooRealVar Kstarmass("Kstarmass",   "M_{K^{*}}", 0,1.5);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    RooRealVar Triggers("Triggers","",0,100);
    
    RooRealVar bkgPeak_c1(TString::Format("bkg%sPeak_c1",decmode),"bkgPeak_c1",readParam(iBin,TString::Format("bkg%sPeak_c1",decmode),0,0.,0.),-10,10);
    RooRealVar bkgPeak_c2(TString::Format("bkg%sPeak_c2",decmode),"bkgPeak_c2",readParam(iBin,TString::Format("bkg%sPeak_c2",decmode),0,0.,0.),-10,10);
    RooRealVar bkgPeak_c3(TString::Format("bkg%sPeak_c3",decmode),"bkgPeak_c3",readParam(iBin,TString::Format("bkg%sPeak_c3",decmode),0,0.,0.),-10,10);
    RooRealVar bkgPeak_c4(TString::Format("bkg%sPeak_c4",decmode),"bkgPeak_c4",readParam(iBin,TString::Format("bkg%sPeak_c4",decmode),0,0.,0.),-10,10);
    RooRealVar bkgPeak_c5(TString::Format("bkg%sPeak_c5",decmode),"bkgPeak_c5",readParam(iBin,TString::Format("bkg%sPeak_c5",decmode),0,0.,0.),-10,10);
    RooRealVar bkgPeak_c6(TString::Format("bkg%sPeak_c6",decmode),"bkgPeak_c6",readParam(iBin,TString::Format("bkg%sPeak_c6",decmode),0,0.,0.),-10,10);
    RooRealVar bkgPeak_c7(TString::Format("bkg%sPeak_c7",decmode),"bkgPeak_c7",readParam(iBin,TString::Format("bkg%sPeak_c7",decmode),0,0.,0.),-10,10);
    RooRealVar bkgPeak_c8(TString::Format("bkg%sPeak_c8",decmode),"bkgPeak_c8",readParam(iBin,TString::Format("bkg%sPeak_c8",decmode),0,0.,0.),-10,10);
    RooArgSet f_bkgPeakA_argset, f_bkgPeakK_argset, f_bkgPeakL_argset;
    f_bkgPeakA_argset.add(RooArgSet(CosThetaL));
    f_bkgPeakA_argset.add(RooArgSet(CosThetaK));
    f_bkgPeakK_argset.add(RooArgSet(CosThetaK,bkgPeak_c1,bkgPeak_c2,bkgPeak_c3,bkgPeak_c4));
    f_bkgPeakL_argset.add(RooArgSet(CosThetaL,bkgPeak_c5,bkgPeak_c6,bkgPeak_c7,bkgPeak_c8));
    f_bkgPeakA_argset.add(RooArgSet(bkgPeak_c1,bkgPeak_c2,bkgPeak_c3,bkgPeak_c4,bkgPeak_c5,bkgPeak_c6,bkgPeak_c7,bkgPeak_c8));
    // Base function+correction, this works without CDF/LHCb cut.
    //TString f_bkgPeakA_format = "(1+\
    //                            bkgPeak_c1*CosThetaK+\
    //                            bkgPeak_c2*CosThetaK**2+\
    //                            bkgPeak_c3*CosThetaL**2+\
    //                            bkgPeak_c4*CosThetaL**4+\
    //                            bkgPeak_c5*CosThetaL**2*CosThetaK+\
    //                            bkgPeak_c6*CosThetaL**4*CosThetaK+\
    //                            bkgPeak_c7*CosThetaL**2*CosThetaK**2+\
    //                            bkgPeak_c8*CosThetaL**4*CosThetaK**2)*\
    //                            exp(-(CosThetaK+0.5*CosThetaL)**4-(CosThetaK-0.5*CosThetaL)**4+1.5*(CosThetaK+0.5*CosThetaL)**2+1.5*(CosThetaK-0.5*CosThetaL)**2)";
    //TString f_bkgPeakK_format = "1+bkgPeak_c1*CosThetaK+bkgPeak_c2*CosThetaK**2+bkgPeak_c3*CosThetaK**3+bkgPeak_c4*CosThetaK**4";
    //TString f_bkgPeakL_format = "1+bkgPeak_c5*CosThetaL+bkgPeak_c6*CosThetaL**2+bkgPeak_c7*CosThetaL**3+bkgPeak_c8*CosThetaL**4";
    TString f_bkgPeakK_format = TString::Format("1+bkg%sPeak_c1*CosThetaK+bkg%sPeak_c2*CosThetaK**2+bkg%sPeak_c3*CosThetaK**3+bkg%sPeak_c4*CosThetaK**4",decmode,decmode,decmode,decmode);
    TString f_bkgPeakL_format = TString::Format("1+bkg%sPeak_c5*CosThetaL+bkg%sPeak_c6*CosThetaL**2+bkg%sPeak_c7*CosThetaL**3+bkg%sPeak_c8*CosThetaL**4",decmode,decmode,decmode,decmode);
    TString f_bkgPeakA_format = TString::Format("(%s)*(%s)",f_bkgPeakK_format.Data(),f_bkgPeakL_format.Data());
    bkgPeak_c3.setVal(0.);
    bkgPeak_c4.setVal(0.);
    bkgPeak_c7.setVal(0.);
    bkgPeak_c8.setVal(0.);
    bkgPeak_c3.setConstant(kTRUE);
    bkgPeak_c4.setConstant(kTRUE);
    bkgPeak_c7.setConstant(kTRUE);
    bkgPeak_c8.setConstant(kTRUE);

    if (strcmp(decmode,"jpsi")==0){
        switch (iBin) {
            case 5:
            case 6:
            case 4:
            case 9:
            case 12:
                return;
            default:
                break;
        }
    }else{
        // Reduce order due to limited statistics
        bkgPeak_c3.setVal(0.);
        bkgPeak_c4.setVal(0.);
        bkgPeak_c3.setConstant(kTRUE);
        bkgPeak_c4.setConstant(kTRUE);
        switch (iBin) {
            case 2:
            case 3:
            case 9:
            case 10:
            case 12:
                return;
            case 6:
                // Linear model for limited statistics
                bkgPeak_c6.setVal(0.);
                bkgPeak_c6.setConstant(kTRUE);
                break;
            default:
                break;
        }
    }
    RooRealVar nbkgPeak(TString::Format("nbkg%sPeak",decmode),"nbkg%sPeak",readParam(iBin,TString::Format("nbkg%sPeak",decmode),0,10),0.1,1E4);
    RooGenericPdf f_bkgPeakK(TString::Format("f_bkg%sPeakK",decmode), f_bkgPeakK_format, f_bkgPeakK_argset);
    RooGenericPdf f_bkgPeakL(TString::Format("f_bkg%sPeakL",decmode), f_bkgPeakL_format, f_bkgPeakL_argset);
    RooGenericPdf f_bkgPeakA(TString::Format("f_bkg%sPeakA",decmode), f_bkgPeakA_format, f_bkgPeakA_argset);
    RooExtendPdf f_bkgPeakK_ext(TString::Format("f_bkg%sPeakK_ext",decmode),"f_bkgPeakK_ext",f_bkgPeakK, nbkgPeak);
    RooExtendPdf f_bkgPeakL_ext(TString::Format("f_bkg%sPeakL_ext",decmode),"f_bkgPeakL_ext",f_bkgPeakL, nbkgPeak);
    RooExtendPdf f_bkgPeakA_ext(TString::Format("f_bkg%sPeakA_ext",decmode),"f_bkgPeakA_ext",f_bkgPeakA, nbkgPeak);

    // Get data
    int mumuMassWindowBin = 1+2*isCDFcut;
    if (iBin==3 || iBin==5) mumuMassWindowBin = 2+2*isCDFcut;
    if (isCDFcut < 0) mumuMassWindowBin = 0;
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(Q2, CosThetaK, Bmass, CosThetaL ,Mumumass, Mumumasserr, Kstarmass, Triggers),TString::Format("(%s) && (%s) && (%s) && (%s)",nTriggeredPath[2], q2range[iBin],mumuMassWindow[mumuMassWindowBin], kstarMassWindow[1]),0);
    f_bkgPeakL.fitTo(*data,Save(kTRUE),Minimizer("Minuit"),Minos(kTRUE));
    f_bkgPeakK.fitTo(*data,Save(kTRUE),Minimizer("Minuit"),Minos(kTRUE));
    bkgPeak_c1.setConstant(kTRUE);
    bkgPeak_c2.setConstant(kTRUE);
    bkgPeak_c3.setConstant(kTRUE);
    bkgPeak_c4.setConstant(kTRUE);
    bkgPeak_c5.setConstant(kTRUE);
    bkgPeak_c6.setConstant(kTRUE);
    bkgPeak_c7.setConstant(kTRUE);
    bkgPeak_c8.setConstant(kTRUE);
    RooFitResult *f_fitresult = f_bkgPeakA_ext.fitTo(*data,Save(kTRUE),Extended(kTRUE),Minimizer("Minuit"),Minos(kTRUE));
    
    // Draw CosThetaK
    TCanvas* c = new TCanvas("c");
    RooPlot* framecosk = CosThetaK.frame();
    data->plotOn(framecosk,Binning(20));
    f_bkgPeakA_ext.plotOn(framecosk);
    framecosk->SetTitle("");
    framecosk->SetMinimum(0);
    framecosk->Draw();
    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = 0;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->SaveSource(TString::Format("%s/%s_%s_cosk_bin%d.cc",plotpath.Data(),outfile,decmode,iBin));
    c->Print(TString::Format("%s/%s_%s_cosk_bin%d.pdf",plotpath.Data(),outfile,decmode,iBin));
    
    // Draw CosThetaL
    RooPlot* framecosl = CosThetaL.frame();
    data->plotOn(framecosl,Binning(20));
    f_bkgPeakA_ext.plotOn(framecosl);
    framecosl->SetTitle("");
    framecosl->SetMinimum(0);
    framecosl->Draw();
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->SaveSource(TString::Format("%s/%s_%s_cosl_bin%d.cc",plotpath.Data(),outfile,decmode,iBin));
    c->Print(TString::Format("%s/%s_%s_cosl_bin%d.pdf",plotpath.Data(),outfile,decmode,iBin));
    
    // Make 2-D plot
    TH1 *h1 = data->createHistogram("CosThetaK,CosThetaL", 5, 6);
    h1->SetMinimum(0);
    h1->Draw("LEGO2");
    TH1 *h1_fit_fine = f_bkgPeakA_ext.createHistogram("CosThetaK,CosThetaL", 100, 120);
    h1_fit_fine->SetMinimum(0);
    h1_fit_fine->Draw("SURF SAME");
    c->Update();
    c->SaveSource(TString::Format("%s/%s_%s_bin%d.cc",plotpath.Data(),outfile,decmode,iBin));
    c->Print(TString::Format("%s/%s_%s_bin%d.pdf",plotpath.Data(),outfile,decmode,iBin));

    TH1 *h1_fit = f_bkgPeakA_ext.createHistogram("CosThetaK,CosThetaL", 5, 6);
    TH2F *h2_compFit = new TH2F("h2_compFit","Yields_{MC}/Yields_{Fit}",5,-1.,1.,6,-1.,1.);
    h2_compFit->SetXTitle("CosThetaK");
    h2_compFit->SetYTitle("CosThetaL");
    for (int i = 1; i <= h2_compFit->GetNbinsX(); i++) {
        for (int j = 1; j <= h2_compFit->GetNbinsY(); j++) {
            h2_compFit->SetBinContent(i,j,h1->GetBinContent(i,j)/h1_fit->GetBinContent(i,j)*h2_compFit->GetNbinsX()*h2_compFit->GetNbinsX()/4.);
        }
    }
    h2_compFit->SetStats(kFALSE);
    h2_compFit->SetMinimum(0.);
    h2_compFit->Draw("LEGO2");
    c->Update();
    c->SaveSource(TString::Format("%s/%s_%s_compFit_2D_bin%d.cc",plotpath.Data(),outfile,decmode,iBin));
    c->Print(TString::Format("%s/%s_%s_compFit_2D_bin%d.pdf",plotpath.Data(),outfile,decmode,iBin));
    h2_compFit->Draw("COL TEXT");
    c->Update();
    c->SaveSource(TString::Format("%s/%s_%s_compFit_2D_TEXT_bin%d.cc",plotpath.Data(),outfile,decmode,iBin));
    c->Print(TString::Format("%s/%s_%s_compFit_2D_TEXT_bin%d.pdf",plotpath.Data(),outfile,decmode,iBin));
    
    // Write to workspace
    if (keepParam){
        //RooWorkspace *wspace = new RooWorkspace("wspace","wspace");
        TFile *f_wspace = new TFile(TString::Format("%s/wspace_PkPl_%s_bin%d.root",owspacepath.Data(),decmode,iBin),"UPDATE");
        RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
        if (!wspace){
            wspace = new RooWorkspace("wspace","wspace");
        }
        bkgPeak_c1.setConstant(kTRUE);
        bkgPeak_c2.setConstant(kTRUE);
        bkgPeak_c3.setConstant(kTRUE);
        bkgPeak_c4.setConstant(kTRUE);
        bkgPeak_c5.setConstant(kTRUE);
        bkgPeak_c6.setConstant(kTRUE);
        bkgPeak_c7.setConstant(kTRUE);
        bkgPeak_c8.setConstant(kTRUE);
        nbkgPeak.setConstant(kTRUE);
        wspace->import(f_bkgPeakK);
        wspace->import(f_bkgPeakL);
        wspace->import(f_bkgPeakA);
        wspace->writeToFile(TString::Format("%s/wspace_PkPl_%s_bin%d.root",owspacepath.Data(),decmode,iBin),true);
    }

    // write datacard
    double val[3] = {0,0,0};
    val[0] = bkgPeak_c1.getVal();val[1] = bkgPeak_c1.getError();
    writeParam(iBin, TString::Format("bkg%sPeak_c1", decmode), val , 2 , keepParam);
    val[0] = bkgPeak_c2.getVal();val[1] = bkgPeak_c2.getError();
    writeParam(iBin, TString::Format("bkg%sPeak_c2", decmode), val , 2 , keepParam);
    val[0] = bkgPeak_c3.getVal();val[1] = bkgPeak_c3.getError();
    writeParam(iBin, TString::Format("bkg%sPeak_c3", decmode), val , 2 , keepParam);
    val[0] = bkgPeak_c4.getVal();val[1] = bkgPeak_c4.getError();
    writeParam(iBin, TString::Format("bkg%sPeak_c4", decmode), val , 2 , keepParam);
    val[0] = bkgPeak_c5.getVal();val[1] = bkgPeak_c5.getError();
    writeParam(iBin, TString::Format("bkg%sPeak_c5", decmode), val , 2 , keepParam);
    val[0] = bkgPeak_c6.getVal();val[1] = bkgPeak_c6.getError();
    writeParam(iBin, TString::Format("bkg%sPeak_c6", decmode), val , 2 , keepParam);
    val[0] = bkgPeak_c7.getVal();val[1] = bkgPeak_c7.getError();
    writeParam(iBin, TString::Format("bkg%sPeak_c7", decmode), val , 2 , keepParam);
    val[0] = bkgPeak_c8.getVal();val[1] = bkgPeak_c8.getError();
    writeParam(iBin, TString::Format("bkg%sPeak_c8", decmode), val , 2 , keepParam);

    // clear
    delete t1;
    delete c;
    delete data;
    
}//}}}

void angular3D_prior(int iBin, const char outfile[] = "angular3D_prior", bool keepParam = false)// DEPRECATED, single PDF for sidebands
{//{{{
    // Fit to signal simulation by YsSm+YcCm to determine Sm
    ch->SetBranchStatus("*",0);
    ch->SetBranchStatus("Bmass"         , 1);
    ch->SetBranchStatus("Mumumass"      , 1);
    ch->SetBranchStatus("Mumumasserr"   , 1);
    ch->SetBranchStatus("Kstarmass"     , 1);
    ch->SetBranchStatus("CosTheta*"     , 1);
    ch->SetBranchStatus("Q2"            , 1);
    ch->SetBranchStatus("Triggers"      , 1);
    RooRealVar Bmass("Bmass","M_{K^{*}#Mu#Mu}",4.5,6.5);// Use wider range for sideband
    RooRealVar CosThetaK("CosThetaK"     , "cos#theta_{K}"       , -1. , 1.   ) ;
    RooRealVar CosThetaL("CosThetaL"     , "cos#theta_{L}"       , -1. , 1.   ) ;
    RooRealVar Mumumass("Mumumass","M_{#mu#mu}",0.,10.);
    RooRealVar Mumumasserr("Mumumasserr","Error of M_{#mu#mu}",0.,10.);
    RooRealVar Kstarmass("Kstarmass","M_{K*}",0.,1.5);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    RooRealVar Triggers("Triggers","",0,100);

    int mumuMassWindowBin = 1+2*isCDFcut;
    if (iBin==3 || iBin==5 || isCDFcut < 0 ) mumuMassWindowBin = 0;
    RooDataSet *dataM = new RooDataSet("dataM","data",ch,RooArgSet(Q2, Bmass, CosThetaK, CosThetaL, Mumumass, Mumumasserr, Kstarmass, Triggers),TString::Format("(%s) && (%s) && (%s) && (%s)",nTriggeredPath[2],q2range[iBin],mumuMassWindow[mumuMassWindowBin],kstarMassWindow[2]),0);
    RooDataSet *dataA = new RooDataSet("dataA","data",ch,RooArgSet(Q2, Bmass, CosThetaK, CosThetaL, Mumumass, Mumumasserr, Kstarmass, Triggers),TString::Format("(%s) && (%s) && (%s) && (%s) && (%s)",nTriggeredPath[2],q2range[iBin],mumuMassWindow[mumuMassWindowBin],kstarMassWindow[1],"Bmass > 5.6 || Bmass < 5.0"),0);
    RooDataSet *dataAHi = new RooDataSet("dataA","data",ch,RooArgSet(Q2, Bmass, CosThetaK, CosThetaL, Mumumass, Mumumasserr, Kstarmass, Triggers),TString::Format("(%s) && (%s) && (%s) && (%s) && (%s)",nTriggeredPath[2],q2range[iBin],mumuMassWindow[mumuMassWindowBin],kstarMassWindow[1],"Bmass > 5.6"),0);
    RooDataSet *dataALo = new RooDataSet("dataA","data",ch,RooArgSet(Q2, Bmass, CosThetaK, CosThetaL, Mumumass, Mumumasserr, Kstarmass, Triggers),TString::Format("(%s) && (%s) && (%s) && (%s) && (%s)",nTriggeredPath[2],q2range[iBin],mumuMassWindow[mumuMassWindowBin],kstarMassWindow[1],"Bmass < 5.0"),0);
    if (dataM->sumEntries() == 0 || dataA->sumEntries() == 0){
        printf("ERROR   : Total number of entries is 0. Exit!\n");
        return;
    }

    // Create peak background distribution(jpsi/psi2s)
    TFile *f_wspace_jpsi_A = new TFile(TString::Format("%s/wspace_PkPl_jpsi_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace_jpsi_A = (RooWorkspace*)f_wspace_jpsi_A->Get("wspace");
    RooGenericPdf *f_bkgjpsiPeakA = 0;
    if (wspace_jpsi_A){
        f_bkgjpsiPeakA = (RooGenericPdf*)wspace_jpsi_A->pdf("f_bkgjpsiPeakA");
    }

    TFile *f_wspace_jpsi_M = new TFile(TString::Format("%s/wspace_YpPm_jpsi_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace_jpsi_M = (RooWorkspace*)f_wspace_jpsi_M->Get("wspace");
    RooRealVar *nbkgjpsiPeak_MC = 0;
    RooAddPdf *f_bkgjpsiPeakM12 = 0;
    if (wspace_jpsi_M){
        nbkgjpsiPeak_MC = (RooRealVar*)wspace_jpsi_M->var("nbkgjpsiPeak_MC");
        f_bkgjpsiPeakM12 = (RooAddPdf*)wspace_jpsi_M->pdf("f_bkgjpsiPeakM12");
    }

    TFile *f_wspace_psi2s_A = new TFile(TString::Format("%s/wspace_PkPl_psi2s_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace_psi2s_A = (RooWorkspace*)f_wspace_psi2s_A->Get("wspace");
    RooGenericPdf *f_bkgpsi2sPeakA = 0;
    if (wspace_psi2s_A){
        f_bkgpsi2sPeakA = (RooGenericPdf*)wspace_psi2s_A->pdf("f_bkgpsi2sPeakA");
    }

    TFile *f_wspace_psi2s_M = new TFile(TString::Format("%s/wspace_YpPm_psi2s_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace_psi2s_M = (RooWorkspace*)f_wspace_psi2s_M->Get("wspace");
    RooRealVar *nbkgpsi2sPeak_MC = 0;
    RooAddPdf *f_bkgpsi2sPeakM12 = 0;
    if (wspace_psi2s_M){
        nbkgpsi2sPeak_MC  = (RooRealVar*)wspace_psi2s_M->var("nbkgpsi2sPeak_MC");
        f_bkgpsi2sPeakM12 = (RooAddPdf*)wspace_psi2s_M->pdf("f_bkgpsi2sPeakM12");
    }
    printf("INFO: f_bkgPeak prepared.\n");
    
    // Create combinatorial background distribution
    TString f_bkgCombM_format = "1+bkgCombM_c1*(Bmass-5)+bkgCombM_c2*exp(bkgCombM_c3*(Bmass-bkgCombM_c4))";
    RooRealVar bkgCombM_c1("bkgCombM_c1","c1",readParam(iBin,"bkgCombM_c1",0, 0),-20,20);
    RooRealVar bkgCombM_c2("bkgCombM_c2","c2",readParam(iBin,"bkgCombM_c2",0,10),-20,20);
    RooRealVar bkgCombM_c3("bkgCombM_c3","c3",readParam(iBin,"bkgCombM_c3",0,-3),-10,10);
    RooRealVar bkgCombM_c4("bkgCombM_c4","c4",readParam(iBin,"bkgCombM_c4",0, 5),-10,10);
    RooArgSet f_bkgCombM_argset;
    f_bkgCombM_argset.add(RooArgSet(bkgCombM_c1,bkgCombM_c2,bkgCombM_c3,bkgCombM_c4));
    f_bkgCombM_argset.add(RooArgSet(Bmass));
    RooGenericPdf f_bkgCombM("f_bkgCombM","f_bkgCombM",f_bkgCombM_format,f_bkgCombM_argset);

    TString f_bkgCombL_format;
    RooRealVar bkgCombL_c1("bkgCombL_c1","c1",readParam(iBin,"bkgCombL_c1",0,0.1),-10,10);
    RooRealVar bkgCombL_c2("bkgCombL_c2","c2",readParam(iBin,"bkgCombL_c2",0,0.1),-10,10);
    RooRealVar bkgCombL_c3("bkgCombL_c3","c3",readParam(iBin,"bkgCombL_c3",0,0.1),-10,10);
    RooRealVar bkgCombL_c4("bkgCombL_c4","c4",readParam(iBin,"bkgCombL_c4",0,0.1),-10,10);
    RooRealVar bkgCombL_c5("bkgCombL_c5","c5",readParam(iBin,"bkgCombL_c5",0,0.5), 0.,10.);
    RooRealVar bkgCombL_c6("bkgCombL_c6","c6",readParam(iBin,"bkgCombL_c6",0,0.1),-10,10.);
    RooRealVar bkgCombL_c7("bkgCombL_c7","c7",readParam(iBin,"bkgCombL_c7",0,0.1),-10,10.);
    RooRealVar bkgCombL_c8("bkgCombL_c8","c8",readParam(iBin,"bkgCombL_c8",0,0.1),-10,10.);
    RooArgSet f_bkgCombL_argset;

    TString f_bkgCombK_format;
    RooRealVar bkgCombK_c1("bkgCombK_c1","c1",readParam(iBin,"bkgCombK_c1",0, 0.1),-10,10);
    RooRealVar bkgCombK_c2("bkgCombK_c2","c2",readParam(iBin,"bkgCombK_c2",0, 0.1),-20,20);
    RooRealVar bkgCombK_c3("bkgCombK_c3","c3",readParam(iBin,"bkgCombK_c3",0, 0.1),-50,50);
    RooRealVar bkgCombK_c4("bkgCombK_c4","c4",readParam(iBin,"bkgCombK_c4",0, 0.1),-100,100);
    RooRealVar bkgCombK_c5("bkgCombK_c5","c5",readParam(iBin,"bkgCombK_c5",0, 0.1),-10,10);
    RooArgSet f_bkgCombK_argset;
    switch (iBin) {
        case 12: //summaryBin[1]
        case 10:
            //f_bkgCombL_format = "exp(-0.5*((CosThetaL-bkgCombL_c1)/bkgCombL_c2)**2)+bkgCombL_c5*exp(-0.5*((CosThetaL-bkgCombL_c3)/bkgCombL_c4)**2)";
            f_bkgCombL_format = "exp(-0.5*((CosThetaL-bkgCombL_c1)/bkgCombL_c2)**2)+bkgCombL_c5*exp(-0.5*((CosThetaL-bkgCombL_c3)/bkgCombL_c4)**2)+bkgCombL_c6*exp(-0.5*((CosThetaL-bkgCombL_c7)/bkgCombL_c8)**2)";
            bkgCombL_c2.setRange(0.,10);
            bkgCombL_c4.setRange(0.,10);
            bkgCombL_c8.setRange(0.,10);
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1));
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c2));
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c3));
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c4));
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c5));
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c6));
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c7));
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c8));
            f_bkgCombK_format = "exp(-0.5*(CosThetaK-bkgCombK_c1)**2/bkgCombK_c2)+bkgCombK_c5*exp(bkgCombK_c3*(bkgCombK_c4-CosThetaK))";
            if (fabs(bkgCombK_c1.getVal())>1) bkgCombK_c1.setVal(-0.5);
            bkgCombK_c1.setRange(-1.,1.);
            if (bkgCombK_c2.getVal()<=0 ) bkgCombK_c2.setVal(0.03);
            bkgCombK_c2.setRange(0.02,0.2);
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c1));
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c2));
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c3));
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c4));
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c5));
            break;
        case 11:
            f_bkgCombL_format = "exp(-0.5*((CosThetaL-bkgCombL_c1)/bkgCombL_c2)**2)+bkgCombL_c5*exp(-0.5*((CosThetaL-bkgCombL_c3)/bkgCombL_c4)**2)+bkgCombL_c6*exp(-0.5*((CosThetaL-bkgCombL_c7)/bkgCombL_c8)**2)";
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1));
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c2));
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c3));
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c4));
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c5));
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c6));
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c7));
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c8));
            f_bkgCombK_format = "exp(bkgCombK_c1*(CosThetaK-bkgCombK_c2))+bkgCombK_c5*exp(bkgCombK_c3*(bkgCombK_c4-CosThetaK))";
            if (fabs(bkgCombK_c1.getVal())>1) bkgCombK_c1.setVal(-0.5);
            bkgCombK_c1.setRange(-1.,1.);
            if (bkgCombK_c2.getVal()<=0 ) bkgCombK_c2.setVal(0.03);
            bkgCombK_c2.setRange(0.02,0.2);
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c1));
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c2));
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c3));
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c4));
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c5));
            break;
        case 4:
            f_bkgCombL_format = "exp(-0.5*((CosThetaL-bkgCombL_c1)/bkgCombL_c2)**2)+bkgCombL_c5*exp(-0.5*((CosThetaL-bkgCombL_c3)/bkgCombL_c4)**2)";
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1));
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c2));
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c3));
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c4));
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c5));
            f_bkgCombK_format = "1+bkgCombK_c1*CosThetaK+bkgCombK_c2*CosThetaK**2+bkgCombK_c3*CosThetaK**3+bkgCombK_c4*CosThetaK**4";
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c1));
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c2));
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c3));
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c4));
            break;
        case 9:
            f_bkgCombL_format = "1+bkgCombL_c1*CosThetaL+bkgCombL_c2*CosThetaL**2+bkgCombL_c3*CosThetaL**3+bkgCombL_c4*CosThetaL**4";
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1));
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c2));
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c3));
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c4));
            f_bkgCombK_format = "1+bkgCombK_c1*CosThetaK+bkgCombK_c2*CosThetaK**2+bkgCombK_c3*CosThetaK**3";
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c1));
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c2));
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c3));
            break;
        default:
            f_bkgCombL_format = "1+bkgCombL_c1*CosThetaL+bkgCombL_c2*CosThetaL**2+bkgCombL_c3*CosThetaL**3+bkgCombL_c4*CosThetaL**4";
            f_bkgCombL_argset.add(RooArgSet(bkgCombL_c1,bkgCombL_c2,bkgCombL_c3,bkgCombL_c4));
            bkgCombL_c1.setVal(0.);
            bkgCombL_c2.setVal(0.1);
            bkgCombL_c3.setVal(0.);
            bkgCombL_c4.setVal(0.1);
            bkgCombL_c5.setVal(0.);
            //bkgCombL_c1.setConstant(kTRUE);
            //bkgCombL_c2.setConstant(kTRUE);
            //bkgCombL_c3.setConstant(kTRUE);
            //bkgCombL_c4.setConstant(kTRUE);
            bkgCombL_c5.setConstant(kTRUE);
            f_bkgCombK_format = "1+bkgCombK_c1*CosThetaK+bkgCombK_c2*CosThetaK**2+bkgCombK_c3*CosThetaK**3";
            f_bkgCombK_argset.add(RooArgSet(bkgCombK_c1,bkgCombK_c2,bkgCombK_c3));
            bkgCombK_c1.setVal(0.);
            bkgCombK_c2.setVal(0.);
            bkgCombK_c3.setVal(0.);
            bkgCombK_c4.setVal(0.);
            bkgCombK_c5.setVal(0.);
            //bkgCombK_c1.setConstant(kTRUE);
            //bkgCombK_c2.setConstant(kTRUE);
            //bkgCombK_c3.setConstant(kTRUE);
            bkgCombK_c4.setConstant(kTRUE);
            bkgCombK_c5.setConstant(kTRUE);
            break;
    }

    f_bkgCombL_argset.add(RooArgSet(CosThetaL));
    RooGenericPdf f_bkgCombL("f_bkgCombL","f_bkgCombL",f_bkgCombL_format, f_bkgCombL_argset);

    f_bkgCombK_argset.add(RooArgSet(CosThetaK));
    RooGenericPdf f_bkgCombK("f_bkgCombK","f_bkgCombK",f_bkgCombK_format,f_bkgCombK_argset);
    
    TString f_bkgCombA_format = TString::Format("(%s)*(%s)",f_bkgCombL_format.Data(),f_bkgCombK_format.Data());
    RooArgSet f_bkgCombA_argset;
    f_bkgCombA_argset.add(f_bkgCombK_argset);
    f_bkgCombA_argset.add(f_bkgCombL_argset);
    RooGenericPdf f_bkgCombA("f_bkgCombA","f_bkgCombA",f_bkgCombA_format,f_bkgCombA_argset);
 
    RooProdPdf f_bkgComb("f_bkgComb", "f_bckComb", RooArgSet(f_bkgCombA, f_bkgCombM));
    RooProdPdf *f_bkgjpsiPeak  = 0;
    RooProdPdf *f_bkgpsi2sPeak = 0;
    switch(iBin){
        case 10:
            f_bkgjpsiPeak  = new RooProdPdf("f_bkgjpsiPeak" , "f_bkgjpsiPeak" ,RooArgSet(*f_bkgjpsiPeakA ,*f_bkgjpsiPeakM12));
            break;
        case 11:
            f_bkgjpsiPeak  = new RooProdPdf("f_bkgjpsiPeak" , "f_bkgjpsiPeak" ,RooArgSet(*f_bkgjpsiPeakA ,*f_bkgjpsiPeakM12 ));
            f_bkgpsi2sPeak = new RooProdPdf("f_bkgpsi2sPeak", "f_bkgpsi2sPeak",RooArgSet(*f_bkgpsi2sPeakA,*f_bkgpsi2sPeakM12));
            break;
        case 4:
            f_bkgpsi2sPeak = new RooProdPdf("f_bkgpsi2sPeak", "f_bkgpsi2sPeak",RooArgSet(*f_bkgpsi2sPeakA,*f_bkgpsi2sPeakM12));
            break;
        case 12:
        default:
            break;
    }
    
    RooRealVar nbkgCombM("nbkgCombM","nbkgCombM",100,0,1E5*scaleFactor);
    RooRealVar nbkgCombA("nbkgCombA","nbkgCombA",100,0,1E5*scaleFactor);
    RooRealVar nbkgComb("nbkgComb","nbkgComb",100,0,1E5*scaleFactor);
    RooRealVar nbkgjpsiPeak("nbkgjpsiPeak","nbkgjpsiPeak",0,0,1E5*scaleFactor);
    RooRealVar nbkgpsi2sPeak("nbkgpsi2sPeak","nbkgpsi2sPeak",0,0,1E5*scaleFactor);
    RooRealVar nbkgComb_dummy("nbkgComb_dummy","nbkgComb_dummy",0);
    nbkgComb_dummy.setConstant(kTRUE);
    RooGenericPdf f_bkgComb_dummy ("f_bkgComb_dummy" ,"f_bkgComb_dummy" ,"1",RooArgSet(Bmass,CosThetaK,CosThetaL));
    RooGenericPdf f_bkgCombA_dummy("f_bkgCombA_dummy","f_bkgCombA_dummy","1",RooArgSet(CosThetaK,CosThetaL));
    RooGenericPdf f_bkgCombM_dummy("f_bkgCombM_dummy","f_bkgCombM_dummy","1",RooArgSet(Bmass));
    RooAddPdf *fM = new RooAddPdf("kernelM","kernelM",RooArgList( f_bkgCombM, f_bkgCombM_dummy),RooArgList(nbkgCombM, nbkgComb_dummy));
    RooAddPdf *fA = 0;
    RooAddPdf *f = 0;
    switch(iBin){
        case 10:
            nbkgjpsiPeak.setVal(nbkgjpsiPeak_MC->getVal()*datasetLumi[0]/datasetLumi[2]*scaleFactor);
            nbkgjpsiPeak.setError(nbkgjpsiPeak_MC->getError()*datasetLumi[0]/datasetLumi[2]*scaleFactor);
            fA = new RooAddPdf("kernelA","kernelA",RooArgList(f_bkgCombA,*f_bkgjpsiPeakA),RooArgList(nbkgCombA,nbkgjpsiPeak));
            f  = new RooAddPdf("kernel" ,"kernel" ,RooArgList(f_bkgComb,*f_bkgjpsiPeak),RooArgList(nbkgComb,nbkgjpsiPeak));
            break;
        case 11: //summaryBin[0]
            nbkgjpsiPeak.setVal(nbkgjpsiPeak_MC->getVal()*datasetLumi[0]/datasetLumi[2]*scaleFactor);
            nbkgjpsiPeak.setError(nbkgjpsiPeak_MC->getError()*datasetLumi[0]/datasetLumi[2]*sqrt(scaleFactor));
            nbkgpsi2sPeak.setVal(nbkgpsi2sPeak_MC->getVal()*datasetLumi[0]/datasetLumi[3]*scaleFactor);
            nbkgpsi2sPeak.setError(nbkgpsi2sPeak_MC->getError()*datasetLumi[0]/datasetLumi[3]*sqrt(scaleFactor));
            fA = new RooAddPdf("kernelA","kernelA",RooArgList(f_bkgCombA,*f_bkgjpsiPeakA,*f_bkgpsi2sPeakA),RooArgList(nbkgCombA,nbkgjpsiPeak,nbkgpsi2sPeak));
            f  = new RooAddPdf("kernel" ,"kernel" ,RooArgList(f_bkgComb ,*f_bkgjpsiPeak,*f_bkgpsi2sPeak),RooArgList(nbkgComb,nbkgjpsiPeak,nbkgpsi2sPeak));
            break;
        case 4:
            nbkgpsi2sPeak.setVal(nbkgpsi2sPeak_MC->getVal()*datasetLumi[0]/datasetLumi[3]*scaleFactor);
            nbkgpsi2sPeak.setError(nbkgpsi2sPeak_MC->getError()*datasetLumi[0]/datasetLumi[3]*sqrt(scaleFactor));
            fA = new RooAddPdf("kernelA","kernelA",RooArgList(f_bkgCombA,*f_bkgpsi2sPeakA),RooArgList(nbkgCombA,nbkgpsi2sPeak));
            f  = new RooAddPdf("kernel" ,"kernel" ,RooArgList(f_bkgComb,*f_bkgpsi2sPeak),RooArgList(nbkgComb,nbkgpsi2sPeak));
            break;
        case 12:// nojpsi for 1<q2<6
        default:
            nbkgjpsiPeak .setVal(0.);
            nbkgpsi2sPeak.setVal(0.);
            fA = new RooAddPdf("kernelA","kernelA",RooArgList( f_bkgCombA,f_bkgCombA_dummy),RooArgList(nbkgCombA,nbkgComb_dummy));
            f  = new RooAddPdf("kernel" ,"kernel" ,RooArgList( f_bkgComb ,f_bkgComb_dummy ),RooArgList(nbkgComb,nbkgComb_dummy));
            break;
    }
    nbkgjpsiPeak .setConstant(kTRUE);
    nbkgpsi2sPeak.setConstant(kTRUE);
    
    // Get sideband data and apply unbinned fit, Minos option doesn't help.
    f_bkgCombL.fitTo(*dataA,Minimizer("Minuit2"));
    f_bkgCombK.fitTo(*dataA,Minimizer("Minuit2"));
    printf("DEBUG\t: Start fitting.\n");
    f_bkgCombA.fitTo(*dataA,Minimizer("Minuit2"));
    printf("DEBUG\t: Start fitting.\n");
    f_bkgCombM.fitTo(*dataM,Minimizer("Minuit2"),Minos(kTRUE));

    // Minuit step with finer precision
    RooAbsReal *nllA = fA->createNLL(*dataA,Extended(kTRUE));// Minos and Save are unknown.
    RooMinuit minuitA(*nllA);
    minuitA.setEps(1e-16);
    minuitA.migrad();
    minuitA.migrad();
    minuitA.save();

    RooAbsReal *nllM = fM->createNLL(*dataM,Extended(kTRUE));// Minos and Save are unknown.
    RooMinuit minuitM(*nllM);
    minuitM.setEps(1e-16);
    minuitM.migrad();
    minuitM.migrad();
    minuitM.save();

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");

    RooPlot* framecosk = CosThetaK.frame(); 
    dataA->plotOn(framecosk,Binning(20)); 
    fA->plotOn(framecosk); 
    fA->plotOn(framecosk,Components(f_bkgCombA),LineColor(2),LineWidth(2),LineStyle(2));
    if ( iBin == 10 || iBin == 11) fA->plotOn(framecosk,Components(*f_bkgjpsiPeakA),LineColor(6),LineWidth(2),LineStyle(2));
    if ( iBin == 4  || iBin == 11) fA->plotOn(framecosk,Components(*f_bkgpsi2sPeakA),LineColor(8),LineWidth(2),LineStyle(2));
    framecosk->SetTitle("");
    framecosk->SetMinimum(0);
    framecosk->Draw();
    
    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = 0;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));

    c->SaveSource(TString::Format("%s/%s_cosk_bin%d.cc",plotpath.Data(),outfile,iBin));
    c->Print(TString::Format("%s/%s_cosk_bin%d.pdf",plotpath.Data(),outfile,iBin));

    RooPlot* framecoskHi = CosThetaK.frame(); 
    dataAHi->plotOn(framecoskHi,Binning(20)); 
    framecoskHi->SetTitle("");
    framecoskHi->SetMinimum(0);
    framecoskHi->Draw();
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->Print(TString::Format("%s/%s_coskHi_bin%d.pdf",plotpath.Data(),outfile,iBin));

    RooPlot* framecoskLo = CosThetaK.frame(); 
    dataALo->plotOn(framecoskLo,Binning(20)); 
    framecoskLo->SetTitle("");
    framecoskLo->SetMinimum(0);
    framecoskLo->Draw();
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->Print(TString::Format("%s/%s_coskLo_bin%d.pdf",plotpath.Data(),outfile,iBin));
    
    
    // 
    RooPlot* framecosl = CosThetaL.frame(); 
    dataA->plotOn(framecosl,Binning(20)); 
    fA->plotOn(framecosl); 
    fA->plotOn(framecosl,Components(f_bkgCombA),LineColor(2),LineWidth(2),LineStyle(2));
    if ( iBin == 10 || iBin == 11) fA->plotOn(framecosl,Components(*f_bkgjpsiPeakA),LineColor(6),LineWidth(2),LineStyle(2));
    if ( iBin == 4  || iBin == 11) fA->plotOn(framecosl,Components(*f_bkgpsi2sPeakA),LineColor(8),LineWidth(2),LineStyle(2));
    framecosl->SetTitle("");
    framecosl->SetMinimum(0);
    framecosl->Draw();

    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->SaveSource(TString::Format("%s/%s_cosl_bin%d.cc",plotpath.Data(),outfile,iBin));
    c->Print(TString::Format("%s/%s_cosl_bin%d.pdf",plotpath.Data(),outfile,iBin));
    
    RooPlot* framecoslHi = CosThetaL.frame(); 
    dataAHi->plotOn(framecoslHi,Binning(20)); 
    framecoslHi->SetTitle("");
    framecoslHi->SetMinimum(0);
    framecoslHi->Draw();
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->Print(TString::Format("%s/%s_coslHi_bin%d.pdf",plotpath.Data(),outfile,iBin));

    RooPlot* framecoslLo = CosThetaL.frame(); 
    dataALo->plotOn(framecoslLo,Binning(20)); 
    framecoslLo->SetTitle("");
    framecoslLo->SetMinimum(0);
    framecoslLo->Draw();
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->Print(TString::Format("%s/%s_coslLo_bin%d.pdf",plotpath.Data(),outfile,iBin));

    //
    RooPlot* framemass = Bmass.frame();
    dataM->plotOn(framemass,Binning(20));
    fM->plotOn(framemass);
    fM->plotOn(framemass,Components(f_bkgCombM),LineColor(2),LineWidth(2),LineStyle(2));
    framemass->SetTitle("");
    framemass->SetMinimum(0);
    framemass->Draw();

    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->SaveSource(TString::Format("%s/%s_bin%d.cc",plotpath.Data(),outfile,iBin));
    c->Print(TString::Format("%s/%s_bin%d.pdf",plotpath.Data(),outfile,iBin));
    
    if (keepParam){
        RooWorkspace *wspace = new RooWorkspace("wspace","wspace");
        Bmass.setRange(5.1,5.6);
        bkgCombL_c1.setConstant(kTRUE);
        bkgCombL_c2.setConstant(kTRUE);
        bkgCombL_c3.setConstant(kTRUE);
        bkgCombL_c4.setConstant(kTRUE);
        bkgCombL_c5.setConstant(kTRUE);
        bkgCombL_c6.setConstant(kTRUE);
        bkgCombL_c7.setConstant(kTRUE);
        bkgCombL_c8.setConstant(kTRUE);
        bkgCombK_c1.setConstant(kTRUE);
        bkgCombK_c2.setConstant(kTRUE);
        bkgCombK_c3.setConstant(kTRUE);
        bkgCombK_c4.setConstant(kTRUE);
        bkgCombK_c5.setConstant(kTRUE);
        bkgCombM_c1.setConstant(kTRUE);
        bkgCombM_c2.setConstant(kTRUE);
        bkgCombM_c3.setConstant(kTRUE);
        bkgCombM_c4.setConstant(kTRUE);
        wspace->import(*f);
        wspace->import(f_bkgCombL);
        wspace->import(f_bkgCombK);

        wspace->writeToFile(TString::Format("%s/wspace_prior_bin%d.root",owspacepath.Data(),iBin),true);

        // Prepare datacard
        double val[3] = {0,0,0};
        val[0] = bkgCombL_c1.getVal();val[1] = bkgCombL_c1.getError();
        writeParam(iBin, "bkgCombL_c1", val , 2 , keepParam);
        val[0] = bkgCombL_c2.getVal();val[1] = bkgCombL_c2.getError();
        writeParam(iBin, "bkgCombL_c2", val , 2 , keepParam);
        val[0] = bkgCombL_c3.getVal();val[1] = bkgCombL_c3.getError();
        writeParam(iBin, "bkgCombL_c3", val , 2 , keepParam);
        val[0] = bkgCombL_c4.getVal();val[1] = bkgCombL_c4.getError();
        writeParam(iBin, "bkgCombL_c4", val , 2 , keepParam);
        val[0] = bkgCombL_c5.getVal();val[1] = bkgCombL_c5.getError();
        writeParam(iBin, "bkgCombL_c5", val , 2 , keepParam);
        val[0] = bkgCombL_c6.getVal();val[1] = bkgCombL_c6.getError();
        writeParam(iBin, "bkgCombL_c6", val , 2 , keepParam);
        val[0] = bkgCombL_c7.getVal();val[1] = bkgCombL_c7.getError();
        writeParam(iBin, "bkgCombL_c7", val , 2 , keepParam);
        val[0] = bkgCombL_c8.getVal();val[1] = bkgCombL_c8.getError();
        writeParam(iBin, "bkgCombL_c8", val , 2 , keepParam);
        val[0] = bkgCombK_c1.getVal();val[1] = bkgCombK_c1.getError();
        writeParam(iBin, "bkgCombK_c1", val , 2 , keepParam);
        val[0] = bkgCombK_c2.getVal();val[1] = bkgCombK_c2.getError();
        writeParam(iBin, "bkgCombK_c2", val , 2 , keepParam);
        val[0] = bkgCombK_c3.getVal();val[1] = bkgCombK_c3.getError();
        writeParam(iBin, "bkgCombK_c3", val , 2 , keepParam);
        val[0] = bkgCombK_c4.getVal();val[1] = bkgCombK_c4.getError();
        writeParam(iBin, "bkgCombK_c4", val , 2 , keepParam);
        val[0] = bkgCombK_c5.getVal();val[1] = bkgCombK_c5.getError();
        writeParam(iBin, "bkgCombK_c5", val , 2 , keepParam);
        val[0] = bkgCombM_c1.getVal();val[1] = bkgCombM_c1.getError();
        writeParam(iBin, "bkgCombM_c1", val , 2 , keepParam);
        val[0] = bkgCombM_c2.getVal();val[1] = bkgCombM_c2.getError();
        writeParam(iBin, "bkgCombM_c2", val , 2 , keepParam);
        val[0] = bkgCombM_c3.getVal();val[1] = bkgCombM_c3.getError();
        writeParam(iBin, "bkgCombM_c3", val , 2 , keepParam);
        val[0] = bkgCombM_c4.getVal();val[1] = bkgCombM_c4.getError();
        writeParam(iBin, "bkgCombM_c4", val , 2 , keepParam);
    }

    // clear
    delete t1;
    delete c;
    delete dataA;
    delete dataM;
}//}}}

void angular3D_prior2(int iBin, const char outfile[] = "angular3D_prior", bool keepParam = false)// DEPRECATED, combination of upper/lower sideband
{//{{{
    // Fit to signal simulation by YsSm+YcCm to determine Sm
    ch->SetBranchStatus("*",0);
    ch->SetBranchStatus("Bmass"         , 1);
    ch->SetBranchStatus("Mumumass"      , 1);
    ch->SetBranchStatus("Mumumasserr"   , 1);
    ch->SetBranchStatus("Kstarmass"     , 1);
    ch->SetBranchStatus("CosTheta*"     , 1);
    ch->SetBranchStatus("Q2"            , 1);
    ch->SetBranchStatus("Triggers"      , 1);
    RooRealVar Bmass("Bmass","M_{K^{*}#Mu#Mu}",4.5,6.5);// Use wider range for sideband
    RooRealVar CosThetaK("CosThetaK"     , "cos#theta_{K}"       , -1. , 1.   ) ;
    RooRealVar CosThetaL("CosThetaL"     , "cos#theta_{L}"       , -1. , 1.   ) ;
    RooRealVar Mumumass("Mumumass","M_{#mu#mu}",0.,10.);
    RooRealVar Mumumasserr("Mumumasserr","Error of M_{#mu#mu}",0.,10.);
    RooRealVar Kstarmass("Kstarmass","M_{K*}",0.,1.5);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    RooRealVar Triggers("Triggers","",0,100);

    int mumuMassWindowBin = 1+2*isCDFcut;
    if (iBin==3 || iBin==5 || isCDFcut < 0 ) mumuMassWindowBin = 0;
    RooDataSet *dataM = new RooDataSet("dataM","data",ch,RooArgSet(Q2, Bmass, CosThetaK, CosThetaL, Mumumass, Mumumasserr, Kstarmass, Triggers),TString::Format("(%s) && (%s) && (%s) && (%s)",nTriggeredPath[2],q2range[iBin],mumuMassWindow[mumuMassWindowBin],kstarMassWindow[2]),0);
    RooDataSet *dataA = new RooDataSet("dataA","data",ch,RooArgSet(Q2, Bmass, CosThetaK, CosThetaL, Mumumass, Mumumasserr, Kstarmass, Triggers),TString::Format("(%s) && (%s) && (%s) && (%s) && (%s)",nTriggeredPath[2],q2range[iBin],mumuMassWindow[mumuMassWindowBin],kstarMassWindow[1],"Bmass > 5.6 || Bmass < 5.0"),0);
    RooDataSet *dataAUp = new RooDataSet("dataAUp","data",ch,RooArgSet(Q2, Bmass, CosThetaK, CosThetaL, Mumumass, Mumumasserr, Kstarmass, Triggers),TString::Format("(%s) && (%s) && (%s) && (%s) && (%s)",nTriggeredPath[2],q2range[iBin],mumuMassWindow[mumuMassWindowBin],kstarMassWindow[1],"Bmass > 5.6"),0);
    RooDataSet *dataALo = new RooDataSet("dataALo","data",ch,RooArgSet(Q2, Bmass, CosThetaK, CosThetaL, Mumumass, Mumumasserr, Kstarmass, Triggers),TString::Format("(%s) && (%s) && (%s) && (%s) && (%s)",nTriggeredPath[2],q2range[iBin],mumuMassWindow[mumuMassWindowBin],kstarMassWindow[1],"Bmass < 5.0"),0);
    if (dataM->sumEntries() == 0 || dataAUp->sumEntries() == 0){
        printf("ERROR   : Total number of entries is 0. Exit!\n");
        return;
    }
    printf("INFO: f_bkgPeak prepared.\n");
    
    // Create combinatorial background distribution
    TString f_bkgCombM_format = "1+bkgCombM_c1*(Bmass-5)+bkgCombM_c2*exp(bkgCombM_c3*(Bmass-bkgCombM_c4))";
    RooRealVar bkgCombM_c1("bkgCombM_c1","c1",readParam(iBin,"bkgCombM_c1",0, 0),-20,20);
    RooRealVar bkgCombM_c2("bkgCombM_c2","c2",readParam(iBin,"bkgCombM_c2",0,10),-20,20);
    RooRealVar bkgCombM_c3("bkgCombM_c3","c3",readParam(iBin,"bkgCombM_c3",0,-3),-10,10);
    RooRealVar bkgCombM_c4("bkgCombM_c4","c4",readParam(iBin,"bkgCombM_c4",0, 5),-10,10);
    RooArgSet f_bkgCombM_argset;
    f_bkgCombM_argset.add(RooArgSet(bkgCombM_c1,bkgCombM_c2,bkgCombM_c3,bkgCombM_c4));
    f_bkgCombM_argset.add(RooArgSet(Bmass));
    RooGenericPdf f_bkgCombM("f_bkgCombM","f_bkgCombM",f_bkgCombM_format,f_bkgCombM_argset);

    TString f_bkgCombLUp_format;
    RooRealVar bkgCombLUp_c1("bkgCombLUp_c1","c1",readParam(iBin,"bkgCombLUp_c1",0,0.1),-10,10);
    RooRealVar bkgCombLUp_c2("bkgCombLUp_c2","c2",readParam(iBin,"bkgCombLUp_c2",0,0.1),-10,10);
    RooRealVar bkgCombLUp_c3("bkgCombLUp_c3","c3",readParam(iBin,"bkgCombLUp_c3",0,0.1),-10,10);
    RooRealVar bkgCombLUp_c4("bkgCombLUp_c4","c4",readParam(iBin,"bkgCombLUp_c4",0,0.1),-10,10);
    RooRealVar bkgCombLUp_c5("bkgCombLUp_c5","c5",readParam(iBin,"bkgCombLUp_c5",0,0.5),-10,10);
    RooRealVar bkgCombLUp_c6("bkgCombLUp_c6","c6",readParam(iBin,"bkgCombLUp_c6",0,0.1),-10,10);
    RooRealVar bkgCombLUp_c7("bkgCombLUp_c7","c7",readParam(iBin,"bkgCombLUp_c7",0,0.1),-10,10);
    RooRealVar bkgCombLUp_c8("bkgCombLUp_c8","c8",readParam(iBin,"bkgCombLUp_c8",0,0.1),-10,10);
    RooArgSet f_bkgCombLUp_argset;

    TString f_bkgCombLLo_format;
    RooRealVar bkgCombLLo_c1("bkgCombLLo_c1","c1",readParam(iBin,"bkgCombLLo_c1",0,0.1),-10,10);
    RooRealVar bkgCombLLo_c2("bkgCombLLo_c2","c2",readParam(iBin,"bkgCombLLo_c2",0,0.1),-10,10);
    RooRealVar bkgCombLLo_c3("bkgCombLLo_c3","c3",readParam(iBin,"bkgCombLLo_c3",0,0.1),-10,10);
    RooRealVar bkgCombLLo_c4("bkgCombLLo_c4","c4",readParam(iBin,"bkgCombLLo_c4",0,0.1),-10,10);
    RooRealVar bkgCombLLo_c5("bkgCombLLo_c5","c5",readParam(iBin,"bkgCombLLo_c5",0,0.5),-10,10);
    RooRealVar bkgCombLLo_c6("bkgCombLLo_c6","c6",readParam(iBin,"bkgCombLLo_c6",0,0.1),-10,10);
    RooRealVar bkgCombLLo_c7("bkgCombLLo_c7","c7",readParam(iBin,"bkgCombLLo_c7",0,0.1),-10,10);
    RooRealVar bkgCombLLo_c8("bkgCombLLo_c8","c8",readParam(iBin,"bkgCombLLo_c8",0,0.1),-10,10);
    RooArgSet f_bkgCombLLo_argset;

    TString f_bkgCombKUp_format;
    RooRealVar bkgCombKUp_c1("bkgCombKUp_c1","c1",readParam(iBin,"bkgCombKUp_c1",0,0.1),-20,20);
    RooRealVar bkgCombKUp_c2("bkgCombKUp_c2","c2",readParam(iBin,"bkgCombKUp_c2",0,0.1),-20,20);
    RooRealVar bkgCombKUp_c3("bkgCombKUp_c3","c3",readParam(iBin,"bkgCombKUp_c3",0,0.1),-20,20);
    RooRealVar bkgCombKUp_c4("bkgCombKUp_c4","c4",readParam(iBin,"bkgCombKUp_c4",0,0.1),-50,50);
    RooRealVar bkgCombKUp_c5("bkgCombKUp_c5","c5",readParam(iBin,"bkgCombKUp_c5",0,0.1),-20,20);
    RooRealVar bkgCombKUp_c6("bkgCombKUp_c6","c6",readParam(iBin,"bkgCombKUp_c6",0,0.1),-20,20);
    RooRealVar bkgCombKUp_c7("bkgCombKUp_c7","c7",readParam(iBin,"bkgCombKUp_c7",0,0.1),-20,20);
    RooRealVar bkgCombKUp_c8("bkgCombKUp_c8","c8",readParam(iBin,"bkgCombKUp_c8",0,0.1),-20,20);
    RooArgSet f_bkgCombKUp_argset;

    TString f_bkgCombKLo_format;
    RooRealVar bkgCombKLo_c1("bkgCombKLo_c1","c1",readParam(iBin,"bkgCombKLo_c1",0,0.1),-20,20);
    RooRealVar bkgCombKLo_c2("bkgCombKLo_c2","c2",readParam(iBin,"bkgCombKLo_c2",0,0.1),-20,20);
    RooRealVar bkgCombKLo_c3("bkgCombKLo_c3","c3",readParam(iBin,"bkgCombKLo_c3",0,0.1),-20,20);
    RooRealVar bkgCombKLo_c4("bkgCombKLo_c4","c4",readParam(iBin,"bkgCombKLo_c4",0,0.1),-50,50);
    RooRealVar bkgCombKLo_c5("bkgCombKLo_c5","c5",readParam(iBin,"bkgCombKLo_c5",0,0.1),-20,20);
    RooRealVar bkgCombKLo_c6("bkgCombKLo_c6","c6",readParam(iBin,"bkgCombKLo_c6",0,0.1),-20,20);
    RooRealVar bkgCombKLo_c7("bkgCombKLo_c7","c7",readParam(iBin,"bkgCombKLo_c7",0,0.1),-20,20);
    RooRealVar bkgCombKLo_c8("bkgCombKLo_c8","c8",readParam(iBin,"bkgCombKLo_c8",0,0.1),-20,20);
    RooArgSet f_bkgCombKLo_argset;

    switch (iBin) {
        case 12: //summaryBin[1]
        case 10:
            f_bkgCombLUp_format = "exp(-0.5*((CosThetaL-bkgCombLUp_c1)/bkgCombLUp_c2)**2)";// ok
            bkgCombLUp_c2.setRange(0.,20);
            f_bkgCombLUp_argset.add(RooArgSet(bkgCombLUp_c1));
            f_bkgCombLUp_argset.add(RooArgSet(bkgCombLUp_c2));
            f_bkgCombLLo_format = "exp(-0.5*((CosThetaL-bkgCombLLo_c1)/bkgCombLLo_c2)**2)+bkgCombLLo_c5*exp(-0.5*((CosThetaL-bkgCombLLo_c3)/bkgCombLLo_c4)**2)+bkgCombLLo_c6+bkgCombLLo_c7*CosThetaL+bkgCombLLo_c8*CosThetaL**2";//ok
            bkgCombLLo_c2.setRange(0.,20);
            bkgCombLLo_c4.setRange(0.,20);
            f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c1));
            f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c2));
            f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c3));
            f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c4));
            f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c5));
            f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c6));
            f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c7));
            f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c8));
            //f_bkgCombK_format = "exp(-0.5*(CosThetaK-bkgCombK_c1)**2/bkgCombK_c2)+bkgCombK_c5*exp(bkgCombK_c3*(bkgCombK_c4-CosThetaK))";
            //f_bkgCombKUp_format = "exp(-0.5*((CosThetaK-bkgCombKUp_c1)/bkgCombKUp_c2)**2)+bkgCombKUp_c3*exp(bkgCombKUp_c4*CosThetaK+bkgCombKUp_c5)+bkgCombKUp_c6+bkgCombKUp_c7*CosThetaK+bkgCombKUp_c8*CosThetaK**2";// not ok
            //f_bkgCombKUp_format = "exp(-0.5*((CosThetaK-bkgCombKUp_c1)/bkgCombKUp_c2)**2)+bkgCombKUp_c6+bkgCombKUp_c7*CosThetaK+bkgCombKUp_c8*CosThetaK**2";// not ok
            //f_bkgCombKUp_format = "exp(-0.5*((CosThetaK-bkgCombKUp_c1)/bkgCombKUp_c2)**2)+bkgCombKUp_c5*exp(bkgCombKUp_c3*CosThetaK+bkgCombKUp_c4)";
            //f_bkgCombKUp_format = "exp(-0.5*((CosThetaK-bkgCombKUp_c1)/bkgCombKUp_c2)**2)+bkgCombKUp_c5*exp(bkgCombKUp_c3*(CosThetaK+bkgCombKUp_c4))";
            f_bkgCombKUp_format = "exp(bkgCombKUp_c2*(CosThetaK-bkgCombKUp_c1))+bkgCombKUp_c5*exp(bkgCombKUp_c3*(CosThetaK-bkgCombKUp_c4))";
            //f_bkgCombKUp_format = "bkgCombKUp_c1+bkgCombKUp_c2*CosThetaK+bkgCombKUp_c3*CosThetaK**2+exp(bkgCombKUp_c4*(CosThetaK+bkgCombKUp_c5))";
            //bkgCombKUp_c3.setRange(0,20);
            f_bkgCombKUp_argset.add(RooArgSet(bkgCombKUp_c1));
            f_bkgCombKUp_argset.add(RooArgSet(bkgCombKUp_c2));
            f_bkgCombKUp_argset.add(RooArgSet(bkgCombKUp_c3));
            f_bkgCombKUp_argset.add(RooArgSet(bkgCombKUp_c4));
            f_bkgCombKUp_argset.add(RooArgSet(bkgCombKUp_c5));
            //f_bkgCombKUp_argset.add(RooArgSet(bkgCombKUp_c6));
            //f_bkgCombKUp_argset.add(RooArgSet(bkgCombKUp_c7));
            //f_bkgCombKUp_argset.add(RooArgSet(bkgCombKUp_c8));
            f_bkgCombKLo_format = "exp(bkgCombKLo_c2*CosThetaK-bkgCombKLo_c1)+bkgCombKLo_c5*exp(bkgCombKLo_c3*CosThetaK+bkgCombKLo_c4)";//ok
            //f_bkgCombKLo_format = "bkgCombKLo_c1+bkgCombKLo_c2*CosThetaK+bkgCombKLo_c3*CosThetaK**2+exp(bkgCombKLo_c4*(CosThetaK+bkgCombKLo_c5))";//ok
            //bkgCombKLo_c5.setRange(0,20);
            f_bkgCombKLo_argset.add(RooArgSet(bkgCombKLo_c1));
            f_bkgCombKLo_argset.add(RooArgSet(bkgCombKLo_c2));
            f_bkgCombKLo_argset.add(RooArgSet(bkgCombKLo_c3));
            f_bkgCombKLo_argset.add(RooArgSet(bkgCombKLo_c4));
            f_bkgCombKLo_argset.add(RooArgSet(bkgCombKLo_c5));
            break;
        case 11:
            //f_bkgCombLUp_format = "exp(-0.5*((CosThetaL-bkgCombLUp_c1)/bkgCombLUp_c2)**2)+bkgCombLUp_c5*exp(-0.5*((CosThetaL-bkgCombLUp_c3)/bkgCombLUp_c4)**2)+bkgCombLUp_c6*exp(-0.5*((CosThetaL-bkgCombLUp_c7)/bkgCombLUp_c8)**2)";
            f_bkgCombLUp_format = "exp(-0.5*((CosThetaL-bkgCombLUp_c1)/bkgCombLUp_c2)**2)";
            if (bkgCombLUp_c5.getVal()<=0 ) bkgCombLUp_c5.setVal(0.1);
            bkgCombLUp_c5.setRange(0,20);
            if (bkgCombLUp_c6.getVal()<=0 ) bkgCombLUp_c6.setVal(0.1);
            bkgCombLUp_c6.setRange(0,20);
            f_bkgCombLUp_argset.add(RooArgSet(bkgCombLUp_c1));
            f_bkgCombLUp_argset.add(RooArgSet(bkgCombLUp_c2));
            f_bkgCombLUp_argset.add(RooArgSet(bkgCombLUp_c3));
            //f_bkgCombLUp_argset.add(RooArgSet(bkgCombLUp_c4));
            //f_bkgCombLUp_argset.add(RooArgSet(bkgCombLUp_c5));
            //f_bkgCombLUp_argset.add(RooArgSet(bkgCombLUp_c6));
            //f_bkgCombLUp_argset.add(RooArgSet(bkgCombLUp_c7));
            //f_bkgCombLUp_argset.add(RooArgSet(bkgCombLUp_c8));
            //f_bkgCombLLo_format = "exp(-0.5*((CosThetaL-bkgCombLLo_c1)/bkgCombLLo_c2)**2)+bkgCombLLo_c5*exp(-0.5*((CosThetaL-bkgCombLLo_c3)/bkgCombLLo_c4)**2)+bkgCombLLo_c6*exp(-0.5*((CosThetaL-bkgCombLLo_c7)/bkgCombLLo_c8)**2)";//ok, but not good
            //f_bkgCombLLo_format = "exp(-0.5*((CosThetaL-bkgCombLLo_c1)/bkgCombLLo_c2)**2)+bkgCombLLo_c3+bkgCombLLo_c4*CosThetaL+bkgCombLLo_c5*CosThetaL**2+bkgCombLLo_c6*CosThetaL**3";// not ok
            f_bkgCombLLo_format = "exp(-0.5*((CosThetaL-bkgCombLLo_c1)/bkgCombLLo_c2)**2)+bkgCombLLo_c3*exp(-0.5*((CosThetaL-bkgCombLLo_c4)/bkgCombLLo_c5)**2)";//
            //f_bkgCombLLo_format = "1+bkgCombLLo_c1*CosThetaL+bkgCombLLo_c2*CosThetaL**2+bkgCombLLo_c3*CosThetaL**3+bkgCombLLo_c4*CosThetaL**4+bkgCombLLo_c5*CosThetaL**5";// not ok at all
            //f_bkgCombLLo_format = "1+bkgCombLLo_c1*CosThetaL+bkgCombLLo_c2*CosThetaL**2+bkgCombLLo_c3*CosThetaL**3+bkgCombLLo_c4*CosThetaL**4+bkgCombLLo_c5*CosThetaL**5";
            if (bkgCombLLo_c3.getVal()<=0 ) bkgCombLLo_c3.setVal(0.1);
            bkgCombLLo_c3.setRange(0,20);
            if (bkgCombLLo_c5.getVal()<=0 ) bkgCombLLo_c5.setVal(0.1);
            bkgCombLLo_c5.setRange(0.05,20);
            f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c1));
            f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c2));
            f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c3));
            f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c4));
            f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c5));
            f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c6));
            //f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c7));
            //f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c8));
            f_bkgCombKUp_format = "exp(bkgCombKUp_c1*(CosThetaK-bkgCombKUp_c2))+bkgCombKUp_c5*exp(bkgCombKUp_c3*(bkgCombKUp_c4-CosThetaK))";
            if (bkgCombKUp_c5.getVal()<=0 ) bkgCombKUp_c5.setVal(0.1);
            bkgCombKUp_c5.setRange(0,20);
            f_bkgCombKUp_argset.add(RooArgSet(bkgCombKUp_c1));
            f_bkgCombKUp_argset.add(RooArgSet(bkgCombKUp_c2));
            f_bkgCombKUp_argset.add(RooArgSet(bkgCombKUp_c3));
            f_bkgCombKUp_argset.add(RooArgSet(bkgCombKUp_c4));
            f_bkgCombKUp_argset.add(RooArgSet(bkgCombKUp_c5));
            f_bkgCombKLo_format = "exp(bkgCombKLo_c1*(CosThetaK-bkgCombKLo_c2))+bkgCombKLo_c5*exp(bkgCombKLo_c3*(bkgCombKLo_c4-CosThetaK))";
            if (bkgCombKLo_c5.getVal()<=0 ) bkgCombKLo_c5.setVal(0.1);
            bkgCombKLo_c5.setRange(0,20);
            f_bkgCombKLo_argset.add(RooArgSet(bkgCombKLo_c1));
            f_bkgCombKLo_argset.add(RooArgSet(bkgCombKLo_c2));
            f_bkgCombKLo_argset.add(RooArgSet(bkgCombKLo_c3));
            f_bkgCombKLo_argset.add(RooArgSet(bkgCombKLo_c4));
            f_bkgCombKLo_argset.add(RooArgSet(bkgCombKLo_c5));
            break;
        case 4:
            f_bkgCombLUp_format = "exp(-0.5*((CosThetaL-bkgCombLUp_c1)/bkgCombLUp_c2)**2)+bkgCombLUp_c5*exp(-0.5*((CosThetaL-bkgCombLUp_c3)/bkgCombLUp_c4)**2)";
            f_bkgCombLUp_argset.add(RooArgSet(bkgCombLUp_c1));
            f_bkgCombLUp_argset.add(RooArgSet(bkgCombLUp_c2));
            f_bkgCombLUp_argset.add(RooArgSet(bkgCombLUp_c3));
            f_bkgCombLUp_argset.add(RooArgSet(bkgCombLUp_c4));
            f_bkgCombLUp_argset.add(RooArgSet(bkgCombLUp_c5));
            f_bkgCombLLo_format = "exp(-0.5*((CosThetaL-bkgCombLLo_c1)/bkgCombLLo_c2)**2)+bkgCombLLo_c5*exp(-0.5*((CosThetaL-bkgCombLLo_c3)/bkgCombLLo_c4)**2)";
            f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c1));
            f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c2));
            f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c3));
            f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c4));
            f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c5));
            f_bkgCombKUp_format = "1+bkgCombKUp_c1*CosThetaK+bkgCombKUp_c2*CosThetaK**2+bkgCombKUp_c3*CosThetaK**3+bkgCombKUp_c4*CosThetaK**4";
            f_bkgCombKUp_argset.add(RooArgSet(bkgCombKUp_c1));
            f_bkgCombKUp_argset.add(RooArgSet(bkgCombKUp_c2));
            f_bkgCombKUp_argset.add(RooArgSet(bkgCombKUp_c3));
            f_bkgCombKUp_argset.add(RooArgSet(bkgCombKUp_c4));
            f_bkgCombKLo_format = "1+bkgCombKLo_c1*CosThetaK+bkgCombKLo_c2*CosThetaK**2+bkgCombKLo_c3*CosThetaK**3+bkgCombKLo_c4*CosThetaK**4";
            f_bkgCombKLo_argset.add(RooArgSet(bkgCombKLo_c1));
            f_bkgCombKLo_argset.add(RooArgSet(bkgCombKLo_c2));
            f_bkgCombKLo_argset.add(RooArgSet(bkgCombKLo_c3));
            f_bkgCombKLo_argset.add(RooArgSet(bkgCombKLo_c4));
            break;
        case 9:
            f_bkgCombLUp_format = "1+bkgCombLUp_c1*CosThetaL+bkgCombLUp_c2*CosThetaL**2+bkgCombLUp_c3*CosThetaL**3+bkgCombLUp_c4*CosThetaL**4";
            f_bkgCombLUp_argset.add(RooArgSet(bkgCombLUp_c1));
            f_bkgCombLUp_argset.add(RooArgSet(bkgCombLUp_c2));
            f_bkgCombLUp_argset.add(RooArgSet(bkgCombLUp_c3));
            f_bkgCombLUp_argset.add(RooArgSet(bkgCombLUp_c4));
            f_bkgCombLLo_format = "1+bkgCombLLo_c1*CosThetaL+bkgCombLLo_c2*CosThetaL**2+bkgCombLLo_c3*CosThetaL**3+bkgCombLLo_c4*CosThetaL**4";
            f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c1));
            f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c2));
            f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c3));
            f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c4));
            f_bkgCombKUp_format = "1+bkgCombKUp_c1*CosThetaK+bkgCombKUp_c2*CosThetaK**2+bkgCombKUp_c3*CosThetaK**3";
            f_bkgCombKUp_argset.add(RooArgSet(bkgCombKUp_c1));
            f_bkgCombKUp_argset.add(RooArgSet(bkgCombKUp_c2));
            f_bkgCombKUp_argset.add(RooArgSet(bkgCombKUp_c3));
            f_bkgCombKLo_format = "1+bkgCombKLo_c1*CosThetaK+bkgCombKLo_c2*CosThetaK**2+bkgCombKLo_c3*CosThetaK**3";
            f_bkgCombKLo_argset.add(RooArgSet(bkgCombKLo_c1));
            f_bkgCombKLo_argset.add(RooArgSet(bkgCombKLo_c2));
            f_bkgCombKLo_argset.add(RooArgSet(bkgCombKLo_c3));
            break;
        default:
            f_bkgCombLUp_format = "1+bkgCombLUp_c1*CosThetaL+bkgCombLUp_c2*CosThetaL**2+bkgCombLUp_c3*CosThetaL**3+bkgCombLUp_c4*CosThetaL**4";
            f_bkgCombLUp_argset.add(RooArgSet(bkgCombLUp_c1,bkgCombLUp_c2,bkgCombLUp_c3,bkgCombLUp_c4));
            bkgCombLUp_c1.setVal(0.);
            bkgCombLUp_c2.setVal(0.1);
            bkgCombLUp_c3.setVal(0.);
            bkgCombLUp_c4.setVal(0.1);
            bkgCombLUp_c5.setVal(0.);
            bkgCombLUp_c5.setConstant(kTRUE);
            f_bkgCombLLo_format = "1+bkgCombLLo_c1*CosThetaL+bkgCombLLo_c2*CosThetaL**2+bkgCombLLo_c3*CosThetaL**3+bkgCombLLo_c4*CosThetaL**4";
            f_bkgCombLLo_argset.add(RooArgSet(bkgCombLLo_c1,bkgCombLLo_c2,bkgCombLLo_c3,bkgCombLLo_c4));
            bkgCombLLo_c1.setVal(0.);
            bkgCombLLo_c2.setVal(0.1);
            bkgCombLLo_c3.setVal(0.);
            bkgCombLLo_c4.setVal(0.1);
            bkgCombLLo_c5.setVal(0.);
            bkgCombLLo_c5.setConstant(kTRUE);
            f_bkgCombKUp_format = "1+bkgCombKUp_c1*CosThetaK+bkgCombKUp_c2*CosThetaK**2+bkgCombKUp_c3*CosThetaK**3";
            f_bkgCombKUp_argset.add(RooArgSet(bkgCombKUp_c1,bkgCombKUp_c2,bkgCombKUp_c3));
            bkgCombKUp_c1.setVal(0.);
            bkgCombKUp_c2.setVal(0.);
            bkgCombKUp_c3.setVal(0.);
            bkgCombKUp_c4.setVal(0.);
            bkgCombKUp_c5.setVal(0.);
            bkgCombKUp_c4.setConstant(kTRUE);
            bkgCombKUp_c5.setConstant(kTRUE);
            f_bkgCombKLo_format = "1+bkgCombKLo_c1*CosThetaK+bkgCombKLo_c2*CosThetaK**2+bkgCombKLo_c3*CosThetaK**3";
            f_bkgCombKLo_argset.add(RooArgSet(bkgCombKLo_c1,bkgCombKLo_c2,bkgCombKLo_c3));
            bkgCombKLo_c1.setVal(0.);
            bkgCombKLo_c2.setVal(0.);
            bkgCombKLo_c3.setVal(0.);
            bkgCombKLo_c4.setVal(0.);
            bkgCombKLo_c5.setVal(0.);
            bkgCombKLo_c4.setConstant(kTRUE);
            bkgCombKLo_c5.setConstant(kTRUE);
            break;
    }

    f_bkgCombLUp_argset.add(RooArgSet(CosThetaL));
    f_bkgCombLLo_argset.add(RooArgSet(CosThetaL));
    RooGenericPdf f_bkgCombLUp("f_bkgCombLUp","f_bkgCombLUp",f_bkgCombLUp_format, f_bkgCombLUp_argset);
    RooGenericPdf f_bkgCombLLo("f_bkgCombLLo","f_bkgCombLLo",f_bkgCombLLo_format, f_bkgCombLLo_argset);

    f_bkgCombKUp_argset.add(RooArgSet(CosThetaK));
    f_bkgCombKLo_argset.add(RooArgSet(CosThetaK));
    RooGenericPdf f_bkgCombKUp("f_bkgCombKUp","f_bkgCombKUp",f_bkgCombKUp_format,f_bkgCombKUp_argset);
    RooGenericPdf f_bkgCombKLo("f_bkgCombKLo","f_bkgCombKLo",f_bkgCombKLo_format,f_bkgCombKLo_argset);
    
    TString f_bkgCombAUp_format = TString::Format("(%s)*(%s)",f_bkgCombLUp_format.Data(),f_bkgCombKUp_format.Data());
    RooArgSet f_bkgCombAUp_argset;
    f_bkgCombAUp_argset.add(f_bkgCombKUp_argset);
    f_bkgCombAUp_argset.add(f_bkgCombLUp_argset);
    RooGenericPdf f_bkgCombAUp("f_bkgCombAUp","f_bkgCombAUp",f_bkgCombAUp_format,f_bkgCombAUp_argset);

    TString f_bkgCombALo_format = TString::Format("(%s)*(%s)",f_bkgCombLLo_format.Data(),f_bkgCombKLo_format.Data());
    RooArgSet f_bkgCombALo_argset;
    f_bkgCombALo_argset.add(f_bkgCombKLo_argset);
    f_bkgCombALo_argset.add(f_bkgCombLLo_argset);
    RooGenericPdf f_bkgCombALo("f_bkgCombALo","f_bkgCombALo",f_bkgCombALo_format,f_bkgCombALo_argset);

    TString f_bkgCombA_format = TString::Format("frac_bkgCombA_UpLo*(%s)+(1-frac_bkgCombA_UpLo)*(%s)",f_bkgCombAUp_format.Data(),f_bkgCombALo_format.Data());
    RooRealVar frac_bkgCombA_UpLo("frac_bkgCombA_UpLo","frac_bkgCombA_UpLo",0.5,0.01,0.99);
    RooArgSet f_bkgCombA_argset;
    f_bkgCombA_argset.add(frac_bkgCombA_UpLo);
    f_bkgCombA_argset.add(f_bkgCombAUp_argset);
    f_bkgCombA_argset.add(f_bkgCombALo_argset);
    RooGenericPdf f_bkgCombA("f_bkgCombA","f_bkgCombA",f_bkgCombA_format,f_bkgCombA_argset);
 
    RooProdPdf f_bkgComb("f_bkgComb", "f_bckComb", RooArgSet(f_bkgCombA, f_bkgCombM));
    RooRealVar nbkgCombM("nbkgCombM","nbkgCombM",100,0,1E5*scaleFactor);
    RooRealVar nbkgCombA("nbkgCombA","nbkgCombA",100,0,1E5*scaleFactor);
    RooRealVar nbkgComb ("nbkgComb" ,"nbkgComb" ,100,0,1E5*scaleFactor);
    RooRealVar nbkgComb_dummy("nbkgComb_dummy","nbkgComb_dummy",0);
    nbkgComb_dummy.setConstant(kTRUE);
    RooGenericPdf f_bkgComb_dummy ("f_bkgComb_dummy" ,"f_bkgComb_dummy" ,"1",RooArgSet(Bmass,CosThetaK,CosThetaL));
    RooGenericPdf f_bkgCombA_dummy("f_bkgCombA_dummy","f_bkgCombA_dummy","1",RooArgSet(CosThetaK,CosThetaL));
    RooGenericPdf f_bkgCombM_dummy("f_bkgCombM_dummy","f_bkgCombM_dummy","1",RooArgSet(Bmass));
    RooAddPdf *fM = new RooAddPdf("kernelM","kernelM",RooArgList( f_bkgCombM,f_bkgCombM_dummy),RooArgList(nbkgCombM,nbkgComb_dummy));
    RooAddPdf *fA = new RooAddPdf("kernelA","kernelA",RooArgList( f_bkgCombA,f_bkgCombA_dummy),RooArgList(nbkgCombA,nbkgComb_dummy));
    RooAddPdf *f  = new RooAddPdf("kernel" ,"kernel" ,RooArgList( f_bkgComb ,f_bkgComb_dummy ),RooArgList(nbkgComb ,nbkgComb_dummy));
    
    // Get sideband data and apply unbinned fit.
    f_bkgCombLUp.fitTo(*dataAUp,Minimizer("Minuit"),Minos(kFALSE));
    f_bkgCombLLo.fitTo(*dataALo,Minimizer("Minuit"),Minos(kFALSE));
    f_bkgCombKUp.fitTo(*dataAUp,Minimizer("Minuit"),Minos(kFALSE));
    f_bkgCombKLo.fitTo(*dataALo,Minimizer("Minuit"),Minos(kFALSE));
    bkgCombLUp_c1.setConstant(kTRUE);
    bkgCombLUp_c2.setConstant(kTRUE);
    bkgCombLUp_c3.setConstant(kTRUE);
    bkgCombLUp_c4.setConstant(kTRUE);
    bkgCombLUp_c5.setConstant(kTRUE);
    bkgCombLUp_c6.setConstant(kTRUE);
    bkgCombLUp_c7.setConstant(kTRUE);
    bkgCombLUp_c8.setConstant(kTRUE);
    bkgCombLLo_c1.setConstant(kTRUE);
    bkgCombLLo_c2.setConstant(kTRUE);
    bkgCombLLo_c3.setConstant(kTRUE);
    bkgCombLLo_c4.setConstant(kTRUE);
    bkgCombLLo_c5.setConstant(kTRUE);
    bkgCombLLo_c6.setConstant(kTRUE);
    bkgCombLLo_c7.setConstant(kTRUE);
    bkgCombLLo_c8.setConstant(kTRUE);
    bkgCombKUp_c1.setConstant(kTRUE);
    bkgCombKUp_c2.setConstant(kTRUE);
    bkgCombKUp_c3.setConstant(kTRUE);
    bkgCombKUp_c4.setConstant(kTRUE);
    bkgCombKUp_c5.setConstant(kTRUE);
    bkgCombKUp_c6.setConstant(kTRUE);
    bkgCombKUp_c7.setConstant(kTRUE);
    bkgCombKUp_c8.setConstant(kTRUE);
    bkgCombKLo_c1.setConstant(kTRUE);
    bkgCombKLo_c2.setConstant(kTRUE);
    bkgCombKLo_c3.setConstant(kTRUE);
    bkgCombKLo_c4.setConstant(kTRUE);
    bkgCombKLo_c5.setConstant(kTRUE);
    bkgCombKLo_c6.setConstant(kTRUE);
    bkgCombKLo_c7.setConstant(kTRUE);
    bkgCombKLo_c8.setConstant(kTRUE);
    RooAbsReal *nllA = fA->createNLL(*dataA,Extended(kTRUE));// Minos and Save are unknown.
    RooMinuit minuitA(*nllA);
    minuitA.setEps(1e-16);
    minuitA.migrad();
    minuitA.migrad();
    minuitA.save();
    
    RooAbsReal *nllM = fM->createNLL(*dataM,Extended(kTRUE));// Minos and Save are unknown.
    RooMinuit minuitM(*nllM);
    minuitM.setEps(1e-16);
    minuitM.migrad();
    minuitM.migrad();
    minuitM.save();
    bkgCombM_c1.setConstant(kTRUE);
    bkgCombM_c2.setConstant(kTRUE);
    bkgCombM_c3.setConstant(kTRUE);
    bkgCombM_c4.setConstant(kTRUE);

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");

    RooPlot* framecosk = CosThetaK.frame(); 
    dataA->plotOn(framecosk,Binning(20)); 
    fA->plotOn(framecosk); 
    framecosk->SetTitle("");
    framecosk->SetMinimum(0);
    framecosk->Draw();
    
    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = 0;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));

    c->SaveSource(TString::Format("%s/%s_cosk_bin%d.cc",plotpath.Data(),outfile,iBin));
    c->Print(TString::Format("%s/%s_cosk_bin%d.pdf",plotpath.Data(),outfile,iBin));

    RooPlot* framecoskUp = CosThetaK.frame(); 
    dataAUp->plotOn(framecoskUp,Binning(20)); 
    f_bkgCombKUp.plotOn(framecoskUp);
    framecoskUp->SetTitle("");
    framecoskUp->SetMinimum(0);
    framecoskUp->Draw();
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->Print(TString::Format("%s/%s_coskUp_bin%d.pdf",plotpath.Data(),outfile,iBin));

    RooPlot* framecoskLo = CosThetaK.frame(); 
    dataALo->plotOn(framecoskLo,Binning(20)); 
    f_bkgCombKLo.plotOn(framecoskLo);
    framecoskLo->SetTitle("");
    framecoskLo->SetMinimum(0);
    framecoskLo->Draw();
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->Print(TString::Format("%s/%s_coskLo_bin%d.pdf",plotpath.Data(),outfile,iBin));
    
    // 
    RooPlot* framecosl = CosThetaL.frame(); 
    dataA->plotOn(framecosl,Binning(20)); 
    fA->plotOn(framecosl); 
    framecosl->SetTitle("");
    framecosl->SetMinimum(0);
    framecosl->Draw();

    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->SaveSource(TString::Format("%s/%s_cosl_bin%d.cc",plotpath.Data(),outfile,iBin));
    c->Print(TString::Format("%s/%s_cosl_bin%d.pdf",plotpath.Data(),outfile,iBin));
    
    RooPlot* framecoslUp = CosThetaL.frame(); 
    dataAUp->plotOn(framecoslUp,Binning(20)); 
    f_bkgCombLUp.plotOn(framecoslUp);
    framecoslUp->SetTitle("");
    framecoslUp->SetMinimum(0);
    framecoslUp->Draw();
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->Print(TString::Format("%s/%s_coslUp_bin%d.pdf",plotpath.Data(),outfile,iBin));

    RooPlot* framecoslLo = CosThetaL.frame(); 
    dataALo->plotOn(framecoslLo,Binning(20)); 
    f_bkgCombLLo.plotOn(framecoslLo);
    framecoslLo->SetTitle("");
    framecoslLo->SetMinimum(0);
    framecoslLo->Draw();
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->Print(TString::Format("%s/%s_coslLo_bin%d.pdf",plotpath.Data(),outfile,iBin));

    //
    RooPlot* framemass = Bmass.frame();
    dataM->plotOn(framemass,Binning(20));
    fM->plotOn(framemass);
    fM->plotOn(framemass,Components(f_bkgCombM),LineColor(2),LineWidth(2),LineStyle(2));
    framemass->SetTitle("");
    framemass->SetMinimum(0);
    framemass->Draw();

    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->SaveSource(TString::Format("%s/%s_bin%d.cc",plotpath.Data(),outfile,iBin));
    c->Print(TString::Format("%s/%s_bin%d.pdf",plotpath.Data(),outfile,iBin));
    
    if (keepParam){
        RooWorkspace *wspace = new RooWorkspace("wspace","wspace");
        Bmass.setRange(5.1,5.6);
        wspace->import(*f);

        wspace->writeToFile(TString::Format("%s/wspace_prior_bin%d.root",owspacepath.Data(),iBin),true);

        // Prepare datacard
        double val[3] = {0,0,0};
        val[0] = bkgCombLUp_c1.getVal();val[1] = bkgCombLUp_c1.getError();
        writeParam(iBin, "bkgCombLUp_c1", val , 2 , keepParam);
        val[0] = bkgCombLUp_c2.getVal();val[1] = bkgCombLUp_c2.getError();
        writeParam(iBin, "bkgCombLUp_c2", val , 2 , keepParam);
        val[0] = bkgCombLUp_c3.getVal();val[1] = bkgCombLUp_c3.getError();
        writeParam(iBin, "bkgCombLUp_c3", val , 2 , keepParam);
        val[0] = bkgCombLUp_c4.getVal();val[1] = bkgCombLUp_c4.getError();
        writeParam(iBin, "bkgCombLUp_c4", val , 2 , keepParam);
        val[0] = bkgCombLUp_c5.getVal();val[1] = bkgCombLUp_c5.getError();
        writeParam(iBin, "bkgCombLUp_c5", val , 2 , keepParam);
        val[0] = bkgCombLUp_c6.getVal();val[1] = bkgCombLUp_c6.getError();
        writeParam(iBin, "bkgCombLUp_c6", val , 2 , keepParam);
        val[0] = bkgCombLUp_c7.getVal();val[1] = bkgCombLUp_c7.getError();
        writeParam(iBin, "bkgCombLUp_c7", val , 2 , keepParam);
        val[0] = bkgCombLUp_c8.getVal();val[1] = bkgCombLUp_c8.getError();
        writeParam(iBin, "bkgCombLUp_c8", val , 2 , keepParam);
        val[0] = bkgCombKUp_c1.getVal();val[1] = bkgCombKUp_c1.getError();
        writeParam(iBin, "bkgCombKUp_c1", val , 2 , keepParam);
        val[0] = bkgCombKUp_c2.getVal();val[1] = bkgCombKUp_c2.getError();
        writeParam(iBin, "bkgCombKUp_c2", val , 2 , keepParam);
        val[0] = bkgCombKUp_c3.getVal();val[1] = bkgCombKUp_c3.getError();
        writeParam(iBin, "bkgCombKUp_c3", val , 2 , keepParam);
        val[0] = bkgCombKUp_c4.getVal();val[1] = bkgCombKUp_c4.getError();
        writeParam(iBin, "bkgCombKUp_c4", val , 2 , keepParam);
        val[0] = bkgCombKUp_c5.getVal();val[1] = bkgCombKUp_c5.getError();
        writeParam(iBin, "bkgCombKUp_c5", val , 2 , keepParam);
        val[0] = bkgCombLLo_c1.getVal();val[1] = bkgCombLLo_c1.getError();
        writeParam(iBin, "bkgCombLLo_c1", val , 2 , keepParam);
        val[0] = bkgCombLLo_c2.getVal();val[1] = bkgCombLLo_c2.getError();
        writeParam(iBin, "bkgCombLLo_c2", val , 2 , keepParam);
        val[0] = bkgCombLLo_c3.getVal();val[1] = bkgCombLLo_c3.getError();
        writeParam(iBin, "bkgCombLLo_c3", val , 2 , keepParam);
        val[0] = bkgCombLLo_c4.getVal();val[1] = bkgCombLLo_c4.getError();
        writeParam(iBin, "bkgCombLLo_c4", val , 2 , keepParam);
        val[0] = bkgCombLLo_c5.getVal();val[1] = bkgCombLLo_c5.getError();
        writeParam(iBin, "bkgCombLLo_c5", val , 2 , keepParam);
        val[0] = bkgCombLLo_c6.getVal();val[1] = bkgCombLLo_c6.getError();
        writeParam(iBin, "bkgCombLLo_c6", val , 2 , keepParam);
        val[0] = bkgCombLLo_c7.getVal();val[1] = bkgCombLLo_c7.getError();
        writeParam(iBin, "bkgCombLLo_c7", val , 2 , keepParam);
        val[0] = bkgCombLLo_c8.getVal();val[1] = bkgCombLLo_c8.getError();
        writeParam(iBin, "bkgCombLLo_c8", val , 2 , keepParam);
        val[0] = bkgCombKLo_c1.getVal();val[1] = bkgCombKLo_c1.getError();
        writeParam(iBin, "bkgCombKLo_c1", val , 2 , keepParam);
        val[0] = bkgCombKLo_c2.getVal();val[1] = bkgCombKLo_c2.getError();
        writeParam(iBin, "bkgCombKLo_c2", val , 2 , keepParam);
        val[0] = bkgCombKLo_c3.getVal();val[1] = bkgCombKLo_c3.getError();
        writeParam(iBin, "bkgCombKLo_c3", val , 2 , keepParam);
        val[0] = bkgCombKLo_c4.getVal();val[1] = bkgCombKLo_c4.getError();
        writeParam(iBin, "bkgCombKLo_c4", val , 2 , keepParam);
        val[0] = bkgCombKLo_c5.getVal();val[1] = bkgCombKLo_c5.getError();
        writeParam(iBin, "bkgCombKLo_c5", val , 2 , keepParam);
        val[0] = bkgCombM_c1.getVal();val[1] = bkgCombM_c1.getError();
        writeParam(iBin, "bkgCombM_c1", val , 2 , keepParam);
        val[0] = bkgCombM_c2.getVal();val[1] = bkgCombM_c2.getError();
        writeParam(iBin, "bkgCombM_c2", val , 2 , keepParam);
        val[0] = bkgCombM_c3.getVal();val[1] = bkgCombM_c3.getError();
        writeParam(iBin, "bkgCombM_c3", val , 2 , keepParam);
        val[0] = bkgCombM_c4.getVal();val[1] = bkgCombM_c4.getError();
        writeParam(iBin, "bkgCombM_c4", val , 2 , keepParam);
    }

    // clear
    delete t1;
    delete c;
    delete dataA;
    delete dataAUp;
    delete dataALo;
    delete dataM;
}//}}}

void angular3D_prior3(int iBin, const char outfile[] = "angular3D_prior", bool keepParam = false)// Linear combination of smooth function from upper/lower sidebands
{//{{{
    // Fit to signal simulation by YsSm+YcCm to determine Sm
    ch->SetBranchStatus("*",0);
    ch->SetBranchStatus("Bmass"         , 1);
    ch->SetBranchStatus("Mumumass"      , 1);
    ch->SetBranchStatus("Mumumasserr"   , 1);
    ch->SetBranchStatus("Kstarmass"     , 1);
    ch->SetBranchStatus("CosTheta*"     , 1);
    ch->SetBranchStatus("Q2"            , 1);
    ch->SetBranchStatus("Triggers"      , 1);
    RooRealVar Bmass("Bmass","M_{K^{*}#Mu#Mu}",4.5,6.5);// Use wider range for sideband
    RooRealVar CosThetaK("CosThetaK"     , "cos#theta_{K}"       , -1. , 1.   ) ;
    RooRealVar CosThetaL("CosThetaL"     , "cos#theta_{L}"       , -1. , 1.   ) ;
    RooRealVar Mumumass("Mumumass","M_{#mu#mu}",0.,10.);
    RooRealVar Mumumasserr("Mumumasserr","Error of M_{#mu#mu}",0.,10.);
    RooRealVar Kstarmass("Kstarmass","M_{K*}",0.,1.5);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    RooRealVar Triggers("Triggers","",0,100);

    int mumuMassWindowBin = 1+2*isCDFcut;
    if (iBin==3 || iBin==5 || isCDFcut < 0 ) mumuMassWindowBin = 0;
    RooDataSet *dataM = new RooDataSet("dataM","data",ch,RooArgSet(Q2, Bmass, CosThetaK, CosThetaL, Mumumass, Mumumasserr, Kstarmass, Triggers),TString::Format("(%s) && (%s) && (%s) && (%s)",nTriggeredPath[2],q2range[iBin],mumuMassWindow[mumuMassWindowBin],kstarMassWindow[2]),0);
    RooDataSet *dataA = new RooDataSet("dataA","data",ch,RooArgSet(Q2, Bmass, CosThetaK, CosThetaL, Mumumass, Mumumasserr, Kstarmass, Triggers),TString::Format("(%s) && (%s) && (%s) && (%s) && (%s)",nTriggeredPath[2],q2range[iBin],mumuMassWindow[mumuMassWindowBin],kstarMassWindow[1],"Bmass > 5.6 || Bmass < 5.0"),0);
    RooDataSet *dataAUp = new RooDataSet("dataAUp","data",ch,RooArgSet(Q2, Bmass, CosThetaK, CosThetaL, Mumumass, Mumumasserr, Kstarmass, Triggers),TString::Format("(%s) && (%s) && (%s) && (%s) && (%s)",nTriggeredPath[2],q2range[iBin],mumuMassWindow[mumuMassWindowBin],kstarMassWindow[1],"Bmass > 5.6"),0);
    RooDataSet *dataALo = new RooDataSet("dataALo","data",ch,RooArgSet(Q2, Bmass, CosThetaK, CosThetaL, Mumumass, Mumumasserr, Kstarmass, Triggers),TString::Format("(%s) && (%s) && (%s) && (%s) && (%s)",nTriggeredPath[2],q2range[iBin],mumuMassWindow[mumuMassWindowBin],kstarMassWindow[1],"Bmass < 5.0"),0);
    if (dataM->sumEntries() == 0 || dataAUp->sumEntries() == 0){
        printf("ERROR   : Total number of entries is 0. Exit!\n");
        return;
    }
    printf("INFO: f_bkgPeak prepared.\n");
    
    // Create combinatorial background distribution
    TString f_bkgCombM_format = "1+bkgCombM_c1*(Bmass-5)+bkgCombM_c2*exp(bkgCombM_c3*(Bmass-bkgCombM_c4))";
    RooRealVar bkgCombM_c1("bkgCombM_c1","c1",readParam(iBin,"bkgCombM_c1",0),-20,40);
    RooRealVar bkgCombM_c2("bkgCombM_c2","c2",readParam(iBin,"bkgCombM_c2",0),-20,20);
    RooRealVar bkgCombM_c3("bkgCombM_c3","c3",readParam(iBin,"bkgCombM_c3",0),-10,10);
    RooRealVar bkgCombM_c4("bkgCombM_c4","c4",readParam(iBin,"bkgCombM_c4",0),-10,10);
    RooArgSet f_bkgCombM_argset;
    f_bkgCombM_argset.add(RooArgSet(bkgCombM_c1,bkgCombM_c2,bkgCombM_c3,bkgCombM_c4));
    f_bkgCombM_argset.add(RooArgSet(Bmass));
    RooGenericPdf f_bkgCombM("f_bkgCombM","f_bkgCombM",f_bkgCombM_format,f_bkgCombM_argset);
    
    RooKeysPdf f_bkgCombKUp("f_bkgCombKUp","f_bkgCombKUp",CosThetaK,*dataAUp,RooKeysPdf::MirrorBoth);
    RooKeysPdf f_bkgCombLUp("f_bkgCombLUp","f_bkgCombLUp",CosThetaL,*dataAUp,RooKeysPdf::MirrorBoth);
    RooKeysPdf f_bkgCombKLo("f_bkgCombKLo","f_bkgCombKLo",CosThetaK,*dataALo,RooKeysPdf::MirrorBoth);
    RooKeysPdf f_bkgCombLLo("f_bkgCombLLo","f_bkgCombLLo",CosThetaL,*dataALo,RooKeysPdf::MirrorBoth);
    RooProdPdf f_bkgCombAUp("f_bkgCombAUp","f_bkgCombAUp",f_bkgCombKUp,f_bkgCombLUp);
    RooProdPdf f_bkgCombALo("f_bkgCombALo","f_bkgCombALo",f_bkgCombKLo,f_bkgCombLLo);

    RooRealVar frac_bkgCombA_UpLo("frac_bkgCombA_UpLo","frac_bkgCombA_UpLo",0.5,0.01,0.99);
    RooAddPdf  f_bkgCombA("f_bkgCombA","f_bkgCombA",f_bkgCombAUp,f_bkgCombALo,frac_bkgCombA_UpLo);
 
    RooProdPdf f_bkgComb("f_bkgComb", "f_bckComb", RooArgSet(f_bkgCombA, f_bkgCombM));
    RooRealVar nbkgCombM("nbkgCombM","nbkgCombM",100,0,1E5*scaleFactor);
    RooRealVar nbkgCombA("nbkgCombA","nbkgCombA",100,0,1E5*scaleFactor);
    RooRealVar nbkgComb ("nbkgComb" ,"nbkgComb" ,100,0,1E5*scaleFactor);
    RooRealVar nbkgComb_dummy("nbkgComb_dummy","nbkgComb_dummy",0);
    nbkgComb_dummy.setConstant(kTRUE);
    RooGenericPdf f_bkgComb_dummy ("f_bkgComb_dummy" ,"f_bkgComb_dummy" ,"1",RooArgSet(Bmass,CosThetaK,CosThetaL));
    RooGenericPdf f_bkgCombA_dummy("f_bkgCombA_dummy","f_bkgCombA_dummy","1",RooArgSet(CosThetaK,CosThetaL));
    RooGenericPdf f_bkgCombM_dummy("f_bkgCombM_dummy","f_bkgCombM_dummy","1",RooArgSet(Bmass));
    RooAddPdf *fM = new RooAddPdf("kernelM","kernelM",RooArgList( f_bkgCombM,f_bkgCombM_dummy),RooArgList(nbkgCombM,nbkgComb_dummy));
    RooAddPdf *fA = new RooAddPdf("kernelA","kernelA",RooArgList( f_bkgCombA,f_bkgCombA_dummy),RooArgList(nbkgCombA,nbkgComb_dummy));
    RooAddPdf *f  = new RooAddPdf("kernel" ,"kernel" ,RooArgList( f_bkgComb ,f_bkgComb_dummy ),RooArgList(nbkgComb ,nbkgComb_dummy));
    
    RooAbsReal *nllA = fA->createNLL(*dataA,Extended(kTRUE));// Minos and Save are unknown.
    RooMinuit minuitA(*nllA);
    minuitA.setEps(1e-16);
    minuitA.migrad();
    minuitA.migrad();
    minuitA.save();

    RooAbsReal *nllM = fM->createNLL(*dataM,Extended(kTRUE));// Minos and Save are unknown.
    RooMinuit minuitM(*nllM);
    minuitM.setEps(1e-16);
    minuitM.migrad();
    minuitM.migrad();
    minuitM.save();

    bkgCombM_c1.setConstant(kTRUE);
    bkgCombM_c2.setConstant(kTRUE);
    bkgCombM_c3.setConstant(kTRUE);
    bkgCombM_c4.setConstant(kTRUE);

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");

    RooPlot* framecosk = CosThetaK.frame(); 
    dataA->plotOn(framecosk,Binning(20)); 
    fA->plotOn(framecosk); 
    framecosk->SetTitle("");
    framecosk->SetMinimum(0);
    framecosk->Draw();
    
    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = 0;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));

    c->SaveSource(TString::Format("%s/%s_cosk_bin%d.cc",plotpath.Data(),outfile,iBin));
    c->Print(TString::Format("%s/%s_cosk_bin%d.pdf",plotpath.Data(),outfile,iBin));

    RooPlot* framecoskUp = CosThetaK.frame(); 
    dataAUp->plotOn(framecoskUp,Binning(20)); 
    f_bkgCombAUp.plotOn(framecoskUp);
    framecoskUp->SetTitle("");
    framecoskUp->SetMinimum(0);
    framecoskUp->Draw();
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->Print(TString::Format("%s/%s_coskUp_bin%d.pdf",plotpath.Data(),outfile,iBin));

    RooPlot* framecoskLo = CosThetaK.frame(); 
    dataALo->plotOn(framecoskLo,Binning(20)); 
    f_bkgCombALo.plotOn(framecoskLo);
    framecoskLo->SetTitle("");
    framecoskLo->SetMinimum(0);
    framecoskLo->Draw();
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->Print(TString::Format("%s/%s_coskLo_bin%d.pdf",plotpath.Data(),outfile,iBin));
    
    // 
    RooPlot* framecosl = CosThetaL.frame(); 
    dataA->plotOn(framecosl,Binning(20)); 
    fA->plotOn(framecosl); 
    framecosl->SetTitle("");
    framecosl->SetMinimum(0);
    framecosl->Draw();

    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->SaveSource(TString::Format("%s/%s_cosl_bin%d.cc",plotpath.Data(),outfile,iBin));
    c->Print(TString::Format("%s/%s_cosl_bin%d.pdf",plotpath.Data(),outfile,iBin));
    
    RooPlot* framecoslUp = CosThetaL.frame(); 
    dataAUp->plotOn(framecoslUp,Binning(20)); 
    f_bkgCombAUp.plotOn(framecoslUp);
    framecoslUp->SetTitle("");
    framecoslUp->SetMinimum(0);
    framecoslUp->Draw();
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->Print(TString::Format("%s/%s_coslUp_bin%d.pdf",plotpath.Data(),outfile,iBin));

    RooPlot* framecoslLo = CosThetaL.frame(); 
    dataALo->plotOn(framecoslLo,Binning(20)); 
    f_bkgCombALo.plotOn(framecoslLo);
    framecoslLo->SetTitle("");
    framecoslLo->SetMinimum(0);
    framecoslLo->Draw();
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->Print(TString::Format("%s/%s_coslLo_bin%d.pdf",plotpath.Data(),outfile,iBin));

    //
    RooPlot* framemass = Bmass.frame();
    dataM->plotOn(framemass,Binning(20));
    fM->plotOn(framemass);
    fM->plotOn(framemass,Components(f_bkgCombM),LineColor(2),LineWidth(2),LineStyle(2));
    framemass->SetTitle("");
    framemass->SetMinimum(0);
    framemass->Draw();

    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    c->Update();
    c->SaveSource(TString::Format("%s/%s_bin%d.cc",plotpath.Data(),outfile,iBin));
    c->Print(TString::Format("%s/%s_bin%d.pdf",plotpath.Data(),outfile,iBin));
    
    if (keepParam){
        RooWorkspace *wspace = new RooWorkspace("wspace","wspace");
        Bmass.setRange(5.1,5.6);
        //wspace->import(f_bkgCombA);
        //wspace->import(f_bkgCombM);
        wspace->import(*f);
        wspace->Print();

        wspace->writeToFile(TString::Format("%s/wspace_prior_bin%d.root",owspacepath.Data(),iBin),true);

        // Prepare datacard
        double val[3] = {0,0,0};
        val[0] = bkgCombM_c1.getVal();val[1] = bkgCombM_c1.getError();
        writeParam(iBin, "bkgCombM_c1", val , 2 , keepParam);
        val[0] = bkgCombM_c2.getVal();val[1] = bkgCombM_c2.getError();
        writeParam(iBin, "bkgCombM_c2", val , 2 , keepParam);
        val[0] = bkgCombM_c3.getVal();val[1] = bkgCombM_c3.getError();
        writeParam(iBin, "bkgCombM_c3", val , 2 , keepParam);
        val[0] = bkgCombM_c4.getVal();val[1] = bkgCombM_c4.getError();
        writeParam(iBin, "bkgCombM_c4", val , 2 , keepParam);
    }

    // clear
    delete t1;
    delete c;
    delete dataA;
    delete dataAUp;
    delete dataALo;
    delete dataM;
}//}}}

void angular2D_data_bin(int iBin, const char outfile[] = "angular2D_data")
{//{{{
    // Same as angular3D_bin, however, it fits to angular axes only.
    //         This function is NOT designed to determine parameters in resonance region.
    if (iBin==3 || iBin==5) return;
    
    // Read data
    ch->SetBranchStatus("*",0);
    ch->SetBranchStatus("Bmass"         , 1);
    ch->SetBranchStatus("Mumumass"      , 1);
    ch->SetBranchStatus("Mumumasserr"   , 1);
    ch->SetBranchStatus("Kstarmass"     , 1);
    ch->SetBranchStatus("CosTheta*"     , 1);
    ch->SetBranchStatus("Q2"            , 1);
    ch->SetBranchStatus("Triggers"      , 1);
    RooRealVar Bmass("Bmass","M_{K^{*}#Mu#Mu}",5.1,5.6);
    RooRealVar CosThetaK("CosThetaK"     , "cos#theta_{K}"       , -1. , 1.   ) ;
    RooRealVar CosThetaL("CosThetaL"     , "cos#theta_{L}"       , -1. , 1.   ) ;
    RooRealVar Mumumass("Mumumass","M^{#mu#mu}",0.,10.);
    RooRealVar Mumumasserr("Mumumasserr","Error of M^{#mu#mu}",0.,10.);
    RooRealVar Kstarmass("Kstarmass", "M_{K%{*}}", 0, 1.5);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    RooRealVar Triggers("Triggers","",0,100);
    RooAddition Bmass_offset("Bmass_offset","Bmass_offset",RooArgSet(Bmass,RooConst(-5)));
    RooProduct Bmass_norm("Bmass_norm","Bmass_norm",RooArgSet(Bmass_offset,RooConst(1./0.5)));
        
    // Create parameters and PDFs
    TFile *f_wspace_sigA = new TFile(TString::Format("%s/wspace_sigA_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace_sigA = (RooWorkspace*)f_wspace_sigA->Get("wspace");
    RooGenericPdf *f_sigA = 0;
    RooRealVar *fl = 0;
    RooRealVar *fs = 0;
    RooRealVar *afb = 0;
    RooRealVar *as = 0;
    if (wspace_sigA){
        f_sigA = (RooGenericPdf*)wspace_sigA->pdf("f_sigA");
        fl = (RooRealVar*)wspace_sigA->var("fl");
        fs = (RooRealVar*)wspace_sigA->var("fs");
        afb = (RooRealVar*)wspace_sigA->var("afb");
        as = (RooRealVar*)wspace_sigA->var("as");
        fl->setRange(-10,10);// unbounded fl
        fl->setVal(0.6);
        afb->setRange(-10,10);// unbounded afb
        afb->setVal(0.9);
    }else{
        printf("ERROR\t\t: Please have wsapce_sigA_bin?.root prepared.\n");
        return;
    }
        
    printf("INFO: f_sigA prepared.\n");
    
    // Create combinatorial background distribution (to be checked)
    TFile *f_wspace_comb_A = new TFile(TString::Format("%s/wspace_prior_bin%d.root",iCombBkgWspacepath.Data(),iBin));
    RooWorkspace *wspace_comb_A = (RooWorkspace*)f_wspace_comb_A->Get("wspace");
    RooGenericPdf *f_bkgCombA = 0;
    if (wspace_comb_A){
        f_bkgCombA = (RooGenericPdf*)wspace_comb_A->pdf("f_bkgCombA");
    }
    printf("INFO: f_bkgCombA prepared.\n");
    
    // Create peak background distribution(jpsi/psi2s)
    TFile *f_wspace_jpsi_M = new TFile(TString::Format("%s/wspace_YpPm_jpsi_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace_jpsi_M = (RooWorkspace*)f_wspace_jpsi_M->Get("wspace");
    RooRealVar *nbkgjpsiPeak_MC = 0;
    if (wspace_jpsi_M){
        nbkgjpsiPeak_MC = (RooRealVar*)wspace_jpsi_M->var("nbkgjpsiPeak_MC");
    }
    TFile *f_wspace_psi2s_M = new TFile(TString::Format("%s/wspace_YpPm_psi2s_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace_psi2s_M = (RooWorkspace*)f_wspace_psi2s_M->Get("wspace");
    RooRealVar *nbkgpsi2sPeak_MC = 0;
    if (wspace_psi2s_M){
        nbkgpsi2sPeak_MC = (RooRealVar*)wspace_psi2s_M->var("nbkgpsi2sPeak_MC");
    }

        // Angular distribution of peaking background
    TFile *f_wspace_jpsi_A = new TFile(TString::Format("%s/wspace_PkPl_jpsi_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace_jpsi_A = (RooWorkspace*)f_wspace_jpsi_A->Get("wspace");
    RooGenericPdf *f_bkgjpsiPeakA = 0;
    if (wspace_jpsi_A){
        f_bkgjpsiPeakA = (RooGenericPdf*)wspace_jpsi_A->pdf("f_bkgjpsiPeakA");
    }
    TFile *f_wspace_psi2s_A = new TFile(TString::Format("%s/wspace_PkPl_psi2s_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace_psi2s_A = (RooWorkspace*)f_wspace_psi2s_A->Get("wspace");
    RooGenericPdf *f_bkgpsi2sPeakA = 0;
    if (wspace_psi2s_A){
        f_bkgpsi2sPeakA = (RooGenericPdf*)wspace_psi2s_A->pdf("f_bkgpsi2sPeakA");
    }
    printf("INFO: f_bkgPeakA prepared.\n");

    // Observed spectrum = model*fullEfficiency
    RooRealVar nsig("nsig","nsig",50,-10,5E3);
    RooRealVar nbkgComb("nbkgComb","nbkgComb",100,-10,1E8);
    RooRealVar nbkgjpsiPeak("nbkgjpsiPeak","nbkgjpsiPeak",10,-10,1E7);
    RooRealVar nbkgpsi2sPeak("nbkgpsi2sPeak","nbkgpsi2sPeak",10,-10,1E7);
    RooAddPdf *fM = 0;
    RooAddPdf *fA = 0;
    RooAddPdf *f = 0;
    switch(iBin){
        case 2:
            nbkgpsi2sPeak.setVal(0.);
            nbkgpsi2sPeak.setConstant(kTRUE);
            nbkgjpsiPeak.setVal(nbkgjpsiPeak_MC->getVal()*datasetLumi[0]/datasetLumi[2]);
            nbkgjpsiPeak.setError(nbkgjpsiPeak_MC->getError()*datasetLumi[0]/datasetLumi[2]);
            fA = new RooAddPdf("kernelA","kernelA",RooArgList(*f_bkgCombA,*f_bkgjpsiPeakA,*f_sigA),RooArgList(nbkgComb,nbkgjpsiPeak,nsig));
            break;
        case 4:
            nbkgjpsiPeak.setVal(nbkgjpsiPeak_MC->getVal()*datasetLumi[0]/datasetLumi[2]);
            nbkgjpsiPeak.setError(nbkgjpsiPeak_MC->getError()*datasetLumi[0]/datasetLumi[2]);
            nbkgpsi2sPeak.setVal(nbkgpsi2sPeak_MC->getVal()*datasetLumi[0]/datasetLumi[3]);
            nbkgpsi2sPeak.setError(nbkgpsi2sPeak_MC->getError()*datasetLumi[0]/datasetLumi[3]);

            fA = new RooAddPdf("kernelA","kernelA",RooArgList(*f_bkgCombA,*f_bkgjpsiPeakA,*f_bkgpsi2sPeakA,*f_sigA),RooArgList(nbkgComb,nbkgjpsiPeak,nbkgpsi2sPeak,nsig));
            break;
        case 6:
        case 9:
            nbkgjpsiPeak.setVal(0.);
            nbkgjpsiPeak.setConstant(kTRUE);
            nbkgpsi2sPeak.setVal(nbkgpsi2sPeak_MC->getVal()*datasetLumi[0]/datasetLumi[3]);
            nbkgpsi2sPeak.setError(nbkgpsi2sPeak_MC->getError()*datasetLumi[0]/datasetLumi[3]);
            fA = new RooAddPdf("kernelA","kernelA",RooArgList(*f_bkgCombA,*f_bkgpsi2sPeakA,*f_sigA),RooArgList(nbkgComb,nbkgpsi2sPeak,nsig));
            break;
        default:
            nbkgjpsiPeak.setVal(0.);
            nbkgpsi2sPeak.setVal(0.);
            nbkgjpsiPeak.setConstant(kTRUE);
            nbkgpsi2sPeak.setConstant(kTRUE);
            fA = new RooAddPdf("kernelA","kernelA",RooArgList(*f_bkgCombA,*f_sigA),RooArgList(nbkgComb,nsig));
            break;
    }
    f = fA;

    // Extra penalty term to confine As, Fs, Fl, Afb.
    //RooRealVar t_penalty("t_penalty","t",0.01);
    //RooGenericPdf f_penaltyAfb("f_penaltyAfb","(1-TMath::Erf((afb-0.75*(1-fl))/(1.5*t_penalty*(1-fl))))*(1-TMath::Erf((-afb-0.75*(1-fl))/(1.5*t_penalty*(1-fl))))",RooArgSet(afb,fl,t_penalty));
    //RooGenericPdf f_penaltyAs("f_penaltyAfb","(1-TMath::Erf((afb-2*(1-fl)/3)/(1.5*t_penalty*(1-fl))))*(1-TMath::Erf((-afb-0.75*(1-fl))/(1.5*t_penalty*(1-fl))))",RooArgSet(afb,fl,t_penalty));
    //RooProdPdf f_penalty("f_penalty","f_penalty",f_penaltyAfb,f_penaltyAs);
    //RooProdPdf f("f","f",f_model,f_penalty);
    printf("INFO: f_penalty NOT applied.\n");

    // Gaussian constraints
    RooGaussian *gaus_nbkgjpsiPeak = 0;
    RooGaussian *gaus_nbkgpsi2sPeak = 0;
    switch(iBin){
        case 2:
            gaus_nbkgjpsiPeak = new RooGaussian("gaus_nbkgjpsiPeak","gaus_nbkgjpsiPeak",nbkgjpsiPeak,RooConst(nbkgjpsiPeak_MC->getVal()*datasetLumi[0]/datasetLumi[2]),RooConst(nbkgjpsiPeak_MC->getError()*datasetLumi[0]/datasetLumi[2]));
            break;
        case 4:
            gaus_nbkgjpsiPeak = new RooGaussian("gaus_nbkgjpsiPeak","gaus_nbkgjpsiPeak",nbkgjpsiPeak,RooConst(nbkgjpsiPeak_MC->getVal()*datasetLumi[0]/datasetLumi[2]),RooConst(nbkgjpsiPeak_MC->getError()*datasetLumi[0]/datasetLumi[2]));
            gaus_nbkgpsi2sPeak = new RooGaussian("gaus_nbkgpsi2sPeak","gaus_nbkgpsi2sPeak",nbkgpsi2sPeak,RooConst(nbkgpsi2sPeak_MC->getVal()*datasetLumi[0]/datasetLumi[3]),RooConst(nbkgpsi2sPeak_MC->getError()*datasetLumi[0]/datasetLumi[3]));
            break;
        case 6:
        case 9:
            gaus_nbkgpsi2sPeak = new RooGaussian("gaus_nbkgpsi2sPeak","gaus_nbkgpsi2sPeak",nbkgpsi2sPeak,RooConst(nbkgpsi2sPeak_MC->getVal()*datasetLumi[0]/datasetLumi[3]),RooConst(nbkgpsi2sPeak_MC->getError()*datasetLumi[0]/datasetLumi[3]));
            break;
        default:
            break;
    }
    //RooBifurGauss gaus_fs("gaus_fs","gaus_fs",fs,RooConst(readParam(iBin,"fs",0,0.00129254)),RooConst(readParam(iBin,"fs",1)),RooConst(readParam(iBin,"fs",1,0.0101371)));
    //RooBifurGauss gaus_as("gaus_as","gaus_as",as,RooConst(readParam(iBin,"as",0,-0.0975917)),RooConst(readParam(iBin,"as",1)),RooConst(readParam(iBin,"as",1,0.0049092)));
//    RooBifurGauss gaus_fs("gaus_fs","gaus_fs",fs,RooConst(0.0129254),RooConst(0.00898344),RooConst(0.0101371));// 2011 result
//    RooBifurGauss gaus_as("gaus_as","gaus_as",as,RooConst(-0.0975919),RooConst(0.00490805),RooConst(0.0049092));
    
    RooArgSet gausConstraints;
    switch(iBin){
        case 2:
            gausConstraints.add(RooArgSet(*gaus_nbkgjpsiPeak));
            break;
        case 4:
            gausConstraints.add(RooArgSet(*gaus_nbkgjpsiPeak,*gaus_nbkgpsi2sPeak));
            break;
        case 6:
        case 9:
            gausConstraints.add(RooArgSet(*gaus_nbkgpsi2sPeak));
            break;
        default:
            break;
    }
    printf("INFO: gausConstraints are settled.\n");
    
    // Get data and apply unbinned fit
    int mumuMassWindowBin = 1+2*isCDFcut;
    if (iBin==3 || iBin==5 || isCDFcut < 0) mumuMassWindowBin = 0; // no cut
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaK, CosThetaL, Bmass, Q2, Mumumass, Mumumasserr, Kstarmass, Triggers),TString::Format("(%s) && (%s) && (%s) && (%s)",nTriggeredPath[2], q2range[iBin],mumuMassWindow[mumuMassWindowBin], kstarMassWindow[1]),0);
    fs  ->setConstant(kTRUE);
    as  ->setConstant(kTRUE);
    f   ->fitTo(*data,Hesse(kFALSE));
    fs  ->setConstant(kFALSE);
    as  ->setConstant(kFALSE);
    fl  ->setConstant(kTRUE);
    afb ->setConstant(kTRUE);
    f   ->fitTo(*data,Hesse(kFALSE));
    fl  ->setConstant(kFALSE);
    afb ->setConstant(kFALSE);
    fs  ->setConstant(kTRUE);
    as  ->setConstant(kTRUE);
    f   ->fitTo(*data,Hesse(kFALSE));
    fs  ->setConstant(kFALSE);
    as  ->setConstant(kFALSE);
    fA  ->fitTo(*data,ExternalConstraints(gausConstraints));
    f   ->fitTo(*data,Extended(kTRUE),ExternalConstraints(gausConstraints),Minos(RooArgSet(*afb,*fl)),Offset(kTRUE));
    if (afb->getErrorLo()+afb->getErrorHi() == 0 || fl->getErrorLo()+fl->getErrorHi() == 0){
        f   ->fitTo(*data,Extended(kTRUE),ExternalConstraints(gausConstraints),Minos(RooArgSet(*afb,*fl)));
    }

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");
    TLatex *t1 = new TLatex();
    t1->SetNDC();
    double fixNDC = 0;
    
    // Draw projection to CosThetaK
    RooPlot* framecosk = CosThetaK.frame(); 
    data->plotOn(framecosk,Binning(20)); 
    f->plotOn(framecosk,LineColor(1)); 
    f->plotOn(framecosk,Components(*f_sigA),LineColor(4),LineWidth(2));
    f->plotOn(framecosk,Components(*f_bkgCombA),LineColor(2),LineWidth(2),LineStyle(2));
    f->plotOn(framecosk,Components(*f_bkgjpsiPeakA),LineColor(6),LineWidth(2),LineStyle(2));
    f->plotOn(framecosk,Components(*f_bkgpsi2sPeakA),LineColor(8),LineWidth(2),LineStyle(2));

    framecosk->SetTitle("");
    framecosk->SetMinimum(0);
    framecosk->Draw();
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    t1->DrawLatex(.65,.79+fixNDC,TString::Format("Y_{Signal}=%.3f",nsig.getVal()));
    t1->DrawLatex(.65,.72+fixNDC,TString::Format("Y_{J/#Psi}=%.3f",nbkgjpsiPeak.getVal()));
    t1->DrawLatex(.65,.65+fixNDC,TString::Format("Y_{#Psi'}=%.3f",nbkgpsi2sPeak.getVal()));
    t1->DrawLatex(.65,.58+fixNDC,TString::Format("Y_{Comb}=%.3f",nbkgComb.getVal()));
    t1->DrawLatex(.35,.79+fixNDC,TString::Format("F_{L}=%.3f#pm%.3f",fl->getVal(),fl->getError()));
    c->Update();
    c->SaveSource(TString::Format("%s/%s_cosk_bin%d.cc",plotpath.Data(),outfile,iBin));
    c->Print(TString::Format("%s/%s_cosk_bin%d.pdf",plotpath.Data(),outfile,iBin));
    printf("Projection on K is created.\n");

    // Draw projection to CosThetaL
    RooPlot* framecosl = CosThetaL.frame(); 
    data->plotOn(framecosl,Binning(20)); 
    f->plotOn(framecosl,LineColor(1)); 
    f->plotOn(framecosl,Components(*f_sigA),LineColor(4),LineWidth(2));
    f->plotOn(framecosl,Components(*f_bkgCombA),LineColor(2),LineWidth(2),LineStyle(2));
    f->plotOn(framecosl,Components(*f_bkgjpsiPeakA),LineColor(6),LineWidth(2),LineStyle(2));
    f->plotOn(framecosl,Components(*f_bkgpsi2sPeakA),LineColor(8),LineWidth(2),LineStyle(2));

    framecosl->SetTitle("");
    framecosl->SetMinimum(0);
    framecosl->Draw();

    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2range[iBin]));
    t1->DrawLatex(.65,.79+fixNDC,TString::Format("Y_{Signal}=%.3f",nsig.getVal()));
    t1->DrawLatex(.65,.72+fixNDC,TString::Format("Y_{J/#Psi}=%.3f",nbkgjpsiPeak.getVal()));
    t1->DrawLatex(.65,.65+fixNDC,TString::Format("Y_{#Psi'}=%.3f",nbkgpsi2sPeak.getVal()));
    t1->DrawLatex(.65,.58+fixNDC,TString::Format("Y_{Comb}=%.3f",nbkgComb.getVal()));
    t1->DrawLatex(.35,.79+fixNDC,TString::Format("A_{FB}=%.3f#pm%.3f",afb->getVal(),afb->getError()));
    c->Update();
    c->SaveSource(TString::Format("%s/%s_cosl_bin%d.cc",plotpath.Data(),outfile,iBin));
    c->Print(TString::Format("%s/%s_cosl_bin%d.pdf",plotpath.Data(),outfile,iBin));
    printf("Projection on L is created.\n");
    
    // Write result
    if (true){
        TFile *f_wspace = new TFile(TString::Format("%s/wspace_angular2D_data_bin%d.root",owspacepath.Data(),iBin),"UPDATE");
        RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
        if (!wspace){
            wspace = new RooWorkspace("wspace","wspace");
        }
        fs->setConstant(kTRUE);
        as->setConstant(kTRUE);
        fl->setConstant(kTRUE);
        afb->setConstant(kTRUE);
        nsig.setConstant(kTRUE);
        nbkgComb.setConstant(kTRUE);
        nbkgjpsiPeak.setConstant(kTRUE);
        nbkgpsi2sPeak.setConstant(kTRUE);
        wspace->import(*f);
        wspace->writeToFile(TString::Format("%s/wspace_angular2D_data_bin%d.root",owspacepath.Data(),iBin),true);
    }

    // write output
    double val[4]={0,0,0,0};
    val[0] = fl->getVal();val[1] = fl->getError();val[2]=fl->getErrorLo();val[3]=fl->getErrorHi();
    writeParam(iBin, "fl", val, 4);
    val[0] = afb->getVal();val[1] = afb->getError();val[2]=afb->getErrorLo();val[3]=afb->getErrorHi();
    writeParam(iBin, "afb",val, 4);
    val[0] = fs->getVal();val[1] = fs->getError();val[2]=fs->getErrorLo();val[3]=fs->getErrorHi();
    writeParam(iBin, "fs", val, 4);
    val[0] = as->getVal();val[1] = as->getError();val[2]=as->getErrorLo();val[3]=as->getErrorHi();
    writeParam(iBin, "as", val, 4);
    
    // clear
    delete t1;
    delete c;
    delete data;

}//}}}

void angular3D_bin(int iBin, const char outfile[] = "angular3D", double dataScaleFactor = 1, double fixFl = -1, double fixAfb = -1, double fixFs=-1, double fixAs=-1)
{//{{{
    // if fixFl or (fixAfb) > -1, the parameter will be fixed to a certain value, and the profile Likelihood is estimated.
    switchRedirectStdio(TString::Format("%s/angular3D_bin_stdout_bin%d.log",odatacardpath.Data(),iBin).Data(),"w",stdout);
    switchRedirectStdio(TString::Format("%s/angular3D_bin_stderr_bin%d.log",odatacardpath.Data(),iBin).Data(),"w",stderr);

    // Remark: You must use RooFit!! It's better in unbinned fit.
    //         Extended ML fit is adopted by Mauro, just follow!!
    //         This function is NOT designed to determine parameters in resonance region.
    printf("INFO\t\t: Processing angular3D_bin(iBin=%d)\n",iBin);
    if (iBin==3 || iBin==5) return;

    // Read data
    ch->SetBranchStatus("*",0);
    ch->SetBranchStatus("Bmass"         , 1);
    ch->SetBranchStatus("Mumumass"      , 1);
    ch->SetBranchStatus("Mumumasserr"   , 1);
    ch->SetBranchStatus("Kstarmass"     , 1);
    ch->SetBranchStatus("CosTheta*"     , 1);
    ch->SetBranchStatus("Q2"            , 1);
    ch->SetBranchStatus("Triggers"      , 1);
    RooRealVar  Bmass("Bmass","M_{K^{*}#Mu#Mu}",5.1,5.6);
    RooRealVar  CosThetaK("CosThetaK"     , "cos#theta_{K}"       , -1. , 1.   ) ;
    RooRealVar  CosThetaL("CosThetaL"     , "cos#theta_{L}"       , -1. , 1.   ) ;
    RooRealVar  Mumumass("Mumumass","M^{#mu#mu}",0.,10.);
    RooRealVar  Mumumasserr("Mumumasserr","Error of M^{#mu#mu}",0.,10.);
    RooRealVar  Kstarmass("Kstarmass", "M_{K^{*}}", 0, 1.5);
    RooRealVar  Q2("Q2","q^{2}",0.5,20.);
    RooRealVar  Triggers("Triggers","",0,100);
    RooAddition Bmass_offset("Bmass_offset","Bmass_offset",RooArgSet(Bmass,RooConst(-5)));
    RooProduct  Bmass_norm("Bmass_norm","Bmass_norm",RooArgSet(Bmass_offset,RooConst(1./0.50)));

    int mumuMassWindowBin = 1+2*isCDFcut;
    if (iBin==3 || iBin==5 || isCDFcut < 0) mumuMassWindowBin = 0; // no cut
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaK, CosThetaL, Bmass, Q2, Mumumass, Mumumasserr, Kstarmass, Triggers),TString::Format("(%s) && (%s) && (%s) && (%s)",nTriggeredPath[2], q2range[iBin],mumuMassWindow[mumuMassWindowBin],kstarMassWindow[1]),0);
    RooDataSet *dataAUp = new RooDataSet("dataAUp","data",ch,RooArgSet(Q2, Bmass, CosThetaK, CosThetaL, Mumumass, Mumumasserr, Kstarmass, Triggers),TString::Format("(%s) && (%s) && (%s) && (%s) && (%s)",nTriggeredPath[2],q2range[iBin],mumuMassWindow[mumuMassWindowBin],kstarMassWindow[1],"Bmass > 5.6"),0);
    RooDataSet *dataALo = new RooDataSet("dataALo","data",ch,RooArgSet(Q2, Bmass, CosThetaK, CosThetaL, Mumumass, Mumumasserr, Kstarmass, Triggers),TString::Format("(%s) && (%s) && (%s) && (%s) && (%s)",nTriggeredPath[2],q2range[iBin],mumuMassWindow[mumuMassWindowBin],kstarMassWindow[1],"Bmass < 5.0"),0);
    if (data->sumEntries() == 0){
        return;
    }

    // Create parameters and PDFs
    TFile *f_wspace_sigA = new TFile(TString::Format("%s/wspace_sigA_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace_sigA = (RooWorkspace*)f_wspace_sigA->Get("wspace");
    RooGenericPdf *f_sigA   = 0;
    RooRealVar *fl          = 0;
    RooRealVar *fs          = 0;
    RooRealVar *afb         = 0;
    RooRealVar *as          = 0;
    if (wspace_sigA){
        f_sigA  = (RooGenericPdf*)wspace_sigA->pdf("f_sigA");
        fl      = (RooRealVar*)wspace_sigA->var("fl");
        fs      = (RooRealVar*)wspace_sigA->var("fs");
        afb     = (RooRealVar*)wspace_sigA->var("afb");
        as      = (RooRealVar*)wspace_sigA->var("as");
        fs      ->setRange(1e-3,1.-1e-3);
        as      ->setRange(-10,10);//transformed as
        fl      ->setRange(-100,100);// unbounded fl
        afb     ->setRange(-100,100);// unbounded afb
        fs      ->setConstant(kFALSE);
        as      ->setConstant(kFALSE);
        fl      ->setConstant(kFALSE);
        afb     ->setConstant(kFALSE);
    }else{
        printf("ERROR\t\t: Please have wsapce_sigA_bin?.root prepared.\n");
        return;
    }
    printf("INFO: f_sigA prepared.\n");

        // Signal double gaussian
    TFile *f_wspace_sigM = new TFile(TString::Format("%s/wspace_Sm_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace_sigM = (RooWorkspace*)f_wspace_sigM->Get("wspace");
    RooGenericPdf *f_sigM       = 0;
    RooRealVar *sigGauss_mean   = 0;
    RooRealVar *sigGauss1_sigma = 0;
    RooRealVar *sigGauss2_sigma = 0;
    RooRealVar *sigM_frac       = 0;
    if (wspace_sigM){
        f_sigM = (RooGenericPdf*)wspace_sigM->pdf("f_sigM");
        sigGauss_mean   = (RooRealVar*)wspace_sigM->var("sigGauss_mean");
        sigGauss1_sigma = (RooRealVar*)wspace_sigM->var("sigGauss1_sigma");
        sigGauss2_sigma = (RooRealVar*)wspace_sigM->var("sigGauss2_sigma");
        sigM_frac       = (RooRealVar*)wspace_sigM->var("sigM_frac");
    }else{
        printf("ERROR\t\t: Please have wsapce_Sm_bin?.root prepared.\n");
        return;
    }
    f_sigA->Print();
    f_sigM->Print();

        // merge mass and angular distro of signal, remark that RooProdPdf is valid for PDFs have no shared variables.
    RooProdPdf f_sig("f_sig","f_sig",RooArgSet(*f_sigM,*f_sigA));
    printf("INFO: f_sig prepared.\n");
    
    // Create peak background distribution(jpsi/psi2s)
    TFile *f_wspace_jpsi_M = new TFile(TString::Format("%s/wspace_YpPm_jpsi_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace_jpsi_M = (RooWorkspace*)f_wspace_jpsi_M->Get("wspace");
    RooRealVar *nbkgjpsiPeak_MC = 0;
    RooAddPdf *f_bkgjpsiPeakM12 = 0;
    if (wspace_jpsi_M){
        nbkgjpsiPeak_MC = (RooRealVar*)wspace_jpsi_M->var("nbkgjpsiPeak_MC");
        f_bkgjpsiPeakM12 = (RooAddPdf*)wspace_jpsi_M->pdf("f_bkgjpsiPeakM12");
    }
    TFile *f_wspace_psi2s_M = new TFile(TString::Format("%s/wspace_YpPm_psi2s_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace_psi2s_M = (RooWorkspace*)f_wspace_psi2s_M->Get("wspace");
    RooRealVar *nbkgpsi2sPeak_MC = 0;
    RooAddPdf *f_bkgpsi2sPeakM12 = 0;
    if (wspace_psi2s_M){
        nbkgpsi2sPeak_MC  = (RooRealVar*)wspace_psi2s_M->var("nbkgpsi2sPeak_MC");
        f_bkgpsi2sPeakM12 = (RooAddPdf*)wspace_psi2s_M->pdf("f_bkgpsi2sPeakM12");
    }

        // Angular distribution of peaking background
    TFile *f_wspace_jpsi_A = new TFile(TString::Format("%s/wspace_PkPl_jpsi_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace_jpsi_A = (RooWorkspace*)f_wspace_jpsi_A->Get("wspace");
    RooGenericPdf *f_bkgjpsiPeakA = 0;
    if (wspace_jpsi_A){
        f_bkgjpsiPeakA = (RooGenericPdf*)wspace_jpsi_A->pdf("f_bkgjpsiPeakA");
    }
    TFile *f_wspace_psi2s_A = new TFile(TString::Format("%s/wspace_PkPl_psi2s_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace_psi2s_A = (RooWorkspace*)f_wspace_psi2s_A->Get("wspace");
    RooGenericPdf *f_bkgpsi2sPeakA = 0;
    if (wspace_psi2s_A){
        f_bkgpsi2sPeakA = (RooGenericPdf*)wspace_psi2s_A->pdf("f_bkgpsi2sPeakA");
    }

        // merge mass with angular term
    RooProdPdf *f_bkgjpsiPeak = 0;
    RooProdPdf *f_bkgpsi2sPeak = 0;
    switch(iBin){
        //case 2:
        case 10:
            f_bkgjpsiPeak  = new RooProdPdf("f_bkgjpsiPeak", "f_bkgjpsiPeak",RooArgSet(*f_bkgjpsiPeakA,*f_bkgjpsiPeakM12));
            break;
        case 11: //summaryBin[0]
            f_bkgjpsiPeak  = new RooProdPdf("f_bkgjpsiPeak", "f_bkgjpsiPeak",RooArgSet(*f_bkgjpsiPeakA,*f_bkgjpsiPeakM12));
            f_bkgpsi2sPeak = new RooProdPdf("f_bkgpsi2sPeak", "f_bkgpsi2sPeak",RooArgSet(*f_bkgpsi2sPeakA,*f_bkgpsi2sPeakM12));
            break;
        //case 6:
        //case 9:
        case 4:
            f_bkgpsi2sPeak = new RooProdPdf("f_bkgpsi2sPeak", "f_bkgpsi2sPeak",RooArgSet(*f_bkgpsi2sPeakA,*f_bkgpsi2sPeakM12));
            break;
        case 12:
        default:
            break;
    }
    printf("INFO: f_bkgPeak prepared.\n");
    
    // Create combinatorial background distribution (to be checked)
    TFile *f_wspace_comb = new TFile(TString::Format("%s/wspace_prior_bin%d.root",iCombBkgWspacepath.Data(),iBin));
    RooWorkspace *wspace_comb = (RooWorkspace*)f_wspace_comb->Get("wspace");
    RooAddPdf  *f_bkgCombA = 0;
    RooGenericPdf *f_bkgCombM = 0;
    if (wspace_comb){
        f_bkgCombA = (RooAddPdf*)wspace_comb->pdf("f_bkgCombA");
        f_bkgCombM = (RooGenericPdf*)wspace_comb->pdf("f_bkgCombM");
    }
    f_bkgCombA->Print();
    f_bkgCombM->Print();
    RooProdPdf f_bkgComb("f_bkgComb","f_bkgComb",RooArgSet(*f_bkgCombA, *f_bkgCombM));
    f_bkgComb.Print();
    printf("INFO: f_bkgComb prepared.\n");
    

    // Observed spectrum = model*fullEfficiency
    const double sigScaleFactor     = dataScaleFactor;
    const double bkgScaleFactor     = dataScaleFactor;
    const double jpsiScaleFactor    = dataScaleFactor;// for validation.
    const double psi2sScaleFactor   = dataScaleFactor;// for validation.
    RooRealVar nsig         ("nsig"          , "nsig"          , 20 *sigScaleFactor   , -10 , 5E3*sigScaleFactor   );
    RooRealVar nbkgComb     ("nbkgComb"      , "nbkgComb"      , 80 *bkgScaleFactor   , -10 , 5E5*bkgScaleFactor   );
    RooRealVar nbkgjpsiPeak ("nbkgjpsiPeak"  , "nbkgjpsiPeak"  , 10 *jpsiScaleFactor  , -10 , 5E5*jpsiScaleFactor  );
    RooRealVar nbkgpsi2sPeak("nbkgpsi2sPeak" , "nbkgpsi2sPeak" , 10 *psi2sScaleFactor , -10 , 5E5*psi2sScaleFactor );
    RooAddPdf *fM = 0;
    RooAddPdf *fA = 0;
    RooAddPdf *f  = 0;
    switch(iBin){
        //case 2:
        //case 10:
        //    nbkgpsi2sPeak.setVal(0.);
        //    nbkgpsi2sPeak.setConstant(kTRUE);
        //    nbkgjpsiPeak.setVal(nbkgjpsiPeak_MC->getVal()*datasetLumi[0]/datasetLumi[2]*jpsiScaleFactor);
        //    nbkgjpsiPeak.setError(nbkgjpsiPeak_MC->getError()*sqrt(datasetLumi[0]/datasetLumi[2]*jpsiScaleFactor));
        //    fM = new RooAddPdf("kernelM","kernelM",RooArgList(*f_bkgCombM,*f_bkgjpsiPeakM12,*f_sigM),RooArgList(nbkgComb,nbkgjpsiPeak,nsig));
        //    fA = new RooAddPdf("kernelA","kernelA",RooArgList(*f_bkgCombA,*f_bkgjpsiPeakA,*f_sigA),RooArgList(nbkgComb,nbkgjpsiPeak,nsig));
        //    f  = new RooAddPdf("kernel","kernel",RooArgList(f_bkgComb,*f_bkgjpsiPeak,f_sig),RooArgList(nbkgComb,nbkgjpsiPeak,nsig));
        //    break;
        //case 11://summaryBin[0]
        //    nbkgjpsiPeak.setVal(nbkgjpsiPeak_MC->getVal()*datasetLumi[0]/datasetLumi[2]*jpsiScaleFactor);
        //    nbkgjpsiPeak.setError(nbkgjpsiPeak_MC->getError()*sqrt(datasetLumi[0]/datasetLumi[2]*jpsiScaleFactor));
        //    nbkgpsi2sPeak.setVal(nbkgpsi2sPeak_MC->getVal()*datasetLumi[0]/datasetLumi[3]*psi2sScaleFactor);
        //    nbkgpsi2sPeak.setError(nbkgpsi2sPeak_MC->getError()*sqrt(datasetLumi[0]/datasetLumi[3]*psi2sScaleFactor));
        //    fM = new RooAddPdf("kernelM","kernelM",RooArgList(*f_bkgCombM,*f_bkgjpsiPeakM12,*f_bkgpsi2sPeakM12,*f_sigM),RooArgList(nbkgComb,nbkgjpsiPeak,nbkgpsi2sPeak,nsig));
        //    fA = new RooAddPdf("kernelA","kernelA",RooArgList(*f_bkgCombA,*f_bkgjpsiPeakA,*f_bkgpsi2sPeakA,*f_sigA),RooArgList(nbkgComb,nbkgjpsiPeak,nbkgpsi2sPeak,nsig));
        //    f  = new RooAddPdf("kernel","kernel",RooArgList(f_bkgComb,*f_bkgjpsiPeak,*f_bkgpsi2sPeak,f_sig),RooArgList(nbkgComb,nbkgjpsiPeak,nbkgpsi2sPeak,nsig));
        //    break;
        ////case 6:
        ////case 9:
        //case 4:
        //    nbkgjpsiPeak.setVal(0.);
        //    nbkgjpsiPeak.setConstant(kTRUE);
        //    nbkgpsi2sPeak.setVal(nbkgpsi2sPeak_MC->getVal()*datasetLumi[0]/datasetLumi[3]*psi2sScaleFactor);
        //    nbkgpsi2sPeak.setError(nbkgpsi2sPeak_MC->getError()*sqrt(datasetLumi[0]/datasetLumi[3]*psi2sScaleFactor));
        //    fM = new RooAddPdf("kernelM","kernelM",RooArgList(*f_bkgCombM,*f_bkgpsi2sPeakM12,*f_sigM),RooArgList(nbkgComb,nbkgpsi2sPeak,nsig));
        //    fA = new RooAddPdf("kernelA","kernelA",RooArgList(*f_bkgCombA,*f_bkgpsi2sPeakA,*f_sigA),RooArgList(nbkgComb,nbkgpsi2sPeak,nsig));
        //    f  = new RooAddPdf("kernel","kernel",RooArgList(f_bkgComb,*f_bkgpsi2sPeak,f_sig),RooArgList(nbkgComb,nbkgpsi2sPeak,nsig));
        //    break;
        default:
            nbkgjpsiPeak .setVal(0.);
            nbkgpsi2sPeak.setVal(0.);
            nbkgjpsiPeak .setConstant(kTRUE);
            nbkgpsi2sPeak.setConstant(kTRUE);
            fM = new RooAddPdf("kernelM","kernelM",RooArgList(*f_bkgCombM,*f_sigM),RooArgList(nbkgComb, nsig));
            fA = new RooAddPdf("kernelA","kernelA",RooArgList(*f_bkgCombA,*f_sigA),RooArgList(nbkgComb, nsig));
            f  = new RooAddPdf("kernel","kernel",RooArgList(f_bkgComb,f_sig),RooArgList(nbkgComb, nsig));
            break;
    }
    fM->Print();
    fA->Print();
    f->Print();

    // Extra penalty term to confine As, Fs, Fl, Afb.
    printf("INFO: f_penalty NOT applied. Instead, unbounded afb/fl is taken into account.\n");

    // Gaussian constraints
    //printf("DEBUG\t\t: preparing gaussian constraints\n");
    RooGaussian gaus_sigGauss_mean("gaus_sigGauss1_mean","gaus_sigGauss1_mean",*sigGauss_mean,RooConst(5.2789),RooConst(sigGauss_mean->getError()));
    RooGaussian gaus_sigGauss1_sigma("gaus_sigGauss1_sigma","gaus_sigGauss1_sigma",*sigGauss1_sigma,RooConst(sigGauss1_sigma->getVal()),RooConst(sigGauss1_sigma->getError()));
    RooGaussian gaus_sigGauss2_sigma("gaus_sigGauss2_sigma","gaus_sigGauss2_sigma",*sigGauss2_sigma,RooConst(sigGauss2_sigma->getVal()),RooConst(sigGauss2_sigma->getError()));
    RooGaussian gaus_sigM_frac("gaus_sigM_frac","gaus_sigM_frac",*sigM_frac,RooConst(sigM_frac->getVal()),RooConst(sigM_frac->getError()));
    RooGaussian *gaus_nbkgjpsiPeak = 0;
    RooGaussian *gaus_nbkgpsi2sPeak = 0;
    switch(iBin){
        //case 2:
        //case 10:
        //    gaus_nbkgjpsiPeak = new RooGaussian("gaus_nbkgjpsiPeak","gaus_nbkgjpsiPeak",nbkgjpsiPeak,RooConst(nbkgjpsiPeak_MC->getVal()*datasetLumi[0]/datasetLumi[2]*jpsiScaleFactor),RooConst(nbkgjpsiPeak_MC->getError()*sqrt(datasetLumi[0]/datasetLumi[2]*jpsiScaleFactor)));
        //    break;
        //case 11: //summaryBin[0]
        ////case 4:
        //    gaus_nbkgjpsiPeak = new RooGaussian("gaus_nbkgjpsiPeak","gaus_nbkgjpsiPeak",nbkgjpsiPeak,RooConst(nbkgjpsiPeak_MC->getVal()*datasetLumi[0]/datasetLumi[2]*jpsiScaleFactor),RooConst(nbkgjpsiPeak_MC->getError()*sqrt(datasetLumi[0]/datasetLumi[2]*jpsiScaleFactor)));
        //    gaus_nbkgpsi2sPeak = new RooGaussian("gaus_nbkgpsi2sPeak","gaus_nbkgpsi2sPeak",nbkgpsi2sPeak,RooConst(nbkgpsi2sPeak_MC->getVal()*datasetLumi[0]/datasetLumi[3]*psi2sScaleFactor),RooConst(nbkgpsi2sPeak_MC->getError()*sqrt(datasetLumi[0]/datasetLumi[3]*psi2sScaleFactor)));
        //    break;
        ////case 6:
        ////case 9:
        //case 4:
        //    gaus_nbkgpsi2sPeak = new RooGaussian("gaus_nbkgpsi2sPeak","gaus_nbkgpsi2sPeak",nbkgpsi2sPeak,RooConst(nbkgpsi2sPeak_MC->getVal()*datasetLumi[0]/datasetLumi[3]*psi2sScaleFactor),RooConst(nbkgpsi2sPeak_MC->getError()*sqrt(datasetLumi[0]/datasetLumi[3]*psi2sScaleFactor)));
        //    break;
        default:
            break;
    }
    
    RooArgSet gausConstraints(gaus_sigGauss_mean,gaus_sigGauss1_sigma,gaus_sigGauss2_sigma,gaus_sigM_frac);
    switch(iBin){
        ////case 2:
        //case 10:
        //    gausConstraints.add(RooArgSet(*gaus_nbkgjpsiPeak));
        //    break;
        //case 11://summaryBin[1]
        //    gausConstraints.add(RooArgSet(*gaus_nbkgjpsiPeak,*gaus_nbkgpsi2sPeak));
        //    break;
        ////case 6:
        ////case 9:
        //case 4:
        //    gausConstraints.add(RooArgSet(*gaus_nbkgpsi2sPeak));
        //    break;
        default:
            break;
    }
    printf("INFO\t: gausConstraints are settled.\n");
    
    // Profile likelihood or fixing parameters
    RooFormulaVar boundedAfb("boundedAfb" , "boundedAfb" , "1.5*(0.5-TMath::ATan(@0)/TMath::Pi())*TMath::ATan(@1)/TMath::Pi()" , RooArgList(*fl    , *afb));
    RooFormulaVar boundedFl("boundedFl"   , "boundedFl"  , "0.5+TMath::ATan(@0)/TMath::Pi()"                                   , RooArgList(*fl));
    RooFormulaVar transAs("transAs"       , "transAs"    , "(@0)/(2*TMath::Sqrt(3*@1*(1-@1)*@2)*0.89)"                         , RooArgList(*as    , *fs     , *fl));
    RooGaussian   gaus_fixUnboundedAfb("gaus_fixUnboundedAfb" , "gaus_fixUnboundedAfb" , boundedAfb , RooConst(fixAfb) , RooConst(0.01));// fix to unbounded value
    RooGaussian   gaus_fixUnboundedFl ("gaus_fixUnboundedFl"  , "gaus_fixUnboundedFl"  , boundedFl  , RooConst(fixFl ) , RooConst(0.01));// fix to unbounded value
    RooGaussian   gaus_fixOrigAs("gaus_fixOrigAs"             , "gaus_fixOrigAs"       , transAs    , RooConst(fixAs ) , RooConst(0.01));// fit to transformed value
    RooGaussian   gaus_fixFs("gaus_fixFs"                     , "gau_fixFs"            , *fs        , RooConst(fixFs)  , RooConst(0.01));
    //RooArgSet transParams(boundedAfb,boundedFl,transAs);
    //fA->addParameters(transParams);
    //fM->addParameters(transParams);
    //f ->addParameters(transParams);
    if (fixFl  != -1){
        gausConstraints.add(RooArgSet(gaus_fixUnboundedFl));
        while ( fabs(toBoundedAfb(afb->getVal(),toUnboundedFl(fixFl))) > (1-fixFl)*0.75 ){
            if (fabs(afb->getVal()) > 0.01){
                afb->setVal(afb->getVal()*0.8);
            }else{
                afb->setVal(0.);
                break;
            }
            printf("DEBUG\t: Reduce initial afb\n");
        }
        printf("DEBUG\t: Fix Fl at %f\n",fixFl);
    }
    if (fixAfb != -1){
        gausConstraints.add(RooArgSet(gaus_fixUnboundedAfb));
        while ( toBoundedFl(fl->getVal()) > 1.-fabs(fixAfb)*4./3. ){
            if (fabs(fl->getVal()) > 0.01){
                fl->setVal(fl->getVal()*0.8);
            }else{
                fl->setVal(0.);
                break;
            }
            printf("DEBUG\t: Reduce initial fl value\n");
        }
        printf("DEBUG\t: Fix Afb at %f\n",fixAfb);
    }
    if (fixAs != -1){
        gausConstraints.add(RooArgSet(gaus_fixOrigAs));
        printf("DEBUG\t: Fix As at %f\n",fixAs);
    }
    if (fixFs != -1){
        gausConstraints.add(RooArgSet(gaus_fixFs));
        printf("DEBUG\t: Fix Fs at %f\n",fixFs);
    } 
    
    // Pro-Fit procedure
    printf("\nINFO\t\t: Enter Pre-fit loop.\n\n");
    fM  ->fitTo(*data,Extended(kTRUE),Hesse(kFALSE),ExternalConstraints(gausConstraints));
    nsig    .setConstant(kTRUE);
    nbkgComb.setConstant(kTRUE);
    fs  ->setConstant(kFALSE);
    as  ->setConstant(kFALSE);
    fl  ->setConstant(kFALSE);
    afb ->setConstant(kFALSE);
    fA  ->fitTo(*data,Extended(kTRUE),Hesse(kFALSE),ExternalConstraints(gausConstraints));
    fA  ->fitTo(*data,Extended(kTRUE),Hesse(kFALSE),ExternalConstraints(gausConstraints));
    for(int iLoop=0; iLoop<5; iLoop++){
        fs  ->setConstant(kFALSE);
        as  ->setConstant(kFALSE);
        fl  ->setConstant(kTRUE);
        afb ->setConstant(kTRUE);
        f   ->fitTo(*data,Extended(kTRUE),Hesse(kFALSE),ExternalConstraints(gausConstraints));
        fl  ->setConstant(kFALSE);
        afb ->setConstant(kFALSE);
        fs  ->setConstant(kTRUE);
        as  ->setConstant(kTRUE);
        f   ->fitTo(*data,Extended(kTRUE),Hesse(kFALSE),ExternalConstraints(gausConstraints));
    }
    if ( fs->getVal() > 0.03 && fs->getVal() < 0.97){
        fs  ->setConstant(kFALSE);
        as  ->setConstant(kFALSE);
    }else{
        fs  ->setVal(0.);
        as  ->setVal(0.);
    }
    printf("\nINFO\t\t: Pre-fit finished.\n\n");

    //// Uncomment this section to drop S-wave contribution.
    //fs  ->setVal(0.);
    //fs  ->setConstant(kTRUE);
    //as  ->setVal(0.);
    //as  ->setConstant(kTRUE);

    // Fitting procedure in TMinuit
    nsig    .setConstant(kFALSE);
    nbkgComb.setConstant(kFALSE);
    sigGauss_mean  ->setConstant(kFALSE);
    sigGauss1_sigma->setConstant(kFALSE);
    sigGauss2_sigma->setConstant(kFALSE);
    sigM_frac      ->setConstant(kFALSE);
    double isMigradConverge[2] = {-1,0};
    double isMinosValid = -1;
    RooAbsReal *nll = f->createNLL(*data,Extended(kTRUE),ExternalConstraints(gausConstraints),Offset(kFALSE),NumCPU(1));// Minos and Save are unknown.
    RooMinuit minuit(*nll);
    minuit.setPrintLevel(0);// -1:None, 0: Reduced, 1: Normal

    printf("INFO\t\t: Start MIGRAD loop\n");
    for(int iLoop = 0; iLoop < 10; iLoop++){
        isMigradConverge[0] = minuit.migrad();
        printf("INFO\t\t: MIGRAD return code=%.0f\n",isMigradConverge[0]);
        if (isMigradConverge[0] == 0) break;
    }
    isMigradConverge[1] = minuit.save()->minNll();
    if (gKeepParam) {
        writeParam(iBin, "migrad", isMigradConverge);
        double val[4]={0,0,0,0};
        val[0] = fl->getVal();val[1] = fl->getError();val[2]=fl->getErrorLo();val[3]=fl->getErrorHi();
        writeParam(iBin, "fl_migrad", val, 4);
        val[0] = afb->getVal();val[1] = afb->getError();val[2]=afb->getErrorLo();val[3]=afb->getErrorHi();
        writeParam(iBin, "afb_migrad",val, 4);
        val[0] = fs->getVal();val[1] = fs->getError();val[2]=fs->getErrorLo();val[3]=fs->getErrorHi();
        writeParam(iBin, "fs_migrad", val, 4);
        val[0] = as->getVal();val[1] = as->getError();val[2]=as->getErrorLo();val[3]=as->getErrorHi();
        writeParam(iBin, "as_migrad", val, 4);
    }
    double isHesseValid = minuit.hesse();
    // Keep HESSE result as preliminary
    if (gKeepParam) {
        writeParam(iBin, "hesse", &isHesseValid, 1);
        minuit.save();
        double val[4]={0,0,0,0};
        val[0] = fl->getVal();val[1] = fl->getError();val[2]=fl->getErrorLo();val[3]=fl->getErrorHi();
        writeParam(iBin, "fl_hesse", val, 4);
        val[0] = afb->getVal();val[1] = afb->getError();val[2]=afb->getErrorLo();val[3]=afb->getErrorHi();
        writeParam(iBin, "afb_hesse",val, 4);
        val[0] = fs->getVal();val[1] = fs->getError();val[2]=fs->getErrorLo();val[3]=fs->getErrorHi();
        writeParam(iBin, "fs_hesse", val, 4);
        val[0] = as->getVal();val[1] = as->getError();val[2]=as->getErrorLo();val[3]=as->getErrorHi();
        writeParam(iBin, "as_hesse", val, 4);
    }
    printf("INFO\t\t: Start MINOS loop\n");
    for(int iLoop = 0; iLoop < 3; iLoop++){
        isMinosValid = minuit.minos(RooArgSet(*afb,*fl,nsig));
        printf("INFO\t\t: MINOS return code=%.0f\n",isMinosValid);
        if (isMinosValid == 0) break;
    }
    if (gKeepParam) {
        writeParam(iBin, "minos", &isMinosValid, 1);
    }
    minuit.save();

    // Draw the frame on the canvas
    TCanvas* c = new TCanvas("c");

    RooPlot* framemass = Bmass.frame();
    data->plotOn(framemass,Binning(20));// TGraphPainter options are allowed by DrawOption()
    f->plotOn(framemass,LineColor(1));
    f->plotOn(framemass,Components(f_sig),LineColor(4),LineWidth(2));
    f->plotOn(framemass,Components(f_bkgComb),LineColor(2),LineWidth(2),LineStyle(2));
    f->plotOn(framemass,Components(f_bkgComb),LineColor(2),LineWidth(2),LineStyle(2));
    f->plotOn(framemass,Components(f_bkgComb),LineColor(2),LineWidth(2),LineStyle(2));
    //if ( iBin == 10 || iBin == summaryBin[0] )f->plotOn(framemass,Components(*f_bkgjpsiPeak),LineColor(6),LineWidth(2),LineStyle(2));
    //if ( iBin == 4  || iBin == summaryBin[1] )f->plotOn(framemass,Components(*f_bkgpsi2sPeak),LineColor(8),LineWidth(2),LineStyle(2));

    framemass->SetTitle("");
    framemass->SetMinimum(0);
    framemass->Draw();

    TLatex *t1 = new TLatex();
    t1->SetNDC();
    t1->SetTextSize(0.04);// default: 0.05
    double fixNDC = 0;
    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2rangeLatex[iBin]));
    t1->DrawLatex(.6,.79+fixNDC,TString::Format("Y_{Signal}=%.3f#pm%.3f",nsig.getVal(),nsig.getError()));
    switch(iBin){
        ////case 2:
        //case 10:
        //    t1->DrawLatex(.6,.72+fixNDC,TString::Format("Y_{J/#Psi}=%.3f#pm%.3f",nbkgjpsiPeak.getVal(),nbkgjpsiPeak.getError()));
        //    t1->DrawLatex(.6,.65+fixNDC,TString::Format("Y_{Comb}=%.3f#pm%.3f",nbkgComb.getVal(),nbkgComb.getError()));
        //    break;
        //case 11: // summaryBin[0]
        //    t1->DrawLatex(.6,.72+fixNDC,TString::Format("Y_{J/#Psi}=%.3f#pm%.3f",nbkgjpsiPeak.getVal(),nbkgjpsiPeak.getError()));
        //    t1->DrawLatex(.6,.65+fixNDC,TString::Format("Y_{#Psi'}=%.3f#pm%.3f",nbkgpsi2sPeak.getVal(),nbkgpsi2sPeak.getError()));
        //    t1->DrawLatex(.6,.58+fixNDC,TString::Format("Y_{Comb}=%.3f#pm%.3f",nbkgComb.getVal(),nbkgComb.getError()));
        //    break;
        ////case 6:
        ////case 9:
        //case 4:
        //    t1->DrawLatex(.6,.72+fixNDC,TString::Format("Y_{#Psi'}=%.3f#pm%.3f",nbkgpsi2sPeak.getVal(),nbkgpsi2sPeak.getError()));
        //    t1->DrawLatex(.6,.65+fixNDC,TString::Format("Y_{Comb}=%.3f#pm%.3f",nbkgComb.getVal(),nbkgComb.getError()));
        //    break;
        default:
            t1->DrawLatex(.6,.72+fixNDC,TString::Format("Y_{Comb}=%.3f#pm%.3f",nbkgComb.getVal(),nbkgComb.getError()));
            break;
    }
    c->SaveSource(TString::Format("%s/%s_bin%d.cc",plotpath.Data(),outfile,iBin));
    c->Print(TString::Format("%s/%s_bin%d.pdf",plotpath.Data(),outfile,iBin));
    printf("\033[1;31mINFO\t\t\033[0m: Projection on mass created.\n");
    
    // Draw projection to CosThetaK
    RooPlot* framecosk = CosThetaK.frame(); 
    data->plotOn(framecosk,Binning(20)); 
    f->plotOn(framecosk,LineColor(1)); 
    f->plotOn(framecosk,Components(f_sig),LineColor(4),LineWidth(2));
    f->plotOn(framecosk,Components(f_bkgComb),LineColor(2),LineWidth(2),LineStyle(2));
    //if (iBin == 10 || iBin == summaryBin[0] )f->plotOn(framecosk,Components(*f_bkgjpsiPeak),LineColor(6),LineWidth(2),LineStyle(2));
    //if (iBin == 4  || iBin == summaryBin[0] )f->plotOn(framecosk,Components(*f_bkgpsi2sPeak),LineColor(8),LineWidth(2),LineStyle(2));

    framecosk->SetTitle("");
    framecosk->SetMinimum(0);
    framecosk->Draw();

    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2rangeLatex[iBin]));
    t1->DrawLatex(.35,.79+fixNDC,TString::Format("F_{L}=%.3f",toBoundedFl(fl->getVal())));
    c->Update();
    c->SaveSource(TString::Format("%s/%s_cosk_bin%d.cc",plotpath.Data(),outfile,iBin));
    c->Print(TString::Format("%s/%s_cosk_bin%d.pdf",plotpath.Data(),outfile,iBin));
    printf("\033[1;31mINFO\t\t\033[0m: Projection on K created.\n");

    // Draw projection to CosThetaL
    RooPlot* framecosl = CosThetaL.frame(); 
    data->plotOn(framecosl,Binning(20)); 
    f->plotOn(framecosl,LineColor(1)); 
    f->plotOn(framecosl,Components(f_sig),LineColor(4),LineWidth(2));
    f->plotOn(framecosl,Components(f_bkgComb),LineColor(2),LineWidth(2),LineStyle(2));
    //if ( iBin == 10 || iBin == summaryBin[0] )f->plotOn(framecosl,Components(*f_bkgjpsiPeak),LineColor(6),LineWidth(2),LineStyle(2));
    //if ( iBin == 4  || iBin == summaryBin[0] )f->plotOn(framecosl,Components(*f_bkgpsi2sPeak),LineColor(8),LineWidth(2),LineStyle(2));

    framecosl->SetTitle("");
    framecosl->SetMinimum(0);
    framecosl->Draw();

    t1->DrawLatex(.35,.86+fixNDC,TString::Format("%s",q2rangeLatex[iBin]));
    t1->DrawLatex(.35,.79+fixNDC,TString::Format("A_{FB}=%.3f",toBoundedAfb(afb->getVal(),fl->getVal())));
    c->Update();
    c->SaveSource(TString::Format("%s/%s_cosl_bin%d.cc",plotpath.Data(),outfile,iBin));
    c->Print(TString::Format("%s/%s_cosl_bin%d.pdf",plotpath.Data(),outfile,iBin));
    printf("\033[1;31mINFO\t\t\033[0m: Projection on L created.\n");
    
    // Write result
    if (gKeepParam){
        RooWorkspace *wspace = new RooWorkspace("wspace","wspace");
        sigGauss_mean  ->setConstant(kTRUE);
        sigGauss1_sigma->setConstant(kTRUE);
        sigGauss2_sigma->setConstant(kTRUE);
        sigM_frac      ->setConstant(kTRUE);
        fs->setConstant(kTRUE);
        as->setConstant(kTRUE);
        fl->setConstant(kTRUE);
        afb->setConstant(kTRUE);
        nsig.setConstant(kTRUE);
        nbkgComb.setConstant(kTRUE);
        nbkgjpsiPeak.setConstant(kTRUE);
        nbkgpsi2sPeak.setConstant(kTRUE);
        wspace->import(*f);
        wspace->import(gausConstraints);
        wspace->writeToFile(TString::Format("%s/wspace_angular3D_bin%d.root",owspacepath.Data(),iBin),true);
        
        double val[4]={0,0,0,0};
        val[0] = fl->getVal();val[1] = fl->getError();val[2]=fl->getErrorLo();val[3]=fl->getErrorHi();
        writeParam(iBin, "fl", val, 4);
        val[0] = afb->getVal();val[1] = afb->getError();val[2]=afb->getErrorLo();val[3]=afb->getErrorHi();
        writeParam(iBin, "afb",val, 4);
        val[0] = fs->getVal();val[1] = fs->getError();val[2]=fs->getErrorLo();val[3]=fs->getErrorHi();
        writeParam(iBin, "fs", val, 4);
        val[0] = as->getVal();val[1] = as->getError();val[2]=as->getErrorLo();val[3]=as->getErrorHi();
        writeParam(iBin, "as", val, 4);
    }

    switchRedirectStdio("_stdio");

    // clear
    delete t1;
    delete c;
    delete data;
    
}//}}} 

void angular3D(const char outfile[] = "angular3D")
{//{{{

    int nWorkBins = 3;
    int workBins[] = {10,4,9};
    double x[nWorkBins];
    double xerr[nWorkBins];
    double yafb[nWorkBins],yerrafbLo[nWorkBins],yerrafbHi[nWorkBins],yfl[nWorkBins],yerrflLo[nWorkBins],yerrflHi[nWorkBins];
    for(int iBin = 0; iBin < nWorkBins; iBin++){
        x[iBin] = (q2rangeup[workBins[iBin]]+q2rangedn[workBins[iBin]])/2;
        xerr[iBin] = (q2rangeup[workBins[iBin]]-q2rangedn[workBins[iBin]])/2;
    }

    // Checkout input data
    for(int ibin = 0; ibin < nWorkBins; ibin++){
        yfl      [ibin] = -100;
        yerrflLo [ibin] = 0;
        yerrflHi [ibin] = 0;
        yafb     [ibin] = -100;
        yerrafbLo[ibin] = 0;
        yerrafbHi[ibin] = 0;
        double ubd_fl   = readParam("fl",TString::Format("%s/wspace_angular3D_bin%d.root",iwspacepath.Data(),workBins[ibin]).Data())->getVal();
        double ubd_afb  = readParam("afb",TString::Format("%s/wspace_angular3D_bin%d.root",iwspacepath.Data(),workBins[ibin]).Data())->getVal();
        yfl[ibin]       = toBoundedFl(ubd_fl);
        yafb[ibin]      = toBoundedAfb(ubd_afb, ubd_fl);

        // errors from Feldman-Cousins method, no transformation needed for direct scanning.
        yerrflLo[ibin] = -1*readParam(workBins[ibin],"FCErrFl",0);
        yerrflHi[ibin] = readParam(workBins[ibin],"FCErrFl",1);
        yerrafbLo[ibin]= -1*readParam(workBins[ibin],"FCErrAfb",0);
        yerrafbHi[ibin]= readParam(workBins[ibin],"FCErrAfb",1);

        //if (yerrflLo[ibin] == 0){// errors from logL scanning
        //    yerrflLo[ibin]      = min(yfl[ibin],fabs(readParam(workBins[ibin],"scanFl",2)));
        //    yerrflHi[ibin]      = min(1-fabs(4./3*yafb[ibin])-yfl[ibin],fabs(readParam(ibin,"scanFl",3)));
        //    yerrafbLo[ibin]     = min(yafb[ibin]-0.75*(yfl[ibin]-1),fabs(readParam(workBins[ibin],"scanAfb",2)));
        //    yerrafbHi[ibin]     = min(-0.75*(yfl[ibin]-1)-yafb[ibin],fabs(readParam(workBins[ibin],"scanAfb",3)));
        //}

        //if (yerrflLo[ibin] == 0){ // errors from HESSE
        //    yerrflLo[ibin]      = fabs(toBoundedFl(readParam(workBins[ibin],"fl_hesse",0)+readParam(workBins[ibin],"fl_hesse",2))-yfl[ibin]);
        //    yerrflHi[ibin]      = fabs(toBoundedFl(readParam(workBins[ibin],"fl_hesse",0)+readParam(workBins[ibin],"fl_hesse",3))-yfl[ibin]);

        //    yerrafbLo[ibin]     = fabs(toBoundedAfb(readParam(workBins[ibin],"afb_hesse",0)+readParam(workBins[ibin],"afb_hesse",2),readParam(workBins[ibin],"fl_hesse",0))-yafb[ibin]);
        //    yerrafbHi[ibin]     = fabs(toBoundedAfb(readParam(workBins[ibin],"afb_hesse",0)+readParam(workBins[ibin],"afb_hesse",3),readParam(workBins[ibin],"fl_hesse",0))-yafb[ibin]);
        //}

        //if (yerrflLo[ibin] == 0){ // errors from MINOS
        //    yerrflLo[ibin]      = fabs(toBoundedFl(readParam(workBins[ibin],"fl",0)+readParam(workBins[ibin],"fl",2))-yfl[ibin]);
        //    yerrflHi[ibin]      = fabs(toBoundedFl(readParam(workBins[ibin],"fl",0)+readParam(workBins[ibin],"fl",3))-yfl[ibin]);
        //    //if (yerrflHi[ibin] == -1){
        //    //    yerrflHi[ibin] = 0;
        //    //}
        //    //if (yerrflLo[ibin] == -1){
        //    //    yerrflLo[ibin] = 0;
        //    //}
        //    //if (yerrflHi[ibin] == yerrflLo[ibin]){
        //    //    yerrflHi[ibin] = 0;
        //    //    yerrflLo[ibin] = 0;
        //    //}

        //    yerrafbLo[ibin]     = fabs(toBoundedAfb(readParam(workBins[ibin],"afb",0)+readParam(workBins[ibin],"afb",2),readParam(workBins[ibin],"fl",0))-yafb[ibin]);
        //    yerrafbHi[ibin]     = fabs(toBoundedAfb(readParam(workBins[ibin],"afb",0)+readParam(workBins[ibin],"afb",3),readParam(workBins[ibin],"fl",0))-yafb[ibin]);
        //    //if (yerrafbHi[ibin] == -1){
        //    //    yerrafbHi[ibin] = 0;
        //    //}
        //    //if (yerrafbLo[ibin] == -1){
        //    //    yerrafbLo[ibin] = 0;
        //    //}
        //    //if (yerrafbHi[ibin] == yerrafbLo[ibin]){
        //    //    yerrafbHi[ibin] = 0;
        //    //    yerrafbLo[ibin] = 0;
        //    //}
        //}
        printf("yafb[%d]=%6.4f + %6.4f - %6.4f\n",workBins[ibin],yafb[ibin],yerrafbHi[ibin],yerrafbLo[ibin]);
        printf("yfl [%d]=%6.4f + %6.4f - %6.4f\n",workBins[ibin],yfl[ibin],yerrflHi[ibin],yerrflLo[ibin]);
    }
    
    // Draw
    TCanvas *c = new TCanvas("c");
    TGraphAsymmErrors *g_fl  = new TGraphAsymmErrors(nWorkBins,x,yfl,xerr,xerr,yerrflLo,yerrflHi);
    g_fl->SetTitle("");
    g_fl->GetXaxis()->SetTitle("q^{2} [(GeV)^{2}]");
    g_fl->GetYaxis()->SetTitle("F_{L}");
    g_fl->GetYaxis()->SetRangeUser(0,1);
    g_fl->SetFillColor(2);
    g_fl->SetFillStyle(3001);
    g_fl->Draw("a2");
    g_fl->Draw("P TEXT");
    double work_genFl[nWorkBins];
    double work_genFlerr[nWorkBins];
    for(int ibin = 0; ibin < nWorkBins; ibin++){
        work_genFl[ibin] = genFl[workBins[ibin]];
        work_genFlerr[ibin] = genFlerr[workBins[ibin]];
    }
    TGraphAsymmErrors *gen_fl  = new TGraphAsymmErrors(nWorkBins,x,work_genFl,xerr,xerr,work_genFlerr,work_genFlerr);
    gen_fl->SetMarkerStyle(21);
    gen_fl->SetFillColor(4);
    gen_fl->SetFillStyle(3001);
    gen_fl->Draw("P2 same");
    c->SaveSource(TString::Format("%s/%s_fl.cc",plotpath.Data(),outfile));
    c->Print(TString::Format("%s/%s_fl.pdf",plotpath.Data(),outfile));
    c->Clear();

    TGraphAsymmErrors *g_afb = new TGraphAsymmErrors(nWorkBins,x,yafb,xerr,xerr,yerrafbLo,yerrafbHi);
    g_afb->SetTitle("");
    g_afb->GetXaxis()->SetTitle("q^{2} [(GeV)^{2}]");
    g_afb->GetYaxis()->SetTitle("A_{FB}");
    g_afb->GetYaxis()->SetRangeUser(-1,1);
    g_afb->SetFillColor(2);
    g_afb->SetFillStyle(3001);
    g_afb->Draw("a2");
    g_afb->Draw("P TEXT");
    double work_genAfb[nWorkBins];
    double work_genAfberr[nWorkBins];
    for(int ibin = 0; ibin < nWorkBins; ibin++){
        work_genAfb[ibin] = genAfb[workBins[ibin]];
        work_genAfberr[ibin] = genAfberr[workBins[ibin]];
    }
    TGraphAsymmErrors *gen_afb = new TGraphAsymmErrors(nWorkBins,x,work_genAfb,xerr,xerr,work_genAfberr,work_genAfberr);
    gen_afb->SetMarkerStyle(21);
    gen_afb->SetFillColor(4);
    gen_afb->SetFillStyle(3001);
    gen_afb->Draw("P2 same");
    c->SaveSource(TString::Format("%s/%s_afb.cc",plotpath.Data(),outfile));
    c->Print(TString::Format("%s/%s_afb.pdf",plotpath.Data(),outfile));

}//}}}

    // Stat Error determination tools
void scanNLL(int iBin, const char outfile[] = "scanNLL")
{//{{{
    
    // Setup data scale for validation
    double dataScaleFactor=1;
    const double sigScaleFactor   = dataScaleFactor;
    const double bkgScaleFactor   = dataScaleFactor;
    const double jpsiScaleFactor  = dataScaleFactor;
    const double psi2sScaleFactor = dataScaleFactor;
    printf("WARNING\t\t: The scale factor for background is %f\n",dataScaleFactor);

    // Get data
    ch->SetBranchStatus("*",0);
    ch->SetBranchStatus("Bmass"         , 1);
    ch->SetBranchStatus("Mumumass"      , 1);
    ch->SetBranchStatus("Mumumasserr"   , 1);
    ch->SetBranchStatus("CosTheta*"     , 1);
    ch->SetBranchStatus("Q2"            , 1);
    ch->SetBranchStatus("Triggers"      , 1);
    RooRealVar Bmass("Bmass","M_{K^{*}#Mu#Mu}",5.1,5.6);
    RooRealVar CosThetaK("CosThetaK"     , "cos#theta_{K}"       , -1. , 1.   ) ;
    RooRealVar CosThetaL("CosThetaL"     , "cos#theta_{L}"       , -1. , 1.   ) ;
    RooRealVar Mumumass("Mumumass","M^{#mu#mu}",0.,10.);
    RooRealVar Mumumasserr("Mumumasserr","Error of M^{#mu#mu}",0.,10.);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    RooRealVar Triggers("Triggers","",0,100);
    
    int mumuMassWindowBin = 1+2*isCDFcut;
    if (iBin==3 || iBin==5 || isCDFcut < 0) mumuMassWindowBin = 0; // no cut
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaK, CosThetaL, Bmass, Q2, Mumumass, Mumumasserr, Triggers),TString::Format("(%s) && (%s) && (%s)",nTriggeredPath[2], q2range[iBin],mumuMassWindow[mumuMassWindowBin]),0);
    if (data->sumEntries() == 0){
        return;
    }

    // Load f from wspace_angular3D_bin and datacard
    double minNll = readParam(iBin,"migrad",1);
    printf("INFO\t\t: minNll in datacard is %f\n",minNll);
    TFile *f_wspace = new TFile(TString::Format("%s/wspace_angular3D_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
    RooAddPdf  *f = 0;
    RooRealVar *fl = 0;
    RooRealVar *afb = 0;
    RooRealVar *bkgCombM_c1 = 0;
    RooRealVar *bkgCombM_c2 = 0;
    RooRealVar *bkgCombM_c3 = 0;
    RooRealVar *bkgCombM_c4 = 0;
    RooRealVar *fs = 0;
    RooRealVar *as = 0;
    RooRealVar *nsig = 0;
    RooRealVar *nbkgComb = 0;
    RooRealVar *nbkgjpsiPeak = 0;
    RooRealVar *nbkgpsi2sPeak = 0;
    RooRealVar *sigGauss_mean  =0;
    RooRealVar *sigGauss1_sigma=0;
    RooRealVar *sigGauss2_sigma=0;
    RooRealVar *sigM_frac      =0;
    RooArgSet  gausConstraints;
    if (wspace){
        f = (RooAddPdf*)wspace->pdf("kernel");
        fl = (RooRealVar*)wspace->var("fl");
        afb = (RooRealVar*)wspace->var("afb");
        bkgCombM_c1 = (RooRealVar*)wspace->var("bkgCombM_c1");
        bkgCombM_c2 = (RooRealVar*)wspace->var("bkgCombM_c2");
        bkgCombM_c3 = (RooRealVar*)wspace->var("bkgCombM_c3");
        bkgCombM_c4 = (RooRealVar*)wspace->var("bkgCombM_c4");
        fs          = (RooRealVar*)wspace->var("fs");
        as          = (RooRealVar*)wspace->var("as");
        nsig        = (RooRealVar*)wspace->var("nsig");
        nbkgComb    = (RooRealVar*)wspace->var("nbkgComb");
        sigGauss_mean    = (RooRealVar*)wspace->var("sigGauss_mean");
        sigGauss1_sigma  = (RooRealVar*)wspace->var("sigGauss1_sigma");
        sigGauss2_sigma  = (RooRealVar*)wspace->var("sigGauss2_sigma");
        sigM_frac        = (RooRealVar*)wspace->var("sigM_frac");
        switch(iBin){
            case 10:
                gausConstraints.add(wspace->argSet("gaus_sigGauss1_mean,gaus_sigGauss1_sigma,gaus_sigGauss2_sigma,gaus_sigM_frac,gaus_nbkgjpsiPeak"));
                nbkgjpsiPeak = (RooRealVar*)wspace->var("nbkgjpsiPeak");
                nbkgjpsiPeak->setVal(nbkgjpsiPeak->getVal()*jpsiScaleFactor);
                nbkgjpsiPeak->setError(nbkgjpsiPeak->getError()*sqrt(jpsiScaleFactor));
                nbkgjpsiPeak    ->setConstant(kFALSE);
                break;
            case 11:
                gausConstraints.add(wspace->argSet("gaus_sigGauss1_mean,gaus_sigGauss1_sigma,gaus_sigGauss2_sigma,gaus_sigM_frac,gaus_nbkgjpsiPeak,gaus_nbkgpsi2sPeak"));
                nbkgjpsiPeak = (RooRealVar*)wspace->var("nbkgjpsiPeak");
                nbkgpsi2sPeak = (RooRealVar*)wspace->var("nbkgpsi2sPeak");
                nbkgjpsiPeak->setVal(nbkgjpsiPeak->getVal()*jpsiScaleFactor);
                nbkgjpsiPeak->setError(nbkgjpsiPeak->getError()*sqrt(jpsiScaleFactor));
                nbkgpsi2sPeak->setVal(nbkgpsi2sPeak->getVal()*psi2sScaleFactor);
                nbkgpsi2sPeak->setError(nbkgpsi2sPeak->getError()*sqrt(psi2sScaleFactor));
                nbkgjpsiPeak    ->setConstant(kFALSE);
                nbkgpsi2sPeak   ->setConstant(kFALSE);
                break;
            case 4:
                gausConstraints.add(wspace->argSet("gaus_sigGauss1_mean,gaus_sigGauss1_sigma,gaus_sigGauss2_sigma,gaus_sigM_frac,gaus_nbkgpsi2sPeak"));
                nbkgpsi2sPeak = (RooRealVar*)wspace->var("nbkgpsi2sPeak");
                nbkgpsi2sPeak->setVal(nbkgpsi2sPeak->getVal()*psi2sScaleFactor);
                nbkgpsi2sPeak->setError(nbkgpsi2sPeak->getError()*sqrt(psi2sScaleFactor));
                nbkgpsi2sPeak   ->setConstant(kFALSE);
                break;
            default:
                gausConstraints = wspace->argSet("gaus_sigGauss1_mean,gaus_sigGauss1_sigma,gaus_sigGauss2_sigma,gaus_sigM_frac");
                break;
        }
        nsig            ->setRange( -10 , 5E3*sigScaleFactor   );
        nbkgComb        ->setRange( -10 , 5E5*bkgScaleFactor   );
        if ( nbkgjpsiPeak  != 0 ) nbkgjpsiPeak    ->setRange( -10 , 5E5*jpsiScaleFactor  );
        if ( nbkgpsi2sPeak != 0 ) nbkgpsi2sPeak   ->setRange( -10 , 5E5*psi2sScaleFactor );
        fl ->setRange(-100,100);// unbounded fl
        afb->setRange(-100,100);// unbounded afb
    }else{
        printf("ERROR\t\t: Please have wsapce_angular3D_bin?.root prepared.\n");
        return;
    }
    bkgCombM_c1     ->setConstant(kFALSE);
    bkgCombM_c2     ->setConstant(kFALSE);
    bkgCombM_c3     ->setConstant(kFALSE);
    bkgCombM_c4     ->setConstant(kFALSE);
    if (fs->getVal() != 0) fs->setConstant(kFALSE);
    if (as->getVal() != 0) as->setConstant(kFALSE);
    nsig            ->setConstant(kFALSE);
    nbkgComb        ->setConstant(kFALSE);
    sigGauss_mean   ->setConstant(kFALSE);
    sigGauss1_sigma ->setConstant(kFALSE);
    sigGauss2_sigma ->setConstant(kFALSE);
    sigM_frac       ->setConstant(kFALSE);
    
    fl              ->setConstant(kTRUE);
    afb             ->setConstant(kTRUE);
    printf("INFO: f prepared.\n");
    
    // Create nll and minuit
    int isMigradConverge = -1;
    int isMinosValid = -1;
    RooAbsReal *nll = f->createNLL(*data,Extended(kTRUE),ExternalConstraints(gausConstraints),NumCPU(8));
    RooMinuit minuit(*nll);
    minuit.setPrintLevel(-1);// -1:None, 0: Reduced, 1: Normal

    // Create test function to find possible domain for fl/afb. CosThetaL as x, CosThetaK as y.
    TString f2_format = "2.*[0]*y**2*(1.-x**2)+0.5*(1-[0])*(1.-y**2)*(1.+x**2)+4./3*[1]*(1-y**2)*x";
    TF2 *f2_model = new TF2("f2_model", f2_format.Data(),-1.,1.,-1.,1.);
    f2_model->SetTitle("PDF value;cos#theta_{L};cos#theta_{K}");
    TCanvas *c = new TCanvas();
    TH2F *h2_minNLLValue        = new TH2F("h2_minNLLValue","",100,-1,1,100,0,1);
    h2_minNLLValue->SetStats(false);
    h2_minNLLValue->SetXTitle("A_{FB}");
    h2_minNLLValue->SetYTitle("F_{L}");

    double minNLLInScan=1e10;
    double maxNLLInScan=-1e10;
    int fatBinWidth = 4;// resolution: *0.02(afb), *0.01(fl)
    int minNLLAfbBin=0;// keep minimum point for second iteration
    int minNLLFlBin =0;// keep minimum value for second iteration
    // remark: the order of afb/fl matters since sometimes it's trapped in local minimum.
    for( int yBin = 1; yBin <= h2_minNLLValue->GetNbinsY(); yBin++){//fl
        for( int xBin = 1; xBin <= h2_minNLLValue->GetNbinsX(); xBin++){//afb
            // Rough scan in first iteration
            if (h2_minNLLValue->GetBinContent(xBin,yBin) != 0) continue;

            // Set Afb/Fl values
            f2_model->SetParameter(0,((double)yBin-0.5)/h2_minNLLValue->GetNbinsY());
            f2_model->SetParameter(1,(2.*xBin-1.)/h2_minNLLValue->GetNbinsX()-1.);
            fl ->setVal(toUnboundedFl(((double)yBin-0.5)/h2_minNLLValue->GetNbinsY()));
            afb->setVal(toUnboundedAfb((2.*xBin-1.)/h2_minNLLValue->GetNbinsX()-1.,((double)yBin-0.5)/h2_minNLLValue->GetNbinsY()));

            // Scan positive-definite PDF
            bool isPositivePDF=true;
            if ( fabs(f2_model->GetParameter(1)) > 3./4 ){
                isPositivePDF = false;
            } else if (f2_model->Eval(1.,0.) < 0 || f2_model->Eval(-1,0) < 0){
                isPositivePDF = false;
            } else {
                for (int i = 0; i <= 100; ++i) {//cosThetaL
                    for (int j = 0; j <= 100; ++j) {//cosThetaK
                        if (f2_model->Eval(0.02*i-1,0.02*j-1) < 0){
                            isPositivePDF = false;
                            break;
                        }
                    }
                    if (!isPositivePDF) break;
                }
            }

            if(isPositivePDF){
                isMigradConverge = minuit.migrad();
                for(int iLoop = 0; iLoop < 5; iLoop++){
                    if (iLoop < 2 || iLoop > 4) isMigradConverge = minuit.minos(RooArgSet(*afb,*fl));
                    if (isMigradConverge != 0) isMigradConverge = minuit.migrad();
                    if (isMigradConverge == 0 && fabs(minuit.save()->minNll()) < fabs(minNll)*10) {// Sometimes FCN give extremely large value
                        double nllValue = minuit.save()->minNll();
                        printf("INFO\t\t: Found NLL=%f at (afb=%f, fl=%f)\n",nllValue,(2.*xBin-1)/h2_minNLLValue->GetNbinsX()-1,((double)yBin-0.5)/h2_minNLLValue->GetNbinsY());
                        for(int xFatBin=0;xFatBin<fatBinWidth;xFatBin++){
                            for(int yFatBin=0;yFatBin<fatBinWidth;yFatBin++){
                                h2_minNLLValue->SetBinContent(xBin+xFatBin,yBin+yFatBin,-nllValue);
                            }
                        }
                        if (nllValue > maxNLLInScan) {
                            maxNLLInScan = nllValue;
                        }
                        if(nllValue < minNLLInScan){
                            minNLLInScan = nllValue;
                            minNLLAfbBin = xBin;
                            minNLLFlBin  = yBin;
                        }
                        break;
                    }else{
                        printf("INFO\t\t: Loop %d, MIGRAD return code=%d at (afb=%f, fl=%f)\n",iLoop,isMigradConverge,(2.*xBin-1)/h2_minNLLValue->GetNbinsX()-1,((double)yBin-0.5)/h2_minNLLValue->GetNbinsY());
                        for(int xFatBin=0;xFatBin<fatBinWidth;xFatBin++){
                            for(int yFatBin=0;yFatBin<fatBinWidth;yFatBin++){
                                h2_minNLLValue->SetBinContent(xBin+xFatBin,yBin+yFatBin,-1);
                            }
                        }
                    }
                }
            }else{
                printf("INFO\t\t: Non-positive-definite PDF at Afb=%.3f, Fl=%.3f\n",(2.*xBin-1)/h2_minNLLValue->GetNbinsX()-1,((double)yBin-0.5)/h2_minNLLValue->GetNbinsY());
                for(int xFatBin=0;xFatBin<fatBinWidth;xFatBin++){
                    for(int yFatBin=0;yFatBin<fatBinWidth;yFatBin++){
                        h2_minNLLValue->SetBinContent(xBin+xFatBin,yBin+yFatBin,-1);
                    }
                }
            }
        }
    }

    // Draw
    h2_minNLLValue->SetMaximum(-(minNLLInScan-0.1));
    h2_minNLLValue->SetMinimum(-(minNLLInScan+2.5));
    h2_minNLLValue->Draw("COLZ0");
    c->SaveSource(TString::Format("%s/%s_bin%d.cc",plotpath.Data(),outfile,iBin));
    c->Print(TString::Format("%s/%s_bin%d.pdf",plotpath.Data(),outfile,iBin));
    h2_minNLLValue->SaveAs(TString::Format("%s/wspace_scanNLL_bin%d.root",owspacepath.Data(),iBin));

    // Scan with fine binning
    TH2F *h2_minNLLValueFine = new TH2F(*h2_minNLLValue);
    h2_minNLLValueFine->SetName("h2_minNLLValueFine");
    h2_minNLLValueFine->SetStats(false);
    h2_minNLLValueFine->SetXTitle("A_{FB}");
    h2_minNLLValueFine->SetYTitle("F_{L}");
    double minNLLInFineScan = 1e10;
    int fineScanBandWidth = 7;
    int fineScanMinXWidth = 7;
    int fineScanMinYWidth = 7;
    for( int yBin = 1; yBin <= h2_minNLLValueFine->GetNbinsY(); yBin++){//fl
        printf("WARNING\t\t: Skip fine scanning for quick test.\n"); continue;// Uncomment this line for quick check
        for( int xBin = 1; xBin <= h2_minNLLValueFine->GetNbinsX(); xBin++){//afb
            if (abs(xBin-minNLLAfbBin) > fineScanBandWidth && abs(yBin-minNLLFlBin) > fineScanBandWidth) continue;
            if (fabs( h2_minNLLValueFine->GetBinContent(xBin,yBin) + minNLLInScan ) > 1 && 
                (abs(xBin-minNLLAfbBin) > fineScanMinXWidth || abs(yBin-minNLLFlBin) > fineScanMinYWidth) ) continue;

            f2_model->SetParameter(0,(yBin-0.5)/h2_minNLLValueFine->GetNbinsY());
            f2_model->SetParameter(1,(2.*xBin-1.)/h2_minNLLValueFine->GetNbinsX()-1.);
            fl ->setVal(toUnboundedFl((yBin-0.5)/h2_minNLLValueFine->GetNbinsY()));
            afb->setVal(toUnboundedAfb((2.*xBin-1.)/h2_minNLLValueFine->GetNbinsX()-1.,((double)yBin-0.5)/h2_minNLLValueFine->GetNbinsY()));

            // Scan positive-definite PDF
            bool isPositivePDF=true;
            if ( fabs(f2_model->GetParameter(1)) > 3./4 ){
                isPositivePDF = false;
            } else if (f2_model->Eval(1.,0.) < 0 || f2_model->Eval(-1,0) < 0){
                isPositivePDF = false;
            }else{
                for (int i = 0; i <= 100; ++i) {//cosThetaL
                    for (int j = 0; j <= 100; ++j) {//cosThetaK
                        if (f2_model->Eval(0.02*i-1,0.02*j-1) < 0){
                            isPositivePDF = false;
                            break;
                        }
                    }
                    if (!isPositivePDF) break;
                }
            }

            if(isPositivePDF){
                isMigradConverge = minuit.migrad();
                for(int iLoop = 0; iLoop < 8; iLoop++){
                    if (iLoop < 2 || iLoop > 6) isMigradConverge = minuit.minos(RooArgSet(*afb,*fl));
                    if (isMigradConverge != 0) isMigradConverge = minuit.migrad();
                    if (isMigradConverge == 0 && fabs(minuit.save()->minNll()) < fabs(minNll)*10){
                        double nllValue = minuit.save()->minNll();
                        printf("INFO\t\t: Found NLL=%f at (afb=%f, fl=%f)\n",nllValue,(2.*xBin-1)/h2_minNLLValueFine->GetNbinsX()-1,(yBin-0.5)/h2_minNLLValueFine->GetNbinsY());
                        h2_minNLLValueFine->SetBinContent(xBin,yBin,-nllValue);
                        if (nllValue < minNLLInFineScan) {
                            minNLLInFineScan = nllValue;
                        }
                        break;
                    }else{
                        printf("INFO\t\t: Loop %d, MIGRAD return code=%d at (afb=%f, fl=%f)\n",iLoop,isMigradConverge,(2.*xBin-1)/h2_minNLLValueFine->GetNbinsX()-1,(yBin-0.5)/h2_minNLLValueFine->GetNbinsY());
                        h2_minNLLValueFine->SetBinContent(xBin,yBin,-1);
                    }
                }
            }else{
                printf("INFO\t\t: Non-positive-definite PDF at Afb=%.3f, Fl=%.3f\n",(2.*xBin-1)/h2_minNLLValueFine->GetNbinsX()-1,(yBin-0.5)/h2_minNLLValueFine->GetNbinsY());
                h2_minNLLValueFine->SetBinContent(xBin,yBin,-1);
            }
        }
    }

    // Draw
    h2_minNLLValueFine->SetMaximum(-(minNLLInScan-0.1));
    h2_minNLLValueFine->SetMinimum(-(minNLLInScan+2.5));
    h2_minNLLValueFine->Draw("COLZ0");

    TLine *line = new TLine();
    line->DrawLine(h2_minNLLValueFine->GetXaxis()->GetBinLowEdge(minNLLAfbBin-fineScanBandWidth),0.,h2_minNLLValueFine->GetXaxis()->GetBinLowEdge(minNLLAfbBin-fineScanBandWidth),1.);
    line->DrawLine(h2_minNLLValueFine->GetXaxis()->GetBinUpEdge(minNLLAfbBin+fineScanBandWidth),0.,h2_minNLLValueFine->GetXaxis()->GetBinUpEdge(minNLLAfbBin+fineScanBandWidth),1.);
    line->DrawLine(-1.,h2_minNLLValueFine->GetYaxis()->GetBinLowEdge(minNLLFlBin-fineScanBandWidth),1.,h2_minNLLValueFine->GetYaxis()->GetBinLowEdge(minNLLFlBin-fineScanBandWidth));
    line->DrawLine(-1.,h2_minNLLValueFine->GetYaxis()->GetBinUpEdge(minNLLFlBin+fineScanBandWidth),1.,h2_minNLLValueFine->GetYaxis()->GetBinUpEdge(minNLLFlBin+fineScanBandWidth));

    c->Update();
    c->SaveSource(TString::Format("%s/%s_fine_bin%d.cc",plotpath.Data(),outfile,iBin));
    c->Print(TString::Format("%s/%s_fine_bin%d.pdf",plotpath.Data(),outfile,iBin));
    h2_minNLLValueFine->SaveAs(TString::Format("%s/wspace_scanNLL_fine_bin%d.root",owspacepath.Data(),iBin));

    return;
}//}}} 

double getNLL(int iBin, double testFl, double testAfb, bool unboundedArg=false, const char outfile[] = "getNLL")
{//{{{

    // This function is aimed for consistancy test of minimum.

    // Setup data scale for validation
    double dataScaleFactor=20;
    const double sigScaleFactor = dataScaleFactor;
    const double bkgScaleFactor = dataScaleFactor;
    const double jpsiScaleFactor = dataScaleFactor;
    const double psi2sScaleFactor = dataScaleFactor;

    // Get data
    ch->SetBranchStatus("*",0);
    ch->SetBranchStatus("Bmass"         , 1);
    ch->SetBranchStatus("Mumumass"      , 1);
    ch->SetBranchStatus("Mumumasserr"   , 1);
    ch->SetBranchStatus("CosTheta*"     , 1);
    ch->SetBranchStatus("Q2"            , 1);
    ch->SetBranchStatus("Triggers"      , 1);
    RooRealVar Bmass("Bmass","M_{K^{*}#Mu#Mu}",5.1,5.6);
    RooRealVar CosThetaK("CosThetaK"     , "cos#theta_{K}"       , -1. , 1.   ) ;
    RooRealVar CosThetaL("CosThetaL"     , "cos#theta_{L}"       , -1. , 1.   ) ;
    RooRealVar Mumumass("Mumumass","M^{#mu#mu}",0.,10.);
    RooRealVar Mumumasserr("Mumumasserr","Error of M^{#mu#mu}",0.,10.);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    RooRealVar Triggers("Triggers","",0,100);
    int mumuMassWindowBin = 1+2*isCDFcut;
    if (iBin==3 || iBin==5 || isCDFcut < 0) mumuMassWindowBin = 0; // no cut
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaK, CosThetaL, Bmass, Q2, Mumumass, Mumumasserr, Triggers),TString::Format("(%s) && (%s) && (%s)",nTriggeredPath[2], q2range[iBin],mumuMassWindow[mumuMassWindowBin]),0);

    // Load f from wspace_angular3D_bin and datacard
    double minNll = readParam(iBin,"migrad",1);
    TFile *f_wspace = new TFile(TString::Format("%s/wspace_angular3D_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
    RooAddPdf  *f = 0;
    RooRealVar *fl = 0;
    RooRealVar *afb = 0;
    RooRealVar *bkgCombM_c1 = 0;
    RooRealVar *fs = 0;
    RooRealVar *as = 0;
    RooRealVar *nsig = 0;
    RooRealVar *nbkgComb = 0;
    RooRealVar *nbkgjpsiPeak = 0;
    RooRealVar *nbkgpsi2sPeak = 0;
    RooRealVar *sigGauss_mean  =0;
    RooRealVar *sigGauss1_sigma=0;
    RooRealVar *sigGauss2_sigma=0;
    RooRealVar *sigM_frac      =0;
    RooArgSet  gausConstraints;
    if (wspace){
        f = (RooAddPdf*)wspace->pdf("kernel");
        fl = (RooRealVar*)wspace->var("fl");
        afb = (RooRealVar*)wspace->var("afb");
        bkgCombM_c1 = (RooRealVar*)wspace->var("bkgCombM_c1");
        fs          = (RooRealVar*)wspace->var("fs");
        as          = (RooRealVar*)wspace->var("as");
        nsig        = (RooRealVar*)wspace->var("nsig");
        nbkgComb    = (RooRealVar*)wspace->var("nbkgComb");
        sigGauss_mean    = (RooRealVar*)wspace->var("sigGauss_mean");
        sigGauss1_sigma  = (RooRealVar*)wspace->var("sigGauss1_sigma");
        sigGauss2_sigma  = (RooRealVar*)wspace->var("sigGauss2_sigma");
        sigM_frac        = (RooRealVar*)wspace->var("sigM_frac");
        switch(iBin){
            case 2:
                gausConstraints.add(wspace->argSet("gaus_sigGauss1_mean,gaus_sigGauss1_sigma,gaus_sigGauss2_sigma,gaus_sigM_frac,gaus_nbkgjpsiPeak"));
                nbkgjpsiPeak = (RooRealVar*)wspace->var("nbkgjpsiPeak");
                nbkgjpsiPeak->setVal(nbkgjpsiPeak->getVal()*jpsiScaleFactor);
                nbkgjpsiPeak->setError(nbkgjpsiPeak->getError()*sqrt(jpsiScaleFactor));
                nbkgjpsiPeak    ->setConstant(kFALSE);
                break;
            case 4:
                gausConstraints.add(wspace->argSet("gaus_sigGauss1_mean,gaus_sigGauss1_sigma,gaus_sigGauss2_sigma,gaus_sigM_frac,gaus_nbkgjpsiPeak,gaus_nbkgpsi2sPeak"));
                nbkgjpsiPeak = (RooRealVar*)wspace->var("nbkgjpsiPeak");
                nbkgpsi2sPeak = (RooRealVar*)wspace->var("nbkgpsi2sPeak");
                nbkgjpsiPeak->setVal(nbkgjpsiPeak->getVal()*jpsiScaleFactor);
                nbkgjpsiPeak->setError(nbkgjpsiPeak->getError()*sqrt(jpsiScaleFactor));
                nbkgpsi2sPeak->setVal(nbkgpsi2sPeak->getVal()*psi2sScaleFactor);
                nbkgpsi2sPeak->setError(nbkgpsi2sPeak->getError()*sqrt(psi2sScaleFactor));
                nbkgjpsiPeak    ->setConstant(kFALSE);
                nbkgpsi2sPeak   ->setConstant(kFALSE);
                break;
            case 6:
            case 9:
                gausConstraints.add(wspace->argSet("gaus_sigGauss1_mean,gaus_sigGauss1_sigma,gaus_sigGauss2_sigma,gaus_sigM_frac,gaus_nbkgpsi2sPeak"));
                nbkgpsi2sPeak = (RooRealVar*)wspace->var("nbkgpsi2sPeak");
                nbkgpsi2sPeak->setVal(nbkgpsi2sPeak->getVal()*psi2sScaleFactor);
                nbkgpsi2sPeak->setError(nbkgpsi2sPeak->getError()*sqrt(psi2sScaleFactor));
                nbkgpsi2sPeak   ->setConstant(kFALSE);
                break;
            default:
                gausConstraints.add(wspace->argSet("gaus_sigGauss1_mean,gaus_sigGauss1_sigma,gaus_sigGauss2_sigma,gaus_sigM_frac"));
                break;
        }
        nsig            ->setRange( -10 , 5E3*sigScaleFactor   );
        nbkgComb        ->setRange( -10 , 5E5*bkgScaleFactor   );
        if ( nbkgjpsiPeak  != 0 ) nbkgjpsiPeak    ->setRange( -10 , 5E5*jpsiScaleFactor  );
        if ( nbkgpsi2sPeak != 0 ) nbkgpsi2sPeak   ->setRange( -10 , 5E5*psi2sScaleFactor );
        fl->setRange(-100,100);// unbounded fl
        afb->setRange(-100,100);// unbounded afb
    }else{
        printf("ERROR\t\t: Please have wsapce_angular3D_bin?.root prepared.\n");
        return -1;
    }
    bkgCombM_c1     ->setConstant(kFALSE);
    if (fs->getVal() != 0) fs              ->setConstant(kFALSE);
    if (as->getVal() != 0) as              ->setConstant(kFALSE);
    nsig            ->setConstant(kFALSE);
    nbkgComb        ->setConstant(kFALSE);
    sigGauss_mean   ->setConstant(kFALSE);
    sigGauss1_sigma ->setConstant(kFALSE);
    sigGauss2_sigma ->setConstant(kFALSE);
    sigM_frac       ->setConstant(kFALSE);

    fl              ->setConstant(kFALSE);
    afb             ->setConstant(kFALSE);
    printf("INFO: f prepared.\n");

    // Create nll and minuit
    int isMigradConverge = -1;
    int isMinosValid = -1;
    RooAbsReal *nll = f->createNLL(*data,Extended(kTRUE),ExternalConstraints(gausConstraints),NumCPU(8));
    RooMinuit minuit(*nll);

    // This section should be the simpliefied version in scanNLL
    // Set Afb/Fl values
    if (unboundedArg){
        fl ->setVal(testFl);
        afb->setVal(testAfb);
    }else{
        fl ->setVal(toUnboundedFl(testFl));
        afb->setVal(toUnboundedAfb(testAfb,testFl));
    }

    // Assuming positive-definite PDF
    if(true){
        isMigradConverge = minuit.migrad();
        for(int iLoop = 0; iLoop < 10; iLoop++){
            if (iLoop < 2 || iLoop > 7) isMigradConverge = minuit.minos(RooArgSet(*afb,*fl));
            if (isMigradConverge != 0) isMigradConverge = minuit.migrad();
            if (isMigradConverge == 0 && fabs(minuit.save()->minNll()) < fabs(minNll)*10) {// Sometimes FCN give extremely large value
                double nllValue = minuit.save()->minNll();
                printf("INFO\t\t: Found NLL=%f at (afb=%f, fl=%f)\n",nllValue,testAfb,testFl);
                if (true){// keep workspace.
                    RooWorkspace *wspace = new RooWorkspace("wspace","wspace");
                    bkgCombM_c1    ->setConstant(kTRUE);
                    sigGauss_mean  ->setConstant(kTRUE);
                    sigGauss1_sigma->setConstant(kTRUE);
                    sigGauss2_sigma->setConstant(kTRUE);
                    sigM_frac      ->setConstant(kTRUE);
                    fs             ->setConstant(kTRUE);
                    as             ->setConstant(kTRUE);
                    fl             ->setConstant(kTRUE);
                    afb            ->setConstant(kTRUE);
                    nsig           ->setConstant(kTRUE);
                    nbkgComb       ->setConstant(kTRUE);
                    nbkgjpsiPeak   ->setConstant(kTRUE);
                    nbkgpsi2sPeak  ->setConstant(kTRUE);
                    wspace         ->import(*f);
                    wspace         ->import(gausConstraints);
                    wspace->writeToFile(TString::Format("%s/wspace_getNLL_bin%d.root",owspacepath.Data(),iBin),true);
                }
                return nllValue;
            }else{
                printf("INFO\t\t: Loop %d, MIGRAD return code=%d at (afb=%f, fl=%f)\n",iLoop,isMigradConverge,testAfb,testFl);
                return 0;
            }
        }
    }else{
        printf("INFO\t\t: Non-positive-definite PDF at Afb=%.3f, Fl=%.3f\n",testAfb,testFl);
        return -2;
    }

    return -1;
}//}}} 

void getErrFromNllScan(int iBin, bool keepParam = false)
{//{{{
    TFile *fin = new TFile(TString::Format("%s/wspace_scanNLL_fine_bin%d.root",iwspacepath.Data(),iBin));
    TH2F  *h2= (TH2F*)fin->Get("h2_minNLLValueFine");

    //int fatBinWidth = 5;

    double maxLLInfo[3]={-1,0,0};//maxLL, afbBin, flBin
    for( int xBin = 1; xBin <= h2->GetNbinsX(); xBin++){//afb
        for( int yBin = 1; yBin <= h2->GetNbinsY(); yBin++){//fl
            if (h2->GetBinContent(xBin,yBin) > maxLLInfo[0]){
                maxLLInfo[0] = h2->GetBinContent(xBin,yBin);
                maxLLInfo[1] = xBin;
                maxLLInfo[2] = yBin;
            }
        }
    }

    double afb[4]={maxLLInfo[1],0,1,(double)h2->GetNbinsX()};// mean symmErr errLo errHi
    for( int xBin = 1; xBin <= h2->GetNbinsX(); xBin++){//afb
        if ( maxLLInfo[0]-h2->GetBinContent(xBin,(int)maxLLInfo[2]) < 0.5 ) continue;
        int categoryAfb= 0;
        if (xBin > maxLLInfo[1]){
            categoryAfb = 1;
        }
        if (abs(xBin-maxLLInfo[1])<abs(afb[2+categoryAfb]-maxLLInfo[1])){
            afb[2+categoryAfb] = xBin;
        }
        if (afb[2] == 1 && h2->GetBinContent(xBin,(int)maxLLInfo[2]) > 0) afb[2] = xBin; 
        if (afb[3] == 1 && h2->GetBinContent(xBin,(int)maxLLInfo[2]) > 0 && h2->GetBinContent(xBin+1,(int)maxLLInfo[2]) <= 0) afb[3] = xBin; 
    }
    
    double fl[4]={maxLLInfo[2],0,1,(double)h2->GetNbinsY()};// mean symmErr errLo errHi
    for( int yBin = 1; yBin <= h2->GetNbinsY(); yBin++){//fl
        if ( maxLLInfo[0]-h2->GetBinContent((int)maxLLInfo[1],yBin) < 0.5 ) continue;
        int categoryFl= 0;
        if (yBin > maxLLInfo[2]){
            categoryFl = 1;
        }
        if (abs(yBin-maxLLInfo[2])<abs(fl[2+categoryFl]-maxLLInfo[2])){
            fl[2+categoryFl] = yBin;
        }
        if (fl[3] == 1 && h2->GetBinContent((int)maxLLInfo[1],yBin) > 0 && h2->GetBinContent((int)maxLLInfo[1],yBin) <= 0) fl[3] = yBin; 
    }
    
    // Handle special case, keep record for hitting boundary.
    if (h2->GetBinContent(afb[2],maxLLInfo[2]) == -1){
        printf("WARNING\t\t: Afb error hit lower boundary.\n");
        afb[1]+=1;
    }
    if (h2->GetBinContent(afb[3],maxLLInfo[2]) == -1){
        printf("WARNING\t\t: Afb error hit upper boundary.\n");
        afb[1]*=-1;
    }
    if (h2->GetBinContent(maxLLInfo[1],fl[2]) == -1){
        printf("WARNING\t\t: Fl error hit lower boundary.\n");
        fl[1]+=1;
    }
    if (h2->GetBinContent(maxLLInfo[1],fl[3]) == -1){
        printf("WARNING\t\t: Fl error hit upper boundary.\n");
        fl[1]*=-1;
    }

    // Transform bin into afb/fl value
    afb[0]=h2->GetXaxis()->GetBinCenter(afb[0]);
    afb[2]=h2->GetXaxis()->GetBinCenter(afb[2])-afb[0];
    afb[3]=h2->GetXaxis()->GetBinCenter(afb[3])-afb[0];
    fl [0]=h2->GetYaxis()->GetBinCenter(fl [0]);
    fl [2]=h2->GetYaxis()->GetBinCenter(fl [2])-fl [0];
    fl [3]=h2->GetYaxis()->GetBinCenter(fl [3])-fl [0];

    printf("INFO\t\t: afb=%f%+f%+f\n",afb[0],afb[2],afb[3]);
    printf("INFO\t\t: fl=%f%+f%+f\n",fl[0],fl[2],fl[3]);

    // write to datacards
    if (true){
        writeParam(iBin,  "scanFl",  fl, 4);
        writeParam(iBin, "scanAfb", afb, 4);
    }

    return;
}//}}}

void harvestToyFitResult(int iBin, int nToy, const char toysDir[])
{//{{{
    // Buffers for checking directory/FILE
    struct stat fiBuff;

    // Setting
    TString idatacardpath_global = idatacardpath;
    TString iwspacepath_global = iwspacepath;
    TString oSummarypath = TString::Format("%s",toysDir);
    iwspacepath=TString::Format("%s",oSummarypath.Data());
    idatacardpath=TString::Format("%s",oSummarypath.Data());
    
    double nllCenter = readParam(iBin,"migrad",0);

    const int nHistBins   = 2000;

    TCanvas *canvas = new TCanvas();
    TFile *fout             = 0;
    TH1F *h_setSummaryAfb   = 0;
    TH1F *h_setSummaryFl    = 0;
    TH1F *h_setSummaryFs    = 0;
    TH1F *h_setSummaryAs    = 0;
    TH1F *h_setSummaryNLL   = 0;
    TH1F *h_setSummaryStatus= 0;

    if (stat(TString::Format("%s/setSummary.root",oSummarypath.Data()).Data(),&fiBuff) == 0){
        printf("INFO\t\t: \"setSummary.root\" exists, replace fit result!\n");
        fout = new TFile(TString::Format("%s/setSummary.root",oSummarypath.Data()),"UPDATE");
        h_setSummaryAfb = (TH1F*)fout->Get("h_setSummaryAfb");
        h_setSummaryFl  = (TH1F*)fout->Get("h_setSummaryFl" );
        h_setSummaryFs  = (TH1F*)fout->Get("h_setSummaryFs" );
        h_setSummaryAs  = (TH1F*)fout->Get("h_setSummaryAs" );
        h_setSummaryNLL = (TH1F*)fout->Get("h_setSummaryNLL");
        h_setSummaryStatus = (TH1F*)fout->Get("h_setSummaryStatus");
    }else{
        fout = new TFile(TString::Format("%s/setSummary.root",oSummarypath.Data()),"RECREATE");
        h_setSummaryAfb = new TH1F("h_setSummaryAfb", TString::Format("h_afb").Data(), nHistBins, -1.,1.);
        h_setSummaryFl  = new TH1F("h_setSummaryFl" , TString::Format("h_fl" ).Data(), nHistBins,  0.,1.);
        h_setSummaryFs  = new TH1F("h_setSummaryFs" , TString::Format("h_fs" ).Data(), nHistBins,  0.,1.);
        h_setSummaryAs  = new TH1F("h_setSummaryAs" , TString::Format("h_as" ).Data(), nHistBins, -1.,1.);
        h_setSummaryNLL = new TH1F("h_setSummaryNLL", TString::Format("h_nll").Data(), 52, nllCenter-26, nllCenter+26);//5 sigma
        h_setSummaryStatus = new TH1F("h_setSummaryStatus",TString::Format("h_status").Data(),11,-10,1);

        for(int iToy = 0; iToy<nToy; iToy++){
            iwspacepath=TString::Format("%s/set%04d",oSummarypath.Data(),iToy+1);
            idatacardpath=TString::Format("%s/set%04d",oSummarypath.Data(),iToy+1);
            double migrad = readParam(iBin,"migrad", 0, -1);
            double fcn    = readParam(iBin,"migrad", 1);
            if (migrad == 0 && fcn > -1e20){
                double afb_ubd = readParam(iBin,"afb_migrad",0);
                double fl_ubd  = readParam(iBin,"fl_migrad",0);
                double fs      = readParam(iBin,"fs_migrad",0);
                double as_tr   = readParam(iBin,"as_migrad",0);
                if (fcn > nllCenter+26){
                    fcn=nllCenter+25.5;
                }else if (fcn<nllCenter-26){
                    fcn=nllCenter-25.5;
                }

                h_setSummaryAfb ->Fill(toBoundedAfb(afb_ubd,fl_ubd));
                h_setSummaryFl  ->Fill(toBoundedFl(fl_ubd));
                h_setSummaryFs  ->Fill(fs);
                h_setSummaryAs  ->Fill(toOriginAs(as_tr,fs,fl_ubd));
                h_setSummaryNLL ->Fill(fcn);
                h_setSummaryStatus->Fill(0);
            }else if(migrad == 0 && fcn < -1e20){
                h_setSummaryStatus->Fill(-9.5);
            }else{
                h_setSummaryStatus->Fill(-1);
            }
        }
    }

    // omit the contribution at boundary.
    TF1 *afbPDFs = new TF1(TString::Format("afbPdf").Data(),"gaus(0)+[3]*exp(-0.5*(x/[4])**2)", -0.74,0.74);
    afbPDFs->SetParameter(2,0.1);
    afbPDFs->SetParameter(3,10);
    afbPDFs->SetParameter(4,0.002);
    afbPDFs->SetParLimits(0,0,100);
    afbPDFs->SetParLimits(2,0.02,1);
    afbPDFs->SetParLimits(3,0.,200);//spike at boundary
    afbPDFs->SetParLimits(4,0.001,0.02);//spike at boundary
    h_setSummaryAfb->Fit(afbPDFs,"RLQ");// Likelihood fit

    TF1 *flPDFs = new TF1(TString::Format("flPdf").Data(),"[0]*exp(-0.5*((x-[1])/[2])**2)+[3]*exp(-0.5*(x-0.005/[4])**2)+[5]*exp(-0.5*((x-0.995)/[6])**2)", 0,1);
    flPDFs->SetParameter(0,20);
    flPDFs->SetParameter(1,0.5);
    flPDFs->SetParameter(2,0.1);
    flPDFs->SetParameter(3,20);
    flPDFs->SetParameter(4,0.02);
    flPDFs->SetParameter(5,20);
    flPDFs->SetParameter(6,0.02);
    flPDFs->SetParLimits(0,5.,100);
    flPDFs->SetParLimits(1,-0.1,1.1);
    flPDFs->SetParLimits(2,0.05,0.2);
    flPDFs->SetParLimits(3,0.,500);//spike at boundary
    flPDFs->SetParLimits(4,0.001,0.02);//spike at boundary
    flPDFs->SetParLimits(5,0.,500);//spike at boundary
    flPDFs->SetParLimits(6,0.001,0.02);//spike at boundary
    h_setSummaryFl->Fit(flPDFs,"RLQ");// Likelihood fit

    h_setSummaryAfb->Draw();
    h_setSummaryAfb->SaveAs(TString::Format("%s/h_setSummaryAfb_fine.cc",oSummarypath.Data()));
    h_setSummaryFl->Draw();
    h_setSummaryFl->SaveAs(TString::Format("%s/h_setSummaryFl_fine.cc",oSummarypath.Data()));
    h_setSummaryFs->Draw();
    h_setSummaryFs->SaveAs(TString::Format("%s/h_setSummaryFs.cc",oSummarypath.Data()));
    h_setSummaryAs->Draw();
    h_setSummaryAs->SaveAs(TString::Format("%s/h_setSummaryAs.cc",oSummarypath.Data()));
    fout->Write(TString::Format("%s/setSummary.root",oSummarypath.Data()).Data(),TObject::kOverwrite);

    //  Rough binning
    h_setSummaryAfb->Rebin(40);
    h_setSummaryAfb->Fit(afbPDFs,"RLQ");// Likelihood fit
    h_setSummaryAfb->Draw();
    h_setSummaryAfb->SaveAs(TString::Format("%s/h_setSummaryAfb.cc",oSummarypath.Data()));
    canvas->Update();
    canvas->Print(TString::Format("%s/h_setSummaryAfb.pdf",oSummarypath.Data()));

    h_setSummaryFl->Rebin(40);
    h_setSummaryFl->Fit(flPDFs,"RLQ");// Likelihood fit
    h_setSummaryFl->Draw();
    h_setSummaryFl->SaveAs(TString::Format("%s/h_setSummaryFl.cc",oSummarypath.Data()));
    canvas->Update();
    canvas->Print(TString::Format("%s/h_setSummaryFl.pdf",oSummarypath.Data()));
    
    fout->Close();
    
    // restore global parameters
    idatacardpath = idatacardpath_global;
    iwspacepath   = iwspacepath_global;

    return;
}//}}}
void harvestFCFitResults(int iBin, int nToy=500)
{//{{{
    // Buffers for checking directory/FILE
    struct stat fiBuff;

    // Setting
    TString otoyspath = TString::Format("./limit/bin%d",iBin);
    const double stepSizeAfb = 0.01;
    const double stepSizeFl  = 0.01;

    // Loop over phase space
        // 1-D loop for errAfb
    for(int iAfb = 0; iAfb*stepSizeAfb < 2.; iAfb++){
        double thisAfb = iAfb*stepSizeAfb-1.+stepSizeAfb/2;
        if (stat(TString::Format("%s/afb%+04.0f",otoyspath.Data(),thisAfb*1000).Data(),&fiBuff) != 0) continue;
        if (stat(TString::Format("%s/afb%+04.0f/set0001/wspace_angular3D_bin%d.root",otoyspath.Data(),thisAfb*1000,iBin).Data(),&fiBuff) != 0) continue;
        if (stat(TString::Format("%s/setSummary.root",owspacepath.Data()).Data(),&fiBuff) == 0) continue;
        harvestToyFitResult(iBin,nToy,TString::Format("%s/afb%+04.0f",otoyspath.Data(),thisAfb*1000));
    }// Afb loop

        // 1-D loop for errFl
    for(int iFl = 0; iFl*stepSizeFl < 1.; iFl++){
        double thisFl = iFl*stepSizeFl+stepSizeFl/2;
        if (stat(TString::Format("%s/fl%+04.0f",otoyspath.Data(),thisFl*1000).Data(),&fiBuff) != 0) continue;
        if (stat(TString::Format("%s/fl%+04.0f/set0001/wspace_angular3D_bin%d.root",otoyspath.Data(),thisFl*1000,iBin).Data(),&fiBuff) != 0) continue;
        if (stat(TString::Format("%s/setSummary.root",owspacepath.Data()).Data(),&fiBuff) == 0) continue;
        harvestToyFitResult(iBin,nToy,TString::Format("%s/fl%+04.0f",otoyspath.Data(),thisFl*1000));
    }// Fl loop

    return;
}//}}}

void getCIFromTH1F(TH1F *hin, double &lowerBd, double &upperBd, double coverage=0.683)
{//{{{
    unsigned int nEntries = 0;
    for(int iBin=1; iBin <= hin->GetNbinsX(); iBin++){nEntries+=hin->GetBinContent(iBin);};

    unsigned int nMinIntervals = 0;
    float lowerBin = 1;
    float upperBin = hin->GetNbinsX();

    unsigned int minBinInterval = hin->GetNbinsX();
    unsigned int buffEntries = 0;
    for(int lBin=1; lBin <= hin->GetNbinsX(); lBin++){
        buffEntries = 0;
        for(int iBin=lBin; iBin <= hin->GetNbinsX(); iBin++){
            buffEntries += hin->GetBinContent(iBin);
            if (buffEntries > coverage*nEntries){
                if (iBin-lBin < minBinInterval){
                    minBinInterval = iBin-lBin;
                    lowerBin = lBin;
                    upperBin = iBin;
                    nMinIntervals = 1;
                }else if (iBin-lBin == minBinInterval){
                    lowerBin = lowerBin*nMinIntervals+lBin;
                    upperBin = upperBin*nMinIntervals+iBin;
                    nMinIntervals++;
                    lowerBin /= nMinIntervals;
                    upperBin /= nMinIntervals;
                }
                break;
            }
        }
    }
    upperBd = hin->GetBinCenter(floor(upperBin));
    lowerBd = hin->GetBinCenter(floor(lowerBin));
    return;
}//}}}
void getFCInterval(int iBin, int nToy=500)
{//{{{
    //  Wrong way to decide the band, should use likelihood ordering instead of minimal interval.
    //  See getFCInterval2
    // Buffers for checking directory/FILE
    struct stat fiBuff;

    // Settings
    TString otoyspath = TString::Format("./limit/bin%d",iBin);
    const double stepSizeAfb = 0.01;
    const double stepSizeFl  = 0.01;

    // Get parameters for the q2 Bin
    TFile *f_wspace = new TFile(TString::Format("%s/wspace_angular3D_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
    if (!wspace) return;
    double  fl = toBoundedFl(readParam("fl",TString::Format("%s/wspace_angular3D_bin%d.root",iwspacepath.Data(),iBin).Data())->getVal());
    double  afb = toBoundedAfb(readParam("afb",TString::Format("%s/wspace_angular3D_bin%d.root",iwspacepath.Data(),iBin).Data())->getVal(),toUnboundedFl(fl));

    // output
    int maxPoints = 1024;
    int nFlPoints = 0;
    int nAfbPoints = 0;
    double flMeasuredHi[maxPoints];
    double flMeasuredLo[maxPoints];
    double flTruth[maxPoints];
    double afbMeasuredHi[maxPoints];
    double afbMeasuredLo[maxPoints];
    double afbTruth[maxPoints];
    double outFCErrFl[2];
    double outFCErrAfb[2];
    for (int i = 0; i < maxPoints; i++) {
        flMeasuredHi [i]=-2;
        flMeasuredLo [i]=-2;
        flTruth      [i]=-2;
        afbMeasuredHi[i]=-2;
        afbMeasuredLo[i]=-2;
        afbTruth     [i]=-2;
    }

    // Loop over phase space
    double thisAfb = afb;
    double thisFl = fl;
    TFile *fin = 0;
    TH1F  *h_setSummaryAfb = 0;
    TH1F  *h_setSummaryFl  = 0;
        // loop for errAfb
    for(int iAfb = 0; iAfb*stepSizeAfb < 2.; iAfb++){
        thisAfb = iAfb*stepSizeAfb-1.+stepSizeAfb/2;
        printf("thisAfb=%+04f\n",thisAfb);

        if (stat(TString::Format("%s/afb%+04.0f_fl%+04.0f",otoyspath.Data(),thisAfb*1000,thisFl*1000).Data(),&fiBuff) != 0) continue;
        fin = new TFile(TString::Format("%s/afb%+04.0f_fl%+04.0f/setSummary.root",otoyspath.Data(),thisAfb*1000,thisFl*1000).Data());
        h_setSummaryAfb = (TH1F*)fin->Get("h_setSummaryAfb");

        getCIFromTH1F(h_setSummaryAfb,afbMeasuredLo[nAfbPoints],afbMeasuredHi[nAfbPoints]);
        afbTruth[nAfbPoints] = thisAfb;
        nAfbPoints++;

        fin->Close();
    }

        // loop for errFl
    thisAfb = afb;
    thisFl = fl;
    for(int iFl = 0; iFl*stepSizeFl < 1.; iFl++){
        thisFl = iFl*stepSizeFl+stepSizeFl/2;

        if (stat(TString::Format("%s/afb%+04.0f_fl%+04.0f",otoyspath.Data(),thisAfb*1000,thisFl*1000).Data(),&fiBuff) != 0) continue;
        fin = new TFile(TString::Format("%s/afb%+04.0f_fl%+04.0f/setSummary.root",otoyspath.Data(),thisAfb*1000,thisFl*1000).Data());
        h_setSummaryFl  = (TH1F*)fin->Get("h_setSummaryFl");

        getCIFromTH1F(h_setSummaryFl,flMeasuredLo[nFlPoints],flMeasuredHi[nFlPoints]);
        flTruth[nFlPoints] = thisFl;
        nFlPoints++;

        fin->Close();
    }// fl loop

    // Make plots, fit, and find confidence interval
    TFile *fout = new TFile(TString::Format("%s/wspace_FCConfInterval_bin%d.root",owspacepath.Data(),iBin),"RECREATE");
    TCanvas *canvas = new TCanvas();
    TLatex *latex = new TLatex();
    latex->SetNDC();
    TLine *line = new TLine();
    TGraph *g_flIntervalHi = new TGraph(nFlPoints,flMeasuredHi,flTruth);
    TGraph *g_flIntervalLo = new TGraph(nFlPoints,flMeasuredLo,flTruth);
    TGraph *g_afbIntervalHi = new TGraph(nAfbPoints,afbMeasuredHi,afbTruth);
    TGraph *g_afbIntervalLo = new TGraph(nAfbPoints,afbMeasuredLo,afbTruth);
    TF1 *f1_flIntervalHi  = new TF1("f1_flIntervalHi",  "[0]+[1]*x+[2]*x**2",0,1);
    TF1 *f1_flIntervalLo  = new TF1("f1_flIntervalLo",  "[0]+[1]*x+[2]*x**2",0,1);
    TF1 *f1_afbIntervalHi = new TF1("f1_afbIntervalHi", "[0]+[1]*x+[2]*x**2",-1,1);
    TF1 *f1_afbIntervalLo = new TF1("f1_afbIntervalLo", "[0]+[1]*x+[2]*x**2",-1,1);
    TFitResultPtr r_flIntervalHi = g_flIntervalHi->Fit("f1_flIntervalHi","S","",max(0.01,fl-0.1),min(0.99,fl+0.3));
    TFitResultPtr r_flIntervalLi = g_flIntervalLo->Fit("f1_flIntervalLo","S","",max(0.01,fl-0.3),min(0.99,fl+0.1));
    TFitResultPtr r_afbIntervalHi= g_afbIntervalHi->Fit("f1_afbIntervalHi","S","",max(0.01,afb-0.1),min(0.99,afb+0.3));
    TFitResultPtr r_afbIntervalLo= g_afbIntervalLo->Fit("f1_afbIntervalLo","S","",max(0.01,afb-0.3),min(0.99,afb+0.1));
        // Calculate the error, you must think about the case in which f1 not defined!
    outFCErrFl[0] = r_flIntervalHi.Get()  != 0 ? max(0.,f1_flIntervalHi->Eval(fl))-fl : 0.-fl;
    outFCErrFl[1] = r_flIntervalLi.Get()  != 0 ? min(flTruth[nFlPoints-1],f1_flIntervalLo->Eval(fl))-fl : flTruth[nFlPoints-1]-fl;
    outFCErrAfb[0]= r_afbIntervalHi.Get() != 0 ? max(-1*afbTruth[nAfbPoints-1],f1_afbIntervalHi->Eval(afb))-afb : -1*afbTruth[nAfbPoints-1]-afb;
    outFCErrAfb[1]= r_afbIntervalLo.Get() != 0 ? min(afbTruth[nAfbPoints-1],f1_afbIntervalLo->Eval(afb))-afb : afbTruth[nAfbPoints-1]-afb;
    if ( outFCErrFl [1] < 0 ){
        outFCErrFl [1] = 0;
    }
    if ( outFCErrAfb[0] > 0 ){
        outFCErrAfb[0] = 0;
        outFCErrFl [1] = 0;
    }
    if ( outFCErrAfb[1] < 0 ){
        outFCErrAfb[1] = 0;
        outFCErrFl [1] = 0;
    }

        // Fl plots
    g_flIntervalHi->SetTitle("");
    g_flIntervalHi->GetXaxis()->SetTitle("Measured F_{L}");
    g_flIntervalHi->GetXaxis()->SetLimits(0.,1.1*max(fl,flTruth[nFlPoints-1]));
    g_flIntervalHi->GetYaxis()->SetTitle("True F_{L}");
    g_flIntervalHi->GetYaxis()->SetRangeUser(0.,1.1*max(fl,flTruth[nFlPoints-1]));
    g_flIntervalHi->SetMarkerSize(2);
    g_flIntervalLo->SetMarkerSize(2);
    g_flIntervalHi->Draw("AP");
    g_flIntervalLo->Draw("P SAME");
    line->SetLineStyle(0.2);
    line->SetLineColor(0.2);
    line->DrawLine(fl,0,fl,flTruth[nFlPoints-1]*1.1);
    line->SetLineStyle(2);
    line->SetLineColor(1);
    line->DrawLine(g_flIntervalHi->GetXaxis()->GetXmin(),flTruth[nFlPoints-1],g_flIntervalHi->GetXaxis()->GetXmax(),flTruth[nFlPoints-1]);
    latex->DrawLatexNDC(0.15,0.75,TString::Format("F_{L}=%.3f^{%+.3f}_{%+.3f}",fl,outFCErrFl[1],outFCErrFl[0]).Data());
    canvas->Update();
    canvas->Print(TString::Format("%s/FCConfInterval_fl_bin%d.pdf",plotpath.Data(),iBin));
    // Afb plots, remark true afb could be negative
    g_afbIntervalHi->SetTitle("");
    g_afbIntervalHi->GetXaxis()->SetTitle("Measured A_{FB}");
    g_afbIntervalHi->GetXaxis()->SetLimits(-1.1*max(afb,afbTruth[nAfbPoints-1]),1.1*max(afb,afbTruth[nAfbPoints-1]));
    g_afbIntervalHi->GetYaxis()->SetTitle("True A_{FB}");
    g_afbIntervalHi->GetYaxis()->SetRangeUser(-1.1*max(afb,afbTruth[nAfbPoints-1]),1.1*max(afb,afbTruth[nAfbPoints-1]));
    g_afbIntervalHi->SetMarkerSize(0.2);
    g_afbIntervalLo->SetMarkerSize(0.2);
    g_afbIntervalHi->Draw("AP");
    g_afbIntervalLo->Draw("P SAME");
    line->SetLineStyle(2);
    line->SetLineColor(2);
    line->DrawLine(afb,-1*afbTruth[nAfbPoints-1]*1.1,afb,afbTruth[nAfbPoints-1]*1.1);
    line->SetLineStyle(2);
    line->SetLineColor(1);
    line->DrawLine(g_afbIntervalHi->GetXaxis()->GetXmin(),afbTruth[nAfbPoints-1],g_afbIntervalHi->GetXaxis()->GetXmax(),afbTruth[nAfbPoints-1]);
    line->DrawLine(g_afbIntervalHi->GetXaxis()->GetXmin(),-1*afbTruth[nAfbPoints-1],g_afbIntervalHi->GetXaxis()->GetXmax(),-1*afbTruth[nAfbPoints-1]);
    latex->DrawLatexNDC(0.15,0.75,TString::Format("A_{FB}=%.3f^{%+.3f}_{%+.3f}",afb,outFCErrAfb[1],outFCErrAfb[0]).Data());
    canvas->Update();
    canvas->Print(TString::Format("%s/FCConfInterval_afb_bin%d.pdf",plotpath.Data(),iBin));

    fout->Close();

    writeParam(iBin,"FCErrFl" ,outFCErrFl );
    writeParam(iBin,"FCErrAfb",outFCErrAfb);

    return;

}//}}}
void getFCInterval2(int iBin, int nToy=500)
{//{{{
    // Buffers for checking directory/FILE
    struct stat fiBuff;

    // Settings
    TString otoyspath = TString::Format("./limit/bin%d",iBin);
    const double targetCoverage = 0.683;
    const double stepSizeAfb = 0.01;
    const double stepSizeFl  = 0.01;

    // Get parameters for the q2 Bin
    TFile *f_wspace = new TFile(TString::Format("%s/wspace_angular3D_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
    if (!wspace) return;
    double  fl = toBoundedFl(readParam("fl",TString::Format("%s/wspace_angular3D_bin%d.root",iwspacepath.Data(),iBin).Data())->getVal());
    double  afb = toBoundedAfb(readParam("afb",TString::Format("%s/wspace_angular3D_bin%d.root",iwspacepath.Data(),iBin).Data())->getVal(),toUnboundedFl(fl));
    printf("INFO\t\t: afb=%f, fl=%f\n",afb,fl);

    // output
    int maxPoints = 1024;//no more than 1024 slices in true value space. (NOW 100,200)
    int nFlPoints = 0;
    int nAfbPoints = 0;
    double flTruth     [maxPoints];// True values for each set
    double flTrueValues[maxPoints];// True values for each toy set
    double flMeasuredHi[maxPoints];// this is the lower bound of confidence belt.
    double flMeasuredLo[maxPoints];// this is the upper bound of confidence belt.
    TF1    *flPDFs     [maxPoints];// P(measured|true);
    double afbTruth     [maxPoints];
    double afbTrueValues[maxPoints];
    double afbMeasuredHi[maxPoints];
    double afbMeasuredLo[maxPoints];
    TF1    *afbPDFs     [maxPoints];
    double outFCErrFl[2];
    double outFCErrAfb[2];
    for (int i = 0; i < maxPoints; i++) {
        flMeasuredHi [i]=-2;
        flMeasuredLo [i]=-2;
        flTruth      [i]=-2;
        flTrueValues [i]=-2;
        afbMeasuredHi[i]=-2;
        afbMeasuredLo[i]=-2;
        afbTruth     [i]=-2;
        afbTrueValues[i]=-2;
    }

    // Loop over phase space
    double thisAfb = afb;
    double thisFl = fl;
    TFile *finAfb[maxPoints];
    TH1F  *h_setSummaryAfb[maxPoints];
        // loop for errAfb
    for(int iTrueAfb = 0; iTrueAfb*stepSizeAfb < 2.; iTrueAfb++){
        thisAfb = iTrueAfb*stepSizeAfb-1.+stepSizeAfb/2;//0.01*n-1+0.005 -> 200 slices
        if ( fabs(thisAfb) < 0.75 ) afbTruth[iTrueAfb] = thisAfb;

        if (stat(TString::Format("%s/afb%+04.0f/setSummary.root",otoyspath.Data(),thisAfb*1000).Data(),&fiBuff) != 0){
            finAfb[iTrueAfb] = 0;
            continue;
        }
        finAfb[iTrueAfb] = new TFile(TString::Format("%s/afb%+04.0f/setSummary.root",otoyspath.Data(),thisAfb*1000).Data(),"READ");
        h_setSummaryAfb[iTrueAfb] = (TH1F*)finAfb[iTrueAfb]->Get("h_setSummaryAfb");
        h_setSummaryAfb[iTrueAfb]->SetName(TString::Format("h_setSummaryAfb%+04.0f",thisAfb*1000).Data());

        afbTrueValues[nAfbPoints]=thisAfb;
        nAfbPoints++;// for TGraph, connecting dots

        // Fit h_setSummaryAfb and keep the TF1
        afbPDFs[iTrueAfb] = new TF1(TString::Format("afbPdf_%+04.0f",thisAfb*1000).Data(),"gaus(0)+[3]*exp(-0.5*(x/[4])**2)", -0.75,0.75);
        afbPDFs[iTrueAfb]->SetParameter(2,0.1);
        afbPDFs[iTrueAfb]->SetParameter(3,10);
        afbPDFs[iTrueAfb]->SetParameter(4,0.002);
        afbPDFs[iTrueAfb]->SetParLimits(0,0,10);
        afbPDFs[iTrueAfb]->SetParLimits(2,0.02,1);
        afbPDFs[iTrueAfb]->SetParLimits(3,0.,200);//spike at boundary
        afbPDFs[iTrueAfb]->SetParLimits(4,0.001,0.02);//spike at boundary

        //h_setSummaryAfb[iTrueAfb]->Rebin(5);
        h_setSummaryAfb[iTrueAfb]->Scale(1./(h_setSummaryAfb[iTrueAfb]->GetSumOfWeights()*h_setSummaryAfb[iTrueAfb]->GetBinWidth(1)));//Normalize
        h_setSummaryAfb[iTrueAfb]->Fit(afbPDFs[iTrueAfb],"RLQ");// Likelihood fit
    }
    thisAfb = afb;
    thisFl = fl;
    TFile *finFl[maxPoints];
    TH1F  *h_setSummaryFl[maxPoints] ;
    for(int iTrueFl = 0; iTrueFl*stepSizeFl < 1.; iTrueFl++){
        thisFl = iTrueFl*stepSizeFl+stepSizeFl/2;
        flTruth[iTrueFl] = thisFl;

        if (stat(TString::Format("%s/fl%+04.0f/setSummary.root",otoyspath.Data(),thisFl*1000).Data(),&fiBuff) != 0){
            finFl[iTrueFl] = 0;
            continue;
        }
        finFl[iTrueFl] = new TFile(TString::Format("%s/fl%+04.0f/setSummary.root",otoyspath.Data(),thisFl*1000).Data(),"READ");
        h_setSummaryFl[iTrueFl] = (TH1F*)finFl[iTrueFl]->Get("h_setSummaryFl");
        h_setSummaryFl[iTrueFl]->SetName(TString::Format("h_setSummaryFl%+04.0f",thisFl*1000).Data());

        flTrueValues[nFlPoints]=thisFl;
        nFlPoints++;// for TGraph, connecting dots

        // Fit h_setSummaryFl and keep the TF1
        flPDFs[iTrueFl] = new TF1(TString::Format("flPdf_%+04.0f",thisFl*1000).Data(),"gaus(0)+[3]*exp(-0.5*(x/[4])**2)+[5]*exp(-0.5*((x-0.995)/[6])**2)", 0,1);
        flPDFs[iTrueFl]->SetParameter(2,0.1);
        flPDFs[iTrueFl]->SetParameter(3,10);
        flPDFs[iTrueFl]->SetParameter(4,0.002);
        flPDFs[iTrueFl]->SetParameter(5,10);
        flPDFs[iTrueFl]->SetParameter(6,0.002);
        flPDFs[iTrueFl]->SetParLimits(0,0,10);
        flPDFs[iTrueFl]->SetParLimits(2,0.02,1);
        flPDFs[iTrueFl]->SetParLimits(3,0.,200);//spike at boundary
        flPDFs[iTrueFl]->SetParLimits(4,0.001,0.02);//spike at boundary
        flPDFs[iTrueFl]->SetParLimits(5,0.,200);//spike at boundary
        flPDFs[iTrueFl]->SetParLimits(6,0.001,0.02);//spike at boundary
        
        h_setSummaryFl[iTrueFl]->Scale(1./(h_setSummaryFl[iTrueFl]->GetSumOfWeights()*h_setSummaryFl[iTrueFl]->GetBinWidth(1)));//Normalize
        h_setSummaryFl[iTrueFl]->Fit(flPDFs[iTrueFl],"RLQ");// Likelihood fit
    }
        
        // As PDFs are settled, build F&C interval by P(x_test|true)/P(x_test|best_x)
    double buffMaxLikelihoodAfb[maxPoints];// find P(x_test|best_x) as a function of measured afb.
    for (int i = 0; i < maxPoints; i++) {
        buffMaxLikelihoodAfb[i]=-2;
    }
    for(int iAfb = 0; iAfb*stepSizeAfb < 2.; iAfb++){
        if (afbTruth[iAfb] < -1) continue;
        //printf("DEBUG\t\t: Looking for P(x,mu_best) for x=thisAfb=%+04f\n",afbTruth[iAfb]);
        for(int iTrueAfb = 0; iTrueAfb*stepSizeAfb < 2.; iTrueAfb++){
            if (finAfb[iTrueAfb] == 0 ) continue;
            if ((*afbPDFs[iTrueAfb])(afbTruth[iAfb])*stepSizeAfb > buffMaxLikelihoodAfb[iAfb]) buffMaxLikelihoodAfb[iAfb] = (*afbPDFs[iTrueAfb])(afbTruth[iAfb])*stepSizeAfb;
        }
    }
    double buffMaxLikelihoodFl[maxPoints];// find P(x_test|best_x) as a function of measured fl.
    for (int i = 0; i < maxPoints; i++) {
        buffMaxLikelihoodFl[i]=-2;
    }
    for(int iFl = 0; iFl*stepSizeFl < 2.; iFl++){
        if (flTruth[iFl] < -1) continue;
        //printf("DEBUG\t\t: Looking for P(x,mu_best) for x=thisFl=%+04f\n",flTruth[iFl]);
        for(int iTrueFl = 0; iTrueFl*stepSizeFl < 1.; iTrueFl++){
            if (finFl[iTrueFl] == 0 ) continue;
            if ((*flPDFs[iTrueFl])(flTruth[iFl])*stepSizeFl > buffMaxLikelihoodFl[iFl]) buffMaxLikelihoodFl[iFl] = (*flPDFs[iTrueFl])(flTruth[iFl])*stepSizeFl;
        }
    }
    
        // Create likelihood ratio ordering and find F-C interval
    TH1F *h_LRAfb[maxPoints];
    TH1F *h_buffProbAsLRAfb = new TH1F("buffProbAsLRAfb","",2000,-1,1);
    int countAfbMeasuredHiLo = 0;
    for(int iTrueAfb = 0; iTrueAfb*stepSizeAfb < 2.; iTrueAfb++){
        if (finAfb[iTrueAfb] == 0) continue;
        //printf("DEBUG\t\t: Looking for P(x,mu)/P(x,mu_best) for mu=thisAfb=%+04f\n",afbTruth[iTrueAfb]);
        
        // Create LR ordering if finAfb exists
        h_LRAfb[iTrueAfb] = new TH1F(TString::Format("h_LRAfb_trueAfb%+04.0f",afbTruth[iTrueAfb]*1000).Data(),"",2/stepSizeAfb,-1,1);
        double buffMaxLikelihoodAfbRatio[maxPoints];
        for (int i = 0; i < maxPoints; i++) {
            buffMaxLikelihoodAfbRatio[i]=-2;
        }
        for(int iAfb = 0; iAfb*stepSizeAfb < 2.; iAfb++){
            if (afbTruth[iAfb]< -1 ) continue;
            h_LRAfb[iTrueAfb]->SetBinContent(iAfb,(*afbPDFs[iTrueAfb])(afbTruth[iAfb])*stepSizeAfb/buffMaxLikelihoodAfb[iAfb]);
        }

        // Find interval
        double buffCov = 0;
        h_buffProbAsLRAfb->Reset();
        for(int iAfb = 0; iAfb*stepSizeAfb < 2.; iAfb++){
            if (afbTruth[iAfb] < -1) continue;
            h_buffProbAsLRAfb->Fill(h_LRAfb[iTrueAfb]->GetBinContent(iAfb),(*afbPDFs[iTrueAfb])(afbTruth[iAfb])*stepSizeAfb);
        }
        for(int iBin = h_buffProbAsLRAfb->GetNbinsX(); iBin > 0; iBin--){
            buffCov += h_buffProbAsLRAfb->GetBinContent(iBin);
            if (buffCov > targetCoverage){
                buffCov = h_buffProbAsLRAfb->GetBinCenter(iBin);// serve as lower bound of LR
                break;
            }
        }
//printf("DEBUG\t\t: lowest LR for iTrueAfb=%d = %f\n", iTrueAfb, buffCov);
        for(int iAfb = 1; iAfb*stepSizeAfb < 2.; iAfb++){
            if (afbTruth[iAfb] < -1) continue;
            if ( h_LRAfb[iTrueAfb]->GetBinContent(iAfb) < buffCov ) continue;
//printf("DEBUG\t\t: h_LRAfb(%d) = %f\n",iAfb,h_LRAfb[iTrueAfb]->GetBinContent(iAfb));
            if ( afbMeasuredLo[countAfbMeasuredHiLo] > -1 && (h_LRAfb[iTrueAfb]->GetBinContent(iAfb)-buffCov)*(h_LRAfb[iTrueAfb]->GetBinContent(iAfb+1)-buffCov) < 0  ){
                if ( fabs(h_LRAfb[iTrueAfb]->GetBinContent(iAfb)-h_LRAfb[iTrueAfb]->GetBinContent(iAfb-1))+fabs(h_LRAfb[iTrueAfb]->GetBinContent(iAfb)-h_LRAfb[iTrueAfb]->GetBinContent(iAfb+1)) > 0.2 ) continue;
                //printf("INFO\t\t: Upper limit found!\n");
                afbMeasuredHi[countAfbMeasuredHiLo] = afbTruth[iAfb];
            }
            if ( afbMeasuredLo[countAfbMeasuredHiLo] < -1 && h_LRAfb[iTrueAfb]->GetBinContent(iAfb-1) > buffCov ){
                //printf("INFO\t\t:   Lower limit found at boundary!\n");
                afbMeasuredLo[countAfbMeasuredHiLo] = afbTruth[iAfb-1];// boundary case
            }
            if ( afbMeasuredLo[countAfbMeasuredHiLo] < -1 && (h_LRAfb[iTrueAfb]->GetBinContent(iAfb)-buffCov)*(h_LRAfb[iTrueAfb]->GetBinContent(iAfb-1)-buffCov) < 0  ){
                //printf("INFO\t\t: Lower limit found!\n");
                afbMeasuredLo[countAfbMeasuredHiLo] = afbTruth[iAfb];
            }
        }
        if (afbMeasuredHi[countAfbMeasuredHiLo] < -1) afbMeasuredHi[countAfbMeasuredHiLo] = 0.75;// Usually not used, simply push to the limit.
        countAfbMeasuredHiLo+=1;
    }
    TH1F *h_LRFl[maxPoints];
    TH1F *h_buffProbAsLRFl = new TH1F("buffProbAsLRFl","",1000,0,1);
    int countFlMeasuredHiLo = 0;
    for(int iTrueFl = 0; iTrueFl*stepSizeFl < 1.; iTrueFl++){
        if (finFl[iTrueFl] == 0) continue;
        //printf("DEBUG\t\t: Looking for P(x,mu)/P(x,mu_best) for mu=thisFl=%+04f\n",flTruth[iTrueFl]);
        
        // Create LR ordering if finFl exists
        h_LRFl[iTrueFl] = new TH1F(TString::Format("h_LRFl_trueFl%+04.0f",flTruth[iTrueFl]*1000).Data(),"",1/stepSizeFl,0,1);
        double buffMaxLikelihoodFlRatio[maxPoints];
        for (int i = 0; i < maxPoints; i++) {
            buffMaxLikelihoodFlRatio[i]=-2;
        }
        for(int iFl = 0; iFl*stepSizeFl < 1.; iFl++){
            if (flTruth[iFl]< -1 ) continue;
            h_LRFl[iTrueFl]->SetBinContent(iFl,(*flPDFs[iTrueFl])(flTruth[iFl])*stepSizeFl/buffMaxLikelihoodFl[iFl]);
        }

        // Find interval
        double buffCov = 0;
        h_buffProbAsLRFl->Reset();
        for(int iFl = 0; iFl*stepSizeFl < 1.; iFl++){
            if (flTruth[iFl] < -1) continue;
            h_buffProbAsLRFl->Fill(h_LRFl[iTrueFl]->GetBinContent(iFl),(*flPDFs[iTrueFl])(flTruth[iFl])*stepSizeFl);
        }
        for(int iBin = h_buffProbAsLRFl->GetNbinsX(); iBin > 0; iBin--){
            buffCov += h_buffProbAsLRFl->GetBinContent(iBin);
            if (buffCov > targetCoverage){
                buffCov = h_buffProbAsLRFl->GetBinCenter(iBin);// serve as lower bound of LR
                break;
            }
        }
        //printf("DEBUG\t\t: lowest LR for iTrueFl=%d = %f\n", iTrueFl, buffCov);
        for(int iFl = 1; iFl*stepSizeFl < 1.; iFl++){
            if (flTruth[iFl] < -1) continue;
            if ( h_LRFl[iTrueFl]->GetBinContent(iFl) < buffCov ) continue;
            //printf("DEBUG\t\t: h_LRFl(%d) = %f\n",iFl,h_LRFl[iTrueFl]->GetBinContent(iFl));
            if ( flMeasuredLo[countFlMeasuredHiLo] > -1 && (h_LRFl[iTrueFl]->GetBinContent(iFl)-buffCov)*(h_LRFl[iTrueFl]->GetBinContent(iFl+1)-buffCov) < 0  ){
                if ( fabs(h_LRFl[iTrueFl]->GetBinContent(iFl)-h_LRFl[iTrueFl]->GetBinContent(iFl-1))+fabs(h_LRFl[iTrueFl]->GetBinContent(iFl)-h_LRFl[iTrueFl]->GetBinContent(iFl+1)) > 0.2 ) continue;
                //printf("INFO\t\t: Upper limit found!\n");
                flMeasuredHi[countFlMeasuredHiLo] = flTruth[iFl];
            }
            if ( flMeasuredLo[countFlMeasuredHiLo] < -1 && h_LRFl[iTrueFl]->GetBinContent(iFl-1) > buffCov ){
                //printf("INFO\t\t:   Lower limit found at boundary!\n");
                flMeasuredLo[countFlMeasuredHiLo] = flTruth[iFl-1];// boundary case
            }
            if ( flMeasuredLo[countFlMeasuredHiLo] < -1 && (h_LRFl[iTrueFl]->GetBinContent(iFl)-buffCov)*(h_LRFl[iTrueFl]->GetBinContent(iFl-1)-buffCov) < 0  ){
                //printf("INFO\t\t: Lower limit found!\n");
                flMeasuredLo[countFlMeasuredHiLo] = flTruth[iFl];
            }
        }
        if (flMeasuredHi[countFlMeasuredHiLo] < -1) flMeasuredHi[countFlMeasuredHiLo] = 1.;// Usually not used, simply push to the limit.
        countFlMeasuredHiLo+=1;
    }

    // Make plots, fit, and find confidence interval
    TFile *fout = new TFile(TString::Format("wspace_FCConfInterval_bin%d.root",iBin).Data(),"RECREATE");
    TCanvas *canvas = new TCanvas();
    TLatex *latex = new TLatex();
    latex->SetNDC();
    TLine *line = new TLine();
    TGraph *g_flIntervalHi = new TGraph(nFlPoints,flMeasuredHi,   flTrueValues);
    TGraph *g_flIntervalLo = new TGraph(nFlPoints,flMeasuredLo,   flTrueValues);
    TGraph *g_afbIntervalHi = new TGraph(nAfbPoints,afbMeasuredHi,afbTrueValues);
    TGraph *g_afbIntervalLo = new TGraph(nAfbPoints,afbMeasuredLo,afbTrueValues);
    TF1 *f1_flIntervalHi  = new TF1("f1_flIntervalHi",  "[0]+[1]*x+[2]*x**3+[3]*x**2",0,1);
    TF1 *f1_flIntervalLo  = new TF1("f1_flIntervalLo",  "[0]+[1]*x+[2]*x**3+[3]*exp([4]*x)",0,1); //
    TF1 *f1_afbIntervalHi = new TF1("f1_afbIntervalHi", "[0]+[1]*x+[2]*x**3+[3]*x**5",-0.75,0.75);
    TF1 *f1_afbIntervalLo = new TF1("f1_afbIntervalLo", "[0]+[1]*x+[2]*x**3+[3]*x**5",-0.75,0.75);
    double flFitInterval [4] = {0.05,0.95,0.05,0.95};
    double afbFitInterval[4] = {-0.7,0.7,-0.7,0.7};
    switch(iBin){
        case 10:
            flFitInterval[2] = 0.05;
            //flFitInterval[3] = 0.5;
            break;
        defalut:
            break;
    }
    //TFitResultPtr r_flIntervalHi = g_flIntervalHi->Fit("f1_flIntervalHi","S","",max(0.01,fl-0.1),min(0.99,fl+0.3));
    //TFitResultPtr r_flIntervalLi = g_flIntervalLo->Fit("f1_flIntervalLo","S","",max(0.01,fl-0.3),min(0.99,fl+0.1));
    //TFitResultPtr r_afbIntervalHi= g_afbIntervalHi->Fit("f1_afbIntervalHi","S","",max(-0.99,afb-0.1),min(0.99,afb+0.3));
    //TFitResultPtr r_afbIntervalLo= g_afbIntervalLo->Fit("f1_afbIntervalLo","S","",max(-0.99,afb-0.3),min(0.99,afb+0.1));
    TFitResultPtr r_flIntervalHi = g_flIntervalHi->Fit("f1_flIntervalHi","SN","",flFitInterval[0],flFitInterval[1]);
    TFitResultPtr r_flIntervalLi = g_flIntervalLo->Fit("f1_flIntervalLo","SN","",flFitInterval[2],flFitInterval[3]);
    TFitResultPtr r_afbIntervalHi= g_afbIntervalHi->Fit("f1_afbIntervalHi","SN","",afbFitInterval[0],afbFitInterval[1]);
    TFitResultPtr r_afbIntervalLo= g_afbIntervalLo->Fit("f1_afbIntervalLo","SN","",afbFitInterval[2],afbFitInterval[3]);
        // Calculate the error, you must think about the case in which f1 not defined!
    outFCErrFl[0] = r_flIntervalHi.Get()  != 0 ? max(0.,f1_flIntervalHi->Eval(fl))-fl : 0.-fl;
    outFCErrFl[1] = r_flIntervalLi.Get()  != 0 ? min(flTrueValues[nFlPoints-1],f1_flIntervalLo->Eval(fl))-fl : flTrueValues[nFlPoints-1]-fl;
    outFCErrAfb[0]= r_afbIntervalHi.Get() != 0 ? max(-1*afbTrueValues[nAfbPoints-1],f1_afbIntervalHi->Eval(afb))-afb : -1*afbTrueValues[nAfbPoints-1]-afb;
    outFCErrAfb[1]= r_afbIntervalLo.Get() != 0 ? min(afbTrueValues[nAfbPoints-1],f1_afbIntervalLo->Eval(afb))-afb : afbTrueValues[nAfbPoints-1]-afb;
    if ( outFCErrFl [1] < 0 ){
        outFCErrFl [1] = 0;
    }
    if ( outFCErrAfb[0] > 0 ){
        outFCErrAfb[0] = 0;
        outFCErrFl [1] = 0;
    }
    if ( outFCErrAfb[1] < 0 ){
        outFCErrAfb[1] = 0;
        outFCErrFl [1] = 0;
    }

        // Fl plots
    g_flIntervalHi->SetTitle("");
    g_flIntervalHi->GetXaxis()->SetTitle("Measured F_{L}");
    g_flIntervalHi->GetXaxis()->SetLimits(0.05,min(0.95,1.1*max(fl,flTrueValues[nFlPoints-1])));
    //g_flIntervalHi->GetXaxis()->SetLimits(0.05,1.1*max(fl,flTrueValues[nFlPoints-1]));
    g_flIntervalHi->GetYaxis()->SetTitle("True F_{L}");
    g_flIntervalHi->GetYaxis()->SetRangeUser(0.,max(fl,flTrueValues[nFlPoints-1]));
    g_flIntervalHi->SetMarkerStyle(7);
    g_flIntervalLo->SetMarkerStyle(7);
    g_flIntervalHi->SetMarkerSize(0.4);
    g_flIntervalLo->SetMarkerSize(0.4);
    g_flIntervalHi->Draw("AP");
    g_flIntervalLo->Draw("P SAME");
    f1_flIntervalHi ->Draw("SAME");
    f1_flIntervalLo ->Draw("SAME"); 
    line->SetLineStyle(2);
    line->SetLineColor(2);
    line->DrawLine(fl,0,fl,flTrueValues[nFlPoints-1]);
    line->SetLineStyle(2);
    line->SetLineColor(1);
    line->DrawLine(g_flIntervalHi->GetXaxis()->GetXmin(),flTrueValues[nFlPoints-1],g_flIntervalHi->GetXaxis()->GetXmax(),flTrueValues[nFlPoints-1]);
    latex->DrawLatexNDC(0.15,0.75,TString::Format("F_{L}=%.3f^{%+.3f}_{%+.3f}",fl,outFCErrFl[1],outFCErrFl[0]).Data());
    canvas->Update();
    canvas->Print(TString::Format("%s/FCConfInterval_fl_bin%d.pdf",plotpath.Data(),iBin));
    // Afb plots, remark true afb could be negative
    g_afbIntervalHi->SetTitle("");
    g_afbIntervalHi->GetXaxis()->SetTitle("Measured A_{FB}");
    g_afbIntervalHi->GetXaxis()->SetLimits(-0.70,0.70);
    //g_afbIntervalHi->GetXaxis()->SetLimits(-0.75,0.75);
    g_afbIntervalHi->GetYaxis()->SetTitle("True A_{FB}");
    g_afbIntervalHi->GetYaxis()->SetRangeUser(-1.1*max(afb,afbTrueValues[nAfbPoints-1]),1.1*max(afb,afbTrueValues[nAfbPoints-1]));
    g_afbIntervalHi->SetMarkerStyle(7);
    g_afbIntervalLo->SetMarkerStyle(7);
    g_afbIntervalHi->SetMarkerSize(0.4);
    g_afbIntervalLo->SetMarkerSize(0.4);
    g_afbIntervalHi->Draw("AP");
    g_afbIntervalLo->Draw("P SAME");
    f1_afbIntervalHi->Draw("SAME"); 
    f1_afbIntervalLo->Draw("SAME"); 
    line->SetLineStyle(2);
    line->SetLineColor(2);
    line->DrawLine(afb,-1*afbTrueValues[nAfbPoints-1]*1.1,afb,afbTrueValues[nAfbPoints-1]*1.1);
    line->SetLineStyle(2);
    line->SetLineColor(1);
    //line->DrawLine(g_afbIntervalHi->GetXaxis()->GetXmin(),afbTrueValues[nAfbPoints-1],g_afbIntervalHi->GetXaxis()->GetXmax(),afbTrueValues[nAfbPoints-1]);
    //line->DrawLine(g_afbIntervalHi->GetXaxis()->GetXmin(),-1*afbTrueValues[nAfbPoints-1],g_afbIntervalHi->GetXaxis()->GetXmax(),-1*afbTrueValues[nAfbPoints-1]);
    latex->DrawLatexNDC(0.15,0.75,TString::Format("A_{FB}=%.3f^{%+.3f}_{%+.3f}",afb,outFCErrAfb[1],outFCErrAfb[0]).Data());
    canvas->Update();
    canvas->Print(TString::Format("%s/FCConfInterval_afb_bin%d.pdf",plotpath.Data(),iBin));

    f1_afbIntervalHi->Write();
    f1_afbIntervalLo->Write();
    for(int iTrueAfb = 0; iTrueAfb*stepSizeAfb < 2.; iTrueAfb++){
        if (finAfb[iTrueAfb] == 0) continue;
        h_setSummaryAfb[iTrueAfb]->Write();
        afbPDFs[iTrueAfb]->Write();
        h_LRAfb[iTrueAfb]->Write();
        finAfb[iTrueAfb]->Close();
    }
    f1_flIntervalHi->Write();
    f1_flIntervalLo->Write();
    for(int iTrueFl = 0; iTrueFl*stepSizeFl < 1.; iTrueFl++){
        if (finFl[iTrueFl] == 0) continue;
        h_setSummaryFl[iTrueFl]->Write();
        flPDFs[iTrueFl]->Write();
        h_LRFl[iTrueFl]->Write();
        finFl[iTrueFl]->Close();
    }
    fout->Write();
    fout->Close();

    writeParam(iBin,"FCErrFl" ,outFCErrFl );
    writeParam(iBin,"FCErrAfb",outFCErrAfb);

    return;

}//}}}
void getFCIntervalRooStat(int iBin)
{//{{{
    // Take advantage of RooStat
    // Ref: https://root.cern.ch/root/html/tutorials/roostats/rs401d_FeldmanCousins.C.html
}//}}}
void getFCResultCoverage(int iBin)
{//{{{// Incorrect definition, use getFCResultCoverage2 for instead
    double afbErrLo = readParam(iBin, "FCErrAfb",0, 1);//lower error is negative
    double afbErrHi = readParam(iBin, "FCErrAfb",1,-1);
    double flErrLo  = readParam(iBin, "FCErrFl" ,0, 1);
    double flErrHi  = readParam(iBin, "FCErrFl" ,1,-1);
    double afb_ubd  = readParam(iBin, "afb"     ,0);
    double fl_ubd   = readParam(iBin, "fl"      ,0);
    double afb      = toBoundedAfb(afb_ubd,fl_ubd);
    double fl       = toBoundedFl(fl_ubd);

    TString toysDir = TString::Format("./limit/bestFit");// FIX ME!

    TFile *f_summary = new TFile(TString::Format("%s/setSummary_bin%d.root",toysDir.Data(),iBin).Data());
        
    TCanvas *canvas = new TCanvas();
    TLine *line = new TLine();
    line->SetLineStyle(2);
    line->SetLineColor(2);
    TLatex *latex = new TLatex();

    TH1F  *h_setSummaryAfb = (TH1F*)f_summary->Get("h_setSummaryAfb");
    double afbCoverage[1]={0};
    for(int jBin=1; jBin<=h_setSummaryAfb->GetNbinsX(); jBin++){
        if (h_setSummaryAfb->GetBinLowEdge(jBin) > afb+afbErrLo && h_setSummaryAfb->GetBinLowEdge(jBin+1) < afb+afbErrHi ){
            afbCoverage[0]+=h_setSummaryAfb->GetBinContent(jBin);
        }
    }
    afbCoverage[0]=afbCoverage[0]/h_setSummaryAfb->GetSumOfWeights();
    
    TH1F  *h_setSummaryFl  = (TH1F*)f_summary->Get("h_setSummaryFl");
    double flCoverage[1]={0};
    for(int jBin=1; jBin<=h_setSummaryFl->GetNbinsX(); jBin++){
        if (h_setSummaryFl->GetBinLowEdge(jBin) > fl+flErrLo && h_setSummaryFl->GetBinLowEdge(jBin+1) < fl+flErrHi ){
            flCoverage[0]+=h_setSummaryFl->GetBinContent(jBin);
        }
    }
    flCoverage[0] =flCoverage[0] /h_setSummaryFl ->GetSumOfWeights();
    
    h_setSummaryAfb->SetTitle("");
    h_setSummaryAfb->Draw();
    line->DrawLine(afb+afbErrLo, h_setSummaryAfb->GetMinimum(), afb+afbErrLo, h_setSummaryAfb->GetMaximum());
    line->DrawLine(afb+afbErrHi, h_setSummaryAfb->GetMinimum(), afb+afbErrHi, h_setSummaryAfb->GetMaximum());
    latex->DrawLatexNDC(0.15,0.8,TString::Format("Coverage=%.2f%%",afbCoverage[0]*100).Data());
    canvas->Update();
    canvas->Print(TString::Format("%s/coverageAfb_bin%d.pdf",plotpath.Data(),iBin).Data());

    h_setSummaryFl->SetTitle("");
    h_setSummaryFl->Draw();
    line->DrawLine(fl+flErrLo, h_setSummaryFl->GetMinimum(), fl+flErrLo, h_setSummaryFl->GetMaximum());
    line->DrawLine(fl+flErrHi, h_setSummaryFl->GetMinimum(), fl+flErrHi, h_setSummaryFl->GetMaximum());
    latex->DrawLatexNDC(0.15,0.8,TString::Format("Coverage=%.2f%%",flCoverage[0]*100));
    canvas->Update();
    canvas->Print(TString::Format("%s/coverageFl_bin%d.pdf",plotpath.Data(),iBin).Data());

    writeParam(iBin, "coverageAfb", afbCoverage, 1);
    writeParam(iBin, "coveragefl" , flCoverage , 1);

    return;
}//}}}
void getFCResultCoverage2(int iBin)
{//{{{
    double afbCoverLo(0.75), afbCoverHi(-0.75);
    double flCoverLo (1.  ), flCoverHi ( 0.  );
    double afb_ubd  = readParam(iBin, "afb"     ,0);
    double fl_ubd   = readParam(iBin, "fl"      ,0);
    double afb      = toBoundedAfb(afb_ubd,fl_ubd);
    double fl       = toBoundedFl(fl_ubd);

    TString toysDir = TString::Format("./limit/bestFit");// FIX ME!

    TFile *f_summary = new TFile(TString::Format("%s/setSummary_bin%d.root",toysDir.Data(),iBin).Data());
    TFile *f_fcresult= new TFile(TString::Format("%s/wspace_FCConfInterval_bin%d.root",iwspacepath.Data(),iBin).Data());
        
    TCanvas *canvas = new TCanvas();
    TLine *line = new TLine();
    line->SetLineStyle(2);
    line->SetLineColor(2);
    TLatex *latex = new TLatex();

    TH1F  *h_setSummaryAfb = (TH1F*)f_summary->Get("h_setSummaryAfb");
    TF1 *f1_afbIntervalLo = (TF1*)f_fcresult->Get("f1_afbIntervalLo");
    TF1 *f1_afbIntervalHi = (TF1*)f_fcresult->Get("f1_afbIntervalHi");
    double afbCoverage[1]={0};
    for(int jBin=1; jBin<=h_setSummaryAfb->GetNbinsX(); jBin++){
        double binCenter = h_setSummaryAfb->GetBinCenter(jBin);
        if (fabs(binCenter)>0.75) {continue;}
        // TODO: FIX HERE!!!
        if (min(0.75,f1_afbIntervalLo->Eval(binCenter)) > afb && max(-0.75,f1_afbIntervalHi->Eval(binCenter)) < afb ){
            if (binCenter > afbCoverHi){
                afbCoverHi = binCenter;
            }
            if (binCenter < afbCoverLo){
                afbCoverLo = binCenter;
            }
            afbCoverage[0]+=h_setSummaryAfb->GetBinContent(jBin);
        }
    }
    afbCoverage[0]=afbCoverage[0]/h_setSummaryAfb->GetSumOfWeights();
    
    TH1F  *h_setSummaryFl  = (TH1F*)f_summary->Get("h_setSummaryFl");
    TF1 *f1_flIntervalLo = (TF1*)f_fcresult->Get("f1_flIntervalLo");
    TF1 *f1_flIntervalHi = (TF1*)f_fcresult->Get("f1_flIntervalHi");
    double flCoverage[1]={0};
    for(int jBin=1; jBin<=h_setSummaryFl->GetNbinsX(); jBin++){
        double binCenter = h_setSummaryFl->GetBinCenter(jBin);
        if (min(1.,f1_flIntervalLo->Eval(binCenter)) > fl && max(0.,f1_flIntervalHi->Eval(binCenter)) < fl ){
            if (binCenter > flCoverHi){
                flCoverHi = binCenter;
            }
            if (binCenter < flCoverLo){
                flCoverLo = binCenter;
            }
            flCoverage[0]+=h_setSummaryFl->GetBinContent(jBin);
        }
    }
    flCoverage[0] =flCoverage[0] /h_setSummaryFl ->GetSumOfWeights();
    
    h_setSummaryAfb->SetTitle("");
    h_setSummaryAfb->Draw();
    //line->DrawLine(afbCoverLo, h_setSummaryAfb->GetMinimum(), afbCoverLo, h_setSummaryAfb->GetMaximum());
    //line->DrawLine(afbCoverHi, h_setSummaryAfb->GetMinimum(), afbCoverHi, h_setSummaryAfb->GetMaximum());
    latex->DrawLatexNDC(0.15,0.8,TString::Format("Coverage=%.2f%%",afbCoverage[0]*100).Data());
    canvas->Update();
    canvas->Print(TString::Format("%s/coverageAfb_bin%d.pdf",plotpath.Data(),iBin).Data());

    h_setSummaryFl->SetTitle("");
    h_setSummaryFl->Draw();
    //line->DrawLine(flCoverLo, h_setSummaryFl->GetMinimum(), flCoverLo, h_setSummaryFl->GetMaximum());
    //line->DrawLine(flCoverHi, h_setSummaryFl->GetMinimum(), flCoverHi, h_setSummaryFl->GetMaximum());
    latex->DrawLatexNDC(0.15,0.8,TString::Format("Coverage=%.2f%%",flCoverage[0]*100));
    canvas->Update();
    canvas->Print(TString::Format("%s/coverageFl_bin%d.pdf",plotpath.Data(),iBin).Data());

    writeParam(iBin, "coverageAfb", afbCoverage, 1);
    writeParam(iBin, "coveragefl" , flCoverage , 1);

    f_summary->Close();
    f_fcresult->Close();

    return;
}//}}}

    // Sys. Error determination tools
bool isEffMapNonNegative(RooGenericPdf *f_sigA)
{//{{{
    return true;
    //return false;
}//}}}
void rndEfficiencyMap(int iBin, RooGenericPdf *f_sigA, TMatrixDSym *errMtx)
{//{{{

    const char varNames[][32]={};
    // Prepare multivariate Gaussian
    RooArgSet *f_sigA_arg = f_sigA->getVariables();

    // randomly generate maps
    do{
    }while(!isEffMapNonNegative(f_sigA));

    printf("INFO: f_sigA transformed randomly.\n");
    return;
}//}}}
void trimErrMtx(TMatrixDSym *iErrMtx, TMatrixDSym *oErrMtx)
{//{{{
    return;
}//}}}
void rndEfficiencyMapTester()
{//{{{
    int iBin = 4;

    // Read data
    ch->SetBranchStatus("*",0);
    ch->SetBranchStatus("Bmass"         , 1);
    ch->SetBranchStatus("Mumumass"      , 1);
    ch->SetBranchStatus("Mumumasserr"   , 1);
    ch->SetBranchStatus("CosTheta*"     , 1);
    ch->SetBranchStatus("Q2"            , 1);
    ch->SetBranchStatus("Triggers"      , 1);
    RooRealVar Bmass("Bmass","M_{K^{*}#Mu#Mu}",5.1,5.6);
    RooRealVar CosThetaK("CosThetaK"     , "cos#theta_{K}"       , -1. , 1.   ) ;
    RooRealVar CosThetaL("CosThetaL"     , "cos#theta_{L}"       , -1. , 1.   ) ;
    RooRealVar Mumumass("Mumumass","M^{#mu#mu}",0.,10.);
    RooRealVar Mumumasserr("Mumumasserr","Error of M^{#mu#mu}",0.,10.);
    RooRealVar Q2("Q2","q^{2}",0.5,20.);
    RooRealVar Triggers("Triggers","",0,100);
    RooAddition Bmass_offset("Bmass_offset","Bmass_offset",RooArgSet(Bmass,RooConst(-5)));
    RooProduct Bmass_norm("Bmass_norm","Bmass_norm",RooArgSet(Bmass_offset,RooConst(1./0.5)));

    int mumuMassWindowBin = 1+2*isCDFcut;
    if (iBin==3 || iBin==5 || isCDFcut < 0) mumuMassWindowBin = 0; // no cut
    RooDataSet *data = new RooDataSet("data","data",ch,RooArgSet(CosThetaK, CosThetaL, Bmass, Q2, Mumumass, Mumumasserr, Triggers),TString::Format("(%s) && (%s) && (%s)",nTriggeredPath[2], q2range[iBin],mumuMassWindow[mumuMassWindowBin]),0);
    if (data->sumEntries() == 0){
        return;
    }

    // Create parameters and PDFs
    TFile *f_wspace_sigA = new TFile(TString::Format("%s/wspace_sigA_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace_sigA = (RooWorkspace*)f_wspace_sigA->Get("wspace");
    RooGenericPdf *f_sigA = 0;
    TMatrixDSym *errMtx = 0;
    if (wspace_sigA){
        f_sigA = (RooGenericPdf*)wspace_sigA->pdf("f_sigA");
        errMtx = (TMatrixDSym*)wspace_sigA->obj("errMtx");
    }else{
        printf("ERROR\t\t: Please have wsapce_sigA_bin?.root prepared.\n");
        return;
    }
    printf("INFO: f_sigA prepared.\n");

    f_sigA->Print();

    // Get variables
    RooArgSet *f_sigA_argset = f_sigA->getVariables();
    f_sigA_argset->Print();
    errMtx->Print();

    double *errMtxArray = errMtx->GetMatrixArray();
    int rowIdx = 1;
    int colIdx = 1;
    printf("%f %f\n",errMtxArray[rowIdx*21+colIdx],errMtxArray[rowIdx*21+colIdx+1]);

    // fix error matrix
    const double maxError = 1e10;
    const double minError = 1e-20;
    const double fixedSmallError = 1e-30;

    //rndEfficiencyMap(f_sigA);
    
    //RooMultiVarGaussian mvaGauss("mvaGauss","mvaGauss",RooArgList(),);

    return;
}//}}}

//_________________________________________________________________________________

// Validation tools
void getToyFromUnfilterGen(int iBin) // For validation, DEPRECATED.
{//{{{
    bool flatModel = false;

    // Generate random toy using unfiltered input and efficiency function
    double gQ2 = 0;
    double gCosThetaK = 0;
    double gCosThetaL = 0;
    double Q2 = 0;
    double CosThetaK = 0;
    double CosThetaL = 0;
    double Mumumass = 0;
    double Mumumasserr = 0;
    int    Triggers = 0;

    //TChain *treein=new TChain("tree");
    //treein->Add("./data/2012/sel_BuToKstarMuMu_NoGenFilter_8TeV_part*_mc.lite.root");//Contact po-hsun.chen@cern.ch to get these files.
    TChain *treein=ch;
    if (treein == NULL) {
        printf("Unfiltered MC sample is missing. Please contact pchen@cern.ch to get it.");
        return;
    }
    treein->SetBranchAddress("genQ2"        , &gQ2);
    treein->SetBranchAddress("genCosThetaK" , &gCosThetaK);
    treein->SetBranchAddress("genCosThetaL" , &gCosThetaL);
    treein->SetBranchAddress("Q2"           , &Q2);
    treein->SetBranchAddress("CosThetaK"    , &CosThetaK);
    treein->SetBranchAddress("CosThetaL"    , &CosThetaL);
    treein->SetBranchAddress("Mumumass"     , &Mumumass);
    treein->SetBranchAddress("Mumumasserr"  , &Mumumasserr);
    treein->SetBranchAddress("Triggers"     , &Triggers);

    // Get efficiency map, y=cosThetaK, x=cosThetaL
    std::string f2_model_format;
    if (flatModel){
        f2_model_format = TString::Format("0.01+0.*x+0.*y");
    }else{
        f2_model_format = TString::Format("%s*(1+([0]+[1]*CosThetaK+[2]*(3*CosThetaK**2-1)/2+[3]*(5*CosThetaK**3-3*CosThetaK)/2)+([4]+[5]*CosThetaK+[6]*(3*CosThetaK**2-1)/2+[7]*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**2+([8]+[9]*CosThetaK+[10]*(3*CosThetaK**2-1)/2+[11]*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**3+([12]+[13]*CosThetaK+[14]*(3*CosThetaK**2-1)/2+[15]*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**4+([16]+[17]*CosThetaK+[18]*(3*CosThetaK**2-1)/2+[19]*(5*CosThetaK**3-3*CosThetaK)/2)*CosThetaL**6)",readParam(iBin,"f_accXrecoEff_ord0","1").c_str());
    }
    //f2_model_format = regex_replace(f2_model_format,regex("CosThetaK"),"y");//use string::replace() for instead
    //f2_model_format = regex_replace(f2_model_format,regex("CosThetaL"),"x");
    while(f2_model_format.find("CosThetaK") !=std::string::npos ){
        f2_model_format.replace(f2_model_format.find("CosThetaK"), 9, "y");
    }
    while(f2_model_format.find("CosThetaL") !=std::string::npos ){
        f2_model_format.replace(f2_model_format.find("CosThetaL"), 9, "x");
    }
    printf("%s\n",f2_model_format.c_str());
    TF2 f2_model("f2_model",f2_model_format.c_str(),-1.,1.,-1.,1.);
    if (flatModel){
        //Set parameters here.
    }else{
        for (int i = 0; i < 20; i++) {
            f2_model.SetParameter(i,readParam(iBin,"accXrecoEff2",i));
            f2_model.SetParError(i,readParam(iBin,"accXrecoEff2Err",i));
        }
    }

    // 
    TFile fout(TString::Format("%s/rndToy_Bin%d.root",owspacepath.Data(),iBin), "RECREATE") ;
    TTree *treeout = treein->CloneTree(0);
    int _count = 0;//number of accepted events
    int _entry = 0;
    TRandom3 *rndGenerator = new TRandom3();
    do {
        treein->GetEntry(_entry); _entry++;
        if (gQ2 > q2rangeup[3] && gQ2 < q2rangedn[3]) continue;//jpsi
        if (gQ2 > q2rangeup[5] && gQ2 < q2rangedn[5]) continue;//psi2s
        if (gQ2 > q2rangedn[iBin] && gQ2 < q2rangeup[iBin]) {
            if (rndGenerator->Rndm() < f2_model.Eval(gCosThetaL,gCosThetaK)) {
                Q2 = gQ2;
                CosThetaL   = gCosThetaL;
                CosThetaK   = gCosThetaK;
                Mumumass    = sqrt(gQ2);
                Mumumasserr = 0.001;
                Triggers    = 1;
                treeout->Fill();
                _count++;
            }
        }
    } while ( _entry < treein->GetEntries() );
    fout.Write();
    fout.Close();
}//}}}

void genToySignal(int iBin, int nEvents = 100, double newAfb = -1., double newFl = -1., double newFs = -1., double newAs=-1.) // For validation.
{//{{{
    // Take Sm, sigA from iwspacepath
    int defNEvents  = 100;
    double defAfb   = genAfb[iBin];
    double defFl    = genFl[iBin];
    double defAfbErr= genAfberr[iBin];
    double defFlErr = genFlerr[iBin];
    double defAs    = -0.1;//-0.1
    double defFs    = 0.05;//0.05

    if (newAfb != -1.) defAfb = newAfb;
    if (newFl  != -1.) defFl  = newFl;
    if (newFs  != -1.) defFs  = newFs;
    if (newAs  != -1.) defAs  = newAs;

    // Get PDF and parameters
    TFile *f_wspace_Sm = new TFile(TString::Format("%s/wspace_Sm_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace_Sm = (RooWorkspace*)f_wspace_Sm->Get("wspace");
    RooAddPdf *f_sigM = 0;
    RooRealVar *Bmass = 0;
    RooRealVar *sigGauss_mean;
    RooRealVar *sigGauss1_sigma;
    RooRealVar *sigGauss2_sigma;
    RooRealVar *sigM_frac;
    if (wspace_Sm){
        f_sigM = (RooAddPdf*)wspace_Sm->pdf("f_sigM");
        Bmass = (RooRealVar*)wspace_Sm->var("Bmass");
        sigGauss_mean   = (RooRealVar*)wspace_Sm->var("sigGauss_mean");
        sigGauss1_sigma = (RooRealVar*)wspace_Sm->var("sigGauss1_sigma");
        sigGauss2_sigma = (RooRealVar*)wspace_Sm->var("sigGauss2_sigma");
        sigM_frac       = (RooRealVar*)wspace_Sm->var("sigM_frac");
    }else{
        printf("ERROR\t\t: Please have wsapce_Sm_bin?.root prepared.\n");
        return;
    }
    TFile *f_wspace_sigA = new TFile(TString::Format("%s/wspace_sigA_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace_sigA = (RooWorkspace*)f_wspace_sigA->Get("wspace");
    RooGenericPdf *f_sigA = 0;
    RooRealVar *CosThetaK = 0;
    RooRealVar *CosThetaL = 0;
    RooRealVar *afb;
    RooRealVar *fl;
    RooRealVar *fs;
    RooRealVar *as;
    if (wspace_sigA){
        f_sigA = (RooGenericPdf*)wspace_sigA->pdf("f_sigA");
        CosThetaK = (RooRealVar*)wspace_sigA->var("CosThetaK");
        CosThetaL = (RooRealVar*)wspace_sigA->var("CosThetaL");
        afb       = (RooRealVar*)wspace_sigA->var("afb");
        as        = (RooRealVar*)wspace_sigA->var("as");
        fl        = (RooRealVar*)wspace_sigA->var("fl");
        fs        = (RooRealVar*)wspace_sigA->var("fs");
        fs        ->setRange(1e-3,1.-1e-3);
        as        ->setRange(-10,10);//transformed as
        fl        ->setRange(-100,100);// unbounded fl
        afb       ->setRange(-100,100);// unbounded afb
    }else{
        printf("ERROR\t\t: Please have wsapce_sigA_bin?.root prepared.\n");
        return;
    }
    // Remark: The 'generate' function doesn't handle the parameter errors!!
    //         In order to do parameter estimation, values must be set manually.
    //fl ->setVal(toUnboundedFl(defFl));
    //fl ->setAsymError(toUnboundedFl(defFl-defFlErr)-toUnboundedFl(defFl),toUnboundedFl(defFl+defFlErr)-toUnboundedFl(defFl));
    //afb->setVal(toUnboundedAfb(defAfb,defFl));
    //afb->setAsymError(toUnboundedAfb(defAfb-defAfbErr,defFl)-toUnboundedAfb(defAfb,defFl),toUnboundedAfb(defAfb+defAfbErr,defFl)-toUnboundedAfb(defAfb,defFl));
    //as ->setVal(toTransformedAs(defFs,defFl,defAs));
    //fs ->setVal(defFs);

    // Random generator for afb/fl values
    RooGaussian gaus_afb("gaus_afb","",*afb,RooConst(defAfb),RooConst(defAfbErr));
    RooGaussian gaus_fl ("gaus_fl" ,"" ,*fl ,RooConst(defFl) ,RooConst(defFlErr));
    gaus_afb.Print();
    gaus_fl.Print();
    if (nEvents == 0) nEvents = defNEvents;

    RooProdPdf *f = new RooProdPdf("f", "f", RooArgSet(*f_sigM,*f_sigA));

    double rndBmass(0);
    double rndCosL(0);
    double rndCosK(0);
    double Q2((q2rangeup[iBin]+q2rangedn[iBin])/2);
    if (iBin == summaryBins[0]){
        Q2 = 2.;// Pick whatever value away from cuts.
    }
    double mumuMass(sqrt(Q2));
    double mumuMasserr(0.01);
    int    Triggers(1);
    double kstarMass(0.892);//signal region

    TFile fout(TString::Format("%s/signalToy_Bin%d.root",owspacepath.Data(),iBin), "RECREATE");
    TTree *tree = new TTree("tree","tree");
    tree->Branch("CosThetaK",&rndCosK,"CosThetaK/D");
    tree->Branch("CosThetaL",&rndCosL,"CosThetaL/D");
    tree->Branch("Bmass",&rndBmass,"Bmass/D");
    tree->Branch("Mumumass",&mumuMass,"Mumumass/D");
    tree->Branch("Mumumasserr",&mumuMasserr,"Mumumasserr/D");
    tree->Branch("Q2",&Q2,"Q2/D");
    tree->Branch("Triggers",&Triggers,"Triggers/I");
    tree->Branch("Kstarmass",&kstarMass,"Kstarmass/D");
    TH1D *h_Afb     =new TH1D("h_Afb","",100,-1,1);
    TH1D *h_Fl      =new TH1D("h_Fl","",100,0,1);
    TH2D *h2_AfbFl  =new TH2D("h2_AfbFl","",100,-1,1,100,0,1);
    TH1D *h_ubdAfb     =new TH1D("h_ubdAfb","",100,-100,100);
    TH1D *h_ubdFl      =new TH1D("h_ubdFl","",100,-100,100);
    TH2D *h2_ubdAfbFl  =new TH2D("h2_ubdAfbFl","",100,-100,100,100,-100,100);

    printf("INFO\t\t: Enter generating loop\n");
    bool needAfbFlError = false; // VERY SLOW if this is true
    if (needAfbFlError){
        const RooArgSet *dataInEntry = 0;
        RooDataSet* data_fl = gaus_fl .generate(*fl ,nEvents);
        RooDataSet* data_afb = gaus_afb.generate(*afb,nEvents);
        for (int iEvt = 0; iEvt < nEvents; iEvt++) {
            double temp_fl  = fabs( ((RooRealVar*)(data_fl->get(iEvt)->find("fl")))->getVal() );
            double temp_afb = ((RooRealVar*)(data_afb->get(iEvt)->find("afb")))->getVal();
            //printf("INFO\t\t: Event%d : temp_fl=%f\t, temp_afb=%f\t\n",iEvt,temp_fl,temp_afb);
            h_Fl        ->Fill(temp_fl);
            h_Afb       ->Fill(temp_afb);
            h2_AfbFl    ->Fill(temp_afb,temp_fl);
            fl ->setVal(toUnboundedFl(temp_fl));
            afb->setVal(toUnboundedAfb(temp_afb,temp_fl));
            as ->setVal(toTransformedAs(defFs,temp_fl,defAs));
            h_ubdFl     ->Fill(fl->getVal());
            h_ubdAfb    ->Fill(afb->getVal());
            h2_ubdAfbFl ->Fill(afb->getVal(),fl->getVal());
            dataInEntry = f->generate(RooArgSet(*Bmass,*CosThetaK,*CosThetaL),1)->get(0);
            rndBmass=((RooRealVar*)dataInEntry->find("Bmass"    ))->getVal();
            rndCosK =((RooRealVar*)dataInEntry->find("CosThetaK"))->getVal();
            rndCosL =((RooRealVar*)dataInEntry->find("CosThetaL"))->getVal();
            tree->Fill();
        }
    }else{
        fl ->setVal(toUnboundedFl(defFl));
        afb->setVal(toUnboundedAfb(defAfb,defFl));
        as ->setVal(toTransformedAs(defFs,defFl,defAs));
        fs ->setVal(defFs);
        printf("INFO\t\t: Generating signal toy with afb=%f , fl=%f , fs=%f , as=%f.\n",defAfb,defFl,defFs,defAs);
        afb->Print();
        fl ->Print();
        const RooDataSet *dataInEntry = 0;
        dataInEntry = f->generate(RooArgSet(*Bmass,*CosThetaK,*CosThetaL),nEvents);
        for (int iEvt = 0; iEvt < nEvents; iEvt++) {
            h_Fl        ->Fill(defFl);
            h_Afb       ->Fill(defAfb);
            h2_AfbFl    ->Fill(defAfb,defFl);
            h_ubdFl     ->Fill(fl->getVal());
            h_ubdAfb    ->Fill(afb->getVal());
            h2_ubdAfbFl ->Fill(afb->getVal(),fl->getVal());
            rndBmass=((RooRealVar*)(dataInEntry->get(iEvt))->find("Bmass"    ))->getVal();
            rndCosK =((RooRealVar*)(dataInEntry->get(iEvt))->find("CosThetaK"))->getVal();
            rndCosL =((RooRealVar*)(dataInEntry->get(iEvt))->find("CosThetaL"))->getVal();
            tree->Fill();
        }
    }
    fout.Write();
    fout.Close();

    f_wspace_Sm->Close();
    f_wspace_sigA->Close();
    delete f_wspace_Sm;
    delete f_wspace_sigA;
    return;
}//}}}

void genToyCombBkg(int iBin, int nEvents = 0, const char outfile[] = "combBkgToy") // For validation.
{//{{{
    // Take angular3D_bin from iwspacepath, prior from iCombBkgWspacepath
    TFile *f_wspace_comb_A = new TFile(TString::Format("%s/wspace_prior_bin%d.root",iCombBkgWspacepath.Data(),iBin));
    RooWorkspace *wspace_comb_A = (RooWorkspace*)f_wspace_comb_A->Get("wspace");
    RooRealVar *CosThetaK = 0;
    RooRealVar *CosThetaL = 0;
    RooGenericPdf *f_bkgCombA = 0;
    RooRealVar *bkgCombK_c1 = 0;
    RooRealVar *bkgCombK_c2 = 0;
    RooRealVar *bkgCombK_c3 = 0;
    RooRealVar *bkgCombK_c4 = 0;
    RooRealVar *bkgCombL_c1 = 0;
    RooRealVar *bkgCombL_c2 = 0;
    RooRealVar *bkgCombL_c3 = 0;
    RooRealVar *bkgCombL_c4 = 0;
    RooRealVar *bkgCombL_c5 = 0;
    if (wspace_comb_A){
        f_bkgCombA = (RooGenericPdf*)wspace_comb_A->pdf("f_bkgCombA");
        CosThetaK = (RooRealVar*)wspace_comb_A->var("CosThetaK");
        CosThetaL = (RooRealVar*)wspace_comb_A->var("CosThetaL");
        bkgCombK_c1 = (RooRealVar*)wspace_comb_A->var("bkgCombK_c1");
        bkgCombK_c2 = (RooRealVar*)wspace_comb_A->var("bkgCombK_c2");
        bkgCombK_c3 = (RooRealVar*)wspace_comb_A->var("bkgCombK_c3");
        bkgCombK_c4 = (RooRealVar*)wspace_comb_A->var("bkgCombK_c4");
        bkgCombL_c1 = (RooRealVar*)wspace_comb_A->var("bkgCombL_c1");
        bkgCombL_c2 = (RooRealVar*)wspace_comb_A->var("bkgCombL_c2");
        bkgCombL_c3 = (RooRealVar*)wspace_comb_A->var("bkgCombL_c3");
        bkgCombL_c4 = (RooRealVar*)wspace_comb_A->var("bkgCombL_c4");
        bkgCombL_c5 = (RooRealVar*)wspace_comb_A->var("bkgCombL_c5");
    }
    TFile *f_wspace_comb_M = new TFile(TString::Format("%s/wspace_angular3D_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace_comb_M = (RooWorkspace*)f_wspace_comb_M->Get("wspace");
    RooRealVar *Bmass = 0;
    RooGenericPdf *f_bkgCombM = 0;
    if (wspace_comb_M){
        Bmass = (RooRealVar*)wspace_comb_M->var("Bmass");
        f_bkgCombM = (RooGenericPdf*)wspace_comb_M->pdf("f_bkgCombM");
    }
    RooProdPdf f("f","f",RooArgSet(*f_bkgCombM,f_bkgCombA));
    
    double rndBmass(0);
    double rndCosL(0);
    double rndCosK(0);
    double Q2((q2rangeup[iBin]+q2rangedn[iBin])/2);
    if (iBin == summaryBins[0]){
        Q2 = 2.;// Pick whatever value away from cuts.
    }
    double mumuMass(sqrt(Q2));
    double mumuMasserr(0.01);
    double kstarMass(0.892);
    int    Triggers(1);
    
    TFile fout(TString::Format("%s/%s_Bin%d.root",owspacepath.Data(),outfile,iBin), "RECREATE");
    TTree *tree = new TTree("tree","tree");
    tree->Branch("CosThetaK",&rndCosK,"CosThetaK/D");
    tree->Branch("CosThetaL",&rndCosL,"CosThetaL/D");
    tree->Branch("Bmass",&rndBmass,"Bmass/D");
    tree->Branch("Mumumass",&mumuMass,"Mumumass/D");
    tree->Branch("Mumumasserr",&mumuMasserr,"Mumumasserr/D");
    tree->Branch("Q2",&Q2,"Q2/D");
    tree->Branch("Triggers",&Triggers,"Triggers/I");
    tree->Branch("Kstarmass",&kstarMass,"Kstarmass/D");

    RooDataSet *outdata = f.generate(RooArgSet(*Bmass,*CosThetaK,*CosThetaL),nEvents);
    const RooArgSet *dataInEntry = 0;
    for (int iEvt = 0; iEvt < nEvents; iEvt++) {
        dataInEntry = outdata->get(iEvt);
        rndBmass=((RooRealVar*)dataInEntry->find("Bmass"    ))->getVal();
        rndCosK =((RooRealVar*)dataInEntry->find("CosThetaK"))->getVal();
        rndCosL =((RooRealVar*)dataInEntry->find("CosThetaL"))->getVal();
        tree->Fill();
    }
    fout.Write();
    fout.Close();

    f_wspace_comb_M->Close();
    f_wspace_comb_A->Close();
    delete f_wspace_comb_A;
    delete f_wspace_comb_M;
}//}}}

void genToyCombBkgFromWspace(int iBin, int nEvents)
{//{{{
    TFile *f_wspace = new TFile(TString::Format("%s/wspace_angular3D_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
    if (!wspace) return;
    RooProdPdf *f_bkgComb = (RooProdPdf*)wspace->pdf("f_bkgComb");
    RooRealVar *Bmass = wspace->var("Bmass");
    RooRealVar *CosThetaK = wspace->var("CosThetaK");
    RooRealVar *CosThetaL = wspace->var("CosThetaL");

    TRandom3 *rndGenerator = new TRandom3();

    double rndBmass(0);
    double rndCosL(0);
    double rndCosK(0);
    double Q2((q2rangeup[iBin]+q2rangedn[iBin])/2);
    if (iBin == summaryBins[0]){
        Q2 = 2.;// Pick whatever value away from cuts.
    }
    double mumuMass(sqrt(Q2));
    double mumuMasserr(0.01);
    double kstarMass(0.892);
    int    Triggers(1);
    
    TFile fout(TString::Format("%s/%s_Bin%d.root",owspacepath.Data(),"combBkgToy",iBin), "RECREATE");
    TTree *tree = new TTree("tree","tree");
    tree->Branch("CosThetaK",&rndCosK,"CosThetaK/D");
    tree->Branch("CosThetaL",&rndCosL,"CosThetaL/D");
    tree->Branch("Bmass",&rndBmass,"Bmass/D");
    tree->Branch("Mumumass",&mumuMass,"Mumumass/D");
    tree->Branch("Mumumasserr",&mumuMasserr,"Mumumasserr/D");
    tree->Branch("Q2",&Q2,"Q2/D");
    tree->Branch("Triggers",&Triggers,"Triggers/I");
    tree->Branch("Kstarmass",&kstarMass,"Kstarmass/D");

    RooDataSet *outdata = f_bkgComb->generate(RooArgSet(*Bmass,*CosThetaK,*CosThetaL),nEvents);
    const RooArgSet *dataInEntry = 0;
    for (int iEvt = 0; iEvt < nEvents; iEvt++) {
        dataInEntry = outdata->get(iEvt);
        rndBmass=((RooRealVar*)dataInEntry->find("Bmass"    ))->getVal();
        rndCosK =((RooRealVar*)dataInEntry->find("CosThetaK"))->getVal();
        rndCosL =((RooRealVar*)dataInEntry->find("CosThetaL"))->getVal();
        tree->Fill();
    }
    fout.Write();
    fout.Close();

    f_wspace->Close();
}//}}}

void genToyPeakBkg(int iBin, int nEventsL = 0, int nEventsR = 0) // For vaidation. Use rndPickMCSamples/genToyPeakBkgFromWspace for instead.
{//{{{
    // Simplify the modeling by Gaussian
    double bkgPeakM_sigma(0.030);
    double bkgPeakM_mean(5.05);
    string f_bkgPeakM_formula = "exp(-0.5*((x-[0])/[1])**2)";
    TF1 *f_bkgPeakM = new TF1("f_bkgPeakM",f_bkgPeakM_formula.c_str(),5,5.56);
    f_bkgPeakM->SetParameter(0,bkgPeakM_mean);
    f_bkgPeakM->SetParameter(1,bkgPeakM_sigma);

    double bkgPeakL_c1(1);//x
    double bkgPeakL_c2(0.05*iBin);// 
    double bkgPeakL_c3(-0.3);// 
    double bkgPeakK_c1(1.);//y
    double bkgPeakK_c2(0.);// mean
    double bkgPeakK_c3(0.5);// 
    string f2_bkgPeakA_formula = "([0]+[2]*(x-[1])**2)*([3]+[5]*(y-[4])**2)";
    TF2 *f2_bkgPeakA = new TF2("f2_bkgPeakA",f2_bkgPeakA_formula.c_str(), -1, 1, -1, 1);
    f2_bkgPeakA->SetParameters(bkgPeakL_c1,bkgPeakL_c2,bkgPeakL_c3,bkgPeakK_c1,bkgPeakK_c2,bkgPeakK_c3);
    
    double rndBmass(0);
    double rndCosL(0);
    double rndCosK(0);
    double Q2((q2rangeup[iBin]+q2rangedn[iBin])/2);
    double mumuMass(sqrt(Q2));
    double mumuMasserr(0.01);
    int    Triggers(1);
    double kstarMass(0.892);
    
    TFile fout(TString::Format("%s/peakBkgToy_Bin%d.root",owspacepath.Data(),iBin), "RECREATE");
    TTree *tree = new TTree("tree","tree");
    tree->Branch("CosThetaK",&rndCosK,"CosThetaK/D");
    tree->Branch("CosThetaL",&rndCosL,"CosThetaL/D");
    tree->Branch("Bmass",&rndBmass,"Bmass/D");
    tree->Branch("Mumumass",&mumuMass,"Mumumass/D");
    tree->Branch("Mumumasserr",&mumuMasserr,"Mumumasserr/D");
    tree->Branch("Q2",&Q2,"Q2/D");
    tree->Branch("Triggers",&Triggers,"Triggers/I");
    tree->Branch("Kstarmass",&kstarMass,"Kstarmass/D");

    // LSB
    for (int iEvt = 0; iEvt < nEventsL; iEvt++) {
        rndBmass = f_bkgPeakM ->GetRandom();
        f2_bkgPeakA->GetRandom2(rndCosL,rndCosK);
        tree->Fill();
    }

    // RSB
    bkgPeakM_sigma = 0.035;
    bkgPeakM_mean = 5.50;
    f_bkgPeakM->SetParameter(0,bkgPeakM_mean);
    f_bkgPeakM->SetParameter(1,bkgPeakM_sigma);
    bkgPeakL_c3 = -0.07;
    bkgPeakK_c3 = -0.07;
    f2_bkgPeakA->SetParameters(bkgPeakL_c1,bkgPeakL_c2,bkgPeakL_c3,bkgPeakK_c1,bkgPeakK_c2,bkgPeakK_c3);
    for (int iEvt = 0; iEvt < nEventsR; iEvt++) {
        rndBmass = f_bkgPeakM ->GetRandom();
        f2_bkgPeakA->GetRandom2(rndCosL,rndCosK);
        tree->Fill();
    }

    fout.Write();
    fout.Close();
}//}}}

void genToyPeakBkgFromWspace(int iBin, double nEvents, const char decmode[]="both") // For error estimation/FC contour
{//{{{ 
    switch(iBin){
        case 2:
        case 4:
        case 6:
        case 9:
        case 10:
        case 11:
            break;
        default:
            return;
    }
    const int nPeakBkgType = 2;
    const char peakBkgType[nPeakBkgType][10] = {"jpsi", "psi2s"};
    int nEventsRound = round(nEvents);

    for(int iType=0; iType < nPeakBkgType; iType++){
        if (strcmp(decmode,peakBkgType[iType])*strcmp(decmode,"both") != 0) continue;
        TFile *f_wspace_peak_A = new TFile(TString::Format("%s/wspace_PkPl_%s_bin%d.root",iwspacepath.Data(),peakBkgType[iType],iBin));
        RooWorkspace *wspace_peak_A = (RooWorkspace*)f_wspace_peak_A->Get("wspace");
        RooRealVar *CosThetaK = 0;
        RooRealVar *CosThetaL = 0;
        RooGenericPdf *f_bkgPeakA = 0;
        if (wspace_peak_A){
            f_bkgPeakA = (RooGenericPdf*)wspace_peak_A->pdf(TString::Format("f_bkg%sPeakA",peakBkgType[iType]).Data());
            CosThetaK = (RooRealVar*)wspace_peak_A->var("CosThetaK");
            CosThetaL = (RooRealVar*)wspace_peak_A->var("CosThetaL");
        }
        TFile *f_wspace_peak_M = new TFile(TString::Format("%s/wspace_YpPm_%s_bin%d.root",iwspacepath.Data(),peakBkgType[iType],iBin));
        RooWorkspace *wspace_peak_M = (RooWorkspace*)f_wspace_peak_M->Get("wspace");
        RooRealVar *Bmass = 0;
        RooAddPdf *f_bkgPeakM = 0;
        if (wspace_peak_M){
            Bmass = (RooRealVar*)wspace_peak_M->var("Bmass");
            f_bkgPeakM = (RooAddPdf*)wspace_peak_M->pdf(TString::Format("f_bkg%sPeakM12",peakBkgType[iType]).Data());
        }
        RooProdPdf f("f","f",RooArgSet(*f_bkgPeakM,f_bkgPeakA));

        double rndBmass(0);
        double rndCosL(0);
        double rndCosK(0);
        double Q2((q2rangeup[iBin]+q2rangedn[iBin])/2);
        double mumuMass(sqrt(Q2));
        double mumuMasserr(0.01);
        int    Triggers(1);
        double kstarMass(0.892);

        TFile fout(TString::Format("%s/%sPeakToy_Bin%d.root",owspacepath.Data(),peakBkgType[iType],iBin), "RECREATE");
        TTree *tree = new TTree("tree","tree");
        tree->Branch("CosThetaK",&rndCosK,"CosThetaK/D");
        tree->Branch("CosThetaL",&rndCosL,"CosThetaL/D");
        tree->Branch("Bmass",&rndBmass,"Bmass/D");
        tree->Branch("Mumumass",&mumuMass,"Mumumass/D");
        tree->Branch("Mumumasserr",&mumuMasserr,"Mumumasserr/D");
        tree->Branch("Q2",&Q2,"Q2/D");
        tree->Branch("Triggers",&Triggers,"Triggers/I");
        tree->Branch("Kstarmass",&kstarMass,"Kstarmass/D");

        RooDataSet *outdata = f.generate(RooArgSet(*Bmass,*CosThetaK,*CosThetaL),nEventsRound);
        const RooArgSet *dataInEntry = 0;
        for (int iEvt = 0; iEvt < nEventsRound; iEvt++) {
            dataInEntry = outdata->get(iEvt);
            rndBmass=((RooRealVar*)dataInEntry->find("Bmass"    ))->getVal();
            rndCosK =((RooRealVar*)dataInEntry->find("CosThetaK"))->getVal();
            rndCosL =((RooRealVar*)dataInEntry->find("CosThetaL"))->getVal();
            tree->Fill();
        }
        fout.Write();
        fout.Close();

        f_wspace_peak_M->Close();
        f_wspace_peak_A->Close();
        delete f_wspace_peak_A;
        delete f_wspace_peak_M;
    }
}//}}}

void genToyPeakBkgFromWspace(int iBin, int scaleFactor = 1)
{//{{{
    TFile *f_wspace = new TFile(TString::Format("%s/wspace_angular3D_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
    if (!wspace) return;

    // Yields of components governed by Poisson and Gaussian
    TRandom3 *rndGenerator = new TRandom3();
    RooGaussian *gaus_nbkgjpsiPeak = 0;
    RooGaussian *gaus_nbkgpsi2sPeak = 0;
    RooRealVar nbkgjpsiPeak("nbkgjpsiPeak","",0,1e8);
    RooRealVar nbkgpsi2sPeak("nbkgpsi2sPeak","",0,1e8);
    RooDataSet *nbkgjpsiPeakDataset = 0;
    RooDataSet *nbkgpsi2sPeakDataset = 0;
    int nbkgjpsiPeakInToy = 0;
    int nbkgpsi2sPeakInToy = 0;
    switch(iBin){
        case 10:
            gaus_nbkgjpsiPeak = (RooGaussian*)wspace->pdf("gaus_nbkgjpsiPeak");
            nbkgjpsiPeakDataset = gaus_nbkgjpsiPeak->generate(nbkgjpsiPeak,1);
            nbkgjpsiPeakInToy=((RooRealVar*)(nbkgjpsiPeakDataset->get(0))->find("nbkgjpsiPeak"))->getVal()*scaleFactor;
            genToyPeakBkgFromWspace(iBin,nbkgjpsiPeakInToy,"jpsi");
            break;
        case 4:
            gaus_nbkgpsi2sPeak = (RooGaussian*)wspace->pdf("gaus_nbkgpsi2sPeak");
            nbkgpsi2sPeakDataset = gaus_nbkgpsi2sPeak->generate(nbkgpsi2sPeak,1);
            nbkgpsi2sPeakInToy=((RooRealVar*)(nbkgpsi2sPeakDataset->get(0))->find("nbkgpsi2sPeak"))->getVal()*scaleFactor;
            genToyPeakBkgFromWspace(iBin,nbkgpsi2sPeakInToy,"psi2s");
            break;
        case 11:
            gaus_nbkgjpsiPeak = (RooGaussian*)wspace->pdf("gaus_nbkgjpsiPeak");
            nbkgjpsiPeakDataset = gaus_nbkgjpsiPeak->generate(nbkgjpsiPeak,1);
            nbkgjpsiPeakInToy=((RooRealVar*)(nbkgjpsiPeakDataset->get(0))->find("nbkgjpsiPeak"))->getVal()*scaleFactor;
            genToyPeakBkgFromWspace(iBin,nbkgjpsiPeakInToy,"jpsi");
            gaus_nbkgpsi2sPeak = (RooGaussian*)wspace->pdf("gaus_nbkgpsi2sPeak");
            nbkgpsi2sPeakDataset = gaus_nbkgpsi2sPeak->generate(nbkgpsi2sPeak,1);
            nbkgpsi2sPeakInToy=((RooRealVar*)(nbkgpsi2sPeakDataset->get(0))->find("nbkgpsi2sPeak"))->getVal()*scaleFactor;
            genToyPeakBkgFromWspace(iBin,nbkgpsi2sPeakInToy,"psi2s");
            break;
        default:
            break;
    }

    return;
}//}}}

void splitMCSamples(double oLumis=datasetLumi[0]) // Split MC samples into small parts
{//{{{
    static char decmode[20];
    while(strcmp(decmode,"jpsi")*strcmp(decmode, "psi2s")*strcmp(decmode, "signal")*strcmp(decmode, "manual") != 0){
        printf("Please insert type [ signal / jpsi / psi2s / manual ]:");
        scanf("%19s",decmode);
    }
    
    double iLumis = 0.;
    if (strcmp(decmode, "signal") == 0){
        iLumis=datasetLumi[1];
    }else if (strcmp(decmode, "jpsi") == 0){
        iLumis = datasetLumi[2];
    }else if (strcmp(decmode, "psi2s") == 0){
        iLumis = datasetLumi[3];
    }else{
        printf("Input luminosity: ");
        scanf("%lf",&iLumis);
    }
    printf("INFO\t\t: Splitting %f/fb into %.0f parts\n",iLumis, floor(iLumis/oLumis));
    
    if (oLumis > iLumis){
        printf("ERROR\t\t: Luminosity of output should be larger than the reservoir. EXIT 1.");
        return;
    }

    int mumuMassWindowBin = 1+2*abs(isCDFcut);
    //if (isCDFcut < 0) mumuMassWindowBin = 0;
    for(int iPart=0; iPart < floor(iLumis/oLumis); iPart++){
        TFile *fout = new TFile(TString::Format("%s/splitMC_%s_part%04d.root",owspacepath.Data(),decmode,iPart+1), "RECREATE");
        TTree *tout = ch->CopyTree("","",(int)floor(oLumis/iLumis*ch->GetEntries()),iPart*(int)floor(oLumis/iLumis*ch->GetEntries()))->CopyTree(mumuMassWindow[mumuMassWindowBin]);
        fout->Write();
        fout->Close();
    }
    return;
}//}}} 

void splitToySamples(const char iFile[], const int nEvtsPerSet, const char oFileLabel[]="splitToy" )
{//{{{
    TFile *fin = new TFile(iFile);
    TTree *tree = (TTree*)fin->Get("tree");
    for(int iPart=0; iPart < ceil(tree->GetEntries())/nEvtsPerSet; iPart++){
        TFile *fout = new TFile(TString::Format("%s/%s_set%04d.root",owspacepath.Data(),oFileLabel,iPart+1), "RECREATE");
        TTree *tout = tout->CopyTree("","",nEvtsPerSet,iPart*nEvtsPerSet);
        fout->Write();
        fout->Close();
    }
    fin->Close();
    return;
}//}}}

void splitToySamples(const char iFile[], const int nSets, const int *nEvtsPerSet, const char oFileLabel[]="splitToy" )
{//{{{
    TFile *fin = new TFile(iFile);
    TTree *tree = (TTree*)fin->Get("tree");
    long int evtCounter = 0;
    for(int iPart=0; iPart < nSets; iPart++){
        TFile *fout = new TFile(TString::Format("%s/%s_set%04d.root",owspacepath.Data(),oFileLabel,iPart+1), "RECREATE");
        TTree *tout = tree->CopyTree("","",nEvtsPerSet[iPart],evtCounter);
        evtCounter+=nEvtsPerSet[iPart];
        fout->Write();
        fout->Close();
    }
    fin->Close();
    return;
}//}}}

void rndPickMCSamples(int nSets = 1, double oLumis=datasetLumi[0]) // Pick distinct events from MC
{//{{{
    static char decmode[20];
    while(strcmp(decmode,"jpsi")*strcmp(decmode, "psi2s")*strcmp(decmode, "signal")*strcmp(decmode, "manual") != 0){
        printf("Please insert type [ signal / jpsi / psi2s / manual ]:");
        scanf("%19s",decmode);
    }
    
    double iLumis = 0.;
    if (strcmp(decmode, "signal") == 0){
        iLumis=datasetLumi[1];
    }else if (strcmp(decmode, "jpsi") == 0){
        iLumis = datasetLumi[2];
    }else if (strcmp(decmode, "psi2s") == 0){
        iLumis = datasetLumi[3];
    }else{
        printf("Input luminosity: ");
        scanf("%lf",&iLumis);
    }
    printf("INFO\t\t: Randomly pick %d sets of events of %f/fb \n", nSets, oLumis);

    if (oLumis > iLumis){
        printf("ERROR\t\t: Luminosity of output should be larger than the reservoir. EXIT 1.");
        return;
    }

    // Initialize random generator
    TRandom3 *random = new TRandom3(time(0));

    int nEntries = ch->GetEntries();
    int nEvtsPerSet = floor(oLumis/iLumis*nEntries);
    int mumuMassWindowBin = 1+2*abs(isCDFcut);
    for(int iSet=0; iSet< nSets; iSet++){
        TFile *fout = new TFile(TString::Format("%s/rndMC_%s_%devts_set%d.root",owspacepath.Data(),decmode,nEvtsPerSet,iSet+1), "RECREATE");
        TTree *tout = ch->CloneTree(0);
        int candEvtId=0;
        std::vector<bool> bits(nEntries,false); // all false;
        for(int iEvt = nEvtsPerSet; iEvt > 0; iEvt--){
            do {
                printf("DEBUG\t\t: iEvt=%d, candEvtId=%d, Picking distict.\n",iEvt,candEvtId);
                candEvtId = random->Rndm()*nEntries;
            }while(bits[candEvtId]);
            bits[candEvtId]=true;
            // Don't access files here, random access really takes time!
        }
        for(int iBit = 0; iBit < nEntries; iBit++){
            if (bits[iBit] == false) continue;
            ch->GetEntry(candEvtId);
            tout->Fill();
        }
        tout = tout->CopyTree("mumuMassWindow[mumuMassWindowBin]");
        tout->FlushBaskets();
        tout->AutoSave();
        fout->Write();
        fout->Close();
    }
}//}}}

void createToy(int iBin, int nToy, double sigACoef[4], double yields[6])
{//{{{
    // sigA  [4]={afb_bdd ,fl_bdd  ,fs     ,as_ori}
    // yields[4]={nSig    ,nBkg    ,nJpsi  ,nJpsiErr ,nPsi2S, nPsi2SErr}

    //printf("INFO\t\t: createToy for Bin%d.\n",iBin);
    struct stat fiBuff;
    FILE *fi = 0;
    TString otoyspath = owspacepath;
    

    // Check output directory
    if (stat(otoyspath,&fiBuff) != 0){
        return;
    }
    
    // Generate nEvts in each toys
    TRandom3 *rndGenerator = new TRandom3();
    int nsigInToy           [nToy];
    int nbkgCombInToy       [nToy];
    int nbkgjpsiPeakInToy   [nToy];
    int nbkgpsi2sPeakInToy  [nToy];
    int nsigInToys          = 0;
    int nbkgCombInToys      = 0;
    int nbkgjpsiPeakInToys  = 0;
    int nbkgpsi2sPeakInToys = 0;
    RooRealVar nbkgjpsiPeak ("nbkgjpsiPeak" ,"",0,1e8);
    RooRealVar nbkgpsi2sPeak("nbkgpsi2sPeak","",0,1e8);
    RooGaussian *gaus_nbkgjpsiPeak  = 0;
    RooGaussian *gaus_nbkgpsi2sPeak = 0;
    RooDataSet  *nbkgjpsiPeakDataset  = 0;
    RooDataSet  *nbkgpsi2sPeakDataset = 0;
    switch(iBin){
        //case 10:
        //    gaus_nbkgjpsiPeak  = new RooGaussian("gaus_nbkgjpsiPeak","",nbkgjpsiPeak,RooConst(yields[2]),RooConst(yields[3]));
        //    nbkgjpsiPeakDataset  = gaus_nbkgjpsiPeak->generate(nbkgjpsiPeak,nToy);
        //    break;
        //case 11:
        //    gaus_nbkgjpsiPeak  = new RooGaussian("gaus_nbkgjpsiPeak","",nbkgjpsiPeak,RooConst(yields[2]),RooConst(yields[3]));
        //    gaus_nbkgpsi2sPeak = new RooGaussian("gaus_nbkgpsi2sPeak","",nbkgpsi2sPeak,RooConst(yields[4]),RooConst(yields[5]));
        //    nbkgjpsiPeakDataset  = gaus_nbkgjpsiPeak->generate(nbkgjpsiPeak,nToy);
        //    nbkgpsi2sPeakDataset = gaus_nbkgpsi2sPeak->generate(nbkgpsi2sPeak,nToy);
        //    break;
        //case 4:
        //    gaus_nbkgpsi2sPeak = new RooGaussian("gaus_nbkgpsi2sPeak","",nbkgpsi2sPeak,RooConst(yields[4]),RooConst(yields[5]));
        //    nbkgpsi2sPeakDataset = gaus_nbkgpsi2sPeak->generate(nbkgpsi2sPeak,nToy);
        //    break;
        default:
            break;
    }
    for(int iToy=0; iToy<nToy; iToy++){
        nsigInToy       [iToy]=rndGenerator->Poisson(yields[0])*scaleFactor;
        nbkgCombInToy   [iToy]=rndGenerator->Poisson(yields[1])*scaleFactor;
        nsigInToys           +=nsigInToy    [iToy];
        nbkgCombInToys       +=nbkgCombInToy[iToy];
        switch(iBin){
            //case 10:
            //    nbkgjpsiPeakInToy[iToy]=((RooRealVar*)(nbkgjpsiPeakDataset->get(iToy))->find("nbkgjpsiPeak"))->getVal()*scaleFactor;
            //    nbkgjpsiPeakInToys+=nbkgjpsiPeakInToy[iToy];
            //    break;
            //case 4:
            //    nbkgpsi2sPeakInToy[iToy]=((RooRealVar*)(nbkgpsi2sPeakDataset->get(iToy))->find("nbkgpsi2sPeak"))->getVal()*scaleFactor;
            //    nbkgpsi2sPeakInToys+=nbkgpsi2sPeakInToy[iToy];
            //    break;
            //case 11:
            //    nbkgjpsiPeakInToy[iToy]=((RooRealVar*)(nbkgjpsiPeakDataset->get(iToy))->find("nbkgjpsiPeak"))->getVal()*scaleFactor;
            //    nbkgjpsiPeakInToys+=nbkgjpsiPeakInToy[iToy];
            //    nbkgpsi2sPeakInToy[iToy]=((RooRealVar*)(nbkgpsi2sPeakDataset->get(iToy))->find("nbkgpsi2sPeak"))->getVal()*scaleFactor;
            //    nbkgpsi2sPeakInToys+=nbkgpsi2sPeakInToy[iToy];
            //    break;
            default:
                nbkgjpsiPeakInToy [iToy] = 0;
                nbkgpsi2sPeakInToy[iToy] = 0;
                break;
        }
    }

    // Generating
    printf("INFO\t\t: Generating signal toys.\n");
    genToySignal(iBin, nsigInToys, sigACoef[0], sigACoef[1], sigACoef[2], sigACoef[3]);
    splitToySamples(TString::Format("%s/signalToy_Bin%d.root" ,otoyspath.Data(),iBin).Data(),nToy,nsigInToy,TString::Format("signalToy_Bin%d",iBin).Data());
    printf("INFO\t\t: Generating combinetorial background toys.\n");
    genToyCombBkgFromWspace(iBin,nbkgCombInToys);
    splitToySamples(TString::Format("%s/combBkgToy_Bin%d.root",otoyspath.Data(),iBin).Data(),nToy,nbkgCombInToy,TString::Format("combBkgToy_Bin%d",iBin).Data());
    printf("INFO\t\t: Generating peaking background toys.\n");
    switch(iBin){
        //case 10:
        //    genToyPeakBkgFromWspace(iBin, nbkgjpsiPeakInToys,"jpsi");
        //    splitToySamples(TString::Format("%s/jpsiPeakToy_Bin%d.root",otoyspath.Data(),iBin).Data(),nToy,nbkgjpsiPeakInToy,TString::Format("jpsiPeakToy_Bin%d",iBin).Data());
        //    break;
        //case 4:
        //    genToyPeakBkgFromWspace(iBin, nbkgpsi2sPeakInToys,"psi2s");
        //    splitToySamples(TString::Format("%s/psi2sPeakToy_Bin%d.root",otoyspath.Data(),iBin).Data(),nToy,nbkgpsi2sPeakInToy,TString::Format("psi2sPeakToy_Bin%d",iBin).Data());
        //    break;
        //case 11:
        //    genToyPeakBkgFromWspace(iBin, nbkgjpsiPeakInToys,"jpsi");
        //    splitToySamples(TString::Format("%s/jpsiPeakToy_Bin%d.root",otoyspath.Data(),iBin).Data(),nToy,nbkgjpsiPeakInToy,TString::Format("jpsiPeakToy_Bin%d",iBin).Data());
        //    genToyPeakBkgFromWspace(iBin, nbkgpsi2sPeakInToys,"psi2s");
        //    splitToySamples(TString::Format("%s/psi2sPeakToy_Bin%d.root",otoyspath.Data(),iBin).Data(),nToy,nbkgpsi2sPeakInToy,TString::Format("psi2sPeakToy_Bin%d",iBin).Data());
        //    break;
        default:
            break;
    }
    for(int iToy = 0; iToy<nToy; iToy++){
        if (stat(TString::Format("%s/set%04d",otoyspath.Data(),iToy+1),&fiBuff) != 0){
            mkdir(TString::Format("%s/set%04d",otoyspath.Data(),iToy+1),0755);
        }
        rename(TString::Format("%s/signalToy_Bin%d_set%04d.root"         ,otoyspath.Data(),       iBin,iToy+1).Data(),
               TString::Format("%s/set%04d/signalToy_Bin%d_set%04d.root" ,otoyspath.Data(),iToy+1,iBin,iToy+1).Data());
        rename(TString::Format("%s/combBkgToy_Bin%d_set%04d.root"        ,otoyspath.Data(),       iBin,iToy+1).Data(),
               TString::Format("%s/set%04d/combBkgToy_Bin%d_set%04d.root",otoyspath.Data(),iToy+1,iBin,iToy+1).Data());
        switch(iBin){
            //case 10:
            //    rename(TString::Format("%s/jpsiPeakToy_Bin%d_set%04d.root"          ,otoyspath.Data(),       iBin,iToy+1).Data(),
            //           TString::Format("%s/set%04d/jpsiPeakToy_Bin%d_set%04d.root"  ,otoyspath.Data(),iToy+1,iBin,iToy+1).Data());
            //    break;
            //case 4:
            //    rename(TString::Format("%s/psi2sPeakToy_Bin%d_set%04d.root"         ,otoyspath.Data(),       iBin,iToy+1).Data(),
            //           TString::Format("%s/set%04d/psi2sPeakToy_Bin%d_set%04d.root" ,otoyspath.Data(),iToy+1,iBin,iToy+1).Data());
            //    break;
            //case 11:
            //    rename(TString::Format("%s/jpsiPeakToy_Bin%d_set%04d.root"          ,otoyspath.Data(),       iBin,iToy+1).Data(),
            //           TString::Format("%s/set%04d/jpsiPeakToy_Bin%d_set%04d.root"  ,otoyspath.Data(),iToy+1,iBin,iToy+1).Data());
            //    rename(TString::Format("%s/psi2sPeakToy_Bin%d_set%04d.root"         ,otoyspath.Data(),       iBin,iToy+1).Data(),
            //           TString::Format("%s/set%04d/psi2sPeakToy_Bin%d_set%04d.root" ,otoyspath.Data(),iToy+1,iBin,iToy+1).Data());
            //    break;
            default:
                break;
        }
    }
}//}}}

void createFCToys(int iBin, int nToy=500)// DEPRECATED, Create toys for contour scanning
{//{{{
    // Remark: keep extension possibility for 2-D contour.

    // Buffers for checking directory/FILE
    bool contourMode=false;// In case contour is needed.
    struct stat fiBuff;
    FILE *fi = 0;

    // Setting
    TString otoyspath = TString::Format("./limit/bin%d",iBin);
    const double stepSizeAfb = 0.01;
    const double stepSizeFl  = 0.01;

    // Create output directory
    if (stat(TString::Format("./limit"),&fiBuff) != 0){
        mkdir(TString::Format("./limit"),0755);
    }
    if (stat(TString::Format("%s",otoyspath.Data()),&fiBuff) != 0){
        mkdir(TString::Format("%s",otoyspath.Data()),0755);
    }

    // Get parameters for the q2 Bin
    TFile *f_wspace = new TFile(TString::Format("%s/wspace_angular3D_bin%d.root",iwspacepath.Data(),iBin));
    RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");
    if (!wspace) return;
    double  fl      = toBoundedFl(wspace->var("fl")->getVal());
    double  afb     = toBoundedAfb(wspace->var("afb")->getVal(),toUnboundedFl(fl));
    double  nsig    = wspace->var("nsig")->getVal();
    double  nbkgComb= wspace->var("nbkgComb")->getVal();

    // Yields of components governed by Poisson and Gaussian
    RooGaussian *gaus_nbkgjpsiPeak = 0;
    RooGaussian *gaus_nbkgpsi2sPeak = 0;
    switch(iBin){
        case 10:
            gaus_nbkgjpsiPeak = (RooGaussian*)wspace->pdf("gaus_nbkgjpsiPeak");
            break;
        case 4:
            gaus_nbkgpsi2sPeak = (RooGaussian*)wspace->pdf("gaus_nbkgpsi2sPeak");
            break;
        case 11:
            gaus_nbkgjpsiPeak = (RooGaussian*)wspace->pdf("gaus_nbkgjpsiPeak");
            gaus_nbkgpsi2sPeak = (RooGaussian*)wspace->pdf("gaus_nbkgpsi2sPeak");
            break;
        default:
            break;
    }
    TRandom3 *rndGenerator = new TRandom3();
    int nsigInToy[nToy];
    int nbkgCombInToy[nToy];
    int nbkgjpsiPeakInToy[nToy];
    int nbkgpsi2sPeakInToy[nToy];
    int nsigInToys = 0;
    int nbkgCombInToys = 0;
    int nbkgjpsiPeakInToys = 0;
    int nbkgpsi2sPeakInToys = 0;
    RooRealVar nbkgjpsiPeak("nbkgjpsiPeak","",0,1e8);
    RooRealVar nbkgpsi2sPeak("nbkgpsi2sPeak","",0,1e8);
    RooDataSet *nbkgjpsiPeakDataset = 0;
    RooDataSet *nbkgpsi2sPeakDataset = 0;
    switch(iBin){
        case 10:
            nbkgjpsiPeakDataset = gaus_nbkgjpsiPeak->generate(nbkgjpsiPeak,nToy);
            break;
        case 4:
            nbkgpsi2sPeakDataset = gaus_nbkgpsi2sPeak->generate(nbkgpsi2sPeak,nToy);
            break;
        case 11:
            nbkgjpsiPeakDataset = gaus_nbkgjpsiPeak->generate(nbkgjpsiPeak,nToy);
            nbkgpsi2sPeakDataset = gaus_nbkgpsi2sPeak->generate(nbkgpsi2sPeak,nToy);
            break;
        default:
            break;
    }
    for(int iToy=0; iToy<nToy; iToy++){
        nsigInToy[iToy]=rndGenerator->Poisson(nsig);
        nbkgCombInToy[iToy]=rndGenerator->Poisson(nbkgComb);
        switch(iBin){
            case 10:
                nbkgjpsiPeakInToy[iToy]=((RooRealVar*)(nbkgjpsiPeakDataset->get(iToy))->find("nbkgjpsiPeak"))->getVal();
                nbkgjpsiPeakInToys+=nbkgjpsiPeakInToy[iToy];
                break;
            case 4:
                nbkgpsi2sPeakInToy[iToy]=((RooRealVar*)(nbkgpsi2sPeakDataset->get(iToy))->find("nbkgpsi2sPeak"))->getVal();
                nbkgpsi2sPeakInToys+=nbkgpsi2sPeakInToy[iToy];
                break;
            case 11:
                nbkgjpsiPeakInToy[iToy]=((RooRealVar*)(nbkgjpsiPeakDataset->get(iToy))->find("nbkgjpsiPeak"))->getVal();
                nbkgjpsiPeakInToys+=nbkgjpsiPeakInToy[iToy];
                nbkgpsi2sPeakInToy[iToy]=((RooRealVar*)(nbkgpsi2sPeakDataset->get(iToy))->find("nbkgpsi2sPeak"))->getVal();
                nbkgpsi2sPeakInToys+=nbkgpsi2sPeakInToy[iToy];
                break;
            default:
                nbkgjpsiPeakInToy[iToy] = 0;
                nbkgpsi2sPeakInToy[iToy] = 0;
                break;
        }
        nsigInToys+=nsigInToy[iToy];
        nbkgCombInToys+=nbkgCombInToy[iToy];
    }

    // Loop over phase space
    double thisAfb = afb;
    double thisFl = fl;
        // loop over afb for errAfb
    for(int iAfb = 0; iAfb*stepSizeAfb < 2.; iAfb++){
        thisAfb = iAfb*stepSizeAfb-1.+stepSizeAfb/2;
        printf("thisAfb=%+04f\n",thisAfb);
        if (!scanAfbFlPositivePdf(thisAfb,thisFl,true)) continue;

        if (stat(TString::Format("%s/afb%+04.0f_fl%+04.0f",otoyspath.Data(),thisAfb*1000,thisFl*1000).Data(),&fiBuff) != 0){
            mkdir(TString::Format("%s/afb%+04.0f_fl%+04.0f",otoyspath.Data(),thisAfb*1000,thisFl*1000).Data(),0755);
        }
        owspacepath=TString::Format("%s/afb%+04.0f_fl%+04.0f",otoyspath.Data(),thisAfb*1000,thisFl*1000);

        genToySignal(iBin, nsigInToys,thisAfb,thisFl);
        splitToySamples(TString::Format("%s/signalToy_Bin%d.root",owspacepath.Data(),iBin).Data(),nToy,nsigInToy,TString::Format("signalToy_Bin%d",iBin).Data());
        genToyCombBkg(iBin,nbkgCombInToys);
        splitToySamples(TString::Format("%s/combBkgToy_Bin%d.root",owspacepath.Data(),iBin).Data(),nToy,nbkgCombInToy,TString::Format("combBkgToy_Bin%d",iBin).Data());
        switch(iBin){
            case 10:
                genToyPeakBkgFromWspace(iBin, nbkgjpsiPeakInToys,"jpsi");
                splitToySamples(TString::Format("%s/jpsiPeakToy_Bin%d.root",owspacepath.Data(),iBin).Data(),nToy,nbkgjpsiPeakInToy,TString::Format("jpsiPeakToy_Bin%d",iBin).Data());
                break;
            case 4:
                genToyPeakBkgFromWspace(iBin, nbkgpsi2sPeakInToys,"psi2s");
                splitToySamples(TString::Format("%s/psi2sPeakToy_Bin%d.root",owspacepath.Data(),iBin).Data(),nToy,nbkgpsi2sPeakInToy,TString::Format("psi2sPeakToy_Bin%d",iBin).Data());
                break;
            case 11:
                genToyPeakBkgFromWspace(iBin, nbkgjpsiPeakInToys,"jpsi");
                splitToySamples(TString::Format("%s/jpsiPeakToy_Bin%d.root",owspacepath.Data(),iBin).Data(),nToy,nbkgjpsiPeakInToy,TString::Format("jpsiPeakToy_Bin%d",iBin).Data());
                genToyPeakBkgFromWspace(iBin, nbkgpsi2sPeakInToys,"psi2s");
                splitToySamples(TString::Format("%s/psi2sPeakToy_Bin%d.root",owspacepath.Data(),iBin).Data(),nToy,nbkgpsi2sPeakInToy,TString::Format("psi2sPeakToy_Bin%d",iBin).Data());
                break;
            default:
                break;
        }
        for(int iToy = 0; iToy<nToy; iToy++){
            if (stat(TString::Format("%s/set%04d",owspacepath.Data(),iToy+1),&fiBuff) != 0){
                mkdir(TString::Format("%s/set%04d",owspacepath.Data(),iToy+1),0755);
            }
            rename(TString::Format("%s/signalToy_Bin%d_set%04d.root"         ,owspacepath.Data(),       iBin,iToy+1).Data(),
                   TString::Format("%s/set%04d/signalToy_Bin%d_set%04d.root" ,owspacepath.Data(),iToy+1,iBin,iToy+1).Data());
            rename(TString::Format("%s/combBkgToy_Bin%d_set%04d.root"        ,owspacepath.Data(),       iBin,iToy+1).Data(),
                   TString::Format("%s/set%04d/combBkgToy_Bin%d_set%04d.root",owspacepath.Data(),iToy+1,iBin,iToy+1).Data());
            switch(iBin){
                case 10:
                    rename(TString::Format("%s/jpsiPeakToy_Bin%d_set%04d.root"          ,owspacepath.Data(),       iBin,iToy+1).Data(),
                           TString::Format("%s/set%04d/jpsiPeakToy_Bin%d_set%04d.root"  ,owspacepath.Data(),iToy+1,iBin,iToy+1).Data());
                    break;
                case 4:
                    rename(TString::Format("%s/psi2sPeakToy_Bin%d_set%04d.root"         ,owspacepath.Data(),       iBin,iToy+1).Data(),
                           TString::Format("%s/set%04d/psi2sPeakToy_Bin%d_set%04d.root" ,owspacepath.Data(),iToy+1,iBin,iToy+1).Data());
                    break;
                case 11:
                    rename(TString::Format("%s/jpsiPeakToy_Bin%d_set%04d.root"          ,owspacepath.Data(),       iBin,iToy+1).Data(),
                           TString::Format("%s/set%04d/jpsiPeakToy_Bin%d_set%04d.root"  ,owspacepath.Data(),iToy+1,iBin,iToy+1).Data());
                    rename(TString::Format("%s/psi2sPeakToy_Bin%d_set%04d.root"         ,owspacepath.Data(),       iBin,iToy+1).Data(),
                           TString::Format("%s/set%04d/psi2sPeakToy_Bin%d_set%04d.root" ,owspacepath.Data(),iToy+1,iBin,iToy+1).Data());
                    break;
                default:
                    break;
            }
        }
    }

        // 2-D loop for errFl and contour
    thisAfb = afb;
    thisFl = fl;
    for(int iAfb = 0; iAfb*stepSizeAfb < 2.; iAfb++){
        if (contourMode) thisAfb = iAfb*stepSizeAfb-1.+stepSizeAfb/2;

        for(int iFl = 0; iFl*stepSizeFl < 1.; iFl++){
            thisFl = iFl*stepSizeFl+stepSizeFl/2;
            if (!scanAfbFlPositivePdf(thisAfb,thisFl,true)) continue;
            if (stat(TString::Format("%s/afb%+04.0f_fl%+04.0f",otoyspath.Data(),thisAfb*1000,thisFl*1000).Data(),&fiBuff) != 0){
                mkdir(TString::Format("%s/afb%+04.0f_fl%+04.0f",otoyspath.Data(),thisAfb*1000,thisFl*1000).Data(),0755);
            }
            owspacepath=TString::Format("%s/afb%+04.0f_fl%+04.0f",otoyspath.Data(),thisAfb*1000,thisFl*1000);

            genToySignal(iBin, nsigInToys,thisAfb,thisFl);
            splitToySamples(TString::Format("%s/signalToy_Bin%d.root",owspacepath.Data(),iBin).Data(),nToy,nsigInToy,TString::Format("signalToy_Bin%d",iBin).Data());
            genToyCombBkg(iBin,nbkgCombInToys);
            splitToySamples(TString::Format("%s/combBkgToy_Bin%d.root",owspacepath.Data(),iBin).Data(),nToy,nbkgCombInToy,TString::Format("combBkgToy_Bin%d",iBin).Data());
            switch(iBin){
                case 10:
                    genToyPeakBkgFromWspace(iBin, nbkgjpsiPeakInToys,"jpsi");
                    splitToySamples(TString::Format("%s/jpsiPeakToy_Bin%d.root",owspacepath.Data(),iBin).Data(),nToy,nbkgjpsiPeakInToy,TString::Format("jpsiPeakToy_Bin%d",iBin).Data());
                    break;
                case 4:
                    genToyPeakBkgFromWspace(iBin, nbkgpsi2sPeakInToys,"psi2s");
                    splitToySamples(TString::Format("%s/psi2sPeakToy_Bin%d.root",owspacepath.Data(),iBin).Data(),nToy,nbkgpsi2sPeakInToy,TString::Format("psi2sPeakToy_Bin%d",iBin).Data());
                    break;
                case 11:
                    genToyPeakBkgFromWspace(iBin, nbkgjpsiPeakInToys,"jpsi");
                    splitToySamples(TString::Format("%s/jpsiPeakToy_Bin%d.root",owspacepath.Data(),iBin).Data(),nToy,nbkgjpsiPeakInToy,TString::Format("jpsiPeakToy_Bin%d",iBin).Data());
                    genToyPeakBkgFromWspace(iBin, nbkgpsi2sPeakInToys,"psi2s");
                    splitToySamples(TString::Format("%s/psi2sPeakToy_Bin%d.root",owspacepath.Data(),iBin).Data(),nToy,nbkgpsi2sPeakInToy,TString::Format("psi2sPeakToy_Bin%d",iBin).Data());
                    break;
                default:
                    break;
            }
            for(int iToy = 0; iToy<nToy; iToy++){
                if (stat(TString::Format("%s/set%04d",owspacepath.Data(),iToy+1),&fiBuff) != 0){
                    mkdir(TString::Format("%s/set%04d",owspacepath.Data(),iToy+1),0755);
                }
                rename(TString::Format("%s/signalToy_Bin%d_set%04d.root"         ,owspacepath.Data(),       iBin,iToy+1).Data(),
                       TString::Format("%s/set%04d/signalToy_Bin%d_set%04d.root" ,owspacepath.Data(),iToy+1,iBin,iToy+1).Data());
                rename(TString::Format("%s/combBkgToy_Bin%d_set%04d.root"        ,owspacepath.Data(),       iBin,iToy+1).Data(),
                       TString::Format("%s/set%04d/combBkgToy_Bin%d_set%04d.root",owspacepath.Data(),iToy+1,iBin,iToy+1).Data());
                switch(iBin){
                    case 10:
                        rename(TString::Format("%s/jpsiPeakToy_Bin%d_set%04d.root"          ,owspacepath.Data(),       iBin,iToy+1).Data(),
                               TString::Format("%s/set%04d/jpsiPeakToy_Bin%d_set%04d.root"  ,owspacepath.Data(),iToy+1,iBin,iToy+1).Data());
                        break;
                    case 4:
                        rename(TString::Format("%s/psi2sPeakToy_Bin%d_set%04d.root"         ,owspacepath.Data(),       iBin,iToy+1).Data(),
                               TString::Format("%s/set%04d/psi2sPeakToy_Bin%d_set%04d.root" ,owspacepath.Data(),iToy+1,iBin,iToy+1).Data());
                        break;
                    case 11:
                        rename(TString::Format("%s/jpsiPeakToy_Bin%d_set%04d.root"          ,owspacepath.Data(),       iBin,iToy+1).Data(),
                               TString::Format("%s/set%04d/jpsiPeakToy_Bin%d_set%04d.root"  ,owspacepath.Data(),iToy+1,iBin,iToy+1).Data());
                        rename(TString::Format("%s/psi2sPeakToy_Bin%d_set%04d.root"         ,owspacepath.Data(),       iBin,iToy+1).Data(),
                               TString::Format("%s/set%04d/psi2sPeakToy_Bin%d_set%04d.root" ,owspacepath.Data(),iToy+1,iBin,iToy+1).Data());
                        break;
                    default:
                        break;
                }
            }
        }// Fl loop
        if (!contourMode) break;
    }// Afb loop

    return;
 }//}}}

void createFCToys2(int iBin, int nToy=500)// Create toys with profiled likelihood considered
{//{{{
    // Buffers for checking directory/FILE
    struct stat fiBuff;
    FILE *fi = 0;

    // Setting
    gKeepParam = true; // It MUST be true to keep profiled fit result.
    redirectStdout = true; // better keep the log to check how profiled fit works.
    TString iwspacepath_global = iwspacepath;
    TString idatacardpath_global = idatacardpath;
    TString otoyspath = TString::Format("./limit/bin%d",iBin);
    const double stepSizeAfb = 0.01;
    const double stepSizeFl  = 0.01;

    // Create output directory
    if (stat(TString::Format("./limit"),&fiBuff) != 0){
        mkdir(TString::Format("./limit"),0755);
    }
    if (stat(TString::Format("%s",otoyspath.Data()),&fiBuff) != 0){
        mkdir(TString::Format("%s",otoyspath.Data()),0755);
    }

    // Get parameters for the q2 Bin
    double  fl_ubd  = readParam(iBin,"fl_migrad",0);
    double  afb_ubd = readParam(iBin,"afb_migrad",0);
    double  fl      = toBoundedFl(fl_ubd);
    double  afb     = toBoundedAfb(afb_ubd,fl_ubd);

    //  Afb
    for(int iAfb = 0; iAfb*stepSizeAfb < 2.; iAfb++){
        double thisAfb = iAfb*stepSizeAfb-1.+stepSizeAfb/2;

        // FIX ME!!!
        // An alternative way needed to confine the range to be scanned...
        if (fabs(thisAfb) > 0.75) continue;
        //if (!scanAfbFlPositivePdf(thisAfb,thisFl,true)) continue;
        printf("thisAfb=%+04f\n",thisAfb);

        if (stat(TString::Format("%s/afb%+04.0f",otoyspath.Data(),thisAfb*1000).Data(),&fiBuff) != 0){
            mkdir(TString::Format("%s/afb%+04.0f",otoyspath.Data(),thisAfb*1000).Data(),0755);
        }else{
            printf(TString::Format("WARNING\t: Skip %s/afb%+04.0f",otoyspath.Data(),thisAfb*1000).Data());
            continue;
        }
        owspacepath  =TString::Format("%s/afb%+04.0f",otoyspath.Data(),thisAfb*1000);
        odatacardpath=TString::Format("%s/afb%+04.0f",otoyspath.Data(),thisAfb*1000);
        plotpath     =TString::Format("%s/afb%+04.0f",otoyspath.Data(),thisAfb*1000);

        TFile *f_wspace = 0;
        RooWorkspace *wspace = 0;
        do{
            printf("INFO\t\t: Fit with trueAfb=%+.4f\n",thisAfb);
            angular3D_bin(iBin,"angular3D",1,-1,thisAfb);// profiling nuisance with fixed bounded afb.
            f_wspace = new TFile(TString::Format("%s/wspace_angular3D_bin%d.root",owspacepath.Data(),iBin));
            wspace = (RooWorkspace*)f_wspace->Get("wspace");
        }while(wspace == 0);
        idatacardpath = odatacardpath;

        if (readParam(iBin,"migrad",0) != 0 || readParam(iBin,"migrad",1) < -1e20){
            f_wspace->Close();
            continue;
        }

        double fl_ubd  = readParam(iBin,"fl_migrad",0); 
        double thisFl  = toBoundedFl(fl_ubd);
        double thisFs  = readParam(iBin,"fs_migrad",0);
        double thisAs  = toOriginAs(readParam(iBin,"as_migrad",0),thisFs,fl_ubd);

        // sigA  [4]={afb_bdd ,fl_bdd  ,fs     ,as_ori}
        // yields[4]={nSig    ,nBkg    ,nJpsi  ,nJpsiErr ,nPsi2S, nPsi2SErr}
        double sigACoef[4] = {thisAfb,thisFl,thisFs,thisAs};
        double yields  [6] = {wspace->var("nsig")->getVal(),wspace->var("nbkgComb")->getVal(),0,0,0,0};
        switch(iBin){
            //case 10:
            //    yields[2] = wspace->var("nbkgjpsiPeak")->getVal();
            //    yields[3] = wspace->var("nbkgjpsiPeak")->getError();
            //    break;
            //case 4:
            //    yields[4] = wspace->var("nbkgpsi2sPeak")->getVal();
            //    yields[5] = wspace->var("nbkgpsi2sPeak")->getError();
            //    break;
            //case 11:
            //    yields[2] = wspace->var("nbkgjpsiPeak")->getVal();
            //    yields[3] = wspace->var("nbkgjpsiPeak")->getError();
            //    yields[4] = wspace->var("nbkgpsi2sPeak")->getVal();
            //    yields[5] = wspace->var("nbkgpsi2sPeak")->getError();
            //    break;
            default:
                break;
        }
        f_wspace->Close();
        createToy(iBin,nToy,sigACoef,yields);
    }

    // Fl
    for(int iFl = 0; iFl*stepSizeFl < 1.; iFl++){
        double thisFl = iFl*stepSizeFl+stepSizeFl/2;

        // FIX ME!!!
        // An alternative way needed to confine the range to be scanned...
        printf("thisFl=%+04f\n",thisFl);

        if (stat(TString::Format("%s/fl%+04.0f",otoyspath.Data(),thisFl*1000).Data(),&fiBuff) != 0){
            mkdir(TString::Format("%s/fl%+04.0f",otoyspath.Data(),thisFl*1000).Data(),0755);
        }else{
            printf(TString::Format("WARNING\t: Skip %s/fl%+04.0f",otoyspath.Data(),thisFl*1000).Data());
            continue;
        }
        owspacepath  =TString::Format("%s/fl%+04.0f",otoyspath.Data(),thisFl*1000);
        odatacardpath=TString::Format("%s/fl%+04.0f",otoyspath.Data(),thisFl*1000);
        plotpath     =TString::Format("%s/fl%+04.0f",otoyspath.Data(),thisFl*1000);
        TFile *f_wspace = 0;
        RooWorkspace *wspace = 0;

        do{
            printf("INFO\t\t: Fit with trueFl=%+.4f\n",thisFl);
            angular3D_bin(iBin,"angular3D",1,thisFl,-1);// profiling nuisance with fixed bounded afb.
            f_wspace = new TFile(TString::Format("%s/wspace_angular3D_bin%d.root",owspacepath.Data(),iBin));
            wspace = (RooWorkspace*)f_wspace->Get("wspace");
        }while(wspace == 0);
        idatacardpath = odatacardpath;

        if (readParam(iBin,"migrad",0) != 0 || readParam(iBin,"migrad",1) < -1e20){
            f_wspace->Close();
            continue;
        }

        double  fl_ubd  = toUnboundedFl(thisFl); 
        double afb_ubd  = readParam(iBin,"afb_migrad",0);
        double thisAfb  = toBoundedAfb(afb_ubd,fl_ubd);
        double thisFs   = readParam(iBin,"fs_migrad",0);
        double thisAs   = toOriginAs(readParam(iBin,"as_migrad",0),thisFs,fl_ubd);

        // sigA  [4]={afb_bdd ,fl_bdd  ,fs     ,as_ori}
        // yields[4]={nSig    ,nBkg    ,nJpsi  ,nJpsiErr ,nPsi2S, nPsi2SErr}
        double sigACoef[4] = {thisAfb,thisFl,thisFs,thisAs};
        double yields  [6] = {wspace->var("nsig")->getVal(),wspace->var("nbkgComb")->getVal(),0,0,0,0};
        switch(iBin){
            //case 10:
            //    yields[2] = wspace->var("nbkgjpsiPeak")->getVal();
            //    yields[3] = wspace->var("nbkgjpsiPeak")->getError();
            //    break;
            //case 4:
            //    yields[4] = wspace->var("nbkgpsi2sPeak")->getVal();
            //    yields[5] = wspace->var("nbkgpsi2sPeak")->getError();
            //    break;
            //case 11:
            //    yields[2] = wspace->var("nbkgjpsiPeak")->getVal();
            //    yields[3] = wspace->var("nbkgjpsiPeak")->getError();
            //    yields[4] = wspace->var("nbkgpsi2sPeak")->getVal();
            //    yields[5] = wspace->var("nbkgpsi2sPeak")->getError();
            //    break;
            default:
                break;
        }
        f_wspace->Close();
        createToy(iBin,nToy,sigACoef,yields);
    }
    return;
}//}}}

void scanProfile(int iBin, const char outfile[]="scanProfile")
{//{{{
    // This is used when we try fixing a value.
    TString oProfilePath = TString::Format("./limit/bin%d",iBin);
    const double stepSizeAfb = 0.01;
    const double stepSizeFl  = 0.01;

    TCanvas *canvas = new TCanvas("canvas");
    TGraph *gr_profiledFl_trueAfb = new TGraph();
    TGraph *gr_profiledFs_trueAfb = new TGraph();
    TGraph *gr_profiledAs_trueAfb = new TGraph();
    TGraph *gr_profiledAfb_trueFl = new TGraph();
    TGraph *gr_profiledFs_trueFl  = new TGraph();
    TGraph *gr_profiledAs_trueFl  = new TGraph();
    int nProfiledFl = 0;
    int nProfiledAfb= 0;

    struct stat fiBuff;
    for(int iAfb = 0; iAfb*stepSizeAfb < 2.; iAfb++){
        double thisAfb = iAfb*stepSizeAfb-1.+stepSizeAfb/2;

        if (stat(TString::Format("%s/afb%+04.0f",oProfilePath.Data(),thisAfb*1000).Data(),&fiBuff) != 0) continue;
        printf("thisAfb=%+04f\n",thisAfb);
        nProfiledFl++;

        idatacardpath  =TString::Format("%s/afb%+04.0f",oProfilePath.Data(),thisAfb*1000);
        if (readParam(iBin,"migrad",0) != 0    ) continue;
        if (readParam(iBin,"migrad",1) <  -1e20) continue;//in case of FAILED
        double fl_ubd  = readParam(iBin,"fl_migrad",0);
        double fs      = readParam(iBin,"fs_migrad",0);
        double as_tr   = readParam(iBin,"as_migrad",0);
        gr_profiledFl_trueAfb->SetPoint(nProfiledFl,thisAfb,toBoundedFl(fl_ubd));
        gr_profiledFs_trueAfb->SetPoint(nProfiledFl,thisAfb,fs);
        gr_profiledAs_trueAfb->SetPoint(nProfiledFl,thisAfb,toOriginAs(as_tr,fs,fl_ubd));
        
    }
    for(int iFl = 0; iFl*stepSizeFl < 1.; iFl++){
        double thisFl = iFl*stepSizeFl+stepSizeFl/2;

        if (stat(TString::Format("%s/fl%+04.0f",oProfilePath.Data(),thisFl*1000).Data(),&fiBuff) != 0) continue;
        printf("thisFl=%+04f\n",thisFl);
        nProfiledAfb++;

        idatacardpath  =TString::Format("%s/fl%+04.0f",oProfilePath.Data(),thisFl*1000);
        if (readParam(iBin,"migrad",0) != 0    ) continue;
        if (readParam(iBin,"migrad",1) <  -1e20) continue;//in case of FAILED
        double afb_ubd = readParam(iBin,"afb_migrad",0);
        double fl_ubd  = readParam(iBin,"fl_migrad",0);
        double fs      = readParam(iBin,"fs_migrad",0);
        double as_tr   = readParam(iBin,"as_migrad",0);
        gr_profiledAfb_trueFl->SetPoint(nProfiledAfb,thisFl,toBoundedAfb(afb_ubd,fl_ubd));
        gr_profiledFs_trueFl ->SetPoint(nProfiledAfb,thisFl,fs);
        gr_profiledAs_trueFl ->SetPoint(nProfiledAfb,thisFl,toOriginAs(as_tr,fs,fl_ubd));
        
    }

    gr_profiledFl_trueAfb->GetXaxis()->SetTitle("True A_{FB}");
    gr_profiledFl_trueAfb->GetYaxis()->SetTitle("Profiled F_{L}");
    gr_profiledFl_trueAfb->GetXaxis()->SetLimits(-0.75,0.75);
    gr_profiledFl_trueAfb->GetYaxis()->SetRangeUser(0,1);
    gr_profiledFl_trueAfb->SetMarkerStyle(7);
    gr_profiledFl_trueAfb->SetMarkerSize(0.4);
    gr_profiledFl_trueAfb->Draw("AP");
    canvas->Print(TString::Format("./plots/%sFl_trueAfb_bin%d.pdf",outfile,iBin).Data());

    gr_profiledFs_trueAfb->GetXaxis()->SetTitle("True A_{FB}");
    gr_profiledFs_trueAfb->GetYaxis()->SetTitle("Profiled F_{s}");
    gr_profiledFs_trueAfb->GetXaxis()->SetLimits(-0.75,0.75);
    gr_profiledFs_trueAfb->GetYaxis()->SetRangeUser(0,1);
    gr_profiledFs_trueAfb->SetMarkerStyle(7);
    gr_profiledFs_trueAfb->SetMarkerSize(0.4);
    gr_profiledFs_trueAfb->Draw("AP");
    canvas->Print(TString::Format("./plots/%sFs_trueAfb_bin%d.pdf",outfile,iBin).Data());

    gr_profiledAs_trueAfb->GetXaxis()->SetTitle("True A_{FB}");
    gr_profiledAs_trueAfb->GetYaxis()->SetTitle("Profiled A_{s}");
    gr_profiledAs_trueAfb->GetXaxis()->SetLimits(-0.75,0.75);
    gr_profiledAs_trueAfb->GetYaxis()->SetRangeUser(-1,1);
    gr_profiledAs_trueAfb->SetMarkerStyle(7);
    gr_profiledAs_trueAfb->SetMarkerSize(0.4);
    gr_profiledAs_trueAfb->Draw("AP");
    canvas->Print(TString::Format("./plots/%sAs_trueAfb_bin%d.pdf",outfile,iBin).Data());

    gr_profiledAfb_trueFl->GetXaxis()->SetTitle("True F_{L}");
    gr_profiledAfb_trueFl->GetYaxis()->SetTitle("Profiled A_{FB}");
    gr_profiledAfb_trueFl->GetXaxis()->SetLimits(0,1.);
    gr_profiledAfb_trueFl->GetYaxis()->SetRangeUser(-0.75,0.75);
    gr_profiledAfb_trueFl->SetMarkerStyle(7);
    gr_profiledAfb_trueFl->SetMarkerSize(0.4);
    gr_profiledAfb_trueFl->Draw("AP");
    canvas->Print(TString::Format("./plots/%sAfb_trueFl_bin%d.pdf",outfile,iBin).Data());

    gr_profiledFs_trueFl->GetXaxis()->SetTitle("True F_{L}");
    gr_profiledFs_trueFl->GetYaxis()->SetTitle("Profiled F_{s}");
    gr_profiledFs_trueFl->GetXaxis()->SetLimits(0,1);
    gr_profiledFs_trueFl->GetYaxis()->SetRangeUser(0,1);
    gr_profiledFs_trueFl->SetMarkerStyle(7);
    gr_profiledFs_trueFl->SetMarkerSize(0.4);
    gr_profiledFs_trueFl->Draw("AP");
    canvas->Print(TString::Format("./plots/%sFs_trueFl_bin%d.pdf",outfile,iBin).Data());

    gr_profiledAs_trueFl->GetXaxis()->SetTitle("True F_{L}");
    gr_profiledAs_trueFl->GetYaxis()->SetTitle("Profiled A_{s}");
    gr_profiledAs_trueFl->GetXaxis()->SetLimits(0,1.);
    gr_profiledAs_trueFl->GetYaxis()->SetRangeUser(-1,1);
    gr_profiledAs_trueFl->SetMarkerStyle(7);
    gr_profiledAs_trueFl->SetMarkerSize(0.4);
    gr_profiledAs_trueFl->Draw("AP");
    canvas->Print(TString::Format("./plots/%sAs_trueFl_bin%d.pdf",outfile,iBin).Data());
    
    return;
}//}}}

void printProfile(int iBin)
{//{{{
    TString oProfilePath = TString::Format("./limit/bin%d",iBin);
    const double stepSizeAfb = 0.01;
    const double stepSizeFl  = 0.01;

    TCanvas *canvas = new TCanvas("canvas");

    struct stat fiBuff;
    for(int iAfb = 0; iAfb*stepSizeAfb < 2.; iAfb++){
        double thisAfb = iAfb*stepSizeAfb-1.+stepSizeAfb/2;

        if (stat(TString::Format("%s/afb%+04.0f",oProfilePath.Data(),thisAfb*1000).Data(),&fiBuff) != 0) continue;
        printf("thisAfb=%+04f\n",thisAfb);

        iwspacepath    = TString::Format("%s/afb%+04.0f",oProfilePath.Data(),thisAfb*1000);
        idatacardpath  = TString::Format("%s/afb%+04.0f",oProfilePath.Data(),thisAfb*1000);
        plotpath       = TString::Format("%s/afb%+04.0f",oProfilePath.Data(),thisAfb*1000);
        if (readParam(iBin,"migrad",0) != 0    ) continue;
        if (readParam(iBin,"migrad",1) <  -1e20) continue;//in case of FAILED

        TFile *f_wspace = new TFile(TString::Format("%s/wspace_angular3D_bin%d.root",iwspacepath.Data(),iBin));
        RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");

        RooAbsPdf *pdf = wspace->pdf("f_sigA");
        RooRealVar *cosK=wspace->var("CosThetaK");
        RooRealVar *cosL=wspace->var("CosThetaL");
        printPdf2D(*pdf,*cosK,*cosL,TString::Format("f_sigA_bin%d",iBin).Data());

        f_wspace->Close();
    }
    for(int iFl = 0; iFl*stepSizeFl < 1.; iFl++){
        double thisFl = iFl*stepSizeFl;

        if (stat(TString::Format("%s/fl%+04.0f",oProfilePath.Data(),thisFl*1000).Data(),&fiBuff) != 0) continue;
        printf("thisFl=%+04f\n",thisFl);

        iwspacepath    = TString::Format("%s/fl%+04.0f",oProfilePath.Data(),thisFl*1000);
        idatacardpath  = TString::Format("%s/fl%+04.0f",oProfilePath.Data(),thisFl*1000);
        plotpath       = TString::Format("%s/fl%+04.0f",oProfilePath.Data(),thisFl*1000);
        if (readParam(iBin,"migrad",0) != 0    ) continue;
        if (readParam(iBin,"migrad",1) <  -1e20) continue;//in case of FAILED

        TFile *f_wspace = new TFile(TString::Format("%s/wspace_angular3D_bin%d.root",iwspacepath.Data(),iBin));
        RooWorkspace *wspace = (RooWorkspace*)f_wspace->Get("wspace");

        RooAbsPdf *pdf = wspace->pdf("f_sigA");
        RooRealVar *cosK=wspace->var("CosThetaK");
        RooRealVar *cosL=wspace->var("CosThetaL");
        printPdf2D(*pdf,*cosK,*cosL,TString::Format("f_sigA_bin%d",iBin).Data());

        f_wspace->Close();
    }
}//}}}

void createProfileToys(int iBin, int nToy=500)
{//{{{
    // Buffers for checking directory/FILE
    struct stat fiBuff;
    FILE *fi = 0;

    // Setting
    gKeepParam     = true; // It MUST be true to keep profiled fit result.
    redirectStdout = true; // better keep the log to check how profiled fit works.
    TString checkTag = "convRateApproach";
    TString otoyspath = TString::Format("./limit/validation/%s_bin%d",checkTag.Data(),iBin);

    // Create output directory
    if (stat(TString::Format("./limit/validation"),&fiBuff) != 0){
        mkdir(TString::Format("./limit/validation"),0755);
    }
    if (stat(TString::Format("%s",otoyspath.Data()),&fiBuff) != 0){
        mkdir(TString::Format("%s",otoyspath.Data()),0755);
    }

    // Get parameters for the q2 Bin
    double  fl_ubd  = readParam(iBin,"fl_migrad",0);
    double  afb_ubd = readParam(iBin,"afb_migrad",0);
    double  fl      = toBoundedFl(fl_ubd);
    double  afb     = toBoundedAfb(afb_ubd,fl_ubd);
    double  as_tr   = readParam(iBin,"as_migrad",0);
    double  fs      = readParam(iBin,"fs_migrad",0);
    double  as      = toOriginAs(as_tr,fs,fl_ubd);

    double sigACoef[4];
    sigACoef[0]=afb;
    sigACoef[1]=fl;
    sigACoef[2]=fs;
    sigACoef[3]=as;

    printf("DEBUG\t\t: afb=%f, fl=%f, fs=%f, as=%f\n",afb,fl,fs,as);

    double yields[6];
    TFile *f_wspace = 0;
    RooWorkspace *wspace = 0;

    f_wspace = new TFile(TString::Format("%s/wspace_angular3D_bin%d.root",iwspacepath.Data(),iBin));
    wspace = (RooWorkspace*)f_wspace->Get("wspace");

    RooRealVar *nsig = wspace->var("nsig");
    RooRealVar *nbkgComb      = wspace->var("nbkgComb");
    RooRealVar *nbkgjpsiPeak  = wspace->var("nbkgjpsiPeak");
    RooRealVar *nbkgpsi2sPeak = wspace->var("nbkgpsi2sPeak");
    yields[0]=nsig->getVal();
    yields[1]=nbkgComb->getVal();
    yields[2]=nbkgjpsiPeak  == NULL ? 0    : nbkgjpsiPeak->getVal();
    yields[3]=nbkgjpsiPeak  == NULL ? 0.01 : nbkgjpsiPeak->getError();
    yields[4]=nbkgpsi2sPeak == NULL ? 0    : nbkgpsi2sPeak->getVal();
    yields[5]=nbkgpsi2sPeak == NULL ? 0.01 : nbkgpsi2sPeak->getError();

    // Start generating
    printf("INFO\t\t: Start Generating!\n");
    for(double iFl=0.005; iFl<1; iFl+=0.1){
        sigACoef[1]=iFl;
        owspacepath = TString::Format("%s/afb%+04.0f_fl%+04.0f",otoyspath.Data(),1000*sigACoef[0],1000*sigACoef[1]);
        if (stat(owspacepath,&fiBuff) != 0){
            mkdir(owspacepath,0755);
        }
        if (!scanAfbFlPositivePdf(sigACoef[0], sigACoef[1], false) ) continue;
        createToy(iBin,500,sigACoef,yields);
    }
    f_wspace->Close();
    return;
}//}}}

void createBestFitToys(int iBin, int nToy=500)
{//{{{
    // Buffers for checking directory/FILE
    struct stat fiBuff;
    FILE *fi = 0;

    // Setting
    gKeepParam     = true; // It MUST be true to keep profiled fit result.
    redirectStdout = true; // better keep the log to check how profiled fit works.
    TString otoyspath = TString::Format("./limit/bestFit");

    // Create output directory
    if (stat(TString::Format("%s",otoyspath.Data()),&fiBuff) != 0){
        mkdir(TString::Format("%s",otoyspath.Data()),0755);
    }

    // Get parameters for the q2 Bin
    double  fl_ubd  = readParam(iBin,"fl_migrad",0);
    double  afb_ubd = readParam(iBin,"afb_migrad",0);
    double  fl      = toBoundedFl(fl_ubd);
    double  afb     = toBoundedAfb(afb_ubd,fl_ubd);
    double  as_tr   = readParam(iBin,"as_migrad",0);
    double  fs      = readParam(iBin,"fs_migrad",0);
    double  as      = toOriginAs(as_tr,fs,fl_ubd);

    double sigACoef[4];
    sigACoef[0]=afb;
    sigACoef[1]=fl;
    sigACoef[2]=fs;
    sigACoef[3]=as;

    printf("DEBUG\t\t: afb=%f, fl=%f, fs=%f, as=%f\n",afb,fl,fs,as);

    double yields[6];
    TFile *f_wspace = 0;
    RooWorkspace *wspace = 0;

    f_wspace = new TFile(TString::Format("%s/wspace_angular3D_bin%d.root",iwspacepath.Data(),iBin));
    wspace = (RooWorkspace*)f_wspace->Get("wspace");

    RooRealVar *nsig = wspace->var("nsig");
    RooRealVar *nbkgComb      = wspace->var("nbkgComb");
    RooRealVar *nbkgjpsiPeak  = wspace->var("nbkgjpsiPeak");
    RooRealVar *nbkgpsi2sPeak = wspace->var("nbkgpsi2sPeak");
    yields[0]=nsig->getVal();
    yields[1]=nbkgComb->getVal();
    yields[2]=nbkgjpsiPeak  == NULL ? 0    : nbkgjpsiPeak->getVal();
    yields[3]=nbkgjpsiPeak  == NULL ? 0.01 : nbkgjpsiPeak->getError();
    yields[4]=nbkgpsi2sPeak == NULL ? 0    : nbkgpsi2sPeak->getVal();
    yields[5]=nbkgpsi2sPeak == NULL ? 0.01 : nbkgpsi2sPeak->getError();

    // Start generating
    printf("INFO\t\t: Start Generating!\n");
    owspacepath = TString::Format("%s",otoyspath.Data());
    if (stat(owspacepath,&fiBuff) != 0){
        mkdir(owspacepath,0755);
    }
    createToy(iBin,nToy,sigACoef,yields);

    f_wspace->Close();
    return;
}//}}}

void createEdgeEffectToys(int iBin, int nToy=500)
{//{{{
    // Buffers for checking directory/FILE
    struct stat fiBuff;
    FILE *fi = 0;

    // Setting
    gKeepParam     = false; // It MUST be true to keep profiled fit result.
    redirectStdout = true; // better keep the log to check how profiled fit works.
    TString otoyspath = TString::Format("./limit/validation/edgeEffect");

    // Create output directory
    if (stat(TString::Format("%s",otoyspath.Data()),&fiBuff) != 0){
        mkdir(TString::Format("%s",otoyspath.Data()),0755);
    }

    // Get parameters for the q2 Bin
    double  fl_ubd  = readParam(iBin,"fl_migrad",0);
    double  afb_ubd = readParam(iBin,"afb_migrad",0);
    double  fl      = toBoundedFl(fl_ubd);
    double  afb     = toBoundedAfb(afb_ubd,fl_ubd);
    double  as_tr   = readParam(iBin,"as_migrad",0);
    double  fs      = readParam(iBin,"fs_migrad",0);
    double  as      = toOriginAs(as_tr,fs,fl_ubd);

    double sigACoef[4];
    sigACoef[0]=afb;
    sigACoef[1]=fl;
    sigACoef[2]=fs;
    sigACoef[3]=as;

    printf("DEBUG\t\t: afb=%f, fl=%f, fs=%f, as=%f\n",afb,fl,fs,as);

    double yields[6];
    TFile *f_wspace = 0;
    RooWorkspace *wspace = 0;

    f_wspace = new TFile(TString::Format("%s/wspace_angular3D_bin%d.root",iwspacepath.Data(),iBin));
    wspace = (RooWorkspace*)f_wspace->Get("wspace");

    RooRealVar *nsig = wspace->var("nsig");
    RooRealVar *nbkgComb      = wspace->var("nbkgComb");
    RooRealVar *nbkgjpsiPeak  = wspace->var("nbkgjpsiPeak");
    RooRealVar *nbkgpsi2sPeak = wspace->var("nbkgpsi2sPeak");
    yields[0]=nsig->getVal();
    yields[1]=nbkgComb->getVal();
    yields[2]=nbkgjpsiPeak  == NULL ? 0    : nbkgjpsiPeak->getVal();
    yields[3]=nbkgjpsiPeak  == NULL ? 0.01 : nbkgjpsiPeak->getError();
    yields[4]=nbkgpsi2sPeak == NULL ? 0    : nbkgpsi2sPeak->getVal();
    yields[5]=nbkgpsi2sPeak == NULL ? 0.01 : nbkgpsi2sPeak->getError();

    // Start generating
    printf("INFO\t\t: Start Generating!\n");
    for(int iScale=1; iScale<4; iScale++){
        owspacepath = TString::Format("%s/yieldsScaleE%d",otoyspath.Data(),iScale);
        if (stat(owspacepath,&fiBuff) != 0){
            mkdir(owspacepath,0755);
        }
        createToy(iBin,nToy,sigACoef,yields);
        for(int iY=0; iY<6; iY++){
            if (iY == 3 || iY == 5) continue;
            yields[iY] = yields[iY] == 0 ? 0 : yields[iY]*10;
        }
    }

    f_wspace->Close();
    return;
}//}}}

// Other tools

//_________________________________________________________________________________

void printListOfTChainElements(TChain *chain){
    TObjArray *fileElements=chain->GetListOfFiles();
    int nFiles = fileElements->GetEntries();
    //printf("DEBUG\t\t: %d files in the chain\n",nFiles);
    TIter next(fileElements);
    TChainElement *chEl=0;
    for( int entry=0; entry < nFiles; entry++ ) {
        chEl=(TChainElement*)next();
        printf("%s\n",chEl->GetTitle());
    }
    printf("DEBUG\t\t: %d files in the chain\n",nFiles);
}

void printHelpMessage(){
    printf("Usage              : ./fit Function infile [--options <arguments>]\n");
    printf("Functions          :\n");
    printf("    bmass                       Fit to mass spectrum using a double Gaussian signal and Chebyshev bkg.\n");
    printf("    fl_gen                      Derive F_{L} from cosThetaK distribution at GEN level.\n");
    printf("    angular_gen                 Derive F_{L} and A_{FB} from angular distribution at GEN level.\n");
    printf("    accXrecoEff                 Get 2D efficiency map from signal simulation.\n");
    printf("    angular2D                   Same as angular_gen, but fit to signal MC at RECO level with efficiency correction, bkg component is NOT considered.\n");
    printf("    angular3D_1a_Sm             Leading step1 to angular3D, determine signal shape from simulation.\n");
    printf("    angular3D_1b_YpPm           Leading step2 to angular3D, determine mass spectrum of peaking bkg from simulation.\n");
    printf("    angular3D_2a_PkPl           Leading step3 to angular3D, determine angular dist. of peaking bkg from simulation.\n");
    printf("    angular3D_prior             Leading step4 to angular3D, fit to data sideband to get initial values of combinatorial bkg.\n");
    //printf("    angular2D_data        Same as 'angular3D', but only fit to angular axes. DEPRECATED for bad distinguish signal/background.\n");
    printf("    angular3D_bins              Derive F_{L} and A_{FB} by fitting to mass and angular distribution for each q2 bin.\n");
    printf("    angular3D_bins_summary      Derive F_{L} and A_{FB} by fitting to mass and angular distribution for each q2 bin.\n");
    printf("    angular3D                   Draw F_{L} and A_{FB} values in q2 bins.\n");
    printf("Expert functions   :\n")  ;
    printf("    createToys                  Generate toys or split MC samples. Please check source code.\n");
    printf("    scanNLL                     In case of low statistics, scan log(L) for better error determination. (VERY SLOW!)\n");
    printf("    createFCToys                Create toys for Feldman-Cousins method (Takes a long while!)\n");
    printf("    createBestFitToys           Create toys with best fit point for coverage test (Takes a long while!)\n");
    printf("    FCScan                      Feldman-Cousins method for better error determination under the constraints of physical boundaries.\n");
    printf("Options     :\n");
    printf("    --help                      Show this help message.\n");
    printf("    --keeplog                   Redirect some messages to odatacardpath/???.log and odatacardpath/???.err.\n");
    printf("    --keepparam                 Keep fit parameters.\n");
    printf("    --plotpath                  Path to output plots.                               [ Default: \"./plots\" ]\n");
    printf("    --iwspacepath               Path to input  wspaces.                             [ Default: \".\" ]\n");
    printf("    --iCombBkgWspacepath        Path to input  wspaces of combinatorial background. [ Default: \".\" ]\n");
    printf("    --owspacepath               Path to output wspaces.                             [ Default: \".\" ]\n");
    printf("    --idatacardpath             Path to input  datacards.                           [ Default: \".\" ]\n");
    printf("    --odatacardpath             Path to output datacards.                           [ Default: \".\" ]\n");
    printf("    --iallpath                  Path to input  datacards/wspaces.                   [ Default: \".\" ]\n");
    printf("    --oallpath                  Path to output plots/datacards/wspaces.             [ Default: \".\" ]\n");
    printf("    --scale                     Modify the scale factor of input. Negative for toys.[ Default: 1 ]\n");
    printf("Remark             :\n");
    printf("    1. Outputs will be by default stored in ./plots, please keep the directory.\n");
    printf("    2. Wildcard is allowed for infile. But you must quote infile like \"inputData_Run*.root\"!\n");
}

int main(int argc, char** argv) {
    srand(time(NULL));
    RooRandom::randomGenerator()->SetSeed(rand());
    //RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    gStyle->SetOptStat(0);

    // Tags
    is7TeVCheck = false;   
    // less than 0 noCut, 1 for CDF . 2 for LHCb . 3 for 16Aug reOptimization . 4 for sel_v3p5
    isToy = false;
    isToy ? isCDFcut*=-1 : isCDFcut=4;
    int nToy=500;
    int nWorkBins = 3;
    int workBins[] = {10,4,9};

    // Help message
    if (argc <= 2) {
        printHelpMessage();
        return 0;
    }

    // main
    TString func    = argv[1];
    TString infile  = argv[2];
    
    // Processing arguments
    char* l_opt_arg;
    const char short_options[] = "hl";
    static struct option long_options[] = {
        {"help"               , 0 , NULL ,'h'} ,
        {"keeplog"            , 0 , NULL ,'l'} ,
        {"keepparam"          , 0 , NULL ,'p'} ,
        {"plotpath"           , 1 , NULL , 1 } ,
        {"idatacardpath"      , 1 , NULL , 2 } ,
        {"odatacardpath"      , 1 , NULL , 3 } ,
        {"iwspacepath"        , 1 , NULL , 4 } ,
        {"iCombBkgWspacepath" , 1 , NULL , 5 } ,
        {"owspacepath"        , 1 , NULL , 6 } ,
        {"iallpath"           , 1 , NULL , 7 } ,
        {"oallpath"           , 1 , NULL , 8 } ,
        {"scale"              , 1 , NULL , 9 } ,
        {NULL                 , 0 , NULL , 0 }
    };
    int arg;
    while ((arg = getopt_long(argc,argv,short_options,long_options,NULL)) != -1){
        switch (arg) {  
            case 'h':   
                printHelpMessage();
                return 0;
                break;  
            case 'l':
                struct stat fiBuff;
                if (stat("/dev/tty",&fiBuff)==0) redirectStdout = true;
                printf("INFO\t\t: Stdout/stderr are redirected to log files...\n");
                break;
            case 'p':
                gKeepParam = true;
                break;
            case 1:   
                l_opt_arg = optarg;
                plotpath=l_opt_arg;
                break;  
            case 2:   
                l_opt_arg = optarg;
                idatacardpath=l_opt_arg;
                break;  
            case 3:   
                l_opt_arg = optarg;
                odatacardpath=l_opt_arg;
                break;  
            case 4:   
                l_opt_arg = optarg;
                iwspacepath=l_opt_arg;
                iCombBkgWspacepath=l_opt_arg;
                break;  
            case 5:
                l_opt_arg = optarg;
                iCombBkgWspacepath=l_opt_arg;
                break;
            case 6:   
                l_opt_arg = optarg;
                owspacepath=l_opt_arg;
                break;  
            case 7:   
                l_opt_arg = optarg;
                idatacardpath=l_opt_arg;
                iwspacepath=l_opt_arg;
                iCombBkgWspacepath=l_opt_arg;
                break;  
            case 8:
                l_opt_arg = optarg;
                plotpath=l_opt_arg;
                odatacardpath=l_opt_arg;
                owspacepath=l_opt_arg;
                break;
            case 9:
                l_opt_arg = optarg;
                scaleFactor=atof(l_opt_arg);
                if (scaleFactor < 0){
                    isToy=true;
                    scaleFactor=fabs(scaleFactor);
                }
                break;
            default:
                printf("WARNING\t\t: %d is not a valid argument. Program terminates!\n",arg);
                abort();
        }
    }
    printf("INFO\t\t: ScaleFactor for input data is %.3f\n",scaleFactor);
    printf("INFO\t\t: Plots will be stored to %s\n",plotpath.Data());
    printf("INFO\t\t: Datacards will be stored to %s\n",odatacardpath.Data());
    printf("INFO\t\t: Workspaces will be stored to %s\n",owspacepath.Data());
  
    if (func == "bmass") {
        if (argc != 4){
            printf("./fit bmass infile binID\n");
            for (int i = 0; i < nQ2Ranges; i++) {
                printf("    Bin %d : %s\n",i,q2range[i]);
            }
            return 0;
        }
        int iBin = atoi(argv[3]);
        ch->Add(infile.Data());
        if (ch == NULL) return 1;
        const char outfile[]="bmass";
        bmass(iBin, outfile); 
    }else if (func == "fl_gen"){
        ch->Add(infile.Data());
        if (ch == NULL) return 1;
        const char outfile[]="fl_gen";
        fl_gen(outfile);
    }else if (func == "angular_gen"){
        ch->Add(infile.Data());
        if (ch == NULL) return 1;
        const char outfile[]="angular_gen";
        angular_gen(outfile);
    }else if (func == "accXrecoEff") {
        // For experts only. Set to true if no given acceptance_8TeV.root and modify the input filepath in the funciton.
        if (false) createAcceptanceHist();

        ch->Add(infile.Data());
        if (ch == NULL) return 1;
        // TODO: accXrecoEff2 returns all parameters 0 for n>=2.
        //int nTempWorkBins = 5;
        //int tempWorkBins[] = {10,4,9,11,12};
        int nTempWorkBins = 1;
        int tempWorkBins[] = {12};
        for (int iBin = 0; iBin < nTempWorkBins; iBin++) {
            accXrecoEff2(tempWorkBins[iBin],gKeepParam);
        }
    }else if (func == "angular2D"){
        static char wantDoFit[10];
        while(strcmp(wantDoFit,"y")*strcmp(wantDoFit, "n") != 0){
            printf("Do you want to redo fitting? [y/n]:");
            scanf("%19s",wantDoFit);
        }
        bool doFit = false;
        if (strcmp(wantDoFit,"y")==0){
            doFit=true;
            ch->Add(infile.Data());
            if (ch == NULL) return 1;
        }
        const char outfile[]="angular2D";
        angular2D(outfile,doFit);
        //for (int iBin = 2; iBin < 9; iBin++) {
        //    if (iBin == 3 || iBin == 5) continue;
        //    angular2D_bin(iBin);
        //}
    }else if (func == "angular3D_1a_Sm" || func == "angular3D_1b_YpPm" || func == "angular3D_2a_PkPl" || func == "angular3D_prior"){
        bool keepParam=false;
        if (!gKeepParam){
            static char wantKeepParam[10];
            while(strcmp(wantKeepParam,"y")*strcmp(wantKeepParam, "n") != 0){
                printf("Do you want to keep fit parameters? [y/n]:");
                scanf("%19s",wantKeepParam);
            }
            if (strcmp(wantKeepParam,"y")==0) keepParam=true;
        }else{
            keepParam=gKeepParam;
        }

        void (*fx)(int, const char*, bool);
        if ( func == "angular3D_1a_Sm" ){
            fx = angular3D_1a_Sm;
        }else if (func == "angular3D_1b_YpPm"){
            fx = angular3D_1b_YpPm;
        }else if (func == "angular3D_2a_PkPl"){
            fx = angular3D_2a_PkPl;
        }else if (func == "angular3D_prior"){
            fx = angular3D_prior3;
        }
        ch->Add(infile.Data());
        if (ch == NULL) return 1;
        for (int iBin = 0; iBin < nWorkBins; iBin++) {
            fx(workBins[iBin],func,keepParam); // By default NOT overwrite exist parameters.
        }
        for (int iBin = 0; iBin < nSummaryBins; iBin++) {
            fx(summaryBins[iBin],func,keepParam); // By default NOT overwrite exist parameters.
        }
    }else if (func == "angular2D_data"){
        ch->Add(infile.Data());
        if (ch == NULL) return 1;
        const char outfile[]="angular2D_data";
        for (int iBin = 0; iBin < nWorkBins; iBin++) {
            angular2D_data_bin(workBins[iBin],outfile);
        }
    }else if (func == "angular3D_bins_summary"){ // same as angular3D, but run summary bins.
        ch->Add(infile.Data());
        if (ch == NULL) return 1;
        const char outfile[]="angular3D";
        for(int iBin=0; iBin<nWorkBins; iBin++){
            angular3D_bin(summaryBins[iBin],outfile,scaleFactor);
        }
    }else if (func == "angular3D_bins"){
        ch->Add(infile.Data());
        if (ch == NULL) return 1;
        const char outfile[]="angular3D";
        for(int iBin=0; iBin<nWorkBins; iBin++){
            angular3D_bin(workBins[iBin],outfile,scaleFactor);
        }
        for (int iBin = 0; iBin < nSummaryBins; iBin++) {
            angular3D_bin(summaryBins[iBin],outfile,scaleFactor); // By default NOT overwrite exist parameters.
        }
    }else if (func == "angular3D_bins_summary"){ // same as angular3D, but run summary bins.
        ch->Add(infile.Data());
        if (ch == NULL) return 1;
        const char outfile[]="angular3D";
        for(int iBin=0; iBin<nSummaryBins; iBin++){
            angular3D_bin(summaryBins[iBin],outfile,scaleFactor);
        }
    }else if (func == "angular3D"){
        const char outfile[]="angular3D";
        angular3D(outfile);
    }else if (func == "scanNLL"){
        int whichBin=-1;
        while(whichBin < 0 || whichBin > 10 || whichBin == 3 || whichBin == 5){
            printf("Which bin do you want to scan? ");
            scanf("%d",&whichBin);
        }
        ch->Add(infile.Data());
        if (ch == NULL) return 1;

        const char outfile[]="scanNLL";
        printf("INFO\t\t: Call scanNLL for bin%d\n",whichBin);
        scanNLL(whichBin,outfile);
        //getErrFromNllScan(whichBin,true);
    }else if (func == "createFCToys"){
        int whichBin=-1;
        while(whichBin < 0 || whichBin == 3 || whichBin == 5){
            printf("INFO\t\t: This function takes a long while. Toys will be stored to ./limit\n");
            printf("    \t\t  Which bin do you need? ");
            scanf("%d",&whichBin);
        }
        ch->Add(infile.Data());
        if (ch == NULL) return 1;

        printf("INFO\t\t: Create toys for bin%d\n",whichBin);
        createFCToys2(whichBin,nToy);// profiled fl/afb value with fixed afb/fl value.
        scanProfile(whichBin);
    }else if (func == "FCScan"){
        int whichBin=-1;
        while(whichBin < 0 || whichBin > 12 || whichBin == 3 || whichBin == 5){
            printf("Which bin do you want to scan? ");
            scanf("%d",&whichBin);
        }

        printf("\t\t  Please run \"./macro/fitFCToys.py\" or \"macro/batchCreate.py\"to fit all toys in ./limit\n");

        printf("INFO\t\t: Press Enter to start scanning confidence interval for bin%d\n",whichBin);
        getchar();
        harvestFCFitResults(whichBin,nToy);
        getFCInterval2(whichBin,nToy);//likelihood ratio ordering

        // Coverage test
        struct stat fiBuff;
        if (stat(TString::Format("./limit/bestFit/setSummary_bin%d.root",whichBin),&fiBuff) != 0){
            harvestToyFitResult(whichBin,2000,"./limit/bestFit");// check "createBestFitToys"
            rename("./limit/bestFit/setSummary.root",TString::Format("./limit/bestFit/setSummary_bin%d.root",whichBin).Data());
        }
        //getFCResultCoverage2(whichBin);// new definition
    }else if (func == "createToys"){
        scaleFactor = 100;
        
        double expNCombBkg[nQ2Ranges] = { 0 , 0 , 0 , 0 , 65.0 , 0 , 0 , 0 , 0 , 55.6 , 136.0 , 0 , 0};
        int    expNSig[nQ2Ranges]     = { 0 , 0 , 0 , 0 , 22   , 0 , 0 , 0 , 0 , 39   , 21    , 0 , 0};
        for (int iBin = 0; iBin < nWorkBins; iBin++) {
            genToySignal(workBins[iBin],expNSig[workBins[iBin]]*scaleFactor);// (iBin, nEvts)
            genToyCombBkg(workBins[iBin],expNCombBkg[workBins[iBin]]*scaleFactor,"combBkgToy");
            genToyPeakBkgFromWspace(workBins[iBin], (int)scaleFactor);// (iBin, scaleFactor)
        }
    }else if (func == "createBestFitToys"){
        scaleFactor = 1;
        for (int iBin = 0; iBin < nWorkBins; iBin++) {
            createBestFitToys(workBins[iBin],2000);
        }
    }else if (func == "test"){
        ch->Add(infile.Data());
        if (ch == NULL) return 1;
        //angular3D_bin(4,"angular3D",scaleFactor);
        //createEdgeEffectToys(4,500);
        //harvestToyFitResult(4,500,"./limit/validation/edgeEffect/yieldsScaleE1/");
        //harvestToyFitResult(4,500,"./limit/validation/edgeEffect/yieldsScaleE2/");
        //harvestToyFitResult(4,500,"./limit/validation/edgeEffect/yieldsScaleE3/");
        //angular3D_prior3(10,"angular3D_prior",true);
        //angular3D_bin(10,"angular3D");
    }else{ 
        cerr << "No function available for: " << func.Data() << endl; 
    }
    //printf("%lld entries processed.\n",ch->GetEntries());

    return 0 ;
}
