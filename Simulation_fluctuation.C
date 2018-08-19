#include <vector>
#include <iostream>
#include <fstream>
#include <Riostream.h>
#include <stdlib.h>
#include <algorithm>
#include <ctime>
#include <TH1D.h>
#include <TH1F.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphAsymmErrors.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLatex.h>
#include <TPad.h>
#include <TFrame.h>
#include <TTree.h>
#include <TSystem.h>
using namespace std;

void Simulation_fluctuation(float f, float ep, float ew, int N)
{ 
    float pi = M_PI ;
    float m = 2.21;
    float s_x , s_y , x=0. , y=0. ; /* s : distance that the sphere travels in dt
                                     x : x coordinate of the sphere
                                     y : y coordinate of the sphere */
    float xtemp , ytemp ;
    float vx=0. , vy=0. ; // v:velocity
    unsigned int seed = time(NULL);
    TRandom *random = new TRandom3(seed);
    float phi_0 = (float) (random->Integer(36000))/100. ; // phi_0 : phase angle of the disk at beginning
    float alpha = pi/18. ; // alpha : tilt angle of the disk
    float omega ; // omega : angular speed of the plate
    float dt = 1./1000. , tsum = 0. , T ;
    float g = 980. ;
    double d[16][3] ; /* d[][0] : x coordinate of pin
                       d[][1] : y coordinate of pin
                       d[][2] : the distance of the pin and sphere */
    float fdotv=0.;
    
    FILE *fptp;
    fptp=fopen("pin_coordinate.txt","r");
    for(int i=0;i<16;i++)
    {
        for(int j=0;j<3;j++)
        {
            fscanf(fptp,"%lf",&d[i][j]);
        }
    }
    fclose(fptp);
    
    cout << "frequency(Hz) : " << f << endl;
    cout << "coefficient of the pin: " << ep << endl;
    cout << "coefficient of the wall: " << ew << endl;
    omega = 2. * pi * f;
    
    char fileout1[60];
    //int N;
    cout << "Simulation duration (in unit of minute): " << N << endl;
    //cin >> N ;
    
    TFile *fo = new TFile(Form("minitree/miniTree_ep%.1f_ew%.1f_fre%.1f_dur%d.root",ep,ew,f,N),"RECREATE");
    //if (!fo || fo->IsZombie())  FATAL("TFile::Open() failed");
    fo->cd();
    TTree *outTree = new TTree("outTree", " ");
    Long64_t event = 0;
    outTree->Branch("event",     &event,     "event/L");
    outTree->Branch("tsum",     &tsum,     "tsum/F");
    outTree->Branch("phi_0", &phi_0, "phi_0/F");
    outTree->Branch("alpha", &alpha, "alpha/F");
    outTree->Branch("f", &f, "f/F");
    outTree->Branch("omega", &omega, "omega/F");
    outTree->Branch("alpha", &alpha, "alpha/F");
    outTree->Branch("g", &g, "g/F");
    outTree->Branch("ep", &ep, "ep/F");
    outTree->Branch("ew", &ew, "ew/F");
    outTree->Branch("x", &x, "x/F");
    outTree->Branch("y", &y, "y/F");
    outTree->Branch("vx", &vx, "vx/F");
    outTree->Branch("vy", &vy, "vy/F");
    outTree->Branch("s_x", &s_x, "s_x/F");
    outTree->Branch("s_y", &s_y, "s_y/F");
    outTree->Branch("fdotv", &fdotv, "fdotv/F");
    
    TH2F *hM_xy = new TH2F("hM_xy"," ",10000,-10,10,10000,-10,10);
    
    for(int i=0;i<N;i++)
    {
        tsum=60.*i;
        
        while(tsum<=60.*(i+1))
        {
            s_x = vx*dt + 0.5*sin(alpha)*cos(phi_0*(pi/180.) + omega*tsum) ;
            s_y = vy*dt + 0.5*dt*dt*g*sin(alpha)*sin(phi_0*(pi/180.) + omega*tsum) ;
            xtemp = x ;
            ytemp = y ;
            x = x + s_x ;
            y = y + s_y ;
            vx = vx + dt*g*sin(alpha)*cos(phi_0*(pi/180.) + omega*tsum) ;
            vy = vy + dt*g*sin(alpha)*sin(phi_0*(pi/180.) + omega*tsum) ;
            float fx = m*g*sin(alpha)*cos(phi_0*(pi/180.) + omega*tsum);
            float fy = m*g*sin(alpha)*sin(phi_0*(pi/180.) + omega*tsum);
            
            for(int a=0;a<=15;a++)
            {
                d[a][2] = sqrt(pow(x-d[a][0],2)+pow(y-d[a][1],2));
                
                if(d[a][2]<=1.04)
                {
                    x=xtemp;
                    y=ytemp;
                    
                    d[a][2] = sqrt(pow(x-d[a][0],2)+pow(y-d[a][1],2));
                    
                    vx = (-ep-1)*((vx*(x-d[a][0])+vy*(y-d[a][1]))/(d[a][2]*d[a][2]))*(x-d[a][0]) + vx;
                    vy = (-ep-1)*((vx*(x-d[a][0])+vy*(y-d[a][1]))/(d[a][2]*d[a][2]))*(y-d[a][1]) + vy;
                    
                    s_x = vx*dt;
                    s_y = vy*dt;
                    x = x + s_x ;
                    y = y + s_y ;
                    
                    break;
                }
            }
            
            if(x<=-7.710000||x>=7.710000)
            {
                x=xtemp;
                y=ytemp;
                vx=(-1)*ew*vx;
                vy=vy;
                
                s_x = vx*dt;
                s_y = vy*dt;
                x = x + s_x ;
                y = y + s_y ;
            }
            if(y<=-7.710000||y>=7.710000)
            {
                x=xtemp;
                y=ytemp;
                vx=vx;
                vy=(-1)*ew*vy;
                
                s_x = vx*dt;
                s_y = vy*dt;
                x = x + s_x ;
                y = y + s_y ;
            }
            
            fdotv = fx*vx + fy*vy;
            
            tsum=tsum+dt;
            
            hM_xy->Fill(x,y);
            
            outTree->Fill();
            event++;
        }
    }
    
    outTree->Write("", TObject::kOverwrite);
    delete outTree;
    
    TCanvas *c = new TCanvas("c"," ",800,800);
    c->cd();
    c->SetLogz();
    hM_xy->GetXaxis()->SetTitle("x(cm)");
    hM_xy->GetYaxis()->SetTitle("y(cm)");
    hM_xy->Draw("colz");
    c->SaveAs("xy_distribution.png");
} 
