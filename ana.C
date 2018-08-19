#include "tdrstyle.C"
void make1Dplot(vector<TH1F*> vhist, vector<const char*> vleg, bool logy, const char* axistitle, const char* plotname)
{
    setTDRStyle();
    gStyle->SetOptStat(0);
    //gStyle->SetPalette(kRainBow);
    
    for(int i=0;i<vhist.size();i++)
    {
        vhist[i]->Sumw2();
        vhist[i]->Scale(1./vhist[i]->Integral(-1,-1));
    }
    
    TCanvas *c = new TCanvas("c", "", 750, 750);
    c->cd();
    if(logy==true) c->SetLogy();
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.07);
    vhist[0]->GetYaxis()->SetTitleOffset(1.3);
    if(logy==true)
        vhist[0]->GetYaxis()->SetRangeUser(0.0001,(vhist[vhist.size()-1]->GetMaximum())*2.);
    else
        vhist[0]->GetYaxis()->SetRangeUser(0.0001,(vhist[vhist.size()-1]->GetMaximum())*1.3);
    vhist[0]->GetXaxis()->SetTitle(axistitle);
    vhist[0]->GetYaxis()->SetTitle("Events");
    vhist[0]->GetXaxis()->SetLabelSize(0.04);
    TLegend* leg = new TLegend(0.65,0.6,0.9,0.9);
    for(int i=0;i<vhist.size();i++)
    {
        vhist[i]->SetLineWidth(2);
        //vhist[i]->SetFillStyle(3004);
        vhist[i]->SetMarkerSize(1.2);
        vhist[i]->SetMarkerStyle(20+i);
        if(i==0)
            vhist[i]->Draw("PE PLC PMC");
        else
            vhist[i]->Draw("same PE PLC PMC");
        leg->AddEntry(vhist[i],vleg[i]);
    }
    leg->Draw("same");
    c->SaveAs(Form("plot/%s%s.png",plotname, (logy==true) ? "_logy" : ""));
    c->SaveAs(Form("plot/%s%s.pdf",plotname, (logy==true) ? "_logy" : ""));
}

void ana(float f, float epin, float ewall, int N)
{
    TFile *fin = new TFile(Form("minitree/miniTree_ep%.1f_ew%.1f_fre%.1f_dur%d.root",epin,ewall,f,N),"READ");
    TTree *t = (TTree*) fin->Get("outTree");
    
    float ep, ew, x, y, vx, vy, s_x, s_y, fdotv;
    t->SetBranchAddress("ep", &ep);
    t->SetBranchAddress("ew", &ew);
    t->SetBranchAddress("x", &x);
    t->SetBranchAddress("y", &y);
    t->SetBranchAddress("vx", &vx);
    t->SetBranchAddress("vy", &vy);
    t->SetBranchAddress("s_x", &s_x);
    t->SetBranchAddress("s_y", &s_y);
    t->SetBranchAddress("fdotv", &fdotv);
    
    int tmp_clock=0;
    float Jtau[10] = {0.0};
    vector<TH1F*> hM_Jtau; hM_Jtau.clear();
    const char* legtex[10] = {"#tau = 0.1s","#tau = 0.2s","#tau = 0.3s","#tau = 0.4s","#tau = 0.5s","#tau = 0.6s","#tau = 0.7s","#tau = 0.8s","#tau = 0.9s","#tau = 1.0s"};
    vector<const char*> v_leg (legtex, legtex + sizeof(legtex) / sizeof(const char*));
    for(int i=0;i<10;i++) hM_Jtau.push_back(new TH1F("hM_Jtau"," ", 115, -8000, 15000));
    
    for(Long64_t ev = 0; ev < t->GetEntriesFast(); ev++)
    {
        t->GetEntry(ev);
        tmp_clock++;
        
        for(int itau=0; itau<10; itau++)
        {
            Jtau[itau] = Jtau[itau]+fdotv;
            
            if(tmp_clock%100*(itau+1) == 0)
            {
                hM_Jtau[itau]->Fill(Jtau[itau]/(100*(itau+1)));
                Jtau[itau] = 0.;
            }
        }
    }
    //make1Dplot(vector<TH1F*> vhist, vector<const char*> vleg, const char* axistitle, const char* plotname)
    make1Dplot(hM_Jtau, v_leg, true, "J_{#tau} (erg)", Form("Jtau_ep%.1f_ew%.1f_fre%.1f_dur%d",epin,ewall,f,N));
    make1Dplot(hM_Jtau, v_leg, false, "J_{#tau} (erg)", Form("Jtau_ep%.1f_ew%.1f_fre%.1f_dur%d",epin,ewall,f,N));
}
