#include <iomanip>

int factor[] = {0,6,18,36,60,90,126,168,216,270,330,396,463,499};
int shift[] = {1,1,1,1,1,1,1,1,1,1,1,1,0,0};
int teltype = 0;  //what is written in the output
string filename = "VTSCamera.txt";
void VTSCameraDraw()
{
  static const int maxpix=2000;
  double r=31.4; //pixel spacing in mm
  double pi=TMath::Pi();
  double deg2rad=pi/180.; //multiply degrees by this for radians
  double xpix[maxpix];
  double ypix[maxpix];
  double plate_radius=760./2.;
  double ext_radius=plate_radius+260./2.;
  //double ext_radius=(760.+180.)/2.;

  for (int i=0;i<maxpix;i++) 
    {
      xpix[i]=0;
      ypix[i]=0;
    } 

  int nrings=20;
  int pixcount=0;
  
  //central pixel
  xpix[pixcount]=0.;
  ypix[pixcount]=0.;
  pixcount+=1;	  

  for (int i=0;i<14;i++) {
    double angle=120;
    double x=i*r;
    double y=0;
    for (int j=0;j<6;j++)
      {
	double vx=r*cos(angle*deg2rad);
	double vy=r*sin(angle*deg2rad);
	for (int k=0;k<i;k++) 
	  {
	    x+=vx;
	    y+=vy;
	    double dist=sqrt(x*x+y*y);
	    double delta=(r/2.);
	    if (dist>plate_radius-delta && dist < plate_radius+delta+15) continue;
	    if (dist>ext_radius-delta) continue;
            //double xrot = x*cos((-60*i)*deg2rad)-y*sin((-60*i)*deg2rad);
            //double yrot = x*sin((-60*i)*deg2rad)+y*cos((-60*i)*deg2rad);
            //xpix[(pixcount)]=x;
            int npix = pixcount>=factor[i] ? pixcount%factor[i]+factor[i-1] : pixcount;
            npix+=shift[i];
	    xpix[npix]=x;
	    ypix[npix]=-1*y;
	    //ypix[(pixcount)]=-1*y;
	    pixcount+=1;	  
	  }
	angle+=60.;
      }
  }

  cout << "Inner plate radius=" << plate_radius << endl;
  cout << "Outer plate radius=" << ext_radius << endl;
  cout << "Number of pixels=" << pixcount << endl;

  gROOT->SetStyle("Plain");

  TH2F *h2=new TH2F("","",1,-550,550,1,-550,550);
  h2->SetStats(0);
  h2->GetXaxis()->SetTitle("X (mm)");
  h2->GetYaxis()->SetTitle("Y (mm)");
  h2->Draw();
  TGraph *g=new TGraph(pixcount,xpix,ypix);
  g->SetTitle("");
  g->GetXaxis()->SetTitle("X (mm)");
  g->GetYaxis()->SetTitle("Y (mm)");
  g->SetMarkerStyle(4);
  g->SetMarkerSize(3);
  g->SetMarkerColor(2);
  g->Draw("P");


ofstream outputf(filename.c_str());
for (int i=0;i<pixcount;i++)
{
  char pix_label[100];
  sprintf(pix_label,"%d",i+1);
  TLatex *l=new TLatex(xpix[i],ypix[i], pix_label);
  l->SetTextSize(0.013);
  l->SetTextAlign(22);
  l->Draw();
  cout<<"* PMPIX "<<teltype<<" 0    "; 
  cout << fixed;
  cout << i << "\t" << setw(10) << setprecision(2) << xpix[i] << "\t" << setw(10) << setprecision(2) << ypix[i] <<" 18.1 0.0"<<endl;
  outputf<<"* PMPIX "<<teltype<<" 0    "; 
  outputf << fixed;
  outputf << i << "\t" << setw(10) << setprecision(2) << xpix[i] << "\t" << setw(10) << setprecision(2) << ypix[i] <<" 18.1 0.0"<<endl;
}

  TEllipse *plate=new TEllipse(0,0,plate_radius);
  plate->SetFillStyle(0);
  plate->Draw("same");

  TEllipse *ext=new TEllipse(0,0,ext_radius);
  ext->SetLineColor(3);
  ext->SetFillStyle(0);
  ext->Draw("same");

  double b1x1=-14.25;
  double b1y1=380;
  double b1x2=14.25;
  double b1y2=480;
  double b1theta=0;

  TPave *bracket1=new TPave(b1x1,b1y1,b1x2,b1y2,1);
  //bracket1->Draw();

  const int npts=7;
  double xpts[npts];
  double ypts[npts];
  xpts[0]=-15;
  xpts[5]=+15;
  xpts[1]=xpts[0];
  xpts[4]=xpts[5];
  xpts[2]=xpts[0]-100;
  //sqrt(100*100+100*100);
  xpts[3]=xpts[5]+100;
  //sqrt(100*100+100*100);
  xpts[6]=xpts[0];
  ypts[0]=380;
  ypts[1]=430;
  ypts[2]=520;
  ypts[3]=520;
  ypts[4]=430;
  ypts[5]=380;
  ypts[6]=ypts[0];

  TGraph *g1=new TGraph(npts,xpts,ypts);
  g1->SetFillStyle(3001);
  g1->Draw("F");
  
  for (int i=0;i<npts;i++)
    {
      double cth=cos(120.*deg2rad);
      double sth=sin(120.*deg2rad);
      double x=xpts[i];
      double y=ypts[i];
      xpts[i]=x*cth-y*sth;
      ypts[i]=x*sth+y*cth;
    }

  TGraph *g2=new TGraph(npts,xpts,ypts);
  g2->SetFillStyle(3001);
  g2->Draw("F");

  for (int i=0;i<npts;i++)
    {
      double cth=cos(120.*deg2rad);
      double sth=sin(120.*deg2rad);
      double x=xpts[i];
      double y=ypts[i];
      xpts[i]=x*cth-y*sth;
      ypts[i]=x*sth+y*cth;
    }
  TGraph *g3=new TGraph(npts,xpts,ypts);
  g3->SetFillStyle(3001);
  g3->Draw("F");

  //TPave *bracket2=new TPave(b2x1,b2y1,b2x2,b2y2,1);
  //bracket2->Draw();*/

  // hexapods extend above and below the bracket - need 10cm to be safe. Start from 2/3rds of the way out along the bracket (~6cm from the plate).


}
