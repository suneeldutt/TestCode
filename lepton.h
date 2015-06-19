#include <iostream>
#include <cmath>

using namespace std;

#define pi 3.1415926539
#define twopi 6.2831853078
class lepton{
double px, py, pz, E;

public:
//constructor function for intialization
lepton(double, double, double, double);
double Pt();
double Eta();
double Phi();
double Theta();
double AngleTheta();
double InvMass();
double DiLeptonPt(lepton l1);
double DiLeptonEta(lepton l1);
double DiLeptonInvMass(lepton l1);
double DiLeptonAngle(lepton l1);
double DiLeptonDeltaEta(lepton l1);
double DiLeptonDeltaPhi(lepton l1);
double DiLeptonDeltaR(lepton l1);

};

double lepton::DiLeptonDeltaR(lepton l1){
return (sqrt( (DiLeptonDeltaEta(l1)*DiLeptonDeltaEta(l1)) + (DiLeptonDeltaPhi(l1)*DiLeptonDeltaPhi(l1)) ));
}

double lepton::DiLeptonDeltaEta(lepton l1){
return (Eta() - l1.Eta());
}

double lepton::DiLeptonDeltaPhi(lepton l1){
double dphi= Phi() - l1.Phi();
if(dphi >= pi) dphi -= twopi;
if(dphi< -pi) dphi += twopi;
return dphi;
}
double lepton::Theta(){
double p = sqrt(px*px+py*py+pz*pz);
return (pz/p);
}
double lepton::AngleTheta(){
double p = sqrt(px*px+py*py+pz*pz);
return (acos(pz/p)*(180./pi));
}
double lepton::Phi(){
double p = sqrt(px*px+py*py+pz*pz);
double theta = acos(pz/p);
double phi = acos(px/(p*sin(theta)));
if(phi >= pi) phi -= twopi;
if(phi< -pi) phi += twopi;
return phi;
}
//definition of intializing constructor
lepton::lepton(double a, double b, double c, double d){
px=a;
py=b;
pz=c;
E=d;
}
double lepton::Pt(){
return (sqrt(px*px+py*py));
}
double lepton::Eta(){
double p = sqrt(px*px+py*py+pz*pz);
double eta = log((p+pz)/(p-pz));
return (eta/2);
}
double lepton::InvMass(){
double p = sqrt(px*px+py*py+pz*pz);
return (sqrt(E*E-p*p));
}
double lepton::DiLeptonPt(lepton l1){
double SumPx = px + l1.px;
double SumPy = py + l1.py;
//double SumPz = pz + l1.pz;
return (sqrt(SumPx*SumPx+SumPy*SumPy));
}
double lepton::DiLeptonEta(lepton l1){
double SumPx = px + l1.px;
double SumPy = py + l1.py;
double SumPz = pz + l1.pz;
double SumP = sqrt(SumPx*SumPx+SumPy*SumPy+SumPz*SumPz);
double eta = log((SumP+SumPz)/(SumP-SumPz));
return (eta/2);
}

double lepton::DiLeptonInvMass(lepton l1){
double SumPx = px + l1.px;
double SumPy = py + l1.py;
double SumPz = pz + l1.pz;
double SumE = pz + l1.E;
double SumP = sqrt(SumPx*SumPx+SumPy*SumPy+SumPz*SumPz);
return (sqrt(SumE*SumE-SumP*SumP));
}
