/*
 Ez a kod az egy komponensu rendszer leirasara lett atalakitva: 
  kikerul belole az f_s(c,T), illetve csak egy komponensnek lesz benne a latens hoje
 
*/


#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <iomanip>

using namespace std;

double qt, qni;
double d0, K0, V_m, R;
double eps2,W;
int K0i=96;
double HT0;

double d,dt,T,gamma_Ni,xi,c_s,c_l,c0;
double dless_M_phi,Q_Cu,Q_Ni,T_Cu,T_Ni,dless_HT_Plapp,dless_HT,A_phase_noise,k,s_0;
double Ds_rat,dless_Mc,E_t;
double dless_M_theta_s_Plapp,dless_M_theta_l_Plapp,dless_M_theta_s,dless_M_theta_l;
  
double x,v,v1,/*dt,*/dx,t,dv,theta,t_r,conc,gradV; 
double theta_mistake, theta_target;
int icount, felt;

double HT1, dHT1;

int num_of_lines;                  // megadja, hogy hany sort irtam ki   
double g(double x)
{
  return 0.25*x*x*(1.0-x)*(1.0-x);
}

double gv(double x)
{
  return 0.5*x*(1.0-x)*(1.0-2.0*x);
}

double gh(double x)
{
  return (7.0*x*x*x-6.0*x*x*x*x)/((1.0-x)*(1.0-x)*(1.0-x));
}

double ghv(double x)
{
  return (21.0*x*x-24.0*x*x*x)/((1.0-x)*(1.0-x)*(1.0-x))+3.0*(7.0*x*x*x-6.0*x*x*x*x)/((1.0-x)*(1.0-x)*(1.0-x)*(1.0-x));
}

double p(double x)
{
  return x*x*x*(10.0-15.0*x+6.0*x*x);
}

double pv(double x) 
{
  30.0*x*x-60.0*x*x*x+30.0*x*x*x*x;
}

double Veff(double x)
{
  return -W*g(x)+p(x)*Q_Ni*(1.0-T/T_Ni)+K0*K0*dless_HT_Plapp/gh(x);  
}

double Veff_inside(double x)
{
  return -W*g(x)+p(x)*Q_Ni*(1.0-T/T_Ni)+K0*K0*HT1/gh(x);  
}


double min(double a, double b)
{
  double retval;
  if (a<b) retval=a;
  else retval=b;
  return retval;
}

void parameters()
{
  
  cout<<"W:\t"<<W<<endl;
  cout<<"c0:\t"<<c0<<endl;
  cout<<"Q_Cu:\t"<<Q_Cu<<endl;
  cout<<"Q_Ni:\t"<<Q_Ni<<endl;
  cout<<"T_Ni:\t"<<T_Ni<<endl;
  cout<<"T_Cu:\t"<<T_Cu<<endl;
  cout<<"T:\t"<<T<<endl;
  cout<<"E_t:\t"<<E_t<<endl;
  cout<<"H:\t"<<dless_HT_Plapp<<endl;
  cout<<"qt:\t"<<qt<<endl;
  cout<<"qni:\t"<<qni<<endl;
  cout<<"K0:\t"<<K0<<endl;
  

}

void cicle(double, double);

double free_energy();

void Find_H();

double Find_K0();

void Main_cicle();

/****************************************************************************************************
                                             MAIN STARTS
*****************************************************************************************************/

int main()
{

 theta_target=0.5;                                  //!!! itt kell allitani a dtheta-t
  

 FILE *f1, *f2, *f3, *f4;
 f1=fopen("other_parameters","r");
 f2=fopen("phase_field_parameters","r");
 f3=fopen("concentration_parameters","r");
 f4=fopen("orientation_parameters","r");
 
 fscanf(f1,"%le %le %le %le %le %le %le %le %le",&d,&dt,&T,&d0,&gamma_Ni,&xi,&c_s,&c_l,&c0);
 fscanf(f2,"%le %le %le %le %le %le %le %le %le %le %le",&dless_M_phi,&W,&Q_Cu,&Q_Ni,&T_Cu,&T_Ni,&dless_HT_Plapp,&dless_HT,&A_phase_noise,&k,&s_0);
 fscanf(f3,"%le %le %le",&Ds_rat,&dless_Mc,&E_t);
 fscanf(f4,"%le %le %le %le",&dless_M_theta_s_Plapp,&dless_M_theta_l_Plapp,&dless_M_theta_s,&dless_M_theta_l);
  
// Itt adok meg egyedul sajat valtozokat, a tobbit beolvasom:
 
 V_m=7.4e-6;
 R=8.3144087;

 eps2=6.0*sqrt(2.0)*gamma_Ni*d0/1728.0;
 qt=(Q_Cu/E_t)*(1.0-T/T_Cu)-(Q_Ni/E_t)*(1.0-T/T_Ni);
 qni=Q_Ni*(1.0-T/T_Ni);
 K0=sqrt(0.96*0.25*W/dless_HT_Plapp);

 HT0=dless_HT_Plapp;
 
 cout<<"Egyutthatok:\t"<<W<<"\t"<<qt<<"\t"<<E_t<<endl;
//  cout<<"Kezdeti H:\t"<<dless_HT_Plapp<<endl;
 
double H_enough=0; 
//      Find_H();                                                     // kiszamolja a H0 értekhez tartozo dtheta-t, ami a theta-kent ter vissza
//      dless_HT_Plapp=dless_HT_Plapp*(theta/0.01)*(theta/0.01); 
//       cout<<"H szorzofaktora:\t"<<theta/0.01<<endl;
 int H_or_K=1;
 while(H_enough!=1)
 {
  num_of_lines=0;
  switch(H_or_K)
  {
    case 0:
      Find_H();                                                     // kiszamolja a H0 értekhez tartozo dtheta-t, ami a theta-kent ter vissza
      dless_HT_Plapp=dless_HT_Plapp*(theta/0.01)*(theta/0.01);
      cout<<"H szorzofaktora:\t"<<theta/0.01<<endl;
      K0=sqrt(0.96*0.25*W/dless_HT_Plapp);
    break;
  
    case 1:
      H_enough=Find_K0();
    break;
  }
 
  cout<<"dless_HT_Plapp ciklusban:\t"<<dless_HT_Plapp<<endl;
 
 Main_cicle();

 ofstream out_free_e("free_energy.dat",ios::app);
 out_free_e<<theta_target<<"\t"<<free_energy()<<endl;
 out_free_e.close();
 
 cout<<"Theta-target:\t"<<theta_target<<"\t"<<"Free energy:\t"<<free_energy()<<endl;
 theta_target>0.2 ? theta_target-=0.05 : theta_target-=0.005;
 }
 
 return 0;
}
/**************************************************************************************************/
/*---------------------------------------------MAIN ENDS------------------------------------------*/
/**************************************************************************************************/

void cicle(double K0, double H)
{
   x=0.9999;  
   v=0.0;  
   v1=0.0;  
   dx=1.0e-9;   
   theta=0.0;  
   t_r=0;
  
   dt=0.00625/128.0*1e-2;
 
  
  while (felt!=3)
  {
   if (felt==0 && x<0.97) {felt=1;}    // eloszor lesz kisebb, mint 0.9
   if (felt==1 && x>0.97) {felt=2;}    // eloszor ter vissza 0.9 fole
   if (x>0.999 && felt==2) felt=3; // eloszor ter vissza 0.9999 fole
   icount++;
   gradV=(Veff_inside(x+dx)-Veff_inside(x-dx))/(2.0*dx);
   dv=-gradV*dt;

   v+=dv;
   x+=(v)*dt;
	v1=v;
   theta+=K0/gh(x)*dt;
 }
  
  
}

double free_energy()
{

  double x_array[num_of_lines];
  double pf_array[num_of_lines];
  double of_array[num_of_lines];
  double grad_pf[num_of_lines];
  double grad_of[num_of_lines];
  double grad_pf_tag[num_of_lines];
  double grad_of_tag[num_of_lines];
  double thermal_tag[num_of_lines];
  
  double epsilon2=6.0*sqrt(2.0)*gamma_Ni*d0/T_Ni;
  double W=6.0*sqrt(2.0)*gamma_Ni/(d0*T_Ni);
  double H=dless_HT_Plapp*epsilon2;
  
  ifstream in_x("x_values.dat");
  ifstream in_pf("pf_values.dat");
  ifstream in_of("of_values.dat");
  
  for(int ibe=0;ibe<num_of_lines;ibe++)
  {
    in_x>>x_array[ibe];
    in_pf>>pf_array[ibe];
    in_of>>of_array[ibe];
  }
  in_x.close();
  in_pf.close();
  in_of.close();
  double dx_free=(x_array[1]-x_array[0])*xi;
  
  for (int i_grad=1;i_grad<(num_of_lines-1);i_grad++)
  {
    grad_pf[i_grad]=(pf_array[i_grad+1]-pf_array[i_grad-1])/(2.0*dx_free);
    grad_of[i_grad]=(of_array[i_grad+1]-of_array[i_grad-1])/(2.0*dx_free); 
  }
  grad_pf[0]=grad_pf[1];  grad_pf[num_of_lines-1]=grad_pf[num_of_lines-2];
  grad_of[0]=grad_of[1];grad_of[num_of_lines-1]=grad_of[num_of_lines-2];
  
  for (int i1=0;i1<num_of_lines;i1++)
  {
    grad_pf_tag[i1]=epsilon2*T/2.0*grad_pf[i1]*grad_pf[i1];
    grad_of_tag[i1]=H*T*gh(pf_array[i1])*grad_of[i1]*grad_of[i1];
    thermal_tag[i1]=p(pf_array[i1])*Q_Ni*(1.0-T/T_Ni);
  }
  
double Free_Energy=0;

// cout<<epsilon2<<"\t"<<T<<"\t"<<W<<"\t"<<H*T<<"\t"<<dless_HT_Plapp<<"\t"<<xi<<endl;
// for (int io=0;io<2;io++) cout<<pf_array[io]<<"\t"<<of_array[io]<<"\t"<<grad_pf_tag[io]<<"\t"<<grad_of_tag[io]<<endl;
  
  for (int i_e=0;i_e<num_of_lines;i_e++)
  {
    Free_Energy+=(grad_pf_tag[i_e]+W*T*g(pf_array[i_e])+thermal_tag[i_e]+grad_of_tag[i_e])*dx_free;
  }
  
  return Free_Energy;
}


void Find_H()
{
   HT1=dless_HT_Plapp;
   K0=sqrt(0.96*0.25*W/HT1);

   icount=0;
   felt=0;
   
   x=0.9999;  
   v=0.0;  
   v1=0.0;  
   dx=1.0e-9;   
   theta=0.0;  
   t_r=0;
  
   dt=0.00625/128.0*1e-2;
 
  
  while (felt!=3)
  {
   if (felt==0 && x<0.97) {felt=1;}    // eloszor lesz kisebb, mint 0.9
   if (felt==1 && x>0.97) {felt=2;}    // eloszor ter vissza 0.9 fole
   if (x>0.999 && felt==2) felt=3; // eloszor ter vissza 0.9999 fole
   icount++;
   gradV=(Veff_inside(x+dx)-Veff_inside(x-dx))/(2.0*dx);
   dv=-gradV*dt;

   v+=dv;
   x+=(v)*dt;
	v1=v;
   theta+=K0/gh(x)*dt;
	
  if (icount%1000000==0) cout<<setprecision(9)<<x<<"\t"<<gradV<<"\t"<<v<<endl;
 }

}

double Find_K0()
{
/*
 Egy iteraciokon keresztul meghatarozza azt a K0 erteket, amit hasznalva az adott es rogzitett HT mellett dtheta=theta_target
 nem lesz. 
 A fuggvenynek nincsen visszateresi erteke, de K0-t megvaltoztatja, amit majd a Main_cicle() at fog venni.
*/  
  HT1=dless_HT_Plapp;
  K0=sqrt(0.96*0.25*W/HT1);
  
  double K1=K0*0.5;
  double theta_mist=0.5;
  int ciklus=1;
  int first=0;
  
  
  
  while(theta_mist>0.001 && K1>1.0e-9)
  {
   if (theta_target==0.3)
   {
       int theta_name=floor(100*theta_target);
  char filename_2[100];
  char filename_3[100];
  sprintf(filename_2,"potential_%d.txt",theta_name); 
  sprintf(filename_3,"potentials_%d.txt",theta_name); 
  ofstream out2(filename_2);
  ofstream out3(filename_3);
  
    for (int i=1;i<100000;i++)
  {
    double i2=i;
    double V0=Veff_inside(1.0);
    i2*=0.00001;
    out2<<i2<<"\t"<<Veff_inside(i2)-V0<<endl;  
    out3<<i2<<"\t"<<-W*g(i2)<<"\t"<<(p(i2)-p(1.0))*Q_Ni*(1.0-T/T_Ni)<<"\t"<<K0*K0*dless_HT_Plapp/gh(i2)<<endl;
 }
  }
    
      felt=0;
   
   x=0.9999;  
   v=0.0;  
   dx=1.0e-9;   
   theta=0.0;  
   dt=0.00625/128.0*1e-2;
 
  
   while (felt!=3)
   {
    if (felt==0 && x<0.97) {felt=1;}    // eloszor lesz kisebb, mint 0.9
    if (felt==1 && x>0.97) {felt=2;}    // eloszor ter vissza 0.9 fole
    if (x>0.999 && felt==2) felt=3; // eloszor ter vissza 0.9999 fole
    icount++;

    gradV=(Veff_inside(x+dx)-Veff_inside(x-dx))/(2.0*dx);
    dv=-gradV*dt;
    v+=dv;
    x+=(v)*dt;
    theta+=K0/gh(x)*dt;
//     if (icount%10000==0) cout<<"x: "<<x<<endl;
   }
   cout<<"theta:\t"<<theta<<"K0:\t"<<K0<<"K1:\t"<<K1<<endl;
  theta_mist=fabs(theta_target-theta);
//  if (ciklus==1 && theta_target<theta){cout<<theta_target<<"\t szog mar nem elerheto ekkora H mellett"<<endl; return 1;}
  while(first==0)
  {
   felt=0;
   
   x=0.9999;  
   v=0.0;  
   dx=1.0e-9;   
   theta=0.0;  
   dt=0.00625/128.0*1e-2;
 
  
   while (felt!=3)
   {
    if (felt==0 && x<0.97) {felt=1;}    // eloszor lesz kisebb, mint 0.9
    if (felt==1 && x>0.97) {felt=2;}    // eloszor ter vissza 0.9 fole
    if (x>0.999 && felt==2) felt=3; // eloszor ter vissza 0.9999 fole
    icount++;

    gradV=(Veff_inside(x+dx)-Veff_inside(x-dx))/(2.0*dx);
    dv=-gradV*dt;
    v+=dv;
    x+=(v)*dt;
    theta+=K0/gh(x)*dt;
    if (icount==floor(1e10)) return 1;
   }
   cout<<"theta:\t"<<theta<<"K0:\t"<<K0<<"K1:\t"<<K1<<endl;
  theta_mist=fabs(theta_target-theta);
    
    
    if (theta_target<theta)
    {
      K1*=2.0;
      K0+=K1;
      cout<<"meg nagyobb"<<endl;
    }
    if (theta_target>theta) {first=1;cout<<"Kisebb lett!"<<endl;}
  }
  
  if (theta_target>theta && theta_mist>0.001 && first==1){K0-=K1;K1*=0.5;}
  if (theta_target<theta && theta_mist>0.001 && first==1){K0+=K1;K1*=0.5;}
 
  ciklus++; 
  if (ciklus>100 || theta_target<=0) return 1;
  }

  return 0;
}


void Main_cicle()
{
  char filename_1[100];
  char filename_2[100];
  char filename_3[100];
  int theta_name=floor(100*theta_target);
  
  sprintf(filename_1,"trajectory_%d.txt",theta_name);
  sprintf(filename_2,"potential_%d.txt",theta_name); 
  sprintf(filename_3,"potentials_%d.txt",theta_name); 
  
  ofstream out(filename_1);
  ofstream out2(filename_2);
  ofstream out3(filename_3);
  ofstream out_x("x_values.dat");
  ofstream out_pf("pf_values.dat");
  ofstream out_of("of_values.dat");
  
  for (int i=1;i<100000;i++)
  {
    double i2=i;
    double V0=Veff(1.0);
    i2*=0.00001;
    out2<<i2<<"\t"<<Veff(i2)-V0<<endl;  
    out3<<i2<<"\t"<<-W*g(i2)<<"\t"<<(p(i2)-p(1.0))*Q_Ni*(1.0-T/T_Ni)<<"\t"<<K0*K0*dless_HT_Plapp/gh(i2)<<endl;
 }
 
 

/*
A szamolt parameterek inicializalasa:
Jeloles a programban   Fizikai jelentes/ jelentes
        x ------------------- phi
        v ------------------- d(phi)/dx
       dx ------------------- gradiens szamolasnal dx
        t ------------------- x (hely koordinata)
       dv ------------------- v+=dv
    theta ------------------- orientacio   
       dt ------------------- dx
*/ 
  x=0.9999;
  v=0.0;
  dx=1.0e-9;
  t=0;
  dv;
  theta=0.0;
  dt=0.00625/128.0*1e-4;
  
  icount=0;
  felt=0;
  
  cout<<"Az orientacio egyutthatoja: "<<K0*K0*HT1<<endl;
  
  out<<t<<"\t"<<x<<"\t"<<theta<<endl;
  out_x<<t<<endl;
  out_pf<<x<<endl;
  out_of<<theta<<endl;
  num_of_lines++;
 while (felt!=3)
 {
   if (felt==0 && x<0.97) felt=1;    // eloszor lesz kisebb, mint 0.9
   if (felt==1 && x>0.97) felt=2;    // eloszor ter vissza 0.9 fole
   if (x>0.999 && felt==2) felt=3;   // eloszor ter vissza 0.9999 fole
   icount++;
   gradV=(Veff(x+dx)-Veff(x-dx))/(2.0*dx);
   dv=-gradV*dt;

   v+=dv;
   x+=(v)*dt;
   theta+=K0/gh(x)*dt;
   t+=dt;

   if (icount%1000==0) 
   {
     out<<setprecision(9)<<t<<"\t"<<x<<"\t"<<theta<<endl;
     out_x<<setprecision(9)<<t<<endl;
     out_pf<<setprecision(9)<<x<<endl;
     out_of<<setprecision(9)<<theta<<endl;
   
     num_of_lines++;
  }
   
 }
 
 cout<<"hossz:\t"<<num_of_lines<<endl;
  out.close();
  out2.close();
  out3.close();
  out_x.close();
  out_pf.close();
  out_of.close();
  
}





