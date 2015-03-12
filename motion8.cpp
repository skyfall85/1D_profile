/*
 Ez a kod az egy komponensu rendszer leirasara lett atalakitva: 
  kikerul belole az f_s(c,T), illetve csak egy komponensnek lesz benne a latens hoje
 
2015.02.12. Csütörtök
- most kezdek el ujra ezzel a koddal foglalkozni
- meg kell erteni, hogy pontosan mi mit is csiinal, es ahol lehet, ott egyszerusiteni kellene
 
 
*/

#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <iomanip>
// #include "functions.h"
#include "find_and_main.h"
#include "free_energy.h"


/*----------------------------------------------------------------
       -------   VARIABLE DECLARATION STARTS     -------
----------------------------------------------------------------*/

double K0, V_m, R;
double eps2,W,dt;

double d,dt0,T,d0,gamma_Ni,xi,c_s,c_l,c0;
double dless_M_phi,Q_Cu,Q_Ni,T_Cu,T_Ni,dless_HT_Plapp,dless_HT,A_phase_noise,k,s_0;
double Ds_rat,dless_Mc,E_t;
double dless_M_theta_s_Plapp,dless_M_theta_l_Plapp,dless_M_theta_s,dless_M_theta_l;
  
double  theta_target;

int num_of_lines;                  // megadja, hogy hany sort irtam ki   
/*----------------------------------------------------------------
        -------   VARIABLE DECLARATION FINISHES   -------
----------------------------------------------------------------*/

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
 
  fscanf(f1,"%le %le %le %le %le %le %le %le %le",&d,&dt0,&T,&d0,&gamma_Ni,&xi,&c_s,&c_l,&c0);
  fscanf(f2,"%le %le %le %le %le %le %le %le %le %le %le",&dless_M_phi,&W,&Q_Cu,&Q_Ni,&T_Cu,&T_Ni,&dless_HT_Plapp,&dless_HT,&A_phase_noise,&k,&s_0);
  fscanf(f3,"%le %le %le",&Ds_rat,&dless_Mc,&E_t);
  fscanf(f4,"%le %le %le %le",&dless_M_theta_s_Plapp,&dless_M_theta_l_Plapp,&dless_M_theta_s,&dless_M_theta_l);
  
// Itt adok meg egyedul sajat valtozokat, a tobbit beolvasom:
 
  calculate_parameters( V_m, R,eps2, K0, gamma_Ni,d0, Q_Cu, Q_Ni, E_t, T, T_Cu,  T_Ni,  W, dless_HT_Plapp);
 
  int H_enough=0; 
  dt=0.00625/128.0*1e-4;

  while(H_enough!=1)
  {
    while(theta_target>=0.02)
    {
    num_of_lines=0;

    H_enough=Find_K0( K0,  theta_target,  dless_HT_Plapp,  W,  Q_Ni, T,  T_Ni);

    cout<<"dless_HT_Plapp ciklusban:\t"<<dless_HT_Plapp<<endl;
 
    Main_cicle(W,Q_Ni,T,T_Ni,K0,dless_HT_Plapp,num_of_lines,dt,theta_target);

    ofstream out_free_e("free_energy.dat",ios::app);
    out_free_e<<theta_target<<"\t"<<free_energy( gamma_Ni, d0, T, T_Ni, Q_Ni, dless_HT_Plapp, xi,num_of_lines)<<endl;
    out_free_e.close();
 
    cout<<"Theta-target:\t"<<theta_target<<"\t"<<"Free energy:\t"<<free_energy( gamma_Ni, d0, T, T_Ni, Q_Ni, dless_HT_Plapp, xi,num_of_lines)<<endl;
    theta_target>0.2 ? theta_target-=0.05 : theta_target-=0.005;
  }
  }
 
 return 0;
}
