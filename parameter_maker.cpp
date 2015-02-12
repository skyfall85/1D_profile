#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
/*------------------Parameterek--------------------*/


using namespace std;

int main(int argc,char *argv[])
{
double H_coeff=atoi(argv[1]);

double d0_coeff=atoi(argv[2]);
/*-----Valtoztathato parameterek-------------------*/
double T=1720;
double  c_s=0.3977252936799706;
double  c_l=0.4649878481647202;
double c0=0;                  // az egy komponensu rendszerben csak az egyik komponens van jelen
int k=4;						//anizotrópiaaban a szimmetria merteke: k fogasu a szimmetria
//-----------Diszkretizalas------------------------
double d=0.00625/4.0;
double dt=1.1875e-6/2.0/(8.0*8.0);



/*-----------Fizikai parameterek-------------------*/
double V_m=7.4e-6;
double dh_f_Ni=2.35e9;		//terfogati olvadasho	
double dh_f_Cu=1.728e9;		//terfogati olvadasho
double T_m_Ni=1728;
double T_m_Cu=1358;
// double d0=4.1577879e-8; 
double d0=1e-9/(2.0*sqrt(2)*log(9))*d0_coeff; 
double gamma_Ni=0.052484417;  
double D_l=1.0e-9;
double D_s=0.0;
double N_A=6.022e23;
double a_0=pow(V_m/N_A,1.0/3.0); //=2.22*10^(-10)
double kszi=5.2e-8*d0_coeff;
/*------------Konstansok--------------------------*/

double R=8.3144087;
double k_B=R/N_A;
double s_0=0.03;
double s_kinetic=0.5;

/*-----Szamitott mennyisegek------------------------*/

double M_phi=1.0/(6.0*sqrt(2.0))*D_l*V_m/(a_0*d0*R*T);
double eps_2=6.0*sqrt(2.0)*gamma_Ni*d0/T_m_Ni;
double tau=kszi*kszi/D_l;
double dim_dx=d*d0;
double dim_dt=dt*tau;

double A_noise=sqrt(2.0*M_phi*k*T/(dim_dx*dim_dx*dim_dx*dim_dt)); 
	


/*-----Dimenziotlan parameterek------------------*/

//	phase-field egyenletben szereplo egyutthatok:
double dless_M_phi=0.9;						//=1/(6*sqrt(2))*(V_m*eps_2/(a_0*d0*R))
double dless_W=6.0*sqrt(2.0)*gamma_Ni/(d0*T_m_Ni)*kszi*kszi/eps_2;		//
double dless_dh_f_Ni=dh_f_Ni*kszi*kszi/(eps_2*T);		// 
double dless_df_f_Cu=dh_f_Cu*kszi*kszi/(eps_2*T);		// 


// koncentracios mezo egyenleteben szereplo egyutthatok:
double dless_E_t_0=(R*T/V_m)*kszi*kszi/(eps_2*T);		// 
double dless_Mc=V_m/(R*T)*eps_2*T/(kszi*kszi);			// 

double dless_E_t_2=dless_E_t_0*dless_Mc;			// 
double dless_df_f_Cu_2=dless_df_f_Cu*dless_Mc;			// 
double dless_dh_f_Ni_2=dless_dh_f_Ni*dless_Mc;			// 
double Ds_rat=0.01;

/*Az orientacios mezohoz kapcsolodo egyutthatok*/
double dless_HT=4.0*gamma_Ni*kszi/(eps_2*T); 			
double dless_M_theta_s=7.2e-4;
double dless_M_theta_l=720.0;
  //Plapp-modell parameterei:
double dless_HT_Plapp=H_coeff*gamma_Ni*d0/(eps_2*T_m_Ni);
double dless_M_theta_l_Plapp=dless_M_theta_l*4.0*d0/kszi;
double dless_M_theta_s_Plapp=dless_M_theta_l_Plapp*1.0e-4;


double A_phase_noise=2.5e-3;
double A_orientation_noise=0.025;

//*************************************************************************
//                                 _    _ _           _            
//                               | | _(_|_)_ __ __ _| |_ __ _ ___ 
//                              | |/ / | | '__/ _` | __/ _` / __|
//                             |   <| | | | | (_| | || (_| \__ \
//                            |_|\_\_|_|_|  \__,_|\__\__,_|___/
//                                  
//*************************************************************************


  
  ofstream out1,out2,out3,out4;
  out1.open("other_parameters");
  out1<<d<<"\t"<<dt<<"\t"<<T<<"\t"<<d0<<"\t"<<gamma_Ni<<"\t"<<kszi<<"\t"<<c_s<<"\t"<<c_l<<"\t"<<c0<<endl;

  out2.open("phase_field_parameters");
  out2<<dless_M_phi<<"\t"<<dless_W<<"\t"<<dless_df_f_Cu<<"\t"<<dless_dh_f_Ni<<"\t"<<T_m_Cu<<"\t"<<T_m_Ni<<"\t"<<dless_HT_Plapp<<"\t"<<dless_HT<<"\t"<<A_phase_noise<<"\t"<<k<<"\t"<<s_0<<endl;

  out3.open("concentration_parameters");
  out3<<Ds_rat<<"\t"<<dless_Mc<<"\t"<<dless_E_t_0<<endl;

  out4.open("orientation_parameters");
  out4<<dless_M_theta_s_Plapp<<"\t"<<dless_M_theta_l_Plapp<<"\t"<<dless_M_theta_s<<"\t"<<dless_M_theta_l<<endl;
  
}