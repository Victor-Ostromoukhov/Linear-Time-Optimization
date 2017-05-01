#define  OMP_ENABLED 1
// #undef  OMP_ENABLED


/* Discrepancy _ Gotta go fast.
© Quentin FERRE ©
@: quentin.ferre@outlook.fr
2015
**************************************************************************************************************************************************************************
Etat de la version actuelle:

1D:ok
2D:ok
3D:ok

sequence:
1D:ok
2D:ok
3D:ok

**************************************************************************************************************************************************************************

Permet, à partir d'une suite de point de calculer la discrepance.
Dimensions implémentées : 1,2 et 3D.
Possibilité d'obtenir le résultat sous forme de séquence.

**************************************************************************************************************************************************************************

Fonctions importantes:
-void quickSort(double tableau[], int debut, int fin):
permet de trier un tableau de dimension 1 dans l'ordre croissant.

-void import_options(int argc, char** argv):
lit les informations passées dans la ligne de commande et les interprète.

-bool import_data_xD():
importe les données de dimension x et retourne true si les données ont bien été chargées.

-double calcul_discrepancy_xD():
retourne la discrepance de la suite, parallélisée au niveau le plus bas(dernière boucle des axes).

-void export_data(double):
exporte le résultat du calcul de discrepance dans le fichier de sortie renseignée dans les paramètres.

-sequence_ ... :
certaines de ces fonctions se retrouvent avec sequence_ suivi du nom de la fonction: ce sont des fonctions avec le même objectif mais
adaptées pour le traitement des données en vue du calcul de la séquence.

**************************************************************************************************************************************************************************

Variables importantes:

-Points_xD *points_tab_xD:
tableau de points de dimension x. Contient tous les points dans l'ordre du fichier d'entrée.

-int pts_number:
renseigne sur le nombre de points de la suite.

-std::string fn_input:
nom du fichier d'entrée.

-std::string fn_output:
nom du fichier de sortie.

-int dimension:
dimension dans laquelle on souhaite que le programme évolue.

-double *sortedX:
tableau 1D contenant la liste des coordonnées de l'axe x triées dans l'ordre croissant.

-double *index_?_sort?:
suite de tableaux permettant à partir de l'index d'une des coordonnées de retrouver les autres corrdonnées du point correspondant.

-int* tmp_ligne:
permet sur le plan 2D <i,j> de retenir le nombre de points dans la boite [i-1,j-1].

-int* tmp_ligne_maj:
permet de mettre à jour tmp_ligne.

-int **tmp_box:
permet, pour une boite dans un repere 3d <i,j,k>, de connaitre le nombre de point dans la boite [i,j,k-1].

-int **tmp_box_maj:
permet de mettre à jour tmp_box.

-bool sequence:
permet de savoir si l'option sequence est activée.

-int nb_threads:
permet de savoir sur combien de threads on souhaite executer le programme.

**************************************************************************************************************************************************************************

METHODE AJOUT DIMENSION :

**************************************************************************************************************************************************************************

*/

#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <limits>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <sys/resource.h>
#include <boost/program_options.hpp>
#include <vector>
#ifdef OMP_ENABLED
  #include <omp.h>
#endif
#include <math.h>

#include "Points_1D.h"
#include "Points_2D.h"
#include "Points_3D.h"

namespace boostPO = boost::program_options;


Points_1D *points_tab_1D ;
Points_2D *points_tab_2D ;
Points_3D *points_tab_3D ;
int pts_number;
std::string fn_input;
std::string fn_output;
int dimension;

time_t td,tf;

double *tab_x;

int* tmp_ligne;
int* tmp_ligne_maj;

int **tmp_box;
int **tmp_box_maj;

bool sequence;
int nb_threads = 1;

int *sortedX;
int *sortedY;
int *sortedZ;

int *index_x_sortY;
int *index_x_sortZ;

int *index_y_sortX;
int *index_y_sortZ;

int *index_z_sortX;
int *index_z_sortY;

///fonctions de tri

bool sortXfunc(int i, int j){

    if(dimension==1){
        return (points_tab_1D[i].get_pos_x() < points_tab_1D[j].get_pos_x());
    }
    else if(dimension==2){
        return (points_tab_2D[i].get_pos_x() < points_tab_2D[j].get_pos_x());
    }
    else {//dim3
        return (points_tab_3D[i].get_pos_x() < points_tab_3D[j].get_pos_x());
    }
}

bool sortYfunc(int i, int j){

    if(dimension==2){
        return (points_tab_2D[i].get_pos_y() < points_tab_2D[j].get_pos_y());
    }
    else {//dim3
        return (points_tab_3D[i].get_pos_y() < points_tab_3D[j].get_pos_y());
    }
}

bool sortZfunc(int i, int j){
    return (points_tab_3D[i].get_pos_z() < points_tab_3D[j].get_pos_z());//dim3
}

///fonctions calcul mesure de lebesgue

double box_area_calcul(double x, double y){
	return x*y;
}

double box_area_calcul(double x, double y,double z){
	return x*y*z;
}

///gestion des options

void import_options(int argc, char** argv){

	boostPO::variables_map vm;
	boostPO::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "makes help message")
		("input,i", boostPO::value(&fn_input)->required(), "input filename")
		("output,o", boostPO::value(&fn_output)->required(), "output filename")
		("dimension,d", boostPO::value(&dimension)->required(), "dimension (presently: 1/2/3)")
		("sequence,s", "sequence enable")
#ifdef OMP_ENABLED
		("threads,t" , boostPO::value(&nb_threads)->default_value(omp_get_max_threads()), "number of threads")
#endif
		;

	boostPO::positional_options_description p;
	p.add("input", -1);

	try
	{
		boostPO::store(
			boostPO::command_line_parser(argc, argv).
			options(desc).positional(p).run(), vm);
		boostPO::notify(vm);
	}
	catch(boost::program_options::error& e)
	{
		std::cerr << desc << std::endl;
		exit(EXIT_FAILURE);
	}

	if(vm.count("help"))
	{
		std::cout << desc << std::endl;
		exit(EXIT_SUCCESS);
	}

	if(vm.count("sequence"))
	{
		std::cout << std::endl<<" Sequence : ON" << std::endl;
		sequence = true;
	}
	else
	{
		std::cout << std::endl<<" Sequence : OFF" << std::endl;
		sequence = false;
	}

#ifdef OMP_ENABLED
	if(nb_threads>omp_get_max_threads())
	{
		std::cerr << "Error Threads : Number of threads available : " <<omp_get_max_threads()<< std::endl;
		exit(EXIT_FAILURE);
	}
	else
		std::cout <<" Threads : " <<nb_threads<< std::endl;
#endif

	std::cout << " Dimension : " <<dimension<< std::endl<< std::endl;
}

///fonctions 1d non sequentielles

bool import_data_1D(){

	std::string line;
	std::ifstream fichier(fn_input.c_str(), std::ios::in);
	if(!fichier.fail())
	{
		std::cout << "  ***************************************"<<std::endl<<std::endl ;
		std::cout << " Chargement des donnees..." <<std::endl;
		pts_number = 0;

		while (std::getline(fichier, line)){//on compte le nombre de points dans le fichier d'entrée
			if(line != "" && line !="\n" && line!=" "){
				pts_number++;
			}
			else continue;
		}

		std::cout <<" "<<pts_number<<" lignes."<<std::endl;
		if(pts_number>0){
			fichier.clear();
			fichier.seekg(0, std::ios::beg);

			double x=0;
			int i=0;

			points_tab_1D = new Points_1D[pts_number];
            sortedX=new int[pts_number];

			while(fichier >> x )
			{
				points_tab_1D[i].set_pos_x(x);//on stocke les points dans le tableau de points
				i++;
			}
			std::cout <<std::endl<< " Data Ok !" <<std::endl<<std::endl;
		}

		fichier.close();
		return true;
	}
	else {
		std::cerr << " Can't open input !"<<std::endl ;
		return false;
	}

}

double calcul_discrepancy_1D(){

	////////////////////////// INITIALISATION /////////////////////////////////

	std::cout << " Calcul discrepance ..."<<std::endl<<std::endl ;
	double discrepancy=0;
	const int const_pt_num = pts_number;

    #pragma omp parallel for num_threads(nb_threads)
    for(int i=0; i<pts_number; i++) sortedX[i] = i;
	std::sort(&sortedX[0],&sortedX[pts_number],sortXfunc);

	////////////////////////// CALCUL DISCREPANCE /////////////////////////////////

	#pragma omp parallel for num_threads(nb_threads)
	for(int j=0;j<const_pt_num;j++)
	{
		int memoire=j;
		double tmp=fabs(points_tab_1D[sortedX[j]].get_pos_x()-((2*(j+1)-1)/(double)(2*const_pt_num)));

		#pragma omp flush(discrepancy)
		if(discrepancy<=tmp)
		{
			#pragma omp crititcal
			if(discrepancy<=tmp)
			{
				discrepancy=tmp;
			}
		}
	}

	discrepancy=discrepancy+((double)1/(double)(2*pts_number));

	////////////////////// FIN DU CALCUL //////////////////////////////
	return discrepancy;
}

void export_data_1D(double discrepancy){

	std::ofstream file(fn_output.c_str());
	if(!file)
		std::cerr << "Can't open ouptut" << std::endl;

    // file << pts_number << "  " << discrepancy << std::endl;
    file << discrepancy << std::endl;
    file.close();
}

///fonctions 2d non sequentielles

bool import_data_2D(){

	std::string line;
	std::ifstream fichier(fn_input.c_str(), std::ios::in);

	if(!fichier.fail())
	{
		std::cout << "  *************************************** " <<std::endl<<std::endl;
		std::cout << " Chargement des donnees..."<<std::endl ;
		pts_number = 0;

		while (std::getline(fichier, line)){//on compte le nombre de points dans le fichier d'entrée
			if(line != "" && line !="\n" && line!=" "){
				pts_number++;
			}
			else continue;
		}

		std::cout <<" "<<pts_number<<" lignes."<<std::endl;
		if(pts_number>0){

			fichier.clear();
			fichier.seekg(0, std::ios::beg);

			double x=0;
			double y=0;
			int i=0;

			points_tab_2D = new Points_2D[pts_number];

			tmp_ligne=new int[pts_number];
			tmp_ligne_maj=new int[pts_number];

			sortedX=new int[pts_number];
            sortedY=new int[pts_number];

            index_x_sortY=new int[pts_number];
            index_y_sortX=new int[pts_number];

			while(fichier >> x >> y)
			{
				points_tab_2D[i].set_pos_x(x);//on stocke les points dans l'ordre de la suite
				points_tab_2D[i].set_pos_y(y);

				i++;
			}
			std::cout << std::endl<<" Data Ok !"<<std::endl<<std::endl ;
		}

		fichier.close();
		return true;
	}
	else {
		std::cerr << " Can't open input !"<<std::endl ;
		return false;
	}

}

double calcul_discrepancy_2D(){

	std::cout << " Calcul discrepance ..."<<std::endl<<std::endl ;
	double discrepancy=0;


	////////////////////////// INITIALISATION /////////////////////////////////

    #pragma omp parallel for num_threads(nb_threads)
     for(int i=0; i<pts_number; i++) sortedX[i] = i;
	std::sort(&sortedX[0],&sortedX[pts_number],sortXfunc);

    #pragma omp parallel for num_threads(nb_threads)
	for(int i=0; i<pts_number; i++) sortedY[i] = i;
	std::sort(&sortedY[0],&sortedY[pts_number],sortYfunc);


     for(int i=0; i<pts_number; i++)
    {
        #pragma omp parallel for num_threads(nb_threads)
         for(int j=0; j<pts_number; j++)
         {
            if(sortedX[i]==sortedY[j]) index_x_sortY[i]=j;
            if(sortedY[i]==sortedX[j]) index_y_sortX[i]=j;
         }
    }

	#pragma omp parallel for num_threads(nb_threads)
	for(int i=0;i<pts_number;i++)
	{
		tmp_ligne[i]=0;
		tmp_ligne_maj[i]=0;
	}

	const int const_pt_num = pts_number;

	////////////////////////// CALCUL DISCREPANCE /////////////////////////////////

	for(int i=0;i<pts_number;i++)
	{
		#pragma omp parallel for num_threads(nb_threads)//mise à jour de tmp_ligne
		for(int k=0;k<const_pt_num;k++)
		{
			tmp_ligne[k]=tmp_ligne_maj[k];
		}

		#pragma omp parallel for firstprivate(points_tab_2D,tmp_ligne,sortedX,sortedY,index_x_sortY,index_y_sortX) num_threads(nb_threads)
		for(int j=0;j<const_pt_num;j++)
		{
			int memoire=0;

			if(j>0)
				memoire=tmp_ligne[j-1];//nb de points ds la boite [i-1][j-1]

			//on regarde si il y a des points sur les axes x et y de la dernière subdivision de la boite qui sont ds la boite
			if(i>index_y_sortX[j])
			{
				memoire++;
			}

			if(j>index_x_sortY[i])
			{
				memoire++;
			}
			//on regarde si il y a un point ds la derniere subdivision de la boîte
			if(i==index_y_sortX[j])
			{
				memoire++;
			}
			//calcul discrepance
			double tmp=0;
			tmp=fabs((((double)memoire))/((double)const_pt_num)-box_area_calcul(points_tab_2D[sortedX[i]].get_pos_x(),points_tab_2D[sortedY[j]].get_pos_y()));
			//on vérifie les autres angles de la boîte

			if(i<const_pt_num-1 && j<const_pt_num-1)
			{
				tmp=std::max(tmp,fabs((((double)memoire))/((double)const_pt_num)-box_area_calcul(points_tab_2D[sortedX[i+1]].get_pos_x(),points_tab_2D[sortedY[j+1]].get_pos_y())));
				tmp=std::max(tmp,fabs((((double)memoire))/((double)const_pt_num)-box_area_calcul(points_tab_2D[sortedX[i+1]].get_pos_x(),points_tab_2D[sortedY[j]].get_pos_y())));
				tmp=std::max(tmp,fabs((((double)memoire))/((double)const_pt_num)-box_area_calcul(points_tab_2D[sortedX[i]].get_pos_x(),points_tab_2D[sortedY[j+1]].get_pos_y())));
			}

			else if (i==const_pt_num-1 && j==const_pt_num-1)
			{
				tmp=std::max(tmp,fabs((((double)memoire))/((double)const_pt_num)-box_area_calcul(1,1)));
			}

			else if(i==const_pt_num-1)
			{
				tmp=std::max(tmp,fabs((((double)memoire))/((double)const_pt_num)-box_area_calcul(1,points_tab_2D[sortedY[j]].get_pos_y())));
                tmp=std::max(tmp,fabs((((double)memoire))/((double)const_pt_num)-box_area_calcul(1,points_tab_2D[sortedY[j+1]].get_pos_y())));
            }

			else if (j==const_pt_num-1)
			{
                tmp=std::max(tmp,fabs((((double)memoire))/((double)const_pt_num)-box_area_calcul(points_tab_2D[sortedX[i]].get_pos_x(),1)));
                tmp=std::max(tmp,fabs((((double)memoire))/((double)const_pt_num)-box_area_calcul(points_tab_2D[sortedX[i+1]].get_pos_x(),1)));
			}

			else
				std::cout<<"impossible"<<std::endl;

			tmp_ligne_maj[j]=memoire;//on prepare la maj de tmp_ligne

			#pragma omp flush(discrepancy)
			if(discrepancy<=tmp)
			{
				#pragma omp crititcal
				if(discrepancy<=tmp)
				{
					discrepancy=tmp;
				}
			}
		}
	}


	/////////// BOUCLES SUPPLEMENTAIRES POUR TESTER LES BOITES VIDES SUR LES AXES X ET Y///////////////////////////////

	#pragma omp parallel for num_threads(nb_threads)
	for(int i=0;i<pts_number;i++)
	{
		double tmp=box_area_calcul(points_tab_2D[sortedX[i]].get_pos_x(),points_tab_2D[sortedY[0]].get_pos_y());
		#pragma omp flush(discrepancy)
		if(discrepancy<=tmp)
		{
			#pragma omp crititcal
			if(discrepancy<=tmp)
			{
				discrepancy=tmp;
			}
		}
	}

	#pragma omp parallel for num_threads(nb_threads)
	for(int i=1;i<pts_number;i++)
	{
		double tmp=box_area_calcul(points_tab_2D[sortedX[0]].get_pos_x(),points_tab_2D[sortedY[i]].get_pos_y());
		#pragma omp flush(discrepancy)
		if(discrepancy<=tmp)
		{
			#pragma omp crititcal
			if(discrepancy<=tmp)
			{
				discrepancy=tmp;
			}
		}
	}

    discrepancy=std::max (discrepancy,fabs(box_area_calcul(points_tab_2D[sortedX[0]].get_pos_x(),1)));
    discrepancy=std::max (discrepancy,fabs(box_area_calcul(1,points_tab_2D[sortedY[0]].get_pos_y())));

	////////////////////// FIN DU CALCUL //////////////////////////////

	return discrepancy;
}

void export_data_2D(double discrepancy){

	std::ofstream file(fn_output.c_str());
	if(!file)
		std::cerr << "Can't open ouptut" << std::endl;

    //file << pts_number << "  " << discrepancy << std::endl;
    file << discrepancy << std::endl;
    file.close();
}

///fonctions 3d non sequentielles

bool import_data_3D(){

	std::string line;

	std::ifstream fichier(fn_input.c_str(), std::ios::in);
	if(!fichier.fail())
	{
		std::cout << "  ***************************************  "<<std::endl<<std::endl ;
		std::cout << " Chargement des donnees..."<<std::endl ;
		pts_number = 0;

		while (std::getline(fichier, line)){
			if(line != "" && line !="\n" && line!=" "){
				pts_number++;
			}
			else continue;
		}

		std::cout <<" "<<pts_number<<" lignes."<<std::endl;
		if(pts_number>0){

			fichier.clear();
			fichier.seekg(0, std::ios::beg);

			double x=0;
			double y=0;
			double z=0;
			int i=0;

			points_tab_3D = new Points_3D[pts_number];

            sortedX=new int[pts_number];
            sortedY=new int[pts_number];
            sortedZ=new int[pts_number];

            index_x_sortY=new int[pts_number];
            index_x_sortZ=new int[pts_number];

            index_y_sortX=new int[pts_number];
            index_y_sortZ=new int[pts_number];

            index_z_sortX=new int[pts_number];
            index_z_sortY=new int[pts_number];

			tmp_ligne=new int[pts_number];
			tmp_ligne_maj=new int[pts_number];

			tmp_box=new int*[pts_number];
			tmp_box_maj=new int*[pts_number];

			for(int i=0;i<pts_number;i++)
			{
				tmp_box[i]=new int [pts_number];
				tmp_box_maj[i]=new int [pts_number];
			}

			while(fichier >> x >> y >> z)
			{
				points_tab_3D[ i ].set_pos_x( x );//on stocke les points dans l'ordre de la suite
				points_tab_3D[ i ].set_pos_y( y );
				points_tab_3D[ i ].set_pos_z( z );
				i++;
			}
			std::cout <<std::endl<< " Data Ok !"<<std::endl<<std::endl ;
		}
		fichier.close();
		return true;
	}
	else {
		std::cerr << " Can't open input !"<<std::endl ;
		return false;
	}

}

double calcul_discrepancy_3D(){

	std::cout << " Calcul discrepance ..."<<std::endl<<std::endl ;
	double discrepancytmp=0;

////////////////////////////////////////////////////////////////////////////////////
    #pragma omp parallel for num_threads(nb_threads)
    for(int i=0; i<pts_number; i++) sortedX[i] = i;
	std::sort(&sortedX[0],&sortedX[pts_number],sortXfunc);

    #pragma omp parallel for num_threads(nb_threads)
	for(int i=0; i<pts_number; i++) sortedY[i] = i;
	std::sort(&sortedY[0],&sortedY[pts_number],sortYfunc);

    #pragma omp parallel for num_threads(nb_threads)
	for(int i=0; i<pts_number; i++) sortedZ[i] = i;
	std::sort(&sortedZ[0],&sortedZ[pts_number],sortZfunc);

    for(int i=0; i<pts_number; i++)
    {
        #pragma omp parallel for num_threads(nb_threads)
         for(int j=0; j<pts_number; j++)
         {
            if(sortedX[i]==sortedY[j]) index_x_sortY[i]=j;
            if(sortedX[i]==sortedZ[j]) index_x_sortZ[i]=j;

            if(sortedY[i]==sortedX[j]) index_y_sortX[i]=j;
            if(sortedY[i]==sortedZ[j]) index_y_sortZ[i]=j;

            if(sortedZ[i]==sortedX[j]) index_z_sortX[i]=j;
            if(sortedZ[i]==sortedY[j]) index_z_sortY[i]=j;
         }
    }

    //////////////////////////////////////////////////////////////////////////////////////////

	#pragma omp parallel for num_threads(nb_threads)
	for(int i=0;i<pts_number;i++)
	{
		tmp_ligne[i]=0;
		tmp_ligne_maj[i]=0;
	}

	/////////////////////////////////////////////////////////////////////////////////////////

	for(int i=0;i<pts_number;i++)
	{
		#pragma omp parallel for num_threads(nb_threads)
		for(int j=0;j<pts_number;j++)
		{
			tmp_box_maj[i][j]=0;
		}
	}

	const int const_pt_num = pts_number;

	////////////////////////// DISCREPANCE //////////////////////////////////////////////////

	for(int k=0; k<pts_number;k++){

		for(int o=0; o<pts_number;o++){//maj de tmp_box
			#pragma omp parallel for
			for(int p=0; p<pts_number;p++){
				tmp_box[o][p]=tmp_box_maj[o][p];

			}
		}

		#pragma omp parallel for num_threads(nb_threads)//raz du tableau de mise à jour
		for(int l=0;l<pts_number;l++)
		{
			tmp_ligne_maj[l]=0;
		}

		for(int i=0;i<pts_number;i++)
		{
			#pragma omp parallel for num_threads(nb_threads)//maj de tmp_ligne
			for(int l=0;l<pts_number;l++)
			{
				tmp_ligne[l]=tmp_ligne_maj[l];
			}

			#pragma omp parallel for firstprivate(points_tab_3D, tmp_ligne,tmp_box,sortedX,sortedY,sortedZ,index_x_sortY,index_x_sortZ,index_y_sortX,index_y_sortZ,index_z_sortX,index_z_sortY) num_threads(nb_threads)
			for(int j=0;j<const_pt_num;j++)
			{
				int memoire=0;
				int etage_inf=0;

				if(j>0)
                memoire=tmp_ligne[j-1];

				etage_inf=tmp_box[i][j];


				if(i>index_y_sortX[j] && index_y_sortZ[j]==k) /////////////////////////////dessus
				memoire++;

                if(j>index_x_sortY[i] && index_x_sortZ[i]==k)  //////////////////////////////cote
				memoire++;

				//on regarde si il y a un point ds la derniere subdivision de la boîte
				if(i==index_y_sortX[j] && index_y_sortZ[j]==k)
				memoire++;


				double tmp=0;
				//calcul discrepance
				tmp=fabs((((double)etage_inf+memoire))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z()));

                double sum=etage_inf+memoire;
				//on vérifie les autres angles de la boite
				if(i<const_pt_num-1 && j<const_pt_num-1 && k<const_pt_num-1)
				{
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j+1]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),points_tab_3D[sortedY[j+1]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k+1]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k+1]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j+1]].get_pos_y(),points_tab_3D[sortedZ[k+1]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),points_tab_3D[sortedY[j+1]].get_pos_y(),points_tab_3D[sortedZ[k+1]].get_pos_z())));
				}

				else if (i==const_pt_num-1 && j==const_pt_num-1 && k== const_pt_num-1)
				{
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(1,1,1)));
				}

				else if(i==const_pt_num-1 && j==const_pt_num-1)
				{
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(1,1,points_tab_3D[sortedZ[k]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k+1]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(1,1,points_tab_3D[sortedZ[k+1]].get_pos_z())));
				}

				else if(k==const_pt_num-1 && j==const_pt_num-1)
				{
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),1,1)));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),1,1)));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z())));
				}

				else if(i==const_pt_num-1 && k==const_pt_num-1)
				{
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(1,points_tab_3D[sortedY[j]].get_pos_y(),1)));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(1,points_tab_3D[sortedY[j+1]].get_pos_y(),1)));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j+1]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z())));
				}

				else if(i==const_pt_num-1)
				{
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(1,points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z())));

					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j+1]].get_pos_y(),points_tab_3D[sortedZ[k+1]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(1,points_tab_3D[sortedY[j+1]].get_pos_y(),points_tab_3D[sortedZ[k+1]].get_pos_z())));

					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k+1]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(1,points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k+1]].get_pos_z())));

					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j+1]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(1,points_tab_3D[sortedY[j+1]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z())));
				}
				else if(j==const_pt_num-1 )
				{
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),1,points_tab_3D[sortedZ[k]].get_pos_z())));

					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k+1]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),1,points_tab_3D[sortedZ[k+1]].get_pos_z())));

					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k+1]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),1,points_tab_3D[sortedZ[k+1]].get_pos_z())));

					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),1,points_tab_3D[sortedZ[k]].get_pos_z())));

				}
				else if(k==const_pt_num-1)
				{
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),1)));

                    tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),points_tab_3D[sortedY[j+1]].get_pos_y(),1)));
                    tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),points_tab_3D[sortedY[j+1]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z())));

                    tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),1)));
                    tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z())));

                    tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j+1]].get_pos_y(),1)));
                    tmp=std::max(tmp,fabs((((double)sum))/((double)const_pt_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j+1]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z())));

				}
				else
					std::cout<<"impossible"<<std::endl;
				//on prépare la maj de tmp_ligne et tmp_box
				tmp_ligne_maj[j]=memoire;
				tmp_box_maj[i][j]=sum;

				//comparer max
				#pragma omp flush(discrepancytmp)
				if(discrepancytmp<=tmp )
				{
					#pragma omp crititcal
					if(discrepancytmp<=tmp  )
					{
						discrepancytmp=tmp;
					}
				}
			}///fin j
		}///fin i
	}///fin k

	/////////// BOUCLES SUPPLEMENTAIRES POUR TESTER LES BOITES VIDES SUR LES AXES X, Y ET Z///////////////////////////////

	#pragma omp parallel for num_threads(nb_threads)
	for(int i=0;i<pts_number;i++)
	{
		double tmp=box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[0]].get_pos_y(),points_tab_3D[sortedZ[0]].get_pos_z());
		#pragma omp flush(discrepancytmp)
		if(discrepancytmp<=tmp)
		{
			#pragma omp crititcal
			if(discrepancytmp<=tmp)
			{
				discrepancytmp=tmp;
			}
		}
	}

	#pragma omp parallel for num_threads(nb_threads)
	for(int i=1;i<pts_number;i++)
	{
		double tmp=box_area_calcul(points_tab_3D[sortedX[0]].get_pos_x(),points_tab_3D[sortedY[i]].get_pos_y(),points_tab_3D[sortedZ[0]].get_pos_z());
		#pragma omp flush(discrepancytmp)
		if(discrepancytmp<=tmp)
		{
			#pragma omp crititcal
			if(discrepancytmp<=tmp)
			{
				discrepancytmp=tmp;
			}
		}
	}

	#pragma omp parallel for num_threads(nb_threads)
	for(int i=1;i<pts_number;i++)
	{
		double tmp=box_area_calcul(points_tab_3D[sortedX[0]].get_pos_x(),points_tab_3D[sortedY[0]].get_pos_y(),points_tab_3D[sortedZ[0]].get_pos_z());
		#pragma omp flush(discrepancytmp)
		if(discrepancytmp<=tmp)
		{
			#pragma omp crititcal
			if(discrepancytmp<=tmp)
			{
				discrepancytmp=tmp;
			}
		}
	}

    discrepancytmp=std::max (discrepancytmp,fabs(box_area_calcul(points_tab_3D[sortedX[0]].get_pos_x(),1,1)));
    discrepancytmp=std::max (discrepancytmp,fabs(box_area_calcul(1,1,points_tab_3D[sortedZ[0]].get_pos_z())));
    discrepancytmp=std::max (discrepancytmp,fabs(box_area_calcul(1,points_tab_3D[sortedY[0]].get_pos_y(),1)));


	////////////////////// FIN DU CALCUL //////////////////////////////

	return discrepancytmp;
}

void export_data_3D(double discrepancy){

	std::ofstream file(fn_output.c_str());
	if(!file)
		std::cerr << "Can't open ouptut" << std::endl;

    //file << pts_number << "  " << discrepancy << std::endl;
    file << discrepancy << std::endl;
    file.close();
}

///fonctions 1d sequentielles

bool sequence_import_data_1D(int palier){

	std::string line;
	std::ifstream fichier(fn_input.c_str(), std::ios::in);
	if(!fichier.fail())
	{
		std::cout << "  *************************************** "<<std::endl<<std::endl ;
		std::cout << " Chargement des donnees..."<<std::endl ;
		pts_number = 0;

		while (std::getline(fichier, line)){
			if(line != "" && line !="\n" && line!=" "){
				pts_number++;
			}
			else continue;
		}

		std::cout <<" "<<pts_number<<" lignes."<<std::endl;
		if(pts_number>0){

			fichier.clear();
			fichier.seekg(0, std::ios::beg);

			double x=0;
			int i=0;

			points_tab_1D = new Points_1D[pts_number];
            sortedX=new int[palier];

			while(fichier >> x )
			{
				points_tab_1D[i].set_pos_x(x);
				i++;
			}
			std::cout << std::endl<<" Data Ok !"<<std::endl<<std::endl ;
		}

		fichier.close();
		return true;
	}
	else {
		std::cerr << " Can't open input !" <<std::endl;
		return false;
	}

}

double sequence_calcul_discrepancy_1D(int palier){

	double discrepancy=0;
	////////////////////////// DISCREPANCE //////////////////////////////

	#pragma omp parallel for num_threads(nb_threads)
	for(int j=0;j<palier;j++)
	{
		int memoire=j;
		double tmp=fabs(points_tab_1D[sortedX[j]].get_pos_x()-((2*(j+1)-1)/(double)(2*palier)));

		#pragma omp flush(discrepancy)
		if(discrepancy<=tmp)
		{
			#pragma omp crititcal
			if(discrepancy<=tmp)
			{
				discrepancy=tmp;
			}
		}
	}

	discrepancy=discrepancy+((double)1/(double)(2*palier));

	////////////////////// FIN DU CALCUL //////////////////////////////

	return discrepancy;
}

///fonctions 2d sequentielles

bool sequence_import_data_2D(int palier){

	std::string line;

	std::ifstream fichier(fn_input.c_str(), std::ios::in);
	if(!fichier.fail())
	{
		std::cout << "  *************************************** " <<std::endl<<std::endl;

		std::cout << " Chargement des donnees..." <<std::endl;
		pts_number = 0;

		while (std::getline(fichier, line)){
			if(line != "" && line !="\n" && line!=" "){
				pts_number++;
			}
			else continue;
		}

		std::cout <<" "<<pts_number<<" lignes."<<std::endl;
		if(pts_number>0){

			fichier.clear();
			fichier.seekg(0, std::ios::beg);

			double x=0;
			double y=0;
			int i=0;

			points_tab_2D = new Points_2D[pts_number];
            sortedX=new int[palier];
            sortedY=new int[palier];

			tmp_ligne=new int[palier];
			tmp_ligne_maj=new int[palier];

			while(fichier >> x >> y)
			{
				points_tab_2D[i].set_pos_x(x);
				points_tab_2D[i].set_pos_y(y);
				i++;
			}
			std::cout << std::endl<<" Data Ok !"<<std::endl<<std::endl ;
		}

		fichier.close();
		return true;
	}
	else {
		std::cerr << " Can't open input !"<<std::endl ;
		return false;
	}

}

double sequence_calcul_discrepancy_2D(int palier){

	double discrepancy=0;

     for(int i=0; i<palier; i++)
    {
        #pragma omp parallel for num_threads(nb_threads)
         for(int j=0; j<palier; j++)
         {
            if(sortedX[i]==sortedY[j]) index_x_sortY[i]=j;
            if(sortedY[i]==sortedX[j]) index_y_sortX[i]=j;
         }
    }
	///////////////////////////////////////////////////////////////////////////

	#pragma omp parallel for num_threads(nb_threads)
	for(int i=0;i<palier;i++)
	{
		tmp_ligne[i]=0;
		tmp_ligne_maj[i]=0;
	}

	const int const_palier_num = palier;

	////////////////////////// DISCREPANCE //////////////////////////////////////

	for(int i=0;i<palier;i++)
	{

		#pragma omp parallel for num_threads(nb_threads)
		for(int k=0;k<palier;k++)
		{
			tmp_ligne[k]=tmp_ligne_maj[k];
		}

		#pragma omp parallel for firstprivate(points_tab_2D,tmp_ligne,sortedX,sortedY,index_x_sortY,index_y_sortX) num_threads(nb_threads)
		for(int j=0;j<const_palier_num;j++)
		{
			int memoire=0;

			if(j>0)
				memoire=tmp_ligne[j-1];

			if(i>index_y_sortX[j])
			{
				memoire++;
			}

			if(j>index_x_sortY[i])
			{
				memoire++;
			}

			if(i==index_y_sortX[j])
			{
				memoire++;
			}

			double tmp=0;
			tmp=fabs((((double)memoire))/((double)const_palier_num)-box_area_calcul(points_tab_2D[sortedX[i]].get_pos_x(),points_tab_2D[sortedY[j]].get_pos_y()));
			//on vérifie les autres angles de la boîte
			if(i<const_palier_num-1 && j<const_palier_num-1)
			{
				tmp=std::max(tmp,fabs((((double)memoire))/((double)const_palier_num)-box_area_calcul(points_tab_2D[sortedX[i+1]].get_pos_x(),points_tab_2D[sortedY[j+1]].get_pos_y())));
				tmp=std::max(tmp,fabs((((double)memoire))/((double)const_palier_num)-box_area_calcul(points_tab_2D[sortedX[i+1]].get_pos_x(),points_tab_2D[sortedY[j]].get_pos_y())));
				tmp=std::max(tmp,fabs((((double)memoire))/((double)const_palier_num)-box_area_calcul(points_tab_2D[sortedX[i]].get_pos_x(),points_tab_2D[sortedY[j+1]].get_pos_y())));
			}

			else if (i==const_palier_num-1 && j==const_palier_num-1)
			{
				tmp=std::max(tmp,fabs((((double)memoire))/((double)const_palier_num)-box_area_calcul(1,1)));
			}

			else if(i==const_palier_num-1)
			{
				tmp=std::max(tmp,fabs((((double)memoire))/((double)const_palier_num)-box_area_calcul(1,points_tab_2D[sortedY[j]].get_pos_y())));
                tmp=std::max(tmp,fabs((((double)memoire))/((double)const_palier_num)-box_area_calcul(1,points_tab_2D[sortedY[j+1]].get_pos_y())));
            }

			else if (j==const_palier_num-1)
			{
                tmp=std::max(tmp,fabs((((double)memoire))/((double)const_palier_num)-box_area_calcul(points_tab_2D[sortedX[i]].get_pos_x(),1)));
                tmp=std::max(tmp,fabs((((double)memoire))/((double)const_palier_num)-box_area_calcul(points_tab_2D[sortedX[i+1]].get_pos_x(),1)));
			}

			else
				std::cout<<"impossible"<<std::endl;

			tmp_ligne_maj[j]=memoire;//on prepare la maj de tmp_ligne

			#pragma omp flush(discrepancy)
			if(discrepancy<=tmp)
			{
				#pragma omp crititcal
				if(discrepancy<=tmp)
				{
					discrepancy=tmp;
				}
			}
		}
	}


	/////////// BOUCLES SUPPLEMENTAIRES POUR TESTER LES BOITES VIDES SUR LES AXES X ET Y///////////////////////////////

	#pragma omp parallel for num_threads(nb_threads)
	for(int i=0;i<palier;i++)
	{
		double tmp=box_area_calcul(points_tab_2D[sortedX[i]].get_pos_x(),points_tab_2D[sortedY[0]].get_pos_y());
		#pragma omp flush(discrepancy)
		if(discrepancy<=tmp)
		{
			#pragma omp crititcal
			if(discrepancy<=tmp)
			{
				discrepancy=tmp;
			}
		}
	}

	#pragma omp parallel for num_threads(nb_threads)
	for(int i=1;i<palier;i++)
	{
		double tmp=box_area_calcul(points_tab_2D[sortedX[0]].get_pos_x(),points_tab_2D[sortedY[i]].get_pos_y());
		#pragma omp flush(discrepancy)
		if(discrepancy<=tmp)
		{
			#pragma omp crititcal
			if(discrepancy<=tmp)
			{
				discrepancy=tmp;
			}
		}
	}

    discrepancy=std::max (discrepancy,fabs(box_area_calcul(points_tab_2D[sortedX[0]].get_pos_x(),1)));
    discrepancy=std::max (discrepancy,fabs(box_area_calcul(1,points_tab_2D[sortedY[0]].get_pos_y())));

	////////////////////// FIN DU CALCUL //////////////////////////////

	return discrepancy;
}

///fonctions 3d sequentielles

bool sequence_import_data_3D(int palier){

	std::string line;

	std::ifstream fichier(fn_input.c_str(), std::ios::in);
	if(!fichier.fail())
	{
		std::cout << "  *************************************** "<<std::endl<<std::endl ;
		std::cout << " Chargement des donnees..."<<std::endl ;
		pts_number = 0;

		while (std::getline(fichier, line)){
			if(line != "" && line !="\n" && line!=" "){
				pts_number++;
			}
			else continue;
		}

		std::cout <<" "<<pts_number<<" lignes."<<std::endl;
		if(pts_number>0){

			fichier.clear();
			fichier.seekg(0, std::ios::beg);

			double x=0;
			double y=0;
			double z=0;
			int i=0;

			points_tab_3D = new Points_3D[pts_number];
            sortedX=new int[palier];
            sortedY=new int[palier];
            sortedZ=new int[palier];

            index_x_sortY=new int[palier];
            index_x_sortZ=new int[palier];

            index_y_sortX=new int[palier];
            index_y_sortZ=new int[palier];

            index_z_sortX=new int[palier];
            index_z_sortY=new int[palier];

			tmp_ligne=new int[palier];
			tmp_ligne_maj=new int[palier];

			tmp_box=new int*[palier];
			tmp_box_maj=new int*[palier];

			for(int i=0;i<palier;i++)
			{
				tmp_box[i]=new int [palier];
				tmp_box_maj[i]=new int [palier];

			}

			while(fichier >> x >> y >> z)
			{
				points_tab_3D[ i ].set_pos_x( x );
				points_tab_3D[ i ].set_pos_y( y );
				points_tab_3D[ i ].set_pos_z( z );
				i++;
			}
			std::cout <<std::endl<< " Data Ok !"<<std::endl<<std::endl ;
		}
		fichier.close();
		return true;
	}
	else {
		std::cerr << " Can't open input !"<<std::endl ;
		return false;
	}

}

double sequence_calcul_discrepancy_3D(int palier){

	double discrepancytmp=0;

    for(int i=0; i<palier; i++)
    {
        #pragma omp parallel for num_threads(nb_threads)
         for(int j=0; j<palier; j++)
         {
            if(sortedX[i]==sortedY[j]) index_x_sortY[i]=j;
            if(sortedX[i]==sortedZ[j]) index_x_sortZ[i]=j;

            if(sortedY[i]==sortedX[j]) index_y_sortX[i]=j;
            if(sortedY[i]==sortedZ[j]) index_y_sortZ[i]=j;

            if(sortedZ[i]==sortedX[j]) index_z_sortX[i]=j;
            if(sortedZ[i]==sortedY[j]) index_z_sortY[i]=j;
         }
    }
	///////////////////////////////////////////////////////////////////////////

	#pragma omp parallel for num_threads(nb_threads)
	for(int i=0;i<palier;i++)
	{
		tmp_ligne[i]=0;
		tmp_ligne_maj[i]=0;
	}
	/////////////////////////////////////////////////////////////////////////////////////////

	for(int i=0;i<palier;i++)
	{
		#pragma omp parallel for num_threads(nb_threads)
		for(int j=0;j<palier;j++)
		{
			tmp_box[i][j]=0;
			tmp_box_maj[i][j]=0;
		}
	}

	const int const_palier_num = palier;

	for(int k=0; k<palier;k++){

		for(int o=0; o<palier;o++){
			#pragma omp parallel for
			for(int p=0; p<palier;p++){
				tmp_box[o][p]=tmp_box_maj[o][p];

			}
		}
		#pragma omp parallel for num_threads(nb_threads)
		for(int l=0;l<palier;l++)
		{
			tmp_ligne_maj[l]=0;
		}

		for(int i=0;i<palier;i++)
		{
			#pragma omp parallel for num_threads(nb_threads)
			for(int l=0;l<palier;l++)
			{
				tmp_ligne[l]=tmp_ligne_maj[l];
			}


			#pragma omp parallel for firstprivate(points_tab_3D,tmp_ligne,tmp_box,sortedX,sortedY,sortedZ,index_x_sortY,index_x_sortZ,index_y_sortX,index_y_sortZ,index_z_sortX,index_z_sortY) num_threads(nb_threads)
			for(int j=0;j<const_palier_num;j++)
			{

				int memoire=0;
				int etage_inf=0;

				if(j>0)
				{
					memoire=tmp_ligne[j-1];
				}

				etage_inf=tmp_box[i][j];

				if(i>index_y_sortX[j] && index_y_sortZ[j]==k) /////////////////////////////dessus
				{
					memoire++;
				}

                if(j>index_x_sortY[i] && index_x_sortZ[i]==k)  //////////////////////////////cote
				{
					memoire++;
				}

				if(i==index_y_sortX[j] && index_y_sortZ[j]==k)
				{
					memoire++;
				}

				double tmp=0;
				///calcul discrepance
				tmp=fabs((((double)etage_inf+memoire))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z()));

                double sum=etage_inf+memoire;
				//on vérifie les autres angles de la boite
				if(i<const_palier_num-1 && j<const_palier_num-1 && k<const_palier_num-1)
				{
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j+1]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),points_tab_3D[sortedY[j+1]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k+1]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k+1]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j+1]].get_pos_y(),points_tab_3D[sortedZ[k+1]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),points_tab_3D[sortedY[j+1]].get_pos_y(),points_tab_3D[sortedZ[k+1]].get_pos_z())));
				}

				else if (i==const_palier_num-1 && j==const_palier_num-1 && k== const_palier_num-1)
				{
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(1,1,1)));
				}

				else if(i==const_palier_num-1 && j==const_palier_num-1)
				{
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(1,1,points_tab_3D[sortedZ[k]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k+1]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(1,1,points_tab_3D[sortedZ[k+1]].get_pos_z())));
				}

				else if(k==const_palier_num-1 && j==const_palier_num-1)
				{
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),1,1)));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),1,1)));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z())));
				}

				else if(i==const_palier_num-1 && k==const_palier_num-1)
				{
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(1,points_tab_3D[sortedY[j]].get_pos_y(),1)));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(1,points_tab_3D[sortedY[j+1]].get_pos_y(),1)));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j+1]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z())));
				}

				else if(i==const_palier_num-1)
				{
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(1,points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z())));

					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j+1]].get_pos_y(),points_tab_3D[sortedZ[k+1]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(1,points_tab_3D[sortedY[j+1]].get_pos_y(),points_tab_3D[sortedZ[k+1]].get_pos_z())));

					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k+1]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(1,points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k+1]].get_pos_z())));

					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j+1]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(1,points_tab_3D[sortedY[j+1]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z())));
				}
				else if(j==const_palier_num-1 )
				{
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),1,points_tab_3D[sortedZ[k]].get_pos_z())));

					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k+1]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),1,points_tab_3D[sortedZ[k+1]].get_pos_z())));

					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k+1]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),1,points_tab_3D[sortedZ[k+1]].get_pos_z())));

					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z())));
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),1,points_tab_3D[sortedZ[k]].get_pos_z())));

				}
				else if(k==const_palier_num-1)
				{
					tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),1)));

                    tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),points_tab_3D[sortedY[j+1]].get_pos_y(),1)));
                    tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),points_tab_3D[sortedY[j+1]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z())));

                    tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),1)));
                    tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i+1]].get_pos_x(),points_tab_3D[sortedY[j]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z())));

                    tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j+1]].get_pos_y(),1)));
                    tmp=std::max(tmp,fabs((((double)sum))/((double)const_palier_num)-box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[j+1]].get_pos_y(),points_tab_3D[sortedZ[k]].get_pos_z())));

				}
				else
					std::cout<<"impossible"<<std::endl;
				//on prépare la maj de tmp_ligne et tmp_box
				tmp_ligne_maj[j]=memoire;
				tmp_box_maj[i][j]=sum;

				//comparer max
				#pragma omp flush(discrepancytmp)
				if(discrepancytmp<=tmp )
				{
					#pragma omp crititcal
					if(discrepancytmp<=tmp  )
					{
						discrepancytmp=tmp;
					}
				}
			}///fin j
		}///fin i
	}///fin k

	/////////// BOUCLES SUPPLEMENTAIRES POUR TESTER LES BOITES VIDES SUR LES AXES X, Y ET Z///////////////////////////////

	#pragma omp parallel for num_threads(nb_threads)
	for(int i=0;i<palier;i++)
	{
		double tmp=box_area_calcul(points_tab_3D[sortedX[i]].get_pos_x(),points_tab_3D[sortedY[0]].get_pos_y(),points_tab_3D[sortedZ[0]].get_pos_z());
		#pragma omp flush(discrepancytmp)
		if(discrepancytmp<=tmp)
		{
			#pragma omp crititcal
			if(discrepancytmp<=tmp)
			{
				discrepancytmp=tmp;
			}
		}
	}

	#pragma omp parallel for num_threads(nb_threads)
	for(int i=1;i<palier;i++)
	{
		double tmp=box_area_calcul(points_tab_3D[sortedX[0]].get_pos_x(),points_tab_3D[sortedY[i]].get_pos_y(),points_tab_3D[sortedZ[0]].get_pos_z());
		#pragma omp flush(discrepancytmp)
		if(discrepancytmp<=tmp)
		{
			#pragma omp crititcal
			if(discrepancytmp<=tmp)
			{
				discrepancytmp=tmp;
			}
		}
	}

	#pragma omp parallel for num_threads(nb_threads)
	for(int i=1;i<palier;i++)
	{
		double tmp=box_area_calcul(points_tab_3D[sortedX[0]].get_pos_x(),points_tab_3D[sortedY[0]].get_pos_y(),points_tab_3D[sortedZ[0]].get_pos_z());
		#pragma omp flush(discrepancytmp)
		if(discrepancytmp<=tmp)
		{
			#pragma omp crititcal
			if(discrepancytmp<=tmp)
			{
				discrepancytmp=tmp;
			}
		}
	}

    discrepancytmp=std::max (discrepancytmp,fabs(box_area_calcul(points_tab_3D[sortedX[0]].get_pos_x(),1,1)));
    discrepancytmp=std::max (discrepancytmp,fabs(box_area_calcul(1,1,points_tab_3D[sortedZ[0]].get_pos_z())));
    discrepancytmp=std::max (discrepancytmp,fabs(box_area_calcul(1,points_tab_3D[sortedY[0]].get_pos_y(),1)));


	return discrepancytmp;
}

///fonction export data sequence

void sequence_export_data(int i,double discrepancy){

	std::ofstream file(fn_output.c_str(),std::ios::app);
	if(!file)
		std::cerr << "Can't open ouptut" << std::endl;

    //file << pts_number << "  " << discrepancy << std::endl;
    file << discrepancy << std::endl;
    file.close();
}

///main

int main(int argc, char** argv){

	double discrepancy=0;

	import_options(argc,argv);

	srand(time(NULL));

	td=time(NULL);

	if(!sequence){

		if(dimension == 1){

			if(import_data_1D()){

				discrepancy=calcul_discrepancy_1D();

				std::cout<<" Discrepancy : "<<discrepancy<<std::endl;

				export_data_1D(discrepancy);

			}
		}

		else if(dimension == 2){

			if(import_data_2D()){


				discrepancy=calcul_discrepancy_2D();

				std::cout<<" Discrepancy : "<<discrepancy<<std::endl;

				export_data_2D(discrepancy);

			}
		}

		else if (dimension == 3){

			if(import_data_3D()){

				discrepancy=calcul_discrepancy_3D();

				std::cout<<" Discrepancy : "<<discrepancy<<std::endl;

				export_data_3D(discrepancy);

			}
		}

		else std::cout<<" Dimension not yet implemented."<<std::endl;

		tf=time(NULL);

		std::cout << " Fin de la recherche en " << difftime(tf, td) << " secondes."<<std::endl ;
	}

	else {

		if(dimension == 1){

			if(sequence_import_data_1D(1)){

				std::ofstream fichier(fn_output.c_str(),std::ios::out);
				if(fichier)
				{
					//fichier<<"Discrepance pour "<<pts_number<<" points."<<std::endl;
					fichier.close();
				}
				else
					std::cerr<<"Erreur à l'ouverture fichier output !"<<std::endl;

				Points_1D* sequence_points_tab_1D = new Points_1D[pts_number];

				#pragma omp parallel for num_threads(nb_threads)
				for(int i=0;i<pts_number;i++)
				{
					sequence_points_tab_1D[i]=points_tab_1D[i];
				}

				delete [] points_tab_1D;

				for(int palier=1;palier<=pts_number;palier++)
				{
					sortedX=new int[palier];
					points_tab_1D=new Points_1D[palier];

                    #pragma omp parallel for num_threads(nb_threads)
                    for(int j=0;j<palier;j++)
                    {
                        points_tab_1D[j].set_pos_x(sequence_points_tab_1D[j].get_pos_x());
                    }

                    #pragma omp parallel for num_threads(nb_threads)
                    for(int i=0; i<palier; i++) sortedX[i] = i;
                    std::sort(&sortedX[0],&sortedX[palier],sortXfunc);
					discrepancy= sequence_calcul_discrepancy_1D(palier);
					sequence_export_data(palier,discrepancy);
				}
			}
        }

		else if(dimension == 2){

			if(sequence_import_data_2D(1)){

				std::ofstream fichier(fn_output.c_str(),std::ios::out);
				if(fichier)
				{
					//fichier<<"Discrepance pour "<<pts_number<<" points."<<std::endl;
					fichier.close();
				}
				else
					std::cerr<<"Erreur à l'ouverture fichier output !"<<std::endl;

				Points_2D* sequence_points_tab_2D = new Points_2D[pts_number];

				#pragma omp parallel for num_threads(nb_threads)
				for(int i=0;i<pts_number;i++)
				{
					sequence_points_tab_2D[i]=points_tab_2D[i];
				}

				delete [] points_tab_2D;

				for(int palier=1;palier<=pts_number;palier++)
				{
					/////////////////////////////////////////////////////////////

					sortedX=new int[palier];
					sortedY=new int[palier];

					index_x_sortY=new int[palier];
                    index_y_sortX=new int[palier];

					tmp_ligne=new int[palier];
					tmp_ligne_maj=new int[palier];

					points_tab_2D=new Points_2D[palier];

					//////////////////////////////////////////////////////////////



                        #pragma omp parallel for num_threads(nb_threads)
						for(int j=0;j<palier;j++)//on ajoute le nv point ds le tab de points
						{
							points_tab_2D[j].set_pos_x(sequence_points_tab_2D[j].get_pos_x());
							points_tab_2D[j].set_pos_y(sequence_points_tab_2D[j].get_pos_y());
						}
                        #pragma omp parallel for num_threads(nb_threads)
                        for(int i=0; i<palier; i++) sortedX[i] = i;
                        std::sort(&sortedX[0],&sortedX[palier],sortXfunc);

                        #pragma omp parallel for num_threads(nb_threads)
                        for(int i=0; i<palier; i++) sortedY[i] = i;
                        std::sort(&sortedY[0],&sortedY[palier],sortYfunc);


					discrepancy= sequence_calcul_discrepancy_2D(palier);
					sequence_export_data(palier,discrepancy);
				}

			}
		}

		else if (dimension == 3){

			if(sequence_import_data_3D(1)){

				std::ofstream fichier(fn_output.c_str(),std::ios::out);
				if(fichier)
				{
					//fichier<<"Discrepance pour "<<pts_number<<" points."<<std::endl;
					fichier.close();
				}
				else
					std::cerr<<"Erreur à l'ouverture fichier output !"<<std::endl;

				Points_3D* sequence_points_tab_3D = new Points_3D[pts_number];

				#pragma omp parallel for num_threads(nb_threads)
				for(int i=0;i<pts_number;i++)
				{
					sequence_points_tab_3D[i]=points_tab_3D[i];
				}

				delete [] points_tab_3D;

				for(int palier=1;palier<=pts_number;palier++)
				{

                    sortedX=new int[palier];
                    sortedY=new int[palier];
                    sortedZ=new int[palier];

                    index_x_sortY=new int[palier];
                    index_x_sortZ=new int[palier];

                    index_y_sortX=new int[palier];
                    index_y_sortZ=new int[palier];

                    index_z_sortX=new int[palier];
                    index_z_sortY=new int[palier];

					tmp_ligne=new int[palier];
					tmp_ligne_maj=new int[palier];

					tmp_box=new int*[palier];
					tmp_box_maj=new int*[palier];

					#pragma omp parallel for num_threads(nb_threads)
					for(int i=0;i<palier;i++)
					{
						tmp_box[i]=new int [palier];
						tmp_box_maj[i]=new int [palier];

					}

					points_tab_3D=new Points_3D[palier];

					//////////////////////////////////////////////////////////////

                        #pragma omp parallel for num_threads(nb_threads)
						for(int j=0;j<palier;j++)
						{
							points_tab_3D[j].set_pos_x(sequence_points_tab_3D[j].get_pos_x());
							points_tab_3D[j].set_pos_y(sequence_points_tab_3D[j].get_pos_y());
							points_tab_3D[j].set_pos_z(sequence_points_tab_3D[j].get_pos_z());
						}


                        #pragma omp parallel for num_threads(nb_threads)
                        for(int i=0; i<palier; i++) sortedX[i] = i;
                        std::sort(&sortedX[0],&sortedX[palier],sortXfunc);

                        #pragma omp parallel for num_threads(nb_threads)
                        for(int i=0; i<palier; i++) sortedY[i] = i;
                        std::sort(&sortedY[0],&sortedY[palier],sortYfunc);

                        #pragma omp parallel for num_threads(nb_threads)
                        for(int i=0; i<palier; i++) sortedZ[i] = i;
                        std::sort(&sortedZ[0],&sortedZ[palier],sortZfunc);


					discrepancy= sequence_calcul_discrepancy_3D(palier);
					sequence_export_data(palier,discrepancy);

                    #pragma omp parallel for num_threads(nb_threads)
					for(int i=0;i<palier;i++)
					{
						delete[] tmp_box[i];
						delete[] tmp_box_maj[i];

					}
				}
			}
		}


	}

    std::cout<<" Traitement termine."<<std::endl<<std::endl;
}











